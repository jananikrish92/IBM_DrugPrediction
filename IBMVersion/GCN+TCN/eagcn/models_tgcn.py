import torch.nn as nn
import torch.nn.functional as F
from layers import *
from utils import *
import torch
from torch.autograd import Variable
from torch.utils.data import Dataset

from models import *
from tcn.tcn import TemporalConvNet
from time import gmtime, strftime

from sklearn import metrics
from sklearn.utils import shuffle, resample
from sklearn.model_selection import train_test_split, KFold

use_cuda = torch.cuda.is_available()
FloatTensor = torch.cuda.FloatTensor if use_cuda else torch.FloatTensor
LongTensor = torch.cuda.LongTensor if use_cuda else torch.LongTensor
IntTensor = torch.cuda.IntTensor if use_cuda else torch.IntTensor
DoubleTensor = torch.cuda.DoubleTensor if use_cuda else torch.DoubleTensor

class dti_pred_net(nn.Module):
    def __init__(self, gcn_params, tcn_params, n_den1=165, n_den2=64, n_den3=32):

        super(dti_pred_net, self).__init__()

        self.gcn = Concate_GCN(n_bfeat=gcn_params.n_bfeat, n_afeat=25,
                            n_sgc1_1=gcn_params.n_sgc1_1, n_sgc1_2=gcn_params.n_sgc1_2, n_sgc1_3=gcn_params.n_sgc1_3, 
                            n_sgc1_4=gcn_params.n_sgc1_4, n_sgc1_5=gcn_params.n_sgc1_5,
                            n_sgc2_1=gcn_params.n_sgc2_1, n_sgc2_2=gcn_params.n_sgc2_2, n_sgc2_3=gcn_params.n_sgc2_3, 
                            n_sgc2_4=gcn_params.n_sgc2_4, n_sgc2_5=gcn_params.n_sgc2_5,
                            n_den1=gcn_params.n_den1, n_den2=gcn_params.n_den2,
                            nclass=gcn_params.nclass, dropout=gcn_params.dropout, use_att = True, molfp_mode = 'sum')
        self.encoder = nn.Embedding(27, tcn_params.input_size)
        self.tcn = TemporalConvNet(tcn_params.input_size, tcn_params.num_channels, kernel_size=tcn_params.kernel_size, dropout=tcn_params.dropout)
        self.den1 = nn.Linear(n_den1, n_den2)
        self.den2 = nn.Linear(n_den2, n_den3)
        self.den3 = nn.Linear(n_den3, 1)
        self.molfp_bn = nn.BatchNorm1d(140).cuda() if use_cuda else nn.BatchNorm1d(140)
        
        self.attn_mat = nn.Parameter(torch.randn(140,25))

    def forward(self, adjs, afms, bfts, OrderAtt, AromAtt, ConjAtt, RingAtt, xp):

        fp = self.gcn(adjs, afms, bfts, OrderAtt, AromAtt, ConjAtt, RingAtt)
        xp = self.encoder(xp)
        pd = self.tcn(xp.transpose(1,2))
        # import pdb
        # pdb.set_trace()
        #pd = pd[:,:,-1]
        fp_temp = torch.mm(fp.view(-1,140), self.attn_mat).view(fp.size(0),fp.size(1),-1)
        algnt_mats = torch.bmm(fp_temp, pd)

        drug_attn_weights,_ = torch.max(algnt_mats, dim=2, keepdim=True)
        prot_attn_weights,_ = torch.max(algnt_mats, dim=1, keepdim=True)
        prot_attn_weights = prot_attn_weights.permute(0,2,1)

        drug_attn_weights = F.softmax(drug_attn_weights, dim=1)
        prot_attn_weights = F.softmax(prot_attn_weights, dim=1)

        fp = self.molfp_bn(torch.sum(fp*drug_attn_weights,dim=1))
        pd = torch.sum(pd.permute(0,2,1)*prot_attn_weights,dim=1)



        x1 = torch.cat((fp, pd), dim=1)
        x1 = F.relu(self.den1(x1))
        x1 = F.relu(self.den2(x1))
        y = self.den3(x1)
        

        return y

class tcn_params():
    def __init__(self, input_size, num_channels, kernel_size, dropout):
        self.input_size = input_size
        self.num_channels = num_channels
        self.kernel_size = kernel_size
        self.dropout = dropout

class gcn_params():
    def __init__(self, n_bfeat, n_sgc1_1, n_sgc1_2, n_sgc1_3, n_sgc1_4, n_sgc1_5, 
                n_sgc2_1, n_sgc2_2, n_sgc2_3, n_sgc2_4, n_sgc2_5, nclass=1, dropout=0.3):
        self.n_bfeat = n_bfeat
        self.n_sgc1_1 = n_sgc1_1
        self.n_sgc1_2 = n_sgc1_2
        self.n_sgc1_3 = n_sgc1_3
        self.n_sgc1_4 = n_sgc1_4
        self.n_sgc1_5 = n_sgc1_5
        self.n_sgc2_1 = n_sgc2_1
        self.n_sgc2_2 = n_sgc2_2
        self.n_sgc2_3 = n_sgc2_3
        self.n_sgc2_4 = n_sgc2_4
        self.n_sgc2_5 = n_sgc2_5
        self.dropout = dropout
        self.n_den1 = n_den1 
        self.n_den2 = n_den2
        self.nclass = nclass

def test_model(loader, model, tasks):
    """
    Help function that tests the model's performance on a dataset
    @param: loader - data loader for the dataset to test against
    """
    true_value = []
    all_out = []
    model.eval()
    out_value_dic = {}
    true_value_dic = {}
    for adj, afm, btf, orderAtt, aromAtt, conjAtt, ringAtt,xp, labels in loader:
        adj_batch, afm_batch, btf_batch, label_batch = Variable(adj), Variable(afm), Variable(btf), Variable(labels)
        orderAtt_batch, aromAtt_batch, conjAtt_batch, ringAtt_batch = Variable(orderAtt), Variable(aromAtt), Variable(
            conjAtt), Variable(ringAtt)
        xp_batch = Variable(xp).long()
        outputs = model(adj_batch, afm_batch, btf_batch, orderAtt_batch, aromAtt_batch, conjAtt_batch, ringAtt_batch,xp_batch)
        probs = F.sigmoid(outputs)

        if use_cuda:
            out_list = probs.cpu().data.view(-1).numpy().tolist()
            all_out.extend(out_list)
            label_list = labels.cpu().numpy().tolist()
            true_value.extend([item for sublist in label_list for item in sublist])
            out_sep_list = probs.cpu().data.view(-1, len(tasks)).numpy().tolist()
        else:
            out_list = probs.data.view(-1).numpy().tolist()
            all_out.extend(out_list)
            label_list = labels.numpy().tolist()
            true_value.extend([item for sublist in label_list for item in sublist])
            out_sep_list = probs.data.view(-1, len(tasks)).numpy().tolist()

        for i in range(0, len(out_sep_list)):
            for j in list(range(0, len(tasks))):
                if label_list[i][j] == -1:
                    #print('Ignore {},{} case: nan'.format(i,j))
                    continue
                if j not in true_value_dic.keys():
                    out_value_dic[j] = [out_sep_list[i][j]]
                    # print(out_value_dic[j])
                    true_value_dic[j] = [int(label_list[i][j])]
                else:
                    out_value_dic[j].extend([out_sep_list[i][j]])
                    # print(out_value_dic[j])
                    true_value_dic[j].extend([int(label_list[i][j])])
    model.train()

    aucs = []
    for key in list(range(0, len(tasks))):
        fpr, tpr, threshold = metrics.roc_curve(true_value_dic[key], out_value_dic[key], pos_label=1)
        auc = metrics.auc(fpr, tpr)
        if math.isnan(auc):
            print('the {}th label has no postive samples, max value {}'.format(key, max(true_value_dic[key])))
        aucs.append(auc)
    return (aucs, sum(aucs)/len(aucs))


def train(tasks, EAGCN_structure, n_den1, n_den2, file_name):


    # x_all,xp_all,y_all, target, sizes = load_data(dataset)
    # max_size = max(sizes)
    # x_all,xp_all, y_all = data_filter(x_all,xp_all, y_all, target, sizes, tasks)
    # x_all, xp_all,y_all = shuffle(x_all,xp_all,y_all, random_state=random_state)

    #IBM data
    x_train,xp_train,y_train, target, sizes = load_data(dataset,path="../data/train/trainData.csv")
   
    max_size = max(sizes)
    #x_train,xp_train,y_train = data_filter(x_train,xp_train, y_train, target, sizes, tasks)
    x_train,xp_train,y_train = shuffle(x_train,xp_train,y_train, random_state=random_state)

    x_test,xp_test,y_test,test_target, sizes_test = load_data(dataset,path="../data/test/testData.csv")
    # max_size = max(sizes)
    #x_test,xp_test,y_test = data_filter(x_test,xp_test,y_test, test_target, sizes_test, tasks)
    x_test,xp_test,y_test = shuffle(x_test,xp_test,y_test, random_state=random_state)



    # X, x_test,Xp,xp_test, y, y_test = train_test_split(x_all,xp_all,y_all, test_size=0.1, random_state=random_state)
    # print(x_all)
    # del x_all,xp_all,y_all
    test_loader = construct_loader(x_test,xp_test,y_test, test_target, batch_size)
    del x_test,xp_test, y_test

    n_bfeat = x_train[0][2].shape[0]

    tcn_param = tcn_params(10, [25]*8, 7, 0.2)
    gcn_param = gcn_params(n_bfeat, 30, 10, 10, 10, 10, 60, 20, 20, 20, 20)

    model = dti_pred_net(gcn_param, tcn_param)

    if use_cuda:
        # lgr.info("Using the GPU")
        model.cuda()

    model.apply(weights_init)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)

    validation_acc_history = []
    stop_training = False
    BCE_weight = set_weight(y_train)

    X = np.array(x_train)
    y = np.array(y_train)
    #x_train, x_val,xp_train,xp_val, y_train, y_val = train_test_split(x_train,xp_train, y_train, test_size=0.1, random_state=random_state)
    #x_val,xp_val,y_val,test_target, sizes_test = load_data(dataset,path="../data/test/testData.csv")
    train_loader = construct_loader(x_train,xp_train, y_train, target, batch_size)
    #validation_loader = construct_loader(x_val, xp_val,y_val, target, batch_size)
    len_train = len(x_train)
    del x_train,xp_train, y_train#, x_val, xp_val,y_val
    for epoch in range(num_epochs):
        print epoch
        for i, (adj, afm, btf, orderAtt, aromAtt, conjAtt, ringAtt,xp, labels) in enumerate(train_loader):
            adj_batch, afm_batch, btf_batch, label_batch = Variable(adj), Variable(afm), Variable(btf), Variable(labels)
            orderAtt_batch, aromAtt_batch, conjAtt_batch, ringAtt_batch = Variable(orderAtt), Variable(
                aromAtt), Variable(
                conjAtt), Variable(ringAtt)
            xp = Variable(xp).long()
            # xp = xp.view(xp.shape[0],-1,xp.shape[1])
            # print xp

            optimizer.zero_grad()
            outputs = model(adj_batch, afm_batch, btf_batch, orderAtt_batch, aromAtt_batch, conjAtt_batch,
                            ringAtt_batch, xp)
            weights = Variable(weight_tensor(BCE_weight, labels=label_batch))
            non_nan_num = Variable(FloatTensor([(labels == 1).sum() + (labels == 0).sum()]))
            loss = F.binary_cross_entropy_with_logits(outputs.view(-1), \
                                                      label_batch.float().view(-1), \
                                                      weight=None, size_average=False) / non_nan_num
            loss.backward()
            optimizer.step()


        # report performance
        if True:
            train_acc_sep, train_acc_tot = test_model(train_loader, model, tasks)
            val_acc_sep, val_acc_tot = (1,1)#test_model(validation_loader, model, tasks)
            print(
                'Epoch: [{}/{}], '
                'Step: [{}/{}], '
                'Loss: {}, \n'
                'Train AUC seperate: {}, \n'
                'Train AUC total: {}, \n'
                'Validation AUC seperate: {}, \n'
                'Validation AUC total: {} \n'.format(
                    epoch + 1, num_epochs, i + 1,
                    math.ceil(len_train / batch_size), loss.data[0], \
                    train_acc_sep, train_acc_tot, val_acc_sep,
                    val_acc_tot))
            if write_file:
                with open(file_name, 'a') as fp:
                    fp.write(
                        'Epoch: [{}/{}], '
                        'Step: [{}/{}], '
                        'Loss: {}, \n'
                        'Train AUC seperate: {}, \n'
                        'Train AUC total: {}, \n'
                        'Validation AUC seperate: {}, \n'
                        'Validation AUC total: {} \n'.format(
                            epoch + 1, num_epochs, i + 1,
                            math.ceil(len_train / batch_size),
                            loss.data[0], \
                            train_acc_sep, train_acc_tot, val_acc_sep,
                            val_acc_tot))
            validation_acc_history.append(val_acc_tot)
            # check if we need to earily stop the model
            # stop_training = earily_stop(validation_acc_history, tasks, early_stop_step_single,
            #                             early_stop_step_multi, early_stop_required_progress) and (train_acc_tot > 0.99)
            # if stop_training:  # early stopping
            #     print("{}th epoch: earily stop triggered".format(epoch))
            #     if write_file:
            #         with open(file_name, 'a') as fp:
            #             fp.write("{}th epoch: earily stop triggered".format(epoch))
            #     break

        # because of the the nested loop
        if stop_training:
            break

    test_auc_sep, test_auc_tot = test_model(test_loader, model, tasks)
    torch.save(model.state_dict(), '{}.pkl'.format(file_name))
    torch.save(model, '{}.pt'.format(file_name))

    print('AUC of the model on the test set for single task: {}\n'
          'AUC of the model on the test set for all tasks: {}'.format(test_auc_sep, test_auc_tot))
    if write_file:
        with open(file_name, 'a') as fp:
            fp.write('AUC of the model on the test set for single task: {}\n'
                     'AUC of the model on the test set for all tasks: {}'.format(test_auc_sep, test_auc_tot))

    return(test_auc_tot)



dataset = 'bindingdb'
EAGCN_structure = 'concate' #  'concate', 'weighted_ave'
write_file = True
n_den1, n_den2= 64, 32
all_tasks = ['label']
batch_size = 64
weight_decay = 0.00001  # L-2 Norm
random_state = 1
num_epochs = 50
learning_rate = 0.0005
tasks = all_tasks # [task]
# print(' learning_rate: {},\n batch_size: {}, \n '
#           'tasks: {},\n random_state: {}, \n EAGCN_structure: {}\n'.format(
#         learning_rate, batch_size, tasks, random_state, EAGCN_structure))
# print('n_sgc1_1, n_sgc1_2, n_sgc1_3, n_sgc1_4, n_sgc1_5, '
#                  'n_sgc2_1, n_sgc2_2, n_sgc2_3, n_sgc2_4, n_sgc2_5: '
#                  '{}, {} {}, {}, {}, / {}, {}, {}, {}, {}'.format(n_sgc1_1, n_sgc1_2, n_sgc1_3, n_sgc1_4, n_sgc1_5,
#                                                                   n_sgc2_1, n_sgc2_2, n_sgc2_3, n_sgc2_4, n_sgc2_5))
# print('n_den1, nden2: {}, {}'.format(n_den1, n_den2))
experiment_date = strftime("%b_%d_%H_%M", gmtime()) +'N'

if use_cuda:
    position = 'server'
else:
    position = 'local'
if len(tasks) == 1:
    directory = '../experiment_result/{}/{}/{}/'.format(position, dataset, tasks)
    if not os.path.exists(directory):
        os.makedirs(directory)
    file_name = '{}{}'.format(directory, experiment_date)
else:
    directory = "../experiment_result/{}/{}/['all_tasks']/".format(position, dataset)
    if not os.path.exists(directory):
        os.makedirs(directory)
    file_name = '{}{}'.format(directory, experiment_date)


if write_file:
    # import pdb
    # pdb.set_trace()
    with open(file_name, 'w') as fp:
        fp.write(' learning_rate: {},\n batch_size: {}, \n '
                     'tasks: {},\n random_state: {} \n,'
                     ' EAGCN_structure: {}\n'.format(learning_rate, batch_size,
                                           tasks, random_state, EAGCN_structure))
        # fp.write('early_stop_step_single: {}, early_stop_step_multi: {}, \n'
        #              'early_stop_required_progress: {},\n early_stop_diff: {}, \n'
        #              'weight_decay: {}\n'.format(early_stop_step_single, early_stop_step_multi,
        #                                                      early_stop_required_progress, early_stop_diff,
        #                                                      weight_decay))
        # fp.write('n_sgc1_1, n_sgc1_2, n_sgc1_3, n_sgc1_4, n_sgc1_5, '
        #          'n_sgc2_1, n_sgc2_2, n_sgc2_3, n_sgc2_4, n_sgc2_5: '
        #          '{}, {} {}, {}, {}, / {}, {}, {}, {}, {}'.format(n_sgc1_1, n_sgc1_2, n_sgc1_3, n_sgc1_4, n_sgc1_5,
        #                                                           n_sgc2_1, n_sgc2_2, n_sgc2_3, n_sgc2_4, n_sgc2_5))


result = train(tasks, EAGCN_structure, n_den1, n_den2,file_name)



