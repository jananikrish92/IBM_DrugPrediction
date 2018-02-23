import csv
import numpy as np
import theano.tensor as T
import theano
from sklearn.cross_validation import train_test_split
def load_labeled_data():
    
    X=np.zeros([21,14564])
    Y=np.zeros([21,])

    with open('newData-new.csv') as tsvfile:
        reader=csv.reader(tsvfile)
        count=0
        for row in reader:
            count+=1  
           # print "bhanu", count
            if(count==1):
                continue
            X[count-2]=np.array(row[2][1:len(row[2])-1].split(',')+row[3][1:len(row[3])-1].split(',')+row[4][1:len(row[4])-1].split(',')+row[5][1:len(row[5])-1].split(',')).astype(np.float)
            Y[count-2]=np.array(row[6]).astype(np.int)
    
    tv_set_x, test_set_x, tv_set_y, test_set_y = train_test_split(X, Y, test_size=0.25, random_state=123)
    train_set_x, valid_set_x, train_set_y, valid_set_y = train_test_split(tv_set_x, tv_set_y, test_size=0.2, random_state=123)
    
    DTI = (
        (
            train_set_x,
            train_set_y
        ),
        (
            valid_set_x,
            valid_set_y
        ),
        (
            test_set_x,
            test_set_y
        )
    )
    
    train_set, valid_set, test_set = DTI
    def shared_dataset(data_xy, borrow=True):
        data_x, data_y = data_xy
        shared_x = theano.shared(np.asarray(data_x,
                                               dtype=theano.config.floatX),
                                 borrow=borrow)
        shared_y = theano.shared(np.asarray(data_y,
                                               dtype=theano.config.floatX),
                                 borrow=borrow)
        # When storing data on the GPU it has to be stored as floats
        # therefore we will store the labels as ``floatX`` as well
        # (``shared_y`` does exactly that). But during our computations
        # we need them as ints (we use labels as index, and if they are
        # floats it doesn't make sense) therefore instead of returning
        # ``shared_y`` we will have to cast it to int. This little hack
        # lets ous get around this issue
        return shared_x, T.cast(shared_y, 'int32')
    test_set_x, test_set_y = shared_dataset(test_set)
    valid_set_x, valid_set_y = shared_dataset(valid_set)
    train_set_x, train_set_y = shared_dataset(train_set)
    
    rval = [(train_set_x, train_set_y), (valid_set_x, valid_set_y), (test_set_x, test_set_y)]
    return rval
         
        
