import csv
from rdkit import Chem, DataStructs
from rdkit.Chem.rdMolDescriptors import _CalcMolWt
from rdkit.Chem import AllChem
import numpy
from read_fasta import calPscIdx


def convToArr(fp):
	arr = numpy.zeros((1,))
  	DataStructs.ConvertToNumpyArray(fp, arr)
  	return arr

count = 0
countp  = 0
countn = 0
pubids = list()
psc_array,prot_id2idx = calPscIdx()

print "collected pds"
# def featureGen(ic50_range):
with open('../BindingDB_All.tsv') as tsvfile,open('../newData.csv','w') as newfile:
	next(tsvfile)
	reader = csv.reader(tsvfile,delimiter='\t')
	writer = csv.writer(newfile,delimiter=',')
	writer.writerow(["PubChem CID","protienId", "fingerprintR1", "fingerprintR2","fingerprintR3","protein descriptor","label"])
	for row in reader:
		# print len(row)
		print count
		count += 1
		pid = row[28]
		if pid!='' and pid not in pubids and row[41]!='' and row[37]!='':
			pubids.append(pid)
			# import pdb
			# pdb.set_trace()
			print "entered"
			proId = row[41].split(',')[0]
			if proId not in prot_id2idx:
				print proId,pid,row[37]
			pd = psc_array[prot_id2idx[proId]]
			try:
				if(row[9]=='' or row[1]=='' or row[1] is None):
					continue
				else:
					md = Chem.MolFromSmiles(row[1])
					if float(row[9])<100 and _CalcMolWt(md)<1000 and md is not None:
						fp1 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,1))
						fp2 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,2))
						fp3 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,3))
						out = pid,proId,fp1.tolist(),fp2.tolist(),fp3.tolist(),pd,1
						writer.writerow(out)
						print "writing done"
						countp += 1
					elif float(row[9])>1000 and _CalcMolWt(md)<100 and md is not None:					
						fp1 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,1))
						fp2 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,2))
						fp3 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,3))
						out = pid,proId,fp1.tolist(),fp2.tolist(),fp3.tolist(),pd,0
						writer.writerow(out)	
						print "writing done"					
						countn += 1
			except ValueError:
				try:
					str = row[9]
					md = Chem.MolFromSmiles(row[1])
					if str[0]=='>' and float(str[1:])>1000 and _CalcMolWt(md)<1000 and md is not None:
						countn += 1
						fp1 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,1))
						fp2 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,2))
						fp3 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,3))
						out = pid,proId,fp1.tolist(),fp2.tolist(),fp3.tolist(),pd,0
						writer.writerow(out)
						print "writing done"
					elif str[0]=='<' and float(str[1:])<100 and _CalcMolWt(md)<1000 and md is not None:
						countp += 1
						fp1 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,1))
						fp2 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,2))
						fp3 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,3))
						out = pid,proId,fp1.tolist(),fp2.tolist(),fp3.tolist(),pd,1
						writer.writerow(out)
						print "writing done"
				except TypeError:
					continue
			except TypeError:
				continue
			except:
				continue


print countn,countp

# featureGen([100,1000])
