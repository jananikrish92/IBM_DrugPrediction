import csv
from rdkit import Chem, DataStructs
from rdkit.Chem.rdMolDescriptors import _CalcMolWt
from rdkit.Chem import AllChem
import numpy


def convToArr(fp):
	arr = numpy.zeros((1,))
  	DataStructs.ConvertToNumpyArray(fp, arr)
  	return arr

count = 0
countp  = 0
countn = 0
with open('../BindingDB_All.tsv') as tsvfile,open('../newData.csv','w') as newfile:
	reader = csv.DictReader(tsvfile,dialect='excel-tab')
	writer = csv.writer(newfile,delimiter=',')
	writer.writerow(["PubChem CID", "fingerprintR1", "fingerprintR2","fingerprintR3"])
	for row in reader:
		pid = row['PubChem CID']
		if pid!='':
			try:
				if(row['IC50 (nM)']=='' or row['Ligand SMILES']=='' or row['Ligand SMILES'] is None):
					continue
				else:
					md = Chem.MolFromSmiles(row['Ligand SMILES'])
					if float(row['IC50 (nM)'])<100 and _CalcMolWt(md)<1000 and md is not None:
						fp1 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,1))
						fp2 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,2))
						fp3 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,3))
						out = pid,fp1.tolist(),fp2.tolist(),fp3.tolist()
						writer.writerow(out)
						countp += 1
					elif float(row['IC50 (nM)'])>1000 and _CalcMolWt(md)<100 and md is not None:
						fp1 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,1))
						fp2 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,2))
						fp3 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,3))
						out = pid,fp1.tolist(),fp2.tolist(),fp3.tolist()
						writer.writerow(out)						
						countn += 1
			except ValueError:
				try:
					str = row['IC50 (nM)']
					md = Chem.MolFromSmiles(row['Ligand SMILES'])
					if str[0]=='>' and float(str[1:])>1000 and _CalcMolWt(md)<1000 and md is not None:
						countn += 1
						fp1 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,1))
						fp2 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,2))
						fp3 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,3))
						out = pid,fp1.tolist(),fp2.tolist(),fp3.tolist()
						writer.writerow(out)
					elif str[0]=='<' and float(str[1:])<100 and _CalcMolWt(md)<1000 and md is not None:
						countp += 1
						fp1 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,1))
						fp2 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,2))
						fp3 = convToArr(AllChem.GetMorganFingerprintAsBitVect(md,3))
						out = pid,fp1.tolist(),fp2.tolist(),fp3.tolist()
						writer.writerow(out)
				except TypeError:
					continue
			except TypeError:
				continue
			except:
				continue

print countn,countp
