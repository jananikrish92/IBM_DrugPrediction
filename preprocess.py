import csv
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import _CalcMolWt
countp = 0
countn = 0
count = 0
with open('../BindingDB_All.tsv') as tsvfile:
	reader = csv.DictReader(tsvfile, dialect='excel-tab')
	for row in reader:
		count = count +1
		# import pdb
		# pdb.set_trace()
		try:
			if(row['IC50 (nM)'] == '' or row['Ligand SMILES']=='' or row['Ligand SMILES'] is None):
				continue
			else:
				if float(row['IC50 (nM)'])<100 and _CalcMolWt(Chem.MolFromSmiles(row['Ligand SMILES']))<1000:
					countp = countp+1
				elif float(row['IC50 (nM)'])>1000 and _CalcMolWt(Chem.MolFromSmiles(row['Ligand SMILES']))<1000:
					countn = countn+1
		except ValueError:
			try:
				str = row['IC50 (nM)']
				if str[0]=='>' and float(str[1:])>1000 and _CalcMolWt(Chem.MolFromSmiles(row['Ligand SMILES']))<1000:
					countn = countn+1
				elif str[0]=='<' and float(str[1:])<100 and _CalcMolWt(Chem.MolFromSmiles(row['Ligand SMILES']))<1000:
					countp = countp+1
			except TypeError:
				continue
		except TypeError:
			continue
		except:
			continue
print "no of positives",countp
print "no of negatives",countn