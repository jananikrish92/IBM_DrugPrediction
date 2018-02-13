import csv
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import _CalcMolWt
countp = 0
countn = 0
count = 0
with open('BindingDB_All.tsv') as tsvfile:
	reader = csv.DictReader(tsvfile, dialect='excel-tab')
	for row in reader:
		try:
			if(row['IC50 (nM)'] == ''):
				continue
			elif float(row['IC50 (nM)'])<100 and _CalcMolWt(Chem.MolFromSmiles(row['Ligand SMILES']))<1000:
				countp = countp+1
			elif float(row['IC50 (nM)'])>1000 and _CalcMolWt(Chem.MolFromSmiles(row['Ligand SMILES']))<1000:
				countn = countn+1
		except ValueError:
			# str = row['IC50 (nM)']
			str = row['IC50 (nM)']
			if str[0]=='>' and float(str[1:])>1000 and _CalcMolWt(Chem.MolFromSmiles(row['Ligand SMILES']))<1000:
				countn = countn+1
			elif str[0]=='<' and float(str[1:])<100 and _CalcMolWt(Chem.MolFromSmiles(row['Ligand SMILES']))<1000:
				countp = countp+1

		except TypeError:
			continue
print "no of positives",countp
print "no of negatives",countn