import csv

header = ['Pubchem CID','proteinId','smile','fastaseq','label']
def addHead(path='./train/train.csv',wpath='./train/trainData.csv'):
	with open(path) as file,open(wpath,'w') as wfile:
		next(file)
		reader = csv.reader(file,delimiter=",")
		writer = csv.writer(wfile,delimiter=",")
		writer.writerow(header)
		for row in reader:
			writer.writerow(row)


addHead('./dev/dev.csv','./dev/devData.csv')