#
# Example run:
# python plot_ami.py corresponding-top-picks-7pathways-6overlap-withcollapsed.txt path/to/output/
#

import sys
sys.path.append("..")
import colors

import os
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics.cluster import adjusted_mutual_info_score
#import glob


def cluster(df,numclusters):
	global model
	X = df.T.values
	link = 'average'
	affin = 'cosine'
	model = AgglomerativeClustering(distance_threshold=None, n_clusters=numclusters,affinity=affin,linkage=link)
	model = model.fit(X)
	mapping = dict(zip(df.columns,model.labels_))
	return model.labels_
	#lets reverse the mapping into a classification
	reversed_dict = dict()
	for k,v in mapping.items():
		try:
			reversed_dict[v].append(k)
		except KeyError:
			reversed_dict[v]=[k]
	return reversed_dict

def plotami(x,y,svname):
	plt.clf()
	plt.scatter(x,y)
	plt.grid()
	plt.xlabel('# clusters')
	plt.ylabel('AMI')
	plt.savefig(svname)
	print(svname.split('/')[-1])

def plot_combined_ami(values,outname):
	max_val = 0
	min_val = 1
	for k1 in values:

		for k2 in values[k1]:
			max_val = max(max_val,max(values[k1][k2][1]))
			min_val = min(min_val,min(values[k1][k2][1]))
			max_x = max(values[k1][k2][0])
			min_x = min(values[k1][k2][0])

	fig = plt.figure(figsize=(6,5))
	i = 1
	for metric in ['gcount','rho']:
		for gt in ['db','pathway']:
			ax = plt.subplot(2,2,i)
			if gt == 'db':
				gt_x = 9
			else:
				gt_x = 6
			ax.plot([gt_x,gt_x],[min_val,max_val],'-r',label=None)

			ax.plot(values[metric][gt][0],values[metric][gt][1],'o-',color=colors.COLORS['lightblue'],ms=5,label='Original')
			ax.plot(values['pred.'+metric][gt][0],values['pred.'+metric][gt][1],'s-',color=colors.COLORS['darkbrown'],ms=5,label='Transformed')
			if i==1:
				ax.set_title('Database Labels',fontsize=10,fontweight='bold')
				ax.set_ylabel('Graphlet Count AMI',fontsize=10,fontweight='bold')
			if i==2:
				ax.set_title('Pathway Labels',fontsize=10,fontweight='bold')
			if i==3:
				ax.set_ylabel('GHuST Coefficient AMI',fontsize=10,fontweight='bold')
				ax.set_xlabel('# Clusters')
			if i == 4:
				ax.set_xlabel('# Clusters')

			plt.legend(loc=1,fontsize=8)
			ax.set_ylim([min_val,max_val])
			ax.set_xlim([1,50])
			ax.set_xticks([min_x,gt_x,max_x])
			ax.set_xticklabels([min_x,gt_x,max_x])
			i+=1
	plt.tight_layout()
	plt.savefig(outname)
	print('Saved to %s' % (outname))


def aggregate_coeff(correspondence_file,path,ext):

	cp = pd.read_csv(correspondence_file,sep='\t',index_col=0)
	dbs = list(cp.columns.values)
	pathways = list(cp.index.values)

	df=None
	for i in pathways:
		for j in dbs:
			if (pd.notnull(cp.loc[i,j])):
				f = path+cp.loc[i,j][:-4]+'-network.txt.orca.'+ext
				df1=pd.read_csv(f,sep =' ', header=None)
				df1=pd.DataFrame(df1[1])
				df1=df1.rename(columns={1:j+'-'+i}).T
				df1['db']=j
				df1['pathway']=i
				if df is not None:
					df = df.append(df1)
				else:
					df = df1
	return df


def main(argv):

	correspondence_file = argv[1]
	path_out = argv[2]

	path_orig = r'../../graphlets/dbs/'
	path_reg = r'../../regression/predicted/dbs/'
	path_cp = r'../../networks/dbs/'

	correspondence_file=path_cp+correspondence_file
	values = {}
	for ext in list(['rho','gcount','pred.rho','pred.gcount']):
		values[ext] = {}
		for gt in list(['db','pathway']):

			#ext = 'pred.gcount'
			if ((ext == 'gcount') | (ext == 'rho')):
				df = aggregate_coeff(correspondence_file,path_orig,ext)
			elif ((ext == 'pred.gcount') | (ext == 'pred.rho')):
				df = aggregate_coeff(correspondence_file,path_reg,ext)

			if (gt == 'db'):
				ground_truth = df['db']
			elif (gt == 'pathway'):
				ground_truth = df['pathway']

			df = df[df.columns[0:len(df.columns)-2]]

			if (ext=='gcount'):
				for i in range(0,len(df.index)):
					s=df.iloc[i,0]
					if s > 0:
						df.iloc[i,0]=df.iloc[i,0]/s

					s=df.iloc[i,1:3].sum()
					if s > 0:
						df.iloc[i,1:3]=df.iloc[i,1:3]/s

					s=df.iloc[i,3:9].sum()
					if s > 0:
						df.iloc[i,3:9]=df.iloc[i,3:9]/s

					s=df.iloc[i,9:30].sum()
					if s > 0:
						df.iloc[i,9:30]=df.iloc[i,9:30]/s
			df=df.T

			mn_nm = 2
			mx_nm = len(df.columns)
			irange = list(range(mn_nm,mx_nm+1))

			amirange = [adjusted_mutual_info_score(cluster(df,i),ground_truth) for i in irange]
			amiplnm = path_out+'by-'+gt+'-ami-'+ext+'.pdf'
			values[ext][gt] = (irange,amirange)
			plotami(irange,amirange,amiplnm)
	plot_combined_ami(values,'%s/combined_ami.pdf' % path_out)

if __name__ == "__main__":
	main(sys.argv)
