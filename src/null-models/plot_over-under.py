# to run:
# python plot_over-under.py gcount
# python plot_over-under.py rho
# assuming all graphlet profiles (for real pathways and random-rewiring model) exist 
# in the appropriate directories

import sys
sys.path.append("..")
import colors

import pandas as pd 
import networkx as nx
import matplotlib.pyplot as plt
#import seaborn as sns
import glob
import os


def plot_barchart(overrep,j,path_out,ext):
	if (ext == 'rho'):
		xl = 'GHuST coefficient'
	else:
		xl = 'Graphlet'

	plt.clf()
	#plt.axes([0, 0, 1.5, 1])
	plt.bar(overrep['graphlet'],overrep['over'],label='over-represented',color=colors.COLORS['darkblue'])
	plt.bar(overrep['graphlet'],overrep['under'],bottom=overrep['over'],label='under-represented',color=colors.COLORS['lightblue'])
	plt.bar(overrep['graphlet'],overrep['not-sig'],bottom=overrep['under']+overrep['over'],label='not-significant',color=colors.COLORS['gray'])
	plt.legend(fontsize=10)
	#plt.legend(fontsize=20,bbox_to_anchor=(1.01, 1), loc='upper left')
	plt.xlabel(xl,fontsize=10)
	plt.ylabel('Fraction of pathways',fontsize=10)
	plt.xticks(fontsize=10,rotation=60)
	plt.yticks(fontsize=10)
	plt.title(j,fontsize=10)
	plt.savefig(path_out+ext+'-over-underrep-'+j+'.pdf',format='pdf',bbox_inches='tight')
	print('--> ',j)



def main(argv):
	ext = argv[1]
	path_orig = r'../../graphlets/dbs/'
	path_rand = r'../../graphlets/null-models/random-rewiring/'
	path_out = r'../../out/'


	if (ext == 'rho'):
		n_char = 20
	else:
		n_char = 23

	d = glob.glob(path_rand+'/*/')
	dbs = [x.split('/')[-2] for x in d]

	for j in dbs:
	    path = path_orig
	    path = path+j+'/'
	    filenames = glob.glob(path + '*-network.txt.orca.'+ext)
	    #print(j, len(filenames))
	    profile=pd.DataFrame()
	    u=0
	    for f in filenames:
	        if os.path.isfile(f.replace(path_orig,path_rand)[:-n_char]+'99randomnetwork-network.txt.orca.'+ext):

	            df=pd.read_table(f,sep=' ',header=None)
	            df=df.rename(columns={0:0,1:"pathway"})

	            for i in range(0,100): #1,101 for 100 random networks
	                df1=pd.read_table(f.replace(path_orig,path_rand)[:-n_char]+str(i)+'randomnetwork-network.txt.orca.'+ext,sep=' ',header=None)
	                df1=df1.rename(columns={0:0,1:i+1})
	                df=pd.merge(df,df1,on=0)
	            
	            df=df.rename(columns={0:"graphlet"})

	            df['mean']=df.iloc[:,2:102].mean(axis=1) #2:102 for 100 networks (last number is not inclusive)
	            df['sd']=df.iloc[:,2:102].std(axis=1)
	            df['z']=(df['pathway']-df['mean'])/df['sd']
	            df['overrep']=0

	            for i in df.index:

	                if df.at[i, 'z'] <= 2.0 and df.at[i, 'z'] >= -2.0:
	                    df.at[i, 'overrep']=0
	                elif df.at[i,'z'] < -2.0:
	                    df.at[i,'overrep'] = -1
	                elif df.at[i,'z'] > 2.0:
	                    df.at[i,'overrep'] = 1
	            profile.insert(u, f.replace(path,'')[:-24], df['overrep'], True)
	            u=u+1
	    overrep=pd.DataFrame(columns=["graphlet","over","under","not-sig"])
	    for u in range(0,len(profile)):
	        if (ext == 'rho'):
	            overrep.loc[u,"graphlet"]=str(u+1)
	        else:
	            overrep.loc[u,"graphlet"]=str(u) # ='G'+str(u) -add prefix G
	        overrep.loc[u,"over"]=(profile.T[u].values==1).sum()
	        overrep.loc[u,"under"]=(profile.T[u].values==-1).sum()
	        overrep.loc[u,"not-sig"]=(profile.T[u].values==0).sum()

	    number_of_pathways=overrep.loc[0][1:4].sum()
	    overrep.iloc[:,1:4]=overrep.iloc[:,1:4].transform(lambda x: x / number_of_pathways)
	    plot_barchart(overrep,j,path_out,ext)

if __name__ == "__main__":
    main(sys.argv)
