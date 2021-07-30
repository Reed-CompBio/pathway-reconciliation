# PCA plot for original rho and transformed rho coefficients
# for corresponding pathways listed in the corresponding pathways file
# to run
# original:
# python plot_pca.py rho
# transformed:
# python plot_pca.py pred.rho
#

import os
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def aggregate_coeff(path_cp,file_cp,path,ext):
    
    cp = pd.read_csv(path_cp+file_cp,sep='\t',index_col=0)
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


def make_plot(path_cp,file_cp,df,ext,path_out):

	cp = pd.read_csv(path_cp+file_cp ,sep='\t',index_col=0)
	pathways = list(cp.index.values)

	Y = df[['db','pathway']]
	#ground_truth = df['pathway']

	df = df[df.columns[0:len(df.columns)-2]]

	if (ext=='gcount'):

	    for i in range(0,len(df.index)):
	        df.iloc[i,0]=1

	        s=df.iloc[i,1:3].sum()
	        if s > 0:
	            df.iloc[i,1:3]=df.iloc[i,1:3]/s

	        s=df.iloc[i,3:9].sum()
	        if s > 0:
	            df.iloc[i,3:9]=df.iloc[i,3:9]/s

	        s=df.iloc[i,9:30].sum()
	        if s > 0:
	            df.iloc[i,9:30]=df.iloc[i,9:30]/s

	Y.index=range(0,len(Y.index))

	x = StandardScaler().fit_transform(df)

	pca = PCA(n_components=2)
	principalComponents = pca.fit_transform(x)
	principalDf = pd.DataFrame(data = principalComponents
	             , columns = ['principal component 1', 'principal component 2'])

	pw_dict = {pathways[i]:i for i in range(len(pathways))}
	markers = {'netpath':'s','kegg_expanded':'P','pid':'d','panther':'v', 'inoh':'*', 'pathbank':'^', 
           'signor_expanded':'>','signor_collapsed':'X','kegg_collapsed':'p'}
	this_cmap = plt.get_cmap('nipy_spectral',len(pathways))

	#plt.scatter(principalDf['principal component 1'],principalDf['principal component 2'],marker=[markers[i] for i in Y['DB']],alpha=0.5)
	for i in range(0,len(x)):
	    plt.scatter(principalDf.loc[i,'principal component 1'],principalDf.loc[i,'principal component 2'],c=[this_cmap(pw_dict[Y.loc[i,'pathway']])],marker=markers[Y.loc[i,'db']],alpha=0.5)
	    #plt.annotate(df.index[i], (principalDf['principal component 1'][i], principalDf['principal component 2'][i]),fontsize=4)

	plt.xlabel("PC-1")
	plt.ylabel("PC-2")
	xmin,xmax = plt.xlim()
	ymin,ymax = plt.ylim()

	for m in markers:
	        plt.scatter(xmin-.5,ymin-.5,c='w',edgecolor='k',lw=0.5,s=50,alpha=1,marker=markers[m],label=m)
	for i in range(len(pathways)):
	        plt.scatter(xmin-.5,ymin-.5,c=[this_cmap(i)],marker='.',s=50,alpha=0.5,label=pathways[i])
	plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=3)

	plt.xlim(xmin,xmax)
	plt.ylim(ymin,ymax)
	plt.savefig(path_out+'PCA-'+ext+'.pdf',format='pdf',bbox_inches='tight')


def main(argv):
	ext = argv[1] #rho, pred.rho, gcount, pred.gcount
	path_orig = r'../../graphlets/dbs/'
	path_reg = r'../../regression/predicted/dbs/'
	path_cp = r'../../networks/dbs/'
	path_out = r'../../out/'

	file_cp='corresponding-top-picks-7pathways-6overlap-withcollapsed.txt'
	
	if ((ext == 'gcount') | (ext == 'rho')):
	    df = aggregate_coeff(path_cp,file_cp,path_orig,ext)
	elif ((ext == 'pred.gcount') | (ext == 'pred.rho')):
	    df = aggregate_coeff(path_cp,file_cp,path_reg,ext)

	make_plot(path_cp,file_cp,df,ext,path_out)

if __name__ == "__main__":
    main(sys.argv)
