# computes transformed rho coefficients from the original using linear regression
# output files are stored as --.pred.rho
# needs three arguments to run
# i.e. paths to the directories containing correponding-paths, original profiles, 
# and to output directory for transformed profiles.
# Example: 
# python regression.py ../../networks/dbs/ ../../graphlets/dbs/ ../../regression/predicted/dbs/
#

import pandas as pd 
import glob
import sys
import os
from scipy import stats
import numpy as np
from sklearn.linear_model import LinearRegression


def transform_profile(path_to_corr: str, inputdir: str,outdir: str, ext: str):

	f = path_to_corr + 'corresponding-top-picks.txt'
	cp = pd.read_csv(f,sep='\t',index_col=0)

	dbs = list(cp.columns.values)
	pathways = list(cp.index.values)	

	df=None
	for i in pathways:
	    for j in dbs:
	        if (pd.notnull(cp.loc[i,j])):
	            f = inputdir+cp.loc[i,j][:-4]+'-network.txt.orca.'+ext
	            df1=pd.read_csv(f,sep =' ', header=None)
	            df1=pd.DataFrame(df1[1])
	            df1=df1.rename(columns={1:j+'-'+i}).T
	            df1['db']=j
	            df1['pathway']=i
	            if df is not None:
	                df = df.append(df1)
	            else:
	                df = df1
	                
	target = df.groupby("pathway", as_index=True)[df.columns[0:len(df.columns)-2]].mean()

	predicted_table = None

	for j in range(0,len(df.columns)-2):
	    df1=None
	    for i in dbs:
	        x = df[df['db']==i][[j]]
	        pathways_in_db = df[df['db']==i][['pathway']]
	        pathways_in_db = list(pathways_in_db['pathway'])
	        y = target.loc[pathways_in_db,j]

	        reg = LinearRegression().fit(x, y)
	        predicted_x = reg.predict(x)
	        
	        data={'db':i,'pathway':pathways_in_db,j:predicted_x}
	        df2=pd.DataFrame(data)
	        if df1 is not None:
	            df1 = df1.append(df2)
	        else:
	            df1 = df2
	        
	    if predicted_table is not None:
	        predicted_table = pd.merge(predicted_table, df1, on = ['db','pathway'])
	    else:
	        predicted_table = df1
	        #print(reg.score(x,y))

	for i in pathways:
	    for j in dbs:
	        if (pd.notnull(cp.loc[i,j])):
	            to_save = predicted_table.loc[(predicted_table["db"]==j) & 
	                    (predicted_table["pathway"]==i),range(0,len(target.columns))].T
	            if ext == 'rho':
	                to_save.index = to_save.index+1
	            f = cp.loc[i,j].split('/')[-1][:-4]+'-network.txt.orca.pred.'+ext
	            fullpath = outdir+j+'/'+f
	            to_save.to_csv(fullpath, sep=' ',header=False, index=True)
	            print(f)

def main(argv):
    path_to_corr = argv[1]
    inputdir = argv[2]
    outdir = argv[3]
    ext='rho' #or gcount

    # transform the coefficients and write results to files
    transform_profile(path_to_corr,inputdir,outdir,ext) 

if __name__ == "__main__":
    main(sys.argv)