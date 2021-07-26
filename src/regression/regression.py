# computes transformed rho coefficients from the original using linear regression
# output files are stored as --.pred.rho
# needs three arguments to run
# i.e. paths to the directories containing correponding-paths, 
# path directories of original profile by database
# and to output 'dbs' directory for transformed profiles.
# An optional fourth argument 'gcount' can be given for graphlet profiles. 
# By default it runs on rho coefficient.
# Transformation normalizes the graphlet profiles for regression.
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

	f = path_to_corr + 'corresponding-top-picks-7pathways-6overlap-withcollapsed.txt'
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

	if ((ext == 'gcount')):
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
    ext='rho' # rho or gcount
    if len(argv) > 4:
        ext = argv[4]

    # transform the coefficients and write results to files
    transform_profile(path_to_corr,inputdir,outdir,ext) 

if __name__ == "__main__":
    main(sys.argv)