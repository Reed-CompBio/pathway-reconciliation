#
# Tobias Rubel | rubelato@reed.edu
# Reed CompBio
#
# This just allows me to easily combine two or more csvs into one from the command line
#
#
import pandas as pd
import sys

def combine(lat):
    return pd.concat([pd.read_csv(x,index_col=0,sep='\s+',names=[x.split('/')[-1].split('.')[0]]) for x in lat], axis=1, sort=True)


def main(argv):
    df = combine(argv[1:-1])
    df.to_csv(argv[-1])

if __name__ == "__main__":
    main(sys.argv)
