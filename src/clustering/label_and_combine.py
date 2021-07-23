#
# in: two formatted csvs
# out: 1 csv of the two combined but with labels of _1 on items from the first and _2 on items from the second.

import pandas as pd 
import os 
import sys


def main(argv):
    df1 = pd.read_csv(argv[1],index_col=0)
    df2 = pd.read_csv(argv[2],index_col=0)
    df1.columns = ['{}_1'.format(x) for x in df1.columns]
    df2.columns = ['{}_2'.format(x) for x in df2.columns]
    df3 = pd.concat([df1,df2],axis=1)
    df3.to_csv(argv[3])


if __name__ == "__main__":
    main(sys.argv)
