import pandas as pd
import numpy as np
import sys, os
import argparse
import csv, sqlite3
import glob

#con = sqlite3.connect(":memory:")
#cur = con.cursor()
#tables=["gene","mechanism","class","group", 
#  ####then per table:
#cur.execute("CREATE TABLE t (col1, col2);") # use your column names here
#
#with open('data.csv','rb') as fin: # `with` statement available in 2.5+
#    #csv.DictReader uses first line in file for column headings by default
#    dr = csv.DictReader(fin) # comma is default delimiter
#    to_db = [(i['col1'], i['col2']) for i in dr]
#
#cur.executemany("INSERT INTO t (col1, col2) VALUES (?, ?);", to_db)
#con.commit()
#con.close()

def main():
    parser = argparse.ArgumentParser(description='getting the missing enzymes of butyrate pathways')
    parser.add_argument('--outFile', '-o', help='output CSV file name encompassing sample names as rows', required=True)
    #to process files as open handels ready to be processed
    parser.add_argument('--inFiles', '-i', type=argparse.FileType('r') , nargs='+', help='input CSV files a file for each sample to be appended', required=True)
    #to process all files given as input file names
    #parser.add_argument('--inFiles', '-i', nargs='*', help='input CSV files a file for each sample; encompassing gene quantification results', required=True)
    args = parser.parse_args()
    inputFhs=args.inFiles
    #snames=[ff.name.split(".tsv")[0] for ff in inputFhs]
    snames=[os.path.basename(ff.name).split(".tsv")[0] for ff in inputFhs]
    print("input file number=",len(inputFhs))
    outputF=args.outFile
    dfs=[pd.read_csv(h,sep="\t") for h in inputFhs]
    print(dfs[0].head())

    for i,df in enumerate(dfs):
        if df.shape[0]>0:
            df.loc[:,'Sample']=snames[i]
            if i==0:
                res=df
            else:
                res=res.append(df)

        else:
            print(df.head())

    print("the combined csv file has a shape= ",res.shape)

    res.to_csv(outputF)

if __name__ == '__main__':
    main()
