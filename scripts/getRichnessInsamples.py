import pandas as pd
import numpy as np
import sys, os
import argparse
import glob

def main():
    analysisDir = sys.argv[1].strip() 
    outputF = sys.argv[2].strip()
    print("gathering kraken label files from "+analysisDir)
    krakenLabelFs = glob.glob(analysisDir+"/*.krakenLabels") #kraken label files
    print(krakenLabelFs[::5])
    print("total kraken label files number=",len(krakenLabelFs))
    print("reading "+analysisDir+"/trackReads.out file")
    trackF=pd.read_csv(analysisDir+"/trackReads.out", sep=" ")[['#sample','#input','#qtrimmed']] #keep only these columns
    trackF.columns=['sname','totalNReads','nQtrimmedReads']
    print(trackF.head())
    snames=[]
    nBacterialReads=[]
    richness=[]
    for f in krakenLabelFs:
       sampleN = os.path.basename(f).split(".krakenLabels")[0]
       with open(f,'rt') as fh:
           if trackF[trackF['sname']==sampleN]['nQtrimmedReads'].values >0 and len(fh.readlines())>0:
               df = pd.read_csv(f,sep="\t",header=None)
           #if there zero qtrimmed reads or zero kraken calssified reads then output empty dataframe
           else:
               df = pd.DataFrame(columns=['sname','taxonomy'])
       snames.append(sampleN)
       nBacterialReads.append(df.shape[0])
       df.columns=['sname','taxonomy']
       richness.append(df['taxonomy'].unique().shape[0])

    df_out=pd.DataFrame({'sname':snames,'nBacterialReads':nBacterialReads,'richness':richness})
    df_out = df_out.merge(trackF, how='left', on='sname')
    df_out.loc[:,'percReadsAreBact']=df_out['nBacterialReads'].values.astype(float)/df_out['totalNReads'].values.astype(float)*100.0
    df_out.to_csv(analysisDir+"/"+outputF)

if __name__ == '__main__':
    df_out=main()
