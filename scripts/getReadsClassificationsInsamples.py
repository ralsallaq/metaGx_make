import pandas as pd
import numpy as np
import sys, os
import argparse
import glob

def main():
    analysisDir = sys.argv[1].strip() 
    outputF = sys.argv[2].strip()
    print("gathering kraken label files from "+analysisDir)
    krakenReportFs = glob.glob(analysisDir+"/*.krakenReport") #kraken2 report files
    print(krakenReportFs[::5])
    print("total kraken report files number=",len(krakenReportFs))
    print("reading "+analysisDir+"/trackReads.out file")
    trackF=pd.read_csv(analysisDir+"/trackReads.out", sep=" ")[['#sample','#input','#qtrimmed']] #keep only these columns
    trackF.columns=['sname','totalNReads','nQtrimmedReads']
    print(trackF.head())
    print("dividing number of reads in track by 2 to get the number of pair-end reads")
    trackF.loc[:,'totalNReads'] = trackF.loc[:,'totalNReads']/2.0 #paired reads in track
    trackF.loc[:,'nQtrimmedReads'] = trackF.loc[:,'nQtrimmedReads']/2.0 #paired reads in track
    print(trackF.head())
    snames=[]
    nBacteriaReads=[]
    nEukaryotaReads=[]
    nArcheaReads=[]
    nVirusReads=[]
    richness=[]
    for f in krakenReportFs:
       sampleN = os.path.basename(f).split(".krakenReport")[0]
       with open(f,'rt') as fh:
           try:
               df = pd.read_csv(f,sep="\t",header=None)
               df.columns=['percReads','numReads','assigned','rankC','tax_id','name']
               ##get rid of the indentations
               df.loc[:,'name'] = df['name'].apply(lambda r: r.strip())
               #add sname as a column
               df.loc[:,'sname']=sampleN
           except:
               df = pd.DataFrame(columns=['percReads','numReads','assigned','rankC','tax_id','name'])
               #add sname as a column
               df.loc[:,'sname']=sampleN

       snames.append(sampleN)
       ranks=['Bacteria','Eukaryota','Archaea','Viruses']
       lists=[nBacteriaReads, nEukaryotaReads, nArcheaReads, nVirusReads]
       for j,r in enumerate(ranks):
           idx = idx= df[(df['rankC']=="D") & (df['name']==r)].index
           try:
               lists[j].append(df.loc[idx,'numReads'].values[0])
           except IndexError:
               lists[j].append(np.nan) #no records found for that rank
               print(sampleN,j,lists[j][-1],idx,df.loc[idx,'numReads'].values)
       ##to get bacterial records in the df
       xx=df[df['rankC']=="D"]
       #index in df for when bacterial records start
       idx1 = xx[xx['name']=="Bacteria"].index.values[0]
       #index in df for when the next to bacterial records start (e.g. Eukaryota)
       idx2 = xx.loc[idx1:,:].iloc[1].name
       #bacterial reads are in between
       df_bacteria = df.loc[idx1:idx2-1,:]
       ## richness is the number of unique tax_ids at the species level
       richness.append(df_bacteria[df_bacteria['rankC']=="S"]['tax_id'].unique().shape[0])

    df_out=pd.DataFrame({'sname':snames,'nBacteriaReads':nBacteriaReads,'nEukaryotaReads':nEukaryotaReads,'nArcheaReads':nArcheaReads,'nVirusReads':nVirusReads,'BactSpeciesRichness':richness})
    df_out = df_out.merge(trackF, how='left', on='sname')
    df_out.loc[:,'percReadsAreBact']=df_out['nBacteriaReads'].values.astype(float)/df_out['nQtrimmedReads'].values.astype(float)*100.0
    df_out.to_csv(analysisDir+"/"+outputF)

if __name__ == '__main__':
    df_out=main()
