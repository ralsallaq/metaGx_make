import pandas as pd
import numpy as np
import sys, os
import argparse

def get_missingECsPerPWY(EC_PWY_FILE, sname):
    egg=pd.read_csv("analysis/eggnog_map."+sname+"/egmap.emapper.annotations.Wcov.csv", index_col=0)
    #ec=pd.read_csv("../../butyrate_EC_pathways.csv")
    ec=pd.read_csv(EC_PWY_FILE)
    ##drop redundent all
    ec = ec[ec['ButyratePwy'] != 'All']
    target_ec=[];temp=ec['EC'].apply(lambda r:target_ec.extend(r.strip().split(",")))
    target_pwy=[];temp=ec['pwy'].apply(lambda r: target_pwy.extend(r.replace("ec","map").strip().split(",")))
    tt=[e.replace("map","ko") for e in target_pwy]
    target_pwy.extend(tt)
    #expand 2.8.3._
    tt=['2.8.3.'+str(i) for i in list(range(1,51))]
    target_ec.remove('2.8.3._')
    target_ec.extend(tt)
    target_pwy=pd.unique(target_pwy)
    target_ec=pd.unique(target_ec)
    ##### trimming the egg dataframe to the maximum depth for each cds######
    gpbycds=egg.groupby("cds")
    mx_ind=gpbycds['depth'].idxmax()
    egg_trimmed=egg.loc[mx_ind,:]
    egg_trimmed.loc[:,'ec'] = egg_trimmed['ec'].apply(lambda r:r if pd.isnull(r) else r.strip().split(","))
    egg_trimmed['kEGG_PWY'] = egg_trimmed['kEGG_PWY'].apply(lambda r: r if pd.isnull(r) else r.strip().split(","))
    temp1 = egg_trimmed['ec'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_ec))
    temp1 = [ e if len(e)>0 else np.nan for i,e in temp1.iteritems() ]
    egg_trimmed.loc[:,'ec_intersxn'] = temp1 
    temp2 = egg_trimmed['kEGG_PWY'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_pwy))
    temp2 = [ e if len(e)>0 else np.nan for i,e in temp2.iteritems()]
    egg_trimmed.loc[:,'pwy_intersxn'] = temp2 
    #trimming away null intersections with target enzymes (keep only related to butyrate production for now)
    egg_trimmed_ecintx=egg_trimmed[~egg_trimmed['ec_intersxn'].isnull()] 
    #because we had the liberty to expand the 2.8.3._ enzymes we might have incorporated pathways that are not related to butyrate, we should remove these as well:
    #one could check these via egg_trimmed_ecintx[ egg_trimmed_ecintx['pwy_intersxn'].isnull() ][['kEGG_PWY','geneOntology','depth','percOfCDSCoveredAtdepth','proteinName','ec_intersxn']]
    egg_trimmed_ecintx = egg_trimmed_ecintx[ ~ egg_trimmed_ecintx['pwy_intersxn'].isnull() ]
    ### to get coverage per butyrate pathway in terms of read depth and in terms of number of involved enzymes I need to get target ec and pwy by Butyrate pathwy
    gbBpwy=ec.groupby("ButyratePwy")
    target_ec_byBpwy={}
    target_pwy_byBpwy={}
    for Bpwy, dfBpwy in list(gbBpwy):
        target_ec_byBpwy[Bpwy]=[]; temp=dfBpwy['EC'].apply(lambda r:target_ec_byBpwy[Bpwy].extend(r.strip().split(",")))
        if Bpwy=='4-aminobutyrate':
            #expand 2.8.3._
            t2=['2.8.3.'+str(i) for i in list(range(1,51))]
            target_ec_byBpwy[Bpwy].remove('2.8.3._')
            target_ec_byBpwy[Bpwy].extend(t2)
        target_ec_byBpwy[Bpwy]=pd.unique(target_ec_byBpwy[Bpwy])
            
        target_pwy_byBpwy[Bpwy]=[]; temp=dfBpwy['pwy'].apply(lambda r:target_pwy_byBpwy[Bpwy].extend(r.replace("ec","map").strip().split(",")))
        tt=[e.replace("map","ko") for e in target_pwy_byBpwy[Bpwy]]
        target_pwy_byBpwy[Bpwy].extend(tt)
        target_pwy_byBpwy[Bpwy]=pd.unique(target_pwy_byBpwy[Bpwy])
    
    ###now we get into each pathway: for a given butyrate pathway (say '4-aminobutyrate') if at least one pathway is detected and out of all ECs a minimum EC-1 is detected (all enzymes but one) 
    
    for Bpwy, dfBpwy in list(gbBpwy):
        temp1 = egg_trimmed_ecintx['pwy_intersxn'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_pwy_byBpwy[Bpwy]))
        temp1 = [ e if len(e)>0 else np.nan for i,e in temp1.iteritems()]
        egg_trimmed_ecintx.loc[:,Bpwy+'_pwyintrsxn']= temp1
    
        temp2 = egg_trimmed_ecintx['ec_intersxn'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_ec_byBpwy[Bpwy]))
        temp2 = [ e if len(e)>0 else np.nan for i,e in temp2.iteritems() ]
        egg_trimmed_ecintx.loc[:,Bpwy+'_ecintrsxn']= temp2
    
    #####To stay away from poor evidence I will trim away genes with <5 depth and %ofbasesAtDepth<0.01:
    egg_trimmed_ecintx = egg_trimmed_ecintx[~((egg_trimmed_ecintx['percOfCDSCoveredAtdepth']<0.01) & (egg_trimmed_ecintx['depth']<5))]

    """One can output here for each sample the following data: egg_trimmed_ecintx[['depth','percOfCDSCoveredAtdepth','4-aminobutyrate_pwyintrsxn','4-aminobutyrate_ecintrsxn','Glutarate_pwyintrsxn','Glutarate_ecintrsxn','Lysine_pwyintrsxn','Lysine_ecintrsxn','Pyruvate_pwyintrsxn', 'Pyruvate_ecintrsxn']] ; which gives depth data and percentage coverage for each found enzyme/pwy"""
    
    #####Now we can quantify which butyrate pathway(s) are present (one enzyme from the pathway can be missing and still count as present pathway)
    Bpways=[b for b, dfBpwy in list(gbBpwy)]
    butyrate_pwys=pd.DataFrame(columns=[i for i in Bpways], index=[sname])
    for p in Bpways:
        pwECs=set()
        #exclude null values
        xx1=egg_trimmed_ecintx[~ egg_trimmed_ecintx[p+'_ecintrsxn'].isnull()]
        for e in xx1[p+'_ecintrsxn'].values:
            pwECs= pwECs.union(e)
        res1=set(target_ec_byBpwy[p])-pwECs
        res2=res1-set(['2.8.3.'+str(i) for i in range(51)])
        if len(res2)>1:
            print("pathways that have more than one missing enzymes", p, res2)
            butyrate_pwys.loc[:,p]=",".join(res2)
        else:
            butyrate_pwys.loc[:,p]=",".join(res2)
    print(butyrate_pwys)
    return butyrate_pwys #dataframe with each pathway as a column and a list of missing enzymes

def main():
    parser = argparse.ArgumentParser(description='getting the missing enzymes of butyrate pathways')
    parser.add_argument('--outFile', '-o', help='output CSV file name encompassing samples as rows, pathways as columns and missing enzymes as entries', required=True)
    #to process files as open handels ready to be processed
    #parser.add_argument('--inFiles', '-i', type=argparse.FileType('r') , nargs='+', help='input CSV files a file for each sample; encompassing gene quantification results', required=True)
    #to process all files given as input file names
    parser.add_argument('--inFiles', '-i', nargs='*', help='input CSV files a file for each sample; encompassing gene quantification results', required=True)
    parser.add_argument('--ECFile', '-ec', nargs='*', help='input CSV file encompassing target pathways and enzyme numbers (ECs)', required=True)
    args = parser.parse_args()
    EC_PWY_FILE=args.ECFile[0]
    print("EC file is set to",EC_PWY_FILE)
    outputFileName=args.outFile
    print("output", outputFileName)
    files=args.inFiles
    print(files)
    for i,ff in enumerate(files):
        sname=os.path.basename(os.path.dirname(ff)).split("eggnog_map.")[1]
        if i==0:
            allsamples_df = get_missingECsPerPWY(EC_PWY_FILE,sname)
        else:
            allsamples_df = allsamples_df.append(get_missingECsPerPWY(EC_PWY_FILE,sname)) 
        #save this step
        allsamples_df.to_csv(outputFileName)
    print(allsamples_df.head())

if __name__ == '__main__':
    main()


