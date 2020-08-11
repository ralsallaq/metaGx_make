import pandas as pd
import numpy as np
import sys, os
import argparse

def get_missingECsPerPWY(EC_PWY_FILE, sname):
    egg=pd.read_csv("analysis/eggnog_map."+sname+"/egmap.emapper.annotations.Wcov.csv", index_col=0)
    #ec=pd.read_csv("../bileSalt_EC_pathways.csv")
    ec=pd.read_csv(EC_PWY_FILE)
    ##drop redundent all
    target_ec=[];temp=ec['EC'].apply(lambda r: target_ec.extend(r.strip().split(",")) if not pd.isna(r) else r)
    target_K=[];temp=ec['KO'].apply(lambda r: target_K.extend(r.strip().split(",")) if not pd.isna(r) else r)
    target_gene=[];temp=ec['ShortName'].apply(lambda r: target_gene.extend(r.lower().strip().split(",")) if not pd.isna(r) else r)
    target_pwy=[];temp=ec['pwy'].apply(lambda r: target_pwy.extend(r.replace("ec","map").strip().split(",")) if not pd.isna(r) else r)
    tt=[e.replace("map","ko") for e in target_pwy]
    target_pwy.extend(tt)
    target_pwy=pd.unique(target_pwy)
    target_ec=pd.unique(target_ec)
    target_K=pd.unique(target_K)
    target_gene=pd.unique(target_gene)
    ##### trimming the egg dataframe to the maximum depth for each cds######
    gpbycds=egg.groupby("cds")
    mx_ind=gpbycds['depth'].idxmax()
    egg_trimmed=egg.loc[mx_ind,:]
    #### splitting pathways and ecs to lists for easy handling later
    egg_trimmed.loc[:,'ec'] = egg_trimmed['ec'].apply(lambda r:r if pd.isnull(r) else r.strip().split(","))
    egg_trimmed['kEGG_PWY'] = egg_trimmed['kEGG_PWY'].apply(lambda r: r if pd.isnull(r) else r.strip().split(","))
    egg_trimmed['kEGG_ko'] = egg_trimmed['kEGG_ko'].apply(lambda r: r if pd.isnull(r) else r.replace("ko:","").strip().split(","))
    egg_trimmed['proteinName'] = egg_trimmed['proteinName'].apply(lambda r: r if pd.isnull(r) else r.replace("ko:","").lower().strip().split(","))

    ###intersecting the ECs detected in each CDS with the target ECs and the PWYs with the target PWYs
    temp1 = egg_trimmed['ec'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_ec))
    temp1 = [ e if len(e)>0 else np.nan for i,e in temp1.iteritems() ]
    egg_trimmed.loc[:,'ec_intersxn'] = temp1 

    temp2 = egg_trimmed['kEGG_PWY'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_pwy))
    temp2 = [ e if len(e)>0 else np.nan for i,e in temp2.iteritems()]
    egg_trimmed.loc[:,'pwy_intersxn'] = temp2 

    temp3 = egg_trimmed['kEGG_ko'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_K))
    temp3 = [ e if len(e)>0 else np.nan for i,e in temp3.iteritems()]
    egg_trimmed.loc[:,'K_intersxn'] = temp3 

    temp4 = egg_trimmed['proteinName'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_gene))
    temp4 = [ e if len(e)>0 else np.nan for i,e in temp4.iteritems()]
    egg_trimmed.loc[:,'gene_intersxn'] = temp4 

    #####trimming away null intersections with target enzymes (keep only related to bile acids); this gives all the CDS that encodes the target enzymes with all the eggnogg annotations for them
    ind1=(~egg_trimmed['ec_intersxn'].isnull())
    ind2=(~egg_trimmed['pwy_intersxn'].isnull())
    ind3=(~egg_trimmed['K_intersxn'].isnull())
    ind4=(~egg_trimmed['gene_intersxn'].isnull())
    #any of these if not null should be accounted for
    ind_trim = ind1 | ind2 | ind3 | ind4
    egg_trimmed_ecintx=egg_trimmed[ind_trim] 
    #####To stay away from poor evidence I will trim away genes with <5 depth and %ofbasesAtDepth<0.01:
    egg_trimmed_ecintx = egg_trimmed_ecintx[~((egg_trimmed_ecintx['percOfCDSCoveredAtdepth']<0.01) & (egg_trimmed_ecintx['depth']<5))]

    ### to get coverage per bile acid pathway in terms of read depth and in terms of number of involved enzymes I need to get target ec and pwy by Butyrate pathwy
    gbBpwy=ec.groupby("BileAcidsPwy")
    target_ec_byBpwy={}
    target_pwy_byBpwy={}
    target_K_byBpwy={}
    target_gene_byBpwy={}
    for Bpwy, dfBpwy in list(gbBpwy):
        target_ec_byBpwy[Bpwy]=[]; temp=dfBpwy['EC'].apply(lambda r:target_ec_byBpwy[Bpwy].extend(r.strip().split(",")) if not pd.isna(r) else r)
        target_ec_byBpwy[Bpwy]=pd.unique(target_ec_byBpwy[Bpwy])
            
        target_pwy_byBpwy[Bpwy]=[]; temp=dfBpwy['pwy'].apply(lambda r:target_pwy_byBpwy[Bpwy].extend(r.replace("ec","map").strip().split(",")) if not pd.isna(r) else r)
        tt=[e.replace("map","ko") for e in target_pwy_byBpwy[Bpwy]]
        target_pwy_byBpwy[Bpwy].extend(tt)
        target_pwy_byBpwy[Bpwy]=pd.unique(target_pwy_byBpwy[Bpwy])

        target_K_byBpwy[Bpwy]=[]; temp=dfBpwy['KO'].apply(lambda r:target_K_byBpwy[Bpwy].extend(r.strip().split(",")) if not pd.isna(r) else r)
        target_K_byBpwy[Bpwy]=pd.unique(target_K_byBpwy[Bpwy])

        target_gene_byBpwy[Bpwy]=[]; temp=dfBpwy['ShortName'].apply(lambda r:target_gene_byBpwy[Bpwy].extend(r.lower().strip().split(",")) if not pd.isna(r) else r)
        target_gene_byBpwy[Bpwy]=pd.unique(target_gene_byBpwy[Bpwy])
    
    ###now we get into each pathway: for a given bile acid pathway (there is only one) if at least one pathway is detected and out of all ECs a minimum EC-1 is detected (all enzymes but one) 
    
    for Bpwy, dfBpwy in list(gbBpwy):
        temp1 = egg_trimmed_ecintx['pwy_intersxn'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_pwy_byBpwy[Bpwy]))
        temp1 = [ e if len(e)>0 else np.nan for i,e in temp1.iteritems()]
        egg_trimmed_ecintx.loc[:,Bpwy+'_pwyintrsxn']= temp1
    
        temp2 = egg_trimmed_ecintx['ec_intersxn'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_ec_byBpwy[Bpwy]))
        temp2 = [ e if len(e)>0 else np.nan for i,e in temp2.iteritems() ]
        egg_trimmed_ecintx.loc[:,Bpwy+'_ecintrsxn']= temp2

        temp3 = egg_trimmed_ecintx['K_intersxn'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_K_byBpwy[Bpwy]))
        temp3 = [ e if len(e)>0 else np.nan for i,e in temp3.iteritems()]
        egg_trimmed_ecintx.loc[:,Bpwy+'_Kintrsxn']= temp3

        temp4 = egg_trimmed_ecintx['gene_intersxn'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_gene_byBpwy[Bpwy]))
        temp4 = [ e if len(e)>0 else np.nan for i,e in temp4.iteritems()]
        egg_trimmed_ecintx.loc[:,Bpwy+'_geneintrsxn']= temp4
    
    
    #####Now we can quantify which bile acid pathway(s) are present (one enzyme from the pathway can be missing and still count as present pathway)
    Bpways=[b for b, dfBpwy in list(gbBpwy)]
    bileacids=pd.DataFrame(columns=[i for i in Bpways], index=[sname])
    for p in Bpways:
        pwIDs=set()
        #exclude null values
        ind1=(~egg_trimmed_ecintx[p+'_ecintrsxn'].isnull())
        ind2=(~egg_trimmed_ecintx[p+'_pwyintrsxn'].isnull())
        ind3=(~egg_trimmed_ecintx[p+'_Kintrsxn'].isnull())
        ind4=(~egg_trimmed_ecintx[p+'_geneintrsxn'].isnull())
        #any of these if not null should be accounted for
        ind_trim = ind1 | ind2 | ind3 | ind4

        xx1=egg_trimmed_ecintx[ind_trim].copy()
        for e in xx1[p+'_ecintrsxn'].values:
            if not pd.isna(e):
                ec_=['ec:'+e.pop()]
                pwIDs= pwIDs.union(ec_)
        for k in xx1[p+'_Kintrsxn'].values:
            if not pd.isna(k):
                ko_=['ko:'+k.pop()]
                pwIDs= pwIDs.union(ko_)
        for g in xx1[p+'_geneintrsxn'].values:
            if not pd.isna(g):
                gene_=['gene:'+g.pop()]
                pwIDs= pwIDs.union(gene_)

        bileacids.loc[:,p]=",".join(pwIDs)
    print(bileacids)
    return bileacids #dataframe with each pathway as a column and a list of identified targets (ecs,kos, genes)

def main():
    parser = argparse.ArgumentParser(description='getting the missing enzymes of butyrate pathways')
    parser.add_argument('--outFile', '-o', help='output CSV file name encompassing samples as rows, pathways as columns and missing enzymes as entries', required=True)
    #to process files as open handels ready to be processed
    #parser.add_argument('--inFiles', '-i', type=argparse.FileType('r') , nargs='+', help='input CSV files a file for each sample; encompassing gene quantification results', required=True)
    #to process all files given as input file names
    parser.add_argument('--inFiles', '-i', nargs='*', help='input CSV files a file for each sample; encompassing gene quantification results', required=True)
    parser.add_argument('--ECFile', '-ec', nargs='*', help='input CSV file encompassing target pathways and enzyme numbers (ECs)', required=True)
    args = parser.parse_args()
    outputFileName=args.outFile
    print("output", outputFileName)
    files=args.inFiles
    print(files)
    EC_PWY_FILE=args.ECFile[0]
    print("EC file is set to",EC_PWY_FILE)
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


