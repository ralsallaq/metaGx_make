import pandas as pd
import numpy as np
import sys, os
import argparse

def get_ECsPWYsFeatures(EC_PWY_FILE,sname,eggNogFile,cdsCovThreshold=0.75):
    egg=pd.read_csv(eggNogFile, index_col=0)
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
    ##### trimming the egg dataframe from zero depth regions and the genome bits######
    ##### add reads that covers various bases of CDS and add the bases covered as well as the percentage of bases covered #####
    egg=egg[egg['cds']!='genome']
    #drop lines corresponding to non-covered portions of the CDS
    egg = egg[egg['depth']>0]
    #save the unique features (one per CDs) from eggNOG
    eggNog_features = egg.loc[:,~egg.columns.isin(['depth', 'num_bases_atDepth', 'sizeOfGene','percOfCDSCoveredAtdepth'])]
    eggNog_features.drop_duplicates(inplace=True)
    gpbycds=egg.groupby("cds")
    egg_trimmed = gpbycds.agg({'depth':'sum','num_bases_atDepth':'sum','sizeOfGene':'first','percOfCDSCoveredAtdepth':'sum'})
    egg_trimmed.reset_index(inplace=True)
    egg_trimmed = egg_trimmed.merge(eggNog_features, how='right', on='cds')

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

    ### drop lines with NaN everything except pwy
    ind=egg_trimmed_ecintx[egg_trimmed_ecintx[['ec_intersxn','gene_intersxn','K_intersxn']].isnull().sum(axis=1)==3].index
    egg_trimmed_ecintx = egg_trimmed_ecintx.loc[~egg_trimmed_ecintx.index.isin(ind),:]

    #No need to proceed if no features were found
    if egg_trimmed_ecintx.shape[0]==0:
        return pd.DataFrame({'kwds':None, 'splitDepth':None,'seedOrth':None},index=[sname])

    #####To stay away from poor evidence I will trim away genes with %ofbasesAtDepth<=cdsCovThreshold:
    #filter CDS that are covered by < threshold coverage (default=0.75)
    egg_trimmed_ecintx = egg_trimmed_ecintx[egg_trimmed_ecintx['percOfCDSCoveredAtdepth']>cdsCovThreshold]

    ### make a column of detected kwds
    egg_trimmed_ecintx.loc[:,'kwds_intersxn'] = egg_trimmed_ecintx['ec_intersxn'].values
    ind = egg_trimmed_ecintx.loc[:,'kwds_intersxn'].isnull()
    egg_trimmed_ecintx.loc[ind,'kwds_intersxn'] = egg_trimmed_ecintx.loc[ind,['gene_intersxn','K_intersxn']].apply(lambda r: r['gene_intersxn'] if pd.isnull(r['K_intersxn']) else r['K_intersxn'], axis=1)

    #split kwds into different columns
    egg_trimmed_ecintx_w = pd.concat([egg_trimmed_ecintx, egg_trimmed_ecintx[['cds', 'depth','kwds_intersxn']].apply(lambda r: pd.Series(list(r['kwds_intersxn'])), axis=1)], axis=1)
    kwds_cols=egg_trimmed_ecintx_w.columns[~egg_trimmed_ecintx_w.columns.isin(egg_trimmed_ecintx.columns.values)].values
    egg_trimmed_ecintx = egg_trimmed_ecintx_w.copy()

    # add a column of distributed depth across ec columns (only non nan ecs will be counted of course)
    egg_trimmed_ecintx.loc[:,'splitDepth'] = egg_trimmed_ecintx.apply(lambda r:r['depth']/(~r.loc[kwds_cols].isnull()).sum(), axis=1)
    total_depth = egg_trimmed_ecintx['depth'].sum()
    egg_trimmed_ecintx = egg_trimmed_ecintx[['cds','percOfCDSCoveredAtdepth','seedOrth']+kwds_cols.tolist()+['splitDepth']]

    # melt the data frame such that the kwds column has one entry:
    egg_trimmed_ecintx_long = pd.melt(egg_trimmed_ecintx, id_vars=['cds','percOfCDSCoveredAtdepth','splitDepth','seedOrth'],value_vars=kwds_cols,value_name="kwds",var_name="kwds_ind")
    egg_trimmed_ecintx_long = egg_trimmed_ecintx_long[~egg_trimmed_ecintx_long['kwds'].isnull()]
    assert(np.isclose(egg_trimmed_ecintx_long['splitDepth'].sum(),total_depth)),"depth differs between long and wide formats, check"
    egg_trimmed_ecintx_long.drop('kwds_ind', axis=1, inplace =True)

    #summarize by simply found kwds and depth per kwd
    gpbykw = egg_trimmed_ecintx_long.groupby('kwds')
    egg_trimmed_ecintx_long = pd.DataFrame([gpbykw.sum()['splitDepth'],gpbykw.agg(lambda r:"|".join(r['seedOrth'].values))['seedOrth']]).T

    #set index to sample name
    egg_trimmed_ecintx_long.reset_index(inplace=True)
    egg_trimmed_ecintx_long.index=[sname]*egg_trimmed_ecintx_long.shape[0]

    return egg_trimmed_ecintx_long #dataframe in long format with each unique kwd found as a row with the corresponding number of reads found



def main():
    parser = argparse.ArgumentParser(description='getting the missing enzymes of butyrate pathways')
    parser.add_argument('--outFile', '-o', help='output CSV file name encompassing samples as rows, pathways as columns and missing enzymes as entries', required=True)
    #to process files as open handels ready to be processed
    #parser.add_argument('--inFiles', '-i', type=argparse.FileType('r') , nargs='+', help='input CSV files a file for each sample; encompassing gene quantification results', required=True)
    #to process all files given as input file names
    parser.add_argument('--inFiles', '-i', nargs='*', help='input CSV files a file for each sample; encompassing gene quantification results', required=True)
    parser.add_argument('--ECFile', '-ec', help='input CSV file encompassing target pathways and enzyme numbers (ECs)', required=True)
    parser.add_argument('--cdsCovThreshold', '-t', help='threshold of CDS coverage by reads to consider the CDS', default=0.75)
    args = parser.parse_args()
    EC_PWY_FILE=args.ECFile
    cdsCovThreshold = float(args.cdsCovThreshold)
    outputFileName=args.outFile
    print("output", outputFileName)
    files=args.inFiles
    print(files)
    print("EC file is set to",EC_PWY_FILE)
    for i,ff in enumerate(files):
        try:
            sname=os.path.basename(os.path.dirname(ff)).split("eggnog_map.")[1]
        except:
            sname="unknown"
        if i==0:
            allsamples_df = get_ECsPWYsFeatures(EC_PWY_FILE,sname,ff,cdsCovThreshold)
        else:
            allsamples_df = allsamples_df.append(get_ECsPWYsFeatures(EC_PWY_FILE,sname,ff,cdsCovThreshold)) 
        #save this step
        allsamples_df.to_csv(outputFileName)
    print(allsamples_df.head())

if __name__ == '__main__':
    main()


