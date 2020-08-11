import pandas as pd
import numpy as np
import sys, os
import argparse

def get_ECsPWYsFeatures(EC_PWY_FILE,sname,eggNogFile,cdsCovThreshold=0.75):
    egg=pd.read_csv(eggNogFile, index_col=0)
    ec=pd.read_csv(EC_PWY_FILE)
    ##drop redundent all
    ec = ec[ec['ButyratePwy'] != 'All']
    target_ec=[];temp=ec['EC'].apply(lambda r:target_ec.extend(r.strip().split(",")))
    target_gene=[];temp=ec['ShortName'].apply(lambda r: target_gene.extend(r.lower().strip().split(",")) if not pd.isna(r) else r)
    target_pwy=[];temp=ec['pwy'].apply(lambda r: target_pwy.extend(r.replace("ec","map").strip().split(",")))
    tt=[e.replace("map","ko") for e in target_pwy]
    target_pwy.extend(tt)
#    #expand 2.8.3.- <--this was an improvisation but I think it is better not to do it
#    tt=['2.8.3.'+str(i) for i in list(range(1,51))]
#    target_ec.remove('2.8.3.-')
#    target_ec.extend(tt)
    target_pwy=pd.unique(target_pwy)
    target_ec=pd.unique(target_ec)
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

    #convert ec and pwy to lists for easy handling later
    egg_trimmed.loc[:,'ec'] = egg_trimmed['ec'].apply(lambda r:r if pd.isnull(r) else r.strip().split(","))
    egg_trimmed.loc[:,'kEGG_PWY'] = egg_trimmed['kEGG_PWY'].apply(lambda r: r if pd.isnull(r) else r.strip().split(","))
    egg_trimmed.loc[:,'proteinName'] = egg_trimmed['proteinName'].apply(lambda r: r if pd.isnull(r) else r.replace("ko:","").lower().strip().split(","))

    #intersecting detected ECs, PWYs, protienNames for each CDs with targets
    temp1 = egg_trimmed['ec'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_ec))
    temp1 = [ e if len(e)>0 else np.nan for i,e in temp1.iteritems() ]
    egg_trimmed.loc[:,'ec_intersxn'] = temp1 

    temp2 = egg_trimmed['kEGG_PWY'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_pwy))
    temp2 = [ e if len(e)>0 else np.nan for i,e in temp2.iteritems()]
    egg_trimmed.loc[:,'pwy_intersxn'] = temp2 

    temp3 = egg_trimmed['proteinName'].apply(lambda r: {} if np.all(pd.isnull(r)) else set(r).intersection(target_gene))
    temp3 = [ e if len(e)>0 else np.nan for i,e in temp3.iteritems()]
    egg_trimmed.loc[:,'gene_intersxn'] = temp3

    #trimming away null intersections with target enzymes (keep only related to butyrate production for now)
    ind1=(~egg_trimmed['ec_intersxn'].isnull())
    ind2=(~egg_trimmed['pwy_intersxn'].isnull())
    ind3=(~egg_trimmed['gene_intersxn'].isnull())
    #any of these if not null should be accounted for
    ind_trim = ind1 | ind2 | ind3
    egg_trimmed_ecintx=egg_trimmed[ind_trim]
    ### drop lines with NaN ec_intersxn and gene_intersxn
    #drop null intersecting ECs so that we only limit the analysis to butyrate enzymes
    ind=egg_trimmed_ecintx[egg_trimmed_ecintx[['ec_intersxn','gene_intersxn']].isnull().sum(axis=1)==2].index
    egg_trimmed_ecintx = egg_trimmed_ecintx.loc[~egg_trimmed_ecintx.index.isin(ind),:]
    #No need to proceed if no features were found
    if egg_trimmed_ecintx.shape[0]==0:
        return pd.DataFrame({'ec':None, 'depth':None,'Bpwy':None},index=[sname])
    #NotApplicable:because we had the liberty to expand the 2.8.3.- enzymes we might have incorporated pathways that are not related to butyrate, we should remove these as well:
    #one could check these via egg_trimmed_ecintx[ egg_trimmed_ecintx['pwy_intersxn'].isnull() ][['kEGG_PWY','geneOntology','depth','percOfCDSCoveredAtdepth','proteinName','ec_intersxn']]
    #####To stay away from poor evidence I will trim away genes with <5 depth and %ofbasesAtDepth<0.01:
    egg_trimmed_ecintx = egg_trimmed_ecintx[~((egg_trimmed_ecintx['percOfCDSCoveredAtdepth']<0.01) & (egg_trimmed_ecintx['depth']<5))]

    ind=egg_trimmed_ecintx[egg_trimmed_ecintx[['ec_intersxn','gene_intersxn']].isnull().sum(axis=1)==2].index
    egg_trimmed_ecintx = egg_trimmed_ecintx.loc[~egg_trimmed_ecintx.index.isin(ind),:]
    
    #replace motifs, short chains ,etc with intersected target protein name by the enzyme number (so to rely on enzyme numbers at the end)
    #first get the index where either the ec or the geneName is missing
    ind=egg_trimmed_ecintx[egg_trimmed_ecintx[['ec_intersxn','gene_intersxn']].isnull().sum(axis=1)==1].index
    #replace the NAN by the EC of detected geneName
    #create translation dictionary
    ec_copy = ec.copy()
    ec_copy.set_index('ShortName', inplace=True)
    ec_copy.index = [e.lower() for e in ec_copy.index.values]
    ec_copy.loc[:,'EC'] = ec_copy['EC'].apply(lambda r:set(r.split(",")))
    trans_nameToEC = ec_copy['EC'].to_dict()
    def replaceNAN(row):
        """replace nan ECs with translation of the found geneName"""
        if pd.isnull(row['gene_intersxn']):
            return row['ec_intersxn']
        elif pd.isnull(row['ec_intersxn']):
            ecs=set()
            for i in row['gene_intersxn']:
                ecs=ecs.union(trans_nameToEC[i])
            return ecs
        else:
            print("both gene_intersxn and ec_intersxn are NaN one should remove this possibility before this step")
            sys.exit(1)

    egg_trimmed_ecintx.loc[ind,'ec_intersxn'] = egg_trimmed_ecintx.loc[ind,:].apply(lambda r: replaceNAN(r), axis=1)

    #split ecs into different columns
    egg_trimmed_ecintx_w = pd.concat([egg_trimmed_ecintx, egg_trimmed_ecintx[['cds', 'depth','ec_intersxn']].apply(lambda r: pd.Series(list(r['ec_intersxn'])), axis=1)], axis=1)
    ec_cols=egg_trimmed_ecintx_w.columns[~egg_trimmed_ecintx_w.columns.isin(egg_trimmed_ecintx.columns.values)].values
    egg_trimmed_ecintx = egg_trimmed_ecintx_w.copy()

    # add a column of distributed depth across ec columns (only non nan ecs will be counted of course)
    egg_trimmed_ecintx.loc[:,'splitDepth'] = egg_trimmed_ecintx.apply(lambda r:r['depth']/(~r.loc[ec_cols].isnull()).sum(), axis=1)
    total_depth = egg_trimmed_ecintx['depth'].sum()
    egg_trimmed_ecintx = egg_trimmed_ecintx[['cds','percOfCDSCoveredAtdepth']+ec_cols.tolist()+['splitDepth']]

    # melt the data frame such that the ec column has one entry:
    egg_trimmed_ecintx_long = pd.melt(egg_trimmed_ecintx, id_vars=['cds','percOfCDSCoveredAtdepth','splitDepth'],value_vars=ec_cols,value_name="ec",var_name="ec_ind")
    egg_trimmed_ecintx_long = egg_trimmed_ecintx_long[~egg_trimmed_ecintx_long['ec'].isnull()]
    assert(np.isclose(egg_trimmed_ecintx_long['splitDepth'].sum(),total_depth)),"depth differs between long and wide formats, check"
    egg_trimmed_ecintx_long.drop('ec_ind', axis=1, inplace =True)

    #filter CDS that are covered by < threshold coverage (default=0.75)
    egg_trimmed_ecintx_long = egg_trimmed_ecintx_long[egg_trimmed_ecintx_long['percOfCDSCoveredAtdepth']>cdsCovThreshold]

    #summarize by simply found ECs and depth per EC
    egg_trimmed_ecintx_long = pd.DataFrame(egg_trimmed_ecintx_long.groupby('ec').sum()['splitDepth'])


    ### to get coverage per butyrate pathway in terms of read depth and in terms of number of involved enzymes I need to get target ec and pwy by Butyrate pathwy
    gbBpwy=ec.groupby("ButyratePwy")
    target_ec_byBpwy={}
    target_pwy_byBpwy={}
    for Bpwy, dfBpwy in list(gbBpwy):
        target_ec_byBpwy[Bpwy]=[]; temp=dfBpwy['EC'].apply(lambda r:target_ec_byBpwy[Bpwy].extend(r.strip().split(",")))
        #if Bpwy=='4-aminobutyrate':
            #expand 2.8.3.-
            #t2=['2.8.3.'+str(i) for i in list(range(1,51))]
            #target_ec_byBpwy[Bpwy].remove('2.8.3.-')
            #target_ec_byBpwy[Bpwy].extend(t2)
        target_ec_byBpwy[Bpwy]=pd.unique(target_ec_byBpwy[Bpwy])
            
        target_pwy_byBpwy[Bpwy]=[]; temp=dfBpwy['pwy'].apply(lambda r:target_pwy_byBpwy[Bpwy].extend(r.replace("ec","map").strip().split(",")))
        tt=[e.replace("map","ko") for e in target_pwy_byBpwy[Bpwy]]
        target_pwy_byBpwy[Bpwy].extend(tt)
        target_pwy_byBpwy[Bpwy]=pd.unique(target_pwy_byBpwy[Bpwy])
    
    def intrsxn(row, Bpwy):
        intersect=set(row).intersection(target_ec_byBpwy[Bpwy])
        if len(intersect)>0:
            return 1
        else:
            return 0

    
    for Bpwy, dfBpwy in list(gbBpwy):
        egg_trimmed_ecintx_long.loc[:,Bpwy] = egg_trimmed_ecintx_long.apply(lambda r:intrsxn([r.name],Bpwy), axis=1)
    #split the depth yet again between multiply represented pwys
    pwy_cols = list(set(egg_trimmed_ecintx_long.columns.values)-{'splitDepth'})
    egg_trimmed_ecintx_long.loc[:,'splitDepth2'] = egg_trimmed_ecintx_long.apply(lambda r:r['splitDepth']/(r.loc[pwy_cols]).sum(), axis=1)
    egg_trimmed_ecintx_long.drop('splitDepth', axis=1, inplace=True)
    egg_trimmed_ecintx_long.rename(columns={'splitDepth2':'depth'}, inplace=True)

    #melt it to long format
    egg_trimmed_ecintx_long.reset_index(inplace=True)
    egg_trimmed_ecintx_long = pd.melt( egg_trimmed_ecintx_long, id_vars=['ec','depth'], value_vars=pwy_cols, var_name='Bpwy', value_name='detected')

    #get rid of zero detected rows
    egg_trimmed_ecintx_long = egg_trimmed_ecintx_long[egg_trimmed_ecintx_long['detected']>0]
    #then remove the detected column
    egg_trimmed_ecintx_long.drop('detected', axis=1, inplace=True)

    #set index to sample name
    egg_trimmed_ecintx_long.index=[sname]*egg_trimmed_ecintx_long.shape[0]
    return egg_trimmed_ecintx_long #dataframe with sample as a row and features are in long format with the number of reads for the features

def main():
    parser = argparse.ArgumentParser(description='getting the missing enzymes of butyrate pathways')
    parser.add_argument('--outFile', '-o', help='output CSV file name encompassing samples as rows, pathway:ec features as columns and depth as entries', required=True)
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
    for i,ff in enumerate(files):
        try:
            sname=os.path.basename(os.path.dirname(ff)).split("eggnog_map.")[1]
        except:
            sname="unknown"
        print("processing sample", sname)
        if i==0:
            allsamples_df = get_ECsPWYsFeatures(EC_PWY_FILE,sname,ff,cdsCovThreshold)
        else:
            allsamples_df = allsamples_df.append(get_ECsPWYsFeatures(EC_PWY_FILE,sname,ff,cdsCovThreshold)) 
        #save this step
        allsamples_df.to_csv(outputFileName)
    print(allsamples_df.head())

if __name__ == '__main__':
    main()


