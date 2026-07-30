#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 08:14:02 2021

@author: saulo
"""
import pathlib, os, sys, csv, glob
import pandas as pd
import numpy as np
from Bio import SeqIO
from functools import reduce
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from ast import literal_eval
import re

pd.set_option('display.max_columns', 10)

def prepare_extraction(path):
    
    directory = os.path.join(path,'ARTS','ARTS_Extractor')
    os.makedirs(directory, exist_ok=True)


def df2xlsx(path, sheet_name, df):
    df = df.copy().reset_index(drop=True)
    writer = pd.ExcelWriter(path, engine='xlsxwriter')

    df.to_excel(writer, sheet_name=sheet_name, startrow=1, header=False, index=False)
    worksheet = writer.sheets[sheet_name]

    (max_row, max_col) = df.shape
    column_settings = [{'header': column} for column in df.columns]
    worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings})
    worksheet.set_column(0, max_col - 1, 15)

    writer._save()


def ARTS_overview(path):

    path = os.path.join(path,'ARTS')

    csv.field_size_limit(sys.maxsize)

    directory=os.path.join(path,'ARTS_Extractor/')    
    print('>>> OVERVIEW <<<\n')
    print('Getting files...')
    
    referencia_fasta = list((pathlib.Path(path).glob("**/*.fna")))
    referencia_clusters = list(pathlib.Path(path).glob("**/*clust.tsv"))
    referencia_knownhits = list(pathlib.Path(path).glob("**/**/knownhits.tsv"))
    referencia_coregenes = list(pathlib.Path(path).glob("**/**/coretable.tsv"))
    referencia_duplic = list(pathlib.Path(path).glob("**/**/duptable.tsv"))
    
    CDS_count = []
    clusters_count = []
    hits_count = []
    core_count = []
    dup_count = []
    bgc_prox_count = []
    
    samples = []
    print('Done. Now, time to get informations.')
    
    for file in referencia_fasta:
        file = str(file)
        if "coregenes" not in file:
            base_fasta = str(os.path.basename(file)[:-4])
            if base_fasta not in samples:
                samples.append(base_fasta)
            recs = list(SeqIO.parse(file, 'fasta'))
            count_cds = [base_fasta, len(recs)]
            CDS_count.append(count_cds)
        
    for file in referencia_clusters:
        base_clust = str(os.path.basename(file)[:-10])
        if base_clust not in samples:
            samples.append(base_clust)
        clust = sum(1 for line in open(file))-1
        count_clust = [base_clust, clust]
        clusters_count.append(count_clust)
    
    for file in referencia_knownhits:
        base_known = str(os.path.abspath(file)).split('/')[-3]
        if base_known not in samples:
            samples.append(base_known)
        hits_file = sum(1 for line in open(file))-1
        count_hits = [base_known, hits_file]
        hits_count.append(count_hits)
    
    for file in referencia_coregenes:
        base_core = str(os.path.abspath(file)).split('/')[-3]
        if base_core not in samples:
            samples.append(base_core)
        hits_core = sum(1 for line in open(file))-1
        countcore = [base_core, hits_core]
        core_count.append(countcore)
        read_core = csv.reader(open(file), delimiter="\t")
        bgc_prox = 0
        for row in read_core:
            if row[4] == 'Yes':
                bgc_prox += 1
        bgc_c = [base_core, bgc_prox]
        bgc_prox_count.append(bgc_c)
        
    for file in referencia_duplic:
        base_dup = str(os.path.abspath(file)).split('/')[-3]
        if base_dup not in samples:
            samples.append(base_dup)
        dup = sum(1 for line in open(file))-1
        countdup = [base_dup, dup]
        dup_count.append(countdup)
    
    CDS_df = pd.DataFrame(CDS_count, columns=['Sample','CDS'])
    clusters_df = pd.DataFrame(clusters_count, columns=['Sample','Clusters'])
    hits_df = pd.DataFrame(hits_count, columns=['Sample','Known Resistence Hits'])
    core_df = pd.DataFrame(core_count, columns=['Sample','Core Genes'])
    dup_df = pd.DataFrame(dup_count, columns=['Sample','Duplicated'])
    bgc_prox_df = pd.DataFrame(bgc_prox_count, columns=['Sample','BGC Proximity'])
    
    print('Total CDSs: ', sum(CDS_df['CDS']), 
          '\nTotal clusters: ', sum(clusters_df['Clusters']),
          '\nTotal Known Resistence Hits: ', sum(hits_df['Known Resistence Hits']),
          '\nTotal Core Genes: ', sum(core_df['Core Genes']),
          '\nTotal Duplicated Genes: ', sum(dup_df['Duplicated']),
          '\nTotal BGC Proximity hits: ', sum(bgc_prox_df['BGC Proximity'],),
          '\n\nAll right. Mounting final table...')
    
    data_frames = [CDS_df, clusters_df, hits_df, core_df, dup_df, bgc_prox_df]

    overview_df = pd.DataFrame(data=samples, columns=['Sample'])
    for data_frame in data_frames:
        if data_frame.empty is False:
            overview_df = overview_df.merge(data_frame, on='Sample')
        else:
            for col in data_frame.columns:
                if col != 'Sample':
                    data_frame[col] = np.NaN
                else:
                    data_frame[col] = samples

            overview_df = overview_df.merge(data_frame, on='Sample')

    
    # data_frames = data_frames[0].join(data_frames[1:])
    
    overview_df.to_csv(os.path.join(directory,'ARTS_overview.tsv'),sep='\t', index=False)
    df2xlsx(os.path.join(directory,'ARTS_overview.xlsx'), 'ARTS_overview', overview_df)
    print('Overview Done.\n')

""" def readTSVKnownHits(tsv):
    df = pd.read_csv(tsv,sep='\t')
    df['Sample'] = str(os.path.basename(os.path.dirname(os.path.dirname(str(tsv)))))
    df.rename(columns={'#Model':'Model'},inplace=True)
    df['Contig'] = df['Sequence description'].str.split('|').str[2].str.split('_').str[1].astype(int)
    df['Contig'] = df['Sample']+'_'+df['Contig'].astype(str)
    df['HitStart'] = df['Sequence description'].str.split('|').str[6].str.split('_').str[0]
    df['HitEnd'] = df['Sequence description'].str.split('|').str[6].str.split('_').str[1]
    df['HitStrand'] = df['Sequence description'].str.split('|').str[6].str.split('_').str[2]
    df.drop(columns=['Sequence description'],inplace=True)
   
    return df """

def readTSVKnownHits(tsv):
    df = pd.read_csv(tsv, sep="\t")
    sample = os.path.basename(os.path.dirname(os.path.dirname(tsv)))
    df["Sample"] = sample
    df.rename(columns={"#Model": "Model"}, inplace=True)

    tsv_abs = os.path.abspath(tsv)
    sample_dir = os.path.dirname(os.path.dirname(tsv_abs))
    arts_dir = os.path.dirname(sample_dir)
    database_dir = os.path.dirname(arts_dir)
    rename_contigs_file = os.path.join(database_dir, "contigsFilesRenamed.tsv")

    # Add conditional structure to change how tables are constucted. If -noprep is activated it will look for contig names inside .fna files
    if os.path.exists(rename_contigs_file):
        # If -noprep is deactivated
        print("\n")
        print(f"Using renamed {rename_contigs_file} in ARTS")

        df["Contig"] = (
            df["Sequence description"]
              .str.split("|").str[2]
              .str.split("_").str[1]
              .astype(int)
        )
        df["Contig"] = df["Sample"] + "_" + df["Contig"].astype(str)
        df["HitStart"] = (
            df["Sequence description"]
              .str.split("|").str[6]
              .str.split("_").str[0]
        )
        df["HitEnd"] = (
            df["Sequence description"]
              .str.split("|").str[6]
              .str.split("_").str[1]
        )
        df["HitStrand"] = (
            df["Sequence description"]
              .str.split("|").str[6]
              .str.split("_").str[2]
        )
        df.drop(columns=["Sequence description"], inplace=True)

    else:
        # if -noprep is activated
        # Build dictionary from the sample .fna
        print("\n")
        print("contigsFilesRenamed.tsv not found during ARTS. Recovering contig names from .fna files.")
        sample_dir = pathlib.Path(tsv).parents[1]
        fna_file = sample_dir / f"{sample}.fna"


        if not fna_file.exists():
            raise FileNotFoundError(
                f"Cannot find FNA file:\n{fna_file}"
            )

        contig_dict = {}

        for record in SeqIO.parse(str(fna_file), "fasta"):
            header = record.description
            try:
                sequence_id = int(header.split("|")[1])
                remainder = header.split(" ", 1)[1]
                contig = re.split(r"[ |]", remainder, maxsplit=1)[0]
                contig_dict[sequence_id] = contig


            except Exception:
                continue

        # Add new column
        df["Contig"] = df["Sequence id"].map(contig_dict)
        # Warn about IDs not found
        missing = df["Contig"].isna().sum()
        if missing > 0:
            print(
                f"WARNING: {missing} Sequence ids could not be matched "
                f"in {fna_file.name}"
            )
        df["HitStart"] = (
            df["Sequence description"]
                .str.split("|").str[6]
                .str.split("_").str[0]
        )
        df["HitEnd"] = (
            df["Sequence description"]
                .str.split("|").str[6]
                .str.split("_").str[1]
        )
        df["HitStrand"] = (
            df["Sequence description"]
                .str.split("|").str[6]
                .str.split("_").str[2]
        )    
        df.drop(columns=["Sequence description"], inplace=True)

    return df
    
def readTSVCoreHits(tsv):
    df = pd.read_csv(tsv,sep='\t')
    df['Sample'] = str(os.path.basename(os.path.dirname(os.path.dirname(str(tsv)))))
    df.rename(columns={'#Core_gene':'Core_gene'},inplace=True)
    df['[Hits_listed]'] = df['[Hits_listed]'].str.replace('[','').str.replace(']','').str.split(';')
    df = df.explode('[Hits_listed]')
    df['Contig'] = df['[Hits_listed]'].str.split('|').str[4].str.split('=').str[-1].str.split('_').str[0].astype(int)
    df['Contig'] = df['Sample']+'_'+df['Contig'].astype(str)
    df['HitStart'] = df['[Hits_listed]'].str.split('|').str[4].str.split(' ').str[0]
    df['HitEnd'] = df['[Hits_listed]'].str.split('|').str[4].str.split(' ').str[1]
    df['HitStrand'] = df['[Hits_listed]'].str.split('|').str[4].str.split(' ').str[2]
    df.drop(columns=['[Hits_listed]'],inplace=True)
    return df
    
def readTSVdupTable(tsv):
    df = pd.read_csv(tsv,sep='\t')
    df['Sample'] = str(os.path.basename(os.path.dirname(os.path.dirname(str(tsv)))))
    return df

def ARTS_Results_Extraction(path):
    path = os.path.join(path,'ARTS')
    os.chdir(path)
    directory=os.path.join(path,'ARTS_Extractor/')

    number_of_samples = 0
    for item in os.listdir(path):
        if os.path.isdir(item) and 'ARTS_Extractor' not in item:
            number_of_samples += 1

    print('>>> CORE GENES INFORMATION <<<\n')
    print('Getting Core Genes info...\n')

    # knownResistenceHits = None
    if glob.glob('**/**/knownhits.tsv') != []:
        knownResistenceHits = pd.concat(map(readTSVKnownHits, glob.glob('**/**/knownhits.tsv')))
        knownResistenceHits.to_csv(os.path.join(directory,'KnownHits.tsv'),sep='\t',index=False)
        df2xlsx(os.path.join(directory,'KnownHits.xlsx'), 'KnownHits', knownResistenceHits)

    # else:
    #     print(number_of_samples)
    #     knownResistenceHits = pd.DataFrame(np.nan, index=[i for i in range(0, number_of_samples)], columns=['Model', 'Description', 'Sequence', 'id', 'evalue' , 'bitscore', 'Sample', 'Contig', 'Contig', 'HitStart', 'HitEnd', 'HitStrand'])
    
    CoreHits = pd.concat(map(readTSVCoreHits, glob.glob('**/**/coretable.tsv')))
    CoreHits.to_csv(os.path.join(directory,'CoreHits.tsv'),sep='\t',index=False)
    df2xlsx(os.path.join(directory,'CoreHits.xlsx'), 'CoreHits', CoreHits)
    
    dupTable = pd.concat(map(readTSVdupTable, glob.glob('**/**/duptable.tsv')))
    dupTable.to_csv(os.path.join(directory,'DupHits.tsv'),sep='\t',index=False)
    df2xlsx(os.path.join(directory,'DupHits.xlsx'), 'DupHits', dupTable)

if __name__ == '__main__':
    # prepare_extraction(os.getcwd())
    ARTS_overview('/media/bioinfo/6tb_hdd/03_ELLEN/krill_runs/NCBI_PROJECTS/AtlanticoSul/')
    # ARTS_Results_Extraction(os.getcwd())
  