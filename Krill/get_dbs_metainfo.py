#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:26:42 2021

@author: saulo
"""
import pandas as pd
import pathlib, os
import ast, numpy as np


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

def get_csvs(path,pattern):
    return list(pathlib.Path(path).glob('**/{}'.format(pattern)))

def map_products_to_category(product):
    replacement_rules = {
        "PKS I": ["t1pks", "T1PKS"],
        "PKS other": ["transatpks", "t2pks", "t3pks", "otherks", "hglks", "transAT-PKS", "transAT-PKS-like", "T2PKS", "T3PKS", "PKS-like", "hglE-KS", "prodigiosin"],
        "NRPS": ["nrps", "NRPS", "NRPS-like", "thioamide-NRP", "NAPAA"],
        "RiPPs": ["lantipeptide", "thiopeptide", "bacteriocin", "linaridin", "cyanobactin", "glycocin", "LAP", "lassopeptide", "sactipeptide", "bottromycin", "head_to_tail", "microcin", "microviridin", "proteusin", "lanthipeptide", "lipolanthine", "RaS-RiPP", "fungal-RiPP","fungalRiPP", "TfuA-related", "guanidinotides", "RiPP-like", "lanthipeptide-class-iii","lanthipeptide-class-i", "lanthipeptide-class-ii","lanthipeptide-class-iv", "lanthipeptide-class-v", "redox-cofactor", "thioamitides", "ranthipeptide",  "epipeptide", "cyclic-lactone-autoinducer", "spliceotide", "RRE-containing", "crocagin"],
        "Saccharides": ["amglyccycl", "oligosaccharide", "cf_saccharide", "saccharide"],
        "Terpene": "terpene",
        "Others": ["acyl_amino_acids", "arylpolyene", "aminocoumarin", "ectoine", "butyrolactone", "nucleoside", "melanin", "phosphoglycolipid", "phenazine", "phosphonate", "other", "cf_putative", "resorcinol", "indole", "ladderane", "PUFA", "furan", "hserlactone", "fused", "cf_fatty_acid", "siderophore", "blactam", "fatty_acid", "PpyS-KS", "CDPS", "betalactone", "PBDE", "tropodithietic-acid", "NAGGN", "halogenated", "pyrrolidine", "mycosporine-like"]
    }

    for category, replacements in replacement_rules.items():
        if product in replacements:
            return category
        elif ',' in product:
            hybrids = product.split(',')
            count = 0
            count_nrps = 0
            count_pksI = 0
            for hybrid in hybrids:
                if hybrid in replacement_rules['NRPS']:
                    count_nrps += 1
                elif hybrid in replacement_rules["PKS I"]:
                    count_pksI += 1
                elif hybrid in replacements:
                    count += 1
                elif category == "PKS other" and hybrid in replacement_rules["PKS I"]:
                    count += 1

            if count_nrps + count_pksI == len(hybrids):
                return "PKS/NRPS Hybrids"
            elif count == len(hybrids):
                return category

    return 'others'


def get(path, ext, root_database):
    ref_rename = get_csvs(path,'fastaFilesRenamed.tsv')
    rename = {}
    
    for ref in ref_rename:
        db = os.path.basename(os.path.dirname(ref))
        df = pd.read_csv(ref,dtype='object',sep='\t')
        df = df.apply(lambda col: col.map(lambda x: x.rstrip(f'.{ext}')))
        rename[db] = pd.Series(df.OriginalName.values, index=df.NewName.values).to_dict()
        
    files = ['BGCs_with_Hits.tsv','Complete_BGCs_with_Hits.tsv','DNABases_and_ORFs_count.tsv']
    
    for f in files:
        tables = get_csvs(path,f)
        print("\n \n tables = get_csv(path,f) for f in files: \n", tables)
        table_clean = pd.DataFrame()
        countORFsAndDNA = pd.DataFrame()
        
        for t in tables:
            # db = os.path.basename(os.path.dirname(t))
            db = t.parents[1].name
            # Create df with BGCs_with_Hits.tsv information
            df = pd.read_csv(t,dtype='object',sep='\t')
            # Add Database name column to the df for DBs_BGCs_with_Hits.tsv table in DBsReportOutput
            df['Database'] = db
            print("\n \n dataframe df depois de pd.read_csv das tabelas de entrada: \n", df)



            if not 'DNABases' in str(t):
                df['OriginalName'] = df['contig'].str.split('_').str[0].replace(rename[db],regex=True)
                table_clean = pd.concat([table_clean, df], ignore_index=True)
                print("\n \n table_clean depois de table_clean = pd.concat([table_clean, df], ignore_index=True): \n", table_clean)

            if 'DNABases' in str(t):
                df['NT (KB)'] = df['NT (KB)'].astype(float)
                df[['CONTIGS', 'ORFS']] = df[['CONTIGS', 'ORFS']].astype(int)
                df[['FILE SIZE (MB)', 'MIN LEN (KB)','MAX LEN (KB)', 'AVG LEN (KB)']] = df[['FILE SIZE (MB)', 'MIN LEN (KB)','MAX LEN (KB)', 'AVG LEN (KB)']].astype(float)
                # df = df.groupby('Database').sum()
                df = df.groupby('Database').agg({'FILE SIZE (MB)': 'sum', 'NT (KB)': 'sum', 'ORFS': 'sum', 'CONTIGS': 'sum', 'MIN LEN (KB)': 'min', 'MAX LEN (KB)': 'max', 'AVG LEN (KB)': 'mean'})
                df = df.reset_index()
                if root_database is False:
                    bgcs_count = len(list(pathlib.Path(os.path.join(path,db)).rglob('*region*.gbk')))
                else:
                    bgcs_count = len(list(pathlib.Path(path).rglob('*region*.gbk')))
                nt = df['NT (KB)']*1000
                df['BGCs/MegaBases'] = bgcs_count/nt*1000000
                df['BGCs/MegaORFs'] = bgcs_count/df['ORFS']*1000000
                # countORFsAndDNA = countORFsAndDNA.append(df)

                countORFsAndDNA = pd.concat([countORFsAndDNA, df], ignore_index=True)

        # Inseri as informações de kind, protocluster e Original_contig na tabela        
        if not 'DNABases' in str(f):
            table_clean = table_clean[["Database","OriginalName","Original_contig","contig","cluster_number","product","kind","protoclusters","completeness","SMILES","Start","End","Size","strand","genes","regulatory_genes",'KnownResistenceHit','CoreHit','BGCs_Hits','BGCs_Hits_Mean_Similarities(%)']].sort_values(by='Database')
            print("\n \n table_clean depois de table_clean = table_clean[[Database, OriginalName, contig, cluster_number, smiles, start, end, etc...]]: \n", table_clean)
            table_clean['product_bigscape'] = table_clean['product'].apply(map_products_to_category)
            table_clean.insert(table_clean.columns.get_loc('product') + 1, 'product_bigscape', table_clean.pop('product_bigscape'))
            print("\n \n table_clean depois de table_clean.insert product_bigscape depois de product: \n", table_clean)
            table_clean.to_csv(os.path.join(path,'DBsReportOutput/DBs_'+f),index=False,sep='\t')
            f = f.replace('tsv', 'xlsx')
            df2xlsx(os.path.join(path,'DBsReportOutput/DBs_'+f), 'DBs Report', table_clean)

        if 'DNABases' in str(f):
            countORFsAndDNA.to_csv(os.path.join(path,'DBsReportOutput/DBs_normalized_info.tsv'),sep='\t')
            df2xlsx(os.path.join(path,'DBsReportOutput/DBs_normalized_info.xlsx'), 'DBs normalized info', countORFsAndDNA)

if __name__ == '__main__':
    get('/media/bioinfo/6tb_hdd/03_ELLEN/krill_runs/NCBI_PROJECTS/','.fasta')