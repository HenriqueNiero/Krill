#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 19:24:23 2021

@author: saulo
"""
import numpy as np
import pandas as pd
import sys, os, ast
from statistics import mean
import pathlib, os
import ast, numpy as np


def df2xlsx(path, sheet_name, df):
    df = df.copy().reset_index()
    writer = pd.ExcelWriter(path, engine='xlsxwriter')

    df.to_excel(writer, sheet_name=sheet_name, startrow=1, header=False, index=False)
    worksheet = writer.sheets[sheet_name]

    (max_row, max_col) = df.shape
    column_settings = [{'header': column} for column in df.columns]
    worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings})
    worksheet.set_column(0, max_col - 1, 15)

    writer._save()


def map_products_to_category(product):
    replacement_rules = {
        "PKS I": ["t1pks", "T1PKS"],
        "PKS other": ["transatpks", "t2pks", "t3pks", "otherks", "hglks", "transAT-PKS", "transAT-PKS-like", "T2PKS", "T3PKS", "PKS-like", "hglE-KS", "prodigiosin"],
        "NRPS": ["nrps", "NRPS", "NRPS-like", "thioamide-NRP", "NAPAA"],
        "RiPPs": ["lantipeptide", "thiopeptide", "bacteriocin", "linaridin", "cyanobactin", "glycocin", "LAP", "lassopeptide", "sactipeptide", "bottromycin", "head_to_tail", "microcin", "microviridin", "proteusin", "lanthipeptide", "lipolanthine", "RaS-RiPP", "fungal-RiPP","fungalRiPP", "TfuA-related", "guanidinotides", "RiPP-like", "lanthipeptide-class-iii","lanthipeptide-class-i", "lanthipeptide-class-ii","lanthipeptide-class-iv", "lanthipeptide-class-v", "redox-cofactor", "thioamitides", "ranthipeptide",  "epipeptide", "cyclic-lactone-autoinducer", "spliceotide", "RRE-containing", "crocagin"],
        "Saccharides": ["amglyccycl", "oligosaccharide", "cf_saccharide", "saccharide"],
        "Terpene": "terpene",
        "Others": ["acyl_amino_acids", "others", "arylpolyene", "aminocoumarin", "ectoine", "butyrolactone", "nucleoside", "melanin", "phosphoglycolipid", "phenazine", "phosphonate", "other", "cf_putative", "resorcinol", "indole", "ladderane", "PUFA", "furan", "hserlactone", "fused", "cf_fatty_acid", "siderophore", "blactam", "fatty_acid", "PpyS-KS", "CDPS", "betalactone", "PBDE", "tropodithietic-acid", "NAGGN", "halogenated", "pyrrolidine", "mycosporine-like"]
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


def get_csvs(path,pattern):
    return list(pathlib.Path(path).glob('**/{}'.format(pattern)))

def run(path):
    if os.path.exists(os.path.join(path,'ARTS','ARTS_Extractor','KnownHits.tsv')):
        KnownHits = pd.read_csv(os.path.join(path,'ARTS','ARTS_Extractor','KnownHits.tsv'),sep='\t',converters={'Sample':str})

    CoreHits = pd.read_csv(os.path.join(path,'ARTS','ARTS_Extractor','CoreHits.tsv'),sep='\t',converters={'Sample':str})
    DupHits = pd.read_csv(os.path.join(path,'ARTS','ARTS_Extractor','DupHits.tsv'),sep='\t',converters={'Sample':str})
    
    try:
        clusters_blast = pd.read_csv(os.path.join(path,'AntiSMASH','AntiSMASH_Extractor','clusters_blast.tsv'),sep='\t')
    except pd.errors.EmptyDataError:
        clusters_blast = pd.DataFrame()
    clusters = pd.read_csv(os.path.join(path,'AntiSMASH','AntiSMASH_Extractor','clusters.tsv'),sep='\t',converters={'file_name':str})


    hit_selected = clusters[['file_name','contig','type','cluster_number','product','contig_edge','kind','protoclusters','SMILES','Start','End','strand','genes']].set_index('contig')
    hit_selected['KnownResistenceHit'] = ''
    hit_selected['BGCs_Hits'] = ''
    hit_selected['BGCs_Hits_Mean_Similarities(%)'] = ''


    for i in clusters.index:
        bgc = clusters.loc[i, 'contig']
        cluster_number = clusters.loc[i,'cluster_number']
        start = int(clusters.loc[i,'Start'])
        end = int(clusters.loc[i,'End'])
        
        if os.path.exists(os.path.join(path,'ARTS','ARTS_Extractor','KnownHits.tsv')):
            tmp = KnownHits[(KnownHits['Contig'] == bgc) & (KnownHits['HitStart'].astype(int) >= start) & (KnownHits['HitEnd'].astype(int) <= end)]
            if not tmp.empty:
                row_list = []
                for row in tmp.itertuples():
                    row_list.append(['KnownResistenceHit',row.Model,row.Description,row.evalue,row.HitStart,row.HitEnd,row.HitStrand])
                hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.Start >= start) & (hit_selected.End <= end)),'KnownResistenceHit'] = str(row_list)
            else:
                hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.Start >= start) & (hit_selected.End <= end)),'KnownResistenceHit'] = 'N/A'
        
        if not clusters_blast.empty:
            tmp_blast = clusters_blast[(clusters_blast['contig'] == bgc) & (clusters_blast['cluster'] == cluster_number)]
        else:
            tmp_blast = pd.DataFrame()
            
        if not tmp_blast.empty:
            row_list = []
            mean_similarities = []
            for row in tmp_blast.itertuples():
                row_list.append(['BGCs_Hits',row.accession,row.cluster_type,row.description,row.similarity])
                mean_similarities.append(row.similarity)
            hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.cluster_number == cluster_number)),'BGCs_Hits'] = str(row_list)
            hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.cluster_number == cluster_number)),'BGCs_Hits_Mean_Similarities(%)'] = round(mean(mean_similarities),2)
        else:
            hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.cluster_number == cluster_number)),'BGCs_Hits'] = 'N/A'
            hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.cluster_number == cluster_number)),'BGCs_Hits_Mean_Similarities(%)'] = 'N/A'

    hit_selected['CoreHit'] = ''
    for i in clusters.index:
        bgc = clusters.loc[i, 'contig']
        start = int(clusters.loc[i,'Start'])
        end = int(clusters.loc[i,'End'])
        tmp = CoreHits[(CoreHits['Contig'] == bgc) & (CoreHits['HitStart'].astype(int) >= start) & (CoreHits['HitEnd'].astype(int) <= end)]
        if not tmp.empty:
            row_list = []
            for row in tmp.itertuples():
                row_list.append(['CoreHit',row.Core_gene,row.Description,row.Function,row.HitStart,row.HitEnd,row.HitStrand])
            hit_selected.loc[((hit_selected.index == bgc) & (hit_selected.Start >= start) & (hit_selected.End <= end)),'CoreHit'] = str(row_list)
        else:
            hit_selected.loc[bgc,'CoreHit'] = 'N/A'
    
    hit_selected = hit_selected.replace('N/A',np.NaN).replace('',np.NaN)
    #hit_selected.dropna(subset=['KnownResistenceHit','CoreHit'],how='all',inplace=True)
    
    
    hit_selected.rename(columns={'contig_edge':'completeness','file_name':'sample'},inplace=True)
    #comment rename contig_edge to completeness and true to fragmented and false to complete
    completeness = {"['True']":"Fragmented","['False']":"Complete"}
    hit_selected['completeness'] = hit_selected['completeness'].replace(completeness)
    hit_selected['Size'] = hit_selected['End']-hit_selected['Start']

    hit_selected['regulatory_genes'] = hit_selected['genes'].str.contains('regulatory',case=False)
    hit_selected = hit_selected[['sample','cluster_number', 'product', 'completeness', 'kind', 'protoclusters', 'SMILES', 'Start', 'End', 'Size', 'strand', 'genes', 'regulatory_genes','KnownResistenceHit','CoreHit','BGCs_Hits','BGCs_Hits_Mean_Similarities(%)']]

    hit_selected['product_bigscape'] = hit_selected['product'].apply(map_products_to_category)
    hit_selected.insert(hit_selected.columns.get_loc('product') + 1, 'product_bigscape', hit_selected.pop('product_bigscape'))


    #Add column original_contig_header to the BGCs_with_Hits.tsv table


    ref_renamecontig = get_csvs(path,'contigsFilesRenamed.tsv')

    renamecontig = {}

    for ref in ref_renamecontig:
        db = os.path.basename(os.path.dirname(ref))
        df = pd.read_csv(ref,dtype='object',sep='\t')
        df = df.apply(lambda col: col.map(lambda x: x.rstrip(f'.fasta')))
        df.rename(columns={'new_contigs':'contig','original_contigs':'Original_contig'},inplace=True)
        hit_selected_originalcontigs = pd.merge(hit_selected, df, how="left", on="contig")
        




    # If the contigs files rename folder exist it will use this file to generate the tables, otherwise
    # it will be skipped and the original hit_selected will be used to generate the tables 
    if os.path.exists(os.path.join(path, 'contigsFilesRenamed.tsv')):
        hit_selected_originalcontigs.to_csv(os.path.join(path,'ReportOutput/BGCs_with_Hits.tsv'),sep='\t')
    
        # hit_selected_originalcontigs.to_csv(os.path.join(path,'BGCs_with_Hits.tsv'),sep='\t')
        df2xlsx(os.path.join(path,'ReportOutput/BGCs_with_Hits.xlsx'), 'BGCs with Hits', hit_selected_originalcontigs)

        hit_selected_originalcontigs = hit_selected_originalcontigs[hit_selected_originalcontigs['completeness'] == "Complete"]
        hit_selected_originalcontigs.to_csv(os.path.join(path,'ReportOutput/Complete_BGCs_with_Hits.tsv'),sep='\t')
        # hit_selected.to_csv(os.path.join(path,'Complete_BGCs_with_Hits.tsv'),sep='\t')
        df2xlsx(os.path.join(path,'ReportOutput/Complete_BGCs_with_Hits.xlsx'), 'Complete BGCs with Hits', hit_selected_originalcontigs)

    else:
        hit_selected.to_csv(os.path.join(path,'ReportOutput/BGCs_with_Hits.tsv'),sep='\t')
    
        # hit_selected.to_csv(os.path.join(path,'BGCs_with_Hits.tsv'),sep='\t')
        df2xlsx(os.path.join(path,'ReportOutput/BGCs_with_Hits.xlsx'), 'BGCs with Hits', hit_selected)

        hit_selected = hit_selected[hit_selected['completeness'] == "Complete"]
        hit_selected.to_csv(os.path.join(path,'ReportOutput/Complete_BGCs_with_Hits.tsv'),sep='\t')
        # hit_selected.to_csv(os.path.join(path,'Complete_BGCs_with_Hits.tsv'),sep='\t')
        df2xlsx(os.path.join(path,'ReportOutput/Complete_BGCs_with_Hits.xlsx'), 'Complete BGCs with Hits', hit_selected)

if __name__ == '__main__':
    run('/media/bioinfo/6tb_hdd/03_ELLEN/krill_runs/01_REPORT_DATABASES/KRILL_RESULTS/BlackSeaSequences')