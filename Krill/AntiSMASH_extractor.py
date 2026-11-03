#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 08:09:19 2021

@author: saulo
"""
from Bio import SeqIO
from tqdm import tqdm
import pathlib, os, sys, json
import pandas as pd
from cprint import *
from pandas import json_normalize
from concurrent.futures import ThreadPoolExecutor, as_completed


def prepare_extraction(path):
    directory = os.path.join(path,'AntiSMASH','AntiSMASH_Extractor/')
    os.makedirs(directory, exist_ok=True)


'''
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
'''

def df2xlsx(path, sheet_name, df):
    """Excel export disabled; TSV/CSV exports remain enabled."""
    print(f"Skipping Excel export: {sheet_name}")
    return

def parse_gbk_cand_cluster(file):

  clusters = []

  base = os.path.basename(file)

  for seq_record in SeqIO.parse(file, "genbank"):
          seqtitle = seq_record.id

          for seq_feature in seq_record.features:

              if 'cand_cluster' in seq_feature.type.lower():
                  note = seq_feature.qualifiers.get('candidate_cluster_number', [])
                  feat_type = seq_feature.type
                  feat_product = seq_feature.qualifiers.get('product', [])
                  edge = seq_feature.qualifiers.get('contig_edge', [])
                  kind = seq_feature.qualifiers.get('kind', [])
                  protos = seq_feature.qualifiers.get('protoclusters', [])
                  SMILES = seq_feature.qualifiers.get('SMILES', [])
                  bpend = int(seq_feature.location.end)
                  bpstart = int(seq_feature.location.start)
                  feat_strand = seq_feature.strand
                  polymer = seq_feature.qualifiers.get('polymer', [])
                  clusters.append([base[:-4], seqtitle, feat_type, note[0], str(",".join(feat_product)), edge, kind, protos, SMILES, bpstart, bpend, polymer, feat_strand])
  return clusters

              
def parse_gbk_cds(file):

  cds_list = []

  base = os.path.basename(file)

  for seq_record in SeqIO.parse(file, "genbank"):
          seqtitle = seq_record.id

          for seq_feature in seq_record.features:

              if 'CDS' == seq_feature.type:

                  feat_type = seq_feature.type
                  loc_start = int(seq_feature.location.start)
                  loc_end = int(seq_feature.location.end)
                  domain = seq_feature.qualifiers.get('NRPS_PKS', [])
                  functions = seq_feature.qualifiers.get('gene_functions', [])
                  kind = seq_feature.qualifiers.get('gene_kind', [])
                  locus_tag = seq_feature.qualifiers.get('locus_tag', [])
                  sec = seq_feature.qualifiers.get('sec_met_domain', [])
                  if functions != []:
                    cds_list.append([base[:-4], seqtitle, feat_type, loc_start, loc_end, domain, functions, kind, locus_tag, sec])
  return cds_list


def protocluster_parse(path, threads):

    path = os.path.join(path,'AntiSMASH')

    directory=os.path.join(path,'AntiSMASH_Extractor/')
    
    # Get list of gbk files (one for each fasta input file)
    reference = [x for x in list(pathlib.Path(path).glob("**/*.gbk")) if not 'region' in str(x)]
    
    candclusters = []  # Candidate clusters
    genes = [] 
    print('>>> CLUSTERS PARSE <<<\n')

    print('Parsing genbank files...')

    db = os.path.basename(path)  # Get database name
    
    # Display a visual progress bar for the region parse
    with tqdm(desc='{db} Region parse'.format(db=os.path.basename(db)),total=len(reference),unit='gbk',colour='red',position=1,leave=False) as pbar:
      # Use multiple threads for the execution of the gbk files
      with ThreadPoolExecutor(max_workers=threads) as executor:
          # Parse gbk candidate clusters
          result_futures = [executor.submit(parse_gbk_cand_cluster, file) for file in reference]

          for future in as_completed(result_futures):
            try:
              if future.result() != None:
                # Add result to the candidate clusters list
                candclusters = candclusters + future.result()
              # Update progress bar
              pbar.update(1)
            except Exception as e:
              cprint.fatal(e,interrupt=False)


    # Display a visual progress bar for the coding sequences (CDS) parse
    with tqdm(desc='{db} CDS parse'.format(db=os.path.basename(db)),total=len(reference),unit='gbk',colour='red',position=1,leave=False) as pbar:
      # Use multiple threads for the execution of the gbk files
      with ThreadPoolExecutor(max_workers=threads) as executor:
          # Parse gbk coding sequences
          result_futures = [executor.submit(parse_gbk_cds, file) for file in reference]
          for future in as_completed(result_futures):
            try:
              if future.result() != None:
                # Add result to the genes list
                genes = genes + future.result()
              # Update progress bar
              pbar.update(1)
            except Exception as e:
              cprint.fatal(e,interrupt=False)

# ==============================================================================
    # Merge candidate clusters and genes (coding sequences) into a dataframe
    print('Building dataframe...')

    candclusters_df = pd.DataFrame(candclusters, columns=['file_name','contig','type', 'cluster_number', 'product', 'contig_edge', 'kind', 'protoclusters', 'SMILES', 'Start', 'End', 'polymer', 'strand'])
    candclusters_df.sort_values(by=['file_name','contig', 'cluster_number'], inplace=True)
    candclusters_df['genes'] = ''
    candclusters_df['genes'].astype(str)
    
    print('Placing ',len(genes),' genes into the right cluster...')
    
    for i in genes:
      for index in candclusters_df.index:
        if candclusters_df.loc[index,'file_name'] == str(i[0]) and candclusters_df.loc[index,'contig'] == str(i[1]) and candclusters_df.loc[index,'Start'] <= float(i[3]) and candclusters_df.loc[index,'End'] >= float(i[4]):
            candclusters_df.loc[index,'genes'] += str(i[1:10]) + ','
    
    print('Saving final dataframe into CSV...')
    
    candclusters_df.to_csv(directory+'clusters.tsv', sep='\t', index=False)
    df2xlsx(directory+'clusters.xlsx', 'clusters', candclusters_df)
# ==============================================================================


       
def parse_blast_results_1(jfile):

        with open(jfile) as json_file:
            data = json.load(json_file)
    
        #json normalize -> busca a sub lista FEATURES, dentro de Records
        #Mantendo as informações input_file, taxo, e o record_id dentro da guia (records, modules, antismash.modules.nrps_pks, record_id)
        return json_normalize(data,
                                    record_path=["records",
                                                 "features"],
        
                                    meta=["input_file",
                                          "taxon",
                                          ["records", "modules", "antismash.modules.nrps_pks", "record_id"]],
                                    meta_prefix="rec_",
                                    errors="ignore")
    
def parse_blast_results_2(jfile):

        with open(jfile) as json_file:
            data = json.load(json_file) 

 
        #Lendo a subtabelas       
        for d in data["records"]:
            
            #Smiles
            try:
                #Criando a tabela de path de busca para as várias possíveis chaves
                for k in d["modules"]["antismash.modules.nrps_pks"]["region_predictions"].keys():
                    return json_normalize(d["modules"],
                                               record_path=["antismash.modules.nrps_pks",
                                                            "region_predictions",
                                                            ]+[k],
                                               meta=[["antismash.modules.nrps_pks", "record_id"]],
                                               meta_prefix="mod_",
                                               errors="ignore")
                    
            except Exception as e:
                pass
            
def parse_blast_results_3(jfile):

        with open(jfile) as json_file:
            data = json.load(json_file) 


        #Lendo a subtabelas       
        for d in data["records"]:
            # Ranking Apenas
            i_file = data["input_file"]
            try:
                #json normalize -> busca a sub lista proteins, dentro de ranking
                #Mas tira a coluna 0 para ficar apenas com os dados da tabela acima
                rkg = json_normalize(d,
                                      record_path=["modules",
                                                  "antismash.modules.clusterblast",
                                                  "knowncluster",
                                                  "results",
                                                  "ranking"],
                                      meta=[["modules",
                                              "antismash.modules.clusterblast",
                                              "knowncluster",
                                              "record_id"],
                                            ["modules",
                                             "antismash.modules.clusterblast",
                                             "knowncluster",
                                             "results",
                                             "region_number"]],
                                      meta_prefix="rkg_",
                                      errors="ignore")  
                rkg["input_file"] = i_file
                return rkg
            except Exception as e:
                pass


def get_clusters_blast(path, threads):

    path = os.path.join(path,'AntiSMASH')

    directory=os.path.join(path,'AntiSMASH_Extractor/')

        #legacy BLAST output format 9 file
    print('\n>>> CLUSTERS BLAST PARSE <<<\n')

    # Initialize dataframes
    recDf = pd.DataFrame()
    modDf = pd.DataFrame()
    rkgDf = pd.DataFrame()
    
    ref = list(pathlib.Path(path).glob("**/*.json"))

    db = os.path.basename(path)
    
    functions = [parse_blast_results_1, parse_blast_results_2, parse_blast_results_3]

    count = 1
    for function in functions:
       with tqdm(desc='BLAST results parse part {db}/3'.format(db=count),total=len(ref),unit='gbk',colour='red',position=1,leave=False) as pbar:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            result_futures = [executor.submit(function, file) for file in ref]
            for future in as_completed(result_futures):

              if function == parse_blast_results_1:
                try:
                  recDf = pd.concat([recDf, future.result()])
                  pbar.update(1)
                except Exception as e:
                  cprint.fatal(e,interrupt=False)
              elif function == parse_blast_results_2:
                try:
                  modDf = pd.concat([modDf, future.result()])
                  pbar.update(1)
                except Exception as e:
                  cprint.fatal(e,interrupt=False)
              elif function == parse_blast_results_3:
                try:
                  rkgDf = pd.concat([rkgDf, future.result()])
                  pbar.update(1)
                except Exception as e:
                  cprint.fatal(e,interrupt=False)
       count+=1
    rkgDf = rkgDf.rename(columns={"rkg_modules.antismash.modules.clusterblast.knowncluster.record_id": "record_id"})
    rkgDf = rkgDf.rename(columns={"rkg_modules.antismash.modules.clusterblast.knowncluster.results.region_number": "region_number"})



    rkgDf.reset_index(drop=True, inplace=True)
    
    rkgDfFinal = pd.DataFrame()
    
    rkg_cols = ["input_file",
                "record_id",
                "accession",
                "region_number",
                "cluster_label",
                "cluster_type",
                "description",
                "blast_score",
                "core_bonus",
                "core_gene_hits",
                "hits",
                "bgc_proteins",
                "synteny_score"]

    for idx in rkgDf.index:
        i_file = rkgDf["input_file"][idx]
        r_id = rkgDf["record_id"][idx]
        aux = pd.concat([json_normalize(rkgDf[0][idx]),
                         json_normalize(rkgDf[1][idx])],
                        axis=1)
        
        
        tags_list = aux["tags"].tolist()
        tags_list2 = tags_list[0]
        tags_len = len(tags_list2)
        reg_num = rkgDf["region_number"][idx]
        
        aux["input_file"] = i_file
        aux["record_id"] = r_id
        aux["bgc_proteins"] = tags_len
        aux["region_number"] = reg_num
        
        rkgDfFinal = pd.concat([rkgDfFinal, aux[rkg_cols]])
    
    rkgDfFinal.reset_index(drop=True, inplace=True)

    if rkgDfFinal.empty is False:
      rkgDfFinal["similarity"] = rkgDfFinal["hits"]*100/rkgDfFinal["bgc_proteins"]
      rkgDfFinal["similarity"] = rkgDfFinal["similarity"].round(decimals=0)
      rkgDfFinal = rkgDfFinal[['input_file','record_id','region_number','accession',
                               'cluster_type','description','blast_score','hits',
                               'bgc_proteins','similarity']]
      rkgDfFinal = rkgDfFinal.rename(columns={'input_file': 'ref','record_id':'contig',
                                  'region_number':'cluster'})
      
      rkgDfFinal = rkgDfFinal.sort_values(by='similarity', ascending=False)
      rkgDfFinal = rkgDfFinal.drop_duplicates()
      
    rkgDfFinal.to_csv(os.path.join(directory,'clusters_blast.tsv'), sep='\t', index=False)
    df2xlsx(os.path.join(directory,'clusters_blast.xlsx'), 'clusters_blast', rkgDfFinal)


if __name__ == '__main__':
  path = '/media/bioinfo/6tb_hdd/03_ELLEN/krill_runs/SRR13515398/'
  threads = 2
  get_clusters_blast(path, threads)