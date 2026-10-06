import sys, os, pathlib, argparse, datetime, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from cprint import *

import prepareFastas
import ARTS_extractor
import AntiSMASH_extractor
import get_screening_results
import count_bases_and_ORFs
import run_AntiSMASH_and_ARTS
import get_dbs_metainfo
import build_charts
import run_bigscape
import run_bigscape_analysis
from bigscape_integration import integrate_bigscape_with_krill

default_threads = len(os.sched_getaffinity(0))

cprint.ok('#################\n### Krill 1.0 ###\n#################')
cprint.ok('\nPlease cite us!\n')

parser = argparse.ArgumentParser()
parser.add_argument('-noprep','--do_not_prepare_fasta_files',help='Rename fasta files, its headers and store changes in a CSV file for control [DEFAULT: TRUE]',action='store_true')
parser.add_argument('-t','--threads',help='Threads to use in analysis [DEFAULT: {}]'.format(default_threads),type=int,default=default_threads)
parser.add_argument("--bigscape", action="store_true", help="Run BiG-SCAPE2 after antiSMASH")
parser.add_argument("--bigscape_cutoff", default="0.3", help="BiG-SCAPE cutoff (default 0.3)")
parser.add_argument('PATH',help='Working path with fasta files',type=str)
parser.add_argument("--pfam-path", dest="pfam", default=None, help="Path to Pfam-A.hmm")
parser.add_argument("--bigscape-analysis", action="store_true", default=False, help="Run the optional BiG-SCAPE category analysis on an existing BiG-SCAPE database")
parser.add_argument("--bigscape-analysis-threshold", type=float, default=0.3, help="Distance threshold used for BiG-SCAPE singleton classification [DEFAULT: 0.3]")
parser.add_argument("--bigscape-analysis-permutations", type=int, default=999, help="Number of permutations used for PERMANOVA/PERMDISP [DEFAULT: 999]")


args = parser.parse_args()

# Fixed BiG-SCAPE settings
args.bigscape_env = "bigscape"
args.bigscape_mix = True
args.bigscape_classify = "category"
args.bigscape_include_singletons = True



root_path = os.path.abspath(args.PATH)
root_database = False
dbs = []
# Find database directories inside the input directory.
# Files such as taxonomy_krill are intentionally ignored.
for x in os.listdir(root_path):
    full_path = os.path.join(root_path, x)

    if os.path.isdir(full_path):
        dbs.append(os.path.abspath(full_path))

# If no database directories exist, the input directory itself is a single database.
if not dbs:
    root_database = True
    dbs.append(root_path)

args.PATH = root_path



# Convert all files extensions (fa, fna, etc) to "fasta"
prepareFastas.convert2fasta(args.PATH)


# Calculate total size of data to be analyzed
fastas=list(pathlib.Path(args.PATH).rglob('*.{}'.format("fasta")))
total_size=int(sum([os.stat(x).st_size for x in fastas])/(10**6))


# Display a visual progress bar for the processing of the whole execution
with tqdm(desc='Total Databases Analysis progress in MB',total=total_size,unit='MB',colour='blue',position=0,leave=False) as dbar:
  for db in dbs:  # Process each database/directory

      os.chdir(db)
      
      if not args.do_not_prepare_fasta_files:
          prepareFastas.run(db,"fasta")  # Rename fasta files, its headers and store changes in a CSV file for control

      os.makedirs(os.path.join(db,'AntiSMASH'),exist_ok=True)  # Create ARTS folder
      os.makedirs(os.path.join(db,'ARTS'),exist_ok=True)  # Create ARTS folder
  
      os.makedirs(os.path.join(db,'ReportOutput'),exist_ok=True)  # Create ARTS folder
    
      db_fastas = list(pathlib.Path(db).glob('*.{}'.format("fasta")))  # List of fasta files

      # Display a visual progress bar for the processing of the current database
      with tqdm(desc='{db} BGCs Analysis'.format(db=os.path.basename(db)),total=len(db_fastas),unit='file',colour='blue',position=1,leave=False) as pbar:
       # Use multiple threads for the execution of the fasta files
       with ThreadPoolExecutor(max_workers=args.threads) as executor:
         # Run AntiSMASH and ARTS
         result_futures = [executor.submit(run_AntiSMASH_and_ARTS.run, fasta_file, db) for fasta_file in db_fastas]
         # Update progress bar
         for future in as_completed(result_futures):
           try:
             pbar.update(1)
             dbar.update(int(future.result()))
           except Exception as e:
             cprint.fatal(e,interrupt=False)

      cprint.info('# Extracting AntiSMASH and ARTS results...')
      AntiSMASH_extractor.prepare_extraction(db)  # Create a folder for the AntiSMASH extraction
      ARTS_extractor.prepare_extraction(db) # Create a folder for the ARTS extraction
      
      AntiSMASH_extractor.protocluster_parse(db, args.threads)  # Parse candidate clusters and coding sequences. output: clusters.tsv
      AntiSMASH_extractor.get_clusters_blast(db, args.threads)
      
      ARTS_extractor.ARTS_overview(db)
      ARTS_extractor.ARTS_Results_Extraction(db)
    
      cprint.info('# Screening AntiSMASH and ARTS results...')
      get_screening_results.run(db)

      cprint.info('# Getting fasta info (Bases and ORFs count)...')
      count_bases_and_ORFs.run(db, "fasta", args.threads)
      
      os.chdir(args.PATH)


os.makedirs(os.path.join(args.PATH,'DBsReportOutput'),exist_ok=True)  # Create ARTS folder

cprint.info('# Getting metainfo (Databases MBases, ORFs, BGCs/Mbases and BGCs/ORFs)...')
get_dbs_metainfo.get(args.PATH, "fasta", root_database)


if args.bigscape:
  cprint.info("# Running global BiG-SCAPE2 analysis...")
  run_bigscape.run(path=args.PATH, threads=args.threads, cutoff=args.bigscape_cutoff, env=args.bigscape_env, pfam=args.pfam)

if args.bigscape:
  cprint.info("# Integrating BiG-SCAPE results with Krill...")
  integrate_bigscape_with_krill(path=args.PATH, cutoff=args.bigscape_cutoff)

if args.bigscape_analysis:
  cprint.info("# Running optional category analysis for BiG-SCAPE results...")
  run_bigscape_analysis.run(path=args.PATH, env=args.bigscape_env, threshold=args.bigscape_analysis_threshold, permutations=args.bigscape_analysis_permutations)

cprint.info('# Building charts...')
build_charts.build_charts(args.PATH)

cprint.ok('\nAll done. Any questions please contact: henrique.niero@lnbio.cnpem.br \nCheers!')
