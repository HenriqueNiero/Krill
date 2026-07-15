import pandas as pd
from pathlib import Path
import sys, os, subprocess
from cprint import *


def convert2fasta(path):
	for item in Path(path).rglob('*'):
		if os.path.isdir(item) is False:
			base, ext = os.path.splitext(item)
			if ext in ['.fa', '.fna', '.fas', '.ffn', '.faa', '.frn']:
				os.rename(item, base + ".fasta")


def run(path,ext):
	cprint.ok('# Preparing input files and its metadata...')
	os.chdir(path)
	if not os.path.isfile(os.path.join(path,'fastaFilesRenamed.tsv')):
		#get folder and files
		fastas = Path(path).glob('*.{}'.format(ext))

		#get info and generate metadata database
		base = []
		df_contigs_header = pd.DataFrame({'original_filename': [], 'new_filename': [], 'original_contigs': [], 'new_contigs': []})
		count = 1
		for f in fastas:			
			# print(f)
			basename = str(os.path.basename(f))
			new_name = f'{count:09d}.{ext}'
			os.rename(basename, new_name)
			base.append([basename,  new_name])
			# Define your pattern and file path
			# Run the awk command: awk '/pattern/ {print}' 
			# Split the output by lines and clean up any trailing empty lines	
			fileresult = subprocess.run(["awk", f"/^>/", new_name], capture_output=True, text=True)
			original_contigs = [line[1:] for line in fileresult.stdout.splitlines() if line]
			# print(original_contigs)
			cmd_rename_fasta = '''gawk -i inplace '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-%i); next} 1' %s''' % (len(ext)+1,new_name)
			cmd_rename_fasta_headers = f"""awk '{{if (/^>/) print ">"substr($0,2)"_"(++i); else print $0;}}' {new_name} > tmp && mv tmp {new_name}"""
			subprocess.run(cmd_rename_fasta, shell=True, executable='/bin/bash')
			subprocess.run(cmd_rename_fasta_headers, shell=True, executable='/bin/bash')
			# Split the output by lines and clean up any trailing empty lines	
			fileresult = subprocess.run(["awk", f"/^>/", new_name], capture_output=True, text=True)
			new_contigs = [line[1:] for line in fileresult.stdout.splitlines() if line]
			# print(new_contigs)
			df = pd.DataFrame({'original_filename': basename, 'new_filename': new_name, 'original_contigs': original_contigs, 'new_contigs': new_contigs})
			df_contigs_header = pd.concat([df_contigs_header, df])
			count += 1



		df = pd.DataFrame(base, columns=['OriginalName', 'NewName'])
		df.to_csv(os.path.join(path,'fastaFilesRenamed.tsv'),sep='\t',index=False)
		cprint.info('# Changed metadata saved in {}'.format(os.path.join(path,'fastaFilesRenamed.tsv')))

		df_contigs_header.to_csv(os.path.join(path,'contigsFilesRenamed.tsv'),sep='\t',index=False)
		cprint.info('# Changed contig metadata saved in {}'.format(os.path.join(path,'contigsFilesRenamed.tsv')))
	else:
		cprint.warn('Files already prepared.')


if __name__ == '__main__':
	path = sys.argv[1]
	# run(path,"fasta")
	convert2fasta(path)