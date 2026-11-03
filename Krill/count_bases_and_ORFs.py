import subprocess, os, tempfile
import pandas as pd
import os


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

def run(path, ext, threads):
    os.chdir(path)
    
    orffinder = os.path.join(
                            os.path.dirname(__file__),
                            'ORF_Finder','ORFfinder'
                            )

    bash="""cp {orffinder} {path} && \
    chmod u+x ORFfinder && \
    echo "FILE,NT,ORFS" && \
    mkdir -p ORFs && \
    parallel -j{threads} 'count=$(grep -v ">" {{}} | tr -d '\\n' | wc -c) && file=$(echo "{{.}}") && ./ORFfinder -in {{}} -g 11 -s 1 -n true -out ORFs/"{{}}" && ORFS=$(grep -c ">" ORFs/"{{}}") && echo -n "$file," && echo -n "$count," && echo $ORFS' ::: *.{ext} && \
    rm ORFfinder""".format(orffinder=orffinder,path=path,ext=ext,threads=threads)
    
    tmp = tempfile.NamedTemporaryFile(suffix='.sh',mode='w+',delete=False)
    tmp.write(bash)
    tmp.seek(0)
    tmp.close()
    
    os.chmod(tmp.name,0o777)

    output, err = subprocess.Popen([tmp.name],
                           shell=True,
                           executable='/bin/bash',
                           universal_newlines=True,
                           text=True,
                           stdout=subprocess.PIPE).communicate()

    os.remove(tmp.name)

    unformatted_tsv = [x.split(',') for x in output.splitlines()]
    df = pd.DataFrame(unformatted_tsv[1:], columns=unformatted_tsv[0])
    
    
    files_list = ''
    files_sizes = {'FILE': [], 'FILE SIZE (MB)': []}
    for file in os.listdir(path):
        file_ext = file.split('.')[-1]
        if ext in file_ext:
            files_list += file + ' '
            files_sizes['FILE'].append(file.replace(f'.{ext}', ''))
            files_sizes['FILE SIZE (MB)'].append(round(os.stat(file).st_size/10**6, 1))

    files_sizes = pd.DataFrame(files_sizes)

    output = subprocess.run(
        f'seqkit stats -a {files_list}',
        shell=True,
        executable='/bin/bash',
        capture_output = True,
        text = True
    )   
    output.stdout = output.stdout.replace(',', '')
    seqkit_lines = output.stdout.split('\n')
    data = []
    import re
    for pos, line in enumerate(seqkit_lines):
        line = re.sub(' +', '\t', line)
        fields = line.split('\t')
        if pos == 0:
            columns = fields
        elif fields != ['']:
            data.append(fields)

    seqkit_results = pd.DataFrame(data, columns=columns)
    seqkit_results['file'] = seqkit_results['file'].str.replace(f'.{ext}', '')
    seqkit_results = seqkit_results.rename(columns={'file': 'FILE', 'num_seqs': 'CONTIGS', 'min_len': 'MIN LEN (KB)', 'avg_len': 'AVG LEN (KB)', 'max_len': 'MAX LEN (KB)'})
    seqkit_results = seqkit_results [['FILE', 'CONTIGS', 'MIN LEN (KB)', 'MAX LEN (KB)', 'AVG LEN (KB)']]
    

    seqkit_results['CONTIGS'] = seqkit_results['CONTIGS'].astype(int)
    seqkit_results['MIN LEN (KB)'] = (seqkit_results['MIN LEN (KB)'].astype(int) / 1000).round(3)
    seqkit_results['MAX LEN (KB)'] = (seqkit_results['MAX LEN (KB)'].astype(int) / 1000).round(3)
    seqkit_results['AVG LEN (KB)'] = (seqkit_results['AVG LEN (KB)'].astype(float) / 1000).round(3)
    
    df = df.merge(files_sizes, on='FILE')
    df = df.merge(seqkit_results, on='FILE')
    
    
    df = df.rename(columns={'NT': 'NT (KB)'})
    df['NT (KB)'] = (df['NT (KB)'].astype(int) / 1000)

    col2move = df.pop('FILE SIZE (MB)') 
  
    df.insert(1, 'FILE SIZE (MB)', col2move) 
    df.to_csv(os.path.join(path,'ReportOutput/DNABases_and_ORFs_count.tsv'),sep='\t',index=False)
    df2xlsx(os.path.join(path,'ReportOutput/DNABases_and_ORFs_count.xlsx'), 'DNABases and ORFs count', df)
    

if __name__ == '__main__':
    run('/media/bioinfo/6tb_hdd/03_ELLEN/krill_runs/01_REPORT_DATABASES/KRILL_RESULTS/DeceptionIslandSequences/', 'fasta', 8)