import subprocess, os, datetime

def run(fasta_file, db):
    project_dir = os.path.dirname(os.path.abspath(__file__))
    name=str(os.path.basename(fasta_file)).rsplit('.',1)[0]
    time=datetime.datetime.now().strftime('%Y-%m-%d %X')
    
    # AntiSMASH
    if not os.path.isdir(os.path.join(db,'AntiSMASH',name)):
        subprocess.run(
            'bash {project_dir}/krill_run_antismash.sh {file} {thr} AntiSMASH AntiSMASH/{file_name}'.format(project_dir=project_dir, file=fasta_file, thr=1, file_name=name),
            shell=True,
            executable='/bin/bash'
        )
    
    #ARTS
    if not os.path.isdir(os.path.join(db,'ARTS',name)):
        # hmm_ref=os.path.join(os.path.dirname(__file__),'db','novos_hmm')
        # hmm_ref=os.path.join(os.path.dirname(__file__),'db','novos_hmm/arylformamidase_new.hmm')
        # hmm_ref=os.path.join(os.path.dirname(__file__),'db','tryptophanDioxygen_blast.hmm')
        # hmm_ref=os.path.join(os.path.dirname(__file__),'db','formate_releasing_enzymes_HMMs')
        hmm_ref=os.path.join(os.path.dirname(__file__),'db','bacterial_metagenomic_resist_models.hmm')
        res_dir=os.path.join(db,'ARTS',name)
        os.makedirs(res_dir,exist_ok=True)
  
        thr=1
        gbk=os.path.join(db,'AntiSMASH',name,name+'.gbk')


        # Block ARTS for showing work process in terminal
        subprocess.run(
            'bash {project_dir}/krill_run_arts.sh {input} {hmm_ref} {res_dir} {thr}'.format(project_dir=project_dir, input=gbk, hmm_ref=hmm_ref, res_dir=res_dir, thr=thr),
            shell=True,
            executable='/bin/bash',
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT
        )


        # Show ARTS working in terminal
        #result = subprocess.run(
        #    'bash {project_dir}/krill_run_arts.sh {input} {hmm_ref} {res_dir} {thr}'.format(project_dir=project_dir, input=gbk, hmm_ref=hmm_ref, res_dir=res_dir, thr=thr),
        #    shell=True,
        #    executable='/bin/bash',
        #    capture_output=True,
        #    text=True
        #)
        #
        #print("Return code:", result.returncode)
        #print("STDOUT:")
        #print(result.stdout)
        #print("STDERR:")
        #print(result.stderr)

    return int(os.stat(fasta_file).st_size/(10**6))
