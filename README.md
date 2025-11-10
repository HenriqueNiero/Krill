# Krill

## :mag_right: SUMMARY
1. :scroll: ABOUT
2. :electric_plug: PRE-REQUISITES
3. :dvd: INSTALLATION
4. :woman_teacher: PREPARING YOUR FILES
5. :woman_technologist: USING

## :scroll: ABOUT

Krill is a tool designed to aid in the analysis of genomic data with a focus on biosynthetic gene clusters (BGCs).

<p align="justify">The package was developed with financial support from the Foundation for Research and Innovation Support of the State of Santa Catarina (FAPESC, process 2020TR1448) and the São Paulo Research Foundation (FAPESP, process 2019/27306-9), as well as the Brazilian National Institute of Science and Technology - INCT-Mar COI (CNPq, Process 400551/2014–4). We also wish to thank the Coordination for the Improvement of Higher Education (CAPES) for the scholarships provided to H.N. (88887.146746/2017–00). A.O.S.L. was further supported by CNPq (312363/2018–4). D.B.B.T. and P.B. also acknowledge funding from the Serrapilheira Institute (Serra-1709–19681 and associated fellowship ref. 3659).</p>

Package Features:
1. Data pre-processing for contigs extraction; 
2. Extraction of biosynthetic gene clusters (BGCs); 
3. Annotation of resistance genes; 
4. Similarity analysis and annotation of biosynthetic families; 
5. Curation of prospected clusters. 

## :electric_plug: PRE-REQUISITES
1. Ubuntu 20.04
2. Python 3.X
3. R
4. Conda
5. AntiSMASH
6. ARTS

Ubuntu/Linux packages:
1. hmmsearch
2. gawk
3. parallel
4. seqkit

Python packages:
1. Pandas
2. Matplotlib
3. CPRINT
4. Numpy
5. Biopython
6. TQDM
7. Xlsxwriter
8. pyScss

R packages:
1. ggplot2
2. dplyr
3. stringr
4. tidyr
5. gridExtra
6. igraph
7. plotly
8. vegan
9. reshape2
10. MASS

## :dvd: INSTALLATION

<details><summary>INSTALLING AntiSMASH via Docker</summary>
<p>
    
1. Pull [AntiSMASH Docker Image](https://hub.docker.com/r/antismash/standalone)
    
```
docker pull antismash/standalone
```
    
2. Prepare AntiSMASH run scripts
    
```
mkdir ~/bin    # not required if you already have that
curl -q https://dl.secondarymetabolites.org/releases/6.1.1/docker-run_antismash-full > ~/bin/run_antismash
chmod a+x ~/bin/run_antismash
export PATH="$HOME/bin:$PATH"
```
    
3. Test installation
    
```
run_antismash . . --version
```
    
or
    
```
run_antismash . . --version
```
    
</p>
</details>

<details><summary>INSTALLING AntiSMASH via Conda</summary>
<p>
    
1. Install AntiSMASH using conda
    
```
conda install -c bioconda antismash
```

</p>
</details>

<details><summary>INSTALLING ARTS via Conda</summary>

<p>
    
1. Download environment spec list file from this repository (arts_specs.txt)
    
2. Create a conda environment for ARTS using the spec list file
    
```
conda create -n "ARTS" --file /path/to/spec-file.txt
```
 
3. Download ARTS project into the conda environment using git
    
```
cd /path/to/ARTS/environment/
git clone https://bitbucket.org/ziemertlab/arts.git
```

</p>
    
</details>

<details><summary>GETTING ARTS REFERENCES</summary>

<p>
    
1. Download [Krill's specific ARTS references](https://bitbucket.org/krill-arts-reference/krillartsreference/downloads/) 
    
2. Unzip the file
3. Replace the "reference" folder in ARTS environment with the new one

</p>
    
</details>



## :woman_teacher: PREPARING YOUR FILES

Krill works in a single folder with fasta files or in a folder with different Projects MAGs (multiple folders). For this second option, it needs to have a [specific folder organization](example/) to start the analysis:

> :warning: DO NOT USE SPACES IN FOLDERS AND FILES NAMES, IT CAN CAUSE ERRORS.

<details><summary>FOLDERS AND FILES STRUCTURE</summary>
<p>
    
#### Flowchart Scheme
```mermaid
flowchart TB
    subgraph A[example/ - Main folder]
        subgraph B[PRJNA602601/ ]
        E[MAG_1.fasta]
        F[MAG_2.fasta]
        G[MAG_3.fasta]
        N[...]
        end
        subgraph C[Project_B/ ]
        H[MAG_1.fasta]
        I[MAG_2.fasta]
        J[MAG_3.fasta]
        O[...]
        end
        subgraph D[Project_C/ ]
        K[MAG_1.fasta]
        L[MAG_2.fasta]
        M[MAG_3.fasta]
        P[...]
        end
    end
```

#### Printscreen Scheme
<p align="center">
    <img src="https://user-images.githubusercontent.com/50638088/184180804-c794655e-3e4c-4509-b38a-3f63eac7c0d5.png"/>
</p>
</p>
</details>
    
## :woman_technologist: USING
```
Krill [OPTIONS] PATH
```

```
usage: Krill [-h] [-noprep] [-t THREADS] [--citation] PATH

positional arguments:
  PATH                  Working path with fasta files

optional arguments:
  -h, --help            show this help message and exit
  -noprep, --do_not_prepare_fasta_files
                        Rename fasta files, its headers and store changes in a CSV file for control [DEFAULT: TRUE]
  -t THREADS, --threads THREADS
                        Threads to use in analysis [DEFAULT: 16]
  --citation            Shows how to cite us
```

## :woman_technologist: :man_technologist: Contributors 
<table>
    <tr>
        <td align="left" valign="top" width="14.28%"><a href="https://github.com/saulobritto"><img src="https://avatars.githubusercontent.com/u/50638088?v=4" width="100px;" alt="Saulo Britto"/><br />Saulo Britto</a></td>
        <td align="left" valign="top" width="14.28%"><a href="https://github.com/ellenjkr"><img src="https://avatars.githubusercontent.com/u/49100458?v=4" width="100px;" alt="Ellen Junker"/><br />Ellen Junker</a></td>
        <td align="left" valign="top" width="14.28%"><a href="https://github.com/machenka-code"><img src="https://avatars.githubusercontent.com/u/195959754?v=4" width="100px;" alt="Maria Nascimento"/><br />Maria Nascimento</a></td>
    </tr>
</table>
