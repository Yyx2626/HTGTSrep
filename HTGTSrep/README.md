# HTGTS-Rep
HTGTS-Rep is a pipeline for comprehensive analysis of HTGTS-Rep-seq.

Author: Zhou Du, Huan Chen, Nia Kyritsis, Jianqiao Hu, Adam Yongxin Ye @ Boston Children's Hospital / Harvard Medical School

## Installation

### Step 1
Download HTGTS-Rep pipeline.

### Step 2
If necessary, install [Miniconda](http://conda.pydata.org/miniconda.html).
Create development conda environment with

    conda env create environment.yml

Activate the environment with

    source activate HTGTSrep

### Step 3
Download proper/latest version of [Igblast](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/), unpack and copy to folder 'external_software/' if your OS is not Linux-based.

### Note
If necessary, manually install some more python and R packages needed for HTGTS-rep:
* Python
  * **Biopython**
  * **itertools**
  * pandas
  * argparse
  * collections
  * multiprocessing
  * scipy
  * random
  * sys
  * os
  * logging
  * re
  * glob
  * csv
  * time
* R
  * **alakazam**
  * **Rcolorbrewer**
  * **Biostrings**
  * tools
  * igraph
  * dplyr
  * plyr


### Demo test

To test whether HTGTSrep pipeline is successfully installed, we provided a demo test dataset in folder 'test/'. To run the demo test, execute:

    python HTGTSrep.py run -m test/metadata.txt -r1 test/r1.fq.gz -r2 test/r2.fq.gz -o test/out

It should generate a folder 'test/out/', in which the output files should be similar to those files in folder 'test/demo_out/'.


## Run pipeline

### Prepare metadata file

Metadata file is a tabular text file, whose format is compatible to our previous [HTGTS pipeline](https://robinmeyers.github.io/transloc_pipeline/). Some columns (Breakseq, Breaksite, Chr, Start, End) in the demo test metadata file are not used in HTGTSrep pipeline, and can be ignored. The metadata file for HTGTSrep pipeline should at least containing these columns: Library, Sequencing, Primer, MID, Adapter.

The sample name (in file name) is determined by Library + "\_" + Sequencing;
for example, if Library=ABC, Sequencing=DEF, then the inferred sample name in filename is ABC\_DEF.


### Run HTGTSrep pipeline

To run the pipeline, you should always activate HTGTS-Rep environment except you have all dependent softwares installed.

    source activate HTGTSrep

Run the pipeline using subcommand "run" for V,D,J usage analysis.

    python HTGTSrep.py run

Run the pipeline using subcommand "mut" for gene mutation profiling.

    python HTGTSrep.py mut

Run the pipeline using subcommand "clonal" for clonotype clustering and phylogenetic tree construction.

    python HTGTSrep.py clonal


### Command-line document

Use option '-h' to see the document (the explanation of each options) in command-line. For example,

    python HTGTSrep.py -h

```
usage: HTGTSrep.py [-h] [--version] {preprocess,run,mut,clonal} ...

HTGTS repertoire sequencing analysis pipeline

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

subcommands:
  {preprocess,run,mut,clonal}
    preprocess          Preprocess sequences files.
    run                 Run preprocess and VDJ assignment from fastq/fasta
                        files.
    mut                 Generate DNA and protein mutation profile along V
                        genes. Required: log file in VDJ assignment analysis
    clonal              Clonal clustering and phylogenetic tree construction.
                        Required: log file in VDJ assignment analysis
```

    python HTGTSrep.py run -h
    
```
usage: HTGTSrep.py run [-h] [-r1 R1] [-r2 R2] [-i INPUT] -m METAFILE -o OUTDIR
                       [--keepunmatch] [--demulti_length DEMULTI_LENGTH]
                       [--demulti_mismatch DEMULTI_MISMATCH]
                       [--overlap OVERLAP] [--diffpercent DIFFPERCENT]
                       [--qscore QSCORE] [--qscore_cov QSCORE_COV]
                       [--fastq_quality_trimmer FASTQ_QUALITY_TRIMMER]
                       [--organism {mouse,human}] [--mousestrain {B6,all}]
                       [--VDJdatabase {IGH,IGK,IGL,TRA,TRB,TRG,TRD}]
                       [--Vdb VDB] [--Ddb DDB] [--Jdb JDB]
                       [--domain_system {kabat,imgt}]
                       [--auxiliary_data AUXILIARY_DATA]
                       [--Vannotation VANNOTATION] [--Vgapseq VGAPSEQ]
                       [--V_score V_SCORE] [--V_identity V_IDENTITY]
                       [--V_coverage V_COVERAGE] [--V_length V_LENGTH]
                       [--J_gene J_GENE] [--J_length J_LENGTH]
                       [--checkProductive] [--skipDedup] [--skipIgBlast]
                       [--skipUnjoined] [--skipDemultiplex] [--D_upstream]
                       [--D_upstream_length D_UPSTREAM_LENGTH]
                       [--D_upstream_stitch_length D_UPSTREAM_STITCH_LENGTH]
                       [--dedup_by_col] [--nproc NPROC]

optional arguments:
  -h, --help            show this help message and exit
  -r1 R1                MiSeq read 1 file in fastq format, mutually exclusive
                        of fasta file
  -r2 R2                MiSeq read 2 file in fastq format
  -i INPUT              Read files in fasta format, seperated multiple files
                        by comma in the same order as in metadata file,
                        mutually exclusive -r1 and -r2 options
  -m METAFILE           Metadata file
  -o OUTDIR             Output directory
  --keepunmatch         Keep unmatched reads from demultiplexing
  --demulti_length DEMULTI_LENGTH
                        Use N (12) bp of barcode+primer to demultiplex
  --demulti_mismatch DEMULTI_MISMATCH
                        N-maximum mismatch (0) to demultiplex
  --overlap OVERLAP     N-minimum overlap to join paired reads (10)
  --diffpercent DIFFPERCENT
                        N-percent maximum difference to join paired reads (8)
  --qscore QSCORE       Mark nucleotide with lower quality score as N (10) in
                        joined sequence
  --qscore_cov QSCORE_COV
                        Minimum percent of bases in a joined sequence have
                        [--qscore] quality score (98)
  --fastq_quality_trimmer FASTQ_QUALITY_TRIMMER
                        Quality threshold - nucleotides with lower quality
                        score will be trimmed from the end of the sequence.
                        (0)
  --organism {mouse,human}
                        Select sample organism, mouse or human (mouse)
  --mousestrain {B6,all}
                        Use mouse C57BL database only or all gene database in
                        the IgBlast searching. NO setting when organism is
                        human (all)
  --VDJdatabase {IGH,IGK,IGL,TRA,TRB,TRG,TRD}
                        Specify database for searching. (IGH)
  --Vdb VDB             Specify IgBlast germline V database.
  --Ddb DDB             Specify IgBlast germline D database.
  --Jdb JDB             Specify IgBlast germline J database.
  --domain_system {kabat,imgt}
                        Domain system kabat or imgt to be used for segment
                        annotation (imgt).
  --auxiliary_data AUXILIARY_DATA
                        Specify file containing the coding frame start
                        positions for sequences in germline J database.
  --Vannotation VANNOTATION
                        Specify V gene annotation file which has three
                        columns: V gene name, coordinates and function.
  --Vgapseq VGAPSEQ     Specify V IMGT-gapped sequences file.
  --V_score V_SCORE     Minimum IgBlast score of V gene (150)
  --V_identity V_IDENTITY
                        Minimum V identity ratio (0.9)
  --V_coverage V_COVERAGE
                        Minimum V coverage ratio (0.1)
  --V_length V_LENGTH   Minimum aligned length of V gene (100)
  --J_gene J_GENE       Specify J gene for alignment. eg: IGHJ4
  --J_length J_LENGTH   Minimum alignment length of J gene (34).
  --checkProductive     Include only records with productive information
  --skipDedup           Specify to skip generate duplicated-removed output
  --skipIgBlast         Specify to skip IgBlast and all upstream steps
  --skipUnjoined        Specify to skip unjoined reads in parsing igblast db
  --skipDemultiplex     Specifiy to skip Demultiplex and all upstream steps
  --D_upstream          Specify to consider D upstream regions in the
                        analysis. The D upstream sequences should be added in
                        the V database and named as 'D gene' plus '_DS'.
  --D_upstream_length D_UPSTREAM_LENGTH
                        Minimum aligned length of D upstream region (100)
  --D_upstream_stitch_length D_UPSTREAM_STITCH_LENGTH
                        Minimum aligned length of D gene D stitched upstream
                        region (70)
  --dedup_by_col        Specify to dedup in the same manner as Version 1
                        pipeline
  --nproc NPROC         The number of simultaneous computational processes to
                        execute (CPU cores to utilized)
```

