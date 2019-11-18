# HTGTS-Rep
HTGTS-Rep is a pipeline for comprehensive analysis of HTGTS-Rep-seq.

Author: Zhou Du, Nia Kyritsis, Adam Yongxin Ye @ Boston Children's Hospital / Harvard Medical School

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
Download proper version of [Igblast](ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.5.0/) and replace programs in igblast_bin folder if your OS is not Linux-based.

### Step 4
If necessary, install some more python and R packages needed for HTGTS-rep:
* Python
  * pandas
  * argparse
  * collections
  * multiprocessing
  * scipy
  * itertools
  * Biopython
  * random
  * sys
  * os
  * logging
  * re
  * glob
  * csv
  * time
* R
  * alakazam
  * tools
  * Rcolorbrewer
  * Biostrings (not available for R 3.6.1)
  * igraph
  * dplyr
  * plyr


## Run pipeline
To run the pipeline, you should always activate HTGTS-Rep environment except you have all dependent softwares installed.

    source activate HTGTSrep

Run the pipeline using subcommand "run" for V,D,J usage analysis.

    python HTGTSrep.py run

Run the pipeline using subcommand "mut" for gene mutation profiling.

    python HTGTSrep.py mut

Run the pipeline using subcommand "clonal" for clonal clustering and phylogenetic tree construction.

    python HTGTSrep.py clonal
