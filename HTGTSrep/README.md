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



## Run pipeline

### Prepare metadata file

Metadata file is a tabular text file, which should at least containing these columns: Library, Sequencing, Primer, MID, Adapter

The sample name (in file name) is determined by Library + "\_" + Sequencing;
for example, if Library=ABC, Sequencing=DEF, then the inferred sample name in filename is ABC\_DEF.


### Run HTGTSrep pipeline

To run the pipeline, you should always activate HTGTS-Rep environment except you have all dependent softwares installed.

    source activate HTGTSrep

Run the pipeline using subcommand "run" for V,D,J usage analysis.

    python HTGTSrep.py run

Run the pipeline using subcommand "mut" for gene mutation profiling.

    python HTGTSrep.py mut

Run the pipeline using subcommand "clonal" for clonal clustering and phylogenetic tree construction.

    python HTGTSrep.py clonal
