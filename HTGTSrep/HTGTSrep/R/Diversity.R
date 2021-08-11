suppressPackageStartupMessages(library(alakazam, quietly=TRUE))

ARGS <- c(
  "datafile","character","file path of data file - rquired columns: CLONE",
  "outprefix","character","prefix of output file"
)

OPTS <- c()

source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local("Rsub.R")
source_local("SHMHelper.R")

parseArgs("Diversity.R", ARGS, OPTS)

############### Start to work ###################

db <- read.csv(datafile, sep="\t")
### 092422020 Lawrence added db$clone_id <- db$CLONE
db$clone_id <- db$CLONE
clones <- estimateAbundance(db, group="SAMPLE", ci=0.95, nboot=200)
pdf(paste(outprefix, ".abundance.pdf", sep=""), width=7, height=5)
plotAbundanceCurve(clones, legend_title="Sample")
graphics.off()

# 'rarefyDiversity' is deprecated. Use 'alphaDiversity' instead. (same arguments; only changed function)
sample_div <- alphaDiversity(db, "SAMPLE", min_q=0, max_q=32, step_q=0.05,
                              ci=0.95, nboot=200)
sample_main <- paste0("Sample diversity (n=", sample_div@n, ")")
pdf(paste(outprefix, ".diversity.pdf", sep=""), width=7, height=5)
# 09302020 JH: Error: Theme element `log_q` is not defined in the element hierarchy. Taken out log_q=TRUE from plotDiversityCurve(arguments)
# 04262021 JH: Error: Theme element `log_d` is not defined in the element hierarchy. taken out lod_d=TRUE too
plotDiversityCurve(sample_div, main_title=sample_main,
                    legend_title="Sample")
graphics.off()
