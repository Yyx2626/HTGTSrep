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

clones <- estimateAbundance(db, group="SAMPLE", ci=0.95, nboot=200)
pdf(paste(outprefix, ".abundance.pdf", sep=""), width=7, height=5)
plotAbundance(clones, legend_title="Sample")
graphics.off()

sample_div <- rarefyDiversity(db, "SAMPLE", min_q=0, max_q=32, step_q=0.05,
                              ci=0.95, nboot=200)
sample_main <- paste0("Sample diversity (n=", sample_div@n, ")")
pdf(paste(outprefix, ".diversity.pdf", sep=""), width=7, height=5)
plotDiversityCurve(sample_div, main_title=sample_main,
                    legend_title="Sample", log_q=TRUE, log_d=TRUE)
graphics.off()
