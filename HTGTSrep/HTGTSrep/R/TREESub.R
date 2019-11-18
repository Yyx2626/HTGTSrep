suppressPackageStartupMessages(library(igraph, quietly=TRUE, verbose=FALSE))

ARGS <- c(
  "graphmlFile","character","file path of graphml file",
  "indexnum","numeric","index of vertex as the root for construct a new subtree"
)

source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}


source_local("Rsub.R")
source_local("SHMHelper.R")
source_local("TREELib.R")

parseArgs("TREESub.R", ARGS)

graph <- read_graph(graphmlFile, "graphml")
a = induced_subgraph(graph, subcomponent(graph, indexnum, "out"))
subtreeConstruct(a, "subset", indexnum)
