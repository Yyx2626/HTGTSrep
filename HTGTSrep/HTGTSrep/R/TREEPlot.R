suppressPackageStartupMessages(library(alakazam, quietly=TRUE))
suppressPackageStartupMessages(library(igraph, quietly=TRUE, verbose=FALSE))
suppressPackageStartupMessages(library(dplyr, quietly=TRUE, verbose=FALSE))
suppressPackageStartupMessages(library(plyr, quietly=TRUE, verbose=FALSE))


ARGS <- c(
  "datafile","character","file path of data file",
  "dnaparspath", "character","file path of dnapars program",
  "readnum","numeric","Minimum # of reads to construct tree"
)

source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}
source_local("Rsub.R")
source_local("SHMHelper.R")
source_local("TREELib.R")

parseArgs("TREEPlot.R", ARGS)

# db<-readChangeoDb('/Users/zhoudu/PycharmProjects/HTGTSrep/clonaltest/tmp_db-treeclonal_germ-collapse.tab')
# sub_db <- subset(db, CLONE == 1355)
# clone <- makeChangeoClone(sub_db, text_fields=c("SHORTCOUNT","SAMPLE","SEQUENCE_ID"), num_fields = "DUPCOUNT")
# dnapars_exec<-'/Users/zhoudu/PycharmProjects/HTGTSrep/igblast_bin/phylip-3.695/exe/dnapars'
# graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp=TRUE)
# treeConstruct(graph)

db<-readChangeoDb(datafile)

bigclone <- c()

for (cloneID in unlist(unique(db[,'CLONE']))){
    sub_db <- subset(db, CLONE == cloneID)
    ### MAYBE SHOULD JUST CHANGE THIS FUNCTION IN ITS PKG???
    sub_db$sequence_id <- sub_db$SEQUENCE_ID
    sub_db$sequence_alignment <- sub_db$SEQUENCE_IMGT
    sub_db$germline_alignment <- sub_db$GERMLINE_IMGT_D_MASK
    sub_db$junction_length <- sub_db$JUNCTION_LENGTH
    sub_db$v_call <- sub_db$V_CALL
    sub_db$j_call <- sub_db$J_CALL
    sub_db$clone_id <- sub_db$CLONE

    cloneRn = dim(sub_db)[1]
    if (cloneRn > 1000){
        bigclone = append(bigclone, cloneID)
    }
    if (cloneRn >= readnum & cloneRn < 1000){
        ### deleted print(cloneID)

        clone <- makeChangeoClone(sub_db, text_fields=c("SEQUENCE_ID","SAMPLE","SHORTCOUNT"), num_fields = "DUPCOUNT")
        dnapars_exec<-dnaparspath
        graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp=TRUE)
        treeConstruct(graph, db)
    }
}
#message("TREEPlot head bigclone: ",head(bigclone))
#message("TREEPlot colnames bigclone: ",colnames(bigclone))
write(bigclone, paste(dirname(datafile), "/bigclone.txt", sep=''))


# print("Loading data")
# db<-readChangeoDb(datafile)
# clones <- db %>%
#     group_by(CLONE) %>%
#     do(CHANGEO=makeChangeoClone(., #text_fields=c("SHORTCOUNT"),
#                                 num_fields="DUPCOUNT"))
# dnapars_exec<-dnaparspath
# print("prepare for pic")
# graphs <- lapply(clones$CHANGEO, buildPhylipLineage,
#                  dnapars_exec=dnapars_exec, rm_temp=TRUE)
# graphs[sapply(graphs, is.null)] <- NULL
# graphs <- graphs[sapply(graphs, vcount) >= readnum]
# print("Beginning construct tree...")
# #lapply(graphs, treeConstruct)
# for (I in 1:length(graphs)){
#     treeConstruct(graphs[[I]], db)
# }
