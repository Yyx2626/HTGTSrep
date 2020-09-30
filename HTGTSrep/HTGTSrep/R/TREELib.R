concatLable <- function(label){
    # use this func to concatenate sample:count
    b = strsplit(gsub("\\|", ",", label), ",")
    c = t(sapply(strsplit(unlist(b), ":"), `[`))
    f = data.frame(name=c[,1], count=as.numeric(c[,2]))
    concatDf = ddply(f, "name", numcolwise(sum))
    #
    q = strsplit(strsplit(label, ",")[[1]][1], "\\|")
    target = sapply(strsplit(unlist(q), ":"), `[`)[1,]

    m = concatDf[match(target, concatDf$name),]
    m = m[m$count>0,]
    outlabel = paste(m$name, m$count, sep=":", collapse=',')
    if(nrow(m)>0){
        return(c(nrow(m), outlabel))
    } else {
        return(c(0, NA))
    }
}

diffInTwoSeq <- function(childseq, parentseq){
    if(nchar(childseq) != nchar(parentseq)){
        stop("Error: Child & parent sequence should be as same length.")
    }
    j = c()
    for (I in 1:nchar(childseq)){
        p = substr(parentseq, I, I)
        c = substr(childseq, I, I)
        if(p != "N" & c != "N" & p != c){
            mismatch = paste(I, p, c, sep="")
            j = append(j, mismatch)
        }
    }
    if (length(j)>0){
        stringdif = paste(j, collapse="")
    } else {
        stringdif = '-'
    }
    return(stringdif)
}

treeConstruct <- function(graph, db){
    V(graph)$label <- c(1:length(V(graph)$SHORTCOUNT))
    t <-V(graph)$SHORTCOUNT
    tmp1 = strsplit(gsub("\\|", ",", t), ",")
    ss = sapply(strsplit(unlist(tmp1), ":"), `[`, 1)
    ss = unique(ss[!is.na(ss)])
    sample_indexes <- vector(mode="list", length=length(ss))
    names(sample_indexes) <- ss

    nl = data.frame("SAMPLENUM"=integer(), "ANNOTATION"=character(),
        stringsAsFactors=FALSE)
    for (I in 1:length(t)){
        if(is.na(t[I])) {
            nl[nrow(nl)+1,] = c(0, NA)
        } else {
            a = concatLable(t[I])
            b = strsplit(a[2], ",")
            c = sapply(strsplit(unlist(b), ":"), `[`, 1)
            for (d in c){
                if(is.null(sample_indexes[[d]])){
                    sample_indexes[[d]][1] = V(graph)$name[I]
                } else {
                    n = length(sample_indexes[[d]])
                    sample_indexes[[d]][n+1] = V(graph)$name[I]
                }
            }
            nl[nrow(nl)+1,] = a
        }
    }

    # Plot graph
    V(graph)$SAMPLENUM <- nl$SAMPLENUM
    V(graph)$ANNOTATION <- nl$ANNOTATION
    V(graph)$color <- "white"
    V(graph)$color[V(graph)$SAMPLENUM == 2] <- "green"
    V(graph)$color[V(graph)$SAMPLENUM == 3] <- "skyblue"
    V(graph)$color[V(graph)$SAMPLENUM >= 4] <- "orange"
    #V(graph)$color[V(graph)$SAMPLENUM > 4] <- "red"
    V(graph)$color[V(graph)$name == "Germline"] <- "black"
    V(graph)$color[grepl("Inferred", V(graph)$name)] <- "gray"


    V(graph)$shape <- "circle"
    #V(graph)$shape[V(graph)$name == "Germline"] <- "circle"
    #V(graph)$shape[grepl("Inferred", V(graph)$name)] <- "circle"


    if(length(t) > 100){
        vsize = 2
        pheight = pwidth = 100
        vcex = 2
        ecex = 2
        lcex = 5
    } else {
        vsize = 10
        pheight = pwidth = 7
        vcex = 1
        ecex = 1
        lcex = .75
    }

    V(graph)$size <- vsize
    #V(graph)$size[V(graph)$name == "Germline"] <- 2
    #V(graph)$size[grepl("Inferred", V(graph)$name)] <- 2

    #E(graph)$label <- ""
    dirpath = paste(dirname(datafile), "/lineageTree/", sep="")
    pdf(paste(dirpath, graph$clone, ".pdf", sep=""), height = pheight, width = pwidth)
    plot(graph, layout=layout_as_tree, edge.arrow.mode=0, edge.label.color="black",
         edge.label.cex = ecex, vertex.frame.color="black", vertex.label.color="black",
         vertex.label.cex=vcex, main=paste("Clone ", graph$clone))
    legend("topright", c("Germline", "Inferred", "SN=1", "SN=2", "SN=3", "SN>=4"),
       fill=c("black", "gray", "white", "green", "skyblue", "orange"), cex=lcex)
    graphics.off()

    # write graph file
    write_graph(graph, paste(dirpath, graph$clone, ".graphml", sep=""), "graphml")
    # write annotation file
    el = as_edgelist(graph)
    stringdifVec = c()
    parentVec = c()
    for (childname in V(graph)$name){
        if (childname == 'Germline'){
            stringdifVec = append(stringdifVec, '-')
            parentVec = append(parentVec, '-')
        } else {
            parentname = el[el[,2]==childname][1]
            childseq = V(graph)$sequence[V(graph)$name == childname]
            parentseq = V(graph)$sequence[V(graph)$name == parentname]
            stringdif = diffInTwoSeq(childseq, parentseq)
            stringdifVec = append(stringdifVec, stringdif)
            parentVec = append(parentVec, V(graph)$label[V(graph)$name == parentname])
        }
    }
    #V(graph)$name[grep("Inferred", V(graph)$name, invert=TRUE)] <- '-'
    outmat <- cbind(V(graph)$label, V(graph)$name, V(graph)$SAMPLENUM, V(graph)$ANNOTATION, parentVec, stringdifVec, V(graph)$sequence)
    colnames(outmat) <- c("INDEX", "NAME", "SAMPLENUM", "SAMPLEDETAIL", "PARENTINDEX", "MUTATION", "SEQUENCE")
    write.table(outmat, paste(dirpath, graph$clone, ".txt", sep=""), sep = "\t",row.names = F, quote = F)

    # func to construct a tree for each individual sample
    treeConstruct_sample(sample_indexes, graph, dirpath)
}

treeConstruct_sample <- function(sample_indexes, graph, dirpath){
    for (sample in names(sample_indexes)){
        l <- unique(sample_indexes[[sample]])
        if (!is.null(l) & length(l) > 1){
            sdb <- db[db$SEQUENCE_ID %in% l,]
            sdb$sequence_id <- sdb$SEQUENCE_ID
            sdb$sequence_alignment <- sdb$SEQUENCE_IMGT
            sdb$germline_alignment <- sdb$GERMLINE_IMGT_D_MASK
            sdb$junction_length <- sdb$JUNCTION_LENGTH
            sdb$v_call <- sdb$V_CALL
            sdb$j_call <- sdb$J_CALL
            sdb$clone_id <- sdb$CLONE
            sclone <- makeChangeoClone(sdb)
            sgraph <- buildPhylipLineage(sclone, dnapars_exec, rm_temp=TRUE)
            V(sgraph)$label[V(sgraph)$name %in% l] <- V(graph)$label[V(graph)$name %in% l]
            subtreeConstruct(sgraph, "sample", paste(graph$clone, sample, sep="."))

            outmat <- cbind(V(sgraph)$label, V(sgraph)$name, V(sgraph)$sequence)
            colnames(outmat) <- c("INDEX", "NAME", "SEQUENCE")
            write.table(outmat, paste(dirpath, graph$clone, ".", sample, ".txt", sep=""), sep = "\t",row.names = F, quote = F)
            # message("TREELib outmat = ", colnames(outmat))
        }
    }
}

subtreeConstruct <- function(graph, type, subname){
    V(graph)$color <- "white"
    V(graph)$color[V(graph)$SAMPLENUM == 2] <- "green"
    V(graph)$color[V(graph)$SAMPLENUM == 3] <- "skyblue"
    V(graph)$color[V(graph)$SAMPLENUM >= 4] <- "orange"
    V(graph)$color[V(graph)$name == "Germline"] <- "black"
    V(graph)$color[grepl("Inferred", V(graph)$name)] <- "gray"


    V(graph)$shape <- "circle"

    if(length(V(graph)) > 100){
        vsize = 2
        pheight = pwidth = 100
        vcex = 2
        ecex = 2
        lcex = 5
    } else {
        vsize = 10
        pheight = pwidth = 7
        vcex = 1
        ecex = 1
        lcex = .75
    }

    V(graph)$size <- vsize

    if(type == "subset"){
        dirpath = dirname(graphmlFile)
        clonenum = strsplit(basename(graphmlFile), "\\.")[[1]][1]
        fn = paste(dirpath, "/", clonenum, ".", indexnum, ".pdf", sep="")
        title = paste("Clone", graph$clone, "Sub ", indexnum)
        pdf(fn, height = pheight, width = pwidth)
        plot(graph, layout=layout_as_tree, edge.arrow.mode=0, edge.label.color="black",
             edge.label.cex = ecex, vertex.frame.color="black", vertex.label.color="black",
             vertex.label.cex=vcex, main=title)
        legend("topright", c("Germline", "Inferred", "SN=1", "SN=2", "SN=3", "SN>=4"),
           fill=c("black", "gray", "white", "green", "skyblue", "orange"), cex=lcex)
        graphics.off()
    } else {
        dirpath = paste(dirname(datafile), "/lineageTree/", sep="")
        fn = paste(dirpath, subname, ".pdf", sep="")
        title = paste("Clone", subname)
        pdf(fn, height = pheight, width = pwidth)
        plot(graph, layout=layout_as_tree, edge.arrow.mode=0, edge.label.color="black",
             edge.label.cex = ecex, vertex.frame.color="black", vertex.label.color="black",
             vertex.label.cex=vcex, main=title)
        legend("topright", c("Germline", "Inferred", "Reads"),
           fill=c("black", "gray", "white"), cex=lcex)
        graphics.off()
    }

}
