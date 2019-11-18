
mutsColClasses <- function() {
  classes <- c("character","character","integer","character","character","character","integer","integer","character")
}

getBases <- function () {
  bases <- c("A","C","G","T","N")
}
getAscii <- function () {
  ascii <- c(65,67,71,84,78)
}
getBasecolors <- function () {
  basecolors <- brewer.pal(7,"Set1")
}

getMutsFromRead <- function(read,muts) {
  return(muts[muts$Expt == read$Expt[1] & muts$Read == read$Read[1],])
}

getReadFromMut <- function(mut,reads) {
  return(reads[reads$Expt == mut$Expt[1] & reads$Read == mut$Read[1],])
}

getMutsFromReads <- function(reads,muts) {
  return(merge(reads[,c("Expt","Read")],muts)) #,by=1:2)[,colnames(muts)])
}

getReadsFromMuts <- function(reads,muts) {
  return(unique(merge(reads,muts[,c("Expt","Read")]))) #,by=1:2)[,colnames(reads)]))
}

convertCoords <- function(coordlist) {
  coordlist <- strsplit(coordlist,",")
  coords <- ldply(lapply(1:length(coordlist),function(i,x) {
                                              df <- ldply(strsplit(x[[i]],"-"))
                                              df <- cbind(i,df)
                                              return(df)},coordlist))
  colnames(coords) <- c("i","start","end")
  coords$start <- as.integer(coords$start)
  coords$end <- as.integer(coords$end)
  return(coords)
}



invertCoords <- function(coordlist,refseq) {
  coords <- convertCoords(coordlist)
  incoords <- coords[c(),]
  refseqlength <- nchar(as.character(refseq))

  # quick reference to check if read
  nextread <- diff(coords$i)

  # Initialize first read
  if (coords$start[1] > 1) incoords[1,] <- c(coords$i[1],1,coords$start[1]-1)

  if (length(nextread) > 0) {
    for (i in 1:length(nextread)) {
      if (nextread[i]) {
        # end last read and start new read
        if (refseqlength > coords$end[i]) incoords[nrow(incoords)+1,] <- c(coords$i[i],coords$end[i]+1,refseqlength)

        if (coords$start[i+1] > 1) incoords[nrow(incoords)+1,] <- c(coords$i[i+1],1,coords$start[i+1]-1)
      } else {
        # take down coords within read
        incoords[nrow(incoords)+1,] <- c(coords$i[i],coords$end[i]+1,coords$start[i+1]-1)
      }
    }
  } else {
    if (refseqlength > coords$end[1]) incoords[nrow(incoords)+1,] <- c(coords$i[1],coords$end[1]+1,refseqlength)
  }

  colnames(incoords) <- c("i","start","end")
  return(incoords)

}

createMutationMatrix <- function(reads,muts,refseq,tstart,tend) {
  mutmat <- matrix(".",ncol=nchar(as.character(refseq)),nrow=nrow(reads))
  incoords <- invertCoords(reads$Coords,refseq)

  for (i in 1:nrow(incoords)) {
    mutmat[incoords[i,"i"],incoords[i,"start"]:incoords[i,"end"]] <- ""
  }


  subs <- muts[muts$Type=="sub",]
  dels <- muts[muts$Type=="del",]

  if (nrow(subs) > 0) {
    for (i in 1:nrow(subs)) {
      sub <- subs[i,]
      read <- getReadFromMut(sub,reads)
      row <- which(reads$Expt == read$Expt[1] & reads$Read == read$Read[1])
      mutmat[row,sub$Pos] <- sub$To
    }
  }

  if (nrow(dels) > 0) {
    for (i in 1:nrow(dels)) {
      del <- dels[i,]
      read <- getReadFromMut(del,reads)
      row <- which(reads$Expt == read$Expt[1] & reads$Read == read$Read[1])
      if (del$Size > 1) {
        mutmat[row,del$Pos] <- "<"
        if (del$Size > 2) mutmat[row,(del$Pos+1):(del$Pos+del$Size-2)] <- "-"
        mutmat[row,del$Pos+del$Size-1] <- ">"
      } else {
        mutmat[row,del$Pos] <- "<>"
      }
    }
  }

  if (tstart > 1) {
    mutmat[,1:(tstart-1)] <- ""
  }

  if (tend < ncol(mutmat)) {
    mutmat[,(tend+1):ncol(mutmat)] <- ""
  }

  return(mutmat)
}

createInsertionMatrix <- function(reads,muts,refseq,tstart,tend) {
  insmat <- matrix(".",ncol=nchar(as.character(refseq)),nrow=nrow(reads))
  incoords <- invertCoords(reads$Coords,refseq)

  for (i in 1:nrow(incoords)) {
    insmat[incoords[i,"i"],incoords[i,"start"]:incoords[i,"end"]] <- ""
  }

  inss <- muts[muts$Type=="ins",]

  if (nrow(inss) > 0) {
    for (i in 1:nrow(inss)) {
      ins <- inss[i,]
      read <- getReadFromMut(ins,reads)
      row <- which(reads$Expt == read$Expt[1] & reads$Read == read$Read[1])
      insmat[row,ins$Pos] <- ins$Ins
    }
  }


  if (tstart > 1) {
    insmat[,1:(tstart-1)] <- ""
  }

  if (tend < ncol(insmat)) {
    insmat[,(tend+1):ncol(insmat)] <- ""
  }

  return(insmat)
}


compareMutations <- function(mutpair) {
  apply(mutpair,2,function(x) {
    mutchars <- c("A","C","G","T")
    delchars <- c("<",">","<>")
    result <- rep(0,3)

    if (all(x %in% c("","."))) return(result)

    if (x[1] %in% mutchars && x[2] != "") result[1] <- result[1] + 1
    if (x[2] %in% mutchars && x[1] != "") result[2] <- result[2] + 1
    if (x[1] %in% mutchars && x[2] %in% mutchars && x[1] == x[2]) result[3] <- result[3] + 1

    if (x[1] %in% delchars && x[2] != "") result[1] <- result[1] + 1
    if (x[2] %in% delchars && x[1] != "") result[2] <- result[2] + 1
    if ( (grepl(delchars[1],x[1]) && grepl(delchars[1],x[2])) ||
           (grepl(delchars[2],x[1]) && grepl(delchars[2],x[2])) ) result[3] <- result[3] + 1

    if (grepl("[ACGTN]",x[3]) && x[4] != "") result[1] <- result[1] + 1
    if (grepl("[ACGTN]",x[4]) && x[3] != "") result[2] <- result[2] + 1
    if (grepl("[ACGTN]",x[3]) && grepl("[ACGTN]",x[4]) && adist(x[3],x[4]) < 2) result[3] <- result[3] + 1
    return(result)
  })
}


calculateProfile <- function(mutmat,refseq) {
  profile <- data.frame(Pos=1:nchar(as.character(refseq)))
  profile$Base <- unlist(strsplit(as.character(refseq),""))
  profile$Reads <- apply(mutmat,2,function(x) {sum(grepl("[.ACGTN]",x))})
  profile$Subs <- apply(mutmat,2,function(x) {sum(grepl("[ACGT]",x))})
  profile$Y <- ifelse(profile$Reads > 0, profile$Subs/profile$Reads, 0)
  return(profile)
}

calculateDeletionProfile <- function(mutmat,refseq) {
  profile <- data.frame(Pos=1:nchar(as.character(refseq)))
  profile$Base <- unlist(strsplit(as.character(refseq),""))
  profile$Reads <- colSums(mutmat != "")
  profile$Dels <- apply(mutmat,2,function(x) {sum(grepl("[<>-]",x))})
  profile$Y <- ifelse(profile$Reads > 0, profile$Dels/profile$Reads, 0)
  return(profile)
}

calculateWeightedDeletionProfile <- function(mutmat,refseq) {
  profile <- data.frame(Pos=1:nchar(as.character(refseq)))

}
