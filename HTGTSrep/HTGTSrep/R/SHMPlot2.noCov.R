ARGS <- c(
  "datafile","character","file path of data file - columns: Pos, Base, Y, Err",
  "output","character","file path of plot pdf"
)

OPTS <- c(
  "tstart","numeric",0,"Start of reference to view",
  "tend","numeric",0,"End of reference to view",
  "plotrows","numeric",4,"Rows on plot",
  "ymax","numeric",0,"Maximum y-axis height",
  "figureheight","numeric",8,"height in inches",
  "showsequence","logical",TRUE,"display sequence on plots",
  "regex1","character","AGCT","",
  "regex2","character","[AGT](?=G[CT][AT])", "",   # this is DGYW/WRCH motif
  "cdr1_start","numeric",0,"Start of cdr1 region",
  "cdr1_end","numeric",0,"End of cdr1 region",
  "cdr2_start","numeric",0,"Start of cdr2 region",
  "cdr2_end","numeric",0,"End of cdr2 region",
  "cdr3_start","numeric",0,"Start of cdr3 region",
  "cdr3_end","numeric",0,"End of cdr3 region",
  "annotation","character","","V allele annotation"
  )

source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local("Rsub.R")
source_local("SHMHelper.R")

parseArgs("SHMPlot2.R", ARGS, OPTS)

suppressPackageStartupMessages(library(RColorBrewer, quietly=TRUE))
suppressPackageStartupMessages(library(Biostrings, quietly=TRUE, verbose=FALSE))
bases <- getBases()
basecolors <- getBasecolors()
ascii <- getAscii()

data <- read.delim(datafile, header=T, as.is=T)
refseq <- DNAString(paste(data$Base,collapse=""))

if (any(diff(data$Pos) != 1)) stop("Pos column must be sequential")


if (tstart < data$Pos[1]) {
  tstart <- data$Pos[1]
}
if (tend == 0 || tend > tail(data$Pos,n=1)) {
  tend <- tail(data$Pos,n=1)
}

data$Style <- match(data$Base,bases)

data$cov <- data$Total/max(data$Total)


rowwidth <- ceiling((tend-tstart+1)/plotrows)

tstarts <- tstart+rowwidth*0:(plotrows-1)
tends <- tstarts+rowwidth-1


plotline <- data.frame(x=c(data$Pos-0.49,data$Pos+0.49),y=c(data$Y,data$Y))
plotline <- plotline[order(plotline$x),]
if (ymax==0) { ymax <- 1.1*max(plotline$y[plotline$x >= tstart & plotline$x <= tend]) }
refy <- -ymax/20

refseq_rc <- reverseComplement(refseq)

regex1_plot <- data.frame(x=NA,y=NA)
regex2_plot <- data.frame(x=NA,y=NA)

if (regex1 != "") {

  regex1_match <- gregexpr(regex1,as.character(refseq))[[1]]
  regex1_match_rc <- gregexpr(regex1,as.character(refseq_rc))[[1]]

  if (regex1_match[1] > 0) {

    for (i in 1:length(regex1_match)) {
      len <- attr(regex1_match,"match.length")[i]
      pos <- data$Pos[regex1_match[i]]
      regex1_plot <- rbind(regex1_plot,c(pos-0.5,2*refy))
      regex1_plot <- rbind(regex1_plot,plotline[plotline$x >= pos - 0.5 & plotline$x <= pos + len - 0.5, ])
      regex1_plot <- rbind(regex1_plot,c(pos + len - 0.5,2*refy))
      regex1_plot <- rbind(regex1_plot,c(NA,NA))
    }
  }

  if (regex1_match_rc[1] > 0) {
    for (i in 1:length(regex1_match_rc)) {
      len <- attr(regex1_match_rc,"match.length")[i]
      pos <- data$Pos[nrow(data) - regex1_match_rc[i] - len + 2]
      regex1_plot <- rbind(regex1_plot,c(pos-0.5,2*refy))
      regex1_plot <- rbind(regex1_plot,plotline[plotline$x >= pos - 0.5 & plotline$x <= pos + len - 0.5, ])
      regex1_plot <- rbind(regex1_plot,c(pos + len - 0.5,2*refy))
      regex1_plot <- rbind(regex1_plot,c(NA,NA))
    }
  }
}

if (regex2 != "") {

  regex2_match <- gregexpr(regex2,as.character(refseq),perl=T)[[1]]
  regex2_match_rc <- gregexpr(regex2,as.character(refseq_rc),perl=T)[[1]]
  if (regex2_match[1] > 0) {

    for (i in 1:length(regex2_match)) {
      len <- 4   # attr(regex2_match,"match.length")[i]
      pos <- data$Pos[regex2_match[i]]
      regex2_plot <- rbind(regex2_plot,c(pos-0.5,2*refy))
      regex2_plot <- rbind(regex2_plot,plotline[plotline$x >= pos - 0.5 & plotline$x <= pos + len - 0.5, ])
      regex2_plot <- rbind(regex2_plot,c(pos + len - 0.5,2*refy))
      regex2_plot <- rbind(regex2_plot,c(NA,NA))
    }
  }

  if (regex2_match_rc[1] > 0) {
    for (i in 1:length(regex2_match_rc)) {
      len <- 4   # attr(regex2_match_rc,"match.length")[i]
      pos <- data$Pos[nrow(data) - regex2_match_rc[i] - len + 2]
      regex2_plot <- rbind(regex2_plot,c(pos-0.5,2*refy))
      regex2_plot <- rbind(regex2_plot,plotline[plotline$x >= pos - 0.5 & plotline$x <= pos + len - 0.5, ])
      regex2_plot <- rbind(regex2_plot,c(pos + len - 0.5,2*refy))
      regex2_plot <- rbind(regex2_plot,c(NA,NA))
    }
  }
}

pdf(output,height=figureheight,width=11)
  par(mai=c(0.2,0.75,0.5,0.75),omi=c(0.5,0,0,0))
  layout(as.matrix(1:plotrows,ncol=1,nrow=plotrows))

  for (i in 1:plotrows) {
    plot(c(),c(),ylab="",xaxt="n",xlab="",xlim=c(tstarts[i]-0.5,tends[i]+0.5),ylim=c(refy,ymax),xaxs="i",bty="o")
    axis(1,lwd=0,lwd.ticks=1)
    if (nrow(regex2_plot) > 1) polygon(regex2_plot, col=rgb(254,196,79,max=255),border=rgb(254,196,79,max=255))
    if (nrow(regex1_plot) > 1) polygon(regex1_plot, col=rgb(217,95,14,max=255),border=rgb(217,95,14,max=255))
    rect(cdr1_start, -0.5, cdr1_end, ymax+0.2, col = rgb(0,0,0,alpha=0.12), border = "NA")
    rect(cdr2_start, -0.5, cdr2_end, ymax+0.2, col = rgb(0,0,0,alpha=0.12), border = "NA")
    rect(cdr3_start, -0.5, cdr3_end, ymax+0.2, col = rgb(0,0,0,alpha=0.12), border = "NA")

    if ("cov" %in% colnames(data)) {
      covpos = c(rep(ymax+0.07, length(data$cov)))
      topcovline <- data.frame(x=c(data$Pos-0.49,data$Pos+0.49),y=c(covpos+data$cov/5,covpos+data$cov/5))
      topcovline <- topcovline[order(topcovline$x),]

      botcovline <- data.frame(x=c(data$Pos-0.49,data$Pos+0.49),y=c(covpos,covpos))

      botcovline <- botcovline[order(botcovline$x),]
      botcovline$y <- unlist(lapply(botcovline$y,function(y) {max(0,y)}))
      #polygon(c(topcovline$x,rev(botcovline$x)),c(topcovline$y,rev(botcovline$y)),col='black',border=T,xpd=TRUE,density=00)

    }

    if ("Err" %in% colnames(data)) {

      toperrline <- data.frame(x=c(data$Pos-0.49,data$Pos+0.49),y=c(data$Y+data$Err,data$Y+data$Err))
      toperrline <- toperrline[order(toperrline$x),]
      boterrline <- data.frame(x=c(data$Pos-0.49,data$Pos+0.49),y=c(data$Y,data$Y))
      #boterrline <- data.frame(x=c(data$Pos-0.49,data$Pos+0.49),y=c(data$Y-data$Err,data$Y-data$Err))
      boterrline <- boterrline[order(boterrline$x),]
      boterrline$y <- unlist(lapply(boterrline$y,function(y) {max(0,y)}))
      polygon(c(toperrline$x,rev(boterrline$x)),c(toperrline$y,rev(boterrline$y)),col="green",border=F)

    }
    #col=grey(0.5,0.5) --> col="green" on 3/18/2015 --> col=rgb(128,255,0,180)

    grid(col=grey(0.5))
    if (showsequence) points(data$Pos[1]:tail(data$Pos,n=1),rep(refy,nrow(data)),col=basecolors[data$Style],pch=ascii[data$Style],cex=0.6)
    lines(plotline$x,plotline$y)
    if(annotation==''){
        genename = tail(strsplit(datafile, "/")[[1]],1)
        genename = strsplit(genename, ".", fixed=TRUE)[[1]][1]
    }else{
        genename = annotation
    }
    legendwords = paste(genename, " TotalReads=", max(data$Total), sep="")
    legend(240, ymax+0.5, legend=legendwords , xpd=TRUE, cex=.6,bty="n")
  }

graphics.off()
