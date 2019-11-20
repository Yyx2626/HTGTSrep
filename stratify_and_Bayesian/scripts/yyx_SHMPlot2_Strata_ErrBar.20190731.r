
#### 2019-07-31, Huan asked me to create a folder for output, to contain the .pdf, .stat.txt, .nuc.txt files for used samples

#### 2019-07-29, Huan asked me to add sample filter on minimum read number

#### 2019-07-26, Huan asked me to re-implement Strata with .nuc.txt files as input

#### 2019-07-17, Huan asked me to also output stat.txt file of the averaged SHM and error bar


# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("Biostrings", version = "3.8")

options(stringsAsFactors=FALSE)
#options(stringsAsFactors=FALSE, error=traceback)

### reference: yyx_echo_string_functions.20160420.r

`%.%` <- function(x,y) paste0(x,y)

join = function(sep, vec, ...) paste(collapse=sep, c(vec, ...))
cat0 = function(...) cat(sep="", ...)

echo_str <- function(x, sep=" =\t", collapse=", ") deparse(substitute(x)) %.% sep %.% join(collapse, x)
echo <- function(x, sep=" =\t", collapse=", ") cat0(deparse(substitute(x)), sep, join(collapse, x), "\n");


str_trim_left = function(str, trimLen){
	strlen = nchar(str)
	start = pmin(strlen, trimLen) + 1
	end = strlen
	substr(str, start, end)
}
str_left = function(str, len){
	strlen = nchar(str)
	start = 1
	end = pmin(len, strlen)
	substr(str, start, end)
}


### ref: yyx_rm_NA.yyx_convert_NA.20141114.r

yyx_NA2FALSE <- function(vec){
	vec[is.na(vec)] <- FALSE
	vec
}

yyx_convert_NA <- function(vec, default=FALSE){
	vec[is.na(vec)] <- default
	vec
}


### reference: SHMHelper.R

getBases <- function () {
	bases <- c("A","C","G","T","N")
}
getAscii <- function () {
	ascii <- c(65,67,71,84,78)
}
getBasecolors <- function () {
	basecolors <- brewer.pal(7,"Set1")
}


### reference: Rsub.R

parseArgs <- function(scriptname, ARGS, OPTS=NULL) {
	if (length(ARGS) %% 3 != 0)
		stop("Error: ARGS must have length divisable by 3")
	if (!is.null(OPTS) && length(OPTS) %% 4 != 0)
		stop("Error: OPTS must have length divisable by 4")
	ARGS <- matrix(ARGS, nrow=length(ARGS)/3, ncol=3, byrow=T, dimnames=list(c(),c("name","type","description")))
	OPTS <- if (!is.null(OPTS)) matrix(OPTS, nrow=length(OPTS)/4, ncol=4, byrow=T, dimnames=list(c(),c("name","type","default","description")))

	
	ARGV <- commandArgs(trailingOnly=TRUE)

	catargs <- paste0("\t",ARGS[,"name"],"\t",ARGS[,"description"], collapse="\n")
	catopts <- if (!is.null(OPTS)) paste0("\t-",OPTS[,"name"],"=",OPTS[,"default"],"\t", OPTS[,"description"], collapse="\n")

	usage <- paste("\nUsage: Rscript ",scriptname," [OPTS] ",paste(ARGS[,"name"],collapse=" "),"\n\nOptions(=defaults):\n", catopts,"\n\nArguments:\n", catargs,"\n\n",sep="")

	types <- c("character","numeric","integer","logical")

	## initialize OPTS
	if (!is.null(OPTS)) for (k in 1:nrow(OPTS)) {
		if (!OPTS[k,"type"] %in% types)
			stop("Unrecognized argument type ",OPTS[k,"type"]," for option ",OPTS[k,"name"],"\n",usage)
		cmd <- paste(OPTS[k,"name"],"<<-as.",OPTS[k,"type"],"(\"",OPTS[k,"default"],"\")",sep="")
		if(is.na(OPTS[k,"default"]) && OPTS[k,"type"]=="numeric"){
			cmd <- OPTS[k,"name"] %.% "<<-NA"
		}
#		cat(paste0("[CMD] ", cmd, "\n"))
		eval(parse(text=cmd))
	}

	if (length(ARGV) < nrow(ARGS)){
		cat(usage)
		stop("Not enough arguments.\n")   # Yyx note 2019-07-29: stop() may has limit on the length of input string
	}
	
	## parse OPTS
	k = 1
	while(str_left(ARGV[k], 1)=="-"){
		now_arg = str_trim_left(ARGV[k], 1)
		now_regexec = regexec("(.*)=(.*)", now_arg)
#		echo(now_arg)
#		print(now_regexec)
		if(now_regexec[[1]][1] >= 0){
			now_matches = regmatches(now_arg, now_regexec)
			now_name = now_matches[[1]][2]
			now_value = now_matches[[1]][3]
			if(now_name %in% OPTS[, "name"]){
				idx <- match(now_name,OPTS[,"name"])
				cmd <- paste(now_name,"<<-as.",OPTS[idx,"type"],"(\"",now_value,"\")",sep="")
#				cat(paste0("[CMD] ", cmd, "\n"))
				eval(parse(text=cmd))
			}else{
				stop("unrecognized option ",now_name," of ",ARGV[k],"\n",usage)
			}
		}else{
			now_name = now_arg
			if(now_name %in% OPTS[, "name"]){
				idx <- match(now_name,OPTS[,"name"])
				if(OPTS[idx,"type"]=="logical"){
					cmd <- paste(OPTS[idx,"name"],"<<-!as.",OPTS[idx,"type"],"(\"",OPTS[idx,"default"],"\")",sep="")
#					cat(paste0("[CMD] ", cmd, "\n"))
					eval(parse(text=cmd))
				}
			}else{
#				stop("option ",now_name," (",ARGV[k],") is not recognized\n",usage)
				stop("unrecognized option or opt=value expression expected instead of ",ARGV[k],"\n",usage)
			}
			
#			cat(paste0("Warning: cannot recognize -", ARGV[k], " as -opt=value\n"))
#			stop("opt=value expression expected instead of ",ARGV[k],"\n",usage)
		}
		k = k + 1
	}
	
	if (length(ARGV) - k + 1 < nrow(ARGS))
		stop("Not enough arguments.\n",usage)


	for (i in 1:nrow(ARGS)) {
		if (! ARGS[i,"type"] %in% types)
			stop("Unrecognized argument type ",ARGS[i,"type"]," for argument ",ARGS[i,"name"],"\n",usage)
		cmd <- paste(ARGS[i,"name"],"<<-as.",ARGS[i,"type"],"(\"",ARGV[k],"\")",sep="")
#		cat(paste0("[CMD] ", cmd, "\n"))
		eval(parse(text=cmd))
		k = k + 1
	}
	
	if (length(ARGV) - k + 1 > 0){
		i = nrow(ARGS)
		tail_ARGV_str = paste0(collapse="\",\"", ARGV[k:length(ARGV)])
		cmd <- paste(ARGS[i,"name"],"<<-c(",ARGS[i,"name"],",as.",ARGS[i,"type"],"(c(\"",tail_ARGV_str,"\")))",sep="")
#		cat(paste0("[CMD] ", cmd, "\n"))
		eval(parse(text=cmd))
	}

#	if (length(ARGV) > i) for (j in (i+1):length(ARGV)) {
#		opt <- unlist(strsplit(ARGV[j],"="))
#		if (length(opt) != 2)
#			stop("opt=value expression expected instead of ",ARGV[j],"\n",usage)
#		if (is.null(OPTS) || ! opt[1] %in% OPTS[,"name"])
#			stop("optional argument ",opt[1]," not recognized\n",usage)
#		idx <- match(opt[1],OPTS[,"name"])
#		cmd <- paste(opt[1],"<<-as.",OPTS[idx,"type"],"(\"",opt[2],"\")",sep="")
#		eval(parse(text=cmd))
#	}
}






### reference: SHMPlot2.noCov.R

ARGS <- c(
  "output","character","file path prefix for output",
  "statfile","character","file path of one guide (same V) .stat.txt file - 1~2th columns: Pos, Base",
  "nucfile","character","file path of .nuc.txt file(s) - columns: Read_ID, base [ACGTN-.] on each position"
)

OPTS <- c(
  "tstart","numeric",0,"Start of reference to view",
  "tend","numeric",0,"End of reference to view",
  "plotrows","numeric",1,"Rows on plot",
  "ymax","numeric",0.75,"Maximum y-axis height",
  "figureheight","numeric",2,"height in inches",
  "showsequence","logical",FALSE,"display sequence on plots",
  "regex1","character","AGCT","",
  "regex2","character","[AGT](?=G[CT][AT])", "",   # this is DGYW/WRCH motif
  "cdr1_start","numeric",0,"Start of cdr1 region",
  "cdr1_end","numeric",0,"End of cdr1 region",
  "cdr2_start","numeric",0,"Start of cdr2 region",
  "cdr2_end","numeric",0,"End of cdr2 region",
  "cdr3_start","numeric",0,"Start of cdr3 region",
  "cdr3_end","numeric",0,"End of cdr3 region",
  "annotation","character","","V allele annotation (default: empty for no annotation)",
  "mutMin","numeric",NA,"Minimum mutation number for stratification",
  "mutMax","numeric",NA,"Maximum mutation number for stratification",
  "mutMinProp","numeric",NA,"Minimum mutation number proportion for stratification",
  "mutMaxProp","numeric",NA,"Maximum mutation number proportion for stratification",
  "minReadNumB","numeric",0,"Minimum read number required for each sample before stratification",
  "minReadNumA","numeric",0,"Minimum read number required for each sample after stratification",
  "fo","logical",FALSE,"force output (default: stop if output file exists)"
  )

#source_local <- function(fname){
#  argv <- commandArgs(trailingOnly = FALSE)
#  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
#  source(paste(base_dir, fname, sep="/"))
#}
#
#source_local("Rsub.R")
#source_local("SHMHelper.R")

parseArgs("yyx_SHMPlot2_Strata_ErrBar.20190731.r", ARGS, OPTS)
cat(paste0("Request to read ", length(nucfile), " nuc.txt files\n"))
for(k in 1:length(nucfile)){
	cat(paste0("\t", k, "\t'", nucfile[k], "'\n"))
}
cat("\n")

output_pdf_filename = output %.% ".pdf"
output_stat_txt_filename = output %.% ".stat.txt"
output_used_nuc_list_filename = output %.% ".used_nuc_list.tsv"
if(!fo){
	if(file.exists(output_pdf_filename)){
		stop(paste0("Output pdf file '", output_pdf_filename, "' already exists.\n\tPlease remove it before run me.\n\tOr you can force output with the option '-fo=T'.\n"))
	}
	if(file.exists(output_stat_txt_filename)){
		stop(paste0("Output stat.txt file '", output_stat_txt_filename, "' already exists.\n\tPlease remove it before run me.\n\tOr you can force output with the option '-fo=T'.\n"))
	}
}

suppressPackageStartupMessages(library(RColorBrewer, quietly=TRUE))
suppressPackageStartupMessages(library(Biostrings, quietly=TRUE, verbose=FALSE))
bases <- getBases()
basecolors <- getBasecolors()
ascii <- getAscii()

cat("\n")
cat(paste0("Now reading guide '", statfile, "' ...\n"))
data <- read.delim(statfile, header=T, as.is=T)
data = data[, c("Pos", "Base")]
data$Pos = as.integer(data$Pos)
data = data[order(data$Pos),]
if (any(diff(data$Pos) != 1)) stop("Pos column must be sequential")
#echo(dim(data))
cat("  input gene number = " %.% nrow(data) %.% "\n")

refseq <- DNAString(paste(data$Base,collapse=""))
#cat(paste0("  initial nrow = ", nrow(data), "\n"))

each_total_Y_colidx = integer(0)
all_sample_prefixes = character(0)
used_nuc_files = character(0)
for(k in 1:length(nucfile)){
	cat("Now reading in " %.% k %.% "-th nuc '" %.% nucfile[k] %.% "' ...\n")
	now_sample_prefix = sub("^.*/", "", nucfile[k])
	now_sample_prefix = sub("[.].*$", "", now_sample_prefix)
	
	now_data <- read.delim(nucfile[k], header=T, as.is=T)
#	echo(dim(now_data))
	rownames(now_data) = now_data[,1]
	now_data = now_data[,-1]
#	echo(dim(now_data))
	cat("  original dim = " %.% join(", ", dim(now_data)) %.% "\n")
	if(ncol(now_data) != nrow(data)){
		cat("Warning: .nuc.txt pos (ncol - 1) is not equal to .stat.txt pos (nrow), so I skip this .nuc.txt file\n")
		next
	}
#	echo(dim(now_data))
	
	## 2019-07-29, add sample filter on minReadNum
	if(nrow(now_data) < minReadNumB){
		cat("Warning: .nuc.txt read number (nrow=" %.% nrow(now_data) %.% ") is less than minReadNumB (" %.% minReadNumB %.% "), so I skip this .nuc.txt file\n")
		next
	}
	
	## 2019-07-26, add stratification
	### ref: intrinsicMutProfileStrata.20190726.py
	eachReadMut = apply(now_data, 1, function(row) sum(row %in% c("A","C","G","T")))
	eachReadWt = apply(now_data, 1, function(row) sum(row == "."))
	eachReadTotal = eachReadMut + eachReadWt
	eachReadMutProp = eachReadMut / eachReadTotal
#	echo(length(eachReadTotal))
	
	eachRead_boolIdx = rep(TRUE, nrow(now_data))
	if(!is.na(mutMin)){
		eachRead_boolIdx = eachRead_boolIdx & yyx_NA2FALSE(eachReadMut >= mutMin)
	}
	if(!is.na(mutMax)){
		eachRead_boolIdx = eachRead_boolIdx & yyx_NA2FALSE(eachReadMut <= mutMax)
	}
	if(!is.na(mutMinProp)){
		eachRead_boolIdx = eachRead_boolIdx & yyx_NA2FALSE(eachReadMutProp >= mutMinProp)
	}
	if(!is.na(mutMaxProp)){
		eachRead_boolIdx = eachRead_boolIdx & yyx_NA2FALSE(eachReadMutProp <= mutMaxProp)
	}
	now_data = now_data[eachRead_boolIdx, ]
#	echo(dim(now_data))
	cat("  stratified dim = " %.% join(", ", dim(now_data)) %.% "\n")
	
	## 2019-07-29, add sample filter on minReadNum
	if(nrow(now_data) < minReadNumA){
		cat("Warning: .nuc.txt read number (nrow=" %.% nrow(now_data) %.% ") is less than minReadNumA (" %.% minReadNumA %.% "), so I skip this .nuc.txt file\n")
		next
	}
	
	eachPosMut = apply(now_data, 2, function(col) sum(col %in% c("A","C","G","T")))
	eachPosWt = apply(now_data, 2, function(col) sum(col == "."))
	eachPosTotal = eachPosMut + eachPosWt
#	echo(length(eachPosTotal))
	
#	data[["Mut." %.% k]] = eachPosMut
	if(now_sample_prefix %in% all_sample_prefixes){
		data[["Total." %.% now_sample_prefix %.% "." %.% (sum(all_sample_prefixes==now_sample_prefix)+1)]] = eachPosTotal
		data[["Y." %.% now_sample_prefix %.% "." %.% (sum(all_sample_prefixes==now_sample_prefix)+1)]] = yyx_convert_NA(eachPosMut / eachPosTotal, 0)
	}else{
		data[["Total." %.% now_sample_prefix]] = eachPosTotal
		data[["Y." %.% now_sample_prefix]] = yyx_convert_NA(eachPosMut / eachPosTotal, 0)
	}
	all_sample_prefixes = c(all_sample_prefixes, now_sample_prefix)
	used_nuc_files = c(used_nuc_files, nucfile[k])
#	cat(paste0("  updated nrow = ", nrow(data), "\n"))
	each_total_Y_colidx = c(each_total_Y_colidx, ncol(data)-1, ncol(data))
}

#each_total_Y_colidx = c((1:length(nucfile))*3, (1:length(nucfile))*3-1)
#each_total_Y_colidx = c((1:length(nucfile))*2, (1:length(nucfile))*2-1)
#each_total_Y_colidx = sort(each_total_Y_colidx) + 2

cat("\nIn total, there are " %.% length(all_sample_prefixes) %.% " valid samples: " %.% join(", ", all_sample_prefixes) %.% "\n")

if(ncol(data) <= 2){   # no valid samples, maybe all filtered by minReadNum
	cat("Note: because there seem to be no valid samples, so I just exit this R script\n")
	q("no")
}
cat("\n")
cat(paste0("Now output used_nuc_list.tsv to '", output_used_nuc_list_filename, "' ...\n"))
write.table(data.frame(sample_name=all_sample_prefixes, nuc_path=used_nuc_files), file=output_used_nuc_list_filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

### reference: intrinsicMutProfileErrBar.20190131.py
data$Total = apply(data[,each_total_Y_colidx], 1, function(row){
	total_vec = row[seq(1, length(row), by=2)]
	sum(total_vec)
})
data$Y = apply(data[,each_total_Y_colidx], 1, function(row){
	total_vec = row[seq(1, length(row), by=2)]
	Y_vec = row[seq(2, length(row), by=2)]
	Y_vec_2 = Y_vec[total_vec > 0]
	mean(Y_vec_2)
})
data$Err = apply(data[,each_total_Y_colidx], 1, function(row){
	total_vec = row[seq(1, length(row), by=2)]
	Y_vec = row[seq(2, length(row), by=2)]
	Y_vec_2 = Y_vec[total_vec > 0]
	sd(Y_vec_2)/sqrt(length(Y_vec_2))
})

cat("\n")
cat(paste0("Now output stat.txt to '", output_stat_txt_filename, "' ...\n"))
write.table(data, file=output_stat_txt_filename, sep="\t", quote=FALSE, row.names=FALSE)

data$Y[is.na(data$Y)] = 0
data$Err[is.na(data$Err)] = 0

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
			len <- 4	 # attr(regex2_match,"match.length")[i]
			pos <- data$Pos[regex2_match[i]]
			regex2_plot <- rbind(regex2_plot,c(pos-0.5,2*refy))
			regex2_plot <- rbind(regex2_plot,plotline[plotline$x >= pos - 0.5 & plotline$x <= pos + len - 0.5, ])
			regex2_plot <- rbind(regex2_plot,c(pos + len - 0.5,2*refy))
			regex2_plot <- rbind(regex2_plot,c(NA,NA))
		}
	}

	if (regex2_match_rc[1] > 0) {
		for (i in 1:length(regex2_match_rc)) {
			len <- 4	 # attr(regex2_match_rc,"match.length")[i]
			pos <- data$Pos[nrow(data) - regex2_match_rc[i] - len + 2]
			regex2_plot <- rbind(regex2_plot,c(pos-0.5,2*refy))
			regex2_plot <- rbind(regex2_plot,plotline[plotline$x >= pos - 0.5 & plotline$x <= pos + len - 0.5, ])
			regex2_plot <- rbind(regex2_plot,c(pos + len - 0.5,2*refy))
			regex2_plot <- rbind(regex2_plot,c(NA,NA))
		}
	}
}

cat("\n")
cat(paste0("Now output pdf to '", output_pdf_filename, "' ...\n"))
pdf(output_pdf_filename, height=figureheight, width=11)
	par(mai=c(0.2,0.75,0.5,0.75),omi=c(0.5,0,0,0))
	if(annotation == ""){
		par(mai=c(0.2,0.75,0.2,0.75))
	}
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
#		if(annotation==''){
#				genename = tail(strsplit(datafile, "/")[[1]],1)
#				genename = strsplit(genename, ".", fixed=TRUE)[[1]][1]
#		}else{
#				genename = annotation
#		}
#		legendwords = paste(genename, " TotalReads=", max(data$Total), sep="")
#		legend(240, ymax+0.5, legend=legendwords , xpd=TRUE, cex=.6,bty="n")
		if(annotation != ""){
			legendwords = paste(annotation, " TotalReads=", max(data$Total), sep="")
			legend(240, ymax+0.5, legend=legendwords , xpd=TRUE, cex=.6,bty="n")
		}
	}

graphics.off()

cat("\nAll done. Congratulations!\n")
