
### ref: yyx_SHMPlot2_Strata_ErrBar.20190729.r

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

#getBases <- function () {
#	bases <- c("A","C","G","T","N")
#}
#getAscii <- function () {
#	ascii <- c(65,67,71,84,78)
#}
#getBasecolors <- function () {
#	basecolors <- brewer.pal(7,"Set1")
#}


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
  "mutMin","numeric",NA,"Minimum mutation number for stratification",
  "mutMax","numeric",NA,"Maximum mutation number for stratification",
  "mutMinProp","numeric",NA,"Minimum mutation number proportion for stratification",
  "mutMaxProp","numeric",NA,"Maximum mutation number proportion for stratification",
  "minReadNumB","numeric",1,"Minimum read number required for each sample before stratification",
  "minReadNumA","numeric",1,"Minimum read number required for each sample after stratification",
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

parseArgs("yyx_mut_rate_per_read.20190904.r", ARGS, OPTS)
cat(paste0("Request to read ", length(nucfile), " nuc.txt files\n"))
for(k in 1:length(nucfile)){
	cat(paste0("\t", k, "\t'", nucfile[k], "'\n"))
}
cat("\n")

output_pdf_filename = output %.% ".pdf"
output_stat_txt_filename = output %.% ".stat.txt"
if(!fo){
#	if(file.exists(output_pdf_filename)){
#		stop(paste0("Output pdf file '", output_pdf_filename, "' already exists.\n\tPlease remove it before run me.\n\tOr you can force output with the option '-fo=T'.\n"))
#	}
	if(file.exists(output_stat_txt_filename)){
		stop(paste0("Output stat.txt file '", output_stat_txt_filename, "' already exists.\n\tPlease remove it before run me.\n\tOr you can force output with the option '-fo=T'.\n"))
	}
}

#suppressPackageStartupMessages(library(RColorBrewer, quietly=TRUE))
#suppressPackageStartupMessages(library(Biostrings, quietly=TRUE, verbose=FALSE))
#bases <- getBases()
#basecolors <- getBasecolors()
#ascii <- getAscii()

cat("\n")
cat(paste0("Now reading guide '", statfile, "' ...\n"))
data <- read.delim(statfile, header=T, as.is=T)
data = data[, c("Pos", "Base")]
data$Pos = as.integer(data$Pos)
data = data[order(data$Pos),]
if (any(diff(data$Pos) != 1)) stop("Pos column must be sequential")
#echo(dim(data))
cat("  input position number = " %.% nrow(data) %.% "\n")

#refseq <- DNAString(paste(data$Base,collapse=""))
#cat(paste0("  initial nrow = ", nrow(data), "\n"))

lowReadTotalColor = hsv(0,1,0.7)
pdf(output_pdf_filename, width=12, height=4)
par(mfrow=c(1,3), mar=c(3,3,3,1)+.1, mgp=c(2,0.8,0))

all_sample_prefixes = character(0)
for(k in 1:length(nucfile)){
	cat("Now reading in " %.% k %.% "-th nuc '" %.% nucfile[k] %.% "' ...\n")
	now_sample_prefix = sub("^.*/", "", nucfile[k])
	now_sample_prefix = sub("[.]nuc.*$", "", now_sample_prefix)
	
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

	## plot mut_rate_per_read distr   ref: yyx_mut_rate_per_read.20190902.r
	eachReadMut = apply(now_data, 1, function(row) sum(row %in% c("A","C","G","T")))
	eachReadWt = apply(now_data, 1, function(row) sum(row == "."))
	eachReadTotal = eachReadMut + eachReadWt
	eachReadMutProp = eachReadMut / eachReadTotal
	
	par(mfrow=c(1,3))
	now_title = "File " %.% k %.% " before stratify : " %.% now_sample_prefix
#	try(silent=TRUE, {
	tmp = hist(eachReadTotal, breaks=seq(-0.5, max(eachReadTotal, na.rm=TRUE)+0.5), xlab="Read length", main=now_title)
	if(sum(eachReadTotal<220.5) > 0){
		hist(eachReadTotal[eachReadTotal<220.5], breaks=seq(-0.5, max(eachReadTotal, na.rm=TRUE)+0.5), xlab="", main="", col=lowReadTotalColor, add=TRUE)
	}
	mtext(nrow(now_data) %.% " reads", cex=0.8)
	abline(v=220.5, col=lowReadTotalColor, lty=2)
	text(220.5, max(tmp$counts), adj=c(1.1, 1), "< 220.5", col=lowReadTotalColor)
	hist(eachReadMut, breaks=seq(-0.5, max(eachReadMut, na.rm=TRUE)+0.5), xlim=c(0,30), xlab="Mutation number per read", main=now_title)
	if(sum(eachReadTotal<220.5) > 0){
		hist(eachReadMut[eachReadTotal<220.5], breaks=seq(-0.5, max(eachReadMut, na.rm=TRUE)+0.5), col=lowReadTotalColor, xlim=c(0,20), xlab="", main="", add=TRUE)
	}
	plot(density(eachReadMutProp, bw=0.003), xlim=c(0,30/300), xlab="Mutation proportion per read", main=now_title)
	if(sum(eachReadTotal<220.5) > 0){
		lines(density(eachReadMutProp[eachReadTotal<220.5], bw=0.003), col=lowReadTotalColor, lty=2)
	}
	rug(eachReadMutProp[eachReadTotal>220.5])
	rug(eachReadMutProp[eachReadTotal<220.5], col=lowReadTotalColor)
#	})
	
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
	
	## plot mut_rate_per_read distr   ref: yyx_mut_rate_per_read.20190902.r
	eachReadMut = apply(now_data, 1, function(row) sum(row %in% c("A","C","G","T")))
	eachReadWt = apply(now_data, 1, function(row) sum(row == "."))
	eachReadTotal = eachReadMut + eachReadWt
	eachReadMutProp = eachReadMut / eachReadTotal
	
	par(mfrow=c(1,3))
	now_title = "File " %.% k %.% " after stratify : " %.% now_sample_prefix
#	try(silent=TRUE, {
	tmp = hist(eachReadTotal, breaks=seq(-0.5, max(eachReadTotal, na.rm=TRUE)+0.5), xlab="Read length", main=now_title)
	if(sum(eachReadTotal<220.5) > 0){
		hist(eachReadTotal[eachReadTotal<220.5], breaks=seq(-0.5, max(eachReadTotal, na.rm=TRUE)+0.5), xlab="", main="", col=lowReadTotalColor, add=TRUE)
	}
	mtext(nrow(now_data) %.% " reads", cex=0.8)
	abline(v=220.5, col=lowReadTotalColor, lty=2)
	text(220.5, max(tmp$counts), adj=c(1.1, 1), "< 220.5", col=lowReadTotalColor)
	hist(eachReadMut, breaks=seq(-0.5, max(eachReadMut, na.rm=TRUE)+0.5), xlim=c(0,30), xlab="Mutation number per read", main=now_title)
	if(sum(eachReadTotal<220.5) > 0){
		hist(eachReadMut[eachReadTotal<220.5], breaks=seq(-0.5, max(eachReadMut, na.rm=TRUE)+0.5), col=lowReadTotalColor, xlim=c(0,20), xlab="", main="", add=TRUE)
	}
	plot(density(eachReadMutProp, bw=0.003), xlim=c(0,30/300), xlab="Mutation proportion per read", main=now_title)
	if(sum(eachReadTotal<220.5) > 0){
		lines(density(eachReadMutProp[eachReadTotal<220.5], bw=0.003), col=lowReadTotalColor, lty=2)
	}
	rug(eachReadMutProp[eachReadTotal>220.5])
	rug(eachReadMutProp[eachReadTotal<220.5], col=lowReadTotalColor)
#	})
	
	
	eachPosMut = apply(now_data, 2, function(col) sum(col %in% c("A","C","G","T")))
	eachPosWt = apply(now_data, 2, function(col) sum(col == "."))
	eachPosTotal = eachPosMut + eachPosWt
#	echo(length(eachPosTotal))
	
#	data[["Mut." %.% k]] = eachPosMut
	if(now_sample_prefix %in% all_sample_prefixes){
		data[["Mut." %.% now_sample_prefix %.% "." %.% (sum(all_sample_prefixes==now_sample_prefix)+1)]] = eachPosMut
		data[["Total." %.% now_sample_prefix %.% "." %.% (sum(all_sample_prefixes==now_sample_prefix)+1)]] = eachPosTotal
		data[["Y." %.% now_sample_prefix %.% "." %.% (sum(all_sample_prefixes==now_sample_prefix)+1)]] = yyx_convert_NA(eachPosMut / eachPosTotal, 0)
	}else{
		data[["Mut." %.% now_sample_prefix]] = eachPosMut
		data[["Total." %.% now_sample_prefix]] = eachPosTotal
		data[["Y." %.% now_sample_prefix]] = yyx_convert_NA(eachPosMut / eachPosTotal, 0)
	}
	all_sample_prefixes = c(all_sample_prefixes, now_sample_prefix)
#	cat(paste0("  updated nrow = ", nrow(data), "\n"))
#	each_total_Y_colidx = c(each_total_Y_colidx, ncol(data)-1, ncol(data))
}

dev.off()

#each_total_Y_colidx = c((1:length(nucfile))*3, (1:length(nucfile))*3-1)
#each_total_Y_colidx = c((1:length(nucfile))*2, (1:length(nucfile))*2-1)
#each_total_Y_colidx = sort(each_total_Y_colidx) + 2

cat("\nIn total, there are " %.% length(all_sample_prefixes) %.% " valid samples: " %.% join(", ", all_sample_prefixes) %.% "\n")

if(ncol(data) <= 2){   # no valid samples, maybe all filtered by minReadNum
	cat("Note: because there seem to be no valid samples, so I just exit this R script\n")
	q("no")
}

#### reference: intrinsicMutProfileErrBar.20190131.py
#data$Total = apply(data[,each_total_Y_colidx], 1, function(row){
#	total_vec = row[seq(1, length(row), by=2)]
#	sum(total_vec)
#})
#data$Y = apply(data[,each_total_Y_colidx], 1, function(row){
#	total_vec = row[seq(1, length(row), by=2)]
#	Y_vec = row[seq(2, length(row), by=2)]
#	Y_vec_2 = Y_vec[total_vec > 0]
#	mean(Y_vec_2)
#})
#data$Err = apply(data[,each_total_Y_colidx], 1, function(row){
#	total_vec = row[seq(1, length(row), by=2)]
#	Y_vec = row[seq(2, length(row), by=2)]
#	Y_vec_2 = Y_vec[total_vec > 0]
#	sd(Y_vec_2)/sqrt(length(Y_vec_2))
#})

cat("\n")
cat(paste0("Now output stat.txt to '", output_stat_txt_filename, "' ...\n"))
write.table(data, file=output_stat_txt_filename, sep="\t", quote=FALSE, row.names=FALSE)

cat("\nAll done. Congratulations!\n")
