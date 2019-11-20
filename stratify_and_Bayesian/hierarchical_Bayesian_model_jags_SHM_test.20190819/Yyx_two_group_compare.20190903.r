## Usage: Rscript this.r <output_prefix> <in1_foreground_prefix> <in2_background_prefix> [pos_range]

options(stringsAsFactors=FALSE)



library(HDInterval)



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
str_trim_right = function(str, trimLen){
	strlen = nchar(str)
	start = 1
	end = pmax(0, strlen - trimLen)
	substr(str, start, end)
}
str_left = function(str, len){
	strlen = nchar(str)
	start = 1
	end = pmin(len, strlen)
	substr(str, start, end)
}
str_right = function(str, len){
	strlen = nchar(str)
	start = pmax(1, strlen - len + 1)
	end = strlen
	substr(str, start, end)
}



### ref: http://varianceexplained.org/r/bayesian_fdr_baseball/
cummean = function(vec){
	cumsum(vec) / (1:length(vec))
}
PEP_to_qvalues = function(PEP_vec){
	order_vec = order(PEP_vec)
	ans = rep(NA, length(PEP_vec))
	PEP_ordered = PEP_vec[order_vec]
	qvalue_ordered = cummean(PEP_ordered)
	ans[order_vec] = qvalue_ordered
	ans
}



if(interactive()==FALSE){
	args = commandArgs(TRUE)
}
if(length(args) < 3){
	stop("Usage: Rscript this.r <output_prefix> <in1_foreground_prefix> <in2_background_prefix> [pos_range]")
}

output_prefix = args[1]
input_foreground_prefix = args[2]
input_background_prefix = args[3]

pos_range = NULL
if(length(args) >= 4){
	pos_range = args[4]
	tmp = as.integer(strsplit(pos_range, "[:-]")[[1]])
	if(length(tmp)==1){
		# good, do nothing
		# 2019-08-21, Yyx debug, pos_range will be character without converting
		pos_range = as.integer(tmp[1])
	}else if(length(tmp)==2){
		pos_range = tmp[1]:tmp[2]
	}else{
		stop("Cannot recognize pos_range = " %.% pos_range)
	}
}



input_foreground_stat_DF = read.delim(input_foreground_prefix %.% ".all.stat.txt")
echo(dim(input_foreground_stat_DF))
input_background_stat_DF = read.delim(input_background_prefix %.% ".all.stat.txt")
echo(dim(input_background_stat_DF))
stopifnot(nrow(input_foreground_stat_DF)==nrow(input_background_stat_DF))

names(input_foreground_stat_DF) = "in1." %.% names(input_foreground_stat_DF)
names(input_background_stat_DF) = "in2." %.% names(input_background_stat_DF)   # there is a bug in previous version, and debugged on 2019-08-29
input_DF = cbind(input_foreground_stat_DF, input_background_stat_DF)
echo(dim(input_DF))

input_DF$in1.simu_num = NA
input_DF$in1.mu_mean = NA
input_DF$in1.mu_LB = NA
input_DF$in1.mu_UB = NA
input_DF$in2.simu_num = NA
input_DF$in2.mu_mean = NA
input_DF$in2.mu_LB = NA
input_DF$in2.mu_UB = NA
input_DF$muDiff.simu_num = NA
input_DF$muDiff_mean = NA
input_DF$muDiff_LB = NA
input_DF$muDiff_UB = NA
input_DF$muDiff_gt0_posterior = NA
input_DF$muDiff_gt0_01_posterior = NA
input_DF$muDiff_gt0_05_posterior = NA
input_DF$muDiff_gt0_1_posterior = NA

real_pos_range = 1:nrow(input_DF)
if(!is.null(pos_range)){
	real_pos_range = pos_range[pos_range %in% 1:nrow(input_DF)]
}

for(i in real_pos_range){
	echo(i)
	try({
	mcmc_1 = read.delim(input_foreground_prefix %.% ".site_" %.% i %.% ".mcmc_rlt.tsv")
	input_DF$in1.simu_num[i] = nrow(mcmc_1)
	input_DF$in1.mu_mean[i] = mean(mcmc_1$mu)
	tmp = hdi(mcmc_1$mu)
	input_DF$in1.mu_LB[i] = tmp[1]
	input_DF$in1.mu_UB[i] = tmp[2]
	})
	
	try({
	mcmc_2 = read.delim(input_background_prefix %.% ".site_" %.% i %.% ".mcmc_rlt.tsv")
	input_DF$in2.simu_num[i] = nrow(mcmc_2)
	input_DF$in2.mu_mean[i] = mean(mcmc_2$mu)
	tmp = hdi(mcmc_2$mu)
	input_DF$in2.mu_LB[i] = tmp[1]
	input_DF$in2.mu_UB[i] = tmp[2]
	})
	
	try({
	NS = min(c(nrow(mcmc_1), nrow(mcmc_2)))
	muDiff_vec = mcmc_1$mu[1:NS] - mcmc_2$mu[1:NS]
	input_DF$muDiff.simu_num[i] = NS
	input_DF$muDiff_mean[i] = mean(muDiff_vec)
	tmp = hdi(muDiff_vec)
	input_DF$muDiff_LB[i] = tmp[1]
	input_DF$muDiff_UB[i] = tmp[2]
	input_DF$muDiff_gt0_posterior[i] = mean(muDiff_vec > 0)
	input_DF$muDiff_gt0_01_posterior[i] = mean(muDiff_vec > 0.01)
	input_DF$muDiff_gt0_05_posterior[i] = mean(muDiff_vec > 0.05)
	input_DF$muDiff_gt0_1_posterior[i] = mean(muDiff_vec > 0.1)
	})
}
write.table(input_DF[real_pos_range,], file=output_prefix %.% ".final_" %.% min(pos_range) %.% "_" %.% max(pos_range) %.% ".stat.txt", row.names=FALSE, sep="\t", quote=FALSE)

for(gt_str in c("0", "0_01", "0_05", "0_1")){
	input_DF[["muDiff_gt" %.% gt_str %.% "_PEP"]] = 1 - input_DF[["muDiff_gt" %.% gt_str %.% "_posterior"]]
	input_DF[["muDiff_gt" %.% gt_str %.% "_FDR"]] = PEP_to_qvalues(input_DF[["muDiff_gt" %.% gt_str %.% "_PEP"]])
}
write.table(input_DF[real_pos_range,], file=output_prefix %.% ".final_" %.% min(pos_range) %.% "_" %.% max(pos_range) %.% ".stat.txt", row.names=FALSE, sep="\t", quote=FALSE)

