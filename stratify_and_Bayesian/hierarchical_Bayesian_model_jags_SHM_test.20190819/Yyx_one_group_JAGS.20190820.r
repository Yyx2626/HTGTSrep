## Usage: Rscript this.r <output_prefix> <merged.stat.txt> [pos_range] [colname_prefix]

options(stringsAsFactors=FALSE)



library(rjags)



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



Yyx_one_group_prop = function(group_DF, n_chains=2, iter_n=100000, burnin=60000, plot_var_names="mu", should_trace=TRUE, should_density=TRUE){
	init_list = NULL
	model = jags.model(file="hierarchical_Bayesian_model_jags_SHM_test.20190819/Yyx_HierarPropBinom_JAGS_single.20190820.txt", data=list(y=group_DF$Mut, n=group_DF$Total, sampleNum=nrow(group_DF)), inits=init_list, n.chains=n_chains, quiet=TRUE, n.adapt=burnin)
	model.samples <- coda.samples(model, c("mu", "eta", "theta"), n.iter=iter_n, progress.bar="none")
	if(should_trace || should_density){
		plot(model.samples[, plot_var_names], trace=should_trace, density=should_density)
	}
	model.samples
}



args = commandArgs(TRUE)
if(length(args) < 2){
	stop("Usage: Rscript this.r <output_prefix> <merged.stat.txt> [pos_range] [colname_prefix]")
}

output_prefix = args[1]
input_filename = args[2]

pos_range = NULL
if(length(args) >= 3){
	pos_range = args[3]
	tmp = as.integer(strsplit(pos_range, "[:-]")[[1]])
	if(length(tmp)==1){
		# good, do nothing
	}else if(length(tmp)==2){
		pos_range = tmp[1]:tmp[2]
	}else{
		stop("Cannot recognize pos_range = " %.% pos_range)
	}
}

colname_prefix = ""
if(length(args) >= 4){
	colname_prefix = args[4]
}
colname_prefix_nchar = nchar(colname_prefix)



cat("Now read in input " %.% input_filename %.% " file ...\n")
all_DF = read.delim(input_filename)
echo(dim(all_DF))
NR = nrow(all_DF)
if(is.null(pos_range)){
	pos_range = 1:NR
}
pos_range = pos_range[pos_range >= 1 & pos_range <= NR]

cat("Now model each site ...\n")
in_Mut_colidx     = which(str_left(names(all_DF), colname_prefix_nchar+4) == colname_prefix %.% "Mut.")
in_Total_colidx = which(str_left(names(all_DF), colname_prefix_nchar+6) == colname_prefix %.% "Total.")

pdf(output_prefix %.% "." %.% min(pos_range) %.% "_" %.% max(pos_range) %.% ".trace_density.pdf", width=9, height=9)

all_DF$mu_mean = NA
all_DF$mu_LB = NA
all_DF$mu_UB = NA
for(i in pos_range){
	echo(i)
	group_DF = data.frame(Total=unlist(all_DF[i, in_Total_colidx]), Mut=unlist(all_DF[i, in_Mut_colidx]))

	try(silent=FALSE, {
		mcmc_rlt = Yyx_one_group_prop(group_DF, plot_var_names=c("mu", "eta"))
		write.table(as.matrix(mcmc_rlt), file=output_prefix %.% ".site_" %.% i %.% ".mcmc_rlt.tsv", row.names=FALSE, sep="\t", quote=FALSE)
	})
}

dev.off()

