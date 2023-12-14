#### Usage: modify the Excel file 'yyx_hideRef_ggseqlogo_inupt.xlsx':
####            I only use the second column, and regard the 1st line as reference sequence
####        install required R packages: tidyverse ggplot2 ggseqlogo readxl
####        change the output_date and working_path variables as follows
####        finally run this R script in R or Rscript


output_date = "20230201"   # TODO: change date
working_path = "/your/path/"   # TODO: change working path

setwd(working_path)



options(stringsAsFactors=FALSE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("Biostrings")
#BiocManager::install("BiocStyle")
#BiocManager::install("Biobase")
#BiocManager::install("seqLogo")
#BiocManager::install("ggseqlogo")

library(tidyverse)
library(ggplot2)
library(ggseqlogo)

library(readxl)



`%.%` = function(x,y) paste0(x,y)



### ref: https://github.com/omarwagih/ggseqlogo/blob/master/R/ggseqlogo.r

matrix_to_heights_rmRef = function(mat, seq_type, decreasing = T, ref_seq){
	mat[is.infinite(mat)] = 0
	if (any(duplicated(rownames(mat)))) 
		stop("Matrix input must have unique row names")
	
	ref_seq_vec = strsplit(ref_seq, "")[[1]]
	for(j in 1:length(ref_seq_vec)){   # Yyx add to rmRef, on 2021-04-14
		mat[ref_seq_vec[j], j] = 0
	}
	
	dat = lapply(1:ncol(mat), function(i) {
		vals = mat[, i]
		pos = sort(vals[vals >= 0], decreasing = decreasing)
		neg = sort(vals[vals < 0], decreasing = !decreasing)
		cs_pos = cumsum(pos)
		cs_neg = cumsum(neg)
		df_pos = df_neg = NULL
		if (length(pos) > 0) 
			df_pos = data.frame(letter = names(pos), position = i, 
				y0 = c(0, cs_pos[-length(cs_pos)]), y1 = cs_pos, 
				stringsAsFactors = F)
		if (length(neg) > 0) 
			df_neg = data.frame(letter = names(neg), position = i, 
				y0 = cs_neg, y1 = c(0, cs_neg[-length(cs_neg)]), 
				stringsAsFactors = F)
		rbind(df_pos, df_neg)
	})
	dat = do.call(rbind, dat)
	space_factor = 0.004
	y_pad = max(dat$y1) * space_factor
	dat$y0 = dat$y0 + y_pad
	dat = subset(dat, dat$y1 > dat$y0)
	dummy = data.frame(letter = dat$letter[1], position = NA, 
		y0 = 0, y1 = 0)
	if (dat$position[1] != 1) {
		dummy$position = 1
		dat = rbind(dummy, dat)
	}
	if (dat$position[nrow(dat)] != ncol(mat)) {
		dummy$position = ncol(mat)
		dat = rbind(dat, dummy)
	}
	rownames(dat) = NULL
	attr(dat, "seq_type") = seq_type
	dat
}

bits_method_rmRef = function (seqs, decreasing, ref_seq, ...){
	pfm = makePFM(seqs, ...)
	ic = attr(pfm, "bits")
	if (all(ic == 0)) {
		warning("All positions have zero information content perhaps due to too few input sequences. Setting all information content to 2.")
		ic = (ic * 0) + 2
	}
	heights = t(t(pfm) * ic)
	seq_type = attr(pfm, "seq_type")
	matrix_to_heights_rmRef(heights, seq_type, decreasing, ref_seq=ref_seq)
}

probability_method_rmRef = function (seqs, decreasing, ref_seq, ...){
	pfm = ggseqlogo:::makePFM(seqs, ...)
	seq_type = attr(pfm, "seq_type")
	matrix_to_heights_rmRef(pfm, seq_type, decreasing, ref_seq=ref_seq)
}

# Generate height data for logo
logo_data_rmRef <- function(seqs, ref_seq, method='bits', stack_width=0.95, rev_stack_order=F, font, seq_group=1, seq_type = 'auto', namespace=NULL ){
	# Get font 
	font_df = ggseqlogo:::get_font(font)
	
	# TODO
	# hh = twosamplelogo_method(seqs, seqs_bg, pval_thresh=0.05, seq_type = seq_type, namespace = namespace)
	
	# Generate heights based on method
	if(method == 'bits'){
		hh = bits_method_rmRef(seqs, decreasing = rev_stack_order, seq_type = seq_type, namespace = namespace, ref_seq=ref_seq)
	}else if(method == 'probability'){
		hh = probability_method_rmRef(seqs, decreasing = rev_stack_order, seq_type = seq_type, namespace = namespace, ref_seq=ref_seq)
	}else if(method == 'custom'){
		if(seq_type == 'auto') seq_type = guessSeqType(rownames(seqs))
		hh = matrix_to_heights_rmRef(seqs, seq_type, decreasing = rev_stack_order, ref_seq=ref_seq)
	}else{
		stop('Invalid method!')
	}
	
	# Merge font df and heights
	ff = merge(font_df, hh, by = 'letter')
	# Scale x and ys to new positions
	x_pad = stack_width/2
	ff$x = ggseqlogo:::newRange(ff$x, ff$position - x_pad, ff$position + x_pad)
	ff$y = ggseqlogo:::newRange(ff$y, ff$y0, ff$y1)
	
	# Rename columns
	ff = as.data.frame(ff)[,c('x', 'y', 'letter', 'position', 'order')]
	ff$seq_group = seq_group
	
	# Set sequence type as attribute, to be used downstream
	attr(ff, 'seq_type') = attr(hh, 'seq_type')
	
	# Return data table
	ff
}

geom_hideRef_logo = function (data = NULL, ref_seq, method = "bits", seq_type = "auto", 
	namespace = NULL, font = "roboto_medium", stack_width = 0.95, 
	rev_stack_order = F, col_scheme = "auto", low_col = "black", 
	high_col = "yellow", na_col = "grey20", plot = T, 
	...) 
{
	if (stack_width > 1 | stack_width <= 0) 
		stop("\"stack_width\" must be between 0 and 1")
	if (is.null(data)) 
		stop("Missing \"data\" parameter!")
	if (!is.null(namespace)) 
		seq_type = "other"
	all_methods = c("bits", "probability", "custom")
	pind = pmatch(method, all_methods)
	method = all_methods[pind]
	if (is.na(method)) 
		stop("method must be one of 'bits' or 'probability', or 'custom'")
	if (is.character(data) | is.matrix(data)) 
		data = list(`1` = data)
	if (is.list(data)) {
		if (is.null(names(data))) 
			names(data) = seq_along(data)
		lvls = names(data)
		data_sp = lapply(names(data), function(n) {
			curr_seqs = data[[n]]
			logo_data_rmRef(seqs = curr_seqs, ref_seq=ref_seq, method = method, stack_width = stack_width, 
				rev_stack_order = rev_stack_order, seq_group = n, 
				seq_type = seq_type, font = font, namespace = namespace)
		})
		data = do.call(rbind, data_sp)
		data$seq_group = factor(data$seq_group, levels = lvls)
	}
	if (!plot) 
		return(data)
	seq_type = attr(data, "seq_type")
	cs = ggseqlogo:::get_col_scheme(col_scheme, seq_type)
	legend_title = attr(cs, "cs_label")
	data = merge(data, cs, by = "letter", all.x = T)
	data = data[order(data$order), ]
	colscale_gradient = is.numeric(cs$group)
	colscale_opts = NULL
	if (colscale_gradient) {
		colscale_opts = scale_fill_gradient(name = legend_title, 
			low = low_col, high = high_col, na.value = na_col)
	}
	else {
		tmp = cs[!duplicated(cs$group) & !is.na(cs$group), ]
		col_map = unlist(split(tmp$col, tmp$group))
		colscale_opts = scale_fill_manual(values = col_map, name = legend_title, 
			na.value = na_col)
	}
	guides_opts = NULL
	if (identical(cs$letter, cs$group)) 
		guides_opts = guides(fill = F)
	y_lim = NULL
	extra_opts = NULL
	if (method == "tsl") {
		y_lab = "Depleted	Enriched"
		tmp = max(abs(data$y))
		row_a = row_b = data[1, ]
		row_a$y = -tmp
		row_b$y = tmp
		data = rbind(data, row_a, row_b)
		data$facet = factor(data$y > 0, c(T, F), c("Enriched", 
			"Depleted"))
		extra_opts = NULL
	}
	else if (method == "custom") {
		y_lab = ""
	}
	else {
		y_lab = method
		substr(y_lab, 1, 1) = toupper(substr(y_lab, 1, 1))
	}
	data$group_by = with(data, interaction(seq_group, letter, 
		position))
	data$x = data$x
	logo_layer = layer(stat = "identity", data = data, 
		mapping = aes_string(x = "x", y = "y", fill = "group", 
			group = "group_by"), geom = "polygon", 
		position = "identity", show.legend = NA, inherit.aes = F, 
		params = list(na.rm = T, ...))
	breaks_fun = function(lim) {
		1:floor(lim[2]/1.05)
	}
	list(logo_layer, scale_x_continuous(breaks = breaks_fun, 
		labels = identity), ylab(y_lab), xlab(""), colscale_opts, 
		guides_opts, coord_cartesian(ylim = y_lim), extra_opts)
}


yyx_hideRef_ggseqlogo = function(ref_seq, aligned_seqs, main="", y_max=1, should_print=TRUE, seq_type="auto", x_label_every = 10){
	R_vec = strsplit(ref_seq, "")[[1]]
	A_list = strsplit(aligned_seqs, "")

	for(k in 1:length(A_list)){
		if(length(A_list[[k]]) < length(R_vec)){
			cat("Warning: shorter length " %.% length(A_list[[k]]) %.% " < ref " %.% length(R_vec) %.% " for " %.% k %.% "-th sequence:\n" %.% paste(collapse="", A_list[[k]]) %.% "\n\n")
		}
	}

	stat_DF = character(0)
	for(i in 1:length(R_vec)){
		now_align_vec = rep(NA, length(A_list))
		for(k in 1:length(A_list)){
			if(length(A_list[[k]]) < i){
				cat("Warning: I will skip " %.% i %.% "-th base of " %.% k %.% "-th sequence\n")
			}
			if(A_list[[k]][i] == "-"){
				A_list[[k]][i] = R_vec[i]
			}
			now_align_vec[k] = A_list[[k]][i]
		}
		now_stat_DF = as.data.frame(table(now_align_vec))
		names(now_stat_DF) = c("AA", "Freq")
		now_stat_DF$Prob = now_stat_DF$Freq / sum(now_stat_DF$Freq)
		now_stat_DF$isRef = now_stat_DF$AA == R_vec[i]
		now_stat_DF = now_stat_DF %>% arrange(-Prob, AA)
		now_stat_DF$Pos = i
		stat_DF = rbind(stat_DF, now_stat_DF)
	}

	original_aligned_seqs = sapply(A_list, function(vec) paste(collapse="", vec))

	#ggseqlogo(original_aligned_seqs, seq_type='aa', method="prob")

	rect_DF = data.frame(x=1:length(R_vec))
	rect_DF$xmin = rect_DF$x - 0.5
	rect_DF$xmax = rect_DF$x + 0.5
	rect_DF$ymin = NA
	rect_DF$ymax = NA
	for(i in 1:length(R_vec)){
		now_stat_DF = stat_DF[stat_DF$Pos==i,]
		now_stat_DF$cumProb = cumsum(now_stat_DF$Prob)
		ref_idx = which(now_stat_DF$isRef)
		if(length(ref_idx) <= 0){
			rect_DF$ymax[i] = 1
			rect_DF$ymin[i] = 1
		}else{
			if(ref_idx==1){
				rect_DF$ymax[i] = 1
			}else{
				rect_DF$ymax[i] = 1 - now_stat_DF$cumProb[ref_idx-1]
			}
			rect_DF$ymin[i] = 1 - now_stat_DF$cumProb[ref_idx]
		}
	}


	y_at = seq(0,1,by=0.1)
	p = ggplot()
	p = p + geom_abline(intercept=y_at, slope=0, linetype="dotted")
	p = p + geom_hideRef_logo(original_aligned_seqs, ref_seq, method="prob", seq_type=seq_type) + theme_logo()
#	p = p + geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), data=rect_DF, fill="white")
#	p = p + geom_abline(intercept=y_at, slope=0, linetype="dotted")
	p = p + coord_cartesian(ylim=c(0,y_max)) + scale_y_continuous(breaks=y_at)
	x_at = 1:length(R_vec)
	x_labels = x_at
	x_labels[x_at %% x_label_every != 0] = ""
	p = p + scale_x_continuous(breaks=x_at, labels=x_labels) + theme(axis.ticks.x = element_line())
	
	if(main != ""){
		p = p + ggtitle(main)
	}

	if(should_print){
		print(p)
	}
}


yyx_convert_base_at = function(vec, base_at=28, base_from="I", base_to="S"){
	no_NA_idxes = !is.na(vec)
	if(sum(no_NA_idxes) > 0){
		vec_no_NA = vec[no_NA_idxes]
		strsplit_list = strsplit(vec_no_NA, "")
		vec[no_NA_idxes] = (sapply(strsplit_list, function(seq_array){
			if(seq_array[base_at]==base_from){
				seq_array[base_at] = base_to
			}
			paste(collapse="", seq_array)
		}))
	}
	return (vec)
}




input = read_excel("yyx_hideRef_ggseqlogo_inupt.xlsx", col_names=FALSE)
dim(input)

## translate DNA/RNA to AA sequence, if needed
# library(Biostrings)
# input[[2]] = input[[2]] %>% DNAStringSet %>% translate %>% as.character

ref_seq = input[[2]][1]
aligned_seqs = input[[2]][-1]



pdf("yyx_hideRef_ggseqlogo_output." %.% output_date %.% ".pdf", width=28, height=8)
nchar(ref_seq)   # 95
aligned_seqs = aligned_seqs[!is.na(aligned_seqs) & aligned_seqs != ""]
aligned_seqs = substr(aligned_seqs, 1, nchar(ref_seq))
now_main = "(" %.% length(aligned_seqs) %.% " aligned sequences)"
yyx_hideRef_ggseqlogo(ref_seq, aligned_seqs, main=now_main, seq_type="aa", x_label_every = 10)
try(silent=TRUE, {
yyx_hideRef_ggseqlogo(ref_seq, aligned_seqs, main=now_main, seq_type="dna", x_label_every = 10)
})
dev.off()


