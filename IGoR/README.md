# IGoR
Some scripts processing IGoR results.

Author: Adam Yongxin Ye @ Boston Children's Hospital / Harvard Medical School

## Installation IGoR

Please refer to [IGoR](https://github.com/qmarcou/IGoR) ([original paper](https://www.nature.com/articles/s41467-018-02832-w)).



## Run IGoR

### Run IGoR-align

#### Prepare "genome" files

Most importantly, prepare three fasta files containing the V/D/J genes/alleles:
* prepared_for_IGoR/reference/genomicVs.fasta
* prepared_for_IGoR/reference/genomicDs.fasta
* prepared_for_IGoR/reference/genomicJs.fasta

#### Prepare input sequence file(s)

For example, use SEQUENCE_INPUT or SEQUENCE_VDJ in .db files of HTGTSrep pipeline (supposing if it is 58-th column)

```bash
cat naiveB.all.db | python3 scripts/yyx_show_or_skip_or_retrieve_columns.20190128.py show

cut -f58 naiveB.all.db | tail -n+2 | perl -ne 'BEGIN{ print "seq_index;sequence\n"; $NR=0; } print join(";", $NR++, $_); ' >naiveB_indexed_sequences.csv
```

#### Run IGoR-align

Make a working sub-directory (e.g. run_igor_20190913), then run IGoR-align for V/D/J gene alignment separately.
According to my practical experience of IGoR program, I set 10 threads for V alignment and only one thread for each of D and J alignment.

```bash
mkdir run_igor_20190913
mkdir run_igor_20190913/aligns
ln -s naiveB_indexed_sequences.csv  run_igor_20190913/aligns/naiveB_indexed_sequences.csv

(date; time igor -threads 10 -set_wd run_igor_20190913 -batch naiveB -species human -chain heavy_naive -align --V  -set_genomic --V prepared_for_IGoR/reference/genomicVs.fasta; date) 2>&1 | tee run_igor_20190913/igor_align_V.20190913.log

(date; time igor -threads 1 -set_wd run_igor_20190913 -batch naiveB -species human -chain heavy_naive -align --D  -set_genomic --D prepared_for_IGoR/reference/genomicDs.fasta; date) 2>&1 | tee run_igor_20190913/igor_align_D.20190913.log

(date; time igor -threads 1 -set_wd run_igor_20190913 -batch naiveB -species human -chain heavy_naive -align --J  -set_genomic --J prepared_for_IGoR/reference/genomicJs.fasta; date) 2>&1 | tee run_igor_20190913/igor_align_D.20190913.log

```


### Run IGoR-infer

#### Split sequences and alignments to about 200,000 sequences per file

My practical experience is to split the input files into small files when running IGoR-infer, due to memory limit. My naiveB_indexed_sequences.csv file contains 2984791 input sequences; therefore, I split it into 15 parts (i=0..14).

```bash
head run_igor_20190913/aligns/naiveB_indexed_sequences.csv
# idx 0.., has headline
tail run_igor_20190913/aligns/naiveB_indexed_sequences.csv
# idx ..2984790

mkdir run_igor_20190913/aligns/split_20190919
rm -f run_igor_20190913/aligns/split_20190919/*

date; for i in {0..9}; do
let shift=-i*200000
part=p0${i}
(date
echo ${part}
echo ${i} idx_seq
echo cat run_igor_20190913/aligns/naiveB_indexed_sequences.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_indexed_sequences.csv
time cat run_igor_20190913/aligns/naiveB_indexed_sequences.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_indexed_sequences.csv
echo ${i} V
echo cat run_igor_20190913/aligns/naiveB_V_alignments.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_V_alignments.csv
time cat run_igor_20190913/aligns/naiveB_V_alignments.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_V_alignments.csv
echo ${i} D
echo cat run_igor_20190913/aligns/naiveB_D_alignments.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_D_alignments.csv
time cat run_igor_20190913/aligns/naiveB_D_alignments.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_D_alignments.csv
echo ${i} J
echo cat run_igor_20190913/aligns/naiveB_J_alignments.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_J_alignments.csv
time cat run_igor_20190913/aligns/naiveB_J_alignments.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_J_alignments.csv
date) 2>&1 | cat >run_igor_20190913/aligns/split_20190919/split_p0${i}.20190919.log &
done
date

date; for i in {10..14}; do
let shift=-i*200000
part=p${i}
(date
echo ${part}
echo ${i} idx_seq
echo cat run_igor_20190913/aligns/naiveB_indexed_sequences.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_indexed_sequences.csv
time cat run_igor_20190913/aligns/naiveB_indexed_sequences.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_indexed_sequences.csv
echo ${i} V
echo cat run_igor_20190913/aligns/naiveB_V_alignments.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_V_alignments.csv
time cat run_igor_20190913/aligns/naiveB_V_alignments.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_V_alignments.csv
echo ${i} D
echo cat run_igor_20190913/aligns/naiveB_D_alignments.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_D_alignments.csv
time cat run_igor_20190913/aligns/naiveB_D_alignments.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_D_alignments.csv
echo ${i} J
echo cat run_igor_20190913/aligns/naiveB_J_alignments.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_J_alignments.csv
time cat run_igor_20190913/aligns/naiveB_J_alignments.csv | python3 scripts/yyx_rearrange_alignment.20181120.py $shift 0 199999 True >run_igor_20190913/aligns/split_20190919/naiveB${part}i0_J_alignments.csv
date) 2>&1 | cat >run_igor_20190913/aligns/split_20190919/split_p${i}.20190919.log &
done
date
```

#### Prepare IGoR initial model file (model_parms.txt)

Select a template model file with similar V(D)J structure with that you are studying.
Then modified it according to the genomicV/D/Js.fasta

For example, I modified
* human BCR heavy (igor/share/igor/models/human/bcr_heavy/models/model_parms.txt) for mouse IgH V-D-J recombination

```bash
cat igor/share/igor/models/human/bcr_heavy/models/model_parms.txt | python3 scripts/yyx_simple_sort_model_parms.20181116.py >tmp.model_parms.txt

cat tmp.model_parms.txt | perl -e '
@V_block = ();
open(IN, "prepared_for_IGoR/reference/genomicVs.fasta") or die;
while(<IN>){
  chomp;
  s/^>//g;
  $title=$_;
  $seq = <IN>;
  chomp($seq);
  push(@V_block, "%".$title.";".$seq.";".scalar(@V_block));
}
close(IN);
$VgeneName = "GeneChoice_V_gene_Undefined_side_prio8_size".scalar(@V_block);

@D_block = ();
open(IN, "prepared_for_IGoR/reference/genomicDs.fasta") or die;
while(<IN>){
  chomp;
  s/^>//g;
  $title=$_;
  $seq = <IN>;
  chomp($seq);
  push(@D_block, "%".$title.";".$seq.";".scalar(@D_block));
}
close(IN);
$DgeneName = "GeneChoice_D_gene_Undefined_side_prio6_size".scalar(@D_block);

@J_block = ();
open(IN, "prepared_for_IGoR/reference/genomicJs.fasta") or die;
while(<IN>){
  chomp;
  s/^>//g;
  $title=$_;
  $seq = <IN>;
  chomp($seq);
  push(@J_block, "%".$title.";".$seq.";".scalar(@J_block));
}
close(IN);
$JgeneName = "GeneChoice_J_gene_Undefined_side_prio7_size".scalar(@J_block);
#print STDERR join("\n", @J_block)."\n";

$skip = "";
while(<STDIN>){
  if($skip ne "" && /^$skip/){  next;  }
  if(/^#/){
    $skip = "";
    if(/^#GeneChoice;V_gene/){
      print;
      print join("\n", @V_block)."\n";
      $skip = "%";
    }elsif(/^#GeneChoice;D_gene/){
      print;
      print join("\n", @D_block)."\n";
      $skip = "%";
    }elsif(/^#GeneChoice;J_gene/){
      print;
      print join("\n", @J_block)."\n";
      $skip = "%";
    }else{
      print;
    }
  }else{
    s/GeneChoice_V_gene_Undefined_side_prio8_size97/$VgeneName/g;
    s/GeneChoice_J_gene_Undefined_side_prio7_size7/$JgeneName/g;
    s/GeneChoice_D_gene_Undefined_side_prio6_size35/$DgeneName/g;
    print;
  }
}
' >prepared_for_IGoR/mouse_bcr_heavy.model_parms.txt
```

* human TCR alpha (igor/share/igor/models/human/tcr_alpha/models/model_parms.txt) for mouse IgL V-J recombination

```bash
cat igor/share/igor/models/human/tcr_alpha/models/model_parms.txt | python3 scripts/yyx_simple_sort_model_parms.20181116.py >tmp.model_parms.txt

cat tmp.model_parms.txt | perl -e '
@V_block = ();
open(IN, "prepared_for_IGoR/reference/genomicVs.fasta") or die;
while(<IN>){
  chomp;
  s/^>//g;
  $title=$_;
  $seq = <IN>;
  chomp($seq);
  push(@V_block, "%".$title.";".$seq.";".scalar(@V_block));
}
close(IN);
$VgeneName = "GeneChoice_V_gene_Undefined_side_prio7_size".scalar(@V_block);

@J_block = ();
open(IN, "prepared_for_IGoR/reference/genomicJs.fasta") or die;
while(<IN>){
  chomp;
  s/^>//g;
  $title=$_;
  $seq = <IN>;
  chomp($seq);
  push(@J_block, "%".$title.";".$seq.";".scalar(@J_block));
}
close(IN);
$JgeneName = "GeneChoice_J_gene_Undefined_side_prio6_size".scalar(@J_block);
#print STDERR join("\n", @J_block)."\n";

$skip = "";
while(<STDIN>){
  if($skip ne "" && /^$skip/){  next;  }
  if(/^#/){
    $skip = "";
    if(/^#GeneChoice;V_gene/){
      print;
      print join("\n", @V_block)."\n";
      $skip = "%";
    }elsif(/^#GeneChoice;J_gene/){
      print;
      print join("\n", @J_block)."\n";
      $skip = "%";
    }else{
      print;
    }
  }else{
    s/GeneChoice_V_gene_Undefined_side_prio7_size103/$VgeneName/g;
    s/GeneChoice_J_gene_Undefined_side_prio6_size68/$JgeneName/g;
    print;
  }
}
' >prepared_for_IGoR/mouse_light.model_parms.txt
```

#### Run IGoR-infer, manually rotating parts for several iterations

I run IGoR-infer, manually rotating 15 splited parts (p00..14) for 6 iterations (i0..5).

```bash
ls run_igor_20190913/aligns/naiveBp*
rm -f run_igor_20190913/aligns/naiveBp*

(date; time for iter in {0..5}; do
 echo Iteration $iter start
 date; time for part in {0..14}; do
  echo iter=$iter part=$part
  let prev=part-1
  if [ $part -eq 0 ]; then
   let prev=14
  fi
  part=`echo ${part}|awk '{printf("%02d\n",$0)}'`
  prev=`echo ${prev}|awk '{printf("%02d\n",$0)}'`
  name=naiveBp${part}i${iter}
  preN=naiveBp${prev}i${iter}
  if [ $part -eq 0 ]; then
   let prevIter=iter-1
   preN=naiveBp${prev}i${prevIter}
  fi
  modelDir=run_igor_20190913/${preN}_inference/
  echo link: ${name}
  ln -s run_igor_20190913/aligns/split_20190919/naiveBp${part}i0_indexed_sequences.csv  run_igor_20190913/aligns/${name}_indexed_sequences.csv
  for VDJ in V D J; do
   ln -s run_igor_20190913/aligns/split_20190919/naiveBp${part}i0_${VDJ}_alignments.csv  run_igor_20190913/aligns/${name}_${VDJ}_alignments.csv
  done
  
  echo output: run_igor_20190913/${name}_inference
  mkdir run_igor_20190913/${name}_inference
  if [ $iter -eq 0 ] && [ $part -eq 0 ]; then
   echo initial model: prepared_for_IGoR/mouse_bcr_heavy.model_parms.txt
   (date; time igor -set_wd run_igor_20190913 -batch ${name} -species human -chain heavy_naive -infer --N_iter 1  -set_custom_model prepared_for_IGoR/mouse_bcr_heavy.model_parms.txt  -set_genomic --V prepared_for_IGoR/reference/genomicVs.fasta --D prepared_for_IGoR/reference/genomicDs.fasta --J prepared_for_IGoR/reference/genomicJs.fasta; date) 2>&1 | tee run_igor_20190913/${name}_inference/igor_infer.20190919.log
  else
   echo initial model: $modelDir/final_
   (date; time igor -set_wd run_igor_20190913 -batch ${name} -species human -chain heavy_naive -infer --N_iter 1  -set_custom_model $modelDir/final_parms.txt $modelDir/final_marginals.txt  -set_genomic --V prepared_for_IGoR/reference/genomicVs.fasta --D prepared_for_IGoR/reference/genomicDs.fasta --J prepared_for_IGoR/reference/genomicJs.fasta; date) 2>&1 | tee run_igor_20190913/${name}_inference/igor_infer.20190919.log
  fi
  
  echo remove: ${name}
  rm -f run_igor_20190913/aligns/${name}_indexed_sequences.csv
  for VDJ in V D J; do
   rm -f run_igor_20190913/aligns/${name}_${VDJ}_alignments.csv
  done
  rmdir run_igor_20190913/${name}_output
  chmod -w run_igor_20190913/${name}_inference
  chmod -w run_igor_20190913/${name}_inference/*
 done; date
 echo Iteration $iter end
done; date) 2>&1 | tee -a run_igor_20190913/igor_infer_iter.20190919.log
```

#### Examine convergence of log-likelihood

```bash
perl -e '
print join("\t", "iter_", 0..5)."\n";
for($p=0; $p<15; $p++){
$sp = $p;
if($sp < 10){ $sp = "0" . $p; }
print "part_$sp";
for($i=0; $i<=5; $i++){
$input_filename = "run_igor_20190913/naiveBp" . $sp . "i". $i . "_inference/likelihoods.out";
open(IN, $input_filename) or die "Error: cannot open $input_filename for input\n";
$headline=<IN>;
$line=<IN>;
chomp($line);
@F=split(/;/, $line);
$out=$F[1]." (".$F[2].")";
print "\t".$out;
close(IN);
}
print "\n";
}
'
```

Since it got convergent fast with so many training input sequences, I just choose p12i5 p04i5 p09i4 p00i4 as final models (and calculate average Pgen generated by them).


### IGoR-evaluate

#### Calculate Pgen of the (training) input sequences based on final four models

```bash
(date; time for part in {0..14}; do
 echo part=$part start
 date; time for modelName in p12i5 p04i5 p09i4 p00i4; do
  modelDir=run_igor_20190913/naiveB${modelName}_inference/
  echo model=$modelName
  part=`echo ${part}|awk '{printf("%02d\n",$0)}'`
  name=naiveBp${part}_${modelName}model
  out=run_igor_20190913/${name}_output/Pgen_counts.csv
  if [ -f $out ] && [ `wc -c $out | awk '{print $1;}'` -gt 10000 ]; then
  	echo output $out exists, so I will skip this sample
  	continue
  fi
  echo link: ${name}
  ln -s /Users/yyx/Documents/work/huan/IGoR/run_igor_20190913/aligns/split_20190919/naiveBp${part}i0_indexed_sequences.csv  /Users/yyx/Documents/work/huan/IGoR/run_igor_20190913/aligns/${name}_indexed_sequences.csv
  for VDJ in V D J; do
   ln -s /Users/yyx/Documents/work/huan/IGoR/run_igor_20190913/aligns/split_20190919/naiveBp${part}i0_${VDJ}_alignments.csv  /Users/yyx/Documents/work/huan/IGoR/run_igor_20190913/aligns/${name}_${VDJ}_alignments.csv
  done
  
  echo use model: $modelDir/final_
  echo output: run_igor_20190913/${name}_evaluate run_igor_20190913/${name}_output
  mkdir run_igor_20190913/${name}_evaluate
  (date; time igor -set_wd run_igor_20190913 -batch ${name} -species human -chain heavy_naive -evaluate -set_custom_model $modelDir/final_parms.txt $modelDir/final_marginals.txt  -output --scenarios 10 --Pgen; date) 2>&1 | tee run_igor_20190913/${name}_evaluate/igor_evaluate.20190923.log
  
  echo remove: ${name}
  rm -f /Users/yyx/Documents/work/huan/IGoR/run_igor_20190913/aligns/${name}_indexed_sequences.csv
  for VDJ in V D J; do
   rm -f /Users/yyx/Documents/work/huan/IGoR/run_igor_20190913/aligns/${name}_${VDJ}_alignments.csv
  done
  chmod -w /Users/yyx/Documents/work/huan/IGoR/run_igor_20190913/${name}_evaluate
  chmod -w /Users/yyx/Documents/work/huan/IGoR/run_igor_20190913/${name}_evaluate/*
  chmod -w /Users/yyx/Documents/work/huan/IGoR/run_igor_20190913/${name}_output
  chmod -w /Users/yyx/Documents/work/huan/IGoR/run_igor_20190913/${name}_output/*
 done; date
 echo part=$part end
done; date) 2>&1 | tee -a run_igor/igor_evaluate_iter.20190923.log
```

#### Merge Pgen results into one tabular file

```bash
mkdir run_igor_20190913/Pgen_merged/
rm -f run_igor_20190913/Pgen_merged/*

for part in {0..14}; do
 if [[ $part -lt 10 ]]; then
  part="0"$part
 fi
 name=naiveBp$part
 echo part=$part start
 time python3 scripts/yyx_merge_Pgen_files.20181124.py \
 p12i5m run_igor_20190913/${name}_p12i5model_output/Pgen_counts.csv \
 p04i5m run_igor_20190913/${name}_p04i5model_output/Pgen_counts.csv \
 p09i4m run_igor_20190913/${name}_p09i4model_output/Pgen_counts.csv \
 p00i4m run_igor_20190913/${name}_p00i4model_output/Pgen_counts.csv \
 >run_igor_20190913/Pgen_merged/p${part}_4models.20190923.tsv
done
```

#### Plot Pgen distribution

Do statistical test and plot figures of violin plot and histogram in R code

```R
options(stringsAsFactors=FALSE)

setwd("run_igor_20190913/Pgen_merged")


`%.%` = function(x,y) paste0(x,y)


input = read.delim("recurrentCDR3_VDJ_4models.20190925.tsv")
input = as.matrix(input[, -1])

mean_log10_Pgen_vec = rowMeans(log10(input))
names(mean_log10_Pgen_vec) = c("GF2", "GF1", "SPF2", "SPF4", "", "NP1", "NP2", "", "", "SPF5", "SPF6", "SPF3", "", "", "SPF1")
mean_log10_Pgen_vec = mean_log10_Pgen_vec[c("GF" %.% 1:2, "SPF" %.% 1:6, "NP" %.% 1:2)]
10^mean_log10_Pgen_vec
'
         GF1          GF2         SPF1         SPF2         SPF3         SPF4         SPF5         SPF6          NP1          NP2 
2.101051e-08 1.916569e-09 4.347344e-09 7.698433e-10 8.182526e-10 9.692503e-09 9.255171e-10 2.973112e-09 4.950832e-08 3.091450e-07 
'

mean_Pgen_vec = rowMeans(input)
names(mean_Pgen_vec) = c("GF2", "GF1", "SPF2", "SPF4", "", "NP1", "NP2", "", "", "SPF5", "SPF6", "SPF3", "", "", "SPF1")
mean_Pgen_vec = mean_Pgen_vec[c("GF" %.% 1:2, "SPF" %.% 1:6, "NP" %.% 1:2)]
mean_Pgen_vec
'
         GF1          GF2         SPF1         SPF2         SPF3         SPF4         SPF5         SPF6          NP1          NP2 
2.162217e-08 1.956763e-09 4.486325e-09 7.788800e-10 8.204955e-10 9.811540e-09 9.344533e-10 3.017750e-09 4.982090e-08 3.102360e-07 
'
range(log10(mean_Pgen_vec), na.rm=TRUE)   # -9.108529 -6.508308
(9.6-6.5)/30   # 0.1033333


all_mean_Pgen_vec = numeric(0)
system.time({
for(part in 0:14){
	print(part)
	if(part < 10){
		part = "0" %.% part
	}
	now_input_filename = "p" %.% part %.% "_4models.20190923.tsv"
	now_input = read.delim(now_input_filename)
	now_mean_Pgen_vec = rowMeans(now_input[, -1])
	all_mean_Pgen_vec = c(all_mean_Pgen_vec, now_mean_Pgen_vec)
}
})   # 27s
length(all_mean_Pgen_vec)   # 2984791
sum(!is.na(all_mean_Pgen_vec))   # 2914266

range(log10(all_mean_Pgen_vec), na.rm=TRUE)   # -52.538420  -6.073687


wilcox.test(mean_Pgen_vec, all_mean_Pgen_vec)
#  p-value = 1.955e-06



## unique on CDR3_SEQ
### ref: work_flow.check_and_run_IGoR_align.20181111.sh

db_input = read.delim("../../naiveB.all.db")
dim(db_input)   # 2984791      67

CDR3_dup_bool_vec = duplicated(db_input$CDR3_SEQ)
table(CDR3_dup_bool_vec)
'
CDR3_dup_bool_vec
  FALSE    TRUE 
 578654 2406137 
'
table(CDR3_dup_bool_vec) / sum(table(CDR3_dup_bool_vec))
'
CDR3_dup_bool_vec
    FALSE      TRUE 
0.1938675 0.8061325 
'


tmpDF = data.frame(log10_Pgen = log10(all_mean_Pgen_vec[!CDR3_dup_bool_vec]), group = "others")
dim(tmpDF)   # 578654      2

for(k in 1:length(mean_Pgen_vec)){
	if(grepl("^NP", names(mean_Pgen_vec)[k])){
		next   # skip NP positive controls
	}
	idx = (which.min(abs(tmpDF$log10_Pgen-log10(mean_Pgen_vec[k]))))[1]
	tmpDF$group[idx] = "recurrent CDR3s"
}
tmpDF$group = factor(tmpDF$group, levels=c("recurrent CDR3s", "others"))

table(tmpDF$group)
'
recurrent CDR3s          others 
              8          578646 
'

with(tmpDF, wilcox.test(log10_Pgen ~ group))
#  p-value = 1.553e-05
pVal = with(tmpDF, wilcox.test(log10_Pgen ~ group))$p.value



### ref: http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
### ref: http://www.sthda.com/english/wiki/ggplot2-dot-plot-quick-start-guide-r-software-and-data-visualization

#install.packages("data.table")
#install.packages("Hmisc")
#library(Hmisc)   # for mean_sdl()
library(ggplot2)

#ggplot(tmpDF, aes(x=group, y=log10_Pgen)) + geom_violin(mapping=aes(fill=group), trim=FALSE, color=NA) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="black"))
#mean_sdl(log10(mean_Pgen_vec))
'
          y     ymin      ymax
1 -8.228504 -9.94939 -6.507618
'

pdf("log10_mean_Pgen_distr.violin.uniq_CDR3_SEQ.20191017.pdf", width=7, height=7)

print(ggplot(subset(tmpDF, "recurrent CDR3s"=="recurrent CDR3s"), aes(x=group, y=log10_Pgen)) + geom_violin(mapping=aes(fill=group), trim=FALSE, color=NA) + geom_boxplot(width=0.04, outlier.shape=NA, fill="gray85") + geom_dotplot(data=subset(tmpDF, group=="recurrent CDR3s"), mapping=aes(x=0.95),  binaxis='y', binwidth=0.05, stackdir='center', dotsize=3, stackratio=2.5, fill=hsv(0,1,0.7), color="gray30") + coord_cartesian(ylim=c(-25, -5)) + scale_fill_manual(breaks=levels(tmpDF$group), values=hsv(c(0,0.66), 0.5, 1)) + theme_bw() + theme(legend.position="none") + ylab("log10(Pgen)") + ggtitle("Mann-Whitney U test p-value = " %.% sprintf("%.2e", pVal)))
#print(ggplot(tmpDF, aes(x=group, y=log10_Pgen)) + geom_violin(mapping=aes(fill=group), trim=FALSE, color=NA) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") + geom_dotplot(data=subset(tmpDF, group=="recurrent CDR3s"), binaxis='y', binwidth=0.1, stackdir='center', dotsize=2, stackratio=2, fill="gray") + coord_cartesian(ylim=c(-25, -5)) + scale_fill_manual(breaks=levels(tmpDF$group), values=hsv(c(0,0.66), 0.5, 1)) + theme_bw() + theme(legend.position="none") + ylab("log10(Pgen)"))
print(ggplot(subset(tmpDF, "recurrent CDR3s"=="recurrent CDR3s"), aes(x=group, y=log10_Pgen)) + geom_violin(mapping=aes(fill=group), trim=FALSE, color=NA) + geom_boxplot(width=0.04, outlier.shape=NA, fill="gray85") + geom_dotplot(data=subset(tmpDF, group=="recurrent CDR3s"), mapping=aes(x=0.95),  binaxis='y', binwidth=0.05, stackdir='center', dotsize=3, stackratio=2.5, fill=hsv(0,1,0.7), color="gray30") + scale_fill_manual(breaks=levels(tmpDF$group), values=hsv(c(0,0.66), 0.5, 1)) + theme_bw() + theme(legend.position="none") + ylab("log10(Pgen)") + ggtitle("Mann-Whitney U test p-value = " %.% sprintf("%.2e", pVal)) + coord_flip(ylim=c(-25, -5)))
#  Warning: Removed 70525 rows containing non-finite values (stat_ydensity or stat_boxplot)

dev.off()


GF_SPF_NP_col = c("red3", "black", "blue")

all_mean_Pgen_vec = all_mean_Pgen_vec[!CDR3_dup_bool_vec]

pdf("log10_mean_Pgen_distr.hist.uniq_CDR3_SEQ.20191017.pdf", width=7, height=7)
#par(mar=c(5,3,1,1)+.1, mgp=c(2,0.8,0))
par(mar=c(3,3,1,1)+.1, mgp=c(2,0.8,0))

hist_rlt = hist(log10(all_mean_Pgen_vec), breaks=seq(-53.5, -5.5, by=1), xlim=c(-25, -5), xlab="log10(Pgen)", ylab="Frequency", main="", axes=TRUE)
#title(xlab="log10(Pgen)", line=3.9)
#axis(2)
#axis(1, line=2)
#rug(log10(mean_Pgen_vec))
#mtext("median = " %.% format(median(log10(all_mean_Pgen_vec), na.rm=TRUE), nsmall=2, digits=4))
median_all_mean_Pgen = median(log10(all_mean_Pgen_vec), na.rm=TRUE)
y_max = max(hist_rlt$counts)
#abline(v = median_all_mean_Pgen, lty=3, col="gray")
segments(median_all_mean_Pgen, -y_max*0.035, median_all_mean_Pgen, y_max, col="gray25", lty=3, xpd=NA)
text(median_all_mean_Pgen, y_max, "median = " %.% sprintf("%.2f", median_all_mean_Pgen), col="gray25", adj=c(0.5, -0.3), cex=0.9)

#rug(log10(mean_Pgen_vec["GF" %.% 1:2]), col=GF_SPF_NP_col[1], lwd=1)
#rug(log10(mean_Pgen_vec["SPF" %.% 1:6]), col=GF_SPF_NP_col[2], lwd=1)
#rug(log10(mean_Pgen_vec["NP" %.% 1:2]), col=GF_SPF_NP_col[3], lwd=1)
for(k in 1:2){
	now_Pgen = log10(mean_Pgen_vec["GF" %.% k])
	segments(now_Pgen, -y_max*0.01, now_Pgen, -y_max*0.035, col=GF_SPF_NP_col[1], xpd=NA)
}
for(k in 1:6){
	now_Pgen = log10(mean_Pgen_vec["SPF" %.% k])
#	segments(now_Pgen, -y_max*0.04, now_Pgen, -y_max*0.07, col=GF_SPF_NP_col[2], xpd=NA)
	segments(now_Pgen, -y_max*0.01, now_Pgen, -y_max*0.035, col=GF_SPF_NP_col[2], xpd=NA)
}
for(k in 1:2){
	now_Pgen = log10(mean_Pgen_vec["NP" %.% k])
#	segments(now_Pgen, -y_max*0.07, now_Pgen, -y_max*0.10, col=GF_SPF_NP_col[3], xpd=NA)
	segments(now_Pgen, -y_max*0.01, now_Pgen, -y_max*0.035, col=GF_SPF_NP_col[3], xpd=NA)
}
#text(-25, -y_max*c(0.025, 0.055, 0.085), c("GF", "SPF", "NP"), adj=c(0, 0.5), xpd=NA, cex=0.9)

#text(-25, y_max, paste0(collapse="\n", names(mean_Pgen_vec) %.% " = " %.% format(mean_Pgen_vec, nsmall=2, digits=4)), adj=c(0,1), cex=0.9)
for(k in 1:2){
	now_name = "GF" %.% k
	now_Pgen = (mean_Pgen_vec[now_name])
	text(-25, y_max*(0.9-0.04*(k-1)), now_name %.% "   = " %.% sprintf("%.1e", now_Pgen) %.% " = 10^" %.% sprintf("%.2f", log10(now_Pgen)), adj=c(0,1), cex=0.9, col=GF_SPF_NP_col[1])
}
for(k in 1:6){
	now_name = "SPF" %.% k
	now_Pgen = (mean_Pgen_vec[now_name])
	text(-25, y_max*(0.9-0.02-0.04*(k-1+2)), now_name %.% " = " %.% sprintf("%.1e", now_Pgen) %.% " = 10^" %.% sprintf("%.2f", log10(now_Pgen)), adj=c(0,1), cex=0.9, col=GF_SPF_NP_col[2])
}
for(k in 1:2){
	now_name = "NP" %.% k
	now_Pgen = (mean_Pgen_vec[now_name])
	text(-25, y_max*(0.9-0.02*2-0.04*(k-1+8)), now_name %.% "   = " %.% sprintf("%.1e", now_Pgen) %.% " = 10^" %.% sprintf("%.2f", log10(now_Pgen)), adj=c(0,1), cex=0.9, col=GF_SPF_NP_col[3])
}

dev.off()
```

#### (Optional) Annotate IGoR bestScenarios in a human-readable way

```
time python3 scripts/yyx_annotate_bestScenarios.20181206.py prepared_for_IGoR/reference/ input_xxx_indexed_sequences.tsv output_xxx_evaluate/initial_model.txt output_xxx_evaluate/initial_marginals.txt output_xxx_output/best_scenarios_counts.csv output_xxx_output/Pgen_counts.csv False | perl -pe 's/;/\t/g' >output_xxx_output/bestScenarios_Pgen_PVDJ.tsv
```



## Back-translate

We back-translated CDR3 amino acid sequence to all possible nucleotide sequences, and then inserted them back to VDJ sequence, and finally examined their Pgen

### Back-translate

```bash
mkdir 20181205_backtranslate

time python3 scripts/yyx_translate.20181205.py F T 20181205_From_codon_usage/Mus_musculus.standard.txt <(tail -n+2 recurrent7CDR3_indexed_sequences.csv | perl -ne 'chomp; @F=split/;/; print $F[1]."\n";') >20181205_backtranslate/recurrent7CDR3_aa_seq.20181205.tsv

time for i in {1..7}; do
echo $i
time python3 scripts/yyx_translate.20181205.py B T 20181205_From_codon_usage/Mus_musculus.standard.txt <(tail -n+$i 20181205_backtranslate/recurrent7CDR3_aa_seq.20181205.tsv | head -n1 | cut -f1) >20181205_backtranslate/backtranslated_nt_seq.$i.20181205.tsv
wc -l 20181205_backtranslate/backtranslated_nt_seq.$i.20181205.tsv
done
```

Similarly, I also manually splited input files into several small files if the number of sequences was much larger than 200,000.

```bash
mkdir run_igor_4/aligns/back_20181205

for i in {1..4} 6; do
echo $i
echo 'seq_index;sequence' >run_igor_4/aligns/back_20181205/back_${i}_indexed_sequences.csv
cat 20181205_backtranslate/backtranslated_nt_seq.$i.20181205.tsv | perl -ne 'BEGIN{$NR=0;} @F=split/\t/; print join(";", $NR++, $F[0])."\n"' >>run_igor_4/aligns/back_20181205/back_${i}_indexed_sequences.csv
done

for j in {1..8}; do
let s=(j-1)*200000+1
echo $j, $s
echo 'seq_index;sequence' >run_igor_4/aligns/back_20181205/back_5${j}_indexed_sequences.csv
tail -n+$s 20181205_backtranslate/backtranslated_nt_seq.5.20181205.tsv | head -n200000 | perl -ne 'BEGIN{$NR=0;} @F=split/\t/; print join(";", $NR++, $F[0])."\n"' >>run_igor_4/aligns/back_20181205/back_5${j}_indexed_sequences.csv
done
```

Check the position of original CDR3 sequence in back-translated files.

```bash
grep -f <(tail -n+2 recurrent7CDR3_indexed_sequences.csv | cut -d';' -f2) run_igor_4/aligns/back_20181205/*_indexed_sequences.csv
```

Manually write into a tabular file.

```bash
vi run_igor_4/aligns/back_20181205/original_CDR3_idx_in_back.txt
'
PP1	1	110307
PP2	2	2161
PP3	3	2075
PP4	4	753
PP5	56	56043
NP1	6	50379
NP2	72	82299
'
```

### Insert back to its original VDJ sequence


### prepare indexed_sequences

```bash
mkdir -p run_igor_PP1_5_NP1_2_20191013/aligns

cat run_igor_4/aligns/back_20181205/original_CDR3_idx_in_back.txt | while read title part idx; do
 echo $title $part $idx
 if [[ -z $part ]]; then
  continue
 fi
 partPre=${part:0:1}
 recIdx=$[partPre-1]
 ls run_igor_4/aligns/back_20181205/back_$partPre*_indexed_sequences.csv | while read inF; do
  f=${inF##*/}
  echo $f
  cat $inF | perl -e '
$refFullVDJ = "";
open(IN, "recurrent7CDR3_indexed_sequences.csv") or die 1;
while(<IN>){
 s/[\r\n]+$//;
 @F = split/;/;
 if($F[0] eq "'$recIdx'"){
  $refFullVDJ = $F[1];
  last;
 }
}
close(IN);

$refCDR3 = "";
open(IN, "run_igor_4/aligns/back_20181205/back_'$part'_indexed_sequences.csv") or die 2;
while(<IN>){
 s/[\r\n]+$//;
 @F = split/;/;
 if($F[0] eq "'$idx'"){
  $refCDR3 = $F[1];
  last;
 }
}
close(IN);

while(<STDIN>){
 s/[\r\n]+$//;
 @F = split/;/;
 $tmp = $refFullVDJ;
 $tmp =~ s/$refCDR3/$F[1]/;
 $F[1] = $tmp;
 print join(";", @F)."\n";
}
' >run_igor_PP1_5_NP1_2_20191013/aligns/$f
 done
done
```

Then, run IGoR-align and IGoR-evaluate similar to previous parts.

And I rearranged IGoR Pgen into tabular file as before.

```bash
mkdir run_igor_PP1_5_NP1_2_20191013/Pgen_merged.20191016

time for part in {1..4} 6 {71..72} {51..58}; do
 part=back_$part
 echo part=$part start
 time python3 scripts/yyx_merge_Pgen_files.20181124.py \
 p12i5m run_igor_PP1_5_NP1_2_20191013/${part}_p12i5model_output/Pgen_counts.csv \
 p04i5m run_igor_PP1_5_NP1_2_20191013/${part}_p04i5model_output/Pgen_counts.csv \
 p09i4m run_igor_PP1_5_NP1_2_20191013/${part}_p09i4model_output/Pgen_counts.csv \
 p00i4m run_igor_PP1_5_NP1_2_20191013/${part}_p00i4model_output/Pgen_counts.csv \
 >run_igor_PP1_5_NP1_2_20191013/Pgen_merged.20191016/${part}_4models.tsv
done
```

Finally, plot histogram to show the position of Pgen of original CDR3 among the Pgen distribution of all back-translated sequences

```bash
mkdir parsed_back_Pgen_merged_20191016

(date; time cat run_igor_4/aligns/back_20181205/original_CDR3_idx_in_back.txt | while read title part idx; do
echo $title $part $idx
seqIdx=${part:0:1}
seqIdx=$[seqIdx-1]
mean_log10_Pgen=`grep "^$seqIdx\t" run_igor_4/Pgen_merged/recurrent7CDR3_VDJ_4models.20190925.tsv | perl -ne 'chomp; @F=split/\t/; $sum=0; for($i=1;$i<=4;$i++){ $sum+=log($F[$i])/log(10); } print ($sum/4)'`
echo -e "$seqIdx\t$mean_log10_Pgen"
time Rscript scripts/yyx_parse_back_Pgen_merged.scenarios_without_mismatches.20191008.r $mean_log10_Pgen $title parsed_back_Pgen_merged_20191016/back_${title}_4models.20191016 run_igor_PP1_5_NP1_2_20191013/Pgen_merged.20191016/back_${part:0:1}*_4models.tsv
done; date) 2>&1 | tee parsed_back_Pgen_merged_20191016/recurrent7CDR3.yyx_parse_back_Pgen_merged.20191016.log
```

