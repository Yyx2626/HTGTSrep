# stratify\_and\_Bayesian

Some scripts for stratification and hierarchical Bayesian model for somatic hypermutation profile

Author: Adam Yongxin Ye @ Boston Children's Hospital / Harvard Medical School


## Pipeline

### Step 0. Prepare .stat.txt and .nuc.txt files from HTGTSrep pipeline

prepare background intrinsic non-productive .nuc.txt for each V gene in folder mut\_profile/nucl\_text/, with filenames \*.${Vgene}-01.NP.nuc.txt

and prepare foreground clonotype .nuc.txt (and one .stat.txt) files in folder clonotype_files/, with filenames ${cln}.\*.nuc.txt (and ${cln}.\*.stat.txt)

* .stat.txt file
  * each row for one nucleotide position in the V gene; the 1st-2nd columns should be Pos, Base

* .nuc.txt file
  * base [ACGTN-.] on each position, each column for one nucleotide position in the V gene, except the 1st column is read ID


### Step 1. Stratify and merge nuc.txt files (output to stat.txt)

traverse each clnV (clonotype\_Vgene; for example, if clnV=clone12323\_IGHV2, then cln=clone12323 and Vgene=IGHV2)

```bash
for mutMaxProp in 0.025; do
basedir=stratify_mutProp_${mutMaxProp/./_}.20190906
echo $basedir

mkdir $basedir
mkdir $basedir/foreground_clone

for cln in clone1711 clone2238 clone3227 clone3621; do
 rm -f $basedir/foreground_clone/$cln*
done

for clnV in clone12323_IGHV2-2 clone13065_IGHV2-2  clone110_IGHV2-9 clone1635_IGHV2-9  clone549_IGHV1-72 clone619_IGHV1-72  clone8844_IGHV9-3  clone776_IGHV1-12  clone1580_IGHV5-17  clone748_IGHV1-11  clone1276_IGHV9-4  clone2644_IGHV6-6  clone2887_IGHV11-2; do
cln=${clnV%_*}
Vgene=${clnV#clone*_}
out=$basedir/foreground_clone/${clnV}_01
statfile=`find clonotype_files | grep $clnV | grep -v __MACOSX | grep stat.txt$ | head -n1`
in_nuc_dir=clonotype_files
(date
echo $cln  $clnV  $Vgene
echo Rscript scripts/yyx_SHMPlot2_Strata_ErrBar.20190731.r -fo  $out.SHMPlot_before_stratify  "$statfile"  $in_nuc_dir/${cln}.*.nuc.txt
time Rscript scripts/yyx_SHMPlot2_Strata_ErrBar.20190731.r -fo  $out.SHMPlot_before_stratify  "$statfile"  $in_nuc_dir/${cln}.*.nuc.txt
echo Rscript scripts/yyx_SHMPlot2_Strata_ErrBar.20190731.r -fo -mutMaxProp=$mutMaxProp  $out.SHMPlot_after_stratify  "$statfile"  $in_nuc_dir/${cln}.*.nuc.txt
time Rscript scripts/yyx_SHMPlot2_Strata_ErrBar.20190731.r -fo -mutMaxProp=$mutMaxProp  $out.SHMPlot_after_stratify  "$statfile"  $in_nuc_dir/${cln}.*.nuc.txt
echo Rscript scripts/yyx_stratify_from_nuc_to_merged_stat.20190906.r -fo -mutMaxProp=$mutMaxProp  $out.all  "$statfile"  $in_nuc_dir/${cln}.*.nuc.txt
time Rscript scripts/yyx_stratify_from_nuc_to_merged_stat.20190906.r -fo -mutMaxProp=$mutMaxProp  $out.all  "$statfile"  $in_nuc_dir/${cln}.*.nuc.txt
date) |& tee $out.yyx_stratify_from_nuc_to_merged_stat.20190911.log
done


mkdir $basedir/background_intrinsic
 
for clnV in clone12323_IGHV2-2 clone13065_IGHV2-2  clone110_IGHV2-9 clone1635_IGHV2-9  clone549_IGHV1-72 clone619_IGHV1-72  clone8844_IGHV9-3  clone776_IGHV1-12  clone1580_IGHV5-17  clone748_IGHV1-11  clone1276_IGHV9-4  clone2644_IGHV6-6  clone2887_IGHV11-2; do
Vgene=${clnV#clone*_}
out=$basedir/background_intrinsic/intrinsic_${Vgene}_01
statfile=`find clonotype_files | grep $clnV | grep -v __MACOSX | grep stat.txt$ | head -n1`
if [[ ! -f $out.all.stat.txt ]]; then
(date
echo Rscript scripts/yyx_SHMPlot2_Strata_ErrBar.20190731.r -fo  $out.SHMPlot_before_stratify  "$statfile"  mut_profile/nucl_text/*.${Vgene}-01.NP.nuc.txt
time Rscript scripts/yyx_SHMPlot2_Strata_ErrBar.20190731.r -fo   $out.SHMPlot_before_stratify  "$statfile"  mut_profile/nucl_text/*.${Vgene}-01.NP.nuc.txt
echo Rscript scripts/yyx_SHMPlot2_Strata_ErrBar.20190731.r -fo -mutMaxProp=$mutMaxProp  $out.SHMPlot_after_stratify  "$statfile"  mut_profile/nucl_text/*.${Vgene}-01.NP.nuc.txt
time Rscript scripts/yyx_SHMPlot2_Strata_ErrBar.20190731.r -fo -mutMaxProp=$mutMaxProp  $out.SHMPlot_after_stratify  "$statfile"  mut_profile/nucl_text/*.${Vgene}-01.NP.nuc.txt
echo Rscript scripts/yyx_stratify_from_nuc_to_merged_stat.20190906.r -fo -mutMaxProp=$mutMaxProp  $out.all  "$statfile"  mut_profile/nucl_text/*.${Vgene}-01.NP.nuc.txt
time Rscript scripts/yyx_stratify_from_nuc_to_merged_stat.20190906.r -fo -mutMaxProp=$mutMaxProp  $out.all  "$statfile"  mut_profile/nucl_text/*.${Vgene}-01.NP.nuc.txt
date) |& tee $out.yyx_stratify_from_nuc_to_merged_stat.20190911.log
fi
done

done   # end for mutMaxProp
```

output to $basedir/foreground\_clone/${clnV}\_01.all.stat.txt and $basedir/background\_intrinsic/intrinsic\_${Vgene}\_01.all.stat.txt

then, symlink to $basedir/mcmc\_rlt.20190911/

```bash
for mutMaxProp in 0.025; do
basedir=stratify_mutProp_${mutMaxProp/./_}.20190906
echo $basedir

mkdir $basedir/mcmc_rlt.20190911
#rm -f $basedir/mcmc_rlt.20190911/*

(date
for clnV in clone12323_IGHV2-2 clone13065_IGHV2-2  clone110_IGHV2-9 clone1635_IGHV2-9  clone549_IGHV1-72 clone619_IGHV1-72  clone8844_IGHV9-3  clone776_IGHV1-12  clone1580_IGHV5-17  clone748_IGHV1-11  clone1276_IGHV9-4  clone2644_IGHV6-6  clone2887_IGHV11-2; do
 Vgene=${clnV#clone*_}
 echo $clnV $Vgene
 echo symlink $basedir/foreground_clone/${clnV}_01.all.stat.txt '->' $basedir/mcmc_rlt.20190911/
 ln -s `realpath $basedir/foreground_clone/${clnV}_01.all.stat.txt` $basedir/mcmc_rlt.20190911/
 if [[ -f $basedir/mcmc_rlt.20190911/intrinsic_${Vgene}_01.all.stat.txt ]]; then
  echo $basedir/mcmc_rlt.20190911/intrinsic_${Vgene}_01.all.stat.txt already exists, so I skip it
 else
  mcmc_intrinsic=`ls $basedir/mcmc_rlt.*/intrinsic_*.all.stat.txt | grep $Vgene | head -n1`
  if [[ ! -z $mcmc_intrinsic ]]; then
   echo $Vgene already has mcmc_rlt
   echo symlink $mcmc_intrinsic '->' $basedir/mcmc_rlt.20190911/
   ln -s `realpath $mcmc_intrinsic`  $basedir/mcmc_rlt.20190911/
   echo symlink ${mcmc_intrinsic%.all.stat.txt}.site_\*.mcmc_rlt.tsv '->' $basedir/mcmc_rlt.20190911/
   ln -s `realpath ${mcmc_intrinsic%.all.stat.txt}.site_*.mcmc_rlt.tsv`  $basedir/mcmc_rlt.20190911/
  else
   echo $Vgene does not mcmc_rlt
   echo symlink $basedir/background_intrinsic/intrinsic_${Vgene}_01.all.stat.txt '->' $basedir/mcmc_rlt.20190911/
   ln -s `realpath $basedir/background_intrinsic/intrinsic_${Vgene}_01.all.stat.txt`  $basedir/mcmc_rlt.20190911/
  fi
 fi
done
date) |& tee $basedir/mcmc_rlt.20190911/symlink.all_stat.20190911.log

done   # end for mutMaxProp

ls stratify_mutProp_*.20190906/mcmc_rlt.20190911/*.all.stat.txt
```

### Step 2. JAGS MCMC

Prerequisites: Install JAGS and R package {rjags}

```bash
time for mutMaxProp in 0.025; do
basedir=stratify_mutProp_${mutMaxProp/./_}.20190906
date
echo $basedir

for clnV in clone12323_IGHV2-2 clone13065_IGHV2-2  clone110_IGHV2-9 clone1635_IGHV2-9  clone549_IGHV1-72 clone619_IGHV1-72  clone8844_IGHV9-3  clone776_IGHV1-12  clone1580_IGHV5-17  clone748_IGHV1-11  clone1276_IGHV9-4  clone2644_IGHV6-6  clone2887_IGHV11-2; do
Vgene=${clnV#clone*_}
echo $clnV  $Vgene

if [[ -f $basedir/mcmc_rlt.20190911/intrinsic_${Vgene}_01.site_150.mcmc_rlt.tsv ]]; then
 echo $Vgene already has mcmc_rlt, so I skip Step2 for it
else
 for ee in {10..300..10}; do
  ss=$((ee-9))
  # TODO: manually adjust above two lines, now for 30 threads in parallel
  echo $ss:$ee
  (date
  echo Rscript hierarchical_Bayesian_model_jags_SHM_test.20190819/Yyx_one_group_JAGS.20190820.r $basedir/mcmc_rlt.20190911/intrinsic_${Vgene}_01 $basedir/mcmc_rlt.20190911/intrinsic_${Vgene}_01.all.stat.txt  $ss:$ee
  time Rscript hierarchical_Bayesian_model_jags_SHM_test.20190819/Yyx_one_group_JAGS.20190820.r $basedir/mcmc_rlt.20190911/intrinsic_${Vgene}_01 $basedir/mcmc_rlt.20190911/intrinsic_${Vgene}_01.all.stat.txt  $ss:$ee
  date) |& cat >$basedir/mcmc_rlt.20190911/intrinsic_${Vgene}_01.Yyx_one_group_JAGS.${ss}_$ee.20190911.log &
 done
 wait
fi

for ee in {10..300..10}; do
 ss=$((ee-9))
 # TODO: manually adjust above two lines, now for 30 threads in parallel
 echo $ss:$ee
 (date
 echo Rscript hierarchical_Bayesian_model_jags_SHM_test.20190819/Yyx_one_group_JAGS.20190820.r $basedir/mcmc_rlt.20190911/${clnV}_01 $basedir/mcmc_rlt.20190911/${clnV}_01.all.stat.txt  $ss:$ee
 time Rscript hierarchical_Bayesian_model_jags_SHM_test.20190819/Yyx_one_group_JAGS.20190820.r $basedir/mcmc_rlt.20190911/${clnV}_01 $basedir/mcmc_rlt.20190911/${clnV}_01.all.stat.txt  $ss:$ee
 date) |& cat >$basedir/mcmc_rlt.20190911/${clnV}_01.Yyx_one_group_JAGS.${ss}_$ee.20190911.log &
done
wait

done   # end for clnV

done   # end for mutMaxProp
```

output to $basedir/mcmc\_rlt.20190911/intrinsic\_${Vgene}\_01.site\_\*.mcmc\_rlt.tsv and $basedir/mcmc\_rlt.20190911/${clnV}\_01.site\_\*.mcmc\_rlt.tsv

### Step 3. Posterior comparison

Prerequisites: Install R package {HDInterval}

```bash
time for mutMaxProp in 0.025; do
basedir=stratify_mutProp_${mutMaxProp/./_}.20190906
date
echo $basedir

mkdir $basedir/two_group_compare.20190911

for clnV in clone12323_IGHV2-2 clone13065_IGHV2-2  clone110_IGHV2-9 clone1635_IGHV2-9  clone549_IGHV1-72 clone619_IGHV1-72  clone8844_IGHV9-3  clone776_IGHV1-12  clone1580_IGHV5-17  clone748_IGHV1-11  clone1276_IGHV9-4  clone2644_IGHV6-6  clone2887_IGHV11-2; do
Vgene=${clnV#clone*_}
echo $clnV  $Vgene

for ee in {10..300..10}; do
 ss=$((ee-9))
 # TODO: manually adjust above two lines, now for 30 threads in parallel
 echo $ss:$ee
 (date
 echo Rscript hierarchical_Bayesian_model_jags_SHM_test.20190819/Yyx_two_group_compare.20190903.r $basedir/two_group_compare.20190911/${clnV}_01 $basedir/mcmc_rlt.20190911/${clnV}_01 $basedir/mcmc_rlt.20190911/intrinsic_${Vgene}_01  $ss:$ee
 time Rscript hierarchical_Bayesian_model_jags_SHM_test.20190819/Yyx_two_group_compare.20190903.r $basedir/two_group_compare.20190911/${clnV}_01 $basedir/mcmc_rlt.20190911/${clnV}_01 $basedir/mcmc_rlt.20190911/intrinsic_${Vgene}_01  $ss:$ee
 date) |& cat >$basedir/two_group_compare.20190911/${clnV}_01.Yyx_two_group_compare.${ss}_$ee.20190911.log &
done
wait

echo Merging  $basedir/two_group_compare.20190911/${clnV}_01.final_*.stat.txt
cp $basedir/two_group_compare.20190911/${clnV}_01.final_1_10.stat.txt  $basedir/two_group_compare.20190911/${clnV}_01.final.stat.txt
for ee in {20..300..10}; do
 ss=$((ee-9))
 # TODO: manually adjust above three lines to match previous setting; before the for loop is the first part, and the for loop should from the second part to the last part
 echo $ss:$ee
 tail -n+2 $basedir/two_group_compare.20190911/${clnV}_01.final_${ss}_$ee.stat.txt  >>$basedir/two_group_compare.20190911/${clnV}_01.final.stat.txt
done

done   # end for clnV
```

intermediate output to $basedir/two\_group\_compare.20190911/${clnV}\_01.final\_${ss}\_$ee.stat.txt

final output to $basedir/two\_group\_compare.20190911/${clnV}\_01.final.stat.txt

Then, check muDiff\_gt0\_1\_PEP column in Excel (muDiff\_gt0\_1\_PEP = PEP(muDiff > 0.1) = 1 - P_posterior(muDiff > 0.1) ). Each row for one nucleotide position in the V gene. Those sites with PEP(muDiff > 0.1) < 0.5 are considered 'significant'.

## Usage prompt

### scripts/yyx\_SHMPlot2\_Strata\_ErrBar.20190731.r 

```
Usage: Rscript yyx_SHMPlot2_Strata_ErrBar.20190731.r [OPTS] output statfile nucfile

Options(=defaults):
	-tstart=0	Start of reference to view
	-tend=0	End of reference to view
	-plotrows=1	Rows on plot
	-ymax=0.75	Maximum y-axis height
	-figureheight=2	height in inches
	-showsequence=FALSE	display sequence on plots
	-regex1=AGCT	
	-regex2=[AGT](?=G[CT][AT])	
	-cdr1_start=0	Start of cdr1 region
	-cdr1_end=0	End of cdr1 region
	-cdr2_start=0	Start of cdr2 region
	-cdr2_end=0	End of cdr2 region
	-cdr3_start=0	Start of cdr3 region
	-cdr3_end=0	End of cdr3 region
	-annotation=	V allele annotation (default: empty for no annotation)
	-mutMin=NA	Minimum mutation number for stratification
	-mutMax=NA	Maximum mutation number for stratification
	-mutMinProp=NA	Minimum mutation number proportion for stratification
	-mutMaxProp=NA	Maximum mutation number proportion for stratification
	-minReadNumB=0	Minimum read number required for each sample before stratification
	-minReadNumA=0	Minimum read number required for each sample after stratification
	-fo=FALSE	force output (default: stop if output file exists)

Arguments:
	output	file path prefix for output
	statfile	file path of one guide (same V) .stat.txt file - 1~2th columns: Pos, Base
	nucfile	file path of .nuc.txt file(s) - columns: Read_ID, base [ACGTN-.] on each position
```

### scripts/yyx\_stratify\_from\_nuc\_to\_merged\_stat.20190906.r

```
Usage: Rscript yyx_mut_rate_per_read.20190904.r [OPTS] output statfile nucfile

Options(=defaults):
	-mutMin=NA	Minimum mutation number for stratification
	-mutMax=NA	Maximum mutation number for stratification
	-mutMinProp=NA	Minimum mutation number proportion for stratification
	-mutMaxProp=NA	Maximum mutation number proportion for stratification
	-minReadNumB=1	Minimum read number required for each sample before stratification
	-minReadNumA=1	Minimum read number required for each sample after stratification
	-fo=FALSE	force output (default: stop if output file exists)

Arguments:
	output	file path prefix for output
	statfile	file path of one guide (same V) .stat.txt file - 1~2th columns: Pos, Base
	nucfile	file path of .nuc.txt file(s) - columns: Read_ID, base [ACGTN-.] on each position
```

### hierarchical\_Bayesian\_model\_jags\_SHM\_test.20190819/Yyx\_one\_group\_JAGS.20190820.r

```
Usage: Rscript Yyx_one_group_JAGS.20190820.r <output_prefix> <merged.stat.txt> [pos_range] [colname_prefix]
```

Note: it will use the Bayesian model setting in hierarchical\_Bayesian\_model\_jags\_SHM\_test.20190819/Yyx\_HierarPropBinom\_JAGS\_single.20190820.txt (in BUGS language)

### hierarchical\_Bayesian\_model\_jags\_SHM\_test.20190819/Yyx\_two\_group\_compare.20190903.r

```
Usage: Rscript Yyx_two_group_compare.20190903.r <output_prefix> <in1_foreground_prefix> <in2_background_prefix> [pos_range]
```
