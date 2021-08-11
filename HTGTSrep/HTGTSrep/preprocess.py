#!/usr/bin/env python3
import os, sys, multiprocessing, logging
from HTGTSrep.lib import loggingRun, reverse_complement, collapse_fasta

def _preprocess_fasta(args):
    ''' create directories and copy fasta into each one'''
    fastafiles = args.input.split(',')
    for i in range(0, len(args.metadict)):
        sample = list(args.__dict__['metadict'].keys())[i]
        eachdir = '%s/%s' % (args.outdir, sample)
        if not os.path.exists(eachdir):
            os.system('mkdir -p %s/reads_fasta' % eachdir)
        os.system('cp %s %s/%s_join.fa' % (fastafiles[i], eachdir, sample))

def _sample_demultiplex(args):
    '''Preprocess
    1) Generate barcode file using first <demulti_length> bps of MID+Primer
    2) demultiplex using fastq-multx program
    '''
    if not args.skipDemultiplex:   # 2021-05-19, Adam_Yyx add this skipDemultiplex condition
        # create barcode file
        fbar = open('%s/barcodes.txt' % args.outdir, 'w')
        for sample in args.metadict:
            primer = args.metadict[sample][1]
            primerSet = args.metadict[sample][1].strip().split(',')
            if len(primerSet) == 1:
                barcode = (args.metadict[sample][0] + \
                           args.metadict[sample][1])[:args.demulti_length]
                fbar.write('%s\t%s\n' % (sample, barcode))
            elif len(primerSet) > 1:
                for i in range(0, len(primerSet)):
                    barcode = (args.metadict[sample][0] + primerSet[i])[:args.demulti_length]
                    fbar.write('%s_%d\t%s\n' % (sample, i, barcode))
        fbar.close()

        # create directories if not exist
        for sample in args.metadict:
            eachdir = '%s/%s' % (args.outdir, sample)
            if not os.path.exists(eachdir):
                os.system('mkdir %s' % eachdir)

        # gunzip raw fastq if gzipped
        if (args.r1).endswith('gz') or (args.r1).endswith('gzip'):
            os.system('gzip -c -d %s > %s/raw.r1.fq' % (args.r1, args.outdir))
            os.system('gzip -c -d %s > %s/raw.r2.fq' % (args.r2, args.outdir))
        else:
            os.system('cp %s %s %s' % (args.r1, args.r2, args.outdir))
        # demultiplex
        try:
            # 02/01/2021 JH: added {2} --demulti_distance
            loggingRun('fastq-multx -m {0} -x -b -d {2} -B {1}/barcodes.txt {1}/raw.r1.fq ' \
                        '{1}/raw.r2.fq -o {1}/%_R1.fq {1}/%_R2.fq'.format(
                        args.demulti_mismatch, args.outdir, args.demulti_distance))
        except:
            logging.error('Cannot run fastq-multx to demultiplex samples\n')
            sys.exit(-1)

    # clean raw data and move files to diff fold
    for sample in args.metadict:
        os.system('mv {0}/{1}_*R* {0}/{1}'.format(args.outdir, sample, args.outdir, sample))
    os.system('rm -rf {0}/raw.r1.fq {0}/raw.r2.fq'.format(args.outdir))
    if args.keepunmatch:
        os.system('mkdir %s/unmatched' % args.outdir)
        os.system('mv {0}/unmatched_* {0}/unmatched'.format(args.outdir))
    else:
        os.system('rm -rf %s/unmatched_*' % args.outdir)

def _primer_process(barcode, primer, sample, subsample, args):
    ''' Process each primer sequence
    '''
    eachdir = '%s/%s' % (args.outdir, sample)
    primertrim = reverse_complement((barcode + primer))[:10]
    adapter = args.metadict[sample][2].strip()
    adapter_rc = reverse_complement(adapter).strip()
    logfile_R1 = '%s/logs/pipeline.preprocess.trimR1.log' % args.outdir
    logfile_R2 = '%s/logs/pipeline.preprocess.trimR2.log' % args.outdir
    # 1. trim adapter at R2 5' end and seq after primer at R2 3' end
    try:
        loggingRun('cutadapt --quiet -g %s -a %s -n 2 -m 50 -o %s/%s_R2.trim.fq %s/%s_R2.fq >> %s' % (
            adapter_rc.strip(), primertrim, eachdir, subsample, eachdir, subsample, logfile_R2))
    except:
        logging.error('Cannot run cutadapt to trim adapter')
        sys.exit(-1)

    # 2. trim barcode at R1 5' end adapter at R1 3' end
    if barcode != '':
        loggingRun('cutadapt --quiet -g ^%s -a %s -n 2 -m 50 -o %s/%s_R1.trim.fq %s/%s_R1.fq >> %s' % (
                barcode, adapter, eachdir, subsample, eachdir, subsample, logfile_R1))
    else:
        loggingRun('cutadapt --quiet -m 50 -a %s -o %s/%s_R1.trim.fq %s/%s_R1.fq >> %s' % (
                adapter, eachdir, subsample, eachdir, subsample, logfile_R1))

    # 3. join R1 and R2 fastq files, and convert to fasta files
    loggingRun("awk '{{print $1}}' {0}/{1}_R1.trim.fq > {0}/{1}_R1.clean.fq".format(
                eachdir, subsample))
    loggingRun("awk '{{print $1}}' {0}/{1}_R2.trim.fq > {0}/{1}_R2.clean.fq".format(
                eachdir, subsample))
    loggingRun("awk '{{ if(NR % 4 == 1) {{print $1}} }}' {0}/{1}_R1.clean.fq |" \
                " sed 's/@//' | sort > {0}/{1}_R1.reads.list".format(
                eachdir, subsample))
    loggingRun("awk '{{ if(NR % 4 == 1) {{print $1}} }}' {0}/{1}_R2.clean.fq |" \
                "sed 's/@//' | sort > {0}/{1}_R2.reads.list".format(
                eachdir, subsample))
    loggingRun('comm -12 {0}/{1}_R1.reads.list {0}/{1}_R2.reads.list > {0}/{1}_comm.reads.list'.format(
                eachdir, subsample))
    loggingRun('seqtk subseq {0}/{1}_R1.clean.fq {0}/{1}_comm.reads.list > {0}/{1}_R1.even.fq'.format(
                eachdir, subsample))
    loggingRun('seqtk subseq {0}/{1}_R2.clean.fq {0}/{1}_comm.reads.list > {0}/{1}_R2.even.fq'.format(
                eachdir, subsample))
    try:
        loggingRun('fastq-join -p {0} -m {1} {2}/{3}_R1.even.fq {2}/{3}_R2.even.fq ' \
                    '-o {2}/{3}_unjoinR1.fq -o {2}/{3}_unjoinR2.fq -o {2}/{3}_join.fq'.format(
                    args.diffpercent, args.overlap, eachdir, subsample))
    except:
        logging.error('Cannot run fastq-join to join read1 and read2')
        sys.exit(-1)

    # 4. Convert fq to fa and keep long reads
    for subset in ['join', 'unjoinR1', 'unjoinR2']:
        prefix = "%s/%s_%s" % (eachdir, subsample, subset)
        try:
            if subset == 'unjoinR1':
                if args.fastq_quality_trimmer > 0:
                    loggingRun("fastq_quality_trimmer -t %d -i %s.fq -Q33 | " \
                                "fastq_to_fasta -Q33 -n | awk '{print $1}' > %s.fa" % (
                                args.fastq_quality_trimmer, prefix, prefix))
                else:
                    loggingRun("fastq_to_fasta -i %s.fq -Q33 -n | awk '{print $1}' > %s.fa" % (
                            prefix, prefix))
            elif subset == 'unjoinR2':
                loggingRun("fastq_to_fasta -i %s.fq -Q33 -n | awk '{print $1}' > %s.fa" % (
                        prefix, prefix))
            else:
                cmdline = "fastq_masker -q %d -Q33 -i %s.fq | " \
                            "fastq_quality_filter -q %d -p %d -Q33 | " \
                            "fastq_to_fasta -Q33 -n | awk '{print $1}' > %s.fa" \
                            % (args.qscore, prefix, args.qscore,
                            args.qscore_cov, prefix)
                loggingRun(cmdline)
        except:
            logging.error('Cannot run FASTX-Toolkit to convert fq to fa format')
            sys.exit(-1)
    #collapse_fasta("%s/%s" % (eachdir, subsample))

def _sample_process(sample, args):
    ''' Perform preprocessing for each sample
    '''
    barcode = args.metadict[sample][0]
    # Use a loop here to process all primer sequences.
    # TO-DO: Need a better way to do this
    primerSet = args.metadict[sample][1].strip().split(',')
    if len(primerSet) == 1:
        _primer_process(barcode, primerSet[0], sample, sample, args)
    elif len(primerSet) > 1:
        subsampleList = ['%s_%d' % (sample, i) for i in range(0, len(primerSet))]
        for i in range(0, len(subsampleList)):
            subsample = subsampleList[i]
            primer = primerSet[i]
            _primer_process(barcode, primer, sample, subsample, args)
        # Pool all sub-samples and run the downstream analysis
        suffixR = ['_unjoinR1.fa', '_unjoinR1.fq', '_unjoinR2.fa',
                    '_unjoinR2.fq', '_comm.reads.list', '_join.fa', '_join.fq']
        for suffix in ['.clean.fq', '.even.fq', '.fq', '.reads.list', '.trim.fq']:
            suffixR.append('_R1' + suffix)
            suffixR.append('_R2' + suffix)

        for suffix in suffixR:
            sl = ['%s/%s/%s%s' % (args.outdir, sample, subsample, suffix) for subsample in subsampleList]
            mergefile = '%s/%s/%s%s' % (args.outdir, sample, sample, suffix)
            os.system('cat %s > %s' % (' '.join(sl), mergefile))
            os.system('rm -rf %s' % ' '.join(sl))

def reads_process(args):
    '''For paired fastq files, will first demultiplex samples using metadata file,
    then remove low quality reads, and join the paired reads. For fasta file, the
    program will copy each fasta file to individual directory.
    '''
    if args.input is None:
        # 1. demultiplex samples
        logging.info('Preprocessing read files.....')
        _sample_demultiplex(args)

        # 2. trim adapter and junk sequences & join R1 and R2 reads
        logging.info('Joining reads......')
        pool = multiprocessing.Pool(processes = len(args.metadict) )
        for sample in args.metadict:                 #_sample_process(sample, args)
            pool.apply_async(_sample_process, (sample, args, ))
        pool.close()
        pool.join()
    else:
        _preprocess_fasta(args)
