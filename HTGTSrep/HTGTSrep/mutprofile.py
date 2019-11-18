#!/usr/bin/env python3
import os, re, sys, logging, csv, multiprocessing, time

import pandas as pd
from glob import glob
from itertools import groupby
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import OrderedDict, Counter

from HTGTSrep.lib import loggingRun, getInputSeq, getCDR, getGermdict, \
                         collapse_db, profile_DNAmut, profile_DNAmut_clonal

def filter_records(records, args):
    # Only keep joined reads with good V coverage
    records = records.loc[records['V_COVERAGE'] > args.min_Vcov, ]

    if not args.includeAllReads:
        records = records.loc[records['V_MUTATION'] > 0, ]

    if args.skipStopCodonRead:
        i = records[(records['STOP']=='T') & (records['IN_FRAME']=='T')].index
        records = records.drop(i)
    return records

def write_Vcov(records, Vcov_file):
    ''' Write V coverage file
    '''
    records_Vcov = records[['SEQUENCE_ID', 'V_ALLELE', 'V_GENE_LEN',
                                'V_COVERAGE', 'STOP', 'IN_FRAME']]
    records_Vcov = records_Vcov.sort_values('V_COVERAGE', ascending=True)
    records_Vcov.to_csv(Vcov_file, sep="\t", index=False)

def gen_mutstat(records, readtype):
    if readtype == 'P':
        records = records[records['PRODUCTIVE']=='T']
    elif readtype == 'Outofframe':
        i = records[(records['STOP']=='T') & (records['IN_FRAME']=='T')].index
        records = records.drop(i)
        records = records[records['PRODUCTIVE']=='F']
    elif readtype == 'IFwithSTOP':
        records = records[(records['STOP']=='T') & (records['IN_FRAME']=='T')]
    elif readtype == 'Total':
        pass
    ALIGN = records['V_ALIGNMENT'].sum()
    MUT = records['V_MUTATION'].sum()
    READS = len(records)
    if ALIGN == 0:
        RATE = "-"
    else:
        RATE = '%.3f' % float(MUT / ALIGN * 100)

    return [MUT, ALIGN, RATE, READS]

def write_mutstat(records, nuc_stat, rates_file):
    ''' Write mutation frequency file
    '''
    records_Vcov = records.loc[records['JOINED']=='T', ]
    nameslist = ["P_MUT_NUM", "P_V_MATCH", "P_MUT_RATE", "P_READS",
                "Outofframe_MUT_NUM", "Outofframe_V_MATCH", "Outofframe_MUT_RATE", "Outofframe_READS",
                "IFwithSTOP_MUT_NUM", "IFwithSTOP_V_MATCH", "IFwithSTOP_RATE", "IFwithSTOP_READS",
                "Total_MUT_NUM", "Total_V_MATCH", "Total_RATE", "Total_READS",
                "ALL_P_MUT_RATE", "ALL_Outofframe_MUT_RATE", "ALL_IFwithSTOP_RATE", "ALL_Total_RATE"]

    V_facts = pd.read_csv(nuc_stat, sep="\t", usecols=[0,1,2], low_memory=False)
    for name in nameslist: V_facts[name] = '-'

    # Count all reads rates
    ALL_P_MUT_RATE = gen_mutstat(records, "P")[2]
    ALL_Outofframe_MUT_RATE = gen_mutstat(records, "Outofframe")[2]
    ALL_IFwithSTOP_RATE = gen_mutstat(records, "IFwithSTOP")[2]
    ALL_Total_RATE = gen_mutstat(records, "Total")[2]

    # Count mut rates for each gene
    for key, group in records.groupby('V_GENE'):
        gene = group["V_GENE"].unique()[0]
        # Assign rates data
        rates = gen_mutstat(group, "P") + gen_mutstat(group, "Outofframe") + \
                gen_mutstat(group, "IFwithSTOP") + gen_mutstat(group, "Total") + \
                [ALL_P_MUT_RATE, ALL_Outofframe_MUT_RATE, ALL_IFwithSTOP_RATE, ALL_Total_RATE]

        V_facts.loc[V_facts["V_GENE"]==gene, nameslist] = rates
    V_facts.to_csv(rates_file, sep="\t", index=False)

def stat_Protein(protein_file, Vgapseq, allele, reads_db, VDJseq=''):
    ''' Prep protein stat file
    '''
    if allele not in Vgapseq:
        return
    else:
        Vseq_ref = Vgapseq[allele]
    # renew ref seq if it is VDJ seq fixed
    if VDJseq != '': Vseq_ref = VDJseq
    stat_file = protein_file.replace('protein.txt', 'pstat.txt')
    # Prep reference protein sequence
    refAA = ''
    pos_listAA = {}
    for i in range(0, len(Vseq_ref), 3):
        triNuc = Vseq_ref[i:i+3]
        if '.' in triNuc or len(triNuc) != 3:
            continue
        refAA = refAA + Seq(Vseq_ref[i:i+3], generic_dna).translate()[0]

    # Load protein file
    # colnames = [i for i in range(1, len(refAA)+1)]
    # reads_db = pd.DataFrame(columns=colnames)
    # for line in open(protein_file):
    #     if not line.startswith('>'):
    #         reads_db.loc[len(reads_db)+1] = [s for s in line.strip()]

    # Gen stat file
    stat_handle = open(stat_file, 'w')
    stat_handle.write('Pos\tMut\tTotal\tBase\tY\tProbability\n')
    readsnum = len(reads_db)
    for i in range(1, len(refAA)+1):
        base = refAA[i-1]
        vcount = reads_db.ix[:,i].value_counts().to_dict()
        misscount = vcount.get('-', 0) + vcount.get('*', 0)
        germcount = vcount.get(base, 0)
        total = readsnum - misscount
        mut = total - germcount
        y = float(mut) / total if total > 0 else 0
        if total > 0:
            probability = ';'.join(['%s:%.2f%%' % (l, float(vcount[l])/total*100) for l in vcount if l not in '-*'])
        else:
            probability = '-'
        stat = '%d\t%d\t%d\t%s\t%.3f\t%s\n' % (i, mut, total, base, y, probability)
        stat_handle.write(stat)
    stat_handle.close()

def profile_Protein(group, protein_file, protein_PDF, Vgapseq, args):
    ''' Prep protein text and PDF file
    '''
    def _parse_peptide_seq(row):
        # sub-function: parse row to obtain protein seq
        readseq = ''
        tripletlist = []
        # Prep triplet array for the seq base on if exists IMGT gaps
        if row['GERMLINE_IMGT_D_MASK'] != '-' and row['SEQUENCE_IMGT'] != '-':
            SEQUENCE_IMGT = row['SEQUENCE_IMGT']
            vlen = int(row['V_SEQ_LENGTH'])
            count = 0
            readseq_gap = ''
            for i in range(0, len(SEQUENCE_IMGT)):
                if count < vlen:
                    readseq_gap += SEQUENCE_IMGT[i]
                    if SEQUENCE_IMGT[i] in 'ATCGN': count += 1
                else:
                    break

            # Trim or complete read V part
            V_GENE_GAP_LEN = int(row['V_GENE_GAP_LEN'])
            if V_GENE_GAP_LEN > len(readseq_gap):
                readseq_gap += '.' * (V_GENE_GAP_LEN - len(readseq_gap))
            else:
                readseq_gap = readseq_gap[0: V_GENE_GAP_LEN]

            for i in range(0, len(readseq_gap), 3):
                # Pass the triplet if exist gap in the germline seq
                if allele in Vgapseq:
                    Vseq_ref = Vgapseq[allele]
                else:
                    Vseq_ref = row['GERMLINE_IMGT_D_MASK']
                # this can be a potential bug as it only trim some nucleic acids
                # at the beginning, luckily it is always triplet within the sequences
                if '.' in Vseq_ref[i:i+3]:
                    continue
                tripletlist.append(readseq_gap[i:i+3])
        else:
            for i in range(0, allele_len):
                p = row['V_ALLELE_NUC'][i]
                if p == '.': readseq += allele_seq[i]
                elif p in 'ATCG': readseq += p
                elif p in 'N-': readseq += '-'
            tripletlist = [readseq[i:i+3] for i in range(0, len(readseq), 3)]

        # Translate the triplet array to protein sequence
        readprotein = ''
        for triplet in tripletlist:
            if len(triplet) != 3: continue
            if '.' in triplet or 'N' in triplet or '-' in triplet:
                readprotein += '-'
            else:
                readprotein += Seq(triplet, generic_dna).translate()[0]
        return readprotein

    allele = group["V_ALLELE"].unique()[0]
    germ_dict = getGermdict(args)
    allele_seq = germ_dict[allele]
    allele_len = len(allele_seq)

    protein_df = pd.DataFrame(group["SEQUENCE_ID"])
    protein_df['READPROTEIN'] = group.apply(_parse_peptide_seq, axis=1)

    # Prep protein file and run weblogo
    protein_handle = open(protein_file, 'w')
    for key, row in protein_df.iterrows():
        protein_handle.write('>%s\n%s\n' % (row["SEQUENCE_ID"], row["READPROTEIN"]))
    protein_handle.close()
    if len(group) >= args.min_profileread:
        os.system('weblogo -f %s -o %s -D fasta --units probability -A protein ' \
                  '--composition none --size large -n 100 --scale-width NO ' \
                  '--color-scheme chemistry -S %f' % (protein_file,
                  protein_PDF, args.ymax_protein))

    # Prep reads_db df for additional stat analysis
    proteinLen = len(protein_df['READPROTEIN'].tolist()[0])
    reads_db = pd.DataFrame(columns=[i for i in range(0, proteinLen)], index=range(0, len(protein_df)))
    reads_db = protein_df['READPROTEIN'].apply(lambda x: pd.Series([p for p in x], index=range(1, proteinLen+1)))
    return reads_db

def run_profile(prod_type, group, sample, Vgapseq, args):
    prod_dict = {'NP':'F', 'P':'T'}
    group = group.loc[group['PRODUCTIVE']==prod_dict[prod_type],]
    if len(group) == 0: return

    allele = group["V_ALLELE"].unique()[0]
    ValleleN = allele.replace('/','-')
    # Generate DNA profiles
    collapse_tag = '_collapse' if (args.collapsePartial or args.collapseIdentical) else ""
    nuc_profile = "%s/%s/mut_profile/nucl_text%s/%s.%s.%s.nuc.txt" % (
                args.dir, sample, collapse_tag, sample, ValleleN, prod_type)
    nuc_stat = "%s/%s/mut_profile/nucl_text%s/%s.%s.%s.stat.txt" % (
                args.dir, sample, collapse_tag, sample, ValleleN, prod_type)
    nuc_PDF = "%s/%s/mut_profile/nucl_profile%s/%s.%s.%s.pdf" % (
                args.dir, sample, collapse_tag, sample, ValleleN, prod_type)
    nuc_profile = re.sub(r'[*)(]', '-', nuc_profile)
    nuc_stat = re.sub(r'[*)(]', '-', nuc_stat)
    nuc_PDF = re.sub(r'[*)(]', '-', nuc_PDF)
    profile_DNAmut(group, nuc_stat, nuc_PDF, nuc_profile, args)

    # Generate Protein profiles
    protein_file = "%s/%s/mut_profile/protein_text%s/%s.%s.%s.protein.txt" % (
                    args.dir, sample, collapse_tag, sample, ValleleN, prod_type)
    protein_PDF = "%s/%s/mut_profile/protein_profile%s/%s.%s.%s.eps" % (
                    args.dir, sample, collapse_tag, sample, ValleleN, prod_type)
    protein_file = re.sub(r'[*)(]', '-', protein_file)
    protein_PDF = re.sub(r'[*)(]', '-', protein_PDF)
    reads_db = profile_Protein(group, protein_file, protein_PDF, Vgapseq, args)

    stat_Protein(protein_file, Vgapseq, allele, reads_db)

def run_profile_fixVDJ(prod_type, group, sample, Vgapseq, args):
    ''' Generate mut profile using concat VDJ seq as reference.
        To note: V+D+J genes are concated directly, so ensure V,D,J genes are
        correct in reference fasta files
    '''
    if prod_type != 'ALL':
        prod_dict = {'NP':'F', 'P':'T'}
        group = group.loc[group['PRODUCTIVE']==prod_dict[prod_type],]

    if len(group) == 0: return

    Vallele = group["V_ALLELE"].unique()[0]
    Dallele = group["D_ALLELE"].unique()[0]
    Jallele = group["J_ALLELE"].unique()[0]
    VDJseq = (Vgapseq[Vallele] + Vgapseq[Dallele] + Vgapseq[Jallele]).replace('.', 'N')

    # Generate DNA profiles
    nuc_profile = "%s/%s/mut_profile/nucl_text_fixedVDJ/%s.%s.%s.nuc.txt" % (
                args.dir, sample, sample, ValleleN, prod_type)
    nuc_stat = "%s/%s/mut_profile/nucl_text_fixedVDJ/%s.%s.%s.stat.txt" % (
                args.dir, sample, sample, ValleleN, prod_type)
    nuc_PDF = "%s/%s/mut_profile/nucl_profile_fixedVDJ/%s.%s.%s.pdf" % (
                args.dir, sample, sample, ValleleN, prod_type)
    nuc_profile = re.sub(r'[*)(]', '-', nuc_profile)
    nuc_stat = re.sub(r'[*)(]', '-', nuc_stat)
    nuc_PDF = re.sub(r'[*)(]', '-', nuc_PDF)
    profile_DNAmut_clonal(VDJseq, group, nuc_stat, nuc_PDF, nuc_profile, args)

    # Generate Protein profiles
    protein_file = "%s/%s/mut_profile/protein_text_fixedVDJ/%s.%s.%s.protein.txt" % (
                    args.dir, sample, sample, ValleleN, prod_type)
    protein_PDF = "%s/%s/mut_profile/protein_profile_fixedVDJ/%s.%s.%s.eps" % (
                    args.dir, sample, sample, ValleleN, prod_type)
    protein_file = re.sub(r'[*)(]', '-', protein_file)
    protein_PDF = re.sub(r'[*)(]', '-', protein_PDF)
    reads_db = profile_Protein(group, protein_file, protein_PDF, Vgapseq, args)

    stat_Protein(protein_file, Vgapseq, Vallele, reads_db, VDJseq)

def mutProfile(args):
    '''
    '''
    logging.info('Parsing DNA & protein mutation profile')

    # Load gapped fasta if existed
    Vgapseq = getInputSeq(args.__dict__['params_dict']['Vgapseq'])

    for sample in args.metadict:
        sampledir = '{0}/{1}'.format(args.dir, sample)

        # Read IgBlast Db and add new sequence feature
        db_path = '{0}/{1}.db.xls'.format(sampledir, sample)
        records = pd.read_csv(db_path, sep="\t", low_memory=False)
        if args.collapsePartial:
            records = collapse_db(records, 'partial', 'F')
            records.to_csv('{0}/mut_profile/{1}.collapse.xls'.format(sampledir, sample), sep="\t", index=False)
        if args.collapseIdentical:
            records = collapse_db(records, 'identical', 'F')
            records.to_csv('{0}/mut_profile/{1}.collapse.xls'.format(sampledir, sample), sep="\t", index=False)
        # Drop D upstream reads
        if 'D_UPSTREAM_MATCH_D_GENE' in records.columns:
            records = records.loc[records['D_UPSTREAM_MATCH_D_GENE'] == '-', ]

        # Write V coverage stat file
        Vcov_file = '{0}/mut_profile/{1}.Vcov.xls'.format(sampledir, sample)
        if (args.collapsePartial or args.collapseIdentical):
            Vcov_file = '{0}/mut_profile/{1}.Vcov.collapse.xls'.format(sampledir, sample)
        write_Vcov(records, Vcov_file)

        # Filter out partial records
        records = filter_records(records, args)
        records['SAMPLE'] = sample

        # Write mutation rate stat file
        statpath = '%s/stat/allsample.stat_join.xls' % args.dir
        rates_file = '{0}/mut_profile/{1}.mutfreq.xls'.format(sampledir, sample)
        if (args.collapsePartial or args.collapseIdentical):
            rates_file = '{0}/mut_profile/{1}.mutfreq.collapse.xls'.format(sampledir, sample)
        write_mutstat(records, statpath, rates_file)

        # DNA & protein profile on V alleles in productive and non-productive seperately
        for key, group in records.groupby('V_ALLELE'):
            run_profile('NP', group, sample, Vgapseq, args)
            run_profile('P', group, sample, Vgapseq, args)

    logging.info('All done. Have a good day!')
