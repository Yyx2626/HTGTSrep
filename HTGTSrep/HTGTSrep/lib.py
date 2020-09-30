#!/usr/bin/env python3
import random, os, sys, logging, re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import numpy as np
from itertools import combinations, product
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.stats import sem

def writeOutputSeq(col_ID, col_seq, output_file, db):
    """
    Generate output fasta file using db dataframe
    """
    output_handle = open(output_file, 'w')
    for key, row in db.iterrows():
        output_handle.write('>%s\n%s\n' % (row[col_ID], row[col_seq]))
    protein_handle.close()

def collapse_fasta(prefix):
    print(prefix)
    # parse joined reads
    seq_dict_joined = getInputSeq(prefix + '_join.raw.fa')
    seq_id_dict = {}
    for id in seq_dict_joined:
        seq = seq_dict_joined[id]
        if seq not in seq_id_dict:
            seq_id_dict[seq] = [id]
        else:
            seq_id_dict[seq].append(id)
    fjoined = open((prefix + '_join.fa'), 'w')
    for seq in seq_id_dict:
        fjoined.write('>%s\n%s\n' % ('_'.join(seq_id_dict[seq]), seq))
    fjoined.close()

    # parse unjoined
    seq_dict_R1 = getInputSeq(prefix + '_unjoinR1.raw.fa')
    seq_dict_R2 = getInputSeq(prefix + '_unjoinR2.raw.fa')
    seq_id_dict_R1R2 = {}  # concated seq: [seq IDs]
    for id in seq_dict_R1:
        concat_seq = seq_dict_R1[id] + seq_dict_R2[id]
        if concat_seq not in seq_id_dict_R1R2:
            seq_id_dict_R1R2[concat_seq] = [id]
        else:
            seq_id_dict_R1R2[concat_seq].append(id)
    fR1 = open(prefix + '_unjoinR1.fa', 'w')
    fR2 = open(prefix + '_unjoinR2.fa', 'w')
    for seq in seq_id_dict_R1R2:
        fR1.write('>%s\n%s\n' % ('_'.join(seq_id_dict_R1R2[seq]), seq_dict_R1[seq_id_dict_R1R2[seq][0]]))
        fR2.write('>%s\n%s\n' % ('_'.join(seq_id_dict_R1R2[seq]), seq_dict_R2[seq_id_dict_R1R2[seq][0]]))
    fR1.close()
    fR2.close()


def getInputSeq(seq_file):
    """
    Arguments:
    seq_file = a fasta file of sequences input

    Returns:
    a dictionary of {ID:Seq}
    """
    ### add print message to warn for empty dict
    if not os.path.exists(seq_file):
        print("%s FAILED TO LOAD. EMPTY DICT IS RETURNED. THIS MAY INFLUENCE YOUR RESULTS" % seq_file, file=sys.stderr)
        return {}

    if seq_file.endswith('.gz'):
        os.system('gunzip %s' % seq_file)
        seq_file_unzip = seq_file.rstrip('.gz')
    else:
        seq_file_unzip = seq_file

    seq_dict = SeqIO.index(seq_file_unzip, "fasta", IUPAC.ambiguous_dna)

    # Create a seq_dict ID translation using IDs truncate up to space or 50 chars
    seqs = {}
    for seq in seq_dict.values():
        seqs.update({seq.description: str(seq.seq).upper()})

    ### .fa files may have a header preceeding each gene. This chunk is added to make sure the header is removed
    ### can't change the brackets, otherwise keyerror
    if ".fa" in seq_file_unzip:
        keys = list(seqs.keys())
        # obtain a list of keys stripped of the header
        for i in range(len(keys)):
            keys[i] = keys[i].replace("lcl|", "", 1)
        seqs = dict(zip(keys, list(seqs.values())))

    if seq_file.endswith('.gz'):
        os.system('gzip %s' % seq_file_unzip)

    return seqs

def getCDR(cdrfile):
    V_CDR = {}
    if not os.path.exists(cdrfile):
        logging.warnings('Cannot find CDR boundary file %s' % os.path.basename(cdrfile))
        return None
    else:
        for line in open(cdrfile):
            l = line.strip().split()
            V_CDR[l[0]] = [int(b) for b in l[1:]]
        return V_CDR

def getCSV(csvfile):
    # Load CSV file by reading by chunks
    tmplist = []
    for chunk in pd.read_csv(csvfile, sep='\t', chunksize=20000):
        tmplist.append(chunk)
    m = pd.concat(tmplist, axis=0)
    del tmplist
    return m

def load_Valign(fname):
    # Load V gene genome alignment position
    V_align = {}
    for line in open(fname):
        l = line.strip().split()
        start_end = '%s_%s' % (l[2], l[3])
        if l[0] not in V_align:
            V_align[l[0]] = {l[1]: [start_end]}
        else:
            if l[1] in V_align[l[0]]:
                V_align[l[0]][l[1]].append(start_end)
            else:
                V_align[l[0]][l[1]] = [start_end]
    return V_align

def CheckAlignOverlap(topinfo, reads_align, Valign, genomealign, hitcheck):
    flag = 'noneed'
    if genomealign == 'T':
        flag = 'unmatch'
        if topinfo[0] not in reads_align:
            flag = 'nohit'
        else:
            for loc in reads_align[topinfo[0]]:
                chrom = loc[0]
                pos = int(loc[1])
                if topinfo[1] not in Valign:
                    flag = 'noVhit'
                    continue
                if chrom in Valign[topinfo[1]]:
                    for start_end in Valign[topinfo[1]][chrom]:
                        start = int(start_end.split('_')[0])
                        end = int(start_end.split('_')[1])
                        # extend 10bp at 5' because V-D or V-J junctions might have matches
                        if (start - 10) <= pos <= end:
                            flag = 'match'
    if flag == 'nohit':
        return 'No_hit_from_genome_alignment'
    elif flag == 'noVhit':
        return 'topVgene_has_no_alignment'
    elif flag == 'unmatch':
        return 'genome_alignment_unmatch_Vgene'
    else:
        return hitcheck

def loggingRun(cmdline):
    logging.info(cmdline)
    os.system(cmdline)

def line_count(fname):
    i = -1
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def fasta_count(fname):
    i = 0
    fin = open(fname)
    for line in fin:
        if line.startswith('>'):
            i += 1
    fin.close()
    return i


def reads_stat(args):
    for sample in sample_info:
        eachdir = '%s/%s' % (args.outdir, sample)
        fstat = open('%s/%s.stat.txt' % (eachdir, sample), 'w')
        # reads count
        total_num = line_count("%s/%s_R1.fq" % (eachdir, sample)) / 4
        join_num = line_count("%s/%s_join.fq" % (eachdir, sample)) / 4
        unjoin_num = total_num - join_num
        fstat.write('Total reads\t%d\nJoined reads\t%d\nUnjoined reads\t%d\n' % (total_num, join_num, unjoin_num))

        # alignment stat
        join_uniq = line_count('%s/%s_join.uniq.xls' % (eachdir, sample))
        R1_uniq = line_count('%s/%s_unjoinR1.uniq.xls' % (eachdir, sample))
        join_NOuniq = line_count('%s/%s_join.NOuniq.xls' % (eachdir, sample))
        R1_NOuniq = line_count('%s/%s_unjoinR1.NOuniq.xls' % (eachdir, sample))
        mergeNum = line_count('%s/%s.IgBlast_merge.xls' % (eachdir, sample))
        fstat.write('# of uniquely/NON-uniquely joined hits\t%d\t%d\n' % (join_uniq, join_NOuniq))
        fstat.write('# of uniquely/NON-uniquely unjoined-R1 hits\t%d\t%d\n' % (R1_uniq, R1_NOuniq))
        fstat.write('# of merged hits\t%d\n' % mergeNum)
        fstat.close()

def random_seq(length):
    ''' Generate random sequnce with input length '''
    seq = ''
    if length == 0:
        return seq
    else:
        seq = ''.join([random.choice('ATCG') for i in range(0, length)])
    return seq

def mutate_seq(orig_string, mutation_rate=0.005):
    ''' Mutate input sequence with point mutations '''
    bases = "ACGT"
    result = []
    mutations = []
    n = 0
    for base in orig_string:
        n += 1
        if random.random() < mutation_rate and base in bases:
            new_base = bases[bases.index(base) - random.randint(1, 3)]  # negatives are OK
            result.append(new_base)
            mutations.append('%s%d%s' % (base, n, new_base))
        else:
            result.append(base)
    return "".join(result), '|'.join(mutations)

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq_rc = "".join(complement.get(base, base) for base in reversed(seq))
    return seq_rc

def fastq_stats(fastqfile):
    # execute fastq-stats
    os.system('fastq-stats %s > %s.fastqstat' % (fastqfile, fastqfile))
    # parse results
    fqstat = {}
    for line in open('%s.fastqstat' % fastqfile):
        l = line.strip().split('\t')
        fqstat[l[0]] = l[1]
    os.system('rm -rf %s.fastqstat' % fastqfile)
    return fqstat

def parsefa_long(file, length):
    id_seq = {}
    id = ''
    for line in open(file):
        if line.startswith('>'):
            id = line.strip()
            id_seq[id] = ''
        else:
            id_seq[id] += seq
    fout = open(file.replace('.fa', '.long.fa'), 'w')
    for id in id_seq:
        if len(id_seq[id]) >= length:
            fout.write('%s\n%s\n' % (id, id_seq[id]))
    fout.close()

def smooth(self, nucMutnum, nucCovnum, genetotalseq, statfilelist):
    # print(nucMutnum['A'], nucMutnum['G'], nucMutnum['C'], nucMutnum['T'])
    # print(nucCovnum['A'], nucCovnum['G'], nucCovnum['C'], nucCovnum['T'])
    nucMucratio = {}
    smoothpower = self.args.additivesmooth
    for nuc in 'AGCT':
        nucMucratio[nuc] = float(nucMutnum[nuc]) / nucCovnum[nuc]
    avecover = sum([nucCovnum[a] for a in 'AGCT']) / len(genetotalseq)
    for gene in statfilelist:
        statfile = statfilelist[gene]
        statnew = statfile.replace('.txt', '.sm%s.txt' % str(smoothpower))
        fnew = open(statnew, 'w')
        for line in open(statfile):
            if line.startswith('Pos'):
                fnew.write(line)
            else:
                l = line.strip().split('\t')
                total_smooth = int(l[2]) + avecover * smoothpower
                mut_smooth = int(l[1]) + nucMucratio[nuc] * avecover * smoothpower
                if total_smooth == 0:
                    l[4] = 0
                else:
                    l[4] = mut_smooth / total_smooth
                l[4] = str(l[4])
                fnew.write('%s\n' % '\t'.join(l))
        fnew.close()
        pdffile = statnew.replace('nucfile', 'profiles').replace('txt', 'pdf')
        ### 09152020: changed showsequence from false to ture
        loggingRun(
            'Rscript scripts/SHMPlot2.R %s %s plotrows=1 figureheight=2 showsequence=TRUE ymax=0.2 cdr1_start=%d cdr1_end=%d cdr2_start=%d cdr2_end=%d cdr3_start=%d cdr3_end=%d' % \
            (statnew, pdffile, self.V_CDR[gene]['CDR1_start'], self.V_CDR[gene]['CDR1_end'], \
             self.V_CDR[gene]['CDR2_start'], self.V_CDR[gene]['CDR2_end'], \
             self.V_CDR[gene]['CDR3_start'], self.V_CDR[gene]['CDR3_end']))

######### this section is for tree construction & file parsing

def mergeSampleCount(shortlist):
    samplelist = [s.split(':')[0] for s in shortlist[0].split('|')]
    sample_count = {}
    for s in samplelist:
        sample_count[s] = 0
    for shortcount in shortlist:
        for oneshort in shortcount.split('|'):
            (a, b) = oneshort.split(':')
            sample_count[a] = sample_count[a] + int(b)
    o = '|'.join(["%s:%d" % (a, sample_count[a]) for a in samplelist])
    return o

def treeCollapseParse(fin, fout):
    db = pd.read_csv(fin, sep="\t", low_memory=False)
    if len(db) < 2: sys.exit('Find no passed read in tmp_db-pass.tab')

    grouped = db.groupby('CLONE')
    idlist = []
    sclist = []
    readlist = []
    fullseqlist = []
    for key, group in grouped:
        seqlist = []
        group = pd.DataFrame(group)
        germseq = list(group['GERMLINE_IMGT_D_MASK'])[0]
        for si in group['SEQUENCE_IMGT']:
            s = []
            for n in range(0, len(si)):
                if si[n] in ['N', '.'] and germseq[n] != 'N':
                    s.append(germseq[n])
                else:
                    s.append(si[n])
            seqlist.append(''.join(s))
        group["FULLSEQ"] = seqlist
        grouped2 = group.groupby("FULLSEQ")
        for subkey, subgroup in grouped2:
            subgroup = pd.DataFrame(subgroup)
            subgroup["trimlen"] = [len(s.replace('.', '').replace('N', '')) for s in subgroup['SEQUENCE_IMGT']]
            subgroup = subgroup.sort_values("trimlen", ascending=False)
            idlist.append(list(subgroup['SEQUENCE_ID'])[0])
            fullseqlist.append(list(subgroup['FULLSEQ'])[0])
            readlist.append('|'.join(list(subgroup['SEQUENCE_ID'])))
            sclist.append(mergeSampleCount(list(subgroup['SHORTCOUNT'])))
    treeCollapse = pd.DataFrame(db.loc[db['SEQUENCE_ID'].isin(idlist),])
    treeCollapse["SHORTCOUNT"] = sclist
    # treeCollapse["SEQUENCE_IMGT"] = fullseqlist
    # treeCollapse["READGROUP"] = readlist
    treeCollapse.to_csv(fout, sep="\t", index=False)


def files_process(args, worktype):
    # IgBlast clean up
    if worktype == 'igblast_clean':
        for sample in args.metadict:
            eachdir = '%s/%s' % (args.outdir, sample)
            dirlist = ['reads_fasta', 'reads_fastq', 'igblast_raw',
                       'igblast_db']  # , 'bowtie_sam']
            for d in dirlist:
                if not os.path.exists('%s/%s' % (eachdir, d)):
                    os.system('mkdir %s/%s' % (eachdir, d))
            os.system('mv {0}/*fa {0}/reads_fasta'.format(eachdir))
            os.system('mv {0}/*.fq {0}/*list {0}/reads_fastq'.format(eachdir))
            os.system('mv {0}/*IgBlast {0}/igblast_raw'.format(eachdir))
            os.system('mv {0}/*IgBlast.db {0}/igblast_db'.format(eachdir))
            # if args.genomealign == 'T':
            #     os.system('mv %s/*.sam %s/bowtie_sam' % (eachdir, eachdir))
            os.system('gzip %s/reads_fast*/*' % (eachdir))
            os.system('gzip %s/igblast*/*' % eachdir)
        if os.path.exists('%s/unmatched/' % args.outdir):
            os.system('gzip %s/unmatched/*' % args.outdir)

def getNmers(sequences, n):
    """
    Breaks input sequences down into n-mers

    Arguments:
      sequences : List of sequences to be broken into n-mers
      n : Length of n-mers to return

    Returns:
      dict : Dictionary mapping sequence to a list of n-mers
    """
    # Add Ns so first nucleotide is center of first n-mer
    sequences_n = ['N' * ((n - 1) // 2) + seq + 'N' * ((n - 1) // 2) for seq in sequences]
    nmers = {}
    for seq, seqn in zip(sequences, sequences_n):
        nmers[seq] = [seqn[i:i + n] for i in range(len(seqn) - n + 1)]
    # nmers = {(seq, [seqn[i:i+n] for i in range(len(seqn)-n+1)]) for seq,seqn in izip(sequences,sequences_n)}

    return nmers

def scoreDNA(a, b, mask_score=None, gap_score=None):
    """
    Returns the score for a pair of IUPAC Ambiguous Nucleotide characters

    Arguments:
      a : First characters
      b : Second character
      n_score : Tuple of length two defining scores for all matches against an N
                character for (a, b), with the score for character (a) taking precedence;
                if None score symmetrically according to IUPAC character identity
      gap_score : Tuple of length two defining score for all matches against a gap (-, .)
                  character for (a, b), with the score for character (a) taking precedence;
                  if None score symmetrically according to IUPAC character identity

    Returns:
      int : Score for the character pair
    """
    # Define ambiguous character translations
    IUPAC_trans = {'AGWSKMBDHV': 'R', 'CTSWKMBDHV': 'Y', 'CGKMBDHV': 'S', 'ATKMBDHV': 'W', 'GTBDHV': 'K',
                   'ACBDHV': 'M', 'CGTDHV': 'B', 'AGTHV': 'D', 'ACTV': 'H', 'ACG': 'V', 'ABCDGHKMRSTVWY': 'N',
                   '-.': '.'}
    # Create list of tuples of synonymous character pairs
    IUPAC_matches = [p for k, v in IUPAC_trans.items() for p in list(product(k, v))]

    # Check gap and N-value conditions, prioritizing score for first character
    if gap_score is not None and a in '-.':
        return gap_score[0]
    elif mask_score is not None and a in 'nN':
        return mask_score[0]
    elif gap_score is not None and b in '-.':
        return gap_score[1]
    elif mask_score is not None and b in 'nN':
        return mask_score[1]

    # Return symmetric and reflexive score for IUPAC match conditions
    if a == b:
        return 1
    elif (a, b) in IUPAC_matches:
        return 1
    elif (b, a) in IUPAC_matches:
        return 1
    else:
        return 0

def getDNADistMatrix(mat=None, mask_dist=0, gap_dist=0):
    """
    Generates a DNA distance matrix

    Arguments:
      mat : Input distance matrix to extend to full alphabet;
            if unspecified, creates Hamming distance matrix that incorporates
            IUPAC equivalencies
      mask_dist : Distance for all matches against an N character
      gap_dist : Distance for all matches against a gap (-, .) character

    Returns:
      DataFrame : pandas.DataFrame of distances
    """
    IUPAC_chars = list('-.ACGTRYSWKMBDHVN')
    mask_char = 'N'

    # Default matrix to inf
    dist_mat = pd.DataFrame(float('inf'), index=IUPAC_chars, columns=IUPAC_chars,
                            dtype=float)
    # Set gap distance
    for c in '-.':
        dist_mat.loc[c] = dist_mat.loc[:, c] = gap_dist

    # Set mask distance
    dist_mat.loc[mask_char] = dist_mat.loc[:, mask_char] = mask_dist

    # Fill in provided distances from input matrix
    if mat is not None:
        for i, j in product(mat.index, mat.columns):
            dist_mat.at[i, j] = mat.at[i, j]
    # If no input matrix, create IUPAC-defined Hamming distance
    else:
        for i, j in product(dist_mat.index, dist_mat.columns):
            dist_mat.at[i, j] = 1 - scoreDNA(i, j,
                                             mask_score=(1 - mask_dist, 1 - mask_dist),
                                             gap_score=(1 - gap_dist, 1 - gap_dist))

    return dist_mat

    pass

def calcDistances(sequences, n, dist_mat, norm, sym):
    """
    Calculate pairwise distances between input sequences

    Arguments:
      sequences : List of sequences for which to calculate pairwise distances
      n : Length of n-mers to be used in calculating distance
      dist_mat : pandas.DataFrame of mutation distances
      norm : Normalization method
      sym : Symmetry method

    Returns:
      ndarray : numpy matrix of pairwise distances between input sequences
    """
    # Initialize output distance matrix
    dists = np.zeros((len(sequences), len(sequences)))
    # Generate dictionary of n-mers from input sequences
    nmers = getNmers(sequences, n)
    # Iterate over combinations of input sequences
    for j, k in combinations(list(range(len(sequences))), 2):
        # Only consider characters and n-mers with mutations
        mutated = [i for i, (c1, c2) in enumerate(zip(sequences[j], sequences[k])) if c1 != c2]
        seq1 = [sequences[j][i] for i in mutated]
        seq2 = [sequences[k][i] for i in mutated]
        nmer1 = [nmers[sequences[j]][i] for i in mutated]
        nmer2 = [nmers[sequences[k]][i] for i in mutated]

        # Determine normalizing factor
        if norm == 'len':
            norm_by = len(sequences[0])
        elif norm == 'mut':
            norm_by = len(mutated)
        else:
            norm_by = 1

        # Determine symmetry function
        if sym == 'avg':
            sym_fun = np.mean
        elif sym == 'min':
            sym_fun = min
        else:
            sym_fun = sum

        # Calculate distances
        try:
            dists[j, k] = dists[k, j] = \
                sum([sym_fun([dist_mat.at[c1, n2], dist_mat.at[c2, n1]]) \
                     for c1, c2, n1, n2 in zip(seq1, seq2, nmer1, nmer2)]) / \
                (norm_by)
        except (KeyError):
            raise KeyError('Unrecognized character in sequence.')
    return dists

def formClusters(dists, link, distance):
    """
    Form clusters based on hierarchical clustering of input distance matrix with
    linkage type and cutoff distance

    Arguments:
      dists : numpy matrix of distances
      link : Linkage type for hierarchical clustering
      distance : Distance at which to cut into clusters

    Returns:
      list : List of cluster assignments
    """
    # Make distance matrix square
    dists = squareform(dists)
    # Compute linkage
    links = linkage(dists, link)
    # Break into clusters based on cutoff
    clusters = fcluster(links, distance, criterion='distance')
    return clusters

def hier_clust(group, distance):
    """
    Form clusters based on hierarchical clustering of input distance matrix with
    linkage type and cutoff distance
    """
    dict_seqID = group.set_index('CDR3_MASK').to_dict()['SEQUENCE_ID']
    seqs = group['CDR3_MASK'].tolist()
    IDs = group['SEQUENCE_ID'].tolist()
    seqs_uniq = list(set(seqs))
    seq_map = {}
    for key, row in group.iterrows():
        seq = row['CDR3_MASK']
        ID = row['SEQUENCE_ID']
        seq_map.setdefault(seq, []).append(ID)

    if len(seqs_uniq) == 1:
        clone_tmp = [IDs[0] for i in range(len(IDs))]
    else:
        dist_mat = getDNADistMatrix(mask_dist=0, gap_dist=0)
        dists = calcDistances(seqs_uniq, 1, dist_mat, 'len', 'avg')
        # Perform hierarchical clustering
        lineage = 'single'
        clusters = formClusters(dists, lineage, distance)

        # Turn clusters into clone dictionary
        clone_dict = {}
        for i, c in enumerate(clusters):
            cdr3seq = seqs_uniq[i]
            for seq_id in seq_map[cdr3seq]:
                clone_dict[seq_id] = c
            # clone_dict.setdefault(c, []).extend(seq_map[seqs_uniq[i]])
        clone_tmp = ['%s_%d' % (IDs[0], clone_dict[seq_id]) for seq_id in IDs]
    return clone_tmp


def getGermdict(args):
    ''' Read VDJ IgBlast database and obtain un-gapped germline sequences
    	2020/09 Lawrence: Add in a check condition to ensure databases are read properly
    '''
    germ_dict = {}
    Vdb = getInputSeq(args.params_dict['Vdb'])
    Ddb = getInputSeq(args.params_dict['Ddb'])
    Jdb = getInputSeq(args.params_dict['Jdb'])
    if not bool(Vdb):
        print('FAILED TO LOAD %s' % args.params_dict['Vdb'], file = sys.stderr)
    else:
        germ_dict.update(Vdb)
    if not bool(Ddb):
        print('FAILED TO LOAD %s' % args.params_dict['Ddb'], file = sys.stderr)
    else:
        germ_dict.update(Ddb)
    if not bool(Jdb):
        print('FAILED TO LOAD %s' % args.params_dict['Jdb'], file = sys.stderr)
    else:
        germ_dict.update(Jdb)
    return germ_dict


def collapse_db(records, collapse_type, N_Diff):
    '''
    Collapse reads Db file

    Input:
        Records: read dataframe
        collaspe_type:
            'identical' -- collpase input sequences are identical
            'partial'  -- collapse shorter sequences to longer ones
            'V1'  -- collapse sequences by multiple columns
        N_Diff:
            'T' consider N as difference
            'F' consider N as not difference
    Output:
    '''

    def _writeOutputSeq(filter_check, col_ID, col_seq, output_file, db):
        """
        Generate output fasta file using db dataframe
        """
        output_handle = open(output_file, 'w')
        for key, row in db.iterrows():
            output_handle.write('>%s\n%s\n' % (row[col_ID], row[col_seq]))
        output_handle.close()

    def _collapse_identical_NasDiff(records):
        def __parse_group(group):
            index_dupreads = ','.join(group['SEQUENCE_ID'])
            ### 20200916 Lawrence: updated .ix to .loc
            top_series = group.loc[group.index[0]]
            top_series['DUPREAD'] = index_dupreads
            return top_series

        grouped = records.groupby('SEQUENCE_INPUT')
        colnames = list(records) + ['DUPREAD']
        records_collapse = pd.DataFrame(columns=colnames, index=range(0, len(grouped)))
        records_collapse = grouped.apply(__parse_group)

        return records_collapse

        # grouped = records.groupby('SEQUENCE_INPUT')
        # index_dupreads = {}
        # indexList = []
        # for key, group in grouped:
        #     idx = group.index[0]
        #     indexList.append(idx)
        #     index_dupreads[idx] = ','.join(group['SEQUENCE_ID'])
        # records_collapse = records.loc[indexList]
        # for idx in index_dupreads:
        #     records_collapse.ix[idx, 'DUPREAD'] = index_dupreads[idx]
        # return records_collapse

    # def _parse_read(row, records_collect):
    #     # Keep read with 'N'
    #     if 'N' in row['SEQUENCE_INPUT']:
    #         records_collect = records_collect.append(row)
    #         return records_collect
    #     else:
    #         records_cdr3 = records_collect[records_collect['CDR3_SEQ']==row['CDR3_SEQ']]
    #         for key,collect in records_cdr3.iterrows():
    #             if row['SEQUENCE_INPUT'] in collect['SEQUENCE_INPUT']:
    #                 records_collect.ix[key, 'DUPREAD'] += ',%s' % row['DUPREAD']
    #                 return records_collect
    #         records_collect = records_collect.append(row)
    #         return records_collect
    #
    # def _collapse_partial_NasDiff(records):
    #     colnames =  list(records) #+ ['N']
    #     records_collect = pd.DataFrame(columns=colnames)
    #     for key,row in records.iterrows():
    #         records_collect = _parse_read(row, records_collect)
    #     return records_collect
    def _collapse_partial_NasDiff(records):
        ''' Collapse shorter reads to longer ones
            Process: check read one by one to see if its input seq is a substring of stored reads
            Need a new method to speed up this
        '''
        records_collect = pd.DataFrame(columns=list(records))
        ### 20200916 Lawrence: updated .ix to .loc
        records_collect.loc[0,] = records.loc[records.index[0]]
        for key, row in records.iterrows():
            if key != records.index[0]:
                inputseq = row['SEQUENCE_INPUT']
                j = pd.Series(records_collect['SEQUENCE_INPUT']).str.contains(inputseq)
                if len(j[j == True]) >= 1:
                    i = j[j == True].index[0]
                    records_collect.loc[i, 'DUPREAD'] += ',%s' % row['DUPREAD']
                elif len(j[j == True]) == 0:
                    records_collect.loc[len(records_collect) + 1,] = row
        return records_collect

    def _parse_SAM(read_readlen, sam_file, collapse_type):
        inputR_refR = {}
        for line in open(sam_file):
            l = line.strip().split()
            if l[5] == '%dM' % read_readlen[l[0]]:
                if collapse_type == 'identical':
                    if l[5] == '%dM' % read_readlen[l[2]]:
                        inputR_refR[l[0]] = l[2]
                else:
                    inputR_refR[l[0]] = l[2]
        return inputR_refR

    def _collapse_V1(records):
        ''' Collapse reads based on various result columns
        '''
        records_new = pd.DataFrame(columns=list(records))
        grouplist = ['V_ALLELE', 'D_ALLELE', 'J_ALLELE', 'STOP', 'IN_FRAME', \
                     'V_END', 'V_D_JUNCTION', 'D_REGION', 'D_J_JUNCTION', \
                     'J_START', 'V_J_JUNCTION']
        for key, group in records.groupby(grouplist):
            dup = ','.join(group['DUPREAD'].values.tolist())
            groupcontent = group.iloc[0]
            groupcontent['DUPREAD'] = dup
            records_new.loc[len(records_new) + 1,] = groupcontent
        return records_new

    def _collapse_NasNoDiff(records, collapse_type):
        ''' Required Bowtie2 software
        '''
        randname = str(random.randint(1, 1000000))

        # Write with/wo 'N' two fa files as input/ref files in bowtie2 searching
        records_woN = records[~records['SEQUENCE_INPUT'].str.contains("N")]
        records_wN = records[records['SEQUENCE_INPUT'].str.contains("N")]
        if len(records_woN) == 0 or len(records_wN) == 0:
            return records
        ref_file = '%s.ref' % randname
        input_file = '%s.input' % randname
        _writeOutputSeq('woN', 'SEQUENCE_ID', 'SEQUENCE_INPUT', ref_file, records_woN)
        _writeOutputSeq('wN', 'SEQUENCE_ID', 'SEQUENCE_INPUT', input_file, records_wN)

        sam_file = '%s.sam' % randname
        os.system('bowtie2-build %s %s -q' % (ref_file, ref_file))
        os.system('bowtie2 -x ./%s -f -U %s --local -S %s --no-head  --np 0 --mp 1000 --rdg 1000,1000 --rfg 1000,1000' % \
                  (ref_file, input_file, sam_file))

        read_readlen = records.set_index('SEQUENCE_ID').to_dict()['INPUT_LEN']
        inputR_refR = _parse_SAM(read_readlen, sam_file, collapse_type)

        records_collapsed = records[~records.SEQUENCE_ID.isin(inputR_refR.keys())].copy()
        records_waitToCollapse = records[records.SEQUENCE_ID.isin(inputR_refR.keys())]
        for inputR in inputR_refR:
            refR = inputR_refR[inputR]
            dup = records_waitToCollapse.loc[records_waitToCollapse['SEQUENCE_ID'] == inputR, 'DUPREAD'].values[0]
            records_collapsed.loc[records_collapsed['SEQUENCE_ID'] == refR, 'DUPREAD'] += ',%s' % dup
        os.system('rm -rf %s %s %s %s*bt2' % (ref_file, input_file, sam_file, ref_file))
        return records_collapsed

    # Main part in this func
    # Collapse identical reads anyway
    records = _collapse_identical_NasDiff(records)
    records['INPUT_LEN'] = records["SEQUENCE_INPUT"].map(len)
    records.sort_values('INPUT_LEN', ascending=False, inplace=True)
    # Collapse identical reads with N as no difference
    if collapse_type == 'identical' and N_Diff == 'F':
        records = _collapse_NasNoDiff(records, 'identical')
    elif collapse_type == 'partial':
        # Collapse shorter reads to longer ones with N as difference
        records = _collapse_partial_NasDiff(records)
        if N_Diff == 'F':
            # Collapse shorter reads to longer ones with N as no difference
            records = _collapse_NasNoDiff(records, 'partial')
            records = records.drop('INPUT_LEN', axis=1)
    elif collapse_type == 'V1':
        # V1 means same way as Version One pipeline
        records = _collapse_NasNoDiff(records, 'identical')
        records = _collapse_V1(records)

    return records

def profile_DNAmut(group, nuc_stat, nuc_PDF, nuc_profile, args):
    ''' Prep DNA mutation profile, text and PDF file
    '''
    
    def _parse_V_ALLELE_NUC(row):
        return pd.Series([row["SEQUENCE_ID"]] + [s for s in row['V_ALLELE_NUC']])

    allele = group["V_ALLELE"].unique()[0]
    germ_dict = getGermdict(args)
    allele_seq = germ_dict[allele]
    allele_len = len(allele_seq)
    colnames = ['ID'] + [l for l in allele_seq]
    allele_mut = pd.DataFrame(columns=colnames, index=range(0, len(group)))
    allele_mut = group.apply(_parse_V_ALLELE_NUC, axis=1)

    statnames = ['Pos', 'Mut', 'Total', 'Base', 'Y', 'A', 'T', 'C', 'G']
    allele_stat = pd.DataFrame(columns=statnames, index=range(1, allele_len + 1))
    allele_stat['Pos'] = range(1, allele_len + 1)
    allele_stat['Base'] = [l for l in allele_seq]
    allele_stat[['Mut', 'Total', 'Y', 'A', 'T', 'C', 'G']] = 0
    for i in range(1, allele_len + 1):
        if len(allele_mut) == 1:
            counts = {}
            counts[allele_mut[[i]].squeeze()] = 1
        else:
            counts = allele_mut[[i]].squeeze().value_counts()
        countA = counts.get('A', 0)
        countT = counts.get('T', 0)
        countC = counts.get('C', 0)
        countG = counts.get('G', 0)
        countMut = countA + countT + countC + countG
        countTotal = countMut + counts.get('.', 0)
        allele_stat.loc[i, 'Mut'] = countMut
        allele_stat.loc[i, 'Total'] = countTotal
        allele_stat.loc[i, 'Y'] = float(countMut / countTotal) if countTotal > 0 else 0
        allele_stat.loc[i, 'A'] = countA
        allele_stat.loc[i, 'T'] = countT
        allele_stat.loc[i, 'C'] = countC
        allele_stat.loc[i, 'G'] = countG

    allele_mut.to_csv(nuc_profile, sep="\t", index=False)
    allele_stat.to_csv(nuc_stat, sep="\t", index=False)

    # run R scripts
    if allele in args.__dict__['V_CDR']:
        cdr = args.__dict__['V_CDR'][allele]
        cdrstring = 'cdr1_start=%s cdr1_end=%s cdr2_start=%s cdr2_end=%s ' \
                    'cdr3_start=%s cdr3_end=%s' % (cdr[0], cdr[1], cdr[2],
                                                   cdr[3], cdr[4], cdr[5])
    else:
        cdrstring = ''

    # Filter group with read number
    if len(group) >= args.min_profileread:
        sample = group["SAMPLE"].unique()[0]
        anno = '_'.join([sample, allele])
    ### 20200915 Lawrence: changed showsequence from false to true
    ### 20200916 Lawrence: changed frpm 'Rscript %s/HTGTSrep/R/SHMPlot2.R %s %s plotrows=1 figureheight=2 showsequence=TRUE ymax=%f %s annotation=%s '
    ### to 'Rscript %s/HTGTSrep/R/SHMPlot2.R \"%s\" \"%s\" plotrows=1 figureheight=2 showsequence=TRUE ymax=%f %s annotation=\"%s\" '
        ### this allows special characters to be in V_allel names
        os.system('Rscript %s/HTGTSrep/R/SHMPlot2.R \"%s\" \"%s\" plotrows=1 figureheight=2 '
                  'showsequence=TRUE ymax=%f %s annotation=\"%s\" ' % (args.scriptdir, nuc_stat,
                                                                    nuc_PDF, args.ymax_DNA, cdrstring, anno))


def getInferSeq(treefile, group):
    ''' Read tree file and get inferred 1 sequence,
        Using germline sequence of V and J to ensure no surprise
    '''
    n = 0
    with open(treefile) as f:
        for line in f:
            n += 1
            l = line.strip().split()
            if n == 2: inferseq = l[-1]
    # get germline parts
    cdr3 = group.iloc[0]['CDR3_MASK']
    sequence_imgt = group.iloc[0]['SEQUENCE_IMGT']
    Vpos = sequence_imgt.find(cdr3)

    germline_imgt_seq = group.iloc[0]['GERMLINE_IMGT_D_MASK']
    seq_V = germline_imgt_seq[:Vpos]
    seq_Junc = inferseq[Vpos:Vpos + len(cdr3)]
    if len(germline_imgt_seq) >= len(inferseq):
        seq_J = germline_imgt_seq[Vpos + len(cdr3): len(inferseq)]
    else:
        seq_J = inferseq[Vpos + len(cdr3): len(inferseq)]

    # Use V and J parts from germline as reference to avoid mismatch at these regions in mut profiling
    newinfer = (seq_V + seq_Junc + seq_J).replace('.', 'N')
    # print(group['CLONE'].tolist()[0], inferseq, newinfer)
    return newinfer

def profile_DNAmut_clonal(inferseq, group, nuc_stat, nuc_PDF, nuc_profile, args):
    ''' Prep DNA mutation profile, text and PDF file
    '''
    allele = group["V_CALL"].unique()[0]

    # Get position list of inferred seq which are not 'N'
    poslist = [i for i in range(0, len(inferseq)) if inferseq[i] != 'N']
    allele_seq = ''.join([inferseq[i] for i in poslist])
    allele_len = len(poslist)

    colnames = ['ID'] + [inferseq[i] for i in poslist]
    allele_mut = pd.DataFrame(columns=colnames)
    for key, row in group.iterrows():
        seq = row['SEQUENCE_IMGT']
        vals = [row["SEQUENCE_ID"]]
        for i in poslist:
            if i >= len(seq):
                vals.append('-')
            else:
                l = seq[i]
                if l == '.': vals.append('-')
                if l == '-': vals.append('N')
                if l in 'ATCGN':
                    if l == inferseq[i]:
                        vals.append('.')
                    else:
                        vals.append(l)
        allele_mut.loc[len(allele_mut) + 1] = vals

    statnames = ['Pos', 'Mut', 'Total', 'Base', 'Y', 'A', 'T', 'C', 'G']
    allele_stat = pd.DataFrame(columns=statnames, index=range(1, allele_len + 1))
    allele_stat['Pos'] = range(1, allele_len + 1)
    allele_stat['Base'] = [l for l in allele_seq]
    allele_stat[['Mut', 'Total', 'Y', 'A', 'T', 'C', 'G']] = 0
    for i in range(1, allele_len + 1):
        if len(allele_mut) == 1:
            counts = {}
            ### 20200916 Lawrence: updated .ix to .iloc
            counts[allele_mut.iloc[:, i].squeeze()] = 1
        else:
            ### 20200916 Lawrence: updated .ix to .iloc
            counts = allele_mut.iloc[:, i].squeeze().value_counts()
        countA = counts.get('A', 0)
        countT = counts.get('T', 0)
        countC = counts.get('C', 0)
        countG = counts.get('G', 0)
        countMut = countA + countT + countC + countG
        countTotal = countMut + counts.get('.', 0)
        allele_stat.loc[i, 'Mut'] = countMut
        allele_stat.loc[i, 'Total'] = countTotal
        allele_stat.loc[i, 'Y'] = float(countMut / countTotal) if countTotal > 0 else 0
        allele_stat.loc[i, 'A'] = countA
        allele_stat.loc[i, 'T'] = countT
        allele_stat.loc[i, 'C'] = countC
        allele_stat.loc[i, 'G'] = countG

    allele_mut.to_csv(nuc_profile, sep="\t", index=False)
    allele_stat.to_csv(nuc_stat, sep="\t", index=False)

    # Obtain CDR1,2 from preprocessed V allele alignment using KABAT definition, and CDR3 from IMGT definition
    cdr3_start = len(inferseq[0:312].replace('N', '')) + 1
    cdr3_end = cdr3_start + group["JUNCTION_LENGTH"].unique()[0]
    if allele in args.__dict__['V_CDR']:
        cdr = args.__dict__['V_CDR'][allele]
        cdrstring = 'cdr1_start=%s cdr1_end=%s cdr2_start=%s cdr2_end=%s ' \
                    'cdr3_start=%d cdr3_end=%d' % (cdr[0], cdr[1], cdr[2],
                                                   cdr[3], cdr3_start, cdr3_end)
    else:
        cdrstring = ''

    CLONE = group["CLONE"].unique()[0]
    J_ALLELE = group["J_CALL"].unique()[0]
    anno = '_'.join(['Clone%d' % CLONE, allele, J_ALLELE])

    # Filter group with read number
    ### 20200915 Lawrence: changed showsequence from false to true
    os.system('Rscript %s/HTGTSrep/R/SHMPlot2.R %s %s plotrows=1 figureheight=2 '
              'showsequence=TRUE ymax=%f %s annotation=%s' % (args.scriptdir, nuc_stat,
                                                               nuc_PDF, args.ymax_DNA, cdrstring, anno))


def profile_DNAmut_clonal_errbar(inferseq, group, nuc_stat_errbar, nuc_PDF_errbar, sample_files, args):
    ''' Prep DNA mutation profile, text and PDF file with error bars
    Input:
        inferseq: inferred sequence as reference
        group: group instance
        nuc_stat_errbar: stat file with sem value
        nuc_PDF_errbar: PDF file with error bars
        samples_files: list of sample files as sample_files[sample] = [root_PDF, root_stat]
    '''
    allele = group["V_CALL"].unique()[0]
    if len(sample_files) > 1:
        stat_list = []
        for sample in sample_files:
            stat = pd.read_csv(sample_files[sample][1], sep="\t")
            stat_list.append(stat)
        stat_all = pd.concat(stat_list, ignore_index=True)

        pos_max = stat_all['Pos'].max()
        statcols = list(stat_all) + ['Err']
        stat_new = pd.DataFrame(columns=statcols, dtype='int64')
        stat_new = stat_new.astype('int64')
        stat_new['Base'] = stat_new['Base'].astype('str')
        stat_new['Y'] = stat_new['Y'].astype('float')
        for i in range(1, pos_max + 1):
            cals = stat_all.loc[stat_all['Pos'] == i]
            Pos = i
            Mut = cals['Mut'].sum()
            Total = cals['Total'].sum()
            Base = list(cals['Base'])[0]
            A = cals['A'].sum()
            T = cals['T'].sum()
            C = cals['C'].sum()
            G = cals['G'].sum()

            cals_Y = cals[cals['Total'] > 0]
            Err = 0
            if len(cals_Y) == 0:
                Y = 0
            else:
                Y = cals_Y['Y'].sum() / len(cals_Y)
            if len(cals_Y) > 1:
                Err = sem(cals_Y['Y'])
            stat_new.loc[i] = [int(Pos), int(Mut), int(Total), str(Base), Y,
                               int(A), int(T), int(C), int(G), Err]
        stat_new.to_csv(nuc_stat_errbar, sep="\t", index=False)

        # Obtain CDR1,2 from preprocessed V allele alignment using KABAT definition, and CDR3 from IMGT definition
        cdr3_start = len(inferseq[0:312].replace('N', '')) + 1
        cdr3_end = cdr3_start + group["JUNCTION_LENGTH"].unique()[0]
        # add this part as some allele names with [] anno
        if allele.split('[')[0] in args.__dict__['V_CDR']:
            cdr = args.__dict__['V_CDR'][allele.split('[')[0]]
            cdrstring = 'cdr1_start=%s cdr1_end=%s cdr2_start=%s cdr2_end=%s ' \
                        'cdr3_start=%d cdr3_end=%d' % (cdr[0], cdr[1], cdr[2],
                                                       cdr[3], cdr3_start, cdr3_end)
        else:
            cdrstring = ''

        CLONE = group["CLONE"].unique()[0]
        J_ALLELE = group["J_CALL"].unique()[0]
        anno = '_'.join(['Clone%d' % CLONE, allele, J_ALLELE])
        # Filter group with read number
        ### 20200915 Lawrence: changed showsequence from false to true
        os.system('Rscript %s/HTGTSrep/R/SHMPlot2.R %s %s plotrows=1 figureheight=2 '
                  'showsequence=TRUE ymax=0.75 %s annotation=%s' % (args.scriptdir,
                                                                     nuc_stat_errbar, nuc_PDF_errbar, cdrstring, anno))
