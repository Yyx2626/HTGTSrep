#!/usr/bin/env python3
import os, re, sys, logging, csv, multiprocessing
import pandas as pd
from itertools import groupby
import itertools, functools
try:
    from Bio.Alphabet import generic_dna, IUPAC
    Bio_Alphabet = True
except ImportError:
    Bio_Alphabet = None
    # usages of generic_dna, IUPAC are not supported in Biopython 1.78 (September 2020).
    print(f"The installed BioPython is a new version that has removed the Alphabet module.",file=sys.stderr)
from Bio.Seq import Seq
from collections import OrderedDict
from HTGTSrep.lib import loggingRun, getInputSeq, getCSV, collapse_db

def mutV(ig_dict):
    ''' Parse V_BTOP string and generate converted V string
    '''
    # parse V_BTOP, do it in a naive way
    V_BTOP = ig_dict['V_BTOP']
    V_mutNum = 0
    V_string = '-' * (ig_dict['V_GERM_START_VDJ'] - 1)
    match_breaks = []
    for m in re.finditer(r'\d+', V_BTOP):
        match_breaks.append(m.start(0))
        match_breaks.append(m.end(0))
    match_breaks.append(len(V_BTOP))

    for i in range(0, len(match_breaks)):
        if i == len(match_breaks)-1: continue
        match = V_BTOP[match_breaks[i]:match_breaks[i+1]]
        if match:
            try:
                match = int(match)
                V_string += '.' * match
            except:
                mismatches = [match[i: i+2] for i in range(0, len(match), 2)]
                for mismatch in mismatches:
                    if mismatch[1] == '-': continue
                    if mismatch[0] == '-':
                        V_string += '-'
                    else:
                        V_string += mismatch[0]
                        if mismatch[0] in 'ATCG': V_mutNum += 1
    if len(V_string) > ig_dict['V_GENE_LEN']:
        logging.error('Error: V mutation string longer than V gene, %s' \
                        % ig_dict['SEQUENCE_ID'])
        sys.exit()
    else:
        V_string += '-' * (ig_dict['V_GENE_LEN'] - len(V_string))
    ig_dict['V_ALLELE_NUC'] = V_string
    ig_dict['V_MUTATION'] = V_mutNum
    return ig_dict

def dmask(ig_dict, repo_dict):
    """
    Join gapped germline sequences aligned with sample sequences

    Arguments:
    ig_dict = iterable yielding dictionaries of sample sequence data
    repo_dict = dictionary of IMGT gapped germline sequences

    Returns:
    dictionary of germline_type: germline_sequence
    """
    align = ig_dict
    vgene = ig_dict['V_ALLELE']
    if vgene in repo_dict:
        vseq = repo_dict[vgene]
        ig_dict['V_GENE_GAP_LEN'] = len(vseq)
        vstart = int(float(align['V_GERM_START_IMGT'])) - 1
        vlen = int(float(align['V_GERM_LENGTH_IMGT']))
        vpad = vlen - len(vseq[vstart:])
        if vpad < 0: vpad = 0
        germ_vseq = vseq[vstart:(vstart + vlen)] + ('N' * vpad)
    else:
        return ig_dict

    dgene = ig_dict['D_ALLELE']
    if dgene != '-':
        if dgene in repo_dict:
            dseq = repo_dict[dgene]
            # Germline start
            try: dstart = int(float(align['D_GERM_START'])) - 1
            except (TypeError, ValueError): dstart = 0
            # Germline length
            try: dlen = int(float(align['D_GERM_LENGTH']))
            except (TypeError, ValueError): dlen = 0
            germ_dseq = repo_dict[dgene][dstart:(dstart + dlen)]
        else:
            # logging.warning('Error: D gene %s not in the IMGT-gapped seq file' % dgene)
            return ig_dict
    else:
        germ_dseq = ''

    jgene = ig_dict['J_ALLELE']
    if jgene != '-':
        if jgene in repo_dict:
            jseq = repo_dict[jgene]
            # Germline start
            try: jstart = int(float(align['J_GERM_START'])) - 1
            except (TypeError, ValueError): jstart = 0
            # Germline length
            try: jlen = int(float(align['J_GERM_LENGTH']))
            except (TypeError, ValueError): jlen = 0
            jpad = jlen - len(jseq[jstart:])
            if jpad < 0: jpad = 0
            germ_jseq = jseq[jstart:(jstart + jlen)] + ('N' * jpad)
        else:
            # logging.warning('Error: J gene %s not in the IMGT-gapped seq file' % jgene)
            return ig_dict
    else:
        try: jlen = int(float(align['J_GERM_LENGTH']))
        except (TypeError, ValueError): jlen = 0
        germ_jseq = 'N' * jlen

    # Assemble pieces starting with V-region
    germ_seq = germ_vseq
    regions = 'V' * len(germ_vseq)

    try:
        np1_len = int(float(align['NP1_LENGTH']))
    except (TypeError, ValueError):
        np1_len = 0

    # PNP nucleotide additions after V
    if np1_len < 0:
        result_log['ERROR'] = 'NP1_LENGTH is negative'
        return result_log, germlines

    regions += 'N' * np1_len

    germ_seq += 'N' * np1_len

    # Add D-region
    germ_seq += germ_dseq
    regions += 'D' * len(germ_dseq)

    # 'VD>', germ_seq, '\nVD>', regions

    try:
        np2_len = int(float(align['NP2_LENGTH']))
    except (TypeError, ValueError):
        np2_len = 0

    # NP nucleotide additions
    if np2_len < 0:
        result_log['ERROR'] = 'NP2_LENGTH is negative'
        return result_log, germlines

    regions += 'N' * np2_len
    germ_seq += 'N' * np2_len

    # Add J-region
    germ_seq += germ_jseq
    regions += 'J' * len(germ_jseq)

    # Full length and regions of germline seq, might be useful in the future
    germlines_full = germ_seq
    germlines_regions = regions

    seq_dmask = germ_seq[:len(germ_vseq)] + \
                'N' * (len(germ_seq) - len(germ_vseq) - len(germ_jseq)) + \
                 germ_seq[-len(germ_jseq):]
    ig_dict['GERMLINE_IMGT_D_MASK'] = seq_dmask

    return ig_dict

def gapV(ig_dict, repo_dict):
    """
    Insert gaps into V region and update alignment information

    Arguments:
      ig_dict : Dictionary of parsed IgBlast output
      repo_dict : Dictionary of IMGT gapped germline sequences

    Returns:
      dict : Updated with SEQUENCE_IMGT, V_GERM_START_IMGT, and V_GERM_LENGTH_IMGT fields
    """

    seq_imgt = '.' * (int(ig_dict['V_GERM_START_VDJ']) - 1) + ig_dict['SEQUENCE_VDJ']
    #print("seq_imgt before gapping in gapV()", seq_imgt, file = sys.stderr)
    # present before gapping
    # Find gapped germline V segment
    vkey = ig_dict['V_ALLELE']
    #print("keys in repo_dict:", repo_dict.keys(), file=sys.stderr)
    #print("vkey of interest:", vkey, file=sys.stderr)
    if vkey in repo_dict:
        vgap = repo_dict[vkey]
        #print("vgap present?", vgap, file=sys.stderr)
        # Iterate over gaps in the germline segment
        gaps = re.finditer(r'\.', vgap)
        #print("vgap present?", gaps, file=sys.stderr)
        gapcount = int(ig_dict['V_GERM_START_VDJ']) - 1
        #print("gapcount present?", gapcount, file=sys.stderr)
        for gap in gaps:
            i = gap.start()
            # Break if gap begins after V region
            if i >= ig_dict['V_GERM_LENGTH_VDJ'] + gapcount:
                break
            # Insert gap into IMGT sequence
            seq_imgt = seq_imgt[:i] + '.' + seq_imgt[i:]
            # Update gap counter
            gapcount += 1
        #print("seq_imgt after gapping", seq_imgt, file=sys.stderr)
        ig_dict['SEQUENCE_IMGT'] = seq_imgt
        # Update IMGT positioning information for V
        ig_dict['V_GERM_START_IMGT'] = 1
        ig_dict['V_GERM_LENGTH_IMGT'] = ig_dict['V_GERM_LENGTH_VDJ'] + gapcount

    # print("seq_imgt keys() check after gapping in Vgap()", ig_dict.keys(), file=sys.stderr)
    return ig_dict

def readOneIgBlastResult(block):
    """
    Parse a single IgBlast query result

    Arguments:
    block =  itertools groupby object of single result

    Returns:
    None if no results, otherwise list of DataFrames for each result block
    """
    results = list()
    i = 0
    for match, subblock in groupby(block, lambda l: l=='\n'):
        if not match:
            # Strip whitespace and comments, and clonal parts
            sub = [s.strip() for s in subblock if not s.startswith((
                    '#', 'Total queries', 'Total identifiable CDR3',
                    'Total unique clonotypes'))]

            # Continue on empty block
            if not sub:  continue
            else:  i += 1

            # Split by tabs
            sub = [s.split('\t') for s in sub]
            # Append list for "V-(D)-J rearrangement summary" (i == 1)
            # And "V-(D)-J junction details" (i == 2)
            # Otherwise append DataFrame of subblock
            if i == 1 or i == 2:
                results.append(sub[0])
            # skip query result if no CDR3 info
            elif i == 3 and sub[0][0] == 'CDR3':
                    results.append(sub[0])
            else:
                df = pd.DataFrame(sub)
                if not df.empty: results.append(df)
    return results if results else None

def writeIgBlastDb(IgBlast_output, IgBlast_db, seq_file, repo_dict, subset, Primer_Jlen, Dupstream=None):
    '''
    Main for IgBlast aligned sample sequences, a parsed file in succinct fashion
    will be generated for each output file.

    Arguments:
    IgBlast_output = IgBlast output file to process
    IgBlast_db = IgBlast database file to write
    seq_file = Query sequence file to read
    repo_dict = dictionary of IMGT gapped germline sequences

    Returns:
    None
    '''
    # Define ordered outputs
    ordered_fields = ['SEQUENCE_ID',
    "V_GENE", "V_ALLELE", "D_ALLELE", "J_ALLELE",
    "STOP", "IN_FRAME", "PRODUCTIVE", "CHAIN_TYPE", "STRAND",
    "V_END", "V_D_JUNCTION", "D_REGION",
    "D_J_JUNCTION", "J_START", "V_J_JUNCTION",
    "V_SCORE", "V_ALIGNMENT", "V_MISMATCH", "V_MUTATION",
    "INDEL", "V_IDENTITY", "V_COVERAGE",
    "V_BTOP", "V_GENE_LEN", "V_GENE_GAP_LEN",
    "D_SCORE", "D_ALIGNMENT",
    "J_SCORE", "J_ALIGNMENT", "J_MUTATION", "J_MUTATION_NOPRIMER",
    "CDR1_SEQ", "CDR1_PEPTIDE",
    "CDR2_SEQ", "CDR2_PEPTIDE",
    "CDR3_SEQ", "CDR3_PEPTIDE",
    "V_CALL", "D_CALL", "J_CALL",
    'D_SEQ_START', 'D_SEQ_LENGTH', 'D_GERM_START', 'D_GERM_LENGTH',
    'J_SEQ_START', 'J_SEQ_LENGTH', 'J_GERM_START', 'J_GERM_LENGTH',
    "NP1_LENGTH", "NP2_LENGTH",
    "V_SEQ_START", "V_SEQ_LENGTH",
    "V_GERM_START_VDJ", "V_GERM_END_VDJ",
    "V_GERM_START_IMGT", "V_GERM_LENGTH_IMGT",
    "V_ALLELE_NUC", "GERMLINE_IMGT_D_MASK",
    "SEQUENCE_INPUT", "SEQUENCE_VDJ", "SEQUENCE_IMGT"]

    if Dupstream:
        ordered_fields.append('D_UPSTREAM_MATCH_D_GENE')
        ordered_fields.append('D_UPSTREAM_STITCH_D_GENE')
    # Initialize IgBlast db
    db_handle = open(IgBlast_db, 'wt')
    db_writer = csv.writer(db_handle, delimiter = "\t")
    # Yyx changed back 2021-06-02, bottleneck is not output
    db_writer.writerow(ordered_fields)
#    print('\t'.join(ordered_fields), file=db_handle)
    # Get input sequence dictionary
    seq_dict = getInputSeq(seq_file)
    # Yyx changed below on 2021-06-02
    def parseOneBlock(block):
#            block = list(block)
#            if k1: continue

            # Initialize db_gen
            db_gen = {}
            for item in ordered_fields:
                db_gen[item] = "-"
            # Extract sequence ID
            query_name = ' '.join(block[0].strip().split(' ')[2:])

            # Parse db_gen to have ID and input sequence
            db_gen['SEQUENCE_ID'] = query_name
            # Yyx add 2021-04-15, should check strand, reverse complement if query is minus strand
            should_reverse_complement = any(('your query represents the minus strand' in x) for x in block)
            # Parse further sub-blocks
            block_list = readOneIgBlastResult(block)
            # Skip read without alignment or V alignment
            if block_list is None: return
            if block_list[0][0] == 'N/A': return
            # Penultimate has to be dataframe of FR&CDR position
            if isinstance(block_list[-2], list): return
            # Parse quality information
            db_gen['STRAND'] = block_list[0][-1]
            db_gen['SEQUENCE_INPUT'] = seq_dict[query_name]
            if db_gen['STRAND'] == '-':
                if Bio_Alphabet:
                    db_gen['SEQUENCE_INPUT'] = str(Seq(db_gen['SEQUENCE_INPUT'],
                                    IUPAC.ambiguous_dna).reverse_complement())
                else:
                    db_gen['SEQUENCE_INPUT'] = str(Seq(db_gen['SEQUENCE_INPUT']).reverse_complement())
            if block_list[0][-2] == 'Yes': db_gen['PRODUCTIVE'] = 'T'
            if block_list[0][-2] == 'No': db_gen['PRODUCTIVE'] = 'F'
            if block_list[0][-3] == 'In-frame': db_gen['IN_FRAME'] = 'T'
            if block_list[0][-3] == 'Out-of-frame': db_gen['IN_FRAME'] = 'F'
            if block_list[0][-4] == 'Yes': db_gen['STOP'] = 'T'
            if block_list[0][-4] == 'No': db_gen['STOP'] = 'F'
            db_gen['CHAIN_TYPE'] = block_list[0][-5]
            # Parse J call
            if block_list[0][-6] != 'N/A': db_gen['J_CALL'] = block_list[0][-6]
            db_gen['J_ALLELE'] = db_gen['J_CALL'].split(',')[0]
            # Parse D call
            if block_list[0][3] == 'VH':
                if block_list[0][1] != 'N/A': db_gen['D_CALL'] = block_list[0][1]
                db_gen['D_ALLELE'] = db_gen['D_CALL'].split(',')[0]
            # Parse V call
            if block_list[0][0] != 'N/A': db_gen['V_CALL'] = block_list[0][0]
            db_gen['V_ALLELE'] = db_gen['V_CALL'].split(',')[0]
            db_gen['V_GENE'] = db_gen['V_ALLELE'].split('*')[0]

            # Parse junction sequence
            if len(block_list[1]) >= 5:
            # ALTERNATIVELY FILTER BY CHAIN_TYPE
            # if db_gen['CHAIN_TYPE'] == 'VH':
                if block_list[1][0] != 'N/A': db_gen['V_END'] = block_list[1][0]
                if block_list[1][1] != 'N/A': db_gen['V_D_JUNCTION'] = block_list[1][1]
                if block_list[1][2] != 'N/A': db_gen['D_REGION'] = block_list[1][2]
                if block_list[1][3] != 'N/A': db_gen['D_J_JUNCTION'] = block_list[1][3]
                if block_list[1][4] != 'N/A': db_gen['J_START'] = block_list[1][4]
            # ALTERNATIVELY FILTER BY CHAIN_TYPE
            # elif db_gen['CHAIN_TYPE'] == 'VK':
            elif len(block_list[1]) == 3:
                if block_list[1][0] != 'N/A': db_gen['V_END'] = block_list[1][0]
                if block_list[1][1] != 'N/A': db_gen['V_J_JUNCTION'] = block_list[1][1]
                if block_list[1][2] != 'N/A': db_gen['J_START'] = block_list[1][2]
            # Parse CDR 1 & 2 sequence, should get reverse complement before extract
            # if subset == 'unjoinR2':
            #     inputseq = seq_dict[query_name]
            # else:

            if Bio_Alphabet:
                inputseq = Seq(seq_dict[query_name], IUPAC.ambiguous_dna)
            else:
                inputseq = Seq(seq_dict[query_name])
            # Yyx add 2021-04-15, should check strand, reverse complement if query is minus strand
            if should_reverse_complement:
                inputseq = inputseq.reverse_complement()
            inputseq = str(inputseq)
            # CDR 1 & 2 block is dataframe
            hit_regions = block_list[-2]
            regions = list(hit_regions.loc[:,0])
            # Parse CDR1/2 only when exist FR1/2
            FR1_EXIST = 0
            FR2_EXIST = 0
            if 'FR1-IMGT' in regions or 'FR1' in regions: FR1_EXIST = 1
            if 'FR2-IMGT' in regions or 'FR2' in regions: FR2_EXIST = 1
            for key, row in hit_regions.iterrows():
                if row[0].startswith('CDR1') and FR1_EXIST == 1:
                    db_gen['CDR1_SEQ'] = inputseq[int(row[1])-1: int(row[2])]
                    if Bio_Alphabet:
                        db_gen['CDR1_PEPTIDE'] = str(Seq(db_gen['CDR1_SEQ'],
                                                generic_dna).translate())
                    else:
                        db_gen['CDR1_PEPTIDE'] = str(Seq(db_gen['CDR1_SEQ']).translate())
                if row[0].startswith('CDR2') and FR2_EXIST == 1:
                    db_gen['CDR2_SEQ'] = inputseq[int(row[1])-1: int(row[2])]
                    if Bio_Alphabet:
                        db_gen['CDR2_PEPTIDE'] = str(Seq(db_gen['CDR2_SEQ'],
                                                generic_dna).translate())
                    else:
                        db_gen['CDR2_PEPTIDE'] = str(Seq(db_gen['CDR2_SEQ']).translate())

            # CDR 1 & 2 block is list
            if isinstance(block_list[2], list) and len(block_list[2]) > 2:
                db_gen['CDR3_SEQ'] = block_list[2][1]
                db_gen['CDR3_PEPTIDE'] = block_list[2][2]
            # Parse segment start and stop positions
            hit_df = block_list[-1]
            seq_vdj = ''

            v_align = hit_df[hit_df[0] == 'V'].iloc[0]
            # Alignment length and mismatch
            db_gen['V_IDENTITY']  = '%.3f' % (float(v_align[3]) / 100.0)
            db_gen['V_ALIGNMENT'] = v_align[4]
            db_gen['V_MISMATCH']  = v_align[5]
            db_gen['INDEL']       = v_align[6]
            db_gen['V_SEQ_START'] = int(v_align[8])
            db_gen['V_SEQ_LENGTH'] = int(v_align[9]) - db_gen['V_SEQ_START'] + 1
            db_gen['V_GERM_START_VDJ']= int(v_align[10])
            db_gen['V_GERM_END_VDJ']= int(v_align[11])
            db_gen['V_GERM_LENGTH_VDJ'] = int(v_align[11]) - db_gen['V_GERM_START_VDJ'] + 1
            db_gen['V_SCORE']     = v_align[13]
            db_gen['V_GENE_LEN']  = int(v_align[15])
            db_gen['V_COVERAGE']  = '%.3f' % (float(v_align[4])/float(v_align[15]))
            db_gen['V_BTOP']      = v_align[16]

            # Update VDJ sequence, removing insertions
            start = 0
            for m in re.finditer(r'-', v_align[18]):
                ins = m.start()
                seq_vdj += v_align[17][start:ins]
                start = ins + 1
            seq_vdj += v_align[17][start:]

            # D alignment results
            if db_gen['D_CALL'] != "-":
                d_align = hit_df[hit_df[0] == 'D'].iloc[0]
                db_gen['D_ALIGNMENT'] = d_align[4]
                db_gen['D_SCORE'] = d_align[13]

                # Determine N-region length and amount of J overlap with V or D alignment
                overlap = 0
                np1_len = int(d_align[8]) - (db_gen['V_SEQ_START'] + db_gen['V_SEQ_LENGTH'])
                if np1_len < 0:
                    db_gen['NP1_LENGTH'] = 0
                    overlap = abs(np1_len)
                else:
                    db_gen['NP1_LENGTH'] = np1_len
                    n1_start = (db_gen['V_SEQ_START'] + db_gen['V_SEQ_LENGTH']-1)
                    n1_end = int(d_align[8])-1
                    seq_vdj += db_gen['SEQUENCE_INPUT'][n1_start:n1_end]

                # Query sequence positions
                db_gen['D_SEQ_START'] = int(d_align[8]) + overlap
                db_gen['D_SEQ_LENGTH'] = max(int(d_align[9]) - db_gen['D_SEQ_START'] + 1, 0)
                # Germline positions
                db_gen['D_GERM_START'] = int(d_align[10]) + overlap
                db_gen['D_GERM_LENGTH'] = max(int(d_align[11]) - db_gen['D_GERM_START'] + 1, 0)

                # Update VDJ sequence, removing insertions
                start = overlap
                for m in re.finditer(r'-', d_align[18]):
                    ins = m.start()
                    seq_vdj += d_align[17][start:ins]
                    start = ins + 1
                seq_vdj += d_align[17][start:]

            # J alignment results
            if db_gen['J_CALL'] != "-":
                j_align = hit_df[hit_df[0] == 'J'].iloc[0]

                db_gen['J_ALIGNMENT'] = j_align[4]
                db_gen['J_MUTATION'] = sum(1 for a, b in zip(j_align[17], j_align[18]) if (a != b and a !='N'))
                vl = Primer_Jlen
                db_gen['J_MUTATION_NOPRIMER'] = sum(1 for a, b in zip(j_align[17][0:-vl], j_align[18][0:-vl])
                                                    if (a != b and a !='N'))
                db_gen['J_SCORE'] =j_align[13]

                # Determine N-region length and amount of J overlap with V or D alignment
                overlap = 0
                if db_gen['D_CALL'] != "-":
                    np2_len = int(j_align[8]) - (db_gen['D_SEQ_START'] + db_gen['D_SEQ_LENGTH'])
                    if np2_len < 0:
                        db_gen['NP2_LENGTH'] = 0
                        overlap = abs(np2_len)
                    else:
                        db_gen['NP2_LENGTH'] = np2_len
                        n2_start = (db_gen['D_SEQ_START']+db_gen['D_SEQ_LENGTH']-1)
                        n2_end = int(j_align[8])-1
                        seq_vdj += db_gen['SEQUENCE_INPUT'][n2_start:n2_end]
                elif db_gen['V_CALL'] != "-":
                    np1_len = int(j_align[8]) - (db_gen['V_SEQ_START'] + db_gen['V_SEQ_LENGTH'])
                    if np1_len < 0:
                        db_gen['NP1_LENGTH'] = 0
                        overlap = abs(np1_len)
                    else:
                        db_gen['NP1_LENGTH'] = np1_len
                        n1_start = (db_gen['V_SEQ_START']+db_gen['V_SEQ_LENGTH']-1)
                        n1_end = int(j_align[8])-1
                        seq_vdj += db_gen['SEQUENCE_INPUT'][n1_start:n1_end]
                else:
                    db_gen['NP1_LENGTH'] = 0

                # Query positions
                db_gen['J_SEQ_START'] = int(j_align[8]) + overlap
                db_gen['J_SEQ_LENGTH'] = max(int(j_align[9]) - db_gen['J_SEQ_START'] + 1, 0)

                # Germline positions
                db_gen['J_GERM_START'] = int(j_align[10]) + overlap
                db_gen['J_GERM_LENGTH'] = max(int(j_align[11]) - db_gen['J_GERM_START'] + 1, 0)

                # Update VDJ sequence, removing insertions
                start = overlap
                for m in re.finditer(r'-', j_align[18]):
                    ins = m.start()
                    seq_vdj += j_align[17][start:ins]
                    start = ins + 1
                seq_vdj += j_align[17][start:]

            db_gen['SEQUENCE_VDJ'] = seq_vdj
            # Create IMGT-gapped sequence and infer IMGT junction
            if not db_gen['V_ALLELE'].endswith('_DS'):
                # print("perform gapV: is db_gen['SEQUENCE_IMGT'] present before gapV?",  db_gen['SEQUENCE_IMGT'], file=sys.stderr)
                # "-" as expected
                db_gen = gapV(db_gen, repo_dict)
                ### initialized upstream: repo_dict = getInputSeq(args.Vgapseq)
                #print("perform gapV: is db_gen['SEQUENCE_IMGT'] present after gapV()?", db_gen['SEQUENCE_IMGT'], file=sys.stderr)
                # "-" empty afterwards
                db_gen = dmask(db_gen, repo_dict)
                db_gen = mutV(db_gen)
            # Update two unique columns for D upstream
            elif Dupstream:
                db_gen['D_UPSTREAM_STITCH_D_GENE'] = 'F'
                if db_gen['V_GENE'].replace('_DS', '') in db_gen['D_CALL'].split(','):
                    db_gen['D_UPSTREAM_MATCH_D_GENE'] = 'T'
                    if db_gen['V_GENE_LEN'] == db_gen['V_GERM_END_VDJ']:
                        db_gen['D_UPSTREAM_STITCH_D_GENE'] = 'T'
                else:
                    db_gen['D_UPSTREAM_MATCH_D_GENE'] = 'F'
            # Yyx changed back 2021-06-02, bottleneck is not output
            db_writer.writerow([db_gen[f] for f in ordered_fields])
#            print('\t'.join([str(db_gen[f]) for f in ordered_fields]), file=db_handle)

    IGBLASTN_pattern = re.compile('# IGBLASTN')
    with open(IgBlast_output) as f:
        # Iterate over individual results (separated by # IgBlastN)
        # Yyx changed back on 2021-06-02, the bottleneck is not input
        for k1, block in groupby(f, lambda x: IGBLASTN_pattern.match(x)):
            if k1: continue
            parseOneBlock(list(block))
#        k1 = ''
#        block = []
#        for x in f:
#            if IGBLASTN_pattern.match(x):
#                if k1 != '':
#                    parseOneBlock(block)
#                k1 = x
#                block = []
#            else:
#                block.append(x)
#        if k1 != '':
#            parseOneBlock(block)

    db_handle.close()

def run_one_IgBlast(sample, args):
    ''' run IgBlast agaist one sample
    '''
    # run IgBlast
    if args.VDJdatabase.startswith('IG'):
        seqtype = 'Ig'
    if args.VDJdatabase.startswith('TR'):
        seqtype = 'TCR'
    for subset in args.readtypes:
        eachdir = '%s/%s' % (args.outdir, sample)
        seq_file = '%s/%s_%s.fa' % (eachdir, sample, subset)
        IgBlast_output = '%s/%s_%s.IgBlast' % (eachdir, sample, subset)
        cmdline  =  '%s/external_software/igblastn -query %s -organism %s' \
                    ' -germline_db_V %s -germline_db_D %s -germline_db_J %s' \
                    ' -auxiliary_data %s -ig_seqtype %s -domain_system %s'  \
                    ' -outfmt "7 std qlen slen btop qseq sseq" -out %s' \
                    ' -num_clonotype 0 -num_threads 2' % (
                    args.scriptdir, seq_file, args.organism,
                    args.Vdb, args.Ddb, args.Jdb,
                    args.auxiliary_data, seqtype, args.domain_system,
                    IgBlast_output)
        os.environ['IGDATA'] = '%s/database/' % args.scriptdir
        os.environ['BLASTDB'] = '%s/database/' % args.scriptdir
        loggingRun(cmdline)

def _write_arguments(logfile, args):
    with open(logfile.replace('.log', '.param'), 'w') as param_file:
        writer = csv.writer(param_file, delimiter="\t")
        for key in args.__dict__:
            writer.writerow([key, args.__dict__[key]])

def run_IgBlast(args):
    ''' Run IgBlast search for each sample
    '''
    logging.info('Runing IgBlast......')

    # perform IgBlast search
    if args.skipIgBlast is False:
        # for sample in metadict: run_one_IgBlast(sample, args)
        pool = multiprocessing.Pool(processes = args.nproc)
        for sample in args.metadict:
            pool.apply_async(run_one_IgBlast, (sample, args, ))
        pool.close()
        pool.join()
    # clean up files if no need to do IgBlast search
    # will never be run if skipIgBlast
    else:
        for sample in args.metadict:
            if os.path.exists('%s/%s/IgBlast_raw' % (args.outdir, sample)):
                os.system('mv {0}/{1}/IgBlast_raw/* {0}/{1}/'.format(args.outdir, sample))
            if os.path.exists('%s/%s/IgBlast_results' % (args.outdir, sample)):
                os.system('mv {0}/{1}/IgBlast_results/* {0}/{1}/'.format(args.outdir, sample))
            if os.path.exists('%s/%s/reads_fasta/' % (args.outdir, sample)):
                os.system('mv {0}/{1}/reads_fasta/* {0}/{1}/'.format(args.outdir, sample))
            # 04262021 JH added -q
            os.system('gunzip -q %s/%s/*.gz' % (args.outdir, sample))

def parse_one_Record(row, args):
    """
    Parse a row of one join read

    Arguments:
    row =  One record in IgBlast dataframe

    Returns:
    A list of tags
    """
    # Fail tags in join reads
    fail_tags = ['V_NO_ALIGNMENT', 'V_LOW_SCORE', 'V_LOW_IDENTITY',
                'V_LOW_COVERAGE', 'V_SHORT_ALIGNMENT', 'V_ALIGN_DIFF_GENE',
                'J_NO_ALIGNMENT', 'J_SHORT_ALIGNMENT', 'J_NOT_MATCH_JGENE',
                'NO_PRODCTIVE_INFO', 'D_STITCH_UPSTREAM_SHORT_ALIGNMENT',
                'D_UPSTREAM_SHORT_ALIGNMENT', 'R1_R2_V_NO_MATCH']
    tags = []
    # Parse D upstream alignment for upstream alignment when required
    if args.D_upstream and row['V_ALLELE'].endswith('_DS'):
        if row['V_GERM_END_VDJ'] == row['V_GENE_LEN']:
            if int(row['V_ALIGNMENT']) < args.D_upstream_stitch_length:
                tags.append('D_STITCH_UPSTREAM_SHORT_ALIGNMENT')
        else:
            if int(row['V_ALIGNMENT']) < args.D_upstream_length:
                tags.append('D_UPSTREAM_SHORT_ALIGNMENT')

    # Parse V gene
    if row['V_ALLELE'] == '-':
        tags.append('V_NO_ALIGNMENT')
    else:
        if float(row['V_SCORE']) < args.V_score:
            tags.append('V_LOW_SCORE')
        if float(row['V_IDENTITY']) < args.V_identity:
            tags.append('V_LOW_IDENTITY')
        if float(row['V_COVERAGE']) < args.V_coverage:
            tags.append('V_LOW_COVERAGE')
        if int(row['V_ALIGNMENT']) < args.V_length:
            tags.append('V_SHORT_ALIGNMENT')
        vlist = [ v.split('*')[0] for v in row['V_CALL'].split(',') ]
        if len(set(vlist)) > 1:
            tags.append('V_ALIGN_DIFF_GENE')
        if args.checkProductive and row['PRODUCTIVE'] == '-':
            tags.append('NO_PRODCTIVE_INFO')

    # Parse J gen=-09
    if not args.skipJAlignmentFilter:
        if row['J_ALLELE'] == '-':
            tags.append('J_NO_ALIGNMENT')
        else:
            if int(row['J_ALIGNMENT']) < args.J_length:
                tags.append('J_SHORT_ALIGNMENT')
            if args.J_gene and args.J_gene not in row['J_ALLELE']:
                tags.append('J_NOT_MATCH_JGENE')

    # Parse R1&R2 V allele
    if 'V_ALLELE_R2' in row:
        if row['V_ALLELE_R2'] != row['V_ALLELE']:
            tags.append('R1_R2_V_NO_MATCH')

    tags = [t for t in tags if t in fail_tags]
    if len(tags) > 0:
        return '|'.join(tags)
    else:
        return '-'

def parse_one_IgBlast(IgBlast_r1, args, IgBlast_r2=None):
    """
    Parse a joined IgBlast (R1 only) or unjoined R1 & R2

    Arguments:
    IgBlast_r1 = IgBlast join database or R1 database (if R2 specified)
    IgBlast_r2 = IgBlast R2 database (optional)
    args = Input arguments

    Returns:
    A pd dataframe containing annotation tags
    """
    if IgBlast_r2 is None:
        # R1 = pd.read_csv(IgBlast_r1, sep="\t", low_memory=False)
        R1 = getCSV(IgBlast_r1)
        if len(R1) == 0:
            return R1
        R1['JOINED'] = 'T'
        R1["NOTE"] = R1.apply(parse_one_Record, args=(args,), axis=1)
        R1["PASS"] = R1["NOTE"].apply(lambda x: 'T' if x == '-' else 'F')
        return R1
    else:
        # Fail tags in R1 & R2 reads
        R1 = getCSV(IgBlast_r1)
        R2 = getCSV(IgBlast_r2)
        if len(R1) == 0 or len(R2) == 0:
            return R1
        R1['JOINED'] = 'F'
        # Replace following fields in R2 to replace in R1
        replace_fields = ['V_SCORE', 'V_ALIGNMENT', 'V_MISMATCH',
                        'INDEL', 'V_MUTATION',
                        'V_IDENTITY', 'V_COVERAGE', 'V_BTOP',
                        'V_GERM_START_VDJ', 'V_GERM_END_VDJ']
        replace_CDR12 = ['CDR1_SEQ', 'CDR1_PEPTIDE', 'CDR2_SEQ', 'CDR2_PEPTIDE']
        # Obtain paired reads
        R1 = R1[R1['SEQUENCE_ID'].isin(R2['SEQUENCE_ID'])].sort_values('SEQUENCE_ID')
        R2 = R2[R2['SEQUENCE_ID'].isin(R1['SEQUENCE_ID'])].sort_values('SEQUENCE_ID')
        R1 = R1.reset_index(drop=True)
        R2 = R2.reset_index(drop=True)
        R1['V_SCORE_R2'] = R2['V_SCORE']
        R1['V_ALLELE_R2'] = R2['V_ALLELE']

        # Use R2 score if higher
        R1[replace_CDR12] = R2[replace_CDR12]
        idx = R1.loc[R1['V_SCORE']<R1['V_SCORE_R2']].index
        R1.loc[idx,replace_fields] = R2.loc[idx,replace_fields]

        R1["NOTE"] = R1.apply(parse_one_Record, args=(args,), axis=1)
        R1["PASS"] = R1["NOTE"].apply(lambda x: 'T' if x == '-' else 'F')
        R1.drop(['V_SCORE_R2', 'V_ALLELE_R2'], inplace=True, axis=1)

        return R1

def summarize_one_IgBlast(sample, records, args):
    ''' Summarize the V usage in one IgBlast task
    '''
    # Prepare V annotation file
    if args.Vannotation:
        annofile = args.Vannotation
    else:
        annofile = '%s/database/annotation/%s_%s_anno.txt' % (
                    args.scriptdir, args.organism, args.VDJdatabase)
        if args.mousestrain in ['B6']:
            annofile = '%s/database/annotation/mouse_%s_B6_anno.txt' % (
                    args.scriptdir, args.VDJdatabase)
    # with no specific organism and VDJdatabase, /usr/pipelines/HTGTSrep_pipeline/database/annotation/mouse_IGH_anno.txt is the default
    # /usr/pipelines/HTGTSrep_pipeline/database/annotation/mouse_IGH_anno.txt has both B6 and mm129
    # Summarize one IgBlast output
    gene_count = {}
    for key, group in records.groupby('V_GENE'):
        gene = group["V_GENE"].unique()[0]
        gene_pass = len(group)
        gene_pass_productive = len(group.loc[group['PRODUCTIVE']=='T'])
        gene_pass_non_productive = len(group.loc[group['PRODUCTIVE']=='F'])

        gene_count[gene] = [gene_pass_productive, gene_pass_non_productive, gene_pass]

    # Do stat summarize
    # First four cols are from db file and used for stat
    nameslist = ["V_GENE", "LOCUS", "FUNCTIONAL",
                "PRODUCTIVE", "NON_PRODCTIVE", "TOTAL"]
    if os.path.exists(annofile):
        V_stat = pd.read_csv(annofile, names=nameslist, sep="\t").fillna(0)
        # Filter out non-mutated reads and drop V_MUTATION col
        for gene in gene_count:
            # update from ix to loc
            gene_isin = V_stat.loc[V_stat["V_GENE"].isin([gene]), ]
            if len(gene_isin) > 0:
                index = gene_isin.index[0]
                # update from ix to iloc
                V_stat.iloc[index, 3:6] = gene_count[gene]
            else:
                newrow = [gene, "-", "-"] + gene_count[gene]
                # update from ix to iloc
                V_stat.loc[len(V_stat)+1] = newrow
                print(newrow, file=sys.stderr)
    else:
        V_stat = pd.DataFrame(columns=nameslist)
        V_stat[["LOCUS","FUNCTIONAL"]] = "-"
        for gene in gene_count:
            V_stat.loc[len(V_stat)+1,0:4] = [gene] + gene_count[gene]
    V_stat = V_stat.drop(V_stat[V_stat['V_GENE']=='V_GENE'].index)
    return V_stat

def summarize_IgBlast(i, V_stat_dict, args):
    ''' Generate a collected V stat report for pass/join reads
        i: stat type
        V_stat_dict: stat data of each type in each sample
        will be [V_stat_pass, V_stat_join, V_stat_joinMut, V_stat_joinNoMut, V_stat_joinMut_VJ, V_stat_joinNoMut_VJ]
    '''
    # Initialize the DataFrame and add column names
    column_gene = ["V_GENE", "LOCUS", "FUNCTIONAL"]
    column_count = ["PRODUCTIVE", "NON_PRODCTIVE", "TOTAL"]
                # "PRODUCTIVE_JOINED", "NON_PRODCTIVE_JOINED", "TOTAL_JOINED"]
    column_all = column_gene
    for sample in V_stat_dict:
        for item in column_count:
            a = '%s_%s' %(sample, item)
            column_all.append(a)
    V_stat_all = pd.DataFrame(columns=column_all)

    # Fill V_GENE columns
    for sample in V_stat_dict:
        for key, row in V_stat_dict[sample][i].iterrows():
            if row['V_GENE'] not in V_stat_all[['V_GENE']].values:
                V_stat_all.loc[len(V_stat_all) + 1, column_gene] = row[column_gene]

    # Fill count columns
    V_stat_all.fillna(0, inplace=True)
    for sample in V_stat_dict:
        for key, row in V_stat_dict[sample][i].iterrows():
            gene = row['V_GENE']
            column_fill = ['%s_%s' % (sample, c) for c in column_count]
            V_stat_all.loc[V_stat_all['V_GENE']==gene, column_fill] = \
            [int(val) for val in row[column_count].values]
    V_stat_all = V_stat_all.drop(V_stat_all[V_stat_all['V_GENE']=='V_GENE'].index)
    return V_stat_all

def summarize_OneSample(dirpath, sample, records_pass, args, dedup=''):
    ''' Stat four types for one sammple
    '''
    statdir = '%s/stat/' % dirpath
    if not os.path.exists(statdir):
        os.system('mkdir %s' % statdir)

    if dedup == '':
        statpath = '%s/stat/%s' % (dirpath, sample)
    else:
        statpath = '%s/stat/%s.dedup' % (dirpath, sample)

    # Stat all pass reads, skip for dedup as all are joined reads
    if dedup == '':
        V_stat_pass_output = '%s.stat.pass.xls' % (statpath)
        V_stat_pass = summarize_one_IgBlast(sample, records_pass, args)
        V_stat_pass.to_csv(V_stat_pass_output, sep="\t", index=False)
    else:
        V_stat_pass = ''

    # Stat pass, joined reads
    records_join = records_pass.loc[records_pass['JOINED']=='T']

    records_join.loc[records_join['V_MUTATION']=="-", 'V_MUTATION'] = 0
    records_join.loc[records_join['J_MUTATION_NOPRIMER']=="-", 'J_MUTATION_NOPRIMER'] = 0
    records_join.loc[:,'V_MUTATION'] = records_join.V_MUTATION.astype(int)
    records_join.loc[:,'J_MUTATION_NOPRIMER'] = records_join.J_MUTATION_NOPRIMER.astype(int)

    V_stat_join_output = '%s.stat.join.xls' % (statpath,)
    V_stat_join = summarize_one_IgBlast(sample, records_join, args)
    V_stat_join.to_csv(V_stat_join_output, sep="\t", index=False)

    records_joinMut = records_join.loc[records_join['V_MUTATION']>0]
    records_joinNoMut = records_join.loc[records_join['V_MUTATION']==0]
    # Stat pass, joined, mutonly reads
    V_stat_joinMut_output = '%s.stat.joinMut_Vonly.xls' % (statpath)
    V_stat_joinMut = summarize_one_IgBlast(sample, records_joinMut, args)
    V_stat_joinMut.to_csv(V_stat_joinMut_output, sep="\t", index=False)

    # Stat pass, joined, no-mut reads
    V_stat_joinNoMut_output = '%s.stat.joinNoMut_Vonly.xls' % (statpath)
    V_stat_joinNoMut = summarize_one_IgBlast(sample, records_joinNoMut, args)
    V_stat_joinNoMut.to_csv(V_stat_joinNoMut_output, sep="\t", index=False)

    records_joinMut_VJ = records_join.loc[ (records_join['V_MUTATION']>0) | (records_join['J_MUTATION_NOPRIMER']>0)]
    records_joinNoMut_VJ = records_join.loc[(records_join['V_MUTATION']==0) & (records_join['J_MUTATION_NOPRIMER']==0)]
    # Stat pass, joined, mutonly reads
    V_stat_joinMut_output_VJ = '%s.stat.joinMut_VJ.xls' % (statpath)
    V_stat_joinMut_VJ = summarize_one_IgBlast(sample, records_joinMut_VJ, args)
    V_stat_joinMut_VJ.to_csv(V_stat_joinMut_output_VJ, sep="\t", index=False)

    # # Stat pass, joined, no-mut reads
    V_stat_joinNoMut_output_VJ = '%s.stat.joinNoMut_VJ.xls' % (statpath)
    V_stat_joinNoMut_VJ = summarize_one_IgBlast(sample, records_joinNoMut_VJ, args)
    V_stat_joinNoMut_VJ.to_csv(V_stat_joinNoMut_output_VJ, sep="\t", index=False)

    return [V_stat_pass, V_stat_join, V_stat_joinMut, V_stat_joinNoMut, V_stat_joinMut_VJ, V_stat_joinNoMut_VJ]

def WriteDbAndStat(records_all, dirpath, sample, args):
    # Summarize V gene usage in IgBlast outputs, the stat order is
    records_pass = records_all[records_all['PASS']=='T'].copy()
    records_jp = records_pass[records_pass['JOINED']=='T'].copy()
    if not args.skipDedup and len(records_jp) > 2:
        if args.dedup_by_col:
            records_dedup = collapse_db(records_jp, 'V1', 'F')
        else:
            records_dedup = collapse_db(records_jp, 'identical', 'F')
    else:
        records_dedup = records_jp.copy()

    V_stat_dict_kedup[sample] = summarize_OneSample(dirpath, sample, records_pass, args)
    V_stat_dict_dedup[sample] = summarize_OneSample(dirpath, sample, records_dedup, args, 'identical')

    # Gen diff database files
    IgBlast_jp_db = '%s/%s.db.xls' % (dirpath, sample)
    records_jp.to_csv(IgBlast_jp_db, sep="\t", index=False)

    # Drop V_GAP and V_BTOP in the output
    droplist = ["V_BTOP", "STRAND",
    'D_SEQ_START', 'D_SEQ_LENGTH', 'D_GERM_START', 'D_GERM_LENGTH',
    'J_SEQ_START', 'J_SEQ_LENGTH', 'J_GERM_START', 'J_GERM_LENGTH',
    "NP1_LENGTH", "NP2_LENGTH",
    "V_SEQ_START", "V_SEQ_LENGTH", "V_GENE_LEN", "V_GENE_GAP_LEN",
    "V_GERM_START_VDJ", "V_GERM_END_VDJ",
    "V_GERM_START_IMGT", "V_GERM_LENGTH_IMGT",
    "V_ALLELE_NUC", "GERMLINE_IMGT_D_MASK",
    "SEQUENCE_INPUT", "SEQUENCE_VDJ", "SEQUENCE_IMGT"]
    records_all.drop(droplist, inplace=True, axis=1)
    records_pass.drop(droplist, inplace=True, axis=1)
    records_dedup.drop(droplist, inplace=True, axis=1)
    # All reads without seqs
    IgBlast_all_output = '%s/%s.all.xls' % (dirpath, sample)
    records_all.to_csv(IgBlast_all_output, sep="\t", index=False)
    # Passed reads without seqs
    IgBlast_pass_output = '%s/%s.pass.xls' % (dirpath, sample)
    records_pass.to_csv(IgBlast_pass_output, sep="\t", index=False)
    # Dedup reads without seqs
    IgBlast_dedup_output = '%s/%s.dedup.xls' % (dirpath, sample)
    records_dedup.to_csv(IgBlast_dedup_output, sep="\t", index=False)

## 2021-05-20, Adam_Yyx add this variable and function writeIgBlastDb_for_oneSampleSubset for parallel
def writeIgBlastDb_for_oneSampleSubset(sample_subset_pair, args, force_output=True):
    sample, subset = sample_subset_pair
    dirpath = '%s/%s' % (args.outdir, sample)
    if args.Vgapseq:
        repo_dict = getInputSeq(args.Vgapseq)

    IgBlast_output = '%s/%s_%s.IgBlast' % (dirpath, sample, subset)
    seq_file = '%s/%s_%s.fa' % (dirpath, sample, subset)
    IgBlast_db = '%s/%s_%s.IgBlast.db' % (dirpath, sample, subset)
    ## 2021-05-20, Adam_Yyx add writeIgBlastDb_oneSample() for parallel
    if os.path.exists(IgBlast_output):
        if force_output or not os.path.exists(IgBlast_db):
            logging.info('start writeIgBlastDb for %s_%s ...' % (sample, subset))
            writeIgBlastDb(IgBlast_output, IgBlast_db, seq_file, repo_dict, subset, args.metadict[sample][3], args.D_upstream)
            logging.info('finish writeIgBlastDb for %s_%s' % (sample, subset))

def writeIgBlastDb_for_oneSample(sample, args, force_output=True):
    for subset in args.readtypes:
        writeIgBlastDb_for_oneSampleSubset((sample, subset), args, force_output)

def writeIgBlastDb_for_oneSubset(subset, args, force_output=True):
    for sample in args.metadict:
        writeIgBlastDb_for_oneSampleSubset((sample, subset), args, force_output)

def parse_IgBlast(args, force_output=True):
    logging.info('Parsing IgBlast......')

    # Initialize summary dict, sample: V_stat_dataframe
    global V_stat_dict_kedup, V_stat_dict_dedup
    V_stat_dict_kedup = OrderedDict()
    V_stat_dict_dedup = OrderedDict()
    if args.Vgapseq:
        repo_dict = getInputSeq(args.Vgapseq)
        # print("if there is args.Vgapseq, args.Vgapseq is:",args.Vgapseq,file=sys.stderr)

    # Parse IgBlast output in each sample

    ## 2021-05-20, Adam_Yyx modified the content of smple loop for parallel
    if args.parallelParseIgBlast.lower() == 'none':
        for sample in args.metadict:
            dirpath = '%s/%s' % (args.outdir, sample)

            # Generate IgBlast db
            for subset in args.readtypes:
                IgBlast_output = '%s/%s_%s.IgBlast' % (dirpath, sample, subset)
                seq_file = '%s/%s_%s.fa' % (dirpath, sample, subset)
                IgBlast_db = '%s/%s_%s.IgBlast.db' % (dirpath, sample, subset)
                ## 2021-05-20, Adam_Yyx add writeIgBlastDb_oneSample() for parallel
                if force_output or not os.path.exists(IgBlast_db):
                    logging.info('start writeIgBlastDb for %s_%s ...' % (sample, subset))
                    writeIgBlastDb(IgBlast_output, IgBlast_db, seq_file, repo_dict, subset, args.metadict[sample][3], args.D_upstream)
                    logging.info('finish writeIgBlastDb for %s_%s' % (sample, subset))
    elif args.parallelParseIgBlast.lower() == 'sample':
        sample_list = list(args.metadict)
        with multiprocessing.Pool(len(sample_list)) as pool:
            pool.map(functools.partial(writeIgBlastDb_for_oneSample, args=args, force_output=force_output), sample_list)
        logging.info('multiprocessing writeIgBlastDb are all done.')
    elif args.parallelParseIgBlast.lower() == 'subset':
        subset_list = list(args.readtypes)
        with multiprocessing.Pool(len(subset_list)) as pool:
            pool.map(functools.partial(writeIgBlastDb_for_oneSubset, args=args, force_output=force_output), subset_list)
        logging.info('multiprocessing writeIgBlastDb are all done.')
    elif args.parallelParseIgBlast.lower() == 'sample_subset':
        sample_subset_pair_list = list(itertools.product(args.metadict, args.readtypes))
        with multiprocessing.Pool(len(sample_subset_pair_list)) as pool:
            pool.map(functools.partial(writeIgBlastDb_for_oneSampleSubset, args=args, force_output=force_output), sample_subset_pair_list)
        logging.info('multiprocessing writeIgBlastDb are all done.')

    ## 2021-05-20, Adam_Yyx separate the content of sample loop for parallel
    for sample in args.metadict:
        dirpath = '%s/%s' % (args.outdir, sample)

        # Parse IgBlast db
        IgBlast_db_join = '%s/%s_join.IgBlast.db' % (dirpath, sample)
        if args.input:
            print('args.input == True, line 907, igblast. This is rarely used. Is this intentional', file = sys.stderr)
            records_all = parse_one_IgBlast(IgBlast_db_join, args)
        else:
            records_join = parse_one_IgBlast(IgBlast_db_join, args)
            # Parse unjoined reads
            if not args.skipUnjoined:
                IgBlast_db_R1 = '%s/%s_unjoinR1.IgBlast.db' % (dirpath, sample)
                IgBlast_db_R2 = '%s/%s_unjoinR2.IgBlast.db' % (dirpath, sample)
                records_unjoin = parse_one_IgBlast(IgBlast_db_R1, args, IgBlast_db_R2)
                records_all = pd.concat([records_join, records_unjoin])
            else:
                records_all = records_join
        #records_all.to_csv('records_all.xls',sep="\t")

        # Write outputs and stats
        if len(records_all) != 0:
            # records_uncollapsed = pd.DataFrame(columns=list(records_all))
            # for index, row in records_all.iterrows():
            #     for id in row['SEQUENCE_ID'].split('_'):
            #         row_add = row
            #         row_add['SEQUENCE_ID'] = id
            #         records_uncollapsed = records_uncollapsed.append(row_add, ignore_index=True)
            WriteDbAndStat(records_all, dirpath, sample, args)

    # Summarize the V usage in all samples, keep-dup and de-dup
    stat_types = ['stat_pass', 'stat_join', 'stat_joinMut_Vonly',
                'stat_joinNoMut_Vonly', 'stat_joinMut_VJ', 'stat_joinNoMut_VJ']
    for i in range(0, len(stat_types)):
        V_stat_all_output_kedup = '%s/stat/allsample.%s.xls' % (args.outdir, stat_types[i])
        V_stat_all_kedup = summarize_IgBlast(i, V_stat_dict_kedup, args)
        V_stat_all_kedup.to_csv(V_stat_all_output_kedup, sep="\t", index=False)
        if i != 0:
            V_stat_all_output_dedup = '%s/stat/allsample.dedup.%s.xls' % (args.outdir, stat_types[i])
            V_stat_all_dedup = summarize_IgBlast(i, V_stat_dict_dedup, args)
            V_stat_all_dedup.to_csv(V_stat_all_output_dedup, sep="\t", index=False)
