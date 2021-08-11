#python3 libConsensus_clonal.py master.xls > allsample.libdetail.xls
import sys
import operator
from Bio.Seq import Seq
try:
    from Bio.Alphabet import generic_dna, IUPAC
    Bio_Alphabet = True
except ImportError:
    Bio_Alphabet = None
    # usages of generic_dna, IUPAC are not supported in Biopython 1.78 (September 2020).
    print(f"The installed BioPython is a new version that has removed the Alphabet module.",file=sys.stderr)
###Read in master.clone.stat.xls

def get_consensus(allCDR):
    nthits = {}
    seq_hits = {}
    pos_base_num = {}
    cdr3_len = len(allCDR[0][0])
    for i in range(0, cdr3_len):
        pos_base_num[i] = {'A':0, 'T':0, 'C':0, 'G':0, 'N':0}
    for seg in allCDR:
        j = seg
        for i in range(0, cdr3_len):
            pos_base_num[i][j[0][i]] += int(j[1])
    consensus = ''
    for i in range(0, cdr3_len):
        consensus += max(pos_base_num[i].items(), key=operator.itemgetter(1))[0]
    return(consensus)

def translate(seq):
    while len(seq)%3 != 0:
        seq += "N"
    if Bio_Alphabet:
        AA = Seq(seq, generic_dna).translate()
    else:
        AA = Seq(seq).translate()
    return AA

def main():

    masterfile = open(sys.argv[1], "r")
    labels = masterfile.readline()[:-1].split("\t")
    libIdx = labels.index("LIB_DETAIL")
    cloneIdx = labels.index("CLONE")

    print("CLONE" + "\t" + "LIBRARY" + "\t" + "CONSENSUS_SEQ" + "\t" + "AA_SEQUENCE")
    for line in masterfile:
        lineItems = line[:-1].split("\t")
        libsline = lineItems[libIdx]
        if libsline[-1] == "|":
            libsline = lineItems[libIdx][:-1]
        clone = lineItems[cloneIdx]
        libs = libsline.split("|")

        libDict = {} ##libname = (seq1,#, seq2:#...)
        for lib in libs:
            libdeets = lib.split(":")
            if libdeets[0] not in libDict:
                libDict[libdeets[0]] = [(libdeets[1], libdeets[2])]
            else:
                libDict[libdeets[0]].append((libdeets[1], libdeets[2]))

        for libname in libDict:
            consensSeq = get_consensus(libDict[libname])
            AAseq = translate(consensSeq)
            print(clone + "\t" + libname + "\t" + consensSeq + "\t" + AAseq)


main()
