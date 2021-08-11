'''import sys
import operator
from Bio.Seq import Seq
try:
    from Bio.Alphabet import generic_dna, IUPAC
    Bio_Alphabet = True
except ImportError:
    Bio_Alphabet = None
    # usages of generic_dna, IUPAC are not supported in Biopython 1.78 (September 2020).
    print(f"The installed BioPython is a new version that has removed the Alphabet module.",file=sys.stderr)

def filter(statFile):
    # statfile = '%s/allsample_clonal/allsample.mix_clone.stat.xls' % (args.outdir)
    ListOfLines = []
    labels = statFile.readline()
    litems = labels[:-1].split("\t")

    ListOfLines.append(litems)
    for line in statFile:
        items = line[:-1].split("\t")
        ListOfLines.append(items)

    return ListOfLines

def consensus(lines, CDR3_AA):
    consensusLines = []
    lines[0].append("CONSENSUS_SEQ")
    jidx = lines[0].index("JUNC_DETAIL")
    consensusLines.append(lines[0])
    for line in lines[1:]:
        # h = 1
        # cdr = ''
        juncDets = line[jidx]
        # commented out useless lines JH 06032021
        # for item in juncDets.split('|'):
        #     i1 = item.split(':')[0]
        #    i2 = int(item.split(':')[1])
        #   if i2 > h:
        #       h = i2
        #       cdr = i1
        if CDR3_AA != "T":
            consensus = get_consensus(juncDets)
        else:
            consensus = get_consensus_AA(juncDets)
        line.append(consensus)
        consensusLines.append(line)

    return consensusLines

def get_consensus_AA(allCDR):
    pos_base_num = {}
    cdr3_len = len(allCDR.split('|')[0].split(':')[0])
    for i in range(0, cdr3_len):
        pos_base_num[i] = {"A": 0, "R": 0, "N": 0, "D": 0, "C": 0, "Q": 0, "E": 0, "G": 0, "H": 0, "I": 0, "L": 0, "K": 0, "M": 0, "F": 0, "P": 0, "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0}
    for seg in allCDR.split('|'):
        j = seg.split(':')
        for i in range(0, cdr3_len):
            pos_base_num[i][j[0][i]] += int(j[1])
    consensus = ''
    for i in range(0, cdr3_len):
        consensus += max(pos_base_num[i].items(), key=operator.itemgetter(1))[0]
    return consensus

def get_consensus(allCDR):
    pos_base_num = {}
    cdr3_len = len(allCDR.split('|')[0].split(':')[0])
    for i in range(0, cdr3_len):
        pos_base_num[i] = {'A':0, 'T':0, 'C':0, 'G':0, "N":0}
    for seg in allCDR.split('|'):
        j = seg.split(':')
        for i in range(0, cdr3_len):
            pos_base_num[i][j[0][i]] += int(j[1])
    consensus = ''
    for i in range(0, cdr3_len):
        consensus += max(pos_base_num[i].items(), key=operator.itemgetter(1))[0]
    return consensus

def translate(listLines, CDR3_AA):
    dnas = []
    listLines[0].append("AA_SEQUENCE")
    dnas.append(listLines[0])
    conSeq = listLines[0].index("CONSENSUS_SEQ")
    # i=0
    if CDR3_AA != "T":
        for line in listLines[1:]:
            seq = line[conSeq]
            while len(seq)%3 != 0:
                seq += "N"
            if Bio_Alphabet:
                AA = Seq(seq, generic_dna).translate()
            else:
                AA = Seq(seq).translate()
            # i+=1
            line.append(str(AA))
            dnas.append(line)
    else:
        for line in listLines[1:]:
            seq = line[conSeq]
            AA = seq
            # i+=1
            line.append(str(AA))
            dnas.append(line)

    return dnas

def main():

    # statfile = '%s/allsample_clonal/allsample.mix_clone.stat.xls' % (args.outdir)
    clonestat = open(sys.argv[1], "r")
    CDR3_AA = sys.argv[2]
    toParse = filter(clonestat)
    toTranslate = consensus(toParse,CDR3_AA)
    translated = translate(toTranslate,CDR3_AA) #already translated allsample.stat file, a list of lines in a file

    libsfile = sys.argv[3:] #read in each library's clonestat file to create lib_detail
    # # useless lines JH 06032021
    # c=0
    # libscloneDict = {} ##lib, clone: cdr3seq, num
    # cdr3dict = {} ###seq, num :  lib, clone
    listOfSingleLibraryDicts = []
    for library in libsfile:
        libstatfile = open(library, "r")
        libDict = {} ##append junction details to the library dictionary, for each clone (cdr3seq,V-allele, J-allele): sample, readnum
        labels = libstatfile.readline()[:-1].split("\t")
        juncIdx = labels.index("JUNC_DETAIL")
        vidx = labels.index("V_ALLELE")
        jidx = labels.index("J_ALLELE")
        sIdx = labels.index("SAMPLE_DETAIL")
        for line in libstatfile: #iterate through all the clones in a single library
            items = line[:-1].split("\t")
            v_allele = items[vidx]
            j_allele = items[jidx]
            juncDetails = items[juncIdx].split("|")
            sample = items[sIdx]
            for j in juncDetails:
                seqSpecs = j.split(":")
                cdr3seq = seqSpecs[0]
                readNum = seqSpecs[1]
                if (cdr3seq, v_allele, j_allele) not in libDict:
                    libDict[(cdr3seq, v_allele, j_allele)] = (sample, readNum)
                #else:
                #    print(cdr3seq, readNum)
        listOfSingleLibraryDicts.append(libDict)

    # print(listOfSingleLibraryDicts)

    translated[0].append("LIB_DETAIL") # will contain all lines of master clone file
    labels = translated[0] # clonestat.readline()[:-1].split("\t")
    idxJD = labels.index("JUNC_DETAIL")
    # idxClone = labels.index("CLONE")
    idxV = labels.index("V_ALLELE")
    idxJ = labels.index("J_ALLELE")
    for line in translated[1:2]: #clonestat:
        v_allele = line[idxV]
        j_allele = line[idxJ]
        juncDetails = line[idxJD].split("|")
        libDetailString = ''
        for j in juncDetails:
            seqSpecs = j.split(":")
            cdr3seqJ = seqSpecs[0]
            # readNum = seqSpecs[1]
            tuple2check = (cdr3seqJ, v_allele, j_allele)
            print(tuple2check,seqSpecs)
            for dict in listOfSingleLibraryDicts:
                if tuple2check in dict:
                    libDetailString += dict[tuple2check][0] + ":" + cdr3seqJ + ":" + dict[tuple2check][1] + "|"
                    print(libDetailString)

        if libDetailString[-1] == "|":
            line.append(libDetailString[:-1])
        else:
            line.append(libDetailString)

    for i in translated:
        print("\t".join(i))


main()
'''

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

def filter(statFile):
    # statfile = '%s/allsample_clonal/allsample.mix_clone.stat.xls' % (args.outdir)
    ListOfLines = []
    labels = statFile.readline()
    # updated to .replace so that the last character of the last line will not be accidentally deleted JH 06042021
    # litems = labels[:-1].split("\t")
    litems = labels.replace("\n", "").split("\t")

    ListOfLines.append(litems)
    for line in statFile:
        # updated to .replace so that the last character of the last line will not be accidentally deleted JH 06042021
        # items = line[:-1].split("\t")
        items = line.replace("\n", "").split("\t")
        ListOfLines.append(items)

    return ListOfLines

def consensus(lines):
    consensusLines = []
    lines[0].append("CONSENSUS_SEQ")
    jidx = lines[0].index("JUNC_DETAIL")
    consensusLines.append(lines[0])
    for line in lines[1:]:
        # h = 1
        # cdr = ''
        juncDets = line[jidx]
        # commented out useless lines JH 06032021
        # for item in juncDets.split('|'):
        #     i1 = item.split(':')[0]
        #    i2 = int(item.split(':')[1])
        #   if i2 > h:
        #       h = i2
        #       cdr = i1
        consensus = get_consensus(juncDets)
        line.append(consensus)
        consensusLines.append(line)

    return consensusLines

def get_consensus(allCDR):
    pos_base_num = {}
    cdr3_len = len(allCDR.split('|')[0].split(':')[0])
    for i in range(0, cdr3_len):
        pos_base_num[i] = {'A':0, 'T':0, 'C':0, 'G':0, "N":0}
    for seg in allCDR.split('|'):
        j = seg.split(':')
        for i in range(0, cdr3_len):
            pos_base_num[i][j[0][i]] += int(j[1])
    consensus = ''
    for i in range(0, cdr3_len):
        consensus += max(pos_base_num[i].items(), key=operator.itemgetter(1))[0]
    return consensus

def translate(listLines):
    dnas = []
    listLines[0].append("AA_SEQUENCE")
    dnas.append(listLines[0])
    conSeq = listLines[0].index("CONSENSUS_SEQ")
    # i=0
    for line in listLines[1:]:
        seq = line[conSeq]
        while len(seq) % 3 != 0:
            seq += "N"
        if Bio_Alphabet:
            AA = Seq(seq, generic_dna).translate()
        else:
            AA = Seq(seq).translate()
        # i+=1
        line.append(str(AA))
        dnas.append(line)
    return dnas

def main():

    # NOTE statfile = '%s/allsample_clonal/allsample.mix_clone.stat.xls' % (args.outdir)
    clonestat = open(sys.argv[1], "r")
    toParse = filter(clonestat)
    toTranslate = consensus(toParse)
    translated = translate(toTranslate) #already translated allsample.stat file, a list of lines in a file

    libsfile = sys.argv[2:] # read in each library's clonestat file to create lib_detail
    # # useless lines JH 06032021
    # c=0
    # libscloneDict = {} ##lib, clone: cdr3seq, num
    # cdr3dict = {} ###seq, num :  lib, clone
    listOfSingleLibraryDicts = []
    for library in libsfile:
        libstatfile = open(library, "r")
        libDict = {} ##append junction details to the library dictionary, for each clone (cdr3seq,V-allele, J-allele): sample, readnum
        # updated to .replace so that the last character of the last line will not be accidentally deleted JH 06042021
        # labels = libstatfile.readline()[:-1].split("\t")
        labels = libstatfile.readline().replace("\n", "").split("\t")
        juncIdx = labels.index("JUNC_DETAIL")
        vidx = labels.index("V_ALLELE")
        jidx = labels.index("J_ALLELE")
        sIdx = labels.index("SAMPLE_DETAIL")
        for line in libstatfile: #iterate through all the clones in a single library
            # updated to .replace so that the last character of the last line will not be accidentally deleted JH 06042021
            # items = line[:-1].split("\t")
            items = line.replace("\n", "").split("\t")
            v_allele = items[vidx]
            j_allele = items[jidx]
            juncDetails = items[juncIdx].split("|")
            sample = items[sIdx]
            for j in juncDetails:
                seqSpecs = j.split(":")
                cdr3seq = seqSpecs[0]
                readNum = seqSpecs[1]
                if (cdr3seq, v_allele, j_allele) not in libDict:
                    libDict[(cdr3seq, v_allele, j_allele)] = (sample, readNum)
                #else:
                #    print(cdr3seq, readNum)
        listOfSingleLibraryDicts.append(libDict)

    # print(listOfSingleLibraryDicts)

    translated[0].append("LIB_DETAIL") # will contain all lines of master clone file
    labels = translated[0] # clonestat.readline()[:-1].split("\t")
    idxJD = labels.index("JUNC_DETAIL")
    # idxClone = labels.index("CLONE")
    idxV = labels.index("V_ALLELE")
    idxJ = labels.index("J_ALLELE")
    for line in translated[1:]: #clonestat:
        v_allele = line[idxV]
        j_allele = line[idxJ]
        juncDetails = line[idxJD].split("|")
        libDetailString = ''
        for j in juncDetails:
            seqSpecs = j.split(":")
            cdr3seqJ = seqSpecs[0]
            # readNum = seqSpecs[1]
            tuple2check = (cdr3seqJ, v_allele, j_allele)
            for dict in listOfSingleLibraryDicts:
                if tuple2check in dict:
                    libDetailString += dict[tuple2check][0] + ":" + cdr3seqJ + ":" + dict[tuple2check][1] + "|"

        if libDetailString[-1] == "|":
            line.append(libDetailString[:-1])
        else:
            line.append(libDetailString)

    for i in translated:
        print("\t".join(i))


main()


