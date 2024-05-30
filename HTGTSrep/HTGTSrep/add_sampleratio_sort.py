##add SAMPLE_RATIO column like allsample stat file between D and E, then sort list by READ_NUM largest to smallest
###python3 scripts/add_sampleratio_sort.py Alt391LC_V2_qs20/clonal_REDO_NEW/allsample_clonal/allsample.mix_clone.stat.xls Alt391LC_V2_qs20/clonal_REDO_NEW/HC547_Alt391_clonal/HC547_Alt391.clone_stat.xls Alt391LC_V2_qs20/clonal_REDO_NEW/HC549_Alt391_clonal/HC549_Alt391.clone_stat.xls
###DO BY COL INSTEAD

import sys, os

def main():


    for smplstat in sys.argv[2:]:
        tempFILE = open(smplstat[:-3]+"_temp.xls", "w")
        
        idx1 = smplstat.index("_clonal/")
        idx2 = smplstat.index(".clone_stat.xls")
        sampleName = smplstat[idx1+8:idx2]#[:smplstat.index("_clonal")]
        ssfile = open(smplstat, "r")
        allLines = ssfile.readlines()
        
        labels = allLines.pop(0)[:-1].split("\t")#ssfile.readline()[:-1].split("\t")
        readidx = labels.index("READ_NUM")
        totalReads = 0
        for l in allLines:
            items = l[:-1].split("\t")
            reads = int(items[readidx])
            totalReads += reads
        
        newLabels = labels[:]
        newLabels.insert(4, "SAMPLE_RATIO")
        readnumDict = {}
        for line in allLines:
            items = line[:-1].split("\t")
            clone_read_num = int(items[readidx])
            ratio = round(clone_read_num/totalReads, 4)
            items.insert(4, str(ratio))

            read_num = items[readidx]
            RDS = ""
            if "." in str(read_num):
                RDS = int(str(read_num)[:str(read_num).index(".")]) ###depends on python version?
            else:
                RDS = int(read_num)
            

            if RDS in readnumDict:
                readnumDict[RDS].append(items)
            else:
                readnumDict[RDS] = [items]

        ##compile list of items to print to a file that is sorted
        lines2write = []
        while len(readnumDict.keys()) != 0:
            #print(max(readnumDict.keys()))
            lines = readnumDict.pop(max(readnumDict.keys()))
            if len(lines) >1:
                for c in lines:
                    lines2write.append(c)
            else:
                lines2write.append(lines[0])

        ls = "\t".join(newLabels)
        tempFILE.write(ls + "\n")
        for p in lines2write:
            tempFILE.write("\t".join(p) + "\n")
            #print(p)

        os.system("mv {0} {1}".format(smplstat[:-3]+"_temp.xls", smplstat))





main()
