import sys


def main():

    masterFile = open(sys.argv[1])
    labels = masterFile.readline()[:-1].split("\t")

    allLines = [] ##print this for output
    allLines.append(labels)
    lidx = labels.index("LIB_DETAIL")
    sidx = labels.index("SAMPLE_DETAIL")
    
    loutput = "\t".join(labels)
    print(loutput)
    for line in masterFile:
        # updated to .replace so that the last character of the last line will not be accidentally deleted JH 06092021
        # items = line[:-1].split("\t")
        items = line.replace("\n", "").split("\t")
        lib_detail = items[lidx].split("|")
        sample_detail= items[sidx].split("|")
        
        
        for s in range(len(sample_detail)): ###order the samples
            idx = sample_detail[s].index(":")
            sample = sample_detail[s][:idx]
            sample_detail[s] = sample

        tempDict = {} #sample_name : all junctions
        for l in range(len(lib_detail)):
            name = lib_detail[l][:lib_detail[l].index(":")]
            if name not in tempDict:
                tempDict[name] = [lib_detail[l]]
            else:
                tempDict[name].append(lib_detail[l])

        newLibDetailString = ""
        for key in tempDict:
            string = "|".join(tempDict[key]) + "|"
            newLibDetailString += string

        #print(newLibDetailString)
        items[lidx] = newLibDetailString


        outputString = "\t".join(items)
        print(outputString)


main()
