###Screen SAMPLE_RATIO with thresholds 0.01, 0.005, 0.003 and keep output in same file, but separate groups with an empty row. Between B and C add threshold used, and between C and D the number of samples passing the screening. SAMPLE_DETAIL and SAMPLE_RATIO should contain only the sampled that passed the screen. Sort each group by number of samples passed screening from largest to smallest. Samples passed should be no larger than SAMPLE_NUM Output file : allsample.master.mix_clone.stat_screen.xls
##python3 scripts/screen_master_stat.py Alt391LC_V2_qs20/clonal_REDO_NEW/allsample_clonal/allsample.master.mix_clone.stat.xls > allsample.master.mix_clone.stat_screen.xls ##don't show clones with only 1 sample passed


import sys

def screen(amfile, threshold, labels):
    thresholhd_List = []
    sample_detal_idx = labels.index("SAMPLE_DETAIL")
    sample_ratio_idx = labels.index("SAMPLE_RATIO")
    sample_num_idx = labels.index("SAMPLE_NUM")
    
    new_sample_detail = ""
    new_sample_ratio = ""
    threshold_used = threshold
    samples_passed_screening = ""
    for line in amfile:
        items=line[:-1].split("\t")
        sampleratios = items[sample_ratio_idx].split("|")
        passed_thresh = [] ##list of sample:ratio, sample2:ratio
        passed_thresh_names = [] ##sample1, sample2
        for s in sampleratios: #s= sample:ratio
            i = s.split(":")
            samplename = i[0]
            ratio = float(i[1])
            if ratio >= threshold:
                passed_thresh.append(s)
                passed_thresh_names.append(samplename)

        passed_details = []
        sampledetails = items[sample_detal_idx].split("|")
        for det in sampledetails:
            d= det.split(":")
            for name in passed_thresh_names:
                if name in det:
                    passed_details.append(det)

        samples_passed_screening = str(len(passed_thresh_names))
        threshold_used = str(threshold)

        items[sample_ratio_idx] = "|".join(passed_thresh) #update sample_ratio
        items[sample_detal_idx] = "|".join(passed_details) #update sample_detail
        items.insert(2, threshold_used)
        items.insert(3, int(samples_passed_screening)) #changed to int

        if len(passed_thresh_names) > 1: ###from 0 originally
            thresholhd_List.append(items)

    return thresholhd_List

def sort_samples(thresholdList):

    sortDict = {}
    for i in thresholdList:
        num_samples_passed = i[3]
        if num_samples_passed in sortDict:
            sortDict[num_samples_passed].append(i)
        else:
            sortDict[num_samples_passed] = [i]

    sortedList = []
    while len(sortDict.keys()) != 0:
        maxItem = sortDict.pop(max(sortDict.keys()))
        if len(maxItem) > 1:
            #maxItem[0][3] = str(maxItem[0][3])
            for m in maxItem:
                m[3] = str(m[3])
                sortedList.append("\t".join(m))
        else:
            maxItem[0][3] = str(maxItem[0][3])
            sortedList.append("\t".join(maxItem[0]))

    return sortedList


def main():

    allsample_master = open(sys.argv[1], "r")
    labels = allsample_master.readline()[:-1].split("\t")
    allLines = allsample_master.readlines()
    newLabels = labels[:]
    newLabels.insert(2, "THRESHOLD")
    newLabels.insert(3, "SAMPLES_PASSED_THRESHOLD")
    list_thresholds = [0.01, 0.005, 0.003]
    list_screened = [] #list of screened groups
    
    for t in list_thresholds:
        screened_output = screen(allLines, t, labels)
        sorted_output = sort_samples(screened_output)
        list_screened.append(sorted_output)

    print("\t".join(newLabels))
    for group in list_screened:
        for l in group:
            print(l)

        print("\n")




main()
