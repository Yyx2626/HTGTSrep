#!/usr/bin/env python3
import os, re, sys, logging, csv, multiprocessing

import pandas as pd
from HTGTSrep.mutprofile import profile_Protein, stat_Protein
from HTGTSrep.lib import getDNADistMatrix, calcDistances, formClusters, \
                         hier_clust, collapse_db, profile_DNAmut, getInputSeq, \
                         getInferSeq, profile_DNAmut_clonal, profile_DNAmut_clonal_errbar
from scipy.stats import sem

def define_clones(records, args):
    '''
    Pre-group records based on V, J and junction length
    '''
    if args.cluster_by_gene:
        records['J_GENE'] = records['J_ALLELE'].str.split('*').str[0]
        ### 09242020 added 'V_ALLELE' to grouped = records.groupby(['V_ALLELE', 'J_GENE', 'JUNCTION_LENGTH'])
        grouped = records.groupby(['V_GENE', 'J_GENE', 'JUNCTION_LENGTH'])
    else:
        ### 09242020 added 'V_ALLELE' to grouped = records.groupby(['V_ALLELE', 'J_GENE', 'JUNCTION_LENGTH'])
        grouped = records.groupby(['V_ALLELE', 'J_ALLELE', 'JUNCTION_LENGTH'])
    records['clonetmp'] = '-'
    records['CLONE'] = '-'
    for key, group in grouped:
        if len(group) == 1:
            ### 09152020 Lawrence: updated from records.ix to records.loc
            records.loc[group.index, 'clonetmp'] = group['SEQUENCE_ID']
        else:
            clonelist = hier_clust(group, args.dist)
            ### 09162020 Lawrence: updated from records.set_value(group.index, 'clonetmp', clonelist)
            ### to records.at[group.index, 'clonetmp'] = clonelist
            records.at[group.index, 'clonetmp'] = clonelist
    clone_num = 1
    for key, group in records.groupby('clonetmp'):
        ### 09152020 Lawrence: updated from records.ix to records.loc
        records.loc[group.index, 'CLONE'] = clone_num
        clone_num += 1
    records.drop('clonetmp', inplace=True, axis=1)
    records.sort_values('CLONE', ascending=True, inplace=True)
    return records

def write_cloneDb(records, outputfile, args):
    droplist = ["V_BTOP", "STRAND",
    'D_SEQ_START', 'D_SEQ_LENGTH', 'D_GERM_START', 'D_GERM_LENGTH',
    'J_SEQ_START', 'J_SEQ_LENGTH', 'J_GERM_START', 'J_GERM_LENGTH',
    "NP1_LENGTH", "NP2_LENGTH",
    "V_SEQ_START", "V_SEQ_LENGTH", "V_GENE_LEN", "V_GENE_GAP_LEN",
    "V_GERM_START_VDJ", "V_GERM_END_VDJ",
    "V_GERM_START_IMGT",
    "V_ALLELE_NUC", "GERMLINE_IMGT_D_MASK",
    "SEQUENCE_INPUT", "SEQUENCE_VDJ"]
    records_output = records.drop(droplist, axis=1)
    records_output.to_csv(outputfile, sep="\t", index=False)

def AAprofile_clones(group, sample, args):
    ''' Do protein profiling
    '''
    clone = group['CLONE'].values[0]
    protein_file = "%s/%s_clonal/protein_text/%s.clone%d.%s.%s.protein.txt" % (
                    args.outdir, sample, sample, clone, args.muttype, args.productivetype)
    protein_file = re.sub(r'[*)(]', '-', protein_file)
    protein_PDF = "%s/%s_clonal/protein_profile/%s.clone%d.%s.%s.eps" % (
                    args.outdir, sample, sample, clone, args.muttype, args.productivetype)
    protein_PDF = re.sub(r'[*)(]', '-', protein_PDF)
    Vgapseq = getInputSeq(args.__dict__['params_dict']['Vgapseq'])
    reads_db = profile_Protein(group, protein_file, protein_PDF, Vgapseq, args)
    stat_file = "%s/%s_clonal/protein_text/%s.clone%d.%s.%s.protein.txt" % (
                    args.outdir, sample, sample, clone, args.muttype, args.productivetype)
    stat_Protein(protein_file, Vgapseq, group['V_ALLELE'].values[0], reads_db)

def DNAprofile_clones(group, sample, args, profileType):
    """
    Generate DNA and Protein mutation profiles for clones
    """
    clone = group['CLONE'].values[0]
    V_ALLELE = group['V_ALLELE'].values[0]
    if profileType == 'sepCluster':
        dirprefix = "%s/%s_clonal" % (args.outdir, sample)
    elif profileType == 'mixCluster':
        dirprefix = "%s/allsample_clonal" % args.outdir

    nuc_profile = "%s/nucl_text/clone%s.%s.%s.%s.nuc.txt" % (
                dirprefix, clone, sample, args.muttype, args.productivetype)
    nuc_stat = "%s/nucl_text/clone%s_%s.%s.%s.%s.stat.txt" % (
                dirprefix, clone, V_ALLELE, sample, args.muttype, args.productivetype)
    nuc_PDF = "%s/nucl_profile/clone%s.%s.%s.%s.pdf" % (
                dirprefix, clone, sample, args.muttype, args.productivetype)

    profile_DNAmut(group, nuc_stat, nuc_PDF, nuc_profile, args)

def DNAprofile_clones_errbar(group, sample_errbar, args):
    ''' Collect stat data from diff samples
    '''
    clone = group['CLONE'].values[0]
    allele = group["V_ALLELE"].values[0]
    dirprefix = "%s/allsample_clonal/" % args.outdir
    sample_files = {}
    for sample in sample_errbar:
        nuc_PDF = "%s/nucl_profile/clone%s.%s.%s.%s.pdf" % (
                dirprefix, clone, sample, args.muttype, args.productivetype)
        nuc_stat = "%s/nucl_text/clone%s_%s.%s.%s.%s.stat.txt" % (
                dirprefix, clone, allele, sample, args.muttype, args.productivetype)
        sample_files[sample] = [nuc_PDF, nuc_stat]
        ### print to find if the correct file paths are generated. trying to solve missing stat file
        # print("finding ", sample_files[sample],file = sys.stderr)
    nuc_PDF_errbar = "%s/nucl_profile_errbar/clone%s.%s.%s.errbar.pdf" % (
            dirprefix, clone, args.muttype, args.productivetype)
    nuc_stat_errbar = "%s/nucl_text_errbar/clone%s_%s.%s.%s.stat.errbar.txt" % (
            dirprefix, clone, allele, args.muttype, args.productivetype)

    # Load CDR3
    if allele in args.__dict__['V_CDR']:
        cdr = args.__dict__['V_CDR'][allele]
        cdrstring = 'cdr1_start=%s cdr1_end=%s cdr2_start=%s cdr2_end=%s ' \
                    'cdr3_start=%s cdr3_end=%s' % (cdr[0], cdr[1], cdr[2],
                    cdr[3], cdr[4], cdr[5])
    else:
        cdrstring = ''

    if len(sample_errbar) == 1:
        nuc_PDF_sample = sample_files[sample_errbar[0]][0]
        print("DNAprofile_clones_errbar did not output a file, because len(sample_errbar) == 1", file = sys.stderr)
        # if os.path.exists(nuc_PDF_sample):
        #     os.system('cp %s %s' % (nuc_PDF_sample, nuc_PDF_errbar))
    else:
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
        for i in range(1, pos_max+1):
            group = stat_all.loc[stat_all['Pos']==i]
            Pos = i
            Mut = group['Mut'].sum()
            Total = group['Total'].sum()
            Base = list(group['Base'])[0]
            A = group['A'].sum()
            T = group['T'].sum()
            C = group['C'].sum()
            G = group['G'].sum()

            group_Y = group[group['Total']>0]
            Err = 0
            if len(group_Y) == 0: Y = 0
            else: Y = group_Y['Y'].sum()/len(group_Y)
            if len(group_Y) > 1:
                Err = sem(group_Y['Y'])
            stat_new.loc[i] = [int(Pos), int(Mut), int(Total), str(Base), Y,
                               int(A), int(T), int(C), int(G), Err]
        stat_new.to_csv(nuc_stat_errbar, sep="\t", index=False)
        if not args.skipTree:
            os.system('Rscript %s/HTGTSrep/R/SHMPlot2.R %s %s plotrows=1 figureheight=2 '
                    'showsequence=FALSE ymax=0.75 %s ' % (args.scriptdir,
                    nuc_stat_errbar, nuc_PDF_errbar, cdrstring))

def Tree_clones(records, sample, args, sample_errbar=[]):
    ''' Parsing records and construct lineage Tree
    '''
    records_Dmask = records.copy()
    records_tree = records.copy()
    records_Dmask['Dmask_N'] = [seq.count('N') for seq in records['GERMLINE_IMGT_D_MASK']]
    records_Dmask['SEQUENCE_IMGT_LENGTH'] = records_Dmask['SEQUENCE_IMGT'].apply(len)
    records_Dmask.sort_values('Dmask_N', ascending=False, inplace=True)

    clone_seqlongest = {}
    for key, group in records_Dmask.groupby('CLONE'):
        # Use GERMLINE_IMGT_D_MASK with most 'N' as consensus one
        seq_Dmask = group['GERMLINE_IMGT_D_MASK'].values[0]
        records_tree.loc[records_tree['CLONE']==key, 'GERMLINE_IMGT_D_MASK'] = seq_Dmask

        # padding shorter SEQUENCE_IMGT to the longest length or trim it by the germline seq length
        seq_IMGT_longest = group['SEQUENCE_IMGT_LENGTH'].max()
        for index, row in group.iterrows():
            padding_length = seq_IMGT_longest - row['SEQUENCE_IMGT_LENGTH']
            SEQUENCE_ID = row['SEQUENCE_ID']
            SEQUENCE_IMGT = row['SEQUENCE_IMGT'] + 'N' * padding_length
            SEQUENCE_IMGT = SEQUENCE_IMGT[0:len(seq_Dmask)]
            records_tree.loc[records_tree['SEQUENCE_ID']==SEQUENCE_ID, 'SEQUENCE_IMGT'] = SEQUENCE_IMGT

    # Get lineage tree from collapsed records
    ### these two statements help find keyError
    # print("sample in treeClone =",sample,file=sys.stderr)
    # print("have allcolumns before collapse?", len(records_tree["SEQUENCE_INPUT"]))
    records_collapse = collapse_db(records_tree, 'partial', 'F')
    records_collapse['DUPCOUNT'] = 1
    for index, row in records_collapse.iterrows():
        readlist = row['DUPREAD'].split(',')
        records_collapse.loc[index, 'DUPCOUNT'] = len(readlist)
        if sample != 'allsample':
            records_collapse.loc[index, 'SHORTCOUNT'] = 's1:%d' % len(readlist)
            records_collapse.loc[index, 'SAMPLECOUNT'] = '%s:%d' % (sample, len(readlist))
        else:
            sample_order = list(args.__dict__['sample_path'].keys())
            sample_num = len(sample_order)
            samples = [records[records['SEQUENCE_ID']==read]['SAMPLE'].item() for read in readlist]
            SAMPLECOUNT = '|'.join(['%s:%d' % (sample, samples.count(sample)) \
                                    for sample in sample_order])
            SHORTCOUNT = '|'.join(['s%d:%d' % (i+1, samples.count(sample_order[i])) \
                                    for i in range(0, sample_num)])
            records_collapse.loc[index, 'SAMPLECOUNT'] = SAMPLECOUNT
            records_collapse.loc[index, 'SHORTCOUNT'] = SHORTCOUNT
    records_collapse.sort_values('CLONE', ascending=True, inplace=True)

    select_Output = ["SEQUENCE_ID", "V_CALL", "J_CALL", "CLONE", "SAMPLE", "DUPCOUNT",
                    "SHORTCOUNT", "SAMPLECOUNT", "JUNCTION_LENGTH", "CDR3_SEQ", "SEQUENCE_IMGT",
                    "GERMLINE_IMGT_D_MASK"]
    records_Output = records_collapse[select_Output]
    file_collapse = "%s/%s_clonal/s_clonal%s.collapse.xls" % (args.outdir, sample, sample)
    records_Output.to_csv(file_collapse, sep="\t", index=False)

    if not args.skipTree:
        os.system('Rscript %s/HTGTSrep/R/TREEPlot.R %s %s/external_software/dnapars %d' % \
              (args.scriptdir, file_collapse, args.scriptdir, args.min_profileread))

    # Generate mutation profile of clones using inferred seq as germline seq
    dirprefix = "%s/%s_clonal" % (args.outdir, sample)
    if sample != 'allsample':
        # For single sample, the inferred sequence is the root sequence in each clone
        for key, group in records_tree.groupby('CLONE'):
            treefile = "%s/%s_clonal/lineageTree/%d.txt" % (args.outdir, sample, key)
            if os.path.exists(treefile):
                inferseq = getInferSeq(treefile, group)
                V_ALLELE = group["V_ALLELE"].unique()[0]
                nuc_stat = "%s/nucl_text_infer/clone%s.%s.%s.%s.%s.stat.txt" % (
                            dirprefix, key, V_ALLELE, sample, args.muttype, args.productivetype)
                nuc_profile = "%s/nucl_text_infer/clone%s.%s.%s.%s.%s.nuc.txt" % (
                            dirprefix, key, V_ALLELE, sample, args.muttype, args.productivetype)
                nuc_PDF = "%s/nucl_profile_infer/clone%s.%s.%s.%s.%s.pdf" % (
                            dirprefix, key, V_ALLELE, sample, args.muttype, args.productivetype)
                profile_DNAmut_clonal(inferseq, group, nuc_stat, nuc_PDF, nuc_profile, args)
    else:
        ''' For all pooled sample
            Folder _profile and _profile_errbar: V gene only profile and with error bar
            Folder _root: tree for each sample, root seq from all sample
            Folder _infer: tree for each sample, root seq from each sample's own seq
            Folder _errbar_infer: the root seq is the root in a clone with all pooled reads
                                  same as _root, so error bar can be added
        '''
        sample_list = list(args.__dict__['sample_path'].keys())
        short_list = ["s%d" % i for i in range(1, len(sample_list)+1)]
        for key, group_clone in records_tree.groupby('CLONE'):
            V_ALLELE = group_clone["V_ALLELE"].unique()[0]
            sample_files = {}
            rootfile = "%s/allsample_clonal/lineageTree/%d.txt" % (args.outdir, key)
            if not os.path.exists(rootfile): continue
            rootseq = getInferSeq(rootfile, group_clone)

            for i in range(0, len(sample_list)):
                sample = sample_list[i]
                short = short_list[i]

                group = group_clone[group_clone['SAMPLE']==sample]
                if len(group) < args.min_profileread_sub: continue

                root_stat = "%s/nucl_text_root/clone%s.%s.%s.%s.%s.stat.txt" % (
                            dirprefix, key, V_ALLELE, sample, args.muttype, args.productivetype)
                root_profile = "%s/nucl_text_root/clone%s.%s.%s.%s.%s.nuc.txt" % (
                            dirprefix, key, V_ALLELE, sample, args.muttype, args.productivetype)
                root_PDF = "%s/nucl_profile_root/clone%s.%s.%s.%s.%s.pdf" % (
                            dirprefix, key, V_ALLELE, sample, args.muttype, args.productivetype)
                profile_DNAmut_clonal(rootseq, group, root_stat, root_PDF, root_profile, args)
                sample_files[sample] = [root_PDF, root_stat]

                treefile = "%s/allsample_clonal/lineageTree/%d.%s.txt" % (args.outdir, key, short)
                if not os.path.exists(treefile): continue
                inferseq = getInferSeq(treefile, group)

                nuc_stat = "%s/nucl_text_infer/clone%s.%s.%s.%s.%s.stat.txt" % (
                            dirprefix, key, V_ALLELE, sample, args.muttype, args.productivetype)
                nuc_profile = "%s/nucl_text_infer/clone%s.%s.%s.%s.%s.nuc.txt" % (
                            dirprefix, key, V_ALLELE, sample, args.muttype, args.productivetype)
                nuc_PDF = "%s/nucl_profile_infer/clone%s.%s.%s.%s.%s.pdf" % (
                            dirprefix, key, V_ALLELE, sample, args.muttype, args.productivetype)
                profile_DNAmut_clonal(inferseq, group, nuc_stat, nuc_PDF, nuc_profile, args)

            # Run for the profile with error bar
            nuc_PDF_errbar = "%s/nucl_profile_errbar_infer/clone%s.%s.%s.%s.errbar.pdf" % (
                    dirprefix, key, V_ALLELE, args.muttype, args.productivetype)
            nuc_stat_errbar = "%s/nucl_text_errbar_infer/clone%s.%s.%s.%s.stat.errbar.txt" % (
                    dirprefix, key, V_ALLELE, args.muttype, args.productivetype)
            if len(sample_files) > 1:
                profile_DNAmut_clonal_errbar(rootseq, group_clone, nuc_stat_errbar, \
                                            nuc_PDF_errbar, sample_files, args)

def series_analyze_onesample(records, sample, args):
    """
    A series of analysis for each clone, including:
    1. Define clones
    2. DNA mutation profile
    3. Amino Acid profile
    4. Summarize clones
    5. Lineage Tree construction
    """
    ### 09162020 Lawrence: added print(...)
    # show the current sample
    print("current sample:", sample, file=sys.stderr)
    # Define clones and write IgBlast Db files
    records = define_clones(records, args)
    outputfile = '%s/%s_clonal/%s.db_clone.xls' % (args.outdir, sample, sample)
    write_cloneDb(records, outputfile, args)

    # Generate mutation profile
    if not args.skipTree:
        Tree_clones(records, sample, args)

    # Generate AA profile
    for key, group in records.groupby('CLONE'):
        if len(group) >= args.min_profileread:
            DNAprofile_clones(group, sample, args, 'sepCluster')
            AAprofile_clones(group, sample, args)

    statfile = '%s/%s_clonal/%s.clone_stat.xls' % (args.outdir, sample, sample)
    clone_stat(records, statfile, 'sep', args)
    return records
'''
    ###ADD CONSENSUS/AA DETAILS HERE!##################################################
    tempFile = '%s/%s_clonal/%s.clone_stat_temp.xls' % (args.outdir, sample, sample)
    scriptLocation = '%s/HTGTSrep/translate_consensus_Clonal.py' % (args.scriptdir) #may need to add %s/HTGTSrep/translate_consensus_clonal
    os.system("python3 {0} {1} > {2}".format(scriptLocation, statfile, tempFile))
    os.system("mv {0} {1}".format(tempFile, statfile))
    
    ###################################################################################
'''


# def selection_analyze():
#     os.system('Rscript %s/external_software/baseline/Baseline_Main_Version1.3.r' \
#               '1 2 1 1 0 0 1:26:38:55:65:104 sample.fasta testfold test' % ())

def diversity_analyze(outputfile, args):
    '''
    Run Diversity.R to analyze diversity and abundance, required alakazam R package
    '''
    if not args.skipDiversity:
        os.system('Rscript %s/HTGTSrep/R/Diversity.R %s %s/allsample_clonal/allsample' % (
                    args.scriptdir, outputfile, args.outdir))

def clone_stat(records, statfile, sampletype, args):
    '''
    Statistically summarize the reads in diff clones
    sampletype:
        mix: clones clustered of all sample reads
        sep: clones clustered of single sample reads
    '''
    statcols = ['CLONE', 'SAMPLE_NUM', 'READ_NUM', 'SAMPLE_DETAIL', 'SAMPLE_RATIO',
                'V_ALLELE', 'D_ALLELE', 'J_ALLELE', 'JUNC_LEN', 'JUNC_NUM', 'JUNC_DETAIL']
    statDf = pd.DataFrame(columns=statcols)
    grouped = records.groupby('CLONE')
    for key, group in grouped:
        CLONE = group['CLONE'].values[0]
        SAMPLE_NUM = len(group['SAMPLE'].unique())
        READ_NUM = len(group)

        if sampletype == 'mix':
            sc = group['SAMPLE'].value_counts()
            SAMPLE_DETAIL = '|'.join(['%s:%d' % (sample, sc.get(sample, 0)) \
                                    for sample in args.__dict__['sample_path']])
            SAMPLE_RATIO = '|'.join(['%s:%.4f' % (sample, float(sc.get(sample, 0))/len(records[records['SAMPLE']==sample])) \
                                    for sample in args.__dict__['sample_path']])
        else:
            SAMPLE_DETAIL = group['SAMPLE'].values[0]
            SAMPLE_RATIO = 1
        V_ALLELE = group["V_ALLELE"].value_counts().keys()[0]
        D_ALLELE = group["D_ALLELE"].value_counts().keys()[0]
        J_ALLELE = group["J_ALLELE"].values[0]
        JUNC_LEN = len(group['CDR3_SEQ'].values[0])
        JUNC_NUM = len(group['CDR3_SEQ'].unique())

        junc_count = group['CDR3_SEQ'].value_counts()
        JUNC_DETAIL = '|'.join(['%s:%d' % (junc, junc_count[junc]) \
                                for junc in group['CDR3_SEQ'].unique()])

        statDf.loc[len(statDf)+1] = [CLONE, SAMPLE_NUM, READ_NUM, SAMPLE_DETAIL, SAMPLE_RATIO,
                                    V_ALLELE, D_ALLELE, J_ALLELE, JUNC_LEN, JUNC_NUM, JUNC_DETAIL]

    if sampletype == 'sep': statDf.drop('SAMPLE_RATIO', inplace=True, axis=1)
    statDf.to_csv(statfile, sep="\t", index=False)

def series_analyze_allsample(records, samplelist, args):
    ''' Series analyze of pooled samples
    '''
    records.CLONE = records.SAMPLE + '.' + records.CLONE.map(str)
    outputfile = '%s/allsample_clonal/allsample.sep_clone.xls' % (args.outdir)
    records.to_csv(outputfile, sep="\t", index=False)

    # Do diversity/abundance analysis using clones clustered within each sample
    diversity_analyze(outputfile, args)

    # Do clonal clustering of all sample reads
    records.drop('CLONE', inplace=True, axis=1)
    records = define_clones(records, args)
    outputfile = '%s/allsample_clonal/allsample.mix_clone.xls' % (args.outdir)
    records.to_csv(outputfile, sep="\t", index=False)
   
    # Generate mutation profile for each sample
    for key, group in records.groupby('CLONE'):
        sample_errbar = []
        # key is each clone
        # group is df
        for sample, subgroup in group.groupby('SAMPLE'):
            ### args.min_profileread_sub = 10
            if len(subgroup) >= args.min_profileread_sub:
                sample_errbar.append(sample)
                DNAprofile_clones(subgroup, sample, args, 'mixCluster')
        if len(group) >= args.min_profileread:
            # Gen all sample profile & with err bar
            DNAprofile_clones(group, "allsample", args, 'mixCluster')
            AAprofile_clones(group, "allsample", args)
            if len(sample_errbar) > 1:
                DNAprofile_clones_errbar(group, sample_errbar, args)

    statfile = '%s/allsample_clonal/allsample.mix_clone.stat.xls' % (args.outdir)
    clone_stat(records, statfile, 'mix', args)
    #########CREATE MASTER TLX HERE ################
    masterstatfile = '%s/allsample_clonal/allsample.master.mix_clone.stat.xls' % (args.outdir)
    #tempFile = '%s/%s_clonal/%s.clone_stat_temp.xls' % (args.outdir, sample, sample)#
    listSamples = ""
    for sample in samplelist:
        clonestatfile = '%s/%s_clonal/%s.clone_stat.xls' % (args.outdir, sample, sample)
        listSamples += clonestatfile + " "
    print("samplelist=",samplelist, file = sys.stderr)

    ##create master xls
    scriptLocation = '%s/HTGTSrep/junctionsPerLibs.py' % (args.scriptdir)
    os.system("python3 {0} {1} {2} > {3}".format(scriptLocation, statfile, listSamples, masterstatfile)) ###junctionsperlibs
    
    ###create lib detail file
    libdetailfile = '%s/allsample_clonal/allsample.lib_detail.xls' % (args.outdir)
    libScript = '%s/HTGTSrep/libConsensus_clonal.py' % (args.scriptdir)
    os.system("python3 {0} {1} > {2}".format(libScript, masterstatfile, libdetailfile))
    
    ###order information in master file
    tempMaster = '%s/allsample_clonal/allsample.master_temp.mix_clone.stat.xls' % (args.outdir)
    orderlibScript = '%s/HTGTSrep/orderLibDetail.py' % (args.scriptdir)
    os.system("python3 {0} {1} > {2}".format(orderlibScript, masterstatfile, tempMaster))
    os.system("mv {0} {1}".format(tempMaster, masterstatfile))
    
    ###add sample ration and sort 6/11
    addRatioScript = '%s/HTGTSrep/add_sampleratio_sort.py' % (args.scriptdir)
    os.system("python3 {0} {1} {2}".format(addRatioScript, statfile, listSamples))
    #print("COMMAND: python3 "+addRatioScript+" "+statfile+" "+listSamples)
    
    ##Screen 6/11
    screenScript = '%s/HTGTSrep/screen_master_stat.py' % (args.scriptdir)
    screenoutputfile = '%s/allsample_clonal/%s.master.mix_clone.stat_screen.xls' % (args.outdir, args.outdir)
    os.system("python3 {0} {1} > {2}".format(screenScript, masterstatfile, screenoutputfile))
    
    
    ###clean up files
    os.system("mv {0} {1}".format('%s/allsample_clonal/allsample.master.mix_clone.stat.xls' % (args.outdir), '%s/allsample_clonal/%s.master.mix_clone.stat.xls' % (args.outdir, args.outdir)))
    os.system("mv {0} {1}".format('%s/allsample_clonal/allsample.lib_detail.xls' % (args.outdir), '%s/allsample_clonal/%s.lib_detail.xls' % (args.outdir, args.outdir)))
    os.system("mv {0} {1}".format('%s/allsample_clonal/allsample.mix_clone.stat.xls' % (args.outdir), '%s/allsample_clonal/%s.mix_clone.stat.xls' % (args.outdir, args.outdir)))
    os.system("mv {0} {1}".format('%s/allsample_clonal/allsample.mix_clone.xls' % (args.outdir), '%s/allsample_clonal/%s.mix_clone.xls' % (args.outdir, args.outdir)))
    os.system("mv {0} {1}".format('%s/allsample_clonal/allsample.sep_clone.xls' % (args.outdir), '%s/allsample_clonal/%s.sep_clone.xls' % (args.outdir, args.outdir)))


    # Generate stat file for shared clones
    ### TEHSE HAbr BEEN GENERATED ABOVE  SO THESE TWO LINES ARE COMMENTED OUT
    # statfile = '%s/allsample_clonal/allsample.mix_clone.stat.xls' % (args.outdir)
    # clone_stat(records, statfile, 'mix', args)

    # Generate lineage Tree
    if not args.skipTree:
        Tree_clones(records, 'allsample', args, sample_errbar)

def clonal_main(args):
    logging.info('Loading reads database')
    # Collapse reads
    sample_path = args.__dict__['sample_path']
    pool = multiprocessing.Pool(processes = len(sample_path))
    results = []
    for sample in sample_path:
        path = sample_path[sample]
        dbpath = path + '/%s.db.xls' % sample
        records = pd.read_csv(dbpath, sep="\t")
        print("initial records file", dbpath, file=sys.stderr)
        # Filter records using V coverage, CDR3 length, whether 'N' in CDR3
        records = records.loc[records['V_COVERAGE'] > args.min_Vcov, ]
        records["JUNCTION_LENGTH"] = records.CDR3_SEQ.map(len)
        records = records.loc[records['JUNCTION_LENGTH'] > 1, ]
        records = records.loc[records['SEQUENCE_IMGT'] != '-', ]

        records['SAMPLE'] = sample
        records['CDR3_MASK'] = [re.sub('[\.-]', 'N', seq) for seq in records.CDR3_SEQ]

        if args.skipCDR3withN:
            records = records[~records['CDR3_MASK'].str.contains("N")]

        if args.muttype == 'MutOnly':
            records = records.loc[records['V_MUTATION'] > 0, ]
        elif args.muttype == 'noMut':
            records = records.loc[records['V_MUTATION'] == 0, ]
        if args.productivetype == 'P':
            records = records.loc[records['PRODUCTIVE'] == 'T', ]
        elif args.productivetype == 'NP':
            records = records.loc[records['PRODUCTIVE'] == 'F', ]
        # Run series analysis with multiprocessing   result = series_analyze_onesample(records, sample, args)
     #   print("mid loop: sample is", sample, file=sys.stderr)
        result = pool.apply_async(series_analyze_onesample, (records, sample, args,))
     #   print("result.get()",result.get(),file=sys.stderr)
        results.append(result)
    pool.close()
    pool.join()

    # Analysis in all samples
    #print("frame: results=",results, file = sys.stderr)
    ### frame: results= [<multiprocessing.pool.ApplyResult object at 0x7fe95d9d13d0>, ....,]
    frames = [result.get() for result in results]
    records_allsample = pd.concat(frames, ignore_index=True)
    series_analyze_allsample(records_allsample, sample_path.keys(), args)
