#!/usr/bin/env python3
__author__    = 'Zhou Du'
__copyright__ = 'Copyright 2017 Alt Lab, Harvard University. All rights reserved.'
__license__   = 'Creative Commons Attribution-NonCommercial-ShareAlike 4.0 Unported'
__version__   = '0.1'
__date__      = '2017.03.02'

import os, re, sys, logging, csv

from pandas import read_csv, DataFrame
from argparse import ArgumentParser
from collections import OrderedDict
import multiprocessing as mp

def check_customdatabase(type, path):
    if type == 'database':
        for suffix in ['.nhr', '.nin', '.nog', '.nsd', '.nsi', '.nsq']:
            filename = path + suffix
            if not os.path.exists(filename):
                logging.error('Cannot find or read database file %s' % filename)
                sys.exit()
    elif type == 'auxiliary':
        if not os.path.exists(path):
            logging.error('Cannot find or read database file %s' % path)
            sys.exit()

def check_metafile(args, metapath):
    # define global variabla, samplename: [barcode/MID, primer, adapter]e
    logging.info('Parsing meta files.....')
    metainfo = read_csv(metapath, sep="\t").fillna("")
    metadict = OrderedDict()
    if 'Primer_Jlen' not in metainfo.keys():
        metainfo['Primer_Jlen'] = 0
    for feature in ['Library', 'Sequencing', 'Primer', 'MID', 'Adapter']:
        if feature not in metainfo.keys():
            logging.error('Cannot find %s feature in metadata file' % feature)
    for i in range(0, len(metainfo)):
        Library = metainfo.loc[i, 'Library'].upper()
        if Library == '': continue
        Primer = metainfo.loc[i, 'Primer'].upper()
        MID = metainfo.loc[i, 'MID'].upper()
        Adapter = metainfo.loc[i, 'Adapter'].upper()
        Primer_Jlen = int(metainfo.loc[i, 'Primer_Jlen'])
        if Primer == '' and not args.input:
            logging.error('Missing primer sequence in library %s' % Library)
            sys.exit()
        if Adapter == '' and not args.input:
            logging.error('Missing adapter sequence in library %s' % Library)
            sys.exit()
        sn = '%s_%s' % (metainfo.loc[i,'Library'], metainfo.loc[i, 'Sequencing'])
        metadict[sn] = [MID, Primer, Adapter, Primer_Jlen]
    args.metadict = metadict
    # Ensure same number of metadata sample and fasta file
    if args.subcmd in ['preprocess', 'run'] and args.input \
        and len(metadict) != len(args.input.split(',')):
        logging.error('Sample records in metadata file not equal to fasta files')
        sys.exit()

def create_dir(args):
    if args.subcmd in ['preprocess', 'run']:
        logdir = '%s/logs' % args.outdir
    if args.subcmd in ['mut', 'clonal']:
        logdir = '%s/logs' % args.dir

    if not os.path.exists(logdir):
        os.system('mkdir -p %s' % logdir)
    if args.subcmd == "preprocess":
        logfile = logdir + "/pipeline.preprocess.log"
    if args.subcmd == "run":
        statdir = '%s/stat' % args.outdir
        if not os.path.exists(statdir):
            os.system('mkdir -p %s' % statdir)
        logfile = logdir + '/pipeline.run.log'
    if args.subcmd == "mut":
        logfile = logdir + '/pipeline.mutprofile.log'
    # if args.subcmd == "compare":
    #     logdir = args.output + '.compare.log'
    if args.subcmd == "clonal":
        logfile = logdir + '/pipeline.clonal.log'
    return logfile

def custom_logging(logfile):
    # init logging
    logging.basicConfig( level = 10,
    format = '%(levelname)-5s @ %(asctime)s: %(message)s ',
    datefmt = '%a, %d %b %Y %H:%M:%S',
    filename = logfile,
    filemode = 'a' )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(levelname)-5s @ %(asctime)s: %(message)s ','%a, %d %b %Y %H:%M:%S')
    #formatter.formatTime('%a, %d %b %Y %H:%M:%S')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

def write_arguments(logfile, args):
    with open(logfile.replace('.log', '.param'), 'w') as param_file:
        writer = csv.writer(param_file, delimiter="\t")
        for key in args.__dict__:
            if key == 'metadict':
                slist = args.__dict__[key].keys()
                rlist = ['{0}/{1}/reads_fasta/{1}_join.fa.gz'.format(
                            args.outdir, s) for s in slist]
                writer.writerow(['sample_name', ','.join(slist)])
                writer.writerow(['sample_joinfa', ','.join(rlist)])
            else:
                writer.writerow([key, args.__dict__[key]])

def check_preprocess(args):
    ''' Check sequence files, metadata file
    '''
    if not (args.r1 and args.r2) and not args.input:
        logging.error('No sequence file is specified')
        sys.exit()
    if args.r1 and not os.path.exists(args.r1):
        logging.error('Cannot find or read read 1 file %s' % args.r1)
        sys.exit()
    if args.r2 and not os.path.exists(args.r2):
        logging.error('Cannot find or read read 2 file %s' % args.r2)
        sys.exit()
    if args.input:
        for inp in args.input.split(','):
            if not os.path.exists(inp):
                logging.error('Cannot find or read fasta file %s' % inp)
                sys.exit()
    if args.fastq_quality_trimmer < 0:
        logging.error('Cannot accept < 0 quality score')

    # Check metadata file
    if not os.path.exists(args.metafile):
        logging.error('Cannot find or read meta file %s' % args.metafile)
        sys.exit()
    else:
        metapath = "%s/metadata.txt" % args.outdir
        # Convert metafile to unix format
        try:
            os.system("perl -p -e 's/\\r/\\n/g' %s > %s" % (
                        args.metafile, metapath))
        except:
            logging.error('Problem converting line feeds to unix format')
            sys.exit()

    # Load metadata info
    check_metafile(args, metapath)

def check_IgBlast(args):
    # check IgBlast program and internal data
    if not os.path.exists('%s/external_software/igblastn' % args.scriptdir):
        logging.error('Cannot find igblastn program at IgBlast_bin/')
        sys.exit()
    if not os.path.exists('%s/database/internal_data' % args.scriptdir):
        logging.error('Cannot find internal_data data database/')
        sys.exit()
    # if args.dedup_by_seq and args.dedup_by_seq < 100:
    #     logging.warning('Over-dedup can happen by only using %d bps' % args.dedup_by_seq)
    if args.auxiliary_data:
        check_customdatabase('auxiliary', args.auxiliary_data)
    if args.Vdb:
        check_customdatabase('database', args.Vdb)
        # if (not args.Valignfile) and args.genomealign == 'T':
        #     logging.error('Must specify V gene alignment file for %s' % args.Vdb)
        #     sys.exit()
    if args.Ddb:
        check_customdatabase('database', args.Ddb)
    if args.Jdb:
        check_customdatabase('database', args.Jdb)
    if args.Vannotation:
        if not os.path.exists(args.Vannotation):
            logging.error('Cannot find or read V gene annotation file %s' % args.Vannotation)
            sys.exit()
    if args.Vgapseq:
        if not os.path.exists(args.Vgapseq):
            logging.error('Cannot find or read V IMGT seq file %s' % args.Vgapseq)
            sys.exit()
    else:
        args.Vgapseq = '%s/database/annotation/%s_%s_GAP.fasta' % (
                    args.scriptdir, args.organism, args.VDJdatabase)
        if not os.path.exists(args.Vgapseq):
            args.Vgapseq = None
            logging.warning('Cannot find or read V IMGT seq file %s' % args.Vgapseq)
    if args.mousestrain == 'B6':
        if args.organism == 'human':
            logging.error('Cannot specify mouse strain in human')
            sys.exit()
        if args.VDJdatabase not in ['IGH', 'IGK']:
            logging.error('Only support IGH or IGK in mouse')
            sys.exit()
        # if args.add_D_upstream == 'mm9_AJ851868ins':
        #     logging.error('Cannot specify mouse strain while including D upstreams')
        #     sys.exit()
    # if args.genomealign == 'T' and (not args.chrseq):
    #     logging.error('Must specify one chromosome to do alignment')
    #     sys.exit()
    # if args.chrseq or args.Valignfile:
    #     args.genomealign = 'T'
    # if args.Valignfile:
    #     if not os.path.exists(args.Valignfile):
    #         logging.error('Cannot find or read V alignment file %s' % args.Valignfile)
    #         sys.exit()

def prepare_IgBlast(args):
    # use prepared VDJ databases
    dpath = '%s/database/database-imgt' % (args.scriptdir)

    if args.VDJdatabase == 'IGALL':
        Vdb_prep = '%s/%s_IGallV_imgt' % (dpath, args.organism)
        Ddb_prep = '%s/%s_IGallD_imgt' % (dpath, args.organism)
        Jdb_prep = '%s/%s_IGallJ_imgt' % (dpath, args.organism)
    elif args.VDJdatabase == 'TRALL':
        Vdb_prep = '%s/%s_TRallV_imgt' % (dpath, args.organism)
        Ddb_prep = '%s/%s_TRallD_imgt' % (dpath, args.organism)
        Jdb_prep = '%s/%s_TRallJ_imgt' % (dpath, args.organism)
    else:
        Vdb_prep = '%s/%s_%sV_imgt' % (dpath, args.organism, args.VDJdatabase)
        Ddb_prep = '%s/%s_%sD_imgt' % (dpath, args.organism, args.VDJdatabase)
        Jdb_prep = '%s/%s_%sJ_imgt' % (dpath, args.organism, args.VDJdatabase)
        if args.mousestrain == 'B6':
            Vdb_prep = Vdb_prep + "_" + args.mousestrain
            Ddb_prep = Ddb_prep + "_" + args.mousestrain

    # user customized databases
    if not args.Vdb: args.Vdb = Vdb_prep
    if not args.Ddb: args.Ddb = Ddb_prep
    if not args.Jdb: args.Jdb = Jdb_prep

    # user customized auxiliary file
    if not args.auxiliary_data:
        args.auxiliary_data = '%s/database/optional_file/%s_gl.aux' % (
                                args.scriptdir, args.organism)
    # Read types need to be parsed
    if args.input or args.skipUnjoined: args.readtypes = ['join']
    else: args.readtypes = ['join', 'unjoinR1', 'unjoinR2']

    # Prep files for skipdemultiplex and skipIgBlast
    if args.skipDemultiplex:
        for sample in args.metadict:
            sampledir = '%s/%s' % (args.outdir, sample)
            os.system('mv %s/reads_fasta/* %s' % (sampledir, sampledir))
            os.system('gunzip %s/*.gz' % sampledir)
    if args.skipIgBlast:
        for sample in args.metadict:
            sampledir = '%s/%s' % (args.outdir, sample)
            os.system('mv %s/reads_fasta/* %s/igblast_raw/* %s' % (
                        sampledir, sampledir, sampledir))
            os.system('gunzip %s/*.gz' % sampledir)

    if not args.nproc: args.nproc = len(args.metadict)

def check_params(args, path=None):
    if path is None:
        # Read paramaters
        paramfile = '%s/logs/pipeline.run.param' % args.dir
        params_dict = {}
        if not os.path.exists(paramfile):
            logging.error('Cannot find param file %s' % paramfile)
            sys.exit()
        else:
            with open(paramfile) as f:
                for line in f:
                    (key, val) = line.split('\t')
                    params_dict[key] = val.strip()
        args.params_dict = params_dict
    else:
        paramfile ='%s/../logs/pipeline.run.param' % path
        params_dict = {}
        if not os.path.exists(paramfile):
            logging.error('Cannot find param file %s' % paramfile)
            sys.exit()
        else:
            with open(paramfile) as f:
                for line in f:
                    (key, val) = line.split('\t')
                    params_dict[key] = val.strip()
        if "params_dict" not in args.__dict__:
            args.__dict__['params_dict'] = params_dict
        else:
            for p in ['Vdb', 'Jdb', 'Ddb', 'domain_system', 'VDJdatabase']:
                if args.__dict__['params_dict'][p] != params_dict[p]:
                    logging.error('The sample diff in %s' % p)
                    sys.exit()

def check_CDRfile(args):
    # Prep CDR boundary file
    if args.genecdr and not os.path.exists(args.genecdr):
        logging.error('Cannot find or read CDR information file %s' % args.genecdr)
        sys.exit()
    else:
        args.genecdr = '%s/database/annotation/%s_CDR.txt' % (args.scriptdir, \
                        args.params_dict['organism'])
    V_CDR = {}
    with open(args.genecdr) as f:
        for line in f:
            l = line.strip().split()
            V_CDR[l[0]] = l[1:]
    args.V_CDR = V_CDR

def check_clonal(args):
    # Check db files
    sample_path = OrderedDict()
    for path in args.dir.strip().split(','):
        sample = path.strip('/').split('/')[-1]
        check_params(args, path)

        dbpath = path + '/%s.db.xls' % sample
        parampath = path + '/../logs/pipeline.run.param'
        if not os.path.exists(dbpath):
            sys.exit('Cannot find database file %s' % dbpath)
        if not os.path.exists(parampath):
            sys.exit('Cannot find param file %s' % parampath)
        sample_path[sample] = path
        if not os.path.exists('%s/%s_clonal' % (args.outdir, sample)):
            os.system('mkdir -p %s/%s_clonal/nucl_text' % (args.outdir, sample))
            os.system('mkdir -p %s/%s_clonal/nucl_profile' % (args.outdir, sample))
            os.system('mkdir -p %s/%s_clonal/nucl_text_infer' % (args.outdir, sample))
            os.system('mkdir -p %s/%s_clonal/nucl_profile_infer' % (args.outdir, sample))
            os.system('mkdir -p %s/%s_clonal/protein_text' % (args.outdir, sample))
            os.system('mkdir -p %s/%s_clonal/protein_profile' % (args.outdir, sample))
            os.system('mkdir -p %s/%s_clonal/lineageTree' % (args.outdir, sample))
    if not os.path.exists('%s/allsample_clonal' % args.outdir):
        os.system('mkdir -p %s/allsample_clonal/nucl_text' % args.outdir)
        os.system('mkdir -p %s/allsample_clonal/nucl_profile' % args.outdir)
        os.system('mkdir -p %s/allsample_clonal/nucl_text_infer' % args.outdir)
        os.system('mkdir -p %s/allsample_clonal/nucl_profile_infer' % args.outdir)
        os.system('mkdir -p %s/allsample_clonal/nucl_text_root' % args.outdir)
        os.system('mkdir -p %s/allsample_clonal/nucl_profile_root' % args.outdir)
        os.system('mkdir -p %s/allsample_clonal/protein_text' % args.outdir)
        os.system('mkdir -p %s/allsample_clonal/protein_profile' % args.outdir)
        os.system('mkdir -p %s/allsample_clonal/nucl_text_errbar' % args.outdir)
        os.system('mkdir -p %s/allsample_clonal/nucl_profile_errbar' % args.outdir)
        os.system('mkdir -p %s/allsample_clonal/nucl_text_errbar_infer' % args.outdir)
        os.system('mkdir -p %s/allsample_clonal/nucl_profile_errbar_infer' % args.outdir)
        os.system('mkdir -p %s/allsample_clonal/lineageTree' % args.outdir)

    args.sample_path = sample_path

    # Read params, CDR, metadata file
    check_CDRfile(args)

def check_mutprofile(args):
    # Read params, CDR, metadata file
    check_params(args)
    check_CDRfile(args)
    metapath = '%s/metadata.txt' % args.dir
    if os.path.exists(metapath): check_metafile(args, metapath)
    else:
        logging.error('Cannot find or read metadata file %s' % metapath)
        sys.exit()
    # check minimum clone read number
    # if args.min_profileread < 10:
    #     logging.error('Need at least 10 reads to do mutation profiling')
    #     sys.exit()
    # Prep VDJ gene seq file
    # if args.geneseq and not os.path.exists(args.geneseq):
    #     logging.error('Cannot find or read gene sequence file %s' % args.geneseq)
    #     sys.exit()
    # else:
    #     args.geneseq = '%s/database/annotation/%s_%s_imgt' % (args.scriptdir, \
    #                     args.params_dict['organism'],
    #                     args.params_dict['VDJdatabase'])
    # Other params
    if args.ymax_DNA <= 0 or args.ymax_protein <= 0:
        logging.error('Cannot use negative value in Y axis')
        sys.exit()
    # if args.min_profileread <= 2:
    #     logging.error('Need at least 3 reads to profile')
    #     sys.exit()
    # Check overall stat file
    statpath = '%s/stat/allsample.stat_join.xls' % args.dir
    if not os.path.exists(statpath):
        logging.error('Cannot find or read stat file %s' % statpath)
        sys.exit()
    # Check for each sample
    for sample in args.metadict:
        sampledir = '%s/%s' % (args.dir, sample)
        p = '%s/reads_fasta/%s_join.fa.gz' % (sampledir, sample)
        if not os.path.exists(p):
            logging.error('Cannot find or read fasta file %s' % p)
            sys.exit()
        t = '%s/%s.pass.xls' % (sampledir, sample)
        if not os.path.exists(t):
            logging.error('Cannot find or read IgBlast file %s' % t)
            sys.exit()
        if not os.path.exists('%s/mut_profile' % sampledir):
            os.system('mkdir -p %s/mut_profile' % sampledir)
        if not os.path.exists('%s/mut_profile/nucl_text' % sampledir):
            os.system('mkdir -p %s/mut_profile/nucl_text' % sampledir)
            os.system('mkdir -p %s/mut_profile/protein_text' % sampledir)
            os.system('mkdir -p %s/mut_profile/nucl_profile' % sampledir)
            os.system('mkdir -p %s/mut_profile/protein_profile' % sampledir)
        if not os.path.exists('%s/mut_profile/nucl_text_collapse' % sampledir) and \
            (args.collapsePartial or args.collapseIdentical):
            os.system('mkdir -p %s/mut_profile/nucl_text_collapse' % sampledir)
            os.system('mkdir -p %s/mut_profile/protein_text_collapse' % sampledir)
            os.system('mkdir -p %s/mut_profile/nucl_profile_collapse' % sampledir)
            os.system('mkdir -p %s/mut_profile/protein_profile_collapse' % sampledir)

        # if not os.path.exists('%s/mutProfile_mutreads' % sampledir):
        #     os.system('mkdir %s/mutProfile_mutreads' % sampledir)
        # if not os.path.exists('%s/clonal' % sampledir):
        #     os.system('mkdir %s/clonal' % sampledir)

def check_args(args):
    # create output fold and define directory for logging
    logfile = create_dir(args)
    custom_logging(logfile)

    # HTGTS script directory
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    scriptdir = '/'.join(scriptdir.split('/')[:-1])
    args.scriptdir = scriptdir

    print('\nWelcome to HTGTSrep pipeline!!! ')
    logging.info('Parameters: '+' '.join(sys.argv))

    # check run subcommand arguments
    if args.subcmd == 'preprocess':
        check_preprocess(args)

    if args.subcmd == 'run':
        check_preprocess(args)
        check_IgBlast(args)
        prepare_IgBlast(args)
        write_arguments(logfile, args)

    if args.subcmd == 'mut':
        check_mutprofile(args)

    # check clonal subcommand arguments
    if args.subcmd == 'clonal':
        check_clonal(args)

def parentParser(subname):
    parser = ArgumentParser(add_help=False)

    # parse args on input fastq or fasta
    if subname == 'preprocess':
        parser.add_argument("-r1", help = '''MiSeq read 1 file in fastq format,
                            mutually exclusive of fasta file''')
        parser.add_argument("-r2", help = "MiSeq read 2 file in fastq format" )
        parser.add_argument("-i", dest = "input", type = str,
                            help = '''Read files in fasta format, seperated multiple
                            files by comma in the same order as in metadata file,
                            mutually exclusive -r1 and -r2 options''' )
        parser.add_argument("-m", required = True, dest = "metafile",
                            help = "Metadata file" )
        parser.add_argument("-o", required = True, dest = "outdir", type = str,
                            help = "Output directory" )
        parser.add_argument("--keepunmatch", action='store_true',
                            dest='keepunmatch', default=False,
                            help = "Keep unmatched reads from demultiplexing")
        parser.add_argument("--demulti_length", type = int, default = 12,
                            help = "Use N (12) bp of barcode+primer to demultiplex" )
        parser.add_argument("--demulti_mismatch", type = int, default = 0,
                            help = "N-maximum mismatch (0) to demultiplex" )
        parser.add_argument("--overlap", type = int, default = 10,
                            help = "N-minimum overlap to join paired reads (10)" )
        parser.add_argument("--diffpercent", type = int, default = 8,
                            help = "N-percent maximum difference to join paired reads (8)" )
        parser.add_argument("--qscore", type = int, default = 10,
                            help = "Mark nucleotide with lower quality score \
                            as N (10) in joined sequence" )
        parser.add_argument("--qscore_cov", type = int, default = 98,
                            help = "Minimum percent of bases in a joined sequence \
                             have [--qscore] quality score (98)" )
        parser.add_argument("--fastq_quality_trimmer", type = int, default = 0,
                            help = '''Quality threshold - nucleotides with lower
                            quality score will be trimmed from the end of the
                            sequence. (0)''' )

    if subname == 'run':
        parser.add_argument("--organism", dest="organism", type = str,
                            choices = ['mouse', 'human'], default = 'mouse',
                            help = "Select sample organism, mouse or human (mouse)" )
        parser.add_argument("--mousestrain", type = str, choices = ['B6', 'all'],
                            help = '''Use mouse C57BL database only or all gene database
                            in the IgBlast searching. NO setting when organism is human (all)''')
        parser.add_argument("--VDJdatabase", type = str, default = 'IGH',
                            choices = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRG', 'TRD'],
                            help = "Specify database for searching. (IGH)" )
        parser.add_argument("--Vdb", type = str, help = "Specify IgBlast germline V database." )
        parser.add_argument("--Ddb", type = str, help = "Specify IgBlast germline D database." )
        parser.add_argument("--Jdb", type = str, help = "Specify IgBlast germline J database." )
        parser.add_argument("--domain_system", choices = ['kabat', 'imgt'], default = 'imgt',
                            help = '''Domain system kabat or imgt to be used for
                            segment annotation (imgt).''' )
        parser.add_argument("--auxiliary_data", type = str,
                            help = '''Specify file containing the coding frame
                            start positions for sequences in germline J database.''' )
        parser.add_argument("--Vannotation", type = str,
                            help = '''Specify V gene annotation file which has three
                            columns: V gene name, coordinates and function. ''' )
        parser.add_argument("--Vgapseq", type = str,
                            help = '''Specify V IMGT-gapped sequences file. ''' )
        parser.add_argument("--V_score", type = int, default = 150,
                            help = "Minimum IgBlast score of V gene (150)" )
        parser.add_argument("--V_identity", type = float, default = 0.9,
                            help = "Minimum V identity ratio (0.9)" )
        parser.add_argument("--V_coverage", type = float, default = 0.1,
                            help = "Minimum V coverage ratio (0.1)" )
        parser.add_argument("--V_length", type = int, default = 100,
                            help = "Minimum aligned length of V gene (100)" )
        parser.add_argument("--J_gene", type = str,
                            help = "Specify J gene for alignment. eg: IGHJ4" )
        parser.add_argument("--J_length", type = int, default = 34,
                            help = "Minimum alignment length of J gene (34)." )
        parser.add_argument("--checkProductive", action='store_true',
                            dest='checkProductive', default=False,
                            help = "Include only records with productive information" )
        parser.add_argument("--skipDedup", action='store_true',
                            dest='skipDedup', default=False,
                            help = "Specify to skip generate duplicated-removed output" )
        parser.add_argument("--skipIgBlast", action='store_true',
                            dest='skipIgBlast', default=False,
                            help = "Specify to skip IgBlast and all upstream steps" )
        parser.add_argument("--skipUnjoined", action='store_true',
                            dest='skipUnjoined', default=False,
                            help = "Specify to skip unjoined reads in parsing igblast db" )
        parser.add_argument("--skipDemultiplex", action='store_true',
                            dest='skipDemultiplex', default=False,
                            help = "Specifiy to skip Demultiplex and all upstream steps" )
        parser.add_argument("--D_upstream", action='store_true',
                            dest='D_upstream', default=False,
                            help = '''Specify to consider D upstream regions
                            in the analysis. The D upstream sequences should be
                            added in the V database and named as 'D gene' plus
                             '_DS'.''')
        parser.add_argument("--D_upstream_length", type = int, default = 100,
                            help = "Minimum aligned length of D upstream region (100)" )
        parser.add_argument("--D_upstream_stitch_length", type = int,
                            default = 70, help = "Minimum aligned length of D \
                            gene D stitched upstream region (70)" )
        parser.add_argument("--dedup_by_col", action='store_true',
                            dest='dedup_by_col', default=False,
                            help = "Specify to dedup in the same manner as Version 1 pipeline" )
        parser.add_argument("--nproc", type = int,
                            help = "The number of simultaneous computational \
                            processes to execute (CPU cores to utilized)" )

    if subname == 'mut':
        parser.add_argument("--dir", "-d", required = True, type = str,
                            help = "HTGTSrep result directory" )
        parser.add_argument("--genecdr", type = str,
                            help = "CDR1, 2, 3 position file" )
        parser.add_argument("--ymax_DNA", type = float, default = 0.75,
                            help = "Maximum of Y axis in DNA profile (0.75)" )
        parser.add_argument("--ymax_protein", type = float, default = 0.3,
                            help = "Maximum of Y axis in protein profile (0.3)" )
        parser.add_argument("--min_Vcov", type = float, default = 0.5,
                            help = "Minimum coverage of V gene by the read (0.5)" )
        parser.add_argument("--min_profileread", type = int, default = 20,
                            help = "Minimum assigned reads of one V allele to \
                            be mutation profiled (20)" )
        parser.add_argument("--skipStopCodonRead", action='store_true',
                            default=False, help = "Specify to exclude reads \
                            of inframe with stop codon")
        parser.add_argument("--includeAllReads", action='store_true',
                            default=False, help = 'Specify to include all passed joined \
                            reads in previous IgBlast search (Default: False, mutation only)')
        parser.add_argument("--collapsePartial", action='store_true',
                            default=False, help = "Collapse short reads to longer reads \
                            that incldue identical sequences")
        parser.add_argument("--collapseIdentical", action='store_true',
                            default=False, help = "Collapse reads to identical ones \
                            that with fewer N")

    if subname == 'clonal':
        parser.add_argument("-d", "--dir", dest="dir", required = True,
                            help = '''HTGTS result folder seperated by comma''' )
        parser.add_argument("--outdir", type = str, required = True,
                            help = "Output directory")
        # parser.add_argument("--outfile", type = str, required = True,
        #                     help = "Output name prefix")

        parser.add_argument("--genecdr", type = str,
                            help = "CDR1, 2, 3 position file" )

        parser.add_argument("--min_Vcov", type = float, default = 0.5,
                            help = "Minimum coverage of V gene by the read (0.5)" )
        parser.add_argument("--muttype", choices = ['MutOnly', 'noMut', 'noMut_filter'],
                            default = 'MutOnly',
                            help = '''Use only reads with mutation (MutOnly),
                            with NO mutation (noMut), or all reads (All)''')
        parser.add_argument("--productivetype", choices = ['P', 'NP', 'noProd_filter'],
                            default = 'noProd_filter',
                            help = """Use Productive/Non-Productive/no Productive filter
                            reads (P/NP/all, default: noProd_filter)""")
        parser.add_argument("--ymax_protein", type = float, default = 0.3,
                            help = "Maximum of Y axis in protein profile (0.3)" )

        parser.add_argument("--min_profileread", type = int, default = 20,
                            help = "The minimum read number (>=10) of one clonal \
                            for mutation profile (0)")
        parser.add_argument("--min_profileread_sub", type = int, default = 10,
                            help = "Minimum assigned reads of one V allele to \
                            be mutation profiled for each sample (10)" )
        parser.add_argument("--dist", type = float, default = 0.16,
                            help = "The distance threshold for clonal grouping (0.16)")
        parser.add_argument("--ymax_DNA", type = float, default = 0.75,
                            help = "Maximum of Y axis in DNA profile image (0.75)" )
        parser.add_argument("--skipDiversity", action='store_true',
                            default=False, help = "Skip the abundance and \
                            diversity analysis. Required: alakazam R package installed")
        parser.add_argument("--skipTree", action='store_true',
                            default=False, help = "Skip the construction of lineage tree.")
        parser.add_argument("--skipCDR3withN", action='store_true',
                            default=False, help = "Skip reads with 'N' in their CDR3 sequences.")
        parser.add_argument("--cluster_by_gene", action='store_true',
                            default=False, help = "Use gene instead of allele in clonal clustering" )

        # parser.add_argument("--norm", type = str,
        #                     choices = ['len', 'mut', 'none'], default = 'len',
        #                     help = "Specifies how to normalize distances. One of none (do not normalize), len (normalize \
        #                     by length), or mut (normalize by number of mutations between sequences). (default: len)")
        # parser.add_argument("--sym", type = str, choices = ['avg', 'min'], default = 'avg',
        #                     help = "Specifies how to combine asymmetric distances. One of \
        #                     bavg (average of A->B and B->A) or min (minimum of A->B and B->A). (default: avg)")
        # parser.add_argument("--link", type = str,  choices = ['single', 'average', 'complete'], default = 'single',
        #                     help = "Type of linkage to use for hierarchical clustering. (default: single)")

        # parser.add_argument("--treeRead", type = int, default = 0,
        #                     help = "The minimum read number (>=3) for each clonal to construct lineage tree (0)")
        # parser.add_argument("--treeSample", type = int, default = 1,
        #                     help = "The minimum shared samples for clonal to construct lineage tree (1)")
        # parser.add_argument("--treeOnly", type = str, choices = ['T', 'F'],
        #                     default = 'F', help = "Only process tree construction")
        #
        # parser.add_argument("--shortname", type = str, default = 'S', help = "Short name for samples in construct lineage tree. \
        #                     Specify one generic or all seperated by comma, eg. S or H1,H2,L1,L2 (S)")
        # parser.add_argument("--model", type = str, choices = ['aa','ham','m1n','hs1f','hs5f'], default = 'm1n',
        #                         help = "Specifies which substitution model to use for calculating distance between sequences. Where m1n is \
        #                         the mouse single nucleotide transition/trasversion model of Smith et al, 1996; hs1f is the human single \
        #                         nucleotide model derived from Yaari et al, 2013; hs5f is the human S5F model of Yaari et al, 2013; ham is \
        #                         nucleotide Hamming distance; and aa is amino acid Hamming distance. The hs5f data should be considered \
        #                         experimental. (default: m1n)")

    return parser

def parse_args():
    '''parse arguments for different sub commands'''
    # Define ArgumentParser
    parser = ArgumentParser(description="HTGTS repertoire sequencing analysis pipeline")
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' Version:%s %s' %(__version__, __date__))

    subparser = parser.add_subparsers(title='subcommands', dest='subcmd')

    # preprocess paired fastq
    subparser.add_parser('preprocess',
        parents=[parentParser(subname='preprocess')],
        help='Preprocess sequences files.')

    # Preprocess and VDj assignment
    subparser.add_parser('run',
        parents=[parentParser(subname='preprocess'), parentParser(subname='run')],
        help='Run preprocess and VDJ assignment from fastq/fasta files.')

    # mutProfile: generate mutation profile along VDJ genes
    subparser.add_parser('mut',
        parents=[parentParser(subname='mut')],
        help='Generate DNA and protein mutation profile along V genes. \
        Required: log file in VDJ assignment analysis')

    # clonal: Clonal analysis and phylogenetic tree construction
    subparser.add_parser('clonal',
        parents=[parentParser(subname='clonal')],
        help='Clonal clustering and phylogenetic tree construction. \
        Required: log file in VDJ assignment analysis')

    args = parser.parse_args()

    if args.subcmd == None:
        parser.print_help()
        sys.exit(0)

    check_args(args)
    return args
