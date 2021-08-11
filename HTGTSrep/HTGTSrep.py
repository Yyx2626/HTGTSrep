#!/usr/bin/env python3

'''
HTGTS repertoire sequencing data analysis pipeline
'''

import os, sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from HTGTSrep.argparser import parse_args
from HTGTSrep.preprocess import reads_process
from HTGTSrep.igblast import run_IgBlast, parse_IgBlast

from HTGTSrep.lib import files_process
from HTGTSrep.mutprofile import mutProfile
from HTGTSrep.clonal import clonal_main

def main():
    args = parse_args()

    if args.subcmd == 'preprocess':
        reads_process(args)

    if args.subcmd == 'run':
        if not args.skipDemultiplex and not args.skipIgBlast:
            reads_process(args)
        if not args.skipIgBlast:
            run_IgBlast(args)
        parse_IgBlast(args)
        # clean up some files from IgBlast
        files_process(args, 'igblast_clean')

    if args.subcmd == 'mut':
        mutProfile(args)

    if args.subcmd == 'clonal':
        clonal_main(args)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) Bye!\n")
        sys.exit(0)
