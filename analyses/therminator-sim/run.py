#!/usr/bin/env python
#
# run.py
#

import sys


def arg_parser():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("file",
                        nargs='+',
                        help='Root files')
    parser.add_argument("-o", "--output",
                        help='output file')
    parser.add_argument("-l", "--limit",
                        type=int,
                        default=-1,
                        help='Limit event count')
    return parser


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = arg_parser().parse_args(argv)

    from ROOT import TChain, gROOT, gSystem, TH1

    gSystem.SetBuildDir(".rootbuild")
    gROOT.SetBatch(True)
    TH1.AddDirectory(False)

    chain = TChain("particles")
    for fname in args.file:
        chain.Add(fname)

    load_ok = gROOT.LoadMacro("RunFemto.C+")
    assert load_ok >= 0, load_ok

    from ROOT import RunFemto
    RunFemto(chain, args.output, args.limit)



if __name__ == "__main__":
    exit(main())
