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
                        default='event000.root',
                        help='Root files')
    return parser


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = arg_parser().parse_args(argv)
    print(args)

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
    RunFemto(chain)



if __name__ == "__main__":
    exit(main())
