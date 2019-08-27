#!/usr/bin/env python
#
# analyses/tpcda/run.py
#

import os
import sys
import asyncio
import tempfile
import subprocess as sp
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from tqdm import tqdm


def arg_parser():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    return parser


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    __dir__ = Path(__file__).parent

    os.chdir(__dir__)

    args = arg_parser().parse_args(argv)

    setup_analysis()
    #run_analysis('foobar.root', ["/alice/data/2015/LHC15o/000245683/pass1/AOD194/3000/AliAOD.root"])
    #return

    files = []
    for i in range(3000, 3040):
    # for i in range(3000, 3003):
        filepath = Path("/alice/data/2015/LHC15o/000245683/pass1/AOD194/%04d/AliAOD.root" % i)
        if not filepath.exists():
            continue
        files.append(filepath)

    loop = asyncio.get_event_loop()
    loop.run_until_complete(run_multiproc_analysis(files))


def setup_analysis(filename='analysis.root'):
    from ROOT import gROOT, TFile, AliAnalysisManager

    tfile = TFile.Open(str(filename), "RECREATE")

    mgr = AliAnalysisManager("mgr")

    gROOT.Macro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C")

    gROOT.Macro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C")
    gROOT.Macro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C")

    gROOT.LoadMacro("$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/macros/Train/PionPionFemto/AddNuTaskPionPionRoot6.C")
    gROOT.ProcessLine("""AddNuTaskPionPionRoot6(
        "container1",
        "",
        "{30:40:50}; (0.2:0.3:0.4:0.5:0.6:0.7:0.8:1.0);"
        "~do_kt_ylm_cf=true; ~do_kt_qinv_cf=true; ~do_kt_pqq3d_cf=true;"
        "~do_sharequality_cf=true;  ~do_avg_sep_cf=true; ~do_detadphistar_cf=true;"
        "~q3d_bin_count=47; ~q3d_maxq = 0.141;"
        "@enable_pair_monitors=false;"
        "@num_events_to_mix=3;"
        "$pion_1_min_tpc_chi_ndof=0.33; $pion_1_max_its_chi_ndof=1.8;"
        "$pion_1_max_tpc_chi_ndof=1.8;"
     )""")

    mgr.Write()
    tfile.Close()


async def run_multiproc_analysis(files, threads=None):
    from pprint import pprint
    from functools import partial
    from itertools import starmap as smap
    from asyncio import wait, FIRST_COMPLETED

    threads = threads or len(os.sched_getaffinity(0))

    executor = ProcessPoolExecutor()
    loop = asyncio.get_running_loop()
    parallel_run = partial(loop.run_in_executor, executor, run_analysis)

    file_list = [files[i::threads] for i in range(threads)]
    file_list = [f for f in file_list if f]
    pprint(file_list)

    with tempfile.TemporaryDirectory() as tdir:
        tdir = Path(tdir)
        print(tdir)

        ofiles = [tdir / ('AnalysisResults%d.root' % i)
                  for i in range(len(file_list))]

        # pending = [parallel_run(*p) for p in zip(ofiles, file_list)]
        pending = list(smap(parallel_run, zip(ofiles, file_list)))

        finished = []

        with tqdm(total=len(pending)) as pbar:

            while pending:
                done, pending = await wait(pending,
                                           loop=loop,
                                           return_when=FIRST_COMPLETED)
                finished.extend([r.result() for r in done])
                pbar.update(len(done))

        #proc = sp.run(["hadd", '-j', str(threads), '-f', 'Results.root', *map(str, ofiles)],
        proc = sp.run(["hadd", '-f', 'Results.root', *map(str, ofiles)],
                      check=True)
        print(proc, proc.returncode)

    return finished


def run_analysis(output_name, files):
    import os
    from pathlib import Path

    # redirect stdout
    sys.stdout.flush()
    newstdout = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 1)
    os.close(devnull)
    sys.stdout = os.fdopen(newstdout, 'w')

    from ROOT import TFile, TChain
    input_files = TChain("aodTree")

    for file in map(str, files):
        input_files.Add(file)

    tfile = TFile.Open("analysis.root")
    mgr = tfile.Get("mgr")

    os.chdir(Path(output_name).parent)
    output_filename = output_name.name

    mgr.SetCommonFileName(output_filename)
    for output in mgr.GetOutputs():
        output.SetFileName(output_filename)

    mgr.InitAnalysis()
    mgr.StartAnalysis("local", input_files)


if __name__ == '__main__':
    exit(main())
