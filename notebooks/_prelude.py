#
# pion-analysis/notebooks/_prelude.py
#

import os
import sys
from os import environ
from pathlib import Path

# directory
_notebook_dir = Path(__file__).parent
_analysis_dir = _notebook_dir.parent

# locate FemtoFitter & insert into sys.path
_ff_path = environ.get('FEMTOFITTER_PATH', _analysis_dir / 'femtofitter')
_ff_path = Path(_ff_path).resolve()

_sys_paths = [Path(p).resolve() for p in sys.path]

if _ff_path.exists() and _ff_path not in _sys_paths:
    sys.path.insert(0, str(_ff_path))

# Add femtofitter build path into LD_LIBRARY_PATH
_libpath = environ.get("LD_LIBRARY_PATH", "").split(':')
_libpath.insert(0, str(_ff_path / 'build'))

environ['LD_LIBRARY_PATH'] = ':'.join(map(str, _libpath))

# Python libraries
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.notebook import tqdm

from itertools import chain, repeat, cycle, islice, product
from functools import partial, reduce
from copy import copy

import json, yaml, feather
import scipy
from scipy.interpolate import interp1d, CubicSpline
import statistics
from statistics import mean
from collections import defaultdict, Counter


try:
    import ROOT
except ImportError:
    import subprocess as _sp
    environ['MODULEPATH'] = '$HOME/alice/sw/MODULES/$ALIBUILD_ARCHITECTURE'
    alice_env = 'AliPhysics/latest-k614-user-root6-cpp14'
    _env_load = _sp.check_output(['modulecmd', 'python', 'load', alice_env],
                                 encoding='utf-8',
                                 stderr=_sp.DEVNULL)
    exec(_env_load)
    sys.path = environ['PYTHONPATH'].split(':') + sys.path
    import ROOT

import femtofitter
from femtofitter import FitResults, PathQuery

# ROOT classes
from ROOT import gROOT, gSystem, gInterpreter, gPad, gStyle, gDirectory
from ROOT import TCanvas, TPad
from ROOT import TLatex, TLine
from ROOT import TFile
from ROOT import TH1, TH3, TF1, TF2
from ROOT import TH1D, TH1F, TH2D, TH2F, TH3D, TH3F
from ROOT import TProfile, TProfile2D

TH1.AddDirectory(False)

# AliROOT Classes
from ROOT import AliFemtoConfigObject

# FemtoFitter Classes
from ROOT import (
    Data1D, Data3D,
    Fitter1DGauss, Fitter1DLevy,
    Fitter1DGaussPolyBg, Fitter1DLevyPolyBg,
    Fitter3DGaussLcms,
    # Fitter3DGaussLcmsFull,
    Fitter3DGaussLcmsOS, Fitter3DGaussLcmsOL,
    Fitter3DGaussLcms, Fitter3DGaussLcmsOS, Fitter3DGaussLcmsOL,

    Fitter3DLevy, Fitter3DLevyFull,

    Mrc1DRatio, Mrc1DRatioMixed,
    Mrc3DRatio, Mrc3DRatioMixed,

    FsiStatic, FsiGamov, FsiKFile,
)

# move out of notebook directory
os.chdir(_analysis_dir)

# load local plotting module
import plotting
from plotting.fitresults import MultiFitResults

# Setup ROOT Styles
gStyle.SetOptStat(0)
# gStyle.SetCanvasDefW(1500)
# gStyle.SetCanvasDefH(900)

# for _ax in 'XYZT':
#     gStyle.SetTitleSize(0.06, _ax)
#     gStyle.SetLabelSize(0.04, _ax)

# for _ax, _offset in (('X', 1.5), ('Y', 1.1), ('Z', 1.3)):
#     gStyle.SetTitleOffset(_offset, _ax)

# gStyle.SetPadTopMargin(0.13)
# gStyle.SetPadBottomMargin(0.25)
# gStyle.SetPadLeftMargin(0.17)
# gStyle.SetPadRightMargin(0.25)

gROOT.ForceStyle()
