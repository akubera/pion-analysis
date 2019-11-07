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
_ff_path = _analysis_dir / 'femtofitter'
_ff_path = environ.get('FEMTOFITTER_PATH', _ff_path)

if _ff_path.exists() and str(_ff_path) not in sys.path:
    sys.path.insert(0, str(_ff_path))

# Add femtofitter build path into LD_LIBRARY_PATH
_libpath = environ.get("LD_LIBRARY_PATH", "").split(':')
_libpath.insert(0, _ff_path / 'build')

environ['LD_LIBRARY_PATH'] = ':'.join(map(str, _libpath))

# Python libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from itertools import chain, repeat, cycle, islice, product
from functools import partial, reduce
from copy import copy

import json, yaml, feather
import scipy
from scipy.interpolate import interp1d, CubicSpline
import statistics
from statistics import mean

import femtofitter
from femtofitter import FitResults, PathQuery

# ROOT classes
import ROOT
from ROOT import gROOT, gSystem, gInterpreter, gPad, gStyle, gDirectory
from ROOT import TCanvas, TPad
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
