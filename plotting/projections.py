#
# pion-analysis/plotting/projections.py
#

from typing import Any, Tuple, Union


from . import PlotData
import re
from functools import lru_cache


def load_datafile(fr):
    from pathlib import Path
    from ROOT import TFile

    path = Path(fr.filename)
    if not path.exists():
        path = 'femtofitter' / path

    if not path.exists():
        raise FileNotFoundError()

    tfile = TFile.Open(str(path))
    assert tfile
    return tfile


def deserialize_args(args: str) -> Tuple:
    """
    Turn argument string into tuple to pass to the thing
    """

    from ast import literal_eval

    if args == '':
        return ()

    res = []
    try:
        args = literal_eval(args)
    except ValueError:
        for arg in args.split(','):
            try:
                a = literal_eval(arg)
            except ValueError:
                a = arg
            res.append(a)
    else:
        if isinstance(args, (tuple, list)):
            res = args
        else:
            res = (args, )

    return tuple(res)


@lru_cache()
def build_fsi_obj(fsi_str: str, shared_ptr: bool=True):
    m = re.match(r'(?P<cls>\w+)\[(?P<args>[^\]]+)\]', fsi_str)
    if not m:
        raise ValueError(f"Could not parse fsi string: {fsi_str!r}")
    import ROOT
    FsiClass = getattr(ROOT, m.group('cls'))
    fsi_args = deserialize_args(m.group('args'))
    if shared_ptr:
        return FsiClass.From(*fsi_args)
    else:
        return FsiClass(*fsi_args)


@lru_cache()
def load_fsi(fitter, fsi_str: str) -> None:
    fitter.fsi = build_fsi_obj(fsi_str)


@lru_cache()
def load_mrc(fitter, mrc_str):
    from pathlib import Path
    import ROOT
    from ROOT import TFile
    regex = re.compile(r"(?P<classname>\w+)\[(?P<filename>[^:]+):(?P<path>[^\]]+)\]")

    if not mrc_str:
        return

    m = regex.match(mrc_str)
    if not m:
        raise ValueError(f"Invalid MRC string {mrc_str}")

    MrcClass = getattr(ROOT, m.group('classname'))
    if not MrcClass:
        raise ValueError(f"Invalid MRC class {m.group('classname')}")

    mrc_filepath_orig = mrc_filepath = Path(m.group('filename'))
    if not mrc_filepath.exists():
        mrc_filepath = 'femtofitter' / mrc_filepath

    if not mrc_filepath.exists():
        raise ValueError(f"Could not find MRC file {mrc_filepath_orig}")

    mrc_file = TFile.Open(str(mrc_filepath))
    mrc_path = mrc_file.Get(m.group('path'))
    fitter.mrc = MrcClass.From(mrc_path)


def plot_fit_projections(fitresult,
                         centrality,
                         cfg=None,
                         kt=None,
                         pair='pip',
                         fclassname=None,
                         ignore_mrc=False,
                         pad=None,
                         canvas_size=(700, 1200),
                         integration=0.05):
    import re
    import ROOT
    from ROOT import TH1, TH3, TCanvas
    from ROOT import TFile, gStyle, gROOT
    from ROOT import TGraphErrors, TLegend, Data3D
    from itertools import product
    from femtofitter import PathQuery

    if not fclassname:
        if 'alpha' in fitresult.df.columns:
            fclassname = 'Fitter3DLevyLcms'
        else:
            fclassname = 'Fitter3DGaussLcms'

    FitterClass = getattr(ROOT, fclassname)

    df = fitresult.df
    df = df[df.cent==centrality]
    df = df[df.pair==pair]

    if df.empty:
        raise ValueError(f"No (centrality, pair) ({centrality!r}, {pair!r}) in fit result")

    if cfg is None:
        cfg = fitresult.config.iloc[0].name
    elif isinstance(cfg, int):
        cfg = fitresult.config.iloc[cfg].name

    df = df[df.cfg==cfg]
    if df.empty:
        raise ValueError(f"No cfg {cfg!r}")

    if isinstance(kt, str):
        kt = (kt, )

    if kt:
        df = df[df.kt.isin(kt)]

    if df.empty:
        raise ValueError(f"No kt in {kt!r}")

    tfile = load_datafile(fitresult)

    if pad is None:
        pad = TCanvas()
        pad.SetCanvasSize(*canvas_size)

    rows = [row for _, row in df.iterrows()]

    plot = PlotData(pad)
    pad.Divide(1, len(rows), 0.01, 0.01)
    pad.SetFillColor(ROOT.kRed)
    int_range = (-integration, integration)

    proj = {
        'x': (TH3.ProjectionX, TH3.GetYaxis, TH3.GetZaxis),
        'y': (TH3.ProjectionY, TH3.GetXaxis, TH3.GetZaxis),
        'z': (TH3.ProjectionZ, TH3.GetXaxis, TH3.GetYaxis),
    }

    for j, series in enumerate(rows, 1):
        pq = PathQuery.From(series)
        tdir = tfile.Get(pq.as_path())

        data = Data3D.From(tdir, series.fitrange)
        fitter = FitterClass(data)

        fit_result = FitterClass.FitResult(series)
        fit_params = fit_result.as_params()

        names = "num", "den"
        num, den = map(tdir.Get, names)

        load_fsi(fitter, series.fsi)
        if series.mrc and not ignore_mrc:
            load_mrc(fitter, series.mrc)

        cf = fitter.get_cf(fit_params)
        if fitter.mrc:
            cf_den = fitter.mrc.GetSmearedDenLike(cf)
        else:
            cf_den = den

        cf.Multiply(cf_den)

        fit_params.Normalize(num)
        fit_params.Normalize(cf)

        def get_range(hist, get_ax):
            return [max(get_ax(hist).FindBin(x), 1)
                    for x in int_range]

        spad = pad.cd(j)
        spad.Divide(3, 1, 0, 0)
        for i, ax in enumerate('xyz'):
            subpad = spad.cd(i+1)

            project, getax1, getax2 = proj[ax]

            range1 = get_range(num, getax1)
            range2 = get_range(num, getax2)

            n = project(num, f"n{ax}", *range1, *range2)
            d = project(den, f"d{ax}", *range1, *range2)

            cf_p = project(cf, f"fit{ax}", *range1, *range2)
            cf_d = project(cf_den, f"fitden{ax}", *range1, *range2)

            ratio = n
            ratio.SetTitle("")
            if ratio.GetSumw2() == 0:
                ratio.Sumw2()
            ratio.Divide(d)
            ratio.SetStats(0)
            ratio.SetLineWidth(2)
            ratio.DrawCopy("HE")

            cf_p.SetLineColor(ROOT.kRed)
            cf_p.SetLineWidth(2)
            cf_p.SetTitle("")
            cf_p.Divide(cf_d)
            cf_p.SetStats(0)
            cf_p.DrawCopy("SAME HIST C")

    return plot

    cfg = centrality
    dfiloc = fitresult.df.iloc[index]

    graphs = plot.graphs = []
    graphs_lam = plot.graphs_lam = []

    leg = plot.legend = TLegend(0.68, 0.6, 0.88, 0.8)
    leg.SetHeader("Centrality", 'C')

    sns_colors = sns.color_palette("colorblind")
    tcolors = plot.tcolors = [ROOT.TColor(*c) for c in sns_colors]
    colors = plot.colors = [c.GetNumber() for c in tcolors]

    print(dfiloc)

    return plot
