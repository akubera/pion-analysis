#
# pion-analysis/plotting/plot1d.py
#

import re
from . import PlotData
from femtofitter import PathQuery
from plotting.fitresults import MultiFitResults
from .projections import load_fsi, load_mrc
import numpy as np
import pandas as pd

_file_cache = {}


def merge_points(df):
    x = []
    y = []
    ye = []

    for kt, subdf in df.sort_values('kT').groupby('kT'):
        # x.append(subdf.kT.mean())
        x.append(subdf.kT.iloc[0])
        y.append(subdf.radius.mean())
        ye.append(subdf.radius_err.mean())

    data = pd.DataFrame([x,y,ye], ['kT', 'radius', 'radius_err']).T
    return data


def build_1D_fit_pair(series, tfile=None, ignore_mrc=False):
    """
    Build a fitter and fit results from
    """
    import ROOT
    from ROOT import TFile, TDirectory, Data1D

    fitter_classname = {
       "Gauss1D": "Fitter1DGauss",
       "Gauss1DPolyBg": "Fitter1DGaussPolyBg",
       "Levy1D": "Fitter1DLevy",
       "Levy1DPolyBg": "Fitter1DLevyPolyBg",
    }[series.fitter]

    FitterClass: type = getattr(ROOT, fitter_classname)

    if isinstance(tfile, str):
        if tfile not in _file_cache:
            _file_cache[tfile] = TFile.Open(tfile)
        tfile = _file_cache[tfile]

    if isinstance(tfile, TFile):
        pq = PathQuery.From(series)
        tdir = tfile.Get(pq.as_path())
    elif isinstance(tfile, TDirectory):
        tdir = tfile

    fit_result = FitterClass.FitResult(series)

    data = Data1D.From(tdir)
    fitter = FitterClass(data)

    load_fsi(fitter, series.fsi)
    if not ignore_mrc:
        load_mrc(fitter, series.mrc)

    return fitter, fit_result


# def load_mrc(mrc_paths):
#     import ROOT
#     from ROOT import TFile
#     mrc_ptrs = {}
#     for mrc in mrc_paths:
#         m = re.match(r'(?P<cls>\w+)\[(?P<file>[^:]+):(?P<path>[^\]]+)\]', mrc)
#         d = m.groupdict()
#         try:
#             mrc_file = _file_cache[d['file']]
#         except KeyError:
#             mrc_file = _file_cache[d['file']] = TFile.Open(d['file'])

#         mrc_cls = getattr(ROOT, d['cls'])
#         mrc_tdir = mrc_file.Get(d['path'])
#         mrc_ptrs[mrc] = mrc_cls.From(mrc_tdir)
#     return mrc_ptrs


# def load_fsi(fsi_paths):
#     import ROOT
#     from ROOT import TFile
#     if isinstance(fsi_paths, str):
#         fsi_paths = (fsi_paths, )

#     fsi_objs = {}
#     for fsi in fsi_paths:
#         m = re.match(r'(?P<cls>\w+)\[(?P<file>[^\]]+)\]', fsi)
#         d = m.groupdict()
#         try:
#             fsi_file = _file_cache[d['file']]
#         except KeyError:
#             fsi_file = _file_cache[d['file']] = TFile.Open(d['file'])

#         fsi_cls = getattr(ROOT, d['cls'])
#         fsi_objs[fsi] = fsi_cls.From(fsi_file)

#     return fsi_objs


def make_1d_plots(df, c=None):
    if c is None:
        from ROOT import TCanvas
        c = TCanvas()
        c.Divide(1, 2)
    import seaborn as sns
    import ROOT
    from ROOT import TGraphErrors, TLegend
    plot = PlotData()
    plot.canvas = c

    graphs = plot.graphs = []
    graphs_lam = plot.graphs_lam = []

    leg = plot.legend = TLegend(0.68, 0.6, 0.88, 0.8)
    leg.SetHeader("Centrality", 'C')

#     colors = [ROOT.kRed, ROOT.kBlue, ROOT.k]
    sns_colors = sns.color_palette("colorblind")
    tcolors = plot.tcolors = [ROOT.TColor(ROOT.TColor.GetFreeColorIndex(), *c) for c in sns_colors]
    colors = plot.colors = [c.GetNumber() for c in tcolors]


    for i, (cent, cdf) in enumerate(df.groupby('cent')):
        color = colors[i]
        data = merge_points(cdf)

        g_rinv = TGraphErrors(data.shape[0])
        np.frombuffer(g_rinv.GetX())[:] = np.array(data.kT)
        np.frombuffer(g_rinv.GetY())[:] = np.array(data.radius)
        np.frombuffer(g_rinv.GetEY())[:] =  np.array(data.radius_err) * 10

        g_lam = TGraphErrors(data.shape[0])
        np.frombuffer(g_lam.GetX())[:] = np.array(data.kT)
        np.frombuffer(g_lam.GetY())[:] = np.array(data.lam)
        np.frombuffer(g_lam.GetEY())[:] =  np.array(data.lam_err) * 10


        for g in (g_rinv, g_lam):
            g.SetMarkerStyle(21)
            g.SetMarkerSize(0.7)
            g.SetMarkerColor(color)
            g.SetLineColor(color)

        graphs.append(g_rinv)
        graphs_lam.append(g_lam)

        c.cd(1)
        if i == 0:
            g_rinv.Draw("APE")
        else:
            g_rinv.Draw("P ")

        c.cd(2)
        if i == 0:
            g_lam.Draw("APE")
        else:
            g_lam.Draw("P ")

        cent_name = cent.replace('_', '-')  + "%"
        leg.AddEntry(g_rinv, cent_name, 'P')


#     for g in graphs:
#         g.Draw('AP')
#     g.Draw('APL')

#     g.SetTitle("Radius")
    leg.Draw()
    return plot


def load_joeys_style():
    from ROOT import gStyle, gROOT
    font_id = 12
    gStyle.SetOptStat(0)
    gStyle.SetTitleFont(font_id, "X") # LaTeX typeface
    gStyle.SetTitleFont(font_id, "Y")
    gStyle.SetTitleFont(font_id, "Z")
    gStyle.SetTitleFont(font_id, "T")
    gStyle.SetLabelFont(font_id, "X") # LaTeX typeface
    gStyle.SetLabelFont(font_id, "Y")
    gStyle.SetLabelFont(font_id, "Z")
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetPadTopMargin(0.13)
    gStyle.SetPadBottomMargin(0.25)
    gStyle.SetPadLeftMargin(0.17)
    gStyle.SetPadRightMargin(0.25)
    gStyle.SetCanvasDefW(1500)
    gStyle.SetCanvasDefH(900)
    gStyle.SetTitleSize(0.06,"X")
    gStyle.SetTitleSize(0.06,"Y")
    gStyle.SetTitleSize(0.06,"Z")
    gStyle.SetTitleSize(0.06,"T") # for title
    gStyle.SetLabelSize(0.04,"X")
    gStyle.SetLabelSize(0.04,"Y")
    gStyle.SetLabelSize(0.04,"Z")
    gStyle.SetTitleOffset(1.5,"X")
    gStyle.SetTitleOffset(1.1,"Y")
    gStyle.SetTitleOffset(1.3,"Z")
    gROOT.ForceStyle()


def make_1d_correlation_function(df, key, tfile=None, size=(800, 600), c=None):
    import ROOT
    from ROOT import TCanvas, TLine
    if c is None:
        c = TCanvas()
        c.SetCanvasSize(*size)

    if tfile is None:
        raise NotImplementedError

    if isinstance(key, int):
        series = df.iloc[key]
    else:
        raise NotImplementedError

    q = PathQuery.From(series)
    tdir = tfile.Get(q.as_path())
    num, den = map(tdir.Get, ("num", "den"))

    norm_rng = [num.GetXaxis().FindBin(l) for l in (0.25, 0.3)]

    ratio = num.Clone("ratio")
    if ratio.GetSumw2N() == 0:
        ratio.Sumw2()
    ratio.Divide(num, den, den.Integral(*norm_rng), num.Integral(*norm_rng))
    ratio.SetStats(0)

    ratio.SetTitle("Correlation Function; q_{inv} (GeV); CF(q_{inv})")
    ratio.SetMarkerStyle(8)
    ratio.SetMarkerSize(0.6)

    xax, yax = ratio.GetXaxis(), ratio.GetYaxis()
    xax.SetRangeUser(0.0, 0.3)
    xax.SetTitleSize(0.040)
    yax.SetTitleOffset(1.0)
    yax.SetTitleOffset(1.0)
    yax.CenterTitle()

    plot = PlotData()
    plot.ratio = ratio
    ratio.Draw('P')
    plot.canvas = c

    plot.unity_line = TLine(0, 1.0, 0.3, 1.0)
    plot.unity_line.Draw()
    plot.unity_line.SetLineStyle(2)
    plot.unity_line.SetLineColor(ROOT.kGray+2)

    # draw fit
    if False:
        fitter, fit_result = build_fitter_pair(series)
        print(fitter.FitResult(series))

    return plot


def plot_ratio_data(df1, df2, data, c=None):
    from ROOT import TH2D
    if c is None:
        from ROOT import TCanvas
        c = TCanvas()

    plot = PlotData()
    plot.canvas = c
    nktbins = 5
    ncentbins = 2
    h = TH2D('hist', 'Hist; kT (GeV); Centrality',
             nktbins, -0.5, nktbins+0.5,
             ncentbins, -0.5, ncentbins+0.5)
    h.SetStats(0)

    cents = []
    kts = []
    ratios = []

    for cent, cdf in df1.groupby('cent'):
        cents.append(cent)
        for kt, ktdf in cdf.groupby('kt'):
            pass

#     h.FillRandom()
#     h2 = h.Clone('h2')
#     h2.Fill

#     h.Divide(h2)
    np.frombuffer(h.GetArray(), count=h.GetNcells())[:] = data.flatten()
    h.SetBinContent(0,0)
    h.Draw('COLZ')

    plot.h = h
    return plot


def plot_projected_fits_abc(fr):
    pass


def plot_projected_fits(fr,
                        centrality,
                        N=4,
                        size=(600, 800),
                        ylim=(0.85, 1.44),
                        c=None,
                        I=1,
                        file_cache={}):
    """
    Plot big matrix of historams and fits
    """
    import ROOT
    from femtofitter import PathQuery
    from ROOT import TFile, gStyle, gROOT
    from ROOT import Data3D, Fitter3DGaussLcms as Fitter

    if c is None:
        from ROOT import TCanvas
        c = TCanvas()
        c.SetCanvasSize(*size)

    plot = PlotData(c)

    tfile = TFile.Open(fr.filename)
    assert tfile

    integration_limits = None

    def make_projections(num, den):
        nonlocal integration_limits
        xax = num.GetXaxis()
        yax = num.GetYaxis()
        zax = num.GetZaxis()

        if not integration_limits:
            xzero = xax.FindBin(0.0)
            yzero = yax.FindBin(0.0)
            zzero = zax.FindBin(0.0)

            xrange = (xzero - N, xzero + N) if xzero > 1 else (1, N)
            yrange = (yzero - N, yzero + N)
            zrange = (zzero - N, zzero + N)

            integration_limits = [
                (xax.GetBinLowEdge(xrange[0]), xax.GetBinUpEdge(xrange[1])),
                (yax.GetBinLowEdge(yrange[0]), yax.GetBinUpEdge(yrange[1])),
                (zax.GetBinLowEdge(zrange[0]), zax.GetBinUpEdge(zrange[1]))]
            # print(integration_limits)
        else:
            xrange, yrange, zrange = integration_limits
            xrange = tuple(map(xax.FindBin, xrange))
            yrange = tuple(map(yax.FindBin, yrange))
            zrange = tuple(map(zax.FindBin, zrange))

        name = 'hist%08x' % (abs(hash(num)))

        rx = num.ProjectionX(f'ratiox{name}', *yrange, *zrange)
        ry = num.ProjectionY(f'ratioy{name}', *xrange, *zrange)
        rz = num.ProjectionZ(f'ratioz{name}', *xrange, *yrange)

        for r in (rx, ry, rz):
            r.SetStats(0)
            r.Sumw2()
            r.SetTitle('')

        dx = den.ProjectionX(f'denx{name}', *yrange, *zrange)
        dy = den.ProjectionY(f'deny{name}', *xrange, *zrange)
        dz = den.ProjectionZ(f'denz{name}', *xrange, *yrange)

        rx.Divide(dx)
        ry.Divide(dy)
        rz.Divide(dz)

        return rx, ry, rz

    df = fr.df[fr.df.cent==centrality]

    mrc_dict = load_mrc(df.mrc.unique())
    fsi_dict = load_fsi(df.fsi.unique())

    hists = plot.hists = []
    fithists = plot.fithists = []
    fits = []
    kts = []

    plot.fit_tcolor = ROOT.TColor(1.0, 0.0, 0.0, 0.5)
    fit_color = plot.fit_tcolor.GetNumber()
    fsi = ROOT.FsiStatic(1.0)
    for kt, ktdf in df.groupby('kt'):
        kts.append(kt)

        series = ktdf.iloc[I]
        query = PathQuery.From(series)
        fitres = Fitter.FitResult(series)
        fits.append(fitres)

        tdir = tfile.Get(query.as_path())
        data = Data3D.From(tdir)
        fitter = Fitter(data)
        fitter.fsi = fsi_dict[series.fsi]
        fitter.mrc = mrc_dict[series.mrc]

        num, den = map(tdir.Get, ("num", "den"))

        rx, ry, rz = make_projections(num, den)
        hists.append((rx, ry, rz))
        qinv = data.src.qinv.get()
#         fsi = fitter.fsi.get()

        fit_params = fitres.as_params()
        fit_params.norm = 1.0

        fitden = fitter.mrc.GetSmearedDen()
#         fitnum = fitter.mrc.GetSmearedFit(fit_params, qinv, fsi)
#         fitnum.Multiply(fitden)

#         np.frombuffer(fitden.GetArray(), np.float32, fitden.GetNcells())[:] = 1.0
        fitnum = fitden.Clone('fitnum')
        fit_params.fill(fitnum, qinv, fsi, 1)

        fitnum.SetLineColor(fit_color)
        fitnum.SetMarkerColor(fit_color)

        fx, fy, fz = make_projections(fitnum, fitden)
        fithists.append((fx, fy, fz))
        for f in fithists[-1]:
            f.SetLineWidth(2)
            f.Scale(1.0/(f.GetBinContent(f.GetNbinsX()) or 1))

    c.Divide(3, len(hists), 0, 0)

    for j, (hlist, flist) in enumerate(zip(hists, fithists)):
        scale = 1.0 / fits[j].norm.value
        for i, (h, f) in enumerate(zip(hlist, flist), 1):
            pad = c.cd(3 * j + i)
            pad.SetTicks(1, 1)
            h.Draw("H")
            h.Scale(scale)
            h.GetXaxis().SetLabelSize(0.12)
            h.GetYaxis().SetRangeUser(*ylim)
            h.GetYaxis().SetLabelSize(0.12)
            h.GetYaxis().SetNdivisions(30205)
            h.GetXaxis().SetNdivisions(5)
            if i == 1:
                pad.SetLeftMargin(0.2)
            elif i == 3:
                pad.SetRightMargin(0.2)
                h.GetYaxis().SetTitle(str(kts[j]).replace('_', ' - '))
                h.GetYaxis().CenterTitle()
                h.GetYaxis().RotateTitle()
                title_size = 0.12
                title_offset = -0.53 / title_size
                h.GetYaxis().SetTitleSize(title_size)
                h.GetYaxis().SetTitleOffset(title_offset)

            if j > len(hists) - 3:
                h.GetXaxis().SetTitleSize(title_size)
#                 pad.SetBottomMargin(0.5)

            f.Draw("SAME HIST L")

    return plot


class PlotResults1D:
    """
    """

    def __init__(self, fr):
        if not isinstance(fr, MultiFitResults):
            fr = MultiFitResults(fr)
        self.fr = fr

    def plot(self,
             centrality,
             df=None,
             kts=None,
             pad=None,
             canvas_size=(500, 900),
             xmax=0.23,
             yrange=(0.8, 1.36),
             ignore_mrc=False,
             ):
        import ROOT
        from ROOT import TPad, TCanvas
        from cppyy.gbl import nullptr

        fitter_classmap = {
            'Gauss1DPolyBg': 'Fitter1DGaussPolyBg',
            'Levy1DPolyBg': 'Fitter1DLevyPolyBg',
        }

        if df is None:
            df = self.fr.df

        df = df[df.cent==centrality]

        if pad is None:
            pad = TCanvas()
            pad.SetCanvasSize(*canvas_size)

        plot = PlotData(pad)

        # pad.SetFillColor(ROOT.kRed)
        # pad.SetFillColor(ROOT.kBlack)

        if kts is None:
            kts = ['0.2_0.3', '0.3_0.4']
            kts = list(df.kt.unique())

        J = len(kts)

        # colors = [ROOT.kBlue, ROOT.kGreen, ROOT.kGreen+2]
        get_color = plot.color_loader('colorblind')

        TOP_BUFFER = 0.02
        BOTTOM_BUFFER = 0.05
        YSPACE = 0.00

        canvas_height_frac = (1.0 - TOP_BUFFER - BOTTOM_BUFFER - (J - 1) * YSPACE) / J

        plot.pads = []
        plot.hists = []
        plot.hists_fit = []
        plot.cf_data_hists = []
        plot.lines = []
        plot.texts = []
        plot.tfiles = {}

        filename = self.fr.frs[0].filename
        if filename not in plot.tfiles:
            plot.tfiles[filename] = ROOT.TFile.Open(filename)

        tfile = plot.tfiles[filename]

        for j, kt in enumerate(kts):
            pad.cd()

            yfrachi = 1.0 - TOP_BUFFER - j * (canvas_height_frac + YSPACE)
            yfraclo = yfrachi - canvas_height_frac

            # extend yfrac for the pad
            if kt == kts[-1]:
                yfraclo -= 0.008

            subpad = TPad(f'pad{j:02d}', '',
                          0.05, yfraclo, 0.96, yfrachi)
            subpad.Draw()
            # subpad.SetFillColor(get_color(kt))
            subpad.cd()
            subpad.SetTopMargin(0.0)
            subpad.SetBottomMargin(0.0)
            subpad.SetTickx(1)
            subpad.SetTicky(1)

            subdf = df[df.kt==kt]

            series = subdf.iloc[0]

            q = PathQuery.From(series)
            tdir = tfile.Get(q.as_path())

            fitter, fit_result = build_1D_fit_pair(series, tdir, ignore_mrc)

            num, den = map(tdir.Get, ("num", "den"))
            ratio = num.Clone('ratio')
            if ratio.GetSumw2N() == 0:
                ratio.Sumw2()
            ratio.Divide(num, den)

            fitcf = fitter.get_cf(fit_result)
            fitcf.SetName("fithist")
            fitcf.SetLineColor(ROOT.kRed)
            fitcf.SetLineWidth(1)

            fit_result.Normalize(ratio)
            fit_result.Normalize(fitcf)
            ratio.SetDirectory(nullptr)
            fitcf.SetDirectory(nullptr)

            # ratio.SetLineWidth(2)

            ratio.Draw('HE')
            fitcf.Draw("SAME HIST L")

            plot.hists.extend([num, den, ratio])
            plot.cf_data_hists.append(ratio)
            plot.hists_fit.append(fitcf)

            kt_lbl_str = 'k_{T} = %s (GeV)' % series.kt.replace('_', '-')
            kt_lbl = ROOT.TLatex(0.01, 1.33, kt_lbl_str)
            # kt_lbl = ROOT.TLatex(0.5, 0.10, kt_lbl_str)
            # kt_lbl.SetNDC()
            kt_lbl.SetTextSize(0.1)
            kt_lbl.SetTextFont(52)

            kt_lbl.Draw('same')
            plot.texts.append(kt_lbl)

            fit_bits = [
                f'R_{{inv}} = {series.radius:0.3f} #pm {series.radius_err:0.3f}',
                f'   #lambda = {series.lam:0.3f} #pm {series.lam_err:0.3f}',
            ]

            if 'alpha' in series:
                fit_bits.append(f'   #alpha = {series.alpha:0.3f}'
                                f' #pm {series.alpha_err:0.3f}')

            if len(fit_bits) == 2:
                fitlines = '#splitline{%s}{%s}' % tuple(fit_bits)
            elif len(fit_bits) == 3:
                fitlines = '#splitline{#splitline{%s}{%s}}{%s}' % tuple(fit_bits)

            fitlbl = ROOT.TLatex(xmax * 0.95, yrange[1] * 0.93, fitlines)
            fitlbl.SetTextAlign(33)
            fitlbl.SetTextSize(.1)
            fitlbl.SetTextFont(42)
            fitlbl.Draw()
            plot.texts.append(fitlbl)
            # subpad.SetLineColor(ROOT.kBlack)
            # subpad.SetBoxlimrderMode(5)
            # subpad.SetBorderSize(2)
            plot.pads.append(subpad)

            plot.hists = []

        for hist in plot.cf_data_hists:
            xax, yax = hist.GetXaxis(), hist.GetYaxis()
            xax.SetRangeUser(0.0, xmax)
            yax.SetRangeUser(*yrange)
            yax.SetLabelSize(0.08)
            yax.SetNdivisions(804)
            xax.SetNdivisions(808)
            xax.SetLabelSize(0.08)

        bottom_pad = plot.pads[-1]
        bottom_pad.SetBottomMargin(0.1)
        x = bottom_pad.FindObject("ratio").GetXaxis()
        x.SetTitle('q_{inv}')
        x.SetTitleSize(0.8)
        x.SetTitleOffset(0.0)

        # # print(bottom_pad.BBoxY2())
        # bbox = bottom_pad.GetBBox()
        # print('>', bbox.fX, bbox.fWidth, bbox.fY, bbox.fHeight)
        # # plot.pads[-1].SetBBoxY2(100)
        # bottom_pad.SetBBoxY2(bbox.fY + 10)
        # bbox = bottom_pad.GetBBox()
        # print('>', bbox.fX, bbox.fWidth, bbox.fY, bbox.fHeight)

        return plot

    def plot_results(self, df=None, **kwargs):
        if df is None:
            df = self.fr.df
        else:
            df = df

        if 'alpha' in df:
            return self.plot_levy_results(df=df, **kwargs)
        else:
            return self.plot_gauss_results(df=df, **kwargs)

    def plot_gauss_results(self,
                           df=None,
                           pad=None,
                           xshift=0.04,
                           palette='colorblind'):
        import ROOT
        from ROOT import TCanvas
        from .gauss import series_to_TGraphErrors
        from .gauss import build_tgraphs_from_data

        if df is None:
            df = self.fr.df.copy()
        else:
            df = df.copy()

        keys = ['radius', 'lam']
        df = self.fr.get_merged_dataframe(keys=keys, df=df)

        if pad is None:
            pad = TCanvas()
            pad.SetCanvasSize(800, 400)

        plot = PlotData(pad)

        pad.Divide(2, 1)

        get_cent_color = plot.color_loader(palette)

        tgraphs = build_tgraphs_from_data(df, keys=keys)
        plot.graphs = tgraphs

        plot.axhists = {
            'radius':  ROOT.TH1C("rinv_axes", "", 1000, 0.2, 1.2),
            'lam':  ROOT.TH1C("lam_axes", "", 1000, 0.2, 1.2),
        }
        plot.axhists['radius'].GetYaxis().SetRangeUser(1.4, 8.4)
        plot.axhists['lam'].GetYaxis().SetRangeUser(0.0, 0.6)

        tgraphs = build_tgraphs_from_data(df, keys=keys)
        plot.tgraphs = tgraphs
        plot.lines = []

        for i, key in enumerate(keys, 1):
            pad.cd(i)
            plot.axhists[key].Draw()

            x_shift = xshift
            for cent, cent_tgraph in tgraphs[key].items():
                color = get_cent_color(cent)

                tgraph, sys_tgraph = cent_tgraph
                tgraph.SetLineColor(color)
                tgraph.SetMarkerColor(color)
                tgraph.SetMarkerStyle(21)
                tgraph.SetMarkerSize(1)

                if sys_tgraph:
                    sys_tgraph.SetLineColor(color)
                    sys_tgraph.SetFillColorAlpha(color, 0.25)
                    sys_tgraph.Draw('SAME E5 P')

                tgraph.Draw('SAME P')

        return plot

    def plot_levy_results(self, df=None, pad=None):
        pass

    def categoryplot(self, key, df=None, kts=None, cents=None):
        import seaborn as sns
        df = df if df is not None else self.fr.df

        if kts:
            df = df[df.kt.isin(kts)]
        if cents:
            df = df[df.cent.isin(cents)]

        opts = dict(x='part:field',
                    col='kt',
                    row='cent',
                    hue='pair',
                   )
        gf = sns.catplot(data=df, y=key, **opts)
        return gf
