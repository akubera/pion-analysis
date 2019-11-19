#
# plotting/levys.py
#


import numpy as np
from . import PlotData
from random import random
from .gauss import series_to_TGraphErrors, merged_dataframe, AddRoutRsideRatio


def plot_levy(df, cfg_filter=None, c=None):
    from femtofitter import FitResults

    if c is None:
        from ROOT import TCanvas
        c = TCanvas()

    if isinstance(df, FitResults):
        df = df.df

    if cfg_filter:
        df = df[df.cfg==cfg_filter]

    for cent, centdf in df.groupby('cent'):
        for kt, ktdf in centdf.groupby('kT'):
            print(kt)
            print(ktdf)

    plot = PlotData(c)

    return plot


def plot_levy3d_df(df,
                   merge=True,
                   cfg_filter=None,
                   c=None,
                   palette='colorblind',
                   marker_size=1.2,
                   marker_style=21,
                   shift=0.03):
    import ROOT
    from ROOT import TCanvas
    import seaborn as sns

    if cfg_filter:
        df = df[df.cfg==cfg_filter].copy()
    else:
        df = df.copy()

    df = AddRoutRsideRatio(df)

    if merge:
        df = merged_dataframe(df, ['Ro', 'Rs', 'Rl', 'alpha', 'lam', 'RoRs'])

    if c is None:
        c = TCanvas()

    plot = PlotData(c)
    cols = 3
    rows = 2
    c.Divide(cols, rows)
    c.SetCanvasSize(1600, 1000)

    tcolors = plot.tcolors = [ROOT.TColor(*rgb)
                              for rgb in sns.color_palette(palette)]
    colors = [tcolor.GetNumber() for tcolor in tcolors]
    color_it = iter(colors)

    centrality_colors = {}

    def centrality_color(cent):
        if cent not in centrality_colors:
            centrality_colors[cent] = next(color_it)
        return centrality_colors[cent]

    plot.axis_hists = []
    plot.graphs = []
    plot.lines = []

    legend = plot.legend = ROOT.TLegend(0.01, 0.7, 0.9, 0.9)
    legend.SetHeader('  #scale[1.2]{#bf{Centrality}}', '')
    legend.SetLineWidth(0)

    legend.SetNColumns(3)
    c.cd(6)
    plot.legend.Draw()

    KT_RANGE = 0.17, 1.04
    TICKCODE_X = 408

    for i in range(1, cols*rows+1):
        pad = c.cd(i)
        pad.SetLeftMargin(0.055)
        pad.SetRightMargin(0.03)
        pad.SetBottomMargin(0.1)
        pad.SetTopMargin(0.1)
        pad.SetTickx(1)
        pad.SetTicky(1)

    YRNG = 1.1, 7.9
    hist_info = [
        ('R_{out}', YRNG),
        ('R_{side}', YRNG),
        ('R_{long}', YRNG),
        ('#lambda', (0.0, 0.8)),
        ('#alpha', (1.0, 3.0)),
        ('R_{out} / R_{side}', (0.5, 1.6))
    ]

    for title, y_range in hist_info:
        axh = ROOT.TH1C("h%s" % random(), f"{title}; k_{{T}} (GeV)", 1000, 0.0, 1.5)
        axh.SetTitleOffset(0.0)
        xax, yax = axh.GetXaxis(), axh.GetYaxis()
        xax.SetTitleOffset(0.89)
        xax.SetTitleSize(0.04)
        xax.SetNdivisions(TICKCODE_X)
        xax.SetRangeUser(*KT_RANGE)
        yax.SetRangeUser(*y_range)
        plot.axis_hists.append(axh)

    shifts = shift * np.linspace(-1.0, 1.0, len(df.cent.unique()))

    for i, key in enumerate(['Ro', 'Rs', 'Rl', 'lam', 'alpha'], 1):
        pad = c.cd(i)

        plot.axis_hists[i-1].Draw()

        for cidx, (cent, cdf) in enumerate(df.groupby('cent')):

            graph = series_to_TGraphErrors(cdf, key)
            if shift:
                x = np.frombuffer(graph.GetX(), np.float64, graph.GetN())
                x += shifts[cidx]

            graph.SetMarkerColor(centrality_color(cent))
            graph.SetLineColor(centrality_color(cent))
            graph.SetMarkerStyle(marker_style)
            graph.SetMarkerSize(marker_size)
            graph.Draw("SAME P")

            plot.graphs.append(graph)

            if i == 1:
                legend.AddEntry(graph, f"{cent.replace('_', '-')}%", 'P')

    # pad = c.cd(5)
    # plot.axis_hists[-1].Draw()

    # for cidx, (cent, cdf) in enumerate(df.groupby('cent')):
    #     graph = series_to_TGraphErrors(cdf, 'RoRs')
    #     plot.graphs.append(graph)

    #     graph.SetMarkerColor(centrality_color(cent))
    #     graph.SetLineColor(centrality_color(cent))
    #     graph.SetMarkerStyle(21)
    #     graph.SetMarkerSize(0.8)

    #     graph.Draw("SAME P")

    # plot.lines.append(ROOT.TLine(KT_RANGE[0], 1.0, KT_RANGE[1], 1.0))
    # plot.lines[-1].SetLineStyle(2)
    # plot.lines[-1].Draw()

    return plot



def plott(df):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 3, figsize=(18, 12))

    import seaborn as sns
    current_palette = sns.color_palette()

    for color, (cent, centdf) in zip(current_palette, df.groupby('cent')):

        for i, key in enumerate(['Ro', 'Rs', 'Rl']):
            centdf.plot.scatter('kT', key, yerr=key + "_err", ax=ax[0,i], color=color)

        for i, key in enumerate(['alphao', 'alphas', 'alphal']):
            centdf.plot.scatter('kT', key, ax=ax[1,i], color=color)


def plot_levy_1d(df):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(18, 7))

    import seaborn as sns
    current_palette = sns.color_palette()

    for color, (cent, centdf) in zip(current_palette, df.groupby('cent')):

        for i, key in enumerate(['radius', 'lam', 'alpha']):
            centdf.plot.scatter('kT', key, yerr=key + "_err", ax=ax[i], color=color)


def plot_levy_projections(tfile, series, c=None, ignore_mrc=False):
    import re
    import matplotlib.pyplot as plt
    from femtofitter import PathQuery

    pq = PathQuery.From(series)
    tdir = tfile.Get(pq.as_path())

    if c is None:
        from ROOT import TCanvas
        c = TCanvas()

    import ROOT
    from ROOT import TH1D

    plot = PlotData(c)
    num, den = map(tdir.Get, ("num", "den"))

    num.Divide(den)
    num.SetStats(0)

    plot.r = num

    Xlim = 0.2

    cf = TH1D()
    num.Copy(cf)

    num.GetXaxis().SetRangeUser(0, Xlim)

    FitterClassname = 'Fitter1DLevyPolyBg'
    FitterClass = getattr(ROOT, FitterClassname)

    params = FitterClass.FitResult(series).as_params()
    print(params)

    params.Normalize(num)

    plot.line_at_one = ROOT.TLine(0.0,1.0, Xlim, 1.0)
    plot.line_at_one.SetLineStyle(2)
    plot.line_at_one.Draw("SAME")

    plot.line_at_one.SetLineStyle(2)
    plot.line_at_one.Draw("SAME")

    fitter = FitterClass(tdir, series.limit)

    cf.SetName(f"FitResult:{series.name}")
    cf.Reset()

    if (fsi := series.fsi):
        print(fsi)
        m = re.match(f"(?P<classname>\w+)\[(?P<filename>[^\]]+)\]", fsi)
        FsiClass = getattr(ROOT, m.group('classname'))
        fitter.fsi = FsiClass.From(m.group('filename'))

    if (mrc := series.mrc) and not ignore_mrc:
        m = re.match(f"(?P<classname>\w+)\[(?P<filename>[^:]+):(?P<path>[^\]]+)\]", mrc)
        MrcClass = getattr(ROOT, m.group('classname'))
        mrc_file = ROOT.TFile.Open(m.group('filename'))
        mrc_path = mrc_file.Get(m.group('path'))
        fitter.mrc = MrcClass.From(mrc_path)
        fitter.fill_smeared_fit(cf, params)
    else:
        fitter.fill(cf, params)

    cf.SetLineColor(ROOT.kRed)
    plot.cf = cf
    cf.Scale(1.0/cf.GetBinContent(cf.GetXaxis().FindBin(Xlim)))

    num.Draw()
    cf.Draw('L SAME')

    return plot


def foo():
    c = TCanvas()
    c.Divide(2, 2, 0, 0)
    c.SetCanvasSize(1000, 800)

    leg = TLegend(.65, .90, 1, 1)
    leg.SetHeader("Centrality", "C")
    leg.SetNColumns(3)
    leg.SetTextSize(0.03)

    for i, key in enumerate(('Ro', 'Rs', 'Rl', 'RoRs'), 1):
        pad = c.cd(i)
    #     pad.SetTopMargin(0)
    #     pad.SetBottomMargin(0)
    #     pad.SetLeftMargin(0.17)
        if i == 4 or i == 3:
            pad.SetBottomMargin(0.18)
    #     if i in (1, 3):
        pad.SetLeftMargin(0.17)
    #     elif i == 1:
    #         pad.SetTopMargin(0.0)

        if i == 3:
            yrange = 0.1, 9.5
        elif i == 4:
            yrange = 0.7, 1.3
        else:
            yrange = 0.1, 7.9

        for g, graph in enumerate(kisiel_graphs[key]):
            graph.GetXaxis().SetRangeUser(0.0, 1.0)
            if g == 0:
                graph.GetYaxis().SetRangeUser(*yrange)
                graph.GetYaxis().SetTitleSize(.07)
                graph.GetYaxis().SetTitleOffset(0.69)
                graph.GetYaxis().CenterTitle()
                if i == 4:
                    graph.GetYaxis().SetNdivisions(404)
                elif i == 3:
                    graph.GetYaxis().SetNdivisions(405)
                if i == 4:
                    graph.GetXaxis().SetTitleSize(0.07)
                    graph.GetXaxis().SetLabelSize(0.07)
    #                 graph.Set

                graph.Draw('AP')
            else:
                graph.Draw('SAME P')

            if i == 1:
                leg.AddEntry(graph, cent_names[g], 'P')

        for g, graph in enumerate(kubera_graphs[key]):
            graph.Draw('P')

    c.cd(0)
    leg.Draw()
    c.Draw()
