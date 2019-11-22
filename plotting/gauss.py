#
# plotting/gauss.py
#


import sys
from . import PlotData
from .systematics import calc_df_systematics
import numpy as np
import pandas as pd
from functools import lru_cache


@lru_cache
def load_kisiel_graphs(cent_names=None):
    import feather
    import ROOT
    from collections import defaultdict
    from itertools import repeat

    kisiel_df_filename = "/home/akubera/Physics/pion-analysis/Run1Data.feather"
    kisiel_df = feather.read_dataframe(kisiel_df_filename)

    colors = [ROOT.kBlue,
              ROOT.kRed,
              ROOT.kGreen+2,
              ROOT.kOrange+2,
              ROOT.kMagenta-3,
              ROOT.kCyan+3]
    markers = repeat(25)
    sizes = repeat(0.8)

    titles = {
        'Ro': 'R_{out} (fm)',
        'Rs': 'R_{side} (fm)',
        'Rl': 'R_{long} (fm)'
    }

    kisiel_graphs = defaultdict(list)
    if cent_names is None:
        cent_names = []

    centgroups = kisiel_df.groupby('cent')

    for key in ('Ro', 'Rs', 'Rl'):
        for (cent, cdf), color, marker, msize in zip(centgroups, colors, markers, sizes):
            graph = series_to_TGraphErrors(cdf, key, ekey=f'{key}_stat_err')
            graph.SetMarkerColor(color)
            graph.SetMarkerStyle(marker)
            graph.SetMarkerSize(msize)
            graph.SetTitle("; k_{T} (GeV/c); %s;" % titles[key])
            kisiel_graphs[key].append(graph)
            if key == 'Ro':
                cent_names.append(cdf.iloc[0].centname)

    for (cent, cdf), color, marker, msize in zip(centgroups, colors, markers, sizes):
        en = np.array(cdf["Ro_stat_err"])
        ed = np.array(cdf["Rs_stat_err"])

        os_ratio = np.array(cdf['Ro'] / cdf['Rs'])
        e = os_ratio * np.sqrt((en / cdf["Ro"])**2 + (ed / cdf["Rs"])**2)

        graph = build_TGraphErrors(cdf.kT, cdf['Ro'] / cdf['Rs'], e)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(msize)
        graph.SetTitle("; k_{T} (GeV/c); R_{out} / R_{side};")

        kisiel_graphs['RoRs'].append(graph)

    return dict(kisiel_graphs), cent_names


def build_TGraphErrors(x, y, ye=None, xe=None):
    from ROOT import TGraphErrors
    if xe is None:
        xe = np.zeros_like(x)
    else:
        xe = np.array(xe)

    if ye is None:
        ye = np.zeros_like(x)
    else:
        ye = np.array(ye)

    x = np.array(x)
    y = np.array(y)
    graph = TGraphErrors(x.size, x, y, xe, ye)
    return graph


def series_to_TGraphErrors(series, ykey, ekey=None, xkey='kT'):
    ekey = ekey or ykey + '_err'

    graph = build_TGraphErrors(series[xkey],
                               series[ykey],
                               series[ekey])
    return graph


def build_TH2(x, y, ye=None, xe=None):
    from ROOT import TH2D
    if xe is None:
        xe = np.zeros_like(x)
    else:
        xe = np.array(xe)

    if ye is None:
        ye = np.zeros_like(x)
    else:
        ye = np.array(ye)

    x = np.array(x)
    y = np.array(y)
    graph = TH2D(f"h{abs(hash(a.tobytes())):016X}", "", x.size, x, y, xe, ye)
    return graph


CENT_TO_COLOR_DICT = {}


def group_df_into_tgraphs(df, colors=None, markers=None, sizes=None, titles=None):
    import ROOT
    from collections import defaultdict
    global CENT_TOC_COLOR_DICT
    graphs = defaultdict(list)

    if colors is None:
        colors = [ROOT.kBlue,
                  ROOT.kRed,
                  # ROOT.kGreen+2,
                  ROOT.kOrange+2,
                  ROOT.kMagenta-3,
                  ROOT.kCyan+3]

    if markers is None:
        markers = [20] * len(colors)

    if sizes is None:
        sizes = [0.8] * len(colors)

    if titles is None:

        titles = {
            'Ro': 'R_{out} (fm)',
            'Rs': 'R_{side} (fm)',
            'Rl': 'R_{long} (fm)'
        }

    cent_dfs = df.groupby('cent')
    centcolor_needs_update = CENT_TO_COLOR_DICT == {}

    for (cent, cdf), color, marker, msize in zip(cent_dfs, colors, markers, sizes):

        if centcolor_needs_update:
            CENT_TO_COLOR_DICT[cent] = color
        else:
#             color = CENT_TO_COLOR_DICT.get(cent, color)
            # print(list(CENT_TO_COLOR_DICT.keys()))
#             color = CENT_TO_COLOR_DICT[cent]
#             color = color
            pass

        for key in ('Ro', 'Rs', 'Rl'):
            graph = series_to_TGraphErrors(cdf, key, key + "_err")
            graph.SetMarkerColor(color)
            graph.SetMarkerStyle(marker)
            graph.SetMarkerSize(msize)
            graph.SetTitle("; k_{T} (GeV/c); " + titles[key])
            graph.GetXaxis().SetTitleSize(0.06)
            graphs[key].append(graph)

        # Handle Ro / Rs
        en = np.array(cdf["Ro_err"])
        ed = np.array(cdf["Rs_err"])
        os_ratio = np.array(cdf['Ro'] / cdf['Rs'])
        e = os_ratio * np.sqrt((en / cdf["Ro"])**2 + (ed / cdf["Rs"])**2)

        graph = build_TGraphErrors(cdf.kT, os_ratio, e)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(msize)

        graphs['RoRs'].append(graph)

    return dict(graphs)


def plot_gauss3d_run1_comparison(df,
                                 c=None,
                                 merge=True,
                                 run1_graphdict=None,
                                 cfg_filter=None,
                                 palette='colorblind'):
    import ROOT
    from ROOT import TCanvas, TLegend

    if run1_graphdict is None:
        run1_graphdict, run1_cents = load_kisiel_graphs()

    if cfg_filter:
        df = df[df.cfg==cfg_filter].copy()
    else:
        df = df.copy()

    df = AddRoutRsideRatio(df)

    if merge:
        df = merged_dataframe(df, ['Ro', 'Rs', 'Rl', 'lam', 'RoRs'])

    if c is None:
        c = TCanvas()

    plot = PlotData(c)

    plot.legend = leg = TLegend()

    get_cent_color = plot.color_loader(palette)

    canvas_div_params = (1, 4, 0, 0)
    canvas_size = (500, 1200)

    canvas_div_params = (2, 2)
    canvas_size = (900, 900)

    c.Divide(*canvas_div_params)
    c.SetCanvasSize(*canvas_size)

    for i in range(4):
        pad = c.cd(i+1)
        pad.SetTickx(1)
        pad.SetTicky(1)

    keys = ['Ro', 'Rs', 'Rl', 'RoRs']
    title_map = {
        'Ro': 'R_{out}',
        'Rs': 'R_{side}',
        'Rl': 'R_{long}',
        'RoRs': 'R_{out} / R_{side}',
    }

    XRANGE = 0.18, 1.02
    YRANGE = 1.23, 9.3

    yranges = [
        YRANGE,
        YRANGE,
        YRANGE,
        (0.45, 1.65),
    ]

    titles = [title_map[key] for key in keys]

    plot.axhists = [ROOT.TH1C('h%d'%i, f'{title}; k_{{T}} GeV; {title}', 1000, *XRANGE)
                    for i, title in enumerate(titles)]

    c.SetFillColor(ROOT.kRed)

    plot.run2_graphs = []

    for i, hist in enumerate(plot.axhists):
        xax, yax = hist.GetXaxis(), hist.GetYaxis()
        yax.SetRangeUser(*yranges[i])
        yax.SetNdivisions(909)
        xax.SetNdivisions(509)

        xax.SetLabelSize(0.06)
        yax.SetLabelSize(0.065)
        xax.SetTitleOffset(1.1)

    for i, key in enumerate(keys, 1):
        pad = c.cd(i)
        hist = plot.axhists[i-1]
        hist.Draw()
        run1_graphs = run1_graphdict[key]

        for graph, cent in zip(run1_graphs, run1_cents):
            graph.SetMarkerColor(get_cent_color(cent))
            graph.SetMarkerSize(1.2)
            graph.Draw("SAME P")

        for cent, cdf in df.groupby('cent'):
            run2_graph = series_to_TGraphErrors(cdf, key)
            run2_graph.SetMarkerColor(get_cent_color(cent))
            run2_graph.SetMarkerSize(1.2)
            run2_graph.SetMarkerStyle(20)
            run2_graph.Draw("SAME P")
            plot.run2_graphs.append(run2_graph)

#         for graph in enumerate(run2_graphs):

#             pass

#         continue

#         print(run1_graphs)
#         continue
#         for g, graph in enumerate(run1_graphs):
#             graph.GetXaxis().SetRangeUser(0.0, 1.0)
#             if g == 0:
#                 graph.GetYaxis().SetRangeUser(*YRANGE)
#                 graph.GetYaxis().SetTitleOffset(0.69)
#                 graph.GetYaxis().CenterTitle()
#                 graph.GetYaxis().SetNdivisions(405)

#                 graph.GetXaxis().SetLabelSize(0.05)

#                 graph.Draw('AP')
#             else:
#                 graph.Draw('SAME P')

#             if i == 1:
#                 leg.AddEntry(graph, cent_names[g], 'P')

    plot.pads = []
    plot.graphs = []

    return plot


def merge_statistical_points(df, key):
    results = []

    for kt, data in df.groupby('kT'):
        errs = np.array(data[key + '_err'])
        weights = errs ** -2
        val = (data[key] * weights).sum() / weights.sum()
        err = np.sqrt(1.0 / weights.sum())
        results.append([kt, val, err])

    return np.array(results).T


def build_tgraphs_from_data(df, keys=None, skip_systematics=False, title_dict=None):
    from collections import defaultdict
    from ROOT import TGraphErrors

    if keys is None:
        keys = ("Ro", "Rs", "Rl")

    if title_dict is None:
        title_dict = {'Ro': "R_{out}",
                      'Rs': "R_{side}",
                      'Rl': "R_{long}",
                      'lam': "#lambda"}

    missing_titles = set(keys) - set(title_dict.keys())
    if missing_titles:
        raise ValueError(f"Title dict missing keys: {missing_titles}")

    def _merge_points(df, key):
        results = []

        for kt, data in df.groupby('kT'):
            errs = np.array(data[key + '_err'])
            weights = errs ** -2
            val = (data[key] * weights).sum() / weights.sum()
            err = np.sqrt(1.0 / weights.sum())
            results.append([kt, val, err])

        return np.array(results).T

    graphs = defaultdict(dict)

    for cent, cdf in df.sort_values('kT').groupby('cent'):

        for r in keys:
            x = []
            y = []
            ye = []
            x, y, ye = _merge_points(cdf, r)

            xe = np.zeros_like(x)
            title = [r]

            gdata = TGraphErrors(x.size)
            np.frombuffer(gdata.GetX(), dtype=np.float64)[:] = x
            np.frombuffer(gdata.GetY(), dtype=np.float64)[:] = y
            np.frombuffer(gdata.GetEY(), dtype=np.float64)[:] = ye

            if skip_systematics:
                gsys = None
            else:
                sys_key = f"{r}_sys_err"
                sys_err = np.array(cdf[sys_key])

                gsys = TGraphErrors(x.size)
                np.frombuffer(gsys.GetX(), dtype=np.float64)[:] = x - 0.004
                np.frombuffer(gsys.GetY(), dtype=np.float64)[:] = y
                np.frombuffer(gsys.GetEY(), dtype=np.float64)[:] = sys_err
                np.frombuffer(gsys.GetEX(), dtype=np.float64)[:] = 0.012

            graphs[r][cent] = (gdata, gsys)

    return dict(graphs)


class TheoryData:

    def __init__(self, datafile=None):
        from pathlib import Path
        import json

        if isinstance(datafile, str):
            datafile = Path(datafile)
        elif datafile is None:
            datafile = Path("Theory-data.json")
            if not datafile.exists():
                datafile = Path('notebooks/Theory-data.json')
            if not datafile.exists():
                datafile = Path('femtofitter/Theory-data.json')

        self.data = json.loads(datafile.read_text())
        self.saved_tlines = {}

    def get_tlines(self, key, centrality, temp):
        from toolz.itertoolz import sliding_window
        from ROOT import TLine
        try:
            return self.saved_tlines[(key, centrality, temp)]
        except KeyError:
            pass

        self.saved_tlines[(key, centrality, temp)] = tlines = []

        line_dict = self.data[key][centrality][temp]
        for p1, p2 in sliding_window(2, zip(*line_dict.values())):
            tlines.append(TLine(*p1, *p2))

#         line_dict = self.data[key][centrality][temp]
#             tlines.extend(TLine(*p1, *p2)
#                           for p1, p2 in sliding_window(2, zip(*line_dict.values())))

        return tlines

    def draw_tlines(self, key, centrality, temp, color, style, width=2):

        tlines = self.get_tlines(key, centrality, temp)

        for tline in tlines:
            tline.Draw()
            tline.SetLineWidth(width)
            tline.SetLineColor(color)
            tline.SetLineStyle(style)


def plot_gauss3d_theory_comparison(df,
                                   c=None,
                                   merge=True,
                                   theory_data=None,
                                   cfg_filter=None,
                                   x_range=(0.19, 1.21),
                                   y_range=(1.3, 8.52),
                                   palette='colorblind'):
    import ROOT
    from ROOT import TCanvas, TLegend, TLine, TH1C

    if c is None:
        c = TCanvas()
        c.Divide(2, 2)
        c.SetCanvasSize(800, 800)
    plot = PlotData(c)

    plot.hists = hists = [TH1C("h%i"%i, "", 1000, *x_range) for i in range(3)]
    for hist in hists:
        hist.GetYaxis().SetRangeUser(*y_range)

    if not isinstance(theory_data, TheoryData):
        theory_data = TheoryData(theory_data)
    plot.theory = theory_data

    centralities = list(plot.theory.data['Ro'].keys())
    assert centralities == ['00_05', '20_30', '40_50']

    df = df[df.cent.isin(centralities)]
    df = plot.df = merged_dataframe(df, theory_data.data.keys())

    linestyle_lotemp = 2
    linestyle_hitemp = 1

    get_cent_color = plot.color_loader(palette)
    get_sys_cent_color = plot.color_loader('muted', 'sys')
    for tcolor in plot.tcolor_dict['sys']:
        tcolor.SetAlpha(0.5)

    data_tgraphs = group_df_into_tgraphs(df)

#     plot.legend_cent = TLegend(0.5, 0.3, 0.745, 0.5)
    plot.legend_cent = TLegend(0.5, 0.4, 0.9, 0.5)
    plot.legend_cent.SetHeader('Centrality', '')
    plot.legend_cent.SetBorderSize(0)
    plot.legend_cent.SetNColumns(3)

#     plot.legend_temp = TLegend(0.75, 0.35, 0.95, 0.5)
    plot.legend_temp = TLegend(0.5, 0.35, 0.95, 0.20)
    plot.legend_temp.SetHeader('T_{Freezeout}', 'C')
    plot.legend_temp.SetBorderSize(0)

    plot.temp_lines = [TLine(), TLine()]
    plot.temp_lines[0].SetLineWidth(2)
    plot.temp_lines[0].SetLineStyle(linestyle_hitemp)
    plot.temp_lines[1].SetLineWidth(2)
    plot.temp_lines[1].SetLineStyle(linestyle_lotemp)

    plot.legend_temp.AddEntry(plot.temp_lines[0], '165 MeV', 'L')
    plot.legend_temp.AddEntry(plot.temp_lines[1], '156 MeV', 'L')

    for i, hist in enumerate(hists, 1):
        c.cd(i)
        hist.Draw()

    for i, key in enumerate(['Ro', 'Rs', 'Rl'], 1):
        c.cd(i)
        for cent in centralities:
            color = get_cent_color(cent)
            theory_data.draw_tlines(key, cent, 'thi', color, linestyle_hitemp)
            theory_data.draw_tlines(key, cent, 'tlo', color, linestyle_lotemp)

    hist_dict = plot.hist_dicts = build_tgraphs_from_data(df)

    for i, key in enumerate(['Ro', 'Rs', 'Rl'], 1):
        c.cd(i)
        for cent in centralities:
            color = get_cent_color(cent)
            tgraph_data, tgraph_sys = hist_dict[key][cent]
            tgraph_data.SetLineColor(color)
            tgraph_data.SetMarkerColor(color)
            tgraph_data.SetMarkerSize(1.25)
            tgraph_data.SetMarkerStyle(8)
            tgraph_data.Draw("SAME P")

            if tgraph_sys:
                sys_color = get_sys_cent_color(cent)
                tgraph_sys.SetFillColor(sys_color)
                tgraph_sys.SetLineColor(color)
                # tgraph_sys.SetLineColor(ROOT.kGray+1)
                # tgraph_sys.SetLineColor(ROOT.kBlack)
                tgraph_sys.SetFillStyle(1001)
                tgraph_sys.Draw('SAME E5')

            if i == 1:
                cent_name = '%2g-%g%%' % tuple(map(int, cent.split('_')))
                plot.legend_cent.AddEntry(tgraph_data, cent_name, 'p')

    c.cd(0)
    plot.legend_cent.Draw()
    plot.legend_temp.Draw()

    return plot

    if cfg_filter:
        df = df[df.cfg==cfg_filter].copy()
    else:
        df = df.copy()

    df = AddRoutRsideRatio(df)

    if merge:
        df = merged_dataframe(df, ['Ro', 'Rs', 'Rl', 'lam', 'RoRs'])

    plot.legend = leg = TLegend()

    get_cent_color = plot.color_loader(palette)

    canvas_div_params = (1, 4, 0, 0)
    canvas_size = (500, 1200)

    canvas_div_params = (2, 2)
    canvas_size = (900, 900)

    c.Divide(*canvas_div_params)
    c.SetCanvasSize(*canvas_size)

    for i in range(4):
        pad = c.cd(i+1)
        pad.SetTickx(1)
        pad.SetTicky(1)

    keys = ['Ro', 'Rs', 'Rl', 'RoRs']
    title_map = {
        'Ro': 'R_{out}',
        'Rs': 'R_{side}',
        'Rl': 'R_{long}',
        'RoRs': 'R_{out} / R_{side}',
    }

    XRANGE = 0.18, 1.02
    YRANGE = 1.23, 9.3

    yranges = [
        YRANGE,
        YRANGE,
        YRANGE,
        (0.45, 1.65),
    ]

    titles = [title_map[key] for key in keys]

    plot.axhists = [ROOT.TH1C('h%d'%i, f'{title}; k_{{T}} GeV; {title}', 1000, *XRANGE)
                    for i, title in enumerate(titles)]

    c.SetFillColor(ROOT.kGreen)

    plot.run2_graphs = []

    for i, hist in enumerate(plot.axhists):
        xax, yax = hist.GetXaxis(), hist.GetYaxis()
        yax.SetRangeUser(*yranges[i])
        yax.SetNdivisions(909)
        xax.SetNdivisions(509)

        xax.SetLabelSize(0.06)
        yax.SetLabelSize(0.065)
        xax.SetTitleOffset(1.1)

    for i, key in enumerate(keys, 1):
        pad = c.cd(i)
        hist = plot.axhists[i-1]
        hist.Draw()
        run1_graphs = run1_graphdict[key]

        for graph, cent in zip(run1_graphs, run1_cents):
            graph.SetMarkerColor(get_cent_color(cent))
            graph.SetMarkerSize(1.2)
            graph.Draw("SAME P")

        for cent, cdf in df.groupby('cent'):
            run2_graph = series_to_TGraphErrors(cdf, key)
            run2_graph.SetMarkerColor(get_cent_color(cent))
            run2_graph.SetMarkerSize(1.2)
            run2_graph.SetMarkerStyle(33)
            run2_graph.Draw("SAME P")
            plot.run2_graphs.append(run2_graph)


def plot_gauss3d_points(df,
                        c=None,
                        kisiel_graphs=None,
                        palette='colorblind'):
    import ROOT
    from ROOT import TCanvas, TLegend

    if c is None:
        c = TCanvas()
        c.Divide(2, 2, 0, 0)
        c.SetCanvasSize(1000, 800)

    if not isinstance(df, dict):
        fit_data = group_df_into_tgraphs(df)
    else:
        fit_data = df

    plot = PlotData(c)

    leg = plot.leg = TLegend(.65, .90, 1, 1)
    leg.SetHeader("Centrality", "C")
    leg.SetNColumns(3)
    leg.SetTextSize(0.03)

    get_color = plot.color_loader(palette)

    if kisiel_graphs is True:
        kisiel_graphs, cent_names = load_kisiel_graphs()

    def plot_kisiel_graphs():
        for g, graph in enumerate(kisiel_graphs[key]):
            graph.GetXaxis().SetRangeUser(0.0, 1.0)
            if g == 0:
                graph.GetYaxis().SetRangeUser(*yrange)
                graph.GetYaxis().SetLabelSize(0.07)
                graph.GetYaxis().SetTitleSize(0.07)
                graph.GetYaxis().SetTitleOffset(0.69)
                graph.GetYaxis().CenterTitle()
                graph.GetXaxis().SetTitleOffset(1.1)
                graph.GetXaxis().SetLabelSize(0.05)
                if i == 4:
                    graph.GetYaxis().SetNdivisions(405)
                else:
                    graph.GetYaxis().SetNdivisions(405)
#                 if i in (3, 4):
                graph.Draw('AP')
            else:
                graph.Draw('SAME P')

            if i == 1:
                leg.AddEntry(graph, cent_names[g], 'P')

    plot.pads = []
    plot.graphs = []

    for i, key in enumerate(('Ro', 'Rs', 'Rl', 'RoRs'), 1):
        pad = c.cd(i)


    #     pad.SetTopMargin(0)
    #     pad.SetBottomMargin(0)
    #     pad.SetLeftMargin(0.17)
#         if i == 4 or i == 3:
#             pad.SetBottomMargin(0.18)
    #     if i in (1, 3):
        pad.SetLeftMargin(0.13)
        pad.SetRightMargin(0.14)
    #     elif i == 1:
    #         pad.SetTopMargin(0.0)

        if i == 3:
            yrange = 0.1, 9.5
        elif i == 4:
            yrange = 0.7, 1.6
        else:
            yrange = 0.1, 7.9

        if kisiel_graphs:
            plot_kisiel_graphs()

        for g, graph in enumerate(fit_data[key]):
            graph.Draw('P')
            plot.graphs.append(graph)

    c.cd(0)
    leg.Draw()
    c.Draw()

    return plot


from . import canvas_divide


def AddRoutRsideRatio(df):
    ratio = df.Ro / df.Rs
    error = ratio * np.hypot(df.Ro_err/df.Ro, df.Rs_err/df.Rs)

    df['RoRs'] = ratio
    df['RoRs_err'] = error
    return df


def calc_weighted_mean(vals, errs, warn=True):
    weights = np.power(errs, -2, where=errs>0, out=np.zeros_like(errs))

    if warn and np.any(weights <= 0.0):
        print(f"Warning: non-positive errors detected", file=sys.stderr)

    if (w := weights.sum()) > 0:
        mean_val = (vals * weights).sum() / w
        err = w ** -0.5
    else:
        mean_val = vals.mean()
        err = 0.0

    return mean_val, err


def merged_dataframe(df, keys, xkey='kT', skip_systematics=False):
    merged_data = []

    merged_keys = ['cent', 'kT']
    for key in keys:
        merged_keys.append(key)
        merged_keys.append(key + '_err')
        if not skip_systematics:
            merged_keys.append(key + '_sys_err')

    for cent, cdf in df.groupby('cent'):

        for kt, kdf in cdf.groupby(xkey):
            values = [cent, kt]

            sys_errors = (None
                          if skip_systematics
                          else calc_df_systematics(kdf, keys))

            for key in keys:
                val, err = calc_weighted_mean(kdf[key], kdf[key + "_err"])

                values.append(val)
                values.append(err)

                if not skip_systematics:
                    values.append(sys_errors[key])

            merged_data.append(values)

    return pd.DataFrame(merged_data, columns=merged_keys)


def plot_gauss3d_df(df,
                    merge=True,
                    cfg_filter=None,
                    c=None,
                    palette='colorblind',
                    systematics=True,
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
        df = merged_dataframe(df, ['Ro', 'Rs', 'Rl', 'lam', 'RoRs'])

    if c is None:
        c = TCanvas()

    plot = PlotData(c)
#     canvas_divide(c)
    cols = 3
    rows = 2
    c.Divide(cols, rows)
    c.SetCanvasSize(1600, 1000)
#     c.SetFillColor(ROOT.kRed)

    from random import random

    get_cent_color = plot.color_loader(palette)

    get_sys_cent_color = plot.color_loader('muted', 'sys')
    for tcolor in plot.tcolor_dict['sys']:
        tcolor.SetAlpha(0.5)

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

    hist_keys = ['Ro', 'Rs', 'Rl', 'lam', 'RoRs']

    YRNG = 1.1, 7.9
    hist_info = {
        'Ro': ('R_{out}', YRNG),
        'Rs': ('R_{side}', YRNG),
        'Rl': ('R_{long}', YRNG),
        'lam': ('#lambda', (0.2, 0.755555)),
        'RoRs': ('R_{out} / R_{side}', (0.5, 1.6)),
    }

    def _build_hist():
        return ROOT.TH1C("h%s" % random(), "; k_{T} (GeV)", 1000, 0.0, 1.5)

    plot.axis_hist_dict = {key: _build_hist() for key in hist_keys}

    for key in hist_keys:
        title, y_range = hist_info[key]
        axh = plot.axis_hist_dict[key]
        axh.SetTitle(title)
        axh.SetTitleOffset(0.0)
        xax, yax = axh.GetXaxis(), axh.GetYaxis()
        xax.SetTitleOffset(0.89)
        xax.SetTitleSize(0.04)
        xax.SetNdivisions(TICKCODE_X)
        xax.SetRangeUser(*KT_RANGE)
        yax.SetRangeUser(*y_range)
        plot.axis_hists.append(axh)

    title_dict = {key: title for key, (title, *_) in hist_info.items()}

    plot.hist_dicts = build_tgraphs_from_data(df,
                                              keys=hist_keys,
                                              title_dict=title_dict)

    shifts = shift * np.linspace(-1.0, 1.0, len(df.cent.unique()))

    ROOT.gStyle.SetErrorX(0)
    ROOT.gROOT.ForceStyle()

    for i, key in enumerate(hist_keys, 1):
        pad = c.cd(i)
        plot.axis_hist_dict[key].Draw()

        for cidx, (cent, cdf) in enumerate(df.groupby('cent')):

            # graph = series_to_TGraphErrors(cdf, key)
            graph, tgraph_sys = plot.hist_dicts[key][cent]
            if shift:
                x = np.frombuffer(graph.GetX(), np.float64, graph.GetN())
                x += shifts[cidx]

            color = get_cent_color(cent)
            graph.SetMarkerColor(color)
            graph.SetLineColor(color)
            graph.SetMarkerStyle(marker_style)
            graph.SetMarkerSize(marker_size)
            graph.Draw("SAME P")

            plot.graphs.append(graph)

            if i == 1:
                legend.AddEntry(graph, f"{cent.replace('_', '-')}%", 'P')

            if tgraph_sys:
                sys_color = get_sys_cent_color(cent)
                tgraph_sys.SetFillColor(sys_color)
                tgraph_sys.SetLineColor(color)
                # tgraph_sys.SetLineColor(ROOT.kGray+1)
                # tgraph_sys.SetLineColor(ROOT.kBlack)
                tgraph_sys.SetFillStyle(1001)
                tgraph_sys.Draw('SAME []')

                if shift:
                    x = np.frombuffer(tgraph_sys.GetX(), np.float64, graph.GetN())
                    x += shifts[cidx]

                    # ex = np.frombuffer(tgraph_sys.GetEX(), np.float64, graph.GetN())
                    # ex = 0.0


    pad = c.cd(5)
    plot.axis_hists[-1].Draw()

    for cidx, (cent, cdf) in enumerate(df.groupby('cent')):
        graph = series_to_TGraphErrors(cdf, 'RoRs')
        if shift:
            x = np.frombuffer(graph.GetX(), np.float64, graph.GetN())
            x += shifts[cidx]
        # graph = series_to_TGraphErrors(cdf, 'RoRs')
        plot.graphs.append(graph)

        color = get_cent_color(cent)
        graph.SetMarkerColor(color)
        graph.SetLineColor(color)
        graph.SetMarkerStyle(marker_style)
        graph.SetMarkerSize(marker_size)

        graph.Draw("SAME P")

    plot.lines.append(ROOT.TLine(KT_RANGE[0], 1.0, KT_RANGE[1], 1.0))
    plot.lines[-1].SetLineStyle(2)
    plot.lines[-1].Draw()

    return plot


def plot_gauss_projections(df, cfg_filter=None, c=None):
    import ROOT
    from ROOT import TCanvas

    if c is None:
        c = TCanvas()

    plot = PlotData(c)
    canvas_divide(c)

    return plot

    colors = ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kOrange

    c.SetFillColor(ROOT.kBlack)

    for i, color in enumerate(colors, 1):
        pad = c.cd(i)
        pad.SetFillColor(color)

#     pad = c.cd(1)
#     pad.SetLeftMargin(3.8)
    return plot
