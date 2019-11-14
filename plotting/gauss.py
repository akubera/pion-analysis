#
# plotting/gauss.py
#


from . import PlotData
import numpy as np
import pandas as pd


KISIEL_GRAPHS = None


def load_kisiel_graphs(cent_names=None):
    import feather
    import ROOT
    from collections import defaultdict
    global KISIEL_GRAPHS
    
#     if KISIEL_GRAPHS:
#         return KISIEL_GRAPHS

    kisiel_df_filename = "/home/akubera/Physics/pion-analysis/Run1Data.feather"
    kisiel_df = feather.read_dataframe(kisiel_df_filename)
    
    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kMagenta-3, ROOT.kCyan+3]
    markers = [25]* len(colors)
    sizes = [0.8]* len(colors)

    titles = {
        'Ro': 'R_{out} (fm)',
        'Rs': 'R_{side} (fm)',
        'Rl': 'R_{long} (fm)'
    }

    KISIEL_GRAPHS = kisiel_graphs = defaultdict(list)
    if cent_names is None:
        cent_names = []
        
    for key in ('Ro', 'Rs', 'Rl'):
        for (cent, cdf), color, marker, msize in zip(kisiel_df.groupby('cent'), colors, markers, sizes):
            graph = series_to_TGraphErrors(cdf, key)
            graph.SetMarkerColor(color)
            graph.SetMarkerStyle(marker)
            graph.SetMarkerSize(msize)
            graph.SetTitle("; k_{T} (GeV/c); %s;" % titles[key])
            kisiel_graphs[key].append(graph)
            if key == 'Ro':
                cent_names.append(cdf.iloc[0].centname)

    for (cent, cdf), color, marker, msize in zip(kisiel_df.groupby('cent'), colors, markers, sizes):
        en = np.array(cdf["Ro_stat_err"])
        ed = np.array(cdf["Rs_stat_err"])

        os_ratio = np.array(cdf['Ro'] / cdf['Rs'])
        e = os_ratio * np.sqrt((en / cdf["Ro"])**2 + (ed / cdf["Rs"])**2)

    #     e = np.zeros_like(en)
        graph = build_TGraphErrors(cdf.kT, cdf['Ro'] / cdf['Rs'], e)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(msize)
        graph.SetTitle("; k_{T} (GeV/c); R_{out} / R_{side};")
    #     cent_names.append(cdf.iloc[0].centname)
        kisiel_graphs['RoRs'].append(graph)
        
    return kisiel_graphs, cent_names


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



def series_to_TH1D(series, ykey, ekey=None, xkey='kT'):
    ekey = ekey or ykey + '_err'
    
    graph = build_TGraphErrors(series[xkey],
                               series[ykey],
                               series[ekey])
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
        
    return graphs


def plot_gauss3d_points(df, c=None, kisiel_graphs=None):
    import ROOT
    from ROOT import TCanvas, TLegend
    
    if c is None:
        c = TCanvas()
        c.Divide(2, 2, 0, 0)
        c.SetCanvasSize(1000, 800)
    
    if not isinstance(df, dict):
        fit_data = group_df_into_tgraphs(df)
#         print(fit_data)
    else:
        fit_data = df
        
    plot = PlotData(c)
    
    leg = plot.leg = TLegend(.65, .90, 1, 1)
    leg.SetHeader("Centrality", "C")
    leg.SetNColumns(3)
    leg.SetTextSize(0.03)
    
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


def plot_gauss_df(df, cfg_filter=None, c=None):
    import ROOT
    from ROOT import TCanvas
    import seaborn as sns

    if c is None:
        c = TCanvas()

    plot = PlotData(c)
#     canvas_divide(c)
    cols = 3
    rows = 2
    c.Divide(cols, rows)
    c.SetCanvasSize(800, 800)
    c.SetFillColor(ROOT.kRed)

#     h = c.h = ROOT.TH1D("h", "Hi", 70, -3, 1)
#     h.FillRandom("gaus")

    tcolors = plot.tcolors = [ROOT.TColor(*rgb) for rgb in sns.color_palette("colorblind")]
    colors = [tcolor.GetNumber() for tcolor in tcolors]
    color_it = iter(colors)

    centrality_colors = {}
    
    def centrality_color(cent):
        if cent not in centrality_colors:
            centrality_colors[cent] = next(color_it)
        return centrality_colors[cent]

    plot.graphs = []
    plot.lines = []
    for i, key in enumerate(['Ro', 'Rs', 'Rl'], 1):
        pad = c.cd(i)
        for cidx, (cent, cdf) in enumerate(df.groupby('cent')):
            
            graph = series_to_TGraphErrors(cdf, key)
            graph.SetMarkerColor(centrality_color(cent))
            graph.SetLineColor(centrality_color(cent))
            graph.SetMarkerStyle(21)
            # graph.SetMarkerSize(msize)
            
            if cidx == 0:
                graph.SetTitle({'Ro': 'R_{out}', 'Rs': 'R_{side}', 'Rl': 'R_{long}'}[key])
                graph.SetTitleOffset(1.0)
                graph.Draw("AP")
                graph.GetXaxis().SetNdivisions(505)
                graph.GetXaxis().SetRangeUser(0.2, 1.25)
                graph.GetYaxis().SetRangeUser(1.0, 8.0)
            else:
                graph.Draw("SAME P")
                
            plot.graphs.append(graph)

    df = AddRoutRsideRatio(df.copy())           

    pad = c.cd(4)
    for cidx, (cent, cdf) in enumerate(df.groupby('cent')):
        graph = series_to_TGraphErrors(cdf, 'RoRs')
        plot.graphs.append(graph)
        
        graph.SetMarkerColor(centrality_color(cent))
        graph.SetLineColor(centrality_color(cent))
        graph.SetMarkerStyle(21)
        graph.SetMarkerSize(0.8)
        
        if cidx == 0:
            graph.Draw("AP")
            graph.GetXaxis().SetNdivisions(505)
            graph.GetXaxis().SetRangeUser(0.2, 1.3)
            graph.GetYaxis().SetRangeUser(0.5, 1.6)
        else:
            graph.Draw("SAME P")

        plot.lines.append(ROOT.TLine(0.2, 1.0, 1.2, 1.0))
        plot.lines[-1].SetLineStyle(2)
        plot.lines[-1].Draw()

    pad = c.cd(1)
    pad.SetLeftMargin(100.01)
    
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
