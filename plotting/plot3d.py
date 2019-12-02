#
# plotting/plot3d.py
#

from plotting.plot1d import build_mrc_dict
from . import PlotData
from .fitresults import MultiFitResults
from .projections import load_datafile, load_fsi, load_mrc
from .projections import build_fsi_obj

_file_cache = {}

def get_fsi_dict(fsi_paths):
    import re
    import ROOT
    from ROOT import TFile
    if isinstance(fsi_paths, str):
        fsi_paths = (fsi_paths, )

    fsi_objs = {}
    for fsi in fsi_paths:
        m = re.match(r'(?P<cls>\w+)\[(?P<file>[^\]]+)\]', fsi)
        d = m.groupdict()
        try:
            fsi_file = _file_cache[d['file']]
        except KeyError:
            fsi_file = _file_cache[d['file']] = TFile.Open(d['file'])

        fsi_cls = getattr(ROOT, d['cls'])
        fsi_objs[fsi] = fsi_cls.From(fsi_file)

    return fsi_objs


class PlotResults3D:
    """
    PlotResults3D
    """

    def __init__(self, fr):
        self.fr = fr

        if isinstance(fr, MultiFitResults):
            for fr in fr.frs:
                if '-0.root' in fr.filename:
                    self.combined_fr = fr
                    break
            else:
                self.combined_fr = None
        else:
            self.combined_fr = fr

    def plot_projections(self,
                         centrality,
                         cfg=None,
                         kt=None,
                         pair='pip',
                         df=None,
                         fclassname=None,
                         ignore_mrc=False,
                         pad=None,
                         canvas_size=(700, 1200),
                         integration=0.05):
        """

        """
        import re
        import ROOT
        from ROOT import TH1, TH3, TCanvas
        from ROOT import TFile, gStyle, gROOT
        from ROOT import TGraphErrors, TLegend, Data3D
        from itertools import product
        from femtofitter import PathQuery

        if self.combined_fr:
            fitresult = self.combined_fr
        else:
            fitresult = self.fr.frs[0]

        df = fitresult.df

        if 'alpha' in df:
            fclassname = 'Fitter3DLevy'
        else:
            fclassname = 'Fitter3DGaussLcms'

        FitterClass = getattr(ROOT, fclassname)

        df = df[df.cent==centrality]
        df = df[df.pair==pair]

        if df.empty:
            raise ValueError(f"No (centrality, pair) ({centrality!r}, {pair!r}) in dataframe")

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


    def plot_projected_fits(self,
                            centrality,
                            fr=None,
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

        if fr is None:
            fr = self.combined_fr or self.fr.frs[0]

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
        if df.empty:
            raise ValueError("bad centrality")

        # if df.mrc.any():
        #     mrc_dict = load_mrc(df.mrc.unique())

        # fsi_dict = load_fsi(df.fsi.unique()[0])
        # load_fsi(fitter, )

        hists = plot.hists = []
        fithists = plot.fithists = []
        fits = []
        kts = []

        plot.fit_tcolor = ROOT.TColor(1.0, 0.0, 0.0, 0.5)
        fit_color = plot.fit_tcolor.GetNumber()
        # fsi = ROOT.FsiStatic.From(1.0)

        for kt, ktdf in df.groupby('kt'):
            kts.append(kt)

            series = ktdf.iloc[I]
            query = PathQuery.From(series)
            fitres = Fitter.FitResult(series)
            fits.append(fitres)
            fit_params = fitres.as_params()

            tdir = tfile.Get(query.as_path())
            data = Data3D.From(tdir)
            fsi = build_fsi_obj(series.fsi, shared_ptr=False)

            # fitter = Fitter(data)
            # fitter.fsi = fsi
            # fitter.fsi = fsi_dict[series.fsi]

            num, den = map(tdir.Get, ("num", "den"))

            rx, ry, rz = make_projections(num, den)
            hists.append((rx, ry, rz))
            qinv = data.src.qinv.get()

            # if series.mrc:
            #     mrc = build_mrc_(series.mrc, shared_ptr=False)
            #     fitden = mrc.GetSmearedDen()
    #         fitnum = fitter.mrc.GetSmearedFit(fit_params, qinv, fsi)
    #         fitnum.Multiply(fitden

    #         np.frombuffer(fitden.GetArray(), np.float32, fitden.GetNcells())[:] = 1.0
            fitden = data.src.den.Clone('fitden')
            fitnum = fitden.Clone('fitnum')
            fit_params.multiply(fitnum, qinv, fsi)

            fitnum.SetLineColor(fit_color)
            fitnum.SetMarkerColor(fit_color)

            fx, fy, fz = make_projections(fitnum, fitden)
            fithists.append((fx, fy, fz))
            for f in fithists[-1]:
                f.SetLineWidth(2)
                # norm_bins = tuple(map(f.GetXaxis().FindBin, (0.07, 0.12)
                # norm_bins = f.GetNbinsX() - 3, f.GetNbinsX()
                # scale = f.Integral(*norm_bins)
                # f.Scale(3.0/scale)
                f.Scale(1.0 / series.norm)

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
