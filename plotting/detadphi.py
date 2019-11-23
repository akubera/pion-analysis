#
# pion-analysis/plotting/detadphi.py
#


from . import PlotData
from pathlib import Path


class DetaDphiPlot:
    """
    Deta Dphi
    """

    def __init__(self, tfile, root='PWG2FEMTO'):
        if isinstance(tfile, (str, Path)):
            tfile = Path(tfile).expanduser()
            if not tfile.exists():
                raise FileNotFoundError(tfile.absolute())
                
            from ROOT import TFile
            tfile = TFile.Open(str(tfile))

        self.tfile = tfile
        self._rootcontainer = tfile.Get(root)

        self.container = self.find_first_subdir(self._rootcontainer)
        self.analysis = self.find_first_subdir(self.container)

    def plot(self, pad=None):
        if pad is None:
            from ROOT import TCanvas
            pad = TCanvas()

        plot = PlotData(pad)
        pad.cd(0)

        n, d = map(self.analysis.Get, ("NumDEtaDPhiStar", "DenDEtaDPhiStar"))
        plot.num = n
        plot.den = d

        norm_bins = 2, 50, 2, 12
        plot.ratio = n.Clone("ratio")
        plot.ratio.Divide(n, d, d.Integral(*norm_bins), n.Integral(*norm_bins))
        plot.ratio.Draw("COL")

        return plot 

    def plot_eta(self, pad=None, N=1, v=None):
        from ROOT import TCanvas, TLine
        
        if pad is None:
            pad = TCanvas()
    
        plot = PlotData(pad)
        pad.cd(0)

        n, d = map(self.analysis.Get, ("NumDEtaDPhiStar", "DenDEtaDPhiStar"))
        if n.GetSumw2N() == 0:
            n.Sumw2()
        
        zero_bin = n.GetYaxis().FindBin(0.0)
        project_bins = zero_bin - N, zero_bin + N
        
        plot.num = nx = n.ProjectionX('nx', *project_bins)
        plot.den = dx = d.ProjectionX('dx', *project_bins)
        
        norm_bins = 2, 20
        plot.ratio = nx.Clone("ratio")
        plot.ratio.SetTitle("#Delta#eta Ratio")
        plot.ratio.Divide(nx, dx, dx.Integral(*norm_bins), nx.Integral(*norm_bins))
        plot.ratio.Draw()
        
        ax = plot.ratio.GetXaxis()
        
        plot.line = TLine(ax.GetXmin(),1.0, ax.GetXmax(), 1.0)
        plot.line.SetLineStyle(2)
        plot.line.Draw()
        
        if v is not None:
            plot.vline = TLine(v, 0.0, v, 1.0)
            plot.vline.SetLineStyle(2)
            plot.vline.Draw()
        
        return plot
        
    def plot_phi(self, pad=None, N=1, v=None):
        from ROOT import TCanvas, TLine
        
        if pad is None:
            pad = TCanvas()
    
        plot = PlotData(pad)
        pad.cd(0)

        n, d = map(self.analysis.Get, ("NumDEtaDPhiStar", "DenDEtaDPhiStar"))
        if n.GetSumw2N() == 0:
            n.Sumw2()
        
        zero_bin = n.GetXaxis().FindBin(0.0)
        project_bins = zero_bin - N, zero_bin + N
        
        plot.num = ny = n.ProjectionY('ny', *project_bins)
        plot.den = dy = d.ProjectionY('dy', *project_bins)
        
        norm_bins = 2, 20
        plot.ratio = ny.Clone("ratio")
        plot.ratio.SetTitle("#Delta#phi Ratio")
        plot.ratio.Divide(ny, dy, dy.Integral(*norm_bins), ny.Integral(*norm_bins))
        plot.ratio.Draw()
        
        ay = plot.ratio.GetXaxis()
        
        plot.line = TLine(ay.GetXmin(),1.0, ay.GetXmax(), 1.0)
        plot.line.SetLineStyle(2)
        plot.line.Draw()
        
        if v is not None:
            plot.vline = TLine(v, 0.0, v, 1.0)
            plot.vline.SetLineStyle(2)
            plot.vline.Draw()
            
        return plot 
        
    @staticmethod
    def find_first_subdir(tdir):
        for key in tdir.GetListOfKeys():
            container = key.ReadObj()
            if container.InheritsFrom("TDirectory"):
                return container
        else:
            raise 
