
# plotting/compare_tof.py

from . import PlotData


def compare_tof(tdir, title=''):
    from ROOT import TCanvas, TPad

#     tdir.Get("Tracks/pass").ls()
    eta_pt = tdir.Get("Tracks/pass/EtaPt")
#     eta_pt = tdir.Get("Tracks/pass/PtPhi")
    events = tdir.Get("Event/pass/cent_mult")
    print(events, events.GetEntries())
    pt = eta_pt.ProjectionY("pt")
    pt.SetTitle(title)
    pt.Scale(1.0 / events.GetEntries())
    c = TCanvas()
    
    plot = PlotData(c)
    pt.Draw('H')
    plot.pt = pt
    
    plot.subpad = TPad('subpad', '', 0.4, 0.4, 0.75, 0.85)
    plot.subpad.cd()
    eta_pt.Draw("COL")
    
    c.cd()
    plot.subpad.Draw()
    
    return plot
