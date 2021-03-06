#
# pion-analysis/plotting/__init__.py
#

import json
import pandas as pd
import numpy as np


class PlotData:

    def __init__(self, c=None):
        self.canvas = c
        self.tcolor_dict = {}
        self.color_dict = {}

    def __getattr__(self, key):
        return getattr(self.canvas, key)

    def color_loader(self, palette, key=''):
        from ROOT import TColor
        import seaborn as sns

        colors = sns.color_palette(palette)

        if key not in self.tcolor_dict:
            self.tcolor_dict[key] = [TColor(*rgb) for rgb in colors]
            colors = self.color_dict[key] = [tcolor.GetNumber() for tcolor in self.tcolor_dict[key]]
        else:
            colors = self.color_dict[key]

        color_it = iter(colors)

        centrality_colors = {}

        def centrality_color(cent):
            cent = cent.rstrip('%')
            centparts = cent.split('-') if '-' in cent else cent.split('_')
            cent_key = tuple('%g' % float(f) for f in centparts)
            if cent_key not in centrality_colors:
                centrality_colors[cent_key] = next(color_it)
            return centrality_colors[cent_key]

        return centrality_color

    @property
    def tcolors(self):
        return self.tcolor_dict.get('', [])

    @property
    def colors(self):
        return self.color_dict.get('', [])

    def set_title_size(self, size):
        self.Update()
        self.GetPrimitive('title').SetTextSize(size)

    @staticmethod
    def get_random_str(N=10, prefix='', suffix=''):
        from random import choices
        from string import ascii_letters
        return prefix + ''.join(choices(ascii_letters, k=N)) + suffix


class FitData:

    def __init__(self, path):
        self.data = json.loads(Path(path).read_text())
        self.df = pd.DataFrame(self.data['df'])


def canvas_divide(c, *, rows=2, cols=3, rmargin=0, lmargin=0, tmargin=0):
    from itertools import product
    import ROOT
    from ROOT import TPad

    rect = c.GetBBox()
    # print(rect.fX, rect.fY, rect.fWidth, rect.fHeight)

    c.Divide(cols, rows)

    pad = c.cd(1)
    rect = pad.GetBBox()
    print(rect.fX, rect.fY, rect.fWidth, rect.fHeight)

#     pad.SetBBoxX1(0)
    pad.SetLeftMargin(3)
    print(rect.fX, rect.fY, rect.fWidth, rect.fHeight)

#     pad = TPad("", "", 0, 0.3, 0.1, 0.2)
    c.SetFillColor(ROOT.kRed)
#     pad.Draw()
#     rect = pad.GetBBox()
#     print(rect.fX, rect.fY, rect.fWidth, rect.fHeight)
    return


#     for i, j in product(range(cols), range(rows)):
    for i in range(cols*rows):
        row, col = i // cols, i % cols
        pad = c.cd(i+1)

        if i == 0:
#             pad =
#             continue
            rect = pad.GetBBox()
            center = pad.GetBBoxCenter()
            print(rect.fX, rect.fY, rect.fWidth, rect.fHeight)
            x1 = 0
            x2 = 26
            y1 = 0
            y2 = 100

#             pad.SetBBoxX1(x1)
            pad.SetBBoxX2(x2)
#             pad.SetBBoxY1(y1)
#             pad.SetBBoxY2(y2)

            rect = pad.GetBBox()
            print(rect.fX, rect.fY, rect.fWidth, rect.fHeight)

#             pad.SetBBoxCenterX(50)
            pad.SetBBoxCenterY(0)

            continue

#         pad.SetBBoxCenter(3 * col + 9)
        x1 = col * 100 // (cols+1)
        x2 = (col+1) * 100 // (cols+1)
        y1 = 0
        y2 = 100

        pad.SetBBoxX1(x1)
        pad.SetBBoxX2(x2)
        pad.SetBBoxY1(y1)
        pad.SetBBoxY2(y2)

#         pad.SetBBoxX2(10 * col + 1)

#         if i > 1:
#             break
