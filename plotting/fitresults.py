#
# pion-analysis/plotting/fitresults.py
#

import re
import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt


from femtofitter import FitResults


class MultiFitResults:


    def __init__(self, pattern, datapath='fitresults'):
        self.datapath = Path(datapath)
        if isinstance(pattern, str):
            self.paths = list(self.datapath.glob(pattern))
        else:
            self.paths = [x for x in (self.datapath.glob(patt) for patt in pattern) for x in x ]

        fr_list = [FitResults(filename) for filename in self.paths]
        dfs = [fr.df for fr in fr_list]
        # for i, df in enumerate(dfs, 1):

        for df, path in zip(dfs, self.paths):

            if (m := re.search(r"(?P<partition>\d+)\.json$", path.name)):
                df['partition'] = int(m.group('partition'))
            else:
                df['partition'] = 1

        df = self.df = pd.concat(dfs).reset_index(drop=True)
        df.loc[df.subset.isnull(), 'subset'] = ''

        if pd.np.any(df.duplicated()):
            raise ValueError("Paths include non-unique fit results")


    def plotstuff(self, df=None, y='Ro'):
        if df is None:
            df = self.df

        plot_opts = dict(col='kT', row='cent')
        fg = sns.catplot(data=df, y=y, hue='partition', **plot_opts)
        return fg

    def otherplot(self, df=None, y='Ro'):
        if df is None:
            df = self.df

        def plotmaker(**kwargs):
            data = kwargs.pop('data')
            ax = plt.gca()

            data = kwargs.pop('data')
            ax = plt.gca()
            plot_opts = {}
            g = sns.catplot(data=data, y=y, x='partition', hue='magfield', ax=ax, **plot_opts)

        fg = sns.FacetGrid(df, col='kT', row='cent', hue='pair')
        fg.map_dataframe(plotmaker)
        return fg
