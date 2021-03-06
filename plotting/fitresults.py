#
# pion-analysis/plotting/fitresults.py
#

from typing import Union, Tuple

import sys
import re
import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt


from femtofitter import FitResults
from .systematics import calc_df_systematics, calc_weighted_mean


class MultiFitResults:

    def __init__(self, pattern, datapath='fitresults'):
        self._datapath = Path(datapath)
        assert self._datapath.exists()

        if isinstance(pattern, (str, Path)):
            self.paths = list(self._datapath.rglob(pattern))
        else:
            self.paths = [filepath
                          for glob in map(self._datapath.glob, pattern)
                          for filepath in glob]

        self.frs = fr_list = [FitResults(filename) for filename in self.paths]
        dfs = [fr.df for fr in fr_list]

        for df, path in zip(dfs, self.paths):

            if (m := re.search(r"(?P<partition>\d+)\.json$", path.name)):
                df['partition'] = int(m.group('partition'))
            else:
                df['partition'] = 1

        df = self.df = pd.concat(dfs).reset_index(drop=True)
        if 'subset' not in df.columns:
            df['subset'] = ''
        else:
            df.loc[df.subset.isnull(), 'subset'] = ''

        df['pair:field'] = df.apply(lambda row: f'{row.pair}:{row.magfield}', axis=1)
        df['part:field'] = df.apply(lambda row: f'{row.partition}:{row.magfield}', axis=1)

        if pd.np.any(df.duplicated()):
            raise ValueError("Paths include non-unique fit results")

        self.orig_df = self.df.copy()

        # combined configurations
        cfg = self.frs[0].config
        for fr in self.frs[1:]:
            new_config_keys = set(fr.config.index) - set(cfg.index)
            if new_config_keys:
                missing_config = fr.config.loc[new_config_keys]
                cfg = pd.concat([cfg, missing_config], sort=True)

        self.cfg_df = cfg

        # find unique configuration properties
        cfg_unique_keys = cfg.columns[cfg.nunique(axis=0)!=1]
        self.feature_df = cfg[[*cfg_unique_keys]]

    def remove_nan_results(self, keys=None, *, inplace=True):
        if keys is None:
            keys = [col for col in self.df.columns if col.endswith('_err')]

        bad_rows = self.df[keys].isnull().any(axis=1)
        df = self.df[~bad_rows].reset_index(drop=True)
        if inplace:
            self.df = df

        return df

    def remove_unfit_lambdas(self, limit=0.95, *, inplace=True):
        if 'lam' not in self.df.columns:
            print("lambda not in this fit result", file=sys.stderr)
            return

        bad_rows = self.df.lam >= limit
        df = self.df[~bad_rows].reset_index(drop=True)
        if inplace:
            self.df = df

        return df

    @classmethod
    def range_reduce(cls, some_dict):
        """
        Searches through dictionary for dicts of form {"MIN": x, "MAX": y}
        and replaces that value with tuple (x, y)
        """
        for key, value in some_dict.items():
            if isinstance(value, dict):
                if {'MIN', 'MAX'} == set(value.keys()):
                    some_dict[key] = (value['MIN'], value['MAX'])
                else:
                    cls.range_reduce(value)

    @classmethod
    def build_feature_string(cls, feature_dict):
        return ', '.join(cls.build_feature_strlist(feature_dict))

    @classmethod
    def build_feature_strlist(cls, feature_dict):
        """
        """
        import re
        fdict = feature_dict.copy()
        keys = '|'.join(feature_dict.keys())

        if (m1 := re.search(r"\b(?P<name>\w+)_bins\b", keys)):
            name = m1.group('name')
            if (m2 := re.search(rf"\b{name}_range\b", keys)):
                val_bins = fdict.pop(m1.group())
                val_range = fdict.pop(m2.group())

                def stringify_number(num):
                    if num % 1:
                        return '%02g' % num
                    else:
                        return '%5.1f' % num

                strval = f'[{val_bins}|%s %s]' % tuple(map(stringify_number, val_range))
                keystr = f"{name}{strval}"

                return [keystr] + cls.build_feature_strlist(fdict)

        parts = []

        for k, v in fdict.items():
            if isinstance(v, dict):
                sublist = cls.build_feature_strlist(v)
                if not k.endswith("_cut"):
                    sublist = [f'{k}/{val}' for val in sublist]
                parts += sublist
            else:
                parts.append('%s:%s' % (k, v))

        return parts

    def add_feature_strings(self, *, inplace=False, sort=True, features=None):
        if 'feature' in self.df.columns:
            return self.df

        feature_dict = {}

        if features is not None:
            if isinstance(features, str):
                features = (features, )

            for fr in self.frs:
                missing_cfgs = set(fr.config.index) - set(feature_dict)
                for cfg in missing_cfgs:
                    series = fr.config.loc[cfg]
                    feature_dict[cfg] = {feature_name: series[feature_name]
                                         for feature_name in features}

        else:

            for cfg, x in self.feature_df.iterrows():
                features = feature_dict[cfg] = {}

                for key, value in x.items():
                    *keypath, keyname = key.split('.')
                    featdict = features
                    for k in keypath:
                        if k not in featdict:
                            featdict[k] = {}
                        featdict = featdict[k]
                    featdict[keyname] = value

                # reduce ranges
                self.range_reduce(features)

        df = self.df if inplace else self.df.copy()

        for cfg in df.cfg.unique():
            feature_str = self.build_feature_string(feature_dict[cfg])
            df.loc[df.cfg==cfg, 'feature'] = feature_str

        if sort:
            df['feature_sort'] = df.feature.apply(self.split_features)
            df.sort_values(['feature_sort', 'cent', 'kt'], inplace=True)
            df.reset_index(drop=True, inplace=True)
            df.drop(columns='feature_sort', inplace=True)

        return df

    @staticmethod
    def split_features(feature: str) -> Union[Tuple[str, float], str]:
        name, _, val = feature.partition(':')
        try:
            return name, float(val)
        except ValueError:
            return feature

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

    def get_features(self):
        raise NotImplementedError

    def get_merged_dataframe(self, keys,
                             df=None,
                             xkey='kT',
                             groups=['cent'],
                             skip_systematics=False):
        """
        Merge dataframe rows based on centrality and k
        """

        if df is None:
            df = self.df

        grouping_keys = [groups] if isinstance(groups, str) else list(groups)
        merged_keys = [*grouping_keys, xkey]

        for key in keys:
            merged_keys.append(key)
            merged_keys.append(key + '_err')
            if not skip_systematics:
                merged_keys.append(key + '_sys_err')

        merged_data = []

        def _recursive_add_groups(subdf, values, groupnames):
            if groupnames:
                for v, vdf in subdf.groupby(groupnames[0]):
                    _recursive_add_groups(vdf, values + [v], groupnames[1:])
                return

            for x, xdf in subdf.groupby(xkey):
                val_list = values + [x]
                sys_errors = (None
                              if skip_systematics
                              else calc_df_systematics(xdf, keys))
                for key in keys:
                    val, err = calc_weighted_mean(xdf[key], xdf[key + "_err"])
                    val_list.append(val)
                    val_list.append(err)
                    if not skip_systematics:
                        val_list.append(sys_errors[key])
                merged_data.append(val_list)


        _recursive_add_groups(df, [], grouping_keys)

        return pd.DataFrame(merged_data, columns=merged_keys)
