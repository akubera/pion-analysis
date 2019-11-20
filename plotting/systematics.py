#
# pion-analysis/plotting/systematics.py
#


import pandas as pd


def calc_df_systematics(df, keys):

    def _calc_sys_err(df, k):
        return calc_weighted_variance(df[k], df[k + "_err"])

    return {k: _calc_sys_err(df, k) for k in keys}


def calc_weighted_variance(vals, errs):
    weights = errs ** -2

    mean = (vals * weights).sum() / weights.sum()

    M = (weights > 0).sum()

    # res_1 = (weights * (vals - mean) ** 2).sum() / ((M - 1) / M * weightssum())
    res = (weights * (vals - mean) ** 2).sum() / (weights.sum() - (weights**2).sum() / weights.sum())
    # print(pd.np.sqrt(res), pd.np.sqrt(res_1), vals.std())
    return pd.np.sqrt(res)
