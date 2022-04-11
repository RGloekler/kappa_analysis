# scaling.py - compute scaling factor for non-identical waveguides
# Ryan Gloekler, 2022
# ECSyD Laboratory

# USAGE:
# import scaling (have this file in directory with file importing from)
#
# IMPORTANT: ENSURE 'scaling_eqn.csv' is in the located at 'data/' or appropriate directory
#
# F, kappa = scaling.get_kappa(kappa_identical, radius, in_wdth, ring_wdth, gap)
# This returns the scaling factor for kappa, and the non-identical kappa value

import pandas as pd
import numpy as np
import sys

# solve the second order polynomial for F
def solve(R, second, first, c):
    return second * (R ** 2) + first * (R) + c

# used to return the proper value of scaling factor for input values, returns kappa
def get_kappa(kappa_id, radius, in_wdth, ring_wdth, gap):
    # open the datafile
    df = open_csv('data/scaling_eqn.csv')

    # scale the values to the closest equation for F
    if in_wdth % 50 != 0: in_wdth = round(in_wdth / 50) * 50
    if ring_wdth % 50 != 0: ring_wdth = round(ring_wdth / 50) * 50
    if gap not in [100, 150, 200]: gap = min([100, 150, 200], key=lambda x:abs(x-gap))

    # for testing
    # print(in_wdth, ring_wdth, gap)

    # access lookup table for these values...
    row = df.loc[(df['in_wdth'] == in_wdth) & (df['ring_wdth'] == ring_wdth) & (df['gaps'] == gap)]

    # pull the equation coefficients from the spreadsheet
    second, first, const = row.iloc[0]['square'], row.iloc[0]['first'], row.iloc[0]['const']

    # compute F and return
    F = solve(float(radius), float(second), float(first), float(const))
    kappa  = kappa_id / F
    return F, kappa

# opens the datafile, returns a pandas dataframe
def open_csv(name):
    df = pd.read_csv(name)
    return df
