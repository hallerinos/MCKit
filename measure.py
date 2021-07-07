import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.lines import Line2D
import numpy as np
import os
import sys
import itertools

import re

import myimporttools as my
import pandas as pd

# plt.style.use('ggplot')

subdir = "./results"
# pwd_paths = [subdir+"M64",subdir+"M128",subdir+"M256"]
pwd_paths = [subdir]
obs_types = ['*']

markers = itertools.cycle(Line2D.filled_markers)

for obs_type in obs_types:
    find_str = '*'+obs_type+'*dat'  # to select a particular substr, use "*<substr>*.dat"
    allData = []
    allFilenames = []
    nf = 0
    for pwd_path in pwd_paths:
        filenames = my.readDir(pwd_path, find_str)
        allFilenames.append(filenames)
        # pnf = nf
        # nf = len(filenames)
        # if (pnf != 0) and (nf != pnf):
        #     raise IndexError("Cannot fully compare data sets - different number of files.")

        datas = []
        for fn in filenames:
            data = pd.read_table(fn, delim_whitespace=True)
            datas.append(data)
        allData.append(datas)

parameters = np.array([re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", fn) for fn in filenames])

magfields = [pp[0] for pp in parameters]
magfields = list(dict.fromkeys(magfields))  # delete duplicates
print(magfields)

local_averages = []
for mf in magfields:
    ns = 0
    for (idx, pp) in enumerate(parameters):
        if pp[0] == mf:  # check which snapshots match the magnetic field
            ns += 1
            if ns == 1:
                tp = np.copy(datas[idx])  # the first match
            else:
                tp += datas[idx]  # add all matches
    local_averages.append(tp/ns)  # divide by the number of matches
[print(la) for la in local_averages]  # show the estimates

local_averages[0]
datas[0].iloc[1,3]

correlations = []
for mf in magfields:
    ns = 0
    for (idx, pp) in enumerate(parameters):
        if pp[0] == mf:  # check which snapshots match the magnetic field
            ns += 1
            data = pd.DataFrame([[datas[idx].iloc[jj,1], datas[idx].iloc[jj,2], datas[idx].iloc[jjpr,1], datas[idx].iloc[jjpr,2], datas[idx].iloc[jj,3]*datas[idx].iloc[jjpr,3], datas[idx].iloc[jj,4]*datas[idx].iloc[jjpr,4], datas[idx].iloc[jj,5]*datas[idx].iloc[jjpr,5]] for jj in range(len(datas[idx])) for jjpr in range(len(datas[idx]))], columns=["x", "y", "xpr", "ypr", "S_{x,r}S_{x,rpr}", "S_{y,r}S_{y,rpr}", "S_{z,r}S_{z,rpr}"])
            if ns == 1:
                tp = np.copy(data)  # the first match
            else:
                tp += data  # add all matches
    correlations.append(tp/ns)  # divide by the number of matches
[print(la) for la in correlations]  # show the estimates

corrected_correlations = []
for (idx, dd) in enumerate(correlations):
    ns = 0
    data = pd.DataFrame([[local_averages[idx].iloc[jj,1], local_averages[idx].iloc[jj,2], local_averages[idx].iloc[jjpr,1], local_averages[idx].iloc[jjpr,2], local_averages[idx].iloc[jj,3]*local_averages[idx].iloc[jjpr,3], local_averages[idx].iloc[jj,4]*local_averages[idx].iloc[jjpr,4], local_averages[idx].iloc[jj,5]*local_averages[idx].iloc[jjpr,5]] for jj in range(len(local_averages[idx])) for jjpr in range(len(local_averages[idx]))], columns=["x", "y", "xpr", "ypr", "S_{x,r}S_{x,rpr}", "S_{y,r}S_{y,rpr}", "S_{z,r}S_{z,rpr}"])
    tp = np.copy(dd)  # the first match
    tp -= data  # add all matches
    corrected_correlations.append(tp)  # divide by the number of matches
[print(la) for la in corrected_correlations]  # show the estimates