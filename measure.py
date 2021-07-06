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

subdir = "/home/andreas/mckit/results"
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
                tp = datas[idx]  # the first match
            else:
                tp += datas[idx]  # add all matches
    local_averages.append(tp/ns)  # divide by the number of matches
[print(la) for la in local_averages]  # show the estimates

correlations = []