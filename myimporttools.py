import sys
import itertools
import numpy as np
import os

spinner = itertools.cycle(['-', '/', '|', '\\'])

def quickRead(fn):
    with open(fn) as data:
        header = data.readline().rsplit()
        print(header)
        entries = []
        for line in data:
            fields = line.rsplit()
            try:
                values = np.array([[float(elem) for elem in fields]])
            except ValueError:
                values = np.array([fields])
            # print(values)
            if entries==[]:
                entries = values
            else:
                entries = np.append(entries, values, axis=0)
            sys.stdout.flush()  # flush stdout buffer
            sys.stdout.write('\b'+next(spinner))  # move spinner
    return entries

def readDir(pwd_path, find_str):
    # pwd_path = "/Users/andreashaller/results/SkyrmionLattice/square_lattice/MPS_results/M256"
    # find_str = "*local*dat"  # to select a subset, use "*<substr>*.dat", e.g. "*B_0*dat"

    # extract all files in the subdirectories of pwd_path
    message = "searching for " + find_str + ": "
    sys.stdout.write(message)
    filenames = []
    for filename in os.popen("find " + str(pwd_path) + " -path "
                            + '"' + str(find_str) + '"').read().split('\n')[0:-1]:
        filenames.append(filename)
    filenames = np.sort(filenames)
    message = str(len(filenames)) + " files found.\n"
    sys.stdout.write(message)

    return filenames

def replaceBraKet(expr):
    return expr.replace("<","$\\langle ").replace(">","\\rangle $")