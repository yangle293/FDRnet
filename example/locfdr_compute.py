from __future__ import print_function
import sys
sys.path.insert(0, "../locfdr-python")
from locfdr import locfdr
from collections import OrderedDict
from os import listdir
from os.path import isfile, join, basename, splitext
import math
import numpy as np
import re
from scipy.stats import norm as norm

def run(pfilename):
    ff = open(join(pfilename))
    gene = []
    pp = []
    zz = []
    for line in ff:
        gene.append(line.strip().split('\t')[0])
        pp.append(float(line.strip().split('\t')[1]))
    ff.close()
    zz = -norm.ppf(pp)
    # eliminate the -infs
    index_noninfs = [i for i in range(len(zz)) if zz[i] != -float('Inf') and zz[i] != float('Inf')]
    tmp = [zz[i] for i in index_noninfs]

    results = locfdr([zz[i] for i in index_noninfs],saveplot=True,saveroot=pfilename,showplot=False)
    fdr_noninfs = results['fdr']
    fdr_dict = OrderedDict(zip(index_noninfs,fdr_noninfs))
    fdr = [1 for i in range(len(zz))]
    for i in range(len(zz)):
        if i in index_noninfs:
            fdr[i] = fdr_dict[i]
        elif zz[i] == float('Inf'):
            fdr[i] = 0

    for i in range(len(fdr)):
        if zz[i] < 0:
            fdr[i] = 1.0
        if (zz[i] > 0) and math.isnan(fdr[i]):
            fdr[i] = 0.0
    data = OrderedDict(zip(gene,fdr))
    output = splitext(pfilename)[0] + '_lfdr.txt'
    with open(output,'w') as fout:
        for (gene,fdr) in data.items():
            fout.write("{}\t{}\n".format(gene,fdr))

if __name__ == "__main__":
    run(sys.argv[1])
