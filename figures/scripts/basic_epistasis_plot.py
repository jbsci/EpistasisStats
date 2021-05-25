#!/usr/bin/env python

'''
Script to generate basic epistasis scatter plots

Written 2020 by Jonathan Barnes

Contact: j@jbsci.dev
'''


#---imports---#

import os,sys,argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr,linregress
from matplotlib.pyplot import cm
import pandas as pd
from matplotlib import rc,rcParams
import matplotlib.patches as mpatches


#---functions---#


def plotter1(fs, out, ptype):
    rcParams['xtick.labelsize'] = fs
    rcParams['ytick.labelsize'] = fs
    rcParams['axes.labelsize'] = fs
    rcParams['legend.fontsize'] = fs
    rcParams['font.family'] = 'serif'
    rc('text', usetex=True)
    plt.figure(figsize=(8,8))
    plt.scatter(null_ep_x, null_ep_y, edgecolor='k', facecolor='None', label = 'No epistasis')
    plt.scatter(pos_ep_x, pos_ep_y, edgecolor='red', facecolor='None', label = 'More stable than additive')
    plt.scatter(neg_ep_x, neg_ep_y, edgecolor='blue', facecolor='None', label = 'Less stable than additive')
    plt.legend(loc=4, prop={'size' : 20}, frameon=False, labelspacing=0, handletextpad=0)
    plt.ylabel('$\Delta\Delta$G$_{1,2}$ (kcal $\cdot$ mol$^{-1}$)')
    plt.xlabel('$\Delta\Delta$G$_{1}$ + $\Delta\Delta$G$_{2}$ (kcal $\cdot$ mol$^{-1}$)')
    if ptype=='bind':
        s1 = [-4,0,4,8,12,16]
        plt.annotate('Binding (SKEMPIv2.0) [572]', xy=(-5.1, 14.8), size=0.9*fs)
    elif ptype=='fold':
        plt.annotate('Folding (ProTherm4) [204]', xy=(-8.7, 10.8), size=0.9*fs)
        s1 = [-8,-4,0,4,8,12]
    plt.xlim(tuple(lims))
    plt.ylim(tuple(lims))
    plt.plot(lims, lims, color='k')
    plt.xticks(s1)
    plt.yticks(s1)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(out+'.pdf', format='pdf', bbox_inches='tight')


#---run---#

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-fs',
            type    =   int,
            help    =   "Font size, defaults to 28",
            default =   28)
    parser.add_argument('-o',
            type=str,
            default='out')
    parser.add_argument('-i',
            type    =   str,
            help    =   "Input datafile")
    parser.add_argument('-ptype',
            type    =   str,
            help    =   "Type of plot: bind or fold")
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    
    try:
    	data = pd.read_csv(args.i, sep='\t')
    except:
        print("Error: Bad file specification")
        sys.exit(1)

    X = data['ddG1'] + data['ddG2']
    y = data['ddG12']

    cutoff = 0.5

    epistasis = y - X

    pos_ep_x = X[epistasis < -cutoff]
    neg_ep_x = X[epistasis > cutoff]
    pos_ep_y = y[epistasis < -cutoff]
    neg_ep_y = y[epistasis > cutoff]
    null_ep_y = y[(epistasis < cutoff)&(epistasis > -cutoff)]
    null_ep_x = X[(epistasis < cutoff)&(epistasis > -cutoff)]

    fitdata = linregress(X,y)

    def fit(x):
        return float(fitdata[0])*x + float(fitdata[1])

    fit_r = pearsonr(X,y)[0]

    if args.ptype == "bind":
        lims = [np.min([X,y])-1, np.max([X,y])+2]
    elif args.ptype == "fold":
        lims = [np.min([X,y])-3, np.max([X,y])+1.2]


    sub_x = np.arange(min(X), max(X), 0.05)

    plotter1(args.fs, args.o, args.ptype)
