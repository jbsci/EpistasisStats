#!/usr/bin/env python


#-# Imports

import numpy as np
import os,sys,argparse
import matplotlib.pyplot as plt
from matplotlib import rc,rcParams

#--# Functions

def plotter1(sim,lm, fs, out, ptype):
    rcParams['xtick.labelsize'] = fs
    rcParams['ytick.labelsize'] = fs
    rcParams['axes.labelsize'] = fs
    rcParams['legend.fontsize'] = fs
    rcParams['font.family'] = 'serif'
    rc('text', usetex=True)
    plt.figure(figsize=(7.5,7.5))
    ultamin = min([loglam, min(simulated)])
    ht = plt.hist(sim, color='k')
    plt.ylabel('Frequency')
    plt.xlabel(r'Log-Likelihood Ratio ($\Lambda$)')
    plt.xlim(ultamin-1,0)
    plt.arrow(lm, 0.08*max(ht[0]), 0, -0.03*max(ht[0]), head_width=0.4, width=0.1, color='k', head_length=0.03*max(ht[0]))
    if ptype == "bind":
        plt.text(lm-0.4, 0.09*max(ht[0]), 'EXP', fontsize=fs)
        plt.text(ultamin-0.8, max(ht[0]), 'Binding (SKEMPIv2.0)', fontsize=fs-2)
    elif ptype=="fold":
        plt.text(lm-0.42, 0.09*max(ht[0]), 'EXP', fontsize=fs)
        plt.text(ultamin-0.8, max(ht[0]), 'Folding (ProTherm4)', fontsize=fs-2)
    plt.tight_layout()
    plt.savefig(out+'.pdf', format='pdf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-sim',
            type    =   str,
            help    =   "Path to textfile with simulated data")
    parser.add_argument("-exp",
            type    =   str,
            help    =   "Path to textfile with the experimental loglikelihood ratio")
    parser.add_argument("-o",
            type    =   str,
            help    =   "Output filename",
            default =   "out")
    parser.add_argument("-ptype",
            type    =   str,
            help    =   "Plot type: bind or fold")
    parser.add_argument("-fs",
            type    =   int,
            help    =   "Font size, defaults to 28",
            default =   28)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    simulated = np.loadtxt(args.sim)
    loglam = float(np.loadtxt(args.exp))
    plotter1(simulated,loglam, args.fs, args.o, args.ptype)

