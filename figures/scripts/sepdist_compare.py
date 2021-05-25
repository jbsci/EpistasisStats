#!/usr/bin/env python



#-# Imports


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import os,sys,argparse
from matplotlib import rc, rcParams
import pandas as pd


#--# Functions 

def pairPropEpSplitter(epistasis_data, aux_property):
    '''
    Reads in epistasis data and whatever property you're comparing.

    Returns 6 arrays, matching positive, negative, and neutral epistasis values
    '''
    neg_ep = epistasis_data[epistasis_data > 0.5]
    neg_ep_aux = aux_property[epistasis_data > 0.5]
    ps_ep = epistasis_data[epistasis_data < -0.5]
    ps_ep_aux = aux_property[epistasis_data < -0.5]
    standardEp = epistasis_data[(epistasis_data > -0.5) & (epistasis_data < 0.5)]
    standardAux = aux_property[(epistasis_data > -0.5) & (epistasis_data < 0.5)]
    return neg_ep, neg_ep_aux, ps_ep, ps_ep_aux, standardEp, standardAux


def plotter1(fs, outfile):
    rcParams['xtick.labelsize'] = fs
    rcParams['ytick.labelsize'] = fs
    rcParams['axes.labelsize'] = fs
    rcParams['legend.fontsize'] = fs
    rcParams['font.family'] = 'serif'
    rc('text', usetex=True)
    bind_neg_y, bind_neg_x, bind_pos_y, bind_pos_x, bind_y, bind_x = pairPropEpSplitter(bindDevData, bindSepData)
    fold_neg_y, fold_neg_x, fold_pos_y, fold_pos_x, fold_y, fold_x = pairPropEpSplitter(foldDevData, foldSepData)
    fig, (ax1, ax2) = plt.subplots(2,1,sharex=True, figsize=(11,9))
    ax1.scatter(bind_x, bind_y, edgecolor='k', facecolor='None', label = 'No epistasis')
    ax1.plot(np.linspace(0,40,40), np.zeros(40), color='k')
    ax2.plot(np.linspace(0,40,40), np.zeros(40), color='k')
    ax1.scatter(bind_pos_x, bind_pos_y, edgecolor='blue', facecolor='None', label = 'More stable than additive')
    ax1.scatter(bind_neg_x, bind_neg_y, edgecolor='red', facecolor='None', label = 'Less stable than additive')
    ax1.legend(loc=4, prop={'size' : 20}, frameon=False, labelspacing=0, handletextpad=0)
    ax2.scatter(fold_x, fold_y, edgecolor='k', facecolor='None')
    ax2.scatter(fold_pos_x, fold_pos_y, edgecolor='blue', facecolor='None')
    ax2.scatter(fold_neg_x, fold_neg_y, edgecolor='red', facecolor='None')
    ax2.set_xlabel(r'$C_{\alpha}$ Separation (\AA)')
    ax1.set_ylabel(r'$\epsilon$ (kcal/mol)')
    ax1.set_yticks([-6, -4, -2, 0,2, 4, 6])
    ax2.set_yticks([-4, -2, 0, 2, 4])
    ax2.set_ylabel(r'$\epsilon$ (kcal/mol)')
    ax1.annotate('Binding (SKEMPIv2.0) [572]', xy=(20.4,4.8), fontsize=fs-3)
    ax2.annotate('Folding (ProTherm4) [204]', xy=(22,ymax-1.1), fontsize=fs-3)
    fig.tight_layout()
    fig.savefig(outfile + '.pdf')

def plotter2(fs, outfile):
    plt.rcParams.update({'font.size' : fs})
    fig, ax = plt.subplots(1,1, figsize=(9,7))
    ax.scatter(bindSepData, bindDevData, c=bindZ, edgecolor='', facecolor='None')
    ax.set_xlabel(r'$C_{\alpha}$ Separation (\AA)')
    ax.set_ylabel(r'$\epsilon$ (kcal$\cdot$mol$^{-1}$)')
    ax.set_yticks(np.round(np.linspace(-6, 5, nticks),2))
    ax.annotate('Binding (SKEMPIv2.0) [572]', xy=(25.5,2))
    fig.tight_layout()
    plt.savefig(outfile + '.pdf', format='pdf')


#---# Run


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fs',
            type=int,
            default=32)
    parser.add_argument('--outfile',
            type    = str,
            default =   "Fig2")
    parser.add_argument('--bind',
            type    =   str,
            help    =   "Path to binding data")
    parser.add_argument("--fold",
            type    =   str,
            help    =   "Path to folding data")
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    bindData = pd.read_csv(args.bind, sep='\t')
    foldData = pd.read_csv(args.fold, sep='\t')

    bindDevData = bindData['epistasis']
    bindSepData = bindData['separation']
    foldDevData = foldData['epistasis']
    foldSepData = foldData['separation']


    ymax = max([bindDevData.max(), foldDevData.max()])
    ymin = min([bindDevData.min(), foldDevData.min()])

    nticks = 9

    plotter1(args.fs, args.outfile)
