#!/usr/bin/env python


#-# Imports

import os,sys,argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams


#--# Constants

a_bind = 1.758976999991862
alpha_bind = 0.03403662808736167
a_fold = 2.7185764312744105
alpha_fold = 0.039461975097656185
bindsepdata = (2.2887, 40)
binddevdata = (-5.77, 4.29425)
foldsepdata = (3.74274, 39.6303)
folddevdata = (-6.4, 12.6)
sigbind = 1.229754

bindingData = (a_bind, alpha_bind)
foldingData = (a_fold, alpha_fold)


#---# Functions


def draw_gaussian_at(support, sd, height, xpos, ypos, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    gaussian = np.exp((-support ** 2.0) / (2 * sd ** 2.0))
    gaussian /= gaussian.max()
    gaussian *= height
    return ax.plot(gaussian + xpos, support + ypos, **kwargs)

def bindplotgenAltNull(npoints, fs, out):
    rcParams['xtick.labelsize'] = fs
    rcParams['ytick.labelsize'] = fs
    rcParams['axes.labelsize'] = fs
    rcParams['legend.fontsize'] = fs
    rcParams['axes.titlesize'] = fs
    rcParams['axes.labelcolor'] = 'black'
    rcParams['font.family'] = 'serif'
    rc('text', usetex=True)
    fig, axes = plt.subplots(2,2, sharex=True, figsize=(7.5,7.5))
    xset = np.linspace(*bindsepdata, npoints)
    yset = np.linspace(*binddevdata, 400)
    sigma = a_bind * np.exp(-alpha_bind * xset)
    axes[1,0].set_xlabel(r'C$_{\alpha}$ Separation (\AA)')
    axes[1,1].set_xlabel(r'C$_{\alpha}$ Separation (\AA)')
    axes[0,0].set_ylabel(r'Distribution of $\epsilon$')
    axes[1,0].set_ylabel(r'$\sigma(\epsilon)$')
    axes[1,1].set_ylim(0.4,1.6)
    axes[0,0].set_title('Alternative Model')
    axes[0,1].set_title('Null Model')
    axes[0,0].yaxis.set_major_locator(plt.MaxNLocator(6))
    axes[0,1].yaxis.set_major_locator(plt.MaxNLocator(6))
    altcolor = (0.0,0.40,0.161)
    for i in range(len(xset)):
        draw_gaussian_at(yset, sd=sigma[i], height=-5.0, xpos=xset[i], ypos=0,ax=axes[0,0], color=altcolor)
        draw_gaussian_at(yset, sd=sigbind, height=-5.0, xpos=xset[i], ypos=0,ax=axes[0,1], color='b')
    axes[1,0].plot(xset, sigma, color=altcolor)
    axes[1,1].plot(xset, [sigbind for i in range(len(xset))], color='b')
    plt.tight_layout()
    fig.savefig(out+'.pdf', format='pdf', bbox_inches='tight')


#----# Run


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-fs',
            type    =   int,
            default =   28)
    parser.add_argument('-npoints',
            type    =   int,
            help    =   "Number points to show",
            default =   20)
    parser.add_argument('-o',
            type    =   str,
            help    =   "Output filename",
            default =   "output")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    bindplotgenAltNull(args.npoints, args.fs, args.o)



    
