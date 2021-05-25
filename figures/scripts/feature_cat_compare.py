#!/usr/bin/env python3

'''
Plot generator for epistasis. Written by Jonathan Barnes 2020.

This is meant to be the final version to handle all presence + distribution plots. 


For questions: j@jbsci.dev
'''


#---imports---#

import numpy as np
import pandas as pd
import os,sys,argparse
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc,rcParams

#---functions---#


def dataload(datafile, separator='\t'):
    '''Loads datafile in specific format'''
    return pd.read_csv(datafile, sep=separator)

def FeatureHandler(datafile, feature):
    '''
    Parses feature data to make plotting easier
    '''
    ep = [item for item in datafile['epistasis']]
    feat = [str(item) for item in datafile[feature]]
    #Make dictionary with categories
    feat_dict = {category : [] for category in np.unique(feat)}
    #Populate dictionary
    for i in range(len(ep)):
        feat_dict[feat[i]].append(ep[i])
    return feat_dict

def draw_gaussian_at(support, sd, height, xpos, ypos,ax=None,**kwargs):
    if ax is None:
        ax = plt.gca()
    gaussian = np.exp((-support ** 2.0) / (2*sd**2.0))
    gaussian /= gaussian.max()
    gaussian *= height
    return ax.plot(gaussian + xpos, support + ypos, **kwargs)


def draw_hist_at(data, scale, height, xpos, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    subhist = np.histogram(data, bins=30)
    binedges = subhist[1]
    bincenters = [np.average([binedges[i],binedges[i+1]]) for i in range(len(binedges)-1)]
    bincolors = ['blue' if bincenters[i] < 0 else 'red' for i in range(len(bincenters))]
    return ax.barh(bincenters, subhist[0]*scale, height, left=xpos, edgecolor=bincolors, **kwargs)
    

def featureErrorBarPlotter(data, feature, system, fs, feature_common_name, output, titleside, figsize, titleshift, rot=0, fsadj=0, stripnum=False):
    '''
    Plots the given features categories wrt. epistasis.
    Errorbars represent the spread of epistasis. Median is the plotted coordinate
    '''
    rcParams['xtick.labelsize'] =fs - fsadj
    rcParams['ytick.labelsize'] =fs
    rcParams['axes.labelsize'] =fs
    rcParams['legend.fontsize'] =fs
    rcParams['font.family'] = 'serif'
    plt.rc('text', usetex=True)
    featDict = FeatureHandler(data, feature)
    numFeats = len(featDict.keys())
    categories = list(featDict.keys())
    if stripnum:
        nonumcats = [cat.strip('0') for cat in categories]
        nonumcatsfeats = { cat.strip('0') : featDict[cat] for cat in categories}
    figsize_split = (float(figsize.split("x")[0]), float(figsize.split("x")[1]))
    fig = plt.figure(figsize=figsize_split)
    ax = fig.add_subplot(111)
    for i in range(len(categories)):
        cat_i = categories[i]
        cat_med = np.mean(featDict[cat_i])
        cat_min = abs(min(featDict[cat_i]) - cat_med)
        cat_max = max(featDict[cat_i]) - cat_med
        ax.errorbar(i, cat_med, yerr=np.array([[cat_min],[cat_max]]), fmt='ko',
                    capsize=5)
        draw_hist_at(featDict[cat_i], 2/(len(featDict[cat_i])), (abs(-cat_min - cat_max))/30, i, ax=ax, color='gray', alpha=0.8)
    ax.plot([-1,numFeats+1],[0,0],c='gray',alpha=0.4)
    ax.set_xlim([-0.5,numFeats-0.5])
    if rot == 0:
        ax.axes.set_xticks([i for i in range(numFeats)])
        if stripnum:
            ax.set_xticklabels(['{:s} ({:d})'.format(category, len(nonumcatsfeats[category])) for category in nonumcats])
        else:
            ax.set_xticklabels(['{:s} ({:d})'.format(category, len(featDict[category])) for category in categories])
    if rot != 0:
        ax.axes.set_xticks([i for i in range(numFeats)])
        if stripnum:
            xlabels = ax.set_xticklabels(['{:s} ({:d})'.format(category, len(nonumcatsfeats[category])) for category in nonumcats], rotation=rot, ha='right')
        else:
            xlabels = ax.set_xticklabels(['{:s} ({:d})'.format(category, len(featDict[category])) for category in categories], rotation=rot, ha='right')
    ax.set_xlabel(feature_common_name)
    if system == 'binding':
        if titleside == 'left':
            ax.set_title("Binding (SKEMPIv2.0) [572]", x=0.62+titleshift, y=0.94, fontdict={'fontsize': fs - 4, 'horizontalalignment' : 'right'})
        elif titleside == 'right':
            ax.set_title("Binding (SKEMPIv2.0) [572]", x=0.38+titleshift, y=0.94, fontdict={'fontsize': fs - 4, 'horizontalalignment' : 'left'})
    if system == 'folding':
        if titleside == 'left':
            ax.set_title("Folding (ProTherm4) [204]", x=0.62+titleshift, y=0.94, fontdict={'fontsize': fs - 4, 'horizontalalignment' : 'right'})
        elif titleside == 'right':
            ax.set_title("Folding (ProTherm4) [204]", x=0.38+titleshift, y=0.94, fontdict={'fontsize': fs - 4, 'horizontalalignment' : 'left'})
    ax.set_ylabel('$\epsilon$ (kcal/mol)')
    plt.tight_layout()
    plt.savefig(output+'.pdf', format='pdf')


def epistasis_split_cutoff(data, cutoff, norm=False):
    ep = data['epistasis']
    total = len(ep)
    additive = data[data['ep_abs'] <= cutoff]
    nonadditive = data[data['ep_abs'] > cutoff]
    pos = len(nonadditive[nonadditive['epistasis'] < 0])
    neg = len(nonadditive[nonadditive['epistasis'] > 0])
    add = len(additive)
    noadd = len(nonadditive)
    if norm:
        return add/total, noadd/total, pos/noadd, neg/noadd
    else:
        return add, noadd, pos, neg



def epistasis_presence_bar(data, cutoff_list, fs, systemlabel, output):
    rcParams['xtick.labelsize'] = fs
    rcParams['ytick.labelsize'] = fs
    rcParams['axes.labelsize'] = fs
    rcParams['legend.fontsize'] = fs
    rcParams['font.family'] = 'serif'
    rc('text', usetex=True)
    width = 2/ len(cutoff_list)
    fig, (ax,ax2) = plt.subplots(2,1,sharex=True,figsize=(15,15))
    ad = []
    no = []
    p = []
    n = []
    for i in range(len(cutoff_list)):
        add, noadd, pos, neg = epistasis_split_cutoff(data, cutoff_list[i], norm=True)
        ad.append(add)
        no.append(noadd)
        p.append(pos)
        n.append(neg)
    x = np.array([i for i in range(len(n))])
    add = np.array(ad)
    noadd = np.array(no)
    pos = np.array(p)
    neg = np.array(n)
    ax.bar(x-width/2, add, width, color='k', label = 'No Epistasis (additive)')
    ax.bar(x+width/2, noadd, width, color='#543A91', label = 'Epistasis (non-additive)')
    ax2.bar(x-width/2, pos, width, color='b', label = 'More stable than predicted by additivity')
    ax2.bar(x+width/2, neg, width, color='r', label = 'Less stable than predicted by additivity')
    if systemlabel == "binding":
        ax.annotate("Binding (SKEMPIv2.0)[572]", xy=(int(len(cutoff_list)/2)-2, 0.91), fontsize=fs)
    elif systemlabel == "folding":
        ax.annotate("Folding (ProTherm4)[204]", xy=(int(len(cutoff_list)/2)-1.75, 0.91), fontsize=fs)
    else:
        ax.annotate(systemlabel, xy=(int(len(cutoff_list)/2)-3.1, 0.905), fontsize=fs)
    ax.set_ylim([0,1])
    ax2.set_ylim([0,1])
    ax.legend(loc=1, prop={'size' : 24}, frameon=False)
    ax2.legend(loc=1, prop={'size' : 24}, frameon=False) 
    ax.set_ylabel('General Breakdown (norm.)')
    ax2.set_ylabel('Sign Breakdown (norm.)')
    ax2.axes.set_xticks([i for i in range(len(cutoff_list))])
    ax2.axes.set_xticklabels(['{:.1f}'.format(cut) for cut in cutoff_list])
    ax2.set_xlabel('Additivity Cutoff (kcal/mol)')
    plt.tight_layout()
    plt.savefig(output + '.pdf', format='pdf')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data',
            type    =   str,
            help    =   "Datafile to load")
    parser.add_argument('--feature',
            type    =   str,
            default =   'cplx_type2',
            help    =   "Feature to generate plot for")
    parser.add_argument('--system',
            type    =   str,
            default =   'binding',
            help    =   "Specify binding or folding, defaults to binding")
    parser.add_argument('--simplename',
            type    =   str,
            default =   'Complex Type',
            help    =   "Name for feature for labels, etc.")
    parser.add_argument('--output',
            type    =   str,
            default =   'output',
            help    =   "Output filename (no extension)")
    parser.add_argument('--presence_bar',
            action  =   'store_true',
            help    =   "Toggles presence bar output")
    parser.add_argument('--fs',
            type    =   int,
            default =   28,
            help    =   "Fontsize")
    parser.add_argument('--rot',
            type    =   int,
            default =   0,
            help    =   "X-axis label rotation")
    parser.add_argument('--fsadj',
            type    =   int,
            default =   0,
            help    =   "Minor font adjustment for titles")
    parser.add_argument('--cutoffmin',
            type    =   float,
            default =   0.1,
            help    =   "Minimum cutoff")
    parser.add_argument('--cutoffmax',
            type    =   float,
            default =   1.0,
            help    =   "Max cutoff")
    parser.add_argument('--titleside',
            type    =   str,
            default =   'left',
            help    =   "Which side for the title to be on")
    parser.add_argument('--cutoffstep',
            type    =   float,
            default =   0.1,
            help    =   "Stepping for cutoff")
    parser.add_argument('--figsize',
            type    =   str,
            default =   '8x8',
            help    =   "Figure size, given as XxY")
    parser.add_argument('--stripnum',
            action  =   'store_true',
            help    =   "Strips leading/trailing number for categories that have one (to force reference)")
    parser.add_argument("--titleshift",
            type    =   float,
            help    =   "Amount to shift title to correct",
            default =   0)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    data = dataload(args.data)
    if args.presence_bar is False:
        featureErrorBarPlotter(data, args.feature, args.system, args.fs, args.simplename, args.output, args.titleside, args.figsize, args.titleshift, args.rot, args.fsadj, args.stripnum)
    elif args.presence_bar:
        cutofflist = list(np.arange(args.cutoffmin, args.cutoffmax+args.cutoffstep, args.cutoffstep))
        epistasis_presence_bar(data, cutofflist, args.fs, args.system, args.output)




