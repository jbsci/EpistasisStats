#!/usr/bin/env python

import os
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc,rcParams

#-# functions #-#

def splitrsq(datafile,typing):
    '''
    Takes full details and strips out the delta R-sq data to a new file
    '''
    with open(datafile, 'r') as f:
        f = f.readlines()
        if typing=="bind":
            f = f[5:f.index("Call:\n")]
        elif typing == "fold":
            f = f[7:f.index("Call:\n")]
    if typing=="bind":
        combined = [elem[5:-3] for elem in f if elem[:3] == "[1]"]
    elif typing=="fold":
        combined = [elem[4:-3] for elem in f if elem[:3] == "[1]"]
    result_dict = { elem.split(':')[0] : float(elem.split(':')[1]) for elem in combined}
    return result_dict

def pull_coeff(datafile):
    '''Pulls coefficients'''
    with open(datafile, 'r') as f:
        f = f.readlines()
        coef = f[f.index("Call:\n")+9:]
    sep = coef[[i for i, elem in enumerate(ceof) if "separation " in elem][0]]


def keyfix(key):
    '''
    Fixes keys for abstracted features
    '''
    keylist = ["charge", "hp", "ss"]
    subkey = key.split("_")
    if subkey[0] in keylist:
        return subkey[0]
    elif len(subkey) > 1:
        if subkey[1] == "sasa":
            return "sasa"
        else:
            return key
    else:
        return key

def rankconvert(instance_dict):
    '''
    Converts list of delta r-square to rank order
    '''
    sorteddrsq = sorted(instance_dict.keys())[::-1][1:]
    asfeat = [instance_dict[elem] for elem in sorteddrsq]
    rankorder = { feat : asfeat.index(feat)+1 for feat in asfeat }
    return rankorder



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-dataroot', help="Data root direcotry")
    parser.add_argument('-outdir', help="output directory, defaults to current", default="./")
    parser.add_argument('-numruns', type=int, help="number of runs to process, defaults to 100", default=100)
    parser.add_argument('-system', help="System, either bind or fold", default="bind")
    parser.add_argument('-summary',help="Will write summary if specified", action="store_true")
    parser.add_argument('-fs', help="Font size", default=28, type=int)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    combodict_raw = { "{:02d}".format(i) : splitrsq("{:s}/run_{:02d}/run_{:02d}.txt".format(args.dataroot,i,i), args.system) for i in range(0,args.numruns)}
    combodict = {}
    for i in range(args.numruns):
        j = "{:02d}".format(i)
        combodict[j] = {keyfix(k) : v for k,v in combodict_raw[j].items()}

    avg_model_size = int(max([len(combodict["{:02d}".format(i)].keys()) - 1 for i in range(args.numruns)]))

    revcombodict = {}
    for i in range(args.numruns):
        j = "{:02d}".format(i)
        revcombodict[j] = {v : k for k,v in combodict[j].items()}

    rankorderdict = {}
    for i in range(100):
        j = "{:02d}".format(i)
        rankorderdict[j] = rankconvert(revcombodict[j])

    allkeys = []
    allkeys_total = {}
    for i in range(0,args.numruns):
        for key in combodict["{:02d}".format(i)]:
            if key not in allkeys:
                allkeys.append(key)
            if key not in allkeys_total.keys():
                allkeys_total[key] = 0
            allkeys_total[key] += 1

    value_dict = {}
    for feature in allkeys:
        value_dict[feature] = []
        for i in range(0,args.numruns):
            try:
                value_dict[feature].append(combodict["{:02d}".format(i)][feature])
            except:
                value_dict[feature].append(0)

    rankavg_dict = {}
    for feature in allkeys:
        rankavg_dict[feature] = []
        for i in range(args.numruns):
            try:
                rankavg_dict[feature].append(rankorderdict["{:02d}".format(i)][feature])
            except:
                pass
    rcParams['xtick.labelsize'] = args.fs
    rcParams['ytick.labelsize'] = args.fs
    rcParams['axes.labelsize'] = args.fs
    rcParams['legend.fontsize'] = args.fs
    rcParams['font.family'] = 'serif'

    f = plt.figure(figsize=(8,8))
    for key in allkeys:
        plt.plot([i for i in range(0,args.numruns)], value_dict[key], label="{:s}: {:d}%, dR^2 avg: {:.3f}".format(key, allkeys_total[key], np.average([elem for elem in value_dict[key] if elem != 0])))
    lgd = plt.legend(bbox_to_anchor=(-0.18,-0.92), loc='lower left')
    plt.xlabel("run")
    plt.ylabel(r"$\Delta R^{2}$")
    if args.system == "bind":
        plt.title("Binding (SKEMPIv2.0)", y=0.94, fontsize=args.fs)
    elif args.system == "fold":
        plt.title("Folding (ProTherm4)", y=0.94, fontsize=args.fs)
    plt.savefig(args.outdir + 'leave_10per_out_{:s}.pdf'.format(args.system), format='pdf', dpi=300, bbox_extra_artists=(lgd,), bbox_inches="tight")

    plt.clf()
    f = plt.figure(figsize=(8,8))
    #plot errorbar plot
    for i in range(len(allkeys)):
        plt.errorbar(i, np.average(value_dict[allkeys[i]]), xerr=0, yerr=np.std(value_dict[allkeys[i]], ddof=1), marker='o', capsize=7, elinewidth=5, capthick=3)
    plt.xticks([i for i in range(len(allkeys))], [allkeys[i] for i in range(len(allkeys))], rotation=75)
    plt.ylabel(r"$\Delta R^{2}$")
    plt.xlabel("Feature")
    plt.savefig(args.outdir + 'leave_10per_out_{:s}_errb.pdf'.format(args.system), format='pdf', dpi=300, bbox_inches="tight")
    #Provides summary data if desired
    if args.summary:
        with open(args.outdir + 'validation_data.txt', 'w') as f:
            f.write("Average Rank\n")
            for key in allkeys:
                f.write("{:s}\t{:.3f}\n".format(key, np.average(rankavg_dict[key])))
            f.write("\n")
            f.write("Average delta Rsq\n")
            for key in allkeys:
                f.write("{:s}\t{:.3f}\n".format(key, np.average([elem for elem in value_dict[key] if elem != 0])))
            f.write("\n")
            f.write("Number of models:\n")
            for key in allkeys:
                f.write("{:s}\t{:.3f}\n".format(key, allkeys_total[key]))
