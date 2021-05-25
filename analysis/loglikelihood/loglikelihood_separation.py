#!/usr/bin/env python3

#---imports---#

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import argparse
import sys

#---functions---#

def dataload(deviationdata, outlier):
    '''Loads the data'''
    rawdata = np.loadtxt(deviationdata, skiprows=1)
    devdata = rawdata[::,0]
    sepdata = rawdata[::,1]
    if outlier:
        sepdatamax = np.where(sepdata == sepdata.max())[0][0]
        #removes outlier, this is specific to binding
        devdata = np.delete(devdata, sepdatamax)
        sepdata = np.delete(sepdata, sepdatamax)
    return devdata, sepdata

def loglikenull(x):
    '''generates log-liklihood for the null model'''
    return np.sum(np.log(norm.pdf(x, loc=np.mean(x), scale=np.std(x))))

def loglikealt1(x,rset, a = 1.7349999, alpha = 0.0327):
    '''generates the log-liklihood for the alternative model'''
    return np.sum(np.log(norm.pdf(x, loc=np.mean(x), scale=f(a,alpha,rset))))


def loglikeall(x,rset, a = 1.73499999, alpha = 0.0327):
    return np.sum(norm.logpdf(x, loc=np.mean(x), scale=f(a,alpha,rset)))

def f(a,alpha,x):
    '''first guess for alternative'''
    return a*np.exp(-alpha*x)


def parameterSearchRecursive(x,rset,init_a,delta_a,init_alpha,delta_alpha, numiter):
    '''
    Returns best-fit parameters
    '''
    if numiter == 0:
        return init_a, init_alpha
    else:
        a_set = np.arange(init_a - delta_a,init_a + delta_a, delta_a/12)
        alpha_set = np.arange(init_alpha - delta_alpha, init_alpha + delta_alpha, delta_alpha/10)
        likelihoodMatrix = np.zeros([a_set.size, alpha_set.size])
        for i in range(a_set.size):
            for j in range(alpha_set.size):
                likelihoodMatrix[i,j] = loglikeall(x,rset,a = a_set[i], alpha = alpha_set[j])
        bestlikelihood = np.where(likelihoodMatrix == np.nanmax(likelihoodMatrix))
        a_best_loc = bestlikelihood[0][0]
        alpha_best_loc = bestlikelihood[1][0]
        if a_best_loc == 0 or a_best_loc == len(a_set) - 1:
            print("WARNING: BEST FIT FOR \"A\" TAKEN FROM END OF ARRAY {:d}".format(numiter))
        if alpha_best_loc == 0 or alpha_best_loc == len(alpha_set) -1:
            print("WARNING: BEST FIT FOR \"ALPHA\" TAKEN FROM EDGE OF ARRAY {:d}".format(numiter))
        a_best = a_set[a_best_loc]
        alpha_best = alpha_set[alpha_best_loc]
        return parameterSearchRecursive(x,rset,a_best,delta_a * 0.5, alpha_best, delta_alpha * 0.5, numiter -1 )

def likelihoodSearchRecursive(x,rset,init_a,delta_a,init_alpha,delta_alpha, numiter, scalefactor=8):
    '''
    Returns maximum likelihood via recursion
    '''
    if numiter == 0:
        return np.sum(norm.logpdf(x,x.mean(),f(init_a, init_alpha,rset)))
    else:
        a_set = np.arange(init_a - delta_a,init_a + delta_a, delta_a/scalefactor)
        alpha_set = np.arange(init_alpha - delta_alpha, init_alpha + delta_alpha, delta_alpha/scalefactor)
        likelihoodMatrix = np.zeros([a_set.size, alpha_set.size])
        for i in range(a_set.size):
            for j in range(alpha_set.size):
                likelihoodMatrix[i,j] = loglikeall(x,rset,a = a_set[i], alpha = alpha_set[j])
        bestlikelihood = np.where(likelihoodMatrix == np.nanmax(likelihoodMatrix))
        a_best_loc = bestlikelihood[0][0]
        alpha_best_loc = bestlikelihood[1][0]
        if a_best_loc == 0 or a_best_loc == len(a_set) - 1:
            print("WARNING: BEST FIT FOR \"A\" TAKEN FROM END OF ARRAY {:d}".format(numiter))
        if alpha_best_loc == 0 or alpha_best_loc == len(alpha_set) -1:
            print("WARNING: BEST FIT FOR \"ALPHA\" TAKEN FROM EDGE OF ARRAY {:d}".format(numiter))
        a_best = a_set[a_best_loc]
        alpha_best = alpha_set[alpha_best_loc]
        return likelihoodSearchRecursive(x,rset,a_best,delta_a/scalefactor, alpha_best, delta_alpha/scalefactor, numiter -1, scalefactor=scalefactor)


def modelgenalt1(devdata, sepdata, a, alpha):
    devmean = devdata.mean()
    devstd = devdata.std()
    alt1 = np.random.normal(loc=devmean, scale=a*np.exp(-alpha*sepdata), size=len(sepdata))
    new_a,new_alpha = parameterSearchRecursive(alt1, sepdata, 1,2,0,0.5,16)
    return new_a, new_alpha

    
def loglikeratiotest(null,alt):
    return null - alt

def nullgen(mean, std, num):
    return np.random.normal(mean, std, num)

def modelgen3(devdata,sepdata,ada,alphadalpha,niter,scalefactor):
    devmean = devdata.mean()
    devstd = devdata.std()
    ai,da = ada.split(',')
    alphai,dalpha = alphadalpha.split(',')
    null1 = np.random.normal(loc=devmean, scale=np.std(devdata,ddof=1), size=len(devdata))
    loglikenull1 = loglikeall(null1, sepdata, a = np.std(null1,ddof=1), alpha = 0)
    loglikealt = likelihoodSearchRecursive(null1, sepdata, float(ai), float(da), float(alphai), float(dalpha), niter,scalefactor=scalefactor)
    lam = loglikenull1 - loglikealt
    return lam

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--datfile',
            type    =   str,
            help    =   "The datafile to load. First colunm should be epistasis and second separation distances")
    parser.add_argument('--nummodel',
            type    =   int,
            default =   1000,
            help    =   "The number of models to simulate")
    parser.add_argument('--out',
            type    =   str,
            default =   'output.txt',
            help    =   "Output file for simulated lambda")
    parser.add_argument('--alphadalpha',
            type    =   str,
            default =   '1,0.5',
            help    =   "Initial estimate for alpha and delta"
    parser.add_argument('--ada',
            type    =   str,
            default =   '1,3',
            help    =   "Initial guess for a and delta")
    parser.add_argument('--numiter',
            type    =   int,
            default =   18,
            help    =   "Number of recursive iterations")
    parser.add_argument('--scale',
            type    =   float,
            default =   8,
            help    =   "The fraction and number of estimates")
    parser.add_argument('--outlier',
            action  = 'store_true',
            help    =   "If outlier exists, will delete from dataset"
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    dev,sep = dataload(args.datfile, outlier=args.outlier)
    ai, da = args.ada.split(',')
    alphai,dalpha = args.alphadalpha.split(',')
    expnull = loglikenull(dev)
    expalt = likelihoodSearchRecursive(dev, sep, float(ai), float(da), float(alphai), float(dalpha), args.numiter, scalefactor=args.scale)
    exp = expnull - expalt
    #Generate set of likelihood ratios. NOTE: THis can return a ton of errors, this is most likely because you've already found the 
    #Optimal likelihood and all adjacent estimates are within error.
    lamset = [ modelgen3(dev, sep, args.ada, args.alphadalpha,args.numiter,scalefactor=args.scale) for i in range(args.nummodel) ]
    print("Exp: {:f}".format(exp))
    print("Min of sims: {:f}".format(min(lamset)))
    with open('exp.txt', 'w') as f:
        f.write("{:f}".format(exp))
    with open(args.out, 'w') as f:
        for lam in lamset:
            f.write("{:f}\n".format(lam))
    
