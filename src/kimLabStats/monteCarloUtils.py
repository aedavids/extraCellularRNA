'''
Created on Sep 9, 2020

@author: Andrew Davidson
aedavids@ucsc.edu

utility functions for use with monteCarloUtils.py

public functions
    kolmogorovSmirnovStat(dist1, dist2)
    kullbackLeiblerStat(dist1, dist2)
    shuffleParallelList(l1, l2)
    
'''

import numpy as np
from   scipy.stats import entropy
from   scipy.stats import ks_2samp

################################################################################
def kolmogorovSmirnovStat(dist1, dist2):
    '''
    MonteCarlo wrapper around ks_2samp 
    arguments
        dist1, and dist2 are 2 'list like' probablity distrbutions
    
    returns 
    '''
    ksStat, pvalue = ks_2samp(dist1, dist2)
    return ksStat

################################################################################
def kullbackLeiblerStat(dist1, dist2):
    '''
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html

    MonteCarlo wrapper around entropy() 
    arguments
        dist1, and dist2 are 2 'list like' probablity distrbutions
    
    returns relative entropy in bits of information
    '''
    ret = entropy(pk=dist1, qk=dist2, base=2)
#     print("p1:{}\np2:{}".format(dist1, e2))

    return ret

# p1 = [0.35, 0.35, 0.10, 0.10, 0.04, 0.04, 0.01, 0.01]
# p2 = p1[::-1]
# p2 = [0.30, 0.30, 0.10, 0.10, 0.04, 0.04, 0.06, 0.06]
# e1 = entropy(p1, base=2)
# e2 = entropy(p2, base=2)
# print("e1:{} e2:{}".format(e1, e2))

# kl1 = entropy(pk=p1, qk=p2, base=2)
# kl2 = entropy(pk=p2, qk=p1, base=2)
# print("kl1:{} kl2:{}".format(kl1, kl2))
################################################################################
def shuffleParallelListDEPRECATED(l1, l2):
    '''
    this is a hack that is probably not correct. see KullbackLeiblerDivergenceMonteCarlo.py
    K-S test assume we have a distribution of a 1 dimentional random variable
    This is a hack so we can calculate biotype p-values using a permutation test
    
    you can consider l1 and l2 as if they are column vectors in a matrix
    will randomly shuffle values in each row
    
    arguments:
        list1 and list2 of equal length
        
    returns tuple with 2 lists
    
    TODO: create a vectorized implementation
    '''
    assert len(l1) == len(l2), "ERROR: l1 and l2 must have same lengths"
    l = len(l1)
    ret1 = np.zeros(l)
    ret2 = np.zeros(l)
    for i in range(l):
        # randomly choose list 1 or 2
        r = np.random.randint(2, size=1)
        if r == 0:
            ret1[i] = l1[i]
            ret2[i] = l2[i]
        else :
            ret2[i] = l1[i]
            ret1[i] = l2[i]
            
    return(ret1, ret2)

# test
# l1 = [0,0,0]
# l2 = [1,1,1]

# for i in range(5):
#     r1, r2 = shuffleParallelList(l1, l2)
#     print("{}\n{}\n".format(r1, r2))
    