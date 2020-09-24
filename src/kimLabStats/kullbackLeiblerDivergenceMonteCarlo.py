'''
Created on Jun 23, 2020

@author: andrewdavidson
Copyright 2020 Santa Cruz Analytics. All rights reserved.
'''

import logging
import numpy as np
import pandas as pd
from   scipy.stats import entropy

################################################################################
class KullbackLeiblerDivergenceMonteCarlo( object ):
    '''
    A monte carlo simulation for estimating the p-value of a the Kullback Leibler Divergence,
     (cross entropy)
    
    public functions:
        __init__( self, numIteration  )
        permutationTest( self, probDist1, probDist2)
    
    over view of algo:
        probDist1 and probDist1 are discrete probability distributions 
        
        Our null hypothesis is that the there is no difference between the 
        2 distributions
        
        they are converted into data sets. i.e. sample counts that match the 
        probability distribution
        
        These data sets are then combined into a single population.
        
        the simulation randomly shuffles the combined data set into 2 samples
        and converts these sample back to probability distributions
        
        The Kullback Leibler Divergence is calculated
        
        If the kl is greater than the kl of the original probability 
        distributions, the count variable is increment
        
        the returned p-value is the count / number of iterations
        
        To accurately recreated the data from a probability distribution, 
        we assume the sample size is equal to 10**max( number of significant digits)
        How ever if you have probabilities that are irrational numbers like 1/3
        it will cause large amounts of memory to be used.
        
        set the numSignifigantDigits to control for this. 
        
    best practices:
        choose numSignifigantDigits wisely. The larger this value the more accurate
        the simulation will be how ever it will use much more memory. This value effects        
    '''

    logger = logging.getLogger( __name__ )

    ################################################################################
    def __init__( self, numIteration, numSignifigantDigits=3 ):
        '''
        
        arguments:
            numIteterations: 
                type int. The number of simulations to run
                
            numSignifigantDigits:
                type int. default = 3. see class doc for details
        '''
        self.numIteration = numIteration
        self.numSignifigantDigits = numSignifigantDigits

    ################################################################################
    def permutationTest( self, probDist1, probDist2):
        '''
        arguments
            x1List, x2List:
                array type objects. They do not have to be the same length. might have to be numpy

        returns:
            (cross entropy, p-value)

        ref:
            Lec 10, BME-263 Data Visualization
            https://en.wikipedia.org/wiki/Resampling_(statistics)#Permutation_tests
        '''

#         x1List = self._clean( probDist1 )
#         x2List = self._clean( probDist2 )

        assert len(probDist1) == len(probDist2), "error: random variables drawn from prob 1 or 2 must have the same number of states"
       
        # calculate the Kullback Leibler Divergence reference statistic
        # if probDist1 == probDist2 
        #    crossEntropy ==  entropy
        # else 
        #    cross entropy > entropy
        # kl = cross entropy - entropy
        
        kl = entropy(probDist1, probDist2, base=2)
        
        # see class doc for rounding explanation
        probDist1 = probDist1.round( self.numSignifigantDigits )
        probDist2 = probDist2.round( self.numSignifigantDigits )
        
        p1Counts = self._probDistributionToCounts(probDist1) 
        p1Items = self._countsToItems(p1Counts)
        
        p2Counts = self._probDistributionToCounts(probDist2) 
        p2Items = self._countsToItems(p2Counts)    
        
        # you can think of combinedList as if it was a bag of
        # colored marbles we want to randomly select from
        combinedList = np.concatenate( (p1Items, p2Items)  )  

        bigger = 0
        probDist1N = len(probDist1)
        for iteration in range( 0, self.numIteration, 1 ):
            newProbDist1, newProbDist2 = self._shuffle(combinedList, probDist1N)

            # cross entropy is equal to zero if the distributions are the same
            # else a negative number
            simulatedKL = entropy(newProbDist1, newProbDist2, base=2)
            if simulatedKL > kl:
                bigger += 1
            
        pValue = bigger / ( float( self.numIteration ) )
#         print("aedwip kl:{}".format(kl))
#         print("aedwip pvalue:{}".format(pValue))
        return (kl, pValue)

    #
    # private functions
    #
    
    ################################################################################
    def _clean(self, x):
        x = x[~np.isnan(x)]
        x = x[ np.logical_and( (x != np.inf), (x != -1* np.inf) )  ]
        return x
        
    ################################################################################    
    def _countsToItems(self, countNP):
        '''
        You can think if this function taking in histogram data and returning a list r
        representing a "bag of marbel" representing the original observed data
        
        assume we started with a random variable that can be in one of 3 states
        
        countNP is a numpy array of length 3. The items are the number of time the random
        variable was observered in a given state
        
        returns
            a numpy array. length == np.sum(countNP)
            
        example :
            countNP = [1,2,3]
            return = [1, 2, 2, 3, 3, 3]
        '''
        ret = np.zeros(countNP[0], dtype=int) + 1 # zero would mean event never occures. by def not valid prob
        for i in range(1, len(countNP)):
            #print("aedwip countNP[{}]:{}".format(i, countNP[i]))
            c = np.zeros(countNP[i], dtype=int) + i + 1
            ret = np.concatenate( (ret, c) )
            
        return ret
    
#     # TODO: AEDWIP: add unit test
#     testProb = np.array([0.1, 0.4, 0.1, 0.1, 0.3])
#     
#     testCounts = probDistributionToCounts(testProb) 
#     print("testCounts:{}".format(testCounts))
#     
#     testResults = countsToItems(testCounts)
#     print("\ntestResults:{}\n".format(testResults))
#     assert len(testResults) == 10
#
#
# as a test see if we can reconstruct the original prob distribution
#
# 
#     # groupby does not work on type int
#     testDF = pd.DataFrame( {"a": testResults} )
#     testDF["b"] = testDF["a"].astype(str)
#     binNP = testDF.groupby("b").count()
#     
#     testNP = binNP.to_numpy().flatten()
#     assert (testProb == (testNP / np.sum(testNP))).all()    
    
    ################################################################################
    def _fixCounts(self, d1, d2):
        '''
        it is possible after we randomly shuffle the combined list and
        split into to 2 sets of observed values, that the sets do not contain all the
        possible observed states of the random variable. This may happen becuase
        The cardinality of the observed values is small or a particular state does not
        have many observations. 
            
        return:
            2 pandas data frames corresponding to d1, and d2
            The indices have been aligned
            psudo count has been added
            
            example of output data frame 
                d1 = [1, 4, 3, 2, 4, 3, 1, 2, 5, 4]
                d2 does not have any observed value for state 4
                d2 = [2, 2, 5, 5, 3, 2, 3, 1, 3, 3]
            
                        count
                state      
                1       2
                2       4
                3       5
                4       1
                5       3
                
        '''
        df1 = pd.DataFrame( {"count":d1} )
        df2 = pd.DataFrame( {"count":d2} )
        
        # groupby does not work on type in
        df1['state'] = df1["count"].astype(str)
        df2['state'] = df2["count"].astype(str)
        
        d1CountsDF = df1.groupby("state").count()
        d2CountsDF = df2.groupby("state").count()
        
        # test to see if any of the possible observed
        # states are missing
        gby1Idx = d1CountsDF.index.to_numpy()
        gby2Idx = d2CountsDF.index.to_numpy()
        
        missingIdxs1 = np.setdiff1d(gby1Idx, gby2Idx)
    #     print("missingIdx1:{}".format(missingIdxs1))
        missingIdxs2 = np.setdiff1d(gby2Idx, gby1Idx)
    #     print("missingIdx2:{}".format(missingIdxs2))
        
        missingIdxs = np.unique( np.concatenate((missingIdxs1, missingIdxs2)) )
    #     print("missingIdxs:{}".format(missingIdxs))
        
        gby1IdxSet = set( gby1Idx )
        gby2IdxSet = set( gby2Idx )
        
    #     print("gby1IdxSet:{}".format(gby1IdxSet))
    #     print("gby2IdxSet:{}".format(gby2IdxSet))
        
        isSorted1 = True
        isSorted2 = True
        for missingIdx in missingIdxs:
            if missingIdx not in gby1IdxSet:
                d1CountsDF.loc[missingIdx] = [0]
                isSorted1 = False
                #print("adding idx:{} to d1".format(missingIdx))
                
            if missingIdx not in gby2IdxSet:
                d2CountsDF.loc[missingIdx] = [0]
                isSorted2 = False
                #print("adding idx:{} to d2".format(missingIdx))
    
                
        # make sure the counts for each state are aligned correctly
        if not isSorted1:
            d1CountsDF.sort_index(inplace=True)
            
        if not isSorted2:
            d2CountsDF.sort_index(inplace=True)        
    
        # add psudo counts
        d1CountsDF = d1CountsDF + 1
        d2CountsDF = d2CountsDF + 1
    
        return (d1CountsDF, d2CountsDF)
    
#     # TODO: AEDWIP add unit test test
#     d1 = [1, 4, 3, 2, 4, 3, 1, 2, 5, 4]
#     # d2 does not have any observation of random varible in state 4
#     d2 = [2, 2, 5, 5, 3, 2, 3, 1, 3, 3]
#     
#     gb1DF, gb2DF = fixCounts(d1, d2)
#     
#     expectedDF = pd.DataFrame( {'state':[ str(i) for i in [1, 2, 3, 4, 5] ], 
#                                  'count':[2, 4, 5, 1, 3]},
#                              )
#     expectedDF.set_index('state', inplace=True)
#     
#     assert expectedDF.equals(gb2DF)    
    
    ################################################################################    
    def _probDistributionToCounts(self, probDist):
        '''
        arguments 
            propDists.
                numpy array like
                
        return an numpy array of same length. where
            probDist = ret/np.sum(ret) 
            
        best practices:
            over view of algo:
                n = max number of signifigant digits
                counts = probDist * n
                
            this accurately recreates a data set matching the probability distribution.
            How ever if you have probabilities that are irrational numbers like 1/3
            it will cause large amounts of memory to be used.
            
            It is recommend that you round the probability distribution to a reasonable
            number of signifigant digits before calling probDistributionToCounts()
        '''
        # find the number of decimal places
        # subtract 2 to compenstate for '0.' 
        p = np.max( [len(s) - 2 for s in probDist.astype(str)] )
        
        # choose number of counts so that fraction are accurate
        # E.G. if prob i specificed to 1/1000 n = 100 would not
        # accuratly represent the original distribution
        # if n = 1000, our counts are larger than needed
        n = 10**p
        ret = (probDist * n).astype(int)
        return ret

#     TODO: AEDWIP add unit test
#     expected = np.array([0.05, 0.001, 0.35, 0.2, 0.25, 0.149])
#     testCounts = probDistributionToCounts(expected) 
#     print(testCounts)
#     testResults = testCounts/np.sum(testCounts)
#     assert (expected == testResults).all()

    ################################################################################    
    def _shuffle(self, combinedList, probDist1N):
        
        
        # randomly select 2 sets of marbels same size as p1Items and p2Items
        np.random.shuffle( combinedList )
        newList1 = combinedList[:probDist1N]
        newList2 = combinedList[probDist1N:]
        
    #     print("\n newilsts")
    #     print("newList1:{}".format(newList1))
    #     print("newList2:{}".format(newList2))
        
        # it is possible that the new lists do not have counts
        # for all possible states. We will fix this bellow
        groupByCountDFList = self._fixCounts(newList1, newList2)
        
        # convert the sets back into probablity distributions
        ret = [0, 1]
    #     groupByCountDFList = [gb1DF, gb2DF]
        for i in range(len(groupByCountDFList)):
    #         print("******* {}".format(i))
            DF = groupByCountDFList[i]
            
    #         DF = pd.DataFrame( {'a':itemsList[i]} )
    #         DF["b"] = DF['a'].astype(str) # group by does not work on int
    #         print("DF:\n{}".format(DF))
    #         #print(DF.info())
            
    #         gCount = DF.groupby('b').count()
    #         #print(type(gSum))
    #         print("\ngCount:{}".format(gCount))
            
            
    # #         tNP = gCount.index.astype(int).to_numpy()
    #         tNP = gCount.to_numpy().flatten()
    #         print()
    
            countsNP = DF.to_numpy().flatten()
    #         print("countsNP:{}".format(countsNP))
            
            retProb = countsNP / np.sum(countsNP)
            ret[i] = retProb
    #         print()
    #         print("aedwip restProb:{}".format(retProb))
    #         print()
            
        return ret
    
    # TODO: AEDWIP add unit test test
    #np.random.seed(42) # 42 trigger need for fixCounts()
#     testP1 = np.array( [0.1, 0.3, 0.4, 0.1, 0.1] )
#     testP2 = np.array( [0.2, 0.2, 0.2, 0.2, 0.2] )
#     result1, result2 = monteCarloProbSuffle( testP1, testP2 )
#     print()
#     print("p1 :{}\nret:{}".format(testP1,result1))
#     print("p2 :{}\nret:{}".format(testP2,result2))
#     
#     assert len(testP1) == len(result1) == len(result2)
#     assert 1.0 == np.sum(result1)
#     assert 1.0 == np.sum(result2)
