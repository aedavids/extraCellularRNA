'''
Created on Jun 23, 2020

@author: andrewdavidson
Copyright 2020 Santa Cruz Analytics. All rights reserved.
'''

import logging
import numpy as np
#from numpy.random.mtrand import shuffle


################################################################################
class MonteCarlo( object ):
    '''
    classdocs

    TODO: AEDWIP
    '''

    logger = logging.getLogger( __name__ )

    ################################################################################
    def __init__( self, numIteration, diffFunc, shuffleFunc=None  ):
        '''
        Constructor
        
        arguments:
            diffFunction:
                example:
                def median2SidedDiff( x1List, x2List ):
                    difference = np.abs( np.median( x1List ) - np.mean( x2List ) )
                    return difference
            
            shuffleFunc 
                default = None
                a function that takes in two list and returns two list
                
                default implementation concatenates xlist1 and xList2 and use 
                np.random.shuffle to split
        '''
        self.numIteration = numIteration
        self.diffFunc = diffFunc
        self.shuffleFunc = shuffleFunc

    ################################################################################
    def permutationTest( self, x1List, x2List):
        '''
        calculates the p-value for a  2 sided  tests

        arguments
            x1List, x2List:
                array type objects. They do not have to be the same length. might have to be numpy

        returns:
            p-value

        ref:
            Lec 10, BME-263 Data Visualization
            https://en.wikipedia.org/wiki/Resampling_(statistics)#Permutation_tests
        '''

        x1List = self._clean( x1List )
        x2List = self._clean( x2List )
        
        # difference = np.abs( np.std( x1List ) - np.std( x2List ) )
        difference = self.diffFunc( x1List, x2List )
        combinedList = np.concatenate( ( x1List, x2List ), axis=None )

        bigger = 0
        for iteration in range( 0, self.numIteration, 1 ):
            if self.shuffleFunc != None:
                newList1, newList2 = self.shuffleFunc(x1List, x2List)
            else :
                np.random.shuffle( combinedList )
                newList1 = combinedList[:len( x1List )]
                newList2 = combinedList[len( x1List ):]

            newDifference = self.diffFunc( newList1, newList2 )

            if newDifference >= difference:
                bigger += 1

        return bigger / ( float( self.numIteration ) )

    #
    # private functions
    #
    ################################################################################
    def _clean(self, x):
        x = x[~np.isnan(x)]
        x = x[ np.logical_and( (x != np.inf), (x != -1* np.inf) )  ]
        return x

        
