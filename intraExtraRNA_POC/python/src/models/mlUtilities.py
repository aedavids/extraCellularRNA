#
# mlUtilities.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

from analysis.utilities import loadDictionary
from analysis.utilities import saveDictionary
import numpy as np
from sklearn.preprocessing import LabelEncoder

################################################################################
def loadEncoder(path: str) -> LabelEncoder:
    '''
    arguments:
        path: file containing labelEncoder values saved as a dictionary
    '''
    encoder = LabelEncoder()
    encoderDict = loadDictionary(path)

    # Manually assign the sorted list of class labels to the classes_ attribute
    # The keys of the dictionary are sorted according to their corresponding values
    # dictionary.get(key) returns the value value
    encoder.classes_ = np.array(sorted(encoderDict, key=encoderDict.get))

    return encoder

################################################################################
def saveLabelEncoder(path : str,
                     encoder : LabelEncoder):
    '''
    saves encoder as a dictionary
    '''

    saveDict = encoder2Dict(encoder)
    saveDictionary(path, saveDict)

################################################################################
def encoder2Dict(encoder : LabelEncoder) -> dict  :
    '''
    key is class
    value is int
    '''
    values = encoder.transform(encoder.classes_)
    retDict = dict(zip(encoder.classes_, values))
    return retDict
