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

    #sklearn label encoder internal types are np.str() and np.int64()
    #convert to python types    
    saveDict = {str(label): int(code) for label, code in zip(encoder.classes_, encoder.transform(encoder.classes_))}
    saveDictionary(path, saveDict)

################################################################################
def saveLabelEncoderDepreciated(path : str,
                     encoder : LabelEncoder):
    '''
    saves encoder as a dictionary
    weird after a python update this stopped working

    $ cat labelEncoder.dict 
{   np.str_('Colorectal Cancer'): np.int64(0),
    np.str_('Esophagus Cancer'): np.int64(1),
    np.str_('Healthy donor'): np.int64(2),
    np.str_('Liver Cancer'): np.int64(3),
    np.str_('Lung Cancer'): np.int64(4),
    np.str_('Stomach Cancer'): np.int64(5)}
    '''

    saveDict = encoder2Dict(encoder)
    saveDictionary(path, saveDict)

################################################################################
def encoder2DictDeprecicated(encoder : LabelEncoder) -> dict  :
    '''
    key is class
    value is int
    '''
    values = encoder.transform(encoder.classes_)
    retDict = dict(zip(encoder.classes_, values))
    return retDict
