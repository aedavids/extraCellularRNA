#
# utilities.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

import pandas as pd
import pathlib as pl
import re
import shutil
###############################################################################
def loadCache(source : str, localCacheDir : str, verbose=False) -> pl.Path:
    '''
    reading large files over a NFS mount is slow. loadCache() will make
    cache directory if needed and copy the source file into the local cache 
    if it does not not already exist

    The gi.ucsc.edu phoenix cluster requires data be cached on local node
    
    arguments:
        localCacheDir : str
            path to local file cache directory
            example : localCacheDir="/scratch/aedavids/tmp"

        source : str
            path of file to copy to local cache

        verbose : bool
            default False. If True will print full path to file in localCache

    returns:
        pathlib.Path : path to local cache directory
    '''
    # we can not join, combine source if it start from the root of the file system
    tmpSource = source
    if source[0] == "/":
        tmpSource = source[1:]        
    
    localTargetPath = pl.Path(localCacheDir,  tmpSource)
    if verbose:
        print("localCachePath:\n{}\n".format(localTargetPath))
            
    localTargetPath.parent.mkdir(parents=True, exist_ok=True)

    if not localTargetPath.exists():
        #print("localTargetPath:{} does not exits".format(localTargetPath))
        #! cp $source $localTargetPath 
        shutil.copy(source, localTargetPath)
        
    return localTargetPath


###############################################################################
def urlify(s):
    '''
    useful function for converting plot titles to strings
    https://stackoverflow.com/a/1007615/4586180

    print(urlify("I can't get no satisfaction!"))
    Prints: I-cant-get-no-satisfaction"
    '''

    # Remove all non-word characters (everything except numbers and letters)
    s = re.sub(r"[^\w\s]", '', s)

    # Replace all runs of whitespace with a single dash
    s = re.sub(r"\s+", '-', s)
    
    return s


