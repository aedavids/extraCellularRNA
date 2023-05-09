'''
aedavids@ucsc.edu
3/16/23

3 positional arguments

argv[1] = string: path to file you want to partition
argv[2] = int : maximum number of columns for a given partition
argv[3] = string: if present file assumed to be in TSV format else CSV
'''

from string import ascii_lowercase
import itertools

import pandas as pd
import sys

# https://stackoverflow.com/a/29351603/4586180
def iter_all_strings():
    for size in itertools.count(1):
        for s in itertools.product(ascii_lowercase, repeat=size):
            yield "".join(s)




def main() :
    dataFile = sys.argv[1]
    numSamplesInPart = int(sys.argv[2])

    sep = ','
    if len(sys.argv) > 3:
        sep = '\t'

    df = pd.read_csv(dataFile, sep=sep, index_col=0)
    #print(f'\n**********\ndf.index : \n{df.index.tolist()}\n')
    #print(f'\n**********\ndf.columns: \n{df.columns}\n');
    # print(f'len(df.columns) : { len(df.columns) }')
    # print('df:')
    # print(df)
    # print("\n\n****************** BEGIN\n\n")

    numCols = len(df.columns)    
    partGenerator = iter_all_strings()
    sortedPartsList = sorted( [ next(partGenerator) for i in range(numCols) ] )

    i = 0
    for start in range(0, numCols, numSamplesInPart) :
        #print("\n**************")

        end = start + numSamplesInPart
        if end > numCols :
            end = start + (numCols - start) 

        #partName = next(partGenerator)
        partName = sortedPartsList[i]
        i += 1
        print(f'start : {start} end: {end} partName : {partName}')
        partDF = df.iloc[:, start:end]
        partFileName = "part_" + partName + ".csv"
        print(f'partFileName : {partFileName}')
        print(partDF.head())
        partDF.to_csv(partFileName, sep=sep)


    #print(df.describe())

    #print("\n*************\n")
    #print(df.info())
        
########################################################################
if __name__ == '__main__':
    main()
