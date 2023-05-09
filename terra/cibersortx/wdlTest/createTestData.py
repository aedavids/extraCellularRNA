'''
aedavids@ucsc.edu
3/16/23
'''

import pandas as pd
import sys

def main() :
    n = int(sys.argv[1])
    idx = list()
    dataDict = dict()
    
    for i in range(n):
        col = [i]*n
        #print(col)
        dataDict['s' + str(i)] = col
        idx.append('g' + str(i))

    #print(dataDict)
    df = pd.DataFrame(dataDict, index=idx)
    df.index.name = "name" 
    print(df.to_csv())
        
########################################################################
if __name__ == '__main__':
    main()
