'''
aedavids@ucsc.edu
3/16/23
'''

import pandas as pd
import sys

# https://stackoverflow.com/a/29351603/4586180
def iter_all_strings():
    for size in itertools.count(1):
        for s in itertools.product(ascii_lowercase, repeat=size):
            yield "".join(s)




def main() :
    csvFile = sys.argv[1]

    df = pd.read_csv(csvFile, index_col=0)
    sumSeries = df.sum()
    sumSeries.index.name = "sampleId"

    print(sumSeries.to_csv())
            
########################################################################
if __name__ == '__main__':
    main()
