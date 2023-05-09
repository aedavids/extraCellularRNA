'''
aedavids@ucsc.edu
3/16/23

partitionDataTask splits data file into parts. we control the file name
and create part file names that will sort. our scatter/gather workflow
pass Array[File] to our gather task. The elements are file paths not file names.
They may have quid. ie do not sort by file names

arguments
argv[1] = string: path to file containing a list of part file paths. one path per line
'''
import os
import sys

def main() :
    filePathFile = sys.argv[1]
    f = open(filePathFile, "r")
    filePathList = f.readlines()

    d = dict()
    for path in filePathList:
        basename = os.path.basename(path)
        d[basename] = path.strip()

    sortedKeys = sorted(d.keys())
    for key in sortedKeys :
        print( d[key] )

########################################################################
if __name__ == '__main__':
    main()
