#
# Andrew E. Davidson, aedavids@ucsc.edu
# 4/18/2020
#

echo To start on a remote server pass arguments --no-browser --port=XXXX

set -x # turn debug on
# set + x # turn debug off

rootDir=`git rev-parse --show-toplevel`
cd ${rootDir}

#
# make script work for both anaconda and minconda
#

# $ conda info | grep -i 'base environment'
# base environment : /Users/andrewdavidson/anaconda3  (writable)
condaBase=`conda info | grep -i 'base environment' | cut -d : -f 2 | cut '-d ' -f 2`
# source ~/anaconda3/etc/profile.d/conda.sh
source ${condaBase}/etc/profile.d/conda.sh

conda activate `basename $PWD` 
export PYTHONPATH="${PYTHONPATH}:`pwd`/src"


# $@ == $1 $2, ...
jupyter notebook "$@"