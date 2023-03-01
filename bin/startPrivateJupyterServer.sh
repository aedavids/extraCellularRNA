#
# Andrew E. Davidson, aedavids@ucsc.edu
# 4/18/2020
#


set -x # turn debug on
# set + x # turn debug off

HOST_PORT=`findUnusedPort.sh`
echo "ssh tunnel port number: " $HOST_PORT

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
#jupyter notebook --no-browser --port $HOST_PORT "$@"
jupyter lab --no-browser --port $HOST_PORT "$@"
