#
# Andrew E. Davidson, aedavids@ucsc.edu
# 4/18/2020
#

set -x # turn debug on
# set + x # turn debug off

rootDir=`git rev-parse --show-toplevel`
cd ${rootDir}
source ~/anaconda3/etc/profile.d/conda.sh
conda activate `basename $PWD` 
export PYTHONPATH="${PYTHONPATH}:`pwd`/src"
#cd jupyterNotebooks
jupyter notebook