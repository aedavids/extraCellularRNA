#
# testUtilities.py
# common unit test functions
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

import pathlib as pl
################################################################################   
def _get1vsAllResults(relativeRootPath : pl.Path) -> list[str]:
        '''
        these are example of DESEq results file produced by 1vsAll.
        they have multiline headers and thousands of rows (they have not been filtered)
        '''
        retList = [
            str(relativeRootPath.joinpath("data/testSignatureGenes/1vsAll/UVM_vs_all.results")),
            str(relativeRootPath.joinpath("data/testSignatureGenes/1vsAll/Vagina_vs_all.results")),
            str(relativeRootPath.joinpath("data/testSignatureGenes/1vsAll/Whole_Blood_vs_all.results")),
        ]
        return retList