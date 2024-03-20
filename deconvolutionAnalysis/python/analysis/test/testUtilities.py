#
# testUtilities.py
# shared test functions
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#


# https://docs.python.org/3/howto/logging.html#changing-the-format-of-displayed-messages
# https://stackoverflow.com/a/48996222
import logging
logger = logging.getLogger(__name__)

import os
import pathlib as pl
from analysis.metrics import CibersortResultsAsKWayClassifier


################################################################################
def getEvaluator( relativeRootPath : pl.Path) -> CibersortResultsAsKWayClassifier:
    
    #dataDir = os.path.join(self.pwd, "../python/analysis/test/data")
    #dataDir = os.path.join(self.pwd, "./data")
    dataDir = relativeRootPath.joinpath("data")
    # expectedFractionsPath = f'{dataDir}/expectedFractions.tsv'
    # resultPath            = f'{dataDir}/results.tsv'
    expectedFractionsPath = dataDir.joinpath('expectedFractions.tsv')
    resultPath            = dataDir.joinpath('results.tsv')

    logger.info(f'pwd: {os.getcwd()}')
    logger.info(f"expectedFractionsPath :{expectedFractionsPath}")

    # evaluator is a object that calculates classification statistics
    evaluator = CibersortResultsAsKWayClassifier(resultPath, expectedFractionsPath, verbose=True)

    return evaluator
