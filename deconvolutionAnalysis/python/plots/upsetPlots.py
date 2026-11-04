#
# upsetPlots.py
#   Creates upsetplots and supporting plot data
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref : extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
#

import logging
from   matplotlib import pyplot as plt
import numpy as np
import os
import pandas as pd
from   pathlib import Path
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration
import pprint as pp
import upsetplot as upsp

###############################################################################
class UpsetPlot( object ):
    '''
    TODO

    public functions: see doc string for details
        __init__()
        configurablePlot()
        findIntersectionElements()
        saveInteresection()        
        savePlot()
    '''

    logger = logging.getLogger(__name__)

   ################################################################################    
    def __init__(self,
                 signatureGeneConfig : SignatureGeneConfiguration,
                 geneSetsDict : dict,
                 geneSetsUpsetPlotData : pd.DataFrame,
                 ):
        '''
        TODO

        arguments:
            signatureGeneConfig:
                contains plot title, base name, ...

            getSetDict:
                data used by findIntersectionElements()

                key : data set . Example "WholeBlood"
                value: list of element ids to plots. Ex. gene names

            geneSetsUpsetPlotData:
                data for drawing upset plot. 
                ref: UpsetPlotDataFactory

                pandas multilevel index data frame. created from a function like 
                upsetPlot from_contents(). 
        '''
        self.signatureGeneConfig = signatureGeneConfig
        self.geneSetsDict = geneSetsDict
        self.geneSetsUpsetPlotData = geneSetsUpsetPlotData
        self.fig = None

    ################################################################################    
    def configurablePlot(self, **kwags) -> tuple[ plt.figure, dict]:
        '''
        return 
            (fig, pltDict)
        
        allows you to pass extra upset plot variables
        examples
            configurablePlot(geneSetsUpsetPlotData, show_counts=True, min_degree=4)
            configurablePlot(geneSetsUpsetPlotData, show_counts=True, min_degree=4, max_degree=7)

        ref: plots.test.testUpsetPlots.testPlot()

        '''
        pltDict = upsp.UpSet(self.geneSetsUpsetPlotData,  **kwags).plot()
        self.fig = plt.gcf()
        title = self.signatureGeneConfig.title
        self.fig.suptitle( title, fontsize=40 ) 
        
        return (self.fig, pltDict)

    ################################################################################  
    def findDegrees(self, intersectionDict) -> list[int]:
        '''
        By defaults the upset plot shows all intersections. This can create very big
        plots that are hard to understand. The upset module has some confusing terms

        count: the number of elements in an intersection. 

        degree: similar. It is the number of sets the intersection is formed from
            example: the degree of A.intersect(B).intersect(C) is 3. The count is the
            number of elements in this intersection

        when you call configurePlot() you can pass arguments optional arguments to
        control which sets are ploted
        
        example:
            configurablePlot(geneSetsUpsetPlotData,  min_degree=4, max_degree=7)
            will create a plot that only contains sets with counts = 4, 5, 6, or 7

        arguments:
            intersectionDict: a dictionary returned by findIntersectionElements()

        returns:
            a list of integers. The integers are the sorted intersection sizes. 
            use these as min_degree and max_degree
        '''
        self.logger.info("BEGIN")

        s = set() 
        for key,items in intersectionDict.items():
            # example key: ('Vagina', 'Whole_Blood')
            l = len(key)
            s.add(l)
        
        ret = sorted( list(s) )


        self.logger.info("END")
        return ret


    ################################################################################    
    def findIntersectionElements(self) -> dict:                                 
        '''
        The upset plot displays set level information about interesection. i.e. the 
        sets that are part of each intersection.

        We need to know what elements are in each intersetion to be able to tune our models

        returns: dict
            retDict
                key: a tuple of sorted  set names
                value: the name of the elements in the intersection

        ref: plots.test.testUpsetPlots.testIntersection()
        '''
        self.logger.info("BEGIN")
        setNamesList = self.geneSetsUpsetPlotData.index.names

        # use numpy fancy indexing
        setNamesNP = np.asarray(setNamesList)
        self.logger.info(f'type(setNamesNP) {type(setNamesNP)}')
        #self.logger.info(f'setNamesNP : {setNamesNP} setNamesNP[True,False, True] :{setNamesNP[[True,False, True]]} ')

        # use groupby level to aggregate by multilevel index instead of columns
        groupByDF = self.geneSetsUpsetPlotData.groupby( level=setNamesList)
        #token = "_XXX_"
        retIntersectionDict = dict()
        for key, item in groupByDF:
            # key is tuple of booleans
            self.logger.debug(f'key: {key}')

            # combine all the intersection set names into a single key
            # use numpy fancy indexing
            # intersectionSetNames = token.join( setNamesNP[ list(key) ] )
            intersectionSetNames = tuple( sorted(setNamesNP[ list(key) ]) )
            self.logger.debug(f'intersectionSetNames: {intersectionSetNames}')
            self.logger.debug(f'get_group({key}) \n{groupByDF.get_group(key)}, "\n\n"')
            ggDF = groupByDF.get_group(key)
            # self.logger.info(f'type(ggDF) : {type(ggDF)} ggDF.columns : {ggDF.columns}')
            # self.logger.info(f'aedwip : {ggDF["id"]} aedwip')
            self.logger.debug(f'ggDF["id"].tolist() : {ggDF["id"].tolist()}')

            retIntersectionDict[intersectionSetNames] = ggDF["id"].tolist()
            
        self.logger.info("END")
        return retIntersectionDict
    
    ################################################################################    
    def saveInteresection(self,  
                            outdir : str,
                            intersectionElementsDict : dict
                            ) -> str :
        '''
        example of use

        ```
        up = UpsetPlotDataFactory(...)
        intersectionElementsDict, retSingleSetDict  = up.findIntersectionElements()

        # update concatenates dictionaries
        # it Adds key:value elements to the dictionary.
        intersectionElementsDict.update(retSingleSetDict)
        up.saveInteresection("outdir", intersectionElementsDict)
        ```
        outdir :
                director to write to. Will create if it does not exist

        returns the path the output file
        '''
        self.logger.info("BEGIN")
        baseName = self._getBaseName( self.signatureGeneConfig )  

        os.makedirs(outdir, exist_ok=True)

        # https://www.adamsmith.haus/python/answers/how-to-read-a-dictionary-from-a-file-in--python
        fileName =  baseName + ".intersection.dict"
        # filePath = aedwip    self.signatureGeneConfig.localCacheRoot + "/" + fileName
        filePath = os.path.join(outdir, fileName)
        with open(filePath,'w') as dataFile: 
            #dataFile.write(str(intersectionElementsDict))
            prettyStr =  pp.pformat(intersectionElementsDict, indent=4, sort_dicts=True)
            dataFile.write(prettyStr)
        
        # intersectionURL = self.signatureGeneConfig.saveGenesOfInterestToBucketURL() + fileName
        # print("save to :\n{}".format(intersectionURL))
        # ! gsutil cp $filePath $intersectionURL 

        self.logger.info(f'save to: {filePath}')
        self.logger.info("END")

        return filePath


    ################################################################################    
    def savePlot(self, 
                 outdir : str,
                 extraFileNameParts : str
                 ) -> str :
        '''

        arguments:
            outdir :
                director to write to. Will create if it does not exist

            extraFileNameParts: 
                a string to base output file name. Use to insure
                filename is unique

                ex.
                extraFileNameParts = 'min_degree=2,max_degree=3'

        returns:
            full path to image in png format

        raises (throws ) FileExistsError

        '''
        baseName = self._getBaseName( self.signatureGeneConfig )
        baseName = baseName + "-" + extraFileNameParts

        os.makedirs(outdir, exist_ok=True)
        
        # save png
        fileName = baseName + ".png"
        #filePath =  self.signatureGeneConfig.localCacheRoot + "/" + fileName
        filePath =  os.path.join(outdir, fileName)

        if os.path.exists(filePath):
            raise FileExistsError(filePath)  

        self.logger.info(f'saving to {filePath}')
        self.fig.savefig(filePath, dpi=300, bbox_inches='tight', facecolor='white')

        # imgURL = self.signatureGeneConfig.saveGenesOfInterestToBucketURL() + fileName
        # print("save to :\n{}".format(imgURL))
        # ! gsutil cp $filePath $imgURL  

        return filePath  

    #
    # private functions
    #
    ################################################################################    
    def _getBaseName(self, signatureGeneConfig):
        # remove chars that can not be used in file names
        baseName =  self.signatureGeneConfig.title.replace(" ", "-")
        baseName = baseName.replace("<", "lt")
        baseName = baseName.replace(">", "gt")
        return baseName
