'''
Created on Jul 12, 2020

@author: aedavids@ucsc.edu
'''

import gseapy as gp
import logging


################################################################################
class GSEA( object ):
    '''
    classdocs
    '''
    logger = logging.getLogger( __name__ )

    ################################################################################
    def __init__( self, params ):
        '''
        Constructor
        '''
        pass
    
    ################################################################################
    def preRank(self):
        '''
        aedwip
            # to speed up debugging do not run if results have already been
        # calculated
     
        returns a tuple (cvsPath, preRes)
            if cvsPath exists preRes == None
        '''
        preRes = gp.prerank( rnk=rnkDf,
                         gene_sets='GO_Biological_Process_2018',
                         processes=4,
                         permutation_num=100,  # reduce number to speed up test
                         outdir=path, format='png' )
        ret = ( csvPath, preRes )

    ################################################################################
    def rankPathWaysForCluster( self, anndata, clusterSigs, clusterId ):
        '''
        aedwip
            # to speed up debugging do not run if results have already been
        # calculated
    
        returns a tuple (cvsPath, preRes)
            if cvsPath exists preRes == None
        '''
        rnkDf = pd.DataFrame( data={'gene':anndata.var['GeneName-0'],
                                   'score':clusterSigs[clusterId]} )
    
        base = './gseapy.out/prerank_report_GO_Biological_Process_2018_clusterId-'
        path = base + clusterId
        # print("path:{}".format(path))
    
        # path to output from previous run
        csvPath = path + "/" + 'gseapy.prerank.gene_sets.report.csv'
        # print("cvsPath:{}".format(csvPath))
    
        # to speed up debugging do not run if results have already been
        # calculated
        ret = None
        exists = os.path.isfile( csvPath )
        if exists:
            # gseapy.gsea.Prerank
            # type(ret):<class 'gseapy.gsea.Prerank'>
            # ret.outdir = csvPath
            ret = ( csvPath, None )
        else:
            preRes = gp.prerank( rnk=rnkDf,
                             gene_sets='GO_Biological_Process_2018',
                             processes=4,
                             permutation_num=100,  # reduce number to speed up test
                             outdir=path, format='png' )
            ret = ( csvPath, preRes )
    
        # print("type(ret):{}".format(type(ret)))
        return ret


    ################################################################################
    def rankPathWays( self, anndata, clusterSigs, topN=5 ):
        '''
        returns a panda dataframe with columns Term, nes, cluster id
            Term: pathway
            nes: normalized enrichment score
            cluster id: an integer
        '''
        retDF = pd.DataFrame()
        for clusterId in clusterSigs.keys():
            csvPath, preRes = rankPathWaysForCluster( anndata,
                                             clusterSigs,
                                             clusterId )
            df = pd.read_csv( csvPath )
            df2 = df.sort_values( by='nes', ascending=False )
            df3 = df2.iloc[0:topN, [0, 2]]
            numRows = df3.shape[0]
            cid = [int( clusterId ) for j in range( numRows )]
            df3['cluster id'] = cid
            retDF = retDF.append( df3 )
    
        retDF = retDF.sort_values( by='cluster id' )
        return retDF
