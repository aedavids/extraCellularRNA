# aedwip

Andrew E. Davidson  
aedavids@ucsc.edu  
3/29/2024  

**createBestCreateSignatureGeneConfig.py**  
Find genes that that are statistically signifigant with  lfc <= -2.0 or >= 2.0
sorted by baseMean

This algorithm selects biomarkers with the strongest 'signal' not biomarkers that
are the most differentially expressed
        
**createBestLFCSignatureGeneConfig.py**  
Find genes that that are statistically signifigant with  lfc <= -2.0 or >= 2.0
and baseMean >= median( baseMean ) sorted by lfc 

This algorithm selects biomarkers that are the most differential expressed
with a strong 'signal' not biomarkers that have the strongest signal
        
**createBestCuratedGeneConfig.py**  
see BestCuratedGeneConfig.py for complete description

select top n genes sorted by base mean from the degree1 intersections

```
psudo code
- load interesectionDictPath
- degree1Dict = findIntersectionsWithDegree(interesectionDict, degree=1


- findGenes(deseqResults, fileName)
  if category in degree1Dict
      genes = degree1Dict[category]
      resultsDF = deseqResuts[genes]
      return resultsDF.sort_values(by'baseMean', accending=False).head(n)
  else 
      return empty dataframe
```
    
```
    Use Case 1: run our training pipeline on hand crafted gene signature matrix
        - interesection dictionary degree1 geneLists define the gene signature matrix
        we want to train with. It was edited by hand or create programmatically. 

        - set n to very large number like 5,000 to cause all the degree1 values to be returns

        - pipeline stage driver shell script config
            + deseqResultsDir contains the results files that will be based to findGenes()
        
            + example: extraCellularRNA/deconvolutionAnalysis/bin/1vsAll-~gender_category/best10CuratedDegree1.sh
                * rootDir="/private/groups/kimlab/GTEx_TCGA"
                * deseqResultsDir="${rootDir}/1vsAll"

    Use Case 2: Automatially add genes to a signature matrix defined by the interesectionDict
        - interesection dictionary degree1 geneLists define the gene signature matrix
        we want to train with. 

        - set n to the number of degree 1 values to return. Typically n is a small number

        example n = 3 will return at most 3 genes from the categories' degree 1 genes

        see pipeine stage dirver instructions above
```

**createBestRemoveHighDegreeSignatureGeneConfig.py**  
genes that are shared between many categories do not make good discrimnators

find best "n" genes. see BestSignatureGeneConfig.findGenes() for details
then remove genes that are elements of intersections with high degree.
The degree of a set is the number of set that have elements in the intersection

**createByDegreeSignatureGeneConfig.py**  
genes that are shared between many categories do not make good discrimnators

selects genes from intersections with degree = x

**createEnrichSignatureGeneConfig.py**  
try to insure each type has a degree1 intersection with a minium of number of genes

**createOptimalSelectiveEnrichSignatureGeneConfig.py**  
    genes that are shared between many categories do not make good discrimnators

find best "n" genes. see BestSignatureGeneConfig.findGenes() for details
then remove genes that are elements of intersections with high degree.
The degree of a set is the number of set that have elements in the intersection

**createRankSignatureGeneConfig.py**  
```
    overview:
    1. concat all deseq results files
    2. select rows that are biologically signigant
    3. sort by base mean
    4. rows with strongest biologically signigant signal will be on top
    5. itterate over the sorted list
        if we find a "new" gene assign it to category with ie highest base mean.
        i.e. signal strength

    Assigning genes to category based on signal strength may still produce poor
    discriminators. It is possible genes discovered this way are members of 
    intersections with high degree
```

**createSelectiveEnrichSignatureGeneConfig.py**  

ensures each category in categories has at least numberOfGenesToAdd degree 1 genes
