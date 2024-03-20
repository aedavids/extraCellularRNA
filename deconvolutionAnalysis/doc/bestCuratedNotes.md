# BestCuratedGeneConfig Notes
Andrew E. Davidson  
aedavids@ucsc.edu  
1/10/24   

**<span style="color:red">See comments in bestCuratedGeneConfig.py. algo may not match discription bellow</span>**

## Big idea
The best discriminators will be genes from degree1 interesections. I.e. genes that are unique to a single category. If we where to search for the degree1 genes in the best 'n' results the larger n the more likely the degree 1 genes is a good discriminator. For example consider n = 1 all 83 genes will be degree 1. There is a good chance in diease cases like cancer that these genes will be differntial expressed in other cancer types.

keep in mind that the top/best 1 vs. all genes are not guaranteed to be unique. They are differentially signifigant comparted to the average of the remain 82 classes. 

our pipeline is as follows 

best500 -> best500FindAllDegree1_wl500 -> best10CuratedDegree1{1,2,3,5,10}

## algorithm overview

- **deconvolutionAnalysis/bin/1vsAll-~gender_category/best500**
for each category finds the best "n" genes. best is defined as genes that are biolgically signifigant. We sort the list in decending order and take the top "n"

- deconvolutionAnalysis/bin/1vsAll-~gender_category/best500FindAllDegree1_wl500.sh
finds all degree1 genes in the first 1000 use the "best" algorithm

  + run OptimalSelectiveEnrichSignatureGeneConfig
    * genes that are shared between many categories do not make good discriminator. 
    Find the best "n" genes see BestSignatureGeneConfig.findGenes() for details.
    Then remove genes that are elements of intersections with degree > 1.
    The degree of a set is the number of set that have elements in the intersection
    
  + deseqResultsDir results will be passed to OptimalSelectiveEnrichSignatureGeneConfig.findGenes()
    * deseqResultsDir=/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best500GTEx_TCGA/training/best500GTEx_TCGA.sh.out/GTEx_TCGA-design-tilda_gender_category-padj-0001-lfc-20-n-500

  + upStreamIntersectionDictionaryPath
  output from upstream stage of pipe. Used to find degree1 genes from previous upstream stages to add
      * upStreamIntersectionDictionaryPath = /private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best500GTEx_TCGA/training/best500GTEx_TCGA.sh.out/upsetPlot.out/best500.intersection.dict
  
  + historic genes was created from deseqResultsDir
  used to make sure we do not add genes that had been removed in previous stage of filtered out from some reason.
  
  + candidate genes are degree1 genes from 
  --resultsDir /private/groups/kimlab/GTEx_TCGA/1vsAll. These files have
  over 77K candidate genes
  
- **deconvolutionAnalysis/bin/1vsAll-~gender_category/best${n}CuratedDegree1**
findGenes() returns the n degree 1 genes from upstream stage

    + runs BestCuratedGeneConfig
      * selects candidate degree1 genes from intersectionDictPath =  /private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best10CuratedDegree1/training/./best10CuratedDegree1.sh.out/upsetPlot.out/best10_from_best500FindAllDegree1_wl500.intersection.dict

## best${n}CuratedDegree1 Tunning Challenges
The version of BestCuratedDegree1 from commit id ce467ff consistently does better than other version. Why?

**overview of findGenes()*  

show git log and diff for commit id
```
 $  git show 4969562
```

show git log
```
 $  git show 4969562 --no-patch
```

check out/resplace teh current version with one from commit id
```
cd extraCellularRNA/deconvolutionAnalysis/python/analysis
cp bestCuratedGeneConfig.py bestCuratedGeneConfig.py.save
git checkout ce467ff -- bestCuratedGeneConfig.py
```

```
(extraCellularRNA) aedavids@mustard $  git log -- bestCuratedGeneConfig.py
commit fbbbd8b8e86aec9ec36f2642f0bd400a247bb3fb (HEAD -> elife)
Author: Andrew E. Davidson <aedavids@ucsc.edu>
Date:   Fri Jan 12 14:41:29 2024 -0800

    return sorted baseMean descending order. seems to under perform

commit a9ce3db8292ae9b5ddaef9d5b1dbf23a4b23d954
Author: Andrew E. Davidson <aedavids@ucsc.edu>
Date:   Thu Jan 11 11:49:20 2024 -0800

    reverted bestCuratedGeneConfig.py to _ce467ff version. it has hte best performance on bestNCuratedDegree1

commit 06f2c4a1cfcf3cd8143ba9356e8cc31863f6184a
Author: Andrew E. Davidson <aedavids@ucsc.edu>
Date:   Tue Jan 9 11:22:44 2024 -0800

    python typo bug fixes.

commit 866b9c19e0a4af13d93eba353d3d1b66a12527b6
Author: Andrew E. Davidson <aedavids@ucsc.edu>
Date:   Mon Jan 8 15:06:01 2024 -0800

    add currated sign gene confi. refactored code to find dict key from file name. All unit tests run. did not add any new tests
```

- commit id ce467ff
  + Date:   Tue Jan 9 11:28:50 2024 -0800
  + if clause: we return the top n from upstream. these are biologically signifignat and sorted in decending order by baseMean. <span style="color:red">we can not gaurantee the ordering</span>. It only works if you choose the upstream stage carefull. It will be easy to misconfigure.
  + else clause: we can not ensure the degree 1 genes are ordered in any way
  ```
  if key not in degree1
      # we pass through
      # deseqResultsDir=/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best500GTEx_TCGA/training/best500GTEx_TCGA.sh.out/GTEx_TCGA-design-tilda_gender_category-padj-0001-lfc-20-n-500
      return deseqDF.head(n=n)
    
  else 
      return n degree1 genes for category sorted by baseMean
  ```
  
- commit id 866b9c
  + git checkout 866b9c -- bestCuratedGeneConfig.py
  + Date:   Mon Jan 8 15:06:01 2024 -0800
  + return first n degree1 genes sorted by baseMean
  + this crashes if key is not in degree1IntersectionDict
  + instersection dict genes created by upset plot. they are not in any partricular order
  + <span style="color:red">return sorted by baseMean in acsending order!</span>

  
- commit id 06f2c4a
  + Date:   Tue Jan 9 11:22:44 2024 -0800
  + else class same as 866b9c
  
- commit id  a9ce3db
  + Date:   Tue Jan 9 11:22:44 2024 -0800
  + base class change  from SignatureGeneConfiguration to BestSignatureGeneConfig
  + <span style="color:red">return sorted by baseMean in acsending order!</span>
  + no structral change. we did not use base class
  + add comments
  
  
- no.artificial enrichement. **we did not commit the code**
  + if category not in degree 1 intersections, return empty dataframe
  
- commit id fbbbd8
  + Date:   Fri Jan 12 14:41:29 2024 -0800
  + return sorted by baseMean decending order


# explore performance difference

what is diff between best10CuratedDegree1 and best10CuratedDegree1_ce467ff same number of genes. I think diff is sort order.

best10CuratedDegree1_ce467ff has sligthly better mean sesitivity, much better LUAD sensitity and slightly better LUSC sensitivity


## <span style="color:red" >weird ce467ffAllGenes - best10AllGenes : set()<span>

best10AllGenes - ce467ffAllGenes : set() == {}

