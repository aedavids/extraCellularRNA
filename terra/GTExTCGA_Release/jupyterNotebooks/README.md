# file summary
Andrew E. Davidson  
aedavids@ucsc.edu  
2/22/24

- boxPlots.ipynb
create a series of box plots. see https://drive.google.com/drive/folders/1MEdIZRbMbMTpRaFY2iad3tS2rGwsT7m3?usp=drive_link


- fixNonDisjointSets.ipynb
looks like when we create the origin train/validate/test data set we lost some samples and some samples appear in more than one of the train/validate/test data sets. Creates usable data sets.

- lumpTrainingSets.ipynb
merge the training, validation, and test groupby counts data files. This will allow us to calculate normalized counts for our data set categorization plots.

- checkColDataForDuplicateSamples.ipynb
looks like when we create the origin train/validate/test data set we lost some samples and some samples appear in more than one of the train/validate/test data sets. This notebook explores what the bug looks like see fixNonDisjointSets.ipynb for resolution

-  normalizeCounts.ipynb 
calculates DESeq2 Normalization on the groupby counts

- calculateScalingFactor.ipynb
uses DESeq2's method for accounting for library size and composition 

- basicSummaryStats.ipynb
this is a format test. count number of male/female, number of samples in each category, ...

- barCharts.ipynb
this is format test. create plots that show proportions + error bar. E.G. for each category plot number of TEs / number of not TEs











