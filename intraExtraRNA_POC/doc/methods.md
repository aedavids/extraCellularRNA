# Intracellular, Extracellular RNA proof of concept methods
Andrew E. Davidson
aedavids@ucsc.edu

for a given type of cancer use the GTEx_TCGA 1vsAll "best" list of genes to train several binary classifiers. Can we demonstrate that we can use intracellular biomarkers to classify extracellal samples

**Types of classifiers**
<span style="color:red">do we have enough data?</span>

1. logistic regression with regulation (linear model)
2. random forest (non linear model)
3. ~~fully connected neural network with regulation~~
   + if we have time might be fun. 
   + we do not have enough data. Might be able to simulate using insilco mixtures

**list of experiments**

a. train logistic regression on intracellular, predict intracellular

b. train logistic regression on intracellular, predict extracellular
    - select a single cancer type
    - select the "best" intracellular gene ids
    - select on the best gene ids from the GTEx_TCGA training set
    - make predictions for control and cancer type extracellular sample
    
c. train logistic regression on extracellular predict extracellular
    - select a single cancer type
    - select the "best" intracellular gene ids
    - selec select "best" gene id from extracellular training set
    - make predictions for control and cancer type extracellular sample
    
d. train logistic regression on extracellular predict intracellular

e. repeat a, b, c, d, using random forest

f. repeat Roman's TE vs naive TE experimetns
