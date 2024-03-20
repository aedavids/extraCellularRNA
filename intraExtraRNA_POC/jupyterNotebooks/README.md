# Table of Contents

- baseSignatureGenesBaseMeanBugFixPOC.ipynb
  test filtering genes using lfc, adj p-value, and baseMean

- testIntraCellularBiomarkersOnExtracellularSamples.ipynb
  Preliminary evidence that intracellular biomarkers can be used on extracellular exomsomal data sets. I trained a logisitc regression classifier on our 32 PANC sample. Select best 20 gene loci. Accuracy: 0.97

- testIntraCellular-TE-BiomarkersOnExtracellularSamples.ipynb
  Preliminary evidence that intracellular biomarkers can be used on extracellular exomsomal data sets. I trained a logisitc regression classifier on our 32 PANC sample. Select best 20 gene loci. Accuracy: 0.81. <span style="color:red">Accuracy is deceiving. Notice false positives > true positives for TE's</span>

- testIntraCellularBiomarkers
  preliminary results
  trained logistic regresion on 309 GTEx and TCGA pacresae and PAAD samples. accuracy 0.984

- testIntraCellular-TE-Biomarkers.ipynb
  preliminary results. trained a binary classifier using "best" 1vsAll PAAD genes. 0.96 accuracty on the training set. I select all the PAAD and panc sample. Todo make predictions on test data set
  
- testIntraCellularModelOnExtracellularSamples
  preliminary results. models trained on intracellular data make poor predictions on extracellular data. Is the evidence of selective packaing? that is say the transcript distributions are different. Accuracy on 32 plasma samples was 0.312. All samples where predicted to be PAAD

- testIntraCellular-TE-ModelOnExtracellularSamples
  never completed code (ran out of time). intracellar model trained using best features did not work. Unlikely TE features will 

  