# How to resolve v35 to v39 HUGO to ENSG mapping errors

we can ignore the ENSG decimal point. It encodes version number

/private/groups/kimlab/genomes.annotations/genomes.annotations/gencode.39/gencode.v39.annotation.expanded.tx.to.gene.tsv

txt2GeneFilePath =  "/private/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv"

```
                gencode.v35.ucsc.rmsk.tx.to.gene.csv          v39 tx2g
	    refSeq	    ENSG	              bioType
97712	AC111149.2	ENSG00000253339.2	  lncRNA        ENSG00000253339.3
172015	AC092140.2	ENSG00000274031.1	  lncRNA        
203144	PCAT19	        ENSG00000267107.8 lncRNA        ENSG00000267107.9
```

