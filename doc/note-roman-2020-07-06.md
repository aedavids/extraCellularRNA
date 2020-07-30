Copied from Slack

Notes from Slack direct msg on 7/6/2020

? looking for exosome data? I though there would data about counts and sized of exosomes. We only have RNA from exosome. The data is really extra cellular RNA it includes what ever was in the dish and the exosome. Acording to Roman  We do not have data where exosomes where collected, purified, then sequenced. I would double check that. Daniel suggested we have 'exosome' data. 

## H358 cell line exosome data
/public/groups/kimlab/exoTIC-biomarkers/data/exotic.v.rneasy.fastq/*.rneasy.*

## AALE cell line exosome data (this is good and under-analyzed) ****
/public/groups/kimlab/aale.luad.exo/aale.exo.data/aale.*

AALE are lung cells

reads are fastQ

star alginment and salmon quantification

Salmon contains information about  transcript counts and lengths

start contains information about what type of transcripts. intron/exon, ...

## Same as ^, but patient data
/public/groups/kimlab/aale.luad.exo/patient.data

## various dirs containing patient exosome and cell line exosome data
/public/groups/kimlab/exoRNA-biomarkers-panc/data/


## Is there evidence of bacteria or viruse DNA in hte patient samples.
- go to ucsc genome browser find data for Gammaproteobacteria. down load gtf or bedfile output

tools --> STAR (custom genome(S)), Kraken2 (I’ve played around with), SILVA database, NCBI genomes

- STAR alignments --> genomic position, distribution, overlap with features (motifs, binding sites, etc)

SALMON quantification --> transcript abundance, variability( important for detection purposes), length, isoform preference. *SOMETHING ELSE I haven’t tried: ‘Differential Isoform Usage’

aale.luad.exo: 3 AALE cell line controls, 3 AALE cell line kras, 3 Lung Cancer Patient samples w/ KRAS mutations ( I have metadata let me know), 3 Healthy Patient


/public/groups/kimlab/aale.luad.exo/aale.exo.data/patient.data

/public/groups/kimlab/aale.luad.exo/patient.data
