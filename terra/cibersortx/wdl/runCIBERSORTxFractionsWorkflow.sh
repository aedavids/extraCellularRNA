set -x

clear

export WDL_TOOLS=/private/home/aedavids/extraCellularRNA/java/bin

rm -rf cromwell-* import.zip;

# -j just store file name ignore directory path
zip -j import.zip cibersortxFractionsTask.wdl ../wdlTest/partitionDataTask.wdl ../wdlTest/mergeTask.wdl ;


java -jar /private/home/aedavids/extraCellularRNA/java/bin/womtool-85.jar validate CIBERSORTxFractionsWorkflow.wdl ;

/private/home/aedavids/extraCellularRNA/bin/runCromwell.sh \
     -Dconfig.file=/private/home/aedavids/extraCellularRNA/terra/wdl/cromwellDebug.conf \
     -jar ${WDL_TOOLS}/cromwell-85.jar \
     run \
     --imports /private/home/aedavids/extraCellularRNA/terra/cibersortx/wdl/import.zip \
     --inputs CIBERSORTxFractionsWorkflow.wdl.input.json \
     CIBERSORTxFractionsWorkflow.wdl

printf "\n\n\n\n*************\n"
find cromwell-executions
