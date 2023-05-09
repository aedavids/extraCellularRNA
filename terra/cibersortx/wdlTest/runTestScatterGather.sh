
set -x

clear;

export WDL_TOOLS=/private/home/aedavids/extraCellularRNA/java/bin

rm -rf cromwell-* import.zip;

java -jar /private/home/aedavids/extraCellularRNA/java/bin/womtool-85.jar validate testScatterGather.wdl;

zip import.zip aggregateTask.wdl  createTestDataTask.wdl  partitionDataTask.wdl mergeTask.wdl;

/private/home/aedavids/extraCellularRNA/bin/runCromwell.sh \
     -Dconfig.file=/private/home/aedavids/extraCellularRNA/terra/wdl/cromwellDebug.conf \
     -jar ${WDL_TOOLS}/cromwell-85.jar \
     run \
     --inputs testScatterGather.wdl.input.json \
     --imports /private/home/aedavids/extraCellularRNA/terra/cibersortx/wdlTest/import.zip \
      testScatterGather.wdl


printf "\n\n\n\n*************\n"
find cromwell-executions
