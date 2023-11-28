set -x

clear

export WDL_TOOLS=/private/home/aedavids/extraCellularRNA/java/bin

rm -rf cromwell-* import.zip;

#zip import.zip cibersortxFractionsTask.wdl ../wdlTest/partitionDataTask.wdl;


java -jar /private/home/aedavids/extraCellularRNA/java/bin/womtool-85.jar validate cibersortxFractionsTask.wdl;

/private/home/aedavids/extraCellularRNA/bin/runCromwell.sh \
     -Dconfig.file=/private/home/aedavids/extraCellularRNA/terra/wdl/cromwellDebug.conf \
     -jar ${WDL_TOOLS}/cromwell-85.jar \
     run \
     --inputs cibersortxFractionsTask.wdl.input.json \
     cibersortxFractionsTask.wdl


printf "\n\n\n\n*************\n"
find cromwell-executions -type f
