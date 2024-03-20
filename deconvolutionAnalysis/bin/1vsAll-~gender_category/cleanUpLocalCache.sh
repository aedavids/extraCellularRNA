
cleanUp(){
    dirsToRemove=$@
    for i in $dirsToRemove;
    do

        printf "\n#####\n"
        du -sh $i

        printf "deleted ${i}?\n"
        select yn in "Yes" "No"; do
            case $yn in
                Yes ) 'rm' -rf "${i}"; break;;
                No )  break;;
            esac
            
        done # end select
    done # end for
}
    
root=/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category

localCacheDirs=`find $root -type d -name private`

printf "\n*********** removing locally cached data files\n"
cleanUp $localCacheDirs

cromwellCacheDirs=`find $root -type d -name cromwell-executions`
cromwellLogs=`find $root -type d -name cromwell-workflow-logs`
printf "\n*********** removing cromwell log files and local cache\n"
cleanUp $cromwellCacheDirs
cleanUp $cromwellLogs


