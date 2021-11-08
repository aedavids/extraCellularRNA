
# java -jar ~/extraCellularRNA/java/bin/womtool-58.jar validate testBashTask.wdl                             
# java -jar ~/extraCellularRNA/java/bin/womtool-58.jar inputs testBashTask.wdl > testBashTask.wdl.input.json
# java -jar ~/extraCellularRNA/java/bin/cromwell-58.jar run --inputs testBashTask.wdl.input.json testBashTask.wdl

# where to find output
# cromwell-executions/testBashTask/61d8ad79-f976-42b4-8a0b-cf2f12594a0f/call-bashTask/execution/stdout
# general form cromwell-execution/workflow name/GUID/call-taskName/execution/
# useful files
#   stdout
#   stderr
#   script
#      cromwell generates a bash script with replacing task input variables with there values
#      useful for debugging
#
#
# the --input file can be found at
# ls cromwell-executions/testBashTask/e23b4923-cc77-4878-ab0c-0f148d26bcf3/call-bashTask/inputs/-1375731305/testBashTask.wdl.input.json 


workflow testBashTask {
    String name
    File fileName
    #File dirName

  call bashTask {
    input :
          name = name,
          fileName = fileName,
          #dirName = dirName
  }
}

task bashTask {
    String name
    File fileName
    #File dirName
    
  command {
    # see bash man page "SHELL BUILTIN COMMANDS" for details
    # set -euxo pipefail
    # set -e Exit immediately if a pipeline see shell builtin command it is more complicated
    # aedwip_set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout; merge these
    # set -o print value of current options

    
    # we do not control the execution of this script
    # cromwell/docker ? puts set -x in stderr. this makes debugging harder
    # https://ops.tips/gists/redirect-all-outputs-of-a-bash-script-to-a-file/#using-exec-to-redirect-stdout-and-stderr
    # Redirect standard error to standard out such that 
    # standard error ends up going to wherever standard
    # out goes (the file).
    exec 2>&1

    set -x
    echo 'argument was specified in task input: hello ${name}!'

    pwd


    ls -l ${fileName}

    cat ${fileName}

    # crash and burn

    # we can not expand varibles if they are not in the input list

    ls .

    # how to local variables work?
    # syntax for local variable is different that syntx for task input variables
    count=`ls . | wc -l`
    echo "the number of files in . is $count"

  }
  output {
    File response = stdout()
  }
  runtime {
   docker: 'ubuntu:latest'
  }
}
