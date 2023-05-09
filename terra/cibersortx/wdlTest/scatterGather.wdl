# https://github.com/openwdl/wdl/blob/main/versions/draft-2/SPEC.md#scatter--gather
# java -jar /private/home/aedavids/extraCellularRNA/java/bin/womtool-85.jar validate scatterGather.wdl;
# rm -rf cromwell*;
# /private/home/aedavids/extraCellularRNA/bin/runCromwell.sh -Dconfig.file=/private/home/aedavids/extraCellularRNA/terra/wdl/cromwellDebug.conf      -jar ${WDL_TOOLS}/cromwell-85.jar      run  scatterGather.wdl 

task inc {
  Int i

  command <<<
  python -c "print(${i} + 1)"
  >>>

  output {
    Int incremented = read_int(stdout())
  }
}

task sum {
  Array[Int] ints

  command <<<
  python -c "print(${sep="+" ints})"
  >>>

  output {
    Int sum = read_int(stdout())
  }
}

workflow wf {
  Array[Int] integers = [1,2,3,4,5]
  scatter (i in integers) {
    call inc {input: i=i}
  }
  call sum {input: ints = inc.incremented}
}

