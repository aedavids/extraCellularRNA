# ref https://api.firecloud.org/ga4gh/v1/tools/aedavids:SamToFastq_copy/versions/1/plain-WDL/descriptor
# Andrew Davidson
# aedavids@ucsc.edu
workflow bamToFastqTest {
    call bamToFastq
}


#
# out of memory error
#
# SamToFastq is a java program
#
# We spent a lot time choosing a runtime memory parameter that allows
# bam2fastq to run most of the time. This value turned out be
#   memoryGb = 20
#   diskSpaceGb=40
#
# we found that for failed GTEx lung, prostate, panc, and colon sample
# all we had to was bump memory to 32Gb. We thought this was over kill
# but in the intesest of time it worked. No need to tune this paramter
#
# GTEx whole blood sample always fail with 32GB.  retry-with-more-memory
# did not resolve this problem
#
# possible explination:
# by default java had a max heap space of  6.97G. Based on based on 'free'
# output our docker img uses 16 gb. (double check? hard to belive it is so big)
# so with memory = 20GB java only had 4Gb heap. SamTofastq can be configure to
# dump memory to a tmp file system
#
# increase memoryGB had no effect above 23 gb.
#
# solution set java max heap size
#
# TODO
# change shell variable syntax use ${foo} for input variable. use "$bar" for shell varibles
# example
#  freeSpace=$(expr 2 '*' 1024 )
#  MAX_HEAP=`echo "$freeSpace"b`
#  ref: https://support.terra.bio/hc/en-us/community/posts/360077548112--task-command-how-to-work-with-local-variables
#

task bamToFastq {

    # If you have hard-coded runtime attributes, like memory, you can change them without breaking call caching. 
    # If you define memory as a variable, they are considered “inputs” and result in a cache miss
    # when you change it even if they aren’t used in the "command" section.    

    # input {
        File inputBam
        String sampleName
        Int memoryGb=32
        Int diskSpaceGb=40

        #
        # out of memory ref
        # https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/2019-02-11-2018-08-12/13259-SamToFastq-on-RNAseq-bams-running-out-of-memory
        #
        
        #
        # https://cloud.google.com/kubernetes-engine/docs/how-to/preemptible-vms
        # instances that last a maximum of 24 hours in general, and provide no availability guarantees.
        # Preemptible VMs are priced lower than standard Compute Engine
        #
        # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#preemptible
        # Take an Int as a value that indicates the maximum number of times Cromwell should request a
        # preemptible machine for this task before defaulting back to a non-preemptible one.
        # default value: 0    
        Int runtime_preemptible = 3

    #}

    command <<<
        # put copy of runtime parameters in output. Makes debug easier
        echo "runtime parameters"
        echo "inputBam   : ${inputBam}"
        echo "sampleName : ${sampleName}"
        echo "memoryGb   : ${memoryGb}"
        echo "diskSpaceGb: ${diskSpaceGb}"

        echo
        echo pwd: `pwd`

        echo
        echo "runtime_preemptible: ${runtime_preemptible}"
        echo
        
        # use cpuinfo to debug preemptiple/cpu/threading quotas and performance
        # cat /proc/cpuinfo; echo
        # example cpuinfo output. you get a block like this for each cpu on the machine
        # processor	: 0
        # vendor_id	: GenuineIntel
        # cpu family	: 6
        # model		: 63
        # model name	: Intel(R) Xeon(R) CPU @ 2.30GHz
        # stepping	: 0
        # microcode	: 0x1
        # cpu MHz		: 2299.998
        # cache size	: 46080 KB
        # physical id	: 0
        # siblings	: 6
        # core id		: 0
        # cpu cores	: 3
        # apicid		: 0
        # initial apicid	: 0
        # fpu		: yes
        # fpu_exception	: yes
        # cpuid level	: 13
        # wp		: yes
        # flags		: fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 ss ht syscall nx pdpe1gb rdtscp lm constant_tsc rep_good nopl xtopology nonstop_tsc cpuid tsc_known_freq pni pclmulqdq ssse3 fma cx16 pcid sse4_1 sse4_2 x2apic movbe popcnt aes xsave avx f16c rdrand hypervisor lahf_lm abm invpcid_single pti ssbd ibrs ibpb stibp fsgsbase tsc_adjust bmi1 avx2 smep bmi2 erms invpcid xsaveopt arat md_clear arch_capabilities
        # bugs		: cpu_meltdown spectre_v1 spectre_v2 spec_store_bypass l1tf mds swapgs
        # bogomips	: 4599.99
        # clflush size	: 64
        # cache_alignment	: 64
        # address sizes	: 46 bits physical, 48 bits virtual
        # power management:

        # is retry with more memory actuall restarting the task with more memory?
        # terra only stdout and stderr from last retry are saved
        free -h; echo 
        #             total       used       free     shared    buffers     cached
        # Mem:           31G        16G        14G       792K        82M        15G
        # -/+ buffers/cache:       1.1G        30G
        # Swap:           0B         0B         0B        

        # free should show swap. double check
        grep Swap /proc/meminfo; echo 
        # SwapCached:            0 kB
        # SwapTotal:             0 kB
        # SwapFree:              0 kB
        
        # file partision sizes
        df -h; echo 
        # Filesystem                         Size  Used Avail Use% Mounted on
        # overlay                             11G  3.2G  7.5G  30% /
        # tmpfs                               64M     0   64M   0% /dev
        # tmpfs                               16G     0   16G   0% /sys/fs/cgroup
        # shm                                 32G     0   32G   0% /dev/shm
        # /dev/sda1                           11G  3.2G  7.5G  30% /google
        # /dev/disk/by-id/google-local-disk   40G   13G   28G  32% /cromwell_root
        # tmpfs                               16G     0   16G   0% /proc/acpi
        # tmpfs                               16G     0   16G   0% /proc/scsi
        # tmpfs                               16G     0   16G   0% /sys/firmware        

        
        # check if java is using all avaiable memory or not
        java -version; echo 

        # https://docs.oracle.com/javase/8/docs/technotes/tools/unix/java.html
        java -XshowSettings:vm; echo 
        # check stderr
        # Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/cromwell_root/tmp.9b5dc172
        # VM settings:
        # Max. Heap Size (Estimated): 6.97G
        # Ergonomics Machine Class: server
        # Using VM: OpenJDK 64-Bit Server VM        

        #
        # find amount of free memory and unused swap
        #
        free --total --bytes | grep Total | tee free.out
        echo
        # Total:  1653340639232 435497533440 59898118144

        # looks like tab is replace with 4 space. Is the a mac terminal thing?
        # replace all white space to make sparcing portabl
        sed -e 's/\s\+/,/g' free.out | tee sed.out
        echo
        # Total:,1653340639232,435497533440,59898118144

        # parse out last numeric value
        freeBytes=`cut -d , -f 4 sed.out`
        printf "freeBytes: $freeBytes \n"
        #59898118144

        # make sure we have multiple of 1024 bytes

        quotient=$(expr $freeBytes / 1024)
        printf "quotient: $quotient \n"
        # 58494256


        # TODO:  make sure we have more than 2 mb

        MAX_HEAP=`echo "$quotient"K`
        printf "MAX_HEAP: $MAX_HEAP \n"
        
	java  \
        -Xmx"$MAX_HEAP" \
        -XshowSettings:vm \
        -jar /usr/gitc/picard.jar SamToFastq \
	TMP_DIR=. \
	INPUT=${inputBam} \
	FASTQ=${sampleName}.1.fastq.gz \
	INTERLEAVE=false \
	SECOND_END_FASTQ=${sampleName}.2.fastq.gz \
	INCLUDE_NON_PF_READS=true \
	CLIPPING_ATTRIBUTE=XT \
	CLIPPING_ACTION=2 \
	UNPAIRED_FASTQ=${sampleName}.unpaired.fastq.gz
    >>>

    output {
	File firstEndFastq = "${sampleName}.1.fastq.gz"
	File secondEndFastq = "${sampleName}.2.fastq.gz"
	File unpairedFastq = "${sampleName}.unpaired.fastq.gz"
    }

    runtime {
	docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
	memory: "${memoryGb} GB"

        #
        # best practice 2 core min: one for os one for work
        # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#cpu
        # In Google Cloud: this is interpreted as "the minimum number of cores to use."
        #
	cpu: "2"
	disks: "local-disk ${diskSpaceGb} SSD"
        preemptible: '${runtime_preemptible}'

        #
        # enable Out of Memory Retry
        # https://support.terra.bio/hc/en-us/articles/4403215299355
        # https://cromwell.readthedocs.io/en/develop/RuntimeAttributes/#maxretries
        # https://cromwell.readthedocs.io/en/develop/cromwell_features/RetryWithMoreMemory/#retry-with-more-memory
        #
        # 8/23/2021
        # retry with more memory does not work
        # all that happens is you burn a lot of compute time
        # maxRetries: 3
    }
}
