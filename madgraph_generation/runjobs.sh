#!/bin/bash

i=0
imax=300

while((i<=imax)); do
    cat > submit.cmd <<EOF
universe = vanilla
Executable = /scratch/anlevin/UserCode/madgraph_generation/make_lhe_weights.sh
Arguments = $i
GetEnv = True
Requirements = (Arch == "X86_64") && (OpSys == "LINUX") && (Disk >= DiskUsage) && ((Memory * 1024) >= ImageSize) && (HasFileTransfer)
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT 
Output = stderr_stdout_${i}.dat
Error =stderr_stdout_${i}.dat
Log = log_${i}.dat
transfer_input_files = /tmp/anlevin/job_files.tar
Queue 1
EOF
    

i=$(($i+1))
condor_submit submit.cmd
rm submit.cmd
  
done
