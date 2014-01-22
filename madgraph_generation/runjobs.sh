#!/bin/bash

i=0
imax=200

output_dir="/scratch/anlevin/data/lhe/qed_2_qcd_99_lm0123_grid/"

gen_command1="'generate p p > w+ w+ p p QED=2 QCD=99, w+ > l+ vl'"
gen_command2="'add process p p > w- w- p p QED=2 QCD=99, w- > l- vl~'"

model="SM_LM0123_UFO"

reweight_file="/scratch/anlevin/UserCode/madgraph_generation/reweight_card_lm0lm1lm2lm3.dat"

if hostname | grep 'lxplus.*\.cern\.ch' >& /dev/null; then 
    echo "running on lxplus"
    while((i<=imax)); do
	bsub -q 1nd "bash /afs/cern.ch/work/a/anlevin/UserCode/madgraph_generation/make_lhe_weights.sh $i"
	i=$(($i+1))
    done
elif hostname | grep 'mit\.edu' &> /dev/null; then
    echo "running at MIT"
    while((i<=imax)); do
	cat > submit.cmd <<EOF
universe = vanilla
Executable = /scratch/anlevin/UserCode/madgraph_generation/make_lhe_weights.sh
Arguments = "$gen_command1 $gen_command2 $i $output_dir $reweight_file $model"
GetEnv = True
Requirements = (Arch == "X86_64") && (OpSys == "LINUX") && (Disk >= DiskUsage) && ((Memory * 1024) >= ImageSize) && (HasFileTransfer)
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT 
Output = stderr_stdout_lm0123_${i}.dat
Error =stderr_stdout_lm0123_${i}.dat
Log = log_lm0123_${i}.dat
Queue 1
EOF
	i=$(($i+1))
	condor_submit submit.cmd
	rm submit.cmd
    done
else
    echo "running on unknown machine, exiting"
    exit
fi
