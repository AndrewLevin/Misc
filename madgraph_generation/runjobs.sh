#!/bin/bash


i=0
imax=1

#output_dir="/scratch/anlevin/data/lhe/qed_4_qcd_99_ls0ls1_grid/"
output_dir="/afs/cern.ch/work/a/anlevin/data/lhe/qed_4_qcd_99_unmerged/"

gen_command1="'generate p p > w+ w+ p p QED=4 QCD=99, w+ > l+ vl'"
gen_command2="'add process p p > w- w- p p QED=4 QCD=99, w- > l- vl~'"

model="SM_LS01_LT012_LM0167"
#model="sm"

#reweight_file="/scratch/anlevin/UserCode/madgraph_generation/reweight_card_lt0lt1lt2.dat"
reweight_file="NONE"

run_card="/scratch/anlevin/UserCode/madgraph_generation/run_card_no_matching.dat"
param_card="/scratch/anlevin/UserCode/madgraph_generation/param_card_lt012.dat"

echo "output directory: "$output_dir
echo "using reweight_card: "$reweight_file
echo "using run_card: "$run_card
echo "using param_card: "$param_card

if ! [ -e $output_dir ]
    then
    echo "output directory does not exist"
    exit
fi

if ((`ls $output_dir | wc -l`))
    then
    echo "output directory is not empty, exiting"
    exit
fi

if [ $reweight_file != "NONE" ] && ! ls $reweight_file >& /dev/null
    then
    echo "reweight file does not exist, exiting"
    exit
fi

if ! ls $run_card >& /dev/null
    then
    echo "runcard does not exist, exiting"
    exit
fi

if ! ls $param_card >& /dev/null
    then
    echo "paramcard does not exist, exiting"
    exit
fi

if ! ls $output_dir >& /dev/null
    then
    echo "output directory does not exist, exiting"
    exit
fi

sleep 15

if hostname | grep 'lxplus.*\.cern\.ch' >& /dev/null; then 
    echo "running on lxplus"
    while((i<=imax)); do
	bsub -q 1nd "bash /afs/cern.ch/work/a/anlevin/UserCode/madgraph_generation/make_lhe_weights.sh $gen_command1 $gen_command2 $i $output_dir $reweight_file $model $run_card $param_card"
	i=$(($i+1))
    done
elif hostname | grep 'mit\.edu' &> /dev/null; then
    echo "running at MIT"
    while((i<=imax)); do
	cat > submit.cmd <<EOF
universe = vanilla
Executable = /scratch/anlevin/UserCode/madgraph_generation/make_lhe_weights.sh
Arguments = "$gen_command1 $gen_command2 $i $output_dir $reweight_file $model $run_card $param_card"
GetEnv = True
Requirements = (Arch == "X86_64") && (OpSys == "LINUX") && (Disk >= DiskUsage) && (HasFileTransfer) && (Machine != "t3btch084.mit.edu") && (Machine != "t3btch085.mit.edu") && (Machine != "t3btch086.mit.edu") && (Machine != "t3btch087.mit.edu") && (Machine != "t3btch088.mit.edu")
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT 
Output = stderr_stdout_ls0ls1_${i}.dat
Error = stderr_stdout_ls0ls1_${i}.dat
Log = log_ls0ls1_${i}.dat
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
