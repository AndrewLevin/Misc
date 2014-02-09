gen_command1=$1
gen_command2=$2
seed=$3
output_dir=$4
reweight_file=$5
model=$6
run_card=$7
param_card=$8

echo "pwd"
pwd

echo ""

echo \$seed
echo $seed

echo \$output_dir
echo $output_dir

echo \$reweight_file
echo $reweight_file

echo \$gen_command1
echo $gen_command1

echo \$gen_command2
echo $gen_command2

echo \$model
echo $model

echo \$run_card
echo $run_card

echo \$param_card
echo $param_card

echo "begin checking if there is a reweight file"
if [ $reweight_file != "NONE" ]; then
    echo "no reweight file"
fi
echo "finished checking if there is a reweight file"

n_events_per_job=1000

if hostname | grep 'lxplus.*\.cern\.ch' >& /dev/null; then
    echo "running on lxplus"
    cp -r /afs/cern.ch/work/a/anlevin/madgraph_beta/2.0.0beta4/ .
elif hostname | grep 'mit\.edu' &> /dev/null; then
    echo "running at MIT"
    cp -r /scratch/anlevin/MG5_aMC_v2_0_0/ .
else
    echo "running on unknown machine, exiting"
    exit
fi

echo "ls"
ls

echo ""

echo "ls MG5_aMC_v2_0_0"
ls MG5_aMC_v2_0_0

echo "" 

cd MG5_aMC_v2_0_0

echo "pwd"
pwd

echo ""

if hostname | grep 'lxplus.*\.cern\.ch' >& /dev/null; then
    echo "running on lxplus"
    cat > commands.mg5 <<EOF
set nb_core 2
import $model
define l+ = e+ mu+ ta+
define l- = e- mu- ta-
$gen_command1
$gen_command2
output mg_dir
launch mg_dir
EOF
if [ $reweight_file != "NONE" ]; then
    echo "reweight=ON" >> commands.mg5 
fi
echo $run_card >> commands.mg5
echo $param_card >> commands.mg5
if [ $reweight_file != "NONE"]; then
    echo $reweight_file >> commands.mg5
fi
cat >> commands.mg5 <<EOF
set clusinfo F
set nevents $n_events_per_job
set iseed $seed
EOF
elif hostname | grep 'mit\.edu' &> /dev/null; then
    echo "running at MIT"
    cat > commands.mg5 <<EOF
set nb_core 2
import $model
define l+ = e+ mu+ ta+
define l- = e- mu- ta-
$gen_command1
$gen_command2
output mg_dir
launch mg_dir
EOF
if [ $reweight_file != "NONE" ]; then
     echo "reweight=ON" >> commands.mg5
fi
echo $run_card >>  commands.mg5
echo $param_card >> commands.mg5
if [ $reweight_file != "NONE" ]; then
    echo $reweight_file >> commands.mg5
fi
cat >> commands.mg5 <<EOF
set clusinfo F
set nevents $n_events_per_job
set iseed $seed
EOF
else
    echo "running on unknown machine, exiting"
    exit
fi

echo "cat commands.mg5"
cat commands.mg5
    
python2.6 ./bin/mg5_aMC commands.mg5

gunzip mg_dir/Events/run_01/unweighted_events.lhe.gz

cp mg_dir/Events/run_01/unweighted_events.lhe ${output_dir}unweighted_events_${seed}.lhe

cd ../

echo "pwd"
pwd

echo ""

echo "ls"
ls

echo ""

echo "removing the madgraph directory"

echo ""

rm -r MG5_aMC_v2_0_0/

echo "ls"
ls
