gen_command1=$1
gen_command2=$2
seed=$3
output_dir=$4
reweight_file=$5
model=$6


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
reweight=ON 
/scratch/anlevin/UserCode/madgraph_generation/run_card.dat
/scratch/anlevin/UserCode/madgraph_generation/param_card.dat
$reweight_file
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
reweight=ON
/scratch/anlevin/UserCode/madgraph_generation/run_card.dat
/scratch/anlevin/UserCode/madgraph_generation/param_card.dat
$reweight_file
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
