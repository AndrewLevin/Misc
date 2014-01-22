seed=$1

echo "pwd"
pwd

echo ""

echo \$seed
echo $seed

if hostname | grep 'lxplus.*\.cern\.ch' >& /dev/null; then
    echo "running on lxplus"
    cp -r /afs/cern.ch/work/a/anlevin/madgraph_beta/2.0.0beta4/ .
elif hostname | grep 'mit\.edu' &> /dev/null; then
    echo "running at MIT"
    cp -r /scratch/anlevin/madgraph_beta/2.0.0beta4/ .
else
    echo "running on unknown machine, exiting"
    exit
fi

echo "ls"
ls

echo ""

echo "ls 2.0.0beta4"
ls 2.0.0beta4

echo "" 

cd 2.0.0beta4

echo "pwd"
pwd

echo ""

if hostname | grep 'lxplus.*\.cern\.ch' >& /dev/null; then
    echo "running on lxplus"
    cat > commands.mg5 <<EOF
    import model SM_LS0_LS1_UFO
    define l+ = e+ mu+ ta+
    define l- = e- mu- ta-
    generate p p > w+ w+ p p QED=4 QCD=99, w+ > l+ vl
    add process p p > w- w- p p QED=4 QCD=99, w- > l- vl~
    output qed_4_qcd_99_ls0ls1_grid
    launch qed_4_qcd_99_ls0ls1_grid
    reweight=ON
    /afs/cern.ch/work/a/anlevin/UserCode/madgraph_generation/run_card.dat
    /afs/cern.ch/work/a/anlevin/UserCode/madgraph_generation/param_card.dat
    /afs/cern.ch/work/a/anlevin/UserCode/madgraph_generation/reweight_card.dat
    set clusinfo F
    set nevents 100
    set iseed $seed
    EOF
elif hostname | grep 'mit\.edu' &> /dev/null; then
    echo "running at MIT"
    cat > commands.mg5 <<EOF
    import model SM_LS0_LS1_UFO
    define l+ = e+ mu+ ta+
    define l- = e- mu- ta-
    generate p p > w+ w+ p p QED=4 QCD=99, w+ > l+ vl
    add process p p > w- w- p p QED=4 QCD=99, w- > l- vl~
    output qed_4_qcd_99_ls0ls1_grid
    launch qed_4_qcd_99_ls0ls1_grid
    reweight=ON
    /scratch/anlevin/UserCode/madgraph_generation/run_card.dat
    /scratch/anlevin/UserCode/madgraph_generation/param_card.dat
    /scratch/anlevin/UserCode/madgraph_generation/reweight_card.dat
    set clusinfo F
    set nevents 100
    set iseed $seed
    EOF
else
    echo "running on unknown machine, exiting"
    exit
fi

echo "cat commands.mg5"
cat commands.mg5
    
python2.6 ./bin/mg5_aMC commands.mg5

gunzip qed_4_qcd_99_ls0ls1_grid/Events/run_01/unweighted_events.lhe.gz

if hostname | grep 'lxplus.*\.cern\.ch' >& /dev/null; then
    echo "running on lxplus"
    cp qed_4_qcd_99_ls0ls1_grid/Events/run_01/unweighted_events.lhe /afs/cern.ch/work/a/anlevin/tmp/unweighted_events_${seed}.lhe
elif hostname | grep 'mit\.edu' &> /dev/null; then
    echo "running at MIT"
    cp qed_4_qcd_99_ls0ls1_grid/Events/run_01/unweighted_events.lhe /scratch/anlevin/tmp/unweighted_events_${seed}.lhe
fi
    
cd ../

echo "pwd"
pwd

echo ""

echo "ls"
ls

echo ""

echo "removing the madgraph directory"

echo ""

rm -r 2.0.0beta4/

echo "ls"
ls
