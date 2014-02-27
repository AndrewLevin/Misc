nevents=$1
seed=$2

echo \$nevents
echo $nevents

echo \$seed
echo $seed

echo "copying gridpack"
echo ""

cp /afs/cern.ch/work/a/anlevin/UserCode/madgraph_generation/run_01_gridpack.tar.gz .

echo""
echo "untaring gridpack"
echo ""

tar -xvf run_01_gridpack.tar.gz

echo ""
echo "compiling"
echo ""

cd madevent
./bin/compile
./bin/clean4grid

cd ../

echo ""
echo "begin generating events" 
echo ""

./run.sh $nevents $seed

echo ""
echo "finished generating events"
echo ""

gunzip events.lhe.gz

echo ""
echo "begin copying output file"
echo ""

cp events.lhe /afs/cern.ch/work/a/anlevin/data/lhe/unmerged_wz/events_${seed}.lhe

echo ""
echo "finished copying output file"
echo ""
