gridpack_filename=$1
nevents=$2
seed=$3
output_dir=$4

echo \$gridpack_filename
echo $gridpack_filename

echo \$nevents
echo $nevents

echo \$seed
echo $seed

echo \$output_dir
echo $output_dir

echo "copying gridpack"
echo ""

cp $gridpack_filename run_01_gridpack.tar.gz

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

cp events.lhe ${output_dir}events_${seed}.lhe

echo ""
echo "finished copying output file"
echo ""
