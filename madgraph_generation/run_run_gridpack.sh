nevents=25
njobs=1000

seed=1

while((seed<=njobs)); do
    bsub -q 8nh "bash run_gridpack.sh $nevents $seed"
    seed=$(($seed+1))
done
