nevents=100
njobs=1000
gridpack_filename=/afs/cern.ch/work/a/anlevin/UserCode/madgraph_generation/wzgamma_ewk_plus_qcd_sm_gridpack.tar.gz
output_dir=/afs/cern.ch/work/a/anlevin/data/lhe/unmerged_wzgamma_qed_5_qcd_99_sm_official/

seed=1

while((seed<=njobs)); do
    bsub -q 8nh "bash run_gridpack.sh $gridpack_filename $nevents $seed $output_dir"
    seed=$(($seed+1))
done
