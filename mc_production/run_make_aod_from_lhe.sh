nevents=`cat /afs/cern.ch/work/a/anlevin/data/lhe/ww_to_ll_same_sign_anom_8_tev_0_tev-4.lhe  | grep "<event>" | wc -l`
echo $nevents
max_events=100000
where_to_start=0
while((where_to_start<=max_events)); do
    echo $where_to_start
    bsub -q 1nd "bash make_aod_from_lhe.sh 300 100 $where_to_start $where_to_start /eos/cms/store/user/anlevin/data/AOD/unmerged/"
    where_to_start=$(($where_to_start+100))
done
