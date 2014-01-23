input_file=/afs/cern.ch/work/a/anlevin/data/lhe/qed_2_qcd_99_lm0123_grid.lhe
output_dir=/eos/cms/store/user/anlevin/data/AOD/qed_2_qcd_99_lm0123/
echo \$input_file
echo $input_file

echo \$output_dir
echo $output_dir

nevents=`cat $input_file  | grep "<event>" | wc -l`

echo \$nevents
echo $nevents

if((nevents==0));
    then
    echo "no events in lhe file, exiting";
    exit;
fi

if ! /afs/cern.ch/project/eos/installation/0.2.31/bin/eos.select ls $output_dir >& /dev/null
    then
    echo "output directory does not exist, exiting"
    exit
fi

outputdirsize=`/afs/cern.ch/project/eos/installation/0.2.31/bin/eos.select ls $output_dir | wc -l`

if((outputdirsize!=0))
    then
    echo "output directory is not empty, exiting"
    exit
fi


max_events=20000
where_to_start=0

if((nevents<max_events))
    then
    echo "nevents < maxevents, exiting"
    exit
fi

sleep 15

while((where_to_start<=max_events)); do
    echo $where_to_start
    bsub -q 1nd "bash make_aod_from_lhe_2.sh $input_file 300 100 $where_to_start $where_to_start $output_dir"
    where_to_start=$(($where_to_start+100))
done
