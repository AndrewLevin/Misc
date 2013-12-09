input_dir=/store/user/anlevin/data/AOD/ww_to_ll_same_sign_sm_qcd_99_qed_4/
output_dir=/eos/cms/store/user/anlevin/data/BAMBU/ww_to_ll_same_sign_sm_qcd_99_qed_4/

echo \$input_dir
echo $input_dir

echo \$output_dir
echo $output_dir

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

sleep 15

for file in `/afs/cern.ch/project/eos/installation/0.2.31/bin/eos.select ls $input_dir`
  do
  bsub -q 1nd "bash make_bambu_from_aod.sh $input_dir/$file $output_dir"
done
