input_dir=/store/user/anlevin/data/AOD/qed_4_qcd_99_lm0123/
output_dir=/eos/cms/store/user/anlevin/data/BAMBU/qed_4_qcd_99_lm0123_unmerged/

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
  bsub -q 8nh "bash make_bambu_from_aod_2.sh $input_dir/$file $output_dir"
done
