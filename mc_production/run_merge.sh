input_dir=$1
output_dir=$2
output_name=$3

if ! [ $input_dir ] || ! [ $output_dir ] || ! [ $output_name ]; then
    echo "usage: bash run_merge.sh input_dir output_dir output_name"
    echo "example: bash run_merge.sh /eos/cms/store/user/anlevin/data/BAMBU/qed_4_qcd_99_sm_unmerged/ /eos/cms/store/user/anlevin/data/BAMBU/qed_4_qcd_99_sm/ qed_4_qcd_99_sm"
    exit
fi

echo \$input_dir
echo $input_dir

echo \$output_dir
echo $output_dir

echo \$output_name
echo $output_name

counter=1
label=1
for file in `/afs/cern.ch/project/eos/installation/0.2.31/bin/eos.select ls $input_dir`
  do

  file_list="\"root://eoscms.cern.ch/${input_dir}${file}\",${file_list}"
  if (( counter == 10 ))
      then
      root  -l -b -q merge.C+\(${file_list}\"${output_name}_${label}\"\)
      counter=0
      /afs/cern.ch/project/eos/installation/0.2.31/bin/eos.select cp ${output_name}_${label}_000.root $output_dir
      rm ${output_name}_${label}_000.root
      label=$(($label+1))
      file_list=""
  fi
  counter=$(($counter+1))
done
