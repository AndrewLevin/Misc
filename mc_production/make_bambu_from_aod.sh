
inputfile=$1
output_dir=$2

echo \$output_dir
echo $output_dir

echo \$inputfile
echo $inputfile

echo pwd
pwd

cd /afs/cern.ch/work/a/anlevin/UserCode/mc_production/CMSSW_5_3_13_patch3/
#cd /afs/cern.ch/work/a/anlevin/UserCode/mc_production/CMSSW_5_3_11/
eval `scramv1 runtime -sh`

cd -

echo pwd
pwd

cp /afs/cern.ch/work/a/anlevin/UserCode/mc_production/BAMBUProd_AODSIM.py delete_this.py


tmpfile=`mktemp`
echo $inputfile | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" >& $tmpfile
inputfile_escapes=`cat $tmpfile`

echo \$inputfile_escapes
echo $inputfile_escapes

cat delete_this.py | sed "s/change_this/$inputfile_escapes/" >& BAMBUProd_AODSIM.py

cmsRun BAMBUProd_AODSIM.py

outputfile=`echo $inputfile | sed "s/\//___/" | sed "s/\//___/" | sed "s/\//___/" | sed "s/\//___/" | sed "s/\//___/" | sed "s/\//___/" | sed "s/\\//___/" | sed "s/\//___/" | sed "s/\//___/" | sed "s/\//___/" | sed "s/\//___/" | sed "s/\//___/" | sed "s/\//___/" | sed "s/\//___/" | sed "s/\//___/" | sed "s/\//___/"` 

echo \$outputfile
echo $outputfile

/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select cp XX-MITDATASET-XX_000.root ${output_dir}$outputfile
