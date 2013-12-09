export SCRAM_ARCH=slc5_amd64_gcc462

scram p CMSSW CMSSW_5_3_11_patch2

cd CMSSW_5_3_11_patch2/src/

input_file=$1
n_events_pileup=$2
n_events=$3
where_to_start=$4
lumi_number=$5
output_dir=$6

echo \$input_file
echo $input_file

echo \$n_events_pileup
echo $n_events_pileup

echo \$n_events
echo $n_events

echo \$output_dir
echo $output_dir

echo \$where_to_start
echo $where_to_start

echo \$lumi_number
echo $lumi_number

echo pwd
pwd

echo date
date

eval `scramv1 runtime -sh`

mkdir -p Configuration/GenProduction/python/EightTeV/

cp /afs/cern.ch/work/a/anlevin/UserCode/mc_production/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_tauola_cff.py Configuration/GenProduction/python/EightTeV/hadronizer.py
cp /afs/cern.ch/work/a/anlevin/UserCode/mc_production/randomizeSeeds.py Configuration/GenProduction/python/

scram b

echo "begin making pileup file"

cmsDriver.py MinBias_8TeV_cfi  --conditions auto:startup -s GEN,SIM --datatier GEN-SIM -n $n_events_pileup --eventcontent RAWSIM --fileout pileup.root --customise Configuration/GenProduction/randomizeSeeds.randomizeSeeds

echo "finished making pileup file"

echo "ls -lh"
ls -lh

echo "begin making GEN-SIM file"

cmsDriver.py Configuration/GenProduction/python/EightTeV/hadronizer.py \
--step GEN,SIM --beamspot Realistic8TeVCollision \
--conditions START53_V7C::All \
--pileup NoPileUp \
--datamix NODATAMIXER \
--eventcontent RAWSIM \
--datatier GEN-SIM \
-n $n_events \
--filein file:$input_file \
--customise Configuration/GenProduction/randomizeSeeds.randomizeSeeds \
--customise_commands process.source.skipEvents\ =\ cms.untracked.uint32\($where_to_start\)\\nprocess.source.firstLuminosityBlock\ =\ cms.untracked.uint32\($lumi_number\)


echo "finished making GEN-SIM file"

echo "ls -lh"
ls -lh

echo "begin making DIGI file"

cmsDriver.py GEN-Fragment \
--step DIGI,L1,DIGI2RAW,HLT:7E33v2 \
--conditions START53_V19::All \
--pileup 2012_Summer_50ns_PoissonOOTPU \
--datamix NODATAMIXER \
--eventcontent RAWSIM \
--datatier GEN-SIM-RAW \
-n -1 \
--filein file:hadronizer_py_GEN_SIM.root \
--customise_commands process.mix.input.fileNames\ =\ cms.untracked.vstring\(\"file:pileup.root\"\) \
--customise Configuration/GenProduction/randomizeSeeds.randomizeSeeds

echo "finished making DIGI file"

echo "ls -lh" 
ls -lh

echo "begin making RECO file"

cmsDriver.py STEP2 \
--step RAW2DIGI,L1Reco,RECO \
--conditions START53_V19::All \
--eventcontent AODSIM \
--datatier AODSIM \
-n -1 \
--filein file:GEN-Fragment_DIGI_L1_DIGI2RAW_HLT_PU.root \
--customise Configuration/GenProduction/randomizeSeeds.randomizeSeeds

echo "finished making RECO file"

echo "ls -lh"
ls -lh

/afs/cern.ch/project/eos/installation/0.2.31/bin/eos.select cp STEP2_RAW2DIGI_L1Reco_RECO.root ${output_dir}STEP2_RAW2DIGI_L1Reco_RECO_$where_to_start.root

echo date
date
