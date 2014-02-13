input_file=$1
n_events=$2

echo \$input_file
echo $input_file

echo \$n_events
echo $n_events

echo pwd
pwd

echo date
date

export SCRAM_ARCH=slc5_amd64_gcc462

cd /afs/cern.ch/work/a/anlevin/UserCode/mc_production/CMSSW_5_3_13_patch3/src/

eval `scramv1 runtime -sh`

cd -;

echo \$CMSSW_BASE
echo $CMSSW_BASE

echo pwd 
pwd

cmsDriver.py Configuration/GenProduction/python/EightTeV/hadronizer_matching.py \
--step GEN --beamspot Realistic8TeVCollision \
--conditions START53_V7C::All \
--pileup NoPileUp \
--datamix NODATAMIXER \
--eventcontent RAWSIM \
--datatier GEN \
-n $n_events \
--filein file:$input_file \
--customise Configuration/GenProduction/randomizeSeeds.randomizeSeeds 

