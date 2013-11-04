for dir in `ls | grep 2019`;
  do echo $dir;
  ls -lh $dir/harvestRun1/PSetDump.py
  cat $dir/harvestRun1/PSetDump.py | grep fileNames | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" >& tmp3.txt
  cat tmp3.txt | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" >& tmp2.txt
  replacement_filename=`cat tmp2.txt`
  echo $replacement_filename 
  cat PSetDump_19.py | sed "s/fileNames = cms.untracked.vstring('file:step3_inDQM.root')/$replacement_filename/" >& tmp.txt
  cat $dir/harvestRun1/PSetDump.py | grep workflow | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLL\
LLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" >& tmp3.txt
  cat tmp3.txt | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//"  >& tmp2.txt
  workflow_name=`cat tmp2.txt`
  echo $workflow_name
  cat tmp.txt | sed "s/workflow = cms.untracked.string('\/Global\/CMSSW_X_Y_Z\/RECO'),/$workflow_name/" >& PSetDump_19_$dir.py
done

for dir in `ls | grep 2017`;
  do echo $dir;
  ls -lh $dir/harvestRun1/PSetDump.py
  cat $dir/harvestRun1/PSetDump.py | grep fileNames | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" >& tmp3.txt
  cat tmp3.txt | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" >& tmp2.txt
  replacement_filename=`cat tmp2.txt`
  echo $replacement_filename 
  cat PSetDump_17.py | sed "s/fileNames = cms.untracked.vstring('file:step3_inDQM.root')/$replacement_filename/" >& tmp.txt
  cat $dir/harvestRun1/PSetDump.py | grep workflow | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLL\
LLLLLLLLLLLLLLLLLLL/" | sed "s/\//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/" >& tmp3.txt
  cat tmp3.txt | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//" | sed "s/LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL/\\\\\//"  >& tmp2.txt
  workflow_name=`cat tmp2.txt`
  echo $workflow_name
  cat tmp.txt | sed "s/workflow = cms.untracked.string('\/Global\/CMSSW_X_Y_Z\/RECO'),/$workflow_name/" >& PSetDump_17_$dir.py
done
