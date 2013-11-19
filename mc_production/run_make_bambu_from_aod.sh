for file in `/afs/cern.ch/project/eos/installation/0.2.31/bin/eos.select ls /eos/cms/store/user/anlevin/data/AOD/unmerged/ | grep 0tev-4`
  do
  bsub -q 8nh "bash make_bambu_from_aod.sh /store/user/anlevin/data/AOD/unmerged/$file /eos/cms/store/user/anlevin/data/BAMBU/"
done
