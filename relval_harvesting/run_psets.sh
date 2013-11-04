for pset in `ls PSet*py | grep -v  PSetDump_23.py`;
  do 
  echo ""
  echo "===============> begin running pset: $pset";
  echo ""

  cmsRun $pset
  
  echo ""
  echo "===============> finished running pset: $pset";
  echo ""
done 