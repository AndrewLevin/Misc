dir=/castor/cern.ch/user/p/paus/filefi/029/s12-zz2e2m-powheg-v7a
for file in `cat /home/mitprod/catalog/cern/filefi/029/s12-zz2e2m-powheg-v7a/Files | awk '{print $2}'`; do
	  staged=`stager_qry -M $dir/$file`
	  staged=`echo $staged | awk '{print $3}'`
	  if [ "$staged" == "STAGED" ]; then    
	      echo "$dir/$file staged" 
		elif [ "$staged" == "STAGEIN" ]; then
	      echo "$dir/$file stagein"
	  else
	      echo -n "            -- staging file: ($dir/$file)"
	      get=`stager_get -M $dir/$file`
	      if echo $get | grep "SUBREQUEST_FAILED"; then echo " SUBREQUEST_FAILED"; else echo ""; fi
	  fi
done
