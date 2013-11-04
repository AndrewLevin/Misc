for file in `ls | grep __DQM`; 
  do echo $file;
  python2.6 DQMUpload_Manual.py $file https://cmsweb.cern.ch/dq\
m/relval
done    

