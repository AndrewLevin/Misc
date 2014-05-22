echo "this script will produce a json file called final.json that can be run on all of the 2012 moriond data except the run 2012 c ecal recovered dataset"
echo 
echo "min and max run numbers for various datasets are found in the process of this and will also be printed out"
echo 

#to get the min and max runs for run 2012c pr v2, do this:
minrun_run2012c_pr_v2=`~/das_client.py --limit=100000 --query="run dataset=/DoubleElectron/Run2012C-PromptReco-v2/AOD | min(run.run_number)" | grep = | awk -F = '{print $2}'`
maxrun_run2012c_pr_v2=`~/das_client.py --limit=100000 --query="run dataset=/DoubleElectron/Run2012C-PromptReco-v2/AOD | max(run.run_number)" | grep = | awk -F = '{print $2}'`

echo minrun_run2012c_pr_v2=$minrun_run2012c_pr_v2
echo maxrun_run2012c_pr_v2=$maxrun_run2012c_pr_v2

python filterJSON.py --min $minrun_run2012c_pr_v2 --max $maxrun_run2012c_pr_v2 /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt --output run2012c_pr_v2_all.json

#to get the min and max runs for run 2012d pr v1, do this:
minrun_run2012d_pr_v1=`~/das_client.py --limit=100000 --query="run dataset=/DoubleElectron/Run2012D-PromptReco-v1/AOD | min(run.run_number)" | grep = | awk -F = '{print $2}'`
maxrun_run2012d_pr_v1=`~/das_client.py --limit=100000 --query="run dataset=/DoubleElectron/Run2012D-PromptReco-v1/AOD | max(run.run_number)" | grep = | awk -F = '{print $2}'`

echo minrun_run2012d_pr_v1=$minrun_run2012d_pr_v1
echo maxrun_run2012d_pr_v1=$maxrun_run2012d_pr_v1

filterJSON.py --min $minrun_run2012d_pr_v1 --max $maxrun_run2012d_pr_v1 /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt --output run2012d_pr_v1.json

python compareJSON.py --sub run2012c_pr_v2_all.json /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_201191-201191_8TeV_11Dec2012ReReco-recover_Collisions12_JSON.txt run2012c_pr_v2.json

#to get the min and max runs for the run 2012c 24 aug rereco, do this: 
minrun_run2012c_24_aug_v1=`~/das_client.py --limit=100000 --query="run dataset=/DoubleElectron/Run2012C-24Aug2012-v1/AOD | min(run.run_number)" | grep = | awk -F = '{print $2}'`
maxrun_run2012c_24_aug_v1=`~/das_client.py --limit=100000 --query="run dataset=/DoubleElectron/Run2012C-24Aug2012-v1/AOD | max(run.run_number)" | grep = | awk -F = '{print $2}'`

echo minrun_run2012c_24_aug_v1=$minrun_run2012c_24_aug_v1
echo maxrun_run2012c_24_aug_v1=$maxrun_run2012c_24_aug_v1

python filterJSON.py --min $minrun_run2012c_24_aug_v1  --max $maxrun_run2012c_24_aug_v1 /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt --output run2012c_24_aug_v1.json

#to get the min and max runs for the run 2012b 13 july rereco, do this:
minrun_run2012b_13_jul_v1=`~/das_client.py --limit=100000 --query="run dataset=/DoubleElectron/Run2012B-13Jul2012-v1/AOD | min(run.run_number)" | grep = | awk -F = '{print $2}'`
maxrun_run2012b_13_jul_v1=`~/das_client.py --limit=100000 --query="run dataset=/DoubleElectron/Run2012B-13Jul2012-v1/AOD | max(run.run_number)" | grep = | awk -F = '{print $2}'`

echo minrun_run2012b_13_jul_v1=$minrun_run2012b_13_jul_v1
echo maxrun_run2012b_13_jul_v1=$maxrun_run2012b_13_jul_v1


python filterJSON.py --min $minrun_run2012b_13_jul_v1  --max $maxrun_run2012b_13_jul_v1 /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt  --output run2012b_13_jul_v1.json

#to get the min and max runs for the run 2012a 13 july rereco, do this:
minrun_run2012a_13_jul_v1=`~/das_client.py --limit=100000 --query="run dataset=/DoubleElectron/Run2012A-13Jul2012-v1/AOD | min(run.run_number)" | grep = | awk -F = '{print $2}'`
maxrun_run2012a_13_jul_v1=`~/das_client.py --limit=100000 --query="run dataset=/DoubleElectron/Run2012A-13Jul2012-v1/AOD | max(run.run_number)" | grep = | awk -F = '{print $2}'`

echo minrun_run2012a_13_jul_v1=$minrun_run2012a_13_jul_v1
echo maxrun_run2012a_13_jul_v1=$maxrun_run2012a_13_jul_v1


python filterJSON.py --min $minrun_run2012a_13_jul_v1  --max $maxrun_run2012a_13_jul_v1 /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt  --output run2012a_13_jul_v1.json

#to get the min and max runs for the run 2012a 13 july rereco, do this:
minrun_run2012a_06_aug_v1=`~/das_client.py --limit=100000 --query="run dataset=/DoubleElectron/Run2012A-recover-06Aug2012-v1/AOD | min(run.run_number)" | grep = | awk -F = '{print $2}'`
maxrun_run2012a_06_aug_v1=`~/das_client.py --limit=100000 --query="run dataset=/DoubleElectron/Run2012A-recover-06Aug2012-v1/AOD | max(run.run_number)" | grep = | awk -F = '{print $2}'`

echo minrun_run2012a_06_aug_v1=$minrun_run2012a_06_aug_v1
echo maxrun_run2012a_06_aug_v1=$maxrun_run2012a_06_aug_v1


python filterJSON.py --min $minrun_run2012a_06_aug_v1  --max $maxrun_run2012a_06_aug_v1  /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt  --output run2012a_06_aug_v1.json

python compareJSON.py --or run2012a_06_aug_v1.json run2012a_13_jul_v1.json temp1.json
python compareJSON.py --or temp1.json run2012b_13_jul_v1.json temp2.json
python compareJSON.py --or temp2.json run2012c_24_aug_v1.json temp3.json
python compareJSON.py --or temp3.json run2012c_pr_v2.json temp4.json
python compareJSON.py --or temp4.json run2012d_pr_v1.json temp5.json

python filterJSON.py --runs 207883,207884,207885,207886,207887,207889,207897,207898,207905,207920,207921,207922,207924,208300,208302,208303,208304,208307 temp5.json --output final.json

#python compareJSON.py --diff final.json /afs/cern.ch/user/b/botta/public/moriond.json




