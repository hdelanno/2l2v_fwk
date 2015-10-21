#! /bin/sh
####################################
#        LaunchOnFarm Script       #
#     Loic.quertenmont@cern.ch     #
#            April 2010            #
####################################

pwd
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc491
export BUILD_ARCH=slc6_amd64_gcc491
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
cd /grid_mnt/storage__localgrid/delannoy/HZZ2L2NU/PUJetID/CompaCinematique/CMSSW_7_4_4/src/UserCode/llvv_fwk/test/CompaCinematique
eval `scramv1 runtime -sh`
ulimit -c 0;

runCompaCinematiqueAnalysis /grid_mnt/storage__localgrid/delannoy/HZZ2L2NU/PUJetID/CompaCinematique/CMSSW_7_4_4/src/UserCode/llvv_fwk/test/CompaCinematique/results/MC8TeV_VBFtoH600toZZto2L2Nu_powheg15_0_cfg.py 
