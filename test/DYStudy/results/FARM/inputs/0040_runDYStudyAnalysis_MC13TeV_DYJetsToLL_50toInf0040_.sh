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
cd /grid_mnt/storage__localgrid/delannoy/HZZ2L2NU/PUJetID/CompaCinematique/CMSSW_7_4_4/src/UserCode/llvv_fwk/test/DYStudy
eval `scramv1 runtime -sh`
ulimit -c 0;

runDYStudyAnalysis /grid_mnt/storage__localgrid/delannoy/HZZ2L2NU/PUJetID/CompaCinematique/CMSSW_7_4_4/src/UserCode/llvv_fwk/test/DYStudy/results/MC13TeV_DYJetsToLL_50toInf_40_cfg.py 
