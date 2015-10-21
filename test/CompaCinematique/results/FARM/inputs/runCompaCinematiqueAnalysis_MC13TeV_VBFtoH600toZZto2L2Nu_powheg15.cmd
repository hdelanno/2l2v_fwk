####################################
#        LaunchOnFarm Script       #
#     Loic.quertenmont@cern.ch     #
#            April 2010            #
####################################

qsub  -q localgrid@cream02 -l walltime=00:30:00 /grid_mnt/storage__localgrid/delannoy/HZZ2L2NU/PUJetID/CompaCinematique/CMSSW_7_4_4/src/UserCode/llvv_fwk/test/CompaCinematique/results/FARM/inputs/0000_runCompaCinematiqueAnalysis_MC13TeV_VBFtoH600toZZto2L2Nu_powheg150000_.sh
