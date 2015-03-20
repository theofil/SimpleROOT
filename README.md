# SimpleROOT

This is a very simple MiniAOD analyzer producing a flat ROOT TTree for physics analysis,
of leptons, photons, jets and MET.

<b> How to install this at CERN:  </b>

ssh username@lxplus.cern.ch  
cd YourFavoriteDir  
setenv SCRAM_ARCH slc6_amd64_gcc491  
cmsrel CMSSW_7_3_2  
cd CMSSW_7_3_2/src  
mkdir Tools  
cd Tools/  
git clone https://github.com/theofil/SimpleROOT.git  
cd SimpleROOT/  
cmsenv  
scramv1 b  

<b> How to produce a test ntuple at CERN: </b>

cd CMSSW_7_3_2/src/Tools/SimpleROOT
cmsRun python/ConfFile_cfg.py


