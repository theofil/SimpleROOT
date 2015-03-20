# SimpleROOT

This is a very simple MiniAOD analyzer producing a flat ROOT TTree for physics analysis,
of leptons, photons, jets and MET.

<h4> How to install this at CERN:  </h4>

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

<h4> The source file: </h4>

https://github.com/theofil/SimpleROOT/blob/master/plugins/SimpleROOT.cc

<h4> The config file: </h4>

https://github.com/theofil/SimpleROOT/blob/master/runme_cfg.py

<h4> How to produce a test ntuple at CERN: </h4>

cd CMSSW_7_3_2/src/Tools/SimpleROOT
cmsRun runme_cfg.py


