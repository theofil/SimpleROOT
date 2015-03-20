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

<h4> Axioms </h4>
<ul>
<li> code should be intuitive and readable by anybody who knows English </li>
<li> skipping events for any reason is not legal, save them all on disk even if ~empty </li>
<li> be flat and use basic types for branches float, (unsigned) int,  (unsigned) short  and STD vectors </li>
<li> be lazy in introducing new structures and classes, use what is already there in edm, reco, std, ROOT namespaces</li>
<li> be laconic and use smart commands when possible C++11 </li>
<li> event selection and logic should be undisturbed by technicalities </li>
<li> using int when short is OK is not OK, double is forbidden</li>
</ul>


