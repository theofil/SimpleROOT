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

<h4> How an ntuple looks like ? </h4>

Download this <a href="http://theofil.web.cern.ch/theofil/get/output.root"> file <a/> 

<pre>
<code>
root -l output.root 
demo->cd() 
events->Draw("l1l2M"," l1l2DR>0.3 && lepPt[1]>20 && abs(lepEta[0])<1.4 && abs(lepEta[1])<1.4 && (lepID[0]*lepID[1] == -11*11 || lepID[0]*lepID[1] == -13*13)","hist") <br>
</code>
</pre>

<h4> Dogma </h4>
<ul>
<li> code should be intuitive and readable by anybody who knows English </li>
<li> skipping events for any reason is not legal, save them all on disk even if ~empty </li>
<li> be flat and use basic types for branches float, (unsigned) int,  (unsigned) short </li>
<li> TTree::Show() should produce easily readable output
<li> be lazy in introducing new structures and classes, avoid over-structuring </li>
<li> use what is already there in edm, reco, std & ROOT namespaces</li>
<li> be laconic, use smart commands when possible C++11 </li>
<li> event selection and logic should be undisturbed by technicalities </li>
<li> using int when short is OK is not OK, double is forbidden</li>
<li> keep number of C++ lines minimum, gather everything in a single file </b>
<li> ntiple file size should be as small as possible </li>
</ul>


