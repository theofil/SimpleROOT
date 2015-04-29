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

Download this <a href="http://theofil.web.cern.ch/theofil/get/output.root">file</a>, which is an example rootuple. 

<pre>
<code>
root -l output.root 
demo->cd() 

TCut sel_SF("(lepID[0]*lepID[1] == -11*11 || lepID[0]*lepID[1] == -13*13)")
TCut sel_lep2020("l1l2DR>0.3 && lepPt[1]>20")
events->Draw("l1l2M",sel_SF && sel_lep2020,"hist") 
events->Show(3)
</code>
</pre>


<pre>
<samp>
======> EVENT:3
 goodVtx         = 1
 nVtx            = 23
 eventNum        = 1001324
 runNum          = 1
 lumi            = 10020
 l1l2M           = 91.1683
 l1l2Pt          = 26.7699
 l1l2Eta         = -2.59355
 l1l2Phi         = 0.80909
 l1l2DPhi        = -2.34752
 l1l2DR          = 2.9592
 nleps           = 2
 lepPt           = 37.2818, 
                  29.2314
 lepEta          = -2.19965, 
                  -0.397985
 lepPhi          = -0.0836522, 
                  2.26387
 lepM            = 0.13957, 
                  0.13957
 lepIso          = 0.00903365, 
                  0
 lepPtRel        = 119.693, 
                  31.4781
 lepID           = 13, 
                  -13
 lepGenMatchIndex = 0, 
                  1
 lepTriggerMatch = 1, 
                  1
 njets           = 1
 jetPt           = 37.2818
 jetEta          = -2.19965
 jetPhi          = -0.0836522
 jetM            = 0.13957
 jetBTag         = 0
 nrjets          = 1
 rjetPt          = 37.2818
 rjetEta         = -2.19965
 rjetPhi         = -0.0836522
 rjetM           = 0.13957
 rjetBTag        = 0
 ngenleps        = 2
 genlepPt        = 36.8366, 
                  30.06
 genlepEta       = -2.19914, 
                  -0.398317
 genlepPhi       = -0.0839442, 
                  2.26368
 genlepM         = 0.105665, 
                  0.105658
 genlepID        = 13, 
                  -13
 genlepMID       = 23, 
                  23
 genlepGMID      = -3, 
                  -3
 genlepGGMID     = 2212, 
                  2212
 genl1l2M        = 91.8685
 genl1l2Pt       = 26.6092
 genl1l2Eta      = -2.58991
 genl1l2Phi      = 0.852789
 genl1l2DPhi     = -2.34762
 genl1l2DR       = 2.95877
 nphos           = 0
 met             = 26.6696
 metPhi          = -1.89484
 genmet          = 1.73535e-06
 t1met           = 28.8569
 t1metPhi        = -1.82337
 sumEt           = 1136.63
 t1sumEt         = 1164.45
 rho             = 12.8935
 nPU             = 30
 nPUTrue         = 20
 vHT             = 11.6015
 t1vHT           = 14.1531
 jvHT            = 60.6142
 HLT_e1e2        = 0
 HLT_mu1mu2      = 1
 HLT_mu1e2       = 0
 HLT_e1mu2       = 0
 HLT_pfmet       = 0
 HLT_pfmetCSV    = 0
 Flag_trackingFailureFilter = 1
 Flag_goodVertices = 1
 Flag_CSCTightHaloFilter = 1
 Flag_trkPOGFilters = 1
 Flag_trkPOG_logErrorTooManyClusters = 1
 Flag_EcalDeadCellTriggerPrimitiveFilter = 1
 Flag_ecalLaserCorrFilter = 1
 Flag_trkPOG_manystripclus53X = 1
 Flag_eeBadScFilter = 1
 Flag_METFilters = 0
 Flag_HBHENoiseFilter = 0
 Flag_trkPOG_toomanystripclus53X = 1
 Flag_hcalLaserEventFilter = 1
 isDYTauTau      = 0
</samp>
</pre>

<h4> Missing elements </h4>
<ul>
<li> trigger bits </li>
<li> event flags </li>
</ul>

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


