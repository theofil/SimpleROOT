# SimpleROOT

This is a MiniAOD analyzer producing a flat ROOT TTree for physics analysis,
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

cd CMSSW_7_3_2/src/Tools/SimpleROOT <br>
cmsRun runme_cfg.py

<h4> How an ntuple looks like ? </h4>

Download this <a href="http://theofil.web.cern.ch/theofil/get/output.root">file</a>, which is an example rootuple. 

<pre>
<code>
root -l output.root 
demo->cd() // The ROOT tree is inside the demo directory 

TCut sel_SF("(lepID[0]*lepID[1] == -11*11 || lepID[0]*lepID[1] == -13*13)")
TCut sel_lep2020("l1l2DR>0.3 && lepPt[1]>20")
events->Draw("l1l2M",sel_SF && sel_lep2020,"hist") 
events->Show(3)
</code>
</pre>


<pre>
<samp>
// ttH event dumb, with h->bb (gen b matched to reco jets) 
// and tt-> WWbb -> tau nu mu nu bb -> el nu nu mu nu bb
======> EVENT:5460
 goodVtx         = 1
 nVtx            = 14
 eventNum        = 8661
 runNum          = 1
 lumi            = 87
 l1l2M           = 104.73
 l1l2Pt          = 66.7505
 l1l2Eta         = 0.0346529
 l1l2Phi         = 2.14853
 l1l2DPhi        = -1.70375
 l1l2DR          = 2.22085
 nleps           = 3
 lepPt           = 59.1267, 
                  39.7924, 20.5614
 lepEta          = -0.573448, 
                  0.851136, -0.962983
 lepPhi          = 1.51639, 
                  -3.06305, -2.98083
 lepM            = -0.0103195, 
                  0.13957, 0.13957
 lepIso          = 0.0268099, 
                  0.0265385, 0.105136
 lepPtRel        = 35.5471, 
                  48.8876, 13.8101
 lepID           = -11, 
                  13, 13
 lepGenMatchIndex = 0, 
                  1, 2
 lepTriggerMatch = 1, 
                  1, 1
 njets           = 4
 jetPt           = 181.486, 
                  151.462, 123.325, 56.1344
 jetEta          = 0.474995, 
                  -0.00419258, -0.386739, 1.6311
 jetPhi          = -0.787108, 
                  1.51442, -2.99017, -0.772743
 jetM            = 14.83, 
                  14.4056, 19.6194, 8.24497
 jetBTag         = 0, 
                  0, 0, 0
 nrjets          = 2
 rjetPt          = 64.4035, 
                  44.9823
 rjetEta         = -0.579714, 
                  0.859564
 rjetPhi         = 1.5235, 
                  -3.05334
 rjetM           = 4.81469, 
                  5.26755
 rjetBTag        = 0, 
                  0
 ngenleps        = 4
 genlepPt        = 60.8975, 
                  40.284, 20.9083, 13.4696
 genlepEta       = -0.573406, 
                  0.85066, -0.962838, -0.0405636
 genlepPhi       = 1.51499, 
                  -3.06309, -2.98075, 1.4617
 genlepM         = 0.000571982, 
                  0.105658, 0.10566, 0.000511
 genlepID        = -11, 
                  13, 13, 11
 genlepMID       = -15, 
                  -24, -521, -511
 genlepGMID      = 24, 
                  -6, -523, -513
 genlepGGMID     = 6, 
                  2212, 5, 5
 genl1l2M        = 106.949
 genl1l2Pt       = 68.3688
 genl1l2Eta      = 0.0246985
 genl1l2Phi      = 2.13852
 genl1l2DPhi     = -1.70511
 genl1l2DR       = 2.22157
 nphos           = 2
 met             = 51.6295
 metPhi          = -0.547856
 genmet          = 53.396
 t1met           = 48.2903
 t1metPhi        = -0.465397
 sumEt           = 1362.14
 t1sumEt         = 1398.48
 rho             = 9.07665
 nPU             = 17
 nPUTrue         = 20
 vHT             = 30.0087
 t1vHT           = 34.8948
 jvHT            = 142.02
 HLT_e1e2        = 0
 HLT_mu1mu2      = 1
 HLT_mu1e2       = 1
 HLT_e1mu2       = 1
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
 Flag_METFilters = 1
 Flag_HBHENoiseFilter = 1
 Flag_trkPOG_toomanystripclus53X = 1
 Flag_hcalLaserEventFilter = 1
 isDYTauTau      = 0
 ngenparts       = 5
 genpartPt       = 228.19, 
                  317.025, 167.13, 153.759, 77.1921
 genpartEta      = 0.892705, 
                  -0.381775, 0.102876, -0.735381, 0.69443
 genpartPhi      = -0.800127, 
                  1.33867, -2.54137, 1.15584, -2.003
 genpartM        = 125, 
                  172.589, 173.982, 78.0383, 80.7928
 genpartID       = 25, 
                  6, -6, 24, -24
 genpartDID1     = 5, 
                  5, -5, 16, -14
 genpartDID2     = -5, 
                  24, -24, -15, 13
 genpartDRMI1    = 3, 
                  1, 2, -1, -1
 genpartDRMI2    = 0, 
                  -1, -1, -1, 1

</samp>
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
<li> ntuple file size should be as small as possible </li>
</ul>


