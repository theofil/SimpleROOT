// -*- C++ -*-
//
// Package:    Tools/SimpleROOT
// Class:      SimpleROOT
//
// Original Author:  Konstantinos Theofilatos
//         Created:  Fri, 20 Feb 2015 16:08:35 GMT
// with many thanks to https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"


#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"
#include "TLorentzVector.h"


#define njetsMax 30 
#define nrjetsMax 10 
#define nlepsMax 10
#define ngenlepsMax 10
#define ngenpartsMax 10

using namespace std;


class SimpleROOT : public edm::EDAnalyzer {
    public:
    	explicit SimpleROOT(const edm::ParameterSet&);
      	~SimpleROOT();

    private:

	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual bool isGoodMuon(const pat::Muon &mu);
        virtual bool isGoodElectron(const pat::Electron &el);
        virtual bool isGoodVertex(const reco::Vertex &PV);
        virtual bool isGoodJet(const pat::Jet &myJet);
        virtual bool isGoodPhoton(const pat::Photon &myPhoton);

        float MuonRelIso(const reco::Candidate *cand);
        float ElectronRelIso(const reco::Candidate *cand);
        float LeptonRelIso(const reco::Candidate *cand){return cand->isElectron() ? ElectronRelIso(cand) : MuonRelIso(cand);}
        bool  isGoodLepton(const reco::Candidate *cand){return cand->isElectron() ? isGoodElectron(*((pat::Electron*)cand)) : isGoodMuon(*((pat::Muon*)cand));}
  
        float PtRel(const reco::Candidate * myLepton, vector<const reco::Candidate *> myJets);
        short getMatchedIndex(const reco::Candidate *, vector<const reco::Candidate *>);
        bool  triggerMatch(const reco::Candidate *, const std::vector<TLorentzVector> &, const std::vector<TLorentzVector> &);
        const reco::Candidate *getGenMother(const reco::Candidate*); 
       
        float getPFIsolation(edm::Handle<pat::PackedCandidateCollection>, const reco::Candidate*, float, float, float, bool, bool); // miniISO 

        TLorentzVector P4(const reco::Candidate* cand){TLorentzVector p4vec; p4vec.SetPxPyPzE( cand->px(), cand->py(), cand->pz(), cand->energy() ); return p4vec;}

    	edm::Handle<reco::VertexCollection> vertices;
   	edm::Handle<pat::MuonCollection> muons;
     	edm::Handle<pat::ElectronCollection> electrons;
        edm::Handle<pat::JetCollection> jets;
        edm::Handle<pat::METCollection> mets;
        edm::Handle<pat::PhotonCollection> photons;
        edm::Handle<double> rhoH;
        edm::Handle<edm::View<PileupSummaryInfo>>  pileup;
        edm::Handle<edm::View<reco::GenParticle>> pruned;
        edm::Handle<edm::View<pat::PackedGenParticle>> packed;
        edm::Handle<pat::PackedCandidateCollection> pfcands;
        edm::Handle<edm::TriggerResults> triggerBits;
        edm::Handle<edm::TriggerResults> filterBits;
        edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
        edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

        
    	edm::EDGetTokenT<reco::VertexCollection> verticesToken;
   	edm::EDGetTokenT<pat::MuonCollection> muonsToken;
     	edm::EDGetTokenT<pat::ElectronCollection> electronsToken;
        edm::EDGetTokenT<pat::JetCollection> jetsToken;
        edm::EDGetTokenT<pat::METCollection> metsToken;
        edm::EDGetTokenT<pat::PhotonCollection> photonsToken;
        edm::EDGetTokenT<double> rhoHToken;
        edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken;
        edm::EDGetTokenT<edm::View<reco::GenParticle>> prunedToken;
        edm::EDGetTokenT<edm::View<pat::PackedGenParticle>> packedToken;
        edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken;
        edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
        edm::EDGetTokenT<edm::TriggerResults> filterBits_;
        edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
        edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

        reco::Vertex vtx; // stores event's primary vertex

	edm::Service<TFileService> fileService_; 
	TTree *events_;
	
        // --- all ntuple vars end with "_" --- 
	bool goodVtx_; 
	unsigned short nVtx_; 

        unsigned long eventNum_;
        unsigned int runNum_;
        unsigned int lumi_;

        float l1l2DPhi_;
        float l1l2DR_;
        float l1l2Pt_;
        float l1l2M_;
        float l1l2Eta_;
        float l1l2Phi_;
      
        unsigned short nleps_;
        float lepPt_                      [nlepsMax]; 
        float lepEta_                     [nlepsMax]; 
        float lepPhi_                     [nlepsMax]; 
        float lepM_                       [nlepsMax];   
        float lepIso_                     [nlepsMax];
        float lepPtRel_                   [nlepsMax];
        short lepID_                      [nlepsMax];
        short lepGenMatchIndex_           [nlepsMax];   // index in the stored genLep[] collection
        bool  lepTriggerMatch_            [nlepsMax];

        unsigned short njets_;
        float jetPt_                      [njetsMax]; 
        float jetEta_                     [njetsMax]; 
        float jetPhi_                     [njetsMax]; 
        float jetM_                       [njetsMax]; 
        float jetBTag_                    [njetsMax];

        unsigned short nrjets_;          
        float rjetPt_                     [nrjetsMax]; 
        float rjetEta_                    [nrjetsMax]; 
        float rjetPhi_                    [nrjetsMax]; 
        float rjetM_                      [nrjetsMax]; 
        float rjetBTag_                   [nrjetsMax];

        unsigned short ngenleps_;
        float genlepPt_                   [ngenlepsMax]; 
        float genlepEta_                  [ngenlepsMax]; 
        float genlepPhi_                  [ngenlepsMax]; 
        float genlepM_                    [ngenlepsMax];   
        short genlepID_                   [ngenlepsMax];
        int   genlepMID_                  [ngenlepsMax]; // mother 
        int   genlepGMID_                 [ngenlepsMax]; // grand mother
        int   genlepGGMID_                [ngenlepsMax]; // grand grand mother

        float genl1l2DPhi_;
        float genl1l2DR_;
        float genl1l2Pt_;
        float genl1l2M_;
        float genl1l2Eta_;
        float genl1l2Phi_;

        unsigned short ngenparts_;
        float genpartPt_                   [ngenpartsMax]; 
        float genpartEta_                  [ngenpartsMax]; 
        float genpartPhi_                  [ngenpartsMax]; 
        float genpartM_                    [ngenpartsMax];   
        int   genpartID_                   [ngenpartsMax]; 
        int   genpartDID1_                 [ngenpartsMax]; // daugther 1
        int   genpartDID2_                 [ngenpartsMax]; // daugther 2
        short genpartDRMI1_                [ngenpartsMax]; // daugther reco match index (genleptons to leptons, genjet to jet)
        short genpartDRMI2_                [ngenpartsMax]; // daugther reco match index 

        float met_;          // raw pf-met
        float metPhi_;
        float t1met_;
        float t1metPhi_;
        float genmet_;          
        float sumEt_;
        float t1sumEt_;

        float rho_;
        float nPU_;
        float nPUTrue_;

        float vHT_;           // pt of the recoil vector of all objects excluding the dilepton [= -met - l1l2] 
        float t1vHT_;         // as vHT but with t1 correction
	float jvHT_;          // recoil of hard jets and subleading leptons

	bool HLT_e1e2_; 
	bool HLT_mu1mu2_; 
	bool HLT_mu1e2_; 
	bool HLT_e1mu2_; 
	bool HLT_pfmet_; 
	bool HLT_pfmetCSV_; 
       
	bool Flag_trackingFailureFilter_;		        
	bool Flag_goodVertices_;			 
	bool Flag_CSCTightHaloFilter_;		 
	bool Flag_trkPOGFilters_;			 
	bool Flag_trkPOG_logErrorTooManyClusters_;	 
	bool Flag_EcalDeadCellTriggerPrimitiveFilter_; 
	bool Flag_ecalLaserCorrFilter_;		 
	bool Flag_trkPOG_manystripclus53X_;		 
	bool Flag_eeBadScFilter_;			 
	bool Flag_METFilters_;			 
	bool Flag_HBHENoiseFilter_;			 
	bool Flag_trkPOG_toomanystripclus53X_;	 
	bool Flag_hcalLaserEventFilter_;		 

	bool isDYTauTau_; 

        unsigned short nphos_;
};

SimpleROOT::SimpleROOT(const edm::ParameterSet& iConfig):  // initialize tokens in the constructor, better performance
verticesToken(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter("offlineSlimmedPrimaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices")))),
muonsToken(consumes<pat::MuonCollection>(iConfig.getUntrackedParameter("slimmedMuons",edm::InputTag("slimmedMuons")))),
electronsToken(consumes<pat::ElectronCollection>(iConfig.getUntrackedParameter("slimmedElectrons",edm::InputTag("slimmedElectrons")))),
jetsToken(consumes<pat::JetCollection>(iConfig.getUntrackedParameter("slimmedJets",edm::InputTag("slimmedJets")))),
metsToken(consumes<pat::METCollection>(iConfig.getUntrackedParameter("slimmedMETs",edm::InputTag("slimmedMETs")))),
photonsToken(consumes<pat::PhotonCollection>(iConfig.getUntrackedParameter("slimmedPhotons",edm::InputTag("slimmedPhotons")))),
rhoHToken(consumes<double>(iConfig.getUntrackedParameter("fixedGridRhoFastjetAll",edm::InputTag("fixedGridRhoFastjetAll")))),
pileupToken(consumes<edm::View<PileupSummaryInfo> >(iConfig.getUntrackedParameter("addPileupInfo", edm::InputTag("addPileupInfo")))),  
prunedToken(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter("prunedGenParticles", edm::InputTag("prunedGenParticles")))),
triggerBits_(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"))),
filterBits_(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","PAT"))),
triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("selectedPatTrigger"))),
triggerPrescales_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger")))
//packedToken(consumes<edm::View<pat::PackedGenParticle> >(edm::InputTag("packedGenParticles"))),
//pfcandsToken(consumes<pat::PackedCandidateCollection> (iConfig.getUntrackedParameter("packedPFCandidates", edm::InputTag("packedPFCandidates"))))
//Token(consumes<>(iConfig.getUntrackedParameter("",edm::InputTag("")))),  // first arg is default the second is used only if is defined in runme_cfg.py
{
    events_ = fileService_->make<TTree>("events","events");
    events_->Branch("goodVtx"          ,&goodVtx_               ,"goodVtx/O");
    events_->Branch("nVtx"             ,&nVtx_                  ,"nVtx/s");
    events_->Branch("eventNum"         ,&eventNum_              ,"eventNum/l");
    events_->Branch("runNum"           ,&runNum_                ,"runNum/i");
    events_->Branch("lumi"             ,&lumi_                  ,"lumi/i");

    events_->Branch("l1l2M"            ,&l1l2M_                 ,"l1l2M/F");
    events_->Branch("l1l2Pt"           ,&l1l2Pt_                ,"l1l2Pt/F");
    events_->Branch("l1l2Eta"          ,&l1l2Eta_               ,"l1l2Eta/F");
    events_->Branch("l1l2Phi"          ,&l1l2Phi_               ,"l1l2Phi/F");
    events_->Branch("l1l2DPhi"         ,&l1l2DPhi_              ,"l1l2DPhi/F");
    events_->Branch("l1l2DR"           ,&l1l2DR_                ,"l1l2DR/F");

    events_->Branch("nleps"            ,&nleps_                 ,"nleps/s");
    events_->Branch("lepPt"            ,lepPt_                  ,"lepPt[nleps]/F");
    events_->Branch("lepEta"           ,lepEta_                 ,"lepEta[nleps]/F");
    events_->Branch("lepPhi"           ,lepPhi_                 ,"lepPhi[nleps]/F");
    events_->Branch("lepM"             ,lepM_                   ,"lepM[nleps]/F");
    events_->Branch("lepIso"           ,lepIso_                 ,"lepIso[nleps]/F");
    events_->Branch("lepPtRel"         ,lepPtRel_               ,"lepPtRel[nleps]/F");
    events_->Branch("lepID"            ,lepID_                  ,"lepID[nleps]/S");
    events_->Branch("lepGenMatchIndex" ,lepGenMatchIndex_       ,"lepGenMatchIndex[nleps]/S");
    events_->Branch("lepTriggerMatch"  ,lepTriggerMatch_        ,"lepTriggerMatch[nleps]/O");

    events_->Branch("njets"            ,&njets_                 ,"njets/s");
    events_->Branch("jetPt"            ,jetPt_                  ,"jetPt[njets]/F");
    events_->Branch("jetEta"           ,jetEta_                 ,"jetEta[njets]/F");
    events_->Branch("jetPhi"           ,jetPhi_                 ,"jetPhi[njets]/F");
    events_->Branch("jetM"             ,jetM_                   ,"jetM[njets]/F");
    events_->Branch("jetBTag"          ,jetBTag_                ,"jetBTag[njets]/F");

    events_->Branch("nrjets"           ,&nrjets_                ,"nrjets/s");
    events_->Branch("rjetPt"           ,rjetPt_                 ,"rjetPt[nrjets]/F");
    events_->Branch("rjetEta"          ,rjetEta_                ,"rjetEta[nrjets]/F");
    events_->Branch("rjetPhi"          ,rjetPhi_                ,"rjetPhi[nrjets]/F");
    events_->Branch("rjetM"            ,rjetM_                  ,"rjetM[nrjets]/F");
    events_->Branch("rjetBTag"         ,rjetBTag_               ,"rjetBTag[nrjets]/F");

    events_->Branch("ngenleps"         ,&ngenleps_              ,"ngenleps/s");
    events_->Branch("genlepPt"         ,genlepPt_               ,"genlepPt[ngenleps]/F");
    events_->Branch("genlepEta"        ,genlepEta_              ,"genlepEta[ngenleps]/F");
    events_->Branch("genlepPhi"        ,genlepPhi_              ,"genlepPhi[ngenleps]/F");
    events_->Branch("genlepM"          ,genlepM_                ,"genlepM[ngenleps]/F");
    events_->Branch("genlepID"         ,genlepID_               ,"genlepID[ngenleps]/S");
    events_->Branch("genlepMID"        ,genlepMID_              ,"genlepMID[ngenleps]/I");
    events_->Branch("genlepGMID"       ,genlepGMID_             ,"genlepGMID[ngenleps]/I");
    events_->Branch("genlepGGMID"      ,genlepGGMID_            ,"genlepGGMID[ngenleps]/I");

    events_->Branch("genl1l2M"         ,&genl1l2M_              ,"genl1l2M/F");
    events_->Branch("genl1l2Pt"        ,&genl1l2Pt_             ,"genl1l2Pt/F");
    events_->Branch("genl1l2Eta"       ,&genl1l2Eta_            ,"genl1l2Eta/F");
    events_->Branch("genl1l2Phi"       ,&genl1l2Phi_            ,"genl1l2Phi/F");
    events_->Branch("genl1l2DPhi"      ,&genl1l2DPhi_           ,"genl1l2DPhi/F");
    events_->Branch("genl1l2DR"        ,&genl1l2DR_             ,"genl1l2DR/F");

    events_->Branch("nphos"            ,&nphos_                 ,"nphos/s");
    events_->Branch("met"              ,&met_                   ,"met/F");
    events_->Branch("metPhi"           ,&metPhi_                ,"metPhi/F");
    events_->Branch("genmet"           ,&genmet_                ,"genmet/F");
    events_->Branch("t1met"            ,&t1met_                 ,"t1met/F");
    events_->Branch("t1metPhi"         ,&t1metPhi_              ,"t1metPhi/F");
    events_->Branch("sumEt"            ,&sumEt_                 ,"sumEt/F");
    events_->Branch("t1sumEt"          ,&t1sumEt_               ,"t1sumEt/F");

    events_->Branch("rho"              ,&rho_                   ,"rho/F");
    events_->Branch("nPU"              ,&nPU_                   ,"nPU/F");
    events_->Branch("nPUTrue"          ,&nPUTrue_               ,"nPUTrue/F");

    events_->Branch("vHT"              ,&vHT_                   ,"vHT/F      ");
    events_->Branch("t1vHT"            ,&t1vHT_                 ,"t1vHT/F    ");
    events_->Branch("jvHT"             ,&jvHT_                  ,"jvHT/F     ");

    events_->Branch("HLT_e1e2"         ,&HLT_e1e2_              ,"HLT_e1e2/O");
    events_->Branch("HLT_mu1mu2"       ,&HLT_mu1mu2_            ,"HLT_mu1mu2/O");
    events_->Branch("HLT_mu1e2"        ,&HLT_mu1e2_             ,"HLT_mu1e2/O");
    events_->Branch("HLT_e1mu2"        ,&HLT_e1mu2_             ,"HLT_e1mu2/O");
    events_->Branch("HLT_pfmet"        ,&HLT_pfmet_             ,"HLT_pfmet/O");
    events_->Branch("HLT_pfmetCSV"     ,&HLT_pfmetCSV_          ,"HLT_pfmetCSV/O");

    events_->Branch("Flag_trackingFailureFilter"                   ,&Flag_trackingFailureFilter_	                ,"Flag_trackingFailureFilter/O");
    events_->Branch("Flag_goodVertices"                            ,&Flag_goodVertices_			                ,"Flag_goodVertices/O");
    events_->Branch("Flag_CSCTightHaloFilter"                      ,&Flag_CSCTightHaloFilter_		                ,"Flag_CSCTightHaloFilter/O");
    events_->Branch("Flag_trkPOGFilters"                           ,&Flag_trkPOGFilters_		                ,"Flag_trkPOGFilters/O");
    events_->Branch("Flag_trkPOG_logErrorTooManyClusters"          ,&Flag_trkPOG_logErrorTooManyClusters_	        ,"Flag_trkPOG_logErrorTooManyClusters/O");
    events_->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter"      ,&Flag_EcalDeadCellTriggerPrimitiveFilter_           ,"Flag_EcalDeadCellTriggerPrimitiveFilter/O");
    events_->Branch("Flag_ecalLaserCorrFilter"                     ,&Flag_ecalLaserCorrFilter_		                ,"Flag_ecalLaserCorrFilter/O");
    events_->Branch("Flag_trkPOG_manystripclus53X"                 ,&Flag_trkPOG_manystripclus53X_		        ,"Flag_trkPOG_manystripclus53X/O");
    events_->Branch("Flag_eeBadScFilter"                           ,&Flag_eeBadScFilter_		                ,"Flag_eeBadScFilter/O");
    events_->Branch("Flag_METFilters"                              ,&Flag_METFilters_			                ,"Flag_METFilters/O");
    events_->Branch("Flag_HBHENoiseFilter"                         ,&Flag_HBHENoiseFilter_			        ,"Flag_HBHENoiseFilter/O");
    events_->Branch("Flag_trkPOG_toomanystripclus53X"              ,&Flag_trkPOG_toomanystripclus53X_	                ,"Flag_trkPOG_toomanystripclus53X/O");
    events_->Branch("Flag_hcalLaserEventFilter"                    ,&Flag_hcalLaserEventFilter_		                ,"Flag_hcalLaserEventFilter/O");

    events_->Branch("isDYTauTau"         ,&isDYTauTau_              ,"isDYTauTau/O");

    events_->Branch("ngenparts"             ,&ngenparts_                ,"ngenparts/s");
    events_->Branch("genpartPt"             ,genpartPt_                 ,"genpartPt[ngenparts]/F");
    events_->Branch("genpartEta"            ,genpartEta_                ,"genpartEta[ngenparts]/F");
    events_->Branch("genpartPhi"            ,genpartPhi_                ,"genpartPhi[ngenparts]/F");
    events_->Branch("genpartM"              ,genpartM_                  ,"genpartM[ngenparts]/F");
    events_->Branch("genpartID"             ,genpartID_                 ,"genpartID[ngenparts]/I");
    events_->Branch("genpartDID1"           ,genpartDID1_               ,"genpartDID1[ngenparts]/I");
    events_->Branch("genpartDID2"           ,genpartDID2_               ,"genpartDID2[ngenparts]/I");
    events_->Branch("genpartDRMI1"          ,genpartDRMI1_              ,"genpartDRMI1[ngenparts]/S");
    events_->Branch("genpartDRMI2"          ,genpartDRMI2_              ,"genpartDRMI2[ngenparts]/S");
    //events_->Branch(""          ,&               ,"");
}


SimpleROOT::~SimpleROOT() {}

// ------------ method called for each event  ------------
void SimpleROOT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    iEvent.getByToken(verticesToken, vertices);  
    iEvent.getByToken(muonsToken, muons);
    iEvent.getByToken(electronsToken, electrons);
    iEvent.getByToken(jetsToken, jets);
    iEvent.getByToken(metsToken, mets);
    iEvent.getByToken(photonsToken, photons);
    iEvent.getByToken(rhoHToken, rhoH);
    iEvent.getByToken(pileupToken, pileup);
    iEvent.getByToken(prunedToken, pruned);
//    iEvent.getByToken(packedToken, packed);   // commented out because not needed at the moment
//    iEvent.getByToken(pfcandsToken, pfcands);
    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(filterBits_, filterBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
    iEvent.getByToken(triggerPrescales_, triggerPrescales);

    vector<const reco::Candidate *> myLeptons; // in this container we will store all selected RECO electrons and RECO muons 
    vector<const reco::Candidate *> myJets; // in this container we will store all prompt jets (PV)
    vector<const reco::Candidate *> myRJets; // in this container we will all rejected prompt jets (PV) -- due to DR matching with myLeptons
    vector<const reco::Candidate *> myPhotons; // in this container we will store all photons

    vector<const reco::Candidate *> myGenLeptons;
    vector<const reco::Candidate *> myGenParticles; // store here interesting particles, like top, W, Z, higgs

    vector<TLorentzVector> hltEgammaCandidates;
    vector<TLorentzVector >hltL3MuonCandidates;

    // --- trigger info
    bool HLT_e1e2(false);
    bool HLT_mu1mu2(false);
    bool HLT_mu1e2(false);
    bool HLT_e1mu2(false);
    bool HLT_pfmet(false);
    bool HLT_pfmetCSV(false);

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) 
    {
	string trigger_name = string(names.triggerName(i));
        trigger_name.pop_back();

        if(trigger_name == string("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v")                       && triggerBits->accept(i)) HLT_e1e2     = true; 
        if(trigger_name == string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")                         && triggerBits->accept(i)) HLT_mu1mu2   = true;         
        if(trigger_name == string("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v") && triggerBits->accept(i)) HLT_mu1e2    = true;         
        if(trigger_name == string("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v")  && triggerBits->accept(i)) HLT_e1mu2    = true; 
        if(trigger_name == string("HLT_PFMET170_NoiseCleaned_v")                                && triggerBits->accept(i)) HLT_pfmet    = true; 
        if(trigger_name == string("HLT_PFMET120_NoiseCleaned_BTagCSV07_v")                      && triggerBits->accept(i)) HLT_pfmetCSV = true; 
    }


   // --- store HLT trigger objects for offline matching
   if(HLT_e1e2 || HLT_mu1mu2 || HLT_mu1e2 || HLT_e1mu2)
   {
	for (pat::TriggerObjectStandAlone obj : *triggerObjects) 
	{
	    obj.unpackPathNames(names);
            if(obj.collection() == "hltL3MuonCandidates::HLT" || obj.collection() == "hltEgammaCandidates::HLT")
            {
		float pt  = obj.pt(); 
                float eta = obj.eta(); 
                float phi = obj.phi();
                TLorentzVector myObj; 
                myObj.SetPtEtaPhiM(pt, eta, phi, 0);
 
	        if(obj.collection() == "hltL3MuonCandidates::HLT") hltL3MuonCandidates.push_back(myObj);
	        if(obj.collection() == "hltEgammaCandidates::HLT") hltEgammaCandidates.push_back(myObj);
	    }
	}
    }


    // --- loop inside the objects
    for (const reco::Vertex &PV : *vertices)
    {
	if( isGoodVertex(PV) ){ vtx = PV; break;} // get first good vertex 
    }

    for (const pat::Muon &mu : *muons) 
    {
	if(isGoodMuon(mu)) myLeptons.push_back(&mu);
    }
    for (const pat::Electron &el : *electrons)
    {
	if(isGoodElectron(el)) myLeptons.push_back(&el);
    }

    for(const pat::Jet &myjet : *jets)
    {
	bool isLeptonMatched = false;
        float DRmax = 0.4;
        for(auto & lep: myLeptons) if( P4(lep).DeltaR( P4(&myjet) ) < DRmax ) isLeptonMatched = true;

	if( isGoodJet(myjet) && !isLeptonMatched ) myJets.push_back(&myjet);
	if( isGoodJet(myjet) && isLeptonMatched ) myRJets.push_back(&myjet);
    }

    for (const pat::Photon &photon : *photons) 
    {
	if( isGoodPhoton(photon) ) myPhotons.push_back(&photon); 
    }

    isDYTauTau_ = false;

    for (const reco::GenParticle & genParticle: *pruned)
    {
	int ID     = genParticle.pdgId();
        float pt   = genParticle.pt() ;
        float eta  = genParticle.eta();
        int status = genParticle.status();
        int nod    = genParticle.numberOfDaughters();


        if(pt > 10 && fabs(eta) < 2.4 && (abs(ID) == 11 || abs(ID) == 13) && status == 1)myGenLeptons.push_back(&genParticle); 
	if(pt > 5 && (abs(ID) == 15) && !isDYTauTau_ && getGenMother(&genParticle)->pdgId() == 23)isDYTauTau_ = true;
 
        if(nod == 2) // save 1->2 decays of interesting particles
        if(ID == 23 || abs(ID)==24 || abs(ID)==25 || abs(ID) == 6)
        {
            myGenParticles.push_back(&genParticle);
	}
    }



    const pat::MET &met = mets->front();
    float rawmet      = met.shiftedPt(pat::MET::NoShift, pat::MET::Raw);
    float rawmetPhi   = met.shiftedPhi(pat::MET::NoShift, pat::MET::Raw);
    float genmet      = met.genMET()->pt();
    float rawmetSumEt = met.shiftedSumEt(pat::MET::NoShift, pat::MET::Raw);
    float t1met       = met.pt();
    float t1metPhi    = met.phi();
    float t1metSumEt  = met.sumEt();

    //  --- sort by pt all objects, the [] is a C++11 lambda func 
    std::sort(myLeptons.begin(), myLeptons.end(), [](const reco::Candidate * a, const reco::Candidate * b){return a->pt() > b->pt();} );
    std::sort(myJets.begin(), myJets.end(), [](const reco::Candidate * a, const reco::Candidate * b){return a->pt() > b->pt();} );
    std::sort(myRJets.begin(), myRJets.end(), [](const reco::Candidate * a, const reco::Candidate * b){return a->pt() > b->pt();} );
    std::sort(myPhotons.begin(), myPhotons.end(), [](const reco::Candidate * a, const reco::Candidate * b){return a->pt() > b->pt();} );
    std::sort(myGenLeptons.begin(), myGenLeptons.end(), [](const reco::Candidate * a, const reco::Candidate * b){return a->pt() > b->pt();} );

    // --- define dilepton pair

    TLorentzVector l1,l2;
    if(myLeptons.size() >=1) l1 = P4(myLeptons[0]);
    if(myLeptons.size() >=2) l2 = P4(myLeptons[1]);

    TLorentzVector genl1,genl2;
    if(myGenLeptons.size() >=1) genl1 = P4(myGenLeptons[0]);
    if(myGenLeptons.size() >=2) genl2 = P4(myGenLeptons[1]);

    TLorentzVector metVector, t1metVector;
    metVector.SetPtEtaPhiE  (rawmet  , 0, rawmetPhi, rawmet );
    t1metVector.SetPtEtaPhiE(t1met   , 0, t1metPhi , t1met  );

    TLorentzVector HTVector, t1HTVector, jHTVector; 
   
    HTVector   = -metVector -l1 -l2;
    t1HTVector = -t1metVector -l1 -l2;
    
    for(auto &myjet : myJets)    jHTVector += P4(myjet);
    for(auto &mylep : myLeptons) jHTVector -= P4(mylep);
 
    // --- get pileup info 
    float nPU = 0;
    float nPUTrue = 0;
    for( auto & pu : *pileup)
    { 
	if(pu.getBunchCrossing() == 0) 
	{
	    nPU      =  pu.getPU_NumInteractions();
            nPUTrue  =  pu.getTrueNumInteractions(); 
	}
    }

    // --- Fill TTree vars
    nVtx_                    = (unsigned short) vertices->size();
    goodVtx_                 = vtx.isValid() ? true:false;

    runNum_                  = iEvent.id().run();
    lumi_                    = iEvent.luminosityBlock();
    eventNum_                = iEvent.id().event();

    nleps_                   = (unsigned short) myLeptons.size();
    for(int ii = 0 ; ii < nlepsMax; ii++) 
    {
       	lepPt_             [ii]  = ii < nleps_ ? myLeptons[ii]->pt()                                                   : 0;
	lepEta_            [ii]  = ii < nleps_ ? myLeptons[ii]->eta()                                                  : 0;
	lepPhi_            [ii]  = ii < nleps_ ? myLeptons[ii]->phi()                                                  : 0;
	lepM_              [ii]  = ii < nleps_ ? myLeptons[ii]->mass()                                                 : 0;
        lepID_             [ii]  = ii < nleps_ ? myLeptons[ii]->pdgId()                                                : 0;
        lepIso_            [ii]  = ii < nleps_ ? LeptonRelIso(myLeptons[ii])                                           : 0;  
        lepPtRel_          [ii]  = ii < nleps_ ? PtRel(myLeptons[ii], myJets)                                          : 1.e+5;  
        lepGenMatchIndex_  [ii]  = ii < nleps_ ? getMatchedIndex(myLeptons[ii], myGenLeptons)                          : -1; //function returns -1 if not matched  
        lepTriggerMatch_   [ii]  = ii < nleps_ ? triggerMatch(myLeptons[ii], hltEgammaCandidates, hltL3MuonCandidates) : 0; // match lepton either to egamma or to muon
    }
    
   
    l1l2DPhi_                = nleps_>=2 ? l1.DeltaPhi(l2) : 0;
    l1l2DR_                  = nleps_>=2 ? l1.DeltaR(l2) : 0;
    l1l2Pt_                  = nleps_>=2 ? (l1+l2).Pt() : 0;
    l1l2M_                   = nleps_>=2 ? (l1+l2).M() : 0;
    l1l2Eta_                 = nleps_>=2 && l1l2Pt_ > 1.e-3 ? (l1+l2).Eta() : 0;
    l1l2Phi_                 = nleps_>=2 ? (l1+l2).Phi() : 0;

    njets_                   = (unsigned short) myJets.size();
    for(int ii = 0 ; ii < njetsMax; ii++) 
    {
       	jetPt_      [ii]  = ii < njets_ ? myJets[ii]->pt()    : 0;
	jetEta_     [ii]  = ii < njets_ ? myJets[ii]->eta()   : 0;
	jetPhi_     [ii]  = ii < njets_ ? myJets[ii]->phi()   : 0;
	jetM_       [ii]  = ii < njets_ ? myJets[ii]->mass()  : 0;
        jetBTag_    [ii]  = 0;  // TBI 
    }

    nrjets_                   = (unsigned short) myRJets.size();
    for(int ii = 0 ; ii < nrjetsMax; ii++) 
    {
       	rjetPt_      [ii]  = ii < nrjets_ ? myRJets[ii]->pt()    : 0;
	rjetEta_     [ii]  = ii < nrjets_ ? myRJets[ii]->eta()   : 0;
	rjetPhi_     [ii]  = ii < nrjets_ ? myRJets[ii]->phi()   : 0;
	rjetM_       [ii]  = ii < nrjets_ ? myRJets[ii]->mass()  : 0;
        rjetBTag_    [ii]  = 0;  // TBI 
    }

    ngenleps_                   = (unsigned short) myGenLeptons.size();
    for(int ii = 0 ; ii < ngenlepsMax; ii++) 
    {
       	genlepPt_      [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->pt()                                                    : 0;
	genlepEta_     [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->eta()                                                   : 0;
	genlepPhi_     [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->phi()                                                   : 0;
	genlepM_       [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->mass()                                                  : 0;
        genlepID_      [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->pdgId()                                                 : 0;
        genlepMID_     [ii]  = ii < ngenleps_ ? (int)getGenMother(myGenLeptons[ii])->pdgId()                              : 0;
        genlepGMID_    [ii]  = ii < ngenleps_ ? (int)getGenMother(getGenMother(myGenLeptons[ii]))->pdgId()                : 0;
        genlepGGMID_   [ii]  = ii < ngenleps_ ? (int)getGenMother(getGenMother(getGenMother(myGenLeptons[ii])))->pdgId()  : 0;
    }

    ngenparts_                 = (unsigned short) myGenParticles.size();
    for(int ii = 0 ; ii < ngenpartsMax; ii++) 
    {
       	genpartPt_      [ii]  = ii < ngenparts_ ? myGenParticles[ii]->pt()                                               : 0;
	genpartEta_     [ii]  = ii < ngenparts_ ? myGenParticles[ii]->eta()                                              : 0;
	genpartPhi_     [ii]  = ii < ngenparts_ ? myGenParticles[ii]->phi()                                              : 0;
	genpartM_       [ii]  = ii < ngenparts_ ? myGenParticles[ii]->mass()                                             : 0;
        genpartID_      [ii]  = ii < ngenparts_ ? myGenParticles[ii]->pdgId()                                            : 0;
        genpartDID1_    [ii]  = ii < ngenparts_ ? myGenParticles[ii]->daughter(0)->pdgId()                               : 0;
        genpartDID2_    [ii]  = ii < ngenparts_ ? myGenParticles[ii]->daughter(1)->pdgId()                               : 0;

        genpartDRMI1_   [ii]  = -1; // default values
        genpartDRMI2_   [ii]  = -1; // default values

        if(abs(genpartDID1_[ii]) == 11 || abs(genpartDID1_[ii])==13) // try to match them to reco leptons
	{
     	    genpartDRMI1_   [ii]  = ii < ngenparts_ ?  getMatchedIndex(myGenParticles[ii]->daughter(0), myLeptons) : -1;
	}

        if(abs(genpartDID2_[ii]) == 11 || abs(genpartDID2_[ii])==13) // try to match them to reco leptons
	{
	    genpartDRMI2_   [ii]  = ii < ngenparts_ ?  getMatchedIndex(myGenParticles[ii]->daughter(1), myLeptons) : -1;
	}

        if(abs(genpartDID1_[ii]) >= 1 && abs(genpartDID1_[ii])<=5) // try to match them to reco jets
	{
     	    genpartDRMI1_   [ii]  = ii < ngenparts_ ?  getMatchedIndex(myGenParticles[ii]->daughter(0), myJets) : -1;
	}

        if(abs(genpartDID2_[ii]) >= 1 && abs(genpartDID2_[ii])<=5) // try to match them to reco jets
	{
	    genpartDRMI2_   [ii]  = ii < ngenparts_ ?  getMatchedIndex(myGenParticles[ii]->daughter(1), myJets) : -1;
        }
    }

    genl1l2DPhi_             = ngenleps_>=2 ? genl1.DeltaPhi(genl2) : 0;
    genl1l2DR_               = ngenleps_>=2 ? genl1.DeltaR(genl2) : 0;
    genl1l2Pt_               = ngenleps_>=2 ? (genl1+genl2).Pt() : 0;
    genl1l2M_                = ngenleps_>=2 ? (genl1+genl2).M() : 0;
    genl1l2Eta_              = ngenleps_>=2 && genl1l2Pt_ > 1.e-3 ? (genl1+genl2).Eta() : 0;
    genl1l2Phi_              = ngenleps_>=2 ? (genl1+genl2).Phi() : 0;

    nphos_                   = (unsigned short) myPhotons.size();
    
    met_                     = rawmet;          
    metPhi_                  = rawmetPhi;
    genmet_                  = genmet;          
    sumEt_                   = rawmetSumEt;
    t1met_                   = t1met;
    t1metPhi_                = t1metPhi;
    t1sumEt_                 = t1metSumEt;

    vHT_                     = HTVector.Pt();            
    t1vHT_                   = t1HTVector.Pt();         
    jvHT_                    = jHTVector.Pt();          
    
    rho_                     = (float) (*rhoH);

    nPU_                     = nPU; 
    nPUTrue_                 = nPUTrue; 

    HLT_e1e2_                = HLT_e1e2;
    HLT_mu1mu2_              = HLT_mu1mu2;
    HLT_mu1e2_               = HLT_mu1e2;
    HLT_e1mu2_               = HLT_e1mu2;
    HLT_pfmet_               = HLT_pfmet;
    HLT_pfmetCSV_            = HLT_pfmetCSV;
    
    // code below was "borrowed" from https://github.com/manuelfs/CfANtupler/blob/master/minicfa/interface/miniAdHocNTupler.h
    const edm::TriggerNames &fnames = iEvent.triggerNames(*filterBits);
    for (unsigned int i = 0, n = filterBits->size(); i < n; ++i) 
    {
      bool filterdecision(true);
      filterdecision = filterBits->accept(i);
      string filterName = fnames.triggerName(i);
      if(filterName=="Flag_trackingFailureFilter")                  Flag_trackingFailureFilter_ = filterdecision; 
      if(filterName=="Flag_goodVertices")			    Flag_goodVertices_ = filterdecision;       
      if(filterName=="Flag_CSCTightHaloFilter")                     Flag_CSCTightHaloFilter_ = filterdecision;                                                         		 
      if(filterName=="Flag_trkPOGFilters")                          Flag_trkPOGFilters_ = filterdecision;                                     			 
      if(filterName=="Flag_trkPOG_logErrorTooManyClusters")         Flag_trkPOG_logErrorTooManyClusters_ = filterdecision;	 
      if(filterName=="Flag_EcalDeadCellTriggerPrimitiveFilter")     Flag_EcalDeadCellTriggerPrimitiveFilter_ = filterdecision;
      if(filterName=="Flag_ecalLaserCorrFilter")                    Flag_ecalLaserCorrFilter_ = filterdecision;                                      		 
      if(filterName=="Flag_trkPOG_manystripclus53X")                Flag_trkPOG_manystripclus53X_ = filterdecision;                                    		 
      if(filterName=="Flag_eeBadScFilter")                          Flag_eeBadScFilter_ = filterdecision;                       			 
      if(filterName=="Flag_METFilters")                             Flag_METFilters_ = filterdecision;                    			
      if(filterName=="Flag_HBHENoiseFilter")                        Flag_HBHENoiseFilter_ = filterdecision;                 			 
      if(filterName=="Flag_trkPOG_toomanystripclus53X")             Flag_trkPOG_toomanystripclus53X_ = filterdecision;                        	 
      if(filterName=="Flag_hcalLaserEventFilter")	            Flag_hcalLaserEventFilter_ = filterdecision;                           	 
    }

    isDYTauTau_            = isDYTauTau_; // variable has been already initialized & previously set, shown here for completeness

    events_->Fill();
}


bool SimpleROOT::isGoodPhoton(const pat::Photon &photon)
{
    bool res = true; // by default is good, unless fails a cut bellow

    if(photon.pt() < 20) res = false; 
    if(fabs(photon.eta()) > 2.4) res = false;

    return res;
}


bool SimpleROOT::isGoodJet(const pat::Jet &myJet)
{
    bool res = true; // by default is good, unless fails a cut bellow

    if(myJet.pt() < 30) res = false; 
    if(fabs(myJet.eta()) > 2.4) res = false;

    return res;
}

bool SimpleROOT::isGoodMuon(const pat::Muon &mu)
{
    // --- https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
    bool res = true; // by default is good, unless fails a cut bellow

    if(mu.pt() < 10) res = false; 
    if(fabs(mu.eta()) > 2.4) res = false;
    if(!mu.isTightMuon(vtx)) res = false;
   
    // --- isolation --- those not used are commented out
    if(res && LeptonRelIso((reco::Candidate*)&mu) > 0.15) res = false;
//    if(getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&mu), 0.05, 0.2, 10., false, false) > 0.15) res = false; // miniISO

    return res;
}

bool SimpleROOT::triggerMatch(const reco::Candidate *myRecoLep, const std::vector<TLorentzVector> & egamma, const std::vector<TLorentzVector> & l3muon)
{
    bool res  = false;
    bool isEl = false;
    bool isMu = false;
    if(abs(myRecoLep->pdgId()) == 11)isEl = true;
    if(abs(myRecoLep->pdgId()) == 13)isMu = true;
   
    if(isEl)
    for(auto & triggerObj: egamma)
    {
	float DR = triggerObj.DeltaR(P4(myRecoLep));
        if(DR < 0.2) res = true;
    }

    if(isMu)
    for(auto & triggerObj: l3muon)
    {
	float DR = triggerObj.DeltaR(P4(myRecoLep));
        if(DR < 0.2) res = true;
    }


    return res;
}


short SimpleROOT::getMatchedIndex(const reco::Candidate * particleToBeMatched, vector<const reco::Candidate *> particleList)
{
    short myIndex = -1; // don't change this
    
    float DR = 1.e+9;

    TLorentzVector particleToBeMatchedP4 = P4(particleToBeMatched);
    int particleToBeMatchedID = particleToBeMatched->pdgId();

    for(unsigned int genIndex = 0; genIndex < particleList.size(); ++genIndex)
    {
	TLorentzVector myRecLepP4 = P4(particleList[genIndex]);
        float myDR = myRecLepP4.DeltaR(particleToBeMatchedP4);
        int particleID = particleList[genIndex]->pdgId();

        if(particleID == particleToBeMatchedID)        // considers matching only of correct charge + ID
	if(myDR < DR && myDR < 0.1) myIndex = genIndex;
    }
    return myIndex;
}


float SimpleROOT::MuonRelIso(const reco::Candidate *cand)
{
    float relIsoWithEA = 0.001;
    pat::Muon mu = *((pat::Muon*)cand);  
    // Effective areas from https://indico.cern.ch/event/367861/contribution/2/material/slides/0.pdf
    const int nEtaBins = 5; 
    const float etaBinLimits[nEtaBins+1] = {0.0, 0.8, 1.3, 2.0, 2.2, 2.5};
    const float effectiveAreaValues[nEtaBins] = {0.0913, 0.0765, 0.0546, 0.0728, 0.1177};

    reco::MuonPFIsolation pfIso =  mu.pfIsolationR03();
    // Find eta bin first. If eta>2.5, the last eta bin is used.
    int etaBin = 0;
    while(etaBin < nEtaBins-1 && abs(mu.eta()) > etaBinLimits[etaBin+1])  ++etaBin;

    float area = effectiveAreaValues[etaBin];
    relIsoWithEA = (float)( pfIso.sumChargedHadronPt + max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - (*rhoH) * area ) )/mu.pt();

    return relIsoWithEA;
}


float SimpleROOT::ElectronRelIso(const reco::Candidate *cand)
{
    float relIsoWithEA = 0;
    pat::Electron el = *((pat::Electron*)cand);  
    // Effective areas from https://indico.cern.ch/event/367861/contribution/2/material/slides/0.pdf
    const int nEtaBins = 5; 
    const float etaBinLimits[nEtaBins+1] = {0.0, 0.8, 1.3, 2.0, 2.2, 2.5};
    const float effectiveAreaValues[nEtaBins] = {0.1013, 0.0988, 0.0572, 0.0842, 0.1530};

    reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
    float etaSC = el.superCluster()->eta();

    // Find eta bin first. If eta>2.5, the last eta bin is used.
    int etaBin = 0;
    while(etaBin < nEtaBins-1 && abs(etaSC) > etaBinLimits[etaBin+1])  ++etaBin;

    float area = effectiveAreaValues[etaBin];
    relIsoWithEA = (float)( pfIso.sumChargedHadronPt + max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - (*rhoH) * area ) )/el.pt();

    return relIsoWithEA;
}

bool SimpleROOT::isGoodElectron(const pat::Electron &el)
{
    bool res = true; // by default is good, unless fails a cut bellow
    bool isEBEEGap = fabs(el.superCluster()->eta()) > 1.4442 && fabs(el.superCluster()->eta()) < 1.5660 ? 1 : 0;

    if(el.pt() < 10) res = false;
    if(fabs(el.eta()) > 2.4 && res == true) res = false;
    if(isEBEEGap && res==true) res=false;

    bool isEB      = fabs(el.superCluster()->eta()) < 1.4442 ? 1 : 0; 
    bool isEE      = fabs(el.superCluster()->eta()) > 1.5660 ? 1 : 0;
    if(res) 
    {
        // --- https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#PHYS14_selection_all_conditions  (v13)
        // --- https://indico.cern.ch/event/292938/  {Loose ID w/o isolation | Phys14}
        // --- http://cmslxr.fnal.gov/source/DataFormats/EgammaCandidates/interface/GsfElectron.h?v=CMSSW_7_3_2
        float  trackMomentumAtVtx               = (float)sqrt(el.trackMomentumAtVtx().mag2());
        float  ecalEnergy                       = (float)el.ecalEnergy();
    
        float  full5x5_sigmaIetaIeta            = (float)el.full5x5_sigmaIetaIeta();
        float  dEtaIn                           = (float)el.deltaEtaSuperClusterTrackAtVtx();
        float  dPhiIn                           = (float)el.deltaPhiSuperClusterTrackAtVtx();
        float  HoE                              = (float)el.hadronicOverEm();
        float  ooEmooP                          = (float)fabs(1/ecalEnergy - 1/trackMomentumAtVtx);
        float  d0                               = (float)el.gsfTrack()->dxy(vtx.position());
        float  dz                               = (float)el.gsfTrack()->dz(vtx.position());
        int    expectedMissingInnerHits         = el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
        bool   passConversionVeto               = el.passConversionVeto(); 
    
        if(isEB)
        {
    	    if(res && full5x5_sigmaIetaIeta         >  0.010557)res=false;    
	    if(res && fabs(dEtaIn)                  >  0.012442)res=false;                   
            if(res && fabs(dPhiIn)                  >  0.072624)res=false;                   
            if(res && HoE                           >  0.121476)res=false; 
            if(res && ooEmooP                       >  0.221803)res=false; 
            if(res && fabs(d0)                      >  0.022664)res=false; 
            if(res && fabs(dz)                      >  0.173670)res=false; 
            if(res && expectedMissingInnerHits      >= 2       )res=false;
            if(res && passConversionVeto            == false   )res=false;
	}

        if(isEE)
        {
    	    if(res && full5x5_sigmaIetaIeta         >  0.032602)res=false;    
	    if(res && fabs(dEtaIn)                  >  0.010654)res=false;                   
            if(res && fabs(dPhiIn)                  >  0.145129)res=false;                   
            if(res && HoE                           >  0.131862)res=false; 
            if(res && ooEmooP                       >  0.142283)res=false; 
            if(res && fabs(d0)                      >  0.097358)res=false; 
            if(res && fabs(dz)                      >  0.198444)res=false; 
            if(res && expectedMissingInnerHits      >= 2       )res=false;
            if(res && passConversionVeto            == false   )res=false;
        }
    }


    // --- isolation -- those not used are commented out
    if(res && LeptonRelIso((reco::Candidate*)&el) > 0.15)res = false;
//    if(res && getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&el), 0.05, 0.2, 10., false, false) > 0.15)res = false;  //miniISO

    return res;
}

float SimpleROOT::PtRel(const reco::Candidate *myLepton, vector<const reco::Candidate *> myJets)
{
    float myPtRel = 1.e+5;
    for(auto & myJet: myJets)
    {
	float ptrel = P4(myLepton).Perp(P4(myJet).Vect());
	myPtRel = ptrel < myPtRel ? ptrel : myPtRel;
    }

    if(myJets.size() == 0) myPtRel = 0;

    return myPtRel;
}

bool SimpleROOT::isGoodVertex(const reco::Vertex &vtx)
{
  // The "good vertex" selection is borrowed Ilya Kravchenko who borrowed from Giovanni Zevi Della Porta
  return ( !(vtx.chi2()==0 && vtx.ndof()==0) && vtx.ndof()>=4. && vtx.position().Rho()<=2.0 && fabs(vtx.position().Z())<=24.0 ) ? true : false; 
}

const reco::Candidate *SimpleROOT::getGenMother(const reco::Candidate* particle) // bet you can't write this in more combact format ;-)
{
  if(particle->numberOfMothers() == 0) return particle; // there are no mothers you are sitting on the proton
  return particle->pdgId() != particle->mother(0)->pdgId() ? particle->mother(0): getGenMother(particle->mother(0));
}


//miniISO from https://github.com/manuelfs/CfANtupler/blob/master/minicfa/interface/miniAdHocNTupler.h
float SimpleROOT::getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
                        const reco::Candidate* ptcl,  
                        float r_iso_min, float r_iso_max, float kt_scale,
                        bool use_pfweight, bool charged_only) {

    if (ptcl->pt()<5.) return 99999.;

    float deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
      if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl->isMuon()) {
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;  
    } else {
      //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    float iso_nh(0.); float iso_ch(0.); 
    float iso_ph(0.); float iso_pu(0.);
    float ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    float r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));
    for (const pat::PackedCandidate &pfc : *pfcands) {
      if (abs(pfc.pdgId())<7) continue;

      float dr = deltaR(pfc, *ptcl);
      if (dr > r_iso) continue;
      
      //////////////////  NEUTRALS  /////////////////////////
      if (pfc.charge()==0){
        if (pfc.pt()>ptThresh) {
          float wpf(1.);
          if (use_pfweight){
            float wpv(0.), wpu(0.);
            for (const pat::PackedCandidate &jpfc : *pfcands) {
              float jdr = deltaR(pfc, jpfc);
              if (pfc.charge()!=0 || jdr<0.00001) continue;
              float jpt = jpfc.pt();
              if (pfc.fromPV()>1) wpv *= jpt/jdr;
              else wpu *= jpt/jdr;
            }
            wpv = log(wpv);
            wpu = log(wpu);
            wpf = wpv/(wpv+wpu);
          }
          /////////// PHOTONS ////////////
          if (abs(pfc.pdgId())==22) {
            if(dr < deadcone_ph) continue;
            iso_ph += wpf*pfc.pt();
	    /////////// NEUTRAL HADRONS ////////////
          } else if (abs(pfc.pdgId())==130) {
            if(dr < deadcone_nh) continue;
            iso_nh += wpf*pfc.pt();
          }
        }
        //////////////////  CHARGED from PV  /////////////////////////
      } else if (pfc.fromPV()>1){
        if (abs(pfc.pdgId())==211) {
          if(dr < deadcone_ch) continue;
          iso_ch += pfc.pt();
        }
        //////////////////  CHARGED from PU  /////////////////////////
      } else {
        if (pfc.pt()>ptThresh){
          if(dr < deadcone_pu) continue;
          iso_pu += pfc.pt();
        }
      }
    }
    float iso(0.);
    if (charged_only){
      iso = iso_ch;
    } else {
      iso = iso_ph + iso_nh;
      if (!use_pfweight) iso -= 0.5*iso_pu;
      if (iso>0) iso += iso_ch;
      else iso = iso_ch;
    }
    iso = iso/ptcl->pt();

    return iso;
  }

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleROOT);
