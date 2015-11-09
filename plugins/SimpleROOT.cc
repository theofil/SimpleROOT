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
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"
#include "TLorentzVector.h"


#define njetsMax 30 
#define njetsFWMax 30 
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
  
        short getMatchedIndex(const reco::Candidate *, vector<const reco::Candidate *>);
        bool  triggerMatch(const reco::Candidate *, const std::vector<TLorentzVector> &, const std::vector<TLorentzVector> &);
        const reco::Candidate *getGenMother(const reco::Candidate*); 
       
        float getPFIsolation(edm::Handle<pat::PackedCandidateCollection>, const reco::Candidate*, float, float, float, bool, bool); // miniISO 

        TLorentzVector P4(const reco::Candidate* cand){TLorentzVector p4vec; p4vec.SetPxPyPzE( cand->px(), cand->py(), cand->pz(), cand->energy() ); return p4vec;}

    	edm::Handle<reco::VertexCollection> vertices;
   	edm::Handle<pat::MuonCollection> muons;
     	edm::Handle<pat::ElectronCollection> electrons;
        edm::Handle<pat::JetCollection> jets;
        edm::Handle<pat::METCollection> puppimets;
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
        edm::Handle<LHEEventProduct> LHEEventInfo;
        edm::Handle<GenEventInfoProduct> genEvtInfo;
        edm::Handle<edm::View<reco::GenJet>> genjets;

        
    	edm::EDGetTokenT<reco::VertexCollection> verticesToken;
   	edm::EDGetTokenT<pat::MuonCollection> muonsToken;
     	edm::EDGetTokenT<pat::ElectronCollection> electronsToken;
        edm::EDGetTokenT<pat::JetCollection> jetsToken;
        edm::EDGetTokenT<pat::METCollection> puppimetsToken;
        edm::EDGetTokenT<pat::METCollection> metsToken;
        edm::EDGetTokenT<pat::PhotonCollection> photonsToken;
        edm::EDGetTokenT<double> rhoHToken;
        edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken;
        edm::EDGetTokenT<edm::View<reco::GenParticle>> prunedToken;
        edm::EDGetTokenT<edm::View<pat::PackedGenParticle>> packedToken;
        edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken;
        edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken;
        edm::EDGetTokenT<edm::TriggerResults> filterBitsToken;
        edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken;
        edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken;
        edm::EDGetTokenT<LHEEventProduct> LHEEventProductToken;
        edm::EDGetTokenT<GenEventInfoProduct> GenEventInfoProductToken;
        edm::EDGetTokenT<edm::View<reco::GenJet>> genjetsToken;

        reco::Vertex vtx; // stores event's primary vertex

	edm::Service<TFileService> fileService_; 
	TTree *events_;
	
        // --- all ntuple vars end with "_" --- 
	bool goodVtx_; 
	unsigned short nVtx_; 

        ULong64_t eventNum_;
        unsigned int runNum_;
        unsigned int lumi_;
        float genWeight_;

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
        short lepID_                      [nlepsMax];
        short lepGenMatchIndex_           [nlepsMax];   // index in the stored genLep[] collection
        bool  lepTriggerMatch_            [nlepsMax];

        unsigned short njets_;
        float jetPt_                      [njetsMax]; 
        float jetEta_                     [njetsMax]; 
        float jetPhi_                     [njetsMax]; 
        float jetM_                       [njetsMax]; 
        float jetBTag_                    [njetsMax];
        float jetGenPt_                   [njetsMax];

        unsigned short njetsFW_;
        float jetFWPt_                    [njetsFWMax]; 
        float jetFWEta_                   [njetsFWMax]; 
        float jetFWPhi_                   [njetsFWMax]; 
        float jetFWM_                     [njetsFWMax]; 
        float jetFWGenPt_                 [njetsFWMax];

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
        bool  genlepIsPrompt_             [ngenlepsMax]; // rivet safe prompt flag

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

        float puppimet_;          
        float met_;          
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

	unsigned short HLT_e1e2_;   // 0 if not fired, otherwise store the prescale (should be 1 for unprescaled paths) 
	unsigned short HLT_mu1mu2_; 
	unsigned short HLT_mu1e2_; 
	unsigned short HLT_e1mu2_; 
	unsigned short HLT_pfmet_; 
	unsigned short HLT_e1_; 
	unsigned short HLT_mu1_; 
	unsigned short HLT_ZeroBias_; 
       
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
        bool isData_;

        unsigned short nphos_;
        TH1I *mcWeights_;
};

SimpleROOT::SimpleROOT(const edm::ParameterSet& iConfig):  // initialize tokens in the constructor, better performance
verticesToken(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter("offlineSlimmedPrimaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices")))),
muonsToken(consumes<pat::MuonCollection>(iConfig.getUntrackedParameter("slimmedMuons",edm::InputTag("slimmedMuons")))),
electronsToken(consumes<pat::ElectronCollection>(iConfig.getUntrackedParameter("slimmedElectrons",edm::InputTag("slimmedElectrons")))),
jetsToken(consumes<pat::JetCollection>(iConfig.getUntrackedParameter("slimmedJets",edm::InputTag("slimmedJets")))),
puppimetsToken(consumes<pat::METCollection>(iConfig.getUntrackedParameter("slimmedMETsPuppi",edm::InputTag("slimmedMETsPuppi")))),
metsToken(consumes<pat::METCollection>(iConfig.getUntrackedParameter("slimmedMETs",edm::InputTag("slimmedMETs")))),
photonsToken(consumes<pat::PhotonCollection>(iConfig.getUntrackedParameter("slimmedPhotons",edm::InputTag("slimmedPhotons")))),
rhoHToken(consumes<double>(iConfig.getUntrackedParameter("fixedGridRhoFastjetAll",edm::InputTag("fixedGridRhoFastjetAll")))),
pileupToken(consumes<edm::View<PileupSummaryInfo> >(iConfig.getUntrackedParameter("slimmedAddPileupInfo", edm::InputTag("slimmedAddPileupInfo")))),  
prunedToken(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter("prunedGenParticles", edm::InputTag("prunedGenParticles")))),
packedToken(consumes<edm::View<pat::PackedGenParticle> >(edm::InputTag("packedGenParticles"))),
pfcandsToken(consumes<pat::PackedCandidateCollection> (iConfig.getUntrackedParameter("packedPFCandidates", edm::InputTag("packedPFCandidates")))),
triggerBitsToken(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"))),
filterBitsToken(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","PAT"))),
triggerObjectsToken(consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("selectedPatTrigger"))),
triggerPrescalesToken(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),
LHEEventProductToken(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"))),
GenEventInfoProductToken(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
genjetsToken(consumes<edm::View<reco::GenJet>>(iConfig.getUntrackedParameter("slimmedGenJets",edm::InputTag("slimmedGenJets"))))
//Token(consumes<>(iConfig.getUntrackedParameter("",edm::InputTag("")))),  // first arg is default the second is used only if is defined in runme_cfg.py
{
    mcWeights_ = fileService_->make<TH1I>("mcWeights","mcWeights", 4, -2, 2);

    events_ = fileService_->make<TTree>("events","events");
    events_->Branch("goodVtx"          ,&goodVtx_               ,"goodVtx/O");
    events_->Branch("nVtx"             ,&nVtx_                  ,"nVtx/s");
    events_->Branch("eventNum"         ,&eventNum_              ,"eventNum/l");
    events_->Branch("runNum"           ,&runNum_                ,"runNum/i");
    events_->Branch("lumi"             ,&lumi_                  ,"lumi/i");
    events_->Branch("genWeight"        ,&genWeight_             ,"genWeight/F");

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
    events_->Branch("lepID"            ,lepID_                  ,"lepID[nleps]/S");
    events_->Branch("lepGenMatchIndex" ,lepGenMatchIndex_       ,"lepGenMatchIndex[nleps]/S");
    events_->Branch("lepTriggerMatch"  ,lepTriggerMatch_        ,"lepTriggerMatch[nleps]/O");

    events_->Branch("njets"            ,&njets_                 ,"njets/s");
    events_->Branch("jetPt"            ,jetPt_                  ,"jetPt[njets]/F");
    events_->Branch("jetEta"           ,jetEta_                 ,"jetEta[njets]/F");
    events_->Branch("jetPhi"           ,jetPhi_                 ,"jetPhi[njets]/F");
    events_->Branch("jetM"             ,jetM_                   ,"jetM[njets]/F");
    events_->Branch("jetBTag"          ,jetBTag_                ,"jetBTag[njets]/F");
    events_->Branch("jetGenPt"         ,jetGenPt_               ,"jetGenPt[njets]/F");

    events_->Branch("njetsFW"          ,&njetsFW_               ,"njetsFW/s");
    events_->Branch("jetFWPt"          ,jetFWPt_                ,"jetFWPt[njetsFW]/F");
    events_->Branch("jetFWEta"         ,jetFWEta_               ,"jetFWEta[njetsFW]/F");
    events_->Branch("jetFWPhi"         ,jetFWPhi_               ,"jetFWPhi[njetsFW]/F");
    events_->Branch("jetFWM"           ,jetFWM_                 ,"jetFWM[njetsFW]/F");
    events_->Branch("jetFWGenPt"       ,jetFWGenPt_             ,"jetFWGenPt[njetsFW]/F");

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
    events_->Branch("genlepIsPrompt"   ,genlepIsPrompt_         ,"genlepIsPrompt[ngenleps]/O");

    events_->Branch("genl1l2M"         ,&genl1l2M_              ,"genl1l2M/F");
    events_->Branch("genl1l2Pt"        ,&genl1l2Pt_             ,"genl1l2Pt/F");
    events_->Branch("genl1l2Eta"       ,&genl1l2Eta_            ,"genl1l2Eta/F");
    events_->Branch("genl1l2Phi"       ,&genl1l2Phi_            ,"genl1l2Phi/F");
    events_->Branch("genl1l2DPhi"      ,&genl1l2DPhi_           ,"genl1l2DPhi/F");
    events_->Branch("genl1l2DR"        ,&genl1l2DR_             ,"genl1l2DR/F");

    events_->Branch("nphos"            ,&nphos_                 ,"nphos/s");
    events_->Branch("puppimet"         ,&puppimet_              ,"puppimet/F");
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

    events_->Branch("HLT_e1e2"         ,&HLT_e1e2_              ,"HLT_e1e2/s");
    events_->Branch("HLT_mu1mu2"       ,&HLT_mu1mu2_            ,"HLT_mu1mu2/s");
    events_->Branch("HLT_mu1e2"        ,&HLT_mu1e2_             ,"HLT_mu1e2/s");
    events_->Branch("HLT_e1mu2"        ,&HLT_e1mu2_             ,"HLT_e1mu2/s");
    events_->Branch("HLT_pfmet"        ,&HLT_pfmet_             ,"HLT_pfmet/s");
    events_->Branch("HLT_e1"           ,&HLT_e1_                ,"HLT_e1/s");
    events_->Branch("HLT_mu1"          ,&HLT_mu1_               ,"HLT_mu1/s");
    events_->Branch("HLT_ZeroBias"     ,&HLT_ZeroBias_          ,"HLT_ZeroBias/s");

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
    events_->Branch("isData"             ,&isData_                  ,"isData/O");

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
    isData_ = iEvent.isRealData();

    iEvent.getByToken(verticesToken, vertices);  
    iEvent.getByToken(muonsToken, muons);
    iEvent.getByToken(electronsToken, electrons);
    iEvent.getByToken(jetsToken, jets);
    iEvent.getByToken(puppimetsToken, puppimets);
    iEvent.getByToken(metsToken, mets);
    iEvent.getByToken(photonsToken, photons);
    iEvent.getByToken(rhoHToken, rhoH);
    iEvent.getByToken(triggerBitsToken, triggerBits);
    iEvent.getByToken(filterBitsToken, filterBits);
    iEvent.getByToken(triggerObjectsToken, triggerObjects);
    iEvent.getByToken(triggerPrescalesToken, triggerPrescales);

    if(!isData_)  // available only for MC
    {
    	iEvent.getByToken(prunedToken, pruned);
        iEvent.getByToken(GenEventInfoProductToken, genEvtInfo);
    	iEvent.getByToken(genjetsToken, genjets);
        iEvent.getByToken(pileupToken, pileup);
//   iEvent.getByToken(LHEEventProductToken, LHEEventInfo); iEvent.getByToken(packedToken, packed);    iEvent.getByToken(pfcandsToken, pfcands);
    }

    vector<const reco::Candidate *> myLeptons; // in this container we will store all selected RECO electrons and RECO muons 
    vector<const reco::Candidate *> myJets; // in this container we will store all prompt jets (PV)
    vector<const reco::Candidate *> myRJets; // in this container we will all rejected prompt jets (PV) -- due to DR matching with myLeptons
    vector<const reco::Candidate *> myJetsFW; // in this container we will store all prompt jets (PV)
    vector<const reco::Candidate *> myPhotons; // in this container we will store all photons

    vector<const reco::Candidate *> myGenLeptons;
    vector<const reco::Candidate *> myGenParticles; // store here interesting particles, like top, W, Z, higgs
    vector<const reco::Candidate *> myGenJets; // in this container we will store all prompt jets (PV)

    vector<TLorentzVector> hltEgammaCandidates;
    vector<TLorentzVector >hltL3MuonCandidates;

    // --- trigger info
    unsigned short HLT_e1e2     = 0;
    unsigned short HLT_mu1mu2   = 0;
    unsigned short HLT_mu1e2    = 0;
    unsigned short HLT_e1mu2    = 0;
    unsigned short HLT_pfmet    = 0;
    unsigned short HLT_e1       = 0;
    unsigned short HLT_mu1      = 0;
    unsigned short HLT_ZeroBias = 0;

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) 
    {
	string trigger_name = string(names.triggerName(i));
        trigger_name.pop_back();

        if(trigger_name == string("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")          && triggerBits->accept(i)) HLT_e1e2     = triggerPrescales->getPrescaleForIndex(i); 
        if(trigger_name == string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")                && triggerBits->accept(i)) HLT_mu1mu2   = triggerPrescales->getPrescaleForIndex(i);
        if(trigger_name == string("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")    && triggerBits->accept(i)) HLT_mu1e2    = triggerPrescales->getPrescaleForIndex(i);  
        if(trigger_name == string("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v")     && triggerBits->accept(i)) HLT_e1mu2    = triggerPrescales->getPrescaleForIndex(i); 
        if(trigger_name == string("HLT_PFMET170_NoiseCleaned_v")                          && triggerBits->accept(i)) HLT_pfmet    = triggerPrescales->getPrescaleForIndex(i); 
        if(trigger_name == string("HLT_Ele22_eta2p1_WPLoose_Gsf_v")                       && triggerBits->accept(i)) HLT_e1       = triggerPrescales->getPrescaleForIndex(i); 
        if(trigger_name == string("HLT_IsoMu17_eta2p1_v")                                 && triggerBits->accept(i)) HLT_mu1 = triggerPrescales->getPrescaleForIndex(i); 
        if(trigger_name == string("HLT_ZeroBias_v")                                       && triggerBits->accept(i)) HLT_ZeroBias = triggerPrescales->getPrescaleForIndex(i); 
    }


   // --- store HLT trigger objects for offline matching
   if(HLT_e1e2 !=0 || HLT_mu1mu2 !=0 || HLT_mu1e2 !=0 || HLT_e1mu2 !=0)
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
        bool isCE = true;

        if(myjet.pt() < 30) continue;  // considers only jets above 30

        if(fabs(myjet.eta()) > 2.4) isCE = false;

        for(auto & lep: myLeptons) if( P4(lep).DeltaR( P4(&myjet) ) < DRmax ) isLeptonMatched = true;

        // save central jets 
	if( isGoodJet(myjet) && !isLeptonMatched && isCE) myJets.push_back(&myjet);
	if( isGoodJet(myjet) && isLeptonMatched  && isCE ) myRJets.push_back(&myjet);

        // save all forward jets
	if( isGoodJet(myjet) && !isCE) myJetsFW.push_back(&myjet);
    }

    if(!isData_)
    for(const reco::GenJet &myjet : *genjets)  // fill-up all gen jets no cut (pdgId = 0 for genjets)
    {
	myGenJets.push_back(&myjet);
    }


    for (const pat::Photon &photon : *photons) 
    {
	if( isGoodPhoton(photon) ) myPhotons.push_back(&photon); 
    }

    isDYTauTau_ = false;

    if(!isData_)
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
    const pat::MET &puppimet = puppimets->front();
    float rawpuppimet        = puppimet.uncorPt();

    const pat::MET &met = mets->front();
    float rawmet      = met.uncorPt();
    float rawmetPhi   = met.uncorPhi();
    float genmet      = !isData_ ? met.genMET()->pt():0;
    float rawmetSumEt = met.shiftedSumEt(pat::MET::NoShift, pat::MET::Raw);
    float t1met       = met.pt();
    float t1metPhi    = met.phi();
    float t1metSumEt  = met.sumEt();

    //  --- sort by pt all objects, the [] is a C++11 lambda func 
    std::sort(myLeptons.begin(), myLeptons.end(), [](const reco::Candidate * a, const reco::Candidate * b){return a->pt() > b->pt();} );
    std::sort(myJets.begin(), myJets.end(), [](const reco::Candidate * a, const reco::Candidate * b){return a->pt() > b->pt();} );
    std::sort(myJetsFW.begin(), myJetsFW.end(), [](const reco::Candidate * a, const reco::Candidate * b){return a->pt() > b->pt();} );
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

    if(!isData_)
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
    isData_                  = isData_; // has been filled upstream
    genWeight_               = !isData_ ? genEvtInfo->weight() : 0 ;

    if(genWeight_<0) mcWeights_->Fill(-0.999);
    if(genWeight_>0) mcWeights_->Fill(+0.999);

    nleps_                   = (unsigned short) myLeptons.size();
    for(int ii = 0 ; ii < nlepsMax; ii++) 
    {
       	lepPt_             [ii]  = ii < nleps_ ? myLeptons[ii]->pt()                                                   : 0;
	lepEta_            [ii]  = ii < nleps_ ? myLeptons[ii]->eta()                                                  : 0;
	lepPhi_            [ii]  = ii < nleps_ ? myLeptons[ii]->phi()                                                  : 0;
	lepM_              [ii]  = ii < nleps_ ? myLeptons[ii]->mass()                                                 : 0;
        lepID_             [ii]  = ii < nleps_ ? myLeptons[ii]->pdgId()                                                : 0;
        lepIso_            [ii]  = ii < nleps_ ? LeptonRelIso(myLeptons[ii])                                           : 0;  
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
        int genjetMatchedIndex = ii < njets_ ? getMatchedIndex(myJets[ii], myGenJets): -1;

       	jetPt_      [ii]  = ii < njets_ ? myJets[ii]->pt()    : 0;
	jetEta_     [ii]  = ii < njets_ ? myJets[ii]->eta()   : 0;
	jetPhi_     [ii]  = ii < njets_ ? myJets[ii]->phi()   : 0;
	jetM_       [ii]  = ii < njets_ ? myJets[ii]->mass()  : 0;
        jetBTag_    [ii]  = ii < njets_ ? ((pat::Jet*) myJets[ii])->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : 0;  //  0.814 nominal cut
        jetGenPt_   [ii]  = ii < njets_ && genjetMatchedIndex >=0 ? myGenJets[genjetMatchedIndex]->pt() : 0;  
    }

    njetsFW_                   = (unsigned short) myJetsFW.size();
    for(int ii = 0 ; ii < njetsFWMax; ii++) 
    {
        int genjetMatchedIndex = ii < njetsFW_ ? getMatchedIndex(myJetsFW[ii], myGenJets): -1;

       	jetFWPt_      [ii]  = ii < njetsFW_ ? myJetsFW[ii]->pt()    : 0;
	jetFWEta_     [ii]  = ii < njetsFW_ ? myJetsFW[ii]->eta()   : 0;
	jetFWPhi_     [ii]  = ii < njetsFW_ ? myJetsFW[ii]->phi()   : 0;
	jetFWM_       [ii]  = ii < njetsFW_ ? myJetsFW[ii]->mass()  : 0;
        jetFWGenPt_   [ii]  = ii < njetsFW_ && genjetMatchedIndex >=0 ? myGenJets[genjetMatchedIndex]->pt() : 0;  
    }

    nrjets_                   = (unsigned short) myRJets.size();
    for(int ii = 0 ; ii < nrjetsMax; ii++) 
    {
       	rjetPt_      [ii]  = ii < nrjets_ ? myRJets[ii]->pt()    : 0;
	rjetEta_     [ii]  = ii < nrjets_ ? myRJets[ii]->eta()   : 0;
	rjetPhi_     [ii]  = ii < nrjets_ ? myRJets[ii]->phi()   : 0;
	rjetM_       [ii]  = ii < nrjets_ ? myRJets[ii]->mass()  : 0;
        rjetBTag_    [ii]  = ii < nrjets_ ? ((pat::Jet*) myRJets[ii])->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")  : 0;  // 0.814 nominal cut
    }

    ngenleps_                   = (unsigned short) myGenLeptons.size();
    for(int ii = 0 ; ii < ngenlepsMax; ii++) 
    {
       	genlepPt_        [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->pt()                                                    : 0;
	genlepEta_       [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->eta()                                                   : 0;
	genlepPhi_       [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->phi()                                                   : 0;
	genlepM_         [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->mass()                                                  : 0;
        genlepID_        [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->pdgId()                                                 : 0;
        genlepMID_       [ii]  = ii < ngenleps_ ? (int)getGenMother(myGenLeptons[ii])->pdgId()                              : 0;
        genlepGMID_      [ii]  = ii < ngenleps_ ? (int)getGenMother(getGenMother(myGenLeptons[ii]))->pdgId()                : 0;
        genlepGGMID_     [ii]  = ii < ngenleps_ ? (int)getGenMother(getGenMother(getGenMother(myGenLeptons[ii])))->pdgId()  : 0;
        genlepIsPrompt_  [ii]  = ii < ngenleps_ ? ((reco::GenParticle*)myGenLeptons[ii])->isPromptFinalState() 	            : 0;
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

        if(abs(genpartDID1_[ii]) >= 1 && abs(genpartDID1_[ii])<=5) // try to match them to reco central jets
	{
     	    genpartDRMI1_   [ii]  = ii < ngenparts_ ?  getMatchedIndex(myGenParticles[ii]->daughter(0), myJets) : -1;
	}

        if(abs(genpartDID2_[ii]) >= 1 && abs(genpartDID2_[ii])<=5) // try to match them to reco central jets
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
    
    puppimet_               =  rawpuppimet;          
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
    HLT_e1_                  = HLT_e1;
    HLT_mu1_                 = HLT_mu1;
    HLT_ZeroBias_            = HLT_ZeroBias;
    
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
    return true; // no jet-id for the moment
}

bool SimpleROOT::isGoodMuon(const pat::Muon &mu)
{
    // use https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon
    bool res = true; // by default is good, unless fails a cut bellow

    if(mu.pt() < 10) res = false; 
    if(fabs(mu.eta()) > 2.4) res = false;
    if(!mu.isTightMuon(vtx)) res = false;

    // --- isolation --- those not used are commented out
    if(res && LeptonRelIso((reco::Candidate*)&mu) > 0.15) res = false;

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
    float relIso = 0.001;
    pat::Muon *mu = ((pat::Muon*)cand);  

    relIso = (mu->pfIsolationR04().sumChargedHadronPt + max(0., mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5*mu->pfIsolationR04().sumPUPt))/mu->pt();
    return relIso;
}


float SimpleROOT::ElectronRelIso(const reco::Candidate *cand)
{
    float relIsoWithEA = 0;
    pat::Electron el = *((pat::Electron*)cand);  
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
    // https://indico.cern.ch/event/370507/contribution/1/attachments/1140657/1633761/Rami_eleCB_ID_25ns.pdf
    // Effective areas from https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
    const int nEtaBins = 7; 
    const float etaBinLimits[nEtaBins+1] = {0.0, 1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
    const float effectiveAreaValues[nEtaBins] = {0.1752 , 0.1862, 0.1411, 0.1534 , 0.1903, 0.2243, 0.2687};

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
        // --- using the EGM loose id for Spring 15 (25 ns): https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
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
    	    if(res && full5x5_sigmaIetaIeta               >  0.0103 )res=false;    
	    if(res && fabs(dEtaIn)                        >  0.0105 )res=false;                   
            if(res && fabs(dPhiIn)                        >  0.115  )res=false;                   
            if(res && HoE                                 >  0.104  )res=false; 
            if(res && LeptonRelIso((reco::Candidate*)&el) >  0.0893 )res=false;
            if(res && ooEmooP                             >  0.102  )res=false; 
            if(res && fabs(d0)                            >  0.0261 )res=false; 
            if(res && fabs(dz)                            >  0.41   )res=false; 
            if(res && expectedMissingInnerHits            >  2      )res=false;
            if(res && passConversionVeto                  == false  )res=false;
	}

        if(isEE)
        {
    	    if(res && full5x5_sigmaIetaIeta               >  0.0301  )res=false;    
	    if(res && fabs(dEtaIn)                        >  0.00814 )res=false;                   
            if(res && fabs(dPhiIn)                        >  0.182   )res=false;                   
            if(res && HoE                                 >  0.0897  )res=false; 
            if(res && LeptonRelIso((reco::Candidate*)&el) >  0.121   )res=false;
            if(res && ooEmooP                             >  0.126   )res=false; 
            if(res && fabs(d0)                            >  0.118   )res=false; 
            if(res && fabs(dz)                            >  0.822   )res=false; 
            if(res && expectedMissingInnerHits            >  1       )res=false;
            if(res && passConversionVeto                  == false   )res=false;
        }
    }

    // --- isolation -- commented out cause is part of Spring15 WP
    //if(res && LeptonRelIso((reco::Candidate*)&el) > 0.15)res = false;

    return res;
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


//define this as a plug-in
DEFINE_FWK_MODULE(SimpleROOT);
