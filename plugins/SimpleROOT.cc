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

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"
#include "TLorentzVector.h"


#define njetsMax 30 
#define nrjetsMax 10 
#define nlepsMax 10
#define ngenlepsMax 10

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
        short getLepGenMatchIndex(const reco::Candidate * myRecoLep, vector<const reco::Candidate *> myGenLeptons);

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

        
    	edm::EDGetTokenT<reco::VertexCollection> verticesToken;
   	edm::EDGetTokenT<pat::MuonCollection> muonsToken;
     	edm::EDGetTokenT<pat::ElectronCollection> electronsToken;
        edm::EDGetTokenT<pat::JetCollection> jetsToken;
        edm::EDGetTokenT<pat::METCollection> metsToken;
        edm::EDGetTokenT<pat::PhotonCollection> photonsToken;
        edm::EDGetTokenT<double> rhoHToken;
        edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken;
        edm::EDGetTokenT<edm::View<reco::GenParticle>> prunedToken;

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

        float genl1l2DPhi_;
        float genl1l2DR_;
        float genl1l2Pt_;
        float genl1l2M_;
        float genl1l2Eta_;
        float genl1l2Phi_;

        float met_;          // raw pf-met
        float metPhi_;
        float t1met_;
        float t1metPhi_;
        float sumEt_;
        float t1sumEt_;

        float rho_;
        float nPU_;
        float nPUTrue_;

        float vHT_;           // pt of the recoil vector of all objects excluding the dilepton [= -met - l1l2] 
        float t1vHT_;         // as vHT but with t1 correction
	float jvHT_;          // recoil of hard jets and subleading leptons


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
prunedToken(consumes<edm::View<reco::GenParticle> >(iConfig.getUntrackedParameter("prunedGenParticles", edm::InputTag("prunedGenParticles"))))
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

    events_->Branch("genl1l2M"         ,&genl1l2M_              ,"genl1l2M/F");
    events_->Branch("genl1l2Pt"        ,&genl1l2Pt_             ,"genl1l2Pt/F");
    events_->Branch("genl1l2Eta"       ,&genl1l2Eta_            ,"genl1l2Eta/F");
    events_->Branch("genl1l2Phi"       ,&genl1l2Phi_            ,"genl1l2Phi/F");
    events_->Branch("genl1l2DPhi"      ,&genl1l2DPhi_           ,"genl1l2DPhi/F");
    events_->Branch("genl1l2DR"        ,&genl1l2DR_             ,"genl1l2DR/F");

    events_->Branch("nphos"            ,&nphos_                 ,"nphos/s");
    events_->Branch("met"              ,&met_                   ,"met/F");
    events_->Branch("metPhi"           ,&metPhi_                ,"metPhi/F");
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

    vector<const reco::Candidate *> myLeptons; // in this container we will store all selected RECO electrons and RECO muons 
    vector<const reco::Candidate *> myJets; // in this container we will store all prompt jets (PV)
    vector<const reco::Candidate *> myRJets; // in this container we will all rejected prompt jets (PV) -- due to DR matching with myLeptons
    vector<const reco::Candidate *> myPhotons; // in this container we will store all photons

    vector<const reco::Candidate *> myGenLeptons;

    // --- loop inside the objects
    for (const reco::Vertex &PV : *vertices)
    {
	if( isGoodVertex(PV) ){ vtx = PV; break;} // get first good vertex 
    }

    for (const pat::Muon &mu : *muons) 
    {
	if( isGoodMuon(mu) )myLeptons.push_back(&mu);
    }
    for (const pat::Electron &el : *electrons)
    {
         if( isGoodElectron(el) ) myLeptons.push_back(&el);
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

    for (const reco::GenParticle & genParticle: *pruned)
    {
	int ID     = genParticle.pdgId();
        float pt   = genParticle.pt() ;
        float eta  = genParticle.eta();
        int status = genParticle.status();
        if(pt > 10 && fabs(eta) < 2.4 && (abs(ID) == 11 || abs(ID) == 13) && status == 1)myGenLeptons.push_back(&genParticle); 
    }

    const pat::MET &met = mets->front();
    float rawmet      = met.shiftedPt(pat::MET::NoShift, pat::MET::Raw);
    float rawmetPhi   = met.shiftedPhi(pat::MET::NoShift, pat::MET::Raw);
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
       	lepPt_             [ii]  = ii < nleps_ ? myLeptons[ii]->pt()                              : 0;
	lepEta_            [ii]  = ii < nleps_ ? myLeptons[ii]->eta()                             : 0;
	lepPhi_            [ii]  = ii < nleps_ ? myLeptons[ii]->phi()                             : 0;
	lepM_              [ii]  = ii < nleps_ ? myLeptons[ii]->mass()                            : 0;
        lepID_             [ii]  = ii < nleps_ ? myLeptons[ii]->pdgId()                           : 0;
        lepIso_            [ii]  = ii < nleps_ ? LeptonRelIso(myLeptons[ii])                      : 0;  
        lepPtRel_          [ii]  = ii < nleps_ ? PtRel(myLeptons[ii], myJets)                     : 1.e+5;  
        lepGenMatchIndex_  [ii]  = ii < nleps_ ? getLepGenMatchIndex(myLeptons[ii], myGenLeptons) : -1; //function returns -1 if not matched  
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
       	jetPt_      [ii]  = ii < nleps_ ? myLeptons[ii]->pt()    : 0;
	jetEta_     [ii]  = ii < nleps_ ? myLeptons[ii]->eta()   : 0;
	jetPhi_     [ii]  = ii < nleps_ ? myLeptons[ii]->phi()   : 0;
	jetM_       [ii]  = ii < nleps_ ? myLeptons[ii]->mass()  : 0;
        jetBTag_    [ii]  = 0;  // TBI 
    }

    nrjets_                   = (unsigned short) myJets.size();
    for(int ii = 0 ; ii < nrjetsMax; ii++) 
    {
       	rjetPt_      [ii]  = ii < nleps_ ? myLeptons[ii]->pt()    : 0;
	rjetEta_     [ii]  = ii < nleps_ ? myLeptons[ii]->eta()   : 0;
	rjetPhi_     [ii]  = ii < nleps_ ? myLeptons[ii]->phi()   : 0;
	rjetM_       [ii]  = ii < nleps_ ? myLeptons[ii]->mass()  : 0;
        rjetBTag_    [ii]  = 0;  // TBI 
    }

    ngenleps_                   = (unsigned short) myGenLeptons.size();
    for(int ii = 0 ; ii < ngenlepsMax; ii++) 
    {
       	genlepPt_      [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->pt()         : 0;
	genlepEta_     [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->eta()        : 0;
	genlepPhi_     [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->phi()        : 0;
	genlepM_       [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->mass()       : 0;
        genlepID_      [ii]  = ii < ngenleps_ ? myGenLeptons[ii]->pdgId()      : 0;
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

    return res;
}

short SimpleROOT::getLepGenMatchIndex(const reco::Candidate * myRecoLep, vector<const reco::Candidate *> myGenLeptons)
{
    short myIndex = -1; // don't change this
    
    float DR = 1.e+9;

    TLorentzVector myRecoLepP4 = P4(myRecoLep);
    for(unsigned int genIndex = 0; genIndex < myGenLeptons.size(); ++genIndex)
    {
	TLorentzVector myGenLepP4 = P4(myGenLeptons[genIndex]);
        float myDR = myGenLepP4.DeltaR(myRecoLepP4);
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

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleROOT);
