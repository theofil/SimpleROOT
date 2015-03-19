// -*- C++ -*-
//
// Package:    Tools/SimpleROOT
// Class:      SimpleROOT
// 
/**\class SimpleROOT SimpleROOT.cc Tools/SimpleROOT/plugins/SimpleROOT.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Konstantinos Theofilatos
//         Created:  Fri, 20 Feb 2015 16:08:35 GMT

// many thanks to https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
//
// Principles
//
// 1) producer's code should be intuitive and readable by anybody who knows English
// 2) don't skip any event for any reason, save them all on disk even if ~empty
// 3) be flat and use basic types for brances float, (unsigned) int,  (unsigned) short  and STD vectors
// 4) don't introduce new structures and classes, unless it's not possible to go ahead without them
// 5) be laconic and use smart commands when possible C++11, exploit ROOT buildin and STD
// 6) event selection and logic should be undisturbed by technicalities
// 7) don't use int when short is OK, double is forbidden

#include <memory>

// user include files
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

#include "TTree.h"
#include "TLorentzVector.h"

using namespace std;

class SimpleROOT : public edm::EDAnalyzer {
    public:
    	explicit SimpleROOT(const edm::ParameterSet&);
      	~SimpleROOT();

    private:
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void reset();
        virtual bool isGoodMuon(const pat::Muon &mu);
        virtual bool isGoodElectron(const pat::Electron &el);
        virtual bool isGoodVertex(const reco::Vertex &PV);
        virtual bool isGoodJet(const pat::Jet &myJet);
        virtual bool isGoodPhoton(const pat::Photon &myPhoton);
        virtual void sortByPt(vector<const reco::Candidate *> myRecoCand);

        TLorentzVector P4(const reco::Candidate* cand){TLorentzVector p4vec; p4vec.SetPxPyPzE( cand->px(), cand->py(), cand->pz(), cand->energy() ); return p4vec;}

	edm::Service<TFileService> fileService_; 
	TTree *events_;
	
        // ---  CP1 --- 
	bool goodVtx_; 
	unsigned short nVtx_; 

        unsigned short nLeps_;
        unsigned short nJets_;
        unsigned short nRjets_;
        unsigned short nPhos_;

        float l1Pt_;         // leading lepton
        float l2Pt_;         // trailing lepton
        float l1Eta_;
        float l2Eta_;
        float l1Phi_;
        float l2Phi_;
        float l1Iso_;
        float l2Iso_;
        float l1l2DPhi_;
        float l1l2DR_;
        float l1l2Pt_;
        float l1l2M_;
        float l1l2Eta_;
        float l1l2Phi_;
     
        float MET_;          // raw pf-met
        float METPhi_;
        float t1MET_;
        float t1METPhi_;
        float sumEt_;
        float t1sumEt_;

        float vHT_;           // pt of the recoil vector of all objects excluding the dilepton [= -MET - l1l2] 
        float t1vHT_;         // as vHT but with t1 correction
	float jvHT_;          // recoil of hard jets and subleading leptons
        
        short l1pdgID_;  // stores the produce of charge (+/-) and id (1,2) for (emu)
        short l2pdgID_;  // i.e. -1 = electron +1 positron, -2 = muon
      
};

SimpleROOT::SimpleROOT(const edm::ParameterSet& iConfig)
{
    // --- CP2 --- 
    events_ = fileService_->make<TTree>("events","events");
    events_->Branch("goodVtx"          ,&goodVtx_               ,"goodVtx/O");
    events_->Branch("nVtx"             ,&nVtx_                  ,"nVtx/s");
    events_->Branch("nLeps"            ,&nLeps_                 ,"nLeps/s");
    events_->Branch("nJets"            ,&nJets_                 ,"nJets/s");
    events_->Branch("nRjets"           ,&nRjets_                ,"nRjets/s");
    events_->Branch("nPhos"            ,&nPhos_                 ,"nPhos/s");
    events_->Branch("l1Pt"             ,&l1Pt_                  ,"l1Pt/F    ");
    events_->Branch("l2Pt"             ,&l2Pt_                  ,"l2Pt/F    ");
    events_->Branch("l1Eta"            ,&l1Eta_                 ,"l1Eta/F   ");
    events_->Branch("l2Eta"            ,&l2Eta_                 ,"l2Eta/F   ");
    events_->Branch("l1Phi"            ,&l1Phi_                 ,"l1Phi/F   ");
    events_->Branch("l2Phi"            ,&l2Phi_                 ,"l2Phi/F   ");
    events_->Branch("l1Iso"            ,&l1Iso_                 ,"l1Iso/F   ");
    events_->Branch("l2Iso"            ,&l2Iso_                 ,"l2Iso/F   ");
    events_->Branch("l1l2DPhi"         ,&l1l2DPhi_              ,"l1l2DPhi/F");
    events_->Branch("l1l2DR"           ,&l1l2DR_                ,"l1l2DR/F  ");
    events_->Branch("l1l2Pt"           ,&l1l2Pt_                ,"l1l2Pt/F  ");
    events_->Branch("l1l2M"            ,&l1l2M_                 ,"l1l2M/F   ");
    events_->Branch("l1l2Eta"          ,&l1l2Eta_               ,"l1l2Eta/F ");
    events_->Branch("l1l2Phi"          ,&l1l2Phi_               ,"l1l2Phi/F ");
    events_->Branch("MET"              ,&MET_                   ,"MET/F     ");
    events_->Branch("METPhi"           ,&METPhi_                ,"METPhi/F  ");
    events_->Branch("t1MET"            ,&t1MET_                 ,"t1MET/F   ");
    events_->Branch("t1METPhi"         ,&t1METPhi_              ,"t1METPhi/F");
    events_->Branch("sumEt"            ,&sumEt_                 ,"sumEt/F   ");
    events_->Branch("t1sumEt"          ,&t1sumEt_               ,"t1sumEt/F ");
    events_->Branch("vHT"              ,&vHT_                   ,"vHT/F      ");
    events_->Branch("t1vHT"            ,&t1vHT_                 ,"t1vHT/F    ");
    events_->Branch("jvHT"             ,&jvHT_                  ,"jvHT/F     ");
    events_->Branch("l1pdgID"          ,&l1pdgID_               ,"l1pdgID/S");
    events_->Branch("l2pdgID"          ,&l2pdgID_               ,"l2pdgID/S");
    //events_->Branch(""          ,&               ,"");
}

void SimpleROOT::reset()
{
    // --- CP3--- 
    goodVtx_                 = 0;
    nVtx_                    = 0;
    nLeps_                   = 0;
    nJets_                   = 0;
    nRjets_                  = 0;
    nPhos_                   = 0;
    l1Pt_                    = 0;         
    l2Pt_                    = 0;         
    l1Eta_                   = 0;
    l2Eta_                   = 0;
    l1Phi_                   = 0;
    l2Phi_                   = 0;
    l1Iso_                   = 0;
    l2Iso_                   = 0;
    l1l2DPhi_                = 0;
    l1l2DR_                  = 0;
    l1l2Pt_                  = 0;
    l1l2M_                   = 0;
    l1l2Eta_                 = 0;
    l1l2Phi_                 = 0;
    
    MET_                     = 0;          
    METPhi_                  = 0;
    t1MET_                   = 0;
    t1METPhi_                = 0;
    sumEt_                   = 0;
    t1sumEt_                 = 0;

    vHT_                     = 0;            
    t1vHT_                   = 0;         
    jvHT_                    = 0;          
    l1pdgID_                 = 0;
    l2pdgID_                 = 0;
}

SimpleROOT::~SimpleROOT() {}

// ------------ method called for each event  ------------
void SimpleROOT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    reset(); // --- initialize all variables to be stored!

    // --- get MiniAOD collections -- names of the collection are hard-coded (no need to read them from python as they don't change so frequently)
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel("offlineSlimmedPrimaryVertices", vertices);  

    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByLabel("slimmedMuons", muons);

    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByLabel("slimmedElectrons", electrons);

    edm::Handle<pat::JetCollection> jets;
    iEvent.getByLabel("slimmedJets", jets);
   
    edm::Handle<pat::METCollection> mets;
    iEvent.getByLabel("slimmedMETs", mets);

    edm::Handle<pat::PhotonCollection> photons;
    iEvent.getByLabel("slimmedPhotons", photons);

    vector<const reco::Candidate *> myLeptons; // in this container we will store all selected RECO electrons and RECO muons 
    vector<const reco::Candidate *> myJets; // in this container we will store all prompt jets (PV)
    vector<const reco::Candidate *> myRjets; // in this container we will all rejected prompt jets (PV) -- due to DR matching with myLeptons
    vector<const reco::Candidate *> myPhotons; // in this container we will store all photons

    // --- loop inside the objects
    reco::Vertex vtx;
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
	if( isGoodJet(myjet) && isLeptonMatched ) myRjets.push_back(&myjet);
    }

    for (const pat::Photon &photon : *photons) 
    {
	if( isGoodPhoton(photon) ) myPhotons.push_back(&photon); 
    }

    const pat::MET &met = mets->front();
    float rawMET      = met.shiftedPt(pat::MET::NoShift, pat::MET::Raw);
    float rawMETPhi   = met.shiftedPhi(pat::MET::NoShift, pat::MET::Raw);
    float rawMETSumEt = met.shiftedSumEt(pat::MET::NoShift, pat::MET::Raw);
    float t1MET       = met.pt();
    float t1METPhi    = met.phi();
    float t1METSumEt  = met.sumEt();

    sortByPt(myLeptons);
    sortByPt(myJets);
    sortByPt(myRjets);
    sortByPt(myPhotons);

//    for(auto &mylep : myLeptons) cout << "mylep->pdgID() = " << mylep->pdgId() << endl;
 
    TLorentzVector l1,l2;
    if(myLeptons.size() >=1) l1 = P4(myLeptons[0]);
    if(myLeptons.size() >=2) l2 = P4(myLeptons[1]);

    TLorentzVector METVector, t1METVector;
    METVector.SetPtEtaPhiE  (rawMET  , 0, rawMETPhi, rawMET );
    t1METVector.SetPtEtaPhiE(t1MET   , 0, t1METPhi , t1MET  );

    TLorentzVector HTVector, t1HTVector, jHTVector; 
   
    HTVector   = -METVector -l1 -l2;
    t1HTVector = -t1METVector -l1 -l2;
    
    for(auto &myjet : myJets)    jHTVector += P4(myjet);
    for(auto &mylep : myLeptons) jHTVector -= P4(mylep);

    // --- Fill branches CP4
    nVtx_                    = (unsigned short) vertices->size();
    goodVtx_                 = vtx.isValid() ? true:false;
    nLeps_                   = (unsigned short) myLeptons.size();
    nJets_                   = (unsigned short) myJets.size();
    nRjets_                  = (unsigned short) myRjets.size();
    nPhos_                   = (unsigned short) myPhotons.size();

    l1Pt_                    = nLeps_>=1 ? l1.Pt() : 0;         
    l2Pt_                    = nLeps_>=2 ? l2.Pt() : 0;         
    l1Eta_                   = nLeps_>=1 && l1Pt_ >1.e-3 ? l1.Eta() : 0;
    l2Eta_                   = nLeps_>=2 && l2Pt_ >1.e-3 ? l2.Eta() : 0;
    l1Phi_                   = nLeps_>=1 ? l1.Phi() : 0;
    l2Phi_                   = nLeps_>=2 ? l2.Phi() : 0;
    l1Iso_                   = 0;
    l2Iso_                   = 0;
    l1l2DPhi_                = nLeps_>=2 ? l1.DeltaPhi(l2) : 0;
    l1l2DR_                  = nLeps_>=2 ? l1.DeltaR(l2) : 0;
    l1l2Pt_                  = nLeps_>=2 ? (l1+l2).Pt() : 0;
    l1l2M_                   = nLeps_>=2 ? (l1+l2).M() : 0;
    l1l2Eta_                 = nLeps_>=2 && l1l2Pt_ > 1.e-3 ? (l1+l2).Eta() : 0;
    l1l2Phi_                 = nLeps_>=2 ? (l1+l2).Phi() : 0;
    
    MET_                     = rawMET;          
    METPhi_                  = rawMETPhi;
    sumEt_                   = rawMETSumEt;
    t1MET_                   = t1MET;
    t1METPhi_                = t1METPhi;
    t1sumEt_                 = t1METSumEt;

    vHT_                     = HTVector.Pt();            
    t1vHT_                   = t1HTVector.Pt();         
    jvHT_                    = jHTVector.Pt();          
    
    l1pdgID_                 = nLeps_ >=1 ? (unsigned short) myLeptons[0]->pdgId() : 0;
    l2pdgID_                 = nLeps_ >=2 ? (unsigned short) myLeptons[1]->pdgId() : 0;
    events_->Fill();
}

void SimpleROOT::sortByPt(vector<const reco::Candidate *> myRecoCand)
{
    //  --- sort by pt all objects, the [] is a C++11 lambda func 
    std::sort(myRecoCand.begin(), myRecoCand.end(), [](const reco::Candidate * a, const reco::Candidate * b){return a->pt() > b->pt();} );
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

    return res;
}

bool SimpleROOT::isGoodElectron(const pat::Electron &el)
{
    bool res = true; // by default is good, unless fails a cut bellow

    if(el.pt() < 10) res = false;
    if(fabs(el.eta()) > 2.4) res = false;

    return res;
}

bool SimpleROOT::isGoodVertex(const reco::Vertex &vtx)
{
  // The "good vertex" selection is borrowed Ilya Kravchenko who borrowed from Giovanni Zevi Della Porta
  return ( !(vtx.chi2()==0 && vtx.ndof()==0) && vtx.ndof()>=4. && vtx.position().Rho()<=2.0 && fabs(vtx.position().Z())<=24.0 ) ? true : false; 
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleROOT);
