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
//
//

// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD

// system include files
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



        TLorentzVector P4(const reco::Candidate* cand){TLorentzVector p4vec; p4vec.SetPxPyPzE( cand->px(), cand->py(), cand->pz(), cand->energy() ); return p4vec;}


	edm::Service<TFileService> fileService_; 
	TTree *events_;
	// ----------member data 
	// to add more variable you should always edit 3 checkpoints of this file 
	// (give typedef, include it in the branch, initialization for each event)
	
		
        // --- CheckPoint 1 --- 
	UChar_t numVtx_; // 8-bit integer can hold up to 255 vertices, denoted as "b" in the leaflist of the branch
	UChar_t flagBit_; // 8-bit integer can hold up to 255 bits, denoted as "b" in the leaflist of the branch


        // --- enumerate possible flags, once and for ever, *never* change the ordering of the enumed flags
        enum flags_{hasGoodPV=0};
      
};

SimpleROOT::SimpleROOT(const edm::ParameterSet& iConfig)
{
    // --- CheckPoint 2 --- 
    events_ = fileService_->make<TTree>("events","events");
    events_->Branch("numVtx",&numVtx_,"numVtx/b");
    events_->Branch("flagBit",&flagBit_,"flagBit/b");
}

void SimpleROOT::reset()
{
    // --- CheckPoint 3 --- 
    numVtx_ = 0;
    flagBit_ = 0;
}

SimpleROOT::~SimpleROOT() {}

// ------------ method called for each event  ------------
void SimpleROOT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    reset(); // --- initialize all variables to be stored!
    //using namespace edm;

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
    printf("MET: pt %5.1f, ptraw %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
        met.pt(), met.shiftedPt(pat::MET::NoShift, pat::MET::Raw), met.phi(), met.sumEt(),
        met.genMET()->pt(),
        met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));
    // raw met use shiftedPt(pat::MET::NoShift, pat::MET::Raw) because of broken uncorrectedPt() method for Phys14

    printf("\n");



    //  --- sort by pt all objects, the [] is a C++11 lambda func 
    std::sort(myLeptons.begin(), myLeptons.end(), [](const reco::Candidate * a, const reco::Candidate * b){return a->pt() > b->pt();} ); 
    std::sort(myJets.begin(), myJets.end(), [](const reco::Candidate * a, const reco::Candidate * b){return a->pt() > b->pt();} ); 
    std::sort(myRjets.begin(), myRjets.end(), [](const reco::Candidate * a, const reco::Candidate * b){return a->pt() > b->pt();} ); 

    cout << "print all sorted objects " << endl;
    for (auto & lep : myLeptons) cout <<"let pt = "<<  lep->pt() << " (eta, phi) = (" << lep->eta() << " , " << lep->phi() << ")" << " lep pdgId = " << lep->pdgId() << endl;
    for(auto &myjet : myJets) cout << "jet pt = " << myjet->pt() << " (eta, phi) = (" << myjet->eta() << " , " << myjet->phi() << ")" << endl;
    for(auto &myjet : myRjets) cout << "rjet pt = " << myjet->pt() << " (eta, phi) = (" << myjet->eta() << " , " << myjet->phi() << ")" << endl;
    for(auto &myphoton : myPhotons) cout << "rjet pt = " << myphoton->pt() << " (eta, phi) = (" << myphoton->eta() << " , " << myphoton->phi() << ")" << endl;
 
    // --- Fill branches
    numVtx_ = UChar_t(vertices->size());
    if( vtx.isValid() ) flagBit_ |= 1 << hasGoodPV;

    events_->Fill();
}

bool SimpleROOT::isGoodPhoton(const pat::Photon &photon)
{
    bool res = true; // by default is good, unless fails a cut bellow

    cout << __LINE__ << endl;
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
