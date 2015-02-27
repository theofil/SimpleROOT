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

#include "TTree.h"

using namespace std;

class SimpleROOT : public edm::EDAnalyzer {
    public:
    	explicit SimpleROOT(const edm::ParameterSet&);
      	~SimpleROOT();

    private:
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void reset();

	edm::Service<TFileService> fileService_; 
	TTree *events_;
	// ----------member data 
	// to add more variable you should always edit 3 checkpoints of this file 
	// (give typedef, include it in the branch, initialization for each event)
	
		
        // --- CheckPoint 1 --- 
	UChar_t numVtx_; // 8-bit integer can hold up to 255 vertices, denoted as "b" in the leaflist of the branch
	UChar_t flagBit_; // 8-bit integer can hold up to 255 bits, denoted as "b" in the leaflist of the branch


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

    // --- loop inside the objects
    const reco::Vertex &PV = vertices->front();  // get event's primary vertex
    if(PV.isValid()) flagBit_ |= 1 << hasGoodPV;
    numVtx_ = UChar_t(vertices->size());
    
    events_->Fill();
}



//define this as a plug-in
DEFINE_FWK_MODULE(SimpleROOT);
