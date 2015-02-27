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
        

//	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
//	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void reset();
//	virtual void endJob() override;

	edm::Service<TFileService> fileService_; 
	TTree *events_;
	// ----------member data 
	// to add more variable you should always edit 3 checkpoints of this file 
	// (give typedef, include it in the branch, initialization for each event)
	
		
        // --- CheckPoint 1 --- 
	UChar_t numVtx_; // 8-bit integer denoted as "b" in the leaflist of the branch
      
};

SimpleROOT::SimpleROOT(const edm::ParameterSet& iConfig)
{
    // --- CheckPoint 2 --- 
    events_ = fileService_->make<TTree>("events","events");
    events_->Branch("numVtx_",&numVtx_,"numVtx_/b");
}

void SimpleROOT::reset()
{
    // --- CheckPoint 3 --- 
    numVtx_ = 0;
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
    bool PVisValid = PV.isValid();
    numVtx_ = UChar_t(vertices->size());
    
    cout << "PVisValid = " << PVisValid << endl;

    events_->Fill();
}


//// ------------ method called once each job just before starting event loop  ------------
//void 
//SimpleROOT::beginJob()
//{
//}
//
//// ------------ method called once each job just after ending the event loop  ------------
//void 
//SimpleROOT::endJob() 
//{
//}
//
//
//// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//void
//SimpleROOT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//  //The following says we do not know what parameters are allowed so do no validation
//  // Please change this to state exactly what you do use, even if it is no parameters
//  edm::ParameterSetDescription desc;
//  desc.setUnknown();
//  descriptions.addDefault(desc);
//}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleROOT);
