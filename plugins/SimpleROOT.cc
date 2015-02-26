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

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      TTree *events;

      // ----------member data ---------------------------
      UShort_t numVtx; // 16-bit integer 
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SimpleROOT::SimpleROOT(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


SimpleROOT::~SimpleROOT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SimpleROOT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

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
    numVtx = UShort_t(vertices->size());
    
    cout << "PVisValid = " << PVisValid << endl;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
SimpleROOT::beginJob()
{
    events = new TTree("events","events");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimpleROOT::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
SimpleROOT::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
SimpleROOT::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
SimpleROOT::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
SimpleROOT::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimpleROOT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleROOT);
