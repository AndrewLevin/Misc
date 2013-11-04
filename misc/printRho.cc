#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#

#include <cmath>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace edm;
using namespace std;
class printRho : public edm::EDAnalyzer {
   public:
      explicit printRho(const edm::ParameterSet&);
      ~printRho();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
       
  ParameterSet conf_;
 
  unsigned int ev;

};

printRho::printRho(const edm::ParameterSet& iConfig):
  conf_(iConfig)
{



}


printRho::~printRho()
{
 
}


void
printRho::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	
  Handle<double> hRho1;
  Handle<double> hRho2;
	
  edm::InputTag tag1("kt6PFJetsForIso1","rho");
  edm::InputTag tag2("kt6PFJetsForIso2","rho");
  iEvent.getByLabel(tag1,hRho1);
  iEvent.getByLabel(tag2,hRho2);
	
	cout << "*hRho1 = " << *hRho1 << endl;
	cout << "*hRho2 = " << *hRho2 << endl;	

}



void 
printRho::beginJob(const edm::EventSetup&)
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
printRho::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(printRho);
