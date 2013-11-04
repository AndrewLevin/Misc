#include "EWKAna/Ntupler/interface/BambuGenDumperMod.hh"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"

#include <iostream>
#include <iomanip>
#include <vector>

using namespace mithep;
using namespace std;

ClassImp(mithep::BambuGenDumperMod)

//--------------------------------------------------------------------------------------------------
BambuGenDumperMod::BambuGenDumperMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPartName(Names::gkMCPartBrn),
  fParticles(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
BambuGenDumperMod::~BambuGenDumperMod(){}

//--------------------------------------------------------------------------------------------------
void BambuGenDumperMod::SlaveBegin()
{
  //
  // Request BAMBU branches
  //
  ReqBranch(fPartName, fParticles);
}

//--------------------------------------------------------------------------------------------------
void BambuGenDumperMod::Process()
{
  LoadBranch(fPartName);

  if(!fParticles) return;

  IncNEventsProcessed(); 

  cout << "Event: " << GetNEventsProcessed() << endl;
  cout << "Particles: " << fParticles->GetEntries() << endl;
  cout << endl;
  cout << "Index";
  cout << setw(8) << "PDG ID";
  cout << setw(8) << "Status";
  cout << setw(8) << "From";
  cout << setw(10) << "E  ";
  cout << setw(10) << "px  ";
  cout << setw(10) << "py  ";
  cout << setw(10) << "pz  ";
  cout << setw(25) << "decay vertex     " << endl;
  cout << "____________________________________________________________________________________________________" << endl;

  //
  // Scan generator level info
  //
  vector<const MCParticle*> index;
    
  for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);
    index.push_back(p);
    
    cout << setw(5) << i+1;
    cout << setw(9) << p->PdgId();
    cout << setw(7) << p->Status();
    
    if(p->HasMother()) {
      for(UInt_t idx=0; idx<index.size(); idx++) {
        const MCParticle *m = p->Mother();
	if(m == index[idx]) {
	  cout << setw(8) << idx+1;
	  break;
	}
      }
    } else {
      cout << setw(8) << " ";
    }
    
    cout << setw(10) << setprecision(3) << p->E();
    cout << setw(10) << setprecision(3) << p->Px();
    cout << setw(10) << setprecision(3) << p->Py();
    cout << setw(10) << setprecision(3) << p->Pz();
    ThreeVector v = p->DecayVertex();
    if( (v.X()!=0) || (v.Y()!=0) || (v.Z()!=0) ) {  // stable particles have decay vertex (0,0,0)
      cout << "   (";
      cout << setw(7) << setprecision(3) << v.X();
      cout << ",";
      cout << setw(7) << setprecision(3) << v.Y();
      cout << ",";
      cout << setw(7) << setprecision(3) << v.Z();
      cout << ")";
    }
    cout << endl;
  }

  cout << "====================================================================================================" << endl;
}

//--------------------------------------------------------------------------------------------------
void BambuGenDumperMod::SlaveTerminate()
{
  cout << endl;
  cout << "Finished with " << GetNEventsProcessed() << " events processed!" << endl;
}

//--------------------------------------------------------------------------------------------------
void BambuGenDumperMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont
  // do anything here.
}

//--------------------------------------------------------------------------------------------------
void BambuGenDumperMod::BeginRun()
{
  // Run startup code on the client machine. For this module, we dont
  // do anything here.
}

//--------------------------------------------------------------------------------------------------
void BambuGenDumperMod::EndRun()
{
  // Run startup code on the client machine. For this module, we dont
  // do anything here.
}

//--------------------------------------------------------------------------------------------------
void BambuGenDumperMod::Terminate()
{
  // Run finishing code on the client computer. 
  // For this module, we dont do anything here.
}
