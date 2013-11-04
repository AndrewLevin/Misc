//--------------------------------------------------------------------------------------------------
// 
// Module that prints out generator level information from BAMBU
//
//==================================================================================================

#ifndef EWKANA_NTUPLER_BAMBUGENDUMPERMOD_HH
#define EWKANA_NTUPLER_BAMBUGENDUMPERMOD_HH

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCParticleFwd.h"

class TTree;

namespace mithep 
{
  class BambuGenDumperMod : public BaseMod 
  {
    public:
      BambuGenDumperMod(const char *name="BambuGenDumperMod", 
                        const char *title="BAMBU Gen Info Dumper");
      ~BambuGenDumperMod();
		      
    protected:
      void                 Begin();
      void                 BeginRun();
      void                 EndRun();
      void                 Process();
      void                 SlaveBegin();
      void                 SlaveTerminate();
      void                 Terminate();

      TString              fPartName;   // branch name of MCParticle collection
      const MCParticleCol *fParticles;  // pointer to generated particle branch

    ClassDef(BambuGenDumperMod, 1)
  };
}
#endif
