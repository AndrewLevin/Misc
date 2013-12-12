// $Id: merge.C,v 1.1 2009/12/09 09:10:03 loizides Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TRandom.h>
#include <TParameter.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/OutputMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void merge(const char *file1, 
	   const char *file2, 
	   const char *file3,
	   const char *file4,
	   const char *file5,
	   const char *file6,
	   const char *file7,
	   const char *file8,
	   const char *file9,
	   const char *file10,
           const char *prefix="mergetest", 
           UInt_t nev=0)
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  OutputMod *omod = new OutputMod;
  omod->SetFileName(prefix);
  omod->Keep("*");

  Analysis *ana = new Analysis;
  ana->SetSuperModule(omod);
  if (nev)
    ana->SetProcessNEvents(nev);


  assert(file1 != "");
  ana->AddFile(file1);

  if(file2 != "")
    ana->AddFile(file2);
  if(file3 != "")
    ana->AddFile(file3);
  if(file4 != "")
    ana->AddFile(file4);
  if(file5!= "")
    ana->AddFile(file5);
  if(file6 != "")
    ana->AddFile(file6);
  if(file7 != "")
    ana->AddFile(file7);
  if(file8 != "")
    ana->AddFile(file8);
  if(file9 != "")
    ana->AddFile(file9);
  if(file10 != "")
    ana->AddFile(file10);
  
  
  // run the analysis after successful initialisation
  ana->Run(0);
}
