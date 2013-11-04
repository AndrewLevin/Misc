//root -l Smurf/ReferenceAnalyses/ZllNormalization.C+\(\"/data/smurf/sixie/data/Run2011_Summer11_EPSHZZV0/mitf-alljets/data_2l.goodlumi1092ipb.root\",\"/data/smurf/sixie/data/Run2011_Summer11_EPSHZZV0/mitf-alljets/zll.root\",1.143\)

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <TClonesArray.h>          
#include "TRandom.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include <iomanip>
#include <TMath.h>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Smurf/Core/SmurfTree.h"
#include "Smurf/Core/LeptonScaleLookup.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TPhoton.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"
#include "MitCommon/MathTools/interface/MathUtils.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 


int    verboseLevel =   0;
const double sigmaB = 0.35;

//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwNtuplerMod");
    if (!dir) {
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}



//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}


//------------------------------------------------------------------------------
// Fill MC PU Distributions
//------------------------------------------------------------------------------
void FillMCPileupDistribution
(
 string mcInputFile    = "/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/dymm.root",
 string outputLabel    = "SmurfV6_dymm"
 )
{
  string outputlabel = ""; if (outputLabel != "") outputlabel = "_" + outputLabel;
  TH1F *NPU = 0;
  TH1F *NPV = 0;

  TFile *infile = new TFile(mcInputFile.c_str(),"READ");
  if (infile->IsOpen()) {
    NPU = (TH1F*)infile->Get("hNPU"); 
    NPV = (TH1F*)infile->Get("hNVtx"); 
  }

  if (!NPU || !NPV) {
    cout << "The NPU and NVtx histograms could not be found from file " << mcInputFile << endl;
    return;
  }

  NormalizeHist(NPU);
  NormalizeHist(NPV);

  //*****************************************************************************
  // Data Loop
  //*****************************************************************************
  TFile *f = new TFile("PileupReweighting.root", "UPDATE");
  f->WriteTObject(NPU, ("NPU"+outputlabel).c_str(), "WriteDelete");
  f->WriteTObject(NPV, ("NPV"+outputlabel).c_str(), "WriteDelete");
  f->Close();
  delete f;

  infile->Close();
  delete infile;
}

//------------------------------------------------------------------------------
// Fill MC PU Distributions
//------------------------------------------------------------------------------
void FillMCPileupDistributionNtuple
(
 string mcInputFile    = "PileupNtuple_DYmm.root",
 string outputLabel    = "SmurfV6_dymm"
 )
{
  string outputlabel = ""; if (outputLabel != "") outputlabel = "_" + outputLabel;

  TFile *infile = new TFile(mcInputFile.c_str(),"READ");
  TH1F *NPU = (TH1F*)infile->Get("hNPU"); 
  TH1F *NPV = (TH1F*)infile->Get("hNVtx"); 

  NormalizeHist(NPU);
  NormalizeHist(NPV);

  //*****************************************************************************
  // Data Loop
  //*****************************************************************************
  TFile *f = new TFile("PileupReweighting.root", "UPDATE");
  f->WriteTObject(NPU, ("NPU"+outputlabel).c_str(), "WriteDelete");
  f->WriteTObject(NPV, ("NPV"+outputlabel).c_str(), "WriteDelete");
  f->Close();
  delete f;

  infile->Close();
  delete infile;
}


//------------------------------------------------------------------------------
// Fill MC PU Distributions
//------------------------------------------------------------------------------
void FillMCPileupDistributionSmurf
(
 string mcInputFile    = "/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/dymm.root",
 string outputLabel    = "SmurfV6_dymm"
 )
{
  string outputlabel = ""; if (outputLabel != "") outputlabel = "_" + outputLabel;

  TChain *chbackground = new TChain("tree");
  chbackground->Add(mcInputFile.c_str());
  TTree *background = (TTree*) chbackground;

  TH1F *NPU = new TH1F(("NPU"+outputlabel).c_str(), ";Number of in-time pileup collisions; NEvents; ", 50, -0.5, 49.5);
  TH1F *NPV = new TH1F(("NPV"+outputlabel).c_str(), ";Number of reconstructed primary vertices; NEvents; ", 50, -0.5, 49.5);

  //----------------------------------------------------------------------------
  UInt_t          cuts;
  UInt_t          dstype;
  UInt_t          nvtx;
  UInt_t          njets;
  UInt_t          event;
  UInt_t          run;
  UInt_t          lumi;
  Float_t         scale1fb;
  LorentzVector*  lep1  = 0;
  LorentzVector*  lep2  = 0;
  LorentzVector*  jet1  = 0;
  LorentzVector*  jet2  = 0;
  Float_t         dPhi;
  Float_t         dR;
  LorentzVector*  dilep = 0;
  UInt_t          type;
  Float_t         pmet;
  Float_t         pTrackMet;
  Float_t         met;
  Float_t         mt;
  Float_t         mt1;
  Float_t         mt2;
  Float_t         dPhiLep1MET;
  Float_t         dPhiLep2MET;
  Float_t         dPhiDiLepMET;
  Float_t         dPhiDiLepJet1;
  Int_t           lq1;
  Int_t           lq2;
  Int_t           lid1;
  Int_t           lid2;
  Int_t           lid3;
  Int_t           processId;
  Float_t         jetLowBtag;
  UInt_t          nSoftMuons;
  Float_t         jet1Btag;
  Float_t         jet2Btag;
  Int_t 	  lep1McId;
  Int_t 	  lep2McId;
  Int_t 	  lep1MotherMcId;
  Int_t 	  lep2MotherMcId;
  UInt_t          npu;



  //*****************************************************************************
  // MC Loop
  //*****************************************************************************

//   background->SetBranchAddress( "cuts"          , &cuts 	  );
//   background->SetBranchAddress( "dstype"        , &dstype	  );
  background->SetBranchAddress( "nvtx"          , &nvtx 	  );
//   background->SetBranchAddress( "njets"         , &njets	  );
//   background->SetBranchAddress( "event"         , &event	  );
//   background->SetBranchAddress( "scale1fb"      , &scale1fb	  );
//   background->SetBranchAddress( "lep1"          , &lep1 	  );
//   background->SetBranchAddress( "lep2"          , &lep2 	  );
//   background->SetBranchAddress( "jet1"          , &jet1 	  );
//   background->SetBranchAddress( "jet2"          , &jet2 	  );
//   background->SetBranchAddress( "dPhi"          , &dPhi 	  );
//   background->SetBranchAddress( "dR"            , &dR		  );
//   background->SetBranchAddress( "dilep"         , &dilep	  );
//   background->SetBranchAddress( "type"          , &type 	  );
//   background->SetBranchAddress( "pmet"          , &pmet 	  );
//   background->SetBranchAddress( "pTrackMet"     , &pTrackMet	  );
//   background->SetBranchAddress( "met"           , &met  	  );
//   background->SetBranchAddress( "mt"            , &mt		  );
//   background->SetBranchAddress( "mt1"           , &mt1  	  );
 //  background->SetBranchAddress( "mt2"           , &mt2  	  );
//   background->SetBranchAddress( "dPhiLep1MET"   , &dPhiLep1MET    );
//   background->SetBranchAddress( "dPhiLep2MET"   , &dPhiLep2MET    );
//   background->SetBranchAddress( "dPhiDiLepMET"  , &dPhiDiLepMET   );
//   background->SetBranchAddress( "dPhiDiLepJet1" , &dPhiDiLepJet1  );
//   background->SetBranchAddress( "lq1"           , &lq1  	  );
//   background->SetBranchAddress( "lq2"           , &lq2  	  );
//   background->SetBranchAddress( "lid1"          , &lid1 	  );
//   background->SetBranchAddress( "lid2"          , &lid2 	  );
//   background->SetBranchAddress( "lid3"          , &lid3 	  );
//   background->SetBranchAddress( "processId"     , &processId	  );
//   background->SetBranchAddress( "jetLowBtag"    , &jetLowBtag	  );
//   background->SetBranchAddress( "nSoftMuons"    , &nSoftMuons	  );
//   background->SetBranchAddress( "jet1Btag"      , &jet1Btag	  );
//   background->SetBranchAddress( "jet2Btag"      , &jet2Btag	  );
//   background->SetBranchAddress( "lep1McId"      , &lep1McId	  );
//   background->SetBranchAddress( "lep2McId"      , &lep2McId	  );
//   background->SetBranchAddress( "lep1MotherMcId", &lep1MotherMcId );
//   background->SetBranchAddress( "lep2MotherMcId", &lep2MotherMcId );
  background->SetBranchAddress( "npu"           , &npu            );

  for (UInt_t i=0; i<background->GetEntries(); i++) {
    
    background->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

    NPU->Fill(npu);
    NPV->Fill(nvtx);
  }
  printf("--- Finished Bgdnal loop\n");


  NormalizeHist(NPU);
  NormalizeHist(NPV);

  //*****************************************************************************
  // Data Loop
  //*****************************************************************************
  TFile *f = new TFile("PileupReweighting.root", "UPDATE");
  f->WriteTObject(NPU, NPU->GetName(), "WriteDelete");
  f->WriteTObject(NPV, NPV->GetName(), "WriteDelete");
  f->Close();

}

void ComputeWeights(string TargetFileName, string SourceHistName, 
                    string Label, Int_t lastPileupBin = 35) {

  TH1F *targetNPU = 0;
  TH1F *sourceNPU = 0;

  TFile *targetfile = new TFile(TargetFileName.c_str(), "READ");
  if (targetfile->IsOpen()) targetNPU = (TH1F*)targetfile->Get("pileup");
  TFile *sourcefile = new TFile("PileupReweighting.root", "UPDATE");
  if (sourcefile->IsOpen()) sourceNPU = (TH1F*)sourcefile->Get(SourceHistName.c_str());

  if (!targetNPU || !sourceNPU) {
    if (!sourceNPU) {
      cout << "Cannot find source hist : " << SourceHistName << "\n"; 
    }
    if (!targetNPU) {
      cout << "Cannot find target hist (\"pileup\") in file : " << TargetFileName << "\n"; 

    }
    return;
  }

  targetNPU->SetDirectory(0);
  NormalizeHist(targetNPU);
  NormalizeHist(sourceNPU);


  TH1D *ReweightingFactor = (TH1D*)targetNPU->Clone("puWeights");
  ReweightingFactor->SetBinContent(0,1.0);

  assert(lastPileupBin <= ReweightingFactor->GetXaxis()->GetNbins());

  for(UInt_t b=1; b < lastPileupBin; ++b) {
    if (sourceNPU->GetBinContent(b)>0) {
      ReweightingFactor->SetBinContent(b, targetNPU->GetBinContent(b) / sourceNPU->GetBinContent(b));
    } else {
      cout << "Error: Bin " << b << " in source histogram has a 0 entry.\n";
      assert(0);
    }
  }

  Double_t sourceLastBin = 0;
  Double_t targetLastBin = 0;
  for(UInt_t b=lastPileupBin; b < ReweightingFactor->GetXaxis()->GetNbins()+2; ++b) {
    sourceLastBin += sourceNPU->GetBinContent(b);
    targetLastBin += targetNPU->GetBinContent(b);
  }
  if (sourceLastBin>0) {
    for (UInt_t b=lastPileupBin; b < ReweightingFactor->GetXaxis()->GetNbins()+2; ++b) {
      ReweightingFactor->SetBinContent(b, targetLastBin / sourceLastBin);
    }
  } else {
    cout << "Error: Last Bin in source histogram is 0.\n";
    assert(0);
  }

  TFile *outputfile = new TFile(("PileupReweighting." + Label + ".root").c_str(), "UPDATE");
  outputfile->WriteTObject(ReweightingFactor, ReweightingFactor->GetName(), "WriteDelete");
  outputfile->Close();
  delete outputfile;
  sourcefile->Close();
  delete sourcefile;
  targetfile->Close();
  delete targetfile;


  sourcefile = new TFile("PileupTargets.root", "UPDATE");
  sourcefile->WriteTObject(targetNPU, ("NPU_Target_" + Label).c_str(), "WriteDelete");
  sourcefile->Close();
  delete sourcefile;

}




//------------------------------------------------------------------------------
// Fill MC PU Distributions
//------------------------------------------------------------------------------
void ValidateReweightingSmurf(

 string mcInputFile       = "/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/dymm.root",
 string dataInputFile     = "/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/data_2l.root",
 string jsonFile          = "/data/smurf/sixie/data/auxiliar/2011Combined.json",
 string PUReweightingFile = "PileupReweighting.42XDYmm_To_Full2011.root",
 string outputLabel       = "SmurfV6DYmm"
 )
{
  string outputlabel = ""; if (outputLabel != "") outputlabel = "_" + outputLabel;

  // ***********************************************************************************************
  // Load Reweight File
  // ***********************************************************************************************
  TFile *fPUS4File = TFile::Open(PUReweightingFile.c_str());
  TH1D *fhDPUS4 = (TH1D*)(fPUS4File->Get("puWeights"));
  assert(fhDPUS4);
  fhDPUS4->SetDirectory(0);
  delete fPUS4File;

  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile(jsonFile.c_str()); 

  // ***********************************************************************************************
  // Load Input files
  // ***********************************************************************************************
  TChain *chbackground = new TChain("tree");
  chbackground->Add(mcInputFile.c_str());
  TTree *background = (TTree*) chbackground;

  TChain *chdata = new TChain("tree");
  chdata->Add(dataInputFile.c_str());
  TTree *data = (TTree*) chdata;

  TH1F *NPV_MC = new TH1F(("NPV_MC"+outputlabel).c_str(), ";Number of reconstructed primary vertices; NEvents; ", 50, -0.5, 49.5);
  TH1F *NPV_Data = new TH1F(("NPV_Data"+outputlabel).c_str(), ";Number of reconstructed primary vertices; NEvents; ", 50, -0.5, 49.5);
  TH1F *Rho_MC = new TH1F(("Rho_MC"+outputlabel).c_str(), ";#rho (Energy density) [GeV]; NEvents; ", 50, -0.5, 49.5);
  TH1F *Rho_Data = new TH1F(("Rho_Data"+outputlabel).c_str(), ";#rho (Energy density) [GeV]; NEvents; ", 50, -0.5, 49.5);
   
  //----------------------------------------------------------------------------
  UInt_t          cuts;
  UInt_t          dstype;
  UInt_t          nvtx;
  UInt_t          njets;
  UInt_t          event;
  UInt_t          run;
  UInt_t          lumi;
  Float_t         scale1fb;
  LorentzVector*  lep1  = 0;
  LorentzVector*  lep2  = 0;
  LorentzVector*  jet1  = 0;
  LorentzVector*  jet2  = 0;
  Float_t         dPhi;
  Float_t         dR;
  LorentzVector*  dilep = 0;
  UInt_t          type;
  Float_t         pmet;
  Float_t         pTrackMet;
  Float_t         met;
  Float_t         mt;
  Float_t         mt1;
  Float_t         mt2;
  Float_t         dPhiLep1MET;
  Float_t         dPhiLep2MET;
  Float_t         dPhiDiLepMET;
  Float_t         dPhiDiLepJet1;
  Int_t           lq1;
  Int_t           lq2;
  Int_t           lid1;
  Int_t           lid2;
  Int_t           lid3;
  Int_t           processId;
  Float_t         jetLowBtag;
  UInt_t          nSoftMuons;
  Float_t         jet1Btag;
  Float_t         jet2Btag;
  Int_t 	  lep1McId;
  Int_t 	  lep2McId;
  Int_t 	  lep1MotherMcId;
  Int_t 	  lep2MotherMcId;
  UInt_t          npu;
  Float_t         rho;



  //*****************************************************************************
  // MC Loop
  //*****************************************************************************

  background->SetBranchAddress( "cuts"          , &cuts 	  );
  background->SetBranchAddress( "nvtx"          , &nvtx 	  );
  background->SetBranchAddress( "lq1"           , &lq1            );
  background->SetBranchAddress( "lq2"           , &lq2            );
  background->SetBranchAddress( "lep1"          , &lep1 	  );
  background->SetBranchAddress( "lep2"          , &lep2 	  );
  background->SetBranchAddress( "dilep"         , &dilep	  );
  background->SetBranchAddress( "type"          , &type 	  );
  background->SetBranchAddress( "npu"           , &npu   	  );
  background->SetBranchAddress( "auxVar0"       , &rho   	  );

  for (UInt_t i=0; i<background->GetEntries(); i++) {
    
    background->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

    //Select Zmumu
    if (!(((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) 
          && ((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))) continue;
    if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
    if( dilep->mass() >= 76 && dilep->mass() <= 106      		 ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
    if( lep2->pt() <= 20	    					 ) continue; // cut on trailing lepton pt
    if(!(type == SmurfTree::mm)                                          ) continue; // mumu only


    //Get PU weight
    //double mynpu = TMath::Min((double)npu,35.499);
    double mynpu = npu;
    Int_t npuxbin = fhDPUS4->GetXaxis()->FindBin(mynpu);
    Double_t weight = fhDPUS4->GetBinContent(npuxbin);

    NPV_MC->Fill(nvtx, weight);
    Rho_MC->Fill(rho, weight);
  }
  printf("--- Finished Bgdnal loop\n");

  //*****************************************************************************
  // Data Loop
  //*****************************************************************************

  data->SetBranchAddress( "cuts"          , &cuts 	  );
  data->SetBranchAddress( "lumi"          , &lumi 	  );
  data->SetBranchAddress( "run"           , &run  	  );
  data->SetBranchAddress( "nvtx"          , &nvtx 	  );
  data->SetBranchAddress( "lq1"           , &lq1          );
  data->SetBranchAddress( "lq2"           , &lq2          );
  data->SetBranchAddress( "lep1"          , &lep1 	  );
  data->SetBranchAddress( "lep2"          , &lep2 	  );
  data->SetBranchAddress( "dilep"         , &dilep	  );
  data->SetBranchAddress( "type"          , &type 	  );
  data->SetBranchAddress( "npu"           , &npu   	  );
  data->SetBranchAddress( "auxVar0"       , &rho   	  );

  for (UInt_t i=0; i<data->GetEntries(); i++) {
    
    data->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)data->GetEntries());

    // not certified run? Skip to next event
    mithep::RunLumiRangeMap::RunLumiPairType rl(run, lumi);      
    if(!rlrm.HasRunLumi(rl)) continue;

    //Select Zmumu
    if (!(((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) 
          && ((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))) continue;
    if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
    if( dilep->mass() >= 76 && dilep->mass() <= 106      		 ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
    if( lep2->pt() <= 20	    					 ) continue; // cut on trailing lepton pt
    if(!(type == SmurfTree::mm)                                          ) continue; // mumu only

    NPV_Data->Fill(nvtx);
    Rho_Data->Fill(rho);
  }
  printf("--- Finished Bgdnal loop\n");


  NormalizeHist(NPV_MC);
  NormalizeHist(NPV_Data);
  NormalizeHist(Rho_MC);
  NormalizeHist(Rho_Data);

  //*****************************************************************************
  // Data Loop
  //*****************************************************************************
  TFile *f = new TFile("PileupReweightingValidation.root", "UPDATE");
  f->WriteTObject(NPV_MC, NPV_MC->GetName(), "WriteDelete");
  f->WriteTObject(NPV_Data, NPV_Data->GetName(), "WriteDelete");
  f->WriteTObject(Rho_MC, Rho_MC->GetName(), "WriteDelete");
  f->WriteTObject(Rho_Data, Rho_Data->GetName(), "WriteDelete");
  f->Close();

}





//------------------------------------------------------------------------------
// Fill MC PU Distributions
//------------------------------------------------------------------------------
void FillPileupNtuple(
 string outputLabel       = "SmurfV6DYmm"
 )
{
  string outputlabel = ""; if (outputLabel != "") outputlabel = "_" + outputLabel;

  // ***********************************************************************************************
  // Ntuples
  // ***********************************************************************************************
  
  UInt_t run,lumi,event;
  Int_t npu,npuMinusOne, npuPlusOne , npv;
  Float_t rho;   
  Bool_t isZMuMu;


  TFile *f = new TFile(("PileupNtuple" +outputlabel + ".root").c_str(), "RECREATE");
  TTree *DataTree = new TTree("DataEvent", "DataEvent");
  DataTree->Branch("run", &run, "run/i");
  DataTree->Branch("lumi", &lumi, "lumi/i");
  DataTree->Branch("event", &event, "event/i");
  DataTree->Branch("npu", &npu, "npu/I");
  DataTree->Branch("npuMinusOne", &npuMinusOne, "npuMinusOne/I");
  DataTree->Branch("npuPlusOne", &npuPlusOne, "npuPlusOne/I");
  DataTree->Branch("npv", &npv, "npv/I");
  DataTree->Branch("rho", &rho, "rho/F");
  DataTree->Branch("isZMuMu", &isZMuMu, "isZMuMu/O");

  TTree *MCTree = new TTree("MCEvent", "MCEvent");
  MCTree->Branch("run", &run, "run/i");
  MCTree->Branch("lumi", &lumi, "lumi/i");
  MCTree->Branch("event", &event, "event/i");
  MCTree->Branch("npu", &npu, "npu/I");
  MCTree->Branch("npuMinusOne", &npuMinusOne, "npuMinusOne/I");
  MCTree->Branch("npuPlusOne", &npuPlusOne, "npuPlusOne/I");
  MCTree->Branch("npv", &npv, "npv/I");
  MCTree->Branch("rho", &rho, "rho/F");
  MCTree->Branch("isZMuMu", &isZMuMu, "isZMuMu/O");



  // ***********************************************************************************************
  // Arrays   
  // ***********************************************************************************************
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  TClonesArray *photonArr = new TClonesArray("mithep::TPhoton");
  Int_t NEvents = 0;


  // ***********************************************************************************************
  // Data
  // ***********************************************************************************************
  vector<string> inputfiles;
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-dmu-m10-v1_LooseLooseSkim.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-dmu-pr-v4_LooseLooseSkim.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-dmu-a05-v1_LooseLooseSkim.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-dmu-pr-v6_LooseLooseSkim.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11b-dmu-pr-v1_LooseLooseSkim.root");
  
  for (UInt_t f = 0; f < inputfiles.size(); ++f) {

    //********************************************************
    // Get Tree
    //********************************************************
    TTree *eventTree = getTreeFromFile(inputfiles[f].c_str(),"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *jetBr;
    TBranch *photonBr;


    //*****************************************************************************************
    //Loop over Data Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
//     eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;
      NEvents++;

      //********************************************************
      // Load the branches
      //********************************************************
//       electronArr->Clear(); 
      muonArr->Clear(); 
//        electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);

      //dilepton preselection
      if (muonArr->GetEntries() < 2) continue;

      //******************************************************************************
      //loop over muon pairs
      //******************************************************************************
      Bool_t goodZMuMuEvent = kFALSE;
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr)[i]);
        Double_t pfRelIso1 = ( mu1->ChargedIso03 + mu1->NeutralIso03_10Threshold ) / mu1->pt;
        for(Int_t j=i+1; j<muonArr->GetEntries(); j++) {
          const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr)[j]);
          Double_t pfRelIso2 = ( mu2->ChargedIso03 + mu2->NeutralIso03_10Threshold) / mu2->pt;
          
          mithep::FourVectorM lepton1;
          mithep::FourVectorM lepton2;
          lepton1.SetCoordinates(mu1->pt, mu1->eta, mu1->phi, 105.658369e-3 );
          lepton2.SetCoordinates(mu2->pt, mu2->eta, mu2->phi, 105.658369e-3 );
          mithep::FourVectorM dilepton = lepton1+lepton2;

          //require loose isolation on at least one of the two leptons to reduce background
          if (!(pfRelIso1 < 0.3 && pfRelIso2 < 0.3)) continue;

          //select Z peak
          if (dilepton.M() > 75.0 && dilepton.M() < 105.0) {
            //Fill These Muons
            goodZMuMuEvent = kTRUE;
            break;

          }
        }
      }
      if (!goodZMuMuEvent) continue;

      run = info->runNum;
      lumi = info->lumiSec;
      event = info->evtNum;
      npu = info->nPUEvents;
      npuMinusOne = 0;
      npuPlusOne = 0;
      npv = info->nPV0;
      rho = double(info->PileupEnergyDensity);
      isZMuMu = goodZMuMuEvent;

      //cout << rho << endl;

      DataTree->Fill();
    }
  }



  // ***********************************************************************************************
  // MC
  // ***********************************************************************************************
  inputfiles.clear();
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0000.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0001.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0002.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0003.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0004.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0005.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0006.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0007.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0008.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0009.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0010.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0011.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0012.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0013.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0014.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0015.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0016.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0017.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0018.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0019.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0020.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0021.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0022.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0023.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0024.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0025.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0026.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0027.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0028.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0029.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0030.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0031.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0032.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0033.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0034.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0035.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0036.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0037.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0038.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0039.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0040.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0041.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0042.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0043.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0044.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0045.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0046.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0047.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0048.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0049.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0050.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0051.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0052.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0053.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0054.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0055.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0056.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0057.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0059.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0061.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0062.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0063.root");

  for (UInt_t f = 0; f < inputfiles.size(); ++f) {

    //********************************************************
    // Get Tree
    //********************************************************
    TTree *eventTree = getTreeFromFile(inputfiles[f].c_str(),"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *jetBr;
    TBranch *photonBr;


    //*****************************************************************************************
    //Loop over MC Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
//     eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;
      NEvents++;

      //********************************************************
      // Load the branches
      //********************************************************
//       electronArr->Clear(); 
      muonArr->Clear(); 
//        electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);

      //dilepton preselection
      if (muonArr->GetEntries() < 2) continue;

      //******************************************************************************
      //loop over muon pairs
      //******************************************************************************
      Bool_t goodZMuMuEvent = kFALSE;
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr)[i]);
        Double_t pfRelIso1 = ( mu1->ChargedIso03 + mu1->NeutralIso03_10Threshold ) / mu1->pt;
        for(Int_t j=i+1; j<muonArr->GetEntries(); j++) {
          const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr)[j]);
          Double_t pfRelIso2 = ( mu2->ChargedIso03 + mu2->NeutralIso03_10Threshold) / mu2->pt;
          
          mithep::FourVectorM lepton1;
          mithep::FourVectorM lepton2;
          lepton1.SetCoordinates(mu1->pt, mu1->eta, mu1->phi, 105.658369e-3 );
          lepton2.SetCoordinates(mu2->pt, mu2->eta, mu2->phi, 105.658369e-3 );
          mithep::FourVectorM dilepton = lepton1+lepton2;

          //require loose isolation on at least one of the two leptons to reduce background
          if (!(pfRelIso1 < 0.3 && pfRelIso2 < 0.3)) continue;

          //select Z peak
          if (dilepton.M() > 75.0 && dilepton.M() < 105.0) {
            //Fill These Muons
            goodZMuMuEvent = kTRUE;
            break;

          }
        }
      }
      if (!goodZMuMuEvent) continue;

      run = info->runNum;
      lumi = info->lumiSec;
      event = info->evtNum;
      npu = info->nPUEvents;
      npuMinusOne = 0;
      npuPlusOne = 0;
      npv = info->nPV0;
      rho = info->PileupEnergyDensity;   
      isZMuMu = goodZMuMuEvent;      
      MCTree->Fill();

    }
  }

  //*****************************************************************************
  // Data Loop
  //*****************************************************************************
  f->WriteTObject(DataTree, "DataTree", "WriteDelete");
  f->WriteTObject(MCTree, "MCTree", "WriteDelete");
  f->Close();

}





//------------------------------------------------------------------------------
// Fill MC PU Distributions
//------------------------------------------------------------------------------
void ValidateReweighting(
 string jsonFile          = "/data/smurf/sixie/data/auxiliar/2011Combined.json",
 string PUReweightingFile = "PileupReweighting.42XDYmm_To_Full2011.root",
 string outputLabel       = "SmurfV6DYmm"
 )
{
  string outputlabel = ""; if (outputLabel != "") outputlabel = "_" + outputLabel;

  // ***********************************************************************************************
  // Load Reweight File
  // ***********************************************************************************************
  TFile *fPUS4File = TFile::Open(PUReweightingFile.c_str());
  TH1D *fhDPUS4 = (TH1D*)(fPUS4File->Get("puWeights"));
  assert(fhDPUS4);
  fhDPUS4->SetDirectory(0);
  delete fPUS4File;

  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile(jsonFile.c_str()); 


  // ***********************************************************************************************
  // Histograms
  // ***********************************************************************************************

  TH1F *NPV_MC = new TH1F(("NPV_MC"+outputlabel).c_str(), ";Number of reconstructed primary vertices; NEvents; ", 50, -0.5, 49.5);
  TH1F *NPV_Data = new TH1F(("NPV_Data"+outputlabel).c_str(), ";Number of reconstructed primary vertices; NEvents; ", 50, -0.5, 49.5);
  TH1F *Rho_MC = new TH1F(("Rho_MC"+outputlabel).c_str(), ";#rho (Energy density) [GeV]; NEvents; ", 50, -0.5, 49.5);
  TH1F *Rho_Data = new TH1F(("Rho_Data"+outputlabel).c_str(), ";#rho (Energy density) [GeV]; NEvents; ", 50, -0.5, 49.5);





  // ***********************************************************************************************
  // Arrays   
  // ***********************************************************************************************
   mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  TClonesArray *photonArr = new TClonesArray("mithep::TPhoton");
  Int_t NEvents = 0;


  // ***********************************************************************************************
  // Data
  // ***********************************************************************************************
  vector<string> inputfiles;
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-dmu-m10-v1_LooseLooseSkim.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-dmu-pr-v4_LooseLooseSkim.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-dmu-a05-v1_LooseLooseSkim.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-dmu-pr-v6_LooseLooseSkim.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11b-dmu-pr-v1_LooseLooseSkim.root");
  
  for (UInt_t f = 0; f < inputfiles.size(); ++f) {

    //********************************************************
    // Get Tree
    //********************************************************
    TTree *eventTree = getTreeFromFile(inputfiles[f].c_str(),"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *jetBr;
    TBranch *photonBr;


    //*****************************************************************************************
    //Loop over Data Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
//     eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;

      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(!rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

      NEvents++;

      //********************************************************
      // Load the branches
      //********************************************************
//       electronArr->Clear(); 
      muonArr->Clear(); 
//        electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);

      //dilepton preselection
      if (muonArr->GetEntries() < 2) continue;

      //******************************************************************************
      //loop over muon pairs
      //******************************************************************************
      Bool_t goodZMuMuEvent = kFALSE;
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr)[i]);
        Double_t pfRelIso1 = ( mu1->ChargedIso03 + mu1->NeutralIso03_10Threshold ) / mu1->pt;
        for(Int_t j=i+1; j<muonArr->GetEntries(); j++) {
          const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr)[j]);
          Double_t pfRelIso2 = ( mu2->ChargedIso03 + mu2->NeutralIso03_10Threshold) / mu2->pt;
          
          mithep::FourVectorM lepton1;
          mithep::FourVectorM lepton2;
          lepton1.SetCoordinates(mu1->pt, mu1->eta, mu1->phi, 105.658369e-3 );
          lepton2.SetCoordinates(mu2->pt, mu2->eta, mu2->phi, 105.658369e-3 );
          mithep::FourVectorM dilepton = lepton1+lepton2;

          //require loose isolation on at least one of the two leptons to reduce background
          if (!(pfRelIso1 < 0.3 && pfRelIso2 < 0.3)) continue;

          //select Z peak
          if (dilepton.M() > 75.0 && dilepton.M() < 105.0) {
            //Fill These Muons
            goodZMuMuEvent = kTRUE;
            break;

          }
        }
      }
      if (!goodZMuMuEvent) continue;

      Rho_Data->Fill(info->PileupEnergyDensity);
      NPV_Data->Fill(info->nPV0); 
    }
  }



  // ***********************************************************************************************
  // MC
  // ***********************************************************************************************
  inputfiles.clear();
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0000.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0001.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0002.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0003.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0004.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0005.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0006.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0007.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0008.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0009.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0010.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0011.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0012.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0013.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0014.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0015.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0016.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0017.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0018.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0019.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0020.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0021.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0022.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0023.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0024.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0025.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0026.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0027.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0028.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0029.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0030.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0031.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0032.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0033.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0034.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0035.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0036.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0037.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0038.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0039.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0040.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0041.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0042.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0043.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0044.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0045.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0046.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0047.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0048.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0049.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0050.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0051.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0052.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0053.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0054.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0055.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0056.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0057.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0059.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0061.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0062.root");
  inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-zmmm20-powheg-v11-pu_noskim_0063.root");

  for (UInt_t f = 0; f < inputfiles.size(); ++f) {

    //********************************************************
    // Get Tree
    //********************************************************
    TTree *eventTree = getTreeFromFile(inputfiles[f].c_str(),"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *jetBr;
    TBranch *photonBr;


    //*****************************************************************************************
    //Loop over MC Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
//     eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;
      NEvents++;

      //********************************************************
      // Load the branches
      //********************************************************
//       electronArr->Clear(); 
      muonArr->Clear(); 
//        electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);

      //dilepton preselection
      if (muonArr->GetEntries() < 2) continue;

      //******************************************************************************
      //loop over muon pairs
      //******************************************************************************
      Bool_t goodZMuMuEvent = kFALSE;
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr)[i]);
        Double_t pfRelIso1 = ( mu1->ChargedIso03 + mu1->NeutralIso03_10Threshold ) / mu1->pt;
        for(Int_t j=i+1; j<muonArr->GetEntries(); j++) {
          const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr)[j]);
          Double_t pfRelIso2 = ( mu2->ChargedIso03 + mu2->NeutralIso03_10Threshold) / mu2->pt;
          
          mithep::FourVectorM lepton1;
          mithep::FourVectorM lepton2;
          lepton1.SetCoordinates(mu1->pt, mu1->eta, mu1->phi, 105.658369e-3 );
          lepton2.SetCoordinates(mu2->pt, mu2->eta, mu2->phi, 105.658369e-3 );
          mithep::FourVectorM dilepton = lepton1+lepton2;

          //require loose isolation on at least one of the two leptons to reduce background
          if (!(pfRelIso1 < 0.3 && pfRelIso2 < 0.3)) continue;

          //select Z peak
          if (dilepton.M() > 75.0 && dilepton.M() < 105.0) {
            //Fill These Muons
            goodZMuMuEvent = kTRUE;
            break;

          }
        }
      }
      if (!goodZMuMuEvent) continue;

       //Get PU weight
      double mynpu = TMath::Min((double)info->nPUEvents,35.499);
      Int_t npuxbin = fhDPUS4->GetXaxis()->FindBin(mynpu);
      Double_t weight = fhDPUS4->GetBinContent(npuxbin);
      
      Rho_MC->Fill(info->PileupEnergyDensity,weight);
      NPV_MC->Fill(info->nPV0,weight); 
    }
  }


   NormalizeHist(NPV_MC);
  NormalizeHist(NPV_Data);
   NormalizeHist(Rho_MC);
  NormalizeHist(Rho_Data);

  //*****************************************************************************
  // Data Loop
  //*****************************************************************************
  TFile *f = new TFile("PileupReweightingValidation.root", "UPDATE");
   f->WriteTObject(NPV_MC, NPV_MC->GetName(), "WriteDelete");
  f->WriteTObject(NPV_Data, NPV_Data->GetName(), "WriteDelete");
   f->WriteTObject(Rho_MC, Rho_MC->GetName(), "WriteDelete");
  f->WriteTObject(Rho_Data, Rho_Data->GetName(), "WriteDelete");
  f->Close();

}





//------------------------------------------------------------------------------
// Fill MC PU Distributions
//------------------------------------------------------------------------------
void ValidateReweightingNtuple(
  string inputFile         = "PileupNtuple_DYmm.root",
  string jsonFile          = "/data/smurf/sixie/data/auxiliar/2011Combined.json",
  string PUReweightingFile = "PileupReweighting.42XDYmm_To_Full2011.root",
  string outputLabel       = "SmurfV6DYmm"
  )
{
  string outputlabel = ""; if (outputLabel != "") outputlabel = "_" + outputLabel;

  // ***********************************************************************************************
  // Load Reweight File
  // ***********************************************************************************************
  TFile *fPUS4File = TFile::Open(PUReweightingFile.c_str());
  TH1D *fhDPUS4 = (TH1D*)(fPUS4File->Get("puWeights"));
  assert(fhDPUS4);
  fhDPUS4->SetDirectory(0);
  delete fPUS4File;

  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile(jsonFile.c_str()); 

  // ***********************************************************************************************
  // Load Input files
  // *********************************************************************************************** 
  TChain *chbackground = new TChain("MCTree");
  chbackground->Add(inputFile.c_str());
  TTree *background = (TTree*) chbackground;

  TChain *chdata = new TChain("DataTree");
  chdata->Add(inputFile.c_str());
  TTree *data = (TTree*) chdata;

  TH1F *NPV_MC = new TH1F(("NPV_MC"+outputlabel).c_str(), ";Number of reconstructed primary vertices; NEvents; ", 50, -0.5, 49.5);
  TH1F *NPV_Data = new TH1F(("NPV_Data"+outputlabel).c_str(), ";Number of reconstructed primary vertices; NEvents; ", 50, -0.5, 49.5);
  TH1F *Rho_MC = new TH1F(("Rho_MC"+outputlabel).c_str(), ";#rho (Energy density) [GeV]; NEvents; ", 50, -0.5, 49.5);
  TH1F *Rho_Data = new TH1F(("Rho_Data"+outputlabel).c_str(), ";#rho (Energy density) [GeV]; NEvents; ", 50, -0.5, 49.5);

  //----------------------------------------------------------------------------
  UInt_t run,lumi,event;
  Int_t npu,npuMinusOne, npuPlusOne , npv;
  Float_t rho;   
  Bool_t isZMuMu;

  //*****************************************************************************
  // MC Loop
  //*****************************************************************************

  background->SetBranchAddress( "run"           , &run 	  );
  background->SetBranchAddress( "lumi"          , &lumi	  );
  background->SetBranchAddress( "event"         , &event  );
  background->SetBranchAddress( "npu"           , &npu	  );
  background->SetBranchAddress( "npuMinusOne"   , &npuMinusOne );
  background->SetBranchAddress( "npuPlusOne"    , &npuPlusOne  );
  background->SetBranchAddress( "npv"           , &npv 	       );
  background->SetBranchAddress( "rho"           , &rho 	       );
  background->SetBranchAddress( "isZMuMu"       , &isZMuMu     );

  for (UInt_t i=0; i<background->GetEntries(); i++) {
    
    background->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)background->GetEntries());

    //Select Zmumu
    if (!isZMuMu) continue;

    //Get PU weight
    double mynpu = TMath::Min((double)npu,25.499);
    Int_t npuxbin = fhDPUS4->GetXaxis()->FindBin(mynpu);
    Double_t weight = fhDPUS4->GetBinContent(npuxbin);

    NPV_MC->Fill(npv, weight);
    Rho_MC->Fill(rho, weight);
  }
  printf("--- Finished Bgdnal loop\n");

  //*****************************************************************************
  // Data Loop
  //*****************************************************************************

  data->SetBranchAddress( "run"           , &run 	  );
  data->SetBranchAddress( "lumi"          , &lumi	  );
  data->SetBranchAddress( "event"         , &event  );
  data->SetBranchAddress( "npu"           , &npu	  );
  data->SetBranchAddress( "npuMinusOne"   , &npuMinusOne );
  data->SetBranchAddress( "npuPlusOne"    , &npuPlusOne  );
  data->SetBranchAddress( "npv"           , &npv 	       );
  data->SetBranchAddress( "rho"           , &rho 	       );
  data->SetBranchAddress( "isZMuMu"       , &isZMuMu     );

  for (UInt_t i=0; i<data->GetEntries(); i++) {
    
    data->GetEntry(i);
    if (i%100000 == 0) printf("--- reading event %5d of %5d\n",i,(int)data->GetEntries());

    // not certified run? Skip to next event
    mithep::RunLumiRangeMap::RunLumiPairType rl(run, lumi);      
    if(!rlrm.HasRunLumi(rl)) continue;

    //Select Zmumu
    if (!isZMuMu) continue;

    NPV_Data->Fill(npv);
    Rho_Data->Fill(rho);
  }
  printf("--- Finished Bgdnal loop\n");


  NormalizeHist(NPV_MC);
  NormalizeHist(NPV_Data);
  NormalizeHist(Rho_MC);
  NormalizeHist(Rho_Data);

  //*****************************************************************************
  // Data Loop
  //*****************************************************************************
  TFile *f = new TFile("PileupReweightingValidation.root", "UPDATE");
  f->WriteTObject(NPV_MC, NPV_MC->GetName(), "WriteDelete");
  f->WriteTObject(NPV_Data, NPV_Data->GetName(), "WriteDelete");
  f->WriteTObject(Rho_MC, Rho_MC->GetName(), "WriteDelete");
  f->WriteTObject(Rho_Data, Rho_Data->GetName(), "WriteDelete");
  f->Close();

}








void DrawValidationPlots(string Label) {
  string label = ""; if (Label != "") label = "_" + Label;

  TFile *f = new TFile("PileupReweightingValidation.root", "READ");
 
  TH1F* NVtx_MC = (TH1F*)f->Get(("NPV_MC"+label).c_str());
  TH1F* NVtx_Data = (TH1F*)f->Get(("NPV_Data"+label).c_str());
  TH1F* Rho_MC = (TH1F*)f->Get(("Rho_MC"+label).c_str());
  TH1F* Rho_Data = (TH1F*)f->Get(("Rho_Data"+label).c_str());


  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(1);


  if (Rho_MC && Rho_Data) {
    tmpLegend->Clear();
    tmpLegend->AddEntry(Rho_MC, "Reweighted MC");  
    tmpLegend->AddEntry(Rho_Data, "Data");
    
    Rho_MC->SetMarkerColor(kRed);
    Rho_MC->SetLineColor(kRed);
    Rho_MC->SetTitle("");
    Rho_MC->SetMaximum(max(double(Rho_MC->GetMaximum()),double(Rho_Data->GetMaximum())) * 1.2);
    Rho_MC->GetYaxis()->SetTitleOffset(1.1);
    Rho_MC->GetXaxis()->SetTitleOffset(1.05);
    
    Rho_MC->Draw("hist");
    Rho_Data->Draw("hist,same");
    tmpLegend->Draw();
    
    cv->SaveAs(("PileupReweightingValidation_Rho" + label + ".png").c_str());
  }

  if (NVtx_MC && NVtx_Data) {
    tmpLegend->Clear();
    tmpLegend->AddEntry(NVtx_MC, "Reweighted MC");  
    tmpLegend->AddEntry(NVtx_Data, "Data");
    
    NVtx_MC->SetMarkerColor(kRed);
    NVtx_MC->SetLineColor(kRed);
    NVtx_MC->SetTitle("");
    NVtx_MC->SetMaximum(max(double(NVtx_MC->GetMaximum()),double(NVtx_Data->GetMaximum())) * 1.2);
    NVtx_MC->GetYaxis()->SetTitleOffset(1.1);
    NVtx_MC->GetXaxis()->SetTitleOffset(1.05);
    
    NVtx_MC->Draw("hist");
    NVtx_Data->Draw("hist,same");
    tmpLegend->Draw();
    
    cv->SaveAs(("PileupReweightingValidation_NVtx" + label + ".png").c_str());
  }



}




void DrawPileupPlots() {

  //*************************************************************************
  //NPU Distributions
  //*************************************************************************
  TFile *sourceFile = new TFile("PileupReweighting.root", "READ");
  TFile *targetFile = new TFile("PileupTargets.root", "READ");
 
  TH1F* NPU_MC = (TH1F*)sourceFile->Get("NPU_Fall11DYmm");
  TH1F* NPU_Target_2011A    = (TH1F*)targetFile->Get("NPU_Target_Fall11DYmm_To_Run2011A");
  TH1F* NPU_Target_2011B    = (TH1F*)targetFile->Get("NPU_Target_Fall11DYmm_To_Run2011B");
  TH1F* NPU_Target_Full2011 = (TH1F*)targetFile->Get("NPU_Target_Fall11DYmm_To_Full2011");

  TCanvas *cv = 0;
  TLegend *tmpLegend = 0;

  if (NPU_MC && NPU_Target_2011A && NPU_Target_2011B && NPU_Target_Full2011) {
    cv = new TCanvas("cv","cv", 800,600);
    tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
    tmpLegend->SetTextSize(0.04);
    tmpLegend->SetBorderSize(0);


    tmpLegend->Clear();
    tmpLegend->AddEntry(NPU_MC,              "Fall 11 MC"   , "L");  
    tmpLegend->AddEntry(NPU_Target_2011A,    "Target Run2011A", "L");  
    tmpLegend->AddEntry(NPU_Target_2011B,    "Target Run2011B", "L");  
    tmpLegend->AddEntry(NPU_Target_Full2011, "Target Full2011", "L");  
  
    NPU_MC->SetMarkerColor(kRed);
    NPU_MC->SetLineColor(kRed);
    NPU_MC->SetTitle("");
    NPU_MC->SetMaximum(max(max(max(double(NPU_MC->GetMaximum()),double(NPU_Target_2011A->GetMaximum())),
                               NPU_Target_2011B->GetMaximum()),NPU_Target_Full2011->GetMaximum()) * 1.2);
    NPU_MC->GetYaxis()->SetTitleOffset(1.1);
    NPU_MC->GetXaxis()->SetTitleOffset(1.05);
    NPU_MC->GetXaxis()->SetRangeUser(-0.5,49.5);
 
    NPU_Target_2011A->SetLineColor(kBlue);
    NPU_Target_2011B->SetLineColor(kMagenta);
    NPU_Target_Full2011->SetLineColor(kBlack);
    NPU_MC->SetLineWidth(2);
    NPU_Target_2011A->SetLineWidth(2);
    NPU_Target_2011B->SetLineWidth(2);
    NPU_Target_Full2011->SetLineWidth(2);
 
 
    NPU_MC->Draw("hist");
    NPU_Target_2011A->Draw("hist,same");
    NPU_Target_2011B->Draw("hist,same");
    NPU_Target_Full2011->Draw("hist,same");
    tmpLegend->Draw(); 

    cv->SaveAs("NPUDistributions.png");
  } else {
    cout << "Could not retrieve all the source or target pileup distributions.\n";
  }


  //*************************************************************************
  //Weights
  //*************************************************************************
  TFile *f = 0;
  TH1F* Weights_Run2011A = 0;
  TH1F* Weights_Run2011B = 0;
  TH1F* Weights_Full2011 = 0;

  f = new TFile("PileupReweighting.Summer11_To_Run2011A.root", "READ");
  Weights_Run2011A = (TH1F*)f->Get("puWeights");
  f = new TFile("PileupReweighting.Summer11_To_Run2011B.root", "READ");
  Weights_Run2011B = (TH1F*)f->Get("puWeights");
  f = new TFile("PileupReweighting.Summer11_To_Full2011.root", "READ");
  Weights_Full2011 = (TH1F*)f->Get("puWeights");

  if (Weights_Run2011A && Weights_Run2011B && Weights_Full2011) { 
    tmpLegend = new TLegend(0.7,0.75,0.93,0.90);   
    tmpLegend->SetTextSize(0.04);
    tmpLegend->SetBorderSize(0);

    tmpLegend->Clear();
    tmpLegend->AddEntry(Weights_Run2011A,    "Run2011A", "L");  
    tmpLegend->AddEntry(Weights_Run2011B,    "Run2011B", "L");  
    tmpLegend->AddEntry(Weights_Full2011,    "Full2011", "L");  
  
    Weights_Run2011A->SetMarkerColor(kBlue);
    Weights_Run2011A->SetLineColor(kBlue);
    Weights_Run2011A->SetTitle("");
    Weights_Run2011A->SetMaximum(max(max(double(Weights_Run2011A->GetMaximum()),double(Weights_Run2011B->GetMaximum())),
                                     Weights_Full2011->GetMaximum()) * 1.2);
    Weights_Run2011A->GetYaxis()->SetTitleOffset(1.1);
    Weights_Run2011A->GetXaxis()->SetTitleOffset(1.05);
    Weights_Run2011A->GetXaxis()->SetRangeUser(-0.5,49.5);
 
    Weights_Run2011B->SetLineColor(kMagenta);
    Weights_Full2011->SetLineColor(kBlack);
    Weights_Run2011A->SetLineWidth(2);
    Weights_Run2011B->SetLineWidth(2);
    Weights_Full2011->SetLineWidth(2);

    Weights_Run2011A->Draw("hist");
    Weights_Run2011B->Draw("hist,same");
    Weights_Full2011->Draw("hist,same");
    tmpLegend->Draw();

    cv->SaveAs("ReweightingFactors_Summer11Source.png");
  } else {
    cout << "Could not find all pileup reweighting histograms for Summer11 -> Run2011. \n";
  }


  f = new TFile("/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedHalfBkgWP/auxiliar/PileupReweighting.Fall11DYmm_To_Run2011A.root", "READ");
  Weights_Run2011A = (TH1F*)f->Get("puWeights");
  f = new TFile("/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedHalfBkgWP/auxiliar/PileupReweighting.Fall11DYmm_To_Run2011B.root", "READ");
  Weights_Run2011B = (TH1F*)f->Get("puWeights");
  f = new TFile("/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedHalfBkgWP/auxiliar/PileupReweighting.Fall11DYmm_To_Full2011.root", "READ");
  Weights_Full2011 = (TH1F*)f->Get("puWeights");

  if (Weights_Run2011A && Weights_Run2011B && Weights_Full2011) { 
    tmpLegend = new TLegend(0.7,0.75,0.93,0.90);   
    tmpLegend->SetTextSize(0.04);
    tmpLegend->SetBorderSize(0);


    tmpLegend->Clear();
    tmpLegend->AddEntry(Weights_Run2011A,    "Run2011A", "L");  
    tmpLegend->AddEntry(Weights_Run2011B,    "Run2011B", "L");  
    tmpLegend->AddEntry(Weights_Full2011,    "Full2011", "L");  
  
    Weights_Run2011A->SetMarkerColor(kBlue);
    Weights_Run2011A->SetLineColor(kBlue);
    Weights_Run2011A->SetTitle("");
    Weights_Run2011A->SetMaximum(max(max(double(Weights_Run2011A->GetMaximum()),double(Weights_Run2011B->GetMaximum())),
                                     Weights_Full2011->GetMaximum()) * 1.2);
    Weights_Run2011A->GetYaxis()->SetTitleOffset(1.1);
    Weights_Run2011A->GetXaxis()->SetTitleOffset(1.05);
    Weights_Run2011A->GetXaxis()->SetRangeUser(-0.5,49.5);
 
    Weights_Run2011B->SetLineColor(kMagenta);
    Weights_Full2011->SetLineColor(kBlack);
    Weights_Run2011A->SetLineWidth(2);
    Weights_Run2011B->SetLineWidth(2);
    Weights_Full2011->SetLineWidth(2);

    Weights_Run2011A->Draw("hist");
    Weights_Run2011B->Draw("hist,same");
    Weights_Full2011->Draw("hist,same");
    tmpLegend->Draw();
    cv->SaveAs("ReweightingFactors_Fall11Source.png");
  } else {
    cout << "Could not find all pileup reweighting histograms for Fall11 -> Run2011. \n";
  }



}



void PileupReweighting(Int_t Option = -1) {

  //**************************************************************************
  // Get Source Distributions from MC
  // Provide input file with the NPU histogram for ALL generated events
  //**************************************************************************
  if (Option == -1 || Option == 0) {
    cout << "**********************************************************" << endl;
    cout << "Retrieving Source Pileup Distributions from Monte Carlo..." << endl;
    FillMCPileupDistribution("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_s11-zmmm20-powheg-v11-pu_noskim.root","Summer11DYmm");
    FillMCPileupDistribution("/data/blue/sixie/Smurf/MVAIDIsoCombinedHalfBkgWP/mc/AllNtuple_f11-zmmm20-powheg-v14b-pu_noskim.root","Fall11DYmm");
  }
  

  //**************************************************************************
  // Load Target Distributions and produce reweighting histograms
  //**************************************************************************
  if (Option == -1 || Option == 1) {

    //Summer11 -> Run2011A
    ComputeWeights("/data/smurf/sixie/Pileup/PUTarget.Run2011A.160404-172619.root","NPU_Summer11DYmm","Summer11DYmm_To_Run2011A", 35);
    
    //Summer11 -> Run2011B
    ComputeWeights("/data/smurf/sixie/Pileup/PUTarget.Run2011B.175832-180252.root","NPU_Summer11DYmm","Summer11DYmm_To_Run2011B", 35);
    
    //Summer11 -> Run2011-Full
    ComputeWeights("/data/smurf/sixie/Pileup/PUTarget.Full2011.160404-180252.root","NPU_Summer11DYmm","Summer11DYmm_To_Full2011", 35);

    //Fall11 -> Run2011A
    ComputeWeights("/data/smurf/sixie/Pileup/PUTarget.Run2011A.160404-172619.root","NPU_Fall11DYmm","Fall11DYmm_To_Run2011A", 35);
    
    //Fall11 -> Run2011B
    ComputeWeights("/data/smurf/sixie/Pileup/PUTarget.Run2011B.175832-180252.root","NPU_Fall11DYmm","Fall11DYmm_To_Run2011B", 35);
    
    //Fall11 -> Run2011-Full
    ComputeWeights("/data/smurf/sixie/Pileup/PUTarget.Full2011.160404-180252.root","NPU_Fall11DYmm","Fall11DYmm_To_Full2011", 35);

  }



  //**************************************************************************
  // Produce Validation plots of reweighting : From SMURF
  //**************************************************************************
  if (Option == -1 || Option == 2) {
    //Run2011-Full
    ValidateReweightingSmurf("/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedHalfBkgWP/mitf-alljets/dymm.root",
                             "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedHalfBkgWP/mitf-alljets/data_2l.root",
                             "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedHalfBkgWP/auxiliar/2011Combined.json",
                             "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedHalfBkgWP/auxiliar/PileupReweighting.Fall11DYmm_To_Full2011.root",
                             "Fall11DYmm_MVAIDIsoCombinedHalfBkgWP_Full2011");
    DrawValidationPlots("Fall11DYmm_MVAIDIsoCombinedHalfBkgWP_Full2011");    
  }
  
  //**************************************************************************
  // Plots for Note
  //**************************************************************************
  if (Option == -1 || Option == 3) {
    DrawPileupPlots();
  }

}
 
