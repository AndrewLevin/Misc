#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "StandardPlotVBSData.C"

void finalPlotVBS(int nsel = 0, int ReBin = 1, char XTitle[300] = "N_{jets}", char units[300] = "", char plotName[300] = "histo_nice.root", char outputName[300] = "njets",
                  bool isLogY = false, double lumi = 4.6) {
  // nsel == 0 --> WW EWK     SM as signal
  // nsel == 1 --> WW EWK non-SM as signal

  gInterpreter->ExecuteMacro("GoodStyle.C");
  gROOT->LoadMacro("StandardPlotVBSData.C");

  TFile* file = new TFile(plotName, "read");

  StandardPlot myPlot;
  myPlot.setLumi(lumi);
  myPlot.setLabel(XTitle);
  if     (lumi ==    4.9)  myPlot.addLabel("#sqrt{s} = 7 TeV");
  else if(lumi ==   19.4)  myPlot.addLabel("#sqrt{s} = 8 TeV");
  else if(lumi ==   24.4)  myPlot.addLabel("#sqrt{s} = 7+8 TeV");
  else                    myPlot.addLabel(""); 
  myPlot.setUnits(units);

  TH1F* hWWEWK   = (TH1F*)file->Get("histo0");
  TH1F* hVV      = (TH1F*)file->Get("histo1");
  TH1F* hWWQCD   = (TH1F*)file->Get("histo2");
  TH1F* hWJets   = (TH1F*)file->Get("histo3");
  TH1F* hWW      = (TH1F*)file->Get("histo4");
  TH1F* hWWEWKNP = (TH1F*)file->Get("histos");
  TH1F *hData    = (TH1F*)file->Get("histo5");

  assert(hWWEWK);
  assert(hVV);
  assert(hWWQCD);
  assert(hWJets);
  assert(hWW);
  assert(hWWEWKNP);
  assert(hData);

  double scale = 1;
  hWWEWK  ->Scale(scale);
  hVV	  ->Scale(scale);
  hWWQCD  ->Scale(scale);
  hWJets  ->Scale(scale);
  hWW	  ->Scale(scale);
  hWWEWKNP->Scale(scale);

  if     (nsel == 1){
    myPlot.setMCHist(iWWEWKNP,(TH1F*)hWWEWKNP->Clone("hWWEWKNP"));
    myPlot.setMCHist(iWWEWK,  (TH1F*)hWWEWK  ->Clone("hWWEWK"));
    myPlot._mass = 1;
    myPlot.setUnits("");
  }
  else if(nsel == 0 || nsel == 10){
    myPlot.setMCHist(iWWEWKNP,(TH1F*)hWWEWK->Clone("hWWEWK"));
    myPlot._mass = 0;
    if(nsel == 10) myPlot.setHWWOverlaid(true);
    myPlot.setUnits("");
  } else assert(0);

  myPlot.setMCHist(iVV,     (TH1F*)hVV     ->Clone("hVV"));
  myPlot.setMCHist(iWWQCD,  (TH1F*)hWWQCD  ->Clone("hWWQCD"));
  myPlot.setMCHist(iWJets,  (TH1F*)hWJets  ->Clone("hWJets")); 
  myPlot.setMCHist(iWW,     (TH1F*)hWW	   ->Clone("hWW"));
  myPlot.setDataHist((TH1F*)hData->Clone("data"));

  printf("%f + %f + %f + %f + %f = %f - %f - sig: %f\n",
          hWWEWK->GetSumOfWeights(),hVV->GetSumOfWeights(),hWWQCD->GetSumOfWeights(),
  	  hWJets->GetSumOfWeights(),hWW->GetSumOfWeights(),
	  hWWEWK->GetSumOfWeights()+hVV->GetSumOfWeights()+hWWQCD->GetSumOfWeights()+
	  hWJets->GetSumOfWeights()+hWW->GetSumOfWeights(),
	  hData->GetSumOfWeights(),hWWEWKNP->GetSumOfWeights());

  TCanvas* c1 = new TCanvas("c1", "c1",5,50,500,500);

  if(isLogY == true) c1->SetLogy();
  myPlot.Draw(ReBin);  // Can pass a rebin 
  c1->GetFrame()->DrawClone();

  char CommandToExec[300];
  sprintf(CommandToExec,"mkdir -p plots");
  gSystem->Exec(CommandToExec);  

  char myOutputFile[300];
  sprintf(myOutputFile,"plots/%s.eps",outputName);
  //c1->SaveAs(myOutputFile);
  sprintf(myOutputFile,"/afs/cern.ch/user/a/anlevin/www/tmp/%s.png",outputName);
  c1->SaveAs(myOutputFile);
  sprintf(myOutputFile,"plots/%s.pdf",outputName);
  //c1->SaveAs(myOutputFile);

  TCanvas* c2 = new TCanvas("c2", "c2",700,50,500,500);
  c2->cd(1);
  if(nsel == 0 || nsel == 10){
    TH1F* hSignal = (TH1F*)file->Get("histo0");
    TH1F* hSumBck = (TH1F*)file->Get("histo1");
    hSumBck->Add(hWWQCD  );
    hSumBck->Add(hWJets  );
    hSumBck->Add(hWW     );
    printf("S/B(%f/%f) = %f\n",hSignal->GetSumOfWeights(),hSumBck->GetSumOfWeights(),hSignal->GetSumOfWeights()/hSumBck->GetSumOfWeights());
    hSignal->Divide(hSumBck);
    hSignal->Draw("e");
  }
}
