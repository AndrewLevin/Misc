#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1D.h>                   // 1D histograms
#include <TH2D.h>                   // 2D histograms
#include <TLatex.h>                 // 
#include <TGaxis.h>                  // 
#include <TPavesText.h>             // 
#include <TBenchmark.h>             // class to track macro running statistics
#include <TEfficiency.h>            // class to handle efficiency calculations
#include <TMath.h>
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "CPlot.hh"          // helper class for plots
#include "MitStyleRemix.hh"  // style settings for drawing
//#include "MyTools.hh"        // miscellaneous helper functions
#include "CEffUser1D.hh"     // class for handling efficiency graphs
#include "CEffUser2D.hh"     // class for handling efficiency tables

//#include "ZSignals.hh"
//#include "ZBackgrounds.hh"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"


#endif

// RooFit headers
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"

#include "../MitHzz4l/macros/EnergyCorrectionData.hh"

using namespace RooFit;

TH1F *puWeights;

float puweight(float npu) {
  int npuxbin = puWeights->GetXaxis()->FindFixBin(TMath::Min(double(npu), 49.499));
  float puwgt = puWeights->GetBinContent(npuxbin);
  return puwgt; 
}

void escaleFit(int cat=2, TString fileDir = "/data/blue/anlevin/") {

  assert(cat == 1 || cat == 2);


  string label = "";
  if(cat == 1) label = "barrel";
  else if (cat == 2) label == "endcap";
  else
    assert(0);
  
  const TString pileupReweightFile = "/data/blue/vdutta/htt/pileup/PUWeights_S12To2012_2968ipb.root";
  TFile *pufile = new TFile(pileupReweightFile.Data());
  puWeights = (TH1F*)pufile->Get("puWeights");

 RooCategory sample("sample","") ;
  sample.defineType("Data", 1) ;
  sample.defineType("MC", 2) ;


  //gROOT->Macro("MitStyleRemix.cc");

  TLatex* text=new TLatex(3.5,23.5,label.c_str());
  text->SetNDC();
  text->SetTextAlign(13);
  text->SetX(0.8);//(0.940);
  text->SetY(0.95);
  text->SetTextFont(42);
  text->SetTextSize(0.035);// dflt=28

  TString fileNameMC;
  TString fileNameData;

  fileNameMC = "electron_energy_all_s12-zllm50-2-v9.root";
  fileNameData = "electron_energy_all_r12a_and_r12b-del-pr-v1.root";

  TFile* mcFile   = new TFile( (fileDir+fileNameMC).Data() );
  TFile* dataFile = new TFile( (fileDir+fileNameData).Data() );

  TTree* mcTree;
  TTree* dataTree;

  mcTree   = (TTree*) mcFile  ->FindObjectAny("Events");
  dataTree = (TTree*) dataFile->FindObjectAny("Events");

  TString basesel = "pt1>10.0 && pt2>10.0";
  TString isebeb = " && abs(eta1) < 1.479 && abs(eta2) < 1.479";
  TString iseeee = " && abs(eta1) > 1.479 && abs(eta2) > 1.479";

  TString catsel;
  if(cat==1) catsel = isebeb;
  else if (cat==2) catsel = iseeee;

  RooRealVar pt1("pt1","pT_1",0.0,10000.0,"GeV");
  RooRealVar pt2("pt2","pT_2",0.0,10000.0,"GeV");
  RooRealVar eta1("eta1","eta1",-1000.0,1000.0);
  RooRealVar eta2("eta2","eta2",-1000.0,1000.0);
  RooRealVar mass("mass","m_{ee}",60.0,120.0,"GeV");
  RooRealVar npu("npu","npu",0,50);
  RooRealVar wgt("wgt","wgt",-5,1000);

  mass.setRange("zrange",60.0,120.0);

  mass.setBins(120);
  mass.setBins(10000,"cache");

  RooArgSet treeSet(RooArgList(mass,pt1,pt2,eta1,eta2,wgt));

  RooFormulaVar puweightf("puweightv","puweightv","puweight(npu)",RooArgList(npu));

  TH1D* hData = new TH1D("hData","",120,60.,120.);
  hData->GetYaxis()->SetTitle("# Events");
  hData->GetXaxis()->SetTitle("m [GeV]");

  TH1D* hMC = new TH1D("hMC","",120,60.,120.);
  hMC->GetYaxis()->SetTitle("# MC Events");
  hMC->GetXaxis()->SetTitle("m [GeV]");

  dataTree->Draw("mass>>hData", (basesel+catsel).Data());
  hData = (TH1D*) gPad->GetPrimitive("hData");
  RooDataSet* dataDS = new RooDataSet("dataDS","",treeSet,"wgt");

  EnergyCorrectionData electron_energy_data;
  dataTree->SetBranchAddress("Events",&electron_energy_data);
  for(UInt_t i = 0; i < dataTree->GetEntries(); i++){
    dataTree->GetEntry(i);
    if(electron_energy_data.mass < 60 || electron_energy_data.mass > 120) continue;
    if(cat == 1 && (fabs(electron_energy_data.eta1) > 1.479 || fabs(electron_energy_data.eta2) > 1.479)) continue;
    if(cat == 2 && (fabs(electron_energy_data.eta1) < 1.479 || fabs(electron_energy_data.eta2) < 1.479)) continue;
    if(electron_energy_data.pt1 < 10 || electron_energy_data.pt2 < 10) continue;

    mass.setVal(electron_energy_data.mass);
    pt1.setVal(electron_energy_data.pt1);
    pt2.setVal(electron_energy_data.pt2);
    eta1.setVal(electron_energy_data.eta1);
    eta2.setVal(electron_energy_data.eta2);
    //wgt.setVal(1);
    dataDS->add(treeSet,1);    
  }

  RooPlot* frameData = mass.frame(Bins(120),Range("zrange"));
  //RooPlot* frameData = mass.frame(Bins(120),Range(87,97));

  mcTree->Draw("mass>>hMC", (basesel+catsel).Data());
  hMC = (TH1D*) gPad->GetPrimitive("hMC");
  hMC->SetFillColor(kYellow-10);
  RooDataSet* mcDS = new RooDataSet("mcDS","",treeSet,"wgt");

  mcTree->SetBranchAddress("Events",&electron_energy_data);

  cout << "mcDS->isWeighted() = " << mcDS->isWeighted() << endl;
  for(UInt_t i = 0; i < mcTree->GetEntries(); i++){
    mcTree->GetEntry(i);
    if(electron_energy_data.mass < 60 || electron_energy_data.mass > 120) continue;
    if(cat == 1 && ((fabs(electron_energy_data.eta1) > 1.479 || fabs(electron_energy_data.eta2) > 1.479))) continue;
    if(cat == 2 && ((fabs(electron_energy_data.eta1) < 1.479 || fabs(electron_energy_data.eta2) < 1.479))) continue;
    if(electron_energy_data.pt1 < 10 || electron_energy_data.pt2 < 10) continue;

    mass.setVal(electron_energy_data.mass);
    pt1.setVal(electron_energy_data.pt1);
    pt2.setVal(electron_energy_data.pt2);
    eta1.setVal(electron_energy_data.eta1);
    eta2.setVal(electron_energy_data.eta2);
    assert(puweight(electron_energy_data.npu) > -5 && puweight(electron_energy_data.npu) < 1000);
    //wgt.setVal(puweight(electron_energy_data.npu));
    mcDS->add(treeSet,puweight(electron_energy_data.npu));
  }
  //for(UInt_t i = 0; i < mcTree->100; i++){
  //  mcTree->GetEntry(i);
  //  cout << "mcDS->weight() = " << mcDS->weight() << endl;
  //}

  //mcDS->addColumn(w::wgt) ; 

  //dataDS = new RooDataSet("dataDS","",dataDS,treeSet,(basesel+catsel).Data(),"wgt");
  //mcDS = new RooDataSet("mcDS","",mcDS,treeSet,(basesel+catsel).Data(),"wgt");
  cout << "mcDS->numEntries() = " << mcDS->numEntries() << endl;
  cout << "mcDS->sumEntries() = " << mcDS->sumEntries() << endl;
  //return;

  //mcDS->addColumn(puweightf);

  //RooDataSet* mcDSW = new RooDataSet("mcDSW","",*(mcDS->get()),RooFit::Import("MC",*mcDS), RooFit::Index(sample),RooFit::WeightVar("puweightv")); 

  RooPlot* frameMC = mass.frame(Bins(120),Range("zrange"));
  //RooPlot* frameMC = mass.frame(Bins(120),Range(88,98));

  RooRealVar numData("numData","", dataDS->sumEntries(),0,10000000.);
  RooRealVar numMC("numMC","", mcDS->sumEntries(),0,10000000.);

  RooRealVar meanMC("meanMC","",0.0,-5.0,5.0);
  meanMC.removeRange();
    
  RooRealVar sigmaMC("sigmaMC","",1.4,0.0,100.0);
  sigmaMC.removeRange();  
    
  RooRealVar alphaMC("alphaMC","",1.0,0.0,10.0);
  alphaMC.removeRange();
    
  RooRealVar nMC("nMC","",1.0,0.0,1000.0);
  nMC.removeRange();   
    
  const double widthzpdg = 2.4952;
  const double masszpdg = 91.1876;

  //when you don't give a range
  //this automatically makes these roo real vars constant
  RooRealVar mz("mz","mz",masszpdg);
  RooRealVar widthz("wz","wz",widthzpdg);
    
  RooBreitWigner zbwMC("zbwMC","zbwMC",mass,mz,widthz);
  RooCBShape zcbMC("zcbMC","zcbMC",mass,meanMC,sigmaMC,alphaMC,nMC);

  RooFFTConvPdf zcbbwMC("zcbbwMC","zcbbwMC",mass,zbwMC,zcbMC);
  zcbbwMC.setBufferFraction(2.5);

  RooRealVar meanData("meanData","",0.0,-5.0,5.0);
  meanData.removeRange();
    
  RooRealVar sigmaData("sigmaData","",1.4,0.0,100.0);
  sigmaData.removeRange();  
    
  RooRealVar alphaData("alphaData","",1.0,0.0,10.0);
  alphaData.removeRange();

  RooRealVar nData("nData","",1.0,0.0,1000.0);
  nData.removeRange();   
    
  RooBreitWigner zbwData("zbwData","zbwData",mass,mz,widthz);
  RooCBShape zcbData("zcbData","zcbData",mass,meanData,sigmaData,alphaData,nData);

  RooFFTConvPdf zcbbwData("zcbbwData","zcbbwData",mass,zbwData,zcbData);
  zcbbwData.setBufferFraction(2.5);

  RooDataSet combinedDS("combinedDS","combinedDS", RooArgList(mass,wgt), RooFit::Index(sample), RooFit::Import("Data",*dataDS), RooFit::Import("MC",*mcDS),RooFit::WeightVar(wgt));

  //RooDataSet combinedDS("combinedDS","combinedDS", RooArgList(mass,npu), RooFit::Index(sample), RooFit::Import("Data",*dataDS), RooFit::Import("MC",*mcDSW),RooFit::WeightVar("puweightv"));
  //RooDataSet combinedDS("combinedDS","combinedDS", RooArgList(mass,npu),  RooFit::Import("MC",*mcDSW),RooFit::WeightVar("puweightv"));
  
  cout << "mcDS->numEntries() = " << mcDS->numEntries() << endl;
  cout << "mcDS->sumEntries() = " << mcDS->sumEntries() << endl;
  cout << "dataDS->numEntries() = " << dataDS->numEntries() << endl;
  cout << "dataDS->sumEntries() = " << dataDS->sumEntries() << endl;
  //RooDataSet * combinedDS = (RooDataSet*)mcDSW->Clone();
  //combinedDS->append(*dataDS);
  cout << "combinedDS.numEntries() = " << combinedDS.numEntries() << endl;
  cout << "combinedDS.sumEntries() = " << combinedDS.sumEntries() << endl;
  //return;
  // PDF for simultaneous fit
  RooSimultaneous totalPdf("totalPdf","totalPdf", sample);
  totalPdf.addPdf(zcbbwData,"Data");
  totalPdf.addPdf(zcbbwMC,"MC");

  RooFormulaVar delE("delE","(meanData-meanMC)/(mz+meanMC)", RooArgList(meanData,meanMC,mz));
  RooFormulaVar delS2("delS2","2*(sigmaData*sigmaData-sigmaMC*sigmaMC)/(mz+meanMC)/(mz+meanMC)", RooArgList(sigmaData,sigmaMC,mz,meanMC));

  RooFitResult* fitR_mc = zcbbwMC.fitTo(*mcDS, RooFit::Save(),RooFit::Strategy(1),RooFit::SumW2Error(kFALSE));
  
  cout << "alphaMC.getVal() = " << alphaMC.getVal() << endl;
  cout << "nMC.getVal() = " << nMC.getVal() << endl;
  alphaData.setVal(alphaMC.getVal());
  nData.setVal(nMC.getVal());

  alphaData.setConstant();
  nData.setConstant();
  //return;

  //RooFitResult* fitR_data = zcbbwData.fitTo(*dataDS, RooFit::Save(),RooFit::Strategy(2),RooFit::SumW2Error(kFALSE));

  RooFitResult* fitR = totalPdf.fitTo(combinedDS, RooFit::Save(true),RooFit::Strategy(2),RooFit::SumW2Error(false));

  
  TCanvas *czfit = new TCanvas("fitcanvas","",1200,600);
  czfit->Divide(2,1);
  czfit->cd(1);
  mcDS->plotOn(frameMC);
  //totalPdf.plotOn(frameMC,Slice(sample,"MC"),ProjWData(sample,combinedDS),RooFit::LineColor(kBlue+2)) ;

  zcbbwMC.plotOn(frameMC,RooFit::LineColor(kBlue+2));
  frameMC->SetTitle("MC");      
  frameMC->GetYaxis()->SetLabelSize(0.03);
  frameMC->GetYaxis()->SetTitleOffset(1.3);
  frameMC->Draw();    
  TLegend *legmc = new TLegend(0.72,0.75,0.92,0.9);
  legmc->AddEntry(frameMC->getObject(0),"MC","LPE");
  legmc->AddEntry(frameMC->getObject(1),"Fit to MC","L");
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);
  legmc->Draw();
  text->Draw();

  double fitalpha = alphaMC.getVal();
  double fitn = nMC.getVal();

  //alphaData.setVal(fitalpha);
  //nData.setVal(fitn);

  //alphaData.setConstant();
  //nData.setConstant();

  //RooFitResult* fitR_data = zcbbwData.fitTo(*dataDS, RooFit::Save(),RooFit::Strategy(2),RooFit::SumW2Error(kFALSE));

  cout << "meanData.getVal() = " << meanData.getVal() << endl;
  cout << "meanMC.getVal() = " << meanMC.getVal()  << endl;

  czfit->cd(2);
  dataDS->plotOn(frameData);
  //totalPdf.plotOn(frameData,Slice(sample,"Data"),ProjWData(sample,combinedDS),RooFit::LineColor(kRed+2)) ;

  zcbbwData.plotOn(frameData,RooFit::LineColor(kRed+2));  
  frameData->GetYaxis()->SetLabelSize(0.03);
  frameData->GetYaxis()->SetTitle("");
  frameData->SetTitle("Data");
  //zcbbwData.paramOn(frameData,Layout(0.2,0.5,0.9));                       
  //frameData->getAttText()->SetTextSize(0.025);
  frameData->Draw();
  TLegend *legdata = new TLegend(0.62,0.75,0.92,0.9);  
  legdata->AddEntry(frameData->getObject(0),"Data, Full2012","LPE");
  legdata->AddEntry(frameData->getObject(1),"Fit to Data","L");
  legdata->SetBorderSize(0);
  legdata->SetFillStyle(0);
  legdata->Draw();    
  text->Draw();
  TString plotname = TString("zfit") + label + TString("_Full2012.png");
  czfit->SaveAs(plotname);

  double sM = delS2.getVal();

  double cbsigmc = sigmaMC.getVal();
  double cbsigdata = sigmaData.getVal();
  double cbsigmcerr = sigmaMC.getError();
  double cbsigdataerr = sigmaData.getError(); 

  double cbsmear = sqrt(cbsigdata*cbsigdata-cbsigmc*cbsigmc);
  double cbsmearerr = (1.0/cbsmear)*sqrt(cbsigdata*cbsigdata*cbsigdataerr*cbsigdataerr + cbsigmc*cbsigmc*cbsigmcerr*cbsigmcerr);

  double emc = masszpdg+meanMC.getVal();
  double emcerr = meanMC.getError();
  double edataerr = meanData.getError();

  double delEerr = sqrt(edataerr*edataerr*emc*emc + emcerr*emcerr*(emc+1)*(emc+1));
  double sMerr = (1.0/cbsmear*emc*emc)*sqrt(2*(emc*emc*(cbsigdata*cbsigdata*cbsigdataerr*cbsigdataerr + cbsigmc*cbsigmc*cbsigmcerr*cbsigmcerr) + 2*cbsmear*cbsmear*cbsmear*cbsmear*emcerr*emcerr));

  cout << "delE.getVal() = (meanData-meanMC)/(mz+meanMC) = " << delE.getVal() << endl;
  cout << "delE.getPropagatedError(*fitR) = " << delE.getPropagatedError(*fitR) << endl;

  double sME = TMath::Sqrt(TMath::Abs(delS2.getPropagatedError(*fitR)));

  cout << "sM =  abs(2*(sigmaData*sigmaData-sigmaMC*sigmaMC)/(mz+meanMC)/(mz+meanMC)) = " << sM << endl;
  cout << "sME = " << sME << endl;
  
  cout << "cbsmear = sqrt(cbsigdata*cbsigdata-cbsigmc*cbsigmc) = " << cbsmear << endl;
  cout << "cmsmearerr = " << cbsmearerr << endl;

  cout << "alphaData.getVal() = " << alphaData.getVal() << endl;
  cout << "alphaMC.getVal() = " << alphaMC.getVal() << endl;

  cout << "nData.getVal() = " << nData.getVal() << endl;
  cout << "nMC.getVal() = " << nMC.getVal()  << endl;


  cout << "sigmaData.getVal() = " << sigmaData.getVal() << endl;
  cout << "sigmaMC.getVal() = " << sigmaMC.getVal() << endl;

  cout << "meanData.getVal() = " << meanData.getVal() << endl;
  cout << "meanMC.getVal() = " << meanMC.getVal()  << endl;

  cout << "mz.getVal() = " << mz.getVal() << endl;
  cout << "widthz.getVal() = " << widthz.getVal()  << endl;

  
}
