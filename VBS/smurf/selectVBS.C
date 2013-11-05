#include "/home/ceballos/releases/CMSSW_5_2_8/src/Smurf/Core/SmurfTree.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include "TLegend.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSystem.h"

void selectVBS
(
 TString bgdInputFile    = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/backgroundA.root",
 TString signalInputFile = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/zz.root",
 TString dataInputFile   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/zz.root"
 )
{

  float lumi = 20;

  TH1F * mlljj_hist = new TH1F("mlljj","mlljj",100,0,3000);

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile,-1);
  bgdEvent.InitTree(0);

  SmurfTree sigEvent;
  sigEvent.LoadTree(signalInputFile,-1);
  sigEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  int nBgd=bgdEvent.tree_->GetEntries();

  for (int i=0; i<nBgd; ++i) {
    if (i%100000 == 0)
      std::cout << "--- reading event " << i << " out of " << nBgd << std::endl;
    bgdEvent.tree_->GetEntry(i);

    if (!(bgdEvent.cuts_ & SmurfTree::Lep1FullSelection))
      continue;

    if (!(bgdEvent.cuts_ & SmurfTree::Lep2FullSelection))
      continue;

    if (!(bgdEvent.cuts_ & SmurfTree::TopVeto))
      continue;

    if(bgdEvent.jet1_.Pt() < 50)
      continue;

    if(bgdEvent.jet2_.Pt() < 50)
      continue;

    if(bgdEvent.njets_ != 2)
      continue;

    if((bgdEvent.jet1_+bgdEvent.jet2_).M() < 600)
      continue;

    if(bgdEvent.lid1_/abs(bgdEvent.lid1_) != bgdEvent.lid2_/abs(bgdEvent.lid2_))
      continue;

    if(bgdEvent.lep1_.Pt() < 25)
      continue;

    if(bgdEvent.lep2_.Pt() < 25)
      continue;

    mlljj_hist->Fill((bgdEvent.lep1_+bgdEvent.lep2_+bgdEvent.jet1_+bgdEvent.jet2_).M(),bgdEvent.scale1fb_*lumi);

    std::cout << "(bgdEvent.lep1_+bgdEvent.lep2_+bgdEvent.jet1_+bgdEvent.jet2_).M() = " << (bgdEvent.lep1_+bgdEvent.lep2_+bgdEvent.jet1_+bgdEvent.jet2_).M() << std::endl; 
    std::cout << "bgdEvent.scale1fb_ = " << bgdEvent.scale1fb_  << std::endl;
    std::cout << "bgdEvent.dstype_ = " << bgdEvent.dstype_ << std::endl;
  }

  mlljj_hist->Draw();

}
