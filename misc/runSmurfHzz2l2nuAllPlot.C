#if !defined(__CINT__) || defined(__MAKECINT__)
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH2D.h"
#include "SmurfHzzPassSelection.hh"
#include "/home/anlevin/cms/cmssw/021/CMSSW_4_2_3_patch2/src/Smurf/Core/LeptonScaleLookup.h"
#include "/home/anlevin/cms/cmssw/021/CMSSW_4_2_3_patch2/src/Smurf/Core/LeptonScaleLookup.cc"
#include "/home/anlevin/cms/cmssw/021/CMSSW_4_2_3_patch2/src/Smurf/Analysis/HZZllvv/PileupReweighting.h"
#include "TLatex.h"
#include <sstream>
#include <map>
#include "tdrStyle.C"
#endif

void LoopOverEventsAndFillHistograms(SmurfTree * tree, std::string process);
void DrawPlot(std::map<std::string, TH1F*> &MAP, std::string TITLE, std::string X_AXIS_TITLE, std::string OUTPUT_NAME);

bool SmurfV6 = false;
bool plotLog = true;
bool full_sel = false;
double luminosity = 1092;
int mH = 300;
std::string mH_s, luminosity_s;

std::map<std::string,TH1F*> njets_map;
std::map<std::string,TH1F*> minmet_0jets_map;
std::map<std::string,TH1F*> minmet_1jet_map;
std::map<std::string,TH1F*> minmet_2jets_map;
std::map<std::string,TH1F*> mt_0jets_map;
std::map<std::string,TH1F*> mt_1jet_map;
std::map<std::string,TH1F*> mt_2jets_map;
std::map<std::string,TH1F*> zpt_0jets_map;
std::map<std::string,TH1F*> zpt_1jet_map;
std::map<std::string,TH1F*> zpt_2jets_map;

double TopAndWWScaleFactor[3] = {1.248, 1.133, 1.181};

//some of the different processes are grouped togethor in the final histograms
//for example, ggww and qqww are grouped togethor into "ww"
//so the group for ggww and qqww would be "ww"
std::vector<std::string> processGroups; 

//retrieves the value of (data efficiency) / (Monte Carlo efficiency) for an electron or a muon with a given pt and eta
double LeptonEfficiencyScaleFactor(double pt, double eta, TH2D *fhDEffMu, TH2D *fhDEffEl, int lid);

void runSmurfHzz2l2nuAllPlot()
{
  std::stringstream mH_ss; mH_ss << mH; mH_s = mH_ss.str();
  
  std::stringstream luminosity_ss; luminosity_ss << luminosity; luminosity_s = luminosity_ss.str();
  
  processGroups.push_back("ww");
  processGroups.push_back("vz");
  processGroups.push_back("top");
  processGroups.push_back("dytt");
  processGroups.push_back("dyeemm");
  processGroups.push_back("data");
  processGroups.push_back("hzz");
  processGroups.push_back("hww");

  setTDRStyle();

  for(unsigned int i = 0; i < processGroups.size(); i++){
    if(full_sel){
      if(mH == 300){
	njets_map[processGroups[i]]=new TH1F(TString("njets_"+processGroups[i]),TString("njets_"+processGroups[i]),3,0,3);
	minmet_0jets_map[processGroups[i]]=new TH1F(TString("minmet_0jets_"+processGroups[i]),TString("minmet_0jets_"+processGroups[i]),25,60,140);
	minmet_1jet_map[processGroups[i]]=new TH1F(TString("minmet_1jet_"+processGroups[i]),TString("minmet_1jet_"+processGroups[i]),25,90,140);
	minmet_2jets_map[processGroups[i]]=new TH1F(TString("minmet_2jets_"+processGroups[i]),TString("minmet_2jets_"+processGroups[i]),25,90,150);
	mt_0jets_map[processGroups[i]]=new TH1F(TString("mt_0jets_"+processGroups[i]),TString("mt_0jets_"+processGroups[i]),25,100,500);
	mt_1jet_map[processGroups[i]]=new TH1F(TString("mt_1jet_"+processGroups[i]),TString("mt_1jet_"+processGroups[i]),25,100,500);
	mt_2jets_map[processGroups[i]]=new TH1F(TString("mt_2jets_"+processGroups[i]),TString("mt_2jets_"+processGroups[i]),25,100,500);
	zpt_0jets_map[processGroups[i]]=new TH1F(TString("zpt_0jets_"+processGroups[i]),TString("mt_0jets_"+processGroups[i]),25,0,200);
	zpt_1jet_map[processGroups[i]]=new TH1F(TString("zpt_1jet_"+processGroups[i]),TString("mt_1jet_"+processGroups[i]),25,0,200);
	zpt_2jets_map[processGroups[i]]=new TH1F(TString("zpt_2jets_"+processGroups[i]),TString("mt_2jets_"+processGroups[i]),25,0,200);
      }
      else if (mH == 250){
	njets_map[processGroups[i]]=new TH1F(TString("njets_"+processGroups[i]),TString("njets_"+processGroups[i]),3,0,3);
	minmet_0jets_map[processGroups[i]]=new TH1F(TString("minmet_0jets_"+processGroups[i]),TString("minmet_0jets_"+processGroups[i]),25,50,120);
	minmet_1jet_map[processGroups[i]]=new TH1F(TString("minmet_1jet_"+processGroups[i]),TString("minmet_1jet_"+processGroups[i]),25,50,140);
	minmet_2jets_map[processGroups[i]]=new TH1F(TString("minmet_2jets_"+processGroups[i]),TString("minmet_2jets_"+processGroups[i]),25,70,150);
	mt_0jets_map[processGroups[i]]=new TH1F(TString("mt_0jets_"+processGroups[i]),TString("mt_0jets_"+processGroups[i]),25,100,500);
	mt_1jet_map[processGroups[i]]=new TH1F(TString("mt_1jet_"+processGroups[i]),TString("mt_1jet_"+processGroups[i]),25,100,500);
	mt_2jets_map[processGroups[i]]=new TH1F(TString("mt_2jets_"+processGroups[i]),TString("mt_2jets_"+processGroups[i]),25,100,500);
	zpt_0jets_map[processGroups[i]]=new TH1F(TString("zpt_0jets_"+processGroups[i]),TString("mt_0jets_"+processGroups[i]),25,0,200);
	zpt_1jet_map[processGroups[i]]=new TH1F(TString("zpt_1jet_"+processGroups[i]),TString("mt_1jet_"+processGroups[i]),25,0,200);
	zpt_2jets_map[processGroups[i]]=new TH1F(TString("zpt_2jets_"+processGroups[i]),TString("mt_2jets_"+processGroups[i]),25,0,200);
      }
      else if (mH == 400){
	njets_map[processGroups[i]]=new TH1F(TString("njets_"+processGroups[i]),TString("njets_"+processGroups[i]),3,0,3);
	minmet_0jets_map[processGroups[i]]=new TH1F(TString("minmet_0jets_"+processGroups[i]),TString("minmet_0jets_"+processGroups[i]),25,90,220);
	minmet_1jet_map[processGroups[i]]=new TH1F(TString("minmet_1jet_"+processGroups[i]),TString("minmet_1jet_"+processGroups[i]),25,90,200);
	minmet_2jets_map[processGroups[i]]=new TH1F(TString("minmet_2jets_"+processGroups[i]),TString("minmet_2jets_"+processGroups[i]),25,90,230);
	mt_0jets_map[processGroups[i]]=new TH1F(TString("mt_0jets_"+processGroups[i]),TString("mt_0jets_"+processGroups[i]),25,100,500);
	mt_1jet_map[processGroups[i]]=new TH1F(TString("mt_1jet_"+processGroups[i]),TString("mt_1jet_"+processGroups[i]),25,100,500);
	mt_2jets_map[processGroups[i]]=new TH1F(TString("mt_2jets_"+processGroups[i]),TString("mt_2jets_"+processGroups[i]),25,100,500);
	zpt_0jets_map[processGroups[i]]=new TH1F(TString("zpt_0jets_"+processGroups[i]),TString("mt_0jets_"+processGroups[i]),25,0,200);
	zpt_1jet_map[processGroups[i]]=new TH1F(TString("zpt_1jet_"+processGroups[i]),TString("mt_1jet_"+processGroups[i]),25,0,200);
	zpt_2jets_map[processGroups[i]]=new TH1F(TString("zpt_2jets_"+processGroups[i]),TString("mt_2jets_"+processGroups[i]),25,0,200);
      }
      else {assert(0);}
    }
    else{
      if(mH == 250 || mH == 300){
	njets_map[processGroups[i]]=new TH1F(TString("njets_"+processGroups[i]),TString("njets_"+processGroups[i]),3,0,3);
	minmet_0jets_map[processGroups[i]]=new TH1F(TString("minmet_0jets_"+processGroups[i]),TString("minmet_0jets_"+processGroups[i]),30,20,200);
	minmet_1jet_map[processGroups[i]]=new TH1F(TString("minmet_1jet_"+processGroups[i]),TString("minmet_1jet_"+processGroups[i]),30,20,200);
	minmet_2jets_map[processGroups[i]]=new TH1F(TString("minmet_2jets_"+processGroups[i]),TString("minmet_2jets_"+processGroups[i]),30,20,200);
	mt_0jets_map[processGroups[i]]=new TH1F(TString("mt_0jets_"+processGroups[i]),TString("mt_0jets_"+processGroups[i]),25,100,500);
	mt_1jet_map[processGroups[i]]=new TH1F(TString("mt_1jet_"+processGroups[i]),TString("mt_1jet_"+processGroups[i]),25,100,500);
	mt_2jets_map[processGroups[i]]=new TH1F(TString("mt_2jets_"+processGroups[i]),TString("mt_2jets_"+processGroups[i]),25,100,500);
	zpt_0jets_map[processGroups[i]]=new TH1F(TString("zpt_0jets_"+processGroups[i]),TString("mt_0jets_"+processGroups[i]),25,0,200);
	zpt_1jet_map[processGroups[i]]=new TH1F(TString("zpt_1jet_"+processGroups[i]),TString("mt_1jet_"+processGroups[i]),25,0,200);
	zpt_2jets_map[processGroups[i]]=new TH1F(TString("zpt_2jets_"+processGroups[i]),TString("mt_2jets_"+processGroups[i]),25,0,200);
      }
      else if (mH == 400){
	njets_map[processGroups[i]]=new TH1F(TString("njets_"+processGroups[i]),TString("njets_"+processGroups[i]),3,0,3);
	minmet_0jets_map[processGroups[i]]=new TH1F(TString("minmet_0jets_"+processGroups[i]),TString("minmet_0jets_"+processGroups[i]),25,40,220);
	minmet_1jet_map[processGroups[i]]=new TH1F(TString("minmet_1jet_"+processGroups[i]),TString("minmet_1jet_"+processGroups[i]),25,40,190);
	minmet_2jets_map[processGroups[i]]=new TH1F(TString("minmet_2jets_"+processGroups[i]),TString("minmet_2jets_"+processGroups[i]),25,40,150);
	mt_0jets_map[processGroups[i]]=new TH1F(TString("mt_0jets_"+processGroups[i]),TString("mt_0jets_"+processGroups[i]),25,160,450);
	mt_1jet_map[processGroups[i]]=new TH1F(TString("mt_1jet_"+processGroups[i]),TString("mt_1jet_"+processGroups[i]),25,150,450);
	mt_2jets_map[processGroups[i]]=new TH1F(TString("mt_2jets_"+processGroups[i]),TString("mt_2jets_"+processGroups[i]),25,150,450);
	zpt_0jets_map[processGroups[i]]=new TH1F(TString("zpt_0jets_"+processGroups[i]),TString("mt_0jets_"+processGroups[i]),25,0,200);
	zpt_1jet_map[processGroups[i]]=new TH1F(TString("zpt_1jet_"+processGroups[i]),TString("mt_1jet_"+processGroups[i]),25,0,200);
	zpt_2jets_map[processGroups[i]]=new TH1F(TString("zpt_2jets_"+processGroups[i]),TString("mt_2jets_"+processGroups[i]),25,0,200);
      }
    }

  }


  TFile *fZZKFactorFile = TFile::Open("/data/smurf/sixie/KFactors/ZZ_KFactor.root");
  TH1D *ZZKFactor;
  ZZKFactor = (TH1D*)(fZZKFactorFile->Get("KFactorZZ_DileptonPt"));
  if (ZZKFactor) {
    ZZKFactor->SetDirectory(0);
  }
  assert(ZZKFactor);
  fZZKFactorFile->Close();
  delete fZZKFactorFile;

  SmurfTree * tree = new SmurfTree;
  
  std::string directory = "/data/smurf/sixie/data/Run2011_Summer11_"+(SmurfV6?std::string("SmurfV6"):std::string("EPSHZZV0"))+"/mitf-alljets/";

//   tree->LoadTree(TString("/data/smurf/data/Run2011_Spring11_SmurfV6/mitf-alljets/wjets.root"));
//   tree->InitTree(0);
//   LoopOverEventsAndFillHistograms(tree, "wjets");

  tree->LoadTree(TString(directory + "ggww2l.root"));
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "ww");
  
  tree->LoadTree(TString(directory + "ww2l.root"));
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "ww");
  
  tree->LoadTree(TString(directory + "wz3l.root"));
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "wz");
  
  tree->LoadTree(TString(directory + "zz2l.root"));
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "zz");
  
  tree->LoadTree(TString(directory + "stop.root"));
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "top");

  tree->LoadTree(TString(directory + "ttop.root"));
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "top");
  
  tree->LoadTree(TString(directory + "wtop.root"));
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "top");
  
  tree->LoadTree(TString(directory + "ttbar.root"));
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "top");  
  
//   tree->LoadTree(TString("/data/smurf/data/Run2011_Spring11_SmurfV6/mitf-alljets/dyee.root"));
//   tree->InitTree(0);
//   LoopOverEventsAndFillHistograms(tree, "dyeemm");
  
//   tree->LoadTree(TString("/data/smurf/data/Run2011_Spring11_SmurfV6/mitf-alljets/dymm.root"));
//   tree->InitTree(0);
//   LoopOverEventsAndFillHistograms(tree, "dyeemm");

  tree->LoadTree(TString(directory + "data_photons.root"),0);
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "dyeemm");
  
  tree->LoadTree(TString("/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/dytt.root"));
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "dytt");

  tree->LoadTree(TString("/data/smurf/data/Run2011_Spring11_SmurfV6/mitf-alljets/hww" + mH_s + ".root"));
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "hww");

  tree->LoadTree(TString(directory + "gfhzz" + mH_s +".root"));
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "hzz");
  
  if(SmurfV6)
    tree->LoadTree(TString(directory + "data_dileptonAll_1092ipb.goodlumiListEPS.root"));
  else
    tree->LoadTree(TString(directory + "data_2l.goodlumi1092ipb.root"));
  tree->InitTree(0);
  LoopOverEventsAndFillHistograms(tree, "data");

  delete tree;

  for(unsigned int i = 0; i < processGroups.size(); i++){
    Color_t color = 0;
    if(processGroups[i] == "wjets")
      color = kCyan;
    else if (processGroups[i] == "ww")
      color = kYellow+2;
    else if (processGroups[i] == "vz")
      color = kGreen;
    else if (processGroups[i] == "top")
      color = kMagenta;
    else if (processGroups[i] == "dytt")
      color = kBlue;
    else
      continue;
    
    njets_map[processGroups[i]]->SetFillColor(color);
    minmet_0jets_map[processGroups[i]]->SetFillColor(color);
    mt_0jets_map[processGroups[i]]->SetFillColor(color);
    zpt_0jets_map[processGroups[i]]->SetFillColor(color);
    minmet_1jet_map[processGroups[i]]->SetFillColor(color);
    mt_1jet_map[processGroups[i]]->SetFillColor(color);
    zpt_1jet_map[processGroups[i]]->SetFillColor(color);
    minmet_2jets_map[processGroups[i]]->SetFillColor(color);
    mt_2jets_map[processGroups[i]]->SetFillColor(color);
    zpt_2jets_map[processGroups[i]]->SetFillColor(color); 
  }
  
  for(unsigned int i = 0; i < processGroups.size(); i++){
    int line_width = 0;
    if(processGroups[i] == "hww")
      line_width = 3;
    else if (processGroups[i] == "hzz")
      line_width = 3;
    else if (processGroups[i] == "dyeemm")
      line_width = 2;
    else if (processGroups[i] == "data")
      line_width = 2;
    else
      continue;

    njets_map[processGroups[i]]->SetLineWidth(line_width);
    minmet_0jets_map[processGroups[i]]->SetLineWidth(line_width);
    mt_0jets_map[processGroups[i]]->SetLineWidth(line_width);
    zpt_0jets_map[processGroups[i]]->SetLineWidth(line_width);
    minmet_1jet_map[processGroups[i]]->SetLineWidth(line_width);
    mt_1jet_map[processGroups[i]]->SetLineWidth(line_width);
    zpt_1jet_map[processGroups[i]]->SetLineWidth(line_width);
    minmet_2jets_map[processGroups[i]]->SetLineWidth(line_width);
    mt_2jets_map[processGroups[i]]->SetLineWidth(line_width);
    zpt_2jets_map[processGroups[i]]->SetLineWidth(line_width);
  }


//   njets_map["dyeemm"]->SetFillColor(38);
//   minmet_0jets_map["dyeemm"]->SetFillColor(38);
//   mt_0jets_map["dyeemm"]->SetFillColor(38);
//   minmet_1jet_map["dyeemm"]->SetFillColor(38);
//   mt_1jet_map["dyeemm"]->SetFillColor(38);
//   minmet_2jets_map["dyeemm"]->SetFillColor(38);
//   mt_2jets_map["dyeemm"]->SetFillColor(38);
   
  njets_map["dyeemm"]->SetLineColor(kRed);
  minmet_0jets_map["dyeemm"]->SetLineColor(kRed);
  mt_0jets_map["dyeemm"]->SetLineColor(kRed);
  zpt_0jets_map["dyeemm"]->SetLineColor(kRed);
  minmet_1jet_map["dyeemm"]->SetLineColor(kRed);
  mt_1jet_map["dyeemm"]->SetLineColor(kRed);
  zpt_1jet_map["dyeemm"]->SetLineColor(kRed);
  minmet_2jets_map["dyeemm"]->SetLineColor(kRed);
  mt_2jets_map["dyeemm"]->SetLineColor(kRed);
  zpt_2jets_map["dyeemm"]->SetLineColor(kRed);

  njets_map["hww"]->SetLineColor(kBlack);
  minmet_0jets_map["hww"]->SetLineColor(kBlack);
  mt_0jets_map["hww"]->SetLineColor(kBlack);
  zpt_0jets_map["hww"]->SetLineColor(kBlack);
  minmet_1jet_map["hww"]->SetLineColor(kBlack);
  mt_1jet_map["hww"]->SetLineColor(kBlack);
  zpt_1jet_map["hww"]->SetLineColor(kBlack);
  minmet_2jets_map["hww"]->SetLineColor(kBlack);
  mt_2jets_map["hww"]->SetLineColor(kBlack);
  zpt_2jets_map["hww"]->SetLineColor(kBlack);

  njets_map["hww"]->SetLineStyle(2);
  minmet_0jets_map["hww"]->SetLineStyle(2);
  mt_0jets_map["hww"]->SetLineStyle(2);
  zpt_0jets_map["hww"]->SetLineStyle(2);
  minmet_1jet_map["hww"]->SetLineStyle(2);
  mt_1jet_map["hww"]->SetLineStyle(2);
  zpt_1jet_map["hww"]->SetLineStyle(2);
  minmet_2jets_map["hww"]->SetLineStyle(2);
  mt_2jets_map["hww"]->SetLineStyle(2);
  zpt_2jets_map["hww"]->SetLineStyle(2);


  std::string selection = full_sel?"fullselection":"preselection";
  DrawPlot(njets_map,"comparison of data and Monte Carlo", "number of jets", "Hm" + mH_s + "_" + selection + "_njets");
  DrawPlot(minmet_0jets_map,"comparison of data and Monte Carlo in the 0 jet bin", "Min(MET, TrkMET) [GeV]", "Hm" + mH_s + "_" +selection + "_0jets_minmet");
  DrawPlot(minmet_1jet_map,"comparison of data and Monte Carlo in the 1 jet bin", "Min(MET, TrkMET) [GeV]", "Hm" + mH_s + "_" +selection + "_1jet_minmet");
  DrawPlot(minmet_2jets_map,"comparison of data and Monte Carlo in the 2 or more jet bin", "Min(MET, TrkMET) [GeV]", "Hm" + mH_s + "_" +selection + "_2jets_minmet");
  DrawPlot(mt_0jets_map,"comparison of data and Monte Carlo in the 0 jet bin", "transverse mass [GeV]", "Hm" + mH_s + "_" +selection + "_0jets_mt");
  DrawPlot(mt_1jet_map,"comparison of data and Monte Carlo in the 1 jet bin", "transverse mass [GeV]", "Hm" + mH_s + "_" +selection + "_1jet_mt");
  DrawPlot(mt_2jets_map,"comparison of data and Monte Carlo in the 2 or more jet bin", "transverse mass [GeV]", "Hm" + mH_s + "_" +selection + "_2jets_mt");
  DrawPlot(zpt_0jets_map,"comparison of data and Monte Carlo in the 0 jet bin", "Z p_t [GeV]", "Hm" + mH_s + "_" +selection + "_0jets_zpt");
  DrawPlot(zpt_1jet_map,"comparison of data and Monte Carlo in the 1 jet bin", "Z p_t [GeV]", "Hm" + mH_s + "_" +selection + "_1jet_zpt");
  DrawPlot(zpt_2jets_map,"comparison of data and Monte Carlo in the 2 or more jet bin", "Z p_t [GeV]", "Hm" + mH_s + "_" +selection + "_2jets_zpt");

}

void DrawPlot(std::map<std::string, TH1F*> &MAP, std::string TITLE, std::string X_AXIS_TITLE, std::string OUTPUT_NAME)
{
  //add the overflow to the last bin
  for(unsigned int i = 0; i < processGroups.size(); i++){
    MAP[processGroups[i]]->SetBinContent(MAP[processGroups[i]]->GetNbinsX() 
				     ,MAP[processGroups[i]]->GetBinContent(MAP[processGroups[i]]->GetNbinsX())+MAP[processGroups[i]]->GetBinContent(MAP[processGroups[i]]->GetNbinsX()+1) );
  }

  TCanvas * c = new TCanvas;

  TLatex *   tex = new TLatex(0.2,0.84,TString("#sqrt{s}=7 TeV, #int Ldt = "+luminosity_s+" pb^{-1}"));
  tex->SetNDC();
  tex->SetTextSize(0.035);
  tex->SetLineWidth(2);
  TLatex *   tex2 = new TLatex(0.2,0.9,"CMS Preliminary 2011");
  tex2->SetNDC();
  tex2->SetTextSize(0.035);
  tex2->SetLineWidth(2);

  TLegend * leg = new TLegend(0.75,0.57,.95,.92);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  for(unsigned int i = 0; i < processGroups.size(); i++){
    if(processGroups[i]!="hzz" && processGroups[i]!="hww"&&processGroups[i]!="data"&&processGroups[i]!="dyeemm")
      leg->AddEntry(MAP[processGroups[i]],TString(processGroups[i]),"f");
    else if(processGroups[i]!="data"&&processGroups[i]!="dyeemm")
      leg->AddEntry(MAP[processGroups[i]],TString(processGroups[i]+mH_s),"l");
    else if(processGroups[i]!="data")
      leg->AddEntry(MAP[processGroups[i]],TString(processGroups[i]),"l");
    else
      leg->AddEntry(MAP[processGroups[i]],TString(processGroups[i]),"pl");
  }

  THStack * hs = new THStack;
  for(unsigned int i = 0; i < processGroups.size(); i++){
    if(processGroups[i]!="hzz" && processGroups[i]!="hww"&&processGroups[i]!="data")
      hs->Add(MAP[processGroups[i]]);
  }

  double maximum = max(max(max(hs->GetMaximum(),MAP["data"]->GetMaximum()),MAP["hzz"]->GetMaximum()),MAP["hww"]->GetMaximum());
  if(plotLog)
    MAP["data"]->GetYaxis()->SetRangeUser(0.0001,100*maximum);
  else
    MAP["data"]->GetYaxis()->SetRangeUser(0,maximum+ 2*min(maximum,sqrt(maximum)));

  std::stringstream gev_per_bin_ss; gev_per_bin_ss << (MAP["data"]->GetXaxis()->GetXmax() - MAP["data"]->GetXaxis()->GetXmin())/double(MAP["data"]->GetNbinsX());
  std::string gev_per_bin_s = gev_per_bin_ss.str();

  MAP["data"]->GetXaxis()->SetTitle(TString(X_AXIS_TITLE));
  MAP["data"]->GetYaxis()->SetTitle(TString("# of events / "+ gev_per_bin_s+" GeV"));
  MAP["data"]->SetStats(kFALSE);
  //MAP["data"]->SetTitle(TString(TITLE));
  MAP["data"]->SetTitle("");
  MAP["data"]->Draw("E1");
  hs->Draw("SAME");
  MAP["hzz"]->Draw("SAME");
  MAP["hww"]->Draw("SAME");
  MAP["data"]->Draw("SAME E1");
  leg->Draw("SAME");
  tex->Draw();
  tex2->Draw();
  c->SetLogy(plotLog);
  gPad->RedrawAxis();
  c->Print(TString(OUTPUT_NAME + ".pdf"));

}

void LoopOverEventsAndFillHistograms(SmurfTree * tree, std::string process)
{

  TFile *fLeptonEffFile = TFile::Open("/data/smurf/data/EPS/auxiliar/efficiency_results_v6.root");
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  delete fLeptonEffFile;

  TFile *fHiggsPtKFactorFile = TFile::Open("/data/smurf/sixie/KFactors/ggHWW_KFactors_PowhegToHQT.root");
  TH1D *HiggsPtKFactor;
  char kfactorHistName[100];
  sprintf(kfactorHistName, "KFactor_PowhegToHQT_mH%d", mH);
  HiggsPtKFactor = (TH1D*)(fHiggsPtKFactorFile->Get(kfactorHistName));
  if (HiggsPtKFactor) {
    HiggsPtKFactor->SetDirectory(0);
  }
  assert(HiggsPtKFactor);
  fHiggsPtKFactorFile->Close();
  delete fHiggsPtKFactorFile;

  TFile *fZZKFactorFile = TFile::Open("/data/smurf/sixie/KFactors/ZZ_KFactor.root");
  TH1D *ZZKFactor;
  ZZKFactor = (TH1D*)(fZZKFactorFile->Get("KFactorZZ_DileptonPt"));
  if (ZZKFactor) {
    ZZKFactor->SetDirectory(0);
  }
  assert(ZZKFactor);
  fZZKFactorFile->Close();
  delete fZZKFactorFile;

  TFile *fGammaKFactorFile = TFile::Open("/home/anlevin/cms/cmssw/021/CMSSW_4_2_3_patch2/src/Smurf/Analysis/HZZllvv/PhotonSampleReweightFactors.root");
  TH2F *GammaKFactor;
  TH1F *GammaZMass;
  GammaKFactor = (TH2F*)(fGammaKFactorFile->Get("h2_weight"));
  GammaZMass = (TH1F*)(fGammaKFactorFile->Get("h1_mass"));
  if (GammaKFactor) {
    GammaKFactor->SetDirectory(0);
  }   
  assert(GammaKFactor);
  if (GammaZMass) {
    GammaZMass->SetDirectory(0);
  }   
  assert(GammaZMass);
  fGammaKFactorFile->Close();
  delete fGammaKFactorFile;

  LeptonScaleLookup trigLookup("/data/smurf/data/EPS/auxiliar/efficiency_results_v6.root");


  mithep::RunLumiRangeMap rlrm;
 if(process=="dyeemm") rlrm.AddJSONFile("Cert_EPSFINAL_May10ReReco_v2_PromptReco_160404_167913_JSON.txt");

  for(int i = 0; i <tree->tree_->GetEntries(); i++){

    tree->tree_->GetEntry(i);

    if(process=="dyeemm"){

	mithep::RunLumiRangeMap::RunLumiPairType rl(tree->run_, tree->lumi_);
	if(!rlrm.HasRunLumi(rl)) continue; 

	double minmet = min(tree->met_,tree->trackMet_);

	if ( tree->njets_ < 2 && (acos(cos(tree->dilep_.Phi() - tree->jet1_.Phi())) 
				  > (165/180.0)*TMath::Pi() && tree->jet1_.Pt() > 15.0)) continue;
	if (!((tree->cuts_ & SmurfTree::TopVeto) == SmurfTree::TopVeto)) continue;
	if ( tree->dilep_.Pt()    < 40.0 ) continue;
	if (fabs(tree->dilep_.Eta()) > 1.5) continue;
	double mZ = GammaZMass->GetRandom();


	if( minmet < 50.0) continue;

	double pxzll = tree->dilep_.Px() + tree->met_*cos( tree->metPhi_);
	double pyzll = tree->dilep_.Py() + tree->met_*sin( tree->metPhi_);
	double mtH = pow(sqrt(tree->dilep_.Pt()*tree->dilep_.Pt()+mZ*mZ)+
			 sqrt(tree->met_*tree->met_      +mZ*mZ),2)
	  -pxzll*pxzll-pyzll*pyzll;
	if(mtH >0) mtH = sqrt(mtH); else mtH = 0.0;
	
	if(full_sel){ 

	  double metLowerLimit = MetLowerLimit(mH,tree->njets_);
	  double mtHLowerLimit =  MtHLowerLimit(mH,tree->njets_);
	  double mtHUpperLimit =  MtHUpperLimit(mH,tree->njets_);
	  
	  if(minmet < metLowerLimit || mtH < mtHLowerLimit || mtH > mtHUpperLimit) continue;

	  double dijet_mass = (tree->jet1_ + tree->jet2_).M();
	  
	  double minEta = min(tree->jet1_.Eta(),tree->jet2_.Eta());
	  double maxEta = max(tree->jet1_.Eta(),tree->jet2_.Eta());
	  
	  if(tree->njets_ >= 2 && ( tree->jet1_.Pt() < 30 || tree->jet2_.Pt() < 30 || dijet_mass < 450 || maxEta - minEta < 3.5 
				    || (tree->jet3_.Pt() > 30 && tree->jet3_.Eta() < maxEta && tree->jet3_.Eta() > minEta)   )) continue;
	}

	float pt = TMath::Min(299.9, tree->dilep_.Pt());
	int njets = TMath::Min(3, int(tree->njets_));
	int bin_pt        = GammaKFactor->GetXaxis()->FindBin(pt);
	int bin_njet      = GammaKFactor->GetYaxis()->FindBin(njets);
	double weight      = GammaKFactor->GetBinContent(bin_pt, bin_njet);

	njets_map[process]->Fill(tree->njets_,weight);
	if(tree->njets_ == 0){
	  minmet_0jets_map[process]->Fill(minmet,weight);
	  mt_0jets_map[process]->Fill(mtH,weight);
	  zpt_0jets_map[process]->Fill(tree->dilep_.Pt(),weight);
	}
	else if (tree->njets_ == 1){
	  minmet_1jet_map[process]->Fill(minmet,weight);
	  mt_1jet_map[process]->Fill(mtH,weight);
	  zpt_1jet_map[process]->Fill(tree->dilep_.Pt(),weight);
	}
	else if (tree->njets_ >= 2){
	  minmet_2jets_map[process]->Fill(minmet,weight);
	  mt_2jets_map[process]->Fill(mtH,weight);
	  zpt_2jets_map[process]->Fill(tree->dilep_.Pt(),weight);
	}

	continue;
    }
      
    if(!PassSelection(tree,mH,full_sel)) continue;

    if(  (abs(tree->lid1_) == 11 && abs(tree->lid2_) == 13) || (abs(tree->lid1_) == 13 && abs(tree->lid2_) == 11)   ) continue;
    
    double mZ = tree->dilep_.M();
    
    double pxzll = tree->dilep_.Px() + tree->met_*cos( tree->metPhi_);
    double pyzll = tree->dilep_.Py() + tree->met_*sin( tree->metPhi_);
    double mtH = pow(sqrt(tree->dilep_.Pt()*tree->dilep_.Pt()+mZ*mZ)+
		     sqrt(tree->met_*tree->met_      +mZ*mZ),2)
      -pxzll*pxzll-pyzll*pyzll;
    if(mtH >0) mtH = sqrt(mtH); else mtH = 0.0;
    

    double minmet = min(tree->trackMet_,tree->met_);
    
    double selectionEfficiencyScaleFactor1 = LeptonEfficiencyScaleFactor(tree->lep1_.Pt(), tree->lep1_.Eta(),fhDEffMu,fhDEffEl,tree->lid1_);
    double selectionEfficiencyScaleFactor2 = LeptonEfficiencyScaleFactor(tree->lep2_.Pt(), tree->lep2_.Eta(),fhDEffMu,fhDEffEl,tree->lid2_);
    double triggerEfficiencyScaleFactor = trigLookup.GetExpectedTriggerEfficiency(fabs(tree->lep1_.eta()), tree->lep1_.pt() , 
    								 fabs(tree->lep2_.eta()), tree->lep2_.pt(), 
								 TMath::Abs(tree->lid1_), TMath::Abs(tree->lid2_));

    
    double weight = 1;
    
    std::string processGroup = process;

    if(process== "hzz" || process=="hww")
      weight *= HiggsPtKFactor->GetBinContent(HiggsPtKFactor->GetXaxis()->FindFixBin(tree->higgsPt_));
    
    if(tree->dstype_ != 0)
      weight *= triggerEfficiencyScaleFactor*selectionEfficiencyScaleFactor1*selectionEfficiencyScaleFactor2*PileupReweightFactor(tree->nvtx_)*tree->scale1fb_*luminosity/1000;
    
    if(process=="zz")
      weight *= ZZKFactor->GetBinContent( ZZKFactor->GetXaxis()->FindFixBin(tree->dilep_.Pt()));
    
    if(process=="top" || process=="ttbar" || process=="ww" )
      weight *= TopAndWWScaleFactor[tree->njets_<2?tree->njets_:2];
    
    if(process=="wz" || process=="zz")
      processGroup ="vz";
    if(process=="ttbar")
      processGroup = "top";


    njets_map[processGroup]->Fill(tree->njets_,weight);
    if(tree->njets_ == 0){
      minmet_0jets_map[processGroup]->Fill(minmet,weight);
      mt_0jets_map[processGroup]->Fill(mtH,weight);
      zpt_0jets_map[processGroup]->Fill(tree->dilep_.Pt(),weight);
    }
    else if (tree->njets_ == 1){
      minmet_1jet_map[processGroup]->Fill(minmet,weight);
      mt_1jet_map[processGroup]->Fill(mtH,weight);
      zpt_1jet_map[processGroup]->Fill(tree->dilep_.Pt(),weight);
    }
    else if (tree->njets_ >= 2){
      minmet_2jets_map[processGroup]->Fill(minmet,weight);
      mt_2jets_map[processGroup]->Fill(mtH,weight);
      zpt_2jets_map[processGroup]->Fill(tree->dilep_.Pt(),weight);
    }
  }
   
}

double LeptonEfficiencyScaleFactor(double pt, double eta, TH2D *fhDEffMu, TH2D *fhDEffEl, int lid){
  // lid == 13 (muon), 11 (electron)
  double mypt   = TMath::Min(pt,49.999);
  double myeta  = TMath::Min(fabs(eta),2.4999);
  double prob = 1.0;
  if     (TMath::Abs(lid) == 13){
    Int_t ptbin = fhDEffMu->GetXaxis()->FindBin(mypt);
    Int_t etabin = fhDEffMu->GetYaxis()->FindBin(myeta);	 
    prob = fhDEffMu->GetBinContent(ptbin,etabin);
  }
  else if(TMath::Abs(lid) == 11){
    Int_t ptbin = fhDEffEl->GetXaxis()->FindBin(mypt);
    Int_t etabin = fhDEffEl->GetYaxis()->FindBin(myeta);	 
    prob = fhDEffEl->GetBinContent(ptbin,etabin);
  }
  return prob;
}

