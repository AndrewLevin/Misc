#include "/home/ceballos/releases/CMSSW_5_2_8/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/releases/CMSSW_5_2_8/src/Smurf/Analysis/HWWlvlv/factors.h"
#include "/home/ceballos/releases/CMSSW_5_2_8/src/Smurf/Core/LeptonScaleLookup.h"
#include "/home/ceballos/releases/CMSSW_5_2_8/src/Ana/nt_scripts/trilepton.h"
#include "/home/ceballos/releases/CMSSW_5_2_8/src/Smurf/Analysis/HWWlvlv/OtherBkgScaleFactors_8TeV.h"
#include "/home/ceballos/releases/CMSSW_5_2_8/src/Ana/nt_scripts/makeSystematicEffects.h"
#include "/home/ceballos/releases/CMSSW_5_2_8/src/Ana/nt_scripts/LeptonEfficiencyZH.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TLegend.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"

void scaleFactor_WS(LorentzVector l,int q, int ld, int mcld, double val[2]);

const int verboseLevel =   1;
bool UseDyttDataDriven = true; // if true, then remove em events in dyll MC
SmurfTree systEvent;
const unsigned int nSelTypes = 3;
const unsigned int nSelTypesSyst = 7;
const bool showSignalOnly = false;

enum selType {WWSEL, BTAGSEL, WZSEL};
TString selTypeName[nSelTypes*2] = {"WWSEL-OS", "BTAGSEL-OS", "WZSEL-OS",
                                    "WWSEL-SS", "BTAGSEL-SS", "WZSEL-SS"};
enum selTypeSyst {JESUP=0, JESDOWN, LEPP, LEPM, MET, EFFP, EFFM};
TString selTypeNameSyst[nSelTypesSyst*2] = {"JESUP-OS", "JESDOWN-OS", "LEPP-OS", "LEPM-OS", "MET-OS", "EFFP-OS", "EFFM-OS",
                                            "JESUP-SS", "JESDOWN-SS", "LEPP-SS", "LEPM-SS", "MET-SS", "EFFP-SS", "EFFM-SS"};

int weight;

std::string filename="/afs/cern.ch/work/a/anlevin/data/lhe/qed_4_qcd_99_ls0ls1_grid.lhe";

vector<pair<float,float> > grid_points;
vector<float> histo_grid;


void parse_grid()
{
  grid_points.push_back(pair<float,float>(0,0));
  histo_grid.push_back(0);

  ifstream infile(filename.c_str());
  assert(infile.is_open());

  while(!infile.eof()){
    std::string line;
    getline(infile,line);

    if(line=="<initrwgt>\0"){
      getline(infile,line);
      assert(line=="<weightgroup type='mg_reweighting'>");
      while(true){
	getline(infile,line);

	if(line=="</initrwgt>\0")
	  return;

	if (line == "</weight>\0" || line=="</weightgroup>\0")
	  continue;

	int param_number1 = 0;
	int param_number2 = 0;
	float param1 = 0;
	float param2 = 0;

	assert(line.find("set param_card anoinputs") != string::npos);
	std::string paraminfo1=line.substr(line.find("set param_card anoinputs ")+std::string("set param_card anoinputs ").size(),line.find("#")-line.find("set param_card anoinputs ")-std::string("set param_card anoinputs ").size());
	stringstream ss1;
	ss1 << paraminfo1;
	ss1 >> param_number1;
	if(param_number1 == 1)
	  ss1 >> param1;
	else if (param_number1==2)
	  ss1 >> param2;
	else
	  assert(0);

	getline(infile,line);

	if (line != "</weight>\0"){

	  assert(line.find("set param_card anoinputs") != string::npos);
	  std::string paraminfo2=line.substr(line.find("set param_card anoinputs ")+std::string("set param_card anoinputs ").size(),line.find("#")-line.find("set param_card anoinputs ")-std::string("set param_card anoinputs ").size());
	  stringstream ss2;
	  ss2 << paraminfo2;
	  ss2 >> param_number2;
	  if(param_number2 == 1)
	    ss2 >> param1;
	  else if (param_number2==2)
	    ss2 >> param2;
	  else
	    assert(0);

	  assert(param_number1 != param_number2);

	}
	
	grid_points.push_back(pair<float,float>(param1,param2));
	histo_grid.push_back(0);
	
      }
    }
  }

  std::cout << "reweight block not found, exiting" << std::endl;
  exit(1);

}


void vbs_ana_anom
(
 int thePlot = 37,
 TString bgdInputFile    = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/backgroundA_skim14.root",
 TString signalInputFile = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/wwss_ewk_ewkdstype.root",
 TString dataInputFile   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/backgroundA_skim14.root",
 TString systInputFile   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/hww_syst_skim14.root",
 int period = 3,
 int lSel = 4
 )
{

  parse_grid();

  std::cout << "grid_points.size() = " << grid_points.size() << std::endl;
  for(int i = 0; i < grid_points.size(); i++){
    std::cout << grid_points[i].first << ", " << grid_points[i].second << std::endl;
  }
  //change to more convenient units  
  for(int i = 0; i < grid_points.size(); i++){
    grid_points[i].first = grid_points[i].first*pow(10.,10);
    grid_points[i].second = grid_points[i].second*pow(10.,10);
  }
  for(int i = 0; i < grid_points.size(); i++){
    std::cout << grid_points[i].first << ", " << grid_points[i].second << std::endl;
  }

  TH1D * data_obs =  new TH1D("data_obs","data_obs",1,0,1);
  TH1D * diboson= new TH1D("diboson","diboson",1,0,1);
  TH1D * background = new TH1D("background","background",1,0,1);
  TH1D * background_ch1boosted_backshapeUp = new TH1D("background_ch1boosted_backshapeUp","background_ch1boosted_backshapeUp",1,0,1);
  TH1D * background_ch1boosted_backshapeDown = new TH1D("background_ch1boosted_backshapeDown","background_ch1boosted_backshapeDown",1,0,1);
  TH1D * background_ch2boosted_backshapeUp = new TH1D("background_ch2boosted_backshapeUp","background_ch2boosted_backshapeUp",1,0,1);
  TH1D * background_ch2boosted_backshapeDown = new TH1D("background_ch2boosted_backshapeDown","background_ch2boosted_backshapeDown",1,0,1);

  TH2D * th2d  = new TH2D("bin_content_lam_dk_1","bin_content_lam_dk_1",11,-2.75,2.75,11,-2.75,2.75);

  double lumi = 1.0;
  double ptJetMin = 30.0;
  double useFullStatTemplates  = true;

  bool fCheckProblem = true;

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile,-1);
  bgdEvent.InitTree(0);

  SmurfTree sigEvent;
  sigEvent.LoadTree(signalInputFile,-1);
  sigEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  if(systInputFile != ""){
    systEvent.LoadTree(systInputFile,-1);
    systEvent.InitTree(0);
  }

  TString ECMsb  = "";
  TString effPath  = "";
  TString fakePath = "";
  TString puPath   = "";
  unsigned int minRun = 0;
  unsigned int maxRun = 999999;
  double lumiE = 1.099; int year = 1;
  if	 (period == 3){ // Full2012-Summer12-V9-19500ipb
    effPath  = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_Moriond_V1.root";
    fakePath = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/summary_fakes_Moriond2012.root";
    puPath   = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/puWeights_Summer12_53x_True_19p5ifb.root";
    lumi     = 19.365;minRun =      0;maxRun = 999999;ECMsb="8TeV";lumiE = 1.026; year = 2012;
  }
  else if(period == 4){ // Full2011-Fall11-V9
    effPath  = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/efficiency_results_Fall11_SmurfV7_Full2011.root";
    fakePath = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/FakeRates_CutBasedMuon_BDTGWithIPInfoElectron.root";
    puPath   = "/data/smurf/data/Run2011_Fall11_SmurfV9_42X/auxiliar/puWeights_Fall11_42x_True.root";
    lumi     = 4.924;minRun =	 0;maxRun = 999999;ECMsb="7TeV"; lumiE = 1.022; year = 2011;
    UseDyttDataDriven = false;fCheckProblem = false;
  }
  else {
    printf("Wrong period(%d)\n",period);
    return;
  }

  //----------------------------------------------------------------------------
  // radio photon to electron
  //----------------------------------------------------------------------------
  TFile *fRatioPhotonElectron = TFile::Open("/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/ratio_photon_electron.root");
  TH1D *fhDRatioPhotonElectron = (TH1D*)(fRatioPhotonElectron->Get("hDRatioPhotonElectron"));
  assert(fhDRatioPhotonElectron);
  fhDRatioPhotonElectron->SetDirectory(0);
  fRatioPhotonElectron->Close();
  delete fRatioPhotonElectron;

  TFile *fLeptonEffFile = TFile::Open(effPath.Data());
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;

  TFile *fLeptonFRFileM = TFile::Open(fakePath.Data());
  TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold30_PtEta"));
  assert(fhDFRMu);
  fhDFRMu->SetDirectory(0);
  fLeptonFRFileM->Close();
  delete fLeptonFRFileM;

  TFile *fLeptonFRFileE = TFile::Open(fakePath.Data());
  TH2D *fhDFREl = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
  assert(fhDFREl);
  fhDFREl->SetDirectory(0);
  fLeptonFRFileE->Close();
  delete fLeptonFRFileE;

  //Fake rate systematics
  TFile *fLeptonFRFileSyst = TFile::Open(fakePath.Data());
  TH2D *fhDFRMuSyst = (TH2D*)(fLeptonFRFileSyst->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
  TH2D *fhDFRElSyst = (TH2D*)(fLeptonFRFileSyst->Get("ElectronFakeRate_V4_ptThreshold20_PtEta"));
  assert(fhDFRMuSyst);
  assert(fhDFRElSyst);
  fhDFRMuSyst->SetDirectory(0);
  fhDFRElSyst->SetDirectory(0);
  fLeptonFRFileSyst->Close();
  delete fLeptonFRFileSyst;
 
  LeptonScaleLookup trigLookup(effPath.Data());

  // useful if using ZZ lepton selection
  LeptonEfficiencyZH theLeptonEfficiencyZH(year);

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  const int nBinMVA = 1;
  Float_t xbins[nBinMVA+1] = {1100, 3000};
  //if(thePlot == 37) {xbins[0] = 500; xbins[1] = 700; xbins[2] = 1100; xbins[3] = 1600; xbins[4] = 2000;}
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();
  TH1D *histo_Data      = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_WWewk_ALT = (TH1D*) histoMVA->Clone("histo_WWewk_ALT");
  TH1D *histo_WWewk     = (TH1D*) histoMVA->Clone("histo_WWewk");
  TH1D *histo_WWqcd     = (TH1D*) histoMVA->Clone("histo_WWqcd");
  TH1D *histo_WZ        = (TH1D*) histoMVA->Clone("histo_WZ");
  TH1D *histo_WS        = (TH1D*) histoMVA->Clone("histo_WS");
  TH1D *histo_VVV       = (TH1D*) histoMVA->Clone("histo_VVV");
  TH1D *histo_Wjets     = (TH1D*) histoMVA->Clone("histo_Wjets");

  char finalStateName[2],effName[10],momName[10];sprintf(effName,"CMS_eff_l");sprintf(momName,"CMS_p_scale_l");
  if     (lSel == 0) {sprintf(finalStateName,"mm");}
  else if(lSel == 1) {sprintf(finalStateName,"me");}
  else if(lSel == 2) {sprintf(finalStateName,"em");}
  else if(lSel == 3) {sprintf(finalStateName,"mm");}
  else if(lSel == 4) {sprintf(finalStateName,"ll");}
  else if(lSel == 5) {sprintf(finalStateName,"sf");}
  else if(lSel == 6) {sprintf(finalStateName,"of");}
  else {printf("Wrong lSel: %d\n",lSel); assert(0);}

  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 400.0;

  if     (thePlot >=  8 && thePlot <=  8) {nBinPlot = 18;  xminPlot =   0.0; xmaxPlot = 900.0;}
  else if(thePlot >=  0 && thePlot <=  0) {nBinPlot = 40;  xminPlot =   0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 13 && thePlot <= 13) {nBinPlot = 40;  xminPlot =   0.0; xmaxPlot = 400.0;}
  else if(thePlot >=  7 && thePlot <=  7) {nBinPlot = 40; xminPlot =    0.0; xmaxPlot = 400.0;} // mll
  else if(thePlot >=  0 && thePlot <= 14) {}
  else if(thePlot >= 15 && thePlot <= 16) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 1.0;}
  else if(thePlot >= 17 && thePlot <= 17) {nBinPlot =  8; xminPlot = -0.5; xmaxPlot =  7.5;}
  else if(thePlot >= 18 && thePlot <= 18) {nBinPlot = 40; xminPlot = -0.5; xmaxPlot = 39.5;}
  else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 20; xminPlot = 0.0; xmaxPlot = 2000.0;} // mlljj
  else if(thePlot >= 20 && thePlot <= 23) {nBinPlot = 18; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 24 && thePlot <= 24) {nBinPlot = 40; xminPlot = 0.0; xmaxPlot = 2.0;}
  else if(thePlot >= 25 && thePlot <= 29) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 5.0;}
  else if(thePlot >= 30 && thePlot <= 30) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 31 && thePlot <= 32) {nBinPlot = 300; xminPlot = 0.0; xmaxPlot = 600.0;}
  else if(thePlot >= 33 && thePlot <= 33) {nBinPlot = 90; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 34 && thePlot <= 34) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  800.0;}
  else if(thePlot >= 35 && thePlot <= 35) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  8.75;}
  else if(thePlot >= 36 && thePlot <= 36) {nBinPlot = 3; xminPlot = -0.5; xmaxPlot =  2.5;}
  else if(thePlot >= 37 && thePlot <= 37) {nBinPlot = 8; xminPlot = 0.0; xmaxPlot =  2000;} // mjj
  else if(thePlot >= 38 && thePlot <= 38) {nBinPlot = 9; xminPlot = 0.0; xmaxPlot =  8.75;} // detajjs
  else if(thePlot >= 39 && thePlot <= 39) {nBinPlot = 50; xminPlot = 0.0; xmaxPlot =  5.0;}
  else if(thePlot >= 40 && thePlot <= 41) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 45 && thePlot <= 46) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;}
  else if(thePlot >= 46 && thePlot <= 46) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;}
  else if(thePlot >= 47 && thePlot <= 47) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot = 400.0;}
  else if(thePlot >= 48 && thePlot <= 48) {nBinPlot = 20; xminPlot = -0.5; xmaxPlot = 19.5;}
  else if(thePlot >= 49 && thePlot <= 52) {nBinPlot = 300; xminPlot = -15.; xmaxPlot = 15.;}
  else if(thePlot >= 53 && thePlot <= 55) {nBinPlot = 36; xminPlot = 0.0; xmaxPlot = 180.0;}
  else if(thePlot >= 56 && thePlot <= 56) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
  else if(thePlot >= 57 && thePlot <= 57) {nBinPlot = 44; xminPlot = 0.0; xmaxPlot = 4.4;}

  TH1D* histos;
  if(thePlot != 19 && thePlot != 37) histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
  else                               histos = new TH1D("histos", "histos", nBinMVA, xbins);  
  histos->Sumw2();
  TH1D* histo0 = (TH1D*) histos->Clone("histo0");
  TH1D* histo1 = (TH1D*) histos->Clone("histo1");
  TH1D* histo2 = (TH1D*) histos->Clone("histo2");
  TH1D* histo3 = (TH1D*) histos->Clone("histo3");
  TH1D* histo4 = (TH1D*) histos->Clone("histo4");
  TH1D* histo5 = (TH1D*) histos->Clone("histo5");
  histos->Scale(0.0);
  histo0->Scale(0.0);
  histo1->Scale(0.0);
  histo2->Scale(0.0);
  histo3->Scale(0.0);
  histo4->Scale(0.0);
  histo5->Scale(0.0);

  TH1D* histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingUp   = new TH1D( Form("histo_WWewk_ALT_CMS_wwss%s_MVAWWewk_ALTStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WWewk_ALT_CMS_wwss%s_MVAWWewk_ALTStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingUp  ->Sumw2();
  TH1D* histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingDown = new TH1D( Form("histo_WWewk_ALT_CMS_wwss%s_MVAWWewk_ALTStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WWewk_ALT_CMS_wwss%s_MVAWWewk_ALTStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingDown->Sumw2();
  TH1D* histo_WWewk_CMS_MVAWWewkStatBoundingUp           = new TH1D( Form("histo_WWewk_CMS_wwss%s__MVAWWewkStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WWewk_CMS_wwss%s__MVAWWewkStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WWewk_CMS_MVAWWewkStatBoundingUp  ->Sumw2();
  TH1D* histo_WWewk_CMS_MVAWWewkStatBoundingDown         = new TH1D( Form("histo_WWewk_CMS_wwss%s_MVAWWewkStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WWewk_CMS_wwss%s__MVAWWewkStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WWewk_CMS_MVAWWewkStatBoundingDown->Sumw2();
  TH1D* histo_WWqcd_CMS_MVAWWqcdStatBoundingUp           = new TH1D( Form("histo_WWqcd_CMS_wwss%s_MVAWWqcdStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WWqcd_CMS_wwss%s_MVAWWqcdStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WWqcd_CMS_MVAWWqcdStatBoundingUp  ->Sumw2();
  TH1D* histo_WWqcd_CMS_MVAWWqcdStatBoundingDown         = new TH1D( Form("histo_WWqcd_CMS_wwss%s_MVAWWqcdStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WWqcd_CMS_wwss%s_MVAWWqcdStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WWqcd_CMS_MVAWWqcdStatBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVAWSStatBoundingUp                 = new TH1D( Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WS_CMS_MVAWSStatBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVAWSStatBoundingDown               = new TH1D( Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WS_CMS_MVAWSStatBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingUp                 = new TH1D( Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingDown               = new TH1D( Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp               = new TH1D( Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown             = new TH1D( Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_Wjets_CMS_MVAWjetsStatBoundingUp           = new TH1D( Form("histo_Wjets_CMS_wwss%s_MVAWjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Wjets_CMS_wwss%s_MVAWjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Wjets_CMS_MVAWjetsStatBoundingUp  ->Sumw2();
  TH1D* histo_Wjets_CMS_MVAWjetsStatBoundingDown         = new TH1D( Form("histo_Wjets_CMS_wwss%s_MVAWjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Wjets_CMS_wwss%s_MVAWjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Wjets_CMS_MVAWjetsStatBoundingDown->Sumw2();

  TH1D* histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingBinUp  [nBinMVA];
  TH1D* histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingBinDown[nBinMVA];
  TH1D* histo_WWewk_CMS_MVAWWewkStatBoundingBinUp          [nBinMVA];
  TH1D* histo_WWewk_CMS_MVAWWewkStatBoundingBinDown        [nBinMVA];
  TH1D* histo_WWqcd_CMS_MVAWWqcdStatBoundingBinUp          [nBinMVA];
  TH1D* histo_WWqcd_CMS_MVAWWqcdStatBoundingBinDown        [nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinUp                [nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinDown              [nBinMVA];
  TH1D* histo_WS_CMS_MVAWSStatBoundingBinUp                [nBinMVA];
  TH1D* histo_WS_CMS_MVAWSStatBoundingBinDown              [nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinUp              [nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinDown            [nBinMVA];
  TH1D* histo_Wjets_CMS_MVAWjetsStatBoundingBinUp          [nBinMVA];
  TH1D* histo_Wjets_CMS_MVAWjetsStatBoundingBinDown        [nBinMVA];
  for(int nb=0; nb<nBinMVA; nb++){
    histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingBinUp[nb]   = new TH1D(Form("histo_WWewk_ALT_CMS_wwss%s_MVAWWewk_ALTStatBounding_%s_Bin%dUp"	     ,finalStateName,ECMsb.Data(),nb), Form("histo_WWewk_ALT_CMS_wwss%s_MVAWWewk_ALTStatBounding_%s_Bin%dUp"	    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingBinUp[nb]	   ->Sumw2();
    histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingBinDown[nb] = new TH1D(Form("histo_WWewk_ALT_CMS_wwss%s_MVAWWewk_ALTStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WWewk_ALT_CMS_wwss%s_MVAWWewk_ALTStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingBinDown[nb]	  ->Sumw2();
    histo_WWewk_CMS_MVAWWewkStatBoundingBinUp[nb]	    = new TH1D(Form("histo_WWewk_CMS_wwss%s_MVAWWewkStatBounding_%s_Bin%dUp"	    ,finalStateName,ECMsb.Data(),nb), Form("histo_WWewk_CMS_wwss%s_MVAWWewkStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WWewk_CMS_MVAWWewkStatBoundingBinUp[nb]        ->Sumw2();
    histo_WWewk_CMS_MVAWWewkStatBoundingBinDown[nb]	    = new TH1D(Form("histo_WWewk_CMS_wwss%s_MVAWWewkStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WWewk_CMS_wwss%s_MVAWWewkStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WWewk_CMS_MVAWWewkStatBoundingBinDown[nb]	 ->Sumw2();
    histo_WWqcd_CMS_MVAWWqcdStatBoundingBinUp[nb]	    = new TH1D(Form("histo_WWqcd_CMS_wwss%s_MVAWWqcdStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_WWqcd_CMS_wwss%s_MVAWWqcdStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WWqcd_CMS_MVAWWqcdStatBoundingBinUp[nb]	  ->Sumw2();
    histo_WWqcd_CMS_MVAWWqcdStatBoundingBinDown[nb]         = new TH1D(Form("histo_WWqcd_CMS_wwss%s_MVAWWqcdStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WWqcd_CMS_wwss%s_MVAWWqcdStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WWqcd_CMS_MVAWWqcdStatBoundingBinDown[nb]	  ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	            = new TH1D(Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	  ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	            = new TH1D(Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	  ->Sumw2();
    histo_WS_CMS_MVAWSStatBoundingBinUp[nb]	            = new TH1D(Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WS_CMS_MVAWSStatBoundingBinUp[nb]	  ->Sumw2();
    histo_WS_CMS_MVAWSStatBoundingBinDown[nb]	            = new TH1D(Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WS_CMS_MVAWSStatBoundingBinDown[nb]	  ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	            = new TH1D(Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	  ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	            = new TH1D(Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_Wjets_CMS_MVAWjetsStatBoundingBinUp[nb]           = new TH1D(Form("histo_Wjets_CMS_wwss%s_MVAWjetsStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb), Form("histo_Wjets_CMS_wwss%s_MVAWjetsStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Wjets_CMS_MVAWjetsStatBoundingBinUp[nb]  ->Sumw2();
    histo_Wjets_CMS_MVAWjetsStatBoundingBinDown[nb]         = new TH1D(Form("histo_Wjets_CMS_wwss%s_MVAWjetsStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb), Form("histo_Wjets_CMS_wwss%s_MVAWjetsStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Wjets_CMS_MVAWjetsStatBoundingBinDown[nb]->Sumw2();
  }

  TH1D* histo_WWewk_ALT_CMS_MVALepEffBoundingUp    = new TH1D( Form("histo_WWewk_ALT_%sUp",effName)  , Form("histo_WWewk_ALT_%sUp",effName)  , nBinMVA, xbins); histo_WWewk_ALT_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_WWewk_ALT_CMS_MVALepEffBoundingDown  = new TH1D( Form("histo_WWewk_ALT_%sDown",effName), Form("histo_WWewk_ALT_%sDown",effName), nBinMVA, xbins); histo_WWewk_ALT_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_WWewk_CMS_MVALepEffBoundingUp        = new TH1D( Form("histo_WWewk_%sUp",effName)  , Form("histo_WWewk_%sUp",effName)  , nBinMVA, xbins); histo_WWewk_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_WWewk_CMS_MVALepEffBoundingDown      = new TH1D( Form("histo_WWewk_%sDown",effName), Form("histo_WWewk_%sDown",effName), nBinMVA, xbins); histo_WWewk_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_WWqcd_CMS_MVALepEffBoundingUp        = new TH1D( Form("histo_WWqcd_%sUp",effName)  , Form("histo_WWqcd_%sUp",effName)  , nBinMVA, xbins); histo_WWqcd_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_WWqcd_CMS_MVALepEffBoundingDown      = new TH1D( Form("histo_WWqcd_%sDown",effName), Form("histo_WWqcd_%sDown",effName), nBinMVA, xbins); histo_WWqcd_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffBoundingUp           = new TH1D( Form("histo_WZ_%sUp",effName)  , Form("histo_WZ_%sUp",effName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffBoundingDown         = new TH1D( Form("histo_WZ_%sDown",effName), Form("histo_WZ_%sDown",effName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVALepEffBoundingUp           = new TH1D( Form("histo_WS_%sUp",effName)  , Form("histo_WS_%sUp",effName)  , nBinMVA, xbins); histo_WS_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVALepEffBoundingDown         = new TH1D( Form("histo_WS_%sDown",effName), Form("histo_WS_%sDown",effName), nBinMVA, xbins); histo_WS_CMS_MVALepEffBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffBoundingUp          = new TH1D( Form("histo_VVV_%sUp",effName)  , Form("histo_VVV_%sUp",effName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffBoundingDown        = new TH1D( Form("histo_VVV_%sDown",effName), Form("histo_VVV_%sDown",effName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffBoundingDown->Sumw2();

  TH1D* histo_WWewk_ALT_CMS_MVALepResBoundingUp    = new TH1D( Form("histo_WWewk_ALT_%sUp",momName)  , Form("histo_WWewk_ALT_%sUp",momName)  , nBinMVA, xbins); histo_WWewk_ALT_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_WWewk_ALT_CMS_MVALepResBoundingDown  = new TH1D( Form("histo_WWewk_ALT_%sDown",momName), Form("histo_WWewk_ALT_%sDown",momName), nBinMVA, xbins); histo_WWewk_ALT_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_WWewk_CMS_MVALepResBoundingUp        = new TH1D( Form("histo_WWewk_%sUp",momName)  , Form("histo_WWewk_%sUp",momName)  , nBinMVA, xbins); histo_WWewk_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_WWewk_CMS_MVALepResBoundingDown      = new TH1D( Form("histo_WWewk_%sDown",momName), Form("histo_WWewk_%sDown",momName), nBinMVA, xbins); histo_WWewk_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_WWqcd_CMS_MVALepResBoundingUp        = new TH1D( Form("histo_WWqcd_%sUp",momName)  , Form("histo_WWqcd_%sUp",momName)  , nBinMVA, xbins); histo_WWqcd_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_WWqcd_CMS_MVALepResBoundingDown      = new TH1D( Form("histo_WWqcd_%sDown",momName), Form("histo_WWqcd_%sDown",momName), nBinMVA, xbins); histo_WWqcd_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepResBoundingUp           = new TH1D( Form("histo_WZ_%sUp",momName)  , Form("histo_WZ_%sUp",momName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepResBoundingDown         = new TH1D( Form("histo_WZ_%sDown",momName), Form("histo_WZ_%sDown",momName), nBinMVA, xbins); histo_WZ_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVALepResBoundingUp           = new TH1D( Form("histo_WS_%sUp",momName)  , Form("histo_WS_%sUp",momName)  , nBinMVA, xbins); histo_WS_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVALepResBoundingDown         = new TH1D( Form("histo_WS_%sDown",momName), Form("histo_WS_%sDown",momName), nBinMVA, xbins); histo_WS_CMS_MVALepResBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepResBoundingUp          = new TH1D( Form("histo_VVV_%sUp",momName)  , Form("histo_VVV_%sUp",momName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepResBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepResBoundingDown        = new TH1D( Form("histo_VVV_%sDown",momName), Form("histo_VVV_%sDown",momName), nBinMVA, xbins); histo_VVV_CMS_MVALepResBoundingDown->Sumw2();

  TH1D* histo_WWewk_ALT_CMS_MVAMETResBoundingUp    = new TH1D( Form("histo_WWewk_ALT_%sUp","CMS_scale_met")  , Form("histo_WWewk_ALT_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_WWewk_ALT_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_WWewk_ALT_CMS_MVAMETResBoundingDown  = new TH1D( Form("histo_WWewk_ALT_%sDown","CMS_scale_met"), Form("histo_WWewk_ALT_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_WWewk_ALT_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_WWewk_CMS_MVAMETResBoundingUp        = new TH1D( Form("histo_WWewk_%sUp","CMS_scale_met")  , Form("histo_WWewk_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_WWewk_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_WWewk_CMS_MVAMETResBoundingDown      = new TH1D( Form("histo_WWewk_%sDown","CMS_scale_met"), Form("histo_WWewk_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_WWewk_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_WWqcd_CMS_MVAMETResBoundingUp        = new TH1D( Form("histo_WWqcd_%sUp","CMS_scale_met")  , Form("histo_WWqcd_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_WWqcd_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_WWqcd_CMS_MVAMETResBoundingDown      = new TH1D( Form("histo_WWqcd_%sDown","CMS_scale_met"), Form("histo_WWqcd_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_WWqcd_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETResBoundingUp           = new TH1D( Form("histo_WZ_%sUp","CMS_scale_met")  , Form("histo_WZ_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_WZ_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETResBoundingDown         = new TH1D( Form("histo_WZ_%sDown","CMS_scale_met"), Form("histo_WZ_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_WZ_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVAMETResBoundingUp           = new TH1D( Form("histo_WS_%sUp","CMS_scale_met")  , Form("histo_WS_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_WS_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVAMETResBoundingDown         = new TH1D( Form("histo_WS_%sDown","CMS_scale_met"), Form("histo_WS_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_WS_CMS_MVAMETResBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETResBoundingUp          = new TH1D( Form("histo_VVV_%sUp","CMS_scale_met")  , Form("histo_VVV_%sUp","CMS_scale_met")  , nBinMVA, xbins); histo_VVV_CMS_MVAMETResBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETResBoundingDown        = new TH1D( Form("histo_VVV_%sDown","CMS_scale_met"), Form("histo_VVV_%sDown","CMS_scale_met"), nBinMVA, xbins); histo_VVV_CMS_MVAMETResBoundingDown->Sumw2();

  TH1D* histo_WWewk_ALT_CMS_MVAJESBoundingUp    = new TH1D( Form("histo_WWewk_ALT_%sUp","CMS_scale_j")  , Form("histo_WWewk_ALT_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_WWewk_ALT_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WWewk_ALT_CMS_MVAJESBoundingDown  = new TH1D( Form("histo_WWewk_ALT_%sDown","CMS_scale_j"), Form("histo_WWewk_ALT_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_WWewk_ALT_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WWewk_CMS_MVAJESBoundingUp        = new TH1D( Form("histo_WWewk_%sUp","CMS_scale_j")  , Form("histo_WWewk_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_WWewk_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WWewk_CMS_MVAJESBoundingDown      = new TH1D( Form("histo_WWewk_%sDown","CMS_scale_j"), Form("histo_WWewk_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_WWewk_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WWqcd_CMS_MVAJESBoundingUp        = new TH1D( Form("histo_WWqcd_%sUp","CMS_scale_j")  , Form("histo_WWqcd_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_WWqcd_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WWqcd_CMS_MVAJESBoundingDown      = new TH1D( Form("histo_WWqcd_%sDown","CMS_scale_j"), Form("histo_WWqcd_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_WWqcd_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingUp           = new TH1D( Form("histo_WZ_%sUp","CMS_scale_j")  , Form("histo_WZ_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingDown         = new TH1D( Form("histo_WZ_%sDown","CMS_scale_j"), Form("histo_WZ_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVAJESBoundingUp           = new TH1D( Form("histo_WS_%sUp","CMS_scale_j")  , Form("histo_WS_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_WS_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVAJESBoundingDown         = new TH1D( Form("histo_WS_%sDown","CMS_scale_j"), Form("histo_WS_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_WS_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingUp          = new TH1D( Form("histo_VVV_%sUp","CMS_scale_j")  , Form("histo_VVV_%sUp","CMS_scale_j")  , nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingDown        = new TH1D( Form("histo_VVV_%sDown","CMS_scale_j"), Form("histo_VVV_%sDown","CMS_scale_j"), nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();

  TH1D* histo_WZ_CMS_WZNLOBoundingUp        = new TH1D( Form("histo_WZ_CMS_wwss_WZNLOBoundingUp"),   Form("histo_WZ_CMS_wwss_WZNLOBoundingUp"),   nBinMVA, xbins); histo_WZ_CMS_WZNLOBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_WZNLOBoundingDown      = new TH1D( Form("histo_WZ_CMS_wwss_WZNLOBoundingDown"), Form("histo_WZ_CMS_wwss_WZNLOBoundingDown"), nBinMVA, xbins); histo_WZ_CMS_WZNLOBoundingDown->Sumw2();

  TH1D* histo_Wjets_CMS_MVAWBoundingUp      = new TH1D( Form("histo_Wjets_CMS_wwss_MVAWBoundingUp"),   Form("histo_Wjets_CMS_wwss_MVAWBoundingUp"),   nBinMVA, xbins); histo_Wjets_CMS_MVAWBoundingUp  ->Sumw2();
  TH1D* histo_Wjets_CMS_MVAWBoundingDown    = new TH1D( Form("histo_Wjets_CMS_wwss_MVAWBoundingDown"), Form("histo_Wjets_CMS_wwss_MVAWBoundingDown"), nBinMVA, xbins); histo_Wjets_CMS_MVAWBoundingDown->Sumw2();

  TH1D* histo_WS_CMS_MVAWSBoundingUp      = new TH1D( Form("histo_WS_CMS_wwss_MVAWSBoundingUp"),   Form("histo_WS_CMS_wwss_MVAWSBoundingUp"),   nBinMVA, xbins); histo_WS_CMS_MVAWSBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVAWSBoundingDown    = new TH1D( Form("histo_WS_CMS_wwss_MVAWSBoundingDown"), Form("histo_WS_CMS_wwss_MVAWSBoundingDown"), nBinMVA, xbins); histo_WS_CMS_MVAWSBoundingDown->Sumw2();

  double nSelectedData[nSelTypes*2];
  double nSigCut[nSelTypes*2],nSigECut[nSelTypes*2];
  double bgdDecay[nSelTypes*2][45],weiDecay[nSelTypes*2][45];
  double nSigCutSyst[nSelTypesSyst*2],nSigECutSyst[nSelTypesSyst*2];
  double bgdDecaySyst[nSelTypesSyst*2][45],weiDecaySyst[nSelTypesSyst*2][45];
  for(unsigned int i=0; i<nSelTypes*2; i++) {
    nSelectedData[i] = 0.0; nSigCut[i] = 0.0; nSigECut[i] = 0.0;
    for(int j=0; j<45; j++) {
      bgdDecay[i][j] = 0.0; weiDecay[i][j] = 0.0; 
    }
  }
  for(unsigned int i=0; i<nSelTypesSyst*2; i++) {
    nSigCutSyst[i] = 0.0; nSigECutSyst[i] = 0.0;
    for(int j=0; j<45; j++) {
      bgdDecaySyst[i][j] = 0.0; weiDecaySyst[i][j] = 0.0; 
    }
  }

  unsigned int patternTopVeto = SmurfTree::TopVeto;

  int nBgd=bgdEvent.tree_->GetEntries();
  for (int evt=0; evt<nBgd; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nBgd);
    bgdEvent.tree_->GetEntry(evt);

    if(!(((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
         ((bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
	  bgdEvent.dstype_ != SmurfTree::data)) continue;
    if(bgdEvent.dstype_ == SmurfTree::data &&
      (bgdEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ <  minRun) continue;
    if(bgdEvent.dstype_ == SmurfTree::data && bgdEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if     (bgdEvent.dstype_ == SmurfTree::data  	   ) fDecay =  1;
    else if(bgdEvent.dstype_ == SmurfTree::wjets 	   ) fDecay =  3;
    else if(bgdEvent.dstype_ == SmurfTree::ttbar 	   ) fDecay =  5;
    else if(bgdEvent.dstype_ == SmurfTree::dyee  	   ) fDecay =  9;
    else if(bgdEvent.dstype_ == SmurfTree::dymm  	   ) fDecay =  9;
    else if(bgdEvent.dstype_ == SmurfTree::dytt  	   ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven  ) fDecay = 10;
    else if(bgdEvent.dstype_ == SmurfTree::tw    	   ) fDecay = 13;
    else if(bgdEvent.dstype_ == SmurfTree::wgamma	   ) fDecay = 19;
    else if(bgdEvent.dstype_ == SmurfTree::wgstar          ) fDecay = 20;
    else if(bgdEvent.dstype_ == SmurfTree::www             ) fDecay = 21;
    else if(bgdEvent.dstype_ == SmurfTree::wz    	   ) fDecay = 27;
    else if(bgdEvent.dstype_ == SmurfTree::zz    	   ) fDecay = 28;
    else if(bgdEvent.dstype_ == SmurfTree::qqww  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::qqww2j  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::qqwwPWG  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::ggzz  	   ) fDecay = 29;
    else if(bgdEvent.dstype_ == SmurfTree::ggww  	   ) fDecay = 30;
    else if(bgdEvent.dstype_ == SmurfTree::wwewk  	   ) fDecay = 31;
    else if(bgdEvent.dstype_ == SmurfTree::other           ) fDecay = 40;
    else if(bgdEvent.processId_==121 ||
            bgdEvent.processId_==122)   fDecay = 41;
    else if(bgdEvent.processId_==24)    fDecay = 42;
    else if(bgdEvent.processId_==26)    fDecay = 43;
    else if(bgdEvent.processId_==10001) fDecay = 44;
    else if(bgdEvent.processId_==10010) fDecay = 44;
    else                                          {fDecay = 0;std::cout << bgdEvent.dstype_ << std::endl;}

    int correctQ = 0;
    if(bgdEvent.lep1McId_ ==  11 && bgdEvent.lq1_ < 0) correctQ++;
    if(bgdEvent.lep1McId_ == -11 && bgdEvent.lq1_ > 0) correctQ++;
    if(bgdEvent.lep1McId_ ==  13 && bgdEvent.lq1_ < 0) correctQ++;
    if(bgdEvent.lep1McId_ == -13 && bgdEvent.lq1_ > 0) correctQ++;
    if(bgdEvent.lep2McId_ ==  11 && bgdEvent.lq2_ < 0) correctQ++;
    if(bgdEvent.lep2McId_ == -11 && bgdEvent.lq2_ > 0) correctQ++;
    if(bgdEvent.lep2McId_ ==  13 && bgdEvent.lq2_ < 0) correctQ++;
    if(bgdEvent.lep2McId_ == -13 && bgdEvent.lq2_ > 0) correctQ++;

    if(fDecay == 29 && correctQ != 2) fDecay = 30;
    if(fDecay == 27 && correctQ  < 2) fDecay = 30;

    bool passSystCuts[2][nSelTypesSyst-2] = {{false, false, false, false, false},
			                     {false, false, false, false, false}};
    bool passCuts[2][nSelTypes] = {{false, false, false},
                                   {false, false, false}};
    bool isRealLepton = false;
    if((TMath::Abs(bgdEvent.lep1McId_) == 11 || TMath::Abs(bgdEvent.lep1McId_) == 13) &&
       (TMath::Abs(bgdEvent.lep2McId_) == 11 || TMath::Abs(bgdEvent.lep2McId_) == 13)) isRealLepton = true;

    double theMET = bgdEvent.met_; double theMETPHI = bgdEvent.metPhi_; 
    
    int lType = 1;
    if     (bgdEvent.lq1_ * bgdEvent.lq2_ < 0) lType = 0;

    int centrality = 0;
    if(((bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() < 0) ||
        (bgdEvent.jet2_.Eta()-bgdEvent.lep1_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep1_.Eta() < 0)) &&
       ((bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() < 0) ||
        (bgdEvent.jet2_.Eta()-bgdEvent.lep2_.Eta() > 0 && bgdEvent.jet1_.Eta()-bgdEvent.lep2_.Eta() < 0))) centrality = 1; 
    double metMin = 30.0; if(bgdEvent.type_ == SmurfTree::ee) metMin = 40.0;
    if(lType == 0) if(bgdEvent.type_ == SmurfTree::mm) metMin = 40.0;

    bool passNjets         = bgdEvent.njets_ >= 2;
    bool passMET           = TMath::Min(bgdEvent.pmet_,bgdEvent.pTrackMet_) > metMin;
    bool preselCuts        = bgdEvent.lep1_.Pt() > 20. && bgdEvent.lep2_.Pt() > 20.;
    bool passBtagVeto      = (bgdEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool pass3rLVeto       = (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
    bool passMass          = bgdEvent.dilep_.M() > 50.0 && (TMath::Abs(bgdEvent.dilep_.M()-91.1876) > 15 || bgdEvent.type_ != SmurfTree::ee);
    bool passVBFSel        = (bgdEvent.jet1_+bgdEvent.jet2_).M() > 500 && TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta()) > 3.5 && centrality == 1;

    if(lType == 0) passMass = passMass && (TMath::Abs(bgdEvent.dilep_.M()-91.1876) > 15 || bgdEvent.type_ != SmurfTree::mm);

    // 0      1      2       3     4   5      6        7           8  9            10            11     12  13    14
    // lep1pt,lep2pt,dilmass,dilpt,met,metPhi,trackMet,trackMetPhi,mt,dPhiDiLepMET,dPhiMETTrkMET,pTFrac,mtZ,mlljj,mjj;
    double outputVarLepP[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 0, outputVarLepP);
    double outputVarLepM[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 1, outputVarLepM);
    double outputVarMET[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 2, outputVarMET);
    double outputVar[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 3, outputVar);
    double outputVarJESP[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 4, outputVarJESP);
    double outputVarJESM[15];
    makeSystematicEffects(bgdEvent.lid1_, bgdEvent.lid2_, bgdEvent.lep1_, bgdEvent.lep2_, bgdEvent.dilep_, 
                         bgdEvent.mt_, theMET, theMETPHI, 
                         bgdEvent.trackMet_, bgdEvent.trackMetPhi_, 
			 bgdEvent.njets_, bgdEvent.jet1_, bgdEvent.jet2_,
			 year, 5, outputVarJESM);
    double MVAVar[6] = {outputVar[13],outputVarJESP[13],outputVarJESM[13],outputVarLepP[13],outputVarLepM[13],outputVarMET[13]};
    if(thePlot == 37) {MVAVar[0]=outputVar[14];MVAVar[1]=outputVarJESP[14];MVAVar[2]=outputVarJESM[14];MVAVar[3]=outputVarLepP[14];MVAVar[4]=outputVarLepM[14];MVAVar[5]=outputVarMET[14];}
    for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
    double addLepEff	 = 1.0; double addLepEffUp   = 1.0; double addLepEffDown = 1.0;
    addLepEff  = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_, 0)*
    		 leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_, 0);
    if(addLepEff > 0) {
      addLepEffUp   = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_, 1)*
        	      leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_, 1);
      addLepEffDown = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_,-1)*
        	      leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_,-1);
    } else {addLepEff = 1.0;}

    double NjetSyst[2] = {0., 0.};
    if(bgdEvent.jet1_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet2_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet3_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet4_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(bgdEvent.jet1_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet2_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet3_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(bgdEvent.jet4_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;

    bool passLSel = false;
    if     (lSel == 0 && bgdEvent.type_ == SmurfTree::mm) passLSel = true;
    else if(lSel == 1 && bgdEvent.type_ == SmurfTree::me) passLSel = true;
    else if(lSel == 2 && bgdEvent.type_ == SmurfTree::em) passLSel = true;
    else if(lSel == 3 && bgdEvent.type_ == SmurfTree::ee) passLSel = true;
    else if(lSel == 4)                                    passLSel = true;
    else if(lSel == 5 && (bgdEvent.type_ == SmurfTree::mm || bgdEvent.type_ == SmurfTree::ee)) passLSel = true;
    else if(lSel == 6 && (bgdEvent.type_ == SmurfTree::me || bgdEvent.type_ == SmurfTree::em)) passLSel = true;

    if(passLSel && NjetSyst[0] >= 2 	      && TMath::Min(outputVarJESP[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarJESP[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin && outputVar[0]     > 20.0 && outputVar[1]     > 20.0 && passBtagVeto && pass3rLVeto && outputVar[2]     > 50.0 && passVBFSel) passSystCuts[lType][JESUP] = true;
    if(passLSel && NjetSyst[1] >= 2 	      && TMath::Min(outputVarJESM[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarJESM[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin && outputVar[0]     > 20.0 && outputVar[1]     > 20.0 && passBtagVeto && pass3rLVeto && outputVar[2]     > 50.0 && passVBFSel) passSystCuts[lType][JESDOWN] = true;
    if(passLSel && bgdEvent.jet2_.Pt() > ptJetMin && TMath::Min(outputVarLepP[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarLepP[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin && outputVarLepP[0] > 20.0 && outputVarLepP[1] > 20.0 && passBtagVeto && pass3rLVeto && outputVarLepP[2] > 50.0 && passVBFSel) passSystCuts[lType][LEPP] = true;
    if(passLSel && bgdEvent.jet2_.Pt() > ptJetMin && TMath::Min(outputVarLepM[4]/bgdEvent.met_*bgdEvent.pmet_,outputVarLepM[6]/bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin && outputVarLepM[0] > 20.0 && outputVarLepM[1] > 20.0 && passBtagVeto && pass3rLVeto && outputVarLepM[2] > 50.0 && passVBFSel) passSystCuts[lType][LEPM] = true;
    if(passLSel && bgdEvent.jet2_.Pt() > ptJetMin && TMath::Min(outputVarMET[4] /bgdEvent.met_*bgdEvent.pmet_,outputVarMET[6] /bgdEvent.trackMet_*bgdEvent.pTrackMet_) > metMin && outputVarMET[0]  > 20.0 && outputVarMET[1]  > 20.0 && passBtagVeto && pass3rLVeto && outputVarMET[2]  > 50.0 && passVBFSel) passSystCuts[lType][MET] = true;

    if(passNjets  == true && passMET == true &&  passLSel == true &&
       preselCuts == true && bgdEvent.dilep_.M() > 15.0) {
       
       if( passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto) passCuts[lType][WWSEL] = true;
       if(!passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto) passCuts[lType][BTAGSEL] = true;
       if( passBtagVeto && passVBFSel == true                      && !pass3rLVeto && bgdEvent.lep3_.Pt() > 20.) passCuts[lType][WZSEL] = true;

      if(isRealLepton == false &&
         (bgdEvent.dstype_ == SmurfTree::ttbar  || bgdEvent.dstype_ == SmurfTree::tw   || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dymm ||
          bgdEvent.dstype_ == SmurfTree::qqww	|| bgdEvent.dstype_ == SmurfTree::ggww || bgdEvent.dstype_ == SmurfTree::wz   || bgdEvent.dstype_ == SmurfTree::zz   ||
          bgdEvent.dstype_ == SmurfTree::wgstar || bgdEvent.dstype_ == SmurfTree::dytt || bgdEvent.dstype_ == SmurfTree::www)) 
        {for(unsigned int i=0; i<nSelTypes; i++) passCuts[lType][i] = false;}
    }

    if(1){
      double theWeight = 0.0;
      double add       = 1.0;
      int nFake = 0;
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2)  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4) && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(nFake < 0) assert(0);
 
      if(nFake > 1){
	add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
											(bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
        add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
											(bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        add = add*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
											(bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
	fDecay = 22;
	theWeight	       = -1.0*add;
      }
      else if(nFake == 1){
        if(bgdEvent.dstype_ == SmurfTree::data){
	  add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          add = add*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
          if(fCheckProblem == true && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)-add)/add>0.0001)
	    printf("PROBLEMA: %f - %f %f %f %f %f = %f\n",add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);
	  // new category, W+jetsM
	  if((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2){
	    fDecay = 23;
	  }
	  else if((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = add*1.0;
	}
	else if(isRealLepton == true || bgdEvent.dstype_ == SmurfTree::wgamma){
          add = add*fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          add = add*fakeRate(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDFRMu, fhDFREl, (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection,
	                                                                                  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection);
	  add = add*nPUScaleFactor2012(fhDPU ,bgdEvent.npu_);
          add = add*leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	  add = add*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
          if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
          add = add*leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_);

          double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt(), 
								   fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
	        						   TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));

          if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) { 
            double trigEff0 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     	 							      fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
        						             TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
      	    double trigEff1 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     	  							      fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt(), 
      	  							     TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid3_));
      	    double trigEff2 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt() , 
     	  							      fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
      	  							     TMath::Abs( bgdEvent.lid3_), TMath::Abs(bgdEvent.lid2_));
      	    trigEff  = 1.0 - ((1.0-trigEff0)*(1.0-trigEff1)*(1.0-trigEff2));
         }
	  
	  add = add*trigEff;
	  if(fCheckProblem == true && (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)+add)/add>0.0001)
	    printf("PROBLEMB: %f - %f %f %f %f %f = %f\n",add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);
	  fDecay                 = 1;

	  if((bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2 ||
	     (bgdEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2){
	    fDecay = 23;
	  }
	  else if((bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 ||
	  	  (bgdEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = -1.0 * bgdEvent.scale1fb_*lumi*add;
	}
	else {
	  theWeight = 0.0;
	}
      }
      else if(bgdEvent.dstype_ == SmurfTree::dyttDataDriven) {
        double sf_trg = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
	        					        fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
							        TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
        double sf_eff = 1.0;
	sf_eff = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_)*
        	 leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);

        theWeight = ZttScaleFactor(period,bgdEvent.scale1fb_,sf_trg,sf_eff)*lumi;
	if(UseDyttDataDriven == false) theWeight = 0.0;
      }
      else if(bgdEvent.dstype_ != SmurfTree::data){

	double add1 = nPUScaleFactor2012(fhDPU,bgdEvent.npu_);
        double add2 = 1.0;
	add2 = leptonEfficiency(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid1_);
	add2 = add2*leptonEfficiency(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid2_);
        if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
        add2 = add2*leptonEfficiency(bgdEvent.lep3_.Pt(), bgdEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, bgdEvent.lid3_);

        double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
								 fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
	        						 TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));

        if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) { 
           double trigEff0 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     								     fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
                						    TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid2_));
      	   double trigEff1 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep1_.Eta()), bgdEvent.lep1_.Pt() , 
     								     fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt(), 
      								    TMath::Abs( bgdEvent.lid1_), TMath::Abs(bgdEvent.lid3_));
      	   double trigEff2 = trigLookup.GetExpectedTriggerEfficiency(fabs(bgdEvent.lep3_.Eta()), bgdEvent.lep3_.Pt() , 
     								     fabs(bgdEvent.lep2_.Eta()), bgdEvent.lep2_.Pt(), 
      								    TMath::Abs( bgdEvent.lid3_), TMath::Abs(bgdEvent.lid2_));
      	   trigEff  = 1.0 - ((1.0-trigEff0)*(1.0-trigEff1)*(1.0-trigEff2));
        }
        add = add1*add2*trigEff;

        if(fCheckProblem == true && (bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) && add != 0 && TMath::Abs((bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_)-add)/add>0.0001)
	 printf("PROBLEMCB(%d): %f %f %f = %f - %f %f %f %f %f = %f\n",bgdEvent.event_,add1,add2,trigEff,add,bgdEvent.sfWeightFR_,bgdEvent.sfWeightPU_,bgdEvent.sfWeightEff_,bgdEvent.sfWeightTrig_,bgdEvent.sfWeightHPt_,bgdEvent.sfWeightFR_*bgdEvent.sfWeightPU_*bgdEvent.sfWeightEff_*bgdEvent.sfWeightTrig_*bgdEvent.sfWeightHPt_);

	if(bgdEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(bgdEvent.type_,theMET);
        // if true, then remove em events in dyll MC
        if(UseDyttDataDriven == true &&
          (bgdEvent.dstype_ == SmurfTree::dymm || bgdEvent.dstype_ == SmurfTree::dyee || bgdEvent.dstype_ == SmurfTree::dytt) &&
          (bgdEvent.type_ == SmurfTree::em || bgdEvent.type_ == SmurfTree::me)) add = 0.0;

	theWeight              = bgdEvent.scale1fb_*lumi*add;
      }

      // uncertainty related to wrong-sign leptons
      double weightWS[2] = {theWeight,theWeight};
      scaleFactor_WS(bgdEvent.lep1_,bgdEvent.lq1_,bgdEvent.lid1_,bgdEvent.lep1McId_,weightWS);
      scaleFactor_WS(bgdEvent.lep2_,bgdEvent.lq2_,bgdEvent.lid2_,bgdEvent.lep2McId_,weightWS);
      if((bgdEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) { 
        scaleFactor_WS(bgdEvent.lep3_,bgdEvent.lq3_,bgdEvent.lid3_,bgdEvent.lep3McId_,weightWS);
      }
      theWeight = weightWS[0];

      if(passCuts[1][WWSEL]){ // begin making plots 
	double myVar = theMET;
	if     (thePlot == 1) myVar = bgdEvent.lep1_.Pt();
	else if(thePlot == 2) myVar = bgdEvent.lep2_.Pt();
	else if(thePlot == 3) myVar = bgdEvent.lep3_.Pt();
	else if(thePlot == 4) myVar = bgdEvent.jet1_.Pt();
	else if(thePlot == 5) myVar = bgdEvent.jet2_.Pt();
	else if(thePlot == 6) myVar = bgdEvent.jet3_.Pt();
	else if(thePlot == 7) myVar = TMath::Min(bgdEvent.dilep_.M(),399.999);
	else if(thePlot == 8) myVar = bgdEvent.mt_;
	else if(thePlot == 9) myVar = bgdEvent.mt1_;
	else if(thePlot ==10) myVar = bgdEvent.mt2_;
	else if(thePlot ==12) myVar = bgdEvent.trackMet_;
	else if(thePlot ==13) myVar = bgdEvent.dilep_.Pt();
	else if(thePlot ==14) myVar = fabs(bgdEvent.dilep_.M()-91.1876);
	else if(thePlot ==15) myVar = fabs(theMET-bgdEvent.dilep_.Pt())/bgdEvent.dilep_.Pt();
	else if(thePlot ==16) myVar = bgdEvent.lep2_.Pt()/bgdEvent.lep1_.Pt();
	else if(thePlot ==17) myVar = bgdEvent.njets_;
	else if(thePlot ==18) myVar = bgdEvent.nvtx_;
	else if(thePlot ==19) myVar = TMath::Max(TMath::Min((bgdEvent.lep1_+bgdEvent.lep2_+bgdEvent.jet1_+bgdEvent.jet2_).M(),2999.999),700.001);
	else if(thePlot ==20) myVar = bgdEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==21) myVar = TMath::Min(bgdEvent.dPhiLep1MET_,bgdEvent.dPhiLep2MET_)*180.0/TMath::Pi();
	else if(thePlot ==22) myVar = DeltaPhi(bgdEvent.dilep_.Phi() ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==23) myVar = DeltaPhi(bgdEvent.trackMetPhi_ ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==24) myVar = bgdEvent.lep1_.Pt()*bgdEvent.lep2_.Pt()/bgdEvent.jet1_.Pt()/bgdEvent.jet2_.Pt();
	else if(thePlot ==25) myVar = fabs(bgdEvent.dilep_.Eta());
	else if(thePlot ==26) myVar = bgdEvent.dymva_;
	else if(thePlot ==27) myVar = TMath::Min(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
	else if(thePlot ==28) myVar = TMath::Max(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
	else if(thePlot ==29) myVar = TMath::Max(fabs(bgdEvent.jet1_.Eta()),fabs(bgdEvent.jet2_.Eta()));
	else if(thePlot ==30) myVar = DeltaPhi(bgdEvent.dilep_.Phi() ,bgdEvent.jet1_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==37) myVar = TMath::Max(TMath::Min((bgdEvent.jet1_+bgdEvent.jet2_).M(),1999.999),500.001);
	else if(thePlot ==38) myVar = TMath::Abs(bgdEvent.jet1_.Eta()-bgdEvent.jet2_.Eta());
	else if(thePlot ==40) myVar = DeltaPhi(bgdEvent.jet1_.Phi() ,bgdEvent.jet2_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==41) myVar = DeltaPhi(bgdEvent.trackMetPhi_,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==44) myVar = bgdEvent.jet1_.Pt()+ bgdEvent.jet2_.Pt()+bgdEvent.jet3_.Pt();
	else if(thePlot ==48) myVar = bgdEvent.type_;
	else if(thePlot ==49) myVar = theMET*cos(theMETPHI);
	else if(thePlot ==50) myVar = theMET*sin(theMETPHI);
	else if(thePlot ==51) myVar = bgdEvent.trackMet_*cos(bgdEvent.trackMetPhi_);
	else if(thePlot ==52) myVar = bgdEvent.trackMet_*sin(bgdEvent.trackMetPhi_);
	else if(thePlot ==53) myVar = DeltaPhi(bgdEvent.jet3_.Phi(),bgdEvent.jet4_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==55) myVar = bgdEvent.dPhiDiLepMET_*180.0/TMath::Pi();
	else if(thePlot ==57) myVar = bgdEvent.dR_;
      	if     (fDecay == 31){
	  histo0->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 21){
      	  histo1->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 27){
      	  histo1->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 29){
      	  histo2->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 1 || fDecay == 23){
      	  histo3->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 30 || fDecay == 28 ||
                fDecay ==  5 || fDecay == 13 || fDecay == 20 || 
		fDecay == 10 || fDecay ==  9 || fDecay == 19){
      	  histo4->Fill(myVar,theWeight);
      	}
      	else if(fDecay == 41 || fDecay == 42 || fDecay == 43){
      	}
      	else {
      	  printf("NOOOOOOOOOOOOOOOOOOOO %d\n",fDecay);
      	}
      } // end making plots
      for(unsigned int i=0; i<nSelTypes; i++) {
        for(int j=0; j<2; j++){
          if(passCuts[j][i]) {
            bgdDecay[i+j*nSelTypes][(int)fDecay] += theWeight;
            weiDecay[i+j*nSelTypes][(int)fDecay] += theWeight*theWeight;
          }
        }
      }
      for(unsigned int i=0; i<5; i++) {
        for(int j=0; j<2; j++){
          if(passSystCuts[j][i]) {
            bgdDecaySyst[i+j*nSelTypesSyst][(int)fDecay] += theWeight;
            weiDecaySyst[i+j*nSelTypesSyst][(int)fDecay] += theWeight*theWeight;
          }
        }
      }
      if(passCuts[lType][WWSEL]) {
        bgdDecaySyst[EFFP+lType*nSelTypesSyst][(int)fDecay] += theWeight          *addLepEffUp  /addLepEff;
        weiDecaySyst[EFFP+lType*nSelTypesSyst][(int)fDecay] += theWeight*theWeight*addLepEffUp  /addLepEff*addLepEffUp  /addLepEff;
        bgdDecaySyst[EFFM+lType*nSelTypesSyst][(int)fDecay] += theWeight          *addLepEffDown/addLepEff;
        weiDecaySyst[EFFM+lType*nSelTypesSyst][(int)fDecay] += theWeight*theWeight*addLepEffDown/addLepEff*addLepEffDown/addLepEff;
      }

      if     (fDecay == 21){
        if(passCuts[1][WWSEL])  	     histo_VVV  			->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_VVV_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_VVV_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_VVV_CMS_MVAJESBoundingUp	->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_VVV_CMS_MVAJESBoundingDown	->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_VVV_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_VVV_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_VVV_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);;
      }
      else if(fDecay == 31){
        if(passCuts[1][WWSEL])  	     histo_WWewk			  ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_WWewk_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_WWewk_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_WWewk_CMS_MVAJESBoundingUp	  ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_WWewk_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_WWewk_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_WWewk_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_WWewk_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);;
      }
      else if(fDecay == 27){
        if(passCuts[1][WWSEL])  	     histo_WZ			       ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_WZ_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_WZ_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_WZ_CMS_MVAJESBoundingUp     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_WZ_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_WZ_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_WZ_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_WZ_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);;
      }
      else if(fDecay == 29){
        if(passCuts[1][WWSEL])  	     histo_WWqcd		          ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_WWqcd_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_WWqcd_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_WWqcd_CMS_MVAJESBoundingUp	  ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_WWqcd_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_WWqcd_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_WWqcd_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_WWqcd_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);;
      }
      else if(fDecay == 30 || fDecay == 28 ||
              fDecay ==  5 || fDecay == 13 || fDecay == 20 || 
	      fDecay == 10 || fDecay ==  9 || fDecay == 19){
        if(passCuts[1][WWSEL])  	     histo_WS		               ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL])  	     histo_WS_CMS_MVALepEffBoundingUp  ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL])  	     histo_WS_CMS_MVALepEffBoundingDown->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_WS_CMS_MVAJESBoundingUp     ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_WS_CMS_MVAJESBoundingDown   ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_WS_CMS_MVALepResBoundingUp  ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_WS_CMS_MVALepResBoundingDown->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_WS_CMS_MVAMETResBoundingUp  ->Fill(MVAVar[5], theWeight);;
        if(passCuts[1][WWSEL])  	     histo_WS_CMS_MVAWSBoundingUp      ->Fill(MVAVar[0], weightWS[1]);
      }
      else if(fDecay == 1 || fDecay == 23){
        double addFR  =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        												     (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
               addFR  =  addFR*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMu    , fhDFREl    , (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        												     (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        double addFRS =        fakeRate(bgdEvent.lep1_.Pt(), bgdEvent.lep1_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        												     (bgdEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
               addFRS = addFRS*fakeRate(bgdEvent.lep2_.Pt(), bgdEvent.lep2_.Eta(), fhDFRMuSyst, fhDFRElSyst, (bgdEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        												     (bgdEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (bgdEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        if(passCuts[1][WWSEL]) 	     histo_Wjets		       ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL]) 	     histo_Wjets_CMS_MVAWBoundingUp    ->Fill(MVAVar[0], theWeight*addFRS/addFR);
      }
      else assert(0);
    } // if passCuts
  } // end background loop
  
  std::cout << "andrew debug 1" << std::endl;

  //for(int a = 0; a < grid_points.size(); a++){
  //  th2d->SetBinContent(th2d->GetXaxis()->FindFixBin(grid_points[a].first), th2d->GetYaxis()->FindFixBin(grid_points[a].second), 1+(grid_points[a].first-1.5)*(grid_points[a].first-1.5) + (grid_points[a].second-1.5)*(grid_points[a].second-1.5));
  //}

  if(systInputFile != ""){
  int nSyst=systEvent.tree_->GetEntries();
  for (int evt=0; evt<nSyst; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nSyst);
    systEvent.tree_->GetEntry(evt);

    if(systEvent.dstype_ == SmurfTree::data &&
      (systEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(systEvent.dstype_ == SmurfTree::data && systEvent.run_ <  minRun) continue;
    if(systEvent.dstype_ == SmurfTree::data && systEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if     (systEvent.dstype_ == SmurfTree::data  	   ) fDecay =  1;
    else if(systEvent.dstype_ == SmurfTree::wjets 	   ) fDecay =  3;
    else if(systEvent.dstype_ == SmurfTree::ttbar 	   ) fDecay =  5;
    else if(systEvent.dstype_ == SmurfTree::dyee  	   ) fDecay =  9;
    else if(systEvent.dstype_ == SmurfTree::dymm  	   ) fDecay =  9;
    else if(systEvent.dstype_ == SmurfTree::dytt  	   ) fDecay = 10;
    else if(systEvent.dstype_ == SmurfTree::dyttDataDriven ) fDecay = 10;
    else if(systEvent.dstype_ == SmurfTree::tw    	   ) fDecay = 13;
    else if(systEvent.dstype_ == SmurfTree::wgamma	   ) fDecay = 19;
    else if(systEvent.dstype_ == SmurfTree::wgstar         ) fDecay = 20;
    else if(systEvent.dstype_ == SmurfTree::www            ) fDecay = 21;
    else if(systEvent.dstype_ == SmurfTree::wz    	   ) fDecay = 27;
    else if(systEvent.dstype_ == SmurfTree::zz    	   ) fDecay = 28;
    else if(systEvent.dstype_ == SmurfTree::qqww  	   ) fDecay = 29;
    else if(systEvent.dstype_ == SmurfTree::qqwwPWG  	   ) fDecay = 29;
    else if(systEvent.dstype_ == SmurfTree::ggzz  	   ) fDecay = 29;
    else if(systEvent.dstype_ == SmurfTree::ggww  	   ) fDecay = 30;
    else if(systEvent.dstype_ == SmurfTree::other          ) fDecay = 40;
    else if(systEvent.processId_==121 ||
            systEvent.processId_==122)   fDecay = 41;
    else if(systEvent.processId_==24)    fDecay = 42;
    else if(systEvent.processId_==26)    fDecay = 43;
    else if(systEvent.processId_==10001) fDecay = 44;
    else if(systEvent.processId_==10010) fDecay = 44;
    else                                          {fDecay = 0;std::cout << systEvent.dstype_ << std::endl;}

    bool passCuts[3][nSelTypes] = {{false, false, false},
                                   {false, false, false}};
    bool isRealLepton = false;
    if((TMath::Abs(systEvent.lep1McId_) == 11 || TMath::Abs(systEvent.lep1McId_) == 13) &&
       (TMath::Abs(systEvent.lep2McId_) == 11 || TMath::Abs(systEvent.lep2McId_) == 13)) isRealLepton = true;

    double theMET = systEvent.met_; double theMETPHI = systEvent.metPhi_; 
    
    int lType = 1;
    if     (systEvent.lq1_ * systEvent.lq2_ < 0) lType = 0;

    int centrality = 0;
    if(((systEvent.jet1_.Eta()-systEvent.lep1_.Eta() > 0 && systEvent.jet2_.Eta()-systEvent.lep1_.Eta() < 0) ||
        (systEvent.jet2_.Eta()-systEvent.lep1_.Eta() > 0 && systEvent.jet1_.Eta()-systEvent.lep1_.Eta() < 0)) &&
       ((systEvent.jet1_.Eta()-systEvent.lep2_.Eta() > 0 && systEvent.jet2_.Eta()-systEvent.lep2_.Eta() < 0) ||
        (systEvent.jet2_.Eta()-systEvent.lep2_.Eta() > 0 && systEvent.jet1_.Eta()-systEvent.lep2_.Eta() < 0))) centrality = 1; 
    double metMin = 30.0; if(systEvent.type_ == SmurfTree::ee) metMin = 40.0;
    if(lType == 0) if(systEvent.type_ == SmurfTree::mm) metMin = 40.0;

    bool passNjets         = systEvent.njets_ >= 2;
    bool passMET           = TMath::Min(systEvent.pmet_,systEvent.pTrackMet_) > metMin;
    bool preselCuts        = systEvent.lep1_.Pt() > 20. && systEvent.lep2_.Pt() > 20.;
    bool passBtagVeto      = (systEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool pass3rLVeto       = (systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
    bool passMass          = systEvent.dilep_.M() > 50.0 && (TMath::Abs(systEvent.dilep_.M()-91.1876) > 15 || systEvent.type_ != SmurfTree::ee);
    bool passVBFSel        = (systEvent.jet1_+systEvent.jet2_).M() > 500 && TMath::Abs(systEvent.jet1_.Eta()-systEvent.jet2_.Eta()) > 3.5 && centrality == 1;

    if(lType == 0) passMass = passMass && (TMath::Abs(systEvent.dilep_.M()-91.1876) > 15 || systEvent.type_ != SmurfTree::mm);

    bool passLSel = false;
    if     (lSel == 0 && systEvent.type_ == SmurfTree::mm) passLSel = true;
    else if(lSel == 1 && systEvent.type_ == SmurfTree::me) passLSel = true;
    else if(lSel == 2 && systEvent.type_ == SmurfTree::em) passLSel = true;
    else if(lSel == 3 && systEvent.type_ == SmurfTree::ee) passLSel = true;
    else if(lSel == 4)                                    passLSel = true;
    else if(lSel == 5 && (systEvent.type_ == SmurfTree::mm || systEvent.type_ == SmurfTree::ee)) passLSel = true;
    else if(lSel == 6 && (systEvent.type_ == SmurfTree::me || systEvent.type_ == SmurfTree::em)) passLSel = true;

    if(passNjets  == true && passMET == true &&  passLSel == true &&
       preselCuts == true && systEvent.dilep_.M() > 15.0) {
      if( passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto) passCuts[lType][WWSEL] = true;

      if(isRealLepton == false &&
         (systEvent.dstype_ == SmurfTree::ttbar  || systEvent.dstype_ == SmurfTree::tw   || systEvent.dstype_ == SmurfTree::dyee || systEvent.dstype_ == SmurfTree::dymm ||
          systEvent.dstype_ == SmurfTree::qqww   || systEvent.dstype_ == SmurfTree::ggww || systEvent.dstype_ == SmurfTree::wz   || systEvent.dstype_ == SmurfTree::zz   ||
          systEvent.dstype_ == SmurfTree::wgstar || systEvent.dstype_ == SmurfTree::dytt || systEvent.dstype_ == SmurfTree::www)) 
        {for(unsigned int i=0; i<nSelTypes; i++) passCuts[lType][i] = false;}

    }

    if(passCuts[lType][WWSEL]){
      double theWeight = 0.0;
      double add       = 1.0;
      int nFake = 0;
      if(((systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep3LooseMuV2)  == SmurfTree::Lep3LooseMuV2)  && (systEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((systEvent.cuts_ & SmurfTree::Lep3LooseEleV4) == SmurfTree::Lep3LooseEleV4) && (systEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection) nFake++;
      if(nFake < 0) assert(0);
 
      if(nFake > 1){
	add = add*fakeRate(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
											  (systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
        add = add*fakeRate(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
											  (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	fDecay = 22;

	theWeight	       = -1.0*add;
      }
      else if(nFake == 1){
        if(systEvent.dstype_ == SmurfTree::data){
	  add = add*fakeRate(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                    (systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                    (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
          if(fCheckProblem == true && TMath::Abs((systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_)-add)/add>0.0001)
	    printf("PROBLEMA: %f - %f %f %f %f %f = %f\n",add,systEvent.sfWeightFR_,systEvent.sfWeightPU_,systEvent.sfWeightEff_,systEvent.sfWeightTrig_,systEvent.sfWeightHPt_,systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_);
	  // new category, W+jetsM
	  if((systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2){
	    fDecay = 23;
	  }
	  else if((systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = add*1.0;
	}
	else if(isRealLepton == true || systEvent.dstype_ == SmurfTree::wgamma){
          add = add*fakeRate(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
	                                                                                  (systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
          add = add*fakeRate(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDFRMu, fhDFREl, (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
	                                                                                  (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (systEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
	  add = add*nPUScaleFactor2012(fhDPU ,systEvent.npu_);
          add = add*leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_);
	  add = add*leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
          if((systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
          add = add*leptonEfficiency(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid3_);

          double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt(), 
								   fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
	        						   TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid2_));
          add = add*trigEff;
	  if(fCheckProblem == true && (systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto && TMath::Abs((systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_)+add)/add>0.0001)
	    printf("PROBLEMBSyst: %f - %f %f %f %f %f = %f\n",add,systEvent.sfWeightFR_,systEvent.sfWeightPU_,systEvent.sfWeightEff_,systEvent.sfWeightTrig_,systEvent.sfWeightHPt_,systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_);
	  fDecay                 = 1;

	  if((systEvent.cuts_ & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2 ||
	     (systEvent.cuts_ & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2){
	    fDecay = 23;
	  }
	  else if((systEvent.cuts_ & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 ||
	  	  (systEvent.cuts_ & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4){
	  }
	  else {
	    assert(0);
	  }
	  theWeight              = -1.0 * systEvent.scale1fb_*lumi*add;
	}
	else {
	  theWeight = 0.0;
	}
      }
      else if(systEvent.dstype_ == SmurfTree::dyttDataDriven) {
        double sf_trg = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt() , 
	        					        fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
							        TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid2_));
        double sf_eff = 1.0;
	sf_eff = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_)*
        	 leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);

        theWeight = ZttScaleFactor(period,systEvent.scale1fb_,sf_trg,sf_eff)*lumi;
	if(UseDyttDataDriven == false) theWeight = 0.0;
      }
      else if(systEvent.dstype_ != SmurfTree::data){

	double add1 = nPUScaleFactor2012(fhDPU,systEvent.npu_);
        double add2 = 1.0;
	add2 = leptonEfficiency(systEvent.lep1_.Pt(), systEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid1_);
	add2 = add2*leptonEfficiency(systEvent.lep2_.Pt(), systEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid2_);
        if((systEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
        add2 = add2*leptonEfficiency(systEvent.lep3_.Pt(), systEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, systEvent.lid3_);

        double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(systEvent.lep1_.Eta()), systEvent.lep1_.Pt() , 
								 fabs(systEvent.lep2_.Eta()), systEvent.lep2_.Pt(), 
	        						 TMath::Abs( systEvent.lid1_), TMath::Abs(systEvent.lid2_));
        add = add1*add2*trigEff;

        if(fCheckProblem == true && add != 0 && TMath::Abs((systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_)-add)/add>0.0001)
	printf("PROBLEMCSy(%d): %f %f %f = %f - %f %f %f %f %f = %f\n",systEvent.event_,add1,add2,trigEff,add,systEvent.sfWeightFR_,systEvent.sfWeightPU_,systEvent.sfWeightEff_,systEvent.sfWeightTrig_,systEvent.sfWeightHPt_,systEvent.sfWeightFR_*systEvent.sfWeightPU_*systEvent.sfWeightEff_*systEvent.sfWeightTrig_*systEvent.sfWeightHPt_);

	if(systEvent.dstype_ == SmurfTree::wgstar) add = add*WGstarScaleFactor(systEvent.type_,theMET);

        // if true, then remove em events in dyll MC
        if(UseDyttDataDriven == true &&
          (systEvent.dstype_ == SmurfTree::dymm || systEvent.dstype_ == SmurfTree::dyee || systEvent.dstype_ == SmurfTree::dytt) &&
          (systEvent.type_ == SmurfTree::em || systEvent.type_ == SmurfTree::me)) add = 0.0;

        //----------------------------------------------------------------------------      
        // Apply weighting factor to wgamma (gamma->electron ratio)
        //----------------------------------------------------------------------------
        if(systEvent.dstype_ == SmurfTree::wgamma) {
          if(!(TMath::Abs(systEvent.lep1McId_) == 11 || TMath::Abs(systEvent.lep1McId_) == 13)) add = add * ratioPhotonElectron(fhDRatioPhotonElectron,TMath::Abs(systEvent.lep1_.Eta()));
          if(!(TMath::Abs(systEvent.lep2McId_) == 11 || TMath::Abs(systEvent.lep2McId_) == 13)) add = add * ratioPhotonElectron(fhDRatioPhotonElectron,TMath::Abs(systEvent.lep2_.Eta()));      
        }

	theWeight              = systEvent.scale1fb_*lumi*add;
      }

      double outputVar[15];
      makeSystematicEffects(systEvent.lid1_, systEvent.lid2_, systEvent.lep1_, systEvent.lep2_, systEvent.dilep_, 
                            systEvent.mt_, theMET, theMETPHI, 
                            systEvent.trackMet_, systEvent.trackMetPhi_, 
			    systEvent.njets_, systEvent.jet1_, systEvent.jet2_,
			    year, 3, outputVar);
      double MVAVar[6] = {outputVar[13],0,0,0,0,0};
      if(thePlot == 37) {MVAVar[0]=outputVar[14];}
      for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
      if(passCuts[1][WWSEL]){
	if     (fDecay == 27){
	  histo_WZ_CMS_WZNLOBoundingUp->Fill(MVAVar[0], theWeight);
        }
      }

    } // if passCuts
  } // end syst loop
  } // if want to use it at all

  int nSig=sigEvent.tree_->GetEntries();
  for (int evt=0; evt<nSig; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
       printf("--- reading Signal event %5d of %5d\n",evt,nSig);
    sigEvent.tree_->GetEntry(evt);

    //some events did not get reweighted
    if (sigEvent.lheWeights_.size() == 1)
      continue;

    //std::cout << (sigEvent.lheWeights_[a]/sigEvent.lheWeights_[0]) << std::endl;


    bool lId = (sigEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && 
               (sigEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;

    if(!lId) continue;

    bool passSystCuts[2][nSelTypesSyst-2] = {{false, false, false, false, false},
			                     {false, false, false, false, false}};
    bool passCuts[2][nSelTypes] = {{false, false, false},
                                   {false, false, false}};
    bool isRealLepton = false;
    if((TMath::Abs(sigEvent.lep1McId_) == 11 || TMath::Abs(sigEvent.lep1McId_) == 13) &&
       (TMath::Abs(sigEvent.lep2McId_) == 11 || TMath::Abs(sigEvent.lep2McId_) == 13)) isRealLepton = true;

    double theMET = sigEvent.met_; double theMETPHI = sigEvent.metPhi_; 
    
    int lType = 1;
    if     (sigEvent.lq1_ * sigEvent.lq2_ < 0) lType = 0;

    int centrality = 0;
    if(((sigEvent.jet1_.Eta()-sigEvent.lep1_.Eta() > 0 && sigEvent.jet2_.Eta()-sigEvent.lep1_.Eta() < 0) ||
        (sigEvent.jet2_.Eta()-sigEvent.lep1_.Eta() > 0 && sigEvent.jet1_.Eta()-sigEvent.lep1_.Eta() < 0)) &&
       ((sigEvent.jet1_.Eta()-sigEvent.lep2_.Eta() > 0 && sigEvent.jet2_.Eta()-sigEvent.lep2_.Eta() < 0) ||
        (sigEvent.jet2_.Eta()-sigEvent.lep2_.Eta() > 0 && sigEvent.jet1_.Eta()-sigEvent.lep2_.Eta() < 0))) centrality = 1; 
    double metMin = 30.0; if(sigEvent.type_ == SmurfTree::ee) metMin = 40.0;
    if(lType == 0) if(sigEvent.type_ == SmurfTree::mm) metMin = 40.0;

    bool passNjets         = sigEvent.njets_ >= 2;
    bool passMET           = TMath::Min(sigEvent.pmet_,sigEvent.pTrackMet_) > metMin;
    bool preselCuts        = sigEvent.lep1_.Pt() > 20. && sigEvent.lep2_.Pt() > 20.;
    bool passBtagVeto      = (sigEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool pass3rLVeto       = (sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
    bool passMass          = sigEvent.dilep_.M() > 50.0 && (TMath::Abs(sigEvent.dilep_.M()-91.1876) > 15 || sigEvent.type_ != SmurfTree::ee);
    bool passVBFSel        = (sigEvent.jet1_+sigEvent.jet2_).M() > 500 && TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta()) > 3.5 && centrality == 1;

    if(lType == 0) passMass = passMass && (TMath::Abs(sigEvent.dilep_.M()-91.1876) > 15 || sigEvent.type_ != SmurfTree::mm);

    double outputVarLepP[15];
    makeSystematicEffects(sigEvent.lid1_, sigEvent.lid2_, sigEvent.lep1_, sigEvent.lep2_, sigEvent.dilep_, 
                         sigEvent.mt_, theMET, theMETPHI, 
                         sigEvent.trackMet_, sigEvent.trackMetPhi_, 
			 sigEvent.njets_, sigEvent.jet1_, sigEvent.jet2_,
			 year, 0, outputVarLepP);
    double outputVarLepM[15];
    makeSystematicEffects(sigEvent.lid1_, sigEvent.lid2_, sigEvent.lep1_, sigEvent.lep2_, sigEvent.dilep_, 
                         sigEvent.mt_, theMET, theMETPHI, 
                         sigEvent.trackMet_, sigEvent.trackMetPhi_, 
			 sigEvent.njets_, sigEvent.jet1_, sigEvent.jet2_,
			 year, 1, outputVarLepM);
    double outputVarMET[15];
    makeSystematicEffects(sigEvent.lid1_, sigEvent.lid2_, sigEvent.lep1_, sigEvent.lep2_, sigEvent.dilep_, 
                         sigEvent.mt_, theMET, theMETPHI, 
                         sigEvent.trackMet_, sigEvent.trackMetPhi_, 
			 sigEvent.njets_, sigEvent.jet1_, sigEvent.jet2_,
			 year, 2, outputVarMET);
    double outputVar[15];
    makeSystematicEffects(sigEvent.lid1_, sigEvent.lid2_, sigEvent.lep1_, sigEvent.lep2_, sigEvent.dilep_, 
                         sigEvent.mt_, theMET, theMETPHI, 
                         sigEvent.trackMet_, sigEvent.trackMetPhi_, 
			 sigEvent.njets_, sigEvent.jet1_, sigEvent.jet2_,
			 year, 3, outputVar);
    double outputVarJESP[15];
    makeSystematicEffects(sigEvent.lid1_, sigEvent.lid2_, sigEvent.lep1_, sigEvent.lep2_, sigEvent.dilep_, 
                         sigEvent.mt_, theMET, theMETPHI, 
                         sigEvent.trackMet_, sigEvent.trackMetPhi_, 
			 sigEvent.njets_, sigEvent.jet1_, sigEvent.jet2_,
			 year, 4, outputVarJESP);
    double outputVarJESM[15];
    makeSystematicEffects(sigEvent.lid1_, sigEvent.lid2_, sigEvent.lep1_, sigEvent.lep2_, sigEvent.dilep_, 
                         sigEvent.mt_, theMET, theMETPHI, 
                         sigEvent.trackMet_, sigEvent.trackMetPhi_, 
			 sigEvent.njets_, sigEvent.jet1_, sigEvent.jet2_,
			 year, 5, outputVarJESM);
    double MVAVar[6] = {outputVar[13],outputVarJESP[13],outputVarJESM[13],outputVarLepP[13],outputVarLepM[13],outputVarMET[13]};
    if(thePlot == 37) {MVAVar[0]=outputVar[14];MVAVar[1]=outputVarJESP[14];MVAVar[2]=outputVarJESM[14];MVAVar[3]=outputVarLepP[14];MVAVar[4]=outputVarLepM[14];MVAVar[5]=outputVarMET[14];}
    for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
    double addLepEff	 = 1.0; double addLepEffUp   = 1.0; double addLepEffDown = 1.0;
    addLepEff  = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_, 0)*
    		 leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_, 0);
    if(addLepEff > 0) {
      addLepEffUp   = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_, 1)*
        	      leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_, 1);
      addLepEffDown = leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_,-1)*
        	      leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_,-1);
    } else {addLepEff = 1.0;}

    double NjetSyst[2] = {0., 0.};
    if(sigEvent.jet1_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(sigEvent.jet2_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(sigEvent.jet3_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(sigEvent.jet4_.Pt()*1.05 > ptJetMin) NjetSyst[0]++;
    if(sigEvent.jet1_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(sigEvent.jet2_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(sigEvent.jet3_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;
    if(sigEvent.jet4_.Pt()*0.95 > ptJetMin) NjetSyst[1]++;

    bool passLSel = false;
    if     (lSel == 0 && sigEvent.type_ == SmurfTree::mm) passLSel = true;
    else if(lSel == 1 && sigEvent.type_ == SmurfTree::me) passLSel = true;
    else if(lSel == 2 && sigEvent.type_ == SmurfTree::em) passLSel = true;
    else if(lSel == 3 && sigEvent.type_ == SmurfTree::ee) passLSel = true;
    else if(lSel == 4)                                    passLSel = true;
    else if(lSel == 5 && (sigEvent.type_ == SmurfTree::mm || sigEvent.type_ == SmurfTree::ee)) passLSel = true;
    else if(lSel == 6 && (sigEvent.type_ == SmurfTree::me || sigEvent.type_ == SmurfTree::em)) passLSel = true;

    if(passLSel && NjetSyst[0] >= 2 	      && TMath::Min(outputVarJESP[4]/sigEvent.met_*sigEvent.pmet_,outputVarJESP[6]/sigEvent.trackMet_*sigEvent.pTrackMet_) > metMin && outputVar[0]	> 20.0 &&  outputVar[1]     > 20.0 && passBtagVeto && pass3rLVeto && outputVar[2]     > 50.0 && passVBFSel) passSystCuts[lType][JESUP] = true;
    if(passLSel && NjetSyst[1] >= 2 	      && TMath::Min(outputVarJESM[4]/sigEvent.met_*sigEvent.pmet_,outputVarJESM[6]/sigEvent.trackMet_*sigEvent.pTrackMet_) > metMin && outputVar[0]	> 20.0 &&  outputVar[1]     > 20.0 && passBtagVeto && pass3rLVeto && outputVar[2]     > 50.0 && passVBFSel) passSystCuts[lType][JESDOWN] = true;
    if(passLSel && sigEvent.jet2_.Pt() > ptJetMin && TMath::Min(outputVarLepP[4]/sigEvent.met_*sigEvent.pmet_,outputVarLepP[6]/sigEvent.trackMet_*sigEvent.pTrackMet_) > metMin && outputVarLepP[0] > 20.0 &&  outputVarLepP[1] > 20.0 && passBtagVeto && pass3rLVeto && outputVarLepP[2] > 50.0 && passVBFSel) passSystCuts[lType][LEPP] = true;
    if(passLSel && sigEvent.jet2_.Pt() > ptJetMin && TMath::Min(outputVarLepM[4]/sigEvent.met_*sigEvent.pmet_,outputVarLepM[6]/sigEvent.trackMet_*sigEvent.pTrackMet_) > metMin && outputVarLepM[0] > 20.0 &&  outputVarLepM[1] > 20.0 && passBtagVeto && pass3rLVeto && outputVarLepM[2] > 50.0 && passVBFSel) passSystCuts[lType][LEPM] = true;
    if(passLSel && sigEvent.jet2_.Pt() > ptJetMin && TMath::Min(outputVarMET[4] /sigEvent.met_*sigEvent.pmet_,outputVarMET[6] /sigEvent.trackMet_*sigEvent.pTrackMet_) > metMin && outputVarMET[0]  > 20.0 &&  outputVarMET[1]  > 20.0 && passBtagVeto && pass3rLVeto && outputVarMET[2]  > 50.0 && passVBFSel) passSystCuts[lType][MET] = true;

    if(passNjets  == true && passMET == true &&  passLSel == true &&
       preselCuts == true && sigEvent.dilep_.M() > 15.0) {
       
       if( passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto) passCuts[lType][WWSEL] = true;
       if(!passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto) passCuts[lType][BTAGSEL] = true;
       if( passBtagVeto && passVBFSel == true                      && !pass3rLVeto && sigEvent.lep3_.Pt() > 20.) passCuts[lType][WZSEL] = true;

      if(isRealLepton == false) {for(unsigned int i=0; i<nSelTypes; i++) passCuts[lType][i] = false;}
    }

    if(isRealLepton == true){
      double add = 1.;
      double addPU = 1.;
      addPU = nPUScaleFactor2012(fhDPU,sigEvent.npu_);
      add = add*leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_);
      add = add*leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_);
      if((sigEvent.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto)
      add = add*leptonEfficiency(sigEvent.lep3_.Pt(), sigEvent.lep3_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid3_);

      double trigEff = trigLookup.GetExpectedTriggerEfficiency(fabs(sigEvent.lep1_.Eta()), sigEvent.lep1_.Pt() , 
  							       fabs(sigEvent.lep2_.Eta()), sigEvent.lep2_.Pt(), 
        						       TMath::Abs( sigEvent.lid1_), TMath::Abs(sigEvent.lid2_));
      add = add*trigEff*addPU;

      double theWeight = sigEvent.scale1fb_*lumi*add;

      if(fCheckProblem == true && TMath::Abs((sigEvent.sfWeightFR_*sigEvent.sfWeightPU_*sigEvent.sfWeightEff_*sigEvent.sfWeightTrig_)-add)/add>0.05 && sigEvent.sfWeightFR_ > 0 && sigEvent.lid3_ != 0) {
	printf("PROBLEM: %f - %f %f %f %f %f = %f\n",add,sigEvent.sfWeightFR_,sigEvent.sfWeightPU_,sigEvent.sfWeightEff_,sigEvent.sfWeightTrig_,sigEvent.sfWeightHPt_,sigEvent.sfWeightFR_*sigEvent.sfWeightPU_*sigEvent.sfWeightEff_*sigEvent.sfWeightTrig_*sigEvent.sfWeightHPt_);
	printf("  		%f %f %f %f\n",1.0,addPU,leptonEfficiency(sigEvent.lep1_.Pt(), sigEvent.lep1_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid1_)*
  	      leptonEfficiency(sigEvent.lep2_.Pt(), sigEvent.lep2_.Eta(), fhDEffMu, fhDEffEl, sigEvent.lid2_), trigLookup.GetExpectedTriggerEfficiency(fabs(sigEvent.lep1_.Eta()), sigEvent.lep1_.Pt() , 
  																			     fabs(sigEvent.lep2_.Eta()), sigEvent.lep2_.Pt(), 
  																			     TMath::Abs( sigEvent.lid1_), TMath::Abs(sigEvent.lid2_)));
      }

      if(passCuts[1][WWSEL] ){ // begin making plots
	double myVar = theMET;
	if     (thePlot == 1) myVar = sigEvent.lep1_.Pt();
	else if(thePlot == 2) myVar = sigEvent.lep2_.Pt();
	else if(thePlot == 3) myVar = sigEvent.lep3_.Pt();
	else if(thePlot == 4) myVar = sigEvent.jet1_.Pt();
	else if(thePlot == 5) myVar = sigEvent.jet2_.Pt();
	else if(thePlot == 6) myVar = sigEvent.jet3_.Pt();
	else if(thePlot == 7) myVar = TMath::Min(sigEvent.dilep_.M(),399.999);
	else if(thePlot == 8) myVar = sigEvent.mt_;
	else if(thePlot == 9) myVar = sigEvent.mt1_;
	else if(thePlot ==10) myVar = sigEvent.mt2_;
	else if(thePlot ==12) myVar = sigEvent.trackMet_;
	else if(thePlot ==13) myVar = sigEvent.dilep_.Pt();
	else if(thePlot ==14) myVar = fabs(sigEvent.dilep_.M()-91.1876);
	else if(thePlot ==15) myVar = fabs(theMET-sigEvent.dilep_.Pt())/sigEvent.dilep_.Pt();
	else if(thePlot ==16) myVar = sigEvent.lep2_.Pt()/sigEvent.lep1_.Pt();
	else if(thePlot ==17) myVar = sigEvent.njets_;
	else if(thePlot ==18) myVar = sigEvent.nvtx_;
	else if(thePlot ==19) myVar = TMath::Max(TMath::Min((sigEvent.lep1_+sigEvent.lep2_+sigEvent.jet1_+sigEvent.jet2_).M(),2999.999),700.001);
	else if(thePlot ==20) myVar = sigEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==21) myVar = TMath::Min(sigEvent.dPhiLep1MET_,sigEvent.dPhiLep2MET_)*180.0/TMath::Pi();
	else if(thePlot ==22) myVar = DeltaPhi(sigEvent.dilep_.Phi() ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==23) myVar = DeltaPhi(sigEvent.trackMetPhi_ ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==24) myVar = sigEvent.lep1_.Pt()*sigEvent.lep2_.Pt()/sigEvent.jet1_.Pt()/sigEvent.jet2_.Pt();
	else if(thePlot ==25) myVar = fabs(sigEvent.dilep_.Eta());
	else if(thePlot ==26) myVar = sigEvent.dymva_;
	else if(thePlot ==27) myVar = TMath::Min(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
	else if(thePlot ==28) myVar = TMath::Max(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
	else if(thePlot ==29) myVar = TMath::Max(fabs(sigEvent.jet1_.Eta()),fabs(sigEvent.jet2_.Eta()));
	else if(thePlot ==30) myVar = DeltaPhi(sigEvent.dilep_.Phi() ,sigEvent.jet1_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==37) myVar = TMath::Max(TMath::Min((sigEvent.jet1_+sigEvent.jet2_).M(),1999.999),500.001);
	else if(thePlot ==38) myVar = TMath::Abs(sigEvent.jet1_.Eta()-sigEvent.jet2_.Eta());
	else if(thePlot ==40) myVar = DeltaPhi(sigEvent.jet1_.Phi() ,sigEvent.jet2_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==41) myVar = DeltaPhi(sigEvent.trackMetPhi_,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==44) myVar = sigEvent.jet1_.Pt()+ sigEvent.jet2_.Pt()+sigEvent.jet3_.Pt();
	else if(thePlot ==48) myVar = sigEvent.type_;
	else if(thePlot ==49) myVar = theMET*cos(theMETPHI);
	else if(thePlot ==50) myVar = theMET*sin(theMETPHI);
	else if(thePlot ==51) myVar = sigEvent.trackMet_*cos(sigEvent.trackMetPhi_);
	else if(thePlot ==52) myVar = sigEvent.trackMet_*sin(sigEvent.trackMetPhi_);
	else if(thePlot ==53) myVar = DeltaPhi(sigEvent.jet3_.Phi(),sigEvent.jet4_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==55) myVar = sigEvent.dPhiDiLepMET_*180.0/TMath::Pi();
	else if(thePlot ==57) myVar = sigEvent.dR_;
      	histos->Fill(myVar,theWeight);
	
	if (sigEvent.lheWeights_.size() != grid_points.size()){
	  std::cout << "theWeight = " << theWeight << std::endl;
	std::cout << "sigEvent.lheWeights_.size() = " << sigEvent.lheWeights_.size() << std::endl;
	std::cout << "grid_points.size() = " << grid_points.size() << std::endl;
	}

	assert(sigEvent.lheWeights_.size() == grid_points.size());
	assert(sigEvent.lheWeights_.size() == histo_grid.size());

	for(int a = 0; a < grid_points.size(); a++){
	  if (theWeight > 0)
	    histo_grid[a]+=theWeight*sigEvent.lheWeights_[a]/sigEvent.lheWeights_[0]*29700./28000.;
	  else 
	    histo_grid[a]+=theWeight*sigEvent.lheWeights_[a]/sigEvent.lheWeights_[0]*16000./15600.;
	}
      } // end making plots

      for(unsigned int i=0; i<nSelTypes; i++) {
        for(int j=0; j<2; j++){
          if(passCuts[j][i]) {
            nSigCut[i+j*nSelTypes]  += theWeight;
            nSigECut[i+j*nSelTypes] += theWeight*theWeight;
          }
        }
      }
      for(unsigned int i=0; i<5; i++) {
        for(int j=0; j<2; j++){
          if(passSystCuts[j][i]) {
            nSigCutSyst[i+j*nSelTypesSyst]  += theWeight;
            nSigECutSyst[i+j*nSelTypesSyst] += theWeight*theWeight;
          }
        }
      }
      if(passCuts[lType][WWSEL]) {
        nSigCutSyst[EFFP+lType*nSelTypesSyst]  += theWeight          *addLepEffUp  /addLepEff;
        nSigECutSyst[EFFP+lType*nSelTypesSyst] += theWeight*theWeight*addLepEffUp  /addLepEff*addLepEffUp  /addLepEff;
        nSigCutSyst[EFFM+lType*nSelTypesSyst]  += theWeight          *addLepEffDown/addLepEff;
        nSigECutSyst[EFFM+lType*nSelTypesSyst] += theWeight*theWeight*addLepEffDown/addLepEff*addLepEffDown/addLepEff;
      }
      if(passLSel){
        if(passCuts[1][WWSEL]) 	             histo_WWewk_ALT  			       ->Fill(MVAVar[0], theWeight);
        if(passCuts[1][WWSEL]) 	             histo_WWewk_ALT_CMS_MVALepEffBoundingUp   ->Fill(MVAVar[0], theWeight*addLepEffUp  /addLepEff);
        if(passCuts[1][WWSEL]) 	             histo_WWewk_ALT_CMS_MVALepEffBoundingDown ->Fill(MVAVar[0], theWeight*addLepEffDown/addLepEff);
        if(passSystCuts[1][JESUP  ] == true) histo_WWewk_ALT_CMS_MVAJESBoundingUp      ->Fill(MVAVar[1], theWeight);
        if(passSystCuts[1][JESDOWN] == true) histo_WWewk_ALT_CMS_MVAJESBoundingDown    ->Fill(MVAVar[2], theWeight);
        if(passSystCuts[1][LEPP]    == true) histo_WWewk_ALT_CMS_MVALepResBoundingUp   ->Fill(MVAVar[3], theWeight);
        if(passSystCuts[1][LEPM]    == true) histo_WWewk_ALT_CMS_MVALepResBoundingDown ->Fill(MVAVar[4], theWeight);
        if(passSystCuts[1][MET]     == true) histo_WWewk_ALT_CMS_MVAMETResBoundingUp   ->Fill(MVAVar[5], theWeight);;
      }
    } // if passCuts
  } // Loop over signal

  for(int a = 0; a < grid_points.size(); a++){
    th2d->SetBinContent(th2d->GetXaxis()->FindFixBin(grid_points[a].first), th2d->GetYaxis()->FindFixBin(grid_points[a].second), histo_grid[a]/histo_grid[0]);
  }

  int nData=dataEvent.tree_->GetEntries();
  for (int evt=0; evt<nData; ++evt) {

    if (evt%100000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",evt,nData);
    dataEvent.tree_->GetEntry(evt);

    bool lId = (dataEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (dataEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;

    if(!lId) continue;

    if((dataEvent.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ <  minRun) continue;
    if(dataEvent.dstype_ == SmurfTree::data && dataEvent.run_ >  maxRun) continue;

    int fDecay = 0;
    if(fDecay == -1 || fDecay > 100) fDecay = 0;

    bool passCuts[3][nSelTypes] = {{false, false, false},
                                   {false, false, false}};

    double theMET = dataEvent.met_; double theMETPHI = dataEvent.metPhi_; 

    int lType = 1;
    if     (dataEvent.lq1_ * dataEvent.lq2_ < 0) lType = 0;

    int centrality = 0;
    if(((dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() < 0) ||
        (dataEvent.jet2_.Eta()-dataEvent.lep1_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep1_.Eta() < 0)) &&
       ((dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() < 0) ||
        (dataEvent.jet2_.Eta()-dataEvent.lep2_.Eta() > 0 && dataEvent.jet1_.Eta()-dataEvent.lep2_.Eta() < 0))) centrality = 1; 
    double metMin = 30.0; if(dataEvent.type_ == SmurfTree::ee) metMin = 40.0;
    if(lType == 0) if(dataEvent.type_ == SmurfTree::mm) metMin = 40.0;

    bool passNjets         = dataEvent.njets_ >= 2;
    bool passMET           = TMath::Min(dataEvent.pmet_,dataEvent.pTrackMet_) > metMin;
    bool preselCuts        = dataEvent.lep1_.Pt() > 20. && dataEvent.lep2_.Pt() > 20.;
    bool passBtagVeto      = (dataEvent.cuts_ & patternTopVeto) == patternTopVeto;
    bool pass3rLVeto       = (dataEvent.cuts_ & SmurfTree::ExtraLeptonVeto) == SmurfTree::ExtraLeptonVeto;
    bool passMass          = dataEvent.dilep_.M() > 50.0 && (TMath::Abs(dataEvent.dilep_.M()-91.1876) > 15 || dataEvent.type_ != SmurfTree::ee);
    bool passVBFSel        = (dataEvent.jet1_+dataEvent.jet2_).M() > 500 && TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta()) > 3.5 && centrality == 1;

    if(lType == 0) passMass = passMass && (TMath::Abs(dataEvent.dilep_.M()-91.1876) > 15 || dataEvent.type_ != SmurfTree::mm);

    bool passLSel = false;
    if     (lSel == 0 && dataEvent.type_ == SmurfTree::mm) passLSel = true;
    else if(lSel == 1 && dataEvent.type_ == SmurfTree::me) passLSel = true;
    else if(lSel == 2 && dataEvent.type_ == SmurfTree::em) passLSel = true;
    else if(lSel == 3 && dataEvent.type_ == SmurfTree::ee) passLSel = true;
    else if(lSel == 4)                                     passLSel = true;
    else if(lSel == 5 && (dataEvent.type_ == SmurfTree::mm || dataEvent.type_ == SmurfTree::ee)) passLSel = true;
    else if(lSel == 6 && (dataEvent.type_ == SmurfTree::me || dataEvent.type_ == SmurfTree::em)) passLSel = true;

    if(passNjets  == true && passMET == true &&  passLSel == true &&
       preselCuts == true && dataEvent.dilep_.M() > 15.0) {
       
       if( passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto) passCuts[lType][WWSEL] = true;
       if(!passBtagVeto && passVBFSel == true && passMass  == true &&  pass3rLVeto) passCuts[lType][BTAGSEL] = true;
       if( passBtagVeto && passVBFSel == true                      && !pass3rLVeto && dataEvent.lep3_.Pt() > 20.) passCuts[lType][WZSEL] = true;

    }

    if(passNjets  == true && passMET == true &&  passLSel == true &&
       preselCuts == true && dataEvent.dilep_.M() > 15.0) {

      if(passCuts[1][WWSEL]){ // begin making plots
	double myVar = theMET;
	if     (thePlot == 1) myVar = dataEvent.lep1_.Pt();
	else if(thePlot == 2) myVar = dataEvent.lep2_.Pt();
	else if(thePlot == 3) myVar = dataEvent.lep3_.Pt();
	else if(thePlot == 4) myVar = dataEvent.jet1_.Pt();
	else if(thePlot == 5) myVar = dataEvent.jet2_.Pt();
	else if(thePlot == 6) myVar = dataEvent.jet3_.Pt();
	else if(thePlot == 7) myVar = TMath::Min(dataEvent.dilep_.M(),399.999);
	else if(thePlot == 8) myVar = dataEvent.mt_;
	else if(thePlot == 9) myVar = dataEvent.mt1_;
	else if(thePlot ==10) myVar = dataEvent.mt2_;
	else if(thePlot ==12) myVar = dataEvent.trackMet_;
	else if(thePlot ==13) myVar = dataEvent.dilep_.Pt();
	else if(thePlot ==14) myVar = fabs(dataEvent.dilep_.M()-91.1876);
	else if(thePlot ==15) myVar = fabs(theMET-dataEvent.dilep_.Pt())/dataEvent.dilep_.Pt();
	else if(thePlot ==16) myVar = dataEvent.lep2_.Pt()/dataEvent.lep1_.Pt();
	else if(thePlot ==17) myVar = dataEvent.njets_;
	else if(thePlot ==18) myVar = dataEvent.nvtx_;
	else if(thePlot ==19) myVar = TMath::Max(TMath::Min((dataEvent.lep1_+dataEvent.lep2_+dataEvent.jet1_+dataEvent.jet2_).M(),2999.999),700.001);
	else if(thePlot ==20) myVar = dataEvent.dPhi_*180.0/TMath::Pi();
	else if(thePlot ==21) myVar = TMath::Min(dataEvent.dPhiLep1MET_,dataEvent.dPhiLep2MET_)*180.0/TMath::Pi();
	else if(thePlot ==22) myVar = DeltaPhi(dataEvent.dilep_.Phi() ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==23) myVar = DeltaPhi(dataEvent.trackMetPhi_ ,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==24) myVar = dataEvent.lep1_.Pt()*dataEvent.lep2_.Pt()/dataEvent.jet1_.Pt()/dataEvent.jet2_.Pt();
	else if(thePlot ==25) myVar = fabs(dataEvent.dilep_.Eta());
	else if(thePlot ==26) myVar = dataEvent.dymva_;
	else if(thePlot ==27) myVar = TMath::Min(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
	else if(thePlot ==28) myVar = TMath::Max(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
	else if(thePlot ==29) myVar = TMath::Max(fabs(dataEvent.jet1_.Eta()),fabs(dataEvent.jet2_.Eta()));
	else if(thePlot ==30) myVar = DeltaPhi(dataEvent.dilep_.Phi() ,dataEvent.jet1_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==37) myVar = TMath::Max(TMath::Min((dataEvent.jet1_+dataEvent.jet2_).M(),1999.999),500.001);
	else if(thePlot ==38) myVar = TMath::Abs(dataEvent.jet1_.Eta()-dataEvent.jet2_.Eta());
	else if(thePlot ==40) myVar = DeltaPhi(dataEvent.jet1_.Phi() ,dataEvent.jet2_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==41) myVar = DeltaPhi(dataEvent.trackMetPhi_,theMETPHI)*180.0/TMath::Pi();
	else if(thePlot ==44) myVar = dataEvent.jet1_.Pt()+ dataEvent.jet2_.Pt()+dataEvent.jet3_.Pt();
	else if(thePlot ==48) myVar = dataEvent.type_;
	else if(thePlot ==49) myVar = theMET*cos(theMETPHI);
	else if(thePlot ==50) myVar = theMET*sin(theMETPHI);
	else if(thePlot ==51) myVar = dataEvent.trackMet_*cos(dataEvent.trackMetPhi_);
	else if(thePlot ==52) myVar = dataEvent.trackMet_*sin(dataEvent.trackMetPhi_);
	else if(thePlot ==53) myVar = DeltaPhi(dataEvent.jet3_.Phi(),dataEvent.jet4_.Phi())*180.0/TMath::Pi();
	else if(thePlot ==55) myVar = dataEvent.dPhiDiLepMET_*180.0/TMath::Pi();
	else if(thePlot ==57) myVar = dataEvent.dR_;
      	histo5->Fill(myVar,1.0);
      } // end making plots

      double outputVar[15];
      makeSystematicEffects(dataEvent.lid1_, dataEvent.lid2_, dataEvent.lep1_, dataEvent.lep2_, dataEvent.dilep_, 
                            dataEvent.mt_, theMET, theMETPHI, 
                            dataEvent.trackMet_, dataEvent.trackMetPhi_, 
			    dataEvent.njets_, dataEvent.jet1_, dataEvent.jet2_, 
			    year, 3, outputVar);
      double MVAVar[6] = {outputVar[13],0,0,0,0,0};
      if(thePlot == 37) {MVAVar[0]=outputVar[14];}
      for(int nv=0; nv<6; nv++) MVAVar[nv] = TMath::Min(TMath::Max(MVAVar[nv],xbins[0]+0.001),xbins[nBinMVA]-0.001);
      if(passCuts[1][WWSEL]){
	histo_Data->Fill(MVAVar[0], 1.0);
      }

      for(unsigned int i=0; i<nSelTypes; i++) {
        for(int j=0; j<2; j++){
          if(passCuts[j][i]) {
            nSelectedData[i+j*nSelTypes]  += 1.0;
          }
        }
      }

    } // if passCuts
  } // End loop data

  char output[200];
  sprintf(output,Form("histo_nice%s.root",ECMsb.Data()));	 
  TFile* outFilePlotsNote = new TFile(output,"recreate");
  outFilePlotsNote->cd();
    double nOldH[6] = {histo0->GetSumOfWeights(),histo1->GetSumOfWeights(),histo2->GetSumOfWeights(),histo3->GetSumOfWeights(),histo4->GetSumOfWeights(),histos->GetSumOfWeights()};
    for(int i=1; i<=histo0->GetNbinsX(); i++){
      if(histo0->GetBinContent(i) < 0) {histo0->SetBinContent(i,0.000001);histo0->SetBinError(i,0.000001);}
      if(histo1->GetBinContent(i) < 0) {histo1->SetBinContent(i,0.000001);histo1->SetBinError(i,0.000001);}
      if(histo2->GetBinContent(i) < 0) {histo2->SetBinContent(i,0.000001);histo2->SetBinError(i,0.000001);}
      if(histo3->GetBinContent(i) < 0) {histo3->SetBinContent(i,0.000001);histo3->SetBinError(i,0.000001);}
      if(histo4->GetBinContent(i) < 0) {histo4->SetBinContent(i,0.000001);histo4->SetBinError(i,0.000001);}
      if(histos->GetBinContent(i) < 0) {histos->SetBinContent(i,0.000001);histos->SetBinError(i,0.000001);}
    }
    if(nOldH[0] > 0) histo0->Scale(nOldH[0]/histo0->GetSumOfWeights());
    if(nOldH[1] > 0) histo1->Scale(nOldH[1]/histo1->GetSumOfWeights());
    if(nOldH[2] > 0) histo2->Scale(nOldH[2]/histo2->GetSumOfWeights());
    if(nOldH[3] > 0) histo3->Scale(nOldH[3]/histo3->GetSumOfWeights());
    if(nOldH[4] > 0) histo4->Scale(nOldH[4]/histo4->GetSumOfWeights());
    if(nOldH[5] > 0) histos->Scale(nOldH[5]/histos->GetSumOfWeights());

    printf("histo -> s: %8.2f d: %8.2f b: %8.2f | %8.2f %8.2f %8.2f %8.2f %8.2f\n",histos->GetSumOfWeights(),histo5->GetSumOfWeights(),
    histo0->GetSumOfWeights()+histo1->GetSumOfWeights()+histo2->GetSumOfWeights()+histo3->GetSumOfWeights()+histo4->GetSumOfWeights(),
    histo0->GetSumOfWeights(),histo1->GetSumOfWeights(),histo2->GetSumOfWeights(),histo3->GetSumOfWeights(),histo4->GetSumOfWeights());

    //histos->Write();
    //histo0->Write();
    //histo1->Write();
    //histo2->Write();
    //histo3->Write();
    //histo4->Write();
    //histo5->Write();

    std::cout << "histo0->GetNbinsX() = " << histo0->GetNbinsX() << std::endl;

    data_obs->Fill(0.5,histo1->Integral()+histo2->Integral()+histo3->Integral()+histo4->Integral()+histo0->Integral());
    diboson->Fill(0.5,histo0->Integral());
    background->Fill(0.5,histo1->Integral()+histo2->Integral()+histo3->Integral()+histo4->Integral()+histo0->Integral());
    background_ch1boosted_backshapeUp->Fill(0.5,histo1->Integral()+histo2->Integral()+histo3->Integral()+histo4->Integral()+histo0->Integral());
    background_ch1boosted_backshapeDown->Fill(0.5,histo1->Integral()+histo2->Integral()+histo3->Integral()+histo4->Integral()+histo0->Integral());
    background_ch2boosted_backshapeUp->Fill(0.5,histo1->Integral()+histo2->Integral()+histo3->Integral()+histo4->Integral()+histo0->Integral());
    background_ch2boosted_backshapeDown->Fill(0.5,histo1->Integral()+histo2->Integral()+histo3->Integral()+histo4->Integral()+histo0->Integral());

    //data_obs->Write();
    //diboson->Write();
    //background->Write();
    //background_ch1boosted_backshapeUp->Write();
    //background_ch1boosted_backshapeDown->Write();
    //background_ch2boosted_backshapeUp->Write();
    //background_ch2boosted_backshapeDown->Write();

    th2d->Write();
    th2d->Clone("bin_content_lam_dg_1")->Write();
    th2d->Clone("bin_content_dk_dg_1")->Write();

  outFilePlotsNote->Close();
  
  const unsigned int nBkg = 6;
  double nTot[nSelTypes*2]; double nETot[nSelTypes*2];
  double bgdCombined[nSelTypes*2][nBkg],bgdCombinedE[nSelTypes*2][nBkg];
  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(unsigned int j=0; j<nBkg; j++) {bgdCombined[i][j] = 0.0; bgdCombinedE[i][j] = 0.0;}
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("selection: %s\n",selTypeName[i].Data());
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("data(%2d): %f\n",i,nSelectedData[i]);
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("nSigCut(%2d): %11.3f +/- %8.3f\n",i,nSigCut[i],sqrt(nSigECut[i]));
    nTot[i] = 0.0; nETot[i] = 0.0;
    for(int j=0; j<45; j++){
      // WWqcd treatment
      if(j == 29 && bgdDecay[i][j] < 0) {printf("negative(29,%d) = %f +/- %f\n",i,bgdDecay[i][j],sqrt(weiDecay[i][j]));bgdDecay[i][j] = 0; weiDecay[i][j] = 0;}

      if(showSignalOnly == false || i%nSelTypes == WWSEL) if(bgdDecay[i][j] != 0) printf("bdg(%2d,%2d) = %11.3f +/- %8.3f\n",i,j,bgdDecay[i][j],sqrt(weiDecay[i][j]));

      nTot[i]  += bgdDecay[i][j];
      nETot[i] += weiDecay[i][j];

      if     (j == 31)                         {bgdCombined[i][0] += bgdDecay[i][j]; bgdCombinedE[i][0] += weiDecay[i][j];}
      else if(j == 29)                         {bgdCombined[i][1] += bgdDecay[i][j]; bgdCombinedE[i][1] += weiDecay[i][j];}
      else if(j == 27)			       {bgdCombined[i][2] += bgdDecay[i][j]; bgdCombinedE[i][2] += weiDecay[i][j];}
      else if(j == 30 || j == 28 ||
              j ==  5 || j == 13 || j == 20 || 
	      j == 10 || j ==  9 || j == 19)   {bgdCombined[i][3] += bgdDecay[i][j]; bgdCombinedE[i][3] += weiDecay[i][j];}
      else if(j == 21)			       {bgdCombined[i][4] += bgdDecay[i][j]; bgdCombinedE[i][4] += weiDecay[i][j];}
      else if(j == 1 || j == 23)               {bgdCombined[i][5] += bgdDecay[i][j]; bgdCombinedE[i][5] += weiDecay[i][j];}
    }
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("nTot(%2d) = %11.3f +/- %8.3f\n",i,nTot[i],sqrt(nETot[i]));
    if(nTot[i] > 0.0 && TMath::Abs(bgdCombined[i][0]+bgdCombined[i][1]+bgdCombined[i][2]+bgdCombined[i][3]+bgdCombined[i][4]+bgdCombined[i][5]-nTot[i])/nTot[i] > 0.00001) 
                    {printf("%f\n",bgdCombined[i][0]+bgdCombined[i][1]+bgdCombined[i][2]+bgdCombined[i][3]+bgdCombined[i][4]+bgdCombined[i][5]);assert(0);}
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("------\n");
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(WWe) = %11.3f +/- %8.3f\n",bgdCombined[i][0],sqrt(bgdCombinedE[i][0]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(WWq) = %11.3f +/- %8.3f\n",bgdCombined[i][1],sqrt(bgdCombinedE[i][1]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(xWZ) = %11.3f +/- %8.3f\n",bgdCombined[i][2],sqrt(bgdCombinedE[i][2]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(xWS) = %11.3f +/- %8.3f\n",bgdCombined[i][3],sqrt(bgdCombinedE[i][3]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(VVV) = %11.3f +/- %8.3f\n",bgdCombined[i][4],sqrt(bgdCombinedE[i][4]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("bgd(xWj) = %11.3f +/- %8.3f\n",bgdCombined[i][5],sqrt(bgdCombinedE[i][5]));
    if(showSignalOnly == false || i%nSelTypes == WWSEL) printf("*******************************\n");
  }

  if(showSignalOnly == false) printf("+++++++++++++++++++++++++++++++\n");
  double nTotSyst[nSelTypesSyst*2]; double nETotSyst[nSelTypesSyst*2];
  double bgdCombinedSyst[nSelTypesSyst*2][nBkg],bgdCombinedESyst[nSelTypesSyst*2][nBkg];
  for(unsigned int i=0; i<nSelTypesSyst*2; i++) {
    if(showSignalOnly == false) printf("selectionSyst: %s\n",selTypeNameSyst[i].Data());
    for(unsigned int j=0; j<nBkg; j++) {bgdCombinedSyst[i][j] = 0.0; bgdCombinedESyst[i][j] = 0.0;}
    if(showSignalOnly == false) printf("nSigCutSyst(%2d): %11.3f +/- %8.3f\n",i,nSigCutSyst[i],sqrt(nSigECutSyst[i]));
    nTotSyst[i] = 0.0; nETotSyst[i] = 0.0;
    for(int j=0; j<45; j++){
      if(showSignalOnly == false) if(bgdDecaySyst[i][j] != 0) printf("bdgSyst(%2d,%2d) = %11.3f +/- %8.3f\n",i,j,bgdDecaySyst[i][j],sqrt(weiDecaySyst[i][j]));
      nTotSyst[i]  += bgdDecaySyst[i][j];
      nETotSyst[i] += weiDecaySyst[i][j];

     if     (j == 31)			      {bgdCombinedSyst[i][0] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][0] += weiDecaySyst[i][j];}
     else if(j == 29)			      {bgdCombinedSyst[i][1] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][1] += weiDecaySyst[i][j];}
     else if(j == 27)			      {bgdCombinedSyst[i][2] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][2] += weiDecaySyst[i][j];}
     else if(j == 30 || j == 28 ||		      
     	     j ==  5 || j == 13 || j == 20 || 
     	     j == 10 || j ==  9 || j == 19)   {bgdCombinedSyst[i][3] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][3] += weiDecaySyst[i][j];}
     else if(j == 21)			      {bgdCombinedSyst[i][4] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][4] += weiDecaySyst[i][j];}
     else if(j == 1 || j == 23) 	      {bgdCombinedSyst[i][5] += bgdDecaySyst[i][j]; bgdCombinedESyst[i][5] += weiDecaySyst[i][j];}
    }
    if(showSignalOnly == false) printf("nTot(%2d) = %11.3f +/- %8.3f\n",i,nTotSyst[i],sqrt(nETotSyst[i]));
    if(nTot[i] > 0.0 && TMath::Abs(bgdCombinedSyst[i][0]+bgdCombinedSyst[i][1]+bgdCombinedSyst[i][2]+bgdCombinedSyst[i][3]+bgdCombinedSyst[i][4]+bgdCombinedSyst[i][5]-nTotSyst[i])/nTotSyst[i] > 0.00001) 
                    {printf("%f\n",bgdCombinedSyst[i][0]+bgdCombinedSyst[i][1]+bgdCombinedSyst[i][2]+bgdCombinedSyst[i][3]+bgdCombinedSyst[i][4]+bgdCombinedSyst[i][5]);assert(0);}
    if(showSignalOnly == false) printf("------\n");
    if(showSignalOnly == false) printf("bgdSyst(WWe) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][0],sqrt(bgdCombinedESyst[i][0]));
    if(showSignalOnly == false) printf("bgdSyst(WWq) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][1],sqrt(bgdCombinedESyst[i][1]));
    if(showSignalOnly == false) printf("bgdSyst(xWZ) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][2],sqrt(bgdCombinedESyst[i][2]));
    if(showSignalOnly == false) printf("bgdSyst(xWS) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][3],sqrt(bgdCombinedESyst[i][3]));
    if(showSignalOnly == false) printf("bgdSyst(VVV) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][4],sqrt(bgdCombinedESyst[i][4]));
    if(showSignalOnly == false) printf("bgdSyst(xWj) = %11.3f +/- %8.3f\n",bgdCombinedSyst[i][5],sqrt(bgdCombinedESyst[i][5]));
    if(showSignalOnly == false) printf("*******************************\n");
  }

  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(unsigned int j=0; j<nBkg; j++) if(bgdCombined[i][j] == 0) {bgdCombined[i][j] = 0.0000000001; bgdCombinedE[i][j] = 0.0;}
  }
  double NFinal[nBkg+1]  = {bgdCombined[WWSEL+nSelTypes][0],bgdCombined[WWSEL+nSelTypes][1],bgdCombined[WWSEL+nSelTypes][2],
                            bgdCombined[WWSEL+nSelTypes][3],bgdCombined[WWSEL+nSelTypes][4],bgdCombined[WWSEL+nSelTypes][5],nSigCut[WWSEL+nSelTypes]};
  //double NFinalE[nBkg+1] = {1.0+sqrt(bgdCombinedE[WWSEL+nSelTypes][0])/bgdCombined[WWSEL+nSelTypes][0],1.0+sqrt(bgdCombinedE[WWSEL+nSelTypes][1])/bgdCombined[WWSEL+nSelTypes][1],
  //                   	    1.0+sqrt(bgdCombinedE[WWSEL+nSelTypes][2])/bgdCombined[WWSEL+nSelTypes][2],1.0+sqrt(bgdCombinedE[WWSEL+nSelTypes][3])/bgdCombined[WWSEL+nSelTypes][3],
  //		            1.0+sqrt(bgdCombinedE[WWSEL+nSelTypes][4])/bgdCombined[WWSEL+nSelTypes][4],1.0+sqrt(nSigECut[WWSEL+nSelTypes])/nSigCut[WWSEL+nSelTypes]};

  double QCDscale_WWewk = 1.05;
  
  double systEffect[nSelTypesSyst][nBkg+1];
  for(unsigned int i=0 ; i<nSelTypesSyst; i++){
    for(unsigned int j=0 ; j<nBkg; j++){
      if(bgdCombinedE[WWSEL+nSelTypes][j] > 0){
        systEffect[i][j] = bgdCombinedSyst[i+nSelTypesSyst][j]/bgdCombined[WWSEL+nSelTypes][j];
        if(systEffect[i][j] < 1) systEffect[i][j] = 1.0/systEffect[i][j];
      } else {systEffect[i][j] = 1.0;}
    }
    systEffect[i][nBkg] = nSigCutSyst[i+nSelTypesSyst]/nSigCut[WWSEL+nSelTypes];    
    if(systEffect[i][nBkg] < 1) systEffect[i][nBkg] = 1.0/systEffect[i][nBkg];
  }
  if(showSignalOnly == false) printf("Syst(sig) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][nBkg]-1,systEffect[JESDOWN][nBkg]-1,systEffect[LEPP][nBkg]-1,systEffect[LEPM][nBkg]-1,systEffect[MET][nBkg]-1,systEffect[EFFP][nBkg]-1,systEffect[EFFM][nBkg]-1);
  if(showSignalOnly == false) printf("Syst(WWe) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][0]-1,systEffect[JESDOWN][0]-1,systEffect[LEPP][0]-1,systEffect[LEPM][0]-1,systEffect[MET][0]-1,systEffect[EFFP][0]-1,systEffect[EFFM][0]-1);
  if(showSignalOnly == false) printf("Syst(WWq) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][1]-1,systEffect[JESDOWN][1]-1,systEffect[LEPP][1]-1,systEffect[LEPM][1]-1,systEffect[MET][1]-1,systEffect[EFFP][1]-1,systEffect[EFFM][1]-1);
  if(showSignalOnly == false) printf("Syst(xWZ) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][2]-1,systEffect[JESDOWN][2]-1,systEffect[LEPP][2]-1,systEffect[LEPM][2]-1,systEffect[MET][2]-1,systEffect[EFFP][2]-1,systEffect[EFFM][2]-1);
  if(showSignalOnly == false) printf("Syst(xWS) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][3]-1,systEffect[JESDOWN][3]-1,systEffect[LEPP][3]-1,systEffect[LEPM][3]-1,systEffect[MET][3]-1,systEffect[EFFP][3]-1,systEffect[EFFM][3]-1);
  if(showSignalOnly == false) printf("Syst(VVV) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][4]-1,systEffect[JESDOWN][4]-1,systEffect[LEPP][4]-1,systEffect[LEPM][4]-1,systEffect[MET][4]-1,systEffect[EFFP][4]-1,systEffect[EFFM][4]-1);
  if(showSignalOnly == false) printf("Syst(xWj) = %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n",systEffect[JESUP][5]-1,systEffect[JESDOWN][5]-1,systEffect[LEPP][5]-1,systEffect[LEPM][5]-1,systEffect[MET][5]-1,systEffect[EFFP][5]-1,systEffect[EFFM][5]-1);

  double WjetsSyst = 1.0;
  if(bgdCombined[WWSEL+nSelTypes][5] > 0){
    WjetsSyst = bgdCombined[WWSEL+nSelTypes][5]*0.36;
    WjetsSyst = 1.0+WjetsSyst/(bgdCombined[WWSEL+nSelTypes][5]);
  }
  if(showSignalOnly == false) printf("WjetsSyst: %f --> %f\n",bgdCombined[WWSEL+nSelTypes][5],WjetsSyst);
  double pdf_qqbar[4] = {1.073,1.073,1.068,1.069};
  double syst_WZ3l = 1.010;

  double nOldWjets = TMath::Max(histo_Wjets->GetSumOfWeights(),0.000001);
  for(int i=1; i<=histo_Wjets->GetNbinsX(); i++){
    if(histo_Wjets->GetBinContent(i)                    < 0) {histo_Wjets		    ->SetBinContent(i,0.000001);histo_Wjets		      ->SetBinError(i,0.000001);}
    if(histo_Wjets_CMS_MVAWBoundingUp->GetBinContent(i) < 0) {histo_Wjets_CMS_MVAWBoundingUp->SetBinContent(i,0.000001);histo_Wjets_CMS_MVAWBoundingUp->SetBinError(i,0.000001);}
  }
  histo_Wjets                   ->Scale(nOldWjets/histo_Wjets                   ->GetSumOfWeights());
  histo_Wjets_CMS_MVAWBoundingUp->Scale(nOldWjets/histo_Wjets_CMS_MVAWBoundingUp->GetSumOfWeights());

  // WWqcd treatment
  if(histo_WWqcd->GetSumOfWeights() > 0){
    double nOldWWQCD[8] = {TMath::Max(histo_WWqcd->GetSumOfWeights(),0.000001),
                           TMath::Max(histo_WWqcd_CMS_MVALepEffBoundingUp->GetSumOfWeights(),0.000001),TMath::Max(histo_WWqcd_CMS_MVALepEffBoundingDown->GetSumOfWeights(),0.000001),
                           TMath::Max(histo_WWqcd_CMS_MVAJESBoundingUp->GetSumOfWeights(),0.000001),   TMath::Max(histo_WWqcd_CMS_MVAJESBoundingDown->GetSumOfWeights(),0.000001),
			   TMath::Max(histo_WWqcd_CMS_MVALepResBoundingUp->GetSumOfWeights(),0.000001),TMath::Max(histo_WWqcd_CMS_MVALepResBoundingDown->GetSumOfWeights(),0.000001),
			   TMath::Max(histo_WWqcd_CMS_MVAMETResBoundingUp->GetSumOfWeights(),0.000001)};
    for(int i=1; i<=histo_WWqcd->GetNbinsX(); i++){
      if(histo_WWqcd			    ->GetBinContent(i) < 0) {histo_WWqcd			  ->SetBinContent(i,0.000001);histo_WWqcd			   ->SetBinError(i,0.000001);}
      if(histo_WWqcd_CMS_MVALepEffBoundingUp  ->GetBinContent(i) < 0) {histo_WWqcd_CMS_MVALepEffBoundingUp  ->SetBinContent(i,0.000001);histo_WWqcd_CMS_MVALepEffBoundingUp  ->SetBinError(i,0.000001);}
      if(histo_WWqcd_CMS_MVALepEffBoundingDown->GetBinContent(i) < 0) {histo_WWqcd_CMS_MVALepEffBoundingDown->SetBinContent(i,0.000001);histo_WWqcd_CMS_MVALepEffBoundingDown->SetBinError(i,0.000001);}
      if(histo_WWqcd_CMS_MVAJESBoundingUp     ->GetBinContent(i) < 0) {histo_WWqcd_CMS_MVAJESBoundingUp	  ->SetBinContent(i,0.000001);histo_WWqcd_CMS_MVAJESBoundingUp     ->SetBinError(i,0.000001);}
      if(histo_WWqcd_CMS_MVAJESBoundingDown   ->GetBinContent(i) < 0) {histo_WWqcd_CMS_MVAJESBoundingDown   ->SetBinContent(i,0.000001);histo_WWqcd_CMS_MVAJESBoundingDown   ->SetBinError(i,0.000001);}
      if(histo_WWqcd_CMS_MVALepResBoundingUp  ->GetBinContent(i) < 0) {histo_WWqcd_CMS_MVALepResBoundingUp  ->SetBinContent(i,0.000001);histo_WWqcd_CMS_MVALepResBoundingUp  ->SetBinError(i,0.000001);}
      if(histo_WWqcd_CMS_MVALepResBoundingDown->GetBinContent(i) < 0) {histo_WWqcd_CMS_MVALepResBoundingDown->SetBinContent(i,0.000001);histo_WWqcd_CMS_MVALepResBoundingDown->SetBinError(i,0.000001);}
      if(histo_WWqcd_CMS_MVAMETResBoundingUp  ->GetBinContent(i) < 0) {histo_WWqcd_CMS_MVAMETResBoundingUp  ->SetBinContent(i,0.000001);histo_WWqcd_CMS_MVAMETResBoundingUp  ->SetBinError(i,0.000001);}
    }
    histo_WWqcd			         ->Scale(nOldWWQCD[0]/histo_WWqcd			         ->GetSumOfWeights());
    histo_WWqcd_CMS_MVALepEffBoundingUp  ->Scale(nOldWWQCD[1]/histo_WWqcd_CMS_MVALepEffBoundingUp  ->GetSumOfWeights());
    histo_WWqcd_CMS_MVALepEffBoundingDown->Scale(nOldWWQCD[2]/histo_WWqcd_CMS_MVALepEffBoundingDown->GetSumOfWeights());
    histo_WWqcd_CMS_MVAJESBoundingUp     ->Scale(nOldWWQCD[3]/histo_WWqcd_CMS_MVAJESBoundingUp     ->GetSumOfWeights());
    histo_WWqcd_CMS_MVAJESBoundingDown   ->Scale(nOldWWQCD[4]/histo_WWqcd_CMS_MVAJESBoundingDown   ->GetSumOfWeights());
    histo_WWqcd_CMS_MVALepResBoundingUp  ->Scale(nOldWWQCD[5]/histo_WWqcd_CMS_MVALepResBoundingUp  ->GetSumOfWeights());
    histo_WWqcd_CMS_MVALepResBoundingDown->Scale(nOldWWQCD[6]/histo_WWqcd_CMS_MVALepResBoundingDown->GetSumOfWeights());
    histo_WWqcd_CMS_MVAMETResBoundingUp  ->Scale(nOldWWQCD[7]/histo_WWqcd_CMS_MVAMETResBoundingUp  ->GetSumOfWeights());
  } else {
    histo_WWqcd			         ->Scale(0.0);
  }

  for(int i=1; i<=histo_WWewk_ALT->GetNbinsX(); i++){
    double factorUp = +1.0; double factorDown = -1.0;
    histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_WWewk_ALT->GetBinContent(i)+factorUp  *histo_WWewk_ALT	 ->GetBinError(i),0.000001));
    histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingDown->SetBinContent(i,TMath::Max(histo_WWewk_ALT->GetBinContent(i)+factorDown*histo_WWewk_ALT	 ->GetBinError(i),0.000001));
    histo_WWewk_CMS_MVAWWewkStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_WWewk    ->GetBinContent(i)+factorUp  *histo_WWewk	 ->GetBinError(i),0.000001));
    histo_WWewk_CMS_MVAWWewkStatBoundingDown        ->SetBinContent(i,TMath::Max(histo_WWewk    ->GetBinContent(i)+factorDown*histo_WWewk	 ->GetBinError(i),0.000001));
    histo_WWqcd_CMS_MVAWWqcdStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_WWqcd    ->GetBinContent(i)+factorUp  *histo_WWqcd	 ->GetBinError(i),0.000001));
    histo_WWqcd_CMS_MVAWWqcdStatBoundingDown        ->SetBinContent(i,TMath::Max(histo_WWqcd    ->GetBinContent(i)+factorDown*histo_WWqcd	 ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingUp	      	    ->SetBinContent(i,TMath::Max(histo_WZ    	->GetBinContent(i)+factorUp  *histo_WZ   	 ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingDown        	    ->SetBinContent(i,TMath::Max(histo_WZ    	->GetBinContent(i)+factorDown*histo_WZ   	 ->GetBinError(i),0.000001));
    histo_WS_CMS_MVAWSStatBoundingUp	      	    ->SetBinContent(i,TMath::Max(histo_WS    	->GetBinContent(i)+factorUp  *histo_WS   	 ->GetBinError(i),0.000001));
    histo_WS_CMS_MVAWSStatBoundingDown        	    ->SetBinContent(i,TMath::Max(histo_WS    	->GetBinContent(i)+factorDown*histo_WS   	 ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingUp        	    ->SetBinContent(i,TMath::Max(histo_VVV   	->GetBinContent(i)+factorUp  *histo_VVV  	 ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown      	    ->SetBinContent(i,TMath::Max(histo_VVV   	->GetBinContent(i)+factorDown*histo_VVV  	 ->GetBinError(i),0.000001));
    histo_Wjets_CMS_MVAWjetsStatBoundingUp    	    ->SetBinContent(i,TMath::Max(histo_Wjets    ->GetBinContent(i)+factorUp  *histo_Wjets        ->GetBinError(i),0.000001));
    histo_Wjets_CMS_MVAWjetsStatBoundingDown  	    ->SetBinContent(i,TMath::Max(histo_Wjets    ->GetBinContent(i)+factorDown*histo_Wjets        ->GetBinError(i),0.000001));

    histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingBinUp[i-1]  ->Add(histo_WWewk_ALT   ); histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingBinUp[i-1]  ->SetBinContent(i,TMath::Max(histo_WWewk_ALT   ->GetBinContent(i)+factorUp  *histo_WWewk_ALT->GetBinError(i),0.000001));
    histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingBinDown[i-1]->Add(histo_WWewk_ALT   ); histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingBinDown[i-1]->SetBinContent(i,TMath::Max(histo_WWewk_ALT   ->GetBinContent(i)+factorDown*histo_WWewk_ALT->GetBinError(i),0.000001));
    histo_WWewk_CMS_MVAWWewkStatBoundingBinUp[i-1]	    ->Add(histo_WWewk       ); histo_WWewk_CMS_MVAWWewkStatBoundingBinUp[i-1]          ->SetBinContent(i,TMath::Max(histo_WWewk       ->GetBinContent(i)+factorUp  *histo_WWewk    ->GetBinError(i),0.000001));
    histo_WWewk_CMS_MVAWWewkStatBoundingBinDown[i-1]	    ->Add(histo_WWewk       ); histo_WWewk_CMS_MVAWWewkStatBoundingBinDown[i-1]        ->SetBinContent(i,TMath::Max(histo_WWewk       ->GetBinContent(i)+factorDown*histo_WWewk    ->GetBinError(i),0.000001));
    histo_WWqcd_CMS_MVAWWqcdStatBoundingBinUp[i-1]	    ->Add(histo_WWqcd       ); histo_WWqcd_CMS_MVAWWqcdStatBoundingBinUp[i-1]          ->SetBinContent(i,TMath::Max(histo_WWqcd       ->GetBinContent(i)+factorUp  *histo_WWqcd    ->GetBinError(i),0.000001));
    histo_WWqcd_CMS_MVAWWqcdStatBoundingBinDown[i-1]	    ->Add(histo_WWqcd       ); histo_WWqcd_CMS_MVAWWqcdStatBoundingBinDown[i-1]        ->SetBinContent(i,TMath::Max(histo_WWqcd       ->GetBinContent(i)+factorDown*histo_WWqcd    ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]	            ->Add(histo_WZ          ); histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]                ->SetBinContent(i,TMath::Max(histo_WZ          ->GetBinContent(i)+factorUp  *histo_WZ       ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]	            ->Add(histo_WZ          ); histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]              ->SetBinContent(i,TMath::Max(histo_WZ          ->GetBinContent(i)+factorDown*histo_WZ       ->GetBinError(i),0.000001));
    histo_WS_CMS_MVAWSStatBoundingBinUp[i-1]	            ->Add(histo_WS          ); histo_WS_CMS_MVAWSStatBoundingBinUp[i-1]                ->SetBinContent(i,TMath::Max(histo_WS          ->GetBinContent(i)+factorUp  *histo_WS       ->GetBinError(i),0.000001));
    histo_WS_CMS_MVAWSStatBoundingBinDown[i-1]              ->Add(histo_WS          ); histo_WS_CMS_MVAWSStatBoundingBinDown[i-1]              ->SetBinContent(i,TMath::Max(histo_WS          ->GetBinContent(i)+factorDown*histo_WS       ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]	            ->Add(histo_VVV         ); histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]              ->SetBinContent(i,TMath::Max(histo_VVV         ->GetBinContent(i)+factorUp  *histo_VVV      ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]            ->Add(histo_VVV         ); histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]            ->SetBinContent(i,TMath::Max(histo_VVV         ->GetBinContent(i)+factorDown*histo_VVV      ->GetBinError(i),0.000001));
    histo_Wjets_CMS_MVAWjetsStatBoundingBinUp[i-1]          ->Add(histo_Wjets       ); histo_Wjets_CMS_MVAWjetsStatBoundingBinUp[i-1]          ->SetBinContent(i,TMath::Max(histo_Wjets       ->GetBinContent(i)+factorUp  *histo_Wjets    ->GetBinError(i),0.000001));
    histo_Wjets_CMS_MVAWjetsStatBoundingBinDown[i-1]        ->Add(histo_Wjets       ); histo_Wjets_CMS_MVAWjetsStatBoundingBinDown[i-1]        ->SetBinContent(i,TMath::Max(histo_Wjets       ->GetBinContent(i)+factorDown*histo_Wjets    ->GetBinError(i),0.000001));
  }
  double mean,up,diff;

  if(showSignalOnly == false) {
    printf("nuisance WZ | Wj: %f/%f | %f/%f\n",histo_WZ   ->GetSumOfWeights(),histo_WZ_CMS_WZNLOBoundingUp  ->GetSumOfWeights(),
                                               histo_Wjets->GetSumOfWeights(),histo_Wjets_CMS_MVAWBoundingUp->GetSumOfWeights());
  }
  histo_Wjets_CMS_MVAWBoundingUp->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_CMS_MVAWBoundingUp->GetSumOfWeights());
  histo_WZ_CMS_WZNLOBoundingUp  ->Scale(histo_WZ   ->GetSumOfWeights()/histo_WZ_CMS_WZNLOBoundingUp  ->GetSumOfWeights());

  for(int i=1; i<=histo_WWewk_ALT->GetNbinsX(); i++){
    
    mean = histo_WWewk_ALT			   ->GetBinContent(i);
    up   = histo_WWewk_ALT_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WWewk_ALT_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WWewk_ALT_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WWewk			   ->GetBinContent(i);
    up   = histo_WWewk_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WWewk_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WWewk_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WWqcd			   ->GetBinContent(i);
    up   = histo_WWqcd_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WWqcd_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WWqcd_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WZ			   ->GetBinContent(i);
    up   = histo_WZ_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_WZ		       ->GetBinContent(i);
    up   = histo_WZ_CMS_WZNLOBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WZ_CMS_WZNLOBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WZ_CMS_WZNLOBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));

    mean = histo_WS			    ->GetBinContent(i);
    up   = histo_WS_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WS_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WS_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_VVV			    ->GetBinContent(i);
    up   = histo_VVV_CMS_MVAMETResBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_VVV_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_VVV_CMS_MVAMETResBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_Wjets 		         ->GetBinContent(i);
    up   = histo_Wjets_CMS_MVAWBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_Wjets_CMS_MVAWBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_Wjets_CMS_MVAWBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  
    mean = histo_WS 		       ->GetBinContent(i);
    up   = histo_WS_CMS_MVAWSBoundingUp->GetBinContent(i);
    diff = TMath::Abs(mean-up);
    if     (mean-up >0) histo_WS_CMS_MVAWSBoundingDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
    else		histo_WS_CMS_MVAWSBoundingDown->SetBinContent(i,TMath::Max(mean-diff,0.000001));
  }
  histo_Wjets_CMS_MVAWBoundingDown->Scale(histo_Wjets->GetSumOfWeights()/histo_Wjets_CMS_MVAWBoundingDown->GetSumOfWeights());
  histo_WZ_CMS_WZNLOBoundingDown  ->Scale(histo_WZ   ->GetSumOfWeights()/histo_WZ_CMS_WZNLOBoundingDown  ->GetSumOfWeights());

  //----------------------------------------------------------------------------
  // Produce output cards for shape-based analyses
  //----------------------------------------------------------------------------
  if(showSignalOnly == false){
  char outputLimits[200];
  sprintf(outputLimits,"wwss%2s.input_%4s.root",finalStateName,ECMsb.Data());
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data	 ->Write();
  histo_WWewk_ALT->Write();
  histo_WWewk	 ->Write();
  histo_WWqcd	 ->Write();
  histo_WZ	 ->Write();
  histo_WS	 ->Write();
  histo_VVV	 ->Write();
  histo_Wjets	 ->Write();

  cout << histo_Data	 ->GetSumOfWeights() << " ";
  cout << histo_WWewk_ALT->GetSumOfWeights() << " ";
  cout << histo_WWewk	 ->GetSumOfWeights() << " ";
  cout << histo_WWqcd	 ->GetSumOfWeights() << " ";
  cout << histo_WZ	 ->GetSumOfWeights() << " ";
  cout << histo_WS	 ->GetSumOfWeights() << " ";
  cout << histo_VVV	 ->GetSumOfWeights() << " ";
  cout << histo_Wjets	 ->GetSumOfWeights() << " ";
  cout << endl;
  printf("uncertainties Stat\n");
  histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingUp  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk_ALT     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingUp 	 ->GetBinContent(i)/histo_WWewk_ALT   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingDown->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk_ALT     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingDown      ->GetBinContent(i)/histo_WWewk_ALT   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_CMS_MVAWWewkStatBoundingUp	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk	->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_CMS_MVAWWewkStatBoundingUp	     ->GetBinContent(i)/histo_WWewk   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_CMS_MVAWWewkStatBoundingDown	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk	->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_CMS_MVAWWewkStatBoundingDown      ->GetBinContent(i)/histo_WWewk   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_CMS_MVAWWqcdStatBoundingUp	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_CMS_MVAWWqcdStatBoundingUp  	->GetBinContent(i)/histo_WWqcd   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_CMS_MVAWWqcdStatBoundingDown	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_CMS_MVAWWqcdStatBoundingDown      ->GetBinContent(i)/histo_WWqcd   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingUp	  	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingUp	->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingDown	 	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingDown	->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_CMS_MVAWSStatBoundingUp	  	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAWSStatBoundingUp	->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_CMS_MVAWSStatBoundingDown	 	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAWSStatBoundingDown	->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingUp	  	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingDown   	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wjets_CMS_MVAWjetsStatBoundingUp  	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_CMS_MVAWjetsStatBoundingUp  ->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wjets_CMS_MVAWjetsStatBoundingDown	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_CMS_MVAWjetsStatBoundingDown->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEff\n");
  histo_WWewk_ALT_CMS_MVALepEffBoundingUp    ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk_ALT  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_ALT_CMS_MVALepEffBoundingUp      ->GetBinContent(i)/histo_WWewk_ALT    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_ALT_CMS_MVALepEffBoundingDown  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk_ALT  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_ALT_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_WWewk_ALT     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_CMS_MVALepEffBoundingUp        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_CMS_MVALepEffBoundingUp	  ->GetBinContent(i)/histo_WWewk	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_CMS_MVALepEffBoundingDown      ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_WWewk	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_CMS_MVALepEffBoundingUp        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_CMS_MVALepEffBoundingUp	     ->GetBinContent(i)/histo_WWqcd	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_CMS_MVALepEffBoundingDown      ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_CMS_MVALepEffBoundingDown	->GetBinContent(i)/histo_WWqcd     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffBoundingUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffBoundingUp     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffBoundingDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_CMS_MVALepEffBoundingUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVALepEffBoundingUp     ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_CMS_MVALepEffBoundingDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVALepEffBoundingDown   ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffBoundingUp          ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffBoundingUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffBoundingDown        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LetRes\n");
  histo_WWewk_ALT_CMS_MVALepResBoundingUp    ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk_ALT  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_ALT_CMS_MVALepResBoundingUp      ->GetBinContent(i)/histo_WWewk_ALT    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_ALT_CMS_MVALepResBoundingDown  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk_ALT  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_ALT_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_WWewk_ALT     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_CMS_MVALepResBoundingUp        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_CMS_MVALepResBoundingUp	  ->GetBinContent(i)/histo_WWewk	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_CMS_MVALepResBoundingDown      ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_WWewk	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_CMS_MVALepResBoundingUp        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_CMS_MVALepResBoundingUp	     ->GetBinContent(i)/histo_WWqcd	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_CMS_MVALepResBoundingDown      ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_CMS_MVALepResBoundingDown	->GetBinContent(i)/histo_WWqcd     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepResBoundingUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepResBoundingUp     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepResBoundingDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_CMS_MVALepResBoundingUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVALepResBoundingUp     ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_CMS_MVALepResBoundingDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVALepResBoundingDown   ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepResBoundingUp          ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepResBoundingUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepResBoundingDown        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepResBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties METRes\n");
  histo_WWewk_ALT_CMS_MVAMETResBoundingUp    ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk_ALT	->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_ALT_CMS_MVAMETResBoundingUp	  ->GetBinContent(i)/histo_WWewk_ALT	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_ALT_CMS_MVAMETResBoundingDown  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk_ALT	->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_ALT_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_WWewk_ALT	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_CMS_MVAMETResBoundingUp        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_CMS_MVAMETResBoundingUp	     ->GetBinContent(i)/histo_WWewk	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_CMS_MVAMETResBoundingDown      ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_CMS_MVAMETResBoundingDown	->GetBinContent(i)/histo_WWewk     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_CMS_MVAMETResBoundingUp        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_CMS_MVAMETResBoundingUp	     ->GetBinContent(i)/histo_WWqcd	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_CMS_MVAMETResBoundingDown      ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_CMS_MVAMETResBoundingDown	->GetBinContent(i)/histo_WWqcd     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETResBoundingUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETResBoundingUp     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETResBoundingDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_CMS_MVAMETResBoundingUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAMETResBoundingUp     ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_CMS_MVAMETResBoundingDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAMETResBoundingDown   ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETResBoundingUp          ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETResBoundingUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETResBoundingDown        ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETResBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties JES\n");
  histo_WWewk_ALT_CMS_MVAJESBoundingUp       ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk_ALT  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_ALT_CMS_MVAJESBoundingUp	      ->GetBinContent(i)/histo_WWewk_ALT    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_ALT_CMS_MVAJESBoundingDown     ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk_ALT  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_ALT_CMS_MVAJESBoundingDown       ->GetBinContent(i)/histo_WWewk_ALT    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_CMS_MVAJESBoundingUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_CMS_MVAJESBoundingUp	  ->GetBinContent(i)/histo_WWewk	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWewk_CMS_MVAJESBoundingDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWewk  ->GetBinContent(i)>0)printf("%5.1f ",histo_WWewk_CMS_MVAJESBoundingDown   ->GetBinContent(i)/histo_WWewk	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_CMS_MVAJESBoundingUp           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_CMS_MVAJESBoundingUp     ->GetBinContent(i)/histo_WWqcd	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WWqcd_CMS_MVAJESBoundingDown         ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WWqcd     ->GetBinContent(i)>0)printf("%5.1f ",histo_WWqcd_CMS_MVAJESBoundingDown   ->GetBinContent(i)/histo_WWqcd	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAJESBoundingUp              ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingUp	     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAJESBoundingDown            ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingDown      ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_CMS_MVAJESBoundingUp              ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAJESBoundingUp	     ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_CMS_MVAJESBoundingDown            ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS	   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAJESBoundingDown      ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingUp             ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingDown           ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingDown	     ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties GEN\n");
  histo_Wjets_CMS_MVAWBoundingUp	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_CMS_MVAWBoundingUp       ->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Wjets_CMS_MVAWBoundingDown	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_Wjets->GetBinContent(i)>0)printf("%5.1f ",histo_Wjets_CMS_MVAWBoundingDown     ->GetBinContent(i)/histo_Wjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_WZNLOBoundingUp            ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_WZNLOBoundingUp	     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_WZNLOBoundingDown          ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_WZNLOBoundingDown       ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_CMS_MVAWSBoundingUp            ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAWSBoundingUp         ->GetBinContent(i)/histo_WS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WS_CMS_MVAWSBoundingDown	  ->Write(); for(int i=1; i<=histo_WWewk->GetNbinsX(); i++) {if(histo_WS   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAWSBoundingDown       ->GetBinContent(i)/histo_WS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");

  for(int nb=0; nb<nBinMVA; nb++){
    histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingBinUp[nb]   ->Write();
    histo_WWewk_ALT_CMS_MVAWWewk_ALTStatBoundingBinDown[nb] ->Write();
    histo_WWewk_CMS_MVAWWewkStatBoundingBinUp[nb]	    ->Write();
    histo_WWewk_CMS_MVAWWewkStatBoundingBinDown[nb]	    ->Write();
    histo_WWqcd_CMS_MVAWWqcdStatBoundingBinUp[nb]	    ->Write();
    histo_WWqcd_CMS_MVAWWqcdStatBoundingBinDown[nb]	    ->Write();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	   	    ->Write();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	   	    ->Write();
    histo_WS_CMS_MVAWSStatBoundingBinUp[nb]	   	    ->Write();
    histo_WS_CMS_MVAWSStatBoundingBinDown[nb]      	    ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	   	    ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	   	    ->Write();
    histo_Wjets_CMS_MVAWjetsStatBoundingBinUp[nb]  	    ->Write();
    histo_Wjets_CMS_MVAWjetsStatBoundingBinDown[nb]         ->Write();
  }

  // HACK!, we don't want alternative models now
  NFinal[nBkg] = 0.0;
  char theWWewk_ALTString[20];
  if(NFinal[nBkg] > 0) sprintf(theWWewk_ALTString,"1.000");
  else                 sprintf(theWWewk_ALTString,"  -  ");

  char outputLimitsShape[200];
  sprintf(outputLimitsShape,"histo_limits_wwss%2s_shape_%4s.txt",finalStateName,ECMsb.Data());
  ofstream newcardShape;
  newcardShape.open(outputLimitsShape);
  newcardShape << Form("imax 1 number of channels\n");
  newcardShape << Form("jmax * number of background\n");
  newcardShape << Form("kmax * number of nuisance parameters\n");
  newcardShape << Form("Observation %d\n",(int)nSelectedData[WWSEL+nSelTypes]);
  newcardShape << Form("shapes *   *   %s  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n",outputLimits);
  newcardShape << Form("shapes data_obs * %s  histo_Data \n",outputLimits);
  newcardShape << Form("bin wwss%2s%4s wwss%2s%4s wwss%2s%4s wwss%2s%4s wwss%2s%4s wwss%2s%4s wwss%2s%4s\n",finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data(),finalStateName,ECMsb.Data());
  newcardShape << Form("process WWewk_ALT WWewk WWqcd WZ WS VVV Wjets\n");
  newcardShape << Form("process -1 0 1 2 3 4 5\n");
  newcardShape << Form("rate %6.3f %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",NFinal[nBkg],NFinal[0],NFinal[1],NFinal[2],NFinal[3],NFinal[4],TMath::Max(NFinal[5],0.0));
  newcardShape << Form("lumi_%4s                               lnN  %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f   -  \n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE,lumiE,lumiE); 		     
  newcardShape << Form("%s                                   shape   %s   1.000 1.000 1.000 1.000 1.000   -  \n",effName,theWWewk_ALTString);
  newcardShape << Form("%s                                   shape   %s   1.000 1.000 1.000 1.000 1.000   -  \n",momName,theWWewk_ALTString);
  newcardShape << Form("CMS_scale_met                        shape   %s   1.000 1.000 1.000 1.000 1.000   -  \n",theWWewk_ALTString);
  newcardShape << Form("CMS_scale_j                          shape   %s   1.000 1.000 1.000 1.000 1.000   -  \n",theWWewk_ALTString);			 
  newcardShape << Form("pdf_qqbar                              lnN  %5.3f %5.3f %5.3f %5.3f   -     -     -  \n",pdf_qqbar[0],pdf_qqbar[1],pdf_qqbar[2],pdf_qqbar[3]);
  newcardShape << Form("QCDscale_WWEwk		               lnN  %5.3f %5.3f   -     -     -     -     -  \n",QCDscale_WWewk,QCDscale_WWewk);	
  newcardShape << Form("QCDscale_VV		               lnN    -     -     -   1.100   -     -     -  \n");		
  newcardShape << Form("CMS_wwss_WZ3l                          lnN    -     -     -   %5.3f   -     -     -  \n",syst_WZ3l);		
  if(NFinal[3] > 0)
  newcardShape << Form("CMS_wwss_WZNLOBounding               shape    -     -     -   1.000   -     -     -  \n");
  if(NFinal[4] > 0)
  newcardShape << Form("CMS_wwss_MVAWSBounding               shape    -     -     -     -   1.000   -     -  \n");		
  if(NFinal[4] > 0)
  newcardShape << Form("QCDscale_VVV		               lnN    -     -     -     -     -   1.500   -  \n");		
  if(NFinal[5] > 0)
  newcardShape << Form("CMS_FakeRate                           lnN    -     -     -     -     -	    -   %5.3f\n",WjetsSyst);  
  if(NFinal[5] > 0)
  newcardShape << Form("CMS_wwss_MVAWBounding                shape    -     -     -	-     -     -   1.000\n");
  if(useFullStatTemplates == false){
    if(NFinal[nBkg] > 0) newcardShape << Form("CMS_wwss%s_MVAWWewk_ALTStatBounding_%s  shape  1.000   -     -	  -	-     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[0]    > 0) newcardShape << Form("CMS_wwss%s_MVAWWewkStatBounding_%s      shape    -   1.000   -	  -	-     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[3]    > 0) newcardShape << Form("CMS_wwss%s_MVAWWqcdStatBounding_%s      shape    -     -   1.000   -	-     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[2]    > 0) newcardShape << Form("CMS_wwss%s_MVAWZStatBounding_%s	       shape    -     -     -	1.000	-     -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[4]    > 0) newcardShape << Form("CMS_wwss%s_MVAWSStatBounding_%s	       shape    -     -     -	  -   1.000   -     -  \n",finalStateName,ECMsb.Data());
    if(NFinal[1]    > 0) newcardShape << Form("CMS_wwss%s_MVAVVVStatBounding_%s        shape    -     -     -	  -	-   1.000   -  \n",finalStateName,ECMsb.Data());
    if(NFinal[5]    > 0) newcardShape << Form("CMS_wwss%s_MVAWjetsStatBounding_%s      shape    -     -     -	  -	-     -   1.000\n",finalStateName,ECMsb.Data());
  } else {
    for(int nb=1; nb<=nBinMVA; nb++){
      if(NFinal[nBkg] > 0 && histo_WWewk_ALT->GetBinContent(nb)    > 0) newcardShape << Form("CMS_wwss%s_MVAWWewk_ALTStatBounding_%s_Bin%d  shape  1.000   -     -     -     -     -	 -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[0]    > 0 && histo_WWewk->GetBinContent(nb)    > 0) newcardShape << Form("CMS_wwss%s_MVAWWewkStatBounding_%s_Bin%d          shape    -   1.000   -     -     -     -	     -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[1]    > 0 && histo_WWqcd->GetBinContent(nb)   > 0) newcardShape << Form("CMS_wwss%s_MVAWWqcdStatBounding_%s_Bin%d           shape    -     -   1.000   -     -     -	  -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[2]    > 0 && histo_WZ->GetBinContent(nb)    > 0) newcardShape << Form("CMS_wwss%s_MVAWZStatBounding_%s_Bin%d                shape    -     -	 -   1.000   -     -	  -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[3]    > 0 && histo_WS->GetBinContent(nb)    > 0) newcardShape << Form("CMS_wwss%s_MVAWSStatBounding_%s_Bin%d                shape    -     -	 -     -   1.000   -	  -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[4]    > 0 && histo_VVV->GetBinContent(nb)    > 0) newcardShape << Form("CMS_wwss%s_MVAVVVStatBounding_%s_Bin%d              shape    -     -	 -     -     -   1.000   -  \n",finalStateName,ECMsb.Data(),nb-1);
      if(NFinal[5]    > 0 && histo_Wjets->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwss%s_MVAWjetsStatBounding_%s_Bin%d             shape    -     -	 -     -     -     -   1.000\n",finalStateName,ECMsb.Data(),nb-1);
    }
  }
  newcardShape.close();
  }

  return;
}

void scaleFactor_WS(LorentzVector l,int lq, int ld, int mcld, double val[2]){
//---------------------------------------------------------------------
// |eta|        data                  mc                    factor
//---------------------------------------------------------------------
//0.0-0.5 0.00043 +/- 0.00002 | 0.00039 +/- 0.00002 ==> 1.098 +/- 0.091
//0.5-1.0 0.00070 +/- 0.00003 | 0.00078 +/- 0.00004 ==> 0.891 +/- 0.059
//1.0-1.5 0.00376 +/- 0.00010 | 0.00309 +/- 0.00009 ==> 1.219 +/- 0.050
//1.5-2.0 0.00869 +/- 0.00021 | 0.00664 +/- 0.00018 ==> 1.309 +/- 0.048
//2.0-2.5 0.01203 +/- 0.00030 | 0.00891 +/- 0.00025 ==> 1.351 +/- 0.051
// additional 10% uncertainty for the overall normalization
  double factor[5]  = {1.098,0.891,1.219,1.309,1.351};
  double factorE[5] = {0.091,0.059,0.050,0.048,0.051};

  if(abs(ld) == 11){
    if((mcld ==  11 && lq > 0) || 
       (mcld == -11 && lq < 0)){ // wrong charge
      if     (abs(l.Eta()) >= 0.0 && abs(l.Eta()) < 0.5) {val[0] = val[0]*factor[0]; val[1] = val[1]*(factor[0]+sqrt(factorE[0]*factorE[0]+0.10*0.10));}
      else if(abs(l.Eta()) >= 0.5 && abs(l.Eta()) < 1.0) {val[0] = val[0]*factor[1]; val[1] = val[1]*(factor[1]+sqrt(factorE[1]*factorE[1]+0.10*0.10));}
      else if(abs(l.Eta()) >= 1.0 && abs(l.Eta()) < 1.5) {val[0] = val[0]*factor[2]; val[1] = val[1]*(factor[2]+sqrt(factorE[2]*factorE[2]+0.10*0.10));}
      else if(abs(l.Eta()) >= 1.5 && abs(l.Eta()) < 2.0) {val[0] = val[0]*factor[3]; val[1] = val[1]*(factor[3]+sqrt(factorE[3]*factorE[3]+0.10*0.10));}
      else if(abs(l.Eta()) >= 2.0 && abs(l.Eta()) < 2.5) {val[0] = val[0]*factor[4]; val[1] = val[1]*(factor[4]+sqrt(factorE[4]*factorE[4]+0.10*0.10));}
    }
  }
}
