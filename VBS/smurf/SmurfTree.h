#ifndef SmurfTree_H
#define SmurfTree_H

#include <set>
#include <vector>
#include <string>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"


//
// Ntuple structure:
//  * plain ntuples without vectors
//  * hypothesis based (not event based). Upto to a user to decide if
//    he allows multiple hypothesis per event or not.
//
// Ntuple content:
//  * event info: run, lumi, event, nVertex, weight (for MC), event type (data, ttbar, wjets etc)
//  * hypothesis type (leading, trailing): ee, em, me, mm
//  * lepton1 (leading)
//  * lepton2 (sub-leading)
//  * jet1 (leading)
//  * jet2 (sub-leading)
//  * jetEventType (0,1,VBF)
//  * nJets
//  * dilepton p4
//  * met
//  * metPhi
//  * sumEt
//  * projectedMet
//  * MT - transverse higgs mass
//  * MT1 - transverse mass of lepton1 and met
//  * MT2 - transverse mass of lepton2 and met
//  * dPhi - delta phi between leptons
//  * dR -   delta R   between leptons   
//  * dPhiLep1Jet1 - delta phi between leading jet and leading lepton
//  * dRLep1Jet1 -   delta R   between leading jet and leading lepton
//  * dPhiLep2Jet1 - delta phi between leading jet and sub-leading lepton
//  * dRLep2Jet1 -   delta R   between leading jet and sub-leading lepton
//  * dPhiLep1Met - delta phi between leading lepton and MET
//  * dPhiLep2Met - delta phi between sub-leading lepton and MET
//  * dPhiDiLepMet - delta phi between di-lepton and MET
//  * dPhiDiLepJet1 - delta phi between di-lepton and leading jet
//  * mc info
//    * lep1 mc_id
//    * lep2 mc_id
//    * jet1 hard-scatter match (flavor)
//    * jet2 hard-scatter match (flavor)
//    * gen met
//
// Lepton variables:
//  * p4
//  * charge
//  * selection type
//  
// Jet variables:
//  * p4
//  * b-tag
// 

class SmurfTree {
 public:
  /// float doesn't have dictionary by default, so use double
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

  /// bit map
  /// DON'T CHANGE ORDER
  enum Selection {
    BaseLine           = 1UL<<0,  // pt(reco)>20/10, acceptance,!STA muon, mll>12
    ChargeMatch        = 1UL<<1,  // q1*q2<0
    Lep1FullSelection  = 1UL<<2,  // full id, isolation, d0, dz etc
    Lep1LooseEleV1     = 1UL<<3,  // electron fakeable object selection is passed V1
    Lep1LooseEleV2     = 1UL<<4,  // electron fakeable object selection is passed V2
    Lep1LooseEleV3     = 1UL<<5,  // electron fakeable object selection is passed V3
    Lep1LooseEleV4     = 1UL<<6,  // electron fakeable object selection is passed V4
    Lep1LooseMuV1      = 1UL<<7,  // muon fakeable object selection (relIso<1.0)
    Lep1LooseMuV2      = 1UL<<8,  // muon fakeable object selection (relIso<0.4)
    Lep2FullSelection  = 1UL<<9,  // full id, isolation, d0, dz etc
    Lep2LooseEleV1     = 1UL<<10, // electron fakeable object selection is passed V1
    Lep2LooseEleV2     = 1UL<<11, // electron fakeable object selection is passed V2
    Lep2LooseEleV3     = 1UL<<12, // electron fakeable object selection is passed V3
    Lep2LooseEleV4     = 1UL<<13, // electron fakeable object selection is passed V4
    Lep2LooseMuV1      = 1UL<<14, // muon fakeable object selection (relIso<1.0)
    Lep2LooseMuV2      = 1UL<<15, // muon fakeable object selection (relIso<0.4)
    FullMET            = 1UL<<16, // full met selection
    ZVeto              = 1UL<<17, // event is not in the Z-mass peak for ee/mm final states
    TopTag             = 1UL<<18, // soft muon and b-jet tagging for the whole event regardless of n-jets (non-zero means tagged)
    TopVeto            = 1UL<<19, // soft muon and b-jet tagging for the whole event regardless of n-jets (zero means tagged)
    OneBJet            = 1UL<<20, // 1-jet events, where the jet is b-tagged (top control sample with one b-quark missing)
    TopTagNotInJets    = 1UL<<21, // soft muon and b-jet tagging for areas outside primary jets (non-zero means tagged)
    ExtraLeptonVeto    = 1UL<<22, // extra lepton veto, DR(muon-electron)>=0.3
    Lep3FullSelection  = 1UL<<23,  // full id, isolation, d0, dz etc
    Lep3LooseEleV1     = 1UL<<24, // electron fakeable object selection is passed V1
    Lep3LooseEleV2     = 1UL<<25, // electron fakeable object selection is passed V2
    Lep3LooseEleV3     = 1UL<<26, // electron fakeable object selection is passed V3
    Lep3LooseEleV4     = 1UL<<27, // electron fakeable object selection is passed V4
    Lep3LooseMuV1      = 1UL<<28, // muon fakeable object selection (relIso<1.0)
    Lep3LooseMuV2      = 1UL<<29, // muon fakeable object selection (relIso<0.4)
    Trigger            = 1UL<<30  // passed a set of triggers
  };

  static const unsigned int FullSelection = BaseLine|ChargeMatch|Lep1FullSelection|Lep2FullSelection|FullMET|ZVeto|TopVeto|ExtraLeptonVeto;

  /// first is leading lepton
  /// DON'T CHANGE ORDER
  enum Type {
    mm, 
    me, 
    em, 
    ee
  };

  /// DON'T CHANGE ORDER
  enum DataType {
    data,
    qqww,
    ggww,
    hww120,
    hww130,
    hww140,
    hww150,
    hww160,
    hww170,
    hww180,
    hww190,
    hww200,
    hww210,
    hww220,
    hww230,
    hww250,
    hww300,
    hww350,
    hww400,
    hww450,
    hww500,
    hww550,
    hww600,
    vbfhww120,
    vbfhww130,
    vbfhww140,
    vbfhww150,
    vbfhww160,
    vbfhww170,
    vbfhww180,
    vbfhww190,
    vbfhww200,
    vbfhww210,
    vbfhww220,
    vbfhww230,
    vbfhww250,
    vbfhww300,
    vbfhww350,
    vbfhww400,
    vbfhww450,
    vbfhww500,
    vbfhww550,
    vbfhww600,
    ttbar,
    tw,
    dyee,
    dymm,
    dytt,
    wjets,
    wz,
    zz,
    wgamma,
    qcd,
    other,
    hww110,
    hww115,
    vbfhww110,
    vbfhww115,
    ggzz,
    www,
    dyttDataDriven,
    wgstar,
    hww118,
    hww122,
    hww124,
    hww126,
    hww128,
    hww135,
    vbfhww118,
    vbfhww122,
    vbfhww124,
    vbfhww126,
    vbfhww128,
    vbfhww135,
    hww125,
    hww145,
    vbfhww125,
    vbfhww145,
    hww155,
    hww700,
    hww800,
    hww900,
    hww1000,
    vbfhww155,
    vbfhww700,
    vbfhww800,
    vbfhww900,
    vbfhww1000,
    qqwwPWG,
    qqww2j,
    qqbarh
  };

  /// variables
  unsigned int   event_;
  unsigned int   run_;
  unsigned int   lumi_;
  unsigned int   nvtx_;
  unsigned int   cuts_;
  float          scale1fb_;
  float          CHSMet_;
  float          NHSMet_;
  float          CHSMetPhi_;
  float          NHSMetPhi_;
  float          met_;
  float          metPhi_;
  float          metMVA_;
  float          metMVAPhi_;
  float          pmetMVA_;
  float          sumet_;
  float          metSig_;
  float          genmet_;
  float          genmetPhi_;
  Type           type_;
  DataType       dstype_;
  LorentzVector  lep1_;
  int            lq1_;
  int            lid1_;
  float          lmva1_;
  LorentzVector  lep2_;
  int            lq2_;
  int            lid2_;
  float          lmva2_;
  LorentzVector  jet1_;
  float          jet1Btag_;
  float          jet1ProbBtag_;
  float          jet1Dz_;
  LorentzVector  jet2_;
  float          jet2Btag_;
  float          jet2ProbBtag_;
  unsigned int   njets_;
  LorentzVector  dilep_;
  LorentzVector  quadlep_;
  float          trackMet_;
  float          trackMetPhi_;
  float          pmet_;
  float          pTrackMet_;
  float          mt_;
  float          mt1_;
  float          mt2_;
  float          dPhi_;
  float          dR_;
  float          dPhiLep1Jet1_;
  float          dRLep1Jet1_;
  float          dPhiLep2Jet1_;
  float          dRLep2Jet1_;
  float          dPhiLep1MET_;
  float          dPhiLep2MET_;
  float          dPhiDiLepMET_;
  float          dPhiDiLepJet1_;
  float          lep1DetEta_;
  float          lep2DetEta_;
  int            lep1McId_;
  int            lep2McId_;
  int            lep1MotherMcId_;
  int            lep2MotherMcId_;
  int            jet1McId_;
  int            jet2McId_;
  TNamed         info_;

  LorentzVector  lep3_;
  int            lq3_;
  int            lid3_;
  float          lmva3_;
  LorentzVector  jet3_;
  float          jet3Btag_;
  float          jet3ProbBtag_;
  LorentzVector  jet4_;
  float          jet4Btag_;
  float          jet4ProbBtag_;
  int            jet3McId_;
  int            jet4McId_;
  float          lep3DetEta_;
  int            lep3McId_;
  int            lep3MotherMcId_;
  float          dPhiLep3Jet1_;
  float          dRLep3Jet1_;
  float          dPhiLep3MET_;
  float          mt3_;
  float          jetLowBtag_;
  unsigned int   nSoftMuons_;
  float          Q_;
  float          id1_;
  float          x1_;
  float          pdf1_;
  float          id2_;
  float          x2_;
  float          pdf2_;
  int            processId_;
  float          higgsPt_;
  float          sfWeightFR_;
  float          sfWeightPU_;
  float          sfWeightTrig_;
  float          sfWeightEff_;
  float          sfWeightHPt_;
  float          dymva_;
  float          npu_;
  float          npuPlusOne_;
  float          npuMinusOne_;

  float          auxVar0_;

 public:
  /// this is the main element
  TTree *tree_;
  TFile *f_;
  
  /// hold the names of variables to facilitate things (filled during Init)
  std::vector<std::string> variables_;

  /// default constructor  
  SmurfTree():info_("info","Smurf ntuple"),
    lepPtr1_(&lep1_),lepPtr2_(&lep2_),jetPtr1_(&jet1_),jetPtr2_(&jet2_),dilepPtr_(&dilep_),quadlepPtr_(&quadlep_),
    lepPtr3_(&lep3_),                 jetPtr3_(&jet3_),jetPtr4_(&jet4_){}
  /// default destructor
  ~SmurfTree(){ 
    if (f_) f_->Close();  
  };

  /// initialize varibles and fill list of available variables
  void InitVariables();

  /// load a SmurfTree
  void LoadTree(const char* file, int type = -1){
    // to load three different ntuples in the same job HwwTree0/1/2/3/4
    // type == 0/1/2/3/4 if all variables was added
    // type = -1 (default) if a minimum set of variables was added with tree as name
    f_ = TFile::Open(file);
    assert(f_);
    if     (type == 0) tree_ = dynamic_cast<TTree*>(f_->Get("HwwTree0"));
    else if(type == 1) tree_ = dynamic_cast<TTree*>(f_->Get("HwwTree1"));
    else if(type == 2) tree_ = dynamic_cast<TTree*>(f_->Get("HwwTree2"));
    else if(type == 3) tree_ = dynamic_cast<TTree*>(f_->Get("HwwTree3"));
    else if(type == 4) tree_ = dynamic_cast<TTree*>(f_->Get("HwwTree4"));
    else               tree_ = dynamic_cast<TTree*>(f_->Get("tree"));
    assert(tree_);
  }

  /// create a SmurfTree
  void CreateTree(int type = -1){
    assert(type==type); // just to suppress warnings
    // to create three different ntuples in the same job HwwTree0/1/2/3/4
    // type == 0/1/2/3/4 add all variables
    // type = -1 (default) add a minimum set of variables with tree as name
    if     (type == 0) tree_ = new TTree("HwwTree0","Smurf ntuples");
    else if(type == 1) tree_ = new TTree("HwwTree1","Smurf ntuples");
    else if(type == 2) tree_ = new TTree("HwwTree2","Smurf ntuples");
    else if(type == 3) tree_ = new TTree("HwwTree3","Smurf ntuples");
    else if(type == 4) tree_ = new TTree("HwwTree4","Smurf ntuples");
    else               tree_ = new TTree("tree","Smurf ntuples");
    f_ = 0;
    InitVariables();
    //book the branches
    tree_->Branch("event"        , &event_        ,   "event/i");
    tree_->Branch("run"          , &run_          ,   "run/i");
    tree_->Branch("lumi"         , &lumi_         ,   "lumi/i");
    tree_->Branch("nvtx"         , &nvtx_         ,   "nvtx/i");
    tree_->Branch("cuts"         , &cuts_         ,   "cuts/i");
    tree_->Branch("scale1fb"     , &scale1fb_     ,   "scale1fb/F");
    tree_->Branch("CHSMet"       , &CHSMet_	  ,   "CHSMet/F");
    tree_->Branch("NHSMet"       , &NHSMet_	  ,   "NHSMet/F");
    tree_->Branch("CHSMetPhi"    , &CHSMetPhi_    ,   "CHSMetPhi/F");
    tree_->Branch("NHSMetPhi"    , &NHSMetPhi_    ,   "NHSMetPhi/F");
    tree_->Branch("met"          , &met_          ,   "met/F");
    tree_->Branch("metPhi"       , &metPhi_       ,   "metPhi/F");
    tree_->Branch("metMVA"       , &metMVA_	  ,   "metMVA/F");
    tree_->Branch("metMVAPhi"    , &metMVAPhi_    ,   "metMVAPhi/F");
    tree_->Branch("pmetMVA"      , &pmetMVA_      ,   "pmetMVA/F");
    tree_->Branch("sumet"        , &sumet_        ,   "sumet/F");
    tree_->Branch("metSig"       , &metSig_       ,   "metSig/F");
    tree_->Branch("genmet"       , &genmet_       ,   "genmet/F");
    tree_->Branch("genmetPhi"    , &genmetPhi_    ,   "genmetPhi/F");
    tree_->Branch("type"         , &type_         ,   "type/I");
    tree_->Branch("dstype"       , &dstype_       ,   "dstype/I");
    tree_->Branch("lep1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr1_);
    tree_->Branch("lq1"          , &lq1_          ,   "lq1/I");
    tree_->Branch("lid1"         , &lid1_         ,   "lid1/I");
    tree_->Branch("lmva1"        , &lmva1_        ,   "lmva1/F");
    tree_->Branch("lep2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr2_);
    tree_->Branch("lq2"          , &lq2_          ,   "lq2/I");
    tree_->Branch("lid2"         , &lid2_         ,   "lid2/I");
    tree_->Branch("lmva2"        , &lmva2_        ,   "lmva2/F");
    tree_->Branch("jet1"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr1_);
    tree_->Branch("jet1Btag"     , &jet1Btag_     ,   "jet1Btag/F");
    tree_->Branch("jet1ProbBtag" , &jet1ProbBtag_ ,   "jet1ProbBtag/F");
    tree_->Branch("jet1Dz"       , &jet1Dz_       ,   "jet1Dz/F");
    tree_->Branch("jet2"         , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr2_);
    tree_->Branch("jet2Btag"     , &jet2Btag_     ,   "jet2Btag/F");
    tree_->Branch("jet2ProbBtag" , &jet2ProbBtag_ ,   "jet2ProbBtag/F");
    tree_->Branch("njets"        , &njets_        ,   "njets/i");
    tree_->Branch("dilep"        , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &dilepPtr_);
    tree_->Branch("quadlep"      , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &quadlepPtr_);
    tree_->Branch("trackMet"     , &trackMet_	  ,   "trackMet/F");
    tree_->Branch("trackMetPhi"  , &trackMetPhi_  ,   "trackMetPhi/F");
    tree_->Branch("pmet"         , &pmet_         ,   "pmet/F");
    tree_->Branch("pTrackMet"    , &pTrackMet_    ,   "pTrackMet/F");
    tree_->Branch("mt"           , &mt_           ,   "mt/F");
    tree_->Branch("mt1"          , &mt1_          ,   "mt1/F");
    tree_->Branch("mt2"          , &mt2_          ,   "mt2/F");
    tree_->Branch("dPhi"         , &dPhi_         ,   "dPhi/F");
    tree_->Branch("dR"           , &dR_           ,   "dR/F");
    tree_->Branch("dPhiLep1Jet1" , &dPhiLep1Jet1_ ,   "dPhiLep1Jet1/F");
    tree_->Branch("dRLep1Jet1"   , &dRLep1Jet1_   ,    "dRLep1Jet1/F");
    tree_->Branch("dPhiLep2Jet1" , &dPhiLep2Jet1_ ,   "dPhiLep2Jet1/F");
    tree_->Branch("dRLep2Jet1"   , &dRLep2Jet1_   ,   "dRLep2Jet1/F");
    tree_->Branch("dPhiDiLepMET" , &dPhiDiLepMET_ ,   "dPhiDiLepMET/F");
    tree_->Branch("dPhiLep1MET"  , &dPhiLep1MET_  ,   "dPhiLep1MET/F");
    tree_->Branch("dPhiLep2MET"  , &dPhiLep2MET_  ,   "dPhiLep2MET/F");
    tree_->Branch("dPhiDiLepJet1", &dPhiDiLepJet1_,   "dPhiDiLepJet1/F");
    tree_->Branch("lep1DetEta"	 , &lep1DetEta_   ,   "lep1DetEta/F");
    tree_->Branch("lep2DetEta"	 , &lep2DetEta_   ,   "lep2DetEta/F");
    tree_->Branch("lep1McId"     , &lep1McId_     ,   "lep1McId/I");
    tree_->Branch("lep2McId"     , &lep2McId_     ,   "lep2McId/I");
    tree_->Branch("lep1MotherMcId",&lep1MotherMcId_,  "lep1MotherMcId/I");
    tree_->Branch("lep2MotherMcId",&lep2MotherMcId_,  "lep2MotherMcId/I");
    tree_->Branch("jet1McId"      ,&jet1McId_	  ,   "jet1McId/I");
    tree_->Branch("jet2McId"      ,&jet2McId_	  ,   "jet2McId/I");
    tree_->Branch("processId",     &processId_    ,   "processId/I");
    tree_->Branch("higgsPt",       &higgsPt_      ,   "higgsPt/F");
    tree_->Branch("sfWeightFR",    &sfWeightFR_   ,   "sfWeightFR/F");
    tree_->Branch("sfWeightPU",    &sfWeightPU_   ,   "sfWeightPU/F");
    tree_->Branch("sfWeightTrig",  &sfWeightTrig_ ,   "sfWeightTrig/F");
    tree_->Branch("sfWeightEff",   &sfWeightEff_  ,   "sfWeightEff/F");
    tree_->Branch("sfWeightHPt",   &sfWeightHPt_  ,   "sfWeightHPt/F");
    tree_->Branch("dymva",         &dymva_        ,   "dymva/F");

    tree_->Branch("lep3", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lepPtr3_);
    tree_->Branch("lq3",           &lq3_,          "lq3/I");
    tree_->Branch("lid3",          &lid3_,         "lid3/I");
    tree_->Branch("lmva3",         &lmva3_,        "lmva3/F");
    tree_->Branch("jet3", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr3_);
    tree_->Branch("jet3Btag",      &jet3Btag_,      "jet3Btag/F");
    tree_->Branch("jet3ProbBtag",  &jet3ProbBtag_,  "jet3ProbBtag/F");
    tree_->Branch("jet4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jetPtr4_);
    tree_->Branch("jet4Btag",      &jet4Btag_,      "jet4Btag/F");
    tree_->Branch("jet4ProbBtag",  &jet4ProbBtag_,  "jet4ProbBtag/F");
    tree_->Branch("lep3DetEta"	 , &lep3DetEta_,    "lep3DetEta/F");
    tree_->Branch("lep3McId",      &lep3McId_,      "lep3McId/I");
    tree_->Branch("lep3MotherMcId",&lep3MotherMcId_,"lep3MotherMcId/I");
    tree_->Branch("jet3McId",      &jet3McId_,      "jet3McId/I");
    tree_->Branch("jet4McId",      &jet4McId_,      "jet4McId/I");
    tree_->Branch("dPhiLep3Jet1",  &dPhiLep3Jet1_,  "dPhiLep3Jet1/F");
    tree_->Branch("dRLep3Jet1",    &dRLep3Jet1_,    "dRLep3Jet1/F");
    tree_->Branch("dPhiLep3MET",   &dPhiLep3MET_,   "dPhiLep3MET/F");
    tree_->Branch("mt3",           &mt3_,           "mt3/F");
    tree_->Branch("jetLowBtag",    &jetLowBtag_,    "jetLowBtag/F");
    tree_->Branch("nSoftMuons",    &nSoftMuons_,    "nSoftMuons/i");
    tree_->Branch("Q",             &Q_	,     "Q/F");
    tree_->Branch("id1",           &id1_	,     "id1/F");
    tree_->Branch("x1",            &x1_	,     "x1/F");
    tree_->Branch("pdf1",          &pdf1_	,     "pdf1/F");
    tree_->Branch("id2",           &id2_	,     "id2/F");
    tree_->Branch("x2",            &x2_	,     "x2/F");
    tree_->Branch("pdf2",          &pdf2_	,     "pdf2/F");
    tree_->Branch("npu",           &npu_,           "npu/F");
    tree_->Branch("npuPlusOne",    &npuPlusOne_,    "npuPlusOne/F");
    tree_->Branch("npuMinusOne",   &npuMinusOne_,   "npuMinusOne/F");

    tree_->Branch("auxVar0",	   &auxVar0_	  ,   "auxVar0/F");
  }

  // initialze a SmurfTree
  void InitTree(int type = -1){
    assert(type==type); // just to suppress warnings
    assert(tree_);
    // don't forget to set pointers to zero before you set address
    // or you will fully appreciate that "ROOT sucks" :)
    InitVariables();
    //Set branch address
    Int_t currentState = gErrorIgnoreLevel;
    // gErrorIgnoreLevel = kError;
    gErrorIgnoreLevel = kBreak;
    tree_->SetBranchAddress("event",         &event_);
    tree_->SetBranchAddress("run",           &run_);
    tree_->SetBranchAddress("lumi",          &lumi_);
    tree_->SetBranchAddress("nvtx",          &nvtx_);
    tree_->SetBranchAddress("cuts",          &cuts_);
    tree_->SetBranchAddress("scale1fb",      &scale1fb_);
    tree_->SetBranchAddress("CHSMet",	     &CHSMet_);
    tree_->SetBranchAddress("NHSMet",	     &NHSMet_);
    tree_->SetBranchAddress("CHSMetPhi",     &CHSMetPhi_);
    tree_->SetBranchAddress("NHSMetPhi",     &NHSMetPhi_);
    tree_->SetBranchAddress("met",           &met_);
    tree_->SetBranchAddress("metPhi",        &metPhi_);
    tree_->SetBranchAddress("metMVA",        &metMVA_);
    tree_->SetBranchAddress("metMVAPhi",     &metMVAPhi_);
    tree_->SetBranchAddress("pmetMVA",       &pmetMVA_);
    tree_->SetBranchAddress("sumet",         &sumet_);
    tree_->SetBranchAddress("metSig",        &metSig_);
    tree_->SetBranchAddress("genmet",        &genmet_);
    tree_->SetBranchAddress("genmetPhi",     &genmetPhi_);
    tree_->SetBranchAddress("type",          &type_);
    tree_->SetBranchAddress("dstype",        &dstype_);
    tree_->SetBranchAddress("lep1",          &lepPtr1_);
    tree_->SetBranchAddress("lq1",           &lq1_);
    tree_->SetBranchAddress("lid1",          &lid1_);
    tree_->SetBranchAddress("lmva1",         &lmva1_);
    tree_->SetBranchAddress("lep2",          &lepPtr2_);
    tree_->SetBranchAddress("lq2",           &lq2_);
    tree_->SetBranchAddress("lid2",          &lid2_);
    tree_->SetBranchAddress("lmva2",         &lmva2_);
    tree_->SetBranchAddress("jet1",          &jetPtr1_);
    tree_->SetBranchAddress("jet1Btag",      &jet1Btag_);
    tree_->SetBranchAddress("jet1ProbBtag",  &jet1ProbBtag_);
    tree_->SetBranchAddress("jet1Dz",        &jet1Dz_);
    tree_->SetBranchAddress("jet2",          &jetPtr2_);
    tree_->SetBranchAddress("jet2Btag",      &jet2Btag_);
    tree_->SetBranchAddress("jet2ProbBtag",  &jet2ProbBtag_);
    tree_->SetBranchAddress("njets",         &njets_);
    tree_->SetBranchAddress("dilep",         &dilepPtr_);
    tree_->SetBranchAddress("quadlep",       &quadlepPtr_);
    tree_->SetBranchAddress("trackMet",      &trackMet_);
    tree_->SetBranchAddress("trackMetPhi",   &trackMetPhi_);
    tree_->SetBranchAddress("pmet",          &pmet_);
    tree_->SetBranchAddress("pTrackMet",     &pTrackMet_);
    tree_->SetBranchAddress("mt",            &mt_);
    tree_->SetBranchAddress("mt1",           &mt1_);
    tree_->SetBranchAddress("mt2",           &mt2_);
    tree_->SetBranchAddress("dPhi",          &dPhi_);
    tree_->SetBranchAddress("dR",            &dR_);
    tree_->SetBranchAddress("dPhiLep1Jet1",  &dPhiLep1Jet1_);
    tree_->SetBranchAddress("dRLep1Jet1",    &dRLep1Jet1_);
    tree_->SetBranchAddress("dPhiLep2Jet1",  &dPhiLep2Jet1_);
    tree_->SetBranchAddress("dRLep2Jet1",    &dRLep2Jet1_);
    tree_->SetBranchAddress("dPhiDiLepMET",  &dPhiDiLepMET_);
    tree_->SetBranchAddress("dPhiLep1MET",   &dPhiLep1MET_);
    tree_->SetBranchAddress("dPhiLep2MET",   &dPhiLep2MET_);
    tree_->SetBranchAddress("dPhiDiLepJet1", &dPhiDiLepJet1_);
    tree_->SetBranchAddress("lep1DetEta",    &lep1DetEta_);
    tree_->SetBranchAddress("lep2DetEta",    &lep2DetEta_);
    tree_->SetBranchAddress("lep1McId",      &lep1McId_);
    tree_->SetBranchAddress("lep2McId",      &lep2McId_);
    tree_->SetBranchAddress("lep1MotherMcId",&lep1MotherMcId_);
    tree_->SetBranchAddress("lep2MotherMcId",&lep2MotherMcId_);
    tree_->SetBranchAddress("jet1McId",      &jet1McId_);
    tree_->SetBranchAddress("jet2McId",      &jet2McId_);
    tree_->SetBranchAddress("processId",     &processId_);
    tree_->SetBranchAddress("higgsPt",       &higgsPt_);
    tree_->SetBranchAddress("sfWeightFR",    &sfWeightFR_);
    tree_->SetBranchAddress("sfWeightPU",    &sfWeightPU_);
    tree_->SetBranchAddress("sfWeightTrig",  &sfWeightTrig_);
    tree_->SetBranchAddress("sfWeightEff",   &sfWeightEff_);
    tree_->SetBranchAddress("sfWeightHPt",   &sfWeightHPt_);
    tree_->SetBranchAddress("dymva",         &dymva_);

    tree_->SetBranchAddress("lep3",	     &lepPtr3_);
    tree_->SetBranchAddress("lq3",	     &lq3_);
    tree_->SetBranchAddress("lid3",	     &lid3_);
    tree_->SetBranchAddress("lmva3",	     &lmva3_);
    tree_->SetBranchAddress("jet3",	     &jetPtr3_);
    tree_->SetBranchAddress("jet3Btag",      &jet3Btag_);
    tree_->SetBranchAddress("jet3ProbBtag",  &jet3ProbBtag_);
    tree_->SetBranchAddress("jet4",	     &jetPtr4_);
    tree_->SetBranchAddress("jet4Btag",      &jet4Btag_);
    tree_->SetBranchAddress("jet4ProbBtag",  &jet4ProbBtag_);
    tree_->SetBranchAddress("jet3McId",      &jet3McId_);
    tree_->SetBranchAddress("jet4McId",      &jet4McId_);
    tree_->SetBranchAddress("lep3DetEta",    &lep3DetEta_);
    tree_->SetBranchAddress("lep3McId",      &lep3McId_);
    tree_->SetBranchAddress("lep3MotherMcId",&lep3MotherMcId_);
    tree_->SetBranchAddress("dPhiLep3Jet1",  &dPhiLep3Jet1_);
    tree_->SetBranchAddress("dRLep3Jet1",    &dRLep3Jet1_);
    tree_->SetBranchAddress("dPhiLep3MET",   &dPhiLep3MET_);
    tree_->SetBranchAddress("mt3",	     &mt3_);
    tree_->SetBranchAddress("jetLowBtag",    &jetLowBtag_);
    tree_->SetBranchAddress("nSoftMuons",    &nSoftMuons_);
    tree_->SetBranchAddress("Q",	     &Q_);
    tree_->SetBranchAddress("id1",	     &id1_);
    tree_->SetBranchAddress("x1",	     &x1_);
    tree_->SetBranchAddress("pdf1",	     &pdf1_);
    tree_->SetBranchAddress("id2",	     &id2_);
    tree_->SetBranchAddress("x2",	     &x2_);
    tree_->SetBranchAddress("pdf2",	     &pdf2_);
    tree_->SetBranchAddress("npu",	     &npu_);
    tree_->SetBranchAddress("npuPlusOne",    &npuPlusOne_);
    tree_->SetBranchAddress("npuMinusOne",   &npuMinusOne_);

    tree_->SetBranchAddress("auxVar0",       &auxVar0_);

    gErrorIgnoreLevel = currentState;
  }

  /// get a built in type variable by name
  double Get(std::string value);
  /// compare two SmurfTrees for a given event on a given level of precision; 
  /// returns the variables that failed the comparison 
  std::vector<std::string> Compare(SmurfTree* value, double prec=0.005);
  /// transform DateType to string
  static std::string name(DataType type){
    switch (type){
    case data: return "data";
    case qqww: return "qqww";
    case ggww: return "ggww";
    case hww110: return "hww110";
    case hww115: return "hww115";
    case hww120: return "hww120";
    case hww130: return "hww130";
    case hww140: return "hww140";
    case hww150: return "hww150";
    case hww160: return "hww160";
    case hww170: return "hww170";
    case hww180: return "hww180";
    case hww190: return "hww190";
    case hww200: return "hww200";
    case hww210: return "hww210";
    case hww220: return "hww220";
    case hww230: return "hww230";
    case hww250: return "hww250";
    case hww300: return "hww300";
    case hww350: return "hww350";
    case hww400: return "hww400";
    case hww450: return "hww450";
    case hww500: return "hww500";
    case hww550: return "hww550";
    case hww600: return "hww600";
    case vbfhww110: return "vbfhww110";
    case vbfhww115: return "vbfhww115";
    case vbfhww120: return "vbfhww120";
    case vbfhww130: return "vbfhww130";
    case vbfhww140: return "vbfhww140";
    case vbfhww150: return "vbfhww150";
    case vbfhww160: return "vbfhww160";
    case vbfhww170: return "vbfhww170";
    case vbfhww180: return "vbfhww180";
    case vbfhww190: return "vbfhww190";
    case vbfhww200: return "vbfhww200";
    case vbfhww210: return "vbfhww210";
    case vbfhww220: return "vbfhww220";
    case vbfhww230: return "vbfhww230";
    case vbfhww250: return "vbfhww250";
    case vbfhww300: return "vbfhww300";
    case vbfhww350: return "vbfhww350";
    case vbfhww400: return "vbfhww400";
    case vbfhww450: return "vbfhww450";
    case vbfhww500: return "vbfhww500";
    case vbfhww550: return "vbfhww550";
    case vbfhww600: return "vbfhww600";
    case ttbar:  return "ttbar";
    case tw:     return "tw";
    case dyee:   return "dyee";
    case dymm:   return "dymm";
    case dytt:   return "dytt";
    case wjets:  return "wjets";
    case wz:     return "wz";
    case zz:     return "zz";
    case wgamma: return "wgamma";
    case qcd:    return "qcd";
    case other:  return "other";
    default:     return "uknown";
    }
  };

  private:

  LorentzVector* lepPtr1_;
  LorentzVector* lepPtr2_;
  LorentzVector* jetPtr1_;
  LorentzVector* jetPtr2_;
  LorentzVector* dilepPtr_;
  LorentzVector* quadlepPtr_;
  LorentzVector* lepPtr3_;
  LorentzVector* jetPtr3_;
  LorentzVector* jetPtr4_;
}; 

inline void 
SmurfTree::InitVariables(){
  // create list of available variables
  if(variables_.empty()){
    //make sure that this is only done once
    variables_.push_back(std::string("event"         ));
    variables_.push_back(std::string("run"           ));
    variables_.push_back(std::string("lumi"          ));
    variables_.push_back(std::string("nvtx"          ));
    variables_.push_back(std::string("scale1fb"      ));
    variables_.push_back(std::string("CHSMet"	     ));
    variables_.push_back(std::string("NHSMet"	     ));
    variables_.push_back(std::string("CHSMetPhi"     ));
    variables_.push_back(std::string("NHSMetPhi"     ));
    variables_.push_back(std::string("met"           ));
    variables_.push_back(std::string("metPhi"        ));
    variables_.push_back(std::string("sumet"         ));
    variables_.push_back(std::string("metSig"        ));
    variables_.push_back(std::string("genmet"        ));
    variables_.push_back(std::string("genmetPhi"     ));
    variables_.push_back(std::string("type"          ));
    variables_.push_back(std::string("dstype"        ));
    variables_.push_back(std::string("lq1"           ));
    variables_.push_back(std::string("lid1"          ));
    variables_.push_back(std::string("lmva1"         ));
    variables_.push_back(std::string("lq2"           ));
    variables_.push_back(std::string("lid2"          ));
    variables_.push_back(std::string("lmva2"         ));
    variables_.push_back(std::string("jet1Btag"      ));
    variables_.push_back(std::string("jet1ProbBtag"  ));
    variables_.push_back(std::string("jet1Dz"        ));
    variables_.push_back(std::string("jet2Btag"      ));
    variables_.push_back(std::string("jet2ProbBtag"  ));
    variables_.push_back(std::string("njets"         ));
    variables_.push_back(std::string("trackMet"      ));
    variables_.push_back(std::string("trackMetPhi"   ));
    variables_.push_back(std::string("pmet"          ));
    variables_.push_back(std::string("pTrackMet"     ));
    variables_.push_back(std::string("mt"            ));
    variables_.push_back(std::string("mt1"           ));
    variables_.push_back(std::string("mt2"           ));
    variables_.push_back(std::string("dPhi"          ));
    variables_.push_back(std::string("dR"            ));
    variables_.push_back(std::string("dPhiLep1Jet1"  ));
    variables_.push_back(std::string("dRLep1Jet1"    ));
    variables_.push_back(std::string("dPhiLep2Jet1"  ));
    variables_.push_back(std::string("dRLep2Jet1"    ));
    variables_.push_back(std::string("dPhiDiLepMET"  ));
    variables_.push_back(std::string("dPhiLep1MET"   ));
    variables_.push_back(std::string("dPhiLep2MET"   ));
    variables_.push_back(std::string("dPhiDiLepJet1" ));
    variables_.push_back(std::string("lep1McId"      ));
    variables_.push_back(std::string("lep2McId"      ));
    variables_.push_back(std::string("lep1MotherMcId"));
    variables_.push_back(std::string("lep2MotherMcId"));
    variables_.push_back(std::string("jet1McId"      ));
    variables_.push_back(std::string("jet2McId"      ));
  }
  // inizialize variables
  event_         = 0;
  run_           = 0;
  lumi_          = 0;
  nvtx_          = 0;
  cuts_          = 0;
  scale1fb_      = 0;
  CHSMet_        = -999.;
  NHSMet_        = -999.;
  CHSMetPhi_     = -999.;
  NHSMetPhi_     = -999.;
  met_           = -999.;
  metPhi_        = -999.;
  metMVA_        = -999.;
  metMVAPhi_     = -999.;
  pmetMVA_       = -999.;
  sumet_         = -999.;
  metSig_        = -999.;
  genmet_        = 0;
  genmetPhi_     = 0;
  type_          = mm;
  dstype_        = data;
  lq1_           = 0;
  lid1_          = 0;
  lmva1_         = 0.0;
  lq2_           = 0;
  lid2_          = 0;
  lmva2_         = 0.0;
  jet1Btag_      = -999.;
  jet1ProbBtag_  = -999.;
  jet1Dz_        = 0.;
  jet2Btag_      = -999.;
  jet2ProbBtag_  = -999.;
  njets_         = 0;
  trackMet_      = -999.;
  trackMetPhi_   = -999.;
  pmet_          = -999.;
  pTrackMet_     = -999.;
  mt_            = -999.;
  mt1_           = -999.;
  mt2_           = -999.;
  dPhi_          = -999.;
  dR_            = -999.;
  dPhiLep1Jet1_  = -999.;
  dRLep1Jet1_    = -999.;
  dPhiLep2Jet1_  = -999.;
  dRLep2Jet1_    = -999.;
  dPhiLep1MET_   = -999.;
  dPhiLep2MET_   = -999.;
  dPhiDiLepMET_  = -999.;
  dPhiDiLepJet1_ = -999.;
  lep1DetEta_    = -999.;
  lep2DetEta_    = -999.;
  lep1McId_   	 = 0;
  lep2McId_   	 = 0;
  lep1MotherMcId_= 0;
  lep2MotherMcId_= 0;
  jet1McId_   	 = 0;
  jet2McId_   	 = 0;
  lep1_       	 = LorentzVector();
  lep2_       	 = LorentzVector();
  jet1_       	 = LorentzVector();
  jet2_       	 = LorentzVector();
  dilep_      	 = LorentzVector();
  quadlep_     	 = LorentzVector();

  lep3_       	 = LorentzVector();
  lq3_  	 = 0;
  lid3_ 	 = 0;
  lmva3_ 	 = 0.0;
  jet3_       	 = LorentzVector();
  jet3Btag_	 = -999.;
  jet3ProbBtag_	 = -999.;
  jet4_       	 = LorentzVector();
  jet4Btag_	 = -999.;
  jet4ProbBtag_	 = -999.;
  jet3McId_	 = 0;
  jet4McId_	 = 0;
  lep3DetEta_    = -999.;
  lep3McId_	 = 0;
  lep3MotherMcId_= 0;
  dPhiLep3Jet1_  = -999.;
  dRLep3Jet1_	 = -999.;
  dPhiLep3MET_   = -999.;
  mt3_  	 = -999.;
  jetLowBtag_	 = -999.;
  nSoftMuons_	 = 0;
  Q_		 = -999.;
  id1_  	 = -999.;
  x1_		 = -999.;
  pdf1_ 	 = -999.;  
  id2_  	 = -999.;  
  x2_		 = -999.;
  pdf2_ 	 = -999.;  
  processId_	 = 0;
  higgsPt_	 = -999;
  sfWeightFR_    = -999;
  sfWeightPU_    = -999;
  sfWeightTrig_  = -999;
  sfWeightEff_   = -999;
  sfWeightHPt_   = -999;
  dymva_         = -999;
  npu_           = -999.;
  npuPlusOne_    = -999.;
  npuMinusOne_   = -999.;

  auxVar0_	 = -999.;
}

inline double
SmurfTree::Get(std::string value)
{
  if(value=="event"         ) { return this->event_;	     }
  if(value=="run"           ) { return this->run_;	     }
  if(value=="lumi"          ) { return this->lumi_;	     }
  if(value=="nvtx"          ) { return this->nvtx_;	     }
  if(value=="cuts"          ) { return this->cuts_;	     }
  if(value=="scale1fb"      ) { return this->scale1fb_;      }
  if(value=="CHSMet"	    ) { return this->CHSMet_; 	     }
  if(value=="NHSMet"	    ) { return this->NHSMet_; 	     }
  if(value=="CHSMetPhi"	    ) { return this->CHSMetPhi_;     }
  if(value=="NHSMetPhi"	    ) { return this->NHSMetPhi_;     }
  if(value=="met"           ) { return this->met_;	     }
  if(value=="metPhi"        ) { return this->metPhi_;	     }
  if(value=="metMVA"        ) { return this->metMVA_;        }
  if(value=="metMVAPhi"     ) { return this->metMVAPhi_;     }
  if(value=="pmetMVA"       ) { return this->pmetMVA_;       }
  if(value=="sumet"         ) { return this->sumet_;	     }
  if(value=="metSig"        ) { return this->metSig_;	     }
  if(value=="genmet"        ) { return this->genmet_;	     }
  if(value=="genmetPhi"     ) { return this->genmetPhi_;     }
  if(value=="type"          ) { return this->type_;	     }
  if(value=="dstype"        ) { return this->dstype_;	     }
  if(value=="lq1"           ) { return this->lq1_;	     }
  if(value=="lid1"          ) { return this->lid1_;	     }
  if(value=="lmva1"         ) { return this->lmva1_;	     }
  if(value=="lq2"           ) { return this->lq2_;	     }
  if(value=="lid2"          ) { return this->lid2_;	     }
  if(value=="lmva2"         ) { return this->lmva2_;	     }
  if(value=="jet1Btag"      ) { return this->jet1Btag_;      }
  if(value=="jet1ProbBtag"  ) { return this->jet1ProbBtag_;  }
  if(value=="jet1Dz"        ) { return this->jet1Dz_;        }
  if(value=="jet2Btag"      ) { return this->jet2Btag_;      }
  if(value=="jet2ProbBtag"  ) { return this->jet2ProbBtag_;  }
  if(value=="njets"         ) { return this->njets_;	     }
  if(value=="trackMet"      ) { return this->trackMet_;      }
  if(value=="trackMetPhi"   ) { return this->trackMetPhi_;   }
  if(value=="pmet"          ) { return this->pmet_;	     }
  if(value=="pTrackMet"     ) { return this->pTrackMet_;     }
  if(value=="mt"            ) { return this->mt_;	     }
  if(value=="mt1"           ) { return this->mt1_;	     }
  if(value=="mt2"           ) { return this->mt2_;	     }
  if(value=="dPhi"          ) { return this->dPhi_;	     }
  if(value=="dR"            ) { return this->dR_;	     }
  if(value=="dPhiLep1Jet1"  ) { return this->dPhiLep1Jet1_;  }
  if(value=="dRLep1Jet1"    ) { return this->dRLep1Jet1_;    }
  if(value=="dPhiLep2Jet1"  ) { return this->dPhiLep2Jet1_;  }
  if(value=="dRLep2Jet1"    ) { return this->dRLep2Jet1_;    }
  if(value=="dPhiDiLepMET"  ) { return this->dPhiDiLepMET_;  }
  if(value=="dPhiLep1MET"   ) { return this->dPhiLep1MET_;   }
  if(value=="dPhiLep2MET"   ) { return this->dPhiLep2MET_;   }
  if(value=="dPhiDiLepJet1" ) { return this->dPhiDiLepJet1_; }
  if(value=="lep1DetEta"    ) { return this->lep1DetEta_; }
  if(value=="lep2DetEta"    ) { return this->lep2DetEta_; }
  if(value=="lep1McId"      ) { return this->lep1McId_;      }
  if(value=="lep2McId"      ) { return this->lep2McId_;      }
  if(value=="lep1MotherMcId") { return this->lep1MotherMcId_;}
  if(value=="lep2MotherMcId") { return this->lep2MotherMcId_;}
  if(value=="jet1McId"      ) { return this->jet1McId_;      }
  if(value=="jet2McId"      ) { return this->jet2McId_;      }

  if(value=="lq3"	    ) { return this->lq3_;	     } 
  if(value=="lid3"	    ) { return this->lid3_;	     } 
  if(value=="lmva3"	    ) { return this->lmva3_;	     } 
  if(value=="jet3Btag"	    ) { return this->jet3Btag_;      } 
  if(value=="jet3ProbBtag"  ) { return this->jet3ProbBtag_;  } 
  if(value=="jet4Btag"	    ) { return this->jet4Btag_;      } 
  if(value=="jet4ProbBtag"  ) { return this->jet4ProbBtag_;  } 
  if(value=="jet3McId"	    ) { return this->jet3McId_;      } 
  if(value=="jet4McId"	    ) { return this->jet4McId_;      } 
  if(value=="lep3DetEta"    ) { return this->lep3DetEta_; }
  if(value=="lep3McId"	    ) { return this->lep3McId_;      } 
  if(value=="lep3MotherMcId") { return this->lep3MotherMcId_;} 
  if(value=="dPhiLep3Jet1"  ) { return this->dPhiLep3Jet1_;  } 
  if(value=="dRLep3Jet1"    ) { return this->dRLep3Jet1_;    } 
  if(value=="dPhiLep3MET"   ) { return this->dPhiLep3MET_;   } 
  if(value=="mt3"	    ) { return this->mt3_;	     } 
  if(value=="jetLowBtag"    ) { return this->jetLowBtag_;    } 
  if(value=="nSoftMuons"    ) { return this->nSoftMuons_ ;   } 
  if(value=="Q"             ) { return this->Q_;	     } 
  if(value=="id1"	    ) { return this->id1_;	     } 
  if(value=="x1"	    ) { return this->x1_;	     } 
  if(value=="pdf1"	    ) { return this->pdf1_;	     } 
  if(value=="id2"           ) { return this->id2_;	     } 
  if(value=="x2"	    ) { return this->x2_;	     } 
  if(value=="pdf2"	    ) { return this->pdf2_;	     } 
  if(value=="processId"	    ) { return this->processId_;     } 
  if(value=="higgsPt"	    ) { return this->higgsPt_;       } 
  if(value=="sfWeightFR"    ) { return this->sfWeightFR_;    } 
  if(value=="sfWeightPU"    ) { return this->sfWeightPU_;    } 
  if(value=="sfWeightTrig"  ) { return this->sfWeightTrig_;  } 
  if(value=="sfWeightEff"   ) { return this->sfWeightEff_;   } 
  if(value=="sfWeightHPt"   ) { return this->sfWeightHPt_;   } 
  if(value=="dymva"         ) { return this->dymva_;   } 
  if(value=="npu"	    ) { return this->npu_;           } 
  if(value=="npuPlusOne"    ) { return this->npuPlusOne_;    } 
  if(value=="npuMinusOne"   ) { return this->npuMinusOne_;   } 

  if(value=="auxVar0"	    ) { return this->auxVar0_;       }

  return -9999.; 
}

inline std::vector<std::string> 
SmurfTree::Compare(SmurfTree* value, double prec){
  std::vector<std::string> fails;
  // this should alway fit with ultimate precision
  if( this->event_ != value->event_ ){ fails.push_back( "event" ); }
  if( this->run_   != value->run_   ){ fails.push_back( "run"   ); }
  if( this->lumi_  != value->lumi_  ){ fails.push_back( "lumi"  ); }

  // check within  (relative) precision
  for(std::vector<std::string>::const_iterator var=variables_.begin(); var!=variables_.end(); ++var){
    if( fabs(Get(*var)-value->Get(*var))/Get(*var)>prec ) fails.push_back(*var);
  }
  return fails;
}


#endif
