#include <iostream>
#include <vector>
using namespace std;
#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include <limits>

#ifndef ROOT_Math_DistFuncMathCore
#include"Math/DistFuncMathCore.h"
#endif


void poisson(int sel = 2, int nevt = 1000000, 
	     double nexp = 1.0, double nsyst = 0., int nobs = 1.0){
 TH1D *hnevt1	= new TH1D("nevt1","nevt1",150,-0.5,149.5);
 TH1D *hnevt2	= new TH1D("nevt2","nevt2",150,-0.5,149.5);
 TH1D *hnevt3	= new TH1D("nevt3","nevt3",150,-0.5,149.5);
 double nevt1,nevt2,nevt3;
 gRandom->SetSeed(666);
 for(int i=0; i<nevt; i++){
   nevt1 = -0.1;
   while(nevt1 < 0){
     nevt1 = gRandom->Gaus(nexp,sqrt(nexp)) + 
             gRandom->Gaus(0.0,nsyst);
   }
   hnevt1->Fill(nevt1);

   nevt2 = -0.1;
   while(nevt2 < 0){
     nevt2 = gRandom->Poisson(nexp) + 
             gRandom->Gaus(0.0,nsyst);
   }
   hnevt2->Fill(nevt2);

   nevt3 = -0.1;
   while(nevt3 < 0){
     nevt3 = gRandom->Gaus(nexp,nsyst);
   }
   nevt3 = gRandom->Poisson(nevt3);
   hnevt3->Fill(nevt3);

   if(i%1000000 == 0) cout << i << endl;
 }
 hnevt1->Scale(1./hnevt1->GetSumOfWeights());
 hnevt2->Scale(1./hnevt2->GetSumOfWeights());
 hnevt3->Scale(1./hnevt3->GetSumOfWeights());
 hnevt1->SetLineColor(kRed);
 hnevt2->SetLineColor(kBlue);
 hnevt3->SetLineColor(kBlack);

 //if     (sel == 0){ // G & P & GP
 //  hnevt1->Draw("");
 //  hnevt2->Draw("same");
 //  hnevt3->Draw("same");
 //}
 //else if(sel == 1){ // P & GP
 //  hnevt2->Draw("");
 //  hnevt3->Draw("same");
 //}
 //else if(sel == 2){ // P
 //  hnevt2->Draw("");
 //}

 assert(hnevt2->GetXaxis()->GetBinCenter(nobs+1)==nobs);

 double probability = hnevt2->Integral(nobs+1,std::numeric_limits<Int_t>::max());
 std::cout << "probability = " << probability << std::endl;
 double significance = ::ROOT::Math::normal_quantile_c(probability,1);
 std::cout << "significance = " << significance << std::endl;

}
