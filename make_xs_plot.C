#include "TGraphErrors.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TStyle.h"
#include "TPad.h"

void make_xs_plot()
{

  float phantom_xs_wminus_wminus = 4.73636178000000037 / 10000;
  float phantom_xs_error_wminus_wminus = 5.08068057350587693 / 10000000;

  float vbfnlo_xs_wminus_wminus = 0.420095092271086 / 1000 ;
  float vbfnlo_xs_error_wminus_wminus = 2.603485234411453 /10000 / 1000;

  float madgraph_xs_wminus_wminus = 0.0009035;
  float madgraph_xs_error_wminus_wminus = 7.1 / 1000000;

  float phantom_xs_wplus_wplus = 7.40094028000000057 / 10000;
  float phantom_xs_error_wplus_wplus = 1.38061691370198703 / 1000000;

  float vbfnlo_xs_wplus_wplus = 1.44889618993773 / 1000;
  float vbfnlo_xs_error_wplus_wplus = 9.496301531743701 /10000 / 1000;

  float madgraph_xs_wplus_wplus = 0.002325;
  float madgraph_xs_error_wplus_wplus = 1.6 /100000;

  float x_phantom_xs_wminus_wminus[1] = {1.1};
  float y_phantom_xs_wminus_wminus[1] = {phantom_xs_wminus_wminus};
  float ex_phantom_xs_wminus_wminus[1] = {0};
  float ey_phantom_xs_wminus_wminus[1] = {phantom_xs_error_wminus_wminus};

  float x_vbfnlo_xs_wminus_wminus[1] = {1.2};
  float y_vbfnlo_xs_wminus_wminus[1] = {vbfnlo_xs_wminus_wminus};
  float ex_vbfnlo_xs_wminus_wminus[1] = {0};
  float ey_vbfnlo_xs_wminus_wminus[1] = {vbfnlo_xs_error_wminus_wminus};

  float x_madgraph_xs_wminus_wminus[1] = {1.3};
  float y_madgraph_xs_wminus_wminus[1] = {madgraph_xs_wminus_wminus};
  float ex_madgraph_xs_wminus_wminus[1] = {0};
  float ey_madgraph_xs_wminus_wminus[1] = {madgraph_xs_error_wminus_wminus};

  float x_phantom_xs_wplus_wplus[1] = {1.1};
  float y_phantom_xs_wplus_wplus[1] = {phantom_xs_wplus_wplus};
  float ex_phantom_xs_wplus_wplus[1] = {0};
  float ey_phantom_xs_wplus_wplus[1] = {phantom_xs_error_wplus_wplus};

  float x_vbfnlo_xs_wplus_wplus[1] = {1.2};
  float y_vbfnlo_xs_wplus_wplus[1] = {vbfnlo_xs_wplus_wplus};
  float ex_vbfnlo_xs_wplus_wplus[1] = {0};
  float ey_vbfnlo_xs_wplus_wplus[1] = {vbfnlo_xs_error_wplus_wplus};

  float x_madgraph_xs_wplus_wplus[1] = {1.3};
  float y_madgraph_xs_wplus_wplus[1] = {madgraph_xs_wplus_wplus};
  float ex_madgraph_xs_wplus_wplus[1] = {0};
  float ey_madgraph_xs_wplus_wplus[1] = {madgraph_xs_error_wplus_wplus};

  TGraphErrors * phantom_wminus_wminus = new TGraphErrors(1,x_phantom_xs_wminus_wminus,y_phantom_xs_wminus_wminus,ex_phantom_xs_wminus_wminus,ey_phantom_xs_wminus_wminus);
  TGraphErrors * vbfnlo_wminus_wminus = new TGraphErrors(1,x_vbfnlo_xs_wminus_wminus,y_vbfnlo_xs_wminus_wminus,ex_vbfnlo_xs_wminus_wminus,ey_vbfnlo_xs_wminus_wminus);
  TGraphErrors * madgraph_wminus_wminus= new TGraphErrors(1,x_madgraph_xs_wminus_wminus,y_madgraph_xs_wminus_wminus,ex_madgraph_xs_wminus_wminus,ey_madgraph_xs_wminus_wminus);

  TGraphErrors *phantom_wplus_wplus= new TGraphErrors(1,x_phantom_xs_wplus_wplus,y_phantom_xs_wplus_wplus,ex_phantom_xs_wplus_wplus,ey_phantom_xs_wplus_wplus);
  TGraphErrors * vbfnlo_wplus_wplus= new TGraphErrors(1,x_vbfnlo_xs_wplus_wplus,y_vbfnlo_xs_wplus_wplus,ex_vbfnlo_xs_wplus_wplus,ey_vbfnlo_xs_wplus_wplus);
  TGraphErrors *madgraph_wplus_wplus= new TGraphErrors(1,x_madgraph_xs_wplus_wplus,y_madgraph_xs_wplus_wplus,ex_madgraph_xs_wplus_wplus,ey_madgraph_xs_wplus_wplus);


  TCanvas * c1 = new TCanvas;
  c1->GetFrame()->SetLineColor(0);
  
  // draw a frame to define the range
  TH1F *hr = c1->DrawFrame(1,0.00035,1.4,1.1*10e-4);
  hr->SetYTitle("cross section (in pb)");
  hr->GetXaxis()->SetAxisColor(0);
  hr->GetXaxis()->SetLabelColor(0);
  c1->GetFrame()->SetLineColor(0);
  c1->SetLeftMargin(0.2);
  //hr->GetYaxis()->SetLabelOffset(0.1);
  hr->GetYaxis()->SetTitleOffset(2);

  phantom_wminus_wminus->SetMarkerStyle(22);
  vbfnlo_wminus_wminus->SetMarkerStyle(22);
  madgraph_wminus_wminus->SetMarkerStyle(22);

  phantom_wminus_wminus->SetMarkerSize(2);
  vbfnlo_wminus_wminus->SetMarkerSize(2);
  madgraph_wminus_wminus->SetMarkerSize(2);


  phantom_wminus_wminus->SetMarkerColor(kRed);
  vbfnlo_wminus_wminus->SetMarkerColor(kBlue);
  phantom_wminus_wminus->SetFillColor(kRed);
  vbfnlo_wminus_wminus->SetFillColor(kBlue);
  phantom_wminus_wminus->SetLineColor(kRed);
  vbfnlo_wminus_wminus->SetLineColor(kBlue);

  phantom_wminus_wminus->Draw("P");
  vbfnlo_wminus_wminus->Draw("P");
  madgraph_wminus_wminus->Draw("P");

  TLegend *legend = new TLegend(.75,.80,.95,.95);
  legend->AddEntry(phantom_wminus_wminus,"phantom");
  legend->AddEntry(vbfnlo_wminus_wminus,"vbfnlo");
  legend->AddEntry(madgraph_wminus_wminus,"madgraph");
  legend->SetFillColor(0);
  legend->Draw();

  c1->Update();
  
  TCanvas * c2 = new TCanvas;
  c2->GetFrame()->SetLineColor(0);
  
  // draw a frame to define the range
  TH1F *hr = c2->DrawFrame(1,0.0001,1.4,0.004);
  hr->SetYTitle("cross section (in pb)");
  hr->GetXaxis()->SetAxisColor(0);
  hr->GetXaxis()->SetLabelColor(0);
  c2->GetFrame()->SetLineColor(0);
  c2->SetLeftMargin(0.2);
  //hr->GetYaxis()->SetLabelOffset(0.1);
  hr->GetYaxis()->SetTitleOffset(2);

  phantom_wplus_wplus->SetMarkerStyle(22);
  vbfnlo_wplus_wplus->SetMarkerStyle(22);
  madgraph_wplus_wplus->SetMarkerStyle(22);

  phantom_wplus_wplus->SetMarkerSize(2);
  vbfnlo_wplus_wplus->SetMarkerSize(2);
  madgraph_wplus_wplus->SetMarkerSize(2);


  phantom_wplus_wplus->SetMarkerColor(kRed);
  vbfnlo_wplus_wplus->SetMarkerColor(kBlue);
  phantom_wplus_wplus->SetFillColor(kRed);
  vbfnlo_wplus_wplus->SetFillColor(kBlue);
  phantom_wplus_wplus->SetLineColor(kRed);
  vbfnlo_wplus_wplus->SetLineColor(kBlue);

  phantom_wplus_wplus->Draw("P");
  vbfnlo_wplus_wplus->Draw("P");
  madgraph_wplus_wplus->Draw("P");

  TLegend *legend = new TLegend(.75,.80,.95,.95);
  legend->AddEntry(phantom_wplus_wplus,"phantom");
  legend->AddEntry(vbfnlo_wplus_wplus,"vbfnlo");
  legend->AddEntry(madgraph_wplus_wplus,"madgraph");
  legend->SetFillColor(0);
  legend->Draw();

  c2->Update();


}
