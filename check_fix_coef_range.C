#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"

void check_fix_coef_range()
{

  RooRealVar x("x","x",-20,20) ;
  RooRealVar m("m","m",0,-10,10) ;
  RooRealVar s("s","s",1,-10,10) ;

  RooGaussian gauss("g","g",x,m,s) ;

  // Construct poly(x,p0)                                                                                                                                           
  RooRealVar p0("p0","p0",0.01,0.,1.) ;
  RooPolynomial poly("p","p",x,p0) ;

  // Construct model = f*gauss(x) + (1-f)*poly(x)                                                                                                                   
  RooRealVar f("f","f",0.5,0.,1.) ;
  RooAddPdf model("model","model",RooArgSet(gauss,poly),f);

RooPlot* frame1 = x.frame() ;
model.plotOn(frame1) ;
 model.plotOn(frame1,RooFit::Components("p"),RooFit::LineStyle(kDashed));

x.setRange(-5,5) ;

RooPlot* frame2 = x.frame() ;
model.plotOn(frame2) ;
 model.plotOn(frame2,RooFit::Components("p"),RooFit::LineStyle(kDashed));

   x.setRange("ref",-20,20);

 model.fixCoefRange("ref");

RooPlot* frame3 = x.frame() ;
model.plotOn(frame3) ;
 model.plotOn(frame3,RooFit::Components("p"),RooFit::LineStyle(kDashed));

 TCanvas * c1 = new TCanvas();

 frame1->Draw();

 c1->SaveAs("/afs/cern.ch/user/a/anlevin/www/tmp/c1.png");

 TCanvas * c2 = new TCanvas();

 frame2->Draw();

 c2->SaveAs("/afs/cern.ch/user/a/anlevin/www/tmp/c2.png");

 TCanvas * c3 = new TCanvas();

 frame3->Draw();

 c3->SaveAs("/afs/cern.ch/user/a/anlevin/www/tmp/c3.png");

}
