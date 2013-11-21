#!/usr/bin/env python

from ROOT import *
from array import array
import sys
import time

gStyle.SetFillStyle(0)
gStyle.SetLegendBorderSize(0); 
gROOT.ForceStyle()

f=TFile("output_distributions_00.root","r")
mjj_hist_1=f.Get("mjj_1")
mjj_hist_2=f.Get("mjj_2")

mjj_hist_1.SetTitle("")
mjj_hist_2.SetTitle("")

mjj_hist_1.Scale(1/mjj_hist_1.Integral())
mjj_hist_2.Scale(1/mjj_hist_2.Integral())
mjj_hist_1.SetLineColor(kRed)
mjj_hist_2.SetLineColor(kBlue)
mjj_hist_1.SetLineWidth(3)
mjj_hist_2.SetLineWidth(3)

mjj_hist_1.GetXaxis().SetTitle("mjj [GeV]")
mjj_hist_2.GetXaxis().SetTitle("mjj [GeV]")

mjj_hist_1.SetStats(0)
mjj_hist_2.SetStats(0)

leg=TLegend(.40,.45,.95,.95)

leg.AddEntry(mjj_hist_1,"FS0 = 0 TeV^-4","l")
leg.AddEntry(mjj_hist_2,"FS0 = 5 TeV^-4","l")
leg.SetFillColor(0)

mjj_hist_1.Draw()
mjj_hist_2.Draw("SAME")
leg.Draw("SAME")

time.sleep(1)                                                                                                                                            

raw_input("hit a key")
