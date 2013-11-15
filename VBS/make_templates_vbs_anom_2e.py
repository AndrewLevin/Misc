#!/usr/bin/env python

from ROOT import *
from array import array
import sys
import time

def fillDownHist(downHist,upHist,unscaledHist):
    if downHist.GetNbinsX() != upHist.GetNbinsX():
        sys.exit(0)
    if downHist.GetNbinsX() != unscaledHist.GetNbinsX():
        sys.exit(0)

    print downHist.GetNbinsX()      

    for bin in range(downHist.GetNbinsX()+1):
       down_hist_bin_content = unscaledHist.GetBinContent(bin) - (upHist.GetBinContent(bin) - unscaledHist.GetBinContent(bin))
       if down_hist_bin_content >= 0:
           downHist.SetBinContent(bin, down_hist_bin_content)
       else:
           downHist.SetBinContent(bin, 0)
    

def countEvents(tree):
    nevents = 0 
    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        if tree.nevents != -1:
            if tree.nevents < 0:
                print "number of events should be greater than 0"
                sys.exit()
            nevents = nevents + tree.nevents
    return nevents        
                
def fillTemplate(filename,xs,hist):
    print filename
    f=TFile(filename,"r")
    vars_tree=f.Get("vars")
    nevents=countEvents(vars_tree)
    for entry in range(vars_tree.GetEntries()):
        vars_tree.GetEntry(entry)
        if vars_tree.weight == -1:
            continue
        weight = vars_tree.weight*xs*luminosity/nevents
        hist.Fill(vars_tree.mlljj,weight)                    

args=sys.argv[1:]
if not len(args)==4:
    print "usage: python make_templates_vbs_anom_2e.py signal1_label signal1_xs signal2_label signal2_xs"
    sys.exit(0)

signal1_label=args[0]
signal2_label=args[2]

#in 1/pb
luminosity = 300000.

#in pb
xs_background_0to300 = 249.97710
xs_background_300to700 = 35.23062
xs_background_700to1300 = 4.13743
xs_background_1300to2100 = 0.41702
xs_background_2100to100000 = 0.04770

xs_signal1 = float(args[1]) 
xs_signal2 = float(args[3])

unscaled_background_filename_0to300="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-0-300-v1510_14TEV_2e_unscaled.root"
unscaled_background_filename_1300to2100="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-1300-2100-v1510_14TEV_2e_unscaled.root"
unscaled_background_filename_2100to100000="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-2100-100000_14TEV_2e_unscaled.root"
unscaled_background_filename_300to700="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-300-700-v1510_14TEV_2e_unscaled.root"
unscaled_background_filename_700to1300="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-700-1300-v1510_14TEV_2e_unscaled.root"

unscaled_signal1_filename="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/signal/ww_to_ll_same_sign_" + signal1_label + "_2e_unscaled.root"
unscaled_signal2_filename="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/signal/ww_to_ll_same_sign_" + signal2_label + "_2e_unscaled.root"

electronscaled_background_filename_0to300="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-0-300-v1510_14TEV_2e_electron_scaled.root"
electronscaled_background_filename_1300to2100="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-1300-2100-v1510_14TEV_2e_electron_scaled.root"
electronscaled_background_filename_2100to100000="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-2100-100000_14TEV_2e_electron_scaled.root"
electronscaled_background_filename_300to700="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-300-700-v1510_14TEV_2e_electron_scaled.root"
electronscaled_background_filename_700to1300="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-700-1300-v1510_14TEV_2e_electron_scaled.root"

electronscaled_signal1_filename="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/signal/ww_to_ll_same_sign_" + signal1_label + "_2e_electron_scaled.root"
electronscaled_signal2_filename="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/signal/ww_to_ll_same_sign_" + signal2_label + "_2e_electron_scaled.root"

jetscaled_background_filename_0to300="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-0-300-v1510_14TEV_2e_jet_scaled.root"
jetscaled_background_filename_1300to2100="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-1300-2100-v1510_14TEV_2e_jet_scaled.root"
jetscaled_background_filename_2100to100000="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-2100-100000_14TEV_2e_jet_scaled.root"
jetscaled_background_filename_300to700="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-300-700-v1510_14TEV_2e_jet_scaled.root"
jetscaled_background_filename_700to1300="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/BB-4p-700-1300-v1510_14TEV_2e_jet_scaled.root"

jetscaled_signal1_filename="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/signal/ww_to_ll_same_sign_" + signal1_label + "_2e_jet_scaled.root"
jetscaled_signal2_filename="/afs/cern.ch/work/a/anlevin/data/select_vbs_output/signal/ww_to_ll_same_sign_"+signal2_label + "_2e_jet_scaled.root"

f = TFile( 'templates.root', 'recreate' )  

jet_scaled_signal1_mlljj_hist = TH1F( 'sig_jet_energy_scaleUp', 'sig_jet_energy_scaleUp', 35, 0., 7000 )
jet_scaled_signal2_mlljj_hist = TH1F( 'sig_ALT_jet_energy_scaleUp', 'sig_ALT_jet_energy_scaleUp', 35, 0., 7000 )
jet_scaled_background_mlljj_hist = TH1F( 'back1_jet_energy_scaleUp', 'back1_jet_energy_scaleUp', 35, 0.,7000. )

electron_scaled_signal1_mlljj_hist = TH1F( 'sig_ele_energy_scaleUp', 'sig_ele_energy_scaleUp', 35, 0., 7000 )
electron_scaled_signal2_mlljj_hist = TH1F( 'sig_ALT_ele_energy_scaleUp', 'sig_ALT_ele_energy_scaleUp', 35, 0., 7000 )
electron_scaled_background_mlljj_hist = TH1F( 'back1_ele_energy_scaleUp', 'back1_ele_energy_scaleUp', 35, 0.,7000. )       

jet_scaled_signal1_mlljj_down_hist = TH1F( 'sig_jet_energy_scaleDown', 'sig_jet_energy_scaleDown', 35, 0., 7000 )
jet_scaled_signal2_mlljj_down_hist = TH1F( 'sig_ALT_jet_energy_scaleDown', 'sig_ALT_jet_energy_scaleDown', 35, 0., 7000 )
jet_scaled_background_mlljj_down_hist = TH1F( 'back1_jet_energy_scaleDown', 'back1_jet_energy_scaleDown', 35, 0.,7000. )

electron_scaled_signal1_mlljj_down_hist = TH1F( 'sig_ele_energy_scaleDown', 'sig_ele_energy_scaleDown', 35, 0., 7000 )
electron_scaled_signal2_mlljj_down_hist = TH1F( 'sig_ALT_ele_energy_scaleDown', 'sig_ALT_ele_energy_scaleDown', 35, 0., 7000 )
electron_scaled_background_mlljj_down_hist = TH1F( 'back1_ele_energy_scaleDown', 'back1_ele_energy_scaleDown', 35, 0.,7000. )       

unscaled_signal1_mlljj_hist = TH1F( 'sig', 'sig', 35, 0., 7000 ) 
unscaled_signal2_mlljj_hist = TH1F( 'sig_ALT', 'sig_ALT', 35, 0., 7000 )
unscaled_background_mlljj_hist = TH1F( 'back1', 'back1', 35, 0.,7000. )

#fill the unscaled signal templates
fillTemplate(unscaled_signal1_filename,xs_signal1,unscaled_signal1_mlljj_hist)
fillTemplate(unscaled_signal2_filename,xs_signal2,unscaled_signal2_mlljj_hist)

#fill the unscaled background template
fillTemplate(unscaled_background_filename_0to300,xs_background_0to300,unscaled_background_mlljj_hist)
fillTemplate(unscaled_background_filename_1300to2100,xs_background_1300to2100,unscaled_background_mlljj_hist)
fillTemplate(unscaled_background_filename_2100to100000,xs_background_2100to100000,unscaled_background_mlljj_hist)
fillTemplate(unscaled_background_filename_300to700,xs_background_300to700,unscaled_background_mlljj_hist)
fillTemplate(unscaled_background_filename_700to1300,xs_background_700to1300,unscaled_background_mlljj_hist)

#fill the jet scaled signal templates
fillTemplate(jetscaled_signal1_filename,xs_signal1,jet_scaled_signal1_mlljj_hist)
fillTemplate(jetscaled_signal2_filename,xs_signal2,jet_scaled_signal2_mlljj_hist)

#fill the jet scaled background template
fillTemplate(jetscaled_background_filename_0to300,xs_background_0to300,jet_scaled_background_mlljj_hist)
fillTemplate(jetscaled_background_filename_1300to2100,xs_background_1300to2100,jet_scaled_background_mlljj_hist)
fillTemplate(jetscaled_background_filename_2100to100000,xs_background_2100to100000,jet_scaled_background_mlljj_hist)
fillTemplate(jetscaled_background_filename_300to700,xs_background_300to700,jet_scaled_background_mlljj_hist)
fillTemplate(jetscaled_background_filename_700to1300,xs_background_700to1300,jet_scaled_background_mlljj_hist)

#fill the lepton scaled signal templates
fillTemplate(electronscaled_signal1_filename,xs_signal1,electron_scaled_signal1_mlljj_hist)
fillTemplate(electronscaled_signal2_filename,xs_signal2,electron_scaled_signal2_mlljj_hist)

#fill the lepton scaled background template
fillTemplate(electronscaled_background_filename_0to300,xs_background_0to300,electron_scaled_background_mlljj_hist)
fillTemplate(electronscaled_background_filename_1300to2100,xs_background_1300to2100,electron_scaled_background_mlljj_hist)
fillTemplate(electronscaled_background_filename_2100to100000,xs_background_2100to100000,electron_scaled_background_mlljj_hist)
fillTemplate(electronscaled_background_filename_300to700,xs_background_300to700,electron_scaled_background_mlljj_hist)
fillTemplate(electronscaled_background_filename_700to1300,xs_background_700to1300,electron_scaled_background_mlljj_hist)

fillDownHist(jet_scaled_signal1_mlljj_down_hist,jet_scaled_signal1_mlljj_hist, unscaled_signal1_mlljj_hist)
fillDownHist(jet_scaled_signal2_mlljj_down_hist,jet_scaled_signal2_mlljj_hist, unscaled_signal2_mlljj_hist)
fillDownHist(jet_scaled_background_mlljj_down_hist,jet_scaled_background_mlljj_hist, unscaled_background_mlljj_hist)

fillDownHist(electron_scaled_signal1_mlljj_down_hist,electron_scaled_signal1_mlljj_hist, unscaled_signal1_mlljj_hist)
fillDownHist(electron_scaled_signal2_mlljj_down_hist,electron_scaled_signal2_mlljj_hist, unscaled_signal2_mlljj_hist)
fillDownHist(electron_scaled_background_mlljj_down_hist,electron_scaled_background_mlljj_hist, unscaled_background_mlljj_hist)

gStyle.SetOptStat(0)
    
c1=TCanvas()

leg1_mlljj = TLegend(.75,.80,.95,.95)
leg1_mlljj.SetFillColor(0)
leg1_mlljj.AddEntry(unscaled_signal1_mlljj_hist,"sm")
leg1_mlljj.AddEntry(unscaled_signal2_mlljj_hist,"anom")
leg1_mlljj.AddEntry(unscaled_background_mlljj_hist,"background")

unscaled_signal1_mlljj_hist.SetLineColor(kRed)
unscaled_signal2_mlljj_hist.SetLineColor(kGreen)
unscaled_background_mlljj_hist.SetLineColor(kBlue)

unscaled_background_mlljj_hist.Draw()
unscaled_signal1_mlljj_hist.Draw("SAME")
unscaled_signal2_mlljj_hist.Draw("SAME")
leg1_mlljj.Draw("SAME")

time.sleep(1)

c2=TCanvas()

leg2_mlljj = TLegend(.75,.80,.95,.95)
leg2_mlljj.SetFillColor(0)
leg2_mlljj.AddEntry(unscaled_signal1_mlljj_hist,"unscaled sm")
leg2_mlljj.AddEntry(jet_scaled_signal1_mlljj_hist,"jet scaled up")
leg2_mlljj.AddEntry(jet_scaled_signal1_mlljj_down_hist,"jet scaled down")

unscaled_signal1_mlljj_hist.SetLineColor(kRed)
jet_scaled_signal1_mlljj_hist.SetLineColor(kGreen)
jet_scaled_signal1_mlljj_down_hist.SetLineColor(kBlue)

jet_scaled_signal1_mlljj_hist.Draw()
unscaled_signal1_mlljj_hist.Draw("SAME")
jet_scaled_signal1_mlljj_down_hist.Draw("SAME")
leg2_mlljj.Draw("SAME")

time.sleep(1)

c3=TCanvas()

leg3_mlljj = TLegend(.75,.80,.95,.95)
leg3_mlljj.SetFillColor(0)
leg3_mlljj.AddEntry(unscaled_signal1_mlljj_hist,"unscaled sm")
leg3_mlljj.AddEntry(electron_scaled_signal1_mlljj_hist,"electron scaled up")
leg3_mlljj.AddEntry(electron_scaled_signal1_mlljj_down_hist,"electron scaled down")

unscaled_signal1_mlljj_hist.SetLineColor(kRed)
electron_scaled_signal1_mlljj_hist.SetLineColor(kGreen)
electron_scaled_signal1_mlljj_down_hist.SetLineColor(kBlue)

electron_scaled_signal1_mlljj_hist.Draw()
unscaled_signal1_mlljj_hist.Draw("SAME")
electron_scaled_signal1_mlljj_down_hist.Draw("SAME")
leg3_mlljj.Draw("SAME")

c4=TCanvas()

leg4_mlljj = TLegend(.75,.80,.95,.95)
leg4_mlljj.SetFillColor(0)
leg4_mlljj.AddEntry(unscaled_signal2_mlljj_hist,"unscaled anom")
leg4_mlljj.AddEntry(jet_scaled_signal2_mlljj_hist,"jet scaled up")
leg4_mlljj.AddEntry(jet_scaled_background_mlljj_down_hist,"jet scaled down")

unscaled_signal2_mlljj_hist.SetLineColor(kRed)
jet_scaled_signal2_mlljj_hist.SetLineColor(kGreen)
jet_scaled_signal2_mlljj_down_hist.SetLineColor(kBlue)

jet_scaled_signal2_mlljj_hist.Draw()
unscaled_signal2_mlljj_hist.Draw("SAME")
jet_scaled_signal2_mlljj_down_hist.Draw("SAME")
leg4_mlljj.Draw("SAME")

time.sleep(1)

c5=TCanvas()

leg5_mlljj = TLegend(.75,.80,.95,.95)
leg5_mlljj.SetFillColor(0)
leg5_mlljj.AddEntry(unscaled_signal2_mlljj_hist,"unscaled anom")
leg5_mlljj.AddEntry(electron_scaled_signal2_mlljj_hist,"electron scaled up")
leg5_mlljj.AddEntry(electron_scaled_background_mlljj_down_hist,"electron scaled down")

unscaled_signal2_mlljj_hist.SetLineColor(kRed)
electron_scaled_signal2_mlljj_hist.SetLineColor(kGreen)
electron_scaled_background_mlljj_down_hist.SetLineColor(kBlue)

electron_scaled_signal2_mlljj_hist.Draw()
unscaled_signal2_mlljj_hist.Draw("SAME")
electron_scaled_signal2_mlljj_down_hist.Draw("SAME")
leg5_mlljj.Draw("SAME")

c6=TCanvas()

leg6_mlljj = TLegend(.75,.80,.95,.95)
leg6_mlljj.SetFillColor(0)
leg6_mlljj.AddEntry(unscaled_background_mlljj_hist,"unscaled background")
leg6_mlljj.AddEntry(jet_scaled_background_mlljj_hist,"jet scaled up")
leg6_mlljj.AddEntry(jet_scaled_background_mlljj_down_hist,"jet scaled down")

unscaled_background_mlljj_hist.SetLineColor(kRed)
jet_scaled_background_mlljj_hist.SetLineColor(kGreen)
jet_scaled_background_mlljj_down_hist.SetLineColor(kBlue)

jet_scaled_background_mlljj_hist.Draw()
unscaled_background_mlljj_hist.Draw("SAME")
jet_scaled_background_mlljj_down_hist.Draw("SAME")
leg6_mlljj.Draw("SAME")

time.sleep(1)

c7=TCanvas()

leg7_mlljj = TLegend(.75,.80,.95,.95)
leg7_mlljj.SetFillColor(0)
leg7_mlljj.AddEntry(unscaled_background_mlljj_hist,"unscaled background")
leg7_mlljj.AddEntry(electron_scaled_background_mlljj_hist,"electron scaled up")
leg7_mlljj.AddEntry(electron_scaled_background_mlljj_down_hist,"electron scaled down")

unscaled_background_mlljj_hist.SetLineColor(kRed)
electron_scaled_background_mlljj_hist.SetLineColor(kGreen)
electron_scaled_background_mlljj_down_hist.SetLineColor(kBlue)

electron_scaled_background_mlljj_hist.Draw()
unscaled_background_mlljj_hist.Draw("SAME")
electron_scaled_background_mlljj_down_hist.Draw("SAME")
leg7_mlljj.Draw("SAME")

print "unscaled_signal2_mlljj_hist.Integral(0,unscaled_signal2_mlljj_hist.GetNbinsX()+1) = " + str(unscaled_signal2_mlljj_hist.Integral(0,unscaled_signal2_mlljj_hist.GetNbinsX()+1))
print "unscaled_signal1_mlljj_hist.Integral(0,unscaled_signal1_mlljj_hist.GetNbinsX()+1) = " + str(unscaled_signal1_mlljj_hist.Integral(0,unscaled_signal1_mlljj_hist.GetNbinsX()+1))
print "unscaled_background_mlljj_hist.Integral(0,unscaled_background_mlljj_hist.GetNbinsX()+1) = " + str(unscaled_background_mlljj_hist.Integral(0,unscaled_background_mlljj_hist.GetNbinsX()+1)) 

time.sleep(1)                                                                                                                                            

#raw_input("hit a key")

f.Write()
f.Close()
