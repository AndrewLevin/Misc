from ROOT import *
import sys
import time

g_higgs = TGraph()
g_no_higgs = TGraph()

g_no_higgs.SetPoint(0,0.2,1183.64)
g_no_higgs.SetPoint(1,0.3,1473)
g_no_higgs.SetPoint(2,0.4,1819)
g_no_higgs.SetPoint(3,0.5,2234)
g_no_higgs.SetPoint(4,0.6,2734)
g_no_higgs.SetPoint(5,0.7,3305)
g_no_higgs.SetPoint(6,0.8,3957)
g_no_higgs.SetPoint(7,0.9,4683)
g_no_higgs.SetPoint(8,1.0,5517)
g_no_higgs.SetPoint(9,1.1,6416)

g_higgs.SetPoint(0,0.2,794.7)
g_higgs.SetPoint(1,0.3,801.1)
g_higgs.SetPoint(2,0.4,802.6)
g_higgs.SetPoint(3,0.5,802.4)
g_higgs.SetPoint(4,0.6,809.5)
g_higgs.SetPoint(5,0.7,813)
g_higgs.SetPoint(6,0.8,816.4029555)
g_higgs.SetPoint(7,0.9,823)
g_higgs.SetPoint(8,1,819.7)
g_higgs.SetPoint(9,1.1,830.3)


g_higgs.SetMarkerColor(kGreen+1)
g_no_higgs.SetMarkerColor(kBlue)

g_higgs.SetLineColor(kGreen+1)
g_no_higgs.SetLineColor(kBlue)

g_higgs.SetLineWidth(3)
g_no_higgs.SetLineWidth(3)

g_no_higgs.GetYaxis().SetTitleOffset(1.2)
g_no_higgs.GetXaxis().SetTitle("$\sqrt{s}$ (TeV)")
g_no_higgs.GetYaxis().SetTitle("cross-section (pb)")

leg=TLegend(0.2,0.6,0.4,0.8)

leg.AddEntry(g_no_higgs,"without Higgs","l")
leg.AddEntry(g_higgs,"with Higgs","l")

g_no_higgs.SetMinimum(0)

g_no_higgs.Draw()
g_higgs.Draw("same")
leg.Draw("same")

gPad.SetLeftMargin(20)

gPad.Update()
#gPad.ForceUpdate()

gPad.SaveAs("/eos/user/a/amlevin/www/tmp/wlwl_scattering_xs.png")
