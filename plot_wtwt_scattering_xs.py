from ROOT import *
import sys
import time

g_higgs = TGraph()
g_no_higgs = TGraph()

g_no_higgs.SetPoint(0,0.2,1804)
g_no_higgs.SetPoint(1,0.3,1809)
g_no_higgs.SetPoint(2,0.4,1819)
g_no_higgs.SetPoint(3,0.5,1815.8)
g_no_higgs.SetPoint(4,0.6,1819)
g_no_higgs.SetPoint(5,0.7,1810.73)
g_no_higgs.SetPoint(6,0.8,1810)
g_no_higgs.SetPoint(7,0.9,1811)
g_no_higgs.SetPoint(8,1.0,1815)
g_no_higgs.SetPoint(9,1.1,1811)

g_higgs.SetPoint(0,0.2,1780)
g_higgs.SetPoint(1,0.3,1799)
g_higgs.SetPoint(2,0.4,1817)
g_higgs.SetPoint(3,0.5,1805)
g_higgs.SetPoint(4,0.6,1808)
g_higgs.SetPoint(5,0.7,1807)
g_higgs.SetPoint(6,0.8,1803)
g_higgs.SetPoint(7,0.9,1812.37057792)
g_higgs.SetPoint(8,1,1813.10038709)
g_higgs.SetPoint(9,1.1,1812)


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
g_no_higgs.SetMaximum(4000)

g_no_higgs.Draw()
g_higgs.Draw("same")
leg.Draw("same")

gPad.SetLeftMargin(20)

gPad.Update()
#gPad.ForceUpdate()

gPad.SaveAs("/eos/user/a/amlevin/www/tmp/wtwt_scattering_xs.png")
