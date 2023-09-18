import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--toysfile',dest='toysfile',default="higgsCombineTest.GoodnessOfFit.mH120.123456.root")
parser.add_argument('--toy',dest='toy',default=1)
parser.add_argument('--channel',dest='channel',default=1)
parser.add_argument('--outfile',dest='outfile',default="toy.png")

args = parser.parse_args()

f=ROOT.TFile.Open(args.toysfile)

d=f.GetDirectory("toys")

toy=d.Get("toy_"+str(args.toy))

cmsth1x=toy.get(1)[0]

frame=cmsth1x.frame()

toy.plotOn(frame,ROOT.RooFit.Cut("CMS_channel == CMS_channel::ch"+str(args.channel)))

c=ROOT.TCanvas()

frame.Draw()

c.SaveAs(args.outfile)

