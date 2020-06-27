import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

data_driven = True
data_driven_correction = True

import json
import sys
import style

import argparse

from math import hypot, pi, sqrt, acos, cos, sin, atan2

#make the Down shape from the Up shape
def makeDownShape(histUp,hist):

    hist_clone=hist.Clone()
    histUp_clone=histUp.Clone()

    hist_clone.Scale(2)
    histUp_clone.Scale(-1)

    hist_clone.Add(histUp_clone)

    return hist_clone

dict_lumi = {"2016" : 35.9, "2017" : 41.5, "2018" : 59.6 , "run2" : 137.1}

parser = argparse.ArgumentParser()

parser.add_argument('--make_unc_table',action='store_true',default=False)
parser.add_argument('--use_yings_electron_hlt_sfs',action='store_true',default=False)
parser.add_argument('--singleproc',action='store_true',default=False)
parser.add_argument('--nthreads',dest='nthreads',default=0) #the argument to EnableImplicitMT
parser.add_argument('--userdir',dest='userdir',default='/afs/cern.ch/user/a/amlevin/') #not used now
parser.add_argument('--workdir',dest='workdir',default='/afs/cern.ch/work/a/amlevin/public/')
parser.add_argument('--draw_non_fid',dest='draw_non_fid',action='store_true')
parser.add_argument('--lep',dest='lep',default='both')
parser.add_argument('--year',dest='year',default='all')
parser.add_argument('--zveto',dest='zveto',action='store_true',default=False)
parser.add_argument('--phoeta',dest='phoeta',default='both')
parser.add_argument('--make_datacard',dest='make_datacard',action='store_true',default=False)
parser.add_argument('--make_cut_and_count_datacard',dest='make_cut_and_count_datacard',action='store_true',default=False)
parser.add_argument('--closure_test',dest='closure_test',action='store_true',default=False)
parser.add_argument('--apply_2018_photon_phi_cut',dest='apply_2018_photon_phi_cut',action='store_true',default=False)
parser.add_argument('--no_wjets_for_2017_and_2018',dest='no_wjets_for_2017_and_2018',action='store_true',default=False)
parser.add_argument('--ewdim6',dest='ewdim6',action='store_true',default=False)
parser.add_argument('--use_wjets_for_fake_photon',dest='use_wjets_for_fake_photon',action='store_true',default=False)
parser.add_argument('--float_sig_fake_cont',dest='float_sig_fake_cont',action='store_true',default=False,help="in the datacard, float the contamination of the fake photon, fake lepton, and double fake backgrounds due to the signal")
parser.add_argument('--draw_ewdim6',dest='draw_ewdim6',action='store_true',default=False)
parser.add_argument('--make_all_plots',dest='make_all_plots',action='store_true',default=False)
parser.add_argument('-o',dest='outputdir',default="/eos/user/a/amlevin/www/tmp/")

args = parser.parse_args()

n_fake_photon_alt = 24
#n_fake_photon_alt = 0

if args.year == "2016":
    years = ["2016"]
    totallumi=dict_lumi["2016"]
elif args.year == "2017":
    years=["2017"]
    totallumi=dict_lumi["2017"]
elif args.year == "2018":
    years=["2018"]
    totallumi=dict_lumi["2018"]
elif args.year == "run2":
    years=["2016","2017","2018"]
    totallumi=dict_lumi["run2"]
else:
    assert(0)

den_pho_sel = 4

sieie_cut_2016_barrel = 0.01022
sieie_cut_2016_endcap = 0.03001
sieie_cut_2017_barrel = 0.01015
sieie_cut_2017_endcap = 0.0272
sieie_cut_2018_barrel = 0.01015
sieie_cut_2018_endcap = 0.0272

chiso_cut_2016_barrel = 1.416
chiso_cut_2016_endcap = 1.012
chiso_cut_2017_barrel = 1.141
chiso_cut_2017_endcap = 1.051
chiso_cut_2018_barrel = 1.141
chiso_cut_2018_endcap = 1.051

if args.lep == "muon":
    lepton_name = "muon"
elif args.lep == "electron":
    lepton_name = "electron"
elif args.lep == "both":
    lepton_name = "both"
else:
    assert(0)

if args.phoeta == "barrel":
    photon_eta_min = 0
    photon_eta_max = 1.5
elif args.phoeta == "endcap":
    photon_eta_min = 1.5
    photon_eta_max = 2.5
elif args.phoeta == "both":
    photon_eta_min = 0
    photon_eta_max = 2.5
else:
    assert(0)

photon_eta_cutstring = "((abs(photon_eta) < "+str(photon_eta_max)+") && (abs(photon_eta) > "+str(photon_eta_min)+"))"

def get_postfilter_selection_string(syst="nominal"):

    assert(syst == "nominal" or syst == "JESUp" or syst == "JERUp")

    if syst == "nominal":
        return "(puppimet > 40)"
    elif syst == "JESUp":
        return "(puppimetJESUp > 40)"
    elif syst == "JERUp":
        return "(puppimetJERUp > 40)"
    else:
        assert(0)

def get_filter_string(year,isdata=True,lep=None):

    if lep == None:
        lep = args.lep

    if not isdata:
        puppimet_cutstring = "(puppimet > 40 || puppimetJESUp > 40 || puppimetJERUp > 40)"
    else:    
        puppimet_cutstring = "(puppimet > 40)"

    if args.zveto:
        zveto_cutstring = "(mlg < 60 || mlg > 120)"
    else:
        zveto_cutstring = "true"

    if args.apply_2018_photon_phi_cut:
        photon_phi_2018_cutstring = "!(photon_phi > 0.55 && photon_phi < 0.7)"
    else:
        photon_phi_2018_cutstring = "true"

    if lep == "muon":
        if year == "2016":
            return "(pass_selection && " + photon_eta_cutstring+" && " + zveto_cutstring + " && " + puppimet_cutstring + " && abs(lepton_pdg_id) == 13 && photon_pt > 25 && lepton_pt > 26)"
        elif year == "2017":
            return "(pass_selection && " + photon_eta_cutstring+" && " + zveto_cutstring + " && " + puppimet_cutstring + " && abs(lepton_pdg_id) == 13 && photon_pt > 25 && lepton_pt > 30)"
        elif year == "2018":
            return "(pass_selection && " + photon_phi_2018_cutstring + " && " + photon_eta_cutstring+" && " + zveto_cutstring + " && " + puppimet_cutstring + " && abs(lepton_pdg_id) == 13 && photon_pt > 25 && lepton_pt > 26)"
        else:
            assert(0)
    elif lep == "electron":                
        if year == "2016":
            return "(pass_selection && " + photon_eta_cutstring+" && " + zveto_cutstring + " && " + puppimet_cutstring + " && abs(lepton_pdg_id) == 11 && photon_pt > 25 && lepton_pt > 30)"
        elif year == "2017":
            return "(pass_selection && " + photon_eta_cutstring+" && " + zveto_cutstring + " && " + puppimet_cutstring + " && abs(lepton_pdg_id) == 11 && photon_pt > 25 && lepton_pt > 35)"
        elif year == "2018":
            return "(pass_selection && " + photon_phi_2018_cutstring + " && " + photon_eta_cutstring+" && " + zveto_cutstring + " && " + puppimet_cutstring + " && abs(lepton_pdg_id) == 11 && photon_pt > 25 && lepton_pt > 35)"
        else:
            assert(0)
    elif lep == "both":    
        if year == "2016":
            return "(pass_selection && " + photon_eta_cutstring+" && " + zveto_cutstring + " && " + puppimet_cutstring + " && ((abs(lepton_pdg_id) == 13 && photon_pt > 25 && lepton_pt > 26) || (abs(lepton_pdg_id) == 11 && photon_pt > 25 && lepton_pt > 30)))"
        elif year == "2017":
            return "(pass_selection && " + photon_eta_cutstring+" && " + zveto_cutstring + " && " + puppimet_cutstring + " && ((abs(lepton_pdg_id) == 13 && photon_pt > 25 && lepton_pt > 30) || (abs(lepton_pdg_id) == 11 && photon_pt > 25 && lepton_pt > 35)))"
        elif year == "2018":
            return "(pass_selection && " + photon_phi_2018_cutstring + " && " + photon_eta_cutstring+" && " + zveto_cutstring + " && " + puppimet_cutstring + " && ((abs(lepton_pdg_id) == 13 && photon_pt > 25 && lepton_pt > 26) || (abs(lepton_pdg_id) == 11 && photon_pt > 25 && lepton_pt > 35)))"
        else:
            assert(0)
    else:
        assert(0)

def pass_json(run,lumi):

    if str(run) not in good_run_lumis.keys():
        return False

    for lumi_pair in good_run_lumis[str(run)]:
        if lumi < lumi_pair[1] and lumi > lumi_pair[0]:
            return True

    return False    

ROOT.gROOT.cd()

ROOT.ROOT.EnableImplicitMT(int(args.nthreads))

#when the TMinuit object is reused, the random seed is not reset after each fit, so the fit result can change when it is run on the same input 
ROOT.TMinuitMinimizer.UseStaticMinuit(False)

import wg_labels
labels = wg_labels.get_labels(args.workdir)

#    from wg_labels_wjets import labels
#from wg_labels_recoil_tree import labels

#mlg_fit_lower_bound = 10
#mlg_fit_upper_bound = 30

mlg_fit_lower_bound = 10
mlg_fit_upper_bound = 250
#mlg_fit_upper_bound = 200
#mlg_fit_upper_bound = 300
mlg_bin_width=2

#the first variable is for the ewdim6 analysis
#variables = ["photon_pt","dphilg","met","lepton_pt","lepton_eta","photon_pt","photon_eta","mlg","lepton_phi","photon_phi","njets40","mt","npvs","drlg"]
#variables_labels = ["ewdim6_photon_pt","dphilg","met","lepton_pt","lepton_eta","photon_pt","photon_eta","mlg","lepton_phi","photon_phi","njets40","mt","npvs","drlg"]

if args.make_all_plots:
    variables = ["detalg","dphilpuppimet","dphilg","puppimet","lepton_pt","lepton_eta","photon_eta","mlg","mlg_overflow","lepton_phi","photon_phi","njets40","puppimt","npvs","drlg","photon_pt","dphigpuppimet","puppimetphi","mlg","mlg","mlg","mlg","mlg","mlg"]
    variables_labels = ["detalg","dphilpuppimet","dphilg","puppimet","lepton_pt","lepton_eta","photon_eta","fit_mlg","mlg","lepton_phi","photon_phi","njets40","puppimt","npvs","drlg","photon_pt","dphigpuppimet","puppimetphi","mlg_large_bins","mlg_3bins","mlg_1bin","mlg_10bins","mlg_15bins","mlg_6bins"]
else:
#    variables = ["mlg_overflow"]
    variables = ["mlg"]
    variables_labels = ["fit_mlg"]

if args.ewdim6:
    variables.append("photon_pt_overflow")
    variables_labels.append("ewdim6_photon_pt")

assert(len(variables) == len(variables_labels))

from array import array

binning_photon_pt = array('f',[400,500,600,900,1500])
#binning_photon_pt = array('f',[300,500,750,1000,1500])
#binning_photon_pt = array('f',[100,200,300,400,500,600])

#binning_mlg = array('f',[10,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,200])
#binning_mlg = array('f',[10,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,200]
binning_mlg = array('f',[60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120])
#binning_mlg = array('f',[10,30,40,50,60,70,80,90,100,125,150,200,300])
#binning_mlg = array('f',[10,50,90,100,125,150,200,300])

n_photon_pt_bins = len(binning_photon_pt)-1

if args.make_all_plots:
    variable_definitions = [
        ["detalg" , "abs(lepton_eta-photon_eta)"],
        ["dphilpuppimet" , "abs(lepton_phi - puppimetphi) > TMath::Pi() ? abs(abs(lepton_phi - puppimetphi) - 2*TMath::Pi()) : abs(lepton_phi - puppimetphi)"],
        ["dphigpuppimet" , "abs(photon_phi - puppimetphi) > TMath::Pi() ? abs(abs(photon_phi - puppimetphi) - 2*TMath::Pi()) : abs(photon_phi - puppimetphi)"],
        ["dphilg" , "abs(lepton_phi - photon_phi) > TMath::Pi() ? abs(abs(lepton_phi - photon_phi) - 2*TMath::Pi()) : abs(lepton_phi - photon_phi)"],
        ["drlg" , "sqrt(dphilg*dphilg+detalg*detalg)" ],
#        ["photon_recoil","cos(photon_phi)*(-lepton_pt*cos(lepton_phi)-puppimet*cos(puppimetphi)) + sin(photon_phi)*(-lepton_pt*sin(lepton_phi) -puppimet*sin(puppimetphi))"],
        ["mlg_overflow","TMath::Min(float(mlg),float("+str(mlg_fit_upper_bound+mlg_bin_width/2.)+"))"]
    ]
else:
    variable_definitions = [ 
#        ["mlg_overflow","TMath::Min(float(mlg),float("+str(mlg_fit_upper_bound+mlg_bin_width/2.)+"))"],
]

if args.ewdim6:
    variable_definitions.append(["photon_pt_overflow","TMath::Min(float(photon_pt),float("+str(   (binning_photon_pt[n_photon_pt_bins] + binning_photon_pt[n_photon_pt_bins-1])/2) +"))"])


if args.make_all_plots:
    histogram_models = [
        ROOT.RDF.TH1DModel('','',50,0,5), #detalg
        ROOT.RDF.TH1DModel('','',48,0,pi), #dphilmet
        ROOT.RDF.TH1DModel('','',12,0,pi), #dphilg
        ROOT.RDF.TH1DModel("met", "", 40, 40., 200 ), 
        ROOT.RDF.TH1DModel('lepton_pt', '', 48, 20., 180 ), 
        ROOT.RDF.TH1DModel('lepton_eta', '', 50, -2.5, 2.5 ),
#        ROOT.RDF.TH1DModel('', '', n_photon_pt_bins, binning_photon_pt ), 
        ROOT.RDF.TH1DModel('photon_eta', '', 50, -2.5, 2.5 ), 
        #ROOT.RDF.TH1DModel("mlg","",mlg_fit_upper_bound/2,0,mlg_fit_upper_bound), 
        ROOT.RDF.TH1DModel("mlg","",(mlg_fit_upper_bound-mlg_fit_lower_bound)/mlg_bin_width,mlg_fit_lower_bound,mlg_fit_upper_bound),  
        ROOT.RDF.TH1DModel("mlg","",(mlg_fit_upper_bound-mlg_fit_lower_bound+mlg_bin_width)/mlg_bin_width,mlg_fit_lower_bound,mlg_fit_upper_bound+mlg_bin_width), 
        ROOT.RDF.TH1DModel("lepton_phi","",56,-3.5,3.5), 
        ROOT.RDF.TH1DModel("photon_phi","",56,-3.5,3.5), 
        ROOT.RDF.TH1DModel("","",7,-0.5,6.5), #njets40
#        ROOT.RDF.TH1DModel("mt","",10,0,200), 
        ROOT.RDF.TH1DModel("puppimt","",40,40,200), 
        ROOT.RDF.TH1DModel("npvs","",51,-0.5,50.5), 
        ROOT.RDF.TH1DModel("","",50,0,5), #drlg
        ROOT.RDF.TH1DModel('photon_pt', '', 48, 20., 180 ),
#        ROOT.RDF.TH1DModel("met", "", 15 , 0., 300 ),
#        ROOT.RDF.TH1DModel('photon_recoil', '', 20, -70., 130 ),
        ROOT.RDF.TH1DModel('','',48,0,pi), #dphigmet
        ROOT.RDF.TH1DModel("","",56,-3.5,3.5), #puppimetphi
        ROOT.RDF.TH1DModel("mlg","",30,0,300), #mlg
        ROOT.RDF.TH1DModel("mlg","",3,0,300), #mlg
        ROOT.RDF.TH1DModel("mlg","",1,0,300), #mlg
        ROOT.RDF.TH1DModel("mlg","",10,0,300), #mlg
        ROOT.RDF.TH1DModel("mlg","",15,0,300), #mlg
        ROOT.RDF.TH1DModel("mlg","",6,0,300), #mlg
#        ROOT.RDF.TH1DModel('', '', len(binning_mlg) - 1, binning_mlg ), #variable mlg binning,
#        ROOT.RDF.TH1DModel('photon_pt', '', 40, 100., 400 ),
#        ROOT.RDF.TH1DModel('photon_pt', '', 10, 300., 400 ),
] 
else:
    histogram_models = [
        ROOT.RDF.TH1DModel("mlg","",(mlg_fit_upper_bound-mlg_fit_lower_bound)/mlg_bin_width,mlg_fit_lower_bound,mlg_fit_upper_bound), 
#        ROOT.RDF.TH1DModel("mlg","",(mlg_fit_upper_bound-mlg_fit_lower_bound+mlg_bin_width)/mlg_bin_width,mlg_fit_lower_bound,mlg_fit_upper_bound+mlg_bin_width), 
    ]

if args.ewdim6:
    histogram_models.append(ROOT.RDF.TH1DModel('', '', n_photon_pt_bins, binning_photon_pt ))

assert(len(variables) == len(histogram_models))

if args.make_all_plots:
    mlg_index = 7
else:
    mlg_index = 0

ewdim6_index = len(histogram_models)-1

ewdim6_samples = {
"2016" : [{"xs" : 0.3097, "filename" : args.workdir+"/data/wg/2017/1June2019/wgjetsewdim6.root"}],
"2017" : [{"xs" : 0.3097, "filename" : args.workdir+"/data/wg/2017/1June2019/wgjetsewdim6.root"}],
"2018" : [{"xs" : 0.3097, "filename" : args.workdir+"/data/wg/2017/1June2019/wgjetsewdim6.root"}]
}

if args.ewdim6:
    for year in years:
        for sample in ewdim6_samples[year]:
            sample["file"] = ROOT.TFile.Open(sample["filename"])
            sample["nweightedevents"] = sample["file"].Get("nEventsGenWeighted").GetBinContent(1)

def getXaxisLabel(varname):
    if varname == "njets40":
        return "number of jets"
    elif varname == "detalg":
        return "#Delta#eta(l,#gamma)"
    elif varname == "dphilpuppimet":
        return "#Delta#phi(l,MET)"
    elif varname == "dphigpuppimet":
        return "#Delta#phi(#gamma,MET)"
    elif varname == "corrdphilpuppimet":
        return "corrected #Delta#phi(l,MET)"
    elif varname == "drlg":
        return "#Delta R(l,#gamma)"
    elif varname == "dphilg":
        return "#Delta#phi(l,#gamma)"
    elif varname == "npvs":
        return "number of PVs"
    elif varname == "mt":
        return "m_{t} (GeV)"
    elif varname == "puppimt":
        return "m_{t} (GeV)"
    elif varname == "corrmt":
        return "corrected m_{t} (GeV)"
    elif varname == "mlg":
        return "m_{lg} (GeV)"
    elif varname == "mlg_overflow":
        return "m_{lg} (GeV)"
    elif varname == "puppimet":
        return "MET (GeV)"
    elif varname == "met":
        return "MET (GeV)"
    elif varname == "corrmet":
        return "corrected MET (GeV)"
    elif varname == "lepton_pt":
        return "lepton p_{T} (GeV)"
    elif varname == "lepton_eta":
        return "lepton #eta"
    elif varname == "lepton_phi":
        return "lepton #phi"
    elif varname == "puppimetphi":
        return "MET #phi"
    elif varname == "photon_pt":
        return "photon p_{T} (GeV)"
    elif varname == "photon_pt_overflow":
        return "photon p_{T} (GeV)"
    elif varname == "photon_eta":
        return "photon #eta"    
    elif varname == "photon_phi":
        return "photon #phi"
    elif varname == "photon_recoil":
        return "photon recoil (GeV)"

    else:
        assert(0)

xoffsetstart = 0.0;
yoffsetstart = 0.0;
xoffset = 0.20;
yoffset = 0.05;

xpositions = [0.65,0.65,0.65,0.65,0.45,0.45,0.45,0.45,0.25,0.25,0.25,0.25]
ypositions = [0,1,2,3,0,1,2,3,0,1,2,3]

#xpositions = [0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.445,0.445,0.445,0.445,0.445,0.445,0.445,0.21,0.21,0.21,0.21,0.21,0.21,0.21]
#ypositions = [0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,6]

style.GoodStyle().cd()

def set_axis_fonts(thstack, coordinate, title):

    if coordinate == "x":
        axis = thstack.GetXaxis();
    elif coordinate == "y":
        axis = thstack.GetYaxis();
    else:
        assert(0)
    
    axis.SetLabelFont  (   42)
    axis.SetLabelOffset(0.015)
    axis.SetLabelSize  (0.050)
    axis.SetNdivisions (  505)
    axis.SetTitleFont  (   42)
    axis.SetTitleOffset(  1.5)
    axis.SetTitleSize  (0.050)
    if (coordinate == "y"):
        axis.SetTitleOffset(1.6)
    axis.SetTitle(title)    

def draw_legend(x1,y1,hist,label,options):

    legend = ROOT.TLegend(x1+xoffsetstart,y1+yoffsetstart,x1+xoffsetstart + xoffset,y1+yoffsetstart + yoffset)

    legend.SetBorderSize(     0)
    legend.SetFillColor (     0)
    legend.SetTextAlign (    12)
    legend.SetTextFont  (    42)
    legend.SetTextSize  ( 0.040)

    legend.AddEntry(hist,label,options)

    legend.Draw("same")

    #otherwise the legend goes out of scope and is deleted once the function finishes
    hist.label = legend

for label in labels.keys():

    labels[label]["hists"] = {}
    for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
        labels[label]["hists-"+unc+"-up"] = {}

    if labels[label]["syst-pdf"]:
        for i in range(0,32):
            labels[label]["hists-pdf-variation"+str(i)] = {}

    if labels[label]["syst-scale"]:
        for i in range(0,8):
            labels[label]["hists-scale-variation"+str(i)] = {}

    for i in range(len(variables)):    
        if labels[label]["color"] == None:
            continue

        labels[label]["hists"][i] = histogram_models[i].GetHistogram()
        labels[label]["hists"][i].SetName(label+" "+variables[i])
        labels[label]["hists"][i].Sumw2()
        for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
            labels[label]["hists-"+unc+"-up"][i] = histogram_models[i].GetHistogram()
            labels[label]["hists-"+unc+"-up"][i].Sumw2()

        if labels[label]["syst-pdf"]:
            for j in range(0,32):
                labels[label]["hists-pdf-variation"+str(j)][i] = histogram_models[i].GetHistogram()
                labels[label]["hists-pdf-variation"+str(j)][i].Sumw2()

        if labels[label]["syst-scale"]:
            for j in range(0,8):
                labels[label]["hists-scale-variation"+str(j)][i] = histogram_models[i].GetHistogram()
                labels[label]["hists-scale-variation"+str(j)][i].Sumw2()


    for year in years:            
        for sample in labels[label]["samples"][year]:
            sample["file"] = ROOT.TFile.Open(sample["filename"])
            sample["tree"] = sample["file"].Get("Events")
            sample["nweightedevents"] = sample["file"].Get("nEventsGenWeighted").GetBinContent(1)

if "w+jets" in labels:
    labels["w+jets"]["hists-prompt-pileup"] = {}
    for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
        if unc == "jes" or unc == "jer":
            continue
        labels["w+jets"]["hists-prompt-pileup-"+unc+"-up"] = {}

    for i in range(len(variables)):    
        labels["w+jets"]["hists-prompt-pileup"][i] = histogram_models[i].GetHistogram()
        labels["w+jets"]["hists-prompt-pileup"][i].Sumw2()
        for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
            if unc == "jes" or unc == "jer":
                continue
            labels["w+jets"]["hists-prompt-pileup-"+unc+"-up"][i] = histogram_models[i].GetHistogram()
            labels["w+jets"]["hists-prompt-pileup-"+unc+"-up"][i].Sumw2()

if "wg+jets" in labels:
    labels["wg+jets"]["hists-pass-fiducial"] = {}
    labels["wg+jets"]["hists-fail-fiducial"] = {}
    for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
        labels["wg+jets"]["hists-pass-fiducial-"+unc+"-up"] = {}
        labels["wg+jets"]["hists-fail-fiducial-"+unc+"-up"] = {}
    for i in range(n_fake_photon_alt):
        labels["wg+jets"]["hists-pass-fiducial-fake-photon-alt"+str(i)] = {}
        labels["wg+jets"]["hists-fail-fiducial-fake-photon-alt"+str(i)] = {}

    labels["wg+jets"]["hists-fake-photon-pass-fiducial"] = {}
    labels["wg+jets"]["hists-fake-photon-fail-fiducial"] = {}
    for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
        labels["wg+jets"]["hists-fake-photon-pass-fiducial-"+unc+"-up"] = {}
        labels["wg+jets"]["hists-fake-photon-fail-fiducial-"+unc+"-up"] = {}
    for i in range(n_fake_photon_alt):
        labels["wg+jets"]["hists-fake-photon-pass-fiducial-fake-photon-alt"+str(i)] = {}
        labels["wg+jets"]["hists-fake-photon-fail-fiducial-fake-photon-alt"+str(i)] = {}

    labels["wg+jets"]["hists-fake-lepton-pass-fiducial"] = {}
    labels["wg+jets"]["hists-fake-lepton-fail-fiducial"] = {}
    for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
        labels["wg+jets"]["hists-fake-lepton-pass-fiducial-"+unc+"-up"] = {}
        labels["wg+jets"]["hists-fake-lepton-fail-fiducial-"+unc+"-up"] = {}

    labels["wg+jets"]["hists-double-fake-pass-fiducial"] = {}
    labels["wg+jets"]["hists-double-fake-fail-fiducial"] = {}
    for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
        labels["wg+jets"]["hists-double-fake-pass-fiducial-"+unc+"-up"] = {}
        labels["wg+jets"]["hists-double-fake-fail-fiducial-"+unc+"-up"] = {}
    for i in range(n_fake_photon_alt):
        labels["wg+jets"]["hists-double-fake-pass-fiducial-fake-photon-alt"+str(i)] = {}
        labels["wg+jets"]["hists-double-fake-fail-fiducial-fake-photon-alt"+str(i)] = {}

    if labels["wg+jets"]["syst-pdf"]:
        for i in range(0,32):
            labels["wg+jets"]["hists-pass-fiducial-pdf-variation"+str(i)] = {}
            labels["wg+jets"]["hists-fail-fiducial-pdf-variation"+str(i)] = {}

    if labels["wg+jets"]["syst-scale"]:
        for i in range(0,8):
            labels["wg+jets"]["hists-pass-fiducial-scale-variation"+str(i)] = {}
            labels["wg+jets"]["hists-fail-fiducial-scale-variation"+str(i)] = {}

    for i in range(len(variables)):    
        labels["wg+jets"]["hists-pass-fiducial"][i] = histogram_models[i].GetHistogram()
        labels["wg+jets"]["hists-pass-fiducial"][i].Sumw2()
        labels["wg+jets"]["hists-fail-fiducial"][i] = histogram_models[i].GetHistogram()
        labels["wg+jets"]["hists-fail-fiducial"][i].Sumw2()
        for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
            labels["wg+jets"]["hists-pass-fiducial-"+unc+"-up"][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-pass-fiducial-"+unc+"-up"][i].Sumw2()
            labels["wg+jets"]["hists-fail-fiducial-"+unc+"-up"][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-fail-fiducial-"+unc+"-up"][i].Sumw2()
        for j in range(n_fake_photon_alt):
            labels["wg+jets"]["hists-pass-fiducial-fake-photon-alt"+str(j)][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-pass-fiducial-fake-photon-alt"+str(j)][i].Sumw2()
            labels["wg+jets"]["hists-fail-fiducial-fake-photon-alt"+str(j)][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-fail-fiducial-fake-photon-alt"+str(j)][i].Sumw2()

        labels["wg+jets"]["hists-fake-photon-pass-fiducial"][i] = histogram_models[i].GetHistogram()
        labels["wg+jets"]["hists-fake-photon-pass-fiducial"][i].Sumw2()
        labels["wg+jets"]["hists-fake-photon-fail-fiducial"][i] = histogram_models[i].GetHistogram()
        labels["wg+jets"]["hists-fake-photon-fail-fiducial"][i].Sumw2()
        for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
            labels["wg+jets"]["hists-fake-photon-pass-fiducial-"+unc+"-up"][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-fake-photon-pass-fiducial-"+unc+"-up"][i].Sumw2()
            labels["wg+jets"]["hists-fake-photon-fail-fiducial-"+unc+"-up"][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-fake-photon-fail-fiducial-"+unc+"-up"][i].Sumw2()
        for j in range(n_fake_photon_alt):
            labels["wg+jets"]["hists-fake-photon-pass-fiducial-fake-photon-alt"+str(j)][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-fake-photon-pass-fiducial-fake-photon-alt"+str(j)][i].Sumw2()
            labels["wg+jets"]["hists-fake-photon-fail-fiducial-fake-photon-alt"+str(j)][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-fake-photon-fail-fiducial-fake-photon-alt"+str(j)][i].Sumw2()

        labels["wg+jets"]["hists-fake-lepton-pass-fiducial"][i] = histogram_models[i].GetHistogram()
        labels["wg+jets"]["hists-fake-lepton-pass-fiducial"][i].Sumw2()
        labels["wg+jets"]["hists-fake-lepton-fail-fiducial"][i] = histogram_models[i].GetHistogram()
        labels["wg+jets"]["hists-fake-lepton-fail-fiducial"][i].Sumw2()
        for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
            labels["wg+jets"]["hists-fake-lepton-pass-fiducial-"+unc+"-up"][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-fake-lepton-pass-fiducial-"+unc+"-up"][i].Sumw2()
            labels["wg+jets"]["hists-fake-lepton-fail-fiducial-"+unc+"-up"][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-fake-lepton-fail-fiducial-"+unc+"-up"][i].Sumw2()

        labels["wg+jets"]["hists-double-fake-pass-fiducial"][i] = histogram_models[i].GetHistogram()
        labels["wg+jets"]["hists-double-fake-pass-fiducial"][i].Sumw2()
        labels["wg+jets"]["hists-double-fake-fail-fiducial"][i] = histogram_models[i].GetHistogram()
        labels["wg+jets"]["hists-double-fake-fail-fiducial"][i].Sumw2()
        for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
            labels["wg+jets"]["hists-double-fake-pass-fiducial-"+unc+"-up"][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-double-fake-pass-fiducial-"+unc+"-up"][i].Sumw2()
            labels["wg+jets"]["hists-double-fake-fail-fiducial-"+unc+"-up"][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-double-fake-fail-fiducial-"+unc+"-up"][i].Sumw2()
        for j in range(n_fake_photon_alt):
            labels["wg+jets"]["hists-double-fake-pass-fiducial-fake-photon-alt"+str(j)][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-double-fake-pass-fiducial-fake-photon-alt"+str(j)][i].Sumw2()
            labels["wg+jets"]["hists-double-fake-fail-fiducial-fake-photon-alt"+str(j)][i] = histogram_models[i].GetHistogram()
            labels["wg+jets"]["hists-double-fake-fail-fiducial-fake-photon-alt"+str(j)][i].Sumw2()

        if labels["wg+jets"]["syst-pdf"]:
            for j in range(0,32):
                labels["wg+jets"]["hists-pass-fiducial-pdf-variation"+str(j)][i] = histogram_models[i].GetHistogram()
                labels["wg+jets"]["hists-pass-fiducial-pdf-variation"+str(j)][i].Sumw2()
                labels["wg+jets"]["hists-fail-fiducial-pdf-variation"+str(j)][i] = histogram_models[i].GetHistogram()
                labels["wg+jets"]["hists-fail-fiducial-pdf-variation"+str(j)][i].Sumw2()
                
        if labels["wg+jets"]["syst-scale"]:
            for j in range(0,8):
                labels["wg+jets"]["hists-pass-fiducial-scale-variation"+str(j)][i] = histogram_models[i].GetHistogram()
                labels["wg+jets"]["hists-pass-fiducial-scale-variation"+str(j)][i].Sumw2()
                labels["wg+jets"]["hists-fail-fiducial-scale-variation"+str(j)][i] = histogram_models[i].GetHistogram()
                labels["wg+jets"]["hists-fail-fiducial-scale-variation"+str(j)][i].Sumw2()

    for year in years:            
        for sample in labels["wg+jets"]["samples"][year]:
            sample["nweightedevents_passfiducial"] = sample["file"].Get("nEventsGenWeighted_PassFidSelection").GetBinContent(1)

        if labels["wg+jets"]["syst-scale"]:
            for i in range(0,8):
                labels["wg+jets"]["samples"][year][0]["nweightedevents_qcdscaleweight"+str(i)]=labels["wg+jets"]["samples"][year][0]["file"].Get("nEventsGenWeighted_PassFidSelection_QCDScaleWeight"+str(i)).GetBinContent(1)

                if labels["wg+jets"]["samples"][year][0]["filename"] == args.workdir+"/data/wg/2016/1June2019/wgjets.root" or labels["wg+jets"]["samples"][year][0]["filename"] == args.workdir+"/data/wg/2016/1June2019jetunc/wgjets.root":
                    labels["wg+jets"]["samples"][year][0]["nweightedevents_qcdscaleweight"+str(i)] *= 2
                    
        if labels["wg+jets"]["syst-pdf"]:
            for i in range(1,32):
                labels["wg+jets"]["samples"][year][0]["nweightedevents_pdfweight"+str(i)]=labels["wg+jets"]["samples"][year][0]["file"].Get("nEventsGenWeighted_PassFidSelection_PDFWeight"+str(i)).GetBinContent(1)

                if labels["wg+jets"]["samples"][year][0]["filename"] == args.workdir+"/data/wg/2016/1June2019/wgjets.root" or labels["wg+jets"]["samples"][year][0]["filename"] == args.workdir+"/data/wg/2016/1June2019jet/wgjets.root":
                    labels["wg+jets"]["samples"][year][0]["nweightedevents_pdfweight"+str(i)] *= 2


#    for year in years:
#        labels["wg+jets"]["samples"][year][0]["nweightedeventspassgenselection"]=labels["wg+jets"]["samples"][year][0]["file"].Get("nWeightedEventsPassGenSelection").GetBinContent(1)
    #labels["wg+jets"]["samples"][year][0]["nweightedeventspassgenselection"]=1

    nweightedeventspassgenselection=0
    nweightedevents = 0
    for year in years:

        lumi = dict_lumi[year]

        nweightedeventspassgenselection+=labels["wg+jets"]["samples"][year][0]["nweightedevents_passfiducial"]*lumi
        nweightedevents+=labels["wg+jets"]["samples"][year][0]["nweightedevents"]*lumi

    fiducial_region_cuts_efficiency = nweightedeventspassgenselection/nweightedevents

#etopbinning = [25,35,45]
etopbinning = [25]

data = {}
fake_signal_contamination = {}

wjets = {}
wjets_fake_photon = {}
wjets_fake_photon_2016 = {}
wjets_fake_photon_chiso_2016 = {}
fake_photon = {}
fake_photon_2016 = {}
wjets_2016 = {}
fake_photon_stat_up = {}
fake_lepton = {}
fake_lepton_stat_down = {}
fake_lepton_stat_up = {}
double_fake = {}
double_fake_stat_up = {}
e_to_p_total = {}
e_to_p = []
for i in range(len(etopbinning)):
    e_to_p.append({})
ewdim6 = {}

data["hists"] = []
fake_signal_contamination["hists"] = []
wjets_fake_photon_2016["hists"] = []
wjets_fake_photon_chiso_2016["hists"] = []
wjets_2016["hists"] = []
fake_photon["hists"] = []
for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
    fake_photon["hists-"+unc+"-up"] = []
for i in range(n_fake_photon_alt):
    fake_photon["hists-alt"+str(i)] = []
fake_photon_2016["hists"] = []
fake_photon_stat_up["hists"] = []
fake_lepton["hists"] = []
for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
    fake_lepton["hists-"+unc+"-up"] = []
fake_lepton_stat_down["hists"] = []
fake_lepton_stat_up["hists"] = []
double_fake["hists"] = []
for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
    double_fake["hists-"+unc+"-up"] = []
for i in range(n_fake_photon_alt):
    double_fake["hists-alt"+str(i)] = []
double_fake_stat_up["hists"] = []
for i in range(len(etopbinning)):
    e_to_p[i]["hists"] = []
    e_to_p[i]["hists-electron-id-sf-up"] = []
    e_to_p[i]["hists-electron-reco-sf-up"] = []
    e_to_p[i]["hists-electron-hlt-sf-up"] = []
    e_to_p[i]["hists-photon-id-sf-up"] = []
    e_to_p[i]["hists-pileup-up"] = []
    e_to_p[i]["hists-prefire-up"] = []
    e_to_p[i]["hists-jes-up"] = []
    e_to_p[i]["hists-jer-up"] = []
e_to_p_total["hists"] = []
e_to_p_total["hists-electron-id-sf-up"] = []
e_to_p_total["hists-electron-reco-sf-up"] = []
e_to_p_total["hists-electron-hlt-sf-up"] = []
e_to_p_total["hists-photon-id-sf-up"] = []
e_to_p_total["hists-pileup-up"] = []
e_to_p_total["hists-prefire-up"] = []
e_to_p_total["hists-jes-up"] = []
e_to_p_total["hists-jer-up"] = []
ewdim6["hists"] = []

for i in range(len(variables)):
    data["hists"].append(histogram_models[i].GetHistogram())
    wjets_fake_photon_2016["hists"].append(histogram_models[i].GetHistogram())
    wjets_fake_photon_chiso_2016["hists"].append(histogram_models[i].GetHistogram())
    wjets_2016["hists"].append(histogram_models[i].GetHistogram())

    fake_photon["hists"].append(histogram_models[i].GetHistogram())
    for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
        fake_photon["hists-"+unc+"-up"].append(histogram_models[i].GetHistogram())
    for j in range(n_fake_photon_alt):
        fake_photon["hists-alt"+str(j)].append(histogram_models[i].GetHistogram())
    fake_photon_2016["hists"].append(histogram_models[i].GetHistogram())
    fake_photon_stat_up["hists"].append(histogram_models[i].GetHistogram())

    fake_lepton["hists"].append(histogram_models[i].GetHistogram())
    for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
        fake_lepton["hists-"+unc+"-up"].append(histogram_models[i].GetHistogram())
    fake_lepton_stat_up["hists"].append(histogram_models[i].GetHistogram())
    fake_lepton_stat_down["hists"].append(histogram_models[i].GetHistogram())

    double_fake["hists"].append(histogram_models[i].GetHistogram())
    for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
        double_fake["hists-"+unc+"-up"].append(histogram_models[i].GetHistogram())
    for j in range(n_fake_photon_alt):
        double_fake["hists-alt"+str(j)].append(histogram_models[i].GetHistogram())
    double_fake_stat_up["hists"].append(histogram_models[i].GetHistogram())

    for j in range(len(etopbinning)):
        e_to_p[j]["hists"].append(histogram_models[i].GetHistogram())
        e_to_p[j]["hists-electron-id-sf-up"].append(histogram_models[i].GetHistogram())
        e_to_p[j]["hists-electron-reco-sf-up"].append(histogram_models[i].GetHistogram())
        e_to_p[j]["hists-electron-hlt-sf-up"].append(histogram_models[i].GetHistogram())
        e_to_p[j]["hists-photon-id-sf-up"].append(histogram_models[i].GetHistogram())
        e_to_p[j]["hists-pileup-up"].append(histogram_models[i].GetHistogram())
        e_to_p[j]["hists-prefire-up"].append(histogram_models[i].GetHistogram())
        e_to_p[j]["hists-jes-up"].append(histogram_models[i].GetHistogram())
        e_to_p[j]["hists-jer-up"].append(histogram_models[i].GetHistogram())
    e_to_p_total["hists"].append(histogram_models[i].GetHistogram())
    e_to_p_total["hists-electron-id-sf-up"].append(histogram_models[i].GetHistogram())
    e_to_p_total["hists-electron-reco-sf-up"].append(histogram_models[i].GetHistogram())
    e_to_p_total["hists-electron-hlt-sf-up"].append(histogram_models[i].GetHistogram())
    e_to_p_total["hists-photon-id-sf-up"].append(histogram_models[i].GetHistogram())
    e_to_p_total["hists-pileup-up"].append(histogram_models[i].GetHistogram())
    e_to_p_total["hists-prefire-up"].append(histogram_models[i].GetHistogram())
    e_to_p_total["hists-jes-up"].append(histogram_models[i].GetHistogram())
    e_to_p_total["hists-jer-up"].append(histogram_models[i].GetHistogram())
    fake_signal_contamination["hists"].append(histogram_models[i].GetHistogram())
    ewdim6["hists"].append(histogram_models[i].GetHistogram())

for i in range(len(variables)):
    data["hists"][i].Sumw2()
    data["hists"][i].SetName("data "+variables[i])
    wjets_fake_photon_2016["hists"][i].Sumw2()
    wjets_fake_photon_chiso_2016["hists"][i].Sumw2()
    wjets_2016["hists"][i].Sumw2()
    fake_photon["hists"][i].Sumw2()
    fake_photon["hists-electron-id-sf-up"][i].Sumw2()
    fake_photon["hists-electron-reco-sf-up"][i].Sumw2()
    fake_photon["hists-electron-hlt-sf-up"][i].Sumw2()
    fake_photon["hists-muon-id-sf-up"][i].Sumw2()
    fake_photon["hists-muon-iso-sf-up"][i].Sumw2()
    fake_photon["hists-muon-hlt-sf-up"][i].Sumw2()
    fake_photon["hists-photon-id-sf-up"][i].Sumw2()
    fake_photon["hists-pileup-up"][i].Sumw2()
    fake_photon["hists-prefire-up"][i].Sumw2()
    fake_photon["hists-jes-up"][i].Sumw2()
    fake_photon["hists-jer-up"][i].Sumw2()
    for j in range(n_fake_photon_alt):
        fake_photon["hists-alt"+str(j)][i].Sumw2()
    fake_photon_2016["hists"][i].Sumw2()
    fake_photon["hists"][i].SetName("fake photon "+variables[i])
    fake_photon_2016["hists"][i].SetName("fake photon 2016 "+variables[i])
    fake_photon_stat_up["hists"][i].Sumw2()
    fake_lepton["hists"][i].Sumw2()
    fake_lepton["hists-electron-id-sf-up"][i].Sumw2()
    fake_lepton["hists-electron-reco-sf-up"][i].Sumw2()
    fake_lepton["hists-electron-hlt-sf-up"][i].Sumw2()
    fake_lepton["hists-muon-id-sf-up"][i].Sumw2()
    fake_lepton["hists-muon-iso-sf-up"][i].Sumw2()
    fake_lepton["hists-muon-hlt-sf-up"][i].Sumw2()
    fake_lepton["hists-photon-id-sf-up"][i].Sumw2()
    fake_lepton["hists-pileup-up"][i].Sumw2()
    fake_lepton["hists-prefire-up"][i].Sumw2()
    fake_lepton["hists-jes-up"][i].Sumw2()
    fake_lepton["hists-jer-up"][i].Sumw2()
    fake_lepton_stat_up["hists"][i].Sumw2()
    fake_lepton_stat_down["hists"][i].Sumw2()
    double_fake["hists"][i].Sumw2()
    double_fake["hists-electron-id-sf-up"][i].Sumw2()
    double_fake["hists-electron-reco-sf-up"][i].Sumw2()
    double_fake["hists-electron-hlt-sf-up"][i].Sumw2()
    double_fake["hists-muon-id-sf-up"][i].Sumw2()
    double_fake["hists-muon-iso-sf-up"][i].Sumw2()
    double_fake["hists-muon-hlt-sf-up"][i].Sumw2()
    double_fake["hists-photon-id-sf-up"][i].Sumw2()
    double_fake["hists-pileup-up"][i].Sumw2()
    double_fake["hists-prefire-up"][i].Sumw2()
    double_fake["hists-jes-up"][i].Sumw2()
    double_fake["hists-jer-up"][i].Sumw2()
    for j in range(n_fake_photon_alt):
        double_fake["hists-alt"+str(j)][i].Sumw2()
    double_fake_stat_up["hists"][i].Sumw2()
    for j in range(len(etopbinning)):
        e_to_p[j]["hists"][i].Sumw2()
        e_to_p[j]["hists-electron-id-sf-up"][i].Sumw2()
        e_to_p[j]["hists-electron-reco-sf-up"][i].Sumw2()
        e_to_p[j]["hists-electron-hlt-sf-up"][i].Sumw2()
        e_to_p[j]["hists-photon-id-sf-up"][i].Sumw2()
        e_to_p[j]["hists-pileup-up"][i].Sumw2()
        e_to_p[j]["hists-prefire-up"][i].Sumw2()
        e_to_p[j]["hists-jes-up"][i].Sumw2()
        e_to_p[j]["hists-jer-up"][i].Sumw2()
    e_to_p_total["hists"][i].Sumw2()
    e_to_p_total["hists-electron-id-sf-up"][i].Sumw2()
    e_to_p_total["hists-electron-reco-sf-up"][i].Sumw2()
    e_to_p_total["hists-electron-hlt-sf-up"][i].Sumw2()
    e_to_p_total["hists-photon-id-sf-up"][i].Sumw2()
    e_to_p_total["hists-pileup-up"][i].Sumw2()
    e_to_p_total["hists-prefire-up"][i].Sumw2()
    e_to_p_total["hists-jes-up"][i].Sumw2()
    e_to_p_total["hists-jer-up"][i].Sumw2()
    ewdim6["hists"][i].Sumw2()
    fake_signal_contamination["hists"][i].Sumw2()

ROOT.gROOT.cd()

include_headers_cpp = '''

#include <iostream>
#include <fstream>
#include <map>
#include "/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/json/3.7.3/include/nlohmann/json.hpp"

'''

wjets_pileup_photons_flags_initialization_cpp = '''

nlohmann::json wjets_2016v1_json;
nlohmann::json wjets_2016v2_json;
nlohmann::json wjets_2017v1_json;
nlohmann::json wjets_2017v2_json;
nlohmann::json wjets_2017v3_json;
nlohmann::json wjets_2018_json;

std::ifstream wjets_2016v1_infile("pileup_photons/wjets_2016v1_events_info.txt");
std::ifstream wjets_2016v2_infile("pileup_photons/wjets_2016v2_events_info.txt");
std::ifstream wjets_2017v1_infile("pileup_photons/wjets_2017v1_events_info.txt");
std::ifstream wjets_2017v2_infile("pileup_photons/wjets_2017v2_events_info.txt");
std::ifstream wjets_2017v3_infile("pileup_photons/wjets_2017v3_events_info.txt");
std::ifstream wjets_2018_infile("pileup_photons/wjets_2018_events_info.txt");

wjets_2016v1_infile >> wjets_2016v1_json;
wjets_2016v2_infile >> wjets_2016v2_json;
wjets_2017v1_infile >> wjets_2017v1_json;
wjets_2017v2_infile >> wjets_2017v2_json;
wjets_2017v3_infile >> wjets_2017v3_json;
wjets_2018_infile >> wjets_2018_json;

std::map<std::pair<int,int>,bool> wjets_2016v1_prompt_pileup_photon_map;
std::map<std::pair<int,int>,bool> wjets_2016v2_prompt_pileup_photon_map;
std::map<std::pair<int,int>,bool> wjets_2017v1_prompt_pileup_photon_map;
std::map<std::pair<int,int>,bool> wjets_2017v2_prompt_pileup_photon_map;
std::map<std::pair<int,int>,bool> wjets_2017v3_prompt_pileup_photon_map;
std::map<std::pair<int,int>,bool> wjets_2018_prompt_pileup_photon_map;

for (int k = 0; k < wjets_2016v1_json.size(); ++k) {
    wjets_2016v1_prompt_pileup_photon_map[std::make_pair(int(wjets_2016v1_json[k]["lumi"]),int(wjets_2016v1_json[k]["event"]))] = wjets_2016v1_json[k]["prompt"];
}

for (int k = 0; k < wjets_2016v2_json.size(); ++k) {
    wjets_2016v2_prompt_pileup_photon_map[std::make_pair(int(wjets_2016v2_json[k]["lumi"]),int(wjets_2016v2_json[k]["event"]))] = wjets_2016v2_json[k]["prompt"];
}

for (int k = 0; k < wjets_2017v1_json.size(); ++k) {
    wjets_2017v1_prompt_pileup_photon_map[std::make_pair(int(wjets_2017v1_json[k]["lumi"]),int(wjets_2017v1_json[k]["event"]))] = wjets_2017v1_json[k]["prompt"];
}

for (int k = 0; k < wjets_2017v2_json.size(); ++k) {
    wjets_2017v2_prompt_pileup_photon_map[std::make_pair(int(wjets_2017v2_json[k]["lumi"]),int(wjets_2017v2_json[k]["event"]))] = wjets_2017v2_json[k]["prompt"];
}

for (int k = 0; k < wjets_2017v3_json.size(); ++k) {
    wjets_2017v3_prompt_pileup_photon_map[std::make_pair(int(wjets_2017v3_json[k]["lumi"]),int(wjets_2017v3_json[k]["event"]))] = wjets_2017v3_json[k]["prompt"];
}

for (int k = 0; k < wjets_2018_json.size(); ++k) {
    wjets_2018_prompt_pileup_photon_map[std::make_pair(int(wjets_2018_json[k]["lumi"]),int(wjets_2018_json[k]["event"]))] = wjets_2018_json[k]["prompt"];
}

'''

wjets_pileup_photons_flags_cpp = '''


bool is_photon_prompt(int lumi,int event, string year, string dsetversion) {

    if (year == "2016" && dsetversion == "v1") {

        if (wjets_2016v1_prompt_pileup_photon_map.find(std::make_pair(lumi,event)) == wjets_2016v1_prompt_pileup_photon_map.end()){
            std::cout << "(lumi,event) = (" << lumi << "," << event << ") not in map" << std::endl;
            exit(1);
        }

        if (wjets_2016v1_prompt_pileup_photon_map[std::make_pair(lumi,event)])
            return true;
    }
    else if (year == "2016" && dsetversion == "v2") {

        if (wjets_2016v2_prompt_pileup_photon_map.find(std::make_pair(lumi,event)) == wjets_2016v2_prompt_pileup_photon_map.end()){
            std::cout << "(lumi,event) = (" << lumi << "," << event << ") not in map" << std::endl;
            exit(1);
        }

        if (wjets_2016v2_prompt_pileup_photon_map[std::make_pair(lumi,event)])
            return true;
    }
    else if (year == "2017" && dsetversion == "v1") {

        if (wjets_2017v1_prompt_pileup_photon_map.find(std::make_pair(lumi,event)) == wjets_2017v1_prompt_pileup_photon_map.end()){
            std::cout << "(lumi,event) = (" << lumi << "," << event << ") not in map" << std::endl;
            exit(1);
        }

        if (wjets_2017v1_prompt_pileup_photon_map[std::make_pair(lumi,event)])
            return true;
    }
    else if (year == "2017" && dsetversion == "v2") {

        if (wjets_2017v2_prompt_pileup_photon_map.find(std::make_pair(lumi,event)) == wjets_2017v2_prompt_pileup_photon_map.end()){
            std::cout << "(lumi,event) = (" << lumi << "," << event << ") not in map" << std::endl;
            exit(1);
        }

        if (wjets_2017v2_prompt_pileup_photon_map[std::make_pair(lumi,event)])
            return true;
    }
    else if (year == "2017" && dsetversion == "v3") {

        if (wjets_2017v3_prompt_pileup_photon_map.find(std::make_pair(lumi,event)) == wjets_2017v3_prompt_pileup_photon_map.end()){
            std::cout << "(lumi,event) = (" << lumi << "," << event << ") not in map" << std::endl;
            exit(1);
        }

        if (wjets_2017v3_prompt_pileup_photon_map[std::make_pair(lumi,event)])
            return true;
    }
    else if (year == "2018" && dsetversion == "") {

        if (wjets_2018_prompt_pileup_photon_map.find(std::make_pair(lumi,event)) == wjets_2018_prompt_pileup_photon_map.end()){
            std::cout << "(lumi,event) = (" << lumi << "," << event << ") not in map" << std::endl;
            exit(1);
        }

        if (wjets_2018_prompt_pileup_photon_map[std::make_pair(lumi,event)])
            return true;
    }
    else
        exit(1);

    return false;
}

'''


photon_eff_scale_factor_cpp = '''

TFile photon_id_2016_sf_file("eff_scale_factors/2016/Fall17V2_2016_Medium_photons.root");
TH2F * photon_id_2016_sf = (TH2F*) photon_id_2016_sf_file.Get("EGamma_SF2D");

TFile photon_id_2017_sf_file("eff_scale_factors/2017/2017_PhotonsMedium.root");
TH2F * photon_id_2017_sf = (TH2F*) photon_id_2017_sf_file.Get("EGamma_SF2D");

TFile photon_id_2018_sf_file("eff_scale_factors/2018/2018_PhotonsMedium.root","read");
TH2F * photon_id_2018_sf = (TH2F*) photon_id_2018_sf_file.Get("EGamma_SF2D");

float photon_efficiency_scale_factor(float pt,float eta,string year,bool err_up=false){

    TH2F * photon_id_sf = 0;

    if (year == "2016") photon_id_sf = photon_id_2016_sf;
    else if (year == "2017") photon_id_sf = photon_id_2017_sf;
    else if (year == "2018") photon_id_sf = photon_id_2018_sf;
    else exit(1);

    float mypt = TMath::Min(pt,float(photon_id_sf->GetYaxis()->GetBinCenter(photon_id_sf->GetNbinsY())));
    float myeta = TMath::Max(TMath::Min(eta,float(photon_id_sf->GetXaxis()->GetBinCenter(photon_id_sf->GetNbinsX()))),float(photon_id_sf->GetXaxis()->GetBinCenter(1)));

    float sf = photon_id_sf->GetBinContent(photon_id_sf->GetXaxis()->FindFixBin(myeta),photon_id_sf->GetYaxis()->FindFixBin(mypt));

    if (err_up) sf += photon_id_sf->GetBinError(photon_id_sf->GetXaxis()->FindFixBin(myeta),photon_id_sf->GetYaxis()->FindFixBin(mypt));

    return sf;
}


'''

if not args.use_yings_electron_hlt_sfs:
    electron_eff_scale_factor_cpp = '''

TFile electron_id_2016_sf_file("eff_scale_factors/2016/2016LegacyReReco_ElectronMedium_Fall17V2.root","read");
TH2F * electron_id_2016_sf = (TH2F*) electron_id_2016_sf_file.Get("EGamma_SF2D");

TFile electron_id_2017_sf_file("eff_scale_factors/2017/2017_ElectronMedium.root","read");
TH2F * electron_id_2017_sf = (TH2F*)electron_id_2017_sf_file.Get("EGamma_SF2D");

TFile electron_id_2018_sf_file("eff_scale_factors/2018/2018_ElectronMedium.root","read");
TH2F * electron_id_2018_sf = (TH2F*)electron_id_2018_sf_file.Get("EGamma_SF2D");

TFile electron_reco_2016_sf_file("eff_scale_factors/2016/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root","read");
TH2F * electron_reco_2016_sf = (TH2F*) electron_reco_2016_sf_file.Get("EGamma_SF2D");

TFile electron_reco_2017_sf_file("eff_scale_factors/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root","read");
TH2F * electron_reco_2017_sf = (TH2F*)electron_reco_2017_sf_file.Get("EGamma_SF2D");

TFile electron_reco_2018_sf_file("eff_scale_factors/2018/egammaEffi.txt_EGM2D_updatedAll.root" ,"read");
TH2F * electron_reco_2018_sf = (TH2F*)electron_reco_2018_sf_file.Get("EGamma_SF2D");

TFile electron_hlt_2016_sf_file("eff_scale_factors/2016/electron_hlt_sfs_2016.root","read");
TH2D * electron_hlt_2016_sf = (TH2D*) electron_hlt_2016_sf_file.Get("hlt_sfs_etapt");

TFile electron_hlt_2017_sf_file("eff_scale_factors/2017/electron_hlt_sfs_2017.root","read");
TH2D * electron_hlt_2017_sf = (TH2D*)electron_hlt_2017_sf_file.Get("hlt_sfs_etapt");

TFile electron_hlt_2018_sf_file("eff_scale_factors/2018/electron_hlt_sfs_2018.root" ,"read");
TH2D * electron_hlt_2018_sf = (TH2D*)electron_hlt_2018_sf_file.Get("hlt_sfs_etapt");

float electron_efficiency_scale_factor(float pt, float eta, string year,bool id_err_up=false, bool reco_err_up=false, bool hlt_err_up=false) {

    TH2F * electron_reco_sf = 0;
    TH2F * electron_id_sf = 0;
    TH2D * electron_hlt_sf = 0;

    if (year == "2016") {
        electron_reco_sf = electron_reco_2016_sf;
        electron_id_sf = electron_id_2016_sf;
        electron_hlt_sf = electron_hlt_2016_sf;
    }
    else if (year == "2017"){
        electron_reco_sf = electron_reco_2017_sf;
        electron_id_sf = electron_id_2017_sf;
        electron_hlt_sf = electron_hlt_2017_sf;
    }
    else if (year == "2018") {
        electron_reco_sf = electron_reco_2018_sf;
        electron_id_sf = electron_id_2018_sf;
        electron_hlt_sf = electron_hlt_2018_sf;
    }
    else
        exit(1);

    int electron_id_sf_xaxisbin = -1;
    int electron_id_sf_yaxisbin = -1;

    if (year == "2016") {    
        electron_id_sf_xaxisbin = electron_id_sf->GetXaxis()->FindFixBin(eta);
        electron_id_sf_yaxisbin = electron_id_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(electron_id_sf->GetYaxis()->GetBinCenter(electron_id_sf->GetNbinsY()))));
    }
    else if (year == "2017") {
        electron_id_sf_xaxisbin = electron_id_sf->GetXaxis()->FindFixBin(eta);
        electron_id_sf_yaxisbin = electron_id_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(electron_id_sf->GetYaxis()->GetBinCenter(electron_id_sf->GetNbinsY()))));
    }
    else if (year == "2018") {
        electron_id_sf_xaxisbin = electron_id_sf->GetXaxis()->FindFixBin(eta);
        electron_id_sf_yaxisbin = electron_id_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(electron_id_sf->GetYaxis()->GetBinCenter(electron_id_sf->GetNbinsY()))));
    }
    else exit(1);

    int electron_hlt_sf_xaxisbin = -1;
    int electron_hlt_sf_yaxisbin = -1;

    if (year == "2016") {    
        electron_hlt_sf_xaxisbin = electron_hlt_sf->GetXaxis()->FindFixBin(eta);
        electron_hlt_sf_yaxisbin = electron_hlt_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(electron_hlt_sf->GetYaxis()->GetBinCenter(electron_hlt_sf->GetNbinsY()))));
    }
    else if (year == "2017") {
        electron_hlt_sf_xaxisbin = electron_hlt_sf->GetXaxis()->FindFixBin(eta);
        electron_hlt_sf_yaxisbin = electron_hlt_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(electron_hlt_sf->GetYaxis()->GetBinCenter(electron_hlt_sf->GetNbinsY()))));
    }
    else if (year == "2018") {
        electron_hlt_sf_xaxisbin = electron_hlt_sf->GetXaxis()->FindFixBin(eta);
        electron_hlt_sf_yaxisbin = electron_hlt_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(electron_hlt_sf->GetYaxis()->GetBinCenter(electron_hlt_sf->GetNbinsY()))));
    }
    else exit(1);


    float sf_id = electron_id_sf->GetBinContent(electron_id_sf_xaxisbin,electron_id_sf_yaxisbin); 
    if (id_err_up) sf_id += electron_id_sf->GetBinError(electron_id_sf_xaxisbin,electron_id_sf_yaxisbin) ;

    //the reco 2D histogram is really a 1D histogram
    float sf_reco=electron_reco_sf->GetBinContent(electron_reco_sf->GetXaxis()->FindFixBin(eta),1);
    if (reco_err_up) sf_reco+=electron_reco_sf->GetBinError(electron_reco_sf->GetXaxis()->FindFixBin(eta),1);

    float sf_hlt = electron_hlt_sf->GetBinContent(electron_hlt_sf_xaxisbin,electron_hlt_sf_yaxisbin); 
    if (hlt_err_up) sf_hlt += electron_hlt_sf->GetBinError(electron_hlt_sf_xaxisbin,electron_hlt_sf_yaxisbin) ;

    return sf_id*sf_reco*sf_hlt;
}


'''
else:
    electron_eff_scale_factor_cpp = '''

TFile electron_id_2016_sf_file("eff_scale_factors/2016/2016LegacyReReco_ElectronMedium_Fall17V2.root","read");
TH2F * electron_id_2016_sf = (TH2F*) electron_id_2016_sf_file.Get("EGamma_SF2D");

TFile electron_id_2017_sf_file("eff_scale_factors/2017/2017_ElectronMedium.root","read");
TH2F * electron_id_2017_sf = (TH2F*)electron_id_2017_sf_file.Get("EGamma_SF2D");

TFile electron_id_2018_sf_file("eff_scale_factors/2018/2018_ElectronMedium.root","read");
TH2F * electron_id_2018_sf = (TH2F*)electron_id_2018_sf_file.Get("EGamma_SF2D");

TFile electron_reco_2016_sf_file("eff_scale_factors/2016/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root","read");
TH2F * electron_reco_2016_sf = (TH2F*) electron_reco_2016_sf_file.Get("EGamma_SF2D");

TFile electron_reco_2017_sf_file("eff_scale_factors/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root","read");
TH2F * electron_reco_2017_sf = (TH2F*)electron_reco_2017_sf_file.Get("EGamma_SF2D");

TFile electron_reco_2018_sf_file("eff_scale_factors/2018/egammaEffi.txt_EGM2D_updatedAll.root" ,"read");
TH2F * electron_reco_2018_sf = (TH2F*)electron_reco_2018_sf_file.Get("EGamma_SF2D");

TFile electron_hlt_2016_sf_file("eff_scale_factors/2016/egammaEffi.txt_EGM2D.root","read");
TH2D * electron_hlt_2016_sf = (TH2D*) electron_hlt_2016_sf_file.Get("EGamma_SF2D");

TFile electron_hlt_2017_sf_file("eff_scale_factors/2017/egammaEffi.txt_EGM2D.root","read");
TH2D * electron_hlt_2017_sf = (TH2D*)electron_hlt_2017_sf_file.Get("EGamma_SF2D");

TFile electron_hlt_2018_sf_file("eff_scale_factors/2018/egammaEffi.txt_EGM2D.root" ,"read");
TH2D * electron_hlt_2018_sf = (TH2D*)electron_hlt_2018_sf_file.Get("EGamma_SF2D");

float electron_efficiency_scale_factor(float pt, float eta, string year,bool id_err_up=false, bool reco_err_up=false, bool hlt_err_up=false) {

    TH2F * electron_reco_sf = 0;
    TH2F * electron_id_sf = 0;
    TH2D * electron_hlt_sf = 0;

    if (year == "2016") {
        electron_reco_sf = electron_reco_2016_sf;
        electron_id_sf = electron_id_2016_sf;
        electron_hlt_sf = electron_hlt_2016_sf;
    }
    else if (year == "2017"){
        electron_reco_sf = electron_reco_2017_sf;
        electron_id_sf = electron_id_2017_sf;
        electron_hlt_sf = electron_hlt_2017_sf;
    }
    else if (year == "2018") {
        electron_reco_sf = electron_reco_2018_sf;
        electron_id_sf = electron_id_2018_sf;
        electron_hlt_sf = electron_hlt_2018_sf;
    }
    else
        exit(1);

    int electron_id_sf_xaxisbin = -1;
    int electron_id_sf_yaxisbin = -1;

    if (year == "2016") {    
        electron_id_sf_xaxisbin = electron_id_sf->GetXaxis()->FindFixBin(eta);
        electron_id_sf_yaxisbin = electron_id_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(electron_id_sf->GetYaxis()->GetBinCenter(electron_id_sf->GetNbinsY()))));
    }
    else if (year == "2017") {
        electron_id_sf_xaxisbin = electron_id_sf->GetXaxis()->FindFixBin(eta);
        electron_id_sf_yaxisbin = electron_id_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(electron_id_sf->GetYaxis()->GetBinCenter(electron_id_sf->GetNbinsY()))));
    }
    else if (year == "2018") {
        electron_id_sf_xaxisbin = electron_id_sf->GetXaxis()->FindFixBin(eta);
        electron_id_sf_yaxisbin = electron_id_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(electron_id_sf->GetYaxis()->GetBinCenter(electron_id_sf->GetNbinsY()))));
    }
    else exit(1);

    int electron_hlt_sf_xaxisbin = -1;
    int electron_hlt_sf_yaxisbin = -1;

    if (year == "2016") {    
        electron_hlt_sf_xaxisbin = electron_hlt_sf->GetXaxis()->FindFixBin(eta);
        electron_hlt_sf_yaxisbin = electron_hlt_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(electron_hlt_sf->GetYaxis()->GetBinCenter(electron_hlt_sf->GetNbinsY()))));
    }
    else if (year == "2017") {
        electron_hlt_sf_xaxisbin = electron_hlt_sf->GetXaxis()->FindFixBin(eta);
        electron_hlt_sf_yaxisbin = electron_hlt_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(electron_hlt_sf->GetYaxis()->GetBinCenter(electron_hlt_sf->GetNbinsY()))));
    }
    else if (year == "2018") {
        electron_hlt_sf_xaxisbin = electron_hlt_sf->GetXaxis()->FindFixBin(eta);
        electron_hlt_sf_yaxisbin = electron_hlt_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(electron_hlt_sf->GetYaxis()->GetBinCenter(electron_hlt_sf->GetNbinsY()))));
    }
    else exit(1);


    float sf_id = electron_id_sf->GetBinContent(electron_id_sf_xaxisbin,electron_id_sf_yaxisbin); 
    if (id_err_up) sf_id += electron_id_sf->GetBinError(electron_id_sf_xaxisbin,electron_id_sf_yaxisbin) ;

    //the reco 2D histogram is really a 1D histogram
    float sf_reco=electron_reco_sf->GetBinContent(electron_reco_sf->GetXaxis()->FindFixBin(eta),1);
    if (reco_err_up) sf_reco+=electron_reco_sf->GetBinError(electron_reco_sf->GetXaxis()->FindFixBin(eta),1);

    float sf_hlt = electron_hlt_sf->GetBinContent(electron_hlt_sf_xaxisbin,electron_hlt_sf_yaxisbin); 
    if (hlt_err_up) sf_hlt += electron_hlt_sf->GetBinError(electron_hlt_sf_xaxisbin,electron_hlt_sf_yaxisbin) ;

    return sf_id*sf_reco*sf_hlt;
}


'''

muon_eff_scale_factor_cpp = '''

TFile muon_iso_2016_sf_file("eff_scale_factors/2016/RunBCDEF_SF_ISO.root","read");
TH2D * muon_iso_2016_sf = (TH2D*) muon_iso_2016_sf_file.Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");

TFile muon_id_2016_sf_file("eff_scale_factors/2016/RunBCDEF_SF_ID.root","read");
TH2D * muon_id_2016_sf = (TH2D*) muon_id_2016_sf_file.Get("NUM_TightID_DEN_genTracks_eta_pt");

TFile muon_iso_2017_sf_file("eff_scale_factors/2017/RunBCDEF_SF_ISO.root","read");
TH2D * muon_iso_2017_sf = (TH2D*) muon_iso_2017_sf_file.Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

TFile muon_id_2017_sf_file("eff_scale_factors/2017/RunBCDEF_SF_ID.root","read");
TH2D * muon_id_2017_sf = (TH2D*) muon_id_2017_sf_file.Get("NUM_TightID_DEN_genTracks_pt_abseta");

TFile muon_iso_2018_sf_file("eff_scale_factors/2018/RunABCD_SF_ISO.root","read");
TH2D * muon_iso_2018_sf = (TH2D*) muon_iso_2018_sf_file.Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

TFile muon_id_2018_sf_file("eff_scale_factors/2018/RunABCD_SF_ID.root","read");
TH2D * muon_id_2018_sf = (TH2D*)muon_id_2018_sf_file.Get("NUM_TightID_DEN_TrackerMuons_pt_abseta");

TFile muon_hlt_2016_sf_file("eff_scale_factors/2016/EfficienciesStudies_2016_trigger_EfficienciesAndSF_RunGtoH.root","read");
TH2F * muon_hlt_2016_sf = (TH2F*) muon_hlt_2016_sf_file.Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");

TFile muon_hlt_2017_sf_file("eff_scale_factors/2017/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root","read");
TH2F * muon_hlt_2017_sf = (TH2F*) muon_hlt_2017_sf_file.Get("IsoMu27_PtEtaBins/abseta_pt_ratio");

TFile muon_hlt_2018_sf_file("eff_scale_factors/2018/EfficienciesStudies_2018_trigger_EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root","read");
TH2F * muon_hlt_2018_sf = (TH2F*) muon_hlt_2018_sf_file.Get("IsoMu24_PtEtaBins/abseta_pt_ratio");

float muon_efficiency_scale_factor(float pt,float eta,string year,bool iso_err_up=false,bool id_err_up=false, bool hlt_err_up=false) {

    TH2D * muon_iso_sf = 0;
    TH2D * muon_id_sf = 0;
    TH2F * muon_hlt_sf = 0;

    if (year == "2016") {
        muon_iso_sf = muon_iso_2016_sf;
        muon_id_sf = muon_id_2016_sf;
        muon_hlt_sf = muon_hlt_2016_sf;
    }
    else if (year == "2017") {
        muon_iso_sf = muon_iso_2017_sf;
        muon_id_sf = muon_id_2017_sf;
        muon_hlt_sf = muon_hlt_2017_sf;
    }
    else if (year == "2018"){
        muon_iso_sf = muon_iso_2018_sf;
        muon_id_sf = muon_id_2018_sf;
        muon_hlt_sf = muon_hlt_2018_sf;
    }
    else exit(1);

    int muon_iso_sf_xaxisbin = -1;
    int muon_iso_sf_yaxisbin = -1;

    if (year == "2016") {   
        muon_iso_sf_xaxisbin = muon_iso_sf->GetXaxis()->FindFixBin(eta);
        muon_iso_sf_yaxisbin = muon_iso_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(muon_iso_sf->GetYaxis()->GetBinCenter(muon_iso_sf->GetNbinsY()))));
    }
    else if (year == "2017") {
        muon_iso_sf_yaxisbin = muon_iso_sf->GetYaxis()->FindFixBin(abs(eta));
        muon_iso_sf_xaxisbin = muon_iso_sf->GetXaxis()->FindFixBin(TMath::Min(pt,float(muon_iso_sf->GetXaxis()->GetBinCenter(muon_iso_sf->GetNbinsX()))));
    }
    else if (year == "2018") {
        muon_iso_sf_yaxisbin = muon_iso_sf->GetYaxis()->FindFixBin(abs(eta));
        muon_iso_sf_xaxisbin = muon_iso_sf->GetXaxis()->FindFixBin(TMath::Min(pt,float(muon_iso_sf->GetXaxis()->GetBinCenter(muon_iso_sf->GetNbinsX()))));
    }
    else exit(1);

    int muon_id_sf_xaxisbin = -1;
    int muon_id_sf_yaxisbin = -1;

    if (year == "2016") {    
        muon_id_sf_xaxisbin = muon_id_sf->GetXaxis()->FindFixBin(eta);
        muon_id_sf_yaxisbin = muon_id_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(muon_id_sf->GetYaxis()->GetBinCenter(muon_id_sf->GetNbinsY()))));
    }
    else if (year == "2017") {
        muon_id_sf_yaxisbin = muon_id_sf->GetYaxis()->FindFixBin(abs(eta));
        muon_id_sf_xaxisbin = muon_id_sf->GetXaxis()->FindFixBin(TMath::Min(pt,float(muon_id_sf->GetXaxis()->GetBinCenter(muon_id_sf->GetNbinsX()))));
    }
    else if (year == "2018") {
        muon_id_sf_yaxisbin = muon_id_sf->GetYaxis()->FindFixBin(abs(eta));
        muon_id_sf_xaxisbin = muon_id_sf->GetXaxis()->FindFixBin(TMath::Min(pt,float(muon_id_sf->GetXaxis()->GetBinCenter(muon_id_sf->GetNbinsX()))));
    }
    else exit(1);

    int muon_hlt_sf_xaxisbin = -1;
    int muon_hlt_sf_yaxisbin = -1;

    if (year == "2016") {    
        muon_hlt_sf_xaxisbin = muon_hlt_sf->GetXaxis()->FindFixBin(abs(eta));
        muon_hlt_sf_yaxisbin = muon_hlt_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(muon_hlt_sf->GetYaxis()->GetBinCenter(muon_hlt_sf->GetNbinsY()))));
    }
    else if (year == "2017") {
        muon_hlt_sf_xaxisbin = muon_hlt_sf->GetXaxis()->FindFixBin(abs(eta));
        muon_hlt_sf_yaxisbin = muon_hlt_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(muon_hlt_sf->GetYaxis()->GetBinCenter(muon_hlt_sf->GetNbinsY()))));
    }
    else if (year == "2018") {
        muon_hlt_sf_xaxisbin = muon_hlt_sf->GetXaxis()->FindFixBin(abs(eta));
        muon_hlt_sf_yaxisbin = muon_hlt_sf->GetYaxis()->FindFixBin(TMath::Min(pt,float(muon_hlt_sf->GetYaxis()->GetBinCenter(muon_hlt_sf->GetNbinsY()))));
    }
    else exit(1);

    float iso_sf = muon_iso_sf->GetBinContent(muon_iso_sf_xaxisbin,muon_iso_sf_yaxisbin);

    if (iso_err_up) iso_sf += muon_iso_sf->GetBinError(muon_iso_sf_xaxisbin,muon_iso_sf_yaxisbin);

    float id_sf = muon_id_sf->GetBinContent(muon_id_sf_xaxisbin,muon_id_sf_yaxisbin); 
    
    if (id_err_up) id_sf += muon_id_sf->GetBinError(muon_id_sf_xaxisbin,muon_id_sf_yaxisbin) ;

    float hlt_sf = muon_hlt_sf->GetBinContent(muon_hlt_sf_xaxisbin,muon_hlt_sf_yaxisbin); 
    
    if (hlt_err_up) hlt_sf += muon_hlt_sf->GetBinError(muon_hlt_sf_xaxisbin,muon_hlt_sf_yaxisbin) ;

    return iso_sf * id_sf * hlt_sf;

}

'''

fake_lepton_weight_cpp = '''

TFile muon_2016_file("fake_lepton_weights/muon_2016_frs.root");
TFile electron_2016_file("fake_lepton_weights/electron_2016_frs.root");

TFile muon_2017_file("fake_lepton_weights/muon_2017_frs.root");
TFile electron_2017_file("fake_lepton_weights/electron_2017_frs.root");

TFile muon_2018_file("fake_lepton_weights/muon_2018_frs.root");
TFile electron_2018_file("fake_lepton_weights/electron_2018_frs.root");

TH2D * muon_2016_fr_hist = (TH2D*)muon_2016_file.Get("muon_frs");
TH2D * electron_2016_fr_hist = (TH2D*)electron_2016_file.Get("electron_frs");
TH2D * muon_2017_fr_hist = (TH2D*)muon_2017_file.Get("muon_frs");
TH2D * electron_2017_fr_hist = (TH2D*)electron_2017_file.Get("electron_frs");
TH2D * muon_2018_fr_hist = (TH2D*)muon_2018_file.Get("muon_frs");
TH2D * electron_2018_fr_hist = (TH2D*)electron_2018_file.Get("electron_frs");

float get_fake_lepton_weight(float eta, float pt, string year, int lepton_pdg_id, string syst = "nominal")
{
    TH2D * fr_hist = 0;

    if (year == "2016" && abs(lepton_pdg_id) == 13) fr_hist = muon_2016_fr_hist;
    else if (year == "2016" && abs(lepton_pdg_id) == 11) fr_hist = electron_2016_fr_hist;
    else if (year == "2017" && abs(lepton_pdg_id) == 13) fr_hist = muon_2017_fr_hist;
    else if (year == "2017" && abs(lepton_pdg_id) == 11) fr_hist = electron_2017_fr_hist;
    else if (year == "2018" && abs(lepton_pdg_id) == 13) fr_hist = muon_2018_fr_hist;
    else if (year == "2018" && abs(lepton_pdg_id) == 11) fr_hist = electron_2018_fr_hist;
    else exit(1);

    float myeta  = TMath::Min(abs(eta),float(2.4999));
    float mypt  = TMath::Min(pt,float(44.999));

    int etabin = fr_hist->GetXaxis()->FindFixBin(myeta);
    int ptbin = fr_hist->GetYaxis()->FindFixBin(mypt);

    float prob = fr_hist->GetBinContent(etabin,ptbin);

    if (syst == "up") prob += fr_hist->GetBinError(etabin,ptbin);
    else if (syst == "down") prob -= fr_hist->GetBinError(etabin,ptbin);
    else assert(syst == "nominal");

    return prob/(1-prob);
}
'''

fake_photon_weight_cpp = '''

float get_fake_photon_weight(float eta, float pt, string year, int lepton_pdg_id, string version = "nominal") {

if (version == "wjets") { 
    assert(year == "2016");
    float fr = 0;
    if (year == "2016") {
        if (abs(eta) < 1.4442) {
            if (pt < 30 && pt > 20) fr = 0.863562091503268;
            else if (pt < 40 && pt > 30) fr = 0.8110749185667753;
            else if (pt < 50 && pt > 40) fr = 0.8356164383561644;
            else if (pt > 50) fr = 0.6576763485477178;
            else exit(1); 
        }
        else if (1.566 < abs(eta) && abs(eta) < 2.5) {
           if (pt < 30 && pt > 20) fr = 0.7814113597246127;
               else if (pt < 40 && pt > 30) fr = 1.2465373961218837;
               else if (pt < 50 && pt > 40) fr = 1.1513157894736843;
               else if (pt > 50) fr = 1.3089430894308942;
               else exit(1); 
           }
        }
        return fr;
    }
else if (version == "wjets_chiso") {
    assert(year == "2016");
    float fr = 0;
    if (year == "2016") {
        if (abs(eta) < 1.4442) {
            if (pt < 30 && pt > 20) fr = 2.0475529330093716;
            else if (pt < 40 && pt > 30) fr = 1.4072790294627382;
            else if (pt < 50 && pt > 40) fr = 1.053030303030303;
            else if (pt > 50) fr = 0.8095238095238096;
            else exit(1); 
    
           }
           else if (1.566 < abs(eta) && abs(eta) < 2.5) {
               if (pt < 30 && pt > 20) fr = 3.556842105263158;
               else if (pt < 40 && pt > 30) fr = 2.292079207920792;
               else if (pt < 50 && pt > 40) fr = 2.1470588235294117;
               else if (pt > 50) fr = 2.5769230769230766;
               else exit(1); 
           }
        }
        return fr;
    }
else if (version == "nominal" || version.substr(0,3) == "alt" || version == "stat_up") { //based on inverting chiso and making the maximum 1.75*chiso_cut
    float fr = 0;
    if (year == "2016") {
       if (abs(eta) < 1.4442) {
          if (pt < 30 && pt > 20) fr = 0.7383849589038014;
          else if (pt < 40 && pt > 30) fr = 0.711832976726155;
          else if (pt < 50 && pt > 40) fr = 0.6096371946961419;
          else if (pt > 50) fr = 0.4248350921227493;
          else exit(1); 
    
       }
       else if (1.566 < abs(eta) && abs(eta) < 2.5) {
          if (pt < 30 && pt > 20) fr = 0.9641149730355477;
          else if (pt < 40 && pt > 30) fr = 0.9881713412122083;
          else if (pt < 50 && pt > 40) fr = 0.9490490494146919;
          else if (pt > 50) fr = 0.9670757894326804;
          else exit(1); 
       }
    }
    else if (year == "2017") {
       if (abs(eta) < 1.4442) {
          if (pt < 30 && pt > 20) fr = 0.74097488174064;
          else if (pt < 40 && pt > 30) fr = 0.7575997778550154;
          else if (pt < 50 && pt > 40) fr = 0.6795843358575392;
          else if (pt > 50) fr = 0.5202412984429878;
          else exit(1); 
    
       }
       else if (1.566 < abs(eta) && abs(eta) < 2.5) {
          if (pt < 30 && pt > 20) fr = 0.3500132190964065;
          else if (pt < 40 && pt > 30) fr = 0.3708241722875352;
          else if (pt < 50 && pt > 40) fr = 0.426602270486085;
          else if (pt > 50) fr = 0.5421704176994768;
          else exit(1); 
       }
    }
    else if (year == "2018") {
       if (abs(eta) < 1.4442) {
          if (pt < 30 && pt > 20) fr = 0.7403414019248873;
          else if (pt < 40 && pt > 30) fr = 0.7790466202005484;
          else if (pt < 50 && pt > 40) fr = 0.71722564116773;
          else if (pt > 50) fr = 0.515278581311353;
          else exit(1); 
    
       }
       else if (1.566 < abs(eta) && abs(eta) < 2.5) {
          if (pt < 30 && pt > 20) fr = 0.3086299852372086;
          else if (pt < 40 && pt > 30) fr = 0.34880515284451147;
          else if (pt < 50 && pt > 40) fr = 0.37604537406362587;
          else if (pt > 50) fr = 0.5203825122902803;
          else exit(1); 
       }
    }

    if (version.substr(0,3) == "alt") {
       if (year == "2016") {
          if (abs(eta) < 1.4442) {
             if (pt < 30 && pt > 20 && version == "alt0") fr += 0.863562091503268 - 0.733248094217973;
             else if (pt < 40 && pt > 30 && version == "alt1") fr += 0.8110749185667753 - 0.7931335006595595;
             else if (pt < 50 && pt > 40 && version == "alt2") fr += 0.8356164383561644 - 0.682480210665109;
             else if (pt > 50 && version == "alt3") fr += 0.6576763485477178 - 0.4236927109550702;
          }
          else if (1.566 < abs(eta) && abs(eta) < 2.5) {
             if (pt < 30 && pt > 20 && version == "alt4") fr += 0.7814113597246127 - 1.0000609792613757;
             else if (pt < 40 && pt > 30 && version == "alt5") fr += 1.2465373961218837 - 0.8299609184544505;
             else if (pt < 50 && pt > 40 && version == "alt6") fr += 1.1513157894736843 - 1.2664123432054573;
             else if (pt > 50 && version == "alt7") fr += 1.3089430894308942 - 0.9603173729944264;
          }
       }
       else if (year == "2017") {
          if (abs(eta) < 1.4442) {
             if (pt < 30 && pt > 20 && version == "alt8") fr += 0.7338530066815144-0.7450558550432999;
             else if (pt < 40 && pt > 30 && version == "alt9") fr += 0.7475409836065574-0.8841120276900355; 
             else if (pt < 50 && pt > 40 && version == "alt10") fr += 1.0245098039215685-0.7939865634894522;
             else if (pt > 50 && version == "alt11") fr += 0.7470817120622568-0.6225170980914996;
          }
          else if (1.566 < abs(eta) && abs(eta) < 2.5) {
             if (pt < 30 && pt > 20 && version == "alt12") fr += 0.6666666666666666-0.3979138739228848;
             else if (pt < 40 && pt > 30 && version == "alt13") fr += 0.7617021276595745-0.5119586040478171;
             else if (pt < 50 && pt > 40 && version == "alt14") fr += 0.8472222222222222-0.5687851363533186;
             else if (pt > 50 && version == "alt15") fr += 1.0857142857142856-0.8108939928786887;
          }
       }
       else if (year == "2018") {
          if (abs(eta) < 1.4442) {
             if (pt < 30 && pt > 20 && version == "alt16") fr += 0.8406708595387841-0.8303084927443096;
             else if (pt < 40 && pt > 30 && version == "alt17") fr += 0.948051948051948-0.9197645786035942;
             else if (pt < 50 && pt > 40 && version == "alt18") fr += 0.7518796992481203-0.7737938893939876;
             else if (pt > 50 && version == "alt19") fr += 0.7315436241610739-0.47164126712162807;
          }
          else if (1.566 < abs(eta) && abs(eta) < 2.5) {
             if (pt < 30 && pt > 20 && version == "alt20") fr += 0.6810344827586207-0.3930515103150424;
             else if (pt < 40 && pt > 30 && version == "alt21") fr += 0.759493670886076-0.47176049745670173;
             else if (pt < 50 && pt > 40 && version == "alt22") fr += 0.6212121212121212-0.43181336073191234;
             else if (pt > 50 && version == "alt23") fr += 0.8333333333333334-0.5139996676356342;
          }
       }
    }

    if (version == "stat_up") {
       if (year == "2016") {
          if (abs(eta) < 1.4442) {
             if (pt < 30 && pt > 20) fr += 0.0026960729583444783;
             else if (pt < 40 && pt > 30) fr += 0.006662460844367098;
             else if (pt < 50 && pt > 40) fr += 0.009055572023850591;
             else if (pt > 50) fr += 0.005570192293811733;
             else exit(1); 
    
          }
          else if (1.566 < abs(eta) && abs(eta) < 2.5) {
             if (pt < 30 && pt > 20) fr += 0.005818625145987401;
             else if (pt < 40 && pt > 30) fr += 0.014817966152871172;
             else if (pt < 50 && pt > 40) fr += 0.024028693336438096;
             else if (pt > 50) fr += 0.02609849029933926;
             else exit(1); 
          }
       }
       else if (year == "2017") {
          if (abs(eta) < 1.4442) {
             if (pt < 30 && pt > 20) fr += 0.0026718396811323355;
             else if (pt < 40 && pt > 30) fr += 0.007173244839053661;
             else if (pt < 50 && pt > 40) fr += 0.009763056116442556;
             else if (pt > 50) fr += 0.006612945121205288;
             else exit(1); 
          }
          else if (1.566 < abs(eta) && abs(eta) < 2.5) {
             if (pt < 30 && pt > 20) fr += 0.002657455941380463;
             else if (pt < 40 && pt > 30) fr += 0.006949455204473568;
             else if (pt < 50 && pt > 40) fr += 0.013323254649466083;
             else if (pt > 50) fr += 0.017033894200413374;
             else exit(1); 
          }
       }
       else if (year == "2018") {
          if (abs(eta) < 1.4442) {
             if (pt < 30 && pt > 20) fr += 0.002262126317855081;
             else if (pt < 40 && pt > 30) fr += 0.006082165017453894;
             else if (pt < 50 && pt > 40) fr += 0.00848747881815101;
             else if (pt > 50) fr += 0.005792952999897535;
             else exit(1); 
          }
          else if (1.566 < abs(eta) && abs(eta) < 2.5) {
             if (pt < 30 && pt > 20) fr += 0.0018299453388268458;
             else if (pt < 40 && pt > 30) fr += 0.005018789997370193;
             else if (pt < 50 && pt > 40) fr += 0.008666577102573342;
             else if (pt > 50) fr += 0.011947487425581427;
             else exit(1); 
          }
       }
    }


    return fr;
} else {

exit(1);

}

return 0;  

}
'''

ROOT.gInterpreter.Declare(include_headers_cpp)
ROOT.gInterpreter.ProcessLine(wjets_pileup_photons_flags_initialization_cpp)
ROOT.gInterpreter.Declare(wjets_pileup_photons_flags_cpp)
ROOT.gInterpreter.Declare(fake_lepton_weight_cpp)
ROOT.gInterpreter.Declare(fake_photon_weight_cpp)
ROOT.gInterpreter.Declare(muon_eff_scale_factor_cpp)
ROOT.gInterpreter.Declare(electron_eff_scale_factor_cpp)
ROOT.gInterpreter.Declare(photon_eff_scale_factor_cpp)

def processMCSample(dummy):

    photon_gen_matching_for_fake_cutstring = "("
    photon_gen_matching_cutstring = "("

    if sample["fsr"]:
        photon_gen_matching_for_fake_cutstring+="photon_gen_matching == 4"
        photon_gen_matching_cutstring+="photon_gen_matching == 4"
    if sample["non_fsr"]:  
        if photon_gen_matching_for_fake_cutstring != "(":
            photon_gen_matching_for_fake_cutstring += " || "
        if photon_gen_matching_cutstring != "(":
            photon_gen_matching_cutstring += " || "
        photon_gen_matching_for_fake_cutstring+="photon_gen_matching == 5 || photon_gen_matching == 6"
        photon_gen_matching_cutstring+="photon_gen_matching == 5 || photon_gen_matching == 6"
    if sample["e_to_p_for_fake"]:
        if photon_gen_matching_for_fake_cutstring != "(":
            photon_gen_matching_for_fake_cutstring += " || "
        photon_gen_matching_for_fake_cutstring+="photon_gen_matching == 1"
    if sample["non-prompt"]:
        pass
        if photon_gen_matching_cutstring != "(":
            photon_gen_matching_cutstring += " || "
        #if we used (!is_photon_prompt && !photon_genjet_matching) instead of (!photon_genjet_matching && !is_photon_prompt), then we would call is_photon_prompt for some events that are not in the std::map
        photon_gen_matching_cutstring+="(!(photon_gen_matching == 1 || photon_gen_matching == 4 || photon_gen_matching == 5 || photon_gen_matching == 6) && photon_genjet_matching) || (!photon_genjet_matching && !is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]))" 
                
    if photon_gen_matching_for_fake_cutstring != "(":    
        photon_gen_matching_for_fake_cutstring+= ")"    
    else:
        photon_gen_matching_for_fake_cutstring= "0"    

    if photon_gen_matching_cutstring != "(":    
        photon_gen_matching_cutstring+= ")"    
    else:
        photon_gen_matching_cutstring= "0"    

    rdf = ROOT.RDataFrame("Events",sample["filename"])

    #the JERUp and JESUp information was not added to the w+jets sample
    if  label == "w+jets" or label == "gg+jets":
        rinterface = rdf.Filter(get_filter_string(year,isdata=True))
    else:    
        rinterface = rdf.Filter(get_filter_string(year,isdata=False))
    if year == "2016" and label == "w+jets":
        rinterface_wjets_2016 = rdf.Filter(get_filter_string(year,isdata=True,lep="both"))

            
    rinterface = rinterface.Define("xs_weight",str(sample["xs"]*1000*lumi/sample["nweightedevents"]) + "*gen_weight/abs(gen_weight)") 

    if year == "2016" or year == "2017":    
        prefire_weight_string = "PrefireWeight"
        prefire_up_weight_string = "PrefireWeight_Up"
    else:    
        prefire_weight_string = "1"
        prefire_up_weight_string = "1"

    rinterface = rinterface.Define("base_weight",get_postfilter_selection_string()+"*xs_weight*puWeight*"+prefire_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\") : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\"))")      
    rinterface = rinterface.Define("prefire_up_base_weight",get_postfilter_selection_string()+"*xs_weight*puWeight*"+prefire_up_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\") : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\"))")    
    rinterface = rinterface.Define("pileup_up_base_weight",get_postfilter_selection_string()+"*xs_weight*puWeightUp*"+prefire_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\") : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\"))")    
    rinterface = rinterface.Define("electron_id_sf_up_base_weight",get_postfilter_selection_string()+"*xs_weight*puWeight*"+prefire_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\") : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\",true))")                  
    rinterface = rinterface.Define("electron_reco_sf_up_base_weight",get_postfilter_selection_string()+"*xs_weight*puWeight*"+prefire_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\") : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\",false,true))")    
    rinterface = rinterface.Define("electron_hlt_sf_up_base_weight",get_postfilter_selection_string()+"*xs_weight*puWeight*"+prefire_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\") : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\",false,false,true))")    
    rinterface = rinterface.Define("muon_id_sf_up_base_weight",get_postfilter_selection_string()+"*xs_weight*puWeight*"+prefire_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\",false,true) : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\"))")                  
    rinterface = rinterface.Define("muon_iso_sf_up_base_weight",get_postfilter_selection_string()+"*xs_weight*puWeight*"+prefire_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\",true) : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\"))")    
    rinterface = rinterface.Define("muon_hlt_sf_up_base_weight",get_postfilter_selection_string()+"*xs_weight*puWeight*"+prefire_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\",false,false,true) : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\"))")    
    rinterface = rinterface.Define("photon_id_sf_up_base_weight",get_postfilter_selection_string()+"*xs_weight*puWeight*"+prefire_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\",true)*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\") : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\"))") 
    if label != "w+jets" and label != "gg+jets":
        rinterface = rinterface.Define("jes_up_base_weight",get_postfilter_selection_string("JESUp")+"*xs_weight*puWeight*"+prefire_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\") : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\"))") 
        rinterface = rinterface.Define("jer_up_base_weight",get_postfilter_selection_string("JERUp")+"*xs_weight*puWeight*"+prefire_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\") : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\"))") 

    if label == "gg+jets":    
        #conversion background
        rinterface = rinterface.Define("pass_lepton_gen","abs(lepton_pdg_id) == 11 &&  lepton_gen_matching == 22 ? 1 : 0")
    else:
        if "is_lepton_real" in rinterface.GetColumnNames():
            rinterface = rinterface.Define("pass_lepton_gen","is_lepton_real == 1")
        else:    
            rinterface = rinterface.Define("pass_lepton_gen","lepton_gen_matching == 1 || lepton_gen_matching == 15")

    rinterface = rinterface.Define("weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*base_weight")

    if label == "w+jets":
        rinterface = rinterface.Define("prompt_pileup_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]))*base_weight")
        rinterface = rinterface.Define("prompt_pileup_pileup_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]))*pileup_up_base_weight")
        rinterface = rinterface.Define("prompt_pileup_prefire_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]))*prefire_up_base_weight")
        rinterface = rinterface.Define("prompt_pileup_photon_id_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]))*photon_id_sf_up_base_weight")
        rinterface = rinterface.Define("prompt_pileup_electron_reco_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]))*electron_reco_sf_up_base_weight")
        rinterface = rinterface.Define("prompt_pileup_electron_id_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]))*electron_id_sf_up_base_weight")
        rinterface = rinterface.Define("prompt_pileup_electron_hlt_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]))*electron_hlt_sf_up_base_weight")
        rinterface = rinterface.Define("prompt_pileup_muon_id_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]))*muon_id_sf_up_base_weight")
        rinterface = rinterface.Define("prompt_pileup_muon_hlt_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]))*muon_hlt_sf_up_base_weight")
        rinterface = rinterface.Define("prompt_pileup_muon_iso_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]))*muon_iso_sf_up_base_weight")

    if label == "wg+jets":
        rinterface = rinterface.Define("fid","pass_fid_selection && fid_met_pt > 0")
        rinterface = rinterface.Define("pass_fiducial_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*base_weight")
        rinterface = rinterface.Define("fail_fiducial_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*base_weight")
        rinterface = rinterface.Define("pass_fiducial_pileup_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*pileup_up_base_weight")
        rinterface = rinterface.Define("fail_fiducial_pileup_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*pileup_up_base_weight")
        rinterface = rinterface.Define("pass_fiducial_prefire_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*prefire_up_base_weight")
        rinterface = rinterface.Define("fail_fiducial_prefire_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*prefire_up_base_weight")
        rinterface = rinterface.Define("pass_fiducial_jes_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*jes_up_base_weight")
        rinterface = rinterface.Define("fail_fiducial_jes_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*jes_up_base_weight")
        rinterface = rinterface.Define("pass_fiducial_jer_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*jer_up_base_weight")
        rinterface = rinterface.Define("fail_fiducial_jer_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*jer_up_base_weight")
        rinterface = rinterface.Define("pass_fiducial_photon_id_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*photon_id_sf_up_base_weight")
        rinterface = rinterface.Define("fail_fiducial_photon_id_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*photon_id_sf_up_base_weight")
        rinterface = rinterface.Define("pass_fiducial_electron_reco_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*electron_reco_sf_up_base_weight")
        rinterface = rinterface.Define("fail_fiducial_electron_reco_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*electron_reco_sf_up_base_weight")
        rinterface = rinterface.Define("pass_fiducial_electron_id_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*electron_id_sf_up_base_weight")
        rinterface = rinterface.Define("fail_fiducial_electron_id_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*electron_id_sf_up_base_weight")
        rinterface = rinterface.Define("pass_fiducial_electron_hlt_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*electron_hlt_sf_up_base_weight")
        rinterface = rinterface.Define("fail_fiducial_electron_hlt_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*electron_hlt_sf_up_base_weight")
        rinterface = rinterface.Define("pass_fiducial_muon_id_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*muon_id_sf_up_base_weight")
        rinterface = rinterface.Define("fail_fiducial_muon_id_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*muon_id_sf_up_base_weight")
        rinterface = rinterface.Define("pass_fiducial_muon_iso_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*muon_iso_sf_up_base_weight")
        rinterface = rinterface.Define("fail_fiducial_muon_iso_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*muon_iso_sf_up_base_weight")
        rinterface = rinterface.Define("pass_fiducial_muon_hlt_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*muon_hlt_sf_up_base_weight")
        rinterface = rinterface.Define("fail_fiducial_muon_hlt_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*muon_hlt_sf_up_base_weight")

        rinterface = rinterface.Define("fake_photon_pass_fiducial_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*base_weight : 0")
        rinterface = rinterface.Define("fake_photon_fail_fiducial_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*base_weight : 0")
        for i in range(n_fake_photon_alt):
            rinterface = rinterface.Define("fake_photon_pass_fiducial_fake_photon_alt"+str(i)+"_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"alt"+str(i)+"\")*base_weight : 0")
            rinterface = rinterface.Define("fake_photon_fail_fiducial_fake_photon_alt"+str(i)+"_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"alt"+str(i)+"\")*base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pass_fiducial_pileup_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*pileup_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_fail_fiducial_pileup_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*pileup_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pass_fiducial_prefire_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*prefire_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_fail_fiducial_prefire_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*prefire_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pass_fiducial_jes_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*jes_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_fail_fiducial_jes_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*jes_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pass_fiducial_jer_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*jer_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_fail_fiducial_jer_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*jer_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pass_fiducial_photon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*photon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_fail_fiducial_photon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*photon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pass_fiducial_electron_reco_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_reco_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_fail_fiducial_electron_reco_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_reco_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pass_fiducial_electron_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_fail_fiducial_electron_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pass_fiducial_electron_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_hlt_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_fail_fiducial_electron_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_hlt_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pass_fiducial_muon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_fail_fiducial_muon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pass_fiducial_muon_iso_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_iso_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_fail_fiducial_muon_iso_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_iso_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pass_fiducial_muon_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_hlt_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_fail_fiducial_muon_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_hlt_sf_up_base_weight : 0")

        rinterface = rinterface.Define("fake_lepton_pass_fiducial_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_fail_fiducial_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_pass_fiducial_pileup_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*pileup_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_fail_fiducial_pileup_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*pileup_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_pass_fiducial_prefire_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*prefire_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_fail_fiducial_prefire_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*prefire_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_pass_fiducial_jes_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*jes_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_fail_fiducial_jes_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*jes_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_pass_fiducial_jer_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*jer_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_fail_fiducial_jer_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*jer_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_pass_fiducial_photon_id_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*photon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_fail_fiducial_photon_id_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*photon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_pass_fiducial_electron_reco_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*electron_reco_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_fail_fiducial_electron_reco_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*electron_reco_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_pass_fiducial_electron_id_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*electron_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_fail_fiducial_electron_id_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*electron_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_pass_fiducial_electron_hlt_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*electron_hlt_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_fail_fiducial_electron_hlt_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*electron_hlt_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_pass_fiducial_muon_id_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*muon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_fail_fiducial_muon_id_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*muon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_pass_fiducial_muon_iso_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*muon_iso_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_fail_fiducial_muon_iso_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*muon_iso_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_pass_fiducial_muon_hlt_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*muon_hlt_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_fail_fiducial_muon_hlt_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*muon_hlt_sf_up_base_weight : 0")

        rinterface = rinterface.Define("double_fake_pass_fiducial_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*base_weight : 0")
        rinterface = rinterface.Define("double_fake_fail_fiducial_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*base_weight : 0") 
        for i in range(n_fake_photon_alt):
            rinterface = rinterface.Define("double_fake_pass_fiducial_fake_photon_alt"+str(i)+"_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"alt"+str(i)+"\")*base_weight : 0")
            rinterface = rinterface.Define("double_fake_fail_fiducial_fake_photon_alt"+str(i)+"_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"alt"+str(i)+"\")*base_weight : 0") 
        rinterface = rinterface.Define("double_fake_pass_fiducial_pileup_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*pileup_up_base_weight : 0")
        rinterface = rinterface.Define("double_fake_fail_fiducial_pileup_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*pileup_up_base_weight : 0") 
        rinterface = rinterface.Define("double_fake_pass_fiducial_prefire_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*prefire_up_base_weight : 0")
        rinterface = rinterface.Define("double_fake_fail_fiducial_prefire_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*prefire_up_base_weight : 0") 
        rinterface = rinterface.Define("double_fake_pass_fiducial_jes_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*jes_up_base_weight : 0")
        rinterface = rinterface.Define("double_fake_fail_fiducial_jes_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*jes_up_base_weight : 0") 
        rinterface = rinterface.Define("double_fake_pass_fiducial_jer_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*jer_up_base_weight : 0")
        rinterface = rinterface.Define("double_fake_fail_fiducial_jer_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*jer_up_base_weight : 0") 
        rinterface = rinterface.Define("double_fake_pass_fiducial_photon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*photon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("double_fake_fail_fiducial_photon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*photon_id_sf_up_base_weight : 0") 
        rinterface = rinterface.Define("double_fake_pass_fiducial_electron_reco_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_reco_sf_up_base_weight : 0")
        rinterface = rinterface.Define("double_fake_fail_fiducial_electron_reco_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_reco_sf_up_base_weight : 0") 
        rinterface = rinterface.Define("double_fake_pass_fiducial_electron_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("double_fake_fail_fiducial_electron_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_id_sf_up_base_weight : 0") 
        rinterface = rinterface.Define("double_fake_pass_fiducial_electron_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_hlt_sf_up_base_weight : 0")
        rinterface = rinterface.Define("double_fake_fail_fiducial_electron_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_hlt_sf_up_base_weight : 0") 
        rinterface = rinterface.Define("double_fake_pass_fiducial_muon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("double_fake_fail_fiducial_muon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_id_sf_up_base_weight : 0") 
        rinterface = rinterface.Define("double_fake_pass_fiducial_muon_iso_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_iso_sf_up_base_weight : 0")
        rinterface = rinterface.Define("double_fake_fail_fiducial_muon_iso_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_iso_sf_up_base_weight : 0") 
        rinterface = rinterface.Define("double_fake_pass_fiducial_muon_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_hlt_sf_up_base_weight : 0")
        rinterface = rinterface.Define("double_fake_fail_fiducial_muon_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" && !fid ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_hlt_sf_up_base_weight : 0") 

        if labels["wg+jets"]["syst-scale"]:
            for i in range(0,8):
                        #this sample has a bug that causes the scale weight to be 1/2 the correct value
                if sample["filename"] == args.workdir+"/data/wg/2016/1June2019/wgjets.root" or sample["filename"] == args.workdir+"/data/wg/2016/1June2019jetunc/wgjets.root":
                    rinterface = rinterface.Define("pass_fiducial_scale"+str(i)+"_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*base_weight*LHEScaleWeight["+str(i)+"]*2")
                    rinterface = rinterface.Define("fail_fiducial_scale"+str(i)+"_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*base_weight*LHEScaleWeight["+str(i)+"]*2")
                else:    
                    rinterface = rinterface.Define("pass_fiducial_scale"+str(i)+"_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*base_weight*LHEScaleWeight["+str(i)+"]")
                    rinterface = rinterface.Define("fail_fiducial_scale"+str(i)+"_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*base_weight*LHEScaleWeight["+str(i)+"]")

        if labels["wg+jets"]["syst-pdf"]:
            for i in range(0,32):
                #this sample has a bug that causes the scale weight to be 1/2 the correct value
                if sample["filename"] == args.workdir+"/data/wg/2016/1June2019/wgjets.root" or sample["filename"] == args.workdir+"/data/wg/2016/1June2019jetunc/wgjets.root":
                    rinterface = rinterface.Define("pass_fiducial_pdf"+str(i)+"_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*base_weight*LHEPdfWeight["+str(i+1)+"]*2")
                    rinterface = rinterface.Define("fail_fiducial_pdf"+str(i)+"_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*base_weight*LHEPdfWeight["+str(i+1)+"]*2")
                else:    
                    rinterface = rinterface.Define("pass_fiducial_pdf"+str(i)+"_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && fid)*base_weight*LHEPdfWeight["+str(i+1)+"]")
                    rinterface = rinterface.Define("fail_fiducial_pdf"+str(i)+"_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && !fid)*base_weight*LHEPdfWeight["+str(i+1)+"]")

    rinterface = rinterface.Define("pileup_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*pileup_up_base_weight")
    rinterface = rinterface.Define("prefire_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*prefire_up_base_weight")
    rinterface = rinterface.Define("electron_id_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*electron_id_sf_up_base_weight")
    rinterface = rinterface.Define("electron_reco_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*electron_reco_sf_up_base_weight")
    rinterface = rinterface.Define("electron_hlt_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*electron_hlt_sf_up_base_weight")
    rinterface = rinterface.Define("muon_id_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*muon_id_sf_up_base_weight")
    rinterface = rinterface.Define("muon_iso_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*muon_iso_sf_up_base_weight")
    rinterface = rinterface.Define("muon_hlt_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*muon_hlt_sf_up_base_weight")
    rinterface = rinterface.Define("photon_id_sf_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*photon_id_sf_up_base_weight")

    if label != "w+jets" and label != "gg+jets":
        rinterface = rinterface.Define("jes_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*jes_up_base_weight")
        rinterface = rinterface.Define("jer_up_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*jer_up_base_weight")

    if labels[label]["syst-scale"]:
        for i in range(0,8):
                     #this sample has a bug that causes the scale weight to be 1/2 the correct value
            if sample["filename"] == args.workdir+"/data/wg/2016/1June2019/wgjets.root" or sample["filename"] == args.workdir+"/data/wg/2016/1June2019jetunc/wgjets.root":
                rinterface = rinterface.Define("scale"+str(i)+"_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*base_weight*LHEScaleWeight["+str(i)+"]*2")
            else:    
                rinterface = rinterface.Define("scale"+str(i)+"_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*base_weight*LHEScaleWeight["+str(i)+"]")

    if labels[label]["syst-pdf"]:
        for i in range(0,32):
            if sample["filename"] == args.workdir+"/data/wg/2016/1June2019/wgjets.root" or sample["filename"] == args.workdir+"/data/wg/2016/1June2019jetunc/wgjets.root":
                rinterface = rinterface.Define("pdf"+str(i)+"_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*base_weight*LHEPdfWeight["+str(i+1)+"]*2")
            else:    
                rinterface = rinterface.Define("pdf"+str(i)+"_weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + ")*base_weight*LHEPdfWeight["+str(i+1)+"]")

#            rinterface = rinterface.Define("fake_lepton_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*base_weight : 0")
#            rinterface = rinterface.Define("fake_photon_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*base_weight : 0")
#            rinterface = rinterface.Define("double_fake_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*base_weight : 0") 

    rinterface = rinterface.Define("fake_lepton_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*xs_weight*puWeight*"+prefire_weight_string+"*" + get_postfilter_selection_string()+" : 0")
    rinterface = rinterface.Define("fake_lepton_pileup_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*pileup_up_base_weight : 0")
    rinterface = rinterface.Define("fake_lepton_prefire_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*prefire_up_base_weight : 0")
    if label != "w+jets" and label != "gg+jets":
        rinterface = rinterface.Define("fake_lepton_jes_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*jes_up_base_weight : 0")
        rinterface = rinterface.Define("fake_lepton_jer_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*jer_up_base_weight : 0")
    rinterface = rinterface.Define("fake_lepton_photon_id_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*photon_id_sf_up_base_weight : 0")
    rinterface = rinterface.Define("fake_lepton_electron_reco_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*electron_reco_sf_up_base_weight : 0")
    rinterface = rinterface.Define("fake_lepton_electron_id_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*electron_id_sf_up_base_weight : 0")
    rinterface = rinterface.Define("fake_lepton_electron_hlt_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*electron_hlt_sf_up_base_weight : 0")
    rinterface = rinterface.Define("fake_lepton_muon_id_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*muon_id_sf_up_base_weight : 0")
    rinterface = rinterface.Define("fake_lepton_muon_iso_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*muon_iso_sf_up_base_weight : 0")
    rinterface = rinterface.Define("fake_lepton_muon_hlt_sf_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*muon_hlt_sf_up_base_weight : 0")
    rinterface = rinterface.Define("fake_lepton_stat_up_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id,\"up\")*xs_weight*puWeight*"+prefire_weight_string + "*" + get_postfilter_selection_string()+" : 0")
    rinterface = rinterface.Define("fake_lepton_stat_down_weight","photon_selection == 0 && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id,\"down\")*xs_weight*puWeight*"+prefire_weight_string+"*" + get_postfilter_selection_string()+" : 0")

    rinterface = rinterface.Define("double_fake_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*xs_weight*puWeight*"+prefire_weight_string+"*" + get_postfilter_selection_string()+" : 0") 
    rinterface = rinterface.Define("double_fake_pileup_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*pileup_up_base_weight : 0")
    rinterface = rinterface.Define("double_fake_prefire_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*prefire_up_base_weight : 0") 
    if label != "w+jets" and label != "gg+jets":
        rinterface = rinterface.Define("double_fake_jer_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*jer_up_base_weight : 0") 
        rinterface = rinterface.Define("double_fake_jes_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*jes_up_base_weight : 0") 
    rinterface = rinterface.Define("double_fake_photon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*photon_id_sf_up_base_weight : 0") 
    rinterface = rinterface.Define("double_fake_electron_reco_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_reco_sf_up_base_weight : 0") 
    rinterface = rinterface.Define("double_fake_electron_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_id_sf_up_base_weight : 0") 
    rinterface = rinterface.Define("double_fake_electron_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_hlt_sf_up_base_weight : 0") 
    rinterface = rinterface.Define("double_fake_muon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_id_sf_up_base_weight : 0") 
    rinterface = rinterface.Define("double_fake_muon_iso_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_iso_sf_up_base_weight : 0") 
    rinterface = rinterface.Define("double_fake_muon_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_hlt_sf_up_base_weight : 0") 
    for i in range(n_fake_photon_alt):
        rinterface = rinterface.Define("double_fake_alt"+str(i)+"_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"alt"+str(i)+"\")*xs_weight*puWeight*"+prefire_weight_string+"*" + get_postfilter_selection_string()+" : 0")
    rinterface = rinterface.Define("double_fake_stat_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"stat_up\")*xs_weight*puWeight*"+prefire_weight_string+"*" + get_postfilter_selection_string()+" : 0") 

    if label == "w+jets":
        rinterface = rinterface.Define("fake_photon_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pileup_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*pileup_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_prefire_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*prefire_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_photon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*photon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_electron_reco_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_reco_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_electron_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_electron_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_hlt_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_muon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_muon_iso_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_iso_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_muon_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_hlt_sf_up_base_weight : 0")
        for j in range(n_fake_photon_alt) :
            rinterface = rinterface.Define("fake_photon_alt"+str(j)+"_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"alt"+str(j)+"\")*base_weight : 0")
        rinterface = rinterface.Define("fake_photon_stat_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !photon_genjet_matching && is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0]) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"stat_up\")*xs_weight*puWeight*"+prefire_weight_string+"*" + get_postfilter_selection_string()+" : 0")
    else:
        rinterface = rinterface.Define("fake_photon_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*base_weight : 0")
        rinterface = rinterface.Define("fake_photon_pileup_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*pileup_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_prefire_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*prefire_up_base_weight : 0")
        if label != "gg+jets":
            rinterface = rinterface.Define("fake_photon_jes_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*jes_up_base_weight : 0")
            rinterface = rinterface.Define("fake_photon_jer_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*jer_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_photon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*photon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_muon_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_muon_iso_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_iso_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_muon_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*muon_hlt_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_electron_reco_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_reco_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_electron_id_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_id_sf_up_base_weight : 0")
        rinterface = rinterface.Define("fake_photon_electron_hlt_sf_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id)*electron_hlt_sf_up_base_weight : 0")
        for j in range(n_fake_photon_alt) :
            rinterface = rinterface.Define("fake_photon_alt"+str(j)+"_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"alt"+str(j)+"\")*base_weight : 0")
        rinterface = rinterface.Define("fake_photon_stat_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_for_fake_cutstring+" ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"stat_up\")*xs_weight*puWeight*"+prefire_weight_string+"*" + get_postfilter_selection_string()+" : 0")



    if label == "w+jets" and year == "2016":
        if "is_lepton_real" in rinterface_wjets_2016.GetColumnNames():
            rinterface_wjets_2016 = rinterface_wjets_2016.Define("pass_lepton_gen","is_lepton_real == 1")
        else:    
            rinterface_wjets_2016 = rinterface_wjets_2016.Define("pass_lepton_gen","lepton_gen_matching == 1 || lepton_gen_matching == 15")
        rinterface_wjets_2016 = rinterface_wjets_2016.Define("xs_weight",str(sample["xs"]*1000*lumi/sample["nweightedevents"]) + "*gen_weight/abs(gen_weight)") 
        rinterface_wjets_2016 = rinterface_wjets_2016.Define("base_weight",get_postfilter_selection_string()+"*xs_weight*puWeight*"+prefire_weight_string+"*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\") : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\"))")      
#        rinterface_wjets_2016 = rinterface_wjets_2016.Define("weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && "+photon_gen_matching_cutstring + " && photon_genjet_matching)*base_weight")
        rinterface_wjets_2016 = rinterface_wjets_2016.Define("weight","(photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && !(photon_gen_matching == 1|| photon_gen_matching == 4 || photon_gen_matching == 5 || photon_gen_matching == 6) && (photon_genjet_matching || !is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0])))*base_weight")
        rinterface_wjets_2016 = rinterface_wjets_2016.Define("wjets_fake_photon_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 && pass_lepton_gen && !(photon_gen_matching == 1|| photon_gen_matching == 4 || photon_gen_matching == 5 || photon_gen_matching == 6) && (photon_genjet_matching || !is_photon_prompt(lumi,event,\""+year+"\",dsetversion[0])) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"wjets\")*xs_weight*puWeight*"+prefire_weight_string+"*" + get_postfilter_selection_string()+": 0")
        rinterface_wjets_2016 = rinterface_wjets_2016.Define("wjets_chiso_fake_photon_weight","photon_selection == 3 && ((abs(photon_eta) < 1.5 && photon_pfRelIso03_chg*photon_pt < "+str(chiso_cut_barrel)+"*1.75) || (abs(photon_eta) > 1.5 && photon_pfRelIso03_chg*photon_pt < "+str(chiso_cut_endcap)+"*1.75)) && is_lepton_tight == 1 && pass_lepton_gen && !(photon_gen_matching == 1|| photon_gen_matching == 4 || photon_gen_matching == 5 || photon_gen_matching == 6) ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"wjets_chiso\")*xs_weight*puWeight*"+prefire_weight_string+"*" + get_postfilter_selection_string()+" : 0")
        for variable_definition in variable_definitions:
            rinterface_wjets_2016 = rinterface_wjets_2016.Define(variable_definition[0],variable_definition[1])


    if sample["e_to_p"] or sample["e_to_p_non_res"]:
        for i in range(len(etopbinning)):
            if i == len(etopbinning) - 1:
                photon_pt = "(photon_pt > "+str(etopbinning[i])+")"
            else:    
                photon_pt = "("+str(etopbinning[i])+" < photon_pt && photon_pt < "+str(etopbinning[i+1])+")"
            rinterface = rinterface.Define("e_to_p_bin"+str(i)+"_weight",photon_pt+" && photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && photon_gen_matching == 1 ? base_weight : 0")
            rinterface = rinterface.Define("e_to_p_bin"+str(i)+"_pileup_up_weight","("+photon_pt+" && photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && photon_gen_matching == 1)*pileup_up_base_weight")
            rinterface = rinterface.Define("e_to_p_bin"+str(i)+"_prefire_up_weight","("+photon_pt+" && photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && photon_gen_matching == 1)*prefire_up_base_weight")
            rinterface = rinterface.Define("e_to_p_bin"+str(i)+"_electron_id_sf_up_weight","("+photon_pt+" && photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && photon_gen_matching == 1)*electron_id_sf_up_base_weight")
            rinterface = rinterface.Define("e_to_p_bin"+str(i)+"_electron_reco_sf_up_weight","("+photon_pt+" && photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && photon_gen_matching == 1)*electron_reco_sf_up_base_weight")
            rinterface = rinterface.Define("e_to_p_bin"+str(i)+"_electron_hlt_sf_up_weight","("+photon_pt+" && photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && photon_gen_matching == 1)*electron_hlt_sf_up_base_weight")
            rinterface = rinterface.Define("e_to_p_bin"+str(i)+"_photon_id_sf_up_weight","("+photon_pt+" && photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && photon_gen_matching == 1)*photon_id_sf_up_base_weight")
            rinterface = rinterface.Define("e_to_p_bin"+str(i)+"_jes_up_weight","("+photon_pt+" && photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && photon_gen_matching == 1)*jes_up_base_weight")
            rinterface = rinterface.Define("e_to_p_bin"+str(i)+"_jer_up_weight","("+photon_pt+" && photon_selection == 0 && is_lepton_tight == 1 && pass_lepton_gen && photon_gen_matching == 1)*jer_up_base_weight")
                
    for variable_definition in variable_definitions:
        rinterface = rinterface.Define(variable_definition[0],variable_definition[1])


    rresultptrs = {}

    for i in range(len(variables)):

        rresultptrs[i] = {}

        if labels[label]["color"] != None:
            rresultptrs[i][""] =rinterface.Histo1D(histogram_models[i],variables[i],"weight")
            rresultptrs[i]["pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"pileup_up_weight")
            rresultptrs[i]["prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"prefire_up_weight")
            rresultptrs[i]["muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"muon_id_sf_up_weight")
            rresultptrs[i]["muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"muon_iso_sf_up_weight")
            rresultptrs[i]["muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"muon_hlt_sf_up_weight")
            rresultptrs[i]["electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"electron_id_sf_up_weight")
            rresultptrs[i]["electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"electron_reco_sf_up_weight")
            rresultptrs[i]["electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"electron_hlt_sf_up_weight")
            rresultptrs[i]["photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"photon_id_sf_up_weight")
            if label == "w+jets":
                rresultptrs[i]["prompt_pileup"] =rinterface.Histo1D(histogram_models[i],variables[i],"prompt_pileup_weight")
                rresultptrs[i]["prompt_pileup_pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"prompt_pileup_pileup_up_weight")
                rresultptrs[i]["prompt_pileup_prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"prompt_pileup_prefire_up_weight")
                rresultptrs[i]["prompt_pileup_photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"prompt_pileup_photon_id_sf_up_weight")
                rresultptrs[i]["prompt_pileup_electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"prompt_pileup_electron_reco_sf_up_weight")
                rresultptrs[i]["prompt_pileup_electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"prompt_pileup_electron_id_sf_up_weight")
                rresultptrs[i]["prompt_pileup_electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"prompt_pileup_electron_hlt_sf_up_weight")
                rresultptrs[i]["prompt_pileup_muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"prompt_pileup_muon_id_sf_up_weight")
                rresultptrs[i]["prompt_pileup_muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"prompt_pileup_muon_iso_sf_up_weight")
                rresultptrs[i]["prompt_pileup_muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"prompt_pileup_muon_hlt_sf_up_weight")
            if label == "wg+jets":
                rresultptrs[i]["pass_fiducial"] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_weight")
                rresultptrs[i]["fail_fiducial"] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_weight")
                rresultptrs[i]["pass_fiducial_pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_pileup_up_weight")
                rresultptrs[i]["fail_fiducial_pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_pileup_up_weight")
                rresultptrs[i]["pass_fiducial_prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_prefire_up_weight")
                rresultptrs[i]["fail_fiducial_prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_prefire_up_weight")
                rresultptrs[i]["pass_fiducial_jes_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_jes_up_weight")
                rresultptrs[i]["fail_fiducial_jes_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_jes_up_weight")
                rresultptrs[i]["pass_fiducial_jer_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_jer_up_weight")
                rresultptrs[i]["fail_fiducial_jer_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_jer_up_weight")
                rresultptrs[i]["pass_fiducial_photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_photon_id_sf_up_weight")
                rresultptrs[i]["fail_fiducial_photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_photon_id_sf_up_weight")
                rresultptrs[i]["pass_fiducial_electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_electron_reco_sf_up_weight")
                rresultptrs[i]["fail_fiducial_electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_electron_reco_sf_up_weight")
                rresultptrs[i]["pass_fiducial_electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_electron_id_sf_up_weight")
                rresultptrs[i]["fail_fiducial_electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_electron_id_sf_up_weight")
                rresultptrs[i]["pass_fiducial_electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_electron_hlt_sf_up_weight")
                rresultptrs[i]["fail_fiducial_electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_electron_hlt_sf_up_weight")
                rresultptrs[i]["pass_fiducial_muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_muon_id_sf_up_weight")
                rresultptrs[i]["fail_fiducial_muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_muon_id_sf_up_weight")
                rresultptrs[i]["pass_fiducial_muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_muon_iso_sf_up_weight")
                rresultptrs[i]["fail_fiducial_muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_muon_iso_sf_up_weight")
                rresultptrs[i]["pass_fiducial_muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_muon_hlt_sf_up_weight")
                rresultptrs[i]["fail_fiducial_muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_muon_hlt_sf_up_weight")
                
                rresultptrs[i]["fake_photon_pass_fiducial"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_weight")
                rresultptrs[i]["fake_photon_fail_fiducial"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_weight")
                rresultptrs[i]["fake_photon_pass_fiducial_pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_pileup_up_weight")
                rresultptrs[i]["fake_photon_fail_fiducial_pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_pileup_up_weight")
                rresultptrs[i]["fake_photon_pass_fiducial_prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_prefire_up_weight")
                rresultptrs[i]["fake_photon_fail_fiducial_prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_prefire_up_weight")
                rresultptrs[i]["fake_photon_pass_fiducial_jes_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_jes_up_weight")
                rresultptrs[i]["fake_photon_fail_fiducial_jes_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_jes_up_weight")
                rresultptrs[i]["fake_photon_pass_fiducial_jer_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_jer_up_weight")
                rresultptrs[i]["fake_photon_fail_fiducial_jer_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_jer_up_weight")
                rresultptrs[i]["fake_photon_pass_fiducial_photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_photon_id_sf_up_weight")
                rresultptrs[i]["fake_photon_fail_fiducial_photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_photon_id_sf_up_weight")
                rresultptrs[i]["fake_photon_pass_fiducial_electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_electron_reco_sf_up_weight")
                rresultptrs[i]["fake_photon_fail_fiducial_electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_electron_reco_sf_up_weight")
                rresultptrs[i]["fake_photon_pass_fiducial_electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_electron_id_sf_up_weight")
                rresultptrs[i]["fake_photon_fail_fiducial_electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_electron_id_sf_up_weight")
                rresultptrs[i]["fake_photon_pass_fiducial_electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_electron_hlt_sf_up_weight")
                rresultptrs[i]["fake_photon_fail_fiducial_electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_electron_hlt_sf_up_weight")
                rresultptrs[i]["fake_photon_pass_fiducial_muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_muon_id_sf_up_weight")
                rresultptrs[i]["fake_photon_fail_fiducial_muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_muon_id_sf_up_weight")
                rresultptrs[i]["fake_photon_pass_fiducial_muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_muon_iso_sf_up_weight")
                rresultptrs[i]["fake_photon_fail_fiducial_muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_muon_iso_sf_up_weight")
                rresultptrs[i]["fake_photon_pass_fiducial_muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_muon_hlt_sf_up_weight")
                rresultptrs[i]["fake_photon_fail_fiducial_muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_muon_hlt_sf_up_weight")
                for j in range(n_fake_photon_alt):
                    rresultptrs[i]["fake_photon_pass_fiducial_fake_photon_alt"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pass_fiducial_fake_photon_alt"+str(j)+"_weight")
                    rresultptrs[i]["fake_photon_fail_fiducial_fake_photon_alt"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_fail_fiducial_fake_photon_alt"+str(j)+"_weight")

                
                rresultptrs[i]["fake_lepton_pass_fiducial"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pass_fiducial_weight")
                rresultptrs[i]["fake_lepton_fail_fiducial"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_fail_fiducial_weight")
                rresultptrs[i]["fake_lepton_pass_fiducial_pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pass_fiducial_pileup_up_weight")
                rresultptrs[i]["fake_lepton_fail_fiducial_pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_fail_fiducial_pileup_up_weight")
                rresultptrs[i]["fake_lepton_pass_fiducial_prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pass_fiducial_prefire_up_weight")
                rresultptrs[i]["fake_lepton_fail_fiducial_prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_fail_fiducial_prefire_up_weight")
                rresultptrs[i]["fake_lepton_pass_fiducial_jes_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pass_fiducial_jes_up_weight")
                rresultptrs[i]["fake_lepton_fail_fiducial_jes_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_fail_fiducial_jes_up_weight")
                rresultptrs[i]["fake_lepton_pass_fiducial_jer_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pass_fiducial_jer_up_weight")
                rresultptrs[i]["fake_lepton_fail_fiducial_jer_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_fail_fiducial_jer_up_weight")
                rresultptrs[i]["fake_lepton_pass_fiducial_photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pass_fiducial_photon_id_sf_up_weight")
                rresultptrs[i]["fake_lepton_fail_fiducial_photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_fail_fiducial_photon_id_sf_up_weight")
                rresultptrs[i]["fake_lepton_pass_fiducial_electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pass_fiducial_electron_reco_sf_up_weight")
                rresultptrs[i]["fake_lepton_fail_fiducial_electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_fail_fiducial_electron_reco_sf_up_weight")
                rresultptrs[i]["fake_lepton_pass_fiducial_electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pass_fiducial_electron_id_sf_up_weight")
                rresultptrs[i]["fake_lepton_fail_fiducial_electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_fail_fiducial_electron_id_sf_up_weight")
                rresultptrs[i]["fake_lepton_pass_fiducial_electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pass_fiducial_electron_hlt_sf_up_weight")
                rresultptrs[i]["fake_lepton_fail_fiducial_electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_fail_fiducial_electron_hlt_sf_up_weight")
                rresultptrs[i]["fake_lepton_pass_fiducial_muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pass_fiducial_muon_id_sf_up_weight")
                rresultptrs[i]["fake_lepton_fail_fiducial_muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_fail_fiducial_muon_id_sf_up_weight")
                rresultptrs[i]["fake_lepton_pass_fiducial_muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pass_fiducial_muon_iso_sf_up_weight")
                rresultptrs[i]["fake_lepton_fail_fiducial_muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_fail_fiducial_muon_iso_sf_up_weight")
                rresultptrs[i]["fake_lepton_pass_fiducial_muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pass_fiducial_muon_hlt_sf_up_weight")
                rresultptrs[i]["fake_lepton_fail_fiducial_muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_fail_fiducial_muon_hlt_sf_up_weight")
                
                
                rresultptrs[i]["double_fake_pass_fiducial"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_weight")
                rresultptrs[i]["double_fake_fail_fiducial"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_weight")
                rresultptrs[i]["double_fake_pass_fiducial_pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_pileup_up_weight")
                rresultptrs[i]["double_fake_fail_fiducial_pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_pileup_up_weight")
                rresultptrs[i]["double_fake_pass_fiducial_prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_prefire_up_weight")
                rresultptrs[i]["double_fake_fail_fiducial_prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_prefire_up_weight")
                rresultptrs[i]["double_fake_pass_fiducial_jes_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_jes_up_weight")
                rresultptrs[i]["double_fake_fail_fiducial_jes_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_jes_up_weight")
                rresultptrs[i]["double_fake_pass_fiducial_jer_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_jer_up_weight")
                rresultptrs[i]["double_fake_fail_fiducial_jer_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_jer_up_weight")
                rresultptrs[i]["double_fake_pass_fiducial_photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_photon_id_sf_up_weight")
                rresultptrs[i]["double_fake_fail_fiducial_photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_photon_id_sf_up_weight")
                rresultptrs[i]["double_fake_pass_fiducial_electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_electron_reco_sf_up_weight")
                rresultptrs[i]["double_fake_fail_fiducial_electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_electron_reco_sf_up_weight")
                rresultptrs[i]["double_fake_pass_fiducial_electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_electron_id_sf_up_weight")
                rresultptrs[i]["double_fake_fail_fiducial_electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_electron_id_sf_up_weight")
                rresultptrs[i]["double_fake_pass_fiducial_electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_electron_hlt_sf_up_weight")
                rresultptrs[i]["double_fake_fail_fiducial_electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_electron_hlt_sf_up_weight")
                rresultptrs[i]["double_fake_pass_fiducial_muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_muon_id_sf_up_weight")
                rresultptrs[i]["double_fake_fail_fiducial_muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_muon_id_sf_up_weight")
                rresultptrs[i]["double_fake_pass_fiducial_muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_muon_iso_sf_up_weight")
                rresultptrs[i]["double_fake_fail_fiducial_muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_muon_iso_sf_up_weight")
                rresultptrs[i]["double_fake_pass_fiducial_muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_muon_hlt_sf_up_weight")
                rresultptrs[i]["double_fake_fail_fiducial_muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_muon_hlt_sf_up_weight")
                for j in range(n_fake_photon_alt):
                    rresultptrs[i]["double_fake_pass_fiducial_fake_photon_alt"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pass_fiducial_fake_photon_alt"+str(j)+"_weight")
                    rresultptrs[i]["double_fake_fail_fiducial_fake_photon_alt"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_fail_fiducial_fake_photon_alt"+str(j)+"_weight")
                
                if labels["wg+jets"]["syst-scale"]:
                    for j in range(0,8):
                        rresultptrs[i]["pass_fiducial_scale"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_scale"+str(j)+"_weight")
                        rresultptrs[i]["fail_fiducial_scale"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_scale"+str(j)+"_weight")
                if labels["wg+jets"]["syst-pdf"]:
                    for j in range(0,32):
                        rresultptrs[i]["pass_fiducial_pdf"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"pass_fiducial_pdf"+str(j)+"_weight")
                        rresultptrs[i]["fail_fiducial_pdf"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"fail_fiducial_pdf"+str(j)+"_weight")

            if label != "w+jets" and label != "gg+jets":
                rresultptrs[i]["jes_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"jes_up_weight")
                rresultptrs[i]["jer_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"jer_up_weight")

        if label == "w+jets" and year == "2016":
            rresultptrs[i]["wjets_fake_photon_chiso"] = rinterface_wjets_2016.Histo1D(histogram_models[i],variables[i],"wjets_chiso_fake_photon_weight")        
            rresultptrs[i]["wjets_fake_photon"] = rinterface_wjets_2016.Histo1D(histogram_models[i],variables[i],"wjets_fake_photon_weight")    
            rresultptrs[i]["wjets"] = rinterface_wjets_2016.Histo1D(histogram_models[i],variables[i],"weight")    

        rresultptrs[i]["fake_photon"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_weight")
        rresultptrs[i]["fake_photon_pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_pileup_up_weight")
        rresultptrs[i]["fake_photon_prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_prefire_up_weight")
        if label != "w+jets" and label != "gg+jets":
            rresultptrs[i]["fake_photon_jes_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_jes_up_weight")
            rresultptrs[i]["fake_photon_jer_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_jer_up_weight")
        rresultptrs[i]["fake_photon_photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_photon_id_sf_up_weight")
        rresultptrs[i]["fake_photon_electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_electron_reco_sf_up_weight")
        rresultptrs[i]["fake_photon_electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_electron_id_sf_up_weight")
        rresultptrs[i]["fake_photon_electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_electron_hlt_sf_up_weight")
        rresultptrs[i]["fake_photon_muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_muon_id_sf_up_weight")
        rresultptrs[i]["fake_photon_muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_muon_iso_sf_up_weight")
        rresultptrs[i]["fake_photon_muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_muon_hlt_sf_up_weight")
        for j in range(n_fake_photon_alt) :
            rresultptrs[i]["fake_photon_alt"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_alt"+str(j)+"_weight")
#            rresultptrs[i]["fake_photon_alt"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_weight")
        rresultptrs[i]["fake_photon_stat_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_stat_up_weight")
        
        rresultptrs[i]["fake_lepton"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_weight")
        rresultptrs[i]["fake_lepton_stat_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_stat_up_weight")
        rresultptrs[i]["fake_lepton_stat_down"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_stat_down_weight")
        rresultptrs[i]["fake_lepton_pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_pileup_up_weight")
        rresultptrs[i]["fake_lepton_prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_prefire_up_weight")
        if label != "w+jets" and label != "gg+jets":
            rresultptrs[i]["fake_lepton_jes_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_jes_up_weight")
            rresultptrs[i]["fake_lepton_jer_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_jer_up_weight")
        rresultptrs[i]["fake_lepton_photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_photon_id_sf_up_weight")
        rresultptrs[i]["fake_lepton_electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_electron_reco_sf_up_weight")
        rresultptrs[i]["fake_lepton_electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_electron_id_sf_up_weight")
        rresultptrs[i]["fake_lepton_electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_electron_hlt_sf_up_weight")
        rresultptrs[i]["fake_lepton_muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_muon_id_sf_up_weight")
        rresultptrs[i]["fake_lepton_muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_muon_iso_sf_up_weight")
        rresultptrs[i]["fake_lepton_muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_muon_hlt_sf_up_weight")

        rresultptrs[i]["double_fake"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_weight")
        for j in range(n_fake_photon_alt):
            rresultptrs[i]["double_fake_alt"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_alt"+str(j)+"_weight")
        rresultptrs[i]["double_fake_stat_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_stat_up_weight")
        rresultptrs[i]["double_fake_pileup_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_pileup_up_weight")
        rresultptrs[i]["double_fake_prefire_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_prefire_up_weight")
        if label != "w+jets" and label != "gg+jets":
            rresultptrs[i]["double_fake_jes_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_jes_up_weight")
            rresultptrs[i]["double_fake_jer_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_jer_up_weight")
        rresultptrs[i]["double_fake_photon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_photon_id_sf_up_weight")
        rresultptrs[i]["double_fake_electron_reco_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_electron_reco_sf_up_weight")
        rresultptrs[i]["double_fake_electron_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_electron_id_sf_up_weight")
        rresultptrs[i]["double_fake_electron_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_electron_hlt_sf_up_weight")
        rresultptrs[i]["double_fake_muon_id_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_muon_id_sf_up_weight")
        rresultptrs[i]["double_fake_muon_iso_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_muon_iso_sf_up_weight")
        rresultptrs[i]["double_fake_muon_hlt_sf_up"] =rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_muon_hlt_sf_up_weight")

        if labels[label]["syst-scale"]:
            for j in range(0,8):
                rresultptrs[i]["scale"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"scale"+str(j)+"_weight")
        if labels[label]["syst-pdf"]:
            for j in range(0,32):
                rresultptrs[i]["pdf"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"pdf"+str(j)+"_weight")

        if sample["e_to_p"] or sample["e_to_p_non_res"]:
            for j in range(len(etopbinning)): 
                rresultptrs[i]["e_to_p"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"e_to_p_bin"+str(j)+"_weight")
                rresultptrs[i]["e_to_p_electron_id_sf_up"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"e_to_p_bin"+str(j)+"_electron_id_sf_up_weight")
                rresultptrs[i]["e_to_p_electron_reco_sf_up"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"e_to_p_bin"+str(j)+"_electron_reco_sf_up_weight")
                rresultptrs[i]["e_to_p_electron_hlt_sf_up"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"e_to_p_bin"+str(j)+"_electron_hlt_sf_up_weight")
                rresultptrs[i]["e_to_p_photon_id_sf_up"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"e_to_p_bin"+str(j)+"_photon_id_sf_up_weight")
                rresultptrs[i]["e_to_p_pileup_up"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"e_to_p_bin"+str(j)+"_pileup_up_weight")
                rresultptrs[i]["e_to_p_prefire_up"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"e_to_p_bin"+str(j)+"_prefire_up_weight")
                rresultptrs[i]["e_to_p_jes_up"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"e_to_p_bin"+str(j)+"_jes_up_weight")
                rresultptrs[i]["e_to_p_jer_up"+str(j)] =rinterface.Histo1D(histogram_models[i],variables[i],"e_to_p_bin"+str(j)+"_jer_up_weight")

    results = {}
    for i in range(len(variables)):
        results[i] = {}
        for key,value in rresultptrs[i].iteritems():
            results[i][key] = ROOT.TH1D(value.GetValue())

    return results
    
def processMCSampleStar(inputs):
    return processMCSample(*inputs)

for year in years:
    for label in labels.keys():

#        if label != "wg+jets":
#            continue

        if label == "w+jets" and (year == "2017" or year == "2018") and args.no_wjets_for_2017_and_2018:
            continue

        if year == "2016":
            lumi=35.9
        elif year == "2017":
            lumi=41.5
        elif year == "2018":
            lumi=59.6
        else:
            assert(0)

        if year == "2016":
            sieie_cut_barrel = sieie_cut_2016_barrel
            sieie_cut_endcap = sieie_cut_2016_endcap
            chiso_cut_barrel = chiso_cut_2016_barrel
            chiso_cut_endcap = chiso_cut_2016_endcap
        elif year == "2017":
            sieie_cut_barrel = sieie_cut_2017_barrel
            sieie_cut_endcap = sieie_cut_2017_endcap
            chiso_cut_barrel = chiso_cut_2017_barrel
            chiso_cut_endcap = chiso_cut_2017_endcap
        elif year == "2018":
            sieie_cut_barrel = sieie_cut_2018_barrel
            sieie_cut_endcap = sieie_cut_2018_endcap
            chiso_cut_barrel = chiso_cut_2018_barrel
            chiso_cut_endcap = chiso_cut_2018_endcap
        else:
            assert(0)

        fake_photon_sieie_cut_barrel = sieie_cut_barrel*1.75
        fake_photon_sieie_cut_endcap = sieie_cut_endcap*1.75
        fake_photon_chiso_cut_barrel = chiso_cut_barrel*1000
        fake_photon_chiso_cut_endcap = chiso_cut_endcap*1000    

        fake_photon_sieie_cut_cutstring = "((abs(photon_eta) < 1.5 && photon_sieie < "+str(fake_photon_sieie_cut_barrel)+ ") || (abs(photon_eta) > 1.5 && photon_sieie < "+str(fake_photon_sieie_cut_endcap)+ "))" 
        fake_photon_chiso_cut_cutstring = "((abs(photon_eta) < 1.5 && photon_pfRelIso03_chg*photon_pt < "+str(fake_photon_chiso_cut_barrel)+ ") || (abs(photon_eta) > 1.5 && photon_pfRelIso03_chg*photon_pt < "+str(fake_photon_chiso_cut_endcap)+ "))" 

        for sample in labels[label]["samples"][year]:

            print "Running over sample " + str(sample["filename"])

            if args.singleproc:
                results = processMCSample(None)
            else:    
                #process each sample in its own thread in order to avoid memory problems
                #don't create any RDataFrames before this, otherwise the RDataFrames in the new processes will only use 1 core
                import multiprocessing
                pool = multiprocessing.Pool(1)
                import itertools

                try:
                    res = pool.map_async(processMCSampleStar, itertools.izip([None]))
                    results = res.get(10000)[0]
                except Exception as e:
                    print "Exception:"
                    print e
                    pool.terminate()
                    sys.exit(0)
                except:
                    e = sys.exc_info()[0]
                    print "exception:"
                    print e
                    pool.terminate()
                    sys.exit(0)
                else:    
                    pool.terminate()

            for i in range(len(variables)):

                if labels[label]["color"] != None:
                    labels[label]["hists"][i].Add(results[i][""])
                    labels[label]["hists-pileup-up"][i].Add(results[i]["pileup_up"])
                    labels[label]["hists-prefire-up"][i].Add(results[i]["prefire_up"])
                    labels[label]["hists-electron-id-sf-up"][i].Add(results[i]["electron_id_sf_up"])
                    labels[label]["hists-electron-reco-sf-up"][i].Add(results[i]["electron_reco_sf_up"])
                    labels[label]["hists-electron-hlt-sf-up"][i].Add(results[i]["electron_hlt_sf_up"])
                    labels[label]["hists-muon-id-sf-up"][i].Add(results[i]["muon_id_sf_up"])
                    labels[label]["hists-muon-iso-sf-up"][i].Add(results[i]["muon_iso_sf_up"])
                    labels[label]["hists-muon-hlt-sf-up"][i].Add(results[i]["muon_hlt_sf_up"])
                    labels[label]["hists-photon-id-sf-up"][i].Add(results[i]["photon_id_sf_up"])
                    if label == "w+jets":
                        labels[label]["hists-prompt-pileup"][i].Add(results[i]["prompt_pileup"])
                        labels[label]["hists-prompt-pileup-pileup-up"][i].Add(results[i]["prompt_pileup_pileup_up"])
                        labels[label]["hists-prompt-pileup-prefire-up"][i].Add(results[i]["prompt_pileup_prefire_up"])
                        labels[label]["hists-prompt-pileup-photon-id-sf-up"][i].Add(results[i]["prompt_pileup_photon_id_sf_up"])
                        labels[label]["hists-prompt-pileup-electron-reco-sf-up"][i].Add(results[i]["prompt_pileup_electron_reco_sf_up"])
                        labels[label]["hists-prompt-pileup-electron-id-sf-up"][i].Add(results[i]["prompt_pileup_electron_id_sf_up"])
                        labels[label]["hists-prompt-pileup-electron-hlt-sf-up"][i].Add(results[i]["prompt_pileup_electron_hlt_sf_up"])
                        labels[label]["hists-prompt-pileup-muon-id-sf-up"][i].Add(results[i]["prompt_pileup_muon_id_sf_up"])
                        labels[label]["hists-prompt-pileup-muon-iso-sf-up"][i].Add(results[i]["prompt_pileup_muon_iso_sf_up"])
                        labels[label]["hists-prompt-pileup-muon-hlt-sf-up"][i].Add(results[i]["prompt_pileup_muon_hlt_sf_up"])
                    if label == "wg+jets":
                        labels[label]["hists-pass-fiducial"][i].Add(results[i]["pass_fiducial"])
                        labels[label]["hists-fail-fiducial"][i].Add(results[i]["fail_fiducial"])
                        labels[label]["hists-pass-fiducial-pileup-up"][i].Add(results[i]["pass_fiducial_pileup_up"])
                        labels[label]["hists-fail-fiducial-pileup-up"][i].Add(results[i]["fail_fiducial_pileup_up"])
                        labels[label]["hists-pass-fiducial-prefire-up"][i].Add(results[i]["pass_fiducial_prefire_up"])
                        labels[label]["hists-fail-fiducial-prefire-up"][i].Add(results[i]["fail_fiducial_prefire_up"])
                        labels[label]["hists-pass-fiducial-jes-up"][i].Add(results[i]["pass_fiducial_jes_up"])
                        labels[label]["hists-fail-fiducial-jes-up"][i].Add(results[i]["fail_fiducial_jes_up"])
                        labels[label]["hists-pass-fiducial-jer-up"][i].Add(results[i]["pass_fiducial_jer_up"])
                        labels[label]["hists-fail-fiducial-jer-up"][i].Add(results[i]["fail_fiducial_jer_up"])
                        labels[label]["hists-pass-fiducial-photon-id-sf-up"][i].Add(results[i]["pass_fiducial_photon_id_sf_up"])
                        labels[label]["hists-fail-fiducial-photon-id-sf-up"][i].Add(results[i]["fail_fiducial_photon_id_sf_up"])
                        labels[label]["hists-pass-fiducial-electron-reco-sf-up"][i].Add(results[i]["pass_fiducial_electron_reco_sf_up"])
                        labels[label]["hists-fail-fiducial-electron-reco-sf-up"][i].Add(results[i]["fail_fiducial_electron_reco_sf_up"])
                        labels[label]["hists-pass-fiducial-electron-id-sf-up"][i].Add(results[i]["pass_fiducial_electron_id_sf_up"])
                        labels[label]["hists-fail-fiducial-electron-id-sf-up"][i].Add(results[i]["fail_fiducial_electron_id_sf_up"])
                        labels[label]["hists-pass-fiducial-electron-hlt-sf-up"][i].Add(results[i]["pass_fiducial_electron_hlt_sf_up"])
                        labels[label]["hists-fail-fiducial-electron-hlt-sf-up"][i].Add(results[i]["fail_fiducial_electron_hlt_sf_up"])
                        labels[label]["hists-pass-fiducial-muon-id-sf-up"][i].Add(results[i]["pass_fiducial_muon_id_sf_up"])
                        labels[label]["hists-fail-fiducial-muon-id-sf-up"][i].Add(results[i]["fail_fiducial_muon_id_sf_up"])
                        labels[label]["hists-pass-fiducial-muon-iso-sf-up"][i].Add(results[i]["pass_fiducial_muon_iso_sf_up"])
                        labels[label]["hists-fail-fiducial-muon-iso-sf-up"][i].Add(results[i]["fail_fiducial_muon_iso_sf_up"])
                        labels[label]["hists-pass-fiducial-muon-hlt-sf-up"][i].Add(results[i]["pass_fiducial_muon_hlt_sf_up"])
                        labels[label]["hists-fail-fiducial-muon-hlt-sf-up"][i].Add(results[i]["fail_fiducial_muon_hlt_sf_up"])

                        if labels["wg+jets"]["syst-scale"]:
                            for j in range(0,8):
                                labels["wg+jets"]["hists-pass-fiducial-scale-variation"+str(j)][i].Add(results[i]["pass_fiducial_scale"+str(j)])
                                labels["wg+jets"]["hists-fail-fiducial-scale-variation"+str(j)][i].Add(results[i]["fail_fiducial_scale"+str(j)])

                        if labels["wg+jets"]["syst-pdf"]:
                            for j in range(0,32):
                                labels["wg+jets"]["hists-pass-fiducial-pdf-variation"+str(j)][i].Add(results[i]["pass_fiducial_pdf"+str(j)])
                                labels["wg+jets"]["hists-fail-fiducial-pdf-variation"+str(j)][i].Add(results[i]["fail_fiducial_pdf"+str(j)])

                    if label != "w+jets" and label != "gg+jets":
                        labels[label]["hists-jes-up"][i].Add(results[i]["jes_up"])
                        labels[label]["hists-jer-up"][i].Add(results[i]["jer_up"])

            for i in range(len(variables)):
                results[i]["fake_photon"].Scale(-1)
                results[i]["fake_photon_pileup_up"].Scale(-1)
                results[i]["fake_photon_prefire_up"].Scale(-1)
                if label != "w+jets" and label != "gg+jets":
                    results[i]["fake_photon_jes_up"].Scale(-1)
                    results[i]["fake_photon_jer_up"].Scale(-1)
                results[i]["fake_photon_photon_id_sf_up"].Scale(-1)
                results[i]["fake_photon_muon_id_sf_up"].Scale(-1)
                results[i]["fake_photon_muon_iso_sf_up"].Scale(-1)
                results[i]["fake_photon_muon_hlt_sf_up"].Scale(-1)
                results[i]["fake_photon_electron_reco_sf_up"].Scale(-1)
                results[i]["fake_photon_electron_id_sf_up"].Scale(-1)
                results[i]["fake_photon_electron_hlt_sf_up"].Scale(-1)
                for j in range(n_fake_photon_alt):
                    results[i]["fake_photon_alt"+str(j)].Scale(-1)
                results[i]["fake_photon_stat_up"].Scale(-1)

                results[i]["fake_lepton"].Scale(-1)
                results[i]["fake_lepton_pileup_up"].Scale(-1)
                results[i]["fake_lepton_prefire_up"].Scale(-1)
                if label != "w+jets" and label != "gg+jets":
                    results[i]["fake_lepton_jes_up"].Scale(-1)
                    results[i]["fake_lepton_jer_up"].Scale(-1)
                results[i]["fake_lepton_photon_id_sf_up"].Scale(-1)
                results[i]["fake_lepton_muon_id_sf_up"].Scale(-1)
                results[i]["fake_lepton_muon_iso_sf_up"].Scale(-1)
                results[i]["fake_lepton_muon_hlt_sf_up"].Scale(-1)
                results[i]["fake_lepton_electron_reco_sf_up"].Scale(-1)
                results[i]["fake_lepton_electron_id_sf_up"].Scale(-1)
                results[i]["fake_lepton_electron_hlt_sf_up"].Scale(-1)

                if label == "wg+jets":
                    results[i]["fake_photon_pass_fiducial"].Scale(-1)
                    results[i]["fake_photon_fail_fiducial"].Scale(-1)
                    results[i]["fake_photon_pass_fiducial_pileup_up"].Scale(-1)
                    results[i]["fake_photon_fail_fiducial_pileup_up"].Scale(-1)
                    results[i]["fake_photon_pass_fiducial_prefire_up"].Scale(-1)
                    results[i]["fake_photon_fail_fiducial_prefire_up"].Scale(-1)
                    results[i]["fake_photon_pass_fiducial_jes_up"].Scale(-1)
                    results[i]["fake_photon_fail_fiducial_jes_up"].Scale(-1)
                    results[i]["fake_photon_pass_fiducial_jer_up"].Scale(-1)
                    results[i]["fake_photon_fail_fiducial_jer_up"].Scale(-1)
                    results[i]["fake_photon_pass_fiducial_photon_id_sf_up"].Scale(-1)
                    results[i]["fake_photon_fail_fiducial_photon_id_sf_up"].Scale(-1)
                    results[i]["fake_photon_pass_fiducial_electron_reco_sf_up"].Scale(-1)
                    results[i]["fake_photon_fail_fiducial_electron_reco_sf_up"].Scale(-1)
                    results[i]["fake_photon_pass_fiducial_electron_id_sf_up"].Scale(-1)
                    results[i]["fake_photon_fail_fiducial_electron_id_sf_up"].Scale(-1)
                    results[i]["fake_photon_pass_fiducial_electron_hlt_sf_up"].Scale(-1)
                    results[i]["fake_photon_fail_fiducial_electron_hlt_sf_up"].Scale(-1)
                    results[i]["fake_photon_pass_fiducial_muon_id_sf_up"].Scale(-1)
                    results[i]["fake_photon_fail_fiducial_muon_id_sf_up"].Scale(-1)
                    results[i]["fake_photon_pass_fiducial_muon_iso_sf_up"].Scale(-1)
                    results[i]["fake_photon_fail_fiducial_muon_iso_sf_up"].Scale(-1)
                    results[i]["fake_photon_pass_fiducial_muon_hlt_sf_up"].Scale(-1)
                    results[i]["fake_photon_fail_fiducial_muon_hlt_sf_up"].Scale(-1)
                    for j in range(n_fake_photon_alt):
                        results[i]["fake_photon_pass_fiducial_fake_photon_alt"+str(j)].Scale(-1)
                        results[i]["fake_photon_fail_fiducial_fake_photon_alt"+str(j)].Scale(-1)

                    results[i]["fake_lepton_pass_fiducial"].Scale(-1)
                    results[i]["fake_lepton_fail_fiducial"].Scale(-1)
                    results[i]["fake_lepton_pass_fiducial_pileup_up"].Scale(-1)
                    results[i]["fake_lepton_fail_fiducial_pileup_up"].Scale(-1)
                    results[i]["fake_lepton_pass_fiducial_prefire_up"].Scale(-1)
                    results[i]["fake_lepton_fail_fiducial_prefire_up"].Scale(-1)
                    results[i]["fake_lepton_pass_fiducial_jes_up"].Scale(-1)
                    results[i]["fake_lepton_fail_fiducial_jes_up"].Scale(-1)
                    results[i]["fake_lepton_pass_fiducial_jer_up"].Scale(-1)
                    results[i]["fake_lepton_fail_fiducial_jer_up"].Scale(-1)
                    results[i]["fake_lepton_pass_fiducial_photon_id_sf_up"].Scale(-1)
                    results[i]["fake_lepton_fail_fiducial_photon_id_sf_up"].Scale(-1)
                    results[i]["fake_lepton_pass_fiducial_electron_reco_sf_up"].Scale(-1)
                    results[i]["fake_lepton_fail_fiducial_electron_reco_sf_up"].Scale(-1)
                    results[i]["fake_lepton_pass_fiducial_electron_id_sf_up"].Scale(-1)
                    results[i]["fake_lepton_fail_fiducial_electron_id_sf_up"].Scale(-1)
                    results[i]["fake_lepton_pass_fiducial_electron_hlt_sf_up"].Scale(-1)
                    results[i]["fake_lepton_fail_fiducial_electron_hlt_sf_up"].Scale(-1)
                    results[i]["fake_lepton_pass_fiducial_muon_id_sf_up"].Scale(-1)
                    results[i]["fake_lepton_fail_fiducial_muon_id_sf_up"].Scale(-1)
                    results[i]["fake_lepton_pass_fiducial_muon_iso_sf_up"].Scale(-1)
                    results[i]["fake_lepton_fail_fiducial_muon_iso_sf_up"].Scale(-1)
                    results[i]["fake_lepton_pass_fiducial_muon_hlt_sf_up"].Scale(-1)
                    results[i]["fake_lepton_fail_fiducial_muon_hlt_sf_up"].Scale(-1)

                if labels[label]["syst-scale"]:
                    for j in range(0,8):
                        labels[label]["hists-scale-variation"+str(j)][i].Add(results[i]["scale"+str(j)])

                if labels[label]["syst-pdf"]:
                    for j in range(0,32):
                        labels[label]["hists-pdf-variation"+str(j)][i].Add(results[i]["pdf"+str(j)])

                if label == "wg+jets":
                    labels[label]["hists-fake-photon-pass-fiducial"][i].Add(results[i]["fake_photon_pass_fiducial"])
                    labels[label]["hists-fake-photon-fail-fiducial"][i].Add(results[i]["fake_photon_fail_fiducial"])
                    labels[label]["hists-fake-photon-pass-fiducial-pileup-up"][i].Add(results[i]["fake_photon_pass_fiducial_pileup_up"])
                    labels[label]["hists-fake-photon-fail-fiducial-pileup-up"][i].Add(results[i]["fake_photon_fail_fiducial_pileup_up"])
                    labels[label]["hists-fake-photon-pass-fiducial-prefire-up"][i].Add(results[i]["fake_photon_pass_fiducial_prefire_up"])
                    labels[label]["hists-fake-photon-fail-fiducial-prefire-up"][i].Add(results[i]["fake_photon_fail_fiducial_prefire_up"])
                    labels[label]["hists-fake-photon-pass-fiducial-jes-up"][i].Add(results[i]["fake_photon_pass_fiducial_jes_up"])
                    labels[label]["hists-fake-photon-fail-fiducial-jes-up"][i].Add(results[i]["fake_photon_fail_fiducial_jes_up"])
                    labels[label]["hists-fake-photon-pass-fiducial-jer-up"][i].Add(results[i]["fake_photon_pass_fiducial_jer_up"])
                    labels[label]["hists-fake-photon-fail-fiducial-jer-up"][i].Add(results[i]["fake_photon_fail_fiducial_jer_up"])
                    labels[label]["hists-fake-photon-pass-fiducial-photon-id-sf-up"][i].Add(results[i]["fake_photon_pass_fiducial_photon_id_sf_up"])
                    labels[label]["hists-fake-photon-fail-fiducial-photon-id-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_photon_id_sf_up"])
                    labels[label]["hists-fake-photon-pass-fiducial-electron-reco-sf-up"][i].Add(results[i]["fake_photon_pass_fiducial_electron_reco_sf_up"])
                    labels[label]["hists-fake-photon-fail-fiducial-electron-reco-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_electron_reco_sf_up"])
                    labels[label]["hists-fake-photon-pass-fiducial-electron-id-sf-up"][i].Add(results[i]["fake_photon_pass_fiducial_electron_id_sf_up"])
                    labels[label]["hists-fake-photon-fail-fiducial-electron-id-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_electron_id_sf_up"])
                    labels[label]["hists-fake-photon-pass-fiducial-electron-hlt-sf-up"][i].Add(results[i]["fake_photon_pass_fiducial_electron_hlt_sf_up"])
                    labels[label]["hists-fake-photon-fail-fiducial-electron-hlt-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_electron_hlt_sf_up"])
                    labels[label]["hists-fake-photon-pass-fiducial-muon-id-sf-up"][i].Add(results[i]["fake_photon_pass_fiducial_muon_id_sf_up"])
                    labels[label]["hists-fake-photon-fail-fiducial-muon-id-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_muon_id_sf_up"])
                    labels[label]["hists-fake-photon-pass-fiducial-muon-iso-sf-up"][i].Add(results[i]["fake_photon_pass_fiducial_muon_iso_sf_up"])
                    labels[label]["hists-fake-photon-fail-fiducial-muon-iso-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_muon_iso_sf_up"])
                    labels[label]["hists-fake-photon-pass-fiducial-muon-hlt-sf-up"][i].Add(results[i]["fake_photon_pass_fiducial_muon_hlt_sf_up"])
                    labels[label]["hists-fake-photon-fail-fiducial-muon-hlt-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_muon_hlt_sf_up"])
                    for j in range(n_fake_photon_alt):
                        labels[label]["hists-fake-photon-pass-fiducial-fake-photon-alt"+str(j)][i].Add(results[i]["fake_photon_pass_fiducial_fake_photon_alt"+str(j)])
                        labels[label]["hists-fake-photon-fail-fiducial-fake-photon-alt"+str(j)][i].Add(results[i]["fake_photon_fail_fiducial_fake_photon_alt"+str(j)])

                    labels[label]["hists-fake-lepton-pass-fiducial"][i].Add(results[i]["fake_lepton_pass_fiducial"])
                    labels[label]["hists-fake-lepton-fail-fiducial"][i].Add(results[i]["fake_lepton_fail_fiducial"])
                    labels[label]["hists-fake-lepton-pass-fiducial-pileup-up"][i].Add(results[i]["fake_lepton_pass_fiducial_pileup_up"])
                    labels[label]["hists-fake-lepton-fail-fiducial-pileup-up"][i].Add(results[i]["fake_lepton_fail_fiducial_pileup_up"])
                    labels[label]["hists-fake-lepton-pass-fiducial-prefire-up"][i].Add(results[i]["fake_lepton_pass_fiducial_prefire_up"])
                    labels[label]["hists-fake-lepton-fail-fiducial-prefire-up"][i].Add(results[i]["fake_lepton_fail_fiducial_prefire_up"])
                    labels[label]["hists-fake-lepton-pass-fiducial-jes-up"][i].Add(results[i]["fake_lepton_pass_fiducial_jes_up"])
                    labels[label]["hists-fake-lepton-fail-fiducial-jes-up"][i].Add(results[i]["fake_lepton_fail_fiducial_jes_up"])
                    labels[label]["hists-fake-lepton-pass-fiducial-jer-up"][i].Add(results[i]["fake_lepton_pass_fiducial_jer_up"])
                    labels[label]["hists-fake-lepton-fail-fiducial-jer-up"][i].Add(results[i]["fake_lepton_fail_fiducial_jer_up"])
                    labels[label]["hists-fake-lepton-pass-fiducial-photon-id-sf-up"][i].Add(results[i]["fake_lepton_pass_fiducial_photon_id_sf_up"])
                    labels[label]["hists-fake-lepton-fail-fiducial-photon-id-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_photon_id_sf_up"])
                    labels[label]["hists-fake-lepton-pass-fiducial-electron-reco-sf-up"][i].Add(results[i]["fake_lepton_pass_fiducial_electron_reco_sf_up"])
                    labels[label]["hists-fake-lepton-fail-fiducial-electron-reco-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_electron_reco_sf_up"])
                    labels[label]["hists-fake-lepton-pass-fiducial-electron-id-sf-up"][i].Add(results[i]["fake_lepton_pass_fiducial_electron_id_sf_up"])
                    labels[label]["hists-fake-lepton-fail-fiducial-electron-id-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_electron_id_sf_up"])
                    labels[label]["hists-fake-lepton-pass-fiducial-electron-hlt-sf-up"][i].Add(results[i]["fake_lepton_pass_fiducial_electron_hlt_sf_up"])
                    labels[label]["hists-fake-lepton-fail-fiducial-electron-hlt-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_electron_hlt_sf_up"])
                    labels[label]["hists-fake-lepton-pass-fiducial-muon-id-sf-up"][i].Add(results[i]["fake_lepton_pass_fiducial_muon_id_sf_up"])
                    labels[label]["hists-fake-lepton-fail-fiducial-muon-id-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_muon_id_sf_up"])
                    labels[label]["hists-fake-lepton-pass-fiducial-muon-iso-sf-up"][i].Add(results[i]["fake_lepton_pass_fiducial_muon_iso_sf_up"])
                    labels[label]["hists-fake-lepton-fail-fiducial-muon-iso-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_muon_iso_sf_up"])
                    labels[label]["hists-fake-lepton-pass-fiducial-muon-hlt-sf-up"][i].Add(results[i]["fake_lepton_pass_fiducial_muon_hlt_sf_up"])
                    labels[label]["hists-fake-lepton-fail-fiducial-muon-hlt-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_muon_hlt_sf_up"])

                    labels[label]["hists-double-fake-pass-fiducial"][i].Add(results[i]["double_fake_pass_fiducial"])
                    labels[label]["hists-double-fake-fail-fiducial"][i].Add(results[i]["double_fake_fail_fiducial"])
                    labels[label]["hists-double-fake-pass-fiducial-pileup-up"][i].Add(results[i]["double_fake_pass_fiducial_pileup_up"])
                    labels[label]["hists-double-fake-fail-fiducial-pileup-up"][i].Add(results[i]["double_fake_fail_fiducial_pileup_up"])
                    labels[label]["hists-double-fake-pass-fiducial-prefire-up"][i].Add(results[i]["double_fake_pass_fiducial_prefire_up"])
                    labels[label]["hists-double-fake-fail-fiducial-prefire-up"][i].Add(results[i]["double_fake_fail_fiducial_prefire_up"])
                    labels[label]["hists-double-fake-pass-fiducial-jes-up"][i].Add(results[i]["double_fake_pass_fiducial_jes_up"])
                    labels[label]["hists-double-fake-fail-fiducial-jes-up"][i].Add(results[i]["double_fake_fail_fiducial_jes_up"])
                    labels[label]["hists-double-fake-pass-fiducial-jer-up"][i].Add(results[i]["double_fake_pass_fiducial_jer_up"])
                    labels[label]["hists-double-fake-fail-fiducial-jer-up"][i].Add(results[i]["double_fake_fail_fiducial_jer_up"])
                    labels[label]["hists-double-fake-pass-fiducial-photon-id-sf-up"][i].Add(results[i]["double_fake_pass_fiducial_photon_id_sf_up"])
                    labels[label]["hists-double-fake-fail-fiducial-photon-id-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_photon_id_sf_up"])
                    labels[label]["hists-double-fake-pass-fiducial-electron-reco-sf-up"][i].Add(results[i]["double_fake_pass_fiducial_electron_reco_sf_up"])
                    labels[label]["hists-double-fake-fail-fiducial-electron-reco-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_electron_reco_sf_up"])
                    labels[label]["hists-double-fake-pass-fiducial-electron-id-sf-up"][i].Add(results[i]["double_fake_pass_fiducial_electron_id_sf_up"])
                    labels[label]["hists-double-fake-fail-fiducial-electron-id-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_electron_id_sf_up"])
                    labels[label]["hists-double-fake-pass-fiducial-electron-hlt-sf-up"][i].Add(results[i]["double_fake_pass_fiducial_electron_hlt_sf_up"])
                    labels[label]["hists-double-fake-fail-fiducial-electron-hlt-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_electron_hlt_sf_up"])
                    labels[label]["hists-double-fake-pass-fiducial-muon-id-sf-up"][i].Add(results[i]["double_fake_pass_fiducial_muon_id_sf_up"])
                    labels[label]["hists-double-fake-fail-fiducial-muon-id-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_muon_id_sf_up"])
                    labels[label]["hists-double-fake-pass-fiducial-muon-iso-sf-up"][i].Add(results[i]["double_fake_pass_fiducial_muon_iso_sf_up"])
                    labels[label]["hists-double-fake-fail-fiducial-muon-iso-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_muon_iso_sf_up"])
                    labels[label]["hists-double-fake-pass-fiducial-muon-hlt-sf-up"][i].Add(results[i]["double_fake_pass_fiducial_muon_hlt_sf_up"])
                    labels[label]["hists-double-fake-fail-fiducial-muon-hlt-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_muon_hlt_sf_up"])
                    for j in range(n_fake_photon_alt):
                        labels[label]["hists-double-fake-pass-fiducial-fake-photon-alt"+str(j)][i].Add(results[i]["double_fake_pass_fiducial_fake_photon_alt"+str(j)])
                        labels[label]["hists-double-fake-fail-fiducial-fake-photon-alt"+str(j)][i].Add(results[i]["double_fake_fail_fiducial_fake_photon_alt"+str(j)])

                    fake_signal_contamination["hists"][i].Add(results[i]["fake_lepton"])
                    fake_signal_contamination["hists"][i].Add(results[i]["fake_photon"])
                    fake_signal_contamination["hists"][i].Add(results[i]["double_fake"])

                if year == "2016":    
                    fake_photon_2016["hists"][i].Add(results[i]["fake_photon"])

                if label == "w+jets" and year == "2016": 
                    wjets_fake_photon_chiso_2016["hists"][i].Add(results[i]["wjets_fake_photon_chiso"])
                    wjets_fake_photon_2016["hists"][i].Add(results[i]["wjets_fake_photon"])
                    wjets_2016["hists"][i].Add(results[i]["wjets"])

                if label != "wg+jets":

                    fake_photon["hists"][i].Add(results[i]["fake_photon"])
                    fake_photon["hists-pileup-up"][i].Add(results[i]["fake_photon_pileup_up"])
                    fake_photon["hists-prefire-up"][i].Add(results[i]["fake_photon_prefire_up"])
                    if label != "w+jets" and label != "gg+jets":
                        fake_photon["hists-jer-up"][i].Add(results[i]["fake_photon_jer_up"])
                        fake_photon["hists-jes-up"][i].Add(results[i]["fake_photon_jes_up"])
                    fake_photon["hists-photon-id-sf-up"][i].Add(results[i]["fake_photon_photon_id_sf_up"])
                    fake_photon["hists-electron-reco-sf-up"][i].Add(results[i]["fake_photon_electron_reco_sf_up"])
                    fake_photon["hists-electron-id-sf-up"][i].Add(results[i]["fake_photon_electron_id_sf_up"])
                    fake_photon["hists-electron-hlt-sf-up"][i].Add(results[i]["fake_photon_electron_hlt_sf_up"])
                    fake_photon["hists-muon-id-sf-up"][i].Add(results[i]["fake_photon_muon_id_sf_up"])
                    fake_photon["hists-muon-iso-sf-up"][i].Add(results[i]["fake_photon_muon_iso_sf_up"])
                    fake_photon["hists-muon-hlt-sf-up"][i].Add(results[i]["fake_photon_muon_hlt_sf_up"])
                    for j in range(n_fake_photon_alt):
                        fake_photon["hists-alt"+str(j)][i].Add(results[i]["fake_photon_alt"+str(j)])
                    fake_photon_stat_up["hists"][i].Add(results[i]["fake_photon_stat_up"])

                    fake_lepton["hists"][i].Add(results[i]["fake_lepton"])
                    fake_lepton["hists-pileup-up"][i].Add(results[i]["fake_lepton_pileup_up"])
                    fake_lepton["hists-prefire-up"][i].Add(results[i]["fake_lepton_prefire_up"])
                    if label != "w+jets" and label != "gg+jets":
                        fake_lepton["hists-jer-up"][i].Add(results[i]["fake_lepton_jer_up"])
                        fake_lepton["hists-jes-up"][i].Add(results[i]["fake_lepton_jes_up"])
                    fake_lepton["hists-photon-id-sf-up"][i].Add(results[i]["fake_lepton_photon_id_sf_up"])
                    fake_lepton["hists-electron-reco-sf-up"][i].Add(results[i]["fake_lepton_electron_reco_sf_up"])
                    fake_lepton["hists-electron-id-sf-up"][i].Add(results[i]["fake_lepton_electron_id_sf_up"])
                    fake_lepton["hists-electron-hlt-sf-up"][i].Add(results[i]["fake_lepton_electron_hlt_sf_up"])
                    fake_lepton["hists-muon-id-sf-up"][i].Add(results[i]["fake_lepton_muon_id_sf_up"])
                    fake_lepton["hists-muon-iso-sf-up"][i].Add(results[i]["fake_lepton_muon_iso_sf_up"])
                    fake_lepton["hists-muon-hlt-sf-up"][i].Add(results[i]["fake_lepton_muon_hlt_sf_up"])
                    fake_lepton_stat_up["hists"][i].Add(results[i]["fake_lepton_stat_up"])
                    fake_lepton_stat_down["hists"][i].Add(results[i]["fake_lepton_stat_down"])

                    double_fake["hists"][i].Add(results[i]["double_fake"])
                    for j in range(n_fake_photon_alt):
                        double_fake["hists-alt"+str(j)][i].Add(results[i]["double_fake_alt"+str(j)])
                    double_fake_stat_up["hists"][i].Add(results[i]["double_fake_stat_up"])
                    double_fake["hists-pileup-up"][i].Add(results[i]["double_fake_pileup_up"])
                    double_fake["hists-prefire-up"][i].Add(results[i]["double_fake_prefire_up"])
                    if label != "w+jets" and label != "gg+jets":
                        double_fake["hists-jer-up"][i].Add(results[i]["double_fake_jer_up"])
                        double_fake["hists-jes-up"][i].Add(results[i]["double_fake_jes_up"])
                    double_fake["hists-photon-id-sf-up"][i].Add(results[i]["double_fake_photon_id_sf_up"])
                    double_fake["hists-electron-reco-sf-up"][i].Add(results[i]["double_fake_electron_reco_sf_up"])
                    double_fake["hists-electron-id-sf-up"][i].Add(results[i]["double_fake_electron_id_sf_up"])
                    double_fake["hists-electron-hlt-sf-up"][i].Add(results[i]["double_fake_electron_hlt_sf_up"])
                    double_fake["hists-muon-id-sf-up"][i].Add(results[i]["double_fake_muon_id_sf_up"])
                    double_fake["hists-muon-iso-sf-up"][i].Add(results[i]["double_fake_muon_iso_sf_up"])
                    double_fake["hists-muon-hlt-sf-up"][i].Add(results[i]["double_fake_muon_hlt_sf_up"])

                if label == "wg+jets":
                    fake_photon["hists"][i].Add(results[i]["fake_photon_fail_fiducial"])
                    fake_photon["hists-pileup-up"][i].Add(results[i]["fake_photon_fail_fiducial_pileup_up"])
                    fake_photon["hists-prefire-up"][i].Add(results[i]["fake_photon_fail_fiducial_prefire_up"])
                    fake_photon["hists-jes-up"][i].Add(results[i]["fake_photon_fail_fiducial_jes_up"])
                    fake_photon["hists-jer-up"][i].Add(results[i]["fake_photon_fail_fiducial_jer_up"])
                    fake_photon["hists-photon-id-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_photon_id_sf_up"])
                    fake_photon["hists-electron-reco-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_electron_reco_sf_up"])
                    fake_photon["hists-electron-id-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_electron_id_sf_up"])
                    fake_photon["hists-electron-hlt-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_electron_hlt_sf_up"])
                    fake_photon["hists-muon-id-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_muon_id_sf_up"])
                    fake_photon["hists-muon-iso-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_muon_iso_sf_up"])
                    fake_photon["hists-muon-hlt-sf-up"][i].Add(results[i]["fake_photon_fail_fiducial_muon_hlt_sf_up"])
                    for j in range(n_fake_photon_alt):
                        fake_photon["hists-alt"+str(j)][i].Add(results[i]["fake_photon_fail_fiducial_fake_photon_alt"+str(j)])

#                    fake_photon_alt["hists"][i].Add(results[i]["fake_photon_alt"])
#                    fake_photon_stat_up["hists"][i].Add(results[i]["fake_photon_stat_up"])
                    fake_lepton["hists"][i].Add(results[i]["fake_lepton_fail_fiducial"])
                    fake_lepton["hists-pileup-up"][i].Add(results[i]["fake_lepton_fail_fiducial_pileup_up"])
                    fake_lepton["hists-prefire-up"][i].Add(results[i]["fake_lepton_fail_fiducial_prefire_up"])
                    fake_lepton["hists-jes-up"][i].Add(results[i]["fake_lepton_fail_fiducial_jes_up"])
                    fake_lepton["hists-jer-up"][i].Add(results[i]["fake_lepton_fail_fiducial_jer_up"])
                    fake_lepton["hists-photon-id-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_photon_id_sf_up"])
                    fake_lepton["hists-electron-reco-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_electron_reco_sf_up"])
                    fake_lepton["hists-electron-id-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_electron_id_sf_up"])
                    fake_lepton["hists-electron-hlt-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_electron_hlt_sf_up"])
                    fake_lepton["hists-muon-id-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_muon_id_sf_up"])
                    fake_lepton["hists-muon-iso-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_muon_iso_sf_up"])
                    fake_lepton["hists-muon-hlt-sf-up"][i].Add(results[i]["fake_lepton_fail_fiducial_muon_hlt_sf_up"])
#                    fake_lepton_stat_up["hists"][i].Add(results[i]["fake_lepton_stat_up"])
#                    fake_lepton_stat_down["hists"][i].Add(results[i]["fake_lepton_stat_down"])

                    double_fake["hists"][i].Add(results[i]["double_fake_fail_fiducial"])
                    double_fake["hists-pileup-up"][i].Add(results[i]["double_fake_fail_fiducial_pileup_up"])
                    double_fake["hists-prefire-up"][i].Add(results[i]["double_fake_fail_fiducial_prefire_up"])
                    double_fake["hists-jes-up"][i].Add(results[i]["double_fake_fail_fiducial_jes_up"])
                    double_fake["hists-jer-up"][i].Add(results[i]["double_fake_fail_fiducial_jer_up"])
                    double_fake["hists-photon-id-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_photon_id_sf_up"])
                    double_fake["hists-electron-reco-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_electron_reco_sf_up"])
                    double_fake["hists-electron-id-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_electron_id_sf_up"])
                    double_fake["hists-electron-hlt-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_electron_hlt_sf_up"])
                    double_fake["hists-muon-id-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_muon_id_sf_up"])
                    double_fake["hists-muon-iso-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_muon_iso_sf_up"])
                    double_fake["hists-muon-hlt-sf-up"][i].Add(results[i]["double_fake_fail_fiducial_muon_hlt_sf_up"])
                    for j in range(n_fake_photon_alt):
                        double_fake["hists-alt"+str(j)][i].Add(results[i]["double_fake_fail_fiducial_fake_photon_alt"+str(j)])
#                    double_fake_alt["hists"][i].Add(results[i]["double_fake_alt"])
#                    double_fake_stat_up["hists"][i].Add(results[i]["double_fake_stat_up"])

                if sample["e_to_p"] or sample["e_to_p_non_res"]:
                    for j in range(len(etopbinning)): 
                        e_to_p[j]["hists"][i].Add(results[i]["e_to_p"+str(j)])
                        e_to_p[j]["hists-electron-id-sf-up"][i].Add(results[i]["e_to_p_electron_id_sf_up"+str(j)])
                        e_to_p[j]["hists-electron-reco-sf-up"][i].Add(results[i]["e_to_p_electron_reco_sf_up"+str(j)])
                        e_to_p[j]["hists-electron-hlt-sf-up"][i].Add(results[i]["e_to_p_electron_hlt_sf_up"+str(j)])
                        e_to_p[j]["hists-photon-id-sf-up"][i].Add(results[i]["e_to_p_photon_id_sf_up"+str(j)])
                        e_to_p[j]["hists-pileup-up"][i].Add(results[i]["e_to_p_pileup_up"+str(j)])
                        e_to_p[j]["hists-prefire-up"][i].Add(results[i]["e_to_p_prefire_up"+str(j)])
                        e_to_p[j]["hists-jes-up"][i].Add(results[i]["e_to_p_jes_up"+str(j)])
                        e_to_p[j]["hists-jer-up"][i].Add(results[i]["e_to_p_jer_up"+str(j)])
                        e_to_p_total["hists"][i].Add(results[i]["e_to_p"+str(j)])
                        e_to_p_total["hists-electron-id-sf-up"][i].Add(results[i]["e_to_p_electron_id_sf_up"+str(j)])
                        e_to_p_total["hists-electron-reco-sf-up"][i].Add(results[i]["e_to_p_electron_reco_sf_up"+str(j)])
                        e_to_p_total["hists-electron-hlt-sf-up"][i].Add(results[i]["e_to_p_electron_hlt_sf_up"+str(j)])
                        e_to_p_total["hists-photon-id-sf-up"][i].Add(results[i]["e_to_p_photon_id_sf_up"+str(j)])
                        e_to_p_total["hists-pileup-up"][i].Add(results[i]["e_to_p_pileup_up"+str(j)])
                        e_to_p_total["hists-prefire-up"][i].Add(results[i]["e_to_p_prefire_up"+str(j)])
                        e_to_p_total["hists-jes-up"][i].Add(results[i]["e_to_p_jes_up"+str(j)])
                        e_to_p_total["hists-jer-up"][i].Add(results[i]["e_to_p_jer_up"+str(j)])

        for i in range(len(variables)):    

            if labels[label]["color"] == None:
                continue

            if label == "w+jets":
                labels[label]["hists-prompt-pileup"][i].SetFillColor(ROOT.kOrange+3)
                labels[label]["hists-prompt-pileup"][i].SetFillStyle(1001)
                labels[label]["hists-prompt-pileup"][i].SetLineColor(ROOT.kOrange+3)

            if label == "wg+jets":
                labels[label]["hists-pass-fiducial"][i].SetFillColor(labels[label]["color-fid"])
                labels[label]["hists-pass-fiducial"][i].SetFillStyle(1001)
                labels[label]["hists-pass-fiducial"][i].SetLineColor(labels[label]["color-fid"])
                labels[label]["hists-fail-fiducial"][i].SetFillColor(labels[label]["color-non-fid"])
                labels[label]["hists-fail-fiducial"][i].SetFillStyle(1001)
                labels[label]["hists-fail-fiducial"][i].SetLineColor(labels[label]["color-non-fid"])

            labels[label]["hists"][i].SetFillColor(labels[label]["color"])
            labels[label]["hists"][i].SetFillStyle(1001)
            labels[label]["hists"][i].SetLineColor(labels[label]["color"])

if args.ewdim6:

    sm_lhe_weight = 0

#    sm_lhe_weight_hist = ROOT.TH1D('', '', n_photon_pt_bins, binning_photon_pt )

#    sm_hist = ROOT.TH1D('', '', n_photon_pt_bins, binning_photon_pt )

    sm_lhe_weight_hist = histogram_models[ewdim6_index].GetHistogram()

    sm_hist = histogram_models[ewdim6_index].GetHistogram()

    cwww_reweights = [0,0+1,0+2,0+3,0+4,0+5,0+6]

    #cwww_coefficients = [0.0, 10.0,-10.0,20.0,-20.0,-30.0,30.0]

    cwww_coefficients = [0.0, 1.0,-1.0,2.0,-2.0,-3.0,3.0]

    cwww_hists = []

    cw_reweights = [0,0+7,0+8,0+9,0+10,0+11,0+12]

    #cw_coefficients = [0.0, 80.0,-80.0,160.0,-160.0,240.0,-240.0]

#    cw_coefficients = [0.0, 17.0,-17.0,34.0,-34.0,51.0,-51.0]

    cw_coefficients = [0.0, 51.0,-51.0,34.0,-34.0,17.0,-17.0]

    cw_hists = []

    cb_reweights = [0,0+13,0+14,0+15,0+16,0+17,0+18]

    #cb_coefficients = [0.0, 80.0,-80.0,160.0,-160.0,240.0,-240.0]

#    cb_coefficients = [0.0, 17.0,-17.0,34.0,-34.0,51.0,-51.0]

    cb_coefficients = [0.0, 51.0,-51.0,34.0,-34.0,17.0,-17.0]

    cb_hists = []

    cpwww_reweights = [0,0+19,0+20,0+21,0+22,0+23,0+24]

    #cpwww_coefficients = [0.0, 4.0,-4.0,8.0,-8.0,12.0,-12.0]

#    cpwww_coefficients = [0.0, 0.5,-0.5,1.0,-1.0,1.5,-1.5]

    cpwww_coefficients = [0.0, 1.5,-1.5,1.0,-1.0,0.5,-0.5]

    cpwww_hists = []

    cpw_reweights = [0,0+25,0+26,0+27,0+28,0+29,0+30]

    #cpw_coefficients = [0.0, 40.0,-40.0,80.0,-80.0,120.0,-120.0]

    cpw_coefficients = [0.0, 24.0,-24.0,16.0,-16.0,8.0,-8.0]

    cpw_hists = []

    for i in range(0,len(cwww_reweights)):
        cwww_hists.append(ROOT.TH1D('', '', n_photon_pt_bins, binning_photon_pt ))

    for i in range(0,len(cw_reweights)):
        cw_hists.append(ROOT.TH1D('', '', n_photon_pt_bins, binning_photon_pt ))

    for i in range(0,len(cb_reweights)):
        cb_hists.append(ROOT.TH1D('', '', n_photon_pt_bins, binning_photon_pt ))

    for i in range(0,len(cpwww_reweights)):
        cpwww_hists.append(ROOT.TH1D('', '', n_photon_pt_bins, binning_photon_pt ))

    for i in range(0,len(cpw_reweights)):
        cpw_hists.append(ROOT.TH1D('', '', n_photon_pt_bins, binning_photon_pt ))

    gen_matching_string = "(pass_lepton_gen && (photon_gen_matching == 4 || photon_gen_matching == 5 || photon_gen_matching == 6))"

    for year in years:    

        print "Running over sample " + str(labels["wg+jets"]["samples"][year][0]["filename"])

        lumi = dict_lumi[year]

        rdf=ROOT.RDataFrame("Events",labels["wg+jets"]["samples"][year][0]["filename"])

        if "is_lepton_real" in rdf.GetColumnNames():
            rinterface = rdf.Define("pass_lepton_gen","is_lepton_real == 1")
        else:    
            rinterface = rdf.Define("pass_lepton_gen","lepton_gen_matching == 1 || lepton_gen_matching == 15")

        rinterface = rinterface.Filter(get_filter_string(year) + " && " + gen_matching_string)

        rinterface = rinterface.Define("xs_weight",str(labels["wg+jets"]["samples"][year][0]["xs"]*1000*lumi/labels["wg+jets"]["samples"][year][0]["nweightedevents"]) + "*gen_weight/abs(gen_weight)")  

        rinterface = rinterface.Define("weight","xs_weight*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\") : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\"))")    

        for variable_definition in variable_definitions:
            rinterface = rinterface.Define(variable_definition[0],variable_definition[1])

        rresultptr = rinterface.Histo1D(histogram_models[ewdim6_index],variables[ewdim6_index],"weight")

        sm_hist.Add(rresultptr.GetValue())

    sm_hist.Print("all")

    gen_matching_string = "((lepton_gen_matching == 1 || lepton_gen_matching == 15) && (photon_gen_matching == 4 || photon_gen_matching == 5 || photon_gen_matching == 6))"

    for year in years:

        print "Running over sample " + str(ewdim6_samples[year][0]["filename"])

        lumi = dict_lumi[year]

        rdf=ROOT.RDataFrame("Events",ewdim6_samples[year][0]["filename"])

        if "is_lepton_real" in rdf.GetColumnNames():
            rinterface = rdf.Define("pass_lepton_gen","is_lepton_real == 1")
        else:    
            rinterface = rdf.Define("pass_lepton_gen","lepton_gen_matching == 1 || lepton_gen_matching == 15")

        rinterface = rinterface.Filter(get_filter_string(year) + " && " + gen_matching_string)

        rinterface = rinterface.Define("xs_weight",str(ewdim6_samples[year][0]["xs"]*1000*lumi/ewdim6_samples[year][0]["nweightedevents"]) + "*gen_weight/abs(gen_weight)")  

        rinterface = rinterface.Define("weight","xs_weight*photon_efficiency_scale_factor(photon_pt,photon_eta,\""+year+"\")*(abs(lepton_pdg_id) == 13 ? muon_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\") : electron_efficiency_scale_factor(lepton_pt,lepton_eta,\""+year+"\"))")    

        for variable_definition in variable_definitions:
            rinterface = rinterface.Define(variable_definition[0],variable_definition[1])

        rresultptrs_cwww = []
        rresultptrs_cw = []
        rresultptrs_cb = []
        rresultptrs_cpwww = []
        rresultptrs_cpw = []

        for i in range(len(cwww_reweights)):
            rinterface = rinterface.Define("cwww_weight_"+str(i),"weight*LHEReweightingWeight["+str(cwww_reweights[i])+"]")
            rresultptrs_cwww.append(rinterface.Histo1D(histogram_models[ewdim6_index],variables[ewdim6_index],"cwww_weight_"+str(i)))
            
        for i in range(len(cw_reweights)):
            rinterface = rinterface.Define("cw_weight_"+str(i),"weight*LHEReweightingWeight["+str(cw_reweights[i])+"]")
            rresultptrs_cw.append(rinterface.Histo1D(histogram_models[ewdim6_index],variables[ewdim6_index],"cw_weight_"+str(i)))

        for i in range(len(cb_reweights)):
            rinterface = rinterface.Define("cb_weight_"+str(i),"weight*LHEReweightingWeight["+str(cb_reweights[i])+"]")
            rresultptrs_cb.append(rinterface.Histo1D(histogram_models[ewdim6_index],variables[ewdim6_index],"cb_weight_"+str(i)))

        for i in range(len(cpwww_reweights)):
            rinterface = rinterface.Define("cpwww_weight_"+str(i),"weight*LHEReweightingWeight["+str(cpwww_reweights[i])+"]")
            rresultptrs_cpwww.append(rinterface.Histo1D(histogram_models[ewdim6_index],variables[ewdim6_index],"cpwww_weight_"+str(i)))

        for i in range(len(cpw_reweights)):
            rinterface = rinterface.Define("cpw_weight_"+str(i),"weight*LHEReweightingWeight["+str(cpw_reweights[i])+"]")
            rresultptrs_cpw.append(rinterface.Histo1D(histogram_models[ewdim6_index],variables[ewdim6_index],"cpw_weight_"+str(i)))

        rinterface = rinterface.Define("sm_weight","weight*LHEReweightingWeight["+str(sm_lhe_weight)+"]")
        rresultptr_sm = rinterface.Histo1D(histogram_models[ewdim6_index],variables[ewdim6_index],"sm_weight")



        for i in range(len(cwww_reweights)):
            cwww_hists[i].Add(rresultptrs_cwww[i].GetValue())

        for i in range(len(cw_reweights)):
            cw_hists[i].Add(rresultptrs_cw[i].GetValue())

        for i in range(len(cb_reweights)):
            cb_hists[i].Add(rresultptrs_cb[i].GetValue())

        for i in range(len(cpwww_reweights)):
            cpwww_hists[i].Add(rresultptrs_cpwww[i].GetValue())

        for i in range(len(cpw_reweights)):
            cpw_hists[i].Add(rresultptrs_cpw[i].GetValue())

        sm_lhe_weight_hist.Add(rresultptr_sm.GetValue())

    cwww_scaling_outfile = ROOT.TFile("cwww_scaling.root",'recreate')
    cw_scaling_outfile = ROOT.TFile("cw_scaling.root",'recreate')
    cb_scaling_outfile = ROOT.TFile("cb_scaling.root",'recreate')
    cpwww_scaling_outfile = ROOT.TFile("cpwww_scaling.root",'recreate')
    cpw_scaling_outfile = ROOT.TFile("cpw_scaling.root",'recreate')

    cwww_hist_max = max(cwww_coefficients) + (max(cwww_coefficients) - min(cwww_coefficients))/(len(cwww_coefficients)-1)/2
    cwww_hist_min = min(cwww_coefficients) - (max(cwww_coefficients) - min(cwww_coefficients))/(len(cwww_coefficients)-1)/2

    cw_hist_max = max(cw_coefficients) + (max(cw_coefficients) - min(cw_coefficients))/(len(cw_coefficients)-1)/2
    cw_hist_min = min(cw_coefficients) - (max(cw_coefficients) - min(cw_coefficients))/(len(cw_coefficients)-1)/2

    cb_hist_max = max(cb_coefficients) + (max(cb_coefficients) - min(cb_coefficients))/(len(cb_coefficients)-1)/2
    cb_hist_min = min(cb_coefficients) - (max(cb_coefficients) - min(cb_coefficients))/(len(cb_coefficients)-1)/2

    cpwww_hist_max = max(cpwww_coefficients) + (max(cpwww_coefficients) - min(cpwww_coefficients))/(len(cpwww_coefficients)-1)/2
    cpwww_hist_min = min(cpwww_coefficients) - (max(cpwww_coefficients) - min(cpwww_coefficients))/(len(cpwww_coefficients)-1)/2

    cpw_hist_max = max(cpw_coefficients) + (max(cpw_coefficients) - min(cpw_coefficients))/(len(cpw_coefficients)-1)/2
    cpw_hist_min = min(cpw_coefficients) - (max(cpw_coefficients) - min(cpw_coefficients))/(len(cpw_coefficients)-1)/2

    sm_lhe_weight_hist.Print("all")

    c = ROOT.TCanvas("c", "c",5,50,500,500)
    sm_hist.SetLineColor(ROOT.kRed)
    sm_lhe_weight_hist.SetLineColor(ROOT.kBlue)
    sm_hist.SetLineWidth(2)
    sm_lhe_weight_hist.SetLineWidth(2)
    sm_hist.SetMaximum(1.55*max(sm_hist.GetMaximum(),sm_lhe_weight_hist.GetMaximum()))
    sm_hist.Draw()
    sm_lhe_weight_hist.Draw("same")
    s=str(totallumi)+" fb^{-1} (13 TeV)"
    lumilabel = ROOT.TLatex (0.95, 0.93, s)
    lumilabel.SetNDC ()
    lumilabel.SetTextAlign (30)
    lumilabel.SetTextFont (42)
    lumilabel.SetTextSize (0.040)
    lumilabel.Draw("same")
    set_axis_fonts(sm_hist,"x",getXaxisLabel(variables[ewdim6_index]))
    j=0
    draw_legend(xpositions[j]-0.05,0.84 - ypositions[j]*yoffset,sm_hist,"SM unweighted","l")
    j=j+1
    draw_legend(xpositions[j]-0.05,0.84 - ypositions[j]*yoffset,sm_lhe_weight_hist,"SM reweighted","l")
    c.SaveAs(args.outputdir + "/" + "sm_reweighting.png")

    cwww_scaling_hists = {}
    cw_scaling_hists = {}
    cb_scaling_hists = {}
    cpw_scaling_hists = {}
    cpwww_scaling_hists = {}

    for i in range(1,cwww_hists[0].GetNbinsX()+1):
        ROOT.gROOT.cd() #so that the histogram created in the next line is not put in a file that is closed
        cwww_scaling_hists[i]=ROOT.TH1D("ewdim6_scaling_bin_"+str(i),"ewdim6_scaling_bin_"+str(i),len(cwww_coefficients),cwww_hist_min,cwww_hist_max)

        for j in range(0,len(cwww_hists)):
            assert(sm_lhe_weight_hist.GetBinContent(i) > 0)

            cwww_scaling_hists[i].SetBinContent(cwww_scaling_hists[i].GetXaxis().FindFixBin(cwww_coefficients[j]), cwww_hists[j].GetBinContent(i)/sm_lhe_weight_hist.GetBinContent(i))
        
        cwww_scaling_outfile.cd()
        cwww_scaling_hists[i].Write()

    cwww_scaling_outfile.Close()

    for i in range(1,cw_hists[0].GetNbinsX()+1):
        cw_scaling_hists[i]=ROOT.TH1D("ewdim6_scaling_bin_"+str(i),"ewdim6_scaling_bin_"+str(i),len(cw_coefficients),cw_hist_min,cw_hist_max)

        for j in range(0,len(cw_hists)):
            assert(sm_lhe_weight_hist.GetBinContent(i) > 0)

            cw_scaling_hists[i].SetBinContent(cw_scaling_hists[i].GetXaxis().FindFixBin(cw_coefficients[j]), cw_hists[j].GetBinContent(i)/sm_lhe_weight_hist.GetBinContent(i))
            
        cw_scaling_outfile.cd()
        cw_scaling_hists[i].Write()

    cw_scaling_outfile.Close()

    for i in range(1,cb_hists[0].GetNbinsX()+1):
        ROOT.gROOT.cd() #so that the histogram created in the next line is not put in a file that is closed
        cb_scaling_hists[i]=ROOT.TH1D("ewdim6_scaling_bin_"+str(i),"ewdim6_scaling_bin_"+str(i),len(cb_coefficients),cb_hist_min,cb_hist_max);

        for j in range(0,len(cb_hists)):
            assert(sm_lhe_weight_hist.GetBinContent(i) > 0)

            cb_scaling_hists[i].SetBinContent(cb_scaling_hists[i].GetXaxis().FindFixBin(cb_coefficients[j]), cb_hists[j].GetBinContent(i)/sm_lhe_weight_hist.GetBinContent(i))
        
        cb_scaling_outfile.cd()
        cb_scaling_hists[i].Write()

    cb_scaling_outfile.Close()

    for i in range(1,cpwww_hists[0].GetNbinsX()+1):
        cpwww_scaling_hists[i]=ROOT.TH1D("ewdim6_scaling_bin_"+str(i),"ewdim6_scaling_bin_"+str(i),len(cpwww_coefficients),cpwww_hist_min,cpwww_hist_max);

        for j in range(0,len(cpwww_hists)):
            assert(sm_lhe_weight_hist.GetBinContent(i) > 0)

            cpwww_scaling_hists[i].SetBinContent(cpwww_scaling_hists[i].GetXaxis().FindFixBin(cpwww_coefficients[j]), cpwww_hists[j].GetBinContent(i)/sm_lhe_weight_hist.GetBinContent(i))
        
        cpwww_scaling_outfile.cd()
        cpwww_scaling_hists[i].Write()

    cpwww_scaling_outfile.Close()

    for i in range(1,cpw_hists[0].GetNbinsX()+1):
        cpw_scaling_hists[i]=ROOT.TH1D("ewdim6_scaling_bin_"+str(i),"ewdim6_scaling_bin_"+str(i),len(cpw_coefficients),cpw_hist_min,cpw_hist_max);

        for j in range(0,len(cpw_hists)):
            assert(sm_lhe_weight_hist.GetBinContent(i) > 0)

            cpw_scaling_hists[i].SetBinContent(cpw_scaling_hists[i].GetXaxis().FindFixBin(cpw_coefficients[j]), cpw_hists[j].GetBinContent(i)/sm_lhe_weight_hist.GetBinContent(i))
        
        cpw_scaling_outfile.cd()
        cpw_scaling_hists[i].Write()

    cpw_scaling_outfile.Close()

data_mlg_tree = ROOT.TTree()

array_data_mlg=array('f',[0])

data_mlg_tree.Branch('m',array_data_mlg,'m/F')

for year in years:

    if year == "2016":
        lumi=35.9
    elif year == "2017":
        lumi=41.5
    elif year == "2018":
        lumi=59.6
    else:
        assert(0)

    if lepton_name == "muon":
        data_filename = args.workdir+"/data/wg/"+year+"/1June2019/single_muon.root"

    elif lepton_name == "electron":
        if year != "2018":
            data_filename = args.workdir+"/data/wg/"+year+"/1June2019/single_electron.root"
        else:    
            data_filename = args.workdir+"/data/wg/"+year+"/1June2019/egamma.root"

    elif lepton_name == "both":
        if year != "2018":
            data_filename = args.workdir+"/data/wg/"+year+"/1June2019/data.root"
        else:
            data_filename = args.workdir+"/data/wg/"+year+"/1June2019/data.root"
    else:
        assert(0)

    if year == "2016":
        sieie_cut_barrel = sieie_cut_2016_barrel
        sieie_cut_endcap = sieie_cut_2016_endcap
        chiso_cut_barrel = chiso_cut_2016_barrel
        chiso_cut_endcap = chiso_cut_2016_endcap
    elif year == "2017":
        sieie_cut_barrel = sieie_cut_2017_barrel
        sieie_cut_endcap = sieie_cut_2017_endcap
        chiso_cut_barrel = chiso_cut_2017_barrel
        chiso_cut_endcap = chiso_cut_2017_endcap
    elif year == "2018":
        sieie_cut_barrel = sieie_cut_2018_barrel
        sieie_cut_endcap = sieie_cut_2018_endcap
        chiso_cut_barrel = chiso_cut_2018_barrel
        chiso_cut_endcap = chiso_cut_2018_endcap
    else:
        assert(0)

    fake_photon_sieie_cut_barrel = sieie_cut_barrel*1.75
    fake_photon_sieie_cut_endcap = sieie_cut_endcap*1.75
    fake_photon_chiso_cut_barrel = chiso_cut_barrel*1000
    fake_photon_chiso_cut_endcap = chiso_cut_endcap*1000    

    print "Running over "+year+" data"

    rdf=ROOT.RDataFrame("Events",data_filename)

    rinterface = rdf.Filter(get_filter_string(year))

    fake_photon_sieie_cut_cutstring = "((abs(photon_eta) < 1.5 && photon_sieie < "+str(fake_photon_sieie_cut_barrel)+ ") || (abs(photon_eta) > 1.5 && photon_sieie < "+str(fake_photon_sieie_cut_endcap)+ "))" 

    fake_photon_chiso_cut_cutstring = "((abs(photon_eta) < 1.5 && photon_pfRelIso03_chg*photon_pt < "+str(fake_photon_chiso_cut_barrel)+ ") || (abs(photon_eta) > 1.5 && photon_pfRelIso03_chg*photon_pt < "+str(fake_photon_chiso_cut_endcap)+ "))" 

    rinterface = rinterface.Define("weight","photon_selection == 0 && is_lepton_tight == 1")
    rinterface = rinterface.Define("fake_lepton_weight","photon_selection == 0 && is_lepton_tight == 0 ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id) : 0")
    rinterface = rinterface.Define("fake_lepton_stat_up_weight","photon_selection == 0 && is_lepton_tight == 0 ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id,\"up\") : 0")
    rinterface = rinterface.Define("fake_lepton_stat_down_weight","photon_selection == 0 && is_lepton_tight == 0 ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id,\"down\") : 0")
    rinterface = rinterface.Define("fake_photon_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id) : 0")

    for i in range(n_fake_photon_alt):
        rinterface = rinterface.Define("fake_photon_alt"+str(i)+"_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"alt"+str(i)+"\") : 0")
    rinterface = rinterface.Define("fake_photon_stat_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 1 ? get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"stat_up\") : 0")
    rinterface = rinterface.Define("double_fake_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id) : 0") 
    for i in range(n_fake_photon_alt):
        rinterface = rinterface.Define("double_fake_alt"+str(i)+"_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"alt"+str(i)+"\") : 0")
    rinterface = rinterface.Define("double_fake_stat_up_weight","photon_selection == 4 && "+fake_photon_sieie_cut_cutstring + " && " + fake_photon_chiso_cut_cutstring+" && is_lepton_tight == 0 ? get_fake_lepton_weight(lepton_eta,lepton_pt,\""+year+"\",lepton_pdg_id)*get_fake_photon_weight(photon_eta,photon_pt,\""+year+"\",lepton_pdg_id,\"stat_up\") : 0") 

    for variable_definition in variable_definitions:
            rinterface = rinterface.Define(variable_definition[0],variable_definition[1])

    rresultptrs = []
    rresultptrs_fake_photon = []
    rresultptrs_fake_photon_stat_up = []
    rresultptrs_fake_photon_alt = []
    for i in range(n_fake_photon_alt):
        rresultptrs_fake_photon_alt.append([])
    rresultptrs_fake_lepton = []
    rresultptrs_fake_lepton_stat_up = []
    rresultptrs_fake_lepton_stat_down = []
    rresultptrs_double_fake = []
    rresultptrs_double_fake_alt = []
    for i in range(n_fake_photon_alt):
        rresultptrs_double_fake_alt.append([])
    rresultptrs_double_fake_stat_up = []

    for i in range(len(variables)):
        rresultptrs_fake_photon.append(rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_weight"))
        rresultptrs.append(rinterface.Histo1D(histogram_models[i],variables[i],"weight"))
        for j in range(n_fake_photon_alt):
            rresultptrs_fake_photon_alt[j].append(rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_alt"+str(j)+"_weight"))
        rresultptrs_fake_photon_stat_up.append(rinterface.Histo1D(histogram_models[i],variables[i],"fake_photon_stat_up_weight"))
        rresultptrs_fake_lepton.append(rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_weight"))
        rresultptrs_fake_lepton_stat_up.append(rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_stat_up_weight"))
        rresultptrs_fake_lepton_stat_down.append(rinterface.Histo1D(histogram_models[i],variables[i],"fake_lepton_stat_down_weight"))
        rresultptrs_double_fake.append(rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_weight"))
        for j in range(n_fake_photon_alt):
            rresultptrs_double_fake_alt[j].append(rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_alt"+str(j)+"_weight"))
        rresultptrs_double_fake_stat_up.append(rinterface.Histo1D(histogram_models[i],variables[i],"double_fake_stat_up_weight"))

    for i in range(len(variables)):
        data["hists"][i].Add(rresultptrs[i].GetValue())
        if year == "2016":    
            fake_photon_2016["hists"][i].Add(rresultptrs_fake_photon[i].GetValue())

        fake_photon["hists"][i].Add(rresultptrs_fake_photon[i].GetValue())
        fake_photon["hists-pileup-up"][i].Add(rresultptrs_fake_photon[i].GetValue())
        fake_photon["hists-prefire-up"][i].Add(rresultptrs_fake_photon[i].GetValue())
        fake_photon["hists-jer-up"][i].Add(rresultptrs_fake_photon[i].GetValue())
        fake_photon["hists-jes-up"][i].Add(rresultptrs_fake_photon[i].GetValue())
        fake_photon["hists-photon-id-sf-up"][i].Add(rresultptrs_fake_photon[i].GetValue())
        fake_photon["hists-electron-reco-sf-up"][i].Add(rresultptrs_fake_photon[i].GetValue())
        fake_photon["hists-electron-id-sf-up"][i].Add(rresultptrs_fake_photon[i].GetValue())
        fake_photon["hists-electron-hlt-sf-up"][i].Add(rresultptrs_fake_photon[i].GetValue())
        fake_photon["hists-muon-id-sf-up"][i].Add(rresultptrs_fake_photon[i].GetValue())
        fake_photon["hists-muon-iso-sf-up"][i].Add(rresultptrs_fake_photon[i].GetValue())
        fake_photon["hists-muon-hlt-sf-up"][i].Add(rresultptrs_fake_photon[i].GetValue())
        for j in range(n_fake_photon_alt):
            fake_photon["hists-alt"+str(j)][i].Add(rresultptrs_fake_photon_alt[j][i].GetValue())
        fake_photon_stat_up["hists"][i].Add(rresultptrs_fake_photon_stat_up[i].GetValue())

        fake_lepton["hists"][i].Add(rresultptrs_fake_lepton[i].GetValue())
        fake_lepton["hists-pileup-up"][i].Add(rresultptrs_fake_lepton[i].GetValue())
        fake_lepton["hists-prefire-up"][i].Add(rresultptrs_fake_lepton[i].GetValue())
        fake_lepton["hists-jer-up"][i].Add(rresultptrs_fake_lepton[i].GetValue())
        fake_lepton["hists-jes-up"][i].Add(rresultptrs_fake_lepton[i].GetValue())
        fake_lepton["hists-photon-id-sf-up"][i].Add(rresultptrs_fake_lepton[i].GetValue())
        fake_lepton["hists-electron-reco-sf-up"][i].Add(rresultptrs_fake_lepton[i].GetValue())
        fake_lepton["hists-electron-id-sf-up"][i].Add(rresultptrs_fake_lepton[i].GetValue())
        fake_lepton["hists-electron-hlt-sf-up"][i].Add(rresultptrs_fake_lepton[i].GetValue())
        fake_lepton["hists-muon-id-sf-up"][i].Add(rresultptrs_fake_lepton[i].GetValue())
        fake_lepton["hists-muon-iso-sf-up"][i].Add(rresultptrs_fake_lepton[i].GetValue())
        fake_lepton["hists-muon-hlt-sf-up"][i].Add(rresultptrs_fake_lepton[i].GetValue())
        fake_lepton_stat_up["hists"][i].Add(rresultptrs_fake_lepton_stat_up[i].GetValue())
        fake_lepton_stat_down["hists"][i].Add(rresultptrs_fake_lepton_stat_down[i].GetValue())

        double_fake["hists"][i].Add(rresultptrs_double_fake[i].GetValue())
        double_fake["hists-pileup-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        double_fake["hists-prefire-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        double_fake["hists-jer-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        double_fake["hists-jes-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        double_fake["hists-photon-id-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        double_fake["hists-electron-reco-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        double_fake["hists-electron-id-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        double_fake["hists-electron-hlt-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        double_fake["hists-muon-id-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        double_fake["hists-muon-iso-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        double_fake["hists-muon-hlt-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        for j in range(n_fake_photon_alt):
            double_fake["hists-alt"+str(j)][i].Add(rresultptrs_double_fake_alt[j][i].GetValue())
        double_fake_stat_up["hists"][i].Add(rresultptrs_double_fake_stat_up[i].GetValue())
        rresultptrs_double_fake[i].GetPtr().Scale(-1)
        for j in range(n_fake_photon_alt):
            rresultptrs_double_fake_alt[j][i].GetPtr().Scale(-1)
        if year == "2016":    
            fake_photon_2016["hists"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_photon["hists"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_photon["hists-pileup-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_photon["hists-prefire-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_photon["hists-jer-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_photon["hists-jes-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_photon["hists-photon-id-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_photon["hists-electron-reco-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_photon["hists-electron-id-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_photon["hists-electron-hlt-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_photon["hists-muon-id-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_photon["hists-muon-iso-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_photon["hists-muon-hlt-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())

        for j in range(n_fake_photon_alt):
            fake_photon["hists-alt"+str(j)][i].Add(rresultptrs_double_fake_alt[j][i].GetValue())

        fake_lepton["hists"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_lepton["hists-pileup-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_lepton["hists-prefire-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_lepton["hists-jer-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_lepton["hists-jes-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_lepton["hists-photon-id-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_lepton["hists-electron-reco-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_lepton["hists-electron-id-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_lepton["hists-electron-hlt-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_lepton["hists-muon-id-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_lepton["hists-muon-iso-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())
        fake_lepton["hists-muon-hlt-sf-up"][i].Add(rresultptrs_double_fake[i].GetValue())

if args.no_wjets_for_2017_and_2018 and "w+jets" in labels:
    for i in range(len(variables)):
        labels["w+jets"]["hists"][i].Scale(fake_photon["hists"][i].Integral()/fake_photon_2016["hists"][i].Integral())
#        wjets_2016["hists"][i].Scale(fake_photon["hists"][i].Integral()/fake_photon_2016["hists"][i].Integral())
#        wjets_fake_photon_2016["hists"][i].Scale(fake_photon["hists"][i].Integral()/fake_photon_2016["hists"][i].Integral())

#wjets_2016["hists"][3].Scale(fake_photon["hists"][mlg_index].Integral()/fake_photon_2016["hists"][mlg_index].Integral())
#wjets_fake_photon["hists"][3].Scale(fake_photon["hists"][mlg_index].Integral()/fake_photon_2016["hists"][mlg_index].Integral())
non_closure = []
for i in range(len(variables)):

    c = ROOT.TCanvas("c", "c",5,50,500,500)

    s="35.9 fb^{-1} (13 TeV)"
    lumilabel = ROOT.TLatex (0.95, 0.93, s)
    lumilabel.SetNDC ()
    lumilabel.SetTextAlign (30)
    lumilabel.SetTextFont (42)
    lumilabel.SetTextSize (0.040)

    minus_one_hist = histogram_models[i].GetHistogram()
    for j in range(1,minus_one_hist.GetNbinsX()+1):
        minus_one_hist.SetBinContent(j,-1)
        minus_one_hist.SetBinError(j,0)
    non_closure.append(histogram_models[i].GetHistogram())


    wjets_2016["hists"][i].SetLineColor(ROOT.kBlue)
    wjets_fake_photon_2016["hists"][i].SetLineColor(ROOT.kRed)
    wjets_fake_photon_2016["hists"][i].SetMinimum(0)
#    if wjets_fake_photon_2016["hists"][i].Integral() > 0:
#        wjets_fake_photon_2016["hists"][i].Scale(wjets_2016["hists"][i].Integral()/wjets_fake_photon_2016["hists"][i].Integral())
    wjets_fake_photon_2016["hists"][i].SetMaximum(1.55*max(wjets_fake_photon_2016["hists"][i].GetMaximum(),wjets_2016["hists"][i].GetMaximum()))
    set_axis_fonts(wjets_fake_photon_2016["hists"][i],"x",getXaxisLabel(variables[i]))
    wjets_fake_photon_2016["hists"][i].Draw()
    wjets_2016["hists"][i].Draw("same")
    lumilabel.Draw("same")


    set_axis_fonts(wjets_fake_photon_2016["hists"][i],"x",getXaxisLabel(variables[i]))
    j=0
    draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,wjets_fake_photon_2016["hists"][i],"fake photon","l")
    j=j+1
    draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,wjets_2016["hists"][i],"w+jets","l")
    j=j+1
    c.SaveAs(args.outputdir + "/" + "closure_test_"+variables_labels[i]+".png")

    wjets_fake_photon_chiso_2016["hists"][i].SetLineColor(ROOT.kRed)
#    if wjets_fake_photon_chiso_2016["hists"][i].Integral() > 0:
#        wjets_fake_photon_chiso_2016["hists"][i].Scale(wjets_2016["hists"][i].Integral()/wjets_fake_photon_chiso_2016["hists"][i].Integral())
    wjets_fake_photon_chiso_2016["hists"][i].SetMaximum(1.55*max(wjets_fake_photon_chiso_2016["hists"][i].GetMaximum(),wjets_2016["hists"][i].GetMaximum()))
    wjets_fake_photon_chiso_2016["hists"][i].SetMinimum(0)
    wjets_fake_photon_chiso_2016["hists"][i].Draw()
    wjets_2016["hists"][i].Draw("same")
    lumilabel.Draw("same")
    set_axis_fonts(wjets_fake_photon_chiso_2016["hists"][i],"x",getXaxisLabel(variables[i]))
    j=0
    draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,wjets_fake_photon_chiso_2016["hists"][i],"fake photon","l")
    j=j+1
    draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,wjets_2016["hists"][i],"w+jets","l")
    j=j+1
    c.SaveAs(args.outputdir + "/" + "closure_test_chiso_"+variables_labels[i]+".png")
    
    window_size=5

    if i == mlg_index:
        for j in range(1,wjets_2016["hists"][i].GetNbinsX()+1):
            rebinnum=ROOT.TH1F("","",2*window_size+1,0,2*window_size+1)
            rebinden=ROOT.TH1F("","",2*window_size+1,0,2*window_size+1)

            rebinnum.Sumw2()
            rebinden.Sumw2()

            if j < window_size+1:

                for k in range(1,2*window_size+2):
                    rebinnum.SetBinContent(k,wjets_2016["hists"][i].GetBinContent(k))
                    rebinnum.SetBinError(k,wjets_2016["hists"][i].GetBinError(k))
                    rebinden.SetBinContent(k,wjets_fake_photon_2016["hists"][i].GetBinContent(k))
                    rebinden.SetBinError(k,wjets_fake_photon_2016["hists"][i].GetBinError(k))
            elif j > wjets_2016["hists"][i].GetNbinsX() - window_size-1:        

                for k in range(wjets_2016["hists"][i].GetNbinsX()+1-2*window_size-1,wjets_2016["hists"][i].GetNbinsX()+1):
                    rebinnum.SetBinContent(k-wjets_2016["hists"][i].GetNbinsX()-1+2*window_size+1+1,wjets_2016["hists"][i].GetBinContent(k))
                    rebinnum.SetBinError(k-wjets_2016["hists"][i].GetNbinsX()-1+2*window_size+1+1,wjets_2016["hists"][i].GetBinError(k))
                    rebinden.SetBinContent(k-wjets_2016["hists"][i].GetNbinsX()-1+2*window_size+1+1,wjets_fake_photon_2016["hists"][i].GetBinContent(k))
                    rebinden.SetBinError(k-wjets_2016["hists"][i].GetNbinsX()-1+2*window_size+1+1,wjets_fake_photon_2016["hists"][i].GetBinError(k))
            else:

                for k in range(j-window_size,j+window_size+1):
                    rebinnum.SetBinContent(k-j+window_size+1,wjets_2016["hists"][i].GetBinContent(k))
                    rebinnum.SetBinError(k-j+window_size+1,wjets_2016["hists"][i].GetBinError(k))
                    rebinden.SetBinContent(k-j+window_size+1,wjets_fake_photon_2016["hists"][i].GetBinContent(k))
                    rebinden.SetBinError(k-j+window_size+1,wjets_fake_photon_2016["hists"][i].GetBinError(k))
                 
#            rebinnum.Print("all")        
#            rebinden.Print("all")        

            rebinnum.Rebin(2*window_size+1)        
            rebinden.Rebin(2*window_size+1)        

            rebinnum.Divide(rebinden)

#            rebinnum.Print("all")        
#            rebinden.Print("all") 

            non_closure[len(non_closure)-1].SetBinContent(j,rebinnum.GetBinContent(1))   
            non_closure[len(non_closure)-1].SetBinError(j,rebinnum.GetBinError(1))   

    else:
        non_closure[len(non_closure)-1].Add(wjets_2016["hists"][i])
#    wjets_fake_photon_2016["hists"][i].Scale(-1)
#    non_closure[len(non_closure)-1].Add(wjets_fake_photon_2016["hists"][i])
#    wjets_fake_photon["hists"][i].Scale(-1)
        non_closure[len(non_closure)-1].Divide(wjets_fake_photon_2016["hists"][i])

    non_closure[len(non_closure)-1].Add(minus_one_hist)
    non_closure[len(non_closure)-1].SetMinimum(-1.5)
    non_closure[len(non_closure)-1].SetMaximum(1.5)
    non_closure[len(non_closure)-1].Draw()
    set_axis_fonts(non_closure[len(non_closure)-1],"x",getXaxisLabel(variables[i]))
    c.SaveAs(args.outputdir + "/" + "non_closure_"+variables_labels[i]+".png")

for i in range(len(variables)):

    c = ROOT.TCanvas("c", "c",5,50,500,500)

    s=str(totallumi)+" fb^{-1} (13 TeV)"
    lumilabel = ROOT.TLatex (0.95, 0.93, s)
    lumilabel.SetNDC ()
    lumilabel.SetTextAlign (30)
    lumilabel.SetTextFont (42)
    lumilabel.SetTextSize (0.040)

    central = labels["wg+jets"]["hists"][i].Clone()
    up = labels["wg+jets"]["hists-pileup-up"][i].Clone()
    down = makeDownShape(labels["wg+jets"]["hists-pileup-up"][i],labels["wg+jets"]["hists"][i]).Clone()

    central.SetLineColor(ROOT.kBlack)
    central.SetLineWidth(2)

    up.SetLineColor(ROOT.kRed)
    up.SetLineWidth(2)

    down.SetLineColor(ROOT.kBlue)
    down.SetLineWidth(2)

    central.Draw("hist")
    up.Draw("hist same")
    down.Draw("hist same")

    lumilabel.Draw("same")

    c.SaveAs(args.outputdir + "/" + "wgjets_pileup_unc_"+variables_labels[i]+".png")

    c.Close()

if "wg+jets" in labels:
    prefire_unc = abs(labels["wg+jets"]["hists-prefire-up"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].Integral())
    pileup_unc = abs(labels["wg+jets"]["hists-pileup-up"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].Integral())
    jes_unc = abs(labels["wg+jets"]["hists-jes-up"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].Integral())
    jer_unc = abs(labels["wg+jets"]["hists-jer-up"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].Integral())
    electron_id_sf_unc = labels["wg+jets"]["hists-electron-id-sf-up"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].Integral()
    electron_reco_sf_unc = labels["wg+jets"]["hists-electron-reco-sf-up"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].Integral()
    electron_hlt_sf_unc = labels["wg+jets"]["hists-electron-hlt-sf-up"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].Integral()
    muon_id_sf_unc = labels["wg+jets"]["hists-muon-id-sf-up"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].Integral()
    muon_iso_sf_unc = labels["wg+jets"]["hists-muon-iso-sf-up"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].Integral()
    muon_hlt_sf_unc = labels["wg+jets"]["hists-muon-hlt-sf-up"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].Integral()
    photon_id_sf_unc = labels["wg+jets"]["hists-photon-id-sf-up"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].Integral()

#    print labels["wg+jets"]["hists-muon-iso-sf-up"][mlg_index].Integral()
#    print labels["wg+jets"]["hists"][mlg_index].Integral()

#    print "(number of wg+jets events run over) = "+str(labels["wg+jets"]["samples"][year][0]["nweightedevents"])

#    print "fiducial_region_cuts_efficiency = "+str(fiducial_region_cuts_efficiency)

if args.draw_ewdim6:
    for i in range(1,n_photon_pt_bins+1):
        #hardcoded to use bin 6 of the scaling histogram for now 
        ewdim6["hists"][ewdim6_index].SetBinContent(i,cb_scaling_hists[i].GetBinContent(7)*labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i))

for i in range(len(variables)):
    if args.float_sig_fake_cont:
        labels["wg+jets"]["hists-pass-fiducial"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial"][i])
        for j in range(n_fake_photon_alt):
            labels["wg+jets"]["hists-pass-fiducial-fake-photon-alt"+str(j)][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-fake-photon-alt"+str(j)][i])
        labels["wg+jets"]["hists-pass-fiducial-pileup-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-pileup-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-prefire-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-prefire-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-jes-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-jes-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-jer-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-jer-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-photon-id-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-photon-id-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-electron-reco-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-electron-reco-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-electron-id-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-electron-id-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-electron-hlt-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-electron-hlt-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-muon-id-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-muon-id-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-muon-iso-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-muon-iso-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-muon-hlt-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-muon-hlt-sf-up"][i])

        labels["wg+jets"]["hists-pass-fiducial"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial"][i])
        labels["wg+jets"]["hists-pass-fiducial-pileup-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-pileup-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-prefire-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-prefire-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-jes-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-jes-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-jer-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-jer-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-photon-id-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-photon-id-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-electron-reco-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-electron-reco-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-electron-id-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-electron-id-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-electron-hlt-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-electron-hlt-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-muon-id-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-muon-id-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-muon-iso-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-muon-iso-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-muon-hlt-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-muon-hlt-sf-up"][i])

        labels["wg+jets"]["hists-pass-fiducial"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial"][i])
        for j in range(n_fake_photon_alt):
            labels["wg+jets"]["hists-pass-fiducial-fake-photon-alt"+str(j)][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-fake-photon-alt"+str(j)][i])
        labels["wg+jets"]["hists-pass-fiducial-pileup-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-pileup-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-prefire-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-prefire-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-jes-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-jes-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-jer-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-jer-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-photon-id-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-photon-id-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-electron-reco-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-electron-reco-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-electron-id-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-electron-id-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-electron-hlt-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-electron-hlt-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-muon-id-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-muon-id-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-muon-iso-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-muon-iso-sf-up"][i])
        labels["wg+jets"]["hists-pass-fiducial-muon-hlt-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-muon-hlt-sf-up"][i])
        pass
    else:
        fake_photon["hists"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial"][i])
        for j in range(n_fake_photon_alt):
            fake_photon["hists-alt"+str(j)][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-fake-photon-alt"+str(j)][i])
        fake_photon["hists-pileup-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-pileup-up"][i])
        fake_photon["hists-prefire-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-prefire-up"][i])
        fake_photon["hists-jes-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-jes-up"][i])
        fake_photon["hists-jer-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-jer-up"][i])
        fake_photon["hists-photon-id-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-photon-id-sf-up"][i])
        fake_photon["hists-electron-reco-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-electron-reco-sf-up"][i])
        fake_photon["hists-electron-id-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-electron-id-sf-up"][i])
        fake_photon["hists-electron-hlt-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-electron-hlt-sf-up"][i])
        fake_photon["hists-muon-id-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-muon-id-sf-up"][i])
        fake_photon["hists-muon-iso-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-muon-iso-sf-up"][i])
        fake_photon["hists-muon-hlt-sf-up"][i].Add(labels["wg+jets"]["hists-fake-photon-pass-fiducial-muon-hlt-sf-up"][i])

        fake_lepton["hists"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial"][i])
        fake_lepton["hists-pileup-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-pileup-up"][i])
        fake_lepton["hists-prefire-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-prefire-up"][i])
        fake_lepton["hists-jes-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-jes-up"][i])
        fake_lepton["hists-jer-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-jer-up"][i])
        fake_lepton["hists-photon-id-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-photon-id-sf-up"][i])
        fake_lepton["hists-electron-reco-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-electron-reco-sf-up"][i])
        fake_lepton["hists-electron-id-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-electron-id-sf-up"][i])
        fake_lepton["hists-electron-hlt-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-electron-hlt-sf-up"][i])
        fake_lepton["hists-muon-id-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-muon-id-sf-up"][i])
        fake_lepton["hists-muon-iso-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-muon-iso-sf-up"][i])
        fake_lepton["hists-muon-hlt-sf-up"][i].Add(labels["wg+jets"]["hists-fake-lepton-pass-fiducial-muon-hlt-sf-up"][i])

        double_fake["hists"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial"][i])
        for j in range(n_fake_photon_alt):
            double_fake["hists-alt"+str(j)][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-fake-photon-alt"+str(j)][i])
        double_fake["hists-pileup-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-pileup-up"][i])
        double_fake["hists-prefire-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-prefire-up"][i])
        double_fake["hists-jes-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-jes-up"][i])
        double_fake["hists-jer-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-jer-up"][i])
        double_fake["hists-photon-id-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-photon-id-sf-up"][i])
        double_fake["hists-electron-reco-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-electron-reco-sf-up"][i])
        double_fake["hists-electron-id-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-electron-id-sf-up"][i])
        double_fake["hists-electron-hlt-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-electron-hlt-sf-up"][i])
        double_fake["hists-muon-id-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-muon-id-sf-up"][i])
        double_fake["hists-muon-iso-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-muon-iso-sf-up"][i])
        double_fake["hists-muon-hlt-sf-up"][i].Add(labels["wg+jets"]["hists-double-fake-pass-fiducial-muon-hlt-sf-up"][i])
        pass

fake_photon_syst2_smoothed = []
non_closure_smoothed = []

for i in range(len(variables)):
    window_size = max(int(wjets_2016["hists"][i].GetNbinsX()/5),1)

    fake_photon_syst2_smoothed.append(histogram_models[i].GetHistogram())
    non_closure_smoothed.append(histogram_models[i].GetHistogram())
    num_smoothed = histogram_models[i].GetHistogram()
    den_smoothed = histogram_models[i].GetHistogram()

    for j in range(1,wjets_2016["hists"][i].GetNbinsX()+1):
        rebinnum=ROOT.TH1F("","",2*window_size+1,0,2*window_size+1)
        rebinden=ROOT.TH1F("","",2*window_size+1,0,2*window_size+1)

        rebinnum.Sumw2()
        rebinden.Sumw2()

        if j < window_size+1:

            for k in range(1,2*window_size+2):
                rebinnum.SetBinContent(k,wjets_2016["hists"][i].GetBinContent(k))
                rebinnum.SetBinError(k,wjets_2016["hists"][i].GetBinError(k))
                rebinden.SetBinContent(k,wjets_fake_photon_2016["hists"][i].GetBinContent(k))
                rebinden.SetBinError(k,wjets_fake_photon_2016["hists"][i].GetBinError(k))
        elif j > wjets_2016["hists"][i].GetNbinsX() - window_size-1:        
            for k in range(wjets_2016["hists"][i].GetNbinsX()+1-2*window_size-1,wjets_2016["hists"][i].GetNbinsX()+1):
                rebinnum.SetBinContent(k-wjets_2016["hists"][i].GetNbinsX()-1+2*window_size+1+1,wjets_2016["hists"][i].GetBinContent(k))
                rebinnum.SetBinError(k-wjets_2016["hists"][i].GetNbinsX()-1+2*window_size+1+1,wjets_2016["hists"][i].GetBinError(k))
                rebinden.SetBinContent(k-wjets_2016["hists"][i].GetNbinsX()-1+2*window_size+1+1,wjets_fake_photon_2016["hists"][i].GetBinContent(k))
                rebinden.SetBinError(k-wjets_2016["hists"][i].GetNbinsX()-1+2*window_size+1+1,wjets_fake_photon_2016["hists"][i].GetBinError(k))
        else:
            for k in range(j-window_size,j+window_size+1):
                rebinnum.SetBinContent(k-j+window_size+1,wjets_2016["hists"][i].GetBinContent(k))
                rebinnum.SetBinError(k-j+window_size+1,wjets_2016["hists"][i].GetBinError(k))
                rebinden.SetBinContent(k-j+window_size+1,wjets_fake_photon_2016["hists"][i].GetBinContent(k))
                rebinden.SetBinError(k-j+window_size+1,wjets_fake_photon_2016["hists"][i].GetBinError(k))

#        rebinden.Scale(0.75)
                 
        rebinnum.Rebin(2*window_size+1)
        rebinden.Rebin(2*window_size+1)

        num_smoothed.SetBinContent(j,rebinnum.GetBinContent(1))
        den_smoothed.SetBinContent(j,rebinden.GetBinContent(1))

        num_smoothed.SetBinError(j,rebinnum.GetBinError(1))
        den_smoothed.SetBinError(j,rebinden.GetBinError(1))

        rebinnum.Divide(rebinden)

        fake_photon_syst2_smoothed[len(fake_photon_syst2_smoothed)-1].SetBinContent(j,fake_photon["hists"][i].GetBinContent(j)*rebinnum.GetBinContent(1))
        fake_photon_syst2_smoothed[len(fake_photon_syst2_smoothed)-1].SetBinError(j,fake_photon["hists"][i].GetBinContent(j)*rebinnum.GetBinError(1))
        non_closure_smoothed[len(non_closure_smoothed)-1].SetBinContent(j,rebinnum.GetBinContent(1))
        non_closure_smoothed[len(non_closure_smoothed)-1].SetBinError(j,rebinnum.GetBinError(1))

    c = ROOT.TCanvas("c", "c",5,50,500,500)
    fake_photon["hists"][i].SetLineColor(ROOT.kRed)
    fake_photon_syst2_smoothed[len(fake_photon_syst2_smoothed)-1].SetLineColor(ROOT.kBlue)
    fake_photon["hists"][i].SetMinimum(0)
    fake_photon["hists"][i].SetMaximum(1.55*max(fake_photon["hists"][i].GetMaximum(),fake_photon_syst2_smoothed[len(fake_photon_syst2_smoothed)-1].GetMaximum()))
    fake_photon["hists"][i].Draw()
    fake_photon_syst2_smoothed[len(fake_photon_syst2_smoothed)-1].Draw("same")
    c.SaveAs(args.outputdir + "/" + "fake_photon_syst2_"+variables_labels[i]+".png")

    c = ROOT.TCanvas("c", "c",5,50,500,500)
    non_closure_smoothed[len(non_closure_smoothed)-1].SetMinimum(0.5)
    non_closure_smoothed[len(non_closure_smoothed)-1].SetMaximum(1.5)
    non_closure_smoothed[len(non_closure_smoothed)-1].Draw()
    c.SaveAs(args.outputdir + "/" + "non_closure_smoothed_"+variables_labels[i]+".png")

    c = ROOT.TCanvas("c", "c",5,50,500,500)
    num_smoothed.SetMinimum(0)
    num_smoothed.SetMaximum(1.55*max(num_smoothed.GetMaximum(),den_smoothed.GetMaximum()))
    num_smoothed.SetLineColor(ROOT.kRed)
    den_smoothed.SetLineColor(ROOT.kBlue)
    num_smoothed.Draw()
    den_smoothed.Draw("same")
    c.SaveAs(args.outputdir + "/" + "closure_smoothed_"+variables_labels[i]+".png")

c1 = ROOT.TCanvas("c1", "c1",5,50,500,500)

for i in range(len(variables)):

#    fake_lepton["hists"][i].Scale(2)

    fake_photon["hists"][i].Scale(1.0)

#    data["hists"][i].Print("all")
#    fake_photon["hists"][i].Print("all")
#    fake_lepton["hists"][i].Print("all")
#    if "wg+jets" in labels:
#        labels["wg+jets"]["hists"][i].Print("all")
#    if "w+jets" in labels:
#       labels["w+jets"]["hists"][i].Print("all")

    data["hists"][i].SetMarkerStyle(ROOT.kFullCircle)
    data["hists"][i].SetLineWidth(3)
    data["hists"][i].SetLineColor(ROOT.kBlack)

    ewdim6["hists"][i].SetLineWidth(3)
    ewdim6["hists"][i].SetLineColor(ROOT.kOrange+3)

    fake_photon["hists"][i].SetFillColor(ROOT.kGray+1)
    fake_lepton["hists"][i].SetFillColor(ROOT.kAzure-1)
    double_fake["hists"][i].SetFillColor(ROOT.kMagenta)
    e_to_p_total["hists"][i].SetFillColor(ROOT.kSpring)

    fake_photon["hists"][i].SetLineColor(ROOT.kGray+1)
    fake_lepton["hists"][i].SetLineColor(ROOT.kAzure-1)
    double_fake["hists"][i].SetLineColor(ROOT.kMagenta)
    e_to_p_total["hists"][i].SetLineColor(ROOT.kSpring)


    fake_photon["hists"][i].SetFillStyle(1001)
    fake_lepton["hists"][i].SetFillStyle(1001)
    double_fake["hists"][i].SetFillStyle(1001)
    e_to_p_total["hists"][i].SetFillStyle(1001)

    s=str(totallumi)+" fb^{-1} (13 TeV)"
    lumilabel = ROOT.TLatex (0.95, 0.93, s)
    lumilabel.SetNDC ()
    lumilabel.SetTextAlign (30)
    lumilabel.SetTextFont (42)
    lumilabel.SetTextSize (0.040)

    if "w+jets" in labels:
        wjets_fake_photon_2016["hists"][i].SetFillColor(labels["w+jets"]["color"])
        wjets_fake_photon_2016["hists"][i].SetFillStyle(1001)
        wjets_fake_photon_2016["hists"][i].SetLineColor(labels["w+jets"]["color"])

#
    hsum = data["hists"][i].Clone()
    hsum.Scale(0.0)

    hstack = ROOT.THStack()

    for label in labels.keys():
        if labels[label]["color"] == None:
            continue

        if label == "w+jets":
            hsum.Add(labels[label]["hists-prompt-pileup"][i])
            hstack.Add(labels[label]["hists-prompt-pileup"][i])
        elif label == "wg+jets":
            if args.draw_non_fid:
                hsum.Add(labels[label]["hists-fail-fiducial"][i])
                hstack.Add(labels[label]["hists-fail-fiducial"][i])
                hsum.Add(labels[label]["hists-pass-fiducial"][i])
                hstack.Add(labels[label]["hists-pass-fiducial"][i])
            else:    
                hsum.Add(labels[label]["hists"][i])
                hstack.Add(labels[label]["hists"][i])
        else:        
            hsum.Add(labels[label]["hists"][i])
            hstack.Add(labels[label]["hists"][i])

    if args.use_wjets_for_fake_photon and "w+jets" in labels:
        hsum.Add(labels["w+jets"]["hists"][i])
        hstack.Add(labels["w+jets"]["hists"][i])
#        hsum.Add(wjets_fake_photon_2016["hists"][i])
#        hstack.Add(wjets_fake_photon_2016["hists"][i])

    if data_driven:
        if not args.use_wjets_for_fake_photon:
            hstack.Add(fake_photon["hists"][i])
        hstack.Add(fake_lepton["hists"][i])
        hstack.Add(double_fake["hists"][i])

    hsum.Add(e_to_p_total["hists"][i])
    hstack.Add(e_to_p_total["hists"][i])

    if data_driven:
        if not args.use_wjets_for_fake_photon:
            hsum.Add(fake_photon["hists"][i])
        hsum.Add(fake_lepton["hists"][i])
        hsum.Add(double_fake["hists"][i])

    if data["hists"][i].GetMaximum() < hsum.GetMaximum():
        hsum.SetMaximum(hsum.GetMaximum()*1.55)
        data["hists"][i].SetMaximum(hsum.GetMaximum())
        hstack.SetMaximum(hsum.GetMaximum())
    else:
        data["hists"][i].SetMaximum(data["hists"][i].GetMaximum()*1.55)
        hstack.SetMaximum(data["hists"][i].GetMaximum())
        hsum.SetMaximum(data["hists"][i].GetMaximum())

    data["hists"][i].SetMinimum(0)
    hstack.SetMinimum(0)
    hsum.SetMinimum(0)

    tpad1=ROOT.TPad("", "", 0, 0.3, 1, 1.0)
    tpad1.SetBottomMargin(0.02)
    tpad1.Draw()
    
    c1.cd()

    tpad2 =ROOT.TPad("", "", 0, 0.0, 1, 0.275)
    tpad2.SetTopMargin(0)
    tpad2.SetBottomMargin(0.35)
    tpad2.Draw()

    tpad1.cd()

#    set_axis_fonts(data["hists"][i],"x",getXaxisLabel(variables[i]))
    hstack.Draw("hist") #you need to draw it before you can get the histogram and therefore the axis
    set_axis_fonts(hstack.GetHistogram(),"y","Events / bin")

    if variables[i] == "photon_pt_overflow":
        hstack.Draw("hist")
    else:
        hstack.Draw("hist")
        data["hists"][i].Draw("same")

    if args.draw_ewdim6:
        ewdim6["hists"][i].Print("all")
        ewdim6["hists"][i].Draw("same")
#wg_qcd.Draw("hist same")
#fake_lepton_hist.Draw("hist same")
#fake_photon_hist.Draw("hist same")

#wg_ewk_hist.Print("all")

    if hstack.GetHistogram().GetMaximum() > 100000:    
        cmslabel = ROOT.TLatex (0.27, 0.93, "#bf{CMS} #it{Preliminary}")
    else:
        cmslabel = ROOT.TLatex (0.18, 0.93, "#bf{CMS} #it{Preliminary}")
    #cmslabel = ROOT.TLatex (0.18, 0.93, "")
    cmslabel.SetNDC ()
    cmslabel.SetTextAlign (10)
    cmslabel.SetTextFont (42)
    cmslabel.SetTextSize (0.040)
    cmslabel.Draw ("same") 
    
    lumilabel.Draw("same")

#wpwpjjewk.Draw("same")

    j=0
    draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,data["hists"][i],"data","lp")

    if  args.use_wjets_for_fake_photon and "w+jets" in labels:
        j=j+1    
        draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,labels["w+jets"]["hists"][i],"W+jets","f")

    if data_driven :
        if not args.use_wjets_for_fake_photon:
            j=j+1
            draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,fake_photon["hists"][i],"nonprompt photon","f")
        j=j+1
        if lepton_name == "muon":
            draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,fake_lepton["hists"][i],"nonprompt muon","f")
        elif lepton_name == "electron":
            draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,fake_lepton["hists"][i],"nonprompt electron","f")
        elif lepton_name == "both":
            draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,fake_lepton["hists"][i],"nonprompt lepton","f")
        else:
            assert(0)
        j=j+1
        draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,double_fake["hists"][i],"double nonprompt","f")

    j=j+1
    draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,e_to_p_total["hists"][i],"e-induced #gamma","f")

    for label in labels.keys():
        if labels[label]["color"] == None:
            continue

        if label == "w+jets":
            j=j+1    
            draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,labels[label]["hists-prompt-pileup"][i],"pileup","f")
        elif label == "wg+jets":
            if args.draw_non_fid:
                j=j+1    
                draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,labels[label]["hists-pass-fiducial"][i],"W#gamma+jets","f")
                j=j+1    
                draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,labels[label]["hists-fail-fiducial"][i],"W#gamma+jets out","f")
            else:    
                j=j+1    
#                draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,labels[label]["hists"][i],label,"f")
                draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,labels[label]["hists"][i],labels[label]["legend"],"f")
        else:        
            j=j+1    
            draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,labels[label]["hists"][i],labels[label]["legend"],"f")

    if args.draw_ewdim6:
        j=j+1
        draw_legend(xpositions[j],0.84 - ypositions[j]*yoffset,ewdim6["hists"][i],"C_{B} = 51","l")

    gstat = ROOT.TGraphAsymmErrors(hsum);

    for j in range(0,gstat.GetN()):
        gstat.SetPointEYlow (j, hsum.GetBinError(j+1));
        gstat.SetPointEYhigh(j, hsum.GetBinError(j+1));

    gstat.SetFillColor(12);
    gstat.SetFillStyle(3345);
    gstat.SetMarkerSize(0);
    gstat.SetLineWidth(0);
    gstat.SetLineColor(ROOT.kWhite);

    gsyst = ROOT.TGraphAsymmErrors(hsum);
    gratiosyst = ROOT.TGraphAsymmErrors(hsum);

    wgjets_pdf_syst=histogram_models[i].GetHistogram()

    for j in range(labels["wg+jets"]["hists-pdf-variation0"][i].GetNbinsX()+1):
        mean_pdf=0

        for k in range(1,32):

            mean_pdf += labels["wg+jets"]["hists-pdf-variation"+str(k)][i].GetBinContent(j)*labels["wg+jets"]["hists"][i].Integral()/labels["wg+jets"]["hists-pdf-variation"+str(k)][i].Integral()

        mean_pdf = mean_pdf/31

        stddev_pdf = 0

        for k in range(1,32):
            stddev_pdf += pow(labels["wg+jets"]["hists-pdf-variation"+str(k)][i].GetBinContent(j)*labels["wg+jets"]["hists"][i].Integral()/labels["wg+jets"]["hists-pdf-variation"+str(k)][i].Integral() - mean_pdf,2)

        stddev_pdf = sqrt(stddev_pdf/(31-1))

        wgjets_pdf_syst.SetBinContent(j,labels["wg+jets"]["hists"][i].GetBinContent(j)+stddev_pdf)

    wgjets_pass_fiducial_pdf_syst=histogram_models[i].GetHistogram()

    for j in range(labels["wg+jets"]["hists-pass-fiducial-pdf-variation0"][i].GetNbinsX()+1):
        mean_pdf=0

        for k in range(1,32):
            mean_pdf += labels["wg+jets"]["hists-pass-fiducial-pdf-variation"+str(k)][i].GetBinContent(j)*labels["wg+jets"]["hists-pass-fiducial"][i].Integral()/labels["wg+jets"]["hists-pass-fiducial-pdf-variation"+str(k)][i].Integral()

        mean_pdf = mean_pdf/31

        stddev_pdf = 0

        for k in range(1,32):
            stddev_pdf += pow(labels["wg+jets"]["hists-pass-fiducial-pdf-variation"+str(k)][i].GetBinContent(j)*labels["wg+jets"]["hists-pass-fiducial"][i].Integral()/labels["wg+jets"]["hists-pass-fiducial-pdf-variation"+str(k)][i].Integral() - mean_pdf,2)

        stddev_pdf = sqrt(stddev_pdf/(31-1))

        wgjets_pass_fiducial_pdf_syst.SetBinContent(j,labels["wg+jets"]["hists-pass-fiducial"][i].GetBinContent(j)+stddev_pdf)

    wgjets_fail_fiducial_pdf_syst=histogram_models[i].GetHistogram()

    for j in range(labels["wg+jets"]["hists-pdf-variation0"][i].GetNbinsX()+1):
        mean_pdf=0

        for k in range(1,32):
            mean_pdf += labels["wg+jets"]["hists-fail-fiducial-pdf-variation"+str(k)][i].GetBinContent(j)

        mean_pdf = mean_pdf/31

        stddev_pdf = 0

        for k in range(1,32):
            stddev_pdf += pow(labels["wg+jets"]["hists-fail-fiducial-pdf-variation"+str(k)][i].GetBinContent(j) - mean_pdf,2)

        stddev_pdf = sqrt(stddev_pdf/(31-1))

        wgjets_fail_fiducial_pdf_syst.SetBinContent(j,labels["wg+jets"]["hists-fail-fiducial"][i].GetBinContent(j)+stddev_pdf)

    wgjets_scale_syst=histogram_models[i].GetHistogram()

    for j in range(labels["wg+jets"]["hists-scale-variation0"][i].GetNbinsX()+1):
        wgjets_scale_syst.SetBinContent(j,labels["wg+jets"]["hists"][i].GetBinContent(j)+labels["wg+jets"]["hists"][i].Integral()*max(
            abs(labels["wg+jets"]["hists-scale-variation0"][i].GetBinContent(j)/labels["wg+jets"]["hists-scale-variation0"][i].Integral() - labels["wg+jets"]["hists"][i].GetBinContent(j)/labels["wg+jets"]["hists"][i].Integral()),
            abs(labels["wg+jets"]["hists-scale-variation1"][i].GetBinContent(j)/labels["wg+jets"]["hists-scale-variation1"][i].Integral() - labels["wg+jets"]["hists"][i].GetBinContent(j)/labels["wg+jets"]["hists"][i].Integral()),
            abs(labels["wg+jets"]["hists-scale-variation3"][i].GetBinContent(j)/labels["wg+jets"]["hists-scale-variation3"][i].Integral() - labels["wg+jets"]["hists"][i].GetBinContent(j)/labels["wg+jets"]["hists"][i].Integral()),
            abs(labels["wg+jets"]["hists-scale-variation4"][i].GetBinContent(j)/labels["wg+jets"]["hists-scale-variation4"][i].Integral() - labels["wg+jets"]["hists"][i].GetBinContent(j)/labels["wg+jets"]["hists"][i].Integral()),
            abs(labels["wg+jets"]["hists-scale-variation5"][i].GetBinContent(j)/labels["wg+jets"]["hists-scale-variation5"][i].Integral() - labels["wg+jets"]["hists"][i].GetBinContent(j)/labels["wg+jets"]["hists"][i].Integral()),
            abs(labels["wg+jets"]["hists-scale-variation6"][i].GetBinContent(j)/labels["wg+jets"]["hists-scale-variation6"][i].Integral() - labels["wg+jets"]["hists"][i].GetBinContent(j)/labels["wg+jets"]["hists"][i].Integral())))

    wgjets_pass_fiducial_scale_syst=histogram_models[i].GetHistogram()

    for j in range(labels["wg+jets"]["hists-pass-fiducial-scale-variation0"][i].GetNbinsX()+1):
        wgjets_pass_fiducial_scale_syst.SetBinContent(j,labels["wg+jets"]["hists-pass-fiducial"][i].GetBinContent(j)+labels["wg+jets"]["hists-pass-fiducial"][i].Integral()*max(
            abs(labels["wg+jets"]["hists-pass-fiducial-scale-variation0"][i].GetBinContent(j)/labels["wg+jets"]["hists-pass-fiducial-scale-variation0"][i].Integral() - labels["wg+jets"]["hists-pass-fiducial"][i].GetBinContent(j)/labels["wg+jets"]["hists-pass-fiducial"][i].Integral()),
            abs(labels["wg+jets"]["hists-pass-fiducial-scale-variation1"][i].GetBinContent(j)/labels["wg+jets"]["hists-pass-fiducial-scale-variation1"][i].Integral() - labels["wg+jets"]["hists-pass-fiducial"][i].GetBinContent(j)/labels["wg+jets"]["hists-pass-fiducial"][i].Integral()),
            abs(labels["wg+jets"]["hists-pass-fiducial-scale-variation3"][i].GetBinContent(j)/labels["wg+jets"]["hists-pass-fiducial-scale-variation3"][i].Integral() - labels["wg+jets"]["hists-pass-fiducial"][i].GetBinContent(j)/labels["wg+jets"]["hists-pass-fiducial"][i].Integral()),
            abs(labels["wg+jets"]["hists-pass-fiducial-scale-variation4"][i].GetBinContent(j)/labels["wg+jets"]["hists-pass-fiducial-scale-variation4"][i].Integral() - labels["wg+jets"]["hists-pass-fiducial"][i].GetBinContent(j)/labels["wg+jets"]["hists-pass-fiducial"][i].Integral()),
            abs(labels["wg+jets"]["hists-pass-fiducial-scale-variation5"][i].GetBinContent(j)/labels["wg+jets"]["hists-pass-fiducial-scale-variation5"][i].Integral() - labels["wg+jets"]["hists-pass-fiducial"][i].GetBinContent(j)/labels["wg+jets"]["hists-pass-fiducial"][i].Integral()),
            abs(labels["wg+jets"]["hists-pass-fiducial-scale-variation6"][i].GetBinContent(j)/labels["wg+jets"]["hists-pass-fiducial-scale-variation6"][i].Integral() - labels["wg+jets"]["hists-pass-fiducial"][i].GetBinContent(j)/labels["wg+jets"]["hists-pass-fiducial"][i].Integral())))

    wgjets_fail_fiducial_scale_syst=histogram_models[i].GetHistogram()

    for j in range(labels["wg+jets"]["hists-fail-fiducial-scale-variation0"][i].GetNbinsX()+1):
        wgjets_fail_fiducial_scale_syst.SetBinContent(j,labels["wg+jets"]["hists-fail-fiducial"][i].GetBinContent(j)+max(
            abs(labels["wg+jets"]["hists-fail-fiducial-scale-variation0"][i].GetBinContent(j) - labels["wg+jets"]["hists-fail-fiducial"][i].GetBinContent(j)),
            abs(labels["wg+jets"]["hists-fail-fiducial-scale-variation1"][i].GetBinContent(j) - labels["wg+jets"]["hists-fail-fiducial"][i].GetBinContent(j)),
            abs(labels["wg+jets"]["hists-fail-fiducial-scale-variation3"][i].GetBinContent(j) - labels["wg+jets"]["hists-fail-fiducial"][i].GetBinContent(j)),
            abs(labels["wg+jets"]["hists-fail-fiducial-scale-variation4"][i].GetBinContent(j) - labels["wg+jets"]["hists-fail-fiducial"][i].GetBinContent(j)),
            abs(labels["wg+jets"]["hists-fail-fiducial-scale-variation5"][i].GetBinContent(j) - labels["wg+jets"]["hists-fail-fiducial"][i].GetBinContent(j)),
            abs(labels["wg+jets"]["hists-fail-fiducial-scale-variation6"][i].GetBinContent(j) - labels["wg+jets"]["hists-fail-fiducial"][i].GetBinContent(j))))

    zgjets_pdf_syst=histogram_models[i].GetHistogram()

    for j in range(labels["zg+jets"]["hists-pdf-variation0"][i].GetNbinsX()+1):
        mean_pdf=0

        for k in range(1,32):
            mean_pdf += labels["zg+jets"]["hists-pdf-variation"+str(k)][i].GetBinContent(j)

        mean_pdf = mean_pdf/31

        stddev_pdf = 0

        for k in range(1,32):
            stddev_pdf += pow(labels["zg+jets"]["hists-pdf-variation"+str(k)][i].GetBinContent(j) - mean_pdf,2)

        stddev_pdf = sqrt(stddev_pdf/(31-1))

        zgjets_pdf_syst.SetBinContent(j,labels["zg+jets"]["hists"][i].GetBinContent(j)+stddev_pdf)

    zgjets_scale_syst=histogram_models[i].GetHistogram()

    for j in range(labels["zg+jets"]["hists-scale-variation0"][i].GetNbinsX()+2):
        zgjets_scale_syst.SetBinContent(j,labels["zg+jets"]["hists"][i].GetBinContent(j)+max(
            abs(labels["zg+jets"]["hists-scale-variation0"][i].GetBinContent(j) - labels["zg+jets"]["hists"][i].GetBinContent(j)),
            abs(labels["zg+jets"]["hists-scale-variation1"][i].GetBinContent(j) - labels["zg+jets"]["hists"][i].GetBinContent(j)),
            abs(labels["zg+jets"]["hists-scale-variation3"][i].GetBinContent(j) - labels["zg+jets"]["hists"][i].GetBinContent(j)),
            abs(labels["zg+jets"]["hists-scale-variation4"][i].GetBinContent(j) - labels["zg+jets"]["hists"][i].GetBinContent(j)),
            abs(labels["zg+jets"]["hists-scale-variation5"][i].GetBinContent(j) - labels["zg+jets"]["hists"][i].GetBinContent(j)),
            abs(labels["zg+jets"]["hists-scale-variation6"][i].GetBinContent(j) - labels["zg+jets"]["hists"][i].GetBinContent(j))))

    for j in range(1,gsyst.GetN()+1):
        total_unc = 0
        
        total_unc += pow(zgjets_scale_syst.GetBinContent(j) - labels["zg+jets"]["hists"][i].GetBinContent(j),2)

        total_unc += pow(zgjets_pdf_syst.GetBinContent(j) - labels["zg+jets"]["hists"][i].GetBinContent(j),2)

        for k in range(n_fake_photon_alt):
            total_unc += pow(fake_photon["hists-alt"+str(k)][i].GetBinContent(j)-fake_photon["hists"][i].GetBinContent(j),2)

        total_unc += pow(fake_photon_syst2_smoothed[i].GetBinContent(j)-fake_photon["hists"][i].GetBinContent(j),2)

        total_unc += pow(0.3*fake_lepton["hists"][i].GetBinContent(j),2)

        for label in ["gg+jets","top+jets","vv+jets","zg+jets"]:
            total_unc += pow(0.018*labels[label]["hists"][i].GetBinContent(j),2)

        if args.draw_non_fid:
            total_unc += pow(0.018*labels["wg+jets"]["hists-pass-fiducial-"+unc+"-up"][i].GetBinContent(j),2)
            total_unc += pow(0.018*labels["wg+jets"]["hists-fail-fiducial-"+unc+"-up"][i].GetBinContent(j),2)
            total_unc += pow(wgjets_pass_fiducial_pdf_syst.GetBinContent(j) - labels["wg+jets"]["hists-pass-fiducial"][i].GetBinContent(j),2)
            total_unc += pow(wgjets_pass_fiducial_scale_syst.GetBinContent(j) - labels["wg+jets"]["hists-pass-fiducial"][i].GetBinContent(j),2)
            total_unc += pow(wgjets_fail_fiducial_pdf_syst.GetBinContent(j) - labels["wg+jets"]["hists-fail-fiducial"][i].GetBinContent(j),2)
            total_unc += pow(wgjets_fail_fiducial_scale_syst.GetBinContent(j) - labels["wg+jets"]["hists-fail-fiducial"][i].GetBinContent(j),2)
        else:
            pass
            total_unc += pow(0.018*labels["wg+jets"]["hists"][i].GetBinContent(j),2)
            total_unc += pow(wgjets_pdf_syst.GetBinContent(j) - labels["wg+jets"]["hists"][i].GetBinContent(j),2)
            total_unc += pow(wgjets_scale_syst.GetBinContent(j) - labels["wg+jets"]["hists"][i].GetBinContent(j),2)

        for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
            if unc != "jer" and unc != "jes":
                total_unc += pow(labels["w+jets"]["hists-prompt-pileup-"+unc+"-up"][i].GetBinContent(j) - labels["w+jets"]["hists-prompt-pileup"][i].GetBinContent(j),2)
                total_unc += pow(labels["gg+jets"]["hists"][i].GetBinContent(j) - labels["gg+jets"]["hists"][i].GetBinContent(j),2)
            for label in ["top+jets","vv+jets","zg+jets"]:
                total_unc += pow(labels[label]["hists-"+unc+"-up"][i].GetBinContent(j) - labels[label]["hists"][i].GetBinContent(j),2)
            if args.draw_non_fid:
                total_unc += pow(labels["wg+jets"]["hists-pass-fiducial-"+unc+"-up"][i].GetBinContent(j) - labels["wg+jets"]["hists-pass-fiducial"][i].GetBinContent(j),2)
                total_unc += pow(labels["wg+jets"]["hists-fail-fiducial-"+unc+"-up"][i].GetBinContent(j) - labels["wg+jets"]["hists-fail-fiducial"][i].GetBinContent(j),2)
            else:
                total_unc += pow(labels["wg+jets"]["hists-"+unc+"-up"][i].GetBinContent(j) - labels["wg+jets"]["hists"][i].GetBinContent(j),2)

        total_unc += pow(hsum.GetBinError(j),2)

        total_unc = sqrt(total_unc)

        gsyst.SetPointEYlow(j-1, total_unc)
        gsyst.SetPointEYlow(j-1, total_unc)

        gratiosyst.SetPoint(j-1, hsum.GetXaxis().GetBinCenter(j), 1)
        if hsum.GetBinContent(j) > 0:
            gratiosyst.SetPointEYlow(j-1, total_unc/hsum.GetBinContent(j))
            gratiosyst.SetPointEYhigh(j-1, total_unc/hsum.GetBinContent(j))
        else:
            gratiosyst.SetPointEYlow(j-1, 0)
            gratiosyst.SetPointEYhigh(j-1, 0)

    gsyst.SetFillColor(12);
    gsyst.SetFillStyle(3345);
    gsyst.SetMarkerSize(0);
    gsyst.SetLineWidth(0);
    gsyst.SetLineColor(ROOT.kWhite);

    gratiosyst.SetFillColor(12);
    gratiosyst.SetFillStyle(3345);
    gratiosyst.SetMarkerSize(0);
    gratiosyst.SetLineWidth(0);
    gratiosyst.SetLineColor(ROOT.kWhite);

#    gstat.Draw("E2same");
    gsyst.Draw("E2same");

    if variables[i] != "photon_pt_overflow":
        data["hists"][i].Draw("same")

    tpad2.cd()

    ratio=histogram_models[i].GetHistogram()

    numerator=data["hists"][i].Clone()
    denominator=hsum.Clone()

    for j in range(1,denominator.GetNbinsX()):
        denominator.SetBinError(j,0)
                                 
    ratio.Add(numerator)    

    ratio.Divide(denominator)


    ratio.GetXaxis().SetLabelFont  (   42)
    ratio.GetXaxis().SetLabelOffset(0.015)
    ratio.GetXaxis().SetLabelSize  (0.150)

    ratio.GetYaxis().SetLabelFont  (   42)
    ratio.GetYaxis().SetLabelOffset(0.015)
    ratio.GetYaxis().SetLabelSize  (0.100)

    ratio.GetYaxis().SetTitleFont  (   42)
    ratio.GetYaxis().SetTitleOffset(0.45)
    ratio.GetYaxis().SetTitleSize  (0.100)
    ratio.GetYaxis().SetTitle("data / pred")

    ratio.GetXaxis().SetTitleFont (42)
    ratio.GetXaxis().SetTitleOffset (0.95)
    ratio.GetXaxis().SetTitleSize (0.15)
    ratio.GetXaxis().SetNdivisions (505)
    ratio.GetXaxis().SetTitle(getXaxisLabel(variables[i]))

    if args.ewdim6 and i == ewdim6_index:
        ratio.Scale(0)

    if i == ewdim6_index:
        ratio.SetMaximum(1.5)
        ratio.SetMinimum(0.5)
    else:
        ratio.SetMaximum(1.3)
        ratio.SetMinimum(0.7)

    ratio.Draw()
    gratiosyst.Draw("E2same")

    c1.cd()

    c1.Update()
    c1.ForceUpdate()
    c1.Modified()

    c1.SaveAs(args.outputdir + "/" + variables_labels[i] + ".png")

c1.Close()

wg_jets_integral_error = ROOT.Double()
wg_jets_integral = labels["wg+jets"]["hists"][mlg_index].IntegralAndError(1,labels["wg+jets"]["hists"][mlg_index].GetXaxis().GetNbins(),wg_jets_integral_error)

wg_jets_fid_integral_error = ROOT.Double()
wg_jets_fid_integral = labels["wg+jets"]["hists-pass-fiducial"][mlg_index].IntegralAndError(1,labels["wg+jets"]["hists-pass-fiducial"][mlg_index].GetXaxis().GetNbins(),wg_jets_fid_integral_error)

wg_jets_nonfid_integral_error = ROOT.Double()
wg_jets_nonfid_integral = labels["wg+jets"]["hists-fail-fiducial"][mlg_index].IntegralAndError(1,labels["wg+jets"]["hists-fail-fiducial"][mlg_index].GetXaxis().GetNbins(),wg_jets_nonfid_integral_error)

gg_jets_integral_error = ROOT.Double()
gg_jets_integral = labels["gg+jets"]["hists"][mlg_index].IntegralAndError(1,labels["gg+jets"]["hists"][mlg_index].GetXaxis().GetNbins(),gg_jets_integral_error)

zg_jets_integral_error = ROOT.Double()
zg_jets_integral = labels["zg+jets"]["hists"][mlg_index].IntegralAndError(1,labels["zg+jets"]["hists"][mlg_index].GetXaxis().GetNbins(),zg_jets_integral_error)

vv_jets_integral_error = ROOT.Double()
vv_jets_integral = labels["vv+jets"]["hists"][mlg_index].IntegralAndError(1,labels["vv+jets"]["hists"][mlg_index].GetXaxis().GetNbins(),vv_jets_integral_error)

top_jets_integral_error = ROOT.Double()
top_jets_integral = labels["top+jets"]["hists"][mlg_index].IntegralAndError(1,labels["top+jets"]["hists"][mlg_index].GetXaxis().GetNbins(),top_jets_integral_error)

fake_lepton_integral_error = ROOT.Double()
fake_lepton_integral = fake_lepton["hists"][mlg_index].IntegralAndError(1,fake_lepton["hists"][mlg_index].GetXaxis().GetNbins(),fake_lepton_integral_error)

fake_photon_integral_error = ROOT.Double()
fake_photon_integral = fake_photon["hists"][mlg_index].IntegralAndError(1,fake_photon["hists"][mlg_index].GetXaxis().GetNbins(),fake_photon_integral_error)

double_fake_integral_error = ROOT.Double()
double_fake_integral = double_fake["hists"][mlg_index].IntegralAndError(1,double_fake["hists"][mlg_index].GetXaxis().GetNbins(),double_fake_integral_error)

data_integral_error = ROOT.Double()
data_integral = data["hists"][mlg_index].IntegralAndError(1,data["hists"][mlg_index].GetXaxis().GetNbins(),data_integral_error)

e_to_p_total_integral_error = ROOT.Double()
e_to_p_total_integral = e_to_p_total["hists"][mlg_index].IntegralAndError(1,e_to_p_total["hists"][mlg_index].GetXaxis().GetNbins(),e_to_p_total_integral_error)

pileup_integral_error = ROOT.Double()
pileup_integral = labels["w+jets"]["hists-prompt-pileup"][mlg_index].IntegralAndError(1,labels["w+jets"]["hists-prompt-pileup"][mlg_index].GetXaxis().GetNbins(),pileup_integral_error)

fake_signal_contamination_integral_error = ROOT.Double()
fake_signal_contamination_integral = fake_signal_contamination["hists"][mlg_index].IntegralAndError(1,fake_signal_contamination["hists"][mlg_index].GetXaxis().GetNbins(),fake_signal_contamination_integral_error)

#print "fake signal contamination = "+str(fake_signal_contamination_integral) + " +/- " +str(fake_signal_contamination_integral_error)

#print "wg+jets = "+str(wg_jets_integral)+" +/- "+str(wg_jets_integral_error)
#print "wg+jets fid = "+str(wg_jets_fid_integral)+" +/- "+str(wg_jets_fid_integral_error)
#print "wg+jets non-fid = "+str(wg_jets_nonfid_integral)+" +/- "+str(wg_jets_nonfid_integral_error)
#print "zg+jets = "+str(zg_jets_integral)+" +/- "+str(zg_jets_integral_error)
#print "gg+jets = "+str(gg_jets_integral)+" +/- "+str(gg_jets_integral_error)
#print "vv+jets = "+str(vv_jets_integral)+" +/- "+str(vv_jets_integral_error)
#print "top+jets = "+str(top_jets_integral)+" +/- "+str(top_jets_integral_error)
#print "fake photon = "+str(fake_photon_integral)+" +/- "+str(fake_photon_integral_error)
#print "fake lepton = "+str(fake_lepton_integral)+" +/- "+str(fake_lepton_integral_error)
#print "double fake = "+str(double_fake_integral)+" +/- "+str(double_fake_integral_error)
#print "data = "+str(data_integral)+" +/- "+str(data_integral_error)
#print "e_to_p_total = "+str(e_to_p_total_integral)+" +/- "+str(e_to_p_total_integral_error)
#print "pileup = "+str(pileup_integral)+" +/- "+str(pileup_integral_error)

print """
\begin{table}[htbp]
\begin{center}
\begin{tabular}{|c|c|}
\hline
process & expected number of events  \\
\hline \hline
$W\gamma$ & $%0.f \pm %0.f$ \\
\hline
$W\gamma$ out & $%0.f \pm %0.f$ \\
\hline
$Z\gamma$ & $%0.f \pm %0.f$ \\
\hline
top & $%0.f \pm %0.f$ \\
\hline
$VV$ & $%0.f \pm %0.f$ \\
\hline
fake photon & $%0.f \pm %0.f$ \\
\hline
fake lepton & $%0.f \pm %0.f$ \\
\hline
double fake & $%0.f \pm %0.f$ \\
\hline
electron-induced photon & $%0.f \pm %0.f$ \\
\hline
pileup & $%0.f \pm %0.f$ \\
\hline
\end{tabular}
\end{center}
\caption{Number of background expected events per category in the %s channel. The uncertainty is statistical. $W\gamma$ and $W\gamma$ out are the contributions to the signal region from the $W\gamma$ process originating from inside and outside the fiducial region, respectively.}
\label{tab:%s_n_background_events}
\end{table}
"""%(
wg_jets_fid_integral,float(wg_jets_fid_integral_error),
wg_jets_nonfid_integral,float(wg_jets_nonfid_integral_error),
zg_jets_integral,float(zg_jets_integral_error),
top_jets_integral,float(top_jets_integral_error),
vv_jets_integral,float(vv_jets_integral_error),
fake_photon_integral,float(fake_photon_integral_error),
fake_lepton_integral,float(fake_lepton_integral_error),
double_fake_integral,float(double_fake_integral_error),
e_to_p_total_integral,float(e_to_p_total_integral_error),
pileup_integral,float(pileup_integral_error),
args.lep,args.lep
)

n_signal = data_integral - double_fake_integral - fake_photon_integral - fake_lepton_integral - top_jets_integral - vv_jets_integral - zg_jets_integral - e_to_p_total_integral

n_signal_error = sqrt(pow(data_integral_error,2) + pow(double_fake_integral_error,2) + pow(fake_lepton_integral_error,2)+ pow(fake_photon_integral_error,2)+pow(top_jets_integral_error,2)+ pow(vv_jets_integral_error,2)+ pow(zg_jets_integral_error,2)+pow(e_to_p_total_integral_error,2))

#print "n_signal = "+str(n_signal) + " +/- " + str(n_signal_error)

#double_fake["hists"][mlg_index].Print("all")
#fake_lepton["hists"][mlg_index].Print("all")
#fake_photon["hists"][mlg_index].Print("all")
#fake_photon_stat_up["hists"][mlg_index].Print("all")

#labels["w+jets"]["hists-prompt-pileup"][mlg_index].Print("all")

zgjets_pdf_syst=histogram_models[mlg_index].GetHistogram()

for i in range(labels["zg+jets"]["hists-pdf-variation0"][mlg_index].GetNbinsX()+1):
    mean_pdf=0

    for j in range(1,32):
        mean_pdf += labels["zg+jets"]["hists-pdf-variation"+str(j)][mlg_index].GetBinContent(i)

    mean_pdf = mean_pdf/31

    stddev_pdf = 0

    for j in range(1,32):
        stddev_pdf += pow(labels["zg+jets"]["hists-pdf-variation"+str(j)][mlg_index].GetBinContent(i) - mean_pdf,2)

    stddev_pdf = sqrt(stddev_pdf/(31-1))

    zgjets_pdf_syst.SetBinContent(i,labels["zg+jets"]["hists"][mlg_index].GetBinContent(i)+stddev_pdf)

zgjets_scale_syst=histogram_models[mlg_index].GetHistogram()

for i in range(labels["zg+jets"]["hists-scale-variation0"][mlg_index].GetNbinsX()+2):
    zgjets_scale_syst.SetBinContent(i,labels["zg+jets"]["hists"][mlg_index].GetBinContent(i)+max(
        abs(labels["zg+jets"]["hists-scale-variation0"][mlg_index].GetBinContent(i) - labels["zg+jets"]["hists"][mlg_index].GetBinContent(i)),
        abs(labels["zg+jets"]["hists-scale-variation1"][mlg_index].GetBinContent(i) - labels["zg+jets"]["hists"][mlg_index].GetBinContent(i)),
        abs(labels["zg+jets"]["hists-scale-variation3"][mlg_index].GetBinContent(i) - labels["zg+jets"]["hists"][mlg_index].GetBinContent(i)),
        abs(labels["zg+jets"]["hists-scale-variation4"][mlg_index].GetBinContent(i) - labels["zg+jets"]["hists"][mlg_index].GetBinContent(i)),
        abs(labels["zg+jets"]["hists-scale-variation5"][mlg_index].GetBinContent(i) - labels["zg+jets"]["hists"][mlg_index].GetBinContent(i)),
        abs(labels["zg+jets"]["hists-scale-variation6"][mlg_index].GetBinContent(i) - labels["zg+jets"]["hists"][mlg_index].GetBinContent(i))))

fake_electron_syst_up=histogram_models[mlg_index].GetHistogram()

fake_electron_syst_up.Add(fake_lepton["hists"][mlg_index])

fake_electron_syst_up.Scale(1.3)

trandom=ROOT.TRandom3()
trandom.SetSeed(0)

fake_photon_syst2_up=[]
fake_photon_syst2_up_relative=[]

#nnuisances = 100
nnuisances = histogram_models[mlg_index].GetHistogram().GetNbinsX()

#for i in range(1,fake_photon["hists"][mlg_index].GetNbinsX()+1):
#for i in range(1,histogram_models[mlg_index].GetHistogram().GetNbinsX()+1):
for i in range(1,nnuisances+1):
#for i in range(1,1000+1):
    fake_photon_syst2_up.append(histogram_models[mlg_index].GetHistogram())
    fake_photon_syst2_up_relative.append(histogram_models[mlg_index].GetHistogram())    

    sequence = []
    old = 0
    for j in range(histogram_models[mlg_index].GetHistogram().GetNbinsX()):
        if j == 0:
            old = trandom.Uniform(-1,1)
            sequence.append(old)
        else:
            found = False
            while (not found):
                new = trandom.Uniform(-1,1)

#                if abs(old - new) > 0.01:
                if False:
                    continue
                else:
                    old = new
                    sequence.append(old)
                    found = True


    for j in range(1,fake_photon["hists"][mlg_index].GetNbinsX()+1):
        if i == j:
#            if fake_photon_syst2_up[len(fake_photon_syst2_up)-1].GetXaxis().GetBinCenter(j) < 60:
#                fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j)*1.2)
#            elif fake_photon_syst2_up[len(fake_photon_syst2_up)-1].GetXaxis().GetBinCenter(j) > 60 and fake_photon_syst2_up[len(fake_photon_syst2_up)-1].GetXaxis().GetBinCenter(j) < 120:
#                fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j)*1.1)
#            elif fake_photon_syst2_up[len(fake_photon_syst2_up)-1].GetXaxis().GetBinCenter(j) > 120:
#                fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j)*1.4)
#            else:
#                assert(0)

#        if True:
#        if j >= 0 and j <= 120 :
#            fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j)+fake_photon["hists"][mlg_index].GetBinError(j))
#            fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j)*1.2)
#            fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j)*(1+sequence[j-1]))

#            non_closure.GetBinContent(non_closure[mlg_index].GetXaxis().FindFixBin(fake_photon["hists"][mlg_index].GetXaxis().GetBinCenter(j))

#            fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j)*(1+non_closure[28].GetBinContent(non_closure[28].GetXaxis().FindFixBin(fake_photon["hists"][mlg_index].GetXaxis().GetBinCenter(j)))))

            bin =  non_closure[mlg_index].GetXaxis().FindFixBin(fake_photon["hists"][mlg_index].GetXaxis().GetBinCenter(j))

#            fake_photon_syst2_up_relative[len(fake_photon_syst2_up_relative)-1].SetBinContent(j,sequence[j-1]*non_closure[mlg_index].GetBinContent(non_closure[mlg_index].GetXaxis().FindFixBin(fake_photon["hists"][mlg_index].GetXaxis().GetBinCenter(j))))

            fake_photon_syst2_up_relative[len(fake_photon_syst2_up_relative)-1].SetBinContent(j,non_closure[mlg_index].GetBinContent(bin)+sequence[j-1]*non_closure[mlg_index].GetBinError(bin))            

#            fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j)*(1+non_closure[mlg_index].GetBinContent(bin)+sequence[j-1]*non_closure[mlg_index].GetBinError(bin)/sqrt(nnuisances)))

#            fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j)*(1+non_closure[mlg_index].GetBinContent(bin)+non_closure[mlg_index].GetBinError(bin)))

            fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j)*(1+abs(non_closure[mlg_index].GetBinContent(bin))+non_closure[mlg_index].GetBinError(bin)))

#            fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j)*(1+non_closure[mlg_index].GetBinContent(bin)))

#            fake_photon["hists"][mlg_index].SetBinError(j,fake_photon["hists"][mlg_index].GetBinContent(j)*(abs(non_closure[mlg_index].GetBinContent(bin))+non_closure[mlg_index].GetBinError(bin)))

#            fake_photon["hists"][mlg_index].SetBinError(j,0)

#            fake_photon["hists"][mlg_index].SetBinError(j,fake_photon["hists"][mlg_index].GetBinContent(j)*2)

        else:
            bin =  non_closure[mlg_index].GetXaxis().FindFixBin(fake_photon["hists"][mlg_index].GetXaxis().GetBinCenter(j))
#            fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j)*(1+non_closure[mlg_index].GetBinContent(bin)))
            fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinContent(j,fake_photon["hists"][mlg_index].GetBinContent(j))

        fake_photon_syst2_up[len(fake_photon_syst2_up)-1].SetBinError(j,0)

#    fake_photon_syst2_up[len(fake_photon_syst2_up)-1].Scale(fake_photon["hists"][mlg_index].Integral()/fake_photon_syst2_up[len(fake_photon_syst2_up)-1].Integral())
#    fake_photon_syst2_up[len(fake_photon_syst2_up)-1].Scale(1/sqrt(histogram_models[mlg_index].GetHistogram().GetNbinsX()))

#for i in range(1,histogram_models[mlg_index].GetHistogram().GetNbinsX()+1):
#
#    if i == 1:
#        fake_photon["hists"][mlg_index].SetBinContent(i,fake_photon["hists"][mlg_index].GetBinContent(i)*(1+non_closure[mlg_index].GetBinContent(i)+non_closure[mlg_index].GetBinError(i)))
#    else:
#        fake_photon["hists"][mlg_index].SetBinContent(i,fake_photon["hists"][mlg_index].GetBinContent(i)*(1+non_closure[mlg_index].GetBinContent(i)))

#    fake_photon["hists"][mlg_index].SetBinError(i,fake_photon["hists"][mlg_index].GetBinContent(i)*(non_closure[mlg_index].GetBinError(i)))



fake_photon["hists"][mlg_index].SetLineWidth(2)

if len(fake_photon_syst2_up) > 4:

    fake_photon_syst2_up[0].SetLineColor(ROOT.kRed)
    fake_photon_syst2_up[1].SetLineColor(ROOT.kBlue)
    fake_photon_syst2_up[2].SetLineColor(ROOT.kGreen)
    fake_photon_syst2_up[3].SetLineColor(ROOT.kMagenta)
    fake_photon_syst2_up[4].SetLineColor(ROOT.kPink)
    
    fake_photon_syst2_up[0].SetLineWidth(2)
    fake_photon_syst2_up[1].SetLineWidth(2)
    fake_photon_syst2_up[2].SetLineWidth(2)
    fake_photon_syst2_up[3].SetLineWidth(2)
    fake_photon_syst2_up[4].SetLineWidth(2)

non_closure[mlg_index].SetLineWidth(2)

if len(fake_photon_syst2_up_relative) > 4:

    fake_photon_syst2_up_relative[0].SetLineColor(ROOT.kRed)
    fake_photon_syst2_up_relative[1].SetLineColor(ROOT.kBlue)
    fake_photon_syst2_up_relative[2].SetLineColor(ROOT.kGreen)
    fake_photon_syst2_up_relative[3].SetLineColor(ROOT.kMagenta)
    fake_photon_syst2_up_relative[4].SetLineColor(ROOT.kPink)

    fake_photon_syst2_up_relative[0].SetLineWidth(2)
    fake_photon_syst2_up_relative[1].SetLineWidth(2)
    fake_photon_syst2_up_relative[2].SetLineWidth(2)
    fake_photon_syst2_up_relative[3].SetLineWidth(2)
    fake_photon_syst2_up_relative[4].SetLineWidth(2)

c = ROOT.TCanvas("","")

if len(fake_photon_syst2_up) > 4:

    fake_photon["hists"][mlg_index].Draw()
    fake_photon_syst2_up[0].Draw("same")
    fake_photon_syst2_up[1].Draw("same")
    fake_photon_syst2_up[2].Draw("same")
    fake_photon_syst2_up[3].Draw("same")
    fake_photon_syst2_up[4].Draw("same")

c.SaveAs(args.outputdir + "/"+ "fake_photon_syst2.png")

if len(fake_photon_syst2_up_relative) > 4:

    non_closure[mlg_index].Draw()
    fake_photon_syst2_up_relative[0].Draw("same")
    fake_photon_syst2_up_relative[1].Draw("same")
    fake_photon_syst2_up_relative[2].Draw("same")
    #fake_photon_syst2_up_relative[3].Draw("same")
    #fake_photon_syst2_up_relative[4].Draw("same")

c.SaveAs(args.outputdir + "/"+ "fake_photon_syst2_relative.png")
c.Close()

fake_lepton_syst_up=[]

for i in range(1,fake_lepton["hists"][mlg_index].GetNbinsX()+1):
    fake_lepton_syst_up.append(histogram_models[mlg_index].GetHistogram())
    for j in range(1,fake_lepton["hists"][mlg_index].GetNbinsX()+1):
        if i == j:
#            fake_lepton_syst_up[len(fake_lepton_syst_up)-1].SetBinContent(j,fake_lepton["hists"][mlg_index].GetBinContent(j)+fake_lepton["hists"][mlg_index].GetBinError(j))
            fake_lepton_syst_up[len(fake_lepton_syst_up)-1].SetBinContent(j,1.3*fake_lepton["hists"][mlg_index].GetBinContent(j))
        else:
            fake_lepton_syst_up[len(fake_lepton_syst_up)-1].SetBinContent(j,fake_lepton["hists"][mlg_index].GetBinContent(j))
        fake_lepton_syst_up[len(fake_lepton_syst_up)-1].SetBinError(j,0)

zgjets_stat_up=[]

for i in range(1,labels["zg+jets"]["hists"][mlg_index].GetNbinsX()+1):
    zgjets_stat_up.append(histogram_models[mlg_index].GetHistogram())
    for j in range(1,labels["zg+jets"]["hists"][mlg_index].GetNbinsX()+1):
        if i == j:
            zgjets_stat_up[len(zgjets_stat_up)-1].SetBinContent(j,labels["zg+jets"]["hists"][mlg_index].GetBinContent(j)+labels["zg+jets"]["hists"][mlg_index].GetBinError(j))
        else:
            zgjets_stat_up[len(zgjets_stat_up)-1].SetBinContent(j,labels["zg+jets"]["hists"][mlg_index].GetBinContent(j))
        zgjets_stat_up[len(zgjets_stat_up)-1].SetBinError(j,0)

wg_stat_up=[]

for i in range(1,labels["wg+jets"]["hists"][mlg_index].GetNbinsX()+1):
    wg_stat_up.append(histogram_models[mlg_index].GetHistogram())
    for j in range(1,labels["wg+jets"]["hists"][mlg_index].GetNbinsX()+1):
        if i == j:
            wg_stat_up[len(wg_stat_up)-1].SetBinContent(j,labels["wg+jets"]["hists"][mlg_index].GetBinContent(j)+labels["wg+jets"]["hists"][mlg_index].GetBinError(j))
        else:    
            wg_stat_up[len(wg_stat_up)-1].SetBinContent(j,labels["wg+jets"]["hists"][mlg_index].GetBinContent(j))
        wg_stat_up[len(wg_stat_up)-1].SetBinError(j,0)

if args.make_cut_and_count_datacard:

    if args.lep == "muon":
        dcard = open("wg_dcard_cut_and_count_mu_chan.txt",'w')
    elif args.lep == "electron":
        dcard = open("wg_dcard_cut_and_count_el_chan.txt",'w')
    else:
        assert(0)

    print >> dcard, "imax 1 number of channels"
    print >> dcard, "jmax * number of background"
    print >> dcard, "kmax * number of nuisance parameters"

    print >> dcard, "Observation "+str(data["hists"][mlg_index].Integral())
    dcard.write("bin")

    if args.lep == "muon":
        dcard.write(" mu_chan")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" mu_chan")

        dcard.write(" mu_chan")
        dcard.write(" mu_chan")
        dcard.write(" mu_chan")
        dcard.write(" mu_chan")
        dcard.write('\n')    
    elif args.lep == "electron":
        dcard.write(" el_chan")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" el_chan")

        dcard.write(" el_chan")
        dcard.write(" el_chan")
        dcard.write(" el_chan")
        dcard.write(" el_chan")
        dcard.write('\n')    
    else:
        assert(0)

    dcard.write("process")
    dcard.write(" Wg")
        
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" " + label.replace("+",""))

    if "w+jets" in labels:
        dcard.write(" pileup")

    dcard.write(" fake_photon")
    dcard.write(" fake_lepton")
    dcard.write(" double_fake")
    dcard.write(" e_to_p_non_res")
    dcard.write('\n')    
    dcard.write("process")
    dcard.write(" 0")

    for j in range(1,len(labels.keys())+2):
        dcard.write(" " + str(j))
    dcard.write('\n')    
    dcard.write('rate')
    dcard.write(' '+str(labels["wg+jets"]["hists"][mlg_index].Integral()))

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" "+ str(labels[label]["hists"][mlg_index].Integral()))

    if "w+jets" in labels:
        dcard.write(" "+str(labels["w+jets"]["hists-prompt-pileup"][mlg_index].Integral())) 
    dcard.write(" "+str(fake_photon["hists"][mlg_index].Integral())) 
    dcard.write(" "+str(fake_lepton["hists"][mlg_index].Integral())) 
    dcard.write(" "+str(double_fake["hists"][mlg_index].Integral())) 
    dcard.write(" "+str(e_to_p_total["hists"][mlg_index].Integral())) 
   
    dcard.write('\n')    

    dcard.write("lumi lnN")
    dcard.write(" 1.018")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" 1.018")

    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" 1.018")
    
    dcard.write('\n')    

    dcard.write("pileup lnN")
    dcard.write(" "+str(labels["wg+jets"]["hists-pileup-up"][mlg_index].Integral()/labels["wg+jets"]["hists"][mlg_index].Integral()))

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" "+str(labels[label]["hists-pileup-up"][mlg_index].Integral()/labels[label]["hists"][mlg_index].Integral()))

    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")

    dcard.write('\n')    

    dcard.write("prefire lnN")
    dcard.write(" "+str(labels["wg+jets"]["hists-prefire-up"][mlg_index].Integral()/labels["wg+jets"]["hists"][mlg_index].Integral()))

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" "+str(labels[label]["hists-prefire-up"][mlg_index].Integral()/labels[label]["hists"][mlg_index].Integral()))

    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")

    dcard.write('\n')    

    dcard.write("jes lnN")
    dcard.write(" "+str(labels["wg+jets"]["hists-jes-up"][mlg_index].Integral()/labels["wg+jets"]["hists"][mlg_index].Integral()))

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" "+str(labels[label]["hists-jes-up"][mlg_index].Integral()/labels[label]["hists"][mlg_index].Integral()))

    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")

    dcard.write('\n')    

    dcard.write("jer lnN")
    dcard.write(" "+str(labels["wg+jets"]["hists-jer-up"][mlg_index].Integral()/labels["wg+jets"]["hists"][mlg_index].Integral()))

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" "+str(labels[label]["hists-jer-up"][mlg_index].Integral()/labels[label]["hists"][mlg_index].Integral()))

    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")

    dcard.write('\n')    

    dcard.write("muonid lnN")
    dcard.write(" "+str(labels["wg+jets"]["hists-muon-id-sf-up"][mlg_index].Integral()/labels["wg+jets"]["hists"][mlg_index].Integral()))

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" "+str(labels[label]["hists-muon-id-sf-up"][mlg_index].Integral()/labels[label]["hists"][mlg_index].Integral()))

    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")

    dcard.write('\n')    

    dcard.write("muoniso lnN")
    dcard.write(" "+str(labels["wg+jets"]["hists-muon-iso-sf-up"][mlg_index].Integral()/labels["wg+jets"]["hists"][mlg_index].Integral()))

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" "+str(labels[label]["hists-muon-iso-sf-up"][mlg_index].Integral()/labels[label]["hists"][mlg_index].Integral()))

    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")

    dcard.write('\n')    

    dcard.write("muonhlt lnN")
    dcard.write(" "+str(labels["wg+jets"]["hists-muon-hlt-sf-up"][mlg_index].Integral()/labels["wg+jets"]["hists"][mlg_index].Integral()))

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" "+str(labels[label]["hists-muon-hlt-sf-up"][mlg_index].Integral()/labels[label]["hists"][mlg_index].Integral()))

    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")

    dcard.write('\n')    

    dcard.write("fakephotonsyst1 lnN")
    dcard.write(" -")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" -")

    dcard.write(" "+str(fake_photon_alt["hists"][mlg_index].Integral()/fake_photon["hists"][mlg_index].Integral()))
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")

    dcard.write('\n')    

    dcard.write("fakephotonsyst2 lnN")
    dcard.write(" -")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" -")

    dcard.write(" "+str(labels["w+jets"]["hists"][mlg_index].Integral()/fake_photon["hists"][mlg_index].Integral()))
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")

    dcard.write('\n')    

    dcard.write("wgscale lnN")
    dcard.write(" "+str(1))

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")

    dcard.write('\n')    

    dcard.write("zgscale lnN")
    dcard.write(" -")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        
        if label == "zg+jets":
            dcard.write(" "+str(zgjets_scale_syst.Integral()/(labels["zg+jets"]["hists"][mlg_index].Integral())))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")

    dcard.write('\n')    

    dcard.write("zgpdf lnN")
    dcard.write(" -")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        
        if label == "zg+jets":
            dcard.write(" "+str(zgjets_pdf_syst.Integral()/(labels["zg+jets"]["hists"][mlg_index].Integral())))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")

    dcard.write('\n')    

    dcard.write("fakemuonsyst lnN")
    dcard.write(" -")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" -")

    dcard.write(" -")
    dcard.write(" 1.3")
    dcard.write(" 1.3")
    dcard.write(" -")
    
    dcard.write('\n')    

if args.make_datacard:

    print "Making datacard"

    data["hists"][mlg_index].Scale(0)
    data["hists"][mlg_index].Add(labels["gg+jets"]["hists"][mlg_index])
    data["hists"][mlg_index].Add(labels["wg+jets"]["hists"][mlg_index])
    data["hists"][mlg_index].Add(labels["top+jets"]["hists"][mlg_index])
    data["hists"][mlg_index].Add(labels["zg+jets"]["hists"][mlg_index])
    data["hists"][mlg_index].Add(labels["vv+jets"]["hists"][mlg_index])
    data["hists"][mlg_index].Add(e_to_p_total["hists"][mlg_index])
#    if args.lep == "electron":
#        data["hists"][mlg_index].Add(e_to_p["hists"][mlg_index])
    data["hists"][mlg_index].Add(fake_photon["hists"][mlg_index])
    data["hists"][mlg_index].Add(fake_lepton["hists"][mlg_index])
    data["hists"][mlg_index].Add(double_fake["hists"][mlg_index])

    tmphist = fake_photon["hists"][mlg_index].Clone("tmphist")

    for i in range(1,histogram_models[mlg_index].GetHistogram().GetNbinsX()+1):
        if i == 1:
            if non_closure[mlg_index].GetBinContent(i) > 0:
                sign = 1
            else:
                sign = -1

#            fake_photon["hists"][mlg_index].SetBinContent(i,fake_photon["hists"][mlg_index].GetBinContent(i)*(1+non_closure[mlg_index].GetBinContent(i)+sign*non_closure[mlg_index].GetBinError(i)))
#        else:
#            fake_photon["hists"][mlg_index].SetBinContent(i,fake_photon["hists"][mlg_index].GetBinContent(i)*(1+non_closure[mlg_index].GetBinContent(i)))


    if args.lep == "muon":
        dcard = open("wg_dcard_mu_chan.txt",'w')
    elif args.lep == "electron":
        dcard = open("wg_dcard_el_chan.txt",'w')
    else:
        assert(0)
        
    print >> dcard, "imax 1 number of channels"
    print >> dcard, "jmax * number of background"
    print >> dcard, "kmax * number of nuisance parameters"

    if args.lep == "muon":
        print >> dcard, "shapes data_obs mu_chan wg_dcard_mu_chan_shapes.root data_obs"
        print >> dcard, "shapes Wg mu_chan wg_dcard_mu_chan_shapes.root wg wg_$SYSTEMATIC" 
    elif args.lep == "electron":
        print >> dcard, "shapes data_obs el_chan wg_dcard_el_chan_shapes.root data_obs"
        print >> dcard, "shapes Wg el_chan wg_dcard_el_chan_shapes.root wg wg_$SYSTEMATIC" 
    else:
        assert(0)    

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if args.lep == "muon":
            print >> dcard, "shapes "+label.replace("+","")+" mu_chan wg_dcard_mu_chan_shapes.root "+label.replace("+","")+ " " +label.replace("+","") + "_$SYSTEMATIC" 
        elif args.lep == "electron":
            print >> dcard, "shapes "+label.replace("+","")+" el_chan wg_dcard_el_chan_shapes.root "+label.replace("+","")+ " " +label.replace("+","") + "_$SYSTEMATIC" 
        else:
            assert(0)    

    if args.lep == "muon":
        print >> dcard, "shapes Wg_out mu_chan wg_dcard_mu_chan_shapes.root wgout wgout_$SYSTEMATIC" 
        print >> dcard, "shapes pileup mu_chan wg_dcard_mu_chan_shapes.root pileup pileup_$SYSTEMATIC" 
        print >> dcard, "shapes fake_photon mu_chan wg_dcard_mu_chan_shapes.root fakephoton fakephoton_$SYSTEMATIC" 
        print >> dcard, "shapes fake_muon mu_chan wg_dcard_mu_chan_shapes.root fakemuon fakemuon_$SYSTEMATIC"
        print >> dcard, "shapes double_fake mu_chan wg_dcard_mu_chan_shapes.root doublefake doublefake_$SYSTEMATIC" 
#        print >> dcard, "shapes e_to_p_non_res mu_chan wg_dcard_mu_chan_shapes.root etopnonres etopnonres_$SYSTEMATIC" 
        for i in range(len(etopbinning)):
            print >> dcard, "shapes e_to_p_bin"+str(i)+" mu_chan wg_dcard_mu_chan_shapes.root etopbin"+str(i)+" etopbin"+str(i)+"_$SYSTEMATIC" 
    elif args.lep == "electron":
        print >> dcard, "shapes Wg_out el_chan wg_dcard_el_chan_shapes.root wgout wgout_$SYSTEMATIC" 
        print >> dcard, "shapes pileup el_chan wg_dcard_el_chan_shapes.root pileup pileup_$SYSTEMATIC" 
        print >> dcard, "shapes fake_photon el_chan wg_dcard_el_chan_shapes.root fakephoton fakephoton_$SYSTEMATIC" 
        print >> dcard, "shapes fake_electron el_chan wg_dcard_el_chan_shapes.root fakeelectron fakeelectron_$SYSTEMATIC"
        print >> dcard, "shapes double_fake el_chan wg_dcard_el_chan_shapes.root doublefake doublefake_$SYSTEMATIC" 
        for i in range(len(etopbinning)):
            print >> dcard, "shapes e_to_p_bin"+str(i)+" el_chan wg_dcard_el_chan_shapes.root etopbin"+str(i)+" etopbin"+str(i)+"_$SYSTEMATIC" 
    else:
        assert(0)    

#    print >> dcard, "shapes data_obs el_chan wg_dcard_el_chan_shapes.root workspace:data_obs"
#    print >> dcard, "shapes Wg el_chan wg_dcard_el_chan_shapes.root workspace:wg workspace:wg_$SYSTEMATIC" 

#    for label in labels.keys():
#        if label == "no label" or label == "wg+jets" or label == "w+jets":
#            continue
#        print >> dcard, "shapes "+label.replace("+","")+" el_chan wg_dcard_el_chan_shapes.root workspace:"+label.replace("+","")+ " workspace:" +label.replace("+","") + "_$SYSTEMATIC" 

#    print >> dcard, "shapes fake_photon el_chan wg_dcard_el_chan_shapes.root workspace:fakephoton workspace:fakephoton_$SYSTEMATIC" 
#    print >> dcard, "shapes fake_electron el_chan wg_dcard_el_chan_shapes.root workspace:fakeelectron workspace:fakeelectron_$SYSTEMATIC"
#    print >> dcard, "shapes double_fake el_chan wg_dcard_el_chan_shapes.root workspace:doublefake workspace:doublefake_$SYSTEMATIC" 
#    print >> dcard, "shapes e_to_p_non_res el_chan wg_dcard_el_chan_shapes.root workspace:etopnonres workspace:etopnonres_$SYSTEMATIC" 
#    print >> dcard, "shapes e_to_p el_chan wg_dcard_el_chan_shapes.root workspace:etop workspace:etop_$SYSTEMATIC" 
    
    print >> dcard, "Observation "+str(data["hists"][mlg_index].Integral())
    dcard.write("bin")
    if args.lep == "muon":
        dcard.write(" mu_chan")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" mu_chan")
        if "w+jets" in labels:
            dcard.write(" mu_chan")    
        dcard.write(" mu_chan")
        dcard.write(" mu_chan")
        dcard.write(" mu_chan")
        dcard.write(" mu_chan")
        for i in range(len(etopbinning)):
            dcard.write(" mu_chan")
        dcard.write('\n')    
    elif args.lep == "electron":
        dcard.write(" el_chan")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" el_chan")
        if "w+jets" in labels:
            dcard.write(" el_chan")        
        dcard.write(" el_chan")
        dcard.write(" el_chan")
        dcard.write(" el_chan")
        dcard.write(" el_chan")
        for i in range(len(etopbinning)):
            dcard.write(" el_chan")
        dcard.write('\n')    
    else:
        assert(0)    

    dcard.write("process")
    sys.stdout.write("process")
    dcard.write(" Wg")
    sys.stdout.write(" Wg")
    dcard.write(" Wg_out")
    sys.stdout.write(" Wg_out")
        
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" " + label.replace("+",""))
        sys.stdout.write(" " + label.replace("+",""))

    if "w+jets" in labels:
        dcard.write(" pileup")
        sys.stdout.write(" pileup")
    dcard.write(" fake_photon")
    sys.stdout.write(" fake_photon")
    if args.lep == "muon":
        dcard.write(" fake_muon")
        sys.stdout.write(" fake_muon")
    elif args.lep == "electron":
        dcard.write(" fake_electron")
        sys.stdout.write(" fake_electron")
    else:
        assert(0)    
    dcard.write(" double_fake")
    sys.stdout.write(" double_fake")
    for i in range(len(etopbinning)):
        dcard.write(" e_to_p_bin"+str(i))
        sys.stdout.write(" e_to_p_bin"+str(i))
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
#        dcard.write(" e_to_p")
        pass
    else:
        assert(0)    
    dcard.write('\n')
    sys.stdout.write('\n')
    dcard.write("process")
#    dcard.write(" 0")

    nprocesses = 0
    nprocesses += 2 #wg+jets and wg+jets out-of-fiducial 
    for label in labels:
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        nprocesses += 1    
    if "w+jets" in labels:
        nprocesses += 1 #pileup
    nprocesses += 3 #fake photon, fake lepton, double fake
    nprocesses += len(etopbinning) 
    for i in range(nprocesses-1):
        dcard.write(" " + str(i))
    dcard.write(" -1")    


    if args.lep == "muon":
#        for j in range(1,len(labels.keys())+3):
#            dcard.write(" " + str(j))
        pass
    elif args.lep == "electron":
#        for j in range(1,len(labels.keys())+4):
#            dcard.write(" " + str(j))
        pass
    else:
        assert(0)    
    dcard.write('\n')
    dcard.write('rate')
    sys.stdout.write('rate')
#    dcard.write(' '+str(labels["wg+jets"]["hists"][mlg_index].Integral()))
    dcard.write(' '+str(labels["wg+jets"]["hists-pass-fiducial"][mlg_index].Integral()))
    sys.stdout.write(' '+str(labels["wg+jets"]["hists-pass-fiducial"][mlg_index].Integral()))
    dcard.write(' '+str(labels["wg+jets"]["hists-fail-fiducial"][mlg_index].Integral()))
    sys.stdout.write(' '+str(labels["wg+jets"]["hists-fail-fiducial"][mlg_index].Integral()))

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" "+ str(labels[label]["hists"][mlg_index].Integral()))
        sys.stdout.write(" "+ str(labels[label]["hists"][mlg_index].Integral()))

    if "w+jets" in labels:
        dcard.write(" "+str(labels["w+jets"]["hists-prompt-pileup"][mlg_index].Integral())) 
        sys.stdout.write(" "+str(labels["w+jets"]["hists-prompt-pileup"][mlg_index].Integral()))
    dcard.write(" "+str(fake_photon["hists"][mlg_index].Integral())) 
    sys.stdout.write(" "+str(fake_photon["hists"][mlg_index].Integral()))
    dcard.write(" "+str(fake_lepton["hists"][mlg_index].Integral())) 
    sys.stdout.write(" "+str(fake_lepton["hists"][mlg_index].Integral()))
    dcard.write(" "+str(double_fake["hists"][mlg_index].Integral())) 
    sys.stdout.write(" "+str(double_fake["hists"][mlg_index].Integral()))
    for i in range(len(etopbinning)):
        dcard.write(" "+str(e_to_p[i]["hists"][mlg_index].Integral())) 
        sys.stdout.write(" "+str(e_to_p[i]["hists"][mlg_index].Integral()))
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
        pass
#        dcard.write(" "+str(e_to_p["hists"][mlg_index].Integral())) 
    else:
        assert(0)    
    dcard.write('\n')
    sys.stdout.write('\n')

    dcard.write("lumi lnN")
    dcard.write(" 1.018")
    dcard.write(" 1.018")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" 1.018")

    if "w+jets" in labels:
        dcard.write(" 1.018")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    for i in range(len(etopbinning)):
        dcard.write(" 1.018")
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
#        dcard.write(" -")
        pass
    else:
        assert(0)    
    
    dcard.write('\n')    

    dcard.write("pileup shape1")
    dcard.write(" 1.0")
    dcard.write(" 1.0")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" 1.0")

    if "w+jets" in labels:
        dcard.write(" 1.0")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    for i in range(len(etopbinning)):
        dcard.write(" -")
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
#        dcard.write(" -")
        pass
    else:
        assert(0)    
    dcard.write('\n')    

    dcard.write("prefire shape1")
    dcard.write(" 1.0")
    dcard.write(" 1.0")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" 1.0")

    if "w+jets" in labels:
        dcard.write(" 1.0")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    for i in range(len(etopbinning)):
        dcard.write(" -")
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
#        dcard.write(" -")
        pass
    else:
        assert(0)    
    dcard.write('\n')    
    
    dcard.write("jes shape1")
    dcard.write(" 1.0")
    dcard.write(" 1.0")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if label != "gg+jets": 
            dcard.write(" 1.0")
        else:
            dcard.write(" -")

    if "w+jets" in labels:
        dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    for i in range(len(etopbinning)):
        dcard.write(" -")
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
#        dcard.write(" -")
        pass
    else:
        assert(0)    
    dcard.write('\n')    
    
    dcard.write("jer shape1")
    dcard.write(" 1.0")
    dcard.write(" 1.0")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if label != "gg+jets":
            dcard.write(" 1.0")
        else:
            dcard.write(" -")

    if "w+jets" in labels:
        dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    for i in range(len(etopbinning)):
        dcard.write(" -")
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
#        dcard.write(" -")
        pass
    else:
        assert(0)    
    dcard.write('\n')    

    dcard.write("photonid shape1")
    dcard.write(" 1.0")
    dcard.write(" 1.0")
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" 1.0")
    if "w+jets" in labels:
        dcard.write(" 1.0")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    for i in range(len(etopbinning)):
        dcard.write(" -")
    dcard.write('\n')    
    
    if args.lep == "muon":
        dcard.write("muonid shape1")
        dcard.write(" 1.0")
        dcard.write(" 1.0")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" 1.0")
        if "w+jets" in labels:
            dcard.write(" 1.0")
        dcard.write(" -")
        dcard.write(" -")
        dcard.write(" -")
        for i in range(len(etopbinning)):
            dcard.write(" -")
        dcard.write('\n')    
        dcard.write("muonhlt shape1")
        dcard.write(" 1.0")
        dcard.write(" 1.0")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" 1.0")
        if "w+jets" in labels:
            dcard.write(" 1.0")
        dcard.write(" -")
        dcard.write(" -")
        dcard.write(" -")
        for i in range(len(etopbinning)):
            dcard.write(" -")
        dcard.write('\n')    
        dcard.write("muoniso shape1")
        dcard.write(" 1.0")
        dcard.write(" 1.0")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" 1.0")
        if "w+jets" in labels:
            dcard.write(" 1.0")
        dcard.write(" -")
        dcard.write(" -")
        dcard.write(" -")
        for i in range(len(etopbinning)):
            dcard.write(" -")
        dcard.write('\n')    
    elif args.lep == "electron":
        dcard.write("electronreco shape1")
        dcard.write(" 1.0")
        dcard.write(" 1.0")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" 1.0")
        if "w+jets" in labels:
            dcard.write(" 1.0")
        dcard.write(" -")
        dcard.write(" -")
        dcard.write(" -")
        for i in range(len(etopbinning)):
            dcard.write(" -")
        dcard.write('\n')    
        dcard.write("electronid shape1")
        dcard.write(" 1.0")
        dcard.write(" 1.0")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" 1.0")
        if "w+jets" in labels:
            dcard.write(" 1.0")
        dcard.write(" -")
        dcard.write(" -")
        dcard.write(" -")
        for i in range(len(etopbinning)):
            dcard.write(" -")
        dcard.write('\n')    
        dcard.write("electronhlt shape1")
        dcard.write(" 1.0")
        dcard.write(" 1.0")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" 1.0")
        if "w+jets" in labels:
            dcard.write(" 1.0")
        dcard.write(" -")
        dcard.write(" -")
        dcard.write(" -")
        for i in range(len(etopbinning)):
            dcard.write(" -")
        dcard.write('\n')    
    else:
        assert(0)    

    for i in range(n_fake_photon_alt):
        dcard.write("fakephotonsyst1var"+str(i)+" shape1")
        dcard.write(" -")
        dcard.write(" -")

        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" -")

        if "w+jets" in labels:
            dcard.write(" -")
        dcard.write(" 1.0")
        dcard.write(" -")
        dcard.write(" -")
        for i in range(len(etopbinning)):
            dcard.write(" -")
        if args.lep == "muon":
            pass
        elif args.lep == "electron":
            pass
#        dcard.write(" -")
        else:
            assert(0)    
        dcard.write('\n')    
    
#    dcard.write("fakephotonsyst2 lnN")
#    dcard.write(" -")

#    for label in labels.keys():
#        if label == "no label" or label == "wg+jets" or label == "w+jets":
#            continue
#        dcard.write(" -")

#    dcard.write(" 1.4")
#    dcard.write(" -")
#    dcard.write(" 1.4")
#    dcard.write(" -")
#    if args.lep == "muon":
#        pass
#    elif args.lep == "electron":
#        dcard.write(" -")
#    else:
#        assert(0)    
#    dcard.write('\n')    

#    dcard.write("fakeelectronsyst lnN")
#    dcard.write(" -")

#    for label in labels.keys():
#        if label == "no label" or label == "wg+jets" or label == "w+jets":
#            continue
#        dcard.write(" -")

#    dcard.write(" -")
#    dcard.write(" 1.3")
#    dcard.write(" 1.3")
#    dcard.write(" -")
#    if args.lep == "muon":
#        pass
#    elif args.lep == "electron":
#        dcard.write(" -")
#    else:
#        assert(0)    
#    dcard.write('\n')    
    
    dcard.write("wgscale shape1")
    dcard.write(" 1.0")
    dcard.write(" -")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" -")

    if "w+jets" in labels:
        dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    for i in range(len(etopbinning)):
        dcard.write(" -")
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
        pass
#        dcard.write(" -")
    else:
        assert(0)    
    dcard.write('\n')    

    dcard.write("wgpdf shape1")
    dcard.write(" 1.0")
    dcard.write(" -")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" -")
    if "w+jets" in labels:
        dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    for i in range(len(etopbinning)):
        dcard.write(" -")
#    dcard.write(" -")
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
        pass
#        dcard.write(" -")
    else:
        assert(0)    
    dcard.write('\n')    

    dcard.write("wgoutscale shape1")
    dcard.write(" -")
    dcard.write(" 1.0")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" -")
    if "w+jets" in labels:
        dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    for i in range(len(etopbinning)):
        dcard.write(" -")
#    dcard.write(" -")
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
#        dcard.write(" -")
        pass
    else:
        assert(0)    
    dcard.write('\n')    

    dcard.write("wgoutpdf shape1")
    dcard.write(" -")
    dcard.write(" 1.0")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" -")
    if "w+jets" in labels:
        dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    for i in range(len(etopbinning)):
        dcard.write(" -")
#    dcard.write(" -")
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
#        dcard.write(" -")
        pass
    else:
        assert(0)    
    dcard.write('\n')    

    dcard.write("zgscale shape1")
    dcard.write(" -")
    dcard.write(" -")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        
        if label == "zg+jets":
            dcard.write(" 1.0")
        else:    
            dcard.write(" -")
    if "w+jets" in labels:
        dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    for i in range(len(etopbinning)):
        dcard.write(" -")
#    dcard.write(" -")
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
#        dcard.write(" -")
        pass
    else:
        assert(0)    
    dcard.write('\n')    

    dcard.write("zgpdf shape1")
    dcard.write(" -")
    dcard.write(" -")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        
        if label == "zg+jets":
            dcard.write(" 1.0")
        else:    
            dcard.write(" -")
    if "w+jets" in labels:
        dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    for i in range(len(etopbinning)):
        dcard.write(" -")
#    dcard.write(" -")
    if args.lep == "muon":
        pass
    elif args.lep == "electron":
#        dcard.write(" -")
        pass
    else:
        assert(0)    
    dcard.write('\n')    

#    dcard.write("sigfakesubtraction shape1")
#    dcard.write(" -")
#    dcard.write(" -")

#    for label in labels.keys():
#        if label == "no label" or label == "wg+jets" or label == "w+jets":
#            continue
#        
#        if label == "zg+jets":
#            dcard.write(" -")
#        else:    
#            dcard.write(" -")
#    if "w+jets" in labels:
#        dcard.write(" -")
#    dcard.write(" -")
#    dcard.write(" 1.0")
#    dcard.write(" -")
#
#    for i in range(len(etopbinning)):
#        dcard.write(" -")
#    if args.lep == "muon":
#        pass
#    elif args.lep == "electron":
#        pass
#    else:
#        assert(0)    
#    dcard.write('\n')    

    dcard.write("* autoMCStats 0\n")

#    if args.lep == "electron":
#        for i in range(len(etopbinning)):
#            dcard.write("etopbin"+str(i)+"norm rateParam el_chan e_to_p_bin"+str(i)+" 2 [0,10]\n")        
#    elif args.lep == "muon":
#        for i in range(len(etopbinning)):
#            dcard.write("etopbin"+str(i)+"norm rateParam mu_chan e_to_p_bin"+str(i)+" 2 [0,10]\n")        
#    else:
#        assert(0)

    if args.lep == "muon":
        pass
    elif args.lep == "electron":
        pass
#        dcard.write("etopnorm rateParam el_chan e_to_p 2 [0,10]\n")
    else:
        assert(0)    


#    for i in range(1,labels["wg+jets"]["hists"][mlg_index].GetNbinsX()+1):
#        dcard.write("wgstatbin"+str(i)+" shape1")
#        dcard.write(" 1.0")

#        for label in labels.keys():
#            if label == "no label" or label == "wg+jets" or label == "w+jets":
#                continue
#            dcard.write(" -")
        
#        dcard.write(" -")
#        dcard.write(" -")
#        dcard.write(" -")
#        dcard.write(" -")
#        dcard.write(" -")

#        dcard.write('\n')    

#    for i in range(1,fake_photon["hists"][mlg_index].GetNbinsX()+1):
    for i in range(1,len(fake_photon_syst2_up)+1):
        dcard.write("fakephotonsyst2var"+str(i)+" shape1")
        dcard.write(" -")
        dcard.write(" -")

        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" -")
        if "w+jets" in labels:
            dcard.write(" -")        
        dcard.write(" 1.0")
        dcard.write(" -")
        dcard.write(" -")
        for i in range(len(etopbinning)):
            dcard.write(" -")
#        dcard.write(" -")
        if args.lep == "muon":
            pass
        elif args.lep == "electron":
            pass
#            dcard.write(" -")
        else:
            assert(0)    

        dcard.write('\n')    

    for i in range(1,fake_photon["hists"][mlg_index].GetNbinsX()+1):
        if args.lep == "muon":
            dcard.write("fakemuonsystbin"+str(i)+" shape1")
        elif args.lep == "electron":
            dcard.write("fakeelectronsystbin"+str(i)+" shape1")
        else:
            assert(0)

        dcard.write(" -")
        dcard.write(" -")

        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" -")
        if "w+jets" in labels:
            dcard.write(" -")        
        dcard.write(" -")
        dcard.write(" 1.0")
        dcard.write(" -")
        for i in range(len(etopbinning)):
            dcard.write(" -")
        if args.lep == "muon":
            pass
        elif args.lep == "electron":
            pass
#            dcard.write(" -")
        else:
            assert(0)    

        dcard.write('\n')

    if args.lep == "muon":
        dcard.write("fakephotonsyst1 group = ")
        for i in range(n_fake_photon_alt):
            dcard.write("fakephotonsyst1var"+str(i)+" ")
        dcard.write(' \n')
        dcard.write("fakephotonsyst2 group = ")
        for i in range(1,len(fake_photon_syst2_up)+1):
            dcard.write("fakephotonsyst2var"+str(i)+" ")
        dcard.write(' \n')
        dcard.write('fakemuonsyst group = ')
        for i in range(1,fake_photon["hists"][mlg_index].GetNbinsX()+1):
            dcard.write("fakemuonsystbin"+str(i)+" ")
        dcard.write(' \n')
    elif args.lep == "electron":
        dcard.write('fakeelectronsyst group = ')
        for i in range(1,fake_photon["hists"][mlg_index].GetNbinsX()+1):
            dcard.write("fakeelectronsystbin"+str(i)+" ")
        dcard.write(' \n')
    else:
        assert(0)

#    for i in range(1,labels["zg+jets"]["hists"][mlg_index].GetNbinsX()+1):
#        dcard.write("zgjetsstatbin"+str(i)+" shape1")
#        dcard.write(" -")

#        for label in labels.keys():
#            if label == "no label" or label == "wg+jets" or label == "w+jets":
#                continue
            
#            if label == "zg+jets":
#                dcard.write(" 1.0")
#            else:    
#                dcard.write(" -")
        
#        dcard.write(" -")
#        dcard.write(" -")
#        dcard.write(" -")
#        dcard.write(" -")
#        dcard.write(" -")

#        dcard.write('\n')    

    dcard.close()

    if args.lep == "muon":
        shapes = ROOT.TFile.Open("wg_dcard_mu_chan_shapes.root","recreate")        
    elif args.lep == "electron":
        shapes = ROOT.TFile.Open("wg_dcard_el_chan_shapes.root","recreate")
    else:
        assert(0)    

    shapes.cd()

    tmphist.Write("tmphist")

    fake_photon["hists"][mlg_index].SetBinContent(fake_photon["hists"][mlg_index].GetNbinsX()+1,0)
    fake_photon["hists"][mlg_index].SetBinError(fake_photon["hists"][mlg_index].GetNbinsX()+1,0)
    fake_photon["hists"][mlg_index].SetBinContent(0,0)
    fake_photon["hists"][mlg_index].SetBinError(0,0)

    data["hists"][mlg_index].Write("data_obs")
#    labels["wg+jets"]["hists"][mlg_index].Write("wg")
    if "w+jets" in labels:
        labels["w+jets"]["hists-prompt-pileup"][mlg_index].Write("pileup")
    labels["wg+jets"]["hists-pass-fiducial"][mlg_index].Write("wg")
    labels["wg+jets"]["hists-fail-fiducial"][mlg_index].Write("wgout")
    labels["top+jets"]["hists"][mlg_index].Write("topjets")
    labels["zg+jets"]["hists"][mlg_index].Write("zgjets")
    labels["vv+jets"]["hists"][mlg_index].Write("vvjets")
    labels["gg+jets"]["hists"][mlg_index].Write("ggjets")
    for i in range(len(etopbinning)):
        e_to_p[i]["hists"][mlg_index].Write("etopbin"+str(i))
    e_to_p_total["hists"][mlg_index].Write("etoptotal")

#    tmphist=labels["wg+jets"]["hists"][mlg_index].Clone("")
#    tmphist.Scale(fake_photon["hists"][mlg_index].Integral()/tmphist.Integral())
#    tmphist.Write("fakephoton")
    fake_photon["hists"][mlg_index].Write("fakephoton")
    if args.lep == "muon":
        fake_lepton["hists"][mlg_index].Write("fakemuon")
    elif args.lep == "electron":
        fake_lepton["hists"][mlg_index].Write("fakeelectron")
    else:
        assert(0)    

    double_fake["hists"][mlg_index].Write("doublefake")

    zgjets_scale_syst.Write("zgjets_zgscaleUp")
    makeDownShape(zgjets_scale_syst,labels["zg+jets"]["hists"][mlg_index]).Write("zgjets_zgscaleDown")

    zgjets_pdf_syst.Write("zgjets_zgpdfUp")
    makeDownShape(zgjets_pdf_syst,labels["zg+jets"]["hists"][mlg_index]).Write("zgjets_zgpdfDown")

    wgjets_pdf_syst=histogram_models[mlg_index].GetHistogram()

    for i in range(labels["wg+jets"]["hists-pdf-variation0"][mlg_index].GetNbinsX()+1):
        mean_pdf=0

        for j in range(1,32):
            mean_pdf += labels["wg+jets"]["hists-pdf-variation"+str(j)][mlg_index].GetBinContent(i)*labels["wg+jets"]["hists"][mlg_index].Integral()/labels["wg+jets"]["hists-pdf-variation"+str(j)][mlg_index].Integral()

        mean_pdf = mean_pdf/31

        stddev_pdf = 0

        for j in range(1,32):
            stddev_pdf += pow(labels["wg+jets"]["hists-pdf-variation"+str(j)][mlg_index].GetBinContent(i)*labels["wg+jets"]["hists"][mlg_index].Integral()/labels["wg+jets"]["hists-pdf-variation"+str(j)][mlg_index].Integral() - mean_pdf,2)

        stddev_pdf = sqrt(stddev_pdf/(31-1))

        wgjets_pdf_syst.SetBinContent(i,labels["wg+jets"]["hists"][mlg_index].GetBinContent(i)+stddev_pdf)

    wgjets_pass_fiducial_pdf_syst=histogram_models[mlg_index].GetHistogram()

    for i in range(labels["wg+jets"]["hists-pass-fiducial-pdf-variation0"][mlg_index].GetNbinsX()+1):
        mean_pdf=0

        for j in range(1,32):
            mean_pdf += labels["wg+jets"]["hists-pass-fiducial-pdf-variation"+str(j)][mlg_index].GetBinContent(i)*labels["wg+jets"]["hists-pass-fiducial"][mlg_index].Integral()/labels["wg+jets"]["hists-pass-fiducial-pdf-variation"+str(j)][mlg_index].Integral()

        mean_pdf = mean_pdf/31

        stddev_pdf = 0

        for j in range(1,32):
            stddev_pdf += pow(labels["wg+jets"]["hists-pass-fiducial-pdf-variation"+str(j)][mlg_index].GetBinContent(i)*labels["wg+jets"]["hists-pass-fiducial"][mlg_index].Integral()/labels["wg+jets"]["hists-pass-fiducial-pdf-variation"+str(j)][mlg_index].Integral() - mean_pdf,2)

        stddev_pdf = sqrt(stddev_pdf/(31-1))

        wgjets_pass_fiducial_pdf_syst.SetBinContent(i,labels["wg+jets"]["hists-pass-fiducial"][mlg_index].GetBinContent(i)+stddev_pdf)

    wgjets_fail_fiducial_pdf_syst=histogram_models[mlg_index].GetHistogram()

    for i in range(labels["wg+jets"]["hists-pdf-variation0"][mlg_index].GetNbinsX()+1):
        mean_pdf=0

        for j in range(1,32):
            mean_pdf += labels["wg+jets"]["hists-fail-fiducial-pdf-variation"+str(j)][mlg_index].GetBinContent(i)

        mean_pdf = mean_pdf/31

        stddev_pdf = 0

        for j in range(1,32):
            stddev_pdf += pow(labels["wg+jets"]["hists-fail-fiducial-pdf-variation"+str(j)][mlg_index].GetBinContent(i) - mean_pdf,2)

        stddev_pdf = sqrt(stddev_pdf/(31-1))

        wgjets_fail_fiducial_pdf_syst.SetBinContent(i,labels["wg+jets"]["hists-fail-fiducial"][mlg_index].GetBinContent(i)+stddev_pdf)

    wgjets_scale_syst=histogram_models[mlg_index].GetHistogram()

    for i in range(labels["wg+jets"]["hists-scale-variation0"][mlg_index].GetNbinsX()+1):
        wgjets_scale_syst.SetBinContent(i,labels["wg+jets"]["hists"][mlg_index].GetBinContent(i)+labels["wg+jets"]["hists"][mlg_index].Integral()*max(
            abs(labels["wg+jets"]["hists-scale-variation0"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-scale-variation0"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists"][mlg_index].Integral()),
            abs(labels["wg+jets"]["hists-scale-variation1"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-scale-variation1"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists"][mlg_index].Integral()),
            abs(labels["wg+jets"]["hists-scale-variation3"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-scale-variation3"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists"][mlg_index].Integral()),
            abs(labels["wg+jets"]["hists-scale-variation4"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-scale-variation4"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists"][mlg_index].Integral()),
            abs(labels["wg+jets"]["hists-scale-variation5"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-scale-variation5"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists"][mlg_index].Integral()),
            abs(labels["wg+jets"]["hists-scale-variation6"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-scale-variation6"][mlg_index].Integral() - labels["wg+jets"]["hists"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists"][mlg_index].Integral())))

    wgjets_pass_fiducial_scale_syst=histogram_models[mlg_index].GetHistogram()

    for i in range(labels["wg+jets"]["hists-pass-fiducial-scale-variation0"][mlg_index].GetNbinsX()+1):
        wgjets_pass_fiducial_scale_syst.SetBinContent(i,labels["wg+jets"]["hists-pass-fiducial"][mlg_index].GetBinContent(i)+labels["wg+jets"]["hists-pass-fiducial"][mlg_index].Integral()*max(
            abs(labels["wg+jets"]["hists-pass-fiducial-scale-variation0"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-pass-fiducial-scale-variation0"][mlg_index].Integral() - labels["wg+jets"]["hists-pass-fiducial"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-pass-fiducial"][mlg_index].Integral()),
            abs(labels["wg+jets"]["hists-pass-fiducial-scale-variation1"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-pass-fiducial-scale-variation1"][mlg_index].Integral() - labels["wg+jets"]["hists-pass-fiducial"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-pass-fiducial"][mlg_index].Integral()),
            abs(labels["wg+jets"]["hists-pass-fiducial-scale-variation3"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-pass-fiducial-scale-variation3"][mlg_index].Integral() - labels["wg+jets"]["hists-pass-fiducial"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-pass-fiducial"][mlg_index].Integral()),
            abs(labels["wg+jets"]["hists-pass-fiducial-scale-variation4"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-pass-fiducial-scale-variation4"][mlg_index].Integral() - labels["wg+jets"]["hists-pass-fiducial"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-pass-fiducial"][mlg_index].Integral()),
            abs(labels["wg+jets"]["hists-pass-fiducial-scale-variation5"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-pass-fiducial-scale-variation5"][mlg_index].Integral() - labels["wg+jets"]["hists-pass-fiducial"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-pass-fiducial"][mlg_index].Integral()),
            abs(labels["wg+jets"]["hists-pass-fiducial-scale-variation6"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-pass-fiducial-scale-variation6"][mlg_index].Integral() - labels["wg+jets"]["hists-pass-fiducial"][mlg_index].GetBinContent(i)/labels["wg+jets"]["hists-pass-fiducial"][mlg_index].Integral())))

    wgjets_fail_fiducial_scale_syst=histogram_models[mlg_index].GetHistogram()

    for i in range(labels["wg+jets"]["hists-fail-fiducial-scale-variation0"][mlg_index].GetNbinsX()+1):
        wgjets_fail_fiducial_scale_syst.SetBinContent(i,labels["wg+jets"]["hists-fail-fiducial"][mlg_index].GetBinContent(i)+max(
            abs(labels["wg+jets"]["hists-fail-fiducial-scale-variation0"][mlg_index].GetBinContent(i) - labels["wg+jets"]["hists-fail-fiducial"][mlg_index].GetBinContent(i)),
            abs(labels["wg+jets"]["hists-fail-fiducial-scale-variation1"][mlg_index].GetBinContent(i) - labels["wg+jets"]["hists-fail-fiducial"][mlg_index].GetBinContent(i)),
            abs(labels["wg+jets"]["hists-fail-fiducial-scale-variation3"][mlg_index].GetBinContent(i) - labels["wg+jets"]["hists-fail-fiducial"][mlg_index].GetBinContent(i)),
            abs(labels["wg+jets"]["hists-fail-fiducial-scale-variation4"][mlg_index].GetBinContent(i) - labels["wg+jets"]["hists-fail-fiducial"][mlg_index].GetBinContent(i)),
            abs(labels["wg+jets"]["hists-fail-fiducial-scale-variation5"][mlg_index].GetBinContent(i) - labels["wg+jets"]["hists-fail-fiducial"][mlg_index].GetBinContent(i)),
            abs(labels["wg+jets"]["hists-fail-fiducial-scale-variation6"][mlg_index].GetBinContent(i) - labels["wg+jets"]["hists-fail-fiducial"][mlg_index].GetBinContent(i))))

    for i in range(1,labels["wg+jets"]["hists"][mlg_index].GetNbinsX()+1):
        wg_stat_up[i-1].Write("wg_wgstatbin"+str(i)+"Up")
        makeDownShape(wg_stat_up[i-1],labels["wg+jets"]["hists"][mlg_index]).Write("wg_wgstatbin"+str(i)+"Down")

    for i in range(1,len(fake_photon_syst2_up)+1):
        fake_photon_syst2_up[i-1].Write("fakephoton_fakephotonsyst2var"+str(i)+"Up")
#        fake_photon_syst2_up[i-1].Write("fakephoton_fakephotonsyst2var"+str(i)+"Down")
        fake_photon["hists"][mlg_index].Write("fakephoton_fakephotonsyst2var"+str(i)+"Down")
#        makeDownShape(fake_photon_syst2_up[i-1],fake_photon["hists"][mlg_index]).Write("fakephoton_fakephotonsyst2var"+str(i)+"Down")

    for i in range(1,fake_lepton["hists"][mlg_index].GetNbinsX()+1):
        if args.lep == "electron":
            fake_lepton_syst_up[i-1].Write("fakeelectron_fakeelectronsystbin"+str(i)+"Up")
            makeDownShape(fake_lepton_syst_up[i-1],fake_lepton["hists"][mlg_index]).Write("fakeelectron_fakeelectronsystbin"+str(i)+"Down")
        elif args.lep == "muon":
            fake_lepton_syst_up[i-1].Write("fakemuon_fakemuonsystbin"+str(i)+"Up")
            makeDownShape(fake_lepton_syst_up[i-1],fake_lepton["hists"][mlg_index]).Write("fakemuon_fakemuonsystbin"+str(i)+"Down")
        else:
            assert(0)

    for i in range(1,labels["zg+jets"]["hists"][mlg_index].GetNbinsX()+1):
        zgjets_stat_up[i-1].Write("zgjets_zgjetsstatbin"+str(i)+"Up")
        makeDownShape(zgjets_stat_up[i-1],labels["zg+jets"]["hists"][mlg_index]).Write("zgjets_zgjetsstatbin"+str(i)+"Down")

#    wgjets_scale_syst.Write("wg_wgscaleUp")
#    makeDownShape(wgjets_scale_syst,labels["wg+jets"]["hists"][mlg_index]).Write("wg_wgscaleDown")

#    wgjets_pdf_syst.Write("wg_wgpdfUp")
#    makeDownShape(wgjets_pdf_syst,labels["wg+jets"]["hists"][mlg_index]).Write("wg_wgpdfDown")

    if not args.float_sig_fake_cont:

        signal_fake_subtraction_up=labels["wg+jets"]["hists-fake-lepton-pass-fiducial"][mlg_index].Clone("signal fake subtraction up")

        signal_fake_subtraction_up.Add(fake_lepton["hists"][mlg_index])

        if args.lep == "muon":
            signal_fake_subtraction_up.Write("fakemuon_sigfakesubtractionUp")
            makeDownShape(signal_fake_subtraction_up,fake_lepton["hists"][mlg_index]).Write("fakemuon_sigfakesubtractionDown")
        elif args.lep == "electron":
            signal_fake_subtraction_up.Write("fakeelectron_sigfakesubtractionUp")
            makeDownShape(signal_fake_subtraction_up,fake_lepton["hists"][mlg_index]).Write("fakeelectron_sigfakesubtractionDown")
        else:
            assert(0)
    
    wgjets_pass_fiducial_scale_syst.Write("wg_wgscaleUp")
    makeDownShape(wgjets_pass_fiducial_scale_syst,labels["wg+jets"]["hists-pass-fiducial"][mlg_index]).Write("wg_wgscaleDown")

    wgjets_pass_fiducial_pdf_syst.Write("wg_wgpdfUp")
    makeDownShape(wgjets_pass_fiducial_pdf_syst,labels["wg+jets"]["hists-pass-fiducial"][mlg_index]).Write("wg_wgpdfDown")

    wgjets_fail_fiducial_scale_syst.Write("wgout_wgoutscaleUp")
    makeDownShape(wgjets_fail_fiducial_scale_syst,labels["wg+jets"]["hists-fail-fiducial"][mlg_index]).Write("wgout_wgoutscaleDown")

    wgjets_fail_fiducial_pdf_syst.Write("wgout_wgoutpdfUp")
    makeDownShape(wgjets_fail_fiducial_pdf_syst,labels["wg+jets"]["hists-fail-fiducial"][mlg_index]).Write("wgout_wgoutpdfDown")

    for i in range(n_fake_photon_alt):
        fake_photon["hists-alt"+str(i)][mlg_index].Write("fakephoton_fakephotonsyst1var"+str(i)+"Up")
        makeDownShape(fake_photon["hists-alt"+str(i)][mlg_index],fake_photon["hists"][mlg_index]).Write("fakephoton_fakephotonsyst1var"+str(i)+"Down")

    labels["w+jets"]["hists"][mlg_index].Write("fakephoton_fakephotonsyst2Up")
    makeDownShape(labels["w+jets"]["hists"][mlg_index],fake_photon["hists"][mlg_index]).Write("fakephoton_fakephotonsyst2Down")

    for unc in ["pileup","prefire","jes","jer","muon-id-sf","muon-iso-sf","muon-hlt-sf","electron-reco-sf","electron-id-sf","electron-hlt-sf","photon-id-sf"]:
        if unc != "jer" and unc != "jes":
            labels["gg+jets"]["hists-"+unc+"-up"][mlg_index].Write("ggjets_"+unc.replace("-sf","").replace("-","")+"Up")
            makeDownShape(labels["gg+jets"]["hists-"+unc+"-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]).Write("ggjets_"+unc.replace("-sf","").replace("-","")+"Down")
            if "w+jets" in labels:
                labels["w+jets"]["hists-prompt-pileup-"+unc+"-up"][mlg_index].Write("pileup_"+unc.replace("-sf","").replace("-","")+"Up")
                makeDownShape(labels["w+jets"]["hists-prompt-pileup-"+unc+"-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]).Write("pileup_"+unc.replace("-sf","").replace("-","")+"Down")

        labels["wg+jets"]["hists-pass-fiducial-"+unc+"-up"][mlg_index].Write("wg_"+unc.replace("-sf","").replace("-","")+"Up")
        makeDownShape(labels["wg+jets"]["hists-pass-fiducial-"+unc+"-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]).Write("wg_"+unc.replace("-sf","").replace("-","")+"Down")
        labels["wg+jets"]["hists-fail-fiducial-"+unc+"-up"][mlg_index].Write("wgout_"+unc.replace("-sf","").replace("-","")+"Up")
        makeDownShape(labels["wg+jets"]["hists-fail-fiducial-"+unc+"-up"][mlg_index],labels["wg+jets"]["hists-fail-fiducial"][mlg_index]).Write("wgout_"+unc.replace("-sf","").replace("-","")+"Down")

        for label in ["top+jets","vv+jets","zg+jets"]:
            labels[label]["hists-"+unc+"-up"][mlg_index].Write(label.replace("+","")+"_"+unc.replace("-sf","").replace("-","")+"Up")
            makeDownShape(labels[label]["hists-"+unc+"-up"][mlg_index],labels[label]["hists"][mlg_index]).Write(label.replace("+","")+"_"+unc.replace("-sf","").replace("-","")+"Down")

    shapes.Close()

if args.make_unc_table:

    goodbins = lambda hist : filter(lambda i : i > 0, [i*int(abs(hist.GetBinContent(i)) > 10) for i in range(1,hist.GetNbinsX()+1)])

    uncmin = lambda up,nom : 100*(min([abs(up.GetBinContent(i)/nom.GetBinContent(i)-1) for i in goodbins(nom)]))

    uncmax = lambda up,nom : 100*(max([abs(up.GetBinContent(i)/nom.GetBinContent(i)-1) for i in goodbins(nom)]))

    statuncmin = lambda nom : 100*min([nom.GetBinError(i)/nom.GetBinContent(i) for i in goodbins(nom)])

    statuncmax = lambda nom : 100*max([nom.GetBinError(i)/nom.GetBinContent(i) for i in goodbins(nom)])

    fakephotonsyst2uncmin = lambda uplist,nom : 100*(min([abs(uplist[i-1].GetBinContent(i)/nom.GetBinContent(i)-1) for i in goodbins(nom)]))

    fakephotonsyst2uncmax = lambda uplist,nom : 100*(max([abs(uplist[i-1].GetBinContent(i)/nom.GetBinContent(i)-1) for i in goodbins(nom)]))

#fakephotonsyst1uncmin([fake_photon["hists-alt1"][mlg_index],fake_photon["hists"][mlg_index]),
#fakephotonsyst1uncmax([fake_photon["hists-alt1"][mlg_index],fake_photon["hists"][mlg_index]),

    print """\\begin{table}[htbp]
\\begin{center}
\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}
\\hline
   & WG & ZG & top & VV & pileup & fake lepton & fake photon & double fake & e to p   \\\\
\\hline
\\hline"""
    for i in range(n_fake_photon_alt):
        print """fake photon comp 1 variation %i & - & - & - & - & - & - & %0.2f-%0.2f & - & - \\\\
\\hline"""%(i,
            uncmin(fake_photon["hists-alt"+str(i)][mlg_index],fake_photon["hists"][mlg_index]),
            uncmax(fake_photon["hists-alt"+str(i)][mlg_index],fake_photon["hists"][mlg_index])
)
    print """\\end{tabular}
\\end{center}
\\caption{}
\\label{tab:fake_photon_comp_1}
\\end{table}"""

    if args.lep == "electron":

        print """\\begin{table}[htbp]
\\begin{center}
\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}
\\hline
   & WG & ZG & top & VV & pileup & fake lepton & fake photon & double fake & e-induced $\gamma$ & conversion    \\\\
\\hline
\\hline
lumi & 1.8 & 1.8 & 1.8 & 1.8 & 1.8 & - & - & - & - & 1.8 \\\\
\\hline
wgscale & %0.2f-%0.2f & - & - & - & - & - & - & - & - & - \\\\
\\hline
wgpdf & %0.2f-%0.2f & - & - & - & - & - & - & - & - & - \\\\
\\hline
zgscale & - & %0.2f-%0.2f & - & - & - & - & - & - & - & - \\\\
\\hline
zgpdf & - & %0.2f-%0.2f & - & - & - & - & - & - & - & - \\\\
\\hline
fake lepton &  - & - & - & - & - & 30 & - & 30 & - & - \\\\
\\hline
fake photon comp 2 &  - & - & - & - & - & - & %0.2f-%0.2f & - & - & - \\\\
\\hline
pileup & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - & %0.2f-%0.2f \\\\
\\hline
prefire & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & -  & %0.2f-%0.2f \\\\
\\hline
JES & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - & %0.2f-%0.2f & %0.2f-%0.2f & - & - \\\\
\\hline
JER & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - & %0.2f-%0.2f & %0.2f-%0.2f & - & - \\\\
\\hline
stat & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f \\\\
\\hline
photon ID and iso & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - & %0.2f-%0.2f \\\\
\\hline
electron reco & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - & %0.2f-%0.2f  \\\\
\\hline
electron ID and iso & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f  & - & %0.2f-%0.2f  \\\\
\\hline
electron HLT & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - & %0.2f-%0.2f  \\\\
\\hline"""%(
uncmin(wgjets_pass_fiducial_scale_syst,labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(wgjets_pass_fiducial_scale_syst,labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),

uncmin(wgjets_pass_fiducial_pdf_syst,labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(wgjets_pass_fiducial_pdf_syst,labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),

uncmin(zgjets_scale_syst,labels["zg+jets"]["hists"][mlg_index]),
uncmax(zgjets_scale_syst,labels["zg+jets"]["hists"][mlg_index]),

uncmin(zgjets_pdf_syst,labels["zg+jets"]["hists"][mlg_index]),
uncmax(zgjets_pdf_syst,labels["zg+jets"]["hists"][mlg_index]),

fakephotonsyst2uncmin(fake_photon_syst2_up,fake_photon["hists"][mlg_index]),
fakephotonsyst2uncmax(fake_photon_syst2_up,fake_photon["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-pileup-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-pileup-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-pileup-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-pileup-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-pileup-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-pileup-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-pileup-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-pileup-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmin(labels["w+jets"]["hists-prompt-pileup-pileup-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmax(labels["w+jets"]["hists-prompt-pileup-pileup-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-pileup-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-pileup-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-pileup-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-pileup-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-pileup-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-pileup-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmin(labels["gg+jets"]["hists-pileup-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]),
uncmax(labels["gg+jets"]["hists-pileup-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-prefire-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-prefire-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-prefire-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-prefire-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-prefire-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-prefire-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-prefire-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-prefire-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmin(labels["w+jets"]["hists-prompt-pileup-prefire-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmax(labels["w+jets"]["hists-prompt-pileup-prefire-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-prefire-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-prefire-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-prefire-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-prefire-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-prefire-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-prefire-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmin(labels["gg+jets"]["hists-prefire-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]),
uncmax(labels["gg+jets"]["hists-prefire-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-jes-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-jes-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-jes-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-jes-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-jes-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-jes-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-jes-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-jes-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
#uncmin(labels["w+jets"]["hists-prompt-pileup-jes-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
#uncmax(labels["w+jets"]["hists-prompt-pileup-jes-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-jes-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-jes-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-jes-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-jes-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-jes-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-jes-up"][mlg_index],double_fake["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-jer-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-jer-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-jer-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-jer-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-jer-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-jer-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-jer-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-jer-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
#uncmin(labels["w+jets"]["hists-prompt-pileup-jer-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
#uncmax(labels["w+jets"]["hists-prompt-pileup-jer-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-jer-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-jer-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-jer-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-jer-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-jer-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-jer-up"][mlg_index],double_fake["hists"][mlg_index]),

statuncmin(labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
statuncmax(labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
statuncmin(labels["zg+jets"]["hists"][mlg_index]),
statuncmax(labels["zg+jets"]["hists"][mlg_index]),
statuncmin(labels["top+jets"]["hists"][mlg_index]),
statuncmax(labels["top+jets"]["hists"][mlg_index]),
statuncmin(labels["vv+jets"]["hists"][mlg_index]),
statuncmax(labels["vv+jets"]["hists"][mlg_index]),
statuncmin(labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
statuncmax(labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
statuncmin(fake_lepton["hists"][mlg_index]),
statuncmax(fake_lepton["hists"][mlg_index]),
statuncmin(fake_photon["hists"][mlg_index]),
statuncmax(fake_photon["hists"][mlg_index]),
statuncmin(double_fake["hists"][mlg_index]),
statuncmax(double_fake["hists"][mlg_index]),
#statuncmin(e_to_p_non_res["hists"][mlg_index]),
#statuncmax(e_to_p_non_res["hists"][mlg_index])
statuncmin(e_to_p_total["hists"][mlg_index]),
statuncmax(e_to_p_total["hists"][mlg_index]),
statuncmin(labels["gg+jets"]["hists"][mlg_index]),
statuncmax(labels["gg+jets"]["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-photon-id-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-photon-id-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-photon-id-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-photon-id-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-photon-id-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-photon-id-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-photon-id-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-photon-id-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmin(labels["w+jets"]["hists-prompt-pileup-photon-id-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmax(labels["w+jets"]["hists-prompt-pileup-photon-id-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-photon-id-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-photon-id-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-photon-id-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-photon-id-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-photon-id-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-photon-id-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmin(labels["gg+jets"]["hists-photon-id-sf-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]),
uncmax(labels["gg+jets"]["hists-photon-id-sf-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-electron-reco-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-electron-reco-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-electron-reco-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-electron-reco-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-electron-reco-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-electron-reco-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-electron-reco-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-electron-reco-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmin(labels["w+jets"]["hists-prompt-pileup-electron-reco-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmax(labels["w+jets"]["hists-prompt-pileup-electron-reco-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-electron-reco-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-electron-reco-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-electron-reco-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-electron-reco-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-electron-reco-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-electron-reco-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmin(labels["gg+jets"]["hists-electron-reco-sf-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]),
uncmax(labels["gg+jets"]["hists-electron-reco-sf-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-electron-id-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-electron-id-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-electron-id-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-electron-id-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-electron-id-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-electron-id-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-electron-id-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-electron-id-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmin(labels["w+jets"]["hists-prompt-pileup-electron-id-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmax(labels["w+jets"]["hists-prompt-pileup-electron-id-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-electron-id-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-electron-id-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-electron-id-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-electron-id-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-electron-id-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-electron-id-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmin(labels["gg+jets"]["hists-electron-id-sf-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]),
uncmax(labels["gg+jets"]["hists-electron-id-sf-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]),
    
uncmin(labels["wg+jets"]["hists-pass-fiducial-electron-hlt-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-electron-hlt-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-electron-hlt-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-electron-hlt-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-electron-hlt-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-electron-hlt-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-electron-hlt-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-electron-hlt-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmin(labels["w+jets"]["hists-prompt-pileup-electron-hlt-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmax(labels["w+jets"]["hists-prompt-pileup-electron-hlt-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-electron-hlt-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-electron-hlt-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-electron-hlt-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-electron-hlt-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-electron-hlt-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-electron-hlt-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmin(labels["gg+jets"]["hists-electron-hlt-sf-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]),
uncmax(labels["gg+jets"]["hists-electron-hlt-sf-up"][mlg_index],labels["gg+jets"]["hists"][mlg_index]))

    elif args.lep == "muon":

        print """\\begin{table}[htbp]
\\begin{center}
\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}
\\hline
   & WG & ZG & top & VV & pileup & fake lepton & fake photon & double fake & e-induced $\gamma$  \\\\
\\hline
\\hline
lumi & 1.8 & 1.8 & 1.8 & 1.8 & 1.8 & - & - & - & - \\\\
\\hline
wgscale & %0.2f-%0.2f & - & - & - & - & - & - & - & - \\\\
\\hline
wgpdf & %0.2f-%0.2f & - & - & - & - & - & - & - & - \\\\
\\hline
zgscale & - & %0.2f-%0.2f & - & - & - & - & - & - & - \\\\
\\hline
zgpdf & - & %0.2f-%0.2f & - & - & - & - & - & - & - \\\\
\\hline
fake lepton &  - & - & - & - & - & 30 & - & 30 & - \\\\
\\hline
fake photon comp 2 &  - & - & - & - & - & - & %0.2f-%0.2f & - & - \\\\
\\hline
pileup & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - \\\\
\\hline
prefire & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & -  \\\\
\\hline
JES & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - & %0.2f-%0.2f & %0.2f-%0.2f & - \\\\
\\hline
JER & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - & %0.2f-%0.2f & %0.2f-%0.2f & - \\\\
\\hline
stat & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f \\\\
\\hline
photon ID and iso & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - \\\\
\\hline
muon ID & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - \\\\
\\hline
muon iso & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - \\\\
\\hline
muon HLT & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & %0.2f-%0.2f & - \\\\
\\hline"""%(
uncmin(wgjets_pass_fiducial_scale_syst,labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(wgjets_pass_fiducial_scale_syst,labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),

uncmin(wgjets_pass_fiducial_pdf_syst,labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(wgjets_pass_fiducial_pdf_syst,labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),

uncmin(zgjets_scale_syst,labels["zg+jets"]["hists"][mlg_index]),
uncmax(zgjets_scale_syst,labels["zg+jets"]["hists"][mlg_index]),

uncmin(zgjets_pdf_syst,labels["zg+jets"]["hists"][mlg_index]),
uncmax(zgjets_pdf_syst,labels["zg+jets"]["hists"][mlg_index]),

fakephotonsyst2uncmin(fake_photon_syst2_up,fake_photon["hists"][mlg_index]),
fakephotonsyst2uncmax(fake_photon_syst2_up,fake_photon["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-pileup-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-pileup-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-pileup-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-pileup-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-pileup-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-pileup-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-pileup-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-pileup-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmin(labels["w+jets"]["hists-prompt-pileup-pileup-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmax(labels["w+jets"]["hists-prompt-pileup-pileup-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-pileup-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-pileup-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-pileup-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-pileup-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-pileup-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-pileup-up"][mlg_index],double_fake["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-prefire-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-prefire-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-prefire-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-prefire-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-prefire-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-prefire-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-prefire-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-prefire-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmin(labels["w+jets"]["hists-prompt-pileup-prefire-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmax(labels["w+jets"]["hists-prompt-pileup-prefire-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-prefire-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-prefire-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-prefire-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-prefire-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-prefire-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-prefire-up"][mlg_index],double_fake["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-jes-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-jes-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-jes-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-jes-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-jes-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-jes-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-jes-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-jes-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
#uncmin(labels["w+jets"]["hists-prompt-pileup-jes-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
#uncmax(labels["w+jets"]["hists-prompt-pileup-jes-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-jes-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-jes-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-jes-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-jes-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-jes-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-jes-up"][mlg_index],double_fake["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-jer-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-jer-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-jer-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-jer-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-jer-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-jer-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-jer-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-jer-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
#uncmin(labels["w+jets"]["hists-prompt-pileup-jer-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
#uncmax(labels["w+jets"]["hists-prompt-pileup-jer-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-jer-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-jer-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-jer-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-jer-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-jer-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-jer-up"][mlg_index],double_fake["hists"][mlg_index]),

statuncmin(labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
statuncmax(labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
statuncmin(labels["zg+jets"]["hists"][mlg_index]),
statuncmax(labels["zg+jets"]["hists"][mlg_index]),
statuncmin(labels["top+jets"]["hists"][mlg_index]),
statuncmax(labels["top+jets"]["hists"][mlg_index]),
statuncmin(labels["vv+jets"]["hists"][mlg_index]),
statuncmax(labels["vv+jets"]["hists"][mlg_index]),
statuncmin(labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
statuncmax(labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
statuncmin(fake_lepton["hists"][mlg_index]),
statuncmax(fake_lepton["hists"][mlg_index]),
statuncmin(fake_photon["hists"][mlg_index]),
statuncmax(fake_photon["hists"][mlg_index]),
statuncmin(double_fake["hists"][mlg_index]),
statuncmax(double_fake["hists"][mlg_index]),
#statuncmin(e_to_p_non_res["hists"][mlg_index]),
#statuncmax(e_to_p_non_res["hists"][mlg_index])
statuncmin(e_to_p_total["hists"][mlg_index]),
statuncmax(e_to_p_total["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-photon-id-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-photon-id-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-photon-id-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-photon-id-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-photon-id-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-photon-id-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-photon-id-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-photon-id-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmin(labels["w+jets"]["hists-prompt-pileup-photon-id-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmax(labels["w+jets"]["hists-prompt-pileup-photon-id-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-photon-id-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-photon-id-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-photon-id-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-photon-id-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-photon-id-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-photon-id-sf-up"][mlg_index],double_fake["hists"][mlg_index]),

uncmin(labels["wg+jets"]["hists-pass-fiducial-muon-id-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-muon-id-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-muon-id-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-muon-id-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-muon-id-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-muon-id-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-muon-id-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-muon-id-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmin(labels["w+jets"]["hists-prompt-pileup-muon-id-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmax(labels["w+jets"]["hists-prompt-pileup-muon-id-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-muon-id-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-muon-id-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-muon-id-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-muon-id-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-muon-id-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-muon-id-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
   
uncmin(labels["wg+jets"]["hists-pass-fiducial-muon-iso-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-muon-iso-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-muon-iso-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-muon-iso-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-muon-iso-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-muon-iso-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-muon-iso-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-muon-iso-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmin(labels["w+jets"]["hists-prompt-pileup-muon-iso-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmax(labels["w+jets"]["hists-prompt-pileup-muon-iso-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-muon-iso-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-muon-iso-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-muon-iso-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-muon-iso-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-muon-iso-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-muon-iso-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
    
uncmin(labels["wg+jets"]["hists-pass-fiducial-muon-hlt-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmax(labels["wg+jets"]["hists-pass-fiducial-muon-hlt-sf-up"][mlg_index],labels["wg+jets"]["hists-pass-fiducial"][mlg_index]),
uncmin(labels["zg+jets"]["hists-muon-hlt-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmax(labels["zg+jets"]["hists-muon-hlt-sf-up"][mlg_index],labels["zg+jets"]["hists"][mlg_index]),
uncmin(labels["top+jets"]["hists-muon-hlt-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmax(labels["top+jets"]["hists-muon-hlt-sf-up"][mlg_index],labels["top+jets"]["hists"][mlg_index]),
uncmin(labels["vv+jets"]["hists-muon-hlt-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmax(labels["vv+jets"]["hists-muon-hlt-sf-up"][mlg_index],labels["vv+jets"]["hists"][mlg_index]),
uncmin(labels["w+jets"]["hists-prompt-pileup-muon-hlt-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmax(labels["w+jets"]["hists-prompt-pileup-muon-hlt-sf-up"][mlg_index],labels["w+jets"]["hists-prompt-pileup"][mlg_index]),
uncmin(fake_lepton["hists-muon-hlt-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmax(fake_lepton["hists-muon-hlt-sf-up"][mlg_index],fake_lepton["hists"][mlg_index]),
uncmin(fake_photon["hists-muon-hlt-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmax(fake_photon["hists-muon-hlt-sf-up"][mlg_index],fake_photon["hists"][mlg_index]),
uncmin(double_fake["hists-muon-hlt-sf-up"][mlg_index],double_fake["hists"][mlg_index]),
uncmax(double_fake["hists-muon-hlt-sf-up"][mlg_index],double_fake["hists"][mlg_index]))

    else:
        assert(0)

    print """\\end{tabular}
\\end{center}
\\caption{}
\\label{tab:wg_n_sig_unc}
\\end{table}"""

if not args.ewdim6:
    sys.exit(0)

wgjets_ewdim6_scale_syst=histogram_models[ewdim6_index].GetHistogram()

for i in range(labels["wg+jets"]["hists-scale-variation0"][ewdim6_index].GetNbinsX()+1):
    wgjets_ewdim6_scale_syst.SetBinContent(i,labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)+max(
        abs(labels["wg+jets"]["hists-scale-variation0"][ewdim6_index].GetBinContent(i)-labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)),
        abs(labels["wg+jets"]["hists-scale-variation1"][ewdim6_index].GetBinContent(i)-labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)),
        abs(labels["wg+jets"]["hists-scale-variation3"][ewdim6_index].GetBinContent(i)-labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)),
        abs(labels["wg+jets"]["hists-scale-variation4"][ewdim6_index].GetBinContent(i)-labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)),
        abs(labels["wg+jets"]["hists-scale-variation5"][ewdim6_index].GetBinContent(i)-labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)),
        abs(labels["wg+jets"]["hists-scale-variation6"][ewdim6_index].GetBinContent(i)-labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i))))

wgjets_ewdim6_pdf_syst=histogram_models[ewdim6_index].GetHistogram()

for i in range(labels["wg+jets"]["hists-pdf-variation0"][ewdim6_index].GetNbinsX()+1):
    mean_pdf=0

    for j in range(1,32):
        mean_pdf += labels["wg+jets"]["hists-pdf-variation"+str(j)][ewdim6_index].GetBinContent(i)

    mean_pdf = mean_pdf/31

    stddev_pdf = 0

    for j in range(1,32):
        stddev_pdf += pow(labels["wg+jets"]["hists-pdf-variation"+str(j)][ewdim6_index].GetBinContent(i) - mean_pdf,2)

    stddev_pdf = sqrt(stddev_pdf/(31-1))

    wgjets_ewdim6_pdf_syst.SetBinContent(i,labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)+stddev_pdf)

zgjets_ewdim6_scale_syst=histogram_models[ewdim6_index].GetHistogram()

for i in range(labels["zg+jets"]["hists-scale-variation0"][ewdim6_index].GetNbinsX()+1):
    zgjets_ewdim6_scale_syst.SetBinContent(i,labels["zg+jets"]["hists"][ewdim6_index].GetBinContent(i)+max(
        abs(labels["zg+jets"]["hists-scale-variation0"][ewdim6_index].GetBinContent(i)-labels["zg+jets"]["hists"][ewdim6_index].GetBinContent(i)),
        abs(labels["zg+jets"]["hists-scale-variation1"][ewdim6_index].GetBinContent(i)-labels["zg+jets"]["hists"][ewdim6_index].GetBinContent(i)),
        abs(labels["zg+jets"]["hists-scale-variation3"][ewdim6_index].GetBinContent(i)-labels["zg+jets"]["hists"][ewdim6_index].GetBinContent(i)),
        abs(labels["zg+jets"]["hists-scale-variation4"][ewdim6_index].GetBinContent(i)-labels["zg+jets"]["hists"][ewdim6_index].GetBinContent(i)),
        abs(labels["zg+jets"]["hists-scale-variation5"][ewdim6_index].GetBinContent(i)-labels["zg+jets"]["hists"][ewdim6_index].GetBinContent(i)),
        abs(labels["zg+jets"]["hists-scale-variation6"][ewdim6_index].GetBinContent(i)-labels["zg+jets"]["hists"][ewdim6_index].GetBinContent(i))))

zgjets_ewdim6_pdf_syst=histogram_models[ewdim6_index].GetHistogram()

for i in range(labels["zg+jets"]["hists-pdf-variation0"][ewdim6_index].GetNbinsX()+1):
    mean_pdf=0

    for j in range(1,32):
        mean_pdf += labels["zg+jets"]["hists-pdf-variation"+str(j)][ewdim6_index].GetBinContent(i)

    mean_pdf = mean_pdf/31

    stddev_pdf = 0

    for j in range(1,32):
        stddev_pdf += pow(labels["zg+jets"]["hists-pdf-variation"+str(j)][ewdim6_index].GetBinContent(i) - mean_pdf,2)

    stddev_pdf = sqrt(stddev_pdf/(31-1))

    zgjets_ewdim6_pdf_syst.SetBinContent(i,labels["zg+jets"]["hists"][ewdim6_index].GetBinContent(i)+stddev_pdf)

for i in range(1,sm_lhe_weight_hist.GetNbinsX()+1):

    dcard = open("wg_dcard_ewdim6_bin"+str(i)+".txt",'w')

    print >> dcard, "imax 1 number of channels"
    print >> dcard, "jmax * number of background"
    print >> dcard, "kmax * number of nuisance parameters"
    print >> dcard, "Observation "+str(data["hists"][ewdim6_index].GetBinContent(i))
    dcard.write("bin")
    dcard.write(" bin1")
    
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" bin1")

    dcard.write(" bin1")    
    dcard.write(" bin1")    
    dcard.write(" bin1")    
    dcard.write(" bin1")    
    dcard.write('\n')    
    
    dcard.write("process")
    dcard.write(" Wg")
        
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" " + label)

    dcard.write(" fake_photon")
    dcard.write(" fake_lepton")
    dcard.write(" double_fake")
    dcard.write(" e_to_p")
    dcard.write('\n')    
    dcard.write("process")
    dcard.write(" 0")
    
    for j in range(1,len(labels.keys())+2):
        dcard.write(" " + str(j))
    dcard.write('\n')    
    dcard.write('rate')
    dcard.write(' '+str(sm_lhe_weight_hist.GetBinContent(i)))
#    dcard.write(' '+str(labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+ str(labels[label]["hists"][ewdim6_index].GetBinContent(i)))
        else:
            dcard.write(" 0.0001") 



    if fake_photon["hists"][ewdim6_index].GetBinContent(i) > 0:        
        dcard.write(" "+str(fake_photon["hists"][ewdim6_index].GetBinContent(i))) 
    else:
        if fake_photon["hists"][ewdim6_index].GetBinContent(i) < 0:
            print "Warning: fake photon estimate is "+str(fake_photon["hists"][ewdim6_index].GetBinContent(i))+ " for bin " + str(i) + ". It will be replaced with 0.0001"
        dcard.write(" 0.0001") 

    if fake_lepton["hists"][ewdim6_index].GetBinContent(i) > 0:        
        dcard.write(" "+str(fake_lepton["hists"][ewdim6_index].GetBinContent(i))) 
    else:
        if fake_lepton["hists"][ewdim6_index].GetBinContent(i) < 0:
            print "Warning: fake lepton estimate is "+str(fake_lepton["hists"][ewdim6_index].GetBinContent(i))+ " for bin " + str(i) + ". It will be replaced with 0.0001"
        dcard.write(" 0.0001") 

    if double_fake["hists"][ewdim6_index].GetBinContent(i) > 0:        
        dcard.write(" "+str(double_fake["hists"][ewdim6_index].GetBinContent(i))) 
    else:
        if double_fake["hists"][ewdim6_index].GetBinContent(i) < 0:
            print "Warning: double fake estimate is "+str(double_fake["hists"][ewdim6_index].GetBinContent(i))+ " for bin " + str(i) + ". It will be replaced with 0.0001"
        dcard.write(" 0.0001") 

    if e_to_p_total["hists"][ewdim6_index].GetBinContent(i) > 0:        
        dcard.write(" "+str(e_to_p_total["hists"][ewdim6_index].GetBinContent(i))) 
    else:
        dcard.write(" 0.0001") 
   
    dcard.write('\n')    

    dcard.write("lumi lnN")
    dcard.write(" 1.018")

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" 1.018")

    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" -")
    dcard.write(" 1.018")

    dcard.write('\n')    

    if sm_lhe_weight_hist.GetBinContent(i) > 0:
        dcard.write("mcstatewdim6bin"+str(i)+" lnN "+str(1+sm_lhe_weight_hist.GetBinError(i)/sm_lhe_weight_hist.GetBinContent(i)))
        
#        dcard.write("mcstat_ewdim6_bin"+str(i)+" lnN "+str(1+labels["wg+jets"]["hists"][ewdim6_index].GetBinError(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" -")

        dcard.write(" -")                
        dcard.write(" -")                
        dcard.write(" -")                
        dcard.write(" -")                
        dcard.write("\n")  

    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue

        if labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write("mcstat"+str(label.replace("+",""))+"bin"+str(i)+" lnN ")
            dcard.write(" -")

            for l in labels.keys():
                if l == "no label" or l == "wg+jets" or label == "w+jets":
                    continue
                if l == label:
                    dcard.write(" "+str(1+labels[label]["hists"][ewdim6_index].GetBinError(i)/labels[label]["hists"][ewdim6_index].GetBinContent(i)))
                else:    
                    dcard.write(" -")

            dcard.write(" -")                
            dcard.write(" -")                
            dcard.write(" -")                
            dcard.write(" -")                
            dcard.write("\n")  

    if fake_lepton["hists"][ewdim6_index].GetBinContent(i) > 0:        
        dcard.write("fakeleptonsyst lnN -")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" -")

        dcard.write(" -")                
        dcard.write(" 1.3")                
        dcard.write(" -")                
        dcard.write(" -")                
        dcard.write("\n")  

    if fake_lepton["hists"][ewdim6_index].GetBinContent(i) > 0:        
        dcard.write("fakeleptonstatbin"+str(i)+" lnN -")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets":
                continue
            dcard.write(" -")

        dcard.write(" -")                
        dcard.write(" "+str(1+fake_lepton["hists"][ewdim6_index].GetBinError(i)/fake_lepton["hists"][ewdim6_index].GetBinContent(i)))
        dcard.write(" -")                
        dcard.write(" -")                
        dcard.write("\n")  

    if fake_photon["hists"][ewdim6_index].GetBinContent(i) > 0:        
        dcard.write("fakephotonstatbin"+str(i)+" lnN -")
        for label in labels.keys():
            if label == "no label" or label == "wg+jets" or label == "w+jets" :
                continue
            dcard.write(" -")

        dcard.write(" "+str(1+fake_photon["hists"][ewdim6_index].GetBinError(i)/fake_photon["hists"][ewdim6_index].GetBinContent(i)))
        dcard.write(" -")                
        dcard.write(" -")                
        dcard.write(" -")                
        dcard.write("\n")  

    if fake_photon["hists"][ewdim6_index].GetBinContent(i) > 0:
        for j in range(n_fake_photon_alt):
            dcard.write("fakephotonsyst1var"+str(j)+" lnN -")
            for label in labels.keys():
                if label == "no label" or label == "wg+jets" or label == "w+jets" :
                    continue
                dcard.write(" -")

            dcard.write(" "+str(fake_photon["hists-alt"+str(j)][ewdim6_index].GetBinContent(i)/fake_photon["hists"][ewdim6_index].GetBinContent(i)))
            dcard.write(" -")
            dcard.write(" -")
            dcard.write(" -")
            dcard.write("\n")

    dcard.write("muonid lnN "+str(labels["wg+jets"]["hists-muon-id-sf-up"][ewdim6_index].GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(labels[label]["hists-muon-id-sf-up"][ewdim6_index].GetBinContent(i)/labels[label]["hists"][ewdim6_index].GetBinContent(i)))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")                
    dcard.write(" -")            
    dcard.write("\n")  

    dcard.write("muonhlt lnN "+str(labels["wg+jets"]["hists-muon-hlt-sf-up"][ewdim6_index].GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(labels[label]["hists-muon-hlt-sf-up"][ewdim6_index].GetBinContent(i)/labels[label]["hists"][ewdim6_index].GetBinContent(i)))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")                
    dcard.write(" -")            
    dcard.write("\n")  

    dcard.write("muoniso lnN "+str(labels["wg+jets"]["hists-muon-iso-sf-up"][ewdim6_index].GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(labels[label]["hists-muon-iso-sf-up"][ewdim6_index].GetBinContent(i)/labels[label]["hists"][ewdim6_index].GetBinContent(i)))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")                
    dcard.write(" -")            
    dcard.write("\n")  

    dcard.write("electronhlt lnN "+str(labels["wg+jets"]["hists-electron-hlt-sf-up"][ewdim6_index].GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(labels[label]["hists-electron-hlt-sf-up"][ewdim6_index].GetBinContent(i)/labels[label]["hists"][ewdim6_index].GetBinContent(i)))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")                
    if e_to_p_total["hists"][ewdim6_index].GetBinContent(i) > 0:
        dcard.write(" "+str(e_to_p_total["hists-electron-hlt-sf-up"][ewdim6_index].GetBinContent(i)/e_to_p_total["hists"][ewdim6_index].GetBinContent(i)))
    else:    
        dcard.write(" -")            
    dcard.write("\n")  

    dcard.write("electronid lnN "+str(labels["wg+jets"]["hists-electron-id-sf-up"][ewdim6_index].GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(labels[label]["hists-electron-id-sf-up"][ewdim6_index].GetBinContent(i)/labels[label]["hists"][ewdim6_index].GetBinContent(i)))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")                
    if e_to_p_total["hists"][ewdim6_index].GetBinContent(i) > 0:
        dcard.write(" "+str(e_to_p_total["hists-electron-id-sf-up"][ewdim6_index].GetBinContent(i)/e_to_p_total["hists"][ewdim6_index].GetBinContent(i)))
    else:    
        dcard.write(" -")            
    dcard.write("\n")  

    dcard.write("electronreco lnN "+str(labels["wg+jets"]["hists-electron-reco-sf-up"][ewdim6_index].GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(labels[label]["hists-electron-reco-sf-up"][ewdim6_index].GetBinContent(i)/labels[label]["hists"][ewdim6_index].GetBinContent(i)))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")                
    if e_to_p_total["hists"][ewdim6_index].GetBinContent(i) > 0:
        dcard.write(" "+str(e_to_p_total["hists-electron-reco-sf-up"][ewdim6_index].GetBinContent(i)/e_to_p_total["hists"][ewdim6_index].GetBinContent(i)))
    else:    
        dcard.write(" -")            
    dcard.write("\n")  

    dcard.write("photonid lnN "+str(labels["wg+jets"]["hists-photon-id-sf-up"][ewdim6_index].GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(labels[label]["hists-photon-id-sf-up"][ewdim6_index].GetBinContent(i)/labels[label]["hists"][ewdim6_index].GetBinContent(i)))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")                
    if e_to_p_total["hists"][ewdim6_index].GetBinContent(i) > 0:
        dcard.write(" "+str(e_to_p_total["hists-photon-id-sf-up"][ewdim6_index].GetBinContent(i)/e_to_p_total["hists"][ewdim6_index].GetBinContent(i)))
    else:    
        dcard.write(" -")            
    dcard.write("\n")  

    dcard.write("jes lnN "+str(labels["wg+jets"]["hists-jes-up"][ewdim6_index].GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if label != "gg+jets" and labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(labels[label]["hists-jes-up"][ewdim6_index].GetBinContent(i)/labels[label]["hists"][ewdim6_index].GetBinContent(i)))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")                
    if e_to_p_total["hists"][ewdim6_index].GetBinContent(i) > 0:
        dcard.write(" "+str(e_to_p_total["hists-jes-up"][ewdim6_index].GetBinContent(i)/e_to_p_total["hists"][ewdim6_index].GetBinContent(i)))
    else:    
        dcard.write(" -")            
    dcard.write("\n")  

    dcard.write("jer lnN "+str(labels["wg+jets"]["hists-jer-up"][ewdim6_index].GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if label != "gg+jets" and labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(labels[label]["hists-jer-up"][ewdim6_index].GetBinContent(i)/labels[label]["hists"][ewdim6_index].GetBinContent(i)))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")                
    if e_to_p_total["hists"][ewdim6_index].GetBinContent(i) > 0:
        dcard.write(" "+str(e_to_p_total["hists-jer-up"][ewdim6_index].GetBinContent(i)/e_to_p_total["hists"][ewdim6_index].GetBinContent(i)))
    else:    
        dcard.write(" -")            
    dcard.write("\n")  

    dcard.write("prefire lnN "+str(labels["wg+jets"]["hists-prefire-up"][ewdim6_index].GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(labels[label]["hists-prefire-up"][ewdim6_index].GetBinContent(i)/labels[label]["hists"][ewdim6_index].GetBinContent(i)))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")            
    if e_to_p_total["hists"][ewdim6_index].GetBinContent(i) > 0:
        dcard.write(" "+str(e_to_p_total["hists-prefire-up"][ewdim6_index].GetBinContent(i)/e_to_p_total["hists"][ewdim6_index].GetBinContent(i)))
    else:    
        dcard.write(" -")            
    dcard.write("\n")  

    dcard.write("pileup lnN "+str(labels["wg+jets"]["hists-pileup-up"][ewdim6_index].GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        if labels[label]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(labels[label]["hists-pileup-up"][ewdim6_index].GetBinContent(i)/labels[label]["hists"][ewdim6_index].GetBinContent(i)))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")        
    if e_to_p_total["hists"][ewdim6_index].GetBinContent(i) > 0:
        dcard.write(" "+str(e_to_p_total["hists-pileup-up"][ewdim6_index].GetBinContent(i)/e_to_p_total["hists"][ewdim6_index].GetBinContent(i)))
    else:    
        dcard.write(" -")        
    dcard.write("\n")  

    dcard.write("wgscale lnN "+str(wgjets_ewdim6_scale_syst.GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")        
    dcard.write(" -")        
    dcard.write("\n")  

    dcard.write("wgpdf lnN "+str(wgjets_ewdim6_pdf_syst.GetBinContent(i)/labels["wg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue
        dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")        
    dcard.write(" -")        
    dcard.write("\n")  

    dcard.write("zgscale lnN -")
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue

        if label == "zg+jets" and labels["zg+jets"]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(zgjets_ewdim6_scale_syst.GetBinContent(i)/labels["zg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
        else:    
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")        
    dcard.write(" -")        
    dcard.write("\n")  

    dcard.write("zgpdf lnN -")
    for label in labels.keys():
        if label == "no label" or label == "wg+jets" or label == "w+jets":
            continue

        if label == "zg+jets" and labels["zg+jets"]["hists"][ewdim6_index].GetBinContent(i) > 0:
            dcard.write(" "+str(zgjets_ewdim6_pdf_syst.GetBinContent(i)/labels["zg+jets"]["hists"][ewdim6_index].GetBinContent(i)))
        else:        
            dcard.write(" -")

    dcard.write(" -")
    dcard.write(" -")                
    dcard.write(" -")        
    dcard.write(" -")        
    dcard.write("\n")  
