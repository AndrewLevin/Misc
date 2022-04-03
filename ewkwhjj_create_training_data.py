from coffea.jetmet_tools import FactorizedJetCorrector, JetCorrectionUncertainty
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory

import awkward as ak
import numpy as np

import uproot as uproot
from coffea.nanoevents import NanoEventsFactory, BaseSchema, NanoAODSchema

import awkward as ak
from coffea import hist, processor

from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

from coffea.lumi_tools import LumiMask

import numba

import math
import numpy as np

import xgboost as xgb

import pandas

import argparse

import pprint

parser = argparse.ArgumentParser()

parser.add_argument('--year',dest='year',default='2016')
parser.add_argument('--nproc',dest='nproc',type=int,default='10')

args = parser.parse_args()

pprint.pprint(vars(args))

assert(args.year == '2016' or args.year == '2017' or args.year == '2018')

year = args.year

from coffea.lookup_tools import extractor

ext = extractor()

if year == '2016':
    ext.add_weight_sets(['electronrecosf EGamma_SF2D /afs/cern.ch/user/a/amlevin/ewkwhjj/data/EGM2D_UL2016preVFP.root'])
    ext.add_weight_sets(['electronrecosfunc EGamma_SF2D_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/EGM2D_UL2016preVFP.root'])
    ext.add_weight_sets(['electronidsf EGamma_SF2D /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Ele_Medium_preVFP_EGM2D.root'])
    ext.add_weight_sets(['electronidsfunc EGamma_SF2D_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Ele_Medium_preVFP_EGM2D.root'])
    ext.add_weight_sets(['muonidsf NUM_TightID_DEN_TrackerMuons_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root'])
    ext.add_weight_sets(['muonidsfunc NUM_TightID_DEN_TrackerMuons_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root'])
    ext.add_weight_sets(['muonisosf NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root'])
    ext.add_weight_sets(['muonisosfunc NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root'])
    ext.add_weight_sets(['muonhltsf NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root'])
    ext.add_weight_sets(['muonhltsfunc NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016pre/Summer19UL16APV_V7_MC_L1FastJet_AK4PFchs.jec.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016pre/Summer19UL16APV_V7_MC_L2Relative_AK4PFchs.jec.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016pre/Summer19UL16APV_V7_MC_L3Absolute_AK4PFchs.jec.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016pre/Summer19UL16APV_V7_MC_L2L3Residual_AK4PFchs.jec.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016pre/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.jr.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016pre/Summer19UL16APV_V7_MC_Uncertainty_AK4PFchs.junc.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016pre/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.jersf.txt'])
#    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016post/Summer19UL16_V7_MC_L1FastJet_AK4PFchs.jec.txt'])
#    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016post/Summer19UL16_V7_MC_L2Relative_AK4PFchs.jec.txt'])
#    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016post/Summer19UL16_V7_MC_L3Absolute_AK4PFchs.jec.txt'])
#    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016post/Summer19UL16_V7_MC_L2L3Residual_AK4PFchs.jec.txt'])
#    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016post/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.jr.txt'])
#    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016post/Summer19UL16_V7_MC_Uncertainty_AK4PFchs.junc.txt'])
#    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2016post/Summer20UL16_JRV3_MC_SF_AK4PFchs.jersf.txt'])

elif year == '2017':
    ext.add_weight_sets(['electronrecosf EGamma_SF2D /afs/cern.ch/user/a/amlevin/ewkwhjj/data/EGM2D_UL2017.root'])
    ext.add_weight_sets(['electronrecosfunc EGamma_SF2D_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/EGM2D_UL2017.root'])
    ext.add_weight_sets(['electronidsf EGamma_SF2D /afs/cern.ch/user/a/amlevin/ewkwhjj/data/EGM2D_Medium_UL17.root'])
    ext.add_weight_sets(['electronidsfunc EGamma_SF2D_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/EGM2D_Medium_UL17.root'])
    ext.add_weight_sets(['muonidsf NUM_TightID_DEN_TrackerMuons_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root'])
    ext.add_weight_sets(['muonidsfunc NUM_TightID_DEN_TrackerMuons_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root'])
    ext.add_weight_sets(['muonisosf NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root'])
    ext.add_weight_sets(['muonisosfunc NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root'])
    ext.add_weight_sets(['muonhltsf NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root'])
    ext.add_weight_sets(['muonhltsfunc NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2017/Summer19UL17_V5_MC_L1FastJet_AK4PFchs.jec.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2017/Summer19UL17_V5_MC_L2Residual_AK4PFchs.jec.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2017/Summer19UL17_V5_MC_L3Absolute_AK4PFchs.jec.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2017/Summer19UL17_V5_MC_L2L3Residual_AK4PFchs.jec.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2017/Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.jr.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2017/Summer19UL17_V5_MC_Uncertainty_AK4PFchs.junc.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2017/Summer19UL17_JRV2_MC_SF_AK4PFchs.jersf.txt'])

elif year == '2018':
    ext.add_weight_sets(['electronrecosf EGamma_SF2D /afs/cern.ch/user/a/amlevin/ewkwhjj/data/EGM2D_UL2018.root'])
    ext.add_weight_sets(['electronrecosfunc EGamma_SF2D_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/EGM2D_UL2018.root'])
    ext.add_weight_sets(['electronidsf EGamma_SF2D /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Ele_Medium_EGM2D.root'])
    ext.add_weight_sets(['electronidsfunc EGamma_SF2D_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Ele_Medium_EGM2D.root'])
    ext.add_weight_sets(['muonidsf NUM_TightID_DEN_TrackerMuons_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root'])
    ext.add_weight_sets(['muonidsfunc NUM_TightID_DEN_TrackerMuons_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root'])
    ext.add_weight_sets(['muonisosf NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root'])
    ext.add_weight_sets(['muonisosfunc NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root'])
    ext.add_weight_sets(['muonhltsf NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root'])
    ext.add_weight_sets(['muonhltsfunc NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2018/Summer19UL18_V5_MC_L1FastJet_AK4PFchs.jec.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2018/Summer19UL18_V5_MC_L2Relative_AK4PFchs.jec.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2018/Summer19UL18_V5_MC_L3Absolute_AK4PFchs.jec.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2018/Summer19UL18_V5_MC_L2L3Residual_AK4PFchs.jec.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2018/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.jr.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2018/Summer19UL18_V5_MC_Uncertainty_AK4PFchs.junc.txt'])
    ext.add_weight_sets(['* * /afs/cern.ch/user/a/amlevin/ewkwhjj/data/2018/Summer19UL18_JRV2_MC_SF_AK4PFchs.jersf.txt'])


ext.add_weight_sets(['pileup ratio_{} /afs/cern.ch/user/a/amlevin/ewkwhjj/data/pileup.root'.format(year)])
ext.add_weight_sets(['pileup_up ratio_{}_up /afs/cern.ch/user/a/amlevin/ewkwhjj/data/pileup.root'.format(year)])
ext.add_weight_sets(['pileup_down ratio_{}_up /afs/cern.ch/user/a/amlevin/ewkwhjj/data/pileup.root'.format(year)])

ext.finalize()

evaluator = ext.make_evaluator()

if year == '2016':
    jec_stack_names = ["Summer19UL16APV_V7_MC_L1FastJet_AK4PFchs","Summer19UL16APV_V7_MC_L2Relative_AK4PFchs","Summer19UL16APV_V7_MC_L3Absolute_AK4PFchs","Summer19UL16APV_V7_MC_L2L3Residual_AK4PFchs","Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs","Summer19UL16APV_V7_MC_Uncertainty_AK4PFchs","Summer20UL16APV_JRV3_MC_SF_AK4PFchs"]
#    jec_stack_names = ["Summer19UL16_V7_MC_L1FastJet_AK4PFchs","Summer19UL16_V7_MC_L2Relative_AK4PFchs","Summer19UL16_V7_MC_L3Absolute_AK4PFchs","Summer19UL16_V7_MC_L2L3Residual_AK4PFchs","Summer20UL16_JRV3_MC_PtResolution_AK4PFchs","Summer19UL16_V7_MC_Uncertainty_AK4PFchs","Summer20UL16_JRV3_MC_SF_AK4PFchs"]
elif year == '2017':
    jec_stack_names = ["Summer19UL17_V5_MC_L1FastJet_AK4PFchs","Summer19UL17_V5_MC_L2Residual_AK4PFchs","Summer19UL17_V5_MC_L3Absolute_AK4PFchs","Summer19UL17_V5_MC_L2L3Residual_AK4PFchs","Summer19UL17_JRV2_MC_PtResolution_AK4PFchs","Summer19UL17_V5_MC_Uncertainty_AK4PFchs","Summer19UL17_JRV2_MC_SF_AK4PFchs"]
elif year == '2018':
    jec_stack_names = ["Summer19UL18_V5_MC_L1FastJet_AK4PFchs","Summer19UL18_V5_MC_L2Relative_AK4PFchs","Summer19UL18_V5_MC_L3Absolute_AK4PFchs","Summer19UL18_V5_MC_L2L3Residual_AK4PFchs","Summer19UL18_JRV2_MC_PtResolution_AK4PFchs","Summer19UL18_V5_MC_Uncertainty_AK4PFchs","Summer19UL18_JRV2_MC_SF_AK4PFchs"]

jec_inputs = {name: evaluator[name] for name in jec_stack_names}


jec_stack = JECStack(jec_inputs)

@numba.njit
def deltar(eta1,phi1,eta2,phi2):
    dphi = (phi1 - phi2 + 3.14) % (2 * 3.14) - 3.14
    return math.sqrt((eta1 - eta2) ** 2 + dphi ** 2)

@numba.njit
def select_event_merged(muon_pt,muon_eta,muon_phi,muon_tightid,muon_pfreliso04all,electron_pt,electron_eta,electron_phi,electron_deltaetasc,electron_dxy,electron_dz,electron_cutbased,jet_pt,jet_eta,jet_phi,jet_btagdeepb,fatjet_pt,fatjet_eta,fatjet_phi,fatjet_msoftdrop,met_pt,dataset,builder):

    builder.begin_list() 

    if met_pt < 30:
        builder.end_list()
        return

    tight_muons = []
    loose_not_tight_muons = []
            
    for i1 in range(len(muon_pt)):
        if muon_tightid[i1]==True and muon_pfreliso04all[i1] < 0.15 and muon_pt[i1] > 26 and abs(muon_eta[i1]) < 2.4:
            tight_muons.append(i1)
        elif muon_tightid[i1]==True and muon_pfreliso04all[i1] < 0.4 and muon_pt[i1] > 26 and abs(muon_eta[i1]) < 2.4:   
            loose_not_tight_muons.append(i1)


    tight_electrons = [] 
            
    for i1 in range(len(electron_pt)):
        if electron_pt[i1] > 30 and abs(electron_eta[i1] + electron_deltaetasc[i1]) < 2.5:
            if (abs(electron_eta[i1] + electron_deltaetasc[i1]) < 1.479 and abs(electron_dz[i1]) < 0.1 and abs(electron_dxy[i1]) < 0.05) or (abs(electron_eta[i1] + electron_deltaetasc[i1]) > 1.479 and abs(electron_dz[i1]) < 0.2 and abs(electron_dxy[i1]) < 0.1):
                if electron_cutbased[i1] >= 3:
                    tight_electrons.append(i1)


    cleaned_fatjets = []

    for i1 in range(len(fatjet_pt)):
        found = False

        for i2 in range(len(tight_muons)):
            if deltar(fatjet_eta[i1],fatjet_phi[i1],muon_eta[tight_muons[i2]],muon_phi[tight_muons[i2]]) < 0.5:
                found = True
                    
        for i2 in range(len(loose_not_tight_muons)):
            if deltar(fatjet_eta[i1],fatjet_phi[i1],muon_eta[loose_not_tight_muons[i2]],muon_phi[loose_not_tight_muons[i2]]) < 0.5:
                found = True

        for i2 in range(len(tight_electrons)):
            if deltar(fatjet_eta[i1],fatjet_phi[i1],electron_eta[tight_electrons[i2]],electron_phi[tight_electrons[i2]]) < 0.5:
                found = True

        if not found:        
            cleaned_fatjets.append(i1)

    cleaned_jets = []                  

    for i1 in range(len(jet_pt)):
        found = False
        
        for i2 in range(len(tight_muons)):
            if deltar(jet_eta[i1],jet_phi[i1],muon_eta[tight_muons[i2]],muon_phi[tight_muons[i2]]) < 0.5:
                found = True
                    
        for i2 in range(len(loose_not_tight_muons)):
            if deltar(jet_eta[i1],jet_phi[i1],muon_eta[loose_not_tight_muons[i2]],muon_phi[loose_not_tight_muons[i2]]) < 0.5:
                found = True

        for i2 in range(len(tight_electrons)):
            if deltar(jet_eta[i1],jet_phi[i1],electron_eta[tight_electrons[i2]],electron_phi[tight_electrons[i2]]) < 0.5:
                found = True

        for i2 in range(len(cleaned_fatjets)):
            if deltar(jet_eta[i1],jet_phi[i1],fatjet_eta[cleaned_fatjets[i2]],fatjet_phi[cleaned_fatjets[i2]]) < 0.5:
                found = True

        if not found:        
            cleaned_jets.append(i1)


    if len(cleaned_jets) < 2 or len(cleaned_fatjets) < 1 or len(tight_muons) + len(tight_electrons) != 1 or len(loose_not_tight_muons) != 0:
        builder.end_list()
        return

    found = False   

    for i1 in range(len(fatjet_pt)):
        
        if found:
            break

        for i2 in range(len(cleaned_jets)):

            if found:
                break

            for i3 in range(len(cleaned_jets)):

                if found:
                    break

                if i2 == i3:
                    continue

                if fatjet_pt[cleaned_fatjets[i1]] < 250 or abs(fatjet_eta[cleaned_fatjets[i1]]) < 2.5 or fatjet_msoftdrop[cleaned_fatjets[i1]] < 50 or fatjet_msoftdrop[cleaned_fatjets[i1]] > 150:
                    continue

                if jet_btagdeepb[cleaned_jets[i2]] > 0.2217 or jet_btagdeepb[cleaned_jets[i3]] > 0.2217 or jet_pt[cleaned_jets[i2]] < 30 or jet_pt[cleaned_jets[i3]] < 30 or abs(jet_eta[cleaned_jets[i2]]) > 4.7 or abs(jet_eta[cleaned_jets[i3]]) > 4.7:
                    continue

            found = True
            
            builder.begin_tuple(7)
            builder.index(0).integer(cleaned_fatjets[i1])
            builder.index(1).integer(cleaned_jets[i2])
            builder.index(2).integer(cleaned_jets[i3])

            if len(tight_muons) > 0:
                builder.index(3).integer(tight_muons[0])
                builder.index(4).integer(-1)
            elif len(tight_electrons) > 0:
                builder.index(3).integer(-1)
                builder.index(4).integer(tight_electrons[0])

            nextrajets=0
            nextrabjets=0
            for i4 in range(len(cleaned_jets)):
                if i4 == i2 or i4 == i3:
                    continue

                if jet_pt[cleaned_jets[i4]] < 30 or abs(jet_eta[cleaned_jets[i4]]) > 4.7:
                    continue

                nextrajets=nextrajets+1

                if jet_btagdeepb[cleaned_jets[i4]] > 0.2217:
                    nextrabjets += 1

            builder.index(5).integer(nextrajets)
            builder.index(6).integer(nextrabjets)

            builder.end_tuple()    
    builder.end_list()                  

@numba.njit
def select_event_resolved(muon_pt,muon_eta,muon_phi,muon_tightid,muon_pfreliso04all,electron_pt,electron_eta,electron_phi,electron_deltaetasc,electron_dxy,electron_dz,electron_cutbased,jet_pt,jet_eta,jet_phi,jet_btagdeepb,met_pt,dataset,builder):

    builder.begin_list() 

    if met_pt < 30:
        builder.end_list()
        return

    found = False

    tight_muons = []
    loose_not_tight_muons = []
            
    for i1 in range(len(muon_pt)):
        if muon_tightid[i1]==True and muon_pfreliso04all[i1] < 0.15 and muon_pt[i1] > 26 and abs(muon_eta[i1]) < 2.4:
            tight_muons.append(i1)
        elif muon_tightid[i1]==True and muon_pfreliso04all[i1] < 0.4 and muon_pt[i1] > 26 and abs(muon_eta[i1]) < 2.4:   
            loose_not_tight_muons.append(i1)
                
    tight_electrons = [] 
            
    for i1 in range(len(electron_pt)):
        if electron_pt[i1] > 30 and abs(electron_eta[i1] + electron_deltaetasc[i1]) < 2.5:
            if (abs(electron_eta[i1] + electron_deltaetasc[i1]) < 1.479 and abs(electron_dz[i1]) < 0.1 and abs(electron_dxy[i1]) < 0.05) or (abs(electron_eta[i1] + electron_deltaetasc[i1]) > 1.479 and abs(electron_dz[i1]) < 0.2 and abs(electron_dxy[i1]) < 0.1):
                if electron_cutbased[i1] >= 3:
                    tight_electrons.append(i1)
                             
    cleaned_jets = []                  

    for i1 in range(len(jet_pt)):

        found = False

        for i2 in range(len(tight_muons)):
            if deltar(jet_eta[i1],jet_phi[i1],muon_eta[tight_muons[i2]],muon_phi[tight_muons[i2]]) < 0.5:
                found = True

        for i2 in range(len(loose_not_tight_muons)):
            if deltar(jet_eta[i1],jet_phi[i1],muon_eta[loose_not_tight_muons[i2]],muon_phi[loose_not_tight_muons[i2]]) < 0.5:
                found = True

        for i2 in range(len(tight_electrons)):
            if deltar(jet_eta[i1],jet_phi[i1],electron_eta[tight_electrons[i2]],electron_phi[tight_electrons[i2]]) < 0.5:
                found = True

        if not found:        
            cleaned_jets.append(i1)



#        if len(cleaned_jets) < 4 or len(tight_muons) + len(loose_not_tight_muons) + len(tight_electrons) != 1:
    if len(cleaned_jets) < 4 or len(tight_muons) + len(tight_electrons) != 1 or len(loose_not_tight_muons) != 0:
        builder.end_list()
        return

    found = False   
            
    for i1 in range(len(cleaned_jets)):

        if found:
            break
                
        for i2 in range(i1+1,len(cleaned_jets)):

            if found:
                break

            for i3 in range(len(cleaned_jets)):

                if found:
                    break

                for i4 in range(i3+1,len(cleaned_jets)):
                        
                    if found:
                        break

                    if i1 == i2 or i1 == i3 or i1 == i4 or i2 == i3 or i2 == i4 or i3 == i4:
                        continue

                    if jet_btagdeepb[cleaned_jets[i1]] < 0.8953 or jet_btagdeepb[cleaned_jets[i2]] < 0.8953 or jet_pt[cleaned_jets[i1]] < 30 or jet_pt[cleaned_jets[i2]] < 30 or abs(jet_eta[cleaned_jets[i1]]) > 2.5 or abs(jet_eta[cleaned_jets[i2]]) > 2.5:
                        continue
                            
                    if jet_btagdeepb[cleaned_jets[i3]] > 0.2217 or jet_btagdeepb[cleaned_jets[i4]] > 0.2217 or jet_pt[cleaned_jets[i3]] < 30 or jet_pt[cleaned_jets[i4]] < 30 or abs(jet_eta[cleaned_jets[i3]]) > 4.7 or abs(jet_eta[cleaned_jets[i4]]) > 4.7:
                        continue
                                
                    found = True
                        
                    builder.begin_tuple(8)
                    builder.index(0).integer(cleaned_jets[i1])
                    builder.index(1).integer(cleaned_jets[i2])
                    builder.index(2).integer(cleaned_jets[i3])
                    builder.index(3).integer(cleaned_jets[i4])

                    if len(tight_muons) > 0:
                        builder.index(4).integer(tight_muons[0])
                        builder.index(5).integer(-1)
                    elif len(tight_electrons) > 0:
                        builder.index(4).integer(-1)
                        builder.index(5).integer(tight_electrons[0])

                    nextrajets=0
                    nextrabjets=0
                    for i5 in range(len(cleaned_jets)):
                        if i5 == i1 or i5 == i2 or i5 == i3 or i5 == i4:
                            continue
                            
                        if jet_pt[cleaned_jets[i5]] < 30 or abs(jet_eta[cleaned_jets[i2]]) > 4.7:
                            continue

                        nextrajets=nextrajets+1

                        if jet_btagdeepb[cleaned_jets[i5]] > 0.2217:
                            nextrabjets += 1

                    builder.index(6).integer(nextrajets)
                    builder.index(7).integer(nextrabjets)
                    
                    builder.end_tuple()    

    builder.end_list()                  
                        
@numba.njit
def select_events(hlt_isotkmu24,hlt_isomu24,hlt_ele27wptightgsf,hlt_isomu27,hlt_ele32wptightgsfl1doubleeg,hlt_ele32wptightgsf,muon_pt,muon_eta,muon_phi,muon_tightid,muon_pfreliso04all,electron_pt,electron_eta,electron_phi,electron_deltaetasc,electron_dxy,electron_dz,electron_cutbased,jet_pt,jet_eta,jet_phi,jet_btagdeepb,fatjet_pt,fatjet_eta,fatjet_phi,fatjet_msoftdrop,met_pt,dataset,builder_merged,builder_resolved):

    for i0 in range(len(muon_pt)):

        pass_hlt = True

        if dataset == 'singlemuon':
            if year == "2016":
                if not hlt_isotkmu24[i0] and not hlt_isomu24[i0]:
                    pass_hlt = False
            elif year == "2017":
                if not hlt_isomu27[i0]:
                    pass_hlt = False
            elif year == "2018":
                if not hlt_isomu24[i0]:
                    pass_hlt = False
        elif dataset == 'singleelectron':
            if year == "2016":
                if hlt_isotkmu24[i0] or hlt_isomu24[i0] or not hlt_ele27wptightgsf[i0]:
                    pass_hlt = False
            elif year == "2017":
                if hlt_isomu27[i0] or not hlt_ele32wptightgsfl1doubleeg[i0]:
                    pass_hlt = False
            elif year == "2018":
                if hlt_isomu24[i0] or not hlt_ele32wptightgsf[i0]:
                    pass_hlt = False
        else: #MC
            if year == "2016":
                if not hlt_isotkmu24[i0] and not hlt_isomu24[i0] and not hlt_ele27wptightgsf[i0]:
                    pass_hlt = False
            elif year == "2017":
                if not hlt_isomu27[i0] and not hlt_ele32wptightgsfl1doubleeg[i0]:
                    pass_hlt = False
            elif year == "2018":
                if not hlt_ele32wptightgsf[i0] and not hlt_isomu24[i0]:
                    pass_hlt = False

        if not pass_hlt:
            builder_merged.begin_list()                     
            builder_merged.end_list()                     
            builder_resolved.begin_list()
            builder_resolved.end_list()
            continue

        select_event_merged(muon_pt[i0],muon_eta[i0],muon_phi[i0],muon_tightid[i0],muon_pfreliso04all[i0],electron_pt[i0],electron_eta[i0],electron_phi[i0],electron_deltaetasc[i0],electron_dxy[i0],electron_dz[i0],electron_cutbased[i0],jet_pt[i0],jet_eta[i0],jet_phi[i0],jet_btagdeepb[i0],fatjet_pt[i0],fatjet_eta[i0],fatjet_phi[i0],fatjet_msoftdrop[i0],met_pt[i0],dataset,builder_merged)

        select_event_resolved(muon_pt[i0],muon_eta[i0],muon_phi[i0],muon_tightid[i0],muon_pfreliso04all[i0],electron_pt[i0],electron_eta[i0],electron_phi[i0],electron_deltaetasc[i0],electron_dxy[i0],electron_dz[i0],electron_cutbased[i0],jet_pt[i0],jet_eta[i0],jet_phi[i0],jet_btagdeepb[i0],met_pt[i0],dataset,builder_resolved)

    return [builder_resolved,builder_merged]
                        
class EwkwhjjProcessor(processor.ProcessorABC):
    def __init__(self):
        self._accumulator = processor.dict_accumulator({
            'sumw': processor.defaultdict_accumulator(float),
            'nevents': processor.defaultdict_accumulator(float),
            'variables': processor.defaultdict_accumulator(processor.column_accumulator(np.transpose(np.array([[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]))).identity),
            'variables_merged':  processor.defaultdict_accumulator(processor.column_accumulator(np.transpose(np.array([[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]))).identity),
            'weights':  processor.defaultdict_accumulator(processor.column_accumulator(np.array([])).identity),
            'weights_merged':  processor.defaultdict_accumulator(processor.column_accumulator(np.array([])).identity),
        })

    @property
    def accumulator(self):
        return self._accumulator

#    @numba.njit
    def process(self, events):

        output = self.accumulator.identity()

        dataset = events.metadata['dataset']

        output['sumw'][dataset] += ak.sum(np.sign(events.Generator.weight))
        output['nevents'][dataset] += len(events)

        name_map = jec_stack.blank_name_map
        name_map['JetPt'] = 'pt'
        name_map['JetMass'] = 'mass'
        name_map['JetEta'] = 'eta'
        name_map['JetA'] = 'area'
        
        jets = events.Jet
        
        jets['pt_raw'] = (1 - jets['rawFactor']) * jets['pt']
        jets['mass_raw'] = (1 - jets['rawFactor']) * jets['mass']
        jets['pt_gen'] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
        jets['rho'] = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, jets.pt)[0]
        name_map['ptGenJet'] = 'pt_gen'
        name_map['ptRaw'] = 'pt_raw'
        name_map['massRaw'] = 'mass_raw'
        name_map['Rho'] = 'rho'
        
        events_cache = events.caches[0]
        
        jet_factory = CorrectedJetsFactory(name_map, jec_stack)
        corrected_jets = jet_factory.build(jets, lazy_cache=events_cache)    

        jet_pt = corrected_jets.pt

        if year == "2016":
            hlt_isotkmu24 = events.HLT.IsoTkMu24
            hlt_isomu24 = events.HLT.IsoMu24
            hlt_ele27wptightgsf = events.HLT.Ele27_WPTight_Gsf
            hlt_isomu27 = ak.Array(len(events)*[False])
            hlt_ele32wptightgsfl1doubleeg = ak.Array(len(events)*[False])
            hlt_ele32wptightgsf = ak.Array(len(events)*[False])
        elif year == "2017":
            hlt_isotkmu24 = ak.Array(len(events)*[False])
            hlt_isomu24 = ak.Array(len(events)*[False])
            hlt_ele27wptightgsf = ak.Array(len(events)*[False])
            hlt_isomu27 = events.HLT.IsoMu27
            hlt_ele32wptightgsfl1doubleeg = events.HLT.Ele32_WPTight_Gsf_L1DoubleEG
            hlt_ele32wptightgsf = ak.Array(len(events)*[False])
        elif year == "2018":
            hlt_isotkmu24 = ak.Array(len(events)*[False])
            hlt_isomu24 = events.HLT.IsoMu24
            hlt_ele27wptightgsf = ak.Array(len(events)*[False])
            hlt_isomu27 = ak.Array(len(events)*[False])
            hlt_ele32wptightgsfl1doubleeg = ak.Array(len(events)*[False])
            hlt_ele32wptightgsf = events.HLT.Ele32_WPTight_Gsf

        builder_resolved,builder_merged = select_events(ak.without_parameters(ak.Array(hlt_isotkmu24)),ak.without_parameters(ak.Array(hlt_isomu24)),ak.without_parameters(ak.Array(hlt_ele27wptightgsf)),ak.without_parameters(ak.Array(hlt_isomu27)),ak.without_parameters(ak.Array(hlt_ele32wptightgsfl1doubleeg)),ak.without_parameters(ak.Array(hlt_ele32wptightgsf)),ak.without_parameters(ak.Array(events.Muon.pt)),ak.without_parameters(ak.Array(events.Muon.eta)),ak.without_parameters(ak.Array(events.Muon.phi)),ak.without_parameters(ak.Array(events.Muon.tightId)),ak.without_parameters(ak.Array(events.Muon.pfRelIso04_all)),ak.without_parameters(ak.Array(events.Electron.pt)),ak.without_parameters(ak.Array(events.Electron.eta)),ak.without_parameters(ak.Array(events.Electron.phi)),ak.without_parameters(ak.Array(events.Electron.deltaEtaSC)),ak.without_parameters(ak.Array(events.Electron.dxy)),ak.without_parameters(ak.Array(events.Electron.dz)),ak.without_parameters(ak.Array(events.Electron.cutBased)),ak.without_parameters(ak.Array(events.Jet.pt)),ak.without_parameters(ak.Array(events.Jet.eta)),ak.without_parameters(ak.Array(events.Jet.phi)),ak.without_parameters(ak.Array(events.Jet.btagDeepB)),ak.without_parameters(ak.Array(events.FatJet.pt)),ak.without_parameters(ak.Array(events.FatJet.eta)),ak.without_parameters(ak.Array(events.FatJet.phi)),ak.without_parameters(ak.Array(events.FatJet.msoftdrop)),ak.without_parameters(ak.Array(events.PuppiMET.pt)),dataset,ak.ArrayBuilder(),ak.ArrayBuilder())

        particleindices = builder_resolved.snapshot()
        particleindices_merged = builder_merged.snapshot()

        basecut = ak.num(particleindices) != 0

        basecut_merged = ak.num(particleindices_merged) != 0

        if ak.any(basecut_merged):
            particleindices_merged = particleindices_merged[basecut_merged]
            events_merged = events[basecut_merged]
            events_merged.FatJet[particleindices_merged['0']]
            jets_merged = [events_merged.Jet[particleindices_merged[idx]] for idx in '12']
            cut7 = ak.firsts((events_merged.FatJet[particleindices_merged['0']].mass > 50) & (events_merged.FatJet[particleindices_merged['0']].mass < 150) & ((jets_merged[0]+jets_merged[1]).mass > 500) & (abs(jets_merged[0].eta - jets_merged[1].eta) > 2.5) & (particleindices_merged['3'] != -1))
            cut8 = ak.firsts((events_merged.FatJet[particleindices_merged['0']].mass > 50) & (events_merged.FatJet[particleindices_merged['0']].mass < 150) & ((jets_merged[0]+jets_merged[1]).mass > 500) & (abs(jets_merged[0].eta - jets_merged[1].eta) > 2.5) & (particleindices_merged['4'] != -1))
#            cut9_merged = cut7_merged | cut8_merged

        if ak.any(basecut):
            particleindices = particleindices[basecut]
            events = events[basecut]
            jets = [events.Jet[particleindices[idx]] for idx in '0123']
            cut1 = ak.firsts(((jets[0]+jets[1]).mass > 50) & ((jets[0]+jets[1]).mass < 150) & ((jets[2]+jets[3]).mass > 500) & (abs(jets[2].eta - jets[3].eta) > 2.5) & (particleindices['4'] != -1))
            cut2 = ak.firsts(((jets[0]+jets[1]).mass > 50) & ((jets[0]+jets[1]).mass < 150) & ((jets[2]+jets[3]).mass > 500) & (abs(jets[2].eta - jets[3].eta) > 2.5) & (particleindices['5'] != -1))
#            cut3 = cut1 | cut2
            cut4 = ak.firsts(((jets[0]+jets[1]).mass > 50) & ((jets[0]+jets[1]).mass < 150) & ((jets[2]+jets[3]).mass > 100) & ((jets[2]+jets[3]).mass < 500) & (particleindices['4'] != -1))
            cut5 = ak.firsts(((jets[0]+jets[1]).mass > 50) & ((jets[0]+jets[1]).mass < 150) & ((jets[2]+jets[3]).mass > 100) & ((jets[2]+jets[3]).mass < 500) & (particleindices['5'] != -1))
#            cut6 = cut4 | cut5

        if ak.any(basecut_merged) and ak.any(cut7):

            sel7_particleindices = particleindices_merged[cut7]
        
            sel7_events = events_merged[cut7]
            
            sel7_fatjets = sel7_events.FatJet[sel7_particleindices['0']]

            sel7_jets = [sel7_events.Jet[sel7_particleindices[idx]] for idx in '12']

            sel7_muons = sel7_events.Muon[sel7_particleindices['3']]

            output["weights_merged"][dataset] += processor.column_accumulator(np.sign(ak.to_numpy(sel7_events.Generator.weight).data))

            output['variables_merged'][dataset] += processor.column_accumulator(np.transpose(np.vstack((
                ak.to_numpy(ak.firsts(sel7_fatjets).pt).data,
                ak.to_numpy(ak.firsts(sel7_fatjets).eta).data,
                ak.to_numpy(ak.firsts(sel7_fatjets).phi).data,
                ak.to_numpy(ak.firsts(sel7_fatjets).btagDeepB).data,
                ak.to_numpy(ak.firsts(sel7_fatjets).btagHbb).data,
                ak.to_numpy(ak.firsts(sel7_fatjets).msoftdrop).data,
                ak.to_numpy(ak.firsts(sel7_particleindices['5'])).data,
                ak.to_numpy(ak.firsts(sel7_particleindices['6'])).data,
                np.zeros(len(sel7_events)),
                np.sign(ak.to_numpy(ak.firsts(sel7_muons).charge).data+1),
                ak.to_numpy(ak.firsts(sel7_muons).pt).data,
                ak.to_numpy(ak.firsts(sel7_muons).eta).data,
                ak.to_numpy(ak.firsts(sel7_muons).phi).data,
                ak.to_numpy(sel7_events.PuppiMET.pt),
                ak.to_numpy(sel7_events.PuppiMET.phi),
                ak.to_numpy(ak.firsts(sel7_jets[0]).pt).data,
                ak.to_numpy(ak.firsts(sel7_jets[1]).pt).data,
                ak.to_numpy(ak.firsts(sel7_jets[0]).eta).data,
                ak.to_numpy(ak.firsts(sel7_jets[1]).eta).data,
                ak.to_numpy(ak.firsts(sel7_jets[0]).phi).data,
                ak.to_numpy(ak.firsts(sel7_jets[1]).phi).data,
                ak.to_numpy(ak.firsts(sel7_jets[0]).btagDeepB).data,
                ak.to_numpy(ak.firsts(sel7_jets[1]).btagDeepB).data,
                ak.to_numpy(ak.firsts((sel7_jets[0]+sel7_jets[1]).mass)).data,
                ak.to_numpy(ak.firsts(sel7_jets[0]).eta - ak.firsts(sel7_jets[1]).eta).data,
                ak.to_numpy(ak.firsts(np.sqrt(2*(sel7_muons+sel7_jets[0]).pt*sel7_events.PuppiMET.pt*(1 - np.cos(sel7_events.PuppiMET.phi - (sel7_muons+sel7_jets[0]).phi))))).data,
                ak.to_numpy(ak.firsts(np.sqrt(2*(sel7_muons+sel7_jets[1]).pt*sel7_events.PuppiMET.pt*(1 - np.cos(sel7_events.PuppiMET.phi - (sel7_muons+sel7_jets[1]).phi))))).data))))

            sel7_muonidsf = ak.firsts(evaluator['muonidsf'](abs(sel7_muons.eta), sel7_muons.pt))
            sel7_muonisosf = ak.firsts(evaluator['muonisosf'](abs(sel7_muons.eta), sel7_muons.pt))
            sel7_muonhltsf = ak.firsts(evaluator['muonhltsf'](abs(sel7_muons.eta), sel7_muons.pt))
            sel7_weight = np.sign(sel7_events.Generator.weight)*sel7_events.L1PreFiringWeight.Nom*sel7_muonidsf*sel7_muonisosf*sel7_muonhltsf

        if ak.any(basecut_merged) and ak.any(cut8):

            sel8_particleindices = particleindices_merged[cut8]
        
            sel8_events = events_merged[cut8]
            
            sel8_fatjets = sel8_events.FatJet[sel8_particleindices['0']]

            sel8_jets = [sel8_events.Jet[sel8_particleindices[idx]] for idx in '12']

            sel8_electrons = sel8_events.Electron[sel8_particleindices['4']]

            output["weights_merged"][dataset] += processor.column_accumulator(np.sign(ak.to_numpy(sel8_events.Generator.weight).data))

            output['variables_merged'][dataset] += processor.column_accumulator(np.transpose(np.vstack((
                ak.to_numpy(ak.firsts(sel8_fatjets).pt).data,
                ak.to_numpy(ak.firsts(sel8_fatjets).eta).data,
                ak.to_numpy(ak.firsts(sel8_fatjets).phi).data,
                ak.to_numpy(ak.firsts(sel8_fatjets).btagDeepB).data,
                ak.to_numpy(ak.firsts(sel8_fatjets).btagHbb).data,
                ak.to_numpy(ak.firsts(sel8_fatjets).msoftdrop).data,
                ak.to_numpy(ak.firsts(sel8_particleindices['5'])).data,
                ak.to_numpy(ak.firsts(sel8_particleindices['6'])).data,
                np.ones(len(sel8_events)),
                np.sign(ak.to_numpy(ak.firsts(sel8_electrons).charge).data+1),
                ak.to_numpy(ak.firsts(sel8_electrons).pt).data,
                ak.to_numpy(ak.firsts(sel8_electrons).eta).data,
                ak.to_numpy(ak.firsts(sel8_electrons).phi).data,
                ak.to_numpy(sel8_events.PuppiMET.pt),
                ak.to_numpy(sel8_events.PuppiMET.phi),
                ak.to_numpy(ak.firsts(sel8_jets[0]).pt).data,
                ak.to_numpy(ak.firsts(sel8_jets[1]).pt).data,
                ak.to_numpy(ak.firsts(sel8_jets[0]).eta).data,
                ak.to_numpy(ak.firsts(sel8_jets[1]).eta).data,
                ak.to_numpy(ak.firsts(sel8_jets[0]).phi).data,
                ak.to_numpy(ak.firsts(sel8_jets[1]).phi).data,
                ak.to_numpy(ak.firsts(sel8_jets[0]).btagDeepB).data,
                ak.to_numpy(ak.firsts(sel8_jets[1]).btagDeepB).data,
                ak.to_numpy(ak.firsts((sel8_jets[0]+sel8_jets[1]).mass)).data,
                ak.to_numpy(ak.firsts(sel8_jets[0]).eta - ak.firsts(sel8_jets[1]).eta).data,
                ak.to_numpy(ak.firsts(np.sqrt(2*(sel8_electrons+sel8_jets[0]).pt*sel8_events.PuppiMET.pt*(1 - np.cos(sel8_events.PuppiMET.phi - (sel8_electrons+sel8_jets[0]).phi))))).data,
                ak.to_numpy(ak.firsts(np.sqrt(2*(sel8_electrons+sel8_jets[1]).pt*sel8_events.PuppiMET.pt*(1 - np.cos(sel8_events.PuppiMET.phi - (sel8_electrons+sel8_jets[1]).phi))))).data))))

            sel8_electronidsf = ak.firsts(evaluator['electronidsf'](sel8_electrons.eta, sel8_electrons.pt))
            sel8_electronrecosf = ak.firsts(evaluator['electronrecosf'](sel8_electrons.eta, sel8_electrons.pt))

            sel8_weight = np.sign(sel8_events.Generator.weight)*sel8_events.L1PreFiringWeight.Nom*sel8_electronidsf*sel8_electronrecosf

        if ak.any(basecut) and ak.any(cut1):

            sel1_particleindices = particleindices[cut1]
        
            sel1_events = events[cut1]
            
            sel1_jets = [sel1_events.Jet[sel1_particleindices[idx]] for idx in '0123']

            sel1_muons = sel1_events.Muon[sel1_particleindices['4']]

            output["weights"][dataset] += processor.column_accumulator(np.sign(ak.to_numpy(sel1_events.Generator.weight).data))

            output['variables'][dataset] += processor.column_accumulator(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel1_particleindices['6'])).data,ak.to_numpy(ak.firsts(sel1_particleindices['7'])).data,np.zeros(len(sel1_events)),np.sign(ak.to_numpy(ak.firsts(sel1_muons).charge).data+1),ak.to_numpy(ak.firsts(sel1_muons).pt).data,ak.to_numpy(ak.firsts(sel1_muons).eta).data,ak.to_numpy(ak.firsts(sel1_muons).phi).data,ak.to_numpy(sel1_events.PuppiMET.pt),ak.to_numpy(sel1_events.PuppiMET.phi),ak.to_numpy(ak.firsts(sel1_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel1_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel1_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel1_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel1_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel1_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel1_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel1_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel1_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel1_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel1_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel1_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel1_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel1_jets[0]+sel1_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel1_jets[2]+sel1_jets[3]).mass)).data, ak.to_numpy(ak.firsts(sel1_jets[2]).eta - ak.firsts(sel1_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_muons+sel1_jets[0]).pt*sel1_events.PuppiMET.pt*(1 - np.cos(sel1_events.PuppiMET.phi - (sel1_muons+sel1_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_muons+sel1_jets[1]).pt*sel1_events.PuppiMET.pt*(1 - np.cos(sel1_events.PuppiMET.phi - (sel1_muons+sel1_jets[1]).phi))))).data))))

            if dataset == 'ewkwhjj_reweighted':
                sel1_pu_weight = evaluator['pileup'](sel1_events.Pileup.nTrueInt)
                sel1_muonidsf = ak.firsts(evaluator['muonidsf'](abs(sel1_muons.eta), sel1_muons.pt))
                sel1_muonisosf = ak.firsts(evaluator['muonisosf'](abs(sel1_muons.eta), sel1_muons.pt))
                sel1_muonhltsf = ak.firsts(evaluator['muonhltsf'](abs(sel1_muons.eta), sel1_muons.pt))
                sel1_weight = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_events.LHEReweightingWeight[:,9]
            else:
                sel1_pu_weight = evaluator['pileup'](sel1_events.Pileup.nTrueInt)
                sel1_muonidsf = ak.firsts(evaluator['muonidsf'](abs(sel1_muons.eta), sel1_muons.pt))
                sel1_muonisosf = ak.firsts(evaluator['muonisosf'](abs(sel1_muons.eta), sel1_muons.pt))
                sel1_muonhltsf = ak.firsts(evaluator['muonhltsf'](abs(sel1_muons.eta), sel1_muons.pt))
                sel1_weight = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf

        if ak.any(basecut) and ak.any(cut2):
                
            sel2_particleindices = particleindices[cut2]
        
            sel2_events = events[cut2]
            
            sel2_jets = [sel2_events.Jet[sel2_particleindices[idx]] for idx in '0123']
            sel2_electrons = sel2_events.Electron[sel2_particleindices['4']]

            output["weights"][dataset] += processor.column_accumulator(np.sign(ak.to_numpy(sel2_events.Generator.weight).data))

            output['variables'][dataset] += processor.column_accumulator(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel2_particleindices['6'])).data,ak.to_numpy(ak.firsts(sel2_particleindices['7'])).data,np.ones(len(sel2_events)),np.sign(ak.to_numpy(ak.firsts(sel2_electrons).charge).data+1),ak.to_numpy(ak.firsts(sel2_electrons).pt).data,ak.to_numpy(ak.firsts(sel2_electrons).eta).data,ak.to_numpy(ak.firsts(sel2_electrons).phi).data,ak.to_numpy(sel2_events.PuppiMET.pt),ak.to_numpy(sel2_events.PuppiMET.phi),ak.to_numpy(ak.firsts(sel2_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel2_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel2_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel2_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel2_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel2_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel2_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel2_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel2_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel2_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel2_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel2_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel2_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel2_jets[0]+sel2_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel2_jets[2]+sel2_jets[3]).mass)).data,ak.to_numpy(ak.firsts(sel2_jets[2]).eta - ak.firsts(sel2_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_electrons+sel2_jets[0]).pt*sel2_events.PuppiMET.pt*(1 - np.cos(sel2_events.PuppiMET.phi - (sel2_electrons+sel2_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_electrons+sel2_jets[1]).pt*sel2_events.PuppiMET.pt*(1 - np.cos(sel2_events.PuppiMET.phi - (sel2_electrons+sel2_jets[1]).phi))))).data))))

            if dataset == 'ewkwhjj_reweighted':
                sel2_pu_weight = evaluator['pileup'](sel2_events.Pileup.nTrueInt)
                sel2_electronidsf = ak.firsts(evaluator['electronidsf'](sel2_electrons.eta, sel2_electrons.pt))
                sel2_electronrecosf = ak.firsts(evaluator['electronrecosf'](sel2_electrons.eta, sel2_electrons.pt))
                sel2_weight = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf*sel2_electronrecosf*sel2_events.LHEReweightingWeight[:,9]
            else:
                sel2_pu_weight = evaluator['pileup'](sel2_events.Pileup.nTrueInt)
                sel2_electronidsf = ak.firsts(evaluator['electronidsf'](sel2_electrons.eta, sel2_electrons.pt))
                sel2_electronrecosf = ak.firsts(evaluator['electronrecosf'](sel2_electrons.eta, sel2_electrons.pt))
                sel2_weight = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf*sel2_electronrecosf

        return output
        
    def postprocess(self, accumulator):
        return accumulator

if year == '2016':
    filelists = {
#        'singlemuon' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/singlemuon.txt',
#        'singleelectron' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/singleelectron.txt',
#        'w' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/w.txt',
#        'ww' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/ww.txt',
        'ewkwhjj': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/ewkwhjj.txt',
        'qcdwhjj': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/qcdwhjj.txt',
        'ttsemi': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/ttsemi.txt',
#        'tthad': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/tthad.txt'
    }
elif year == '2017':
    filelists = {
#        'singlemuon' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2017/singlemuon.txt',
#        'singleelectron' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2017/singleelectron.txt',
        'ttsemi' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2017/ttsemi.txt'
    }
elif year == '2018':
    filelists = {
#        'singlemuon' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/singlemuon.txt',
#        'egamma' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/egamma.txt',
#        'ewkwhjj_reweighted': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/ewkwhjj_reweighted.txt',
        'ewkwhjj': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/ewkwhjj_reweighted.txt',
        'ttsemi': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/ttsemi.txt',
#        'tthad': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/tthad.txt',
#        'qcdwphjj': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/qcdwphjj.txt',
#        'qcdwmhjj': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/qcdwmhjj.txt',
#        'qcdwph': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/qcdwph.txt',
#        'qcdwmh': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/qcdwmh.txt',
#        'wlep' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/wlep.txt',
#        'wlep2j': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/wlep2j.txt',

    }
else:
    assert(0)

samples = {}

for filelist in filelists:
    f = open(filelists[filelist])
    samples[filelist] = f.read().rstrip('\n').split('\n')

"""
    
result = processor.run_uproot_job(
    samples,
    'Events',
    EwkwhjjProcessor(),
    processor.iterative_executor,
    {'schema': NanoAODSchema},
    chunksize=10000000,
)

"""

result = processor.run_uproot_job(
    samples,
    'Events',
    EwkwhjjProcessor(),
    processor.futures_executor,
    {'schema': NanoAODSchema, 'workers': args.nproc},
    chunksize=10000000,
)

for key in result['nevents'].keys():
    print('result[\'nevents\'][\'{}\'] = {}'.format(key,result['nevents'][key]))

from coffea.util import save

#save(result,'training_data_{}'.format(args.year))

#[result['variables'][key] for key in result['variables'].keys()]

#save(pandas.DataFrame(np.hstack(

#['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt']
#key="ww"

#pandas.concat([pandas.DataFrame(f["variables"][key].value.data,columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','met','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt']),pandas.DataFrame(f["weights"][key].value,columns=["weight"]),

#pandas.DataFrame(len(f["variables"][key].value.data)*[key],columns=["label"])

#],axis=1)

save(pandas.concat(
    [
        pandas.concat(
            [
                pandas.DataFrame(result["variables"][key].value.data,columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','met','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt']),
                pandas.DataFrame(result["weights"][key].value,columns=["weight"]),
                pandas.DataFrame(len(result["variables"][key].value.data)*[key],columns=["label"])
         ],axis=1)
        for key in result['variables'].keys()
    ]
),"ewkwhjj_training_data_{}".format(args.year))

save(pandas.concat(
    [
        pandas.concat(
            [
                pandas.DataFrame(result["variables_merged"][key].value.data,columns=['higgsjetpt','higgsjeteta','higgsjetphi','higgsjetbtag1','higgsjetbtag2','higgsjetsoftdropmass','nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','met','metphi','vbsjet1pt','vbsjet2pt','vbsjet1eta','vbsjet2eta','vbsjet1phi','vbsjet2phi','vbsjet1btag','vbsjet2btag','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt']),
                pandas.DataFrame(result["weights_merged"][key].value,columns=["weight"]),
                pandas.DataFrame(len(result["variables_merged"][key].value.data)*[key],columns=["label"])
            ],axis=1)
        for key in result['variables'].keys()
    ]
),"ewkwhjj_merged_training_data_{}".format(args.year))
