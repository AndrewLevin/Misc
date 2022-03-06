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

bst = xgb.Booster({'nthread': 1})

bst.load_model('/afs/cern.ch/user/a/amlevin/ewkwhjj/model.bin')

bst_merged = xgb.Booster({'nthread': 1})

bst_merged.load_model('/afs/cern.ch/user/a/amlevin/ewkwhjj/merged.model')

lumimask = LumiMask("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt")

from coffea.lookup_tools import extractor

ext = extractor()
ext.add_weight_sets(["muonidsf NUM_TightID_DEN_TrackerMuons_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root"])
ext.add_weight_sets(["muonidsfunc NUM_TightID_DEN_TrackerMuons_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root"])
ext.add_weight_sets(["muonisosf NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root"])
ext.add_weight_sets(["muonisosfunc NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root"])
ext.add_weight_sets(["muonhltsf NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root"])
ext.add_weight_sets(["muonhltsfunc NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root"])
ext.add_weight_sets(["electronidsf EGamma_SF2D /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Ele_Medium_preVFP_EGM2D.root"])
ext.add_weight_sets(["electronidsfunc EGamma_SF2D_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Ele_Medium_preVFP_EGM2D.root"])
ext.finalize()

evaluator = ext.make_evaluator()

@numba.njit
def delta_r(a,b):
    dphi = (a.phi - b.phi + 3.14) % (2 * 3.14) - 3.14
    return math.sqrt((a.eta - b.eta) ** 2 + dphi ** 2)

@numba.njit
def select_events_merged (events,dataset,builder,syst="nominal"):

    for i0 in range(len(events)):
        
        builder.begin_list() 
        
        if dataset == "SingleMuon":
            if not events[i0].HLT.IsoTkMu24 and not events[i0].HLT.IsoMu24:
                builder.end_list()
                continue
        elif dataset == "SingleElectron":
            if events[i0].HLT.IsoTkMu24 or events[i0].HLT.IsoMu24 or not events[i0].HLT.Ele27_WPTight_Gsf :
                builder.end_list()
                continue
        else: #MC
            if not events[i0].HLT.IsoTkMu24 and not events[i0].HLT.IsoMu24 and not events[i0].HLT.Ele27_WPTight_Gsf:
                builder.end_list()
                continue
        
        if syst == "nominal":        
            if events[i0].PuppiMET.pt < 30:
                builder.end_list()
                continue
        elif syst == "JESUp":
            if events[i0].PuppiMET.ptJESUp < 30:
                builder.end_list()
                continue
        elif syst == "JERUp":
            if events[i0].PuppiMET.ptJERUp < 30:
                builder.end_list()
                continue
        else:
            builder.end_list()
            continue

        found = False

        tight_muons = []
        loose_not_tight_muons = []
            
        for i1 in range(len(events[i0].Muon)):
            if events[i0].Muon[i1].tightId==True and events[i0].Muon[i1].pfRelIso04_all < 0.15 and events[i0].Muon[i1].pt > 26 and abs(events[i0].Muon[i1].eta) < 2.4:
                tight_muons.append(i1)
            elif events[i0].Muon[i1].tightId==True and events[i0].Muon[i1].pfRelIso04_all < 0.4 and events[i0].Muon[i1].pt > 26 and abs(events[i0].Muon[i1].eta) < 2.4:   
                loose_not_tight_muons.append(i1)
                
        tight_electrons = [] 
            
        for i1 in range(len(events[i0].Electron)):
            if events[i0].Electron[i1].pt > 30 and abs(events[i0].Electron[i1].eta + events[i0].Electron[i1].deltaEtaSC) < 2.5:
                if (abs(events[i0].Electron[i1].eta + events[i0].Electron[i1].deltaEtaSC) < 1.479 and abs(events[i0].Electron[i1].dz) < 0.1 and abs(events[i0].Electron[i1].dxy) < 0.05) or (abs(events[i0].Electron[i1].eta + events[i0].Electron[i1].deltaEtaSC) > 1.479 and abs(events[i0].Electron[i1].dz) < 0.2 and abs(events[i0].Electron[i1].dxy) < 0.1):
                    if events[i0].Electron[i1].cutBased >= 3:
                        tight_electrons.append(i1)

        cleaned_fatjets = []

        for i1 in range(len(events[i0].FatJet)):
            found = False

            for i2 in range(len(tight_muons)):
                if delta_r(events[i0].FatJet[i1],events[i0].Muon[tight_muons[i2]]) < 0.5:
                    found = True
                    
            for i2 in range(len(loose_not_tight_muons)):
                if delta_r(events[i0].FatJet[i1],events[i0].Muon[loose_not_tight_muons[i2]]) < 0.5:
                    found = True

            for i2 in range(len(tight_electrons)):
                if delta_r(events[i0].FatJet[i1],events[i0].Electron[tight_electrons[i2]]) < 0.5:
                    found = True

            if not found:        
                cleaned_fatjets.append(i1)

        cleaned_jets = []                  

        for i1 in range(len(events[i0].Jet)):
            found = False

            for i2 in range(len(tight_muons)):
                if delta_r(events[i0].Jet[i1],events[i0].Muon[tight_muons[i2]]) < 0.5:
                    found = True
                    
            for i2 in range(len(loose_not_tight_muons)):
                if delta_r(events[i0].Jet[i1],events[i0].Muon[loose_not_tight_muons[i2]]) < 0.5:
                    found = True

            for i2 in range(len(tight_electrons)):
                if delta_r(events[i0].Jet[i1],events[i0].Electron[tight_electrons[i2]]) < 0.5:
                    found = True

            for i2 in range(len(cleaned_fatjets)):
                if delta_r(events[i0].Jet[i1],events[i0].FatJet[cleaned_fatjets[i2]]) < 0.5:
                    found = True

            if not found:        
                cleaned_jets.append(i1)

#        if len(cleaned_jets) < 4 or len(tight_muons) + len(loose_not_tight_muons) + len(tight_electrons) != 1:
        if len(cleaned_jets) < 2 or len(cleaned_fatjets) < 1 or len(tight_muons) + len(tight_electrons) != 1 or len(loose_not_tight_muons) != 0:
            builder.end_list()
            continue
            
        found = False   

        for i1 in range(len(events[i0].FatJet)):

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

                    if events[i0].FatJet[cleaned_fatjets[i1]].pt < 250 or abs(events[i0].FatJet[cleaned_fatjets[i1]].eta) < 2.5 or events[i0].FatJet[cleaned_fatjets[i1]].msoftdrop < 50 or events[i0].FatJet[cleaned_fatjets[i1]].msoftdrop > 150:
                        continue

                    if events[i0].Jet[cleaned_jets[i2]].btagDeepB > 0.2217 or events[i0].Jet[cleaned_jets[i3]].btagDeepB > 0.2217 or events[i0].Jet[cleaned_jets[i2]].pt < 30 or events[i0].Jet[cleaned_jets[i3]].pt < 30 or abs(events[i0].Jet[cleaned_jets[i2]].eta) > 4.7 or abs(events[i0].Jet[cleaned_jets[i3]].eta) > 4.7:
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

                    if events[i0].Jet[cleaned_jets[i4]].pt < 30 or abs(events[i0].Jet[cleaned_jets[i4]].eta) > 4.7:
                        continue

                    nextrajets=nextrajets+1

                    if events[i0].Jet[cleaned_jets[i4]].btagDeepB > 0.2217:
                        nextrabjets += 1

                builder.index(5).integer(nextrajets)
                builder.index(6).integer(nextrabjets)

                builder.end_tuple()    
                builder.end_list()                  
                        
    return builder

@numba.njit
def select_events_resolved (events,dataset,builder,syst="nominal"):
    
    for i0 in range(len(events)):
        
        builder.begin_list() 
        
        if dataset == "SingleMuon":
            if not events[i0].HLT.IsoTkMu24 and not events[i0].HLT.IsoMu24:
                builder.end_list()
                continue
        elif dataset == "SingleElectron":
            if events[i0].HLT.IsoTkMu24 or events[i0].HLT.IsoMu24 or not events[i0].HLT.Ele27_WPTight_Gsf :
                builder.end_list()
                continue
        else: #MC
            if not events[i0].HLT.IsoTkMu24 and not events[i0].HLT.IsoMu24 and not events[i0].HLT.Ele27_WPTight_Gsf:
                builder.end_list()
                continue
        
        if syst == "nominal":        
            if events[i0].PuppiMET.pt < 30:
                builder.end_list()
                continue
        elif syst == "JESUp":
            if events[i0].PuppiMET.ptJESUp < 30:
                builder.end_list()
                continue
        elif syst == "JERUp":
            if events[i0].PuppiMET.ptJERUp < 30:
                builder.end_list()
                continue
        else:
            builder.end_list()
            continue

        found = False

        tight_muons = []
        loose_not_tight_muons = []
            
        for i1 in range(len(events[i0].Muon)):
            if events[i0].Muon[i1].tightId==True and events[i0].Muon[i1].pfRelIso04_all < 0.15 and events[i0].Muon[i1].pt > 26 and abs(events[i0].Muon[i1].eta) < 2.4:
                tight_muons.append(i1)
            elif events[i0].Muon[i1].tightId==True and events[i0].Muon[i1].pfRelIso04_all < 0.4 and events[i0].Muon[i1].pt > 26 and abs(events[i0].Muon[i1].eta) < 2.4:   
                loose_not_tight_muons.append(i1)
                
        tight_electrons = [] 
            
        for i1 in range(len(events[i0].Electron)):
            if events[i0].Electron[i1].pt > 30 and abs(events[i0].Electron[i1].eta + events[i0].Electron[i1].deltaEtaSC) < 2.5:
                if (abs(events[i0].Electron[i1].eta + events[i0].Electron[i1].deltaEtaSC) < 1.479 and abs(events[i0].Electron[i1].dz) < 0.1 and abs(events[i0].Electron[i1].dxy) < 0.05) or (abs(events[i0].Electron[i1].eta + events[i0].Electron[i1].deltaEtaSC) > 1.479 and abs(events[i0].Electron[i1].dz) < 0.2 and abs(events[i0].Electron[i1].dxy) < 0.1):
                    if events[i0].Electron[i1].cutBased >= 3:
                        tight_electrons.append(i1)
                             
        cleaned_jets = []                  

        for i1 in range(len(events[i0].Jet)):
            found = False

            for i2 in range(len(tight_muons)):
                if delta_r(events[i0].Jet[i1],events[i0].Muon[tight_muons[i2]]) < 0.5:
                    found = True
                    
            for i2 in range(len(loose_not_tight_muons)):
                if delta_r(events[i0].Jet[i1],events[i0].Muon[loose_not_tight_muons[i2]]) < 0.5:
                    found = True

            for i2 in range(len(tight_electrons)):
                if delta_r(events[i0].Jet[i1],events[i0].Electron[tight_electrons[i2]]) < 0.5:
                    found = True

            if not found:        
                cleaned_jets.append(i1)

#        if len(cleaned_jets) < 4 or len(tight_muons) + len(loose_not_tight_muons) + len(tight_electrons) != 1:
        if len(cleaned_jets) < 4 or len(tight_muons) + len(tight_electrons) != 1 or len(loose_not_tight_muons) != 0:
            builder.end_list()
            continue
            
        found = False   
            
        for i1 in range(len(cleaned_jets)):

            if found:
                break

            for i2 in range(len(cleaned_jets)):

                if found:
                    break

                for i3 in range(len(cleaned_jets)):

                    if found:
                        break

                    for i4 in range(len(cleaned_jets)):
                        
                        if found:
                            break

                        if i1 == i2 or i1 == i3 or i1 == i4 or i2 == i3 or i2 == i4 or i3 == i4:
                            continue

                        if events[i0].Jet[cleaned_jets[i1]].btagDeepB < 0.8953 or events[i0].Jet[cleaned_jets[i2]].btagDeepB < 0.8953 or events[i0].Jet[cleaned_jets[i1]].pt < 30 or events[i0].Jet[cleaned_jets[i2]].pt < 30 or abs(events[i0].Jet[cleaned_jets[i1]].eta) > 2.5 or abs(events[i0].Jet[cleaned_jets[i2]].eta) > 2.5:
                            continue
                            
                        if events[i0].Jet[cleaned_jets[i3]].btagDeepB > 0.2217 or events[i0].Jet[cleaned_jets[i4]].btagDeepB > 0.2217 or events[i0].Jet[cleaned_jets[i3]].pt < 30 or events[i0].Jet[cleaned_jets[i4]].pt < 30 or abs(events[i0].Jet[cleaned_jets[i3]].eta) > 4.7 or abs(events[i0].Jet[cleaned_jets[i4]].eta) > 4.7:
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

                            if events[i0].Jet[cleaned_jets[i5]].pt < 30 or abs(events[i0].Jet[cleaned_jets[i2]].eta) > 4.7:
                                continue

                            nextrajets=nextrajets+1

                            if events[i0].Jet[cleaned_jets[i5]].btagDeepB > 0.2217:
                                nextrabjets += 1

                        builder.index(6).integer(nextrajets)
                        builder.index(7).integer(nextrabjets)

                        builder.end_tuple()    

        builder.end_list()                  
                        
    return builder                            
                        
class EwkwhjjProcessor(processor.ProcessorABC):
    def __init__(self):
        self._accumulator = processor.dict_accumulator({
            "sumw": processor.defaultdict_accumulator(float),
            "nevents": processor.defaultdict_accumulator(float),
            "sel1_bdtscore_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel1_bdtscore_binning1_pileupUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel1_bdtscore_binning1_prefireUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel1_bdtscore_binning1_electronidsfUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel1_bdtscore_binning1_muonidsfUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel1_bdtscore_binning1_muonisosfUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel1_bdtscore_binning1_muonhltsfUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel1_bdtscore_binning1_JESUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel1_bdtscore_binning1_JERUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel1_bdtscore_binning2": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 40, -0.5, 1.5),
            ),
            "sel1_bdtscore_binning3": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 21, 0, 1.05),
            ),
            "sel1_higgsdijetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsdijetmass", "Higgs dijet mass [GeV]", 75, 0, 300),
            ),
            "sel1_higgsdijetpt_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsdijetpt", "Higgs dijet pt [GeV]", 20, 0, 200),
            ),
            "sel1_vbsdijetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetmass", "VBS dijet mass [GeV]", 19, 100, 2000),
            ),
            "sel1_vbsdijetmass_binning2": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetmass", "VBS dijet mass [GeV]", 39, 100, 4000),
            ),
            "sel1_vbsdijetabsdeta_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetabsdeta", "VBS dijet $\Delta \eta$", 55, 2.5, 8),
            ),
            "sel1_leptonpt_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("leptonpt", "Lepton pt [GeV]", 19, 20, 200),
            ),
            "sel1_met_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("met", "MET [GeV]", 20, 0, 200),
            ),
            "sel2_bdtscore_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel2_bdtscore_binning1_pileupUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel2_bdtscore_binning1_prefireUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel2_bdtscore_binning1_electronidsfUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel2_bdtscore_binning1_JESUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel2_bdtscore_binning1_JERUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel2_bdtscore_binning2": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 40, -0.5, 1.5),
            ),
            "sel2_bdtscore_binning3": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 21, 0, 1.05),
            ),
            "sel2_higgsdijetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsdijetmass", "Higgs dijet mass [GeV]", 75, 0, 300),
            ),
            "sel2_higgsdijetpt_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsdijetpt", "Higgs dijet pt [GeV]", 20, 0, 200),
            ),
            "sel2_vbsdijetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetmass", "VBS dijet mass [GeV]", 19, 100, 2000),
            ),
            "sel2_vbsdijetmass_binning2": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetmass", "VBS dijet mass [GeV]", 39, 100, 4000),
            ),
            "sel2_vbsdijetabsdeta_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetabsdeta", "VBS dijet $\Delta \eta$", 55, 2.5, 8),
            ),
            "sel2_leptonpt_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("leptonpt", "Lepton pt [GeV]", 19, 20, 200),
            ),
            "sel2_met_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("met", "MET [GeV]", 20, 0, 200),
            ),
            "sel3_bdtscore_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel3_bdtscore_binning1_pileupUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel3_bdtscore_binning1_prefireUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel3_bdtscore_binning1_electronidsfUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel3_bdtscore_binning1_muonidsfUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel3_bdtscore_binning1_muonisosfUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel3_bdtscore_binning1_muonhltsfUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel3_bdtscore_binning1_JESUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel3_bdtscore_binning1_JERUp": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1),
            ),
            "sel3_bdtscore_binning2": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 40, -0.5, 1.5),
            ),
            "sel3_bdtscore_binning3": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 21, 0, 1.05),
            ),
            "sel3_higgsdijetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsdijetmass", "Higgs dijet mass [GeV]", 75, 0, 300),
            ),
            "sel3_higgsdijetpt_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsdijetpt", "Higgs dijet pt [GeV]", 20, 0, 200),
            ),
            "sel3_vbsdijetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetmass", "VBS dijet mass [GeV]", 19, 100, 2000),
            ),
            "sel3_vbsdijetmass_binning2": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetmass", "VBS dijet mass [GeV]", 39, 100, 4000),
            ),
            "sel3_vbsdijetabsdeta_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetabsdeta", "VBS dijet $\Delta \eta$", 55, 2.5, 8),
            ),
            "sel3_leptonpt_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("leptonpt", "Lepton pt [GeV]", 19, 20, 200),
            ),
            "sel3_met_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("met", "MET [GeV]", 20, 0, 200),
            ),
            "sel4_higgsdijetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsdijetmass", "Higgs dijet mass [GeV]", 75, 0, 300),
            ),
            "sel4_higgsdijetpt_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsdijetpt", "Higgs dijet pt [GeV]", 20, 0, 200),
            ),
            "sel4_vbsdijetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetmass", "VBS dijet mass [GeV]", 19, 100, 2000),
            ),
            "sel4_vbsdijetmass_binning2": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetmass", "VBS dijet mass [GeV]", 4, 100, 500),
            ),
            "sel4_leptonpt_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("leptonpt", "Lepton pt [GeV]", 19, 20, 200),
            ),
            "sel4_met_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("met", "MET [GeV]", 20, 0, 200),
            ),
            "sel5_higgsdijetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsdijetmass", "Higgs dijet mass [GeV]", 75, 0, 300),
            ),
            "sel5_higgsdijetpt_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsdijetpt", "Higgs dijet pt [GeV]", 20, 0, 200),
            ),
            "sel5_vbsdijetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetmass", "VBS dijet mass [GeV]", 19, 100, 2000),
            ),
            "sel5_vbsdijetmass_binning2": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetmass", "VBS dijet mass [GeV]", 4, 100, 500),
            ),
            "sel5_leptonpt_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("leptonpt", "Lepton pt [GeV]", 19, 20, 200),
            ),
            "sel5_met_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("met", "MET [GeV]", 20, 0, 200),
            ),
            "sel6_higgsdijetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsdijetmass", "Higgs dijet mass [GeV]", 75, 0, 300),
            ),
            "sel6_higgsdijetpt_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsdijetpt", "Higgs dijet pt [GeV]", 20, 0, 200),
            ),
            "sel6_vbsdijetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetmass", "VBS dijet mass [GeV]", 19, 100, 2000),
            ),
            "sel6_vbsdijetmass_binning2": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("vbsdijetmass", "VBS dijet mass [GeV]", 4, 100, 500),
            ),
            "sel6_leptonpt_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("leptonpt", "Lepton pt [GeV]", 19, 20, 200),
            ),
            "sel6_met_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("met", "MET [GeV]", 20, 0, 200),
            ),
            "sel7_higgsjetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsjetmass", "Higgs jet mass [GeV]", 75, 0, 300),
            ),
            "sel7_higgsjetsoftdropmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsjetsoftdropmass", "Higgs jet soft drop mass [GeV]", 75, 0, 300),
            ),
            "sel7_bdtscore_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1.0),
            ),
            "sel7_bdtscore_binning2": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 40, -0.5, 1.5),
            ),
            "sel7_bdtscore_binning3": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 21, 0, 1.05),
            ),
            "sel8_higgsjetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsjetmass", "Higgs jet mass [GeV]", 75, 0, 300),
            ),
            "sel8_higgsjetsoftdropmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsjetsoftdropmass", "Higgs jet soft drop mass [GeV]", 75, 0, 300),
            ),
            "sel8_bdtscore_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1.0),
            ),
            "sel8_bdtscore_binning2": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 40, -0.5, 1.5),
            ),
            "sel8_bdtscore_binning3": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 21, 0, 1.05),
            ),
            "sel9_higgsjetmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsjetmass", "Higgs jet mass [GeV]", 75, 0, 300),
            ),
            "sel9_higgsjetsoftdropmass_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("higgsjetsoftdropmass", "Higgs jet soft drop mass [GeV]", 75, 0, 300),
            ),
            "sel9_bdtscore_binning1": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 20, 0, 1.0),
            ),
            "sel9_bdtscore_binning2": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 40, -0.5, 1.5),
            ),
            "sel9_bdtscore_binning3": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("bdtscore", "BDT score", 21, 0, 1.05),
            ),

        })

    @property
    def accumulator(self):
        return self._accumulator

#    @numba.njit
    def process(self, events):

        output = self.accumulator.identity()

        dataset = events.metadata['dataset']

        if "Single" not in dataset:
            output["sumw"][dataset] += ak.sum(np.sign(events.Generator.weight))
        output["nevents"][dataset] += len(events)

        if "Single" in dataset :
            events = events[lumimask(events.run,events.luminosityBlock)]
   
        particleindices = select_events_resolved(events,dataset,ak.ArrayBuilder()).snapshot()

        if "Single" not in dataset:
            particleindices_JESUp = select_events_resolved(events,dataset,ak.ArrayBuilder(),syst="JESUp").snapshot()
            particleindices_JERUp = select_events_resolved(events,dataset,ak.ArrayBuilder(),syst="JERUp").snapshot()

        particleindices_merged = select_events_merged(events,dataset,ak.ArrayBuilder()).snapshot()
        
        basecut = ak.num(particleindices) != 0

        if "Single" not in dataset:
            basecut_JESUp = ak.num(particleindices_JESUp) != 0
            basecut_JERUp = ak.num(particleindices_JERUp) != 0

        basecut_merged = ak.num(particleindices_merged) != 0

        if "Single" in dataset:
            dataset = "Data"

        if ak.any(basecut_merged):
            particleindices_merged = particleindices_merged[basecut_merged]
            events_merged = events[basecut_merged]
            events_merged.FatJet[particleindices_merged["0"]]
            jets_merged = [events_merged.Jet[particleindices_merged[idx]] for idx in "12"]
            cut7 = ak.firsts((events_merged.FatJet[particleindices_merged["0"]].mass > 50) & (events_merged.FatJet[particleindices_merged["0"]].mass < 150) & ((jets_merged[0]+jets_merged[1]).mass > 500) & (abs(jets_merged[0].eta - jets_merged[1].eta) > 2.5) & (particleindices_merged["3"] != -1))
            cut8 = ak.firsts((events_merged.FatJet[particleindices_merged["0"]].mass > 50) & (events_merged.FatJet[particleindices_merged["0"]].mass < 150) & ((jets_merged[0]+jets_merged[1]).mass > 500) & (abs(jets_merged[0].eta - jets_merged[1].eta) > 2.5) & (particleindices_merged["4"] != -1))
#            cut9_merged = cut7_merged | cut8_merged

        if dataset != "Data" and ak.any(basecut_JESUp):
            particleindices_JESUp = particleindices_JESUp[basecut_JESUp]
            events_JESUp = events[basecut_JESUp]
            jets_JESUp = [events_JESUp.Jet[particleindices_JESUp[idx]] for idx in "0123"]
            cut1_JESUp = ak.firsts(((jets_JESUp[0]+jets_JESUp[1]).mass > 50) & ((jets_JESUp[0]+jets_JESUp[1]).mass < 150) & ((jets_JESUp[2]+jets_JESUp[3]).mass > 500) & (abs(jets_JESUp[2].eta - jets_JESUp[3].eta) > 2.5) & (particleindices_JESUp["4"] != -1))
            cut2_JESUp = ak.firsts(((jets_JESUp[0]+jets_JESUp[1]).mass > 50) & ((jets_JESUp[0]+jets_JESUp[1]).mass < 150) & ((jets_JESUp[2]+jets_JESUp[3]).mass > 500) & (abs(jets_JESUp[2].eta - jets_JESUp[3].eta) > 2.5) & (particleindices_JESUp["5"] != -1))
#            cut3_JESUp = cut1_JESUp | cut2_JESUp

        if dataset != "Data" and ak.any(basecut_JERUp):
            particleindices_JERUp = particleindices_JERUp[basecut_JERUp]
            events_JERUp = events[basecut_JERUp]
            jets_JERUp= [events_JERUp.Jet[particleindices_JERUp[idx]] for idx in "0123"]        
            cut1_JERUp = ak.firsts(((jets_JERUp[0]+jets_JERUp[1]).mass > 50) & ((jets_JERUp[0]+jets_JERUp[1]).mass < 150) & ((jets_JERUp[2]+jets_JERUp[3]).mass > 500) & (abs(jets_JERUp[2].eta - jets_JERUp[3].eta) > 2.5) & (particleindices_JERUp["4"] != -1))
            cut2_JERUp = ak.firsts(((jets_JERUp[0]+jets_JERUp[1]).mass > 50) & ((jets_JERUp[0]+jets_JERUp[1]).mass < 150) & ((jets_JERUp[2]+jets_JERUp[3]).mass > 500) & (abs(jets_JERUp[2].eta - jets_JERUp[3].eta) > 2.5) & (particleindices_JERUp["5"] != -1))
#            cut3_JERUp = cut1_JERUp | cut2_JERUp

        if ak.any(basecut):
            particleindices = particleindices[basecut]
            events = events[basecut]
            jets = [events.Jet[particleindices[idx]] for idx in "0123"]
            cut1 = ak.firsts(((jets[0]+jets[1]).mass > 50) & ((jets[0]+jets[1]).mass < 150) & ((jets[2]+jets[3]).mass > 500) & (abs(jets[2].eta - jets[3].eta) > 2.5) & (particleindices["4"] != -1))
            cut2 = ak.firsts(((jets[0]+jets[1]).mass > 50) & ((jets[0]+jets[1]).mass < 150) & ((jets[2]+jets[3]).mass > 500) & (abs(jets[2].eta - jets[3].eta) > 2.5) & (particleindices["5"] != -1))
#            cut3 = cut1 | cut2
            cut4 = ak.firsts(((jets[0]+jets[1]).mass > 50) & ((jets[0]+jets[1]).mass < 150) & ((jets[2]+jets[3]).mass > 100) & ((jets[2]+jets[3]).mass < 500) & (particleindices["4"] != -1))
            cut5 = ak.firsts(((jets[0]+jets[1]).mass > 50) & ((jets[0]+jets[1]).mass < 150) & ((jets[2]+jets[3]).mass > 100) & ((jets[2]+jets[3]).mass < 500) & (particleindices["5"] != -1))
#            cut6 = cut4 | cut5

        if ak.any(basecut_merged) and ak.any(cut7):

            sel7_particleindices = particleindices_merged[cut7]
        
            sel7_events = events_merged[cut7]
            
            sel7_fatjets = sel7_events.FatJet[sel7_particleindices["0"]]

            sel7_jets = [sel7_events.Jet[sel7_particleindices[idx]] for idx in "12"]

            sel7_muons = sel7_events.Muon[sel7_particleindices["3"]]

            sel7_muonidsf = ak.firsts(evaluator["muonidsf"](abs(sel7_muons.eta), sel7_muons.pt))
            sel7_muonisosf = ak.firsts(evaluator["muonisosf"](abs(sel7_muons.eta), sel7_muons.pt))
            sel7_muonhltsf = ak.firsts(evaluator["muonhltsf"](abs(sel7_muons.eta), sel7_muons.pt))

            sel7_X = pandas.DataFrame(np.transpose(np.vstack((
                ak.to_numpy(ak.firsts(sel7_fatjets).pt).data,
                ak.to_numpy(ak.firsts(sel7_fatjets).eta).data,
                ak.to_numpy(ak.firsts(sel7_fatjets).phi).data,
                ak.to_numpy(ak.firsts(sel7_fatjets).btagDeepB).data,
                ak.to_numpy(ak.firsts(sel7_fatjets).btagHbb).data,
                ak.to_numpy(ak.firsts(sel7_fatjets).msoftdrop).data,
                ak.to_numpy(ak.firsts(sel7_particleindices["5"])).data,
                ak.to_numpy(ak.firsts(sel7_particleindices["6"])).data,
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

            sel7_d = xgb.DMatrix(sel7_X)

            sel7_bdtscore = bst_merged.predict(sel7_d)


            if dataset == "Data":
                sel7_weight = np.ones(len(sel7_events))
            else:    
                sel7_weight = np.sign(sel7_events.Generator.weight)*sel7_events.L1PreFiringWeight.Nom*sel7_muonidsf*sel7_muonisosf*sel7_muonhltsf

            output["sel7_higgsjetmass_binning1"].fill(
                dataset=dataset,
                higgsjetmass=ak.firsts(sel7_events.FatJet[sel7_particleindices["0"]].mass),
                weight=sel7_weight
            )

            output["sel7_higgsjetsoftdropmass_binning1"].fill(
                dataset=dataset,
                higgsjetsoftdropmass=ak.firsts(sel7_events.FatJet[sel7_particleindices["0"]].msoftdrop),
                weight=sel7_weight
            )

            output["sel7_bdtscore_binning1"].fill(
                dataset=dataset,
                bdtscore=sel7_bdtscore,
                weight=sel7_weight
            )

            output["sel7_bdtscore_binning2"].fill(
                dataset=dataset,
                bdtscore=sel7_bdtscore,
                weight=sel7_weight
            )

            output["sel7_bdtscore_binning3"].fill(
                dataset=dataset,
                bdtscore=sel7_bdtscore,
                weight=sel7_weight
            )

            output["sel9_higgsjetmass_binning1"].fill(
                dataset=dataset,
                higgsjetmass=ak.firsts(sel7_events.FatJet[sel7_particleindices["0"]].mass),
                weight=sel7_weight
            )

            output["sel9_higgsjetsoftdropmass_binning1"].fill(
                dataset=dataset,
                higgsjetsoftdropmass=ak.firsts(sel7_events.FatJet[sel7_particleindices["0"]].msoftdrop),
                weight=sel7_weight
            )

            output["sel9_bdtscore_binning1"].fill(
                dataset=dataset,
                bdtscore=sel7_bdtscore,
                weight=sel7_weight
            )

            output["sel9_bdtscore_binning2"].fill(
                dataset=dataset,
                bdtscore=sel7_bdtscore,
                weight=sel7_weight
            )

            output["sel9_bdtscore_binning3"].fill(
                dataset=dataset,
                bdtscore=sel7_bdtscore,
                weight=sel7_weight
            )

        if ak.any(basecut_merged) and ak.any(cut8):

            sel8_particleindices = particleindices_merged[cut8]
        
            sel8_events = events_merged[cut8]
            
            sel8_fatjets = sel8_events.FatJet[sel8_particleindices["0"]]

            sel8_jets = [sel8_events.Jet[sel8_particleindices[idx]] for idx in "12"]

            sel8_electrons = sel8_events.Electron[sel8_particleindices["4"]]

            sel8_electronidsf = ak.firsts(evaluator["electronidsf"](sel8_electrons.eta, sel8_electrons.pt))

            sel8_X = pandas.DataFrame(np.transpose(np.vstack((
                ak.to_numpy(ak.firsts(sel8_fatjets).pt).data,
                ak.to_numpy(ak.firsts(sel8_fatjets).eta).data,
                ak.to_numpy(ak.firsts(sel8_fatjets).phi).data,
                ak.to_numpy(ak.firsts(sel8_fatjets).btagDeepB).data,
                ak.to_numpy(ak.firsts(sel8_fatjets).btagHbb).data,
                ak.to_numpy(ak.firsts(sel8_fatjets).msoftdrop).data,
                ak.to_numpy(ak.firsts(sel8_particleindices["5"])).data,
                ak.to_numpy(ak.firsts(sel8_particleindices["6"])).data,
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

            sel8_d = xgb.DMatrix(sel8_X)

            sel8_bdtscore = bst_merged.predict(sel8_d)

            if dataset == "Data":
                sel8_weight = np.ones(len(sel8_events))
            else:    
                sel8_weight = np.sign(sel8_events.Generator.weight)*sel8_events.L1PreFiringWeight.Nom*sel8_electronidsf

            output["sel8_higgsjetmass_binning1"].fill(
                dataset=dataset,
                higgsjetmass=ak.firsts(sel8_events.FatJet[sel8_particleindices["0"]].mass),
                weight=sel8_weight
            )

            output["sel8_higgsjetsoftdropmass_binning1"].fill(
                dataset=dataset,
                higgsjetsoftdropmass=ak.firsts(sel8_events.FatJet[sel8_particleindices["0"]].msoftdrop),
                weight=sel8_weight
            )

            output["sel8_bdtscore_binning1"].fill(
                dataset=dataset,
                bdtscore=sel8_bdtscore,
                weight=sel8_weight
            )

            output["sel8_bdtscore_binning2"].fill(
                dataset=dataset,
                bdtscore=sel8_bdtscore,
                weight=sel8_weight
            )

            output["sel8_bdtscore_binning3"].fill(
                dataset=dataset,
                bdtscore=sel8_bdtscore,
                weight=sel8_weight
            )

            output["sel9_higgsjetmass_binning1"].fill(
                dataset=dataset,
                higgsjetmass=ak.firsts(sel8_events.FatJet[sel8_particleindices["0"]].mass),
                weight=sel8_weight
            )

            output["sel9_higgsjetsoftdropmass_binning1"].fill(
                dataset=dataset,
                higgsjetsoftdropmass=ak.firsts(sel8_events.FatJet[sel8_particleindices["0"]].msoftdrop),
                weight=sel8_weight
            )

            output["sel9_bdtscore_binning1"].fill(
                dataset=dataset,
                bdtscore=sel8_bdtscore,
                weight=sel8_weight
            )

            output["sel9_bdtscore_binning2"].fill(
                dataset=dataset,
                bdtscore=sel8_bdtscore,
                weight=sel8_weight
            )

            output["sel9_bdtscore_binning3"].fill(
                dataset=dataset,
                bdtscore=sel8_bdtscore,
                weight=sel8_weight
            )

        if dataset != "Data" and ak.any(basecut_JESUp) and ak.any(cut1_JESUp):

            sel1_JESUp_particleindices = particleindices_JESUp[cut1_JESUp]
        
            sel1_JESUp_events = events_JESUp[cut1_JESUp]
            
            sel1_JESUp_jets = [sel1_JESUp_events.Jet[sel1_JESUp_particleindices[idx]] for idx in "0123"]

            sel1_JESUp_muons = sel1_JESUp_events.Muon[sel1_JESUp_particleindices["4"]]

            sel1_JESUp_muonidsf = ak.firsts(evaluator["muonidsf"](abs(sel1_JESUp_muons.eta), sel1_JESUp_muons.pt))
            sel1_JESUp_muonisosf = ak.firsts(evaluator["muonisosf"](abs(sel1_JESUp_muons.eta), sel1_JESUp_muons.pt))
            sel1_JESUp_muonhltsf = ak.firsts(evaluator["muonhltsf"](abs(sel1_JESUp_muons.eta), sel1_JESUp_muons.pt))

            sel1_JESUp_X = pandas.DataFrame(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel1_JESUp_particleindices["6"])).data,ak.to_numpy(ak.firsts(sel1_JESUp_particleindices["7"])).data,np.zeros(len(sel1_JESUp_events)),np.sign(ak.to_numpy(ak.firsts(sel1_JESUp_muons).charge).data+1),ak.to_numpy(ak.firsts(sel1_JESUp_muons).pt).data,ak.to_numpy(ak.firsts(sel1_JESUp_muons).eta).data,ak.to_numpy(ak.firsts(sel1_JESUp_muons).phi).data,ak.to_numpy(sel1_JESUp_events.PuppiMET.ptJESUp),ak.to_numpy(sel1_JESUp_events.PuppiMET.phiJESUp),ak.to_numpy(ak.firsts(sel1_JESUp_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel1_JESUp_jets[0]+sel1_JESUp_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel1_JESUp_jets[2]+sel1_JESUp_jets[3]).mass)).data, ak.to_numpy(ak.firsts(sel1_JESUp_jets[2]).eta - ak.firsts(sel1_JESUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_JESUp_muons+sel1_JESUp_jets[0]).pt*sel1_JESUp_events.PuppiMET.ptJESUp*(1 - np.cos(sel1_JESUp_events.PuppiMET.phiJESUp - (sel1_JESUp_muons+sel1_JESUp_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_JESUp_muons+sel1_JESUp_jets[1]).pt*sel1_JESUp_events.PuppiMET.ptJESUp*(1 - np.cos(sel1_JESUp_events.PuppiMET.phiJESUp - (sel1_JESUp_muons+sel1_JESUp_jets[1]).phi))))).data))),columns=["nextrajets","nextrabjets","leptonflavor","leptoncharge","leptonpt","leptoneta","leptonphi",'metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel1_JESUp_d = xgb.DMatrix(sel1_JESUp_X)

            sel1_JESUp_bdtscore = bst.predict(sel1_JESUp_d)

            sel1_JESUp_weight = np.sign(sel1_JESUp_events.Generator.weight)*sel1_JESUp_events.L1PreFiringWeight.Nom*sel1_JESUp_muonidsf*sel1_JESUp_muonisosf*sel1_JESUp_muonhltsf

            output["sel1_bdtscore_binning1_JESUp"].fill(
                dataset=dataset,
                bdtscore=sel1_JESUp_bdtscore,
                weight=sel1_JESUp_weight
            )

            output["sel3_bdtscore_binning1_JESUp"].fill(
                dataset=dataset,
                bdtscore=sel1_JESUp_bdtscore,
                weight=sel1_JESUp_weight
            )

        if dataset != "Data" and ak.any(basecut_JERUp) and ak.any(cut1_JERUp):

            sel1_JERUp_particleindices = particleindices_JERUp[cut1_JERUp]
        
            sel1_JERUp_events = events_JERUp[cut1_JERUp]
            
            sel1_JERUp_jets = [sel1_JERUp_events.Jet[sel1_JERUp_particleindices[idx]] for idx in "0123"]

            sel1_JERUp_muons = sel1_JERUp_events.Muon[sel1_JERUp_particleindices["4"]]

            sel1_JERUp_muonidsf = ak.firsts(evaluator["muonidsf"](abs(sel1_JERUp_muons.eta), sel1_JERUp_muons.pt))
            sel1_JERUp_muonisosf = ak.firsts(evaluator["muonisosf"](abs(sel1_JERUp_muons.eta), sel1_JERUp_muons.pt))
            sel1_JERUp_muonhltsf = ak.firsts(evaluator["muonhltsf"](abs(sel1_JERUp_muons.eta), sel1_JERUp_muons.pt))

            sel1_JERUp_X = pandas.DataFrame(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel1_JERUp_particleindices["6"])).data,ak.to_numpy(ak.firsts(sel1_JERUp_particleindices["7"])).data,np.zeros(len(sel1_JERUp_events)),np.sign(ak.to_numpy(ak.firsts(sel1_JERUp_muons).charge).data+1),ak.to_numpy(ak.firsts(sel1_JERUp_muons).pt).data,ak.to_numpy(ak.firsts(sel1_JERUp_muons).eta).data,ak.to_numpy(ak.firsts(sel1_JERUp_muons).phi).data,ak.to_numpy(sel1_JERUp_events.PuppiMET.ptJERUp),ak.to_numpy(sel1_JERUp_events.PuppiMET.phiJERUp),ak.to_numpy(ak.firsts(sel1_JERUp_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel1_JERUp_jets[0]+sel1_JERUp_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel1_JERUp_jets[2]+sel1_JERUp_jets[3]).mass)).data, ak.to_numpy(ak.firsts(sel1_JERUp_jets[2]).eta - ak.firsts(sel1_JERUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_JERUp_muons+sel1_JERUp_jets[0]).pt*sel1_JERUp_events.PuppiMET.ptJERUp*(1 - np.cos(sel1_JERUp_events.PuppiMET.phiJERUp - (sel1_JERUp_muons+sel1_JERUp_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_JERUp_muons+sel1_JERUp_jets[1]).pt*sel1_JERUp_events.PuppiMET.ptJERUp*(1 - np.cos(sel1_JERUp_events.PuppiMET.phiJERUp - (sel1_JERUp_muons+sel1_JERUp_jets[1]).phi))))).data))),columns=["nextrajets","nextrabjets","leptonflavor","leptoncharge","leptonpt","leptoneta","leptonphi",'metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel1_JERUp_d = xgb.DMatrix(sel1_JERUp_X)

            sel1_JERUp_bdtscore = bst.predict(sel1_JERUp_d)

            sel1_JERUp_weight = np.sign(sel1_JERUp_events.Generator.weight)*sel1_JERUp_events.L1PreFiringWeight.Nom*sel1_JERUp_muonidsf*sel1_JERUp_muonisosf*sel1_JERUp_muonhltsf

            output["sel1_bdtscore_binning1_JERUp"].fill(
                dataset=dataset,
                bdtscore=sel1_JERUp_bdtscore,
                weight=sel1_JERUp_weight
            )

            output["sel3_bdtscore_binning1_JERUp"].fill(
                dataset=dataset,
                bdtscore=sel1_JERUp_bdtscore,
                weight=sel1_JERUp_weight
            )

        if ak.any(basecut) and ak.any(cut1):

            sel1_particleindices = particleindices[cut1]
        
            sel1_events = events[cut1]
            
            sel1_jets = [sel1_events.Jet[sel1_particleindices[idx]] for idx in "0123"]

            sel1_muons = sel1_events.Muon[sel1_particleindices["4"]]

            sel1_muonidsf = ak.firsts(evaluator["muonidsf"](abs(sel1_muons.eta), sel1_muons.pt))
            sel1_muonisosf = ak.firsts(evaluator["muonisosf"](abs(sel1_muons.eta), sel1_muons.pt))
            sel1_muonhltsf = ak.firsts(evaluator["muonhltsf"](abs(sel1_muons.eta), sel1_muons.pt))
            sel1_muonidsfUp = ak.firsts(evaluator["muonidsfunc"](abs(sel1_muons.eta), sel1_muons.pt))+sel1_muonidsf
            sel1_muonisosfUp = ak.firsts(evaluator["muonisosfunc"](abs(sel1_muons.eta), sel1_muons.pt))+sel1_muonisosf
            sel1_muonhltsfUp = ak.firsts(evaluator["muonhltsfunc"](abs(sel1_muons.eta), sel1_muons.pt))+sel1_muonhltsf

            sel1_X = pandas.DataFrame(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel1_particleindices["6"])).data,ak.to_numpy(ak.firsts(sel1_particleindices["7"])).data,np.zeros(len(sel1_events)),np.sign(ak.to_numpy(ak.firsts(sel1_muons).charge).data+1),ak.to_numpy(ak.firsts(sel1_muons).pt).data,ak.to_numpy(ak.firsts(sel1_muons).eta).data,ak.to_numpy(ak.firsts(sel1_muons).phi).data,ak.to_numpy(sel1_events.PuppiMET.pt),ak.to_numpy(sel1_events.PuppiMET.phi),ak.to_numpy(ak.firsts(sel1_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel1_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel1_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel1_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel1_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel1_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel1_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel1_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel1_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel1_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel1_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel1_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel1_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel1_jets[0]+sel1_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel1_jets[2]+sel1_jets[3]).mass)).data, ak.to_numpy(ak.firsts(sel1_jets[2]).eta - ak.firsts(sel1_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_muons+sel1_jets[0]).pt*sel1_events.PuppiMET.pt*(1 - np.cos(sel1_events.PuppiMET.phi - (sel1_muons+sel1_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_muons+sel1_jets[1]).pt*sel1_events.PuppiMET.pt*(1 - np.cos(sel1_events.PuppiMET.phi - (sel1_muons+sel1_jets[1]).phi))))).data))),columns=["nextrajets","nextrabjets","leptonflavor","leptoncharge","leptonpt","leptoneta","leptonphi",'metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel1_d = xgb.DMatrix(sel1_X)

            sel1_bdtscore = bst.predict(sel1_d)    

            if dataset == "Data":
                sel1_weight = np.ones(len(sel1_events))
            else:
                sel1_weight = np.sign(sel1_events.Generator.weight)*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf
                sel1_weight_prefireUp = np.sign(sel1_events.Generator.weight)*sel1_events.L1PreFiringWeight.Up*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf
                sel1_weight_muonidsfUp = np.sign(sel1_events.Generator.weight)*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsfUp*sel1_muonisosf*sel1_muonhltsf
                sel1_weight_muonisosfUp = np.sign(sel1_events.Generator.weight)*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosfUp*sel1_muonhltsf
                sel1_weight_muonhltsfUp = np.sign(sel1_events.Generator.weight)*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsfUp

                
            output["sel1_bdtscore_binning1"].fill(
                dataset=dataset,
                bdtscore=sel1_bdtscore,
                weight=sel1_weight
            )

            if dataset != "Data":
                output["sel1_bdtscore_binning1_prefireUp"].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_prefireUp
                )

                output["sel1_bdtscore_binning1_muonidsfUp"].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonidsfUp
                )
                
                output["sel1_bdtscore_binning1_muonisosfUp"].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonisosfUp
                )

                output["sel1_bdtscore_binning1_muonhltsfUp"].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonhltsfUp
                )

            output["sel1_bdtscore_binning2"].fill(
                dataset=dataset,
                bdtscore=sel1_bdtscore,
                weight=sel1_weight
            )

            output["sel1_bdtscore_binning3"].fill(
                dataset=dataset,
                bdtscore=sel1_bdtscore,
                weight=sel1_weight
            )

            output["sel1_higgsdijetmass_binning1"].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel1_jets[0]+sel1_jets[1]).mass),
                weight=sel1_weight
            )

            output["sel1_higgsdijetpt_binning1"].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel1_jets[0]+sel1_jets[1]).pt),
                weight=sel1_weight
            )
        
            output["sel1_vbsdijetmass_binning1"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel1_jets[2]+sel1_jets[3]).mass),
                weight=sel1_weight
            )

            output["sel1_vbsdijetmass_binning2"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel1_jets[2]+sel1_jets[3]).mass),
                weight=sel1_weight
            )

            output["sel1_vbsdijetabsdeta_binning1"].fill(
                dataset=dataset,
                vbsdijetabsdeta=ak.firsts(sel1_jets[2]).eta - ak.firsts(sel1_jets[3]).eta,
                weight=sel1_weight
            )

            output["sel1_leptonpt_binning1"].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel1_muons.pt),
                weight=sel1_weight
            )

            output["sel1_met_binning1"].fill(
                dataset=dataset,
                met=sel1_events.PuppiMET.pt,
                weight=sel1_weight
            )

            output["sel3_bdtscore_binning1"].fill(
                dataset=dataset,
                bdtscore = sel1_bdtscore,
                weight=sel1_weight
            )

            if dataset != "Data":

                output["sel3_bdtscore_binning1_prefireUp"].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_prefireUp
                )
                
                output["sel3_bdtscore_binning1_electronidsfUp"].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight
                )
                
                output["sel3_bdtscore_binning1_muonidsfUp"].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonidsfUp
                )
                
                output["sel3_bdtscore_binning1_muonisosfUp"].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonisosfUp
                )
                
                output["sel3_bdtscore_binning1_muonhltsfUp"].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonhltsfUp
                )

            output["sel3_bdtscore_binning2"].fill(
                dataset=dataset,
                bdtscore = sel1_bdtscore,
                weight=sel1_weight
            )

            output["sel3_bdtscore_binning3"].fill(
                dataset=dataset,
                bdtscore=sel1_bdtscore,
                weight=sel1_weight
            )

            output["sel3_higgsdijetmass_binning1"].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel1_jets[0]+sel1_jets[1]).mass),
                weight=sel1_weight
            )

            output["sel3_higgsdijetpt_binning1"].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel1_jets[0]+sel1_jets[1]).pt),
                weight=sel1_weight
            )
        
            output["sel3_vbsdijetmass_binning1"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel1_jets[2]+sel1_jets[3]).mass),
                weight=sel1_weight
            )

            output["sel3_vbsdijetmass_binning2"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel1_jets[2]+sel1_jets[3]).mass),
                weight=sel1_weight
            )

            output["sel3_vbsdijetabsdeta_binning1"].fill(
                dataset=dataset,
                vbsdijetabsdeta=ak.firsts(sel1_jets[2]).eta - ak.firsts(sel1_jets[3]).eta,
                weight=sel1_weight
            )

            output["sel3_leptonpt_binning1"].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel1_muons.pt),
                weight=sel1_weight
            )

            output["sel3_met_binning1"].fill(
                dataset=dataset,
                met=sel1_events.PuppiMET.pt,
                weight=sel1_weight
            )

        if dataset != "Data" and ak.any(basecut_JESUp) and ak.any(cut2_JESUp):

            sel2_JESUp_particleindices = particleindices_JESUp[cut2_JESUp]
        
            sel2_JESUp_events = events_JESUp[cut2_JESUp]
            
            sel2_JESUp_jets = [sel2_JESUp_events.Jet[sel2_JESUp_particleindices[idx]] for idx in "0123"]

            sel2_JESUp_electrons = sel2_JESUp_events.Electron[sel2_JESUp_particleindices["4"]]

            sel2_JESUp_electronidsf = ak.firsts(evaluator["electronidsf"](sel2_JESUp_electrons.eta, sel2_JESUp_electrons.pt))

            sel2_JESUp_X = pandas.DataFrame(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel2_JESUp_particleindices["6"])).data,ak.to_numpy(ak.firsts(sel2_JESUp_particleindices["7"])).data,np.zeros(len(sel2_JESUp_events)),np.sign(ak.to_numpy(ak.firsts(sel2_JESUp_electrons).charge).data+1),ak.to_numpy(ak.firsts(sel2_JESUp_electrons).pt).data,ak.to_numpy(ak.firsts(sel2_JESUp_electrons).eta).data,ak.to_numpy(ak.firsts(sel2_JESUp_electrons).phi).data,ak.to_numpy(sel2_JESUp_events.PuppiMET.ptJESUp),ak.to_numpy(sel2_JESUp_events.PuppiMET.phiJESUp),ak.to_numpy(ak.firsts(sel2_JESUp_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel2_JESUp_jets[0]+sel2_JESUp_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel2_JESUp_jets[2]+sel2_JESUp_jets[3]).mass)).data, ak.to_numpy(ak.firsts(sel2_JESUp_jets[2]).eta - ak.firsts(sel2_JESUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_JESUp_electrons+sel2_JESUp_jets[0]).pt*sel2_JESUp_events.PuppiMET.ptJESUp*(1 - np.cos(sel2_JESUp_events.PuppiMET.phiJESUp - (sel2_JESUp_electrons+sel2_JESUp_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_JESUp_electrons+sel2_JESUp_jets[1]).pt*sel2_JESUp_events.PuppiMET.ptJESUp*(1 - np.cos(sel2_JESUp_events.PuppiMET.phiJESUp - (sel2_JESUp_electrons+sel2_JESUp_jets[1]).phi))))).data))),columns=["nextrajets","nextrabjets","leptonflavor","leptoncharge","leptonpt","leptoneta","leptonphi",'metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel2_JESUp_d = xgb.DMatrix(sel2_JESUp_X)

            sel2_JESUp_bdtscore = bst.predict(sel2_JESUp_d)

            sel2_JESUp_weight = np.sign(sel2_JESUp_events.Generator.weight)*sel2_JESUp_events.L1PreFiringWeight.Nom*sel2_JESUp_electronidsf

            output["sel2_bdtscore_binning1_JESUp"].fill(
                dataset=dataset,
                bdtscore=sel2_JESUp_bdtscore,
                weight=sel2_JESUp_weight
            )

            output["sel3_bdtscore_binning1_JESUp"].fill(
                dataset=dataset,
                bdtscore=sel2_JESUp_bdtscore,
                weight=sel2_JESUp_weight
            )

        if dataset != "Data" and ak.any(basecut_JERUp) and ak.any(cut2_JERUp):

            sel2_JERUp_particleindices = particleindices_JERUp[cut2_JERUp]
        
            sel2_JERUp_events = events_JERUp[cut2_JERUp]
            
            sel2_JERUp_jets = [sel2_JERUp_events.Jet[sel2_JERUp_particleindices[idx]] for idx in "0123"]

            sel2_JERUp_electrons = sel2_JERUp_events.Electron[sel2_JERUp_particleindices["4"]]

            sel2_JERUp_electronidsf = ak.firsts(evaluator["electronidsf"](sel2_JERUp_electrons.eta, sel2_JERUp_electrons.pt))

            sel2_JERUp_X = pandas.DataFrame(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel2_JERUp_particleindices["6"])).data,ak.to_numpy(ak.firsts(sel2_JERUp_particleindices["7"])).data,np.zeros(len(sel2_JERUp_events)),np.sign(ak.to_numpy(ak.firsts(sel2_JERUp_electrons).charge).data+1),ak.to_numpy(ak.firsts(sel2_JERUp_electrons).pt).data,ak.to_numpy(ak.firsts(sel2_JERUp_electrons).eta).data,ak.to_numpy(ak.firsts(sel2_JERUp_electrons).phi).data,ak.to_numpy(sel2_JERUp_events.PuppiMET.ptJERUp),ak.to_numpy(sel2_JERUp_events.PuppiMET.phiJERUp),ak.to_numpy(ak.firsts(sel2_JERUp_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel2_JERUp_jets[0]+sel2_JERUp_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel2_JERUp_jets[2]+sel2_JERUp_jets[3]).mass)).data, ak.to_numpy(ak.firsts(sel2_JERUp_jets[2]).eta - ak.firsts(sel2_JERUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_JERUp_electrons+sel2_JERUp_jets[0]).pt*sel2_JERUp_events.PuppiMET.ptJERUp*(1 - np.cos(sel2_JERUp_events.PuppiMET.phiJERUp - (sel2_JERUp_electrons+sel2_JERUp_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_JERUp_electrons+sel2_JERUp_jets[1]).pt*sel2_JERUp_events.PuppiMET.ptJERUp*(1 - np.cos(sel2_JERUp_events.PuppiMET.phiJERUp - (sel2_JERUp_electrons+sel2_JERUp_jets[1]).phi))))).data))),columns=["nextrajets","nextrabjets","leptonflavor","leptoncharge","leptonpt","leptoneta","leptonphi",'metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel2_JERUp_d = xgb.DMatrix(sel2_JERUp_X)

            sel2_JERUp_bdtscore = bst.predict(sel2_JERUp_d)

            sel2_JERUp_weight = np.sign(sel2_JERUp_events.Generator.weight)*sel2_JERUp_events.L1PreFiringWeight.Nom*sel2_JERUp_electronidsf

            output["sel2_bdtscore_binning1_JERUp"].fill(
                dataset=dataset,
                bdtscore=sel2_JERUp_bdtscore,
                weight=sel2_JERUp_weight
            )

            output["sel3_bdtscore_binning1_JERUp"].fill(
                dataset=dataset,
                bdtscore=sel2_JERUp_bdtscore,
                weight=sel2_JERUp_weight
            )

        if ak.any(basecut) and ak.any(cut2):
                
            sel2_particleindices = particleindices[cut2]
        
            sel2_events = events[cut2]
            
            sel2_jets = [sel2_events.Jet[sel2_particleindices[idx]] for idx in "0123"]
            sel2_electrons = sel2_events.Electron[sel2_particleindices["4"]]

            sel2_electronidsf = ak.firsts(evaluator["electronidsf"](sel2_electrons.eta, sel2_electrons.pt))
            sel2_electronidsfUp = ak.firsts(evaluator["electronidsfunc"](sel2_electrons.eta, sel2_electrons.pt))+sel2_electronidsf
            
            sel2_X = pandas.DataFrame(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel2_particleindices["6"])).data,ak.to_numpy(ak.firsts(sel2_particleindices["7"])).data,np.ones(len(sel2_events)),np.sign(ak.to_numpy(ak.firsts(sel2_electrons).charge).data+1),ak.to_numpy(ak.firsts(sel2_electrons).pt).data,ak.to_numpy(ak.firsts(sel2_electrons).eta).data,ak.to_numpy(ak.firsts(sel2_electrons).phi).data,ak.to_numpy(sel2_events.PuppiMET.pt),ak.to_numpy(sel2_events.PuppiMET.phi),ak.to_numpy(ak.firsts(sel2_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel2_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel2_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel2_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel2_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel2_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel2_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel2_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel2_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel2_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel2_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel2_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel2_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel2_jets[0]+sel2_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel2_jets[2]+sel2_jets[3]).mass)).data,ak.to_numpy(ak.firsts(sel2_jets[2]).eta - ak.firsts(sel2_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_electrons+sel2_jets[0]).pt*sel2_events.PuppiMET.pt*(1 - np.cos(sel2_events.PuppiMET.phi - (sel2_electrons+sel2_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_electrons+sel2_jets[1]).pt*sel2_events.PuppiMET.pt*(1 - np.cos(sel2_events.PuppiMET.phi - (sel2_electrons+sel2_jets[1]).phi))))).data))),columns=["nextrajets","nextrabjets","leptonflavor","leptoncharge","leptonpt","leptoneta","leptonphi",'metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel2_d = xgb.DMatrix(sel2_X)

            sel2_bdtscore = bst.predict(sel2_d)    

            if dataset == "Data":
                sel2_weight = np.ones(len(sel2_events))
            else:
                sel2_weight = np.sign(sel2_events.Generator.weight)*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf
                sel2_weight_prefireUp = np.sign(sel2_events.Generator.weight)*sel2_events.L1PreFiringWeight.Up*sel2_electronidsf
                sel2_weight_electronidsfUp = np.sign(sel2_events.Generator.weight)*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsfUp

            output["sel2_bdtscore_binning1"].fill(
                dataset=dataset,
                bdtscore = sel2_bdtscore,
                weight=sel2_weight
            )

            if dataset != "Data":

                output["sel2_bdtscore_binning1_prefireUp"].fill(
                    dataset=dataset,
                    bdtscore = sel2_bdtscore,
                    weight=sel2_weight_prefireUp
                )

                output["sel2_bdtscore_binning1_electronidsfUp"].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight_electronidsfUp
                )

            output["sel2_bdtscore_binning2"].fill(
                dataset=dataset,
                bdtscore = sel2_bdtscore,
                weight=sel2_weight
            )

            output["sel2_bdtscore_binning3"].fill(
                dataset=dataset,
                bdtscore=sel2_bdtscore,
                weight=sel2_weight
            )

            output["sel2_higgsdijetmass_binning1"].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel2_jets[0]+sel2_jets[1]).mass),
                weight=sel2_weight
            )

            output["sel2_higgsdijetpt_binning1"].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel2_jets[0]+sel2_jets[1]).pt),
                weight=sel2_weight
            )
        
            output["sel2_vbsdijetmass_binning1"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel2_jets[2]+sel2_jets[3]).mass),
                weight=sel2_weight
            )

            output["sel2_vbsdijetmass_binning2"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel2_jets[2]+sel2_jets[3]).mass),
                weight=sel2_weight
            )

            output["sel2_vbsdijetabsdeta_binning1"].fill(
                dataset=dataset,
                vbsdijetabsdeta=ak.firsts(sel2_jets[2]).eta - ak.firsts(sel2_jets[3]).eta,
                weight=sel2_weight
            )

            output["sel2_leptonpt_binning1"].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel2_electrons.pt),
                weight=sel2_weight
            )

            output["sel2_met_binning1"].fill(
                dataset=dataset,
                met=sel2_events.PuppiMET.pt,
                weight=sel2_weight
            )

            output["sel3_bdtscore_binning1"].fill(
                dataset=dataset,
                bdtscore = sel2_bdtscore,
                weight=sel2_weight
            )

            if dataset != "Data":

                output["sel3_bdtscore_binning1_prefireUp"].fill(
                    dataset=dataset,
                    bdtscore = sel2_bdtscore,
                    weight=sel2_weight_prefireUp
                )
                
                output["sel3_bdtscore_binning1_electronidsfUp"].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight_electronidsfUp
                )
                
                output["sel3_bdtscore_binning1_muonidsfUp"].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight
                )
                
                output["sel3_bdtscore_binning1_muonisosfUp"].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight
                )
                
                output["sel3_bdtscore_binning1_muonhltsfUp"].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight
                )

            output["sel3_bdtscore_binning2"].fill(
                dataset=dataset,
                bdtscore = sel2_bdtscore,
                weight=sel2_weight
            )

            output["sel3_bdtscore_binning3"].fill(
                dataset=dataset,
                bdtscore=sel2_bdtscore,
                weight=sel2_weight
            )

            output["sel3_higgsdijetmass_binning1"].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel2_jets[0]+sel2_jets[1]).mass),
                weight=sel2_weight
            )

            output["sel3_higgsdijetpt_binning1"].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel2_jets[0]+sel2_jets[1]).pt),
                weight=sel2_weight
            )
        
            output["sel3_vbsdijetmass_binning1"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel2_jets[2]+sel2_jets[3]).mass),
                weight=sel2_weight
            )

            output["sel3_vbsdijetmass_binning2"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel2_jets[2]+sel2_jets[3]).mass),
                weight=sel2_weight
            )

            output["sel3_vbsdijetabsdeta_binning1"].fill(
                dataset=dataset,
                vbsdijetabsdeta=ak.firsts(sel2_jets[2]).eta - ak.firsts(sel2_jets[3]).eta,
                weight=sel2_weight
            )

            output["sel3_leptonpt_binning1"].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel2_electrons.pt),
                weight=sel2_weight
            )

            output["sel3_met_binning1"].fill(
                dataset=dataset,
                met=sel2_events.PuppiMET.pt,
                weight=sel2_weight
            )

        if ak.any(basecut) and ak.any(cut4):
                
            sel4_particleindices = particleindices[cut4]
        
            sel4_events = events[cut4]
            
            sel4_jets = [sel4_events.Jet[sel4_particleindices[idx]] for idx in "0123"]
            sel4_muons = sel4_events.Muon[sel4_particleindices["4"]]
            
            if dataset == "Data":
                sel4_weight = np.ones(len(sel4_events))
            else:
                sel4_weight = np.sign(sel4_events.Generator.weight)

            output["sel4_higgsdijetmass_binning1"].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel4_jets[0]+sel4_jets[1]).mass),
                weight=sel4_weight
            )

            output["sel4_higgsdijetpt_binning1"].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel4_jets[0]+sel4_jets[1]).pt),
                weight=sel4_weight
            )
        
            output["sel4_vbsdijetmass_binning1"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel4_jets[2]+sel4_jets[3]).mass),
                weight=sel4_weight
            )

            output["sel4_vbsdijetmass_binning2"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel4_jets[2]+sel4_jets[3]).mass),
                weight=sel4_weight
            )

            output["sel4_leptonpt_binning1"].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel4_muons.pt),
                weight=sel4_weight
            )

            output["sel4_met_binning1"].fill(
                dataset=dataset,
                met=sel4_events.PuppiMET.pt,
                weight=sel4_weight
            )

            output["sel6_higgsdijetmass_binning1"].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel4_jets[0]+sel4_jets[1]).mass),
                weight=sel4_weight
            )

            output["sel6_higgsdijetpt_binning1"].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel4_jets[0]+sel4_jets[1]).pt),
                weight=sel4_weight
            )
        
            output["sel6_vbsdijetmass_binning1"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel4_jets[2]+sel4_jets[3]).mass),
                weight=sel4_weight
            )

            output["sel6_vbsdijetmass_binning2"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel4_jets[2]+sel4_jets[3]).mass),
                weight=sel4_weight
            )

            output["sel6_leptonpt_binning1"].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel4_muons.pt),
                weight=sel4_weight
            )

            output["sel6_met_binning1"].fill(
                dataset=dataset,
                met=sel4_events.PuppiMET.pt,
                weight=sel4_weight
            )

        if ak.any(basecut) and ak.any(cut5):
                
            sel5_particleindices = particleindices[cut5]
        
            sel5_events = events[cut5]
            
            sel5_jets = [sel5_events.Jet[sel5_particleindices[idx]] for idx in "0123"]
            sel5_electrons = sel5_events.Electron[sel5_particleindices["4"]]
            
            if dataset == "Data":
                sel5_weight = np.ones(len(sel5_events))
            else:
                sel5_weight = np.sign(sel5_events.Generator.weight)

            output["sel5_higgsdijetmass_binning1"].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel5_jets[0]+sel5_jets[1]).mass),
                weight=sel5_weight
            )

            output["sel5_higgsdijetpt_binning1"].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel5_jets[0]+sel5_jets[1]).pt),
                weight=sel5_weight
            )
        
            output["sel5_vbsdijetmass_binning1"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel5_jets[2]+sel5_jets[3]).mass),
                weight=sel5_weight
            )

            output["sel5_vbsdijetmass_binning2"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel5_jets[2]+sel5_jets[3]).mass),
                weight=sel5_weight
            )

            output["sel5_leptonpt_binning1"].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel5_electrons.pt),
                weight=sel5_weight
            )

            output["sel5_met_binning1"].fill(
                dataset=dataset,
                met=sel5_events.PuppiMET.pt,
                weight=sel5_weight
            )

            output["sel6_higgsdijetmass_binning1"].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel5_jets[0]+sel5_jets[1]).mass),
                weight=sel5_weight
            )

            output["sel6_higgsdijetpt_binning1"].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel5_jets[0]+sel5_jets[1]).pt),
                weight=sel5_weight
            )
        
            output["sel6_vbsdijetmass_binning1"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel5_jets[2]+sel5_jets[3]).mass),
                weight=sel5_weight
            )

            output["sel6_vbsdijetmass_binning2"].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel5_jets[2]+sel5_jets[3]).mass),
                weight=sel5_weight
            )

            output["sel6_leptonpt_binning1"].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel5_electrons.pt),
                weight=sel5_weight
            )

            output["sel6_met_binning1"].fill(
                dataset=dataset,
                met=sel5_events.PuppiMET.pt,
                weight=sel5_weight
            )

        return output
        
    def postprocess(self, accumulator):
        return accumulator

samples_local = {
"SingleElectron" : [
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/00406F80-7083-D74E-8795-49F431270448.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/00518AB9-C001-FB4F-BBD8-8CB818642CFD.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/005712D8-A755-DE41-8140-E772DE0AECFA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0093A31E-1326-3145-926E-B7C57A48E0E8.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0156340F-C336-1241-B2A8-A7ED20B2A043.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/01EAD3EB-94BF-2D44-B6FF-7661E7503B77.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0201470A-A7EA-764F-9BA2-4DEF0CC37E42.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/02A9B576-790A-9141-9762-AA2831FABF58.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/03472ECC-0FBC-A845-87F6-85BA3D045B89.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0391B391-CD8F-0E4A-872C-E0B00FCA5FD2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/040C6059-EFC7-E249-AC65-5E2F6FE388E8.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0456095A-061E-874D-BFE2-FBC1578F643F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0560B858-391D-1242-B3FB-3649B576FA85.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/059A9B2B-480D-964F-9B62-6B9500FDCA61.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/05BF8BE9-A1CE-F842-AD24-141FC2A1A0C5.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0600EBC7-D6C8-1946-A2B8-65D3D0997034.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/068FE8AA-A338-4543-8181-4687B3D703D5.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/07368BF9-C206-914E-8700-AA03337A049B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0745C857-9E1F-DA45-9A13-96A8D4100EFC.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/074FFE1D-6F75-354F-B664-52989853B986.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/07C71732-FF2B-D84C-A51E-DFC68A827B19.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/085F2824-0564-CD44-B713-E2205F6A3A85.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/08695D89-1DCA-9146-B360-E8A9D3E2ED91.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0A61CBE2-CB6E-774D-956A-29D4E5F4A7EE.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0B0A7F66-16A8-E642-A9A5-FD164C9DDBF2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0B80B1F1-7A7F-AC40-846D-CFCDDD05DE01.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0BE4DB7C-2398-E442-B3AB-AC4AC3DFBE2C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0BEA9DDE-EE29-6C4F-9181-F892B34B7603.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0D19EF2A-904C-5347-A8FA-8DCA87868F37.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0D3360E0-CCE2-3543-A953-620C362B10C2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0DEB37B4-BB35-054E-9C07-EC4221E6B173.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0E0F0650-D922-584E-9F3E-C4E992A8EE47.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0E9DDCD3-9C7F-214E-9B9B-32033D2A2F37.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/0F8521B0-5770-9244-915E-F643EDAA32D2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/11768EFC-1382-1B41-8113-1D18FC71B750.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/11B65D6A-EE0E-EC4B-8392-9B5B9F30948E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/122AB677-34D9-B744-B7EF-65AFA7DDAAFF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/12340C6E-9F9B-D548-ABCD-31D7D4DD019D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1354D433-BF36-084E-B92C-544161DA667E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1498486F-78DC-EB47-88C6-CB9731F2386D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/14CCEB5D-D3A7-8F4E-87B2-C05703CAD405.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/156943D8-0B36-DB42-AD64-7B5DA38DB5DF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1586BC26-EF29-034B-AC43-A4902F266A71.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/15C5ACB6-9C35-C548-9BFB-486375AC715F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/16756593-AF66-CA45-B9F6-677D4559877D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/16A9A2A2-DFBF-2E44-8D6C-D69F3FF9CCFD.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/172C3AB6-A094-EA48-82C2-312C587F810B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/17A63A4D-7D3B-9B4D-8826-80CD5DB51631.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/18541411-4F4B-334B-89B2-FD6E60C5E5B9.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1906B16F-F286-6E4D-8772-7FF81EF94805.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/19096B97-3AD2-FE45-B316-1991C4634BF7.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/19B10F2C-40CC-C245-A1C8-9A722EA672B6.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/19E16619-BDB5-FD46-9F16-8273ACFD5C0C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1A0481F9-CBE9-094E-A7EC-71720E338FC9.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1AE20B09-CC0F-774D-9AE6-71A9B14CC34B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1B1FC377-5B36-ED46-B468-88AF8A8AF30E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1BB5575E-4545-FE43-A3BB-37032DFA3EBD.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1BD912DB-A9DC-244B-AF20-F7406D805B5A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1C622A93-6528-CA43-B7CD-6C48F0081580.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1C9FE2CF-36A1-AD49-938D-26F24CA7D045.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1D1E459B-A4CF-D840-A207-42A894EEAEEF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1D6931B7-A652-4A45-8984-3BFD9EE16D84.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1E9D0884-D34C-AC48-B0CF-8800AF91C1D2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1E9E9AE1-D585-7240-8195-E050F34A7866.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1EA55646-335F-754C-BA97-D730EE291731.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1EB301B8-F0DF-7340-987B-58641312B381.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1F403F75-76E0-3C48-BB90-189FF77A17E5.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/1F60AD8E-9AB6-534D-9E3A-54C4954E90C3.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/202ADFD3-8A93-484F-8416-18D9BED9BF4D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/2104C1B4-DE48-704F-A32E-24BB26447D34.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/24D98199-2C12-2649-8475-888C0C793350.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/25CEB250-9B73-1644-BE9C-D638EC630653.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/26E40CC2-412B-824D-ABAC-EF03F6B58939.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/288F554A-12A5-544D-B7E4-29782FDA25EC.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/289B7B64-29E7-9C4D-9E2F-9322EDB66A7D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/299F9360-40E9-5941-BD14-374289B0500F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/29A176D4-51B6-364A-B916-3057A5063712.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/29EFFA6A-006E-3E45-A632-DAECB235D015.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/2A05AC1C-E7D9-BF48-ACE6-B35F9C8BCD08.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/2AB2E3B2-47F0-6142-8F93-6AFDE18B18EA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/2AC8AF01-D299-EC47-B6F2-490F32F2219B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/2AE354D8-757A-5B40-91DC-A497FF4EB495.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/2B36CBCC-5E7A-744C-B785-0AF275AAFCCA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/2B9CDE6B-BF57-A943-824B-35CB02C2D245.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/2C477666-3961-EE40-9ABB-9413C2BA6219.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/2EA54E12-8CC5-B842-8133-062E5B8E5180.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/2F32E3AD-6290-7C45-8DA9-DE53B01D73B3.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/2FD00ED1-01B2-CF4B-809A-DCF028B780FF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/30583D9F-D553-964E-BEE9-F98C24AB89F8.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/30598339-090E-034F-B8DC-8FC3E737934F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/312E5E53-7D22-6846-9690-F09AB1C5A1A3.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/31A241F3-D4C9-D94E-8B34-B426CA8E0649.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/31B8E512-A204-4145-A37F-A7517C927203.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/31BADB98-BCD8-4B40-9825-7709B08299BD.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/329E327C-3B0C-334E-92C2-1EB909137EDD.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/34CB386E-BBBF-ED46-87EA-103E6CC185BB.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/34EAF75C-DA7A-A843-A867-3431B97FA39F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3534209E-4967-D84F-AB0F-15BAFF357FF4.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/353EBFA2-20B4-364C-A605-07079C8C3F57.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/35C60E3C-E08A-544C-9912-908C1C597E4C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/365492D2-7E94-8E47-86A5-C9490E1C9EB7.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3661903C-9BE9-044B-91F9-6849FF39B6CA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3668E5C5-1E9A-6849-BA2C-07717C0A2ABA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/36D35EAC-3BDC-5044-BAE7-E79CB3791441.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/374C9420-29A1-744C-93D9-9FDC59A6CDC8.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/38A4D060-A677-1A4E-A51B-E7D9ECDAD04D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/38EF97ED-4BD2-9244-B2DF-D5F67A1A1DC6.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/396C0A74-ABA7-134E-AEA5-64CF1630506C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3998437C-F335-9D49-A58A-D77FAAA88A41.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3A3C3FFD-A94D-CA4C-9331-D69A8BBD4139.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3A6980CF-C706-954B-8F43-713570908903.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3BCFFD67-6F81-C34F-AB41-847027DF4211.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3BFC4811-566B-CA4F-B19B-C4B38BBF9639.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3C0D7450-430F-2342-BFCF-619B4F1D03B5.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3DC05B89-E5C2-9B42-A9D5-D313C1CEF278.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3E658EF7-F6B2-8542-AF10-304B350655F1.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3ED994B3-7933-7C40-AF29-847FAC984D19.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3EDC5465-057A-094C-B9D4-465D8443A5B9.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/3F6A495C-AC15-8148-8876-FEB0A1A71665.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/402A8D1A-A37D-9447-ACB8-411D1471991D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4095B633-B907-114E-82DE-A386EBA97082.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/41309394-7723-2540-B5A3-812B6393A43A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/41A092F9-9A52-7849-BE17-D643C70B5557.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/422A6C58-F739-284B-8B08-B46B29D7E0CD.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/426D6D9F-F08F-5749-AA5D-7B31EEEC67C9.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4299B418-248E-E648-A1C9-EC5309A43205.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/42FA2CFB-9CAB-2244-8CFD-A572DD73B43A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/43B00DA0-41AF-D042-9F41-7F95FBE5E59F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/43C65197-79F1-ED46-BA93-2C4E09DD5D45.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/43E392AF-541C-0240-BB31-62C8E6CB4E7F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/44F6400C-4CAD-E34F-9390-BF0F3FBB748D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/452CB7D6-C3C5-BE43-BE0A-E08152AA912F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/455EE3F7-D21E-9B44-BE5E-C732BF0ADD56.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4638A641-9114-344D-B0BD-118CB973190D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/468D599C-03EB-A64B-BC67-89B3B552A31D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/471C9C80-BD3B-CF45-BA88-8599C822A8B5.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/47D03AC6-462D-8045-944B-CAE235BA0193.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/47E91A34-7E96-8A41-B6D3-FB2FD9DA8849.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/47FD57C5-246A-724E-84ED-65B13FBDA314.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4800493A-3630-054A-A470-6FFE9C7B7EF9.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/49058605-6AD6-734F-AA8F-CC2E9350937A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/493981A8-92AF-7A41-8BF3-5107039184D1.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/49FAAA2C-84ED-864E-80EA-55A175C5D62D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4AE1DA71-05D5-E044-ACEF-5CAEDF8F0C0B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4AEA3462-C718-1B42-BC44-D5FD319816A7.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4B763B6C-6D4D-9644-B017-FE6B42FF71CA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4C3A1566-39B4-C844-9129-463D6D885117.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4C567043-645C-8744-9E57-CA08A9581A21.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4CCC9B83-5578-3647-896F-F50B13697072.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4D8B76BC-6460-4349-A8E6-716737F68F80.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4E5993AD-4D25-8947-9727-FF41361D3A52.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4EE67019-119E-6B46-8B15-B7D79EBB9D56.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/4F1BF128-FAF9-DD4D-922D-75780BD4D321.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/501A4354-5EAA-354A-8EA3-A172BF270326.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/50BBFFEB-AFDC-D948-8247-2019BD69B4C4.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/50FB010A-EA41-DC48-BA1C-03BA96249863.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5101C75D-475F-3849-A359-B01FC494B798.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/51036EF0-D064-C74D-BAA9-0015598DD7DF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/51A0A604-31F5-3E48-B8FB-6C8B3F680ACF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/51B1401A-433A-5C4D-AB45-75DB4F5A4843.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5261DCA4-D296-DD42-BB28-5A6DFC2F9F97.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5369CE02-A723-4349-94EA-DA2F86C22C9E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/54067BDD-326D-DC4D-99FF-C62B368AA847.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5416D54B-F2A9-7042-8107-2FFE8609DC36.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/54226DD5-CD7A-9A43-8B79-C6F6389D0C7F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5482B2CE-D532-074E-AAE9-D21372FA3A56.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/55167C12-9643-AC46-8834-F8A83C5A6BF6.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5627014D-82E3-BD47-BC68-8F13134B9177.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/565B1708-EE04-D849-A300-1A6EC9FFB9DC.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/57391810-B5CA-5148-8CF9-FB7318A851B3.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/588B1F7B-3E24-6947-A3DE-ED0831FD6C37.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/58B3B2AF-11C7-694E-A39D-29172B4971EB.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/59193DBF-CEF8-064F-B44A-987DA423590B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/592F9AD5-7F3F-A64C-A2DD-DA5D4C69CA92.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/59722D78-EED8-4F4E-B7F7-30C971A69A80.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/599FAF9B-1A8C-834B-BA3A-5EC74104F334.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5A7C9DDD-89A9-884F-821B-B0B2D10FE399.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5ADE31DD-0F1E-4147-B739-47A4983B6DF9.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5B394585-E5F0-044D-92D5-511AFA1C290F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5C6FE773-FB0F-DF4E-AA00-26F23FC726B0.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5CC83044-4544-B244-8F48-34900C6B05B8.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5CF86C8F-A4D3-AF4E-BB14-F98978D0E0C2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5D546170-4CC2-8047-B7A2-DFA78CD0B796.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5DC146D2-8EE0-734B-A98C-FF046331BCB3.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5DDDD6F3-3A7A-5042-B125-1E7293C64CE0.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5E8D1B09-FE30-214A-85FB-DC380DFD8A0E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5FC6D4C8-7C25-0D4B-9802-C89AA5F298AA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5FD0FB93-9596-0D4E-80C0-E14E93F4AE02.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5FF5DFDF-3A06-9449-AC7B-40D1A65ACF6C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/5FF6F6D1-89D2-2B41-8D13-9161F4D969CD.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/61163327-B48D-B84C-B9B6-718E787F508A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/61B4057A-C0AD-1E49-9687-B302FACAB4E0.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/62DD0A02-B6A4-CC4E-8114-D9C80879D9CE.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/630D9647-7E3A-0C4D-836F-D7069BFE98CB.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/636D149E-9B84-8040-A662-62D7323AF84C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/63AA3C41-D10E-EA42-B56E-A153B031F06A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/64735898-795A-DA48-B1C8-E5973308DCFF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/648A8B49-B1AA-AF49-A689-190D492AD117.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/65275DC0-CC85-A745-8A86-554631746543.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6548AC94-1D4E-1A40-85FA-7250A162F400.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6577E09B-FFE6-B44A-B89C-8ABE6DC5245E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/65C93760-0935-B641-98F7-63F88BB46F4C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/65C9F1B6-EFA8-724F-9165-C5B909BF8C8F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6676EF67-5115-D04B-A941-B27E10F4F56F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/668B508C-9E6D-0A4D-AA2F-9F0F87C79F0A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/670E43C7-359C-FA45-84A5-43F07F1F0F9C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/674CD296-D0D9-FA49-B7D3-FC0EFCA8DDDE.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6805CC1A-6B1A-5044-83F6-9658DE602A85.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/69F2EFAB-0664-EA47-BB96-FC1AAFD1E9D6.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6A2102DB-55A4-6B43-891F-6B7AB6DE88E2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6A5190BA-7564-6644-B56D-EAF5E321DCC1.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6B336E06-5895-C64B-92BB-D6AA3B88F7E9.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6B988E2E-9E9B-5D41-ABF2-0DA28EFCE109.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6C322BC9-180D-654E-93A8-92BE581B7DDE.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6D5E86E1-9B77-594A-B4F8-DBAEF00E25CE.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6E12FF60-020D-3B47-A391-04EF83FFEF02.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6EE5F862-32BE-174B-ACAF-7658B830A0EB.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6F79B8D8-AB14-9C4F-A3F8-8C3420592E52.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6FC72B0D-9AC5-7445-BB81-9AE2207E2ADE.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/6FC85334-058B-FD44-A465-61D337F0F63D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7099C3E4-CB90-734F-8C13-8DFC95D5BCA3.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/71957B52-D338-2345-9ED3-77B224E5C32C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7222BC00-0952-3D47-A5EE-2478418AF001.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/734A8849-21FF-844F-871D-41B537E6CB05.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/73CA981F-03FB-6B47-9078-D44CA909BD0C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7405CFA4-B92C-EB4E-B236-C82AC31FCF0A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/74E71684-C239-D64C-ACC5-DED88659AC42.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/75CBFF47-8E27-084F-A0B7-510D3B37E0E3.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/770028FA-6698-D04B-BE84-CEFCEFA83041.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/77747B18-CF3B-FA49-AA93-2080C8455582.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/777945B0-97DB-A340-9FB9-E45DC3DC2F69.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/779C9B59-F473-2143-8CA2-5750D1ACDA68.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/77FBAFF1-C7E9-6541-9B0B-E2A76A1EA0D2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/78BAFCBE-346E-424F-AAA2-B3312D364FDE.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/78D25636-AAF5-CE40-8C58-E5964258A524.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/78E08653-AC9A-7A4F-B58C-B45688AB55CE.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7945B514-DC21-6E40-8002-583C273ED61B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7AF670B3-BB25-B94E-B725-90494FF51C1A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7B003B6E-E1CA-6D4F-9BCD-BDC44C7BBF11.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7C407625-D83B-0E4C-8590-7D9BEACC7C65.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7D5B1FF9-B956-3B4E-A793-A7E73DBA4E65.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7E1813EC-2659-F340-B99F-411724CB9160.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7E3D6EEC-EFE7-924A-9C2E-4BEF9BDF9E2B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7E676600-5E77-EB47-A85E-1B1F6707A890.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7E86C060-BE21-CB46-A9C3-C3A75779DDCF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7FEC9187-64EF-AF40-BBF3-1D94F9645304.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/7FF8956B-8048-1849-B974-0BB3CEB6A950.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8089153C-B503-4548-8F80-06F55977910B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8091C050-3B27-1945-84D7-85BED88B8FFC.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/80A07488-D292-A041-B236-C1E918358086.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8211630D-1E4D-3540-B931-9858AC128E1F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8246C19C-FEC7-A948-AF9E-A7ECE45AA21B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/831B428D-3C88-2046-AB9A-A71FB706A9DA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/832A95C1-777A-C94F-88ED-8DA9DE901754.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/834A7A94-CF92-2B45-AF3C-E748FF887AF3.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/85165689-FC86-7F4F-B71B-E90C9A06BCE0.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8521F0B3-6CFA-7C4D-B9C9-72148C7FAB2E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/854400C6-7EA3-7F45-AF46-A1D382E1FFA2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/85CF736C-46B9-0547-B980-B4F53BD64405.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/86FCCAC0-E439-2C4A-BFC9-13029BDD8431.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/882BE27B-4C53-D64C-A663-1A0A5390D879.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/887B7CC3-DB70-B64D-A43E-5434D6A5BB35.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8885BED7-15B2-B945-914C-C95FFA67BA6E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/88D7296A-853F-E440-B1CA-F5F2E5B8EB15.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/895A7F3B-9E1C-0349-BDBA-BABEA391F09E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/89F85466-311F-E048-A950-F9582099E303.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8B085BFD-9F70-B244-8E24-0CAA6BF4AECF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8B5B43FA-0C26-0649-B03C-EDBD3FE74012.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8BA8A2F6-5E13-AC4E-A98E-EB69DF70E392.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8C1B7286-8537-0C4B-BA54-8C7651303409.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8D1C2BAD-670F-914E-985D-4CE2D69DFB08.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8D614E8A-021E-E040-AEBA-5FB59FF9385F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8D8BC583-9112-2E4E-BEB7-CCCA131FED95.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8D8ECD9B-89CA-4145-9660-4EE036876C70.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8DF45CB4-BB08-054E-BBEA-B264423E6AB0.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8E39FF54-BB03-0F49-A74B-5A52533E7ACB.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8E4E41F9-FCAB-A34E-8D79-F3C2E76759FA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8EE1AAA9-A485-F24B-BD04-1100E8275C91.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8F0889E1-C79B-4344-85FF-07EDCDE5CA83.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8F883EE8-B7D4-904E-A511-6536F689A700.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/8FFD944D-35CC-AA43-AFCF-51DE5D5264D2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/907406FD-AB7D-7A4B-8C42-4BC3B85849C2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/91A8661B-37AC-4840-98D6-216F50D70387.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/93C0A632-793A-FA4C-9CA7-21EBD9202703.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/947C376C-218A-7B49-AA38-BBA98DDF40C6.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/94DF495A-ACD6-1E4E-BB99-A191538DB0A9.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/952DAF45-AFFF-1641-A82B-D4DC33CCA526.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/954D99E8-1411-764A-B443-CF0CE05172DC.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/956AE607-98F0-A54A-B19F-DFC390D76EF5.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/95C9BA58-E9BA-F243-845A-6B835CC97BD8.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/98BBB482-211A-4F43-9670-D60A697478B1.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/98F2E8AE-DD19-484B-8B0A-AE0569F7B171.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/99E44A75-8BBB-FA4B-B4C9-09E3927BA891.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/99F6E8C8-4166-5647-81F5-444043F24923.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9ABE5270-E3B7-BE4A-A892-D1C24D401C21.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9B386D86-8FE6-2B4F-B832-A30CADD528E4.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9B589371-E753-1241-91B9-E0F856AD2CE9.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9B8C8771-4326-8E4C-B770-AC7CAEED7972.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9C509FC5-FB3C-4C41-A9C0-6EB7DB7EBA5D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9C58016B-336D-D343-ACB3-140764486CB5.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9CC99A0B-0618-A044-A5AD-23A46DD2D39E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9D36BBC1-48C6-9E4C-92BE-899F87090BBE.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9D6E2D81-1647-2747-9369-8F9BF60B36E4.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9D743DE7-24B3-2B48-8E82-232958A62E8C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9DB2F132-AB01-7040-9F10-F8D294026436.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9E7515D7-26DB-3249-AD58-F2AD6462D5E4.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9EE3333E-1B40-5147-BF55-50A6B30F714D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9EF53D6A-5841-7144-9BB8-A2B01A24D0D4.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9F50ADDF-3DBD-9642-B1E6-65F20AD66598.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9F78B65F-A706-7745-878F-08E4B83BA3EF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9FCC4F6A-006A-B844-9DFE-81D637BDAF4D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/9FE22EC1-66C3-DF4F-8B94-6946501D9F7A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A000B86E-A75F-7549-A06F-7C19B723ED14.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A0179389-0560-C147-986F-2A0AFBF85175.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A060ECF3-44B5-0744-99F2-10919CB66BB1.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A072DDBB-BD72-4748-B679-DF4A94FC2731.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A087F674-68BF-AF45-8E15-CCC6028B2576.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A0A4CB42-5B86-5E46-9D69-3CC6BDD144CA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A10A1F56-D828-C343-AAE2-CFABDC2AB48E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A10CE2DF-E000-134B-8891-E144290CBC74.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A195C2FB-0498-1348-A3D8-A4DCCB246E5D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A31FFD55-DC21-C541-90F1-CAA9C69951EB.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A3FDBA64-4103-1E42-AAB5-F9A8A9882F9E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A41C89B0-F9BC-7344-96E6-69FB691CC65F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A492637D-E7EA-0449-9C59-A173F215D367.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A49BDE29-E5F4-154C-8AE2-0E1A61238D60.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A4FD2136-E33A-1343-9B10-FB349F76AF93.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A54D3203-4296-A444-90C9-3222D94DFF1E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A5558E22-AAE5-4142-B682-78429571052A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A577D301-DA2C-5E43-B632-439A1EABA3C5.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A5DCF875-C8C7-1649-8527-21FBA592F93C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A63D4B67-4DCB-604E-B209-EED86D2CE696.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A6C78DE2-B20E-164D-B872-8AD76E834DE0.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A7261CAC-EAE2-9A47-AA8B-C51420837582.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A7A6DBCC-4B21-F241-91CB-BDD056CF9ED7.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A7DA6992-50AB-8645-85C4-24976D4ACAFF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A88562C6-0E25-A54D-ADCC-98AB99C62365.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A96C2EF5-C5E9-084E-AC3F-D36D39EDA96A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/A999573C-FBD6-5843-93E2-7BD1595DCBB6.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/AA9D0D9A-F249-BD47-93EA-45E35DF69653.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/ABDFFB20-B3E7-C544-A70A-F81B561235AC.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/AC9A488C-E402-D646-86A8-6A087CBD4140.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/ACFAC4F3-D8E2-F943-8E4C-9737F393B90C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/ACFBFDBB-DE20-A84F-87FC-FF7015BC2403.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/AD446F6A-A99C-4C46-9989-8781AAA537DE.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/AE83E80C-23B4-0E4A-B997-5C35145FE004.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/AF46DDEF-B00F-7D4C-9B17-9A0F9F11043A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/AF4B15F4-A11E-4B4B-9AA0-93DB28827E57.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/AF81D2E6-342E-7147-AC1D-ADAFF63B583A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B151FD17-7F84-D946-96A1-2C31AF639B8E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B16E1841-1340-7C47-B2F4-59431DDA2A55.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B1775694-F485-7A4D-8029-A57C5C0FB041.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B3E8DB25-D5B4-514C-8301-C35CF69E61B0.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B4F234BD-C3A7-2948-9C5F-F8C509BDC4D7.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B538D4E4-9E8C-1C40-BADC-43E3425DBB8C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B5C22796-A504-3947-AB48-897B1A09C88B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B689989B-E126-3141-AE79-D39FC85F8489.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B6C987DC-3945-D04A-8732-80A898272C29.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B70EFC66-C6F3-2F44-BE3E-43AB354A7880.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B7AAD791-1DB1-EA43-BEA6-E3609D4170CE.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B7B3D772-C182-AE41-9955-2B197B87AFF8.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B7F2548D-7DB4-5D41-830C-26D13199E07E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B86C69E3-1889-5042-9053-CA7CE8105181.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/B8F2411B-A1EB-A146-89EF-FC0DF23E5C71.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/BA260F73-22EC-6042-BE70-20D6AA3C716F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/BA6CC564-204D-F646-A098-48582B490B6A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/BA99198E-DBC8-814F-AEA2-04F77BBF26AF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/BB6CC8E2-C16E-1045-A3C9-78BACD26BFF9.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/BB6D24A4-9686-AC49-A999-368731C7CBCA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/BC2C1B66-83B7-0E4B-AC75-CA33FC89DDCB.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/BE31B8A6-86D5-BB4D-B063-2BB5B6F70894.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/BEA4C7C9-A3AE-A84F-9A0D-2B673E09B9A4.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/BEAEFAFC-8D68-1B48-AC43-104702D74E38.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/BF02A9B5-63C3-9042-A7E4-61CD050B94A6.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/BF866CC1-98D5-C74D-A209-D509D7CA309D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/BF96FE11-3194-1B4F-B2B1-4B6DD0191A7C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C00CB48E-952D-F445-B5A0-F8A52C1740C0.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C1006A47-9519-7C42-B3D3-A140704D155C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C15AF7C7-B1E8-874A-B47E-D6C8691CCB07.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C25FF17A-B64D-5E4A-8AC1-15187DA11A1F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C2AF04CB-F2CD-B947-9B26-68A05373ECAD.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C330644D-F164-4448-AD70-634508AB65FF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C34A36A0-C476-0646-88FB-EBE4025D3189.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C43EBD69-0017-7F43-9431-8E65A3AFB542.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C4D31476-4702-CB43-A8B3-8C4AD9FD7E34.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C532A1EA-655F-CA4B-B0E7-A6AD14D08B01.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C546E75E-3669-F441-AF25-2CE0188EA6B2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C93D0569-0832-2849-92E2-EBEE3F05A75F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C988C8DB-CE37-E540-9E03-C5E5BB4429A5.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/C9D9BADA-386E-514D-BDDB-D618F4F187B7.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CA4B5810-8D1A-EC47-B3C5-466274C16645.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CA95BEBB-BD20-A345-B517-6D32E443F129.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CB4CCA7A-F632-1F40-8937-6E8FBF17203E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CB4E572E-22C1-9443-BD9F-3DCD79908FD4.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CB7BF1F8-DBBD-874D-BC11-D0F7871CAF66.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CCC4D6FC-1C8D-EC47-B5B4-28423B1803DF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CCC5823F-3F0B-8649-AE10-B7AC896008E1.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CCC5B162-06E3-2A42-81FE-92E3A6004CDB.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CD7F433E-F782-864E-B8ED-8A3E76950164.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CDDB3A4D-910E-E146-BE8C-AF9878D9B5B8.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CE187FC0-27D1-764B-BB24-5C7F8418EDD6.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CF9724FA-9B71-1744-AD92-A6C7539819FA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CFBCA453-B91E-6241-AFBF-6CF183C096ED.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/CFF5AAAF-AED1-2649-AE6E-2A4550F94673.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D089EB55-5567-DA40-9CB0-7E58441E7B2B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D11737FA-B488-7D41-8D55-B2020C615568.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D1AC23F3-50D1-544B-A953-582474C1C129.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D23F11A1-5EBD-944C-91E5-061B6C12E4D0.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D27C0D05-06EE-1D41-A387-4936D5C0EB75.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D2A7C4F3-BFCE-654D-9452-521071EDDCE4.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D2D2CEAE-A5CE-A04B-B81B-536ABB14BFFA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D3816C9A-9CF3-5B4B-874B-099007494AD3.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D3F6924E-0A11-1B48-B19F-53F696CC8D65.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D4604979-85DE-454C-B1BA-C12DD9CB8885.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D4AD64B9-E1C0-7A47-A019-177A80536947.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D50F5D26-8161-2C49-829C-801DD131BA7A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D582F0FF-D728-804C-A6D0-C41D7B9B1491.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D5DD2A55-CAD5-1344-94F6-A7E167B6F39C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D621F4C4-ED52-9F48-938B-46641315AB8F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D69F6362-8722-524B-AFAD-EAC281B32330.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D6B06099-81BB-DE4C-8314-A73D619C3232.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D7D5BB71-C4B6-5E46-9182-06BC71B4E585.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D8363897-8101-1F4F-9E7C-296BB2A7B8D4.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D8458C85-92B1-F646-958C-40D78648F51E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D8650C7E-FD08-F94F-A7C2-AB11A0042F6F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D8C82F21-CD49-3C4D-BD90-33F1613369A4.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D8FAD8A5-4065-4B48-8C39-76C9E217132C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D9AE3B22-4A97-164E-8635-D16400A2F32E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/D9F74687-1F0F-F147-8C62-D53EE93F7B56.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DA0A7930-DE7B-6141-AF4E-F9D06865E07A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DABA9A9C-9058-0541-B079-46ADCBAFFB9E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DB1C9C47-C0C3-5A4A-A1A1-E0212A6DAAA0.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DBBF126A-EBEE-8A41-B918-C0756E263781.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DC663292-716F-2A49-AFD5-D73F20DC0AE4.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DC75BF22-1FF1-6641-951C-D493F8033E4A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DCB2EA4A-2822-2B41-9EE7-8F0FF60DE61B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DDB4F7E2-C332-764F-9DCD-57FE35214C51.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DE481FAE-961C-D442-A37C-3AB41456AC66.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DF8D44E4-BF83-EE4D-B81C-6127A086C3C3.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DFBC9D7A-6E2A-6647-9853-AE3E560E3E53.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DFE2DCA5-A771-7448-8506-0D2740C41E39.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DFF3409E-2275-704F-A18A-5E9E41FD4776.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/DFF57A2F-0223-A74E-8095-257CE9755AB7.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E10D0CF9-1A6A-BE49-A769-6674FA825C77.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E12CFB28-D365-6547-A8B5-B4177C189DBA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E1CFE337-9EA5-AF46-A15D-D898630E5564.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E338C697-ABF9-0340-BCBB-9E682476FF45.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E3B70EA3-4108-8946-A018-DD42B727BC81.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E42593E7-1168-A540-B0DA-F125DD5CB0D7.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E486DBDE-35B3-C242-B7C2-897A6ED8572C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E4BEA898-1092-D34B-B82A-E70804FA4D95.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E4EC6CCB-3782-E145-AA5F-21CC8B18CA61.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E4F15409-AA04-7744-A7E2-DAFB531CE641.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E4F2AAB1-0763-DC45-966C-3387251D2815.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E60D2890-C62D-B940-A811-0026B3A7C57D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E61922D6-9ADF-A24F-B5C4-A2DC99610E7B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E680A96B-D4AE-354D-B7E1-BE908FF501AB.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E854F380-6F39-DF4C-8459-CC2CCB0A0D64.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E85626D2-633D-0645-9779-D9228A67F20E.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/E8E8522A-2D9E-C946-B064-2AE77206E7A2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/EA370243-CD3D-A249-9DA3-8ED6B1CAA6BD.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/EAC30FBC-777A-724C-A1BD-907BD842DC1A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/EB250A66-8F9B-DD4A-B3F9-0E61677286D2.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/EBBB61E2-3368-AE46-B245-446F8CCFC0B6.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/EC4E2828-1392-924D-9F01-0D77647873FC.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/EE4B4D9F-C51E-9945-81A7-7932462A6144.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/EE6D5BA0-C698-464C-AFF6-606E54749C35.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/EEDA2F2C-185F-9048-AC62-5F69DF139770.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/EEECA562-0093-6B4C-AA11-6875D594A852.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/EFC8894F-D06E-6F4A-AE0E-500B693E7FA3.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F012F7AD-FFA2-4342-B811-2FBB538D4E6A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F0356F98-9FEB-F042-BD9E-252ED253A349.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F1685621-2D00-5949-84C5-673476883AA8.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F17689C1-08B7-8D4C-81A1-F06E4EE3B2B1.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F19774F9-5135-7642-A893-AE666A0CF394.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F1A06213-2A8F-144F-878F-FD1C58CD7F4B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F41A5313-B380-2D41-89CD-489CD8375B1B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F42D380D-96D3-6345-A7CF-A564CE26C2B7.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F433AB65-CB8A-0C4E-85BC-CD68C12EDB5F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F474DFE4-8DFB-6E4D-ADB1-2E2A1C74E404.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F476C7F1-5EEB-A44D-B31E-D04A9B168713.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F4CF4708-8D00-8043-A750-3B96BD1D5599.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F50B322A-29F0-C84E-82D2-AABDFAAA9496.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F53FD642-9BD0-EA4C-A19D-4BE13B23FA8F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F54A04DC-57C2-0A4A-BC08-BE67B7D15D0C.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F59C65FC-FABA-5E4E-ABB3-C6AED5045CFE.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F5F8AADD-09B7-A04C-874A-F20805252EAF.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F681E14B-7A56-DD4C-AA49-18FFD1E6DE8B.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F72E59C9-C8CB-104F-9160-7D777568F8F0.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F753A434-7079-584B-94BE-E854B3D437D3.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F7DB5190-D9BD-FC40-B83B-89CC35647BA0.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/F8CCE6A9-219D-D74E-B04E-A8A0BE6C683D.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FA1A4E59-4350-5040-8AB7-7EF830AFA78A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FA9A8BD7-1754-B14B-B08C-A4A8A41F91F8.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FAAA3332-67B4-4E49-AC47-1DD2FB2653D8.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FAF3CE2D-C56F-8A49-989F-45E0B547FA00.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FB0C9CBA-6846-CE46-B2A5-414CD9C159CC.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FB585F64-DF88-6247-A3F6-D5CF704FA023.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FC83A351-A31F-A448-B497-8117D8F2CEDB.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FCD22D30-BC0A-804F-87A1-09F8D5D0C7C8.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FCD59491-8863-0A4D-9F63-6E7457CACE01.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FD812374-389A-B745-91A5-DB00A00BBAAA.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FDBFAC07-DB59-2647-BFC3-BAAAC291223F.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FDD8276A-7968-9C45-9D38-280DE9DDF650.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FE6975B4-F048-0F47-A784-6206B98359DB.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FF06EEDC-7CF4-C74D-9A03-8D475F04E08A.root",
"/eos/user/a/aejung/andrew/data/nano/singleelectron/2016/FFC104D0-5F59-A347-8AD8-FF445E556C4E.root"],
"Signal" : [
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1758.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1699.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1722.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1629.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1418.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1135.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1328.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1458.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1662.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1597.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1636.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1833.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1913.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1251.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1278.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1774.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1503.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1556.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1978.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1992.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1134.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1675.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1740.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1718.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1815.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1777.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1546.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1138.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1893.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1532.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1725.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1345.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1344.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1642.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1267.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1190.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1579.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1637.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2018.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1767.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1789.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1442.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1943.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1879.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1596.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1847.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1958.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1272.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1504.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1179.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1160.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1163.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1346.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1181.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1529.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1902.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1875.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1153.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1956.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1338.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2075.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1628.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2099.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1850.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1414.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1343.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1178.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1429.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1268.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1854.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1150.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1680.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1895.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1177.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1733.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1426.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1439.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1105.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1479.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1492.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2003.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2039.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1422.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1271.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1340.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1543.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1672.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1342.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1233.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1207.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1322.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1889.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1238.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1242.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1535.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1436.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1991.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1215.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1955.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1708.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1214.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1337.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1928.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1878.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1244.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1234.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1360.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1147.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1679.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1404.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1555.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1169.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2077.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1240.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1441.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1424.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1352.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1856.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1187.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1141.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1145.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1469.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1674.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1136.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1247.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1901.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1876.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2009.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1806.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1144.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1919.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1972.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1195.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1417.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1170.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1625.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1581.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1766.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1639.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1741.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1312.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1127.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1405.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1110.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1967.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1668.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1761.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1132.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1111.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1447.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1926.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1388.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1624.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2016.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1123.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1476.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1120.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1601.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1489.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1949.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1142.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1953.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1845.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1334.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1444.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1670.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1482.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1182.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1318.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1329.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1916.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1550.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1803.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1372.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1186.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1677.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1885.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1689.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1291.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1133.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1721.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1496.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1865.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1369.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1243.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1859.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1952.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1844.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1869.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2022.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1813.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1763.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1923.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1109.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1402.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2058.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1918.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1219.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1461.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1622.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1456.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1748.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1892.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1164.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2076.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1408.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1729.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1537.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1306.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1401.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1216.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1633.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2053.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1805.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2006.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1656.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1295.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1498.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1600.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1365.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1829.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1823.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1791.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1920.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1523.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1632.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1775.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1807.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1698.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1925.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2094.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1986.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2038.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1987.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2026.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1696.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1985.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1113.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1704.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1950.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1488.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2037.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1659.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1778.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1522.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1694.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2010.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1302.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1592.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1896.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2044.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1410.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1914.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2047.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1941.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2040.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1661.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2004.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1922.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1685.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1969.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1936.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1420.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1904.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1799.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1994.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2036.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1225.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2042.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1654.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1827.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2028.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1536.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1880.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2031.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1782.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1907.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1783.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1339.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1931.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1979.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1121.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1517.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1692.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1982.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1245.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1882.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1676.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1248.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1801.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2008.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1580.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1413.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1871.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2027.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1933.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1386.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1712.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1849.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2090.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1239.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1130.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1831.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1818.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1650.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1795.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1486.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2021.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1154.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1316.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1828.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1921.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1808.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1749.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1824.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1298.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1519.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1608.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1938.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1392.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2050.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1825.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1470.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1968.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1964.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1915.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1793.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1792.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1518.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1487.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1996.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1886.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1836.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1948.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1669.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1771.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1617.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1924.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2057.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1411.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1905.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1327.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2049.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2011.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1598.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2014.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1500.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2032.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1226.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1714.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1412.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1166.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1375.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1820.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1937.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1997.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1379.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1619.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1810.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1116.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1995.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1321.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1483.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1317.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1690.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1653.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1809.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1647.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1491.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2048.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1115.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1573.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1236.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1638.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1509.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1516.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1703.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1940.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1657.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1927.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1648.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1508.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1311.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1570.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1484.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2054.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1314.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1645.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1819.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1957.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1300.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1607.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1750.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1332.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1762.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2025.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1319.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1303.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1156.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1765.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1746.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2013.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1755.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1802.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1881.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1465.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1382.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1811.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1129.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1594.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1984.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1155.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1839.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1609.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1814.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2059.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1735.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1506.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1935.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1292.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1620.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1117.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1700.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1588.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1159.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1126.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2046.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1151.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1122.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1394.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1222.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1511.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2052.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1203.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1784.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1932.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1383.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1593.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1988.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1962.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1211.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1963.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1971.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1945.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1897.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1899.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1641.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1157.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2055.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1521.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1205.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2035.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1326.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1884.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1325.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1180.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1683.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1530.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1947.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1912.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1569.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1976.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2034.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1898.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1751.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1304.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1834.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1727.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1227.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1780.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1599.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2023.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1837.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1840.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1860.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2056.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1217.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1906.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1167.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1398.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1146.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1687.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1841.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1772.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1946.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1158.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1589.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1224.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1713.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1660.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1396.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1883.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1575.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1728.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1939.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1705.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1864.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1688.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1323.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1644.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1551.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1998.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1531.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1485.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1357.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1453.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1684.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1241.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1583.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1450.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1843.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1255.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1908.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1139.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1513.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1591.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1826.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1400.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1104.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1610.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1942.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1309.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1393.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1313.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1838.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1279.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1747.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2002.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1395.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1230.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1611.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1961.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1397.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1125.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1989.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1384.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1773.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1737.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1148.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2001.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1835.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1162.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1406.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1756.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1331.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1615.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1663.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2051.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1798.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1911.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1367.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1618.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1752.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1586.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1266.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1520.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1499.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1315.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1821.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1299.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2029.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1846.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1744.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1960.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2020.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1448.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1137.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1655.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1228.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1466.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1320.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1526.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1606.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1754.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1501.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1868.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1634.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1861.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1612.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1124.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1757.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1709.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2045.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1296.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1128.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1736.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1867.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1977.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1232.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1734.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2062.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1745.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1390.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1682.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1707.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1576.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1490.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1212.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1980.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1475.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1407.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1389.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1776.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1463.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2033.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1917.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2024.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1930.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1223.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1770.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1785.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1894.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1858.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1419.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1460.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1459.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1446.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1261.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1830.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1888.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1515.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2030.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2007.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1753.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1566.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1804.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1206.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1524.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1107.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1161.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1308.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1464.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1781.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1817.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2068.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1788.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1587.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1235.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1387.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1563.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1229.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1437.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1800.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1590.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1454.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1118.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2070.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1209.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1443.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1702.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2012.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1627.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1165.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1231.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1726.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1288.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1285.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1686.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1377.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1794.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1112.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2067.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1143.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1678.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1204.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1635.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1260.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1467.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2017.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1764.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1545.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1478.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1063.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2005.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1822.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1693.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1416.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1082.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1585.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2072.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2089.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1552.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1577.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1560.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1305.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1539.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1069.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2095.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1064.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1701.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1790.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1072.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1253.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1114.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.2043.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1934.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1873.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1101.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1098.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1103.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1097.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1096.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1099.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1102.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1100.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1095.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1970.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1505.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1614.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1090.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1067.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1068.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1073.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1070.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1076.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1066.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1071.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1079.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1086.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1075.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1083.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1084.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1092.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1088.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1093.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1094.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1062.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1087.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1091.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1056.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1077.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1060.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1054.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1052.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1078.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1058.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1055.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1080.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1081.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1057.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1074.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1059.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1053.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1051.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1085.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1743.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1061.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1089.root",
"/eos/user/a/amlevin/data/nano/ewkwhjj/2016/ewkwhjj.1065.root"],
"WW" : [
"/eos/user/a/amlevin/data/nano/ww/2016/004E4C32-9FC4-A24C-AF08-374E95DF0B4C.root",
"/eos/user/a/amlevin/data/nano/ww/2016/044FC323-464D-464F-A7B3-D2BEA2EC53C7.root",
"/eos/user/a/amlevin/data/nano/ww/2016/10DA70C0-DD84-2049-82FB-DCC828A74F5D.root",
"/eos/user/a/amlevin/data/nano/ww/2016/13A6B559-A234-3446-97C0-1AD9550E912B.root",
"/eos/user/a/amlevin/data/nano/ww/2016/1D6F3E2F-7A4E-6841-926F-800FDA9C16E7.root",
"/eos/user/a/amlevin/data/nano/ww/2016/21F66F87-19EE-404F-93C5-17526EA1C0ED.root",
"/eos/user/a/amlevin/data/nano/ww/2016/237E819E-2AFE-A445-98E4-1FC0CB5E1D61.root",
"/eos/user/a/amlevin/data/nano/ww/2016/309C4E04-65DC-9B42-995B-5D45916B2D06.root",
"/eos/user/a/amlevin/data/nano/ww/2016/37C01D05-8328-9A4F-AD4B-1AD09B4288AA.root",
"/eos/user/a/amlevin/data/nano/ww/2016/3E60A247-EEB8-D54C-98DC-1AE489D38F19.root",
"/eos/user/a/amlevin/data/nano/ww/2016/4A13C66D-3984-7148-BAA2-D01EA8AA2685.root",
"/eos/user/a/amlevin/data/nano/ww/2016/5093F917-D106-F342-AA0E-9BCA519D87C0.root",
"/eos/user/a/amlevin/data/nano/ww/2016/5925A9F5-0C75-B54F-950A-501D77409970.root",
"/eos/user/a/amlevin/data/nano/ww/2016/6DD0AF9E-308B-1C40-B1E3-EDD5BF96DE91.root",
"/eos/user/a/amlevin/data/nano/ww/2016/6F56FDE1-5554-4B49-8B64-A841D42BD8C3.root",
"/eos/user/a/amlevin/data/nano/ww/2016/76541735-6144-DC45-8D92-A68A9A8A51FE.root",
"/eos/user/a/amlevin/data/nano/ww/2016/85D3A6CC-3F6A-B946-A843-1853EE2B96D0.root",
"/eos/user/a/amlevin/data/nano/ww/2016/864423BD-4FA5-FE4E-A67E-44B94331731C.root",
"/eos/user/a/amlevin/data/nano/ww/2016/914978FA-D96B-E04C-9D19-3849F81B3EE0.root",
"/eos/user/a/amlevin/data/nano/ww/2016/93D0C189-5EE1-0140-9307-9F90661FFF90.root",
"/eos/user/a/amlevin/data/nano/ww/2016/9B8578C7-49D7-8844-B082-8172C1FB534D.root",
"/eos/user/a/amlevin/data/nano/ww/2016/A19D9628-DA0B-0B44-ACCC-5597E7F8DF9C.root",
"/eos/user/a/amlevin/data/nano/ww/2016/A1BCCAA8-14FA-0845-BAEE-031B0EDA1D19.root",
"/eos/user/a/amlevin/data/nano/ww/2016/A4640085-D831-AF4F-8661-70A108B00584.root",
"/eos/user/a/amlevin/data/nano/ww/2016/A85AC2BA-F1DD-FB4B-ADAB-A9C4C8AE35E6.root",
"/eos/user/a/amlevin/data/nano/ww/2016/C28F9BCD-5A01-D645-8759-EA075D16E17D.root",
"/eos/user/a/amlevin/data/nano/ww/2016/C8277FEB-771C-B34D-B638-335566291F6A.root",
"/eos/user/a/amlevin/data/nano/ww/2016/D5C526BF-95A2-F248-BBA6-82AAC42CD06B.root",
"/eos/user/a/amlevin/data/nano/ww/2016/D6C803BD-AF23-BB4F-8D52-A248F2186127.root",
"/eos/user/a/amlevin/data/nano/ww/2016/D86B5C0C-B2FA-9A45-8CFD-F6BDFE6260EB.root",
"/eos/user/a/amlevin/data/nano/ww/2016/DA100535-8EE3-7749-BBDC-9F75EBA8D7F5.root",
"/eos/user/a/amlevin/data/nano/ww/2016/E32B21CA-53DA-814F-9F2B-2149A0E76712.root",
"/eos/user/a/amlevin/data/nano/ww/2016/E4130C64-2E68-A649-80FC-2161759A5174.root",
"/eos/user/a/amlevin/data/nano/ww/2016/F44E06AC-57B3-2E42-9B5D-6E82C4443387.root",
"/eos/user/a/amlevin/data/nano/ww/2016/FAF7F5F5-8750-8347-A890-F6106D178BB6.root",
"/eos/user/a/amlevin/data/nano/ww/2016/FB4B351A-BD0F-CF4A-9FAF-100A2A3A86BC.root"],
"QCDWHJJ": [
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.100.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.104.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.105.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.10.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.11.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.15.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.16.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.18.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.20.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.21.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.22.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.23.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.24.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.25.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.26.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.27.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.28.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.29.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.2.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.31.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.32.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.33.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.34.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.35.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.37.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.38.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.39.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.40.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.42.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.43.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.44.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.45.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.46.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.50.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.51.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.56.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.57.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.58.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.5.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.60.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.61.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.63.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.64.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.66.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.6.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.70.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.71.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.73.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.74.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.76.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.79.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.81.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.83.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.84.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.86.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.87.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.88.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.89.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.94.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.95.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.96.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.97.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.98.root",
"/eos/user/a/amlevin/data/nano/qcdwhjj/2016/qcdwhjj.9.root"],
"W" : [
"/eos/user/a/amlevin/data/nano/wjets/2016/06D35860-716B-E84F-96EA-DBFD9A0AC14A.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/07C70462-2333-8B40-8977-02C350035338.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/0F75C859-C4F1-4249-9328-E7BBCE7860F5.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/10A236B0-E23B-F34A-AF83-6EFED2748771.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/12C7D7AD-ACFA-9C43-8A36-55BA69ADEF50.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/1AF787D9-106E-4E46-8C54-3F9090A708E0.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/1D79C4F2-8F87-9B4D-9660-0646A27F1025.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/1D8749B0-6514-F24F-A42C-369253FF8F4A.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/1FA1242B-66A9-D74A-9D22-AE0549FEFCED.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/23842545-1A1D-284C-8D46-46303BAD75F2.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/26E3E1F6-18BF-664C-84F1-FB18D69679DF.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/2B0DCBCE-8852-2B4E-A993-1398BF2E3392.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/3120BC9D-D7C6-1A43-8AFB-D3E8423B79BE.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/35ECA093-07B3-BA4A-A18C-75098C8EECFF.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/3C87B3D0-A038-5247-A39A-94AB8A443336.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/45E1165D-FBB0-6448-9157-2788BACCF9A8.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/483C6FA1-D83F-D349-BC48-FE99C71A3F96.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/4A3D4710-698B-424F-9EAE-C289222249E4.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/4E411982-8D43-394A-9A7B-D4769C8ACAE1.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/4F38E1D4-09B5-3442-B1B3-0A35D9E8E8BA.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/51BB7415-AC5A-0946-BBBA-5E88D0D5E029.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/5D0C1B02-06F5-E041-AC72-D26EF75DCBD7.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/6ACF2FB0-C11E-6C48-A1B9-642FAF1D5411.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/6E6D002B-1330-074B-B506-8281F9DA5434.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/70C436CE-7FB4-4D46-9070-651C53839065.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/768AFAF6-0083-AF4C-930B-7CA94884F860.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/78DE8C69-9824-5846-8920-3AD7322F19F8.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/7B19847D-A4A2-6E4E-BC7F-225047C861B9.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/7FA64119-B032-4840-B225-D6F842991206.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/809BF93D-EE1F-2148-882A-3DCC332BAFFA.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/8167B3FD-F6CC-3C44-9BF0-C1D64EA8183F.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/81CB2648-DEC4-4846-AA4C-7E37ADE90DC8.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/8511AE5A-6258-6E44-ACC1-58443EF920F2.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/98937CEC-4F07-FC4F-8048-BA0D5F50F8C2.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/9FF62D73-FA09-AD4B-B1DD-863E006F8927.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/A2E11814-5D2F-2E49-BA7F-7D28D6D90733.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/A8E108BE-AE45-4842-AEBD-302BABF69F8F.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/A8EE4133-2D84-464C-ABFA-C4C4918F4DC6.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/B34F7845-2C9B-7245-82A8-5E3A7785FA75.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/BBC37465-E2A8-0A44-92EE-B13FA9D781D7.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/BE05B6E6-B207-1946-9FC2-82E35E43464F.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/C5F0F6DA-3F77-2644-B362-7C461FB8DC3A.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/C6514BB9-3259-D047-BC77-1DC88093BD00.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/E68DD6C1-8FC3-CC4B-A831-C533F39ECF5B.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/E703D4B2-AA80-C74E-B51C-E0BD114D2DC8.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/E8E4CF17-4CAB-F843-8C7C-559F53821D3E.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/EAAC026C-178E-784E-B1E3-7009701F21B4.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/EFD50059-A086-714A-AFCD-CF3D374D75AA.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/F9334797-7706-DE49-A76D-88EF9A67176B.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/FA8E115E-70DF-A644-A0AA-404F70DB4A27.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/FBCFCB4F-62E3-4E42-924A-AE76897CBBD0.root",
"/eos/user/a/amlevin/data/nano/wjets/2016/FD4FBAD9-DA92-2648-B6A3-E0032410F13E.root"],
"SingleMuon" : [
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0107961B-4308-F845-8F96-E14622BBA484.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/01479010-3B51-B04A-9C5C-7800085D11B8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0264F565-9A1D-BA44-A1CF-E3EDAF5E134D.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0292AEDD-4149-7D4F-9F08-4B8D0C751BD0.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/02C01469-552A-7E45-9894-6406F0030DB3.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/038B9744-94C1-024F-B78F-E42B8AE52F50.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/03F57A33-8FF6-6F49-81CE-FAE32DC3A8DC.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/044B8A93-571B-0F46-BBC4-FFEDD253984D.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/05964001-92E0-8C4F-AF6C-9A24C11D1311.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/06831354-D977-364D-B76A-B0410B1C56E4.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/06AF8B5B-B7CC-6442-A94F-AD4D9A7306E7.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/070EC05A-1BD0-B940-A8D2-AD6650144A72.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/08211FAE-1506-CC41-9034-D0E7E4BE3DC1.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/08B84717-5E23-2044-BE33-840A2BD7ED41.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/08BEF135-37C0-F145-88E9-63F7A10D3BC7.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0978ED93-2E13-8B4B-869C-589FEDBE8BFB.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0993AF2B-E117-854B-B51F-FC3A4C3C8899.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0A30ED25-D2D2-2143-A4AD-250FFF2E9F57.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0A4230E2-0C75-604D-890F-A4CE5E5C164E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0AFA8A0B-AC2C-3941-9274-B691E0D0FC76.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0B54974C-D6CA-0748-9ED9-03BB3DDD82C8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0DC685D9-6630-1A4D-8508-7F3EBEA4CF46.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0DEE1709-0416-F24B-ACB2-C68997CB6465.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0F03FF10-7F20-8A4E-A19A-418C090A3F97.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0F6506E1-94CC-E14B-AE40-F72BD37BF8D4.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/0FE0804D-33EA-7C41-A968-132B08494EA6.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/10573E58-8061-F349-B708-275CD2609AF7.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/10887F9B-5FA6-C349-84A7-DAE5B1D6BAEF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/120B99EA-7A22-094A-889D-E9049DB1A21C.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/122D6DE9-04AE-D14B-8663-8447B5AD17BD.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/136621ED-C32B-B347-9FC1-1BE4AF4F8EF3.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/138679D1-FE57-EA44-B64F-99BDD543B901.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/13EF306B-E9EA-B24F-8094-DAC8A7BC4AE0.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/143F9660-3FE1-3E4E-A58D-FD56967BA9B8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/148777C2-F638-714C-AF76-00A32CD6F32A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/14EA6AFB-593D-074F-AF31-BF588987BE99.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/1534AD93-BAF5-014E-8BF5-5C1CAFF1077D.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/15A37A29-EEDC-D342-94BA-57BE290A5B7E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/15A7AA88-B1A1-384F-8AA5-13031C1A5E87.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/16CB6056-D4D5-2046-9D30-82D8985F831F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/16F5A02F-B9CD-504B-AB4F-E07DCD182818.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/171D273F-B79C-984A-8531-74411579001B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/172A58F3-41E6-7441-B0BB-AAEC6476F392.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/17F6E03B-6B69-4341-863B-538EED26B89B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/17FE90D3-5ABF-0D44-AF65-88388A6FF0E9.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/19C28AC5-2AC9-834E-B9C3-D03B7342EE54.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/1A3CC77C-70B2-0842-8536-8D66639E306C.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/1B8651A5-21DD-9F4D-BD58-4D20FB836DDF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/1BE80946-8275-7A46-842F-6DBDB9193FAD.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/1C08614E-0C0E-6044-966A-CAF630CAEF8F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/1D280BA6-2EBB-A544-91EB-E821F3CF0B06.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/1D772C99-43DE-2749-9EFF-0A5DFFC5EFD8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/1D87B4FB-E31C-9F43-AC21-C32469DE9FC6.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/1E594251-5C49-3649-9413-14663553AE7D.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/1E728DC9-D98B-A94C-8556-218835F47DAD.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/1EA71AED-29E5-7F45-8646-BFFDB586A5E9.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/1EB443F2-1230-8042-B8AE-FD50329CA59B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2045F967-9F0A-7C46-9946-787B27D56E88.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/20ACE5E9-4359-FB44-9E8A-56CBB5742839.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/20CCC9D4-0AD5-9E42-B642-DE86F25153ED.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2139E3EF-ACE7-A841-8FB5-5B4F2E8848F6.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/235310C5-0AE2-5447-AA04-49442DDB4005.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/236A04EE-C105-D947-8A2E-F8CC6731644F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/23734531-415C-BD43-83DE-94E787D3F57C.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/23953892-6295-B341-9D0B-6892EBDD67FC.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2417FC8B-0363-4D48-9843-E42B14C8FF53.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/247E64DA-FE25-DB4E-ACFE-034455CFD3EB.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2489751B-A7BA-644B-8FE7-BFF0073158B0.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/25C2350E-9C24-D845-835C-30393D17D35A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/25DEFBD5-45E9-544B-9E88-B1D3ED80EFFA.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/261B08A0-22B5-5149-B00E-F9A774477005.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/26657EBD-54DA-AE48-922E-FE8C516349BB.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/266B2CD6-17D4-C04F-B1FD-445C371464E4.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/27F10BCC-A23F-5B47-9800-DF87D00DD0C8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/28415B3F-DFB1-5542-9217-7D134447AC48.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2A6FBB1A-7918-C548-BE13-F3BA8A463323.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2AABE6AB-784E-1C40-AE82-634DCFBC2ECA.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2AF26192-4B3A-324D-B14C-CFFFB1F69F4C.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2B51593C-FFCB-034F-9127-15A9BFCCD882.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2C22A425-8A80-704D-8800-F654BE0D2EA7.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2C7EC1F2-BE79-AB40-A7B1-58740E6440B4.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2CA777F4-83CF-0748-BC6A-E05FD63DE257.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2CBD5330-87CB-EF4F-B3BF-7B132A90B72A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2D385C50-1179-624D-B03D-E3F26649355C.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2D400D95-5C56-5843-B127-AEFD1CFBF08A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/2EE0F7C4-0464-A848-8B5E-E8902FE75AB2.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3006FC71-3B8B-E444-BA65-5B9B5A624427.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/30799F04-9B6A-9E44-9F8F-BBBEC2334324.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/31B08676-9EBB-A94F-9000-825AB7A15CB6.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/31BF084F-D414-1542-B61D-0F58520A17E8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/321A93D2-8910-2E47-8CBB-1751B1E90B5F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/34393CF3-C0AB-DC42-B730-5A21648842BD.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/34709CBD-22FA-A344-9DA4-336CB6FA3F21.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/35B1E60C-FD5E-1F47-953A-9256C1488490.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/36804A9C-3C31-E345-9E96-37B6E764E986.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/370BE877-DA24-DB41-A875-07A86EAB6852.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/37B993E8-7A11-6141-B42F-C20F3F6BB2D0.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/388AB3E1-8708-7D42-91BA-83E52373E808.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/393744CB-F2F6-554E-B61E-0E33DB7A2436.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/39B0B7DF-2C71-9E4E-952F-833475062880.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/39C22886-986F-934C-A8F6-C1E0E11B5F2B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3A5970C2-4B39-8147-B5F5-D13256BA4A0B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3A8A8FD7-1FB0-FC49-93FD-20378419A49A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3ADD6FD9-5B2A-1441-BC2E-D68FDA486B4B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3B6DFB90-29EB-0641-A32D-01A8EB3ECCE3.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3BCDC9CF-57F5-104F-8C20-2BA300439111.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3BEA5C1B-C344-D74F-937C-D6FB3D241731.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3C47E0A9-A390-284F-B560-9C5737BA8254.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3D66ACA7-CB10-E940-8FB0-3E0478BD4923.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3DF19DD0-6B23-0548-94AB-A20C605C1C5B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3E76909B-34D1-DE46-8355-51C2ADD587D5.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3F457462-2AEB-E247-9A27-5AEF9384A933.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/3FCF8DC3-CFC0-3543-A351-CF663AE50BBF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/40EA32B0-1284-7246-83A3-A703D908F9FC.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/4175B42A-D042-F246-9F10-24B7E7A22ABA.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/41A65DB5-8F3D-E049-A007-12A007B89033.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/43C4322E-4EF0-F046-AEEA-A38587E36062.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/44AC3FD4-94DD-FB43-A7C7-C102155521AF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/46E26356-FF36-2941-B03A-972CA50D266E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/473F6C9A-F7AF-2640-886D-D814F98886D8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/489586A9-8DDB-174C-B88E-E6C82CED947F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/49036A83-AC2C-5D4D-83DB-DECD692E1708.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/4989DEE1-8D09-824C-AC86-2E191393EBCF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/49E2BB63-31FC-7247-A176-B67D0C0F1924.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/4A97CBF1-A25D-294D-99DB-D874853F9C64.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/4AD09F91-D53E-E345-817A-BBC670875478.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/4AF3D049-8CA4-6D4C-AB8D-F8DC97553F8A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/4B985344-7CD7-6142-BCB6-B95E47279B9A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/4C688018-90E1-5847-AEA3-D5F5DEB330EA.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/4CE0AEF4-0128-1947-9CBB-6459DF4138AE.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/50ADF677-F507-D044-BCF2-40392F4A89DB.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/51E02A2D-02CE-F845-BC7E-A9AB5D3FB0A8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/51E1A1D8-FBF4-0E45-9F37-FBFD4B5A7A09.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/53945FFD-11B3-4240-A82C-976090AB37E5.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/53A51E60-A5C9-1E45-AE99-CE8D7A539E09.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/54C96DFE-FE30-E649-B71C-58C35FB9A273.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5553FFE8-D1DF-9842-9F49-16BB8EFFD4A3.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/55B39F65-9230-1D4B-80B9-C74894153453.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/55BB0652-F28A-8C44-8EF5-C8DC05A1FCE8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5709DEB4-AAE5-7949-975C-119CF30ABA70.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/576759DA-4A35-534B-B926-2A9E4A5A7268.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/579013D1-626D-D943-A8B0-A1A558C54F33.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/594F3D70-ED17-3D42-A322-269D0A074F9A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/59905258-7CCF-BB41-AF01-D9FCBBECC03D.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/59EC82CD-D050-5F46-88F7-526776849F77.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5B0F93E6-6976-D445-8700-255064AC7C78.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5B25D22F-F2CB-884E-931A-420B7E75E7B3.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5B9B9627-5FB6-9941-8638-51E5E18A07D4.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5BF83F4E-336D-1641-9576-53FC198765B5.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5CA1752A-6DE7-6042-B320-6AF27E4E1557.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5CA4CE73-629C-2F48-939C-4274B369F112.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5CE8C621-31EB-1749-A93E-426EB6863338.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5E060727-1BCF-1F40-BBD8-D6841127F7BA.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5EC1EF0A-B8A8-5F42-920B-AE253638717B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5F5362B6-CEC3-9641-A110-9191AB6720F1.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/5F5E0726-4074-DA46-9AD9-C39389E885EC.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/60E957CB-2729-5D4C-B261-43406A639E8F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/61FC1E38-F75C-6B44-AD19-A9894155874E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/626C68F8-9DBD-194B-993E-8672B8E73BAE.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/6322429F-B620-854D-B4AC-6892359CAB2C.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/6351C999-8EF7-6345-822D-E2EF1C203603.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/6539EE36-6427-5249-8178-20DB40114D75.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/65EE57DC-3578-6246-BD7C-81F6B6EF62B4.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/66A58A5A-4A7F-0546-B370-BABDDE02ADFD.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/66A7A94B-7B8B-154F-AB44-8F6552EDBF42.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/6789206B-7C28-5C44-B975-A5EF9D7BC6A0.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/680D1B84-37F9-E547-9244-E5833F3395D6.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/6C8BBA90-79DC-CA45-9E0F-4AB683E91D8F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/6CBB587F-1366-6F40-973C-7F98EA9C8335.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/6F0D8ECB-94C4-AD48-878A-519E5B4BA158.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/6FFED094-8391-E342-AA75-5B7F8379978D.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/7063C5EB-0B17-C747-BEFF-E49AA9230484.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/710A09CB-D2D9-0948-B4EC-6D94FD9F612B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/715BF52A-7FE3-D848-A45F-99A0DDA80205.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/716B3E8D-FE51-C44E-A522-185F9614E18C.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/71A0EA85-9A1A-984C-9736-EA2B289815F5.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/73EFC7F5-FE87-B241-AF88-A0B4F2557AFD.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/755A1FC4-28F1-5E49-8466-8C3B34733F1B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/773694FA-CCEA-A94E-B460-065761E6DC20.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/77C51055-617E-FF4E-A832-BD87CC8DC3D4.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/788CA5CB-13E0-544A-8427-33E4591CF20B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/78995FCF-4A0B-364B-8875-9207B74685CF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/791500CE-B4AB-7C46-91D7-4E689B94F55A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/794286F9-A4E6-2449-9729-DFD9F44DEF2D.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/795E7335-9CF7-5F45-800E-E61BA4F46B54.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/79F5CA2C-DD37-C742-A19B-8EAEEB531072.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/7A2C148A-EDA6-6341-85BC-936039FE6D6B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/7A9FF3A3-9B82-204A-8FCE-1CB9ECFC65FB.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/7BA44B30-9EBA-6E4A-8DD6-2A239BDDBA50.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/7BCD0C03-8572-9B40-84A6-527976F62CE9.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/7C32E1BE-A26A-9B48-8D42-3110D13440DA.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/7C47FA39-A93D-BB47-A2D5-902701077DC6.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/7D65C758-5003-DB4C-95AF-5B324BE45358.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/7DD9EB04-8421-C343-B1FC-6D3533566A6E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/7DE32971-D86B-704A-B2D5-4CE16229F84F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/7EC9A11A-ACD9-4A46-8770-B7A48670CBBF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/7F5791FE-262B-EA45-BE19-04E4F21FEE1E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/809D6534-33BF-244B-B689-8C63F763959A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/81CAFEE3-F676-9848-91F1-A698B854EE9B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/82A53FB2-CE9D-744A-AD5B-E05213C353A5.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/84E83B69-73A2-BF4C-A59C-C3C10F71C76E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/8531DBAE-434E-FF45-9647-1FFA5B964660.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/874B6634-10CD-4541-97E6-7A97BA44F4C6.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/88C9D533-7FC0-1A4A-A43E-ECF3270D2200.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/89029286-01FF-E643-A9B3-D35C09F54C1B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/8A517AC7-989B-1A41-A0AB-3EF7E5DFEDE0.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/8E33310F-1983-EE41-A97A-7D6853AF1BAD.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/8F131E28-A74E-9D45-B799-2A8CFBC3E45D.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/8F3012A6-72CA-6C4F-81AA-64EEB640EA84.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/8F48F391-94FB-524F-82B8-D03C17C3BBBA.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/8F669943-6CC9-7248-9ECA-858EC1ACD9DE.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/8F897B77-7EDA-B54B-812C-2185345CE97F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/8FB46D51-A84B-8043-A86B-25743CE20170.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/9110E45A-9003-CB4F-9220-18B2EBFA2B42.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/912D8BEA-4A7F-1243-B1F9-79C53F8602CB.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/92238EFF-D3BC-B745-B707-BB7096790F83.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/93E71464-3924-9548-BFEF-C590B902DAD8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/94CF3C61-CC22-FA43-BDB6-34609AD9B155.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/95AF7A2E-F973-5642-96C6-B37EA7C96FF7.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/961E3D28-0C03-C944-B0AD-5C9BE7DF8419.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/9713F817-B8FF-C048-94FF-9F45A5750A1E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/977A32B2-1672-0646-8029-C12B37014B64.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/97A94352-7E24-FE4C-83A8-071AB6D56AE6.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/99824E36-CD3E-0445-85BF-5B18A4C4126C.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/9C99B9F7-7F24-FE4A-A28D-E31D70DBDF66.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/9D4060EE-88CD-1942-AECB-A85E7634F5A1.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/9EB5A008-3D02-6A45-817F-71888DFFA070.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/9F77D7D1-C075-C44F-B2B4-FAFEC0126484.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/9FD1BF22-3392-004F-A9AE-4FA8A02AD545.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A11CF168-BAF6-C04E-92FC-67ECFE105A24.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A1DF7082-FF8C-114D-BA6C-FC0EDFDF61B2.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A20283AC-8E40-104E-A56F-0450F7FFF3C8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A283D96B-1BCF-7042-A31C-79200BA2DDEC.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A319F9B7-FB48-834A-9F56-BBBDBC4EF77A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A3279183-BDA6-F44C-8EF3-C063BC598B68.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A35F36C7-8A32-094D-BE5A-022DF8A7E031.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A3664F2C-DE66-C64E-8471-88F341DBC4BF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A3B52E09-7E69-1646-B3BC-0192BB4559E7.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A40B199F-D758-5640-BDFB-DD724205E5A9.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A41FFD8B-188C-5540-8BA7-79F5B9A42FFA.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A4428D33-117D-5B4E-9641-4BED58AD15F3.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A468306F-D7AD-004D-ABE1-7BF9E8647E48.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A5B416ED-2D1C-4846-92FF-CA7E5C6DBA61.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A605975B-76C4-F54B-9DD3-15A14FB2D850.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A6585C88-B184-6440-AD3F-FCF1EF930B4C.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A6CA373A-583D-BD42-9556-0351CF148FDC.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A7256324-3DFC-2243-9FB6-FBB68C93E24A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A7ADD16A-5EA4-BD40-AB01-9A3F9E8BB161.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A8CDFFDE-BE32-8741-A0EF-E2B1A983CD4D.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A93A0756-9C67-864C-AF43-C48069B6B7C3.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A97F82B8-F32C-7C46-8128-7D6578798F0F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/A9BD13D1-CFE2-CF46-A483-2557A80E06F7.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/AAC7B364-12FC-9147-926E-91B98AAB39FA.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/ABAF299E-EAF9-5043-B89B-1D13A9475BBC.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/ABC5FF26-39EB-5543-835F-4C1B5981D052.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/ABF5CD64-0947-964E-8104-84D880B737CD.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/ABF89EFE-CE22-184B-979B-A53A75039EFB.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/ACF68248-100D-D345-80C9-AECE7BCB73AF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/AD99EB27-9D69-B74E-981C-F1869AED6582.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B02F3578-EDAA-D649-BF16-B0090E99B34A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B0438805-AA68-EB43-AE71-4E68A36554CB.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B0858743-46D3-7C41-BF23-0EDBF408FC25.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B0CFD000-E84C-344D-BED7-8AB35F3DBA36.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B14C8B3C-8299-3245-B2F6-06AB1352CA4F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B15FFF7B-30C2-2449-915B-A172881A22B2.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B16FC28D-21F4-7640-8E4C-7097B3CB591D.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B31B5932-416C-C04A-9CBE-F419FB2C35D3.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B4E0C1C3-E86E-3F4B-9B1A-D17E3ADB0908.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B5840BE5-911B-9542-8583-E6C51AB2A6E6.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B59B700B-BF01-7B42-BC16-292C09285BC2.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B61CB2F3-3F7B-C143-B349-C9FE421103AC.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B7EEE475-047D-CA4D-B680-617FE16540AC.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B875B8F0-1932-E34A-8E54-76EC840A9A8D.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/B8DE09E5-093E-2343-AE22-AFB36AAA6D5C.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/BC07CDBA-E54A-AD4A-BC24-960043783955.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/BC6A0B90-68D1-F842-B46D-7ED1A4597C81.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/BCBAEA3A-2460-A34F-9CD8-A8ED2680A645.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/BCE9916A-1796-5843-A0C4-BE23ADF386DE.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/BD979CEA-AA60-104B-B656-A2A7ACE54A4B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/BE433087-4E55-6F44-8EBF-EDD13E745BAE.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/BE5B1EAA-9245-5346-9311-14DDF222C580.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/BF01A0D3-6171-E64D-A28E-48C587B0942E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/BF3B627A-F37C-3F4F-97EA-153B4B08E240.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/C08430A5-EB2C-EB4D-A653-7CA621B11BFC.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/C10D77FF-118D-3C45-A997-B586C4324D31.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/C1175B23-982A-EF46-8420-4BBA32D41163.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/C144B7AE-E8C0-9D40-BECF-D027FA604946.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/C250DA77-36D8-2440-8F71-4A913F18EA6C.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/C2A7D179-B263-704C-A8EB-B75F4D82F7B5.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/C40154F9-3E74-D34A-BC9E-0F7B5655468F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/C4CBCCB7-BF30-9847-9136-309EC3E7311D.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/C653EBCD-C75E-8441-BB1D-12A90F76D597.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/C715AA35-8956-354D-B32F-7396713FA3F3.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/C760C7A6-E4AB-E447-93BB-1EAD0177DB41.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/C8653579-E6F0-5D44-ACEA-154A482DF4AD.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/CA9164D3-3CB4-F84C-98B3-C93E46450BF3.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/CBFBEA23-7C46-624E-98B1-2CAEE35458C8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/CC4EEB3C-0B58-9445-BC41-D00F916E17C5.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/CC55E1DA-ADBB-0342-A53C-C7B7CFA6AD85.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/CD1172F5-A00E-0041-9A9C-5E1862026357.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/CDB1F701-F889-854D-B8DC-7C7F6378AA66.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/CE1F1724-2933-C44F-829E-0FA4CB49D49A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/CEB3C22E-30A5-704F-BD07-BD6D0838DAFA.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/CEE8B116-DBE4-E54F-B038-A81090B2A147.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/CF7476C7-BFC1-1F4C-BB99-CF90CB99A6B4.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D0387C2C-B889-DE4E-B69B-6C85C9DF5DBE.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D05C985B-A5E8-8E4E-9A65-C756D09C5723.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D115B24D-D04F-4E42-B5C0-4DDF94B4EC51.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D2392663-56C1-5849-A224-28CD5251166A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D2C497A2-F632-AC4B-8655-1ACA9F2A43DE.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D2D73C3D-3578-9445-B131-F0CE82C9C278.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D2F76766-1F43-264A-A39F-2DE68BF68752.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D31162AA-AC66-3F4F-97D5-1DDE64E93D99.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D3E1C7D0-4615-C249-92C7-103153B85522.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D548310D-36A5-DE46-97C5-15933EAF7DC9.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D5642602-6BA0-E843-B8CC-03B0072E1796.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D56ED17F-D34E-A14C-BA80-F7A28A294B87.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D65EB6D0-8B85-154A-BC19-59EBF79E98E8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D86602E7-FD64-F144-813F-0BC7E92E4801.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/D87A4CC6-6CF7-6943-838E-475A7533EAA9.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/DAF58E02-6236-A143-9D45-A8A3B15355FA.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/DBA71BF9-1B93-EC42-A96C-BD2B35E3538E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/DCA6F24D-0006-634B-9964-65CE4B8E70D5.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/DCC13781-ED79-3D48-A43A-91553C70581E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/DD3399EB-8498-E848-BB6A-D0BC92D966AA.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/DE630F0B-48FC-AC49-9BD8-6951A585BA50.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/E037AB14-FEA9-4D40-AA30-9ADF56BD6270.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/E32CD72B-0233-E547-9C3C-4786ADBFDAB1.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/E52E3913-C4BD-714A-9A24-793D29214038.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/E5645EEB-122B-734F-9F04-704C91CE8586.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/E5B3AF0A-A0F0-594B-9E72-B9EA1859D549.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/E5CC21BE-26FC-334B-971E-B6B390EBB693.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/E6C5C5AA-946B-A34C-ABB5-4B14AA0E9ED7.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/E801FC2F-0877-9F4C-978D-FE9C94B92CAF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/E8E0A6FB-E2D7-BE40-B098-7FEABFB340FF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/E93D6A9B-27A7-A64A-9CC4-D0B1E190BE0A.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/EA9DD56A-11DE-AE48-8E43-07852D3E3694.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/EAABF77C-A79B-FB4E-B307-7B2F9CFC6C9E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/EB4F3DFA-05D3-D04E-BFA7-A241D08EC34E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/EC12C4E4-4F1F-4441-86B8-5F15678EDC79.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/EC4283CE-BB96-B845-8012-6ED0CC0FA7A2.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/EC4D9A98-59D5-DC42-ADA2-B63FC447E2A9.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/EE1BAEDB-88FC-674E-A969-707818596CA4.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/EE7627B5-D250-D740-9898-ECCB3B1E59C8.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/EF40339D-452F-DD48-80B2-B8BE511913E2.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F032F4EF-B33C-5F4C-957F-15D6F2ADCAF0.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F24742D6-8BF8-F640-AD2F-063EE84A759F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F250A551-A8DC-BC47-895C-3BE19459F5DC.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F302A865-17B0-064B-8154-41526BB38244.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F3289382-F0E6-5B43-83A8-B7D08209FFAF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F361D15A-F30D-4D45-8874-446B7450D4BF.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F3F7C96C-CF45-EB45-918C-9E5A05F8D899.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F3FD0787-A332-1B45-951F-4CD540A72ACB.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F4101B00-7984-C54F-98E9-2EF72AD66F43.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F4BB5330-5236-5940-BA0B-F13088D94F5E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F5ED054F-7091-7243-B66B-47770BCFBA97.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F7349EE3-AB18-774B-9B14-9FBBA2350E0B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F748EA94-609F-0C40-8D34-6F686FF6C6AD.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F773638E-CE20-944A-A120-755C3E071B2F.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F794187B-BAF4-C840-A9C9-93709CD6087E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/F9ECEECC-1226-DC4C-BD06-0C218B250D60.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/FAB9736A-50B6-9C4A-ABC8-09738672BA66.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/FBD685D1-F6C9-8946-BA61-14F8DB0B1B82.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/FC3AB504-0AB2-3145-A883-5BB67213097E.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/FD262A89-0289-F941-9494-5FD7556833B4.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/FE81EE6E-0E21-1847-9105-D700F2E762DB.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/FF8230B6-5C45-0F43-B910-7C2912F2A74B.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/FF9B5EE5-4C8E-3F41-B608-AC8DD47E51F2.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/FFD6CCA5-95F8-BE42-B747-5DF0E011D4AB.root",
"/eos/user/a/aera/andrew/data/nano/singlemuon/2016/FFED1DB7-3626-A04C-8438-BF055BD16672.root"],
"TTBarSemi" : [
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/00AAB546-6769-3345-A154-8B41C01307B5.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/00F8ACF6-2C09-A144-BA86-AFE99B2AEF33.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/04DB6487-69D2-8148-989A-2F184E3432B3.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/057ECF7E-3304-CE4C-9570-35E9253768F2.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/05C36F92-AF8E-0843-B76C-A8028F37F00B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/064FBC40-1FC1-9C42-9575-838F41D6242F.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/07FB0797-A762-724F-9ECF-F3D44B36C2DA.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/08D594BF-1983-8945-9ECA-8DC836CBFBDF.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/09EA2721-660A-0F49-80CC-E8AB8E038BFA.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/0A18E677-F3C5-1446-9C46-B499FE4B7A7B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/0AC371A4-D247-7349-A4DD-C7FE00A48D2D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/0CFE4E5D-6128-7047-A021-229946BE1578.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/0D520B44-5A58-9746-A85A-3B9A3A9126DB.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/0D943C9B-8B3C-A54B-9354-77967309063C.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/0EE97716-16BB-724B-904E-C92091989D26.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/0EFF6AA7-78FA-2448-AB33-AD3602D05EA3.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/0FCDE449-0BBF-6047-946F-3D6EEC9F771B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/1066DEF0-5770-A446-9BBF-621EA19C5D87.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/116DFE58-E199-514F-AC2A-1A322430AA17.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/1183F600-D74B-444B-A2BD-1F17E290CF33.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/11AD1CCC-8E67-E343-A014-77E672050C51.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/12D382E7-232B-CD4D-888B-78784C3F3292.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/1323D181-40E7-5E40-A8EF-2AB13075DEBD.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/13B83EBD-9741-924A-90EC-EB9EA771DCDD.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/14D103D0-4DB4-3F46-B37B-63B25808E6B6.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/14E0325E-24B0-0046-95E1-61C3770D04E4.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/14FCFAA2-C801-6548-8964-1FF6C4A54965.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/1763EF73-6904-074D-AE07-5758A5B953FF.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/1803C9C4-AC2E-364A-844F-4621039F1394.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/18CDD58B-CAF6-FD41-805C-C7E1CF34A0BA.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/1A951BEB-6A18-3D48-8A07-D5DE70FC5F88.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/1B77AC8C-4212-894B-A376-8E1BBDD8F3CE.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/1C818548-FBC8-D744-94BC-B69CB67B9490.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/1F70F02D-36D1-BD4A-ACBE-C797B81405F9.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/1F770C42-DEB3-ED44-8C2D-EDDCD2CF6D68.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/208910E3-3AE6-4041-9065-E4E6FF80AF35.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/218A36E7-5615-0540-8FA0-2B0B25CC34E1.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/23C4890A-CB7E-2D4A-B6F9-EEFE0190E6D9.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/244BF8B9-4430-4942-A39E-F2D62218239A.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/245CC65B-624C-1544-8B51-CDEAF235682B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/24C56E50-6CFF-5B4F-BFFF-2DD9CC389C6A.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/266C29DC-E1AF-1948-BFFD-F4E34F97EFED.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/27F3AA6C-9A57-3246-AF7B-F19045EA6189.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/29B28726-B537-CE47-944A-94BD682A1CA3.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/2BC0E066-0A59-8140-8302-F553B5AF40F9.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/2BF76729-6867-E844-ACD9-8836BA757C80.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/2BFA2CFA-E97B-EE45-BE6B-BB1F60D7468A.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/2C91BDC8-EBDC-824A-A358-9CD5CCB2A264.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/2D725446-D07E-BE4A-8DB3-E1E7AE60076E.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/2DA0D5B6-0BB7-0D4E-9888-90B5BDF747DA.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/2EB00EE4-6BEA-664E-B679-B208BB3BA100.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/2EB5DBB5-BDEF-F949-BA74-88E392C966AF.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/2EC0708C-7684-B144-953F-4E37D7482DB7.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/2F9E6E5B-B73D-0B45-B921-43AC625D0DE1.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/30F6D697-D63D-1843-971E-01075D9D46D2.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/32557AA4-C71F-614B-931B-9A61C52BC190.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/3315DF43-6EAF-F147-8B3D-BCBADFDEAF94.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/3432CF59-5280-054C-82EC-01132DD8E196.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/365B2653-D634-1C4E-B4B3-6F6A12DC2D9E.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/373A5DA0-03CE-514C-88B3-621F763C069A.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/37652545-830A-6C41-A7B6-1A511510E28D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/395ACBEA-1F5F-4E45-9C1A-8A3FE82A154D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/3A50DB1B-526E-F84B-B423-C25886D8B9E2.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/3BA80D09-FE89-5C4B-92B3-DB59E085155F.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/3BB7F920-2C92-1E4C-AF34-39ED0D0A42A3.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/3D5010B0-0238-984D-946C-6E89DA586AA9.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/3E41D18F-3814-FC48-8DAE-2B6B7F998923.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/3F23AC07-5316-DB41-9D55-C79A3F7BFF6F.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/40352B25-8DCF-A44D-9E67-D6F8C98F1B5F.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/40FFB613-6499-6D4B-9AEB-FF81DA6F4987.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/436E39C2-4D96-2B48-B7A8-DA2AEBD7F83A.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/439526D3-0EB1-3346-9658-8CE3693240A7.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/440B77E9-4741-4C43-A46B-9FD1BF0D09AE.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/44357630-287B-7D4A-B7B5-710ECF2CBBE6.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/45F2F0C0-8DA3-494D-B4F9-F0CA15746068.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/471A9040-AB3A-D744-9CF6-4D17C6C36BAD.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/47C164F8-349A-924B-A94C-A49AC4C3667E.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/4A14D069-CDAE-8645-BE0A-2B4835A9C309.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/4B62D5B3-24C1-C941-AF5C-82FDF5F0FAED.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/4BBBFF2C-6EA3-8141-A267-A6AEB7E8FA9E.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/4C139123-3067-DB41-BECB-AD8D2E49B38B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/4CA26608-6F0B-F044-AD29-EF3EB5B14D28.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/4CC4FBB8-0DAD-1F49-903A-DA5DCB865BFF.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/4CDA4B9B-1890-224C-ADFA-323781C3A716.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/4E9F9B08-B8E3-6349-B56D-4FFB02B72537.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/4F75A76B-811A-6B4C-B882-685AC409AF1A.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/4FA668B7-E166-584E-A0C8-0A1F36E4BD09.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/508E87E6-01CD-9F46-9921-DAE7F0CFAF2B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/51452850-C2FF-0640-815B-4F11464A5750.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/51D99FAD-7411-3845-84C8-C1737E8E4234.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/5244DF32-7C19-AD46-842B-1B95066C88BC.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/529F2EDE-2C63-A846-A6F0-7804D4E4AFF9.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/531BC3DE-480B-754A-8276-927FC1675136.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/533019FB-B6BC-844A-8765-CD9353CCAFAA.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/5361ABE6-5C42-044E-9B06-6C37EF0B2214.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/56A4989D-EFC6-E34B-ADE1-DDE4B1001C78.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/5716CD88-B306-7D4F-8F70-4FDCC6133E9B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/599BCD72-A0A5-FB41-9AA2-73A156C3236E.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/5C750601-4EDC-6F46-9B2A-D3B60F4BE9DF.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/5D38355A-B1E5-1445-B3E5-5BC0742CD211.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/5D5EF4CE-4D7E-574A-AF90-7400766B6CB9.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/5E3D3778-8D5F-2442-913A-3EA219284611.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/5EF886FD-62B4-7249-B686-7E0C8DD88515.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/5F0F6CBB-9724-484A-BCC4-B06E3EF242BC.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/5F5FBB5C-170A-5E4E-A370-4EA224E273F3.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/624BF3E8-CE3F-AF43-8A09-49C42D9C420F.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/6546823D-A1FF-ED4F-BD37-17BA3BB84942.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/65C5763B-875B-2944-BFCF-439C041083DA.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/66EAAD64-051E-C145-B2B7-A0F4514C171D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/671C3A1D-3BF1-D041-8CE7-2C57F24C3CAC.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/67CECA18-E87E-524C-B356-EC9151A056AD.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/67EA9919-C62C-AA4A-8203-BCE945537787.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/6858170A-7730-6443-8F89-5366592DEA24.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/68EA23B4-F40F-CC43-8001-6DF063B788F5.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/693D5AE7-DA09-A041-91D0-F00A9635EC17.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/6C52FC70-C68A-F349-B2E3-970ECB7B6554.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/6D7A9A97-E4D0-0D4C-9EE1-3745E0F010DF.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/6E243349-628A-194C-B995-59FE18A614EF.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/6EACF54C-673C-DB44-83C6-4A7FB9D5F9F5.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/6F7529A8-DE38-984A-94C8-BC5962D538F0.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/7080FCCB-134F-3D48-AA0F-4859884DB679.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/709C01EF-E138-FF41-9AA7-C419B2BC7287.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/710AB23F-A2EE-214A-9B8C-6FE1215D6E78.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/7135B247-03F1-D949-AC58-8361D0627DF5.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/74643869-4F32-4741-A961-DA3D998E5F6D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/74FADE0E-E80D-CA41-80BB-D155D0A43ED6.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/75129056-24AB-3046-A7AE-1441851052DC.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/7827F181-DECD-144A-A242-BB68EFB57EDC.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/78C383C6-B23D-014E-8FD3-6A724DCB8BBD.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/79ECF990-A30B-4146-A3DE-BC0245233C18.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/7A0CFDFF-B350-0F49-9A2F-4BD536FDDF9B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/7B00BBFA-4E9E-D143-87E1-3398DD41635D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/7C647ED4-F8ED-B64B-B951-DB3D5C664B9A.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/7CD7649C-2BDB-6447-B4CE-6E09285B4C6B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/7EDDC59F-2427-9A44-8B14-4F840971A47A.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/7EFD6B73-443E-594C-A7C4-0A7D74C5F2FB.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/7F4B43E0-BE1B-0448-85CD-BA5DAE32000A.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/7FA74DCA-F28B-9745-8C42-DEC40A5C1EEF.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/816F182E-1B87-D642-B992-9FE579792702.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/81E86B24-56F5-D340-9745-C0E3E8BFBADD.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/8331A3A4-7E76-EA43-A9B6-260A2BB2A75A.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/8465E4CF-3304-A94B-940F-964D0E4E85D8.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/854F1E7A-AF22-7142-95E3-CCAF68083579.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/861F48B6-DCE9-0D47-8A67-0FB0CCCF4337.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/86424EF5-BDA1-0C46-B926-7A0501D4E650.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/8670A0B6-E7B4-6F44-B9D2-5F25C6E59F65.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/87B4ABB4-CFF0-E443-896E-2CEDA9297C9E.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/87CE44C7-CD54-9E43-97BB-ED5F67410B92.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/888AF839-450F-9A43-8F89-EB9BB04BB76B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/88A75C2C-2538-A24A-B917-3968FE8AE7FA.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/88DC1030-AB84-4B43-A0F8-4705431B2B6C.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/891C8ECF-BF5C-DD4B-85EB-72B2E946293B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/8924BA72-6D4E-B84F-9531-AD65D9F80636.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/89A0B075-AEC1-F54D-AA8C-2DF61D0DAD86.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/8B00731A-0BAE-4641-9A56-45D2D3F5F8F2.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/8BAF562E-0AA2-F74E-B390-14EDA2AE4902.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/8C0F25CB-E2B2-344E-BAB3-2A3DD5F402AD.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/8E45A7D5-55E4-4F4C-B771-BCF51D8CF0EB.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/90BDD5B7-E753-AF40-8014-9EF9FF30566D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/91AAF311-DDB4-C44A-9C55-2A7E8A05E0ED.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/9255B0C9-7BE4-304C-9EBC-5416E9C6196E.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/92864F41-14D8-9249-A2BE-AEF5691EC01B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/9483B8FC-2845-7A4B-93FE-2DB403538EE1.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/94E16651-ABFF-D441-8D33-46C6C2CED96D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/94E90D21-998B-3941-A0C5-51171E8618A5.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/94F81334-3677-5F48-971A-CC47DB478710.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/963D970C-3EE8-0E49-88CE-0FFC09312BD0.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/97E91EDC-6320-8544-80BB-83D9C4FCD812.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/98F148C6-FD40-6A4C-A704-5E3398A792C6.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/99B0F95F-022D-514A-8149-97711021E9D5.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/99D53430-ADA9-0D4B-B305-354BE04AEDA4.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/9BF51CE0-A75F-F048-A413-12B2565172D2.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/9CA16A2D-B4C7-8C4D-BF4D-E5A2197060DF.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/9E0D145D-7926-614C-B737-876DA6320927.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/9F33EEFE-EBD3-CA47-BED9-0FB415BC0826.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/9F637E6E-2137-B242-9A4B-E78ED895A733.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/A010C6EB-0BD7-6146-9A2E-A20D85B6C701.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/A155F415-1738-0240-871C-739178C64D71.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/A2676F32-7167-0A4A-B9FB-604938751AAD.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/A5FBF355-99B0-034A-801D-AA5A4A1FAA8C.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/A62F70D8-3E4F-D14A-9AA7-4D51D17D4502.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/AA726D80-DB34-BB46-9354-A895B05F1E7F.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/AACA0570-DAA3-C845-AE2C-D3F6FD69D041.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/AC929C47-AB76-604D-8741-1D33D3B64284.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/AD2D16EE-B5E0-1B4B-8768-8C204D6284CB.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/AD2D6458-8688-A641-9665-159E1BA85B6B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/B06A8CD4-BAA1-834F-B1B0-BBEF5F928371.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/B09BD29A-A667-5C4C-94CA-126E497E5089.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/B2840F6E-B45F-C544-88C5-1B86E45E0AB4.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/B29E2A00-8B60-514D-B1BE-8EF19E2FB9B0.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/B3F43573-3175-914F-8FE6-72819D5036B2.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/B425EB28-90CF-2A48-B262-35C4D2D40F52.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/B5090F68-DE36-E345-97E4-C48233BBCECC.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/B527F684-9AFD-7E4F-9C83-C945C306CCC5.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/B63C52EB-65EA-8E43-A146-90E10ED5A541.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/B7E1EC36-AE45-5747-B0F2-12F30E240C15.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/B99D4D47-37DA-C749-8E26-E6C277C4EF10.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/B9D0EAE9-52B4-B64B-BA46-0C710CA05EE6.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/BA2B9E76-28BE-F34C-8B83-1670677266A1.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/BAA7792E-4861-0240-B2E8-8ADCFC0D75D4.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/BAACDC1B-7D41-2E43-A0A8-B9CC17BD58FA.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/BD44AC20-9752-1842-8BDD-11B3D9A7F959.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C131B97A-C8ED-2B41-B188-F5D8E5EDDCD8.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C289EC9F-6137-B449-9EAA-C9F7C62D48DC.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C2EA7AA0-DD62-EC4D-A97D-A3E1770C4AC5.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C340D20B-32AE-6C47-B38F-633BB7970E9A.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C34AAAF3-E400-1C46-B778-03F3C8B953A7.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C3856C01-9E0C-FB4C-BC13-CF333A32A963.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C581EEDA-365D-8243-AA9B-487474E836E3.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C70DEEFB-5CAA-D941-A44C-D40811A51A43.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C71B4285-6607-7C48-AEAF-0D1E01FE80B9.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C75DCDA5-71ED-9D4E-9C26-D33D502E96DA.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C806AA99-562F-7446-84D5-6D1C8CC5A96C.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C8736A1E-7E5C-F641-BFA4-E542027A4A65.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C8F4B3BC-30D7-9E40-99FB-1D179C1A72F0.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/C99FAE2F-A81F-494E-8906-13C6C6B57C9A.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/CB90BCA9-C3A6-3343-8C5D-3B9CD2012D77.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/CC0D0A0C-F7EF-D943-BF99-1A5449E681EA.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/CDD4AEEA-C13C-9343-87B3-D3901D086AB8.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/D19CAD37-DC4B-9D43-845A-022F0EC49425.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/D215AC3B-4815-9443-8971-D487E7367E2C.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/D3D79859-F7A4-7F43-9CB9-C7F2AF82AD18.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/D6D09A3D-FCA5-0F4F-96D7-B7C6D60CE1F2.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/D82E0EA0-5CF3-4B4A-A25F-A68F7B4880D7.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/D8C1771F-6731-AC43-A443-736E66AF6E8B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/DA8B29C6-106D-8E49-A6EE-41EE67F1EBB6.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/DB8D3D2E-4390-B542-89AA-3CCF530FB307.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/DC966EF5-4969-1C41-A427-C24B7BFD0125.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/DC985944-2EBF-5344-A1D7-AA5B26156DAD.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/DCFC538F-FF51-5A4D-AF54-39D0C159271D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/DE487056-664D-BC4A-B2B7-3D18FDD801DC.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/E086BC9D-79A9-1943-9BDD-78FDF5326C7E.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/E12BFE3D-1266-0540-B5FE-2E03CB055D72.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/E1DB956C-38DF-8F46-A231-6801FF3BC549.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/E2084E7A-98CB-C542-BAC3-FFC22ED4C571.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/E24E86CE-0B2D-854D-A42D-EA1AA389316D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/E371A98B-387E-7B40-928E-A66D704FF034.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/E3F90874-EBC4-1944-BF0F-9A571307B515.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/E42CF18F-07E9-E444-AB5E-A5D3DAA83425.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/E4797596-5308-E740-A731-779A412BC5AB.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/E4A0848C-DF7B-EB4F-A4AE-FEE77FCAF780.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/E5EABD17-62DC-3746-A9C4-1F0FFA7606EB.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/E62CBB8D-97FF-1741-BCF2-BB266E024275.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/EA29823C-B8A1-364D-A318-6ECDB2250EA7.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/EBE2860D-B054-444F-A344-0732E9FFD954.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/EBFA7E89-F442-0B42-B5B5-01A4A708DC19.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/ECC000DF-4905-D24B-927C-CB5D01FADC0D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/EDE0846D-2A74-0B44-B1F6-EBF792DC5490.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/EE6E4C51-A469-7741-94F6-9F8EF6E874E9.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/EF0BA10E-6EAA-ED45-A20D-55626BB9EED2.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/EFA324AE-C520-1F4F-8FA7-794CE841282D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/EFAB90DE-FE83-3B4B-9D5C-6FF9311FF7B5.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/F039F6C0-39E4-574D-BB19-7BF8158BC218.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/F28D58B7-7D20-D849-B1AB-8E0C641F34F3.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/F332D3E1-8858-DC4D-A72D-49745F1F4946.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/F3BE4A93-6499-EF48-BBF8-766EB134849D.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/F7958108-A321-2043-95DD-40E57C5BA5C1.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/F9F7AC89-8467-CC4F-BF1C-2D2A3EBA6AE1.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/FA12F918-8DD1-D847-AE38-24087289ACC0.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/FA2B7671-BE11-964F-8A3A-E9B5FBD3EB3B.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/FA5DAC57-9463-CA4C-ABE3-83E86E94D6BB.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/FB4B38B3-9391-504F-ACBF-9EA5EA4B6664.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/FC0857D5-8883-3642-8342-756D5398D8FE.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/FD1E47F2-94D9-5F47-A455-8A02AC398D36.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/FD71F2AB-C482-364A-9DFC-35015F1C48AD.root",
"/eos/user/a/aera/andrew/data/nano/ttsemi/2016/FDE701D2-9FDA-534B-AD08-0423A452F98A.root"],
"TTBarHad": [
"/eos/user/a/amlevin/data/nano/tthad/2016/02D3C041-82C6-BC49-BE63-FA1FC6FEA67F.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/02FFB3FD-8340-BB4A-A7F6-817AD991B569.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/04D88F01-CC1E-7D41-A559-DE6FE6DCE503.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/0533C1CB-0FDC-A14F-A957-4D21F41C5D6D.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/07F3F40B-7DCF-154F-A28F-28D9BBC41838.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/0B885551-5C5B-F84F-8C05-3127670817A1.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/1690A2E7-2093-B042-84C3-B7F31A93BAE0.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/18BF8110-D450-5D4F-A7F8-41AB9E673500.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/241F72CE-7D3B-0447-8BB3-4F5EC23762B1.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/257A66C9-549D-1D4D-BB76-D10D7C4F8261.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/27CC3A71-D800-BE4A-A8EA-9ED1475B5DE4.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/28B99BD1-A138-0D41-87E1-6766C8FC9F62.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/29B45CE6-9F66-484D-92C4-1BA046E6C789.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/30501AA6-9996-A04F-B5C1-D7D0CB4731EE.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/3B829CE8-F37F-BB4D-8249-40300D34CB25.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/3BEB557E-B073-DB40-A6DD-FE1339A7D475.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/4A86C3CF-E357-CE43-B58D-988969908CA7.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/575DF761-312E-9E42-8A91-95BF8B05A906.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/5A6BA99F-8F8B-EA42-B8D2-C4E56DE3493D.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/5BE6DD20-F4D2-4F42-A2EE-8167401ED77A.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/5F60CEDF-CDF2-D740-84D9-CBE7E471565C.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/60BB9E52-D45F-F245-9DBE-5B975741B63E.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/68429ECB-961A-1E45-9956-EDD931CA9356.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/6950E54A-1ED4-0B4B-AB5B-5DD2898931E0.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/6DB41A0D-0CF4-1E4A-93D3-E74EB7DEDF51.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/6EB6A93D-F0BA-4341-B37A-0DB9CFE7D841.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/70A416BD-91E6-3748-88AF-1AB3541AF586.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/751B4F9C-190D-A245-895A-BC713FD8273B.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/7A5A4F23-541D-3149-8596-867F773DF518.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/7EA34091-6754-2446-992D-9AEE392B9CF5.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/83299B8F-6ECC-1A49-A01E-8400B3A639FD.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/851C2D8A-58CC-1343-9FE4-D59218080406.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/8608255F-1A40-A04F-BEC7-C11C33D7FC45.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/86331BBA-CAF2-774F-9066-0BD5D64C5668.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/88C6FAA5-FD65-764D-9FB7-DADB2825716D.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/8E578FC4-03CE-F445-B8CE-337B16F20A3B.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/8FAD84A6-8941-834A-9F5D-208C557F89B9.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/9453AE91-696D-7449-B73A-3654ACF38252.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/9F998C26-1E31-7944-9383-5583D1BDB2A7.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/A54524C1-7F4E-6546-A3BA-7D2487CD817C.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/A634885A-E9A4-5A42-8C61-9A37D3E9F7D4.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/A9C58A7A-9CEC-1843-A135-64F2F343428B.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/B084A8F1-956B-304C-9BC9-D576E9A0F395.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/B514CE93-05EE-0346-99F9-F4F19303EDE9.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/BA9580A1-0A57-714F-BF50-E793D3B64410.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/C3EC8C42-C48A-1A47-9ADB-E347703F282C.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/C91CBF95-0732-3342-AB41-D0BC62BD90DD.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/CC0566AE-6E5D-9D4D-A098-B8ED48186FE4.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/CE0F4D29-8629-F145-BF1B-AB6FFDA418A4.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/DD92E735-B4E5-9842-A0B5-8DA3DD369247.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/E2C38295-3CA1-3A4C-864C-2D9FF9ED9555.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/E9246021-CDAA-FC4E-9DDB-6E1AFAC50F9F.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/EB5D2B94-D1DA-3C45-A09C-31D4FC1FA12E.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/ECA26E33-5625-334E-9982-71F60514937C.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/F176A370-DB7B-5140-B7C3-B24149C078E4.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/F2F89FBE-9906-3243-A48C-08FEB5C5C19A.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/F536211E-3833-1C46-BC71-DAA4FC91CAF7.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/F5F4B458-5790-A949-AB66-CF0632D120B5.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/F62287E5-0083-BA47-BD8A-CDF7161C328D.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/F9B7CF60-DDF9-F44D-8FD3-44DB365F36B2.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/FACDA9D8-924B-C044-8D35-56CB2A517CAF.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/FCDD3859-6C76-2148-9710-3DCED5CF81E3.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/004BAF76-21C7-2345-AE72-91576585A54F.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/0AB6839C-1BE6-A643-B26A-59CE8CF6C336.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/0BF74053-5B53-7744-BA31-7D08CE71A3CF.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/175EF5FF-1C48-E547-98E4-D9745DB1804A.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/1CE557BC-1E7C-F547-A379-BA0B83A90B4E.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/1DA07572-2E05-0C48-BC9E-FFE51014CA9E.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/22101332-6274-6B43-BBE6-6879367E8849.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/24F53C0F-FC77-2742-B4F8-B9CB24916C3B.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/2605CE21-AB56-2240-B40F-3277E2A963EC.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/2B3F5867-1C26-D94B-844E-4D7E72B86594.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/314EA76F-1A84-A64C-9115-86A23AD93A50.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/37D0EDAE-27A1-3B4C-86B2-0E01C5CF08CC.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/4358D38F-BBD7-8544-BE99-495F32C03853.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/516D7828-DD76-9444-911D-65FF913D5853.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/54AEAC26-2CB2-414A-A606-B26F43041E4E.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/6CBB0B2D-917B-0343-B323-15026F6B2941.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/737B5A6C-CAFE-8849-9B18-9E3F355F91E3.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/7FE7FB74-71A8-6B41-8140-6AD179FD0065.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/802378C2-A418-FF42-9CBE-7B9EFCC25416.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/821CC7A5-40C9-FE49-B762-FD4C250D1B92.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/82C0A405-EB70-0740-ABAF-3E7F2DB5C56D.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/866C548E-77D0-FD42-B1C7-6EB5E8046299.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/88E8D1A6-7E79-644B-A87E-27E249F0E83B.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/9A9E6F98-F2FD-AF4D-BB8A-E47212C3BAC1.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/9CA0344A-0478-8D4A-8CDD-90AF91AB9043.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/A3D430C6-1980-504C-AF7C-62980AF73C2E.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/B4B00882-648B-6649-83D7-319F47583383.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/C4DC9C0A-64B3-454C-877E-F393E56B037A.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/C968C36B-C2D8-9B47-8582-B65D283ABBB7.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/CD923292-5E72-9145-A95F-E3ED7408AE1D.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/CEBD58F1-AD07-2E49-A0A6-499C2DF6DA33.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/D073105B-7730-154F-B319-E912DA3BF7BD.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/E681C0E3-3F30-F349-A54D-89DED0946267.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/ECF0CC2F-B8CA-E94F-A6CA-616139C9D554.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/EEC7D98C-CEDF-A244-9830-3EB49B9C73A5.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/F14CFB5B-6C53-AC47-A6C8-311F6F71814E.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/F44BE850-88C4-4A4C-A5D4-948F3E382F5B.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/F5124192-2E6F-8D45-BBAB-A4DE17873C75.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/14D1B390-C6AB-7F4E-B780-0BDD579E2CD0.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/1AC9A787-4585-BE4A-B74D-B95863FCBC02.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/1E20CD7D-30D5-C949-8B2F-1823788BA729.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/29740D81-4B3F-8D4C-A44D-49EE77CD0633.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/2F0D05E1-510B-6745-81FB-F5903A08A916.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/40D601FE-FDCE-A147-995A-CDE429035908.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/4185235E-5605-A04C-A0F4-DA44BBF9048E.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/530BE54B-F74E-4B45-9EC4-C378DF5F8E89.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/5D29A844-191B-2648-AC78-3EE299F03412.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/6359C85E-2F52-914B-B9D5-7E0715BFCF5C.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/683A729C-7BBB-0B4F-BE4F-09C230760097.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/74F51E18-1512-954F-AC56-569593A6E125.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/760EBFC2-11DA-8446-9C0E-FB1EF89C203E.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/912AD694-B04A-1F49-98B3-F43C41BF6565.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/91746986-7A00-4F42-A436-7741849A56AC.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/9A759338-51C9-5F44-B7CD-20AD76E2846F.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/CD77C409-4020-E64E-9525-13ED2B516FBD.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/1A9CB2EB-3907-BA49-9E9B-45572FB51C72.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/27D3E0CE-9FC5-1D4C-BC07-989FEB5B2B68.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/2AD1BD66-077E-3C4A-9C5D-CFC0A83A54D9.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/2F1DC8F2-E831-5440-A790-EABA21F5B39A.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/37757CDA-FC83-D04E-9329-9C54CD5BEB62.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/3B24D6D2-01C6-4D48-B3D2-0E6B4ECE4628.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/3E8269B0-48DF-C74B-A6CC-C692F7B2FAD4.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/42E12133-3B2F-8F46-A3FE-5DAD9F0DC73B.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/480C3721-B4F4-7346-9B06-A96D5B7D9B40.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/4972F8F7-47D1-AB44-8B1F-2DA3CA614062.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/4C631976-7D9D-5847-B559-9F7742F50AD1.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/4DB56BBF-05D1-D64A-AD82-710F325E4CB9.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/52B3B9DC-8E26-DD43-8FDE-CAC746726BEC.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/56C480C6-65C3-CB47-B1B5-92B3985ADEC5.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/58C43E46-845C-7B41-9012-51E05DEA4619.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/5D65561D-C351-6141-8A5C-01997F6F54FF.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/60692F68-BEEF-6A46-AA27-AA05489674E9.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/62AE6EFF-2C0E-8041-9F58-BE17D8CFE67E.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/6388BB0C-E3D8-6446-8126-6AFDBF3ADEDE.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/66660185-ED3F-1A48-9D11-3248964EEBE3.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/6F6E8296-2988-C04C-95BA-61FE35065ACA.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/7572F666-7524-DF4A-932E-CD55BEFAF164.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/81F3F474-A787-004B-B712-9EE90D8E2DA6.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/A1C37D60-8A3F-E241-B2FE-4B9621F3C12A.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/A4D67F9D-CC6E-4048-B7F7-99E7E51F67B7.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/A6F2D061-2B7C-F945-A98B-1A850EE3CE36.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/B3E4B6CA-45B2-6A42-82E2-58A1EFAB2FE3.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/BD5FA02D-DD1F-1C45-A2A4-855A782E6460.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/BE0C6A06-C03E-C144-BB22-9CCDC8CFA018.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/BE9F55E2-9AAC-B24A-BA8D-B601FDBB9EFC.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/BEE3F240-861F-6442-82BE-F7431D9E9A71.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/C0DA5779-427D-6442-A448-CE1981CBD773.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/C1895281-4008-4147-93E2-7C9C80A1FE7B.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/C7A26094-7A22-1C44-9FFE-6E6EE6E39DCF.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/CA3CD3F5-785D-D846-9B7A-9362B1DA0D8C.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/CEE0EA97-77BA-D84F-980B-466D52FF8F57.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/D5FFC8D6-347B-0848-9426-1367CF23A7A4.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/D9172EA5-9F68-6A4A-875D-BB0CE213A113.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/DFE61B6D-8179-E143-A5FF-1A1D977D7CAE.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/E8240ED8-D6E8-3547-BBC4-96007F02A3DA.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/F01F108F-B011-DF4A-AF94-3A64CFF48F0E.root",
"/eos/user/a/amlevin/data/nano/tthad/2016/EE2AC44A-69E7-1E40-BAAD-C3CEDA720813.root"]
}

result = processor.run_uproot_job(
    samples_local,
    "Events",
    EwkwhjjProcessor(),
    processor.futures_executor,
    {"schema": NanoAODSchema, "workers": 10},
    chunksize=10000000,
)

if "Signal" in result["nevents"]:
    print("result[\"nevents\"][\"Signal\"] = {}".format(result["nevents"]["Signal"]))

if "TTBarSemi" in result["nevents"]:
    print("result[\"nevents\"][\"TTBarSemi\"] = {}".format(result["nevents"]["TTBarSemi"]))

if "QCDWHJJ" in result["nevents"]:
    print("result[\"nevents\"][\"QCDWHJJ\"] = {}".format(result["nevents"]["QCDWHJJ"]))

if "W" in result["nevents"]:
    print("result[\"nevents\"][\"W\"] = {}".format(result["nevents"]["W"]))

if "WW" in result["nevents"]:
    print("result[\"nevents\"][\"WW\"] = {}".format(result["nevents"]["WW"]))

if "TTBarHad" in result["nevents"]:
    print("result[\"nevents\"][\"TTBarHad\"] = {}".format(result["nevents"]["TTBarHad"]))

from coffea.util import save

save(result,"outfile")
