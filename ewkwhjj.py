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

parser = argparse.ArgumentParser()

parser.add_argument('--year',dest='year',default='2016')
parser.add_argument('--nprocesses',dest='nprocesses',type=int,default='10')

args = parser.parse_args()

assert(args.year == '2016' or args.year == '2017' or args.year == '2018')

year = args.year

bst = xgb.Booster({'nthread': 1})

bst.load_model('/afs/cern.ch/user/a/amlevin/ewkwhjj/model.bin')

bst_merged = xgb.Booster({'nthread': 1})

bst_merged.load_model('/afs/cern.ch/user/a/amlevin/ewkwhjj/merged.model')

if year == '2016':
    lumimask = LumiMask('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt')
elif year == '2017':
    lumimask = LumiMask('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt')
elif year == '2018':
    lumimask = LumiMask('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt')
else:
    assert(0)

from coffea.lookup_tools import extractor

ext = extractor()
ext.add_weight_sets(['muonidsf NUM_TightID_DEN_TrackerMuons_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root'])
ext.add_weight_sets(['muonidsfunc NUM_TightID_DEN_TrackerMuons_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root'])
ext.add_weight_sets(['muonisosf NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root'])
ext.add_weight_sets(['muonisosfunc NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root'])
ext.add_weight_sets(['muonhltsf NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root'])
ext.add_weight_sets(['muonhltsfunc NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root'])
ext.add_weight_sets(['electronidsf EGamma_SF2D /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Ele_Medium_preVFP_EGM2D.root'])
ext.add_weight_sets(['electronidsfunc EGamma_SF2D_error /afs/cern.ch/user/a/amlevin/ewkwhjj/data/Ele_Medium_preVFP_EGM2D.root'])
ext.add_weight_sets(['pileup ratio_2016 /afs/cern.ch/user/a/amlevin/ewkwhjj/data/pileup.root'])
ext.add_weight_sets(['pileup_up ratio_2016_up /afs/cern.ch/user/a/amlevin/ewkwhjj/data/pileup.root'])
ext.finalize()

evaluator = ext.make_evaluator()

@numba.njit
def delta_r(a,b):
    dphi = (a.phi - b.phi + 3.14) % (2 * 3.14) - 3.14
    return math.sqrt((a.eta - b.eta) ** 2 + dphi ** 2)

@numba.njit
def select_events_merged (events,dataset,builder,syst='nominal'):

    for i0 in range(len(events)):
        
        builder.begin_list() 
        
        if dataset == 'singlemuon':
            if year == "2016":
                if not events[i0].HLT.IsoTkMu24 and not events[i0].HLT.IsoMu24:
                    builder.end_list()
                    continue
            elif year == "2017":
                if not events[i0].HLT.IsoMu27:
                    builder.end_list()
                    continue
            elif year == "2018":
                if not events[i0].HLT.IsoMu24:
                    builder.end_list()
                    continue
        elif dataset == 'singleelectron':
            if year == "2016":
                if events[i0].HLT.IsoTkMu24 or events[i0].HLT.IsoMu24 or not events[i0].HLT.Ele27_WPTight_Gsf :
                    builder.end_list()
                    continue
            elif year == "2017":
                if events[i0].HLT.IsoMu27 or not events[i0].HLT.Ele32_WPTight_Gsf_L1DoubleEG:
                    builder.end_list()
                    continue
            elif year == "2018":
                if events[i0].HLT.IsoMu24 or not events[i0].HLT.Ele32_WPTight_Gsf:
                    builder.end_list()
                    continue
        else: #MC
            if year == "2016":
                if not events[i0].HLT.IsoTkMu24 and not events[i0].HLT.IsoMu24 and not events[i0].HLT.Ele27_WPTight_Gsf:
                    builder.end_list()
                    continue
            elif year == "2017":
                if not events[i0].HLT.IsoMu27 and not events[i0].HLT.Ele32_WPTight_Gsf_L1DoubleEG:
                    builder.end_list()
                    continue
            elif year == "2018":
                if not events[i0].HLT.Ele32_WPTight_Gsf and not events[i0].HLT.IsoMu24:
                    builder.end_list()
                    continue
        
        if syst == 'nominal':        
            if events[i0].PuppiMET.pt < 30:
                builder.end_list()
                continue
        elif syst == 'JESUp':
            if events[i0].PuppiMET.ptJESUp < 30:
                builder.end_list()
                continue
        elif syst == 'JERUp':
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
def select_events_resolved (events,dataset,builder,syst='nominal'):
    
    for i0 in range(len(events)):
        
        builder.begin_list() 

        if dataset == 'singlemuon':
            if year == "2016":
                if not events[i0].HLT.IsoTkMu24 and not events[i0].HLT.IsoMu24:
                    builder.end_list()
                    continue
            elif year == "2017":
                if not events[i0].HLT.IsoMu27:
                    builder.end_list()
                    continue
            elif year == "2018":
                if not events[i0].HLT.IsoMu24:
                    builder.end_list()
                    continue
        elif dataset == 'singleelectron':
            if year == "2016":
                if events[i0].HLT.IsoTkMu24 or events[i0].HLT.IsoMu24 or not events[i0].HLT.Ele27_WPTight_Gsf :
                    builder.end_list()
                    continue
            elif year == "2017":
                if events[i0].HLT.IsoMu27 or not events[i0].HLT.Ele32_WPTight_Gsf_L1DoubleEG:
                    builder.end_list()
                    continue
            elif year == "2018":
                if events[i0].HLT.IsoMu24 or not events[i0].HLT.Ele32_WPTight_Gsf:
                    builder.end_list()
                    continue
        else: #MC
            if year == "2016":
                if not events[i0].HLT.IsoTkMu24 and not events[i0].HLT.IsoMu24 and not events[i0].HLT.Ele27_WPTight_Gsf:
                    builder.end_list()
                    continue
            elif year == "2017":
                if not events[i0].HLT.IsoMu27 and not events[i0].HLT.Ele32_WPTight_Gsf_L1DoubleEG:
                    builder.end_list()
                    continue
            elif year == "2018":
                if not events[i0].HLT.Ele32_WPTight_Gsf and not events[i0].HLT.IsoMu24:
                    builder.end_list()
                    continue
        
        if syst == 'nominal':        
            if events[i0].PuppiMET.pt < 30:
                builder.end_list()
                continue
        elif syst == 'JESUp':
            if events[i0].PuppiMET.ptJESUp < 30:
                builder.end_list()
                continue
        elif syst == 'JERUp':
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
            'sumw': processor.defaultdict_accumulator(float),
            'nevents': processor.defaultdict_accumulator(float),
            'sel1_bdtscore_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_pileupUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_prefireUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_electronidsfUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_muonidsfUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_muonisosfUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_muonhltsfUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_JESUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_JERUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning2': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 40, -0.5, 1.5),
            ),
            'sel1_bdtscore_binning3': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 21, 0, 1.05),
            ),
            'sel1_higgsdijetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsdijetmass', 'Higgs dijet mass [GeV]', 75, 0, 300),
            ),
            'sel1_higgsdijetpt_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsdijetpt', 'Higgs dijet pt [GeV]', 20, 0, 200),
            ),
            'sel1_vbsdijetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetmass', 'VBS dijet mass [GeV]', 19, 100, 2000),
            ),
            'sel1_vbsdijetmass_binning2': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetmass', 'VBS dijet mass [GeV]', 39, 100, 4000),
            ),
            'sel1_vbsdijetabsdeta_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetabsdeta', 'VBS dijet $\Delta \eta$', 55, 2.5, 8),
            ),
            'sel1_leptonpt_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('leptonpt', 'Lepton pt [GeV]', 19, 20, 200),
            ),
            'sel1_met_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('met', 'MET [GeV]', 20, 0, 200),
            ),
            'sel2_bdtscore_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_pileupUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_prefireUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_electronidsfUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_JESUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_JERUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning2': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 40, -0.5, 1.5),
            ),
            'sel2_bdtscore_binning3': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 21, 0, 1.05),
            ),
            'sel2_higgsdijetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsdijetmass', 'Higgs dijet mass [GeV]', 75, 0, 300),
            ),
            'sel2_higgsdijetpt_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsdijetpt', 'Higgs dijet pt [GeV]', 20, 0, 200),
            ),
            'sel2_vbsdijetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetmass', 'VBS dijet mass [GeV]', 19, 100, 2000),
            ),
            'sel2_vbsdijetmass_binning2': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetmass', 'VBS dijet mass [GeV]', 39, 100, 4000),
            ),
            'sel2_vbsdijetabsdeta_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetabsdeta', 'VBS dijet $\Delta \eta$', 55, 2.5, 8),
            ),
            'sel2_leptonpt_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('leptonpt', 'Lepton pt [GeV]', 19, 20, 200),
            ),
            'sel2_met_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('met', 'MET [GeV]', 20, 0, 200),
            ),
            'sel3_bdtscore_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_pileupUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_prefireUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_electronidsfUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_muonidsfUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_muonisosfUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_muonhltsfUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_JESUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_JERUp': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning2': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 40, -0.5, 1.5),
            ),
            'sel3_bdtscore_binning3': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 21, 0, 1.05),
            ),
            'sel3_higgsdijetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsdijetmass', 'Higgs dijet mass [GeV]', 75, 0, 300),
            ),
            'sel3_higgsdijetpt_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsdijetpt', 'Higgs dijet pt [GeV]', 20, 0, 200),
            ),
            'sel3_vbsdijetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetmass', 'VBS dijet mass [GeV]', 19, 100, 2000),
            ),
            'sel3_vbsdijetmass_binning2': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetmass', 'VBS dijet mass [GeV]', 39, 100, 4000),
            ),
            'sel3_vbsdijetabsdeta_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetabsdeta', 'VBS dijet $\Delta \eta$', 55, 2.5, 8),
            ),
            'sel3_leptonpt_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('leptonpt', 'Lepton pt [GeV]', 19, 20, 200),
            ),
            'sel3_met_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('met', 'MET [GeV]', 20, 0, 200),
            ),
            'sel4_higgsdijetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsdijetmass', 'Higgs dijet mass [GeV]', 75, 0, 300),
            ),
            'sel4_higgsdijetpt_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsdijetpt', 'Higgs dijet pt [GeV]', 20, 0, 200),
            ),
            'sel4_vbsdijetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetmass', 'VBS dijet mass [GeV]', 19, 100, 2000),
            ),
            'sel4_vbsdijetmass_binning2': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetmass', 'VBS dijet mass [GeV]', 4, 100, 500),
            ),
            'sel4_leptonpt_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('leptonpt', 'Lepton pt [GeV]', 19, 20, 200),
            ),
            'sel4_met_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('met', 'MET [GeV]', 20, 0, 200),
            ),
            'sel5_higgsdijetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsdijetmass', 'Higgs dijet mass [GeV]', 75, 0, 300),
            ),
            'sel5_higgsdijetpt_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsdijetpt', 'Higgs dijet pt [GeV]', 20, 0, 200),
            ),
            'sel5_vbsdijetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetmass', 'VBS dijet mass [GeV]', 19, 100, 2000),
            ),
            'sel5_vbsdijetmass_binning2': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetmass', 'VBS dijet mass [GeV]', 4, 100, 500),
            ),
            'sel5_leptonpt_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('leptonpt', 'Lepton pt [GeV]', 19, 20, 200),
            ),
            'sel5_met_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('met', 'MET [GeV]', 20, 0, 200),
            ),
            'sel6_higgsdijetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsdijetmass', 'Higgs dijet mass [GeV]', 75, 0, 300),
            ),
            'sel6_higgsdijetpt_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsdijetpt', 'Higgs dijet pt [GeV]', 20, 0, 200),
            ),
            'sel6_vbsdijetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetmass', 'VBS dijet mass [GeV]', 19, 100, 2000),
            ),
            'sel6_vbsdijetmass_binning2': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('vbsdijetmass', 'VBS dijet mass [GeV]', 4, 100, 500),
            ),
            'sel6_leptonpt_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('leptonpt', 'Lepton pt [GeV]', 19, 20, 200),
            ),
            'sel6_met_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('met', 'MET [GeV]', 20, 0, 200),
            ),
            'sel7_higgsjetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsjetmass', 'Higgs jet mass [GeV]', 75, 0, 300),
            ),
            'sel7_higgsjetsoftdropmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsjetsoftdropmass', 'Higgs jet soft drop mass [GeV]', 75, 0, 300),
            ),
            'sel7_bdtscore_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1.0),
            ),
            'sel7_bdtscore_binning2': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 40, -0.5, 1.5),
            ),
            'sel7_bdtscore_binning3': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 21, 0, 1.05),
            ),
            'sel8_higgsjetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsjetmass', 'Higgs jet mass [GeV]', 75, 0, 300),
            ),
            'sel8_higgsjetsoftdropmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsjetsoftdropmass', 'Higgs jet soft drop mass [GeV]', 75, 0, 300),
            ),
            'sel8_bdtscore_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1.0),
            ),
            'sel8_bdtscore_binning2': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 40, -0.5, 1.5),
            ),
            'sel8_bdtscore_binning3': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 21, 0, 1.05),
            ),
            'sel9_higgsjetmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsjetmass', 'Higgs jet mass [GeV]', 75, 0, 300),
            ),
            'sel9_higgsjetsoftdropmass_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('higgsjetsoftdropmass', 'Higgs jet soft drop mass [GeV]', 75, 0, 300),
            ),
            'sel9_bdtscore_binning1': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 20, 0, 1.0),
            ),
            'sel9_bdtscore_binning2': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 40, -0.5, 1.5),
            ),
            'sel9_bdtscore_binning3': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('bdtscore', 'BDT score', 21, 0, 1.05),
            ),

        })

    @property
    def accumulator(self):
        return self._accumulator

#    @numba.njit
    def process(self, events):

        output = self.accumulator.identity()

        dataset = events.metadata['dataset']

        if 'single' not in dataset:
            output['sumw'][dataset] += ak.sum(np.sign(events.Generator.weight))
        output['nevents'][dataset] += len(events)

        if 'single' in dataset :
            events = events[lumimask(events.run,events.luminosityBlock)]
   
        particleindices = select_events_resolved(events,dataset,ak.ArrayBuilder()).snapshot()

        if 'single' not in dataset:
            particleindices_JESUp = select_events_resolved(events,dataset,ak.ArrayBuilder(),syst='JESUp').snapshot()
            particleindices_JERUp = select_events_resolved(events,dataset,ak.ArrayBuilder(),syst='JERUp').snapshot()

        particleindices_merged = select_events_merged(events,dataset,ak.ArrayBuilder()).snapshot()
        
        basecut = ak.num(particleindices) != 0

        if 'single' not in dataset:
            basecut_JESUp = ak.num(particleindices_JESUp) != 0
            basecut_JERUp = ak.num(particleindices_JERUp) != 0

        basecut_merged = ak.num(particleindices_merged) != 0

        if 'single' in dataset:
            dataset = 'Data'

        if ak.any(basecut_merged):
            particleindices_merged = particleindices_merged[basecut_merged]
            events_merged = events[basecut_merged]
            events_merged.FatJet[particleindices_merged['0']]
            jets_merged = [events_merged.Jet[particleindices_merged[idx]] for idx in '12']
            cut7 = ak.firsts((events_merged.FatJet[particleindices_merged['0']].mass > 50) & (events_merged.FatJet[particleindices_merged['0']].mass < 150) & ((jets_merged[0]+jets_merged[1]).mass > 500) & (abs(jets_merged[0].eta - jets_merged[1].eta) > 2.5) & (particleindices_merged['3'] != -1))
            cut8 = ak.firsts((events_merged.FatJet[particleindices_merged['0']].mass > 50) & (events_merged.FatJet[particleindices_merged['0']].mass < 150) & ((jets_merged[0]+jets_merged[1]).mass > 500) & (abs(jets_merged[0].eta - jets_merged[1].eta) > 2.5) & (particleindices_merged['4'] != -1))
#            cut9_merged = cut7_merged | cut8_merged

        if dataset != 'Data' and ak.any(basecut_JESUp):
            particleindices_JESUp = particleindices_JESUp[basecut_JESUp]
            events_JESUp = events[basecut_JESUp]
            jets_JESUp = [events_JESUp.Jet[particleindices_JESUp[idx]] for idx in '0123']
            cut1_JESUp = ak.firsts(((jets_JESUp[0]+jets_JESUp[1]).mass > 50) & ((jets_JESUp[0]+jets_JESUp[1]).mass < 150) & ((jets_JESUp[2]+jets_JESUp[3]).mass > 500) & (abs(jets_JESUp[2].eta - jets_JESUp[3].eta) > 2.5) & (particleindices_JESUp['4'] != -1))
            cut2_JESUp = ak.firsts(((jets_JESUp[0]+jets_JESUp[1]).mass > 50) & ((jets_JESUp[0]+jets_JESUp[1]).mass < 150) & ((jets_JESUp[2]+jets_JESUp[3]).mass > 500) & (abs(jets_JESUp[2].eta - jets_JESUp[3].eta) > 2.5) & (particleindices_JESUp['5'] != -1))
#            cut3_JESUp = cut1_JESUp | cut2_JESUp

        if dataset != 'Data' and ak.any(basecut_JERUp):
            particleindices_JERUp = particleindices_JERUp[basecut_JERUp]
            events_JERUp = events[basecut_JERUp]
            jets_JERUp= [events_JERUp.Jet[particleindices_JERUp[idx]] for idx in '0123']        
            cut1_JERUp = ak.firsts(((jets_JERUp[0]+jets_JERUp[1]).mass > 50) & ((jets_JERUp[0]+jets_JERUp[1]).mass < 150) & ((jets_JERUp[2]+jets_JERUp[3]).mass > 500) & (abs(jets_JERUp[2].eta - jets_JERUp[3].eta) > 2.5) & (particleindices_JERUp['4'] != -1))
            cut2_JERUp = ak.firsts(((jets_JERUp[0]+jets_JERUp[1]).mass > 50) & ((jets_JERUp[0]+jets_JERUp[1]).mass < 150) & ((jets_JERUp[2]+jets_JERUp[3]).mass > 500) & (abs(jets_JERUp[2].eta - jets_JERUp[3].eta) > 2.5) & (particleindices_JERUp['5'] != -1))
#            cut3_JERUp = cut1_JERUp | cut2_JERUp

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

            sel7_muonidsf = ak.firsts(evaluator['muonidsf'](abs(sel7_muons.eta), sel7_muons.pt))
            sel7_muonisosf = ak.firsts(evaluator['muonisosf'](abs(sel7_muons.eta), sel7_muons.pt))
            sel7_muonhltsf = ak.firsts(evaluator['muonhltsf'](abs(sel7_muons.eta), sel7_muons.pt))

            sel7_X = pandas.DataFrame(np.transpose(np.vstack((
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

            sel7_d = xgb.DMatrix(sel7_X)

            sel7_bdtscore = bst_merged.predict(sel7_d)


            if dataset == 'Data':
                sel7_weight = np.ones(len(sel7_events))
            else:    
                sel7_weight = np.sign(sel7_events.Generator.weight)*sel7_events.L1PreFiringWeight.Nom*sel7_muonidsf*sel7_muonisosf*sel7_muonhltsf

            output['sel7_higgsjetmass_binning1'].fill(
                dataset=dataset,
                higgsjetmass=ak.firsts(sel7_events.FatJet[sel7_particleindices['0']].mass),
                weight=sel7_weight
            )

            output['sel7_higgsjetsoftdropmass_binning1'].fill(
                dataset=dataset,
                higgsjetsoftdropmass=ak.firsts(sel7_events.FatJet[sel7_particleindices['0']].msoftdrop),
                weight=sel7_weight
            )

            output['sel7_bdtscore_binning1'].fill(
                dataset=dataset,
                bdtscore=sel7_bdtscore,
                weight=sel7_weight
            )

            output['sel7_bdtscore_binning2'].fill(
                dataset=dataset,
                bdtscore=sel7_bdtscore,
                weight=sel7_weight
            )

            output['sel7_bdtscore_binning3'].fill(
                dataset=dataset,
                bdtscore=sel7_bdtscore,
                weight=sel7_weight
            )

            output['sel9_higgsjetmass_binning1'].fill(
                dataset=dataset,
                higgsjetmass=ak.firsts(sel7_events.FatJet[sel7_particleindices['0']].mass),
                weight=sel7_weight
            )

            output['sel9_higgsjetsoftdropmass_binning1'].fill(
                dataset=dataset,
                higgsjetsoftdropmass=ak.firsts(sel7_events.FatJet[sel7_particleindices['0']].msoftdrop),
                weight=sel7_weight
            )

            output['sel9_bdtscore_binning1'].fill(
                dataset=dataset,
                bdtscore=sel7_bdtscore,
                weight=sel7_weight
            )

            output['sel9_bdtscore_binning2'].fill(
                dataset=dataset,
                bdtscore=sel7_bdtscore,
                weight=sel7_weight
            )

            output['sel9_bdtscore_binning3'].fill(
                dataset=dataset,
                bdtscore=sel7_bdtscore,
                weight=sel7_weight
            )

        if ak.any(basecut_merged) and ak.any(cut8):

            sel8_particleindices = particleindices_merged[cut8]
        
            sel8_events = events_merged[cut8]
            
            sel8_fatjets = sel8_events.FatJet[sel8_particleindices['0']]

            sel8_jets = [sel8_events.Jet[sel8_particleindices[idx]] for idx in '12']

            sel8_electrons = sel8_events.Electron[sel8_particleindices['4']]

            sel8_electronidsf = ak.firsts(evaluator['electronidsf'](sel8_electrons.eta, sel8_electrons.pt))

            sel8_X = pandas.DataFrame(np.transpose(np.vstack((
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

            sel8_d = xgb.DMatrix(sel8_X)

            sel8_bdtscore = bst_merged.predict(sel8_d)

            if dataset == 'Data':
                sel8_weight = np.ones(len(sel8_events))
            else:    
                sel8_weight = np.sign(sel8_events.Generator.weight)*sel8_events.L1PreFiringWeight.Nom*sel8_electronidsf

            output['sel8_higgsjetmass_binning1'].fill(
                dataset=dataset,
                higgsjetmass=ak.firsts(sel8_events.FatJet[sel8_particleindices['0']].mass),
                weight=sel8_weight
            )

            output['sel8_higgsjetsoftdropmass_binning1'].fill(
                dataset=dataset,
                higgsjetsoftdropmass=ak.firsts(sel8_events.FatJet[sel8_particleindices['0']].msoftdrop),
                weight=sel8_weight
            )

            output['sel8_bdtscore_binning1'].fill(
                dataset=dataset,
                bdtscore=sel8_bdtscore,
                weight=sel8_weight
            )

            output['sel8_bdtscore_binning2'].fill(
                dataset=dataset,
                bdtscore=sel8_bdtscore,
                weight=sel8_weight
            )

            output['sel8_bdtscore_binning3'].fill(
                dataset=dataset,
                bdtscore=sel8_bdtscore,
                weight=sel8_weight
            )

            output['sel9_higgsjetmass_binning1'].fill(
                dataset=dataset,
                higgsjetmass=ak.firsts(sel8_events.FatJet[sel8_particleindices['0']].mass),
                weight=sel8_weight
            )

            output['sel9_higgsjetsoftdropmass_binning1'].fill(
                dataset=dataset,
                higgsjetsoftdropmass=ak.firsts(sel8_events.FatJet[sel8_particleindices['0']].msoftdrop),
                weight=sel8_weight
            )

            output['sel9_bdtscore_binning1'].fill(
                dataset=dataset,
                bdtscore=sel8_bdtscore,
                weight=sel8_weight
            )

            output['sel9_bdtscore_binning2'].fill(
                dataset=dataset,
                bdtscore=sel8_bdtscore,
                weight=sel8_weight
            )

            output['sel9_bdtscore_binning3'].fill(
                dataset=dataset,
                bdtscore=sel8_bdtscore,
                weight=sel8_weight
            )

        if dataset != 'Data' and ak.any(basecut_JESUp) and ak.any(cut1_JESUp):

            sel1_JESUp_particleindices = particleindices_JESUp[cut1_JESUp]
        
            sel1_JESUp_events = events_JESUp[cut1_JESUp]
            
            sel1_JESUp_jets = [sel1_JESUp_events.Jet[sel1_JESUp_particleindices[idx]] for idx in '0123']

            sel1_JESUp_muons = sel1_JESUp_events.Muon[sel1_JESUp_particleindices['4']]

            sel1_JESUp_pu_weight = evaluator['pileup'](sel1_JESUp_events.Pileup.nTrueInt)
            sel1_JESUp_muonidsf = ak.firsts(evaluator['muonidsf'](abs(sel1_JESUp_muons.eta), sel1_JESUp_muons.pt))
            sel1_JESUp_muonisosf = ak.firsts(evaluator['muonisosf'](abs(sel1_JESUp_muons.eta), sel1_JESUp_muons.pt))
            sel1_JESUp_muonhltsf = ak.firsts(evaluator['muonhltsf'](abs(sel1_JESUp_muons.eta), sel1_JESUp_muons.pt))

            sel1_JESUp_X = pandas.DataFrame(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel1_JESUp_particleindices['6'])).data,ak.to_numpy(ak.firsts(sel1_JESUp_particleindices['7'])).data,np.zeros(len(sel1_JESUp_events)),np.sign(ak.to_numpy(ak.firsts(sel1_JESUp_muons).charge).data+1),ak.to_numpy(ak.firsts(sel1_JESUp_muons).pt).data,ak.to_numpy(ak.firsts(sel1_JESUp_muons).eta).data,ak.to_numpy(ak.firsts(sel1_JESUp_muons).phi).data,ak.to_numpy(sel1_JESUp_events.PuppiMET.ptJESUp),ak.to_numpy(sel1_JESUp_events.PuppiMET.phiJESUp),ak.to_numpy(ak.firsts(sel1_JESUp_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_JESUp_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel1_JESUp_jets[0]+sel1_JESUp_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel1_JESUp_jets[2]+sel1_JESUp_jets[3]).mass)).data, ak.to_numpy(ak.firsts(sel1_JESUp_jets[2]).eta - ak.firsts(sel1_JESUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_JESUp_muons+sel1_JESUp_jets[0]).pt*sel1_JESUp_events.PuppiMET.ptJESUp*(1 - np.cos(sel1_JESUp_events.PuppiMET.phiJESUp - (sel1_JESUp_muons+sel1_JESUp_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_JESUp_muons+sel1_JESUp_jets[1]).pt*sel1_JESUp_events.PuppiMET.ptJESUp*(1 - np.cos(sel1_JESUp_events.PuppiMET.phiJESUp - (sel1_JESUp_muons+sel1_JESUp_jets[1]).phi))))).data))),columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel1_JESUp_d = xgb.DMatrix(sel1_JESUp_X)

            sel1_JESUp_bdtscore = bst.predict(sel1_JESUp_d)

            sel1_JESUp_weight = np.sign(sel1_JESUp_events.Generator.weight)*sel1_JESUp_pu_weight*sel1_JESUp_events.L1PreFiringWeight.Nom*sel1_JESUp_muonidsf*sel1_JESUp_muonisosf*sel1_JESUp_muonhltsf

            output['sel1_bdtscore_binning1_JESUp'].fill(
                dataset=dataset,
                bdtscore=sel1_JESUp_bdtscore,
                weight=sel1_JESUp_weight
            )

            output['sel3_bdtscore_binning1_JESUp'].fill(
                dataset=dataset,
                bdtscore=sel1_JESUp_bdtscore,
                weight=sel1_JESUp_weight
            )

        if dataset != 'Data' and ak.any(basecut_JERUp) and ak.any(cut1_JERUp):

            sel1_JERUp_particleindices = particleindices_JERUp[cut1_JERUp]
        
            sel1_JERUp_events = events_JERUp[cut1_JERUp]
            
            sel1_JERUp_jets = [sel1_JERUp_events.Jet[sel1_JERUp_particleindices[idx]] for idx in '0123']

            sel1_JERUp_muons = sel1_JERUp_events.Muon[sel1_JERUp_particleindices['4']]

            sel1_JERUp_pu_weight = evaluator['pileup'](sel1_JERUp_events.Pileup.nTrueInt)
            sel1_JERUp_muonidsf = ak.firsts(evaluator['muonidsf'](abs(sel1_JERUp_muons.eta), sel1_JERUp_muons.pt))
            sel1_JERUp_muonisosf = ak.firsts(evaluator['muonisosf'](abs(sel1_JERUp_muons.eta), sel1_JERUp_muons.pt))
            sel1_JERUp_muonhltsf = ak.firsts(evaluator['muonhltsf'](abs(sel1_JERUp_muons.eta), sel1_JERUp_muons.pt))

            sel1_JERUp_X = pandas.DataFrame(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel1_JERUp_particleindices['6'])).data,ak.to_numpy(ak.firsts(sel1_JERUp_particleindices['7'])).data,np.zeros(len(sel1_JERUp_events)),np.sign(ak.to_numpy(ak.firsts(sel1_JERUp_muons).charge).data+1),ak.to_numpy(ak.firsts(sel1_JERUp_muons).pt).data,ak.to_numpy(ak.firsts(sel1_JERUp_muons).eta).data,ak.to_numpy(ak.firsts(sel1_JERUp_muons).phi).data,ak.to_numpy(sel1_JERUp_events.PuppiMET.ptJERUp),ak.to_numpy(sel1_JERUp_events.PuppiMET.phiJERUp),ak.to_numpy(ak.firsts(sel1_JERUp_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_JERUp_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel1_JERUp_jets[0]+sel1_JERUp_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel1_JERUp_jets[2]+sel1_JERUp_jets[3]).mass)).data, ak.to_numpy(ak.firsts(sel1_JERUp_jets[2]).eta - ak.firsts(sel1_JERUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_JERUp_muons+sel1_JERUp_jets[0]).pt*sel1_JERUp_events.PuppiMET.ptJERUp*(1 - np.cos(sel1_JERUp_events.PuppiMET.phiJERUp - (sel1_JERUp_muons+sel1_JERUp_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_JERUp_muons+sel1_JERUp_jets[1]).pt*sel1_JERUp_events.PuppiMET.ptJERUp*(1 - np.cos(sel1_JERUp_events.PuppiMET.phiJERUp - (sel1_JERUp_muons+sel1_JERUp_jets[1]).phi))))).data))),columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel1_JERUp_d = xgb.DMatrix(sel1_JERUp_X)

            sel1_JERUp_bdtscore = bst.predict(sel1_JERUp_d)

            sel1_JERUp_weight = np.sign(sel1_JERUp_events.Generator.weight)*sel1_JERUp_pu_weight*sel1_JERUp_events.L1PreFiringWeight.Nom*sel1_JERUp_muonidsf*sel1_JERUp_muonisosf*sel1_JERUp_muonhltsf

            output['sel1_bdtscore_binning1_JERUp'].fill(
                dataset=dataset,
                bdtscore=sel1_JERUp_bdtscore,
                weight=sel1_JERUp_weight
            )

            output['sel3_bdtscore_binning1_JERUp'].fill(
                dataset=dataset,
                bdtscore=sel1_JERUp_bdtscore,
                weight=sel1_JERUp_weight
            )

        if ak.any(basecut) and ak.any(cut1):

            sel1_particleindices = particleindices[cut1]
        
            sel1_events = events[cut1]
            
            sel1_jets = [sel1_events.Jet[sel1_particleindices[idx]] for idx in '0123']

            sel1_muons = sel1_events.Muon[sel1_particleindices['4']]

            sel1_X = pandas.DataFrame(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel1_particleindices['6'])).data,ak.to_numpy(ak.firsts(sel1_particleindices['7'])).data,np.zeros(len(sel1_events)),np.sign(ak.to_numpy(ak.firsts(sel1_muons).charge).data+1),ak.to_numpy(ak.firsts(sel1_muons).pt).data,ak.to_numpy(ak.firsts(sel1_muons).eta).data,ak.to_numpy(ak.firsts(sel1_muons).phi).data,ak.to_numpy(sel1_events.PuppiMET.pt),ak.to_numpy(sel1_events.PuppiMET.phi),ak.to_numpy(ak.firsts(sel1_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel1_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel1_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel1_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel1_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel1_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel1_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel1_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel1_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel1_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel1_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel1_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel1_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel1_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel1_jets[0]+sel1_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel1_jets[2]+sel1_jets[3]).mass)).data, ak.to_numpy(ak.firsts(sel1_jets[2]).eta - ak.firsts(sel1_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_muons+sel1_jets[0]).pt*sel1_events.PuppiMET.pt*(1 - np.cos(sel1_events.PuppiMET.phi - (sel1_muons+sel1_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel1_muons+sel1_jets[1]).pt*sel1_events.PuppiMET.pt*(1 - np.cos(sel1_events.PuppiMET.phi - (sel1_muons+sel1_jets[1]).phi))))).data))),columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel1_d = xgb.DMatrix(sel1_X)

            sel1_bdtscore = bst.predict(sel1_d)    

            if dataset == 'Data':
                sel1_weight = np.ones(len(sel1_events))
            else:
                sel1_pu_weight = evaluator['pileup'](sel1_events.Pileup.nTrueInt)
                sel1_puUp_weight = evaluator['pileup_up'](sel1_events.Pileup.nTrueInt)
                sel1_muonidsf = ak.firsts(evaluator['muonidsf'](abs(sel1_muons.eta), sel1_muons.pt))
                sel1_muonisosf = ak.firsts(evaluator['muonisosf'](abs(sel1_muons.eta), sel1_muons.pt))
                sel1_muonhltsf = ak.firsts(evaluator['muonhltsf'](abs(sel1_muons.eta), sel1_muons.pt))
                sel1_muonidsfUp = ak.firsts(evaluator['muonidsfunc'](abs(sel1_muons.eta), sel1_muons.pt))+sel1_muonidsf
                sel1_muonisosfUp = ak.firsts(evaluator['muonisosfunc'](abs(sel1_muons.eta), sel1_muons.pt))+sel1_muonisosf
                sel1_muonhltsfUp = ak.firsts(evaluator['muonhltsfunc'](abs(sel1_muons.eta), sel1_muons.pt))+sel1_muonhltsf

                sel1_weight = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf
                sel1_weight_pileupUp = np.sign(sel1_events.Generator.weight)*sel1_puUp_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf
                sel1_weight_prefireUp = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Up*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf
                sel1_weight_muonidsfUp = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsfUp*sel1_muonisosf*sel1_muonhltsf
                sel1_weight_muonisosfUp = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosfUp*sel1_muonhltsf
                sel1_weight_muonhltsfUp = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsfUp

                
            output['sel1_bdtscore_binning1'].fill(
                dataset=dataset,
                bdtscore=sel1_bdtscore,
                weight=sel1_weight
            )

            if dataset != 'Data':
                output['sel1_bdtscore_binning1_pileupUp'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_pileupUp
                )

                output['sel1_bdtscore_binning1_prefireUp'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_prefireUp
                )

                output['sel1_bdtscore_binning1_muonidsfUp'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonidsfUp
                )
                
                output['sel1_bdtscore_binning1_muonisosfUp'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonisosfUp
                )

                output['sel1_bdtscore_binning1_muonhltsfUp'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonhltsfUp
                )

            output['sel1_bdtscore_binning2'].fill(
                dataset=dataset,
                bdtscore=sel1_bdtscore,
                weight=sel1_weight
            )

            output['sel1_bdtscore_binning3'].fill(
                dataset=dataset,
                bdtscore=sel1_bdtscore,
                weight=sel1_weight
            )

            output['sel1_higgsdijetmass_binning1'].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel1_jets[0]+sel1_jets[1]).mass),
                weight=sel1_weight
            )

            output['sel1_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel1_jets[0]+sel1_jets[1]).pt),
                weight=sel1_weight
            )
        
            output['sel1_vbsdijetmass_binning1'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel1_jets[2]+sel1_jets[3]).mass),
                weight=sel1_weight
            )

            output['sel1_vbsdijetmass_binning2'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel1_jets[2]+sel1_jets[3]).mass),
                weight=sel1_weight
            )

            output['sel1_vbsdijetabsdeta_binning1'].fill(
                dataset=dataset,
                vbsdijetabsdeta=ak.firsts(sel1_jets[2]).eta - ak.firsts(sel1_jets[3]).eta,
                weight=sel1_weight
            )

            output['sel1_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel1_muons.pt),
                weight=sel1_weight
            )

            output['sel1_met_binning1'].fill(
                dataset=dataset,
                met=sel1_events.PuppiMET.pt,
                weight=sel1_weight
            )

            output['sel3_bdtscore_binning1'].fill(
                dataset=dataset,
                bdtscore = sel1_bdtscore,
                weight=sel1_weight
            )

            if dataset != 'Data':

                output['sel3_bdtscore_binning1_pileupUp'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_pileupUp
                )

                output['sel3_bdtscore_binning1_prefireUp'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_prefireUp
                )
                
                output['sel3_bdtscore_binning1_electronidsfUp'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight
                )
                
                output['sel3_bdtscore_binning1_muonidsfUp'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonidsfUp
                )
                
                output['sel3_bdtscore_binning1_muonisosfUp'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonisosfUp
                )
                
                output['sel3_bdtscore_binning1_muonhltsfUp'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonhltsfUp
                )

            output['sel3_bdtscore_binning2'].fill(
                dataset=dataset,
                bdtscore = sel1_bdtscore,
                weight=sel1_weight
            )

            output['sel3_bdtscore_binning3'].fill(
                dataset=dataset,
                bdtscore=sel1_bdtscore,
                weight=sel1_weight
            )

            output['sel3_higgsdijetmass_binning1'].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel1_jets[0]+sel1_jets[1]).mass),
                weight=sel1_weight
            )

            output['sel3_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel1_jets[0]+sel1_jets[1]).pt),
                weight=sel1_weight
            )
        
            output['sel3_vbsdijetmass_binning1'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel1_jets[2]+sel1_jets[3]).mass),
                weight=sel1_weight
            )

            output['sel3_vbsdijetmass_binning2'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel1_jets[2]+sel1_jets[3]).mass),
                weight=sel1_weight
            )

            output['sel3_vbsdijetabsdeta_binning1'].fill(
                dataset=dataset,
                vbsdijetabsdeta=ak.firsts(sel1_jets[2]).eta - ak.firsts(sel1_jets[3]).eta,
                weight=sel1_weight
            )

            output['sel3_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel1_muons.pt),
                weight=sel1_weight
            )

            output['sel3_met_binning1'].fill(
                dataset=dataset,
                met=sel1_events.PuppiMET.pt,
                weight=sel1_weight
            )

        if dataset != 'Data' and ak.any(basecut_JESUp) and ak.any(cut2_JESUp):

            sel2_JESUp_particleindices = particleindices_JESUp[cut2_JESUp]
        
            sel2_JESUp_events = events_JESUp[cut2_JESUp]
            
            sel2_JESUp_jets = [sel2_JESUp_events.Jet[sel2_JESUp_particleindices[idx]] for idx in '0123']

            sel2_JESUp_electrons = sel2_JESUp_events.Electron[sel2_JESUp_particleindices['4']]

            sel2_JESUp_pu_weight = evaluator['pileup'](sel2_JESUp_events.Pileup.nTrueInt)
            sel2_JESUp_electronidsf = ak.firsts(evaluator['electronidsf'](sel2_JESUp_electrons.eta, sel2_JESUp_electrons.pt))

            sel2_JESUp_X = pandas.DataFrame(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel2_JESUp_particleindices['6'])).data,ak.to_numpy(ak.firsts(sel2_JESUp_particleindices['7'])).data,np.zeros(len(sel2_JESUp_events)),np.sign(ak.to_numpy(ak.firsts(sel2_JESUp_electrons).charge).data+1),ak.to_numpy(ak.firsts(sel2_JESUp_electrons).pt).data,ak.to_numpy(ak.firsts(sel2_JESUp_electrons).eta).data,ak.to_numpy(ak.firsts(sel2_JESUp_electrons).phi).data,ak.to_numpy(sel2_JESUp_events.PuppiMET.ptJESUp),ak.to_numpy(sel2_JESUp_events.PuppiMET.phiJESUp),ak.to_numpy(ak.firsts(sel2_JESUp_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_JESUp_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel2_JESUp_jets[0]+sel2_JESUp_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel2_JESUp_jets[2]+sel2_JESUp_jets[3]).mass)).data, ak.to_numpy(ak.firsts(sel2_JESUp_jets[2]).eta - ak.firsts(sel2_JESUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_JESUp_electrons+sel2_JESUp_jets[0]).pt*sel2_JESUp_events.PuppiMET.ptJESUp*(1 - np.cos(sel2_JESUp_events.PuppiMET.phiJESUp - (sel2_JESUp_electrons+sel2_JESUp_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_JESUp_electrons+sel2_JESUp_jets[1]).pt*sel2_JESUp_events.PuppiMET.ptJESUp*(1 - np.cos(sel2_JESUp_events.PuppiMET.phiJESUp - (sel2_JESUp_electrons+sel2_JESUp_jets[1]).phi))))).data))),columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel2_JESUp_d = xgb.DMatrix(sel2_JESUp_X)

            sel2_JESUp_bdtscore = bst.predict(sel2_JESUp_d)

            sel2_JESUp_weight = np.sign(sel2_JESUp_events.Generator.weight)*sel2_JESUp_pu_weight*sel2_JESUp_events.L1PreFiringWeight.Nom*sel2_JESUp_electronidsf

            output['sel2_bdtscore_binning1_JESUp'].fill(
                dataset=dataset,
                bdtscore=sel2_JESUp_bdtscore,
                weight=sel2_JESUp_weight
            )

            output['sel3_bdtscore_binning1_JESUp'].fill(
                dataset=dataset,
                bdtscore=sel2_JESUp_bdtscore,
                weight=sel2_JESUp_weight
            )

        if dataset != 'Data' and ak.any(basecut_JERUp) and ak.any(cut2_JERUp):

            sel2_JERUp_particleindices = particleindices_JERUp[cut2_JERUp]
        
            sel2_JERUp_events = events_JERUp[cut2_JERUp]
            
            sel2_JERUp_jets = [sel2_JERUp_events.Jet[sel2_JERUp_particleindices[idx]] for idx in '0123']

            sel2_JERUp_electrons = sel2_JERUp_events.Electron[sel2_JERUp_particleindices['4']]

            sel2_JERUp_pu_weight = evaluator['pileup'](sel2_JERUp_events.Pileup.nTrueInt)
            sel2_JERUp_electronidsf = ak.firsts(evaluator['electronidsf'](sel2_JERUp_electrons.eta, sel2_JERUp_electrons.pt))

            sel2_JERUp_X = pandas.DataFrame(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel2_JERUp_particleindices['6'])).data,ak.to_numpy(ak.firsts(sel2_JERUp_particleindices['7'])).data,np.zeros(len(sel2_JERUp_events)),np.sign(ak.to_numpy(ak.firsts(sel2_JERUp_electrons).charge).data+1),ak.to_numpy(ak.firsts(sel2_JERUp_electrons).pt).data,ak.to_numpy(ak.firsts(sel2_JERUp_electrons).eta).data,ak.to_numpy(ak.firsts(sel2_JERUp_electrons).phi).data,ak.to_numpy(sel2_JERUp_events.PuppiMET.ptJERUp),ak.to_numpy(sel2_JERUp_events.PuppiMET.phiJERUp),ak.to_numpy(ak.firsts(sel2_JERUp_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_JERUp_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel2_JERUp_jets[0]+sel2_JERUp_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel2_JERUp_jets[2]+sel2_JERUp_jets[3]).mass)).data, ak.to_numpy(ak.firsts(sel2_JERUp_jets[2]).eta - ak.firsts(sel2_JERUp_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_JERUp_electrons+sel2_JERUp_jets[0]).pt*sel2_JERUp_events.PuppiMET.ptJERUp*(1 - np.cos(sel2_JERUp_events.PuppiMET.phiJERUp - (sel2_JERUp_electrons+sel2_JERUp_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_JERUp_electrons+sel2_JERUp_jets[1]).pt*sel2_JERUp_events.PuppiMET.ptJERUp*(1 - np.cos(sel2_JERUp_events.PuppiMET.phiJERUp - (sel2_JERUp_electrons+sel2_JERUp_jets[1]).phi))))).data))),columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel2_JERUp_d = xgb.DMatrix(sel2_JERUp_X)

            sel2_JERUp_bdtscore = bst.predict(sel2_JERUp_d)

            sel2_JERUp_weight = np.sign(sel2_JERUp_events.Generator.weight)*sel2_JERUp_pu_weight*sel2_JERUp_events.L1PreFiringWeight.Nom*sel2_JERUp_electronidsf

            output['sel2_bdtscore_binning1_JERUp'].fill(
                dataset=dataset,
                bdtscore=sel2_JERUp_bdtscore,
                weight=sel2_JERUp_weight
            )

            output['sel3_bdtscore_binning1_JERUp'].fill(
                dataset=dataset,
                bdtscore=sel2_JERUp_bdtscore,
                weight=sel2_JERUp_weight
            )

        if ak.any(basecut) and ak.any(cut2):
                
            sel2_particleindices = particleindices[cut2]
        
            sel2_events = events[cut2]
            
            sel2_jets = [sel2_events.Jet[sel2_particleindices[idx]] for idx in '0123']
            sel2_electrons = sel2_events.Electron[sel2_particleindices['4']]

            sel2_X = pandas.DataFrame(np.transpose(np.vstack((ak.to_numpy(ak.firsts(sel2_particleindices['6'])).data,ak.to_numpy(ak.firsts(sel2_particleindices['7'])).data,np.ones(len(sel2_events)),np.sign(ak.to_numpy(ak.firsts(sel2_electrons).charge).data+1),ak.to_numpy(ak.firsts(sel2_electrons).pt).data,ak.to_numpy(ak.firsts(sel2_electrons).eta).data,ak.to_numpy(ak.firsts(sel2_electrons).phi).data,ak.to_numpy(sel2_events.PuppiMET.pt),ak.to_numpy(sel2_events.PuppiMET.phi),ak.to_numpy(ak.firsts(sel2_jets[0]).pt).data,ak.to_numpy(ak.firsts(sel2_jets[1]).pt).data,ak.to_numpy(ak.firsts(sel2_jets[2]).pt).data,ak.to_numpy(ak.firsts(sel2_jets[3]).pt).data,ak.to_numpy(ak.firsts(sel2_jets[0]).eta).data,ak.to_numpy(ak.firsts(sel2_jets[1]).eta).data,ak.to_numpy(ak.firsts(sel2_jets[2]).eta).data,ak.to_numpy(ak.firsts(sel2_jets[3]).eta).data,ak.to_numpy(ak.firsts(sel2_jets[0]).phi).data,ak.to_numpy(ak.firsts(sel2_jets[1]).phi).data,ak.to_numpy(ak.firsts(sel2_jets[2]).phi).data,ak.to_numpy(ak.firsts(sel2_jets[3]).phi).data,ak.to_numpy(ak.firsts(sel2_jets[0]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_jets[1]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_jets[2]).btagDeepB).data,ak.to_numpy(ak.firsts(sel2_jets[3]).btagDeepB).data,ak.to_numpy(ak.firsts((sel2_jets[0]+sel2_jets[1]).mass)).data,ak.to_numpy(ak.firsts((sel2_jets[2]+sel2_jets[3]).mass)).data,ak.to_numpy(ak.firsts(sel2_jets[2]).eta - ak.firsts(sel2_jets[3]).eta).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_electrons+sel2_jets[0]).pt*sel2_events.PuppiMET.pt*(1 - np.cos(sel2_events.PuppiMET.phi - (sel2_electrons+sel2_jets[0]).phi))))).data,ak.to_numpy(ak.firsts(np.sqrt(2*(sel2_electrons+sel2_jets[1]).pt*sel2_events.PuppiMET.pt*(1 - np.cos(sel2_events.PuppiMET.phi - (sel2_electrons+sel2_jets[1]).phi))))).data))),columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbsjet1pt','vbsjet2pt','higgsjet1eta','higgsjet2eta','vbsjet1eta','vbsjet2eta','higgsjet1phi','higgsjet2phi','vbsjet1phi','vbsjet2phi','higgsjet1btag','higgsjet2btag','vbsjet1btag','vbsjet2btag','higgsdijetmass','vbsdijetmass','vbsdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel2_d = xgb.DMatrix(sel2_X)

            sel2_bdtscore = bst.predict(sel2_d)    

            if dataset == 'Data':
                sel2_weight = np.ones(len(sel2_events))
            else:
                sel2_pu_weight = evaluator['pileup'](sel2_events.Pileup.nTrueInt)
                sel2_puUp_weight = evaluator['pileup_up'](sel2_events.Pileup.nTrueInt)
                sel2_electronidsf = ak.firsts(evaluator['electronidsf'](sel2_electrons.eta, sel2_electrons.pt))
                sel2_electronidsfUp = ak.firsts(evaluator['electronidsfunc'](sel2_electrons.eta, sel2_electrons.pt))+sel2_electronidsf

                sel2_weight = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf
                sel2_weight_pileupUp = np.sign(sel2_events.Generator.weight)*sel2_puUp_weight*sel2_events.L1PreFiringWeight.Up*sel2_electronidsf
                sel2_weight_prefireUp = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Up*sel2_electronidsf
                sel2_weight_electronidsfUp = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsfUp

            output['sel2_bdtscore_binning1'].fill(
                dataset=dataset,
                bdtscore = sel2_bdtscore,
                weight=sel2_weight
            )

            if dataset != 'Data':

                output['sel2_bdtscore_binning1_pileupUp'].fill(
                    dataset=dataset,
                    bdtscore = sel2_bdtscore,
                    weight=sel2_weight_pileupUp
                )

                output['sel2_bdtscore_binning1_prefireUp'].fill(
                    dataset=dataset,
                    bdtscore = sel2_bdtscore,
                    weight=sel2_weight_prefireUp
                )

                output['sel2_bdtscore_binning1_electronidsfUp'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight_electronidsfUp
                )

            output['sel2_bdtscore_binning2'].fill(
                dataset=dataset,
                bdtscore = sel2_bdtscore,
                weight=sel2_weight
            )

            output['sel2_bdtscore_binning3'].fill(
                dataset=dataset,
                bdtscore=sel2_bdtscore,
                weight=sel2_weight
            )

            output['sel2_higgsdijetmass_binning1'].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel2_jets[0]+sel2_jets[1]).mass),
                weight=sel2_weight
            )

            output['sel2_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel2_jets[0]+sel2_jets[1]).pt),
                weight=sel2_weight
            )
        
            output['sel2_vbsdijetmass_binning1'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel2_jets[2]+sel2_jets[3]).mass),
                weight=sel2_weight
            )

            output['sel2_vbsdijetmass_binning2'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel2_jets[2]+sel2_jets[3]).mass),
                weight=sel2_weight
            )

            output['sel2_vbsdijetabsdeta_binning1'].fill(
                dataset=dataset,
                vbsdijetabsdeta=ak.firsts(sel2_jets[2]).eta - ak.firsts(sel2_jets[3]).eta,
                weight=sel2_weight
            )

            output['sel2_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel2_electrons.pt),
                weight=sel2_weight
            )

            output['sel2_met_binning1'].fill(
                dataset=dataset,
                met=sel2_events.PuppiMET.pt,
                weight=sel2_weight
            )

            output['sel3_bdtscore_binning1'].fill(
                dataset=dataset,
                bdtscore = sel2_bdtscore,
                weight=sel2_weight
            )

            if dataset != 'Data':

                output['sel3_bdtscore_binning1_pileupUp'].fill(
                    dataset=dataset,
                    bdtscore = sel2_bdtscore,
                    weight=sel2_weight_pileupUp
                )

                output['sel3_bdtscore_binning1_prefireUp'].fill(
                    dataset=dataset,
                    bdtscore = sel2_bdtscore,
                    weight=sel2_weight_prefireUp
                )
                
                output['sel3_bdtscore_binning1_electronidsfUp'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight_electronidsfUp
                )
                
                output['sel3_bdtscore_binning1_muonidsfUp'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight
                )
                
                output['sel3_bdtscore_binning1_muonisosfUp'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight
                )
                
                output['sel3_bdtscore_binning1_muonhltsfUp'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight
                )

            output['sel3_bdtscore_binning2'].fill(
                dataset=dataset,
                bdtscore = sel2_bdtscore,
                weight=sel2_weight
            )

            output['sel3_bdtscore_binning3'].fill(
                dataset=dataset,
                bdtscore=sel2_bdtscore,
                weight=sel2_weight
            )

            output['sel3_higgsdijetmass_binning1'].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel2_jets[0]+sel2_jets[1]).mass),
                weight=sel2_weight
            )

            output['sel3_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel2_jets[0]+sel2_jets[1]).pt),
                weight=sel2_weight
            )
        
            output['sel3_vbsdijetmass_binning1'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel2_jets[2]+sel2_jets[3]).mass),
                weight=sel2_weight
            )

            output['sel3_vbsdijetmass_binning2'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel2_jets[2]+sel2_jets[3]).mass),
                weight=sel2_weight
            )

            output['sel3_vbsdijetabsdeta_binning1'].fill(
                dataset=dataset,
                vbsdijetabsdeta=ak.firsts(sel2_jets[2]).eta - ak.firsts(sel2_jets[3]).eta,
                weight=sel2_weight
            )

            output['sel3_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel2_electrons.pt),
                weight=sel2_weight
            )

            output['sel3_met_binning1'].fill(
                dataset=dataset,
                met=sel2_events.PuppiMET.pt,
                weight=sel2_weight
            )

        if ak.any(basecut) and ak.any(cut4):
                
            sel4_particleindices = particleindices[cut4]
        
            sel4_events = events[cut4]
            
            sel4_jets = [sel4_events.Jet[sel4_particleindices[idx]] for idx in '0123']
            sel4_muons = sel4_events.Muon[sel4_particleindices['4']]
            
            if dataset == 'Data':
                sel4_weight = np.ones(len(sel4_events))
            else:
                sel4_weight = np.sign(sel4_events.Generator.weight)

            output['sel4_higgsdijetmass_binning1'].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel4_jets[0]+sel4_jets[1]).mass),
                weight=sel4_weight
            )

            output['sel4_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel4_jets[0]+sel4_jets[1]).pt),
                weight=sel4_weight
            )
        
            output['sel4_vbsdijetmass_binning1'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel4_jets[2]+sel4_jets[3]).mass),
                weight=sel4_weight
            )

            output['sel4_vbsdijetmass_binning2'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel4_jets[2]+sel4_jets[3]).mass),
                weight=sel4_weight
            )

            output['sel4_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel4_muons.pt),
                weight=sel4_weight
            )

            output['sel4_met_binning1'].fill(
                dataset=dataset,
                met=sel4_events.PuppiMET.pt,
                weight=sel4_weight
            )

            output['sel6_higgsdijetmass_binning1'].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel4_jets[0]+sel4_jets[1]).mass),
                weight=sel4_weight
            )

            output['sel6_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel4_jets[0]+sel4_jets[1]).pt),
                weight=sel4_weight
            )
        
            output['sel6_vbsdijetmass_binning1'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel4_jets[2]+sel4_jets[3]).mass),
                weight=sel4_weight
            )

            output['sel6_vbsdijetmass_binning2'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel4_jets[2]+sel4_jets[3]).mass),
                weight=sel4_weight
            )

            output['sel6_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel4_muons.pt),
                weight=sel4_weight
            )

            output['sel6_met_binning1'].fill(
                dataset=dataset,
                met=sel4_events.PuppiMET.pt,
                weight=sel4_weight
            )

        if ak.any(basecut) and ak.any(cut5):
                
            sel5_particleindices = particleindices[cut5]
        
            sel5_events = events[cut5]
            
            sel5_jets = [sel5_events.Jet[sel5_particleindices[idx]] for idx in '0123']
            sel5_electrons = sel5_events.Electron[sel5_particleindices['4']]
            
            if dataset == 'Data':
                sel5_weight = np.ones(len(sel5_events))
            else:
                sel5_weight = np.sign(sel5_events.Generator.weight)

            output['sel5_higgsdijetmass_binning1'].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel5_jets[0]+sel5_jets[1]).mass),
                weight=sel5_weight
            )

            output['sel5_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel5_jets[0]+sel5_jets[1]).pt),
                weight=sel5_weight
            )
        
            output['sel5_vbsdijetmass_binning1'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel5_jets[2]+sel5_jets[3]).mass),
                weight=sel5_weight
            )

            output['sel5_vbsdijetmass_binning2'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel5_jets[2]+sel5_jets[3]).mass),
                weight=sel5_weight
            )

            output['sel5_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel5_electrons.pt),
                weight=sel5_weight
            )

            output['sel5_met_binning1'].fill(
                dataset=dataset,
                met=sel5_events.PuppiMET.pt,
                weight=sel5_weight
            )

            output['sel6_higgsdijetmass_binning1'].fill(
                dataset=dataset,
                higgsdijetmass=ak.firsts((sel5_jets[0]+sel5_jets[1]).mass),
                weight=sel5_weight
            )

            output['sel6_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=ak.firsts((sel5_jets[0]+sel5_jets[1]).pt),
                weight=sel5_weight
            )
        
            output['sel6_vbsdijetmass_binning1'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel5_jets[2]+sel5_jets[3]).mass),
                weight=sel5_weight
            )

            output['sel6_vbsdijetmass_binning2'].fill(
                dataset=dataset,
                vbsdijetmass=ak.firsts((sel5_jets[2]+sel5_jets[3]).mass),
                weight=sel5_weight
            )

            output['sel6_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=ak.firsts(sel5_electrons.pt),
                weight=sel5_weight
            )

            output['sel6_met_binning1'].fill(
                dataset=dataset,
                met=sel5_events.PuppiMET.pt,
                weight=sel5_weight
            )

        return output
        
    def postprocess(self, accumulator):
        return accumulator

if year == '2016':
    filelists = {
        'singlemuon' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/singlemuon.txt',
        'singleelectron' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/singleelectron.txt',
        'w' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/w.txt',
        'ww' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/ww.txt',
        'ewkwhjj': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/ewkwhjj.txt',
        'qcdwhjj': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/qcdwhjj.txt',
        'ttsemi': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/ttsemi.txt',
        'tthad': '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/tthad.txt'
    }
elif year == '2017':
    filelists = {
        'singlemuon' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2017/singlemuon.txt'
    }
elif year == '2018':
    filelists = {
        'singlemuon' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/singlemuon.txt'
    }
else:
    assert(0)

samples = {}

for filelist in filelists:
    f = open(filelists[filelist])
    samples[filelist] = f.read().rstrip('\n').split('\n')
    
result = processor.run_uproot_job(
    samples,
    'Events',
    EwkwhjjProcessor(),
    processor.futures_executor,
    {'schema': NanoAODSchema, 'workers': args.nprocesses},
    chunksize=10000000,
)

for key in result['nevents'].keys():
    print('result[\'nevents\'][\'{}\'] = {}'.format(key,result['nevents'][key]))

from coffea.util import save

save(result,'outfile_{}'.format(year))
