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

        output = self.accumulator.identity()

        dataset = events.metadata['dataset']

        if dataset not in ['singleelectron','singlemuon','egamma']:
            output['sumw'][dataset] += ak.sum(np.sign(events.Generator.weight))
        output['nevents'][dataset] += len(events)

        if dataset in ['singleelectron','singlemuon','egamma']:
            events = events[lumimask(events.run,events.luminosityBlock)]

        events = events[(events.PuppiMET.pt > 30) | (events.PuppiMET.ptJERUp > 30) | (events.PuppiMET.ptJESUp > 30)]    

        if year == "2016":
            if dataset == 'singlemuon':
                events = events[events.HLT.IsoTkMu24 | events.HLT.IsoMu24]
            elif dataset == 'singleelectron':
                events = events[vents.HLT.IsoTkMu24 | events.HLT.IsoMu24 | events.HLT.Ele27_WPTight_Gsf]
            else:    
                events = events[events.HLT.IsoTkMu24 | events.HLT.IsoMu24 | events.HLT.Ele27_WPTight_Gsf]
        elif year == "2017":
            if dataset == 'singlemuon':
                events = events[events.HLT.IsoMu27]
            elif dataset == 'singleelectron':
                events = events[events.HLT.Ele32_WPTight_Gsf_L1DoubleEG]    
            else:
                events = events[events.HLT.IsoMu27 | events.HLT.Ele32_WPTight_Gsf_L1DoubleEG]        
        elif year == "2018":
            if dataset == 'singlemuon':
                events = events[events.HLT.IsoMu24]
            elif dataset == 'egamma':    
                events = events[events.HLT.Ele32_WPTight_Gsf]
            else:
                events = events[events.HLT.IsoMu24 |events.HLT.Ele32_WPTight_Gsf]

        events = events[(ak.num(events.Jet) > 3) | ((ak.num(events.Jet) > 1) & (ak.num(events.FatJet) > 0))]

        events = events[(ak.num(events.Electron) > 0) | (ak.num(events.Muon) > 0)]

        tight_muons = events.Muon[events.Muon.tightId & (events.Muon.pfRelIso04_all < 0.15) & (events.Muon.pt > 26) & (abs(events.Muon.eta) < 2.4)]

        loose_not_tight_muons = events.Muon[events.Muon.tightId & (events.Muon.pfRelIso04_all < 0.4) & (events.Muon.pfRelIso04_all > 0.15) & (events.Muon.pt > 20) & (abs(events.Muon.eta) < 2.4)]

        tight_electrons = events.Electron[(events.Electron.pt > 30) & (events.Electron.cutBased >= 3) & (events.Electron.eta + events.Electron.deltaEtaSC < 2.5) & ((abs(events.Electron.dz) < 0.1) & (abs(events.Electron.dxy) < 0.05) & (events.Electron.eta + events.Electron.deltaEtaSC < 1.479)) | ((abs(events.Electron.dz) < 0.2) & (abs(events.Electron.dxy) < 0.1) & (events.Electron.eta + events.Electron.deltaEtaSC > 1.479))]

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
        jet_pt_jesup = corrected_jets.JES_jes.up.pt
        jet_pt_jerup = corrected_jets.JER.up.pt

        corrected_jets = ak.zip({
            "pt": corrected_jets.pt,
            "eta": corrected_jets.eta,
            "phi": corrected_jets.phi,
            "mass": corrected_jets.mass,
            "charge": np.ones(len(corrected_jets.pt)),
            "btagDeepB": corrected_jets.btagDeepB
        }, with_name="PtEtaPhiMCandidate")    

        fatjets = events.FatJet[(events.FatJet.pt > 250) & (abs(events.FatJet.eta) < 2.5) & (events.FatJet.msoftdrop > 50) & (events.FatJet.msoftdrop < 150)]    
        b_jets = corrected_jets[(events.Jet.cleanmask == 1) & (jet_pt > 30) & (abs(events.Jet.eta) < 2.5) & (events.Jet.btagDeepB > 0.8953)]
        vbf_jets = corrected_jets[(events.Jet.cleanmask == 1) & (jet_pt > 30) & (abs(events.Jet.eta) < 4.7) & (events.Jet.btagDeepB < 0.2217)]
        nextrajets = ak.num(events.Jet[(events.Jet.cleanmask == 1) & (jet_pt > 30) & (abs(events.Jet.eta) < 4.7)]) - 4
        nextrabjets = ak.num(events.Jet[(events.Jet.cleanmask == 1) & (jet_pt > 30) & (abs(events.Jet.eta) < 4.7) & (events.Jet.btagDeepB > 0.2217)]) - 2

        basecut_merged = (ak.num(fatjets) > 0) & (ak.num(vbf_jets) > 1) & (ak.num(tight_muons) + ak.num(tight_electrons) == 1) & (ak.num(loose_not_tight_muons) == 0) & (events.PuppiMET.pt > 30)
        events_merged = events[basecut_merged]
        fatjets_merged = fatjets[basecut_merged]
        vbf_jets_merged = vbf_jets[basecut_merged]
        tight_muons_merged = tight_muons[basecut_merged]
        tight_electrons_merged = tight_electrons[basecut_merged]
        nextrajets_merged = nextrajets[basecut_merged]
        nextrabjets_merged = nextrabjets[basecut_merged]

        basecut = (ak.num(b_jets) > 1) & (ak.num(vbf_jets) > 1) & (ak.num(tight_muons) + ak.num(tight_electrons) == 1) & (ak.num(loose_not_tight_muons) == 0) & (events.PuppiMET.pt > 30)
        events = events[basecut]
        b_jets = b_jets[basecut]
        vbf_jets = vbf_jets[basecut]
        tight_muons = tight_muons[basecut]
        tight_electrons = tight_electrons[basecut]
        nextrajets = nextrajets[basecut]
        nextrabjets = nextrabjets[basecut]

        if dataset in ['singleelectron','singlemuon','egamma']:
            dataset = 'data'

        if ak.any(basecut_merged):
            cut7 = (fatjets_merged[:,0].mass > 50) & (fatjets_merged[:,0].mass < 150) & ((vbf_jets_merged[:,0]+vbf_jets_merged[:,1]).mass > 500) & (abs(vbf_jets_merged[:,0].eta - vbf_jets_merged[:,1].eta) > 2.5) & (ak.num(tight_muons_merged) > 0)
            cut8 = (fatjets_merged[:,0].mass > 50) & (fatjets_merged[:,0].mass < 150) & ((vbf_jets_merged[:,0]+vbf_jets_merged[:,1]).mass > 500) & (abs(vbf_jets_merged[:,0].eta - vbf_jets_merged[:,1].eta) > 2.5) & (ak.num(tight_electrons_merged) > 0)
#            cut9 = cut7 | cut8

        cut1 = ((b_jets[:,0] + b_jets[:,1]).mass > 50) & ((b_jets[:,0] + b_jets[:,1]).mass < 150) & ((vbf_jets[:,0] + vbf_jets[:,1]).mass > 500) & (abs(vbf_jets[:,0].eta - vbf_jets[:,1].eta) > 2.5) & (ak.num(tight_muons) > 0)
        cut2 = ((b_jets[:,0] + b_jets[:,1]).mass > 50) & ((b_jets[:,0] + b_jets[:,1]).mass < 150) & ((vbf_jets[:,0] + vbf_jets[:,1]).mass > 500) & (abs(vbf_jets[:,0].eta - vbf_jets[:,1].eta) > 2.5) & (ak.num(tight_electrons) > 0)
#            cut3 = cut1 | cut2

        if ak.any(basecut_merged) and ak.any(cut7):

            sel7_events = events_merged[cut7]
            sel7_fatjets = fatjets_merged[cut7]
            sel7_vbf_jets = vbf_jets_merged[cut7]
            sel7_muons = tight_muons_merged[cut7][:,0]
            sel7_nextrajets = nextrajets_merged[cut7]
            sel7_nextrabjets = nextrabjets_merged[cut7]

            output["weights_merged"][dataset] += processor.column_accumulator(np.sign(ak.to_numpy(sel7_events.Generator.weight).data))

            output['variables_merged'][dataset] += processor.column_accumulator(np.transpose(np.vstack((
                ak.to_numpy(sel7_fatjets[:,0].pt),
                ak.to_numpy(sel7_fatjets[:,0].eta),
                ak.to_numpy(sel7_fatjets[:,0].phi),
                ak.to_numpy(sel7_fatjets[:,0].btagDeepB),
                ak.to_numpy(sel7_fatjets[:,0].btagHbb),
                ak.to_numpy(sel7_fatjets[:,0].msoftdrop),
                ak.to_numpy(sel7_nextrajets),
                ak.to_numpy(sel7_nextrabjets),
                np.zeros(len(sel7_events)),
                np.sign(ak.to_numpy(sel7_muons.charge)+1),
                ak.to_numpy(sel7_muons.pt),
                ak.to_numpy(sel7_muons.eta),
                ak.to_numpy(sel7_muons.phi),
                ak.to_numpy(sel7_events.PuppiMET.pt),
                ak.to_numpy(sel7_events.PuppiMET.phi),
                ak.to_numpy(sel7_vbf_jets[:,0].pt),
                ak.to_numpy(sel7_vbf_jets[:,1].pt),
                ak.to_numpy(sel7_vbf_jets[:,0].eta),
                ak.to_numpy(sel7_vbf_jets[:,1].eta),
                ak.to_numpy(sel7_vbf_jets[:,0].phi),
                ak.to_numpy(sel7_vbf_jets[:,1].phi),
                ak.to_numpy(sel7_vbf_jets[:,0].btagDeepB),
                ak.to_numpy(sel7_vbf_jets[:,1].btagDeepB),
                ak.to_numpy((sel7_vbf_jets[:,0]+sel7_vbf_jets[:,1]).mass),
                ak.to_numpy(sel7_vbf_jets[:,0].eta - sel7_vbf_jets[:,1].eta),
                ak.to_numpy(np.sqrt(2*(sel7_muons+sel7_vbf_jets[:,0]).pt*sel7_events.PuppiMET.pt*(1 - np.cos(sel7_events.PuppiMET.phi - (sel7_muons+sel7_vbf_jets[:,0]).phi)))),
                ak.to_numpy(np.sqrt(2*(sel7_muons+sel7_vbf_jets[:,1]).pt*sel7_events.PuppiMET.pt*(1 - np.cos(sel7_events.PuppiMET.phi - (sel7_muons+sel7_vbf_jets[:,1]).phi))))))))

            sel7_muonidsf = evaluator['muonidsf'](abs(sel7_muons.eta), sel7_muons.pt)
            sel7_muonisosf = evaluator['muonisosf'](abs(sel7_muons.eta), sel7_muons.pt)
            sel7_muonhltsf = evaluator['muonhltsf'](abs(sel7_muons.eta), sel7_muons.pt)
            sel7_weight = np.sign(sel7_events.Generator.weight)*sel7_events.L1PreFiringWeight.Nom*sel7_muonidsf*sel7_muonisosf*sel7_muonhltsf

        if ak.any(basecut_merged) and ak.any(cut8):

            sel8_events = events_merged[cut8]
            sel8_fatjets = fatjets_merged[cut8]
            sel8_vbf_jets = vbf_jets_merged[cut8]
            sel8_electrons = tight_electrons_merged[cut8][:,0]
            sel8_nextrajets = nextrajets_merged[cut8]
            sel8_nextrabjets = nextrabjets_merged[cut8]

            output["weights_merged"][dataset] += processor.column_accumulator(np.sign(ak.to_numpy(sel8_events.Generator.weight).data))

            output['variables_merged'][dataset] += processor.column_accumulator(np.transpose(np.vstack((
                ak.to_numpy(sel8_fatjets[:,0].pt),
                ak.to_numpy(sel8_fatjets[:,0].eta),
                ak.to_numpy(sel8_fatjets[:,0].phi),
                ak.to_numpy(sel8_fatjets[:,0].btagDeepB),
                ak.to_numpy(sel8_fatjets[:,0].btagHbb),
                ak.to_numpy(sel8_fatjets[:,0].msoftdrop),
                ak.to_numpy(sel8_nextrajets),
                ak.to_numpy(sel8_nextrabjets),
                np.ones(len(sel8_events)),
                np.sign(ak.to_numpy(sel8_electrons.charge)+1),
                ak.to_numpy(sel8_electrons.pt),
                ak.to_numpy(sel8_electrons.eta),
                ak.to_numpy(sel8_electrons.phi),
                ak.to_numpy(sel8_events.PuppiMET.pt),
                ak.to_numpy(sel8_events.PuppiMET.phi),
                ak.to_numpy(sel8_vbf_jets[:,0].pt),
                ak.to_numpy(sel8_vbf_jets[:,1].pt),
                ak.to_numpy(sel8_vbf_jets[:,0].eta),
                ak.to_numpy(sel8_vbf_jets[:,1].eta),
                ak.to_numpy(sel8_vbf_jets[:,0].phi),
                ak.to_numpy(sel8_vbf_jets[:,1].phi),
                ak.to_numpy(sel8_vbf_jets[:,0].btagDeepB),
                ak.to_numpy(sel8_vbf_jets[:,1].btagDeepB),
                ak.to_numpy((sel8_vbf_jets[:,0]+sel8_vbf_jets[:,1]).mass),
                ak.to_numpy(sel8_vbf_jets[:,0].eta - sel8_vbf_jets[:,1].eta),
                ak.to_numpy(np.sqrt(2*(sel8_electrons+sel8_vbf_jets[:,0]).pt*sel8_events.PuppiMET.pt*(1 - np.cos(sel8_events.PuppiMET.phi - (sel8_electrons+sel8_vbf_jets[:,0]).phi)))),
                ak.to_numpy(np.sqrt(2*(sel8_electrons+sel8_vbf_jets[:,1]).pt*sel8_events.PuppiMET.pt*(1 - np.cos(sel8_events.PuppiMET.phi - (sel8_electrons+sel8_vbf_jets[:,1]).phi))))))))

            sel8_electronidsf = evaluator['electronidsf'](sel8_electrons.eta, sel8_electrons.pt)
            sel8_electronrecosf = evaluator['electronrecosf'](sel8_electrons.eta, sel8_electrons.pt)
            sel8_weight = np.sign(sel8_events.Generator.weight)*sel8_events.L1PreFiringWeight.Nom*sel8_electronidsf*sel8_electronrecosf

        if ak.any(basecut) and ak.any(cut1):

            sel1_events = events[cut1]
            
            sel1_b_jets = b_jets[cut1]

            sel1_vbf_jets = vbf_jets[cut1]

            sel1_muons = tight_muons[cut1][:,0]

            sel1_nextrajets = nextrajets[cut1]

            sel1_nextrabjets = nextrabjets[cut1]

            output["weights"][dataset] += processor.column_accumulator(np.sign(ak.to_numpy(sel1_events.Generator.weight).data))

            output['variables'][dataset] += processor.column_accumulator(np.transpose(np.vstack((
                ak.to_numpy(sel1_nextrajets),
                ak.to_numpy(sel1_nextrabjets),
                np.zeros(len(sel1_events)),
                np.sign(ak.to_numpy(sel1_muons.charge)+1),
                ak.to_numpy(sel1_muons.pt),
                ak.to_numpy(sel1_muons.eta),
                ak.to_numpy(sel1_muons.phi),
                ak.to_numpy(sel1_events.PuppiMET.pt),
                ak.to_numpy(sel1_events.PuppiMET.phi),
                ak.to_numpy(sel1_b_jets[:,0].pt),
                ak.to_numpy(sel1_b_jets[:,1].pt),
                ak.to_numpy(sel1_vbf_jets[:,0].pt),
                ak.to_numpy(sel1_vbf_jets[:,1].pt),
                ak.to_numpy(sel1_b_jets[:,0].eta),
                ak.to_numpy(sel1_b_jets[:,1].eta),
                ak.to_numpy(sel1_vbf_jets[:,0].eta),
                ak.to_numpy(sel1_vbf_jets[:,1].eta),
                ak.to_numpy(sel1_b_jets[:,0].phi),
                ak.to_numpy(sel1_b_jets[:,1].phi),
                ak.to_numpy(sel1_vbf_jets[:,0].phi),
                ak.to_numpy(sel1_vbf_jets[:,1].phi),
                ak.to_numpy(sel1_b_jets[:,0].btagDeepB),
                ak.to_numpy(sel1_b_jets[:,1].btagDeepB),
                ak.to_numpy(sel1_vbf_jets[:,0].btagDeepB),
                ak.to_numpy(sel1_vbf_jets[:,1].btagDeepB),
                ak.to_numpy((sel1_b_jets[:,0]+sel1_b_jets[:,1]).mass),
                ak.to_numpy((sel1_vbf_jets[:,0]+sel1_vbf_jets[:,1]).mass),
                ak.to_numpy(sel1_vbf_jets[:,0].eta - sel1_vbf_jets[:,1].eta),
                ak.to_numpy(np.sqrt(2*(sel1_muons+sel1_b_jets[:,0]).pt*sel1_events.PuppiMET.pt*(1 - np.cos(sel1_events.PuppiMET.phi - (sel1_muons+sel1_b_jets[:,0]).phi)))),ak.to_numpy(np.sqrt(2*(sel1_muons+sel1_b_jets[:,1]).pt*sel1_events.PuppiMET.pt*(1 - np.cos(sel1_events.PuppiMET.phi - (sel1_muons+sel1_b_jets[:,1]).phi))))))))

            sel1_pu_weight = evaluator['pileup'](sel1_events.Pileup.nTrueInt)
            sel1_muonidsf = evaluator['muonidsf'](abs(sel1_muons.eta), sel1_muons.pt)
            sel1_muonisosf = evaluator['muonisosf'](abs(sel1_muons.eta), sel1_muons.pt)
            sel1_muonhltsf = evaluator['muonhltsf'](abs(sel1_muons.eta), sel1_muons.pt)
            sel1_weight = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf

        if ak.any(basecut) and ak.any(cut2):
       
            sel2_events = events[cut2]
            sel2_b_jets = b_jets[cut2]
            sel2_vbf_jets = vbf_jets[cut2]
            sel2_electrons = tight_electrons[cut2][:,0]
            sel2_nextrajets = nextrajets[cut2]
            sel2_nextrabjets = nextrabjets[cut2]
         
            output["weights"][dataset] += processor.column_accumulator(np.sign(ak.to_numpy(sel2_events.Generator.weight).data))

            output['variables'][dataset] += processor.column_accumulator(np.transpose(np.vstack((
                ak.to_numpy(sel2_nextrajets),
                ak.to_numpy(sel2_nextrabjets),
                np.ones(len(sel2_events)),
                np.sign(ak.to_numpy(sel2_electrons.charge)+1),
                ak.to_numpy(sel2_electrons.pt),
                ak.to_numpy(sel2_electrons.eta),
                ak.to_numpy(sel2_electrons.phi),
                ak.to_numpy(sel2_events.PuppiMET.pt),
                ak.to_numpy(sel2_events.PuppiMET.phi),
                ak.to_numpy(sel2_b_jets[:,0].pt),
                ak.to_numpy(sel2_b_jets[:,1].pt),
                ak.to_numpy(sel2_vbf_jets[:,0].pt),
                ak.to_numpy(sel2_vbf_jets[:,1].pt),
                ak.to_numpy(sel2_b_jets[:,0].eta),
                ak.to_numpy(sel2_b_jets[:,1].eta),
                ak.to_numpy(sel2_vbf_jets[:,0].eta),
                ak.to_numpy(sel2_vbf_jets[:,1].eta),
                ak.to_numpy(sel2_b_jets[:,0].phi),
                ak.to_numpy(sel2_b_jets[:,1].phi),
                ak.to_numpy(sel2_vbf_jets[:,0].phi),
                ak.to_numpy(sel2_vbf_jets[:,1].phi),
                ak.to_numpy(sel2_b_jets[:,0].btagDeepB),
                ak.to_numpy(sel2_b_jets[:,1].btagDeepB),
                ak.to_numpy(sel2_vbf_jets[:,0].btagDeepB),
                ak.to_numpy(sel2_vbf_jets[:,1].btagDeepB),
                ak.to_numpy((sel2_b_jets[:,0]+sel2_b_jets[:,1]).mass),
                ak.to_numpy((sel2_vbf_jets[:,0]+sel2_vbf_jets[:,1]).mass),
                ak.to_numpy(sel2_vbf_jets[:,0].eta - sel2_vbf_jets[:,1].eta),
                ak.to_numpy(np.sqrt(2*(sel2_electrons+sel2_b_jets[:,0]).pt*sel2_events.PuppiMET.pt*(1 - np.cos(sel2_events.PuppiMET.phi - (sel2_electrons+sel2_b_jets[:,0]).phi)))),ak.to_numpy(np.sqrt(2*(sel2_electrons+sel2_b_jets[:,1]).pt*sel2_events.PuppiMET.pt*(1 - np.cos(sel2_events.PuppiMET.phi - (sel2_electrons+sel2_b_jets[:,1]).phi))))))))

            sel2_pu_weight = evaluator['pileup'](sel2_events.Pileup.nTrueInt)

            sel2_electronidsf = evaluator['electronidsf'](sel2_electrons.eta, sel2_electrons.pt)
            sel2_electronrecosf = evaluator['electronrecosf'](sel2_electrons.eta, sel2_electrons.pt)

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
        'stoptchan' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/stoptchan.txt',
        'santitoptchan' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/santitoptchan.txt',
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

save(result,"training_data_{}".format(year))

"""

save(pandas.concat(
    [
        pandas.concat(
            [
                pandas.DataFrame(result["variables"][key].value,columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','met','metphi','higgsjet1pt','higgsjet2pt','vbfjet1pt','vbfjet2pt','higgsjet1eta','higgsjet2eta','vbfjet1eta','vbfjet2eta','higgsjet1phi','higgsjet2phi','vbfjet1phi','vbfjet2phi','higgsjet1btag','higgsjet2btag','vbfjet1btag','vbfjet2btag','higgsdijetmass','vbfdijetmass','vbfdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt']),
                pandas.DataFrame(result["weights"][key].value,columns=["weight"]),
                pandas.DataFrame(len(result["variables"][key].value)*[key],columns=["label"])
         ],axis=1)
        for key in result['variables'].keys()
    ],
ignore_index=True),"ewkwhjj_training_data_{}".format(args.year))

save(pandas.concat(
    [
        pandas.concat(
            [
                pandas.DataFrame(result["variables_merged"][key].value,columns=['higgsjetpt','higgsjeteta','higgsjetphi','higgsjetbtag1','higgsjetbtag2','higgsjetsoftdropmass','nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','met','metphi','vbfjet1pt','vbfjet2pt','vbfjet1eta','vbfjet2eta','vbfjet1phi','vbfjet2phi','vbfjet1btag','vbfjet2btag','vbfdijetmass','vbfdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt']),
                pandas.DataFrame(result["weights_merged"][key].value,columns=["weight"]),
                pandas.DataFrame(len(result["variables_merged"][key].value)*[key],columns=["label"])
            ],axis=1)
        for key in result['variables'].keys()
    ],
ignore_index=True),"ewkwhjj_merged_training_data_{}".format(args.year))

"""
