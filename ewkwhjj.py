from coffea.jetmet_tools import FactorizedJetCorrector, JetCorrectionUncertainty
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory

import awkward as ak
import numpy as np

import uproot as uproot
from coffea.nanoevents import NanoEventsFactory, BaseSchema, NanoAODSchema

import awkward as ak
import coffea.hist
from coffea import processor


from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

from coffea.lumi_tools import LumiMask

import numba

import hist

import math
import numpy as np

import xgboost as xgb

import pandas

import argparse

import pprint

from coffea.lookup_tools.dense_lookup import dense_lookup

import correctionlib

cset = correctionlib.CorrectionSet.from_file("/afs/cern.ch/user/a/amlevin/ewkwhjj/data/2018/btagging.json.gz")
#print([c for c in cset])

def lighttagSF(j, syst="central"):
    # until correctionlib handles jagged data natively we have to flatten and unflatten
    j, nj = ak.flatten(j), ak.num(j)
    sf = cset["deepCSV_incl"].evaluate(syst, "M", np.array(j.hadronFlavour), np.array(abs(j.eta)), np.array(j.pt))
    return ak.unflatten(sf, nj)


def btagSF(j, syst="central"):
    # until correctionlib handles jagged data natively we have to flatten and unflatten
    j, nj = ak.flatten(j), ak.num(j)
    sf = cset["deepCSV_comb"].evaluate(syst, "T", np.array(j.hadronFlavour), np.array(abs(j.eta)), np.array(j.pt))
    return ak.unflatten(sf, nj)

def combine(eff, sf, passbtag):
    # tagged SF = SF*eff / eff = SF
    tagged_sf = ak.prod(sf[passbtag], axis=-1)
    # untagged SF = (1 - SF*eff) / (1 - eff)
    # untagged_sf = ak.prod(((1 - sf*eff) / (1 - eff))[~passbtag], axis=-1)
    untagged_sf = ak.prod((((1 - sf*eff)[~passbtag]) / ((1 - eff)[~passbtag])), axis=-1)
    return tagged_sf * untagged_sf

parser = argparse.ArgumentParser()

parser.add_argument('--year',dest='year',default='2016')
parser.add_argument('--nproc',dest='nproc',type=int,default='10')
parser.add_argument('--nreweights',dest='nreweights',type=int,default='10')

args = parser.parse_args()

pprint.pprint(vars(args))

assert(args.year == '2016' or args.year == '2017' or args.year == '2018')

year = args.year

bst = xgb.Booster({'nthread': 1})

bst.load_model('/afs/cern.ch/user/a/amlevin/ewkwhjj/models/resolved.model')

bst_merged = xgb.Booster({'nthread': 1})

bst_merged.load_model('/afs/cern.ch/user/a/amlevin/ewkwhjj/models/merged.model')

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
            'sel1_bdtscore_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_pileupup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_pileupdown': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_prefireup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_electronidsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_electronrecosfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_muonidsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_muonisosfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_muonhltsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_bcbtagsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_lightbtagsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_jesup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning1_jerup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel1_bdtscore_binning2': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 40, -0.5, 1.5),
            ),
            'sel1_bdtscore_binning3': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 21, 0, 1.05),
            ),
            'sel1_higgsdijetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsdijetmass', 'Higgs dijet mass [GeV]', 75, 0, 300),
            ),
            'sel1_higgsdijetpt_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsdijetpt', 'Higgs dijet pt [GeV]', 20, 0, 200),
            ),
            'sel1_vbfdijetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetmass', 'VBF dijet mass [GeV]', 19, 100, 2000),
            ),
            'sel1_vbfdijetmass_binning2': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetmass', 'VBF dijet mass [GeV]', 39, 100, 4000),
            ),
            'sel1_vbfdijetabsdeta_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetabsdeta', 'VBF dijet $\Delta \eta$', 55, 2.5, 8),
            ),
            'sel1_leptonpt_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('leptonpt', 'Lepton pt [GeV]', 19, 20, 200),
            ),
            'sel1_leptonabseta_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('leptonabseta', 'Lepton |eta|', 25, 0, 2.5),
            ),
            'sel1_met_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('met', 'MET [GeV]', 20, 0, 200),
            ),
            'sel2_bdtscore_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_pileupup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_pileupdown': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_prefireup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_electronidsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_electronrecosfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_bcbtagsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_lightbtagsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_jesup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning1_jerup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel2_bdtscore_binning2': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 40, -0.5, 1.5),
            ),
            'sel2_bdtscore_binning3': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 21, 0, 1.05),
            ),
            'sel2_higgsdijetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsdijetmass', 'Higgs dijet mass [GeV]', 75, 0, 300),
            ),
            'sel2_higgsdijetpt_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsdijetpt', 'Higgs dijet pt [GeV]', 20, 0, 200),
            ),
            'sel2_vbfdijetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetmass', 'VBF dijet mass [GeV]', 19, 100, 2000),
            ),
            'sel2_vbfdijetmass_binning2': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetmass', 'VBF dijet mass [GeV]', 39, 100, 4000),
            ),
            'sel2_vbfdijetabsdeta_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetabsdeta', 'VBF dijet $\Delta \eta$', 55, 2.5, 8),
            ),
            'sel2_leptonpt_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('leptonpt', 'Lepton pt [GeV]', 19, 20, 200),
            ),
            'sel2_leptonabseta_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('leptonabseta', 'Lepton |eta|', 25, 0, 2.5),
            ),
            'sel2_met_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('met', 'MET [GeV]', 20, 0, 200),
            ),
            'sel3_bdtscore_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_pileupup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_pileupdown': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_prefireup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_electronidsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_electronrecosfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_muonidsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_muonisosfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_muonhltsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_bcbtagsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_lightbtagsfup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_jesup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning1_jerup': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1),
            ),
            'sel3_bdtscore_binning2': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 40, -0.5, 1.5),
            ),
            'sel3_bdtscore_binning3': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 21, 0, 1.05),
            ),
            'sel3_higgsdijetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsdijetmass', 'Higgs dijet mass [GeV]', 75, 0, 300),
            ),
            'sel3_higgsdijetpt_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsdijetpt', 'Higgs dijet pt [GeV]', 20, 0, 200),
            ),
            'sel3_vbfdijetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetmass', 'VBF dijet mass [GeV]', 19, 100, 2000),
            ),
            'sel3_vbfdijetmass_binning2': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetmass', 'VBF dijet mass [GeV]', 39, 100, 4000),
            ),
            'sel3_vbfdijetabsdeta_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetabsdeta', 'VBF dijet $\Delta \eta$', 55, 2.5, 8),
            ),
            'sel3_leptonpt_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('leptonpt', 'Lepton pt [GeV]', 19, 20, 200),
            ),
            'sel3_leptonabseta_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('leptonabseta', 'Lepton |eta|', 25, 0, 2.5),
            ),
            'sel3_met_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('met', 'MET [GeV]', 20, 0, 200),
            ),
            'sel4_higgsdijetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsdijetmass', 'Higgs dijet mass [GeV]', 75, 0, 300),
            ),
            'sel4_higgsdijetpt_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsdijetpt', 'Higgs dijet pt [GeV]', 20, 0, 200),
            ),
            'sel4_vbfdijetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetmass', 'VBF dijet mass [GeV]', 19, 100, 2000),
            ),
            'sel4_vbfdijetmass_binning2': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetmass', 'VBF dijet mass [GeV]', 4, 100, 500),
            ),
            'sel4_leptonpt_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('leptonpt', 'Lepton pt [GeV]', 19, 20, 200),
            ),
            'sel4_leptonabseta_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('leptonabseta', 'Lepton |eta|', 25, 0, 2.5),
            ),
            'sel4_met_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('met', 'MET [GeV]', 20, 0, 200),
            ),
            'sel5_higgsdijetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsdijetmass', 'Higgs dijet mass [GeV]', 75, 0, 300),
            ),
            'sel5_higgsdijetpt_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsdijetpt', 'Higgs dijet pt [GeV]', 20, 0, 200),
            ),
            'sel5_vbfdijetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetmass', 'VBF dijet mass [GeV]', 19, 100, 2000),
            ),
            'sel5_vbfdijetmass_binning2': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetmass', 'VBF dijet mass [GeV]', 4, 100, 500),
            ),
            'sel5_leptonpt_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('leptonpt', 'Lepton pt [GeV]', 19, 20, 200),
            ),
            'sel5_leptonabseta_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('leptonabseta', 'Lepton |eta|', 25, 0, 2.5),
            ),
            'sel5_met_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('met', 'MET [GeV]', 20, 0, 200),
            ),
            'sel6_higgsdijetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsdijetmass', 'Higgs dijet mass [GeV]', 75, 0, 300),
            ),
            'sel6_higgsdijetpt_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsdijetpt', 'Higgs dijet pt [GeV]', 20, 0, 200),
            ),
            'sel6_vbfdijetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetmass', 'VBF dijet mass [GeV]', 19, 100, 2000),
            ),
            'sel6_vbfdijetmass_binning2': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('vbfdijetmass', 'VBF dijet mass [GeV]', 4, 100, 500),
            ),
            'sel6_leptonpt_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('leptonpt', 'Lepton pt [GeV]', 19, 20, 200),
            ),
            'sel6_leptonabseta_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('leptonabseta', 'Lepton |eta|', 25, 0, 2.5),
            ),
            'sel6_met_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('met', 'MET [GeV]', 20, 0, 200),
            ),
            'sel7_higgsjetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsjetmass', 'Higgs jet mass [GeV]', 75, 0, 300),
            ),
            'sel7_higgsjetsoftdropmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsjetsoftdropmass', 'Higgs jet soft drop mass [GeV]', 75, 0, 300),
            ),
            'sel7_bdtscore_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1.0),
            ),
            'sel7_bdtscore_binning2': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 40, -0.5, 1.5),
            ),
            'sel7_bdtscore_binning3': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 21, 0, 1.05),
            ),
            'sel8_higgsjetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsjetmass', 'Higgs jet mass [GeV]', 75, 0, 300),
            ),
            'sel8_higgsjetsoftdropmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsjetsoftdropmass', 'Higgs jet soft drop mass [GeV]', 75, 0, 300),
            ),
            'sel8_bdtscore_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1.0),
            ),
            'sel8_bdtscore_binning2': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 40, -0.5, 1.5),
            ),
            'sel8_bdtscore_binning3': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 21, 0, 1.05),
            ),
            'sel9_higgsjetmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsjetmass', 'Higgs jet mass [GeV]', 75, 0, 300),
            ),
            'sel9_higgsjetsoftdropmass_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('higgsjetsoftdropmass', 'Higgs jet soft drop mass [GeV]', 75, 0, 300),
            ),
            'sel9_bdtscore_binning1': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 20, 0, 1.0),
            ),
            'sel9_bdtscore_binning2': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 40, -0.5, 1.5),
            ),
            'sel9_bdtscore_binning3': coffea.hist.Hist(
                'Events',
                coffea.hist.Cat('dataset', 'Dataset'),
                coffea.hist.Bin('bdtscore', 'BDT score', 21, 0, 1.05),
            ),

        })

    @property
    def accumulator(self):
        return self._accumulator

#    @numba.njit
    def process(self, events):

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
                events = events[events.HLT.IsoTkMu24 | events.HLT.IsoMu24 | events.HLT.Ele27_WPTight_Gsf]
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

        if dataset not in ['singleelectron','singlemuon','egamma']:

            phasespace_cuts = (
                (events.Jet.cleanmask == 1) 
                & (abs(events.Jet.eta) < 2.5)
                & (events.Jet.pt > 30.)
            )
            
            jets = ak.flatten(events.Jet[phasespace_cuts])

            efficiencyinfo = (
                hist.Hist.new
                .Reg(10, 30, 300, name="pt")
                .Reg(4, 0, 2.5, name="abseta")
                .IntCat([0, 4, 5], name="flavor")
                .Bool(name="passWP")
                .Double()
                .fill(
                    pt=np.minimum(jets.pt,299),
                    abseta=np.minimum(abs(jets.eta),2.49),
                    flavor=jets.hadronFlavour,
                    passWP=jets.btagDeepB > 0.8953, # UL 2018 medium WP
                )
            )

            eff = efficiencyinfo[{"passWP": True}] / efficiencyinfo[{"passWP": sum}]

            # note this seems to turn 0,4,5 into 0,1,2
            # efflookup = dense_lookup(eff.values(), [ax.edges for ax in eff.axes])

            edges = [ax.edges for ax in eff.axes]

            assert(np.all(edges[2] == np.array([0.,1.,2.,3.])))
            
            edges = edges[0:2]

            edges.append(np.array([0.,4.,5.]))

            efflookup = dense_lookup(eff.values(), edges)

            lightJets = events.Jet[phasespace_cuts & (events.Jet.hadronFlavour == 0)]
            bcJets = events.Jet[phasespace_cuts & (events.Jet.hadronFlavour > 0)]
            
            lightEff = efflookup(lightJets.pt, abs(lightJets.eta), lightJets.hadronFlavour)
            bcEff = efflookup(bcJets.pt, abs(bcJets.eta), bcJets.hadronFlavour)

            lightweight = combine(
                lightEff,
                lighttagSF(lightJets),
                lightJets.btagDeepB > 0.8953,
            )
            bcweight = combine(
                bcEff,
                btagSF(bcJets),
                bcJets.btagDeepB > 0.8953,
            )
            btagsf = lightweight * bcweight

            btagsf_lightup = combine(
                lightEff,
                lighttagSF(lightJets, "up"),
                lightJets.btagDeepB > 0.8953,
            ) * bcweight

            btagsf_bcup = combine(
                bcEff,
                btagSF(bcJets, "up"),
                bcJets.btagDeepB > 0.8953,
            ) * lightweight

        if dataset not in ['singleelectron','singlemuon','egamma']:

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

        else:    

            corrected_jets = events.Jet

            jet_pt = events.Jet.pt
            jet_pt_jesup = events.Jet.pt
            jet_pt_jerup = events.Jet.pt


        if dataset not in ['singleelectron','singlemuon','egamma']:
            corrected_jets_jesup = ak.zip({
                "pt": corrected_jets.JES_jes.up.pt,
                "eta": corrected_jets.JES_jes.up.eta,
                "phi": corrected_jets.JES_jes.up.phi,
                "mass": corrected_jets.JES_jes.up.mass,
                "charge": np.ones(len(corrected_jets.pt)),
                "btagDeepB": corrected_jets.btagDeepB
            }, with_name="PtEtaPhiMCandidate")    
            
            corrected_jets_jerup = ak.zip({
                "pt": corrected_jets.JER.up.pt,
                "eta": corrected_jets.JER.up.eta,
                "phi": corrected_jets.JER.up.phi,
                "mass": corrected_jets.JER
                .up.mass,
                "charge": np.ones(len(corrected_jets.pt)),
                "btagDeepB": corrected_jets.btagDeepB
            }, with_name="PtEtaPhiMCandidate")    

        corrected_jets = ak.zip({
            "pt": corrected_jets.pt,
            "eta": corrected_jets.eta,
            "phi": corrected_jets.phi,
            "mass": corrected_jets.mass,
            "charge": np.ones(len(corrected_jets.pt)),
            "btagDeepB": corrected_jets.btagDeepB
        }, with_name="PtEtaPhiMCandidate")    

        if dataset not in ['singleelectron','singlemuon','egamma']:
            b_jets_jesup = corrected_jets_jesup[(events.Jet.cleanmask == 1) & (jet_pt_jesup > 30) & (abs(events.Jet.eta) < 2.5) & (events.Jet.btagDeepB > 0.8953)]
            vbf_jets_jesup = corrected_jets_jesup[(events.Jet.cleanmask == 1) & (jet_pt_jesup > 30) & (abs(events.Jet.eta) < 4.7) & (events.Jet.btagDeepB < 0.2217)]
            nextrajets_jesup = ak.num(events.Jet[(events.Jet.cleanmask == 1) & (jet_pt_jesup > 30) & (abs(events.Jet.eta) < 4.7)]) - 4
            nextrabjets_jesup = ak.num(events.Jet[(events.Jet.cleanmask == 1) & (jet_pt_jesup > 30) & (abs(events.Jet.eta) < 4.7) & (events.Jet.btagDeepB > 0.2217)]) - 2

            b_jets_jerup = corrected_jets_jerup[(events.Jet.cleanmask == 1) & (jet_pt_jerup > 30) & (abs(events.Jet.eta) < 2.5) & (events.Jet.btagDeepB > 0.8953)]
            vbf_jets_jerup = corrected_jets_jerup[(events.Jet.cleanmask == 1) & (jet_pt_jerup > 30) & (abs(events.Jet.eta) < 4.7) & (events.Jet.btagDeepB < 0.2217)]
            nextrajets_jerup = ak.num(events.Jet[(events.Jet.cleanmask == 1) & (jet_pt_jerup > 30) & (abs(events.Jet.eta) < 4.7)]) - 4
            nextrabjets_jerup = ak.num(events.Jet[(events.Jet.cleanmask == 1) & (jet_pt_jerup > 30) & (abs(events.Jet.eta) < 4.7) & (events.Jet.btagDeepB > 0.2217)]) - 2

            basecut_jesup = (ak.num(b_jets_jesup) > 1) & (ak.num(vbf_jets_jesup) > 1) & (ak.num(tight_muons) + ak.num(tight_electrons) == 1) & (ak.num(loose_not_tight_muons) == 0) & (events.PuppiMET.ptJESUp > 30)
            events_jesup = events[basecut_jesup]
            b_jets_jesup = b_jets_jesup[basecut_jesup]
            vbf_jets_jesup = vbf_jets_jesup[basecut_jesup]
            tight_muons_jesup = tight_muons[basecut_jesup]
            loose_not_tight_muons_jesup = loose_not_tight_muons[basecut_jesup]
            tight_electrons_jesup = tight_electrons[basecut_jesup]
            nextrajets_jesup = nextrajets_jesup[basecut_jesup]
            nextrabjets_jesup = nextrabjets_jesup[basecut_jesup]
            btagsf_jesup = btagsf[basecut_jesup]

            basecut_jerup = (ak.num(b_jets_jerup) > 1) & (ak.num(vbf_jets_jerup) > 1) & (ak.num(tight_muons) + ak.num(tight_electrons) == 1) & (ak.num(loose_not_tight_muons) == 0) & (events.PuppiMET.ptJESUp > 30)
            events_jerup = events[basecut_jerup]
            b_jets_jerup = b_jets_jerup[basecut_jerup]
            vbf_jets_jerup = vbf_jets_jerup[basecut_jerup]
            tight_muons_jerup = tight_muons[basecut_jerup]
            loose_not_tight_muons_jerup = loose_not_tight_muons[basecut_jerup]
            tight_electrons_jerup = tight_electrons[basecut_jerup]
            nextrajets_jerup = nextrajets_jerup[basecut_jerup]
            nextrabjets_jerup = nextrabjets_jerup[basecut_jerup]
            btagsf_jerup = btagsf[basecut_jerup]



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
        if dataset not in ['singleelectron','singlemuon','egamma']:
            btagsf = btagsf[basecut]
            btagsf_lightup = btagsf_lightup[basecut]
            btagsf_bcup = btagsf_bcup[basecut]


        if dataset in ['singleelectron','singlemuon','egamma']:
            dataset = 'data'

        if ak.any(basecut_merged):
            cut7 = (fatjets_merged[:,0].mass > 50) & (fatjets_merged[:,0].mass < 150) & ((vbf_jets_merged[:,0]+vbf_jets_merged[:,1]).mass > 500) & (abs(vbf_jets_merged[:,0].eta - vbf_jets_merged[:,1].eta) > 2.5) & (ak.num(tight_muons_merged) > 0)
            cut8 = (fatjets_merged[:,0].mass > 50) & (fatjets_merged[:,0].mass < 150) & ((vbf_jets_merged[:,0]+vbf_jets_merged[:,1]).mass > 500) & (abs(vbf_jets_merged[:,0].eta - vbf_jets_merged[:,1].eta) > 2.5) & (ak.num(tight_electrons_merged) > 0)
#            cut9 = cut7 | cut8

        if dataset != 'data' and ak.any(basecut_jesup):
            cut1_jesup = ((b_jets_jesup[:,0]+b_jets_jesup[:,1]).mass > 50) & ((b_jets_jesup[:,0]+b_jets_jesup[:,1]).mass < 150) & ((vbf_jets_jesup[:,0]+vbf_jets_jesup[:,1]).mass > 500) & (abs(vbf_jets_jesup[:,0].eta - vbf_jets_jesup[:,1].eta) > 2.5) & (ak.num(tight_muons_jesup) > 0)
            cut2_jesup = ((b_jets_jesup[:,0]+b_jets_jesup[:,1]).mass > 50) & ((b_jets_jesup[:,0]+b_jets_jesup[:,1]).mass < 150) & ((vbf_jets_jesup[:,0]+vbf_jets_jesup[:,1]).mass > 500) & (abs(vbf_jets_jesup[:,0].eta - vbf_jets_jesup[:,1].eta) > 2.5) & (ak.num(tight_electrons_jesup) > 0)
#            cut3_jesup = cut1_jesup | cut2_jesup

        if dataset != 'data' and ak.any(basecut_jerup):
            cut1_jerup = ((b_jets_jerup[:,0]+b_jets_jerup[:,1]).mass > 50) & ((b_jets_jerup[:,0]+b_jets_jerup[:,1]).mass < 150) & ((vbf_jets_jerup[:,0]+vbf_jets_jerup[:,1]).mass > 500) & (abs(vbf_jets_jerup[:,0].eta - vbf_jets_jerup[:,1].eta) > 2.5) & (ak.num(tight_muons_jerup) > 0)
            cut2_jerup = ((b_jets_jerup[:,0]+b_jets_jerup[:,1]).mass > 50) & ((b_jets_jerup[:,0]+b_jets_jerup[:,1]).mass < 150) & ((vbf_jets_jerup[:,0]+vbf_jets_jerup[:,1]).mass > 500) & (abs(vbf_jets_jerup[:,0].eta - vbf_jets_jerup[:,1].eta) > 2.5) & (ak.num(tight_electrons_jerup) > 0)
#            cut3_jerup = cut1_jerup | cut2_jerup

        cut1 = ((b_jets[:,0] + b_jets[:,1]).mass > 50) & ((b_jets[:,0] + b_jets[:,1]).mass < 150) & ((vbf_jets[:,0] + vbf_jets[:,1]).mass > 500) & (abs(vbf_jets[:,0].eta - vbf_jets[:,1].eta) > 2.5) & (ak.num(tight_muons) > 0)
        cut2 = ((b_jets[:,0] + b_jets[:,1]).mass > 50) & ((b_jets[:,0] + b_jets[:,1]).mass < 150) & ((vbf_jets[:,0] + vbf_jets[:,1]).mass > 500) & (abs(vbf_jets[:,0].eta - vbf_jets[:,1].eta) > 2.5) & (ak.num(tight_electrons) > 0)
#            cut3 = cut1 | cut2
        cut4 = ((b_jets[:,0] + b_jets[:,1]).mass > 50) & ((b_jets[:,0] + b_jets[:,1]).mass < 150) & ((vbf_jets[:,0] + vbf_jets[:,1]).mass > 100) & ((vbf_jets[:,0] + vbf_jets[:,1]).mass < 500) & (ak.num(tight_muons) > 0)
        cut5 = ((b_jets[:,0] + b_jets[:,1]).mass > 50) & ((b_jets[:,0] + b_jets[:,1]).mass < 150) & ((vbf_jets[:,0] + vbf_jets[:,1]).mass > 100) & ((vbf_jets[:,0] + vbf_jets[:,1]).mass < 500) & (ak.num(tight_electrons) > 0)
#            cut6 = cut4 | cut5

        if ak.any(basecut_merged) and ak.any(cut7):

            sel7_events = events_merged[cut7]
            sel7_fatjets = fatjets_merged[cut7]
            sel7_vbf_jets = vbf_jets_merged[cut7]
            sel7_muons = tight_muons_merged[cut7][:,0]
            sel7_nextrajets = nextrajets_merged[cut7]
            sel7_nextrabjets = nextrabjets_merged[cut7]

            sel7_X = pandas.DataFrame(np.transpose(np.vstack((
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

            sel7_d = xgb.DMatrix(sel7_X)

            sel7_bdtscore = bst_merged.predict(sel7_d)

            if dataset == 'data':
                sel7_weight = np.ones(len(sel7_events))
            else:    
                sel7_muonidsf = evaluator['muonidsf'](abs(sel7_muons.eta), sel7_muons.pt)
                sel7_muonisosf = evaluator['muonisosf'](abs(sel7_muons.eta), sel7_muons.pt)
                sel7_muonhltsf = evaluator['muonhltsf'](abs(sel7_muons.eta), sel7_muons.pt)
                sel7_weight = np.sign(sel7_events.Generator.weight)*sel7_events.L1PreFiringWeight.Nom*sel7_muonidsf*sel7_muonisosf*sel7_muonhltsf

            output['sel7_higgsjetmass_binning1'].fill(
                dataset=dataset,
                higgsjetmass=sel7_fatjets[:,0].mass,
                weight=sel7_weight
            )

            output['sel7_higgsjetsoftdropmass_binning1'].fill(
                dataset=dataset,
                higgsjetsoftdropmass=sel7_fatjets[:,0].msoftdrop,
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
                higgsjetmass=sel7_fatjets[:,0].mass,
                weight=sel7_weight
            )

            output['sel9_higgsjetsoftdropmass_binning1'].fill(
                dataset=dataset,
                higgsjetsoftdropmass=sel7_fatjets[:,0].msoftdrop,
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

            sel8_events = events_merged[cut8]
            sel8_fatjets = fatjets_merged[cut8]
            sel8_vbf_jets = vbf_jets_merged[cut8]
            sel8_electrons = tight_electrons_merged[cut8][:,0]
            sel8_nextrajets = nextrajets_merged[cut8]
            sel8_nextrabjets = nextrabjets_merged[cut8]

            sel8_X = pandas.DataFrame(np.transpose(np.vstack((
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

            sel8_d = xgb.DMatrix(sel8_X)

            sel8_bdtscore = bst_merged.predict(sel8_d)

            if dataset == 'data':
                sel8_weight = np.ones(len(sel8_events))
            else:    
                sel8_electronidsf = evaluator['electronidsf'](sel8_electrons.eta, sel8_electrons.pt)
                sel8_electronrecosf = evaluator['electronrecosf'](sel8_electrons.eta, sel8_electrons.pt)
                sel8_weight = np.sign(sel8_events.Generator.weight)*sel8_events.L1PreFiringWeight.Nom*sel8_electronidsf*sel8_electronrecosf

            output['sel8_higgsjetmass_binning1'].fill(
                dataset=dataset,
                higgsjetmass=sel8_fatjets[:,0].mass,
                weight=sel8_weight
            )

            output['sel8_higgsjetsoftdropmass_binning1'].fill(
                dataset=dataset,
                higgsjetsoftdropmass=sel8_fatjets[:,0].msoftdrop,
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
                higgsjetmass=sel8_fatjets[:,0].mass,
                weight=sel8_weight
            )

            output['sel9_higgsjetsoftdropmass_binning1'].fill(
                dataset=dataset,
                higgsjetsoftdropmass=sel8_fatjets[:,0].msoftdrop,
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

        if dataset != 'data' and ak.any(basecut_jesup) and ak.any(cut1_jesup):

            sel1_jesup_events = events_jesup[cut1_jesup]
            sel1_jesup_b_jets = b_jets_jesup[cut1_jesup]
            sel1_jesup_vbf_jets = vbf_jets_jesup[cut1_jesup]
            sel1_jesup_muons = tight_muons_jesup[cut1_jesup][:,0]
            sel1_jesup_nextrajets = nextrajets_jesup[cut1_jesup]
            sel1_jesup_nextrabjets = nextrabjets_jesup[cut1_jesup]

            sel1_jesup_pu_weight = evaluator['pileup'](sel1_jesup_events.Pileup.nTrueInt)
            sel1_jesup_muonidsf = evaluator['muonidsf'](abs(sel1_jesup_muons.eta), sel1_jesup_muons.pt)
            sel1_jesup_muonisosf = evaluator['muonisosf'](abs(sel1_jesup_muons.eta), sel1_jesup_muons.pt)
            sel1_jesup_muonhltsf = evaluator['muonhltsf'](abs(sel1_jesup_muons.eta), sel1_jesup_muons.pt)
            sel1_jesup_btagsf = btagsf_jesup[cut1_jesup]

            sel1_jesup_X = pandas.DataFrame(np.transpose(np.vstack((
                ak.to_numpy(sel1_jesup_nextrajets),
                ak.to_numpy(sel1_jesup_nextrabjets),
                np.zeros(len(sel1_jesup_events)),
                np.sign(ak.to_numpy(sel1_jesup_muons.charge)+1),
                ak.to_numpy(sel1_jesup_muons.pt),
                ak.to_numpy(sel1_jesup_muons.eta),
                ak.to_numpy(sel1_jesup_muons.phi),
                ak.to_numpy(sel1_jesup_events.PuppiMET.ptJESUp),
                ak.to_numpy(sel1_jesup_events.PuppiMET.phiJESUp),
                ak.to_numpy(sel1_jesup_b_jets[:,0].pt),
                ak.to_numpy(sel1_jesup_b_jets[:,1].pt),
                ak.to_numpy(sel1_jesup_vbf_jets[:,0].pt),
                ak.to_numpy(sel1_jesup_vbf_jets[:,1].pt),
                ak.to_numpy(sel1_jesup_b_jets[:,0].eta),
                ak.to_numpy(sel1_jesup_b_jets[:,1].eta),
                ak.to_numpy(sel1_jesup_vbf_jets[:,0].eta),
                ak.to_numpy(sel1_jesup_vbf_jets[:,1].eta),
                ak.to_numpy(sel1_jesup_b_jets[:,0].phi),
                ak.to_numpy(sel1_jesup_b_jets[:,1].phi),
                ak.to_numpy(sel1_jesup_vbf_jets[:,0].phi),
                ak.to_numpy(sel1_jesup_vbf_jets[:,1].phi),
                ak.to_numpy(sel1_jesup_b_jets[:,0].btagDeepB),
                ak.to_numpy(sel1_jesup_b_jets[:,1].btagDeepB),
                ak.to_numpy(sel1_jesup_vbf_jets[:,0].btagDeepB),
                ak.to_numpy(sel1_jesup_vbf_jets[:,1].btagDeepB),
                ak.to_numpy((sel1_jesup_b_jets[:,0]+sel1_jesup_b_jets[:,1]).mass),
                ak.to_numpy((sel1_jesup_vbf_jets[:,0]+sel1_jesup_vbf_jets[:,1]).mass),
                ak.to_numpy(sel1_jesup_vbf_jets[:,0].eta - sel1_jesup_vbf_jets[:,1].eta),
                ak.to_numpy(np.sqrt(2*(sel1_jesup_muons+sel1_jesup_b_jets[:,0]).pt*sel1_jesup_events.PuppiMET.ptJESUp*(1 - np.cos(sel1_jesup_events.PuppiMET.phi - (sel1_jesup_muons+sel1_jesup_b_jets[:,0]).phi)))),
                ak.to_numpy(np.sqrt(2*(sel1_jesup_muons+sel1_jesup_b_jets[:,1]).pt*sel1_jesup_events.PuppiMET.ptJESUp*(1 - np.cos(sel1_jesup_events.PuppiMET.phiJESUp - (sel1_jesup_muons+sel1_jesup_b_jets[:,1]).phi))))))),
                columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbfjet1pt','vbfjet2pt','higgsjet1eta','higgsjet2eta','vbfjet1eta','vbfjet2eta','higgsjet1phi','higgsjet2phi','vbfjet1phi','vbfjet2phi','higgsjet1btag','higgsjet2btag','vbfjet1btag','vbfjet2btag','higgsdijetmass','vbfdijetmass','vbfdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel1_jesup_d = xgb.DMatrix(sel1_jesup_X)

            sel1_jesup_bdtscore = bst.predict(sel1_jesup_d)

            if dataset == 'ewkwhjj_reweighted':
                sel1_jesup_weight = np.sign(sel1_jesup_events.Generator.weight)*sel1_jesup_pu_weight*sel1_jesup_events.L1PreFiringWeight.Nom*sel1_jesup_muonidsf*sel1_jesup_muonisosf*sel1_jesup_muonhltsf*sel1_jesup_btagsf*sel1_jesup_events.LHEReweightingWeight[:,9]
            else:    
                sel1_jesup_weight = np.sign(sel1_jesup_events.Generator.weight)*sel1_jesup_pu_weight*sel1_jesup_events.L1PreFiringWeight.Nom*sel1_jesup_muonidsf*sel1_jesup_muonisosf*sel1_jesup_muonhltsf*sel1_jesup_btagsf

            output['sel1_bdtscore_binning1_jesup'].fill(
                dataset=dataset,
                bdtscore=sel1_jesup_bdtscore,
                weight=sel1_jesup_weight
            )

            output['sel3_bdtscore_binning1_jesup'].fill(
                dataset=dataset,
                bdtscore=sel1_jesup_bdtscore,
                weight=sel1_jesup_weight
            )

        if dataset != 'data' and ak.any(basecut_jerup) and ak.any(cut1_jerup):

            sel1_jerup_events = events_jerup[cut1_jerup]
            sel1_jerup_b_jets = b_jets_jerup[cut1_jerup]
            sel1_jerup_vbf_jets = vbf_jets_jerup[cut1_jerup]
            sel1_jerup_muons = tight_muons_jerup[cut1_jerup][:,0]
            sel1_jerup_nextrajets = nextrajets_jerup[cut1_jerup]
            sel1_jerup_nextrabjets = nextrabjets_jerup[cut1_jerup]

            sel1_jerup_pu_weight = evaluator['pileup'](sel1_jerup_events.Pileup.nTrueInt)
            sel1_jerup_muonidsf = evaluator['muonidsf'](abs(sel1_jerup_muons.eta), sel1_jerup_muons.pt)
            sel1_jerup_muonisosf = evaluator['muonisosf'](abs(sel1_jerup_muons.eta), sel1_jerup_muons.pt)
            sel1_jerup_muonhltsf = evaluator['muonhltsf'](abs(sel1_jerup_muons.eta), sel1_jerup_muons.pt)
            sel1_jerup_btagsf = btagsf_jerup[cut1_jerup]

            sel1_jerup_X = pandas.DataFrame(np.transpose(np.vstack((
                ak.to_numpy(sel1_jerup_nextrajets),
                ak.to_numpy(sel1_jerup_nextrabjets),
                np.zeros(len(sel1_jerup_events)),
                np.sign(ak.to_numpy(sel1_jerup_muons.charge)+1),
                ak.to_numpy(sel1_jerup_muons.pt),
                ak.to_numpy(sel1_jerup_muons.eta),
                ak.to_numpy(sel1_jerup_muons.phi),
                ak.to_numpy(sel1_jerup_events.PuppiMET.ptJERUp),
                ak.to_numpy(sel1_jerup_events.PuppiMET.phiJERUp),
                ak.to_numpy(sel1_jerup_b_jets[:,0].pt),
                ak.to_numpy(sel1_jerup_b_jets[:,1].pt),
                ak.to_numpy(sel1_jerup_vbf_jets[:,0].pt),
                ak.to_numpy(sel1_jerup_vbf_jets[:,1].pt),
                ak.to_numpy(sel1_jerup_b_jets[:,0].eta),
                ak.to_numpy(sel1_jerup_b_jets[:,1].eta),
                ak.to_numpy(sel1_jerup_vbf_jets[:,0].eta),
                ak.to_numpy(sel1_jerup_vbf_jets[:,1].eta),
                ak.to_numpy(sel1_jerup_b_jets[:,0].phi),
                ak.to_numpy(sel1_jerup_b_jets[:,1].phi),
                ak.to_numpy(sel1_jerup_vbf_jets[:,0].phi),
                ak.to_numpy(sel1_jerup_vbf_jets[:,1].phi),
                ak.to_numpy(sel1_jerup_b_jets[:,0].btagDeepB),
                ak.to_numpy(sel1_jerup_b_jets[:,1].btagDeepB),
                ak.to_numpy(sel1_jerup_vbf_jets[:,0].btagDeepB),
                ak.to_numpy(sel1_jerup_vbf_jets[:,1].btagDeepB),
                ak.to_numpy((sel1_jerup_b_jets[:,0]+sel1_jerup_b_jets[:,1]).mass),
                ak.to_numpy((sel1_jerup_vbf_jets[:,0]+sel1_jerup_vbf_jets[:,1]).mass),
                ak.to_numpy(sel1_jerup_vbf_jets[:,0].eta - sel1_jerup_vbf_jets[:,1].eta),
                ak.to_numpy(np.sqrt(2*(sel1_jerup_muons+sel1_jerup_b_jets[:,0]).pt*sel1_jerup_events.PuppiMET.ptJERUp*(1 - np.cos(sel1_jerup_events.PuppiMET.phi - (sel1_jerup_muons+sel1_jerup_b_jets[:,0]).phi)))),
                ak.to_numpy(np.sqrt(2*(sel1_jerup_muons+sel1_jerup_b_jets[:,1]).pt*sel1_jerup_events.PuppiMET.ptJERUp*(1 - np.cos(sel1_jerup_events.PuppiMET.phiJERUp - (sel1_jerup_muons+sel1_jerup_b_jets[:,1]).phi))))))),
                columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbfjet1pt','vbfjet2pt','higgsjet1eta','higgsjet2eta','vbfjet1eta','vbfjet2eta','higgsjet1phi','higgsjet2phi','vbfjet1phi','vbfjet2phi','higgsjet1btag','higgsjet2btag','vbfjet1btag','vbfjet2btag','higgsdijetmass','vbfdijetmass','vbfdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel1_jerup_d = xgb.DMatrix(sel1_jerup_X)

            sel1_jerup_bdtscore = bst.predict(sel1_jerup_d)

            if dataset == 'ewkwhjj_reweighted':
                sel1_jerup_weight = np.sign(sel1_jerup_events.Generator.weight)*sel1_jerup_pu_weight*sel1_jerup_events.L1PreFiringWeight.Nom*sel1_jerup_muonidsf*sel1_jerup_muonisosf*sel1_jerup_muonhltsf*sel1_jerup_btagsf*sel1_jerup_events.LHEReweightingWeight[:,9]
            else:
                sel1_jerup_weight = np.sign(sel1_jerup_events.Generator.weight)*sel1_jerup_pu_weight*sel1_jerup_events.L1PreFiringWeight.Nom*sel1_jerup_muonidsf*sel1_jerup_muonisosf*sel1_jerup_muonhltsf*sel1_jerup_btagsf

            output['sel1_bdtscore_binning1_jerup'].fill(
                dataset=dataset,
                bdtscore=sel1_jerup_bdtscore,
                weight=sel1_jerup_weight
            )

            output['sel3_bdtscore_binning1_jerup'].fill(
                dataset=dataset,
                bdtscore=sel1_jerup_bdtscore,
                weight=sel1_jerup_weight
            )

        if ak.any(basecut) and ak.any(cut1):

            sel1_events = events[cut1]
            sel1_b_jets = b_jets[cut1]
            sel1_vbf_jets = vbf_jets[cut1]
            sel1_muons = tight_muons[cut1][:,0]
            sel1_nextrajets = nextrajets[cut1]
            sel1_nextrabjets = nextrabjets[cut1]
            if dataset != 'data':
                sel1_btagsf = btagsf[cut1]
                sel1_btagsf_lightup = btagsf_lightup[cut1]
                sel1_btagsf_bcup = btagsf_bcup[cut1]

            sel1_X = pandas.DataFrame(np.transpose(np.vstack((
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
                ak.to_numpy(np.sqrt(2*(sel1_muons+sel1_b_jets[:,0]).pt*sel1_events.PuppiMET.pt*(1 - np.cos(sel1_events.PuppiMET.phi - (sel1_muons+sel1_b_jets[:,0]).phi)))),ak.to_numpy(np.sqrt(2*(sel1_muons+sel1_b_jets[:,1]).pt*sel1_events.PuppiMET.pt*(1 - np.cos(sel1_events.PuppiMET.phi - (sel1_muons+sel1_b_jets[:,1]).phi))))))),
                columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbfjet1pt','vbfjet2pt','higgsjet1eta','higgsjet2eta','vbfjet1eta','vbfjet2eta','higgsjet1phi','higgsjet2phi','vbfjet1phi','vbfjet2phi','higgsjet1btag','higgsjet2btag','vbfjet1btag','vbfjet2btag','higgsdijetmass','vbfdijetmass','vbfdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel1_d = xgb.DMatrix(sel1_X)

            sel1_bdtscore = bst.predict(sel1_d)    

            if dataset == 'data':
                sel1_weight = np.ones(len(sel1_events))
            elif dataset == 'ewkwhjj_reweighted':
                sel1_pu_weight = evaluator['pileup'](sel1_events.Pileup.nTrueInt)
                sel1_puUp_weight = evaluator['pileup_up'](sel1_events.Pileup.nTrueInt)
                sel1_puDown_weight = evaluator['pileup_down'](sel1_events.Pileup.nTrueInt)
                sel1_muonidsf = evaluator['muonidsf'](abs(sel1_muons.eta), sel1_muons.pt)
                sel1_muonisosf = evaluator['muonisosf'](abs(sel1_muons.eta), sel1_muons.pt)
                sel1_muonhltsf = evaluator['muonhltsf'](abs(sel1_muons.eta), sel1_muons.pt)
                sel1_muonidsfUp = evaluator['muonidsfunc'](abs(sel1_muons.eta), sel1_muons.pt)+sel1_muonidsf
                sel1_muonisosfUp = evaluator['muonisosfunc'](abs(sel1_muons.eta), sel1_muons.pt)+sel1_muonisosf
                sel1_muonhltsfUp = evaluator['muonhltsfunc'](abs(sel1_muons.eta), sel1_muons.pt)+sel1_muonhltsf

                sel1_weight = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_events.LHEReweightingWeight[:,9]*sel1_btagsf
                sel1_weight_reweighted = []
                for i in range(args.nreweights):
                    sel1_weight_reweighted.append(np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_events.LHEReweightingWeight[:,i]*sel1_btagsf)

                sel1_weight_pileupup = np.sign(sel1_events.Generator.weight)*sel1_puUp_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_events.LHEReweightingWeight[:,9]*sel1_btagsf
                sel1_weight_pileupdown = np.sign(sel1_events.Generator.weight)*sel1_puDown_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_events.LHEReweightingWeight[:,9]*sel1_btagsf
                sel1_weight_prefireup = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Up*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_events.LHEReweightingWeight[:,9]*sel1_btagsf
                sel1_weight_muonidsfup = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsfUp*sel1_muonisosf*sel1_muonhltsf*sel1_events.LHEReweightingWeight[:,9]*sel1_btagsf
                sel1_weight_muonisosfup = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosfUp*sel1_muonhltsf*sel1_events.LHEReweightingWeight[:,9]*sel1_btagsf
                sel1_weight_muonhltsfup = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsfUp*sel1_events.LHEReweightingWeight[:,9]*sel1_btagsf
                sel1_weight_bcbtagsfup = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_events.LHEReweightingWeight[:,9]*sel1_btagsf_bcup
                sel1_weight_lightbtagsfup = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_events.LHEReweightingWeight[:,9]*sel1_btagsf_lightup

            else:
                sel1_pu_weight = evaluator['pileup'](sel1_events.Pileup.nTrueInt)
                sel1_puup_weight = evaluator['pileup_up'](sel1_events.Pileup.nTrueInt)
                sel1_pudown_weight = evaluator['pileup_down'](sel1_events.Pileup.nTrueInt)
                sel1_muonidsf = evaluator['muonidsf'](abs(sel1_muons.eta), sel1_muons.pt)
                sel1_muonisosf = evaluator['muonisosf'](abs(sel1_muons.eta), sel1_muons.pt)
                sel1_muonhltsf = evaluator['muonhltsf'](abs(sel1_muons.eta), sel1_muons.pt)
                sel1_muonidsfUp = evaluator['muonidsfunc'](abs(sel1_muons.eta), sel1_muons.pt)+sel1_muonidsf
                sel1_muonisosfUp = evaluator['muonisosfunc'](abs(sel1_muons.eta), sel1_muons.pt)+sel1_muonisosf
                sel1_muonhltsfUp = evaluator['muonhltsfunc'](abs(sel1_muons.eta), sel1_muons.pt)+sel1_muonhltsf

                sel1_weight = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_btagsf
                sel1_weight_pileupup = np.sign(sel1_events.Generator.weight)*sel1_puup_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_btagsf
                sel1_weight_pileupdown = np.sign(sel1_events.Generator.weight)*sel1_pudown_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_btagsf
                sel1_weight_prefireup = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Up*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_btagsf
                sel1_weight_muonidsfup = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsfUp*sel1_muonisosf*sel1_muonhltsf*sel1_btagsf
                sel1_weight_muonisosfup = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosfUp*sel1_muonhltsf*sel1_btagsf
                sel1_weight_muonhltsfup = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsfUp*sel1_btagsf
                sel1_weight_bcbtagsfup = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_btagsf_bcup
                sel1_weight_lightbtagsfup = np.sign(sel1_events.Generator.weight)*sel1_pu_weight*sel1_events.L1PreFiringWeight.Nom*sel1_muonidsf*sel1_muonisosf*sel1_muonhltsf*sel1_btagsf_lightup
                
            output['sel1_bdtscore_binning1'].fill(
                dataset=dataset,
                bdtscore=sel1_bdtscore,
                weight=sel1_weight
            )

            if dataset == 'ewkwhjj_reweighted':
                for i in range(args.nreweights):
                    output['sel1_bdtscore_binning1'.format(i)].fill(
                        dataset='{}_reweightingweight{}'.format(dataset,i),
                        bdtscore = sel1_bdtscore,
                        weight=sel1_weight_reweighted[i]
                    )

            if dataset != 'data':
                output['sel1_bdtscore_binning1_pileupup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_pileupup
                )

                output['sel1_bdtscore_binning1_pileupdown'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_pileupdown
                )

                output['sel1_bdtscore_binning1_prefireup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_prefireup
                )

                output['sel1_bdtscore_binning1_muonidsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonidsfup
                )
                
                output['sel1_bdtscore_binning1_muonisosfup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonisosfup
                )

                output['sel1_bdtscore_binning1_muonhltsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonhltsfup
                )

                output['sel1_bdtscore_binning1_bcbtagsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_bcbtagsfup
                )

                output['sel1_bdtscore_binning1_lightbtagsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_lightbtagsfup
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
                higgsdijetmass=(sel1_b_jets[:,0]+sel1_b_jets[:,1]).mass,
                weight=sel1_weight
            )

            output['sel1_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=(sel1_b_jets[:,0]+sel1_b_jets[:,1]).pt,
                weight=sel1_weight
            )
        
            output['sel1_vbfdijetmass_binning1'].fill(
                dataset=dataset,
                vbfdijetmass=(sel1_vbf_jets[:,0]+sel1_vbf_jets[:,1]).mass,
                weight=sel1_weight
            )

            output['sel1_vbfdijetmass_binning2'].fill(
                dataset=dataset,
                vbfdijetmass=(sel1_vbf_jets[:,0]+sel1_vbf_jets[:,1]).mass,
                weight=sel1_weight
            )

            output['sel1_vbfdijetabsdeta_binning1'].fill(
                dataset=dataset,
                vbfdijetabsdeta=abs(sel1_vbf_jets[:,0].eta - sel1_vbf_jets[:,1].eta),
                weight=sel1_weight
            )

            output['sel1_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=sel1_muons.pt,
                weight=sel1_weight
            )

            output['sel1_leptonabseta_binning1'].fill(
                dataset=dataset,
                leptonabseta=abs(sel1_muons.eta),
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

            if dataset == "ewkwhjj_reweighted":
                for i in range(args.nreweights):
                    output['sel3_bdtscore_binning1'].fill(
                        dataset='{}_reweightingweight{}'.format(dataset,i),
                        bdtscore = sel1_bdtscore,
                        weight=sel1_weight_reweighted[i]
                    )

            if dataset != 'data':

                output['sel3_bdtscore_binning1_pileupup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_pileupup
                )

                output['sel3_bdtscore_binning1_pileupdown'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_pileupdown
                )

                output['sel3_bdtscore_binning1_prefireup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_prefireup
                )
                
                output['sel3_bdtscore_binning1_electronidsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight
                )

                output['sel3_bdtscore_binning1_electronrecosfup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight
                )
                
                output['sel3_bdtscore_binning1_muonidsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonidsfup
                )
                
                output['sel3_bdtscore_binning1_muonisosfup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonisosfup
                )
                
                output['sel3_bdtscore_binning1_muonhltsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_muonhltsfup
                )

                output['sel3_bdtscore_binning1_bcbtagsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_bcbtagsfup
                )

                output['sel3_bdtscore_binning1_lightbtagsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel1_bdtscore,
                    weight=sel1_weight_lightbtagsfup
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
                higgsdijetmass=(sel1_b_jets[:,0]+sel1_b_jets[:,1]).mass,
                weight=sel1_weight
            )

            output['sel3_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=(sel1_b_jets[:,0]+sel1_b_jets[:,1]).mass,
                weight=sel1_weight
            )
        
            output['sel3_vbfdijetmass_binning1'].fill(
                dataset=dataset,
                vbfdijetmass=(sel1_vbf_jets[:,0]+sel1_vbf_jets[:,1]).mass,
                weight=sel1_weight
            )

            output['sel3_vbfdijetmass_binning2'].fill(
                dataset=dataset,
                vbfdijetmass=(sel1_vbf_jets[:,0]+sel1_vbf_jets[:,1]).mass,
                weight=sel1_weight
            )

            output['sel3_vbfdijetabsdeta_binning1'].fill(
                dataset=dataset,
                vbfdijetabsdeta=abs(sel1_vbf_jets[:,0].eta - sel1_vbf_jets[:,1].eta),
                weight=sel1_weight
            )

            output['sel3_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=sel1_muons.pt,
                weight=sel1_weight
            )

            output['sel3_leptonabseta_binning1'].fill(
                dataset=dataset,
                leptonabseta=abs(sel1_muons.eta),
                weight=sel1_weight
            )

            output['sel3_met_binning1'].fill(
                dataset=dataset,
                met=sel1_events.PuppiMET.pt,
                weight=sel1_weight
            )

        if dataset != 'data' and ak.any(basecut_jesup) and ak.any(cut2_jesup):

            sel2_jesup_events = events_jesup[cut2_jesup]
            sel2_jesup_b_jets = b_jets_jesup[cut2_jesup]
            sel2_jesup_vbf_jets = vbf_jets_jesup[cut2_jesup]
            sel2_jesup_electrons = tight_electrons_jesup[cut2_jesup][:,0]
            sel2_jesup_nextrajets = nextrajets_jesup[cut2_jesup]
            sel2_jesup_nextrabjets = nextrabjets_jesup[cut2_jesup]

            sel2_jesup_pu_weight = evaluator['pileup'](sel2_jesup_events.Pileup.nTrueInt)
            sel2_jesup_electronidsf = evaluator['electronidsf'](sel2_jesup_electrons.eta, sel2_jesup_electrons.pt)
            sel2_jesup_electronrecosf = evaluator['electronrecosf'](sel2_jesup_electrons.eta, sel2_jesup_electrons.pt)
            sel2_jesup_btagsf = btagsf_jesup[cut2_jesup]

            sel2_jesup_X = pandas.DataFrame(np.transpose(np.vstack((
                ak.to_numpy(sel2_jesup_nextrajets),
                ak.to_numpy(sel2_jesup_nextrabjets),
                np.ones(len(sel2_jesup_events)),
                np.sign(ak.to_numpy(sel2_jesup_electrons.charge)+1),
                ak.to_numpy(sel2_jesup_electrons.pt),
                ak.to_numpy(sel2_jesup_electrons.eta),
                ak.to_numpy(sel2_jesup_electrons.phi),
                ak.to_numpy(sel2_jesup_events.PuppiMET.ptJESUp),
                ak.to_numpy(sel2_jesup_events.PuppiMET.phiJESUp),
                ak.to_numpy(sel2_jesup_b_jets[:,0].pt),
                ak.to_numpy(sel2_jesup_b_jets[:,1].pt),
                ak.to_numpy(sel2_jesup_vbf_jets[:,0].pt),
                ak.to_numpy(sel2_jesup_vbf_jets[:,1].pt),
                ak.to_numpy(sel2_jesup_b_jets[:,0].eta),
                ak.to_numpy(sel2_jesup_b_jets[:,1].eta),
                ak.to_numpy(sel2_jesup_vbf_jets[:,0].eta),
                ak.to_numpy(sel2_jesup_vbf_jets[:,1].eta),
                ak.to_numpy(sel2_jesup_b_jets[:,0].phi),
                ak.to_numpy(sel2_jesup_b_jets[:,1].phi),
                ak.to_numpy(sel2_jesup_vbf_jets[:,0].phi),
                ak.to_numpy(sel2_jesup_vbf_jets[:,1].phi),
                ak.to_numpy(sel2_jesup_b_jets[:,0].btagDeepB),
                ak.to_numpy(sel2_jesup_b_jets[:,1].btagDeepB),
                ak.to_numpy(sel2_jesup_vbf_jets[:,0].btagDeepB),
                ak.to_numpy(sel2_jesup_vbf_jets[:,1].btagDeepB),
                ak.to_numpy((sel2_jesup_b_jets[:,0]+sel2_jesup_b_jets[:,1]).mass),
                ak.to_numpy((sel2_jesup_vbf_jets[:,0]+sel2_jesup_vbf_jets[:,1]).mass),
                ak.to_numpy(sel2_jesup_vbf_jets[:,0].eta - sel2_jesup_vbf_jets[:,1].eta),
                ak.to_numpy(np.sqrt(2*(sel2_jesup_electrons+sel2_jesup_b_jets[:,0]).pt*sel2_jesup_events.PuppiMET.ptJESUp*(1 - np.cos(sel2_jesup_events.PuppiMET.phi - (sel2_jesup_electrons+sel2_jesup_b_jets[:,0]).phi)))),
                ak.to_numpy(np.sqrt(2*(sel2_jesup_electrons+sel2_jesup_b_jets[:,1]).pt*sel2_jesup_events.PuppiMET.ptJESUp*(1 - np.cos(sel2_jesup_events.PuppiMET.phiJESUp - (sel2_jesup_electrons+sel2_jesup_b_jets[:,1]).phi))))))),
                columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbfjet1pt','vbfjet2pt','higgsjet1eta','higgsjet2eta','vbfjet1eta','vbfjet2eta','higgsjet1phi','higgsjet2phi','vbfjet1phi','vbfjet2phi','higgsjet1btag','higgsjet2btag','vbfjet1btag','vbfjet2btag','higgsdijetmass','vbfdijetmass','vbfdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel2_jesup_d = xgb.DMatrix(sel2_jesup_X)

            sel2_jesup_bdtscore = bst.predict(sel2_jesup_d)

            if dataset == 'ewkwhjj_reweighted':
                sel2_jesup_weight = np.sign(sel2_jesup_events.Generator.weight)*sel2_jesup_pu_weight*sel2_jesup_events.L1PreFiringWeight.Nom*sel2_jesup_electronidsf*sel2_jesup_electronrecosf*sel2_jesup_btagsf*sel2_jesup_events.LHEReweightingWeight[:,9]
            else:    
                sel2_jesup_weight = np.sign(sel2_jesup_events.Generator.weight)*sel2_jesup_pu_weight*sel2_jesup_events.L1PreFiringWeight.Nom*sel2_jesup_electronidsf*sel2_jesup_electronrecosf*sel2_jesup_btagsf

            output['sel2_bdtscore_binning1_jesup'].fill(
                dataset=dataset,
                bdtscore=sel2_jesup_bdtscore,
                weight=sel2_jesup_weight
            )

            output['sel3_bdtscore_binning1_jesup'].fill(
                dataset=dataset,
                bdtscore=sel2_jesup_bdtscore,
                weight=sel2_jesup_weight
            )

        if dataset != 'data' and ak.any(basecut_jerup) and ak.any(cut2_jerup):

            sel2_jerup_events = events_jerup[cut2_jerup]
            sel2_jerup_b_jets = b_jets_jerup[cut2_jerup]
            sel2_jerup_vbf_jets = vbf_jets_jerup[cut2_jerup]
            sel2_jerup_electrons = tight_electrons_jerup[cut2_jerup][:,0]
            sel2_jerup_nextrajets = nextrajets_jerup[cut2_jerup]
            sel2_jerup_nextrabjets = nextrabjets_jerup[cut2_jerup]

            sel2_jerup_pu_weight = evaluator['pileup'](sel2_jerup_events.Pileup.nTrueInt)
            sel2_jerup_electronidsf = evaluator['electronidsf'](sel2_jerup_electrons.eta, sel2_jerup_electrons.pt)
            sel2_jerup_electronrecosf = evaluator['electronrecosf'](sel2_jerup_electrons.eta, sel2_jerup_electrons.pt)
            sel2_jerup_btagsf = btagsf_jerup[cut2_jerup]

            sel2_jerup_X = pandas.DataFrame(np.transpose(np.vstack((
                ak.to_numpy(sel2_jerup_nextrajets),
                ak.to_numpy(sel2_jerup_nextrabjets),
                np.ones(len(sel2_jerup_events)),
                np.sign(ak.to_numpy(sel2_jerup_electrons.charge)+1),
                ak.to_numpy(sel2_jerup_electrons.pt),
                ak.to_numpy(sel2_jerup_electrons.eta),
                ak.to_numpy(sel2_jerup_electrons.phi),
                ak.to_numpy(sel2_jerup_events.PuppiMET.ptJERUp),
                ak.to_numpy(sel2_jerup_events.PuppiMET.phiJERUp),
                ak.to_numpy(sel2_jerup_b_jets[:,0].pt),
                ak.to_numpy(sel2_jerup_b_jets[:,1].pt),
                ak.to_numpy(sel2_jerup_vbf_jets[:,0].pt),
                ak.to_numpy(sel2_jerup_vbf_jets[:,1].pt),
                ak.to_numpy(sel2_jerup_b_jets[:,0].eta),
                ak.to_numpy(sel2_jerup_b_jets[:,1].eta),
                ak.to_numpy(sel2_jerup_vbf_jets[:,0].eta),
                ak.to_numpy(sel2_jerup_vbf_jets[:,1].eta),
                ak.to_numpy(sel2_jerup_b_jets[:,0].phi),
                ak.to_numpy(sel2_jerup_b_jets[:,1].phi),
                ak.to_numpy(sel2_jerup_vbf_jets[:,0].phi),
                ak.to_numpy(sel2_jerup_vbf_jets[:,1].phi),
                ak.to_numpy(sel2_jerup_b_jets[:,0].btagDeepB),
                ak.to_numpy(sel2_jerup_b_jets[:,1].btagDeepB),
                ak.to_numpy(sel2_jerup_vbf_jets[:,0].btagDeepB),
                ak.to_numpy(sel2_jerup_vbf_jets[:,1].btagDeepB),
                ak.to_numpy((sel2_jerup_b_jets[:,0]+sel2_jerup_b_jets[:,1]).mass),
                ak.to_numpy((sel2_jerup_vbf_jets[:,0]+sel2_jerup_vbf_jets[:,1]).mass),
                ak.to_numpy(sel2_jerup_vbf_jets[:,0].eta - sel2_jerup_vbf_jets[:,1].eta),
                ak.to_numpy(np.sqrt(2*(sel2_jerup_electrons+sel2_jerup_b_jets[:,0]).pt*sel2_jerup_events.PuppiMET.ptJERUp*(1 - np.cos(sel2_jerup_events.PuppiMET.phi - (sel2_jerup_electrons+sel2_jerup_b_jets[:,0]).phi)))),
                ak.to_numpy(np.sqrt(2*(sel2_jerup_electrons+sel2_jerup_b_jets[:,1]).pt*sel2_jerup_events.PuppiMET.ptJERUp*(1 - np.cos(sel2_jerup_events.PuppiMET.phiJERUp - (sel2_jerup_electrons+sel2_jerup_b_jets[:,1]).phi))))))),
                columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbfjet1pt','vbfjet2pt','higgsjet1eta','higgsjet2eta','vbfjet1eta','vbfjet2eta','higgsjet1phi','higgsjet2phi','vbfjet1phi','vbfjet2phi','higgsjet1btag','higgsjet2btag','vbfjet1btag','vbfjet2btag','higgsdijetmass','vbfdijetmass','vbfdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel2_jerup_d = xgb.DMatrix(sel2_jerup_X)

            sel2_jerup_bdtscore = bst.predict(sel2_jerup_d)

            if dataset == 'ewkwhjj_reweighted':
                sel2_jerup_weight = np.sign(sel2_jerup_events.Generator.weight)*sel2_jerup_pu_weight*sel2_jerup_events.L1PreFiringWeight.Nom*sel2_jerup_electronidsf*sel2_jerup_electronrecosf*sel2_jerup_btagsf*sel2_jerup_events.LHEReweightingWeight[:,9]
            else:    
                sel2_jerup_weight = np.sign(sel2_jerup_events.Generator.weight)*sel2_jerup_pu_weight*sel2_jerup_events.L1PreFiringWeight.Nom*sel2_jerup_electronidsf*sel2_jerup_electronrecosf*sel2_jerup_btagsf

            output['sel2_bdtscore_binning1_jerup'].fill(
                dataset=dataset,
                bdtscore=sel2_jerup_bdtscore,
                weight=sel2_jerup_weight
            )

            output['sel3_bdtscore_binning1_jerup'].fill(
                dataset=dataset,
                bdtscore=sel2_jerup_bdtscore,
                weight=sel2_jerup_weight
            )

        if ak.any(basecut) and ak.any(cut2):
                
            sel2_events = events[cut2]
            sel2_b_jets = b_jets[cut2]
            sel2_vbf_jets = vbf_jets[cut2]
            sel2_electrons = tight_electrons[cut2][:,0]
            sel2_nextrajets = nextrajets[cut2]
            sel2_nextrabjets = nextrabjets[cut2]
            if dataset != 'data':
                sel2_btagsf = btagsf[cut2]
                sel2_btagsf_bcup = btagsf_bcup[cut2]
                sel2_btagsf_lightup = btagsf_lightup[cut2]

            sel2_X = pandas.DataFrame(np.transpose(np.vstack((
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
                ak.to_numpy(np.sqrt(2*(sel2_electrons+sel2_b_jets[:,0]).pt*sel2_events.PuppiMET.pt*(1 - np.cos(sel2_events.PuppiMET.phi - (sel2_electrons+sel2_b_jets[:,0]).phi)))),ak.to_numpy(np.sqrt(2*(sel2_electrons+sel2_b_jets[:,1]).pt*sel2_events.PuppiMET.pt*(1 - np.cos(sel2_events.PuppiMET.phi - (sel2_electrons+sel2_b_jets[:,1]).phi))))))),
                columns=['nextrajets','nextrabjets','leptonflavor','leptoncharge','leptonpt','leptoneta','leptonphi','metpt','metphi','higgsjet1pt','higgsjet2pt','vbfjet1pt','vbfjet2pt','higgsjet1eta','higgsjet2eta','vbfjet1eta','vbfjet2eta','higgsjet1phi','higgsjet2phi','vbfjet1phi','vbfjet2phi','higgsjet1btag','higgsjet2btag','vbfjet1btag','vbfjet2btag','higgsdijetmass','vbfdijetmass','vbfdijetabsdeta','leptonhiggsjet1mt','leptonhiggsjet2mt'])

            sel2_d = xgb.DMatrix(sel2_X)

            sel2_bdtscore = bst.predict(sel2_d)    

            if dataset == 'data':
                sel2_weight = np.ones(len(sel2_events))
            elif dataset == 'ewkwhjj_reweighted':
                sel2_pu_weight = evaluator['pileup'](sel2_events.Pileup.nTrueInt)
                sel2_puUp_weight = evaluator['pileup_up'](sel2_events.Pileup.nTrueInt)
                sel2_puDown_weight = evaluator['pileup_down'](sel2_events.Pileup.nTrueInt)
                sel2_electronidsf = evaluator['electronidsf'](sel2_electrons.eta, sel2_electrons.pt)
                sel2_electronidsfUp = evaluator['electronidsfunc'](sel2_electrons.eta, sel2_electrons.pt)+sel2_electronidsf
                sel2_electronrecosf = evaluator['electronrecosf'](sel2_electrons.eta, sel2_electrons.pt)
                sel2_electronrecosfUp = evaluator['electronrecosfunc'](sel2_electrons.eta, sel2_electrons.pt)+sel2_electronrecosf

                sel2_weight = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf*sel2_electronrecosf*sel2_events.LHEReweightingWeight[:,9]*sel2_btagsf
                sel2_weight_reweighted = []
                for i in range(args.nreweights):
                    sel2_weight_reweighted.append(np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf*sel2_electronrecosf*sel2_events.LHEReweightingWeight[:,i]*sel2_btagsf)
                sel2_weight_pileupup = np.sign(sel2_events.Generator.weight)*sel2_puUp_weight*sel2_events.L1PreFiringWeight.Up*sel2_electronidsf*sel2_electronrecosf*sel2_events.LHEReweightingWeight[:,9]*sel2_btagsf
                sel2_weight_pileupdown = np.sign(sel2_events.Generator.weight)*sel2_puDown_weight*sel2_events.L1PreFiringWeight.Up*sel2_electronidsf*sel2_electronrecosf*sel2_events.LHEReweightingWeight[:,9]*sel2_btagsf
                sel2_weight_prefireup = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Up*sel2_electronidsf*sel2_electronrecosf*sel2_events.LHEReweightingWeight[:,9]*sel2_btagsf
                sel2_weight_electronidsfup = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsfUp*sel2_electronrecosf*sel2_events.LHEReweightingWeight[:,9]*sel2_btagsf
                sel2_weight_electronrecosfup = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf*sel2_electronrecosfUp*sel2_events.LHEReweightingWeight[:,9]*sel2_btagsf
                sel2_weight_bcbtagsfup = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf*sel2_electronrecosf*sel2_events.LHEReweightingWeight[:,9]*sel2_btagsf_bcup
                sel2_weight_lightbtagsfup = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf*sel2_electronrecosf*sel2_events.LHEReweightingWeight[:,9]*sel2_btagsf_lightup
            else:
                sel2_pu_weight = evaluator['pileup'](sel2_events.Pileup.nTrueInt)
                sel2_puup_weight = evaluator['pileup_up'](sel2_events.Pileup.nTrueInt)
                sel2_pudown_weight = evaluator['pileup_down'](sel2_events.Pileup.nTrueInt)
                sel2_electronidsf = evaluator['electronidsf'](sel2_electrons.eta, sel2_electrons.pt)
                sel2_electronidsfUp = evaluator['electronidsfunc'](sel2_electrons.eta, sel2_electrons.pt)+sel2_electronidsf
                sel2_electronrecosf = evaluator['electronrecosf'](sel2_electrons.eta, sel2_electrons.pt)
                sel2_electronrecosfUp = evaluator['electronrecosfunc'](sel2_electrons.eta, sel2_electrons.pt)+sel2_electronrecosf

                sel2_weight = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf*sel2_electronrecosf*sel2_btagsf
                sel2_weight_pileupup = np.sign(sel2_events.Generator.weight)*sel2_puup_weight*sel2_events.L1PreFiringWeight.Up*sel2_electronidsf*sel2_electronrecosf*sel2_btagsf
                sel2_weight_pileupdown = np.sign(sel2_events.Generator.weight)*sel2_pudown_weight*sel2_events.L1PreFiringWeight.Up*sel2_electronidsf*sel2_electronrecosf*sel2_btagsf
                sel2_weight_prefireup = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Up*sel2_electronidsf*sel2_electronrecosf*sel2_btagsf
                sel2_weight_electronidsfup = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsfUp*sel2_electronrecosf*sel2_btagsf
                sel2_weight_electronrecosfup = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf*sel2_electronrecosfUp*sel2_btagsf
                sel2_weight_bcbtagsfup = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf*sel2_electronrecosf*sel2_btagsf_bcup
                sel2_weight_lightbtagsfup = np.sign(sel2_events.Generator.weight)*sel2_pu_weight*sel2_events.L1PreFiringWeight.Nom*sel2_electronidsf*sel2_electronrecosf*sel2_btagsf_lightup

            output['sel2_bdtscore_binning1'].fill(
                dataset=dataset,
                bdtscore = sel2_bdtscore,
                weight=sel2_weight
            )


            if dataset == "ewkwhjj_reweighted":
                for i in range(args.nreweights):
                    output['sel2_bdtscore_binning1'].fill(
                        dataset='{}_reweightingweight{}'.format(dataset,i),
                        bdtscore = sel2_bdtscore,
                        weight=sel2_weight_reweighted[i]
                    )

            if dataset != 'data':

                output['sel2_bdtscore_binning1_pileupup'].fill(
                    dataset=dataset,
                    bdtscore = sel2_bdtscore,
                    weight=sel2_weight_pileupup
                )

                output['sel2_bdtscore_binning1_pileupdown'].fill(
                    dataset=dataset,
                    bdtscore = sel2_bdtscore,
                    weight=sel2_weight_pileupdown
                )

                output['sel2_bdtscore_binning1_prefireup'].fill(
                    dataset=dataset,
                    bdtscore = sel2_bdtscore,
                    weight=sel2_weight_prefireup
                )

                output['sel2_bdtscore_binning1_electronidsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight_electronidsfup
                )

                output['sel2_bdtscore_binning1_electronrecosfup'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight_electronrecosfup
                )

                output['sel2_bdtscore_binning1_bcbtagsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight_bcbtagsfup
                )

                output['sel2_bdtscore_binning1_lightbtagsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight_lightbtagsfup
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
                higgsdijetmass=(sel2_b_jets[:,0]+sel2_b_jets[:,1]).mass,
                weight=sel2_weight
            )

            output['sel2_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=(sel2_b_jets[:,0]+sel2_b_jets[:,1]).pt,
                weight=sel2_weight
            )
        
            output['sel2_vbfdijetmass_binning1'].fill(
                dataset=dataset,
                vbfdijetmass=(sel2_vbf_jets[:,0]+sel2_vbf_jets[:,1]).mass,
                weight=sel2_weight
            )

            output['sel2_vbfdijetmass_binning2'].fill(
                dataset=dataset,
                vbfdijetmass=(sel2_vbf_jets[:,0]+sel2_vbf_jets[:,1]).mass,
                weight=sel2_weight
            )

            output['sel2_vbfdijetabsdeta_binning1'].fill(
                dataset=dataset,
                vbfdijetabsdeta=abs(sel2_vbf_jets[:,0].eta - sel2_vbf_jets[:,1].eta),
                weight=sel2_weight
            )

            output['sel2_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=sel2_electrons.pt,
                weight=sel2_weight
            )

            output['sel2_leptonabseta_binning1'].fill(
                dataset=dataset,
                leptonabseta=abs(sel2_electrons.eta),
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

            if dataset == "ewkwhjj_reweighted":
                for i in range(args.nreweights):
                    output['sel3_bdtscore_binning1'].fill(
                        dataset='{}_reweightingweight{}'.format(dataset,i),
                        bdtscore = sel2_bdtscore,
                        weight=sel2_weight_reweighted[i]
                    )

            if dataset != 'data':

                output['sel3_bdtscore_binning1_pileupup'].fill(
                    dataset=dataset,
                    bdtscore = sel2_bdtscore,
                    weight=sel2_weight_pileupup
                )

                output['sel3_bdtscore_binning1_pileupdown'].fill(
                    dataset=dataset,
                    bdtscore = sel2_bdtscore,
                    weight=sel2_weight_pileupdown
                )

                output['sel3_bdtscore_binning1_prefireup'].fill(
                    dataset=dataset,
                    bdtscore = sel2_bdtscore,
                    weight=sel2_weight_prefireup
                )
                
                output['sel3_bdtscore_binning1_electronidsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight_electronidsfup
                )

                output['sel3_bdtscore_binning1_electronrecosfup'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight_electronrecosfup
                )
                
                output['sel3_bdtscore_binning1_muonidsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight
                )
                
                output['sel3_bdtscore_binning1_muonisosfup'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight
                )
                
                output['sel3_bdtscore_binning1_muonhltsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight
                )

                output['sel3_bdtscore_binning1_bcbtagsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight_bcbtagsfup
                )

                output['sel3_bdtscore_binning1_lightbtagsfup'].fill(
                    dataset=dataset,
                    bdtscore=sel2_bdtscore,
                    weight=sel2_weight_lightbtagsfup
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
                higgsdijetmass=(sel2_b_jets[:,0]+sel2_b_jets[:,1]).mass,
                weight=sel2_weight
            )

            output['sel3_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=(sel2_b_jets[:,0]+sel2_b_jets[:,1]).pt,
                weight=sel2_weight
            )
        
            output['sel3_vbfdijetmass_binning1'].fill(
                dataset=dataset,
                vbfdijetmass=(sel2_vbf_jets[:,0]+sel2_vbf_jets[:,1]).mass,
                weight=sel2_weight
            )

            output['sel3_vbfdijetmass_binning2'].fill(
                dataset=dataset,
                vbfdijetmass=(sel2_vbf_jets[:,0]+sel2_vbf_jets[:,1]).mass,
                weight=sel2_weight
            )

            output['sel3_vbfdijetabsdeta_binning1'].fill(
                dataset=dataset,
                vbfdijetabsdeta=abs(sel2_vbf_jets[:,0].eta - sel2_vbf_jets[:,1].eta),
                weight=sel2_weight
            )

            output['sel3_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=sel2_electrons.pt,
                weight=sel2_weight
            )

            output['sel3_leptonabseta_binning1'].fill(
                dataset=dataset,
                leptonabseta=abs(sel2_electrons.eta),
                weight=sel2_weight
            )

            output['sel3_met_binning1'].fill(
                dataset=dataset,
                met=sel2_events.PuppiMET.pt,
                weight=sel2_weight
            )

        if ak.any(basecut) and ak.any(cut4):
                
            sel4_events = events[cut4]
            sel4_b_jets = b_jets[cut4]
            sel4_vbf_jets = vbf_jets[cut4]
            sel4_muons = tight_muons[cut4][:,0]
            
            if dataset == 'data':
                sel4_weight = np.ones(len(sel4_events))
            else:
                sel4_weight = np.sign(sel4_events.Generator.weight)

            output['sel4_higgsdijetmass_binning1'].fill(
                dataset=dataset,
                higgsdijetmass=(sel4_b_jets[:,0]+sel4_b_jets[:,1]).mass,
                weight=sel4_weight
            )

            output['sel4_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=(sel4_b_jets[:,0]+sel4_b_jets[:,1]).pt,
                weight=sel4_weight
            )
        
            output['sel4_vbfdijetmass_binning1'].fill(
                dataset=dataset,
                vbfdijetmass=(sel4_vbf_jets[:,0]+sel4_vbf_jets[:,1]).mass,
                weight=sel4_weight
            )

            output['sel4_vbfdijetmass_binning2'].fill(
                dataset=dataset,
                vbfdijetmass=(sel4_vbf_jets[:,0]+sel4_vbf_jets[:,1]).mass,
                weight=sel4_weight
            )

            output['sel4_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=sel4_muons.pt,
                weight=sel4_weight
            )

            output['sel4_leptonabseta_binning1'].fill(
                dataset=dataset,
                leptonabseta=abs(sel4_muons.eta),
                weight=sel4_weight
            )

            output['sel4_met_binning1'].fill(
                dataset=dataset,
                met=sel4_events.PuppiMET.pt,
                weight=sel4_weight
            )

            output['sel6_higgsdijetmass_binning1'].fill(
                dataset=dataset,
                higgsdijetmass=(sel4_b_jets[:,0]+sel4_b_jets[:,1]).mass,
                weight=sel4_weight
            )

            output['sel6_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=(sel4_b_jets[:,0]+sel4_b_jets[:,1]).pt,
                weight=sel4_weight
            )
        
            output['sel6_vbfdijetmass_binning1'].fill(
                dataset=dataset,
                vbfdijetmass=(sel4_vbf_jets[:,0]+sel4_vbf_jets[:,1]).mass,
                weight=sel4_weight
            )

            output['sel6_vbfdijetmass_binning2'].fill(
                dataset=dataset,
                vbfdijetmass=(sel4_vbf_jets[:,0]+sel4_vbf_jets[:,1]).mass,
                weight=sel4_weight
            )

            output['sel6_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=sel4_muons.pt,
                weight=sel4_weight
            )

            output['sel6_leptonabseta_binning1'].fill(
                dataset=dataset,
                leptonabseta=abs(sel4_muons.eta),
                weight=sel4_weight
            )

            output['sel6_met_binning1'].fill(
                dataset=dataset,
                met=sel4_events.PuppiMET.pt,
                weight=sel4_weight
            )

        if ak.any(basecut) and ak.any(cut5):
                
            sel5_events = events[cut5]
            
            sel5_b_jets = b_jets[cut5]
            sel5_vbf_jets = vbf_jets[cut5]
            sel5_electrons = tight_electrons[cut5][:,0]
            
            if dataset == 'data':
                sel5_weight = np.ones(len(sel5_events))
            else:
                sel5_weight = np.sign(sel5_events.Generator.weight)

            output['sel5_higgsdijetmass_binning1'].fill(
                dataset=dataset,
                higgsdijetmass=(sel5_b_jets[:,0]+sel5_b_jets[:,1]).mass,
                weight=sel5_weight
            )

            output['sel5_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=(sel5_b_jets[:,0]+sel5_b_jets[:,1]).pt,
                weight=sel5_weight
            )
        
            output['sel5_vbfdijetmass_binning1'].fill(
                dataset=dataset,
                vbfdijetmass=(sel5_vbf_jets[:,0]+sel5_vbf_jets[:,1]).mass,
                weight=sel5_weight
            )

            output['sel5_vbfdijetmass_binning2'].fill(
                dataset=dataset,
                vbfdijetmass=(sel5_vbf_jets[:,0]+sel5_vbf_jets[:,1]).mass,
                weight=sel5_weight
            )

            output['sel5_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=sel5_electrons.pt,
                weight=sel5_weight
            )

            output['sel5_leptonabseta_binning1'].fill(
                dataset=dataset,
                leptonabseta=abs(sel5_electrons.eta),
                weight=sel5_weight
            )

            output['sel5_met_binning1'].fill(
                dataset=dataset,
                met=sel5_events.PuppiMET.pt,
                weight=sel5_weight
            )

            output['sel6_higgsdijetmass_binning1'].fill(
                dataset=dataset,
                higgsdijetmass=(sel5_b_jets[:,0]+sel5_b_jets[:,1]).mass,
                weight=sel5_weight
            )

            output['sel6_higgsdijetpt_binning1'].fill(
                dataset=dataset,
                higgsdijetpt=(sel5_b_jets[:,0]+sel5_b_jets[:,1]).pt,
                weight=sel5_weight
            )
        
            output['sel6_vbfdijetmass_binning1'].fill(
                dataset=dataset,
                vbfdijetmass=(sel5_vbf_jets[:,0]+sel5_vbf_jets[:,1]).mass,
                weight=sel5_weight
            )

            output['sel6_vbfdijetmass_binning2'].fill(
                dataset=dataset,
                vbfdijetmass=(sel5_vbf_jets[:,0]+sel5_vbf_jets[:,1]).mass,
                weight=sel5_weight
            )

            output['sel6_leptonpt_binning1'].fill(
                dataset=dataset,
                leptonpt=sel5_electrons.pt,
                weight=sel5_weight
            )

            output['sel6_leptonabseta_binning1'].fill(
                dataset=dataset,
                leptonabseta=abs(sel5_electrons.eta),
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
    samples = {
        'singlemuon' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/singlemuon.txt', 'frac' : 1},
        'singleelectron' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/singleelectron.txt', 'frac' : 1},
        'w' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/w.txt', 'frac' : 1},
        'ww' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/ww.txt', 'frac' : 1},
        'ewkwhjj': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/ewkwhjj.txt', 'frac' : 1},
        'qcdwhjj': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/qcdwhjj.txt', 'frac' : 1},
        'ttsemi': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/ttsemi.txt', 'frac' : 1},
        'tthad': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2016/tthad.txt', 'frac' : 1},
    }
elif year == '2017':
    samples = {
        'singlemuon' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2017/singlemuon.txt', 'frac' : 1},
        'singleelectron' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2017/singleelectron.txt', 'frac' : 1},
        'ttsemi' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2017/ttsemi.txt', 'frac' : 1},
    }
elif year == '2018':
    samples = {
        'singlemuon' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/singlemuon.txt', 'frac' : 1},
        'egamma' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/egamma.txt', 'frac' : 1},
        'ewkwhjj_reweighted': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/ewkwhjj_reweighted.txt', 'frac' : 1},
        'ewkwhjj': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/ewkwhjj_reweighted.txt', 'frac' : 1},
        'ttsemi': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/ttsemi.txt', 'frac' : 1},
        'tthad': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/tthad.txt', 'frac' : 1},
        #        'qcdwphjj': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/qcdwphjj.txt', 'frac' : 1},
        #        'qcdwmhjj': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/qcdwmhjj.txt', 'frac' : 1},
        'qcdwph': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/qcdwph.txt', 'frac' : 1},
        'qcdwmh': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/qcdwmh.txt', 'frac' : 1},
        #        'wlepbb' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/wbb.txt', 'frac' : 1},
        'wlep' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/wlep.txt', 'frac' : 1},
        #        'wlep2j': {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/wlep2j.txt', 'frac' : 1},
        'stoptchan' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/stoptchan.txt', 'frac' : 1},
        'santitoptchan' : {'filelist' : '/afs/cern.ch/user/a/amlevin/ewkwhjj/filelists/2018/santitoptchan.txt', 'frac' : 1},

    }
else:
    assert(0)

for sample in samples:
    f = open(samples[sample]['filelist'])
    frac = samples[sample]['frac']
    samples[sample] = f.read().rstrip('\n').split('\n')
    samples[sample] = samples[sample][0:int(frac*len(samples[sample]))]

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

save(result,'outfile_{}'.format(year))
