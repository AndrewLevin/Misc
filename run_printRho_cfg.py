import FWCore.ParameterSet.Config as cms

process = cms.Process("DEMO")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START44_V1::All'   # CMSSW_4XY
##process.GlobalTag.globaltag = 'START52_V4::All'    # CMSSW_52Y
print "hi"
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    #'rfio:/castor/cern.ch/user/b/benedet/Fall11_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_0.root'  #CMSSW_4XY
    #'root://eoscms//eos/cms//store/relval/CMSSW_5_2_2/RelValZEE/GEN-SIM-RECO//START52_V4-v3/0004/CA611CFC-9074-E111-A583-003048CF94A6.root'   # CMSSW_52Y
     'file:/data/blue/khahn/AOD/HToZZTo4L_M-120_Fall11S6.00215E21D5C4.root'
    ),
    secondaryFileNames = cms.untracked.vstring(),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
    )
print "hi2"

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )

process.kt6PFJetsForIso1 = process.kt6PFJets.clone( Rho_EtaMax = cms.double(2.5), Ghost_EtaMax = cms.double(3.1) )
process.kt6PFJetsForIso2 = process.kt6PFJets.clone( Rho_EtaMax = cms.double(2.5), Ghost_EtaMax = cms.double(2.5) )

process.demo = cms.EDAnalyzer("printRho")

process.pAna = cms.Path(process.kt6PFJetsForIso1 * process.kt6PFJetsForIso2 * process.demo)
print "hi3"

process.schedule = cms.Schedule(process.pAna)
