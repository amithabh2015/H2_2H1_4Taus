import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeProducer")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.GlobalTag.globaltag = 'MCRUN2_74_V9'
process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v2'
#process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v1'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

# MINIAOD

# ZZ_TuneCUETP8M1
#'/store/mc/RunIISpring15DR74/ZZ_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0A03D5CB-1609-E511-9D94-0025904CF710.root',
#'/store/mc/RunIISpring15DR74/ZZ_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0C1BB328-0D09-E511-8029-0025905C2CBE.root',
#'/store/mc/RunIISpring15DR74/ZZ_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0C479546-7209-E511-BA6A-3417EBE8862E.root'

# WZTo3LNu
#'/store/mc/RunIISpring15DR74/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/008E7FBF-9218-E511-81E0-001E675A5244.root',
#'/store/mc/RunIISpring15DR74/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/0473AA1C-AE18-E511-A22A-A0040420FE80.root',
#'/store/mc/RunIISpring15DR74/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/0C4172B8-9218-E511-B9C9-001E675A58D9.root',
#'/store/mc/RunIISpring15DR74/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/1299A85E-9E18-E511-AE6C-A0040420FE80.root'

# WJetsToLNu
#        '/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02215B44-2D70-E411-90A3-0025905A60B8.root',
#        '/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0603D444-2D70-E411-AF03-002618943922.root',
#        '/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/08947C88-3570-E411-974E-002618FDA26D.root',
#        '/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0E8B81D9-2D70-E411-94AB-0025905A4964.root',
#        '/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/1225D443-2D70-E411-9D85-0025905B85F6.root'

# DYJetsToLL_M-50
#'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/00C4781D-6B08-E511-8A0A-0025905A6084.root',
#'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/02DE3B74-6C08-E511-ABE3-0025905A60D0.root',

# TTJets
#        '/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root',
#        '/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00D3EAF1-3174-E411-A5B2-0025904B144E.root',
#        '/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02EF3EFC-0475-E411-A9DB-002590DB9166.root'
#
# QCD_Pt-50to80  
#        '/store/mc/Phys14DR/QCD_Pt-50to80_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/00000/164247BB-6577-E411-A063-00259073E50A.root',
#        '/store/mc/Phys14DR/QCD_Pt-50to80_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/00000/1EF05EBB-6577-E411-BD3A-20CF305616EF.root',
#        '/store/mc/Phys14DR/QCD_Pt-50to80_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/00000/24BF07BB-6577-E411-89E4-00259073E3AE.root',
#        '/store/mc/Phys14DR/QCD_Pt-50to80_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/00000/303020BC-6577-E411-A7A7-00259073E516.root'

# TT_TuneCUETP8M1_13TeV
#'/store/mc/RunIISpring15DR74/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v4/10000/00D2A247-2910-E511-9F3D-0CC47A4DEDD2.root',
#'/store/mc/RunIISpring15DR74/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v4/10000/024F2046-2910-E511-A1A5-0CC47A4DEE2A.root',
#'/store/mc/RunIISpring15DR74/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v4/10000/02690835-3310-E511-8FDC-0025905A60AA.root'

# WW_TuneCUETP8M1_13TeV

# Signal Samples
#'file:Signal/H2ToH1H1_To4Tau_mH2_125_mH1_8_13TeV_MINIAODSIM_0.root',
#'file:Signal/H2ToH1H1_To4Tau_mH2_125_mH1_8_13TeV_MINIAODSIM_1.root',

#Data-2015B-DoubleMuon-MINIAOD-PromptReco-v1
#'/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/244/00000/E42FEF61-6E27-E511-B93A-02163E0143C0.root',
#'/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/252/00000/9ADEE140-9C27-E511-919A-02163E011D23.root'

# Data-2015C-DoubleMuon-MINIAOD-PromptReco-v1
#'/store/data/Run2015C/DoubleMuon/MINIAOD/PromptReco-v1/000/254/227/00000/78B24362-2A46-E511-AEFB-02163E014330.root'

# Data-2015D-DoubleMuon-MINIAOD-05Oct2015-v1
'/store/data/Run2015D/DoubleMuon/MINIAOD/05Oct2015-v1/30000/067DA5FD-876F-E511-B2A1-0025905A60FE.root'
#    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
##    inputCommands=cms.untracked.vstring(
##        'keep *',
##        'drop patTaus_*_*_*',
##        'drop *PFTau*_*_*_*'
##    )
),
)
 
#####################################################
  
#--------------------------------------------------------------------------------
# produce PAT-tuple
process.load("PhysicsTools/PatAlgos/patSequences_cff")
# configure pat::Jet production
# (enable L2L3Residual corrections in case running on Data)
#jetCorrections = ( 'L1FastJet', 'L2Relative', 'L3Absolute')
#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJetCollection(
#    process,
#    jetSource = cms.InputTag('ak4PFJetsCHS'),
#    jetCorrections = ( 'AK4PFchs', jetCorrections, "" ),
#    outputModules = []
#)

#process.patJets.addTagInfos = cms.bool(True)
#process.patJets.tagInfoSources = cms.VInputTag("impactParameterTagInfosAOD","secondaryVertexTagInfosAOD","softMuonTagInfosAOD")
#process.patJets.tagInfoSources = cms.VInputTag("secondaryVertexTagInfosEI")


#--------------------------------------------------------------------------------
# switch to HPS PFTaus (and disable all "cleaning" cuts)
#from PhysicsTools.PatAlgos.tools.tauTools import *
#switchToPFTauHPS(process)

##################################################
# Main
process.NTupleProducer = cms.EDAnalyzer("NTupleProducerFromMiniAOD",
# data, year, period, skim
IsData = cms.untracked.bool(True),
# collections
MuonCollectionTag = cms.InputTag("slimmedMuons"), 
JetCollectionTag = cms.InputTag("slimmedJets"),
MetCollectionTag = cms.InputTag("slimmedMETs"),
TrackCollectionTag = cms.InputTag("packedPFCandidates"),
GenParticleCollectionTag = cms.InputTag("prunedGenParticles"),
TriggerObjectCollectionTag = cms.InputTag("selectedPatTrigger"),
BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
# trigger info
HLTriggerPaths = cms.untracked.vstring(
'HLT_Mu24_v',
'HLT_Mu24_eta2p1_v',
'HLT_Mu27_v',
'HLT_Mu34_v',
'HLT_Mu45_eta2p1_v',
'HLT_Mu17_Mu8_v',
'HLT_Mu17_Mu8_DZ_v',
'HLT_Mu17_Mu8_SameSign_v'
'HLT_Mu17_Mu8_SameSign_DZ_v',
'HLT_IsoMu18_v',
'HLT_IsoMu20_v',
'HLT_IsoMu17_eta2p1_v'
'HLT_IsoMu24_eta2p1_v'
),
TriggerProcess = cms.untracked.string("HLT"),
# tracks
RecTrackPtMin = cms.untracked.double(0.5),
RecTrackEtaMax = cms.untracked.double(2.4),
RecTrackNum = cms.untracked.int32(0),
# muons
RecMuonPtMin = cms.untracked.double(10.),
RecMuonEtaMax = cms.untracked.double(2.4),
RecMuonHLTriggerMatching = cms.untracked.vstring(
'HLT_Mu24_v.*:hltL3fL1sMu16L1f0L2f16L3Filtered24',
'HLT_Mu24_eta2p1_v.*:hltL3fL1sMu20Eta2p1L1f0L2f10QL3Filtered24Q',
'HLT_Mu27_v.*:hltL3fL1sMu25L1f0L2f10QL3Filtered27Q',
'HLT_Mu34_v.*:hltL3fL1sMu20L1f0L2f20L3Filtered34',
'HLT_Mu45_eta2p1_v.*:hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered45e2p1Q',
#'HLT_Mu45_eta2p1_v.*:hltL3fL1sMu16orMu25L1f0L2f16QL3Filtered45e2p1Q',
'HLT_Mu17_Mu8_v.*:hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_Mu8_v.*:hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_Mu8_DZ_v.*:hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_Mu8_DZ_v.*:hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_Mu8_DZ_v.*:hltDiMuonGlb17Glb8DzFiltered0p2',
'HLT_Mu17_Mu8_SameSign_v.*:hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_Mu8_SameSign_v.*:hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_Mu8_SameSign_v.*:hltDiMuonGlb17Glb8DzFiltered0p2',
'HLT_Mu17_Mu8_SameSign_v.*:hltDiMuonGlb17Glb8DzFiltered0p2SameSign' 
'HLT_Mu17_Mu8_SameSign_DZ_v.*:hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8',
'HLT_Mu17_Mu8_SameSign_DZ_v.*:hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17',
'HLT_Mu17_Mu8_SameSign_DZ_v.*:hltDiMuonGlb17Glb8DzFiltered0p2',
'HLT_Mu17_Mu8_SameSign_DZ_v.*:hltDiMuonGlb17Glb8DzFiltered0p2SameSign',
'HLT_IsoMu18_v.*:hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09',
'HLT_IsoMu20_v.*:hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09',
'HLT_IsoMu17_eta2p1_v.*:hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09'
'HLT_IsoMu24_eta2p1_v.*:hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09'
),
RecMuonNum = cms.untracked.int32(0),
# jets
RecJetPtMin = cms.untracked.double(30.),
RecJetEtaMax = cms.untracked.double(2.4),
RecJetHLTriggerMatching = cms.untracked.vstring(),
RecJetBtagDiscriminators = cms.untracked.vstring(
'jetBProbabilityBJetTags',
'jetProbabilityBJetTags',
'trackCountingHighPurBJetTags',
'trackCountingHighEffBJetTags',
'simpleSecondaryVertexHighEffBJetTags',
'simpleSecondaryVertexHighPurBJetTags',
'combinedInclusiveSecondaryVertexV2BJetTags',
'pfCombinedSecondaryVertexBJetTags',
'combinedMVABJetTags'
),
RecJetNum = cms.untracked.int32(0),
# sample name (Monte Carlo)
SampleName = cms.untracked.string("Data") 
)

#process.patJets.addBTagInfo = cms.bool(True)

process.p = cms.Path(process.NTupleProducer)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("Data2015D_DoubleMu_ntuple.root"
				  # fileName = cms.string("WZTo3LNu_ntuple_100000.root"
	)
                                   )

#processDumpFile = open('MyRootMaker.dump', 'w')
#print >> processDumpFile, process.dumpPython()

