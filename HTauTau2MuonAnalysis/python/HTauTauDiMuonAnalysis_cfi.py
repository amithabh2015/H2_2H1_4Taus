import FWCore.ParameterSet.Config as cms
import RecoMET.METProducers.METSigParams_cfi as met_config

HTauTauDiMuonDefaultsBlock = cms.PSet(
	datasetID = cms.int32(0),
	datasetNickname = cms.string(""),
	year = cms.int32(2014),
	period = cms.string("2014"),
	MuonSource = cms.InputTag("selectedPatMuons"),
	PatMETSource = cms.InputTag("patMETs"),
	#PatMETMVASource = cms.InputTag("patPFMetByMVA"),
	JetSource = cms.InputTag("selectedPatJets"),
	#PuJetIdSource = cms.InputTag("pileupJetIdProducer"),
	#PuJetMvaSource = cms.InputTag("pileupJetIdProducer"),
	RecoJetSource = cms.InputTag("ak4PFJets"),
	PhotonSource = cms.InputTag("selectedPatPhotons"),
	PrimaryVertexSource = cms.InputTag("offlinePrimaryVertices"),
	GenParticleSource = cms.InputTag("genParticles"),
	##   TriggerSource = cms.InputTag("patTriggerPFlow"),
	TriggerSource = cms.InputTag("TriggerResults"), # correct process name is automatically determined
	Triggers = cms.vstring(
		"Mu9",
		"Mu11",
		"Mu12",
		"Mu13",
		"Mu15",

		"Mu20",
		"Mu24",
		"Mu24_eta2p1",
		"Mu30",
		"Mu30_eta2p1",
		"Mu40",
		"Mu40_eta2p1",
		"Mu60_eta2p1",

		"IsoMu9",
		"IsoMu12",
		"IsoMu13",
		"IsoMu15",
		"IsoMu15_eta2p1",
		"IsoMu17",
		"IsoMu20",
		"IsoMu24",
		"IsoMu24_eta2p1",
		"IsoMu30_eta2p1",
		"IsoMu34_eta2p1",

		"DoubleMu3",
		"DoubleMu4_Acoplanarity03",
		"DoubleMu5",
		"DoubleMu5_Acoplanarity03",
#		"DoubleMu5_IsoMu5",
		"DoubleMu6",
		"DoubleMu6_Acoplanarity03",
		"DoubleMu7",
		"DoubleMu45",

		"Mu13_Mu8",
		"Mu17_Mu8",
                "Mu17_TkMu8",
	),
	TriggerMatching = cms.PSet(
		Objects = cms.vstring("PosMuon", "NegMuon"),
		Triggers = cms.VPSet(
			cms.PSet(
				TriggerName = cms.string("IsoMu17"),
				FilterModules = cms.VPSet(
					cms.PSet(
						Nick = cms.string("Leg17"),
						Modules = cms.vstring("hltSingleMuIsoL3IsoFiltered17") # 2011 menu
					)
				)
			),
			cms.PSet(
				TriggerName = cms.string("Mu24"),
				FilterModules = cms.VPSet(
					cms.PSet(
						Nick = cms.string("Leg24"),
						Modules = cms.vstring("hltSingleMu24L2QualL3Filtered24", # 2011 menu
						                      "hltL3fL1sMu16L1f0L2f16QL3Filtered24Q") # 2012 menu
					)
				)
			),
			cms.PSet(
				TriggerName = cms.string("Mu24_eta2p1"),
				FilterModules = cms.VPSet(
					cms.PSet(
						Nick = cms.string("Leg24"),
						Modules = cms.vstring("hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered24Q") # 2012 menu
					)
				)
			),
			cms.PSet(
				TriggerName = cms.string("DoubleMu7"),
				FilterModules = cms.VPSet(
					cms.PSet(
						Nick = cms.string("Leg7"),
						Modules = cms.vstring("hltDiMuonL3PreFiltered7") # 2011 menu
					)
				)
			),
			cms.PSet(
				TriggerName = cms.string("Mu13_Mu8"),
				FilterModules = cms.VPSet(
					cms.PSet(
						Nick = cms.string("Leg8"),
						Modules = cms.vstring("hltDiMuonL3PreFiltered8", # early 2011 menu
						                      "hltDiMuonL3p5PreFiltered8", # late 2011 menu
						                      "hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8", # early 2012 menu
						                      "hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8") # late 2012 menu
					),
					cms.PSet(
						Nick = cms.string("Leg13"),
						Modules = cms.vstring("hltSingleMu13L3Filtered13", # 2011 menu
						                      "hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered13", # early 2012 menu
						                      "hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered13") # late 2012 menu
					),
					cms.PSet(
						Nick = cms.string("DzFilter"),
						AllowNonExisting = cms.untracked.bool(True),
						Modules = cms.vstring("hltDiMuonMu13Mu8DzFiltered0p2", # early 2012 menu
						                      "hltDiMuonGlb13Glb8DzFiltered0p2") # late 2012 menu
					)
				)
			),
			cms.PSet(
				TriggerName = cms.string("Mu17_Mu8"),
				FilterModules = cms.VPSet(
					cms.PSet(
						Nick = cms.string("Leg8"),
						Modules = cms.vstring("hltDiMuonL3PreFiltered8",
						                      "hltDiMuonL3p5PreFiltered8",
						                      "hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8",
						                      "hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8"),
					),
					cms.PSet(
						Nick = cms.string("Leg17"),
						Modules = cms.vstring("hltSingleMu13L3Filtered17",
						                      "hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17",
						                      "hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17"),
					),
					cms.PSet(
						Nick = cms.string("DzFilter"),
						AllowNonExisting = cms.untracked.bool(True),
						Modules = cms.vstring("hltDiMuonMu17Mu8DzFiltered0p2",
						                      "hltDiMuonGlb17Glb8DzFiltered0p2")
					)
				)
			),
		)
	),
	PFParticleSource = cms.InputTag("particleFlow"),
	pfPileUpSource = cms.InputTag("pfPileUp"),
	pfNoPileUpSource = cms.InputTag("pfNoPileUp"),
	PtMinCut = cms.double(10.0),
	EtaMaxCut = cms.double(2.5),
	PtMinHardLeptonCut = cms.double(10.0),
	IsolationCut = cms.double(9999999.0),
	ZVertexCut = cms.double(24.0),
	DVertexCut = cms.double(3.0),
	ProbVertexCut = cms.double(0.005),
	NumberOfJetsCut = cms.int32(999999),
	JetEnergyMinCut = cms.double(10.0),
	METCut = cms.double(999999.0),
	DiLeptonDPhiCut =cms.double(-1.0),
	IsoDeltaR = cms.double(0.4),
	PrintOut = cms.int32(0),
	UseTriggerResults = cms.int32(0),
   L1MuonSource = cms.InputTag("l1extraParticles"),
)
