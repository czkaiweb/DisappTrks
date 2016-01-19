import FWCore.ParameterSet.Config as cms
from OSUT3Analysis.Configuration.processingUtilities import *
import math
import os

################################################################################
##### Set up the 'process' object ##############################################
################################################################################

process = cms.Process ('OSUAnalysis')

# how often to print a log message
process.load ('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# input source when running interactively
# ---------------------------------------
process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring (
        "/store/user/ahart/AMSB_chargino500GeV_ctau100cm_step4_User.root",  # signal
        # "file:condor/metMinimalSkim/AMSB_chargino_500GeV_100cm/metMinimalSkim/skim_0.root",
        # "file:condor/candTrkSelection/WJetsToLNu/CandTrkSelection/merged/skim_merged.root",
        # "/store/user/ahart/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-DisappearingTracks-v1/151009_193826/0001/miniAODWithCandidateTracks_1001.root", # bkgd MC
        # "/store/user/ahart/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-DisappearingTracks-v1/151009_193736/0000/miniAODWithCandidateTracks_1.root", # bkgd MC
        # '/store/user/wulsin/MET/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-DisappearingTracks-v1/151011_155240/0000/miniAODWithCandidateTracks_10.root', # data
    ),
)

# process.source.eventsToProcess = cms.untracked.VEventRange (
#     "1:60:53",
# )



# FIXME:  set_input does not work (because of error with /usr/bin/file) in CMSSW_7_4_5_ROOT5
# argument can be a ROOT file, directory, or dataset name*
# *registered dataset names are listed in 'datasets' in:
#    https://github.com/OSU-CMS/OSUT3Analysis/blob/master/Configuration/python/configurationOptions.py

# sample direcotory
# set_input(process, "/store/user/ahart/BN_stopToBottom_M_800_10mm_Tune4C_8TeV_pythia8_lantonel-Summer12_DR53X-PU_S10_START53_V19-v1-ab45720b22c4f98257a2f100c39d504b_USER_1/")

# sample ROOT file
#set_input(process, "/store/user/ahart/BN_stopToBottom_M_800_10mm_Tune4C_8TeV_pythia8_lantonel-Summer12_DR53X-PU_S10_START53_V19-v1-ab45720b22c4f98257a2f100c39d504b_USER_1/stopToBottom_M_800_10mm_Tune4C_8TeV_pythia8_lantonel-Summer12_DR53X-PU_S10_START53_V19-v1-ab45720b22c4f98257a2f100c39d504b_USER_10_2_Dzw.root")

# sample dataset nickname
#set_input(process, "DYToTauTau_20")
#set_input(process, "DYToMuMu_20")

# output histogram file name when running interactively
process.TFileService = cms.Service ('TFileService',
    fileName = cms.string ('hist.root')
)

# number of events to process when running interactively
process.maxEvents = cms.untracked.PSet (
    input = cms.untracked.int32 (100)
)

################################################################################
##### Set up the 'collections' map #############################################
################################################################################

from OSUT3Analysis.AnaTools.osuAnalysis_cfi import collectionMap  # miniAOD

# collectionMap.tracks = cms.InputTag ('candidateDisappearingTracks')  # For signal (old version)
collectionMap.tracks = cms.InputTag ('candidateTrackProducer')       # For data skim and new samples

################################################################################
##### Set up weights to be used in plotting and cutflows  ######################
################################################################################
weights = cms.VPSet (
    # cms.PSet (
    #     inputCollections = cms.vstring("muons"),
    #     inputVariable = cms.string("pt")
    # ),
)


################################################################################
##### Set up any user-defined variable producers ###############################
################################################################################

variableProducers = []
#variableProducers.append("MyVariableProducer")

################################################################################
##### Import the channels to be run ############################################
################################################################################

from DisappTrks.StandardAnalysis.EventSelections import *
from DisappTrks.StandardAnalysis.ElectronTagProbeSelections import *
from DisappTrks.StandardAnalysis.MuonTagProbeSelections import *

################################################################################
##### Import the histograms to be plotted ######################################
################################################################################

from DisappTrks.StandardAnalysis.HistogramDefinitions import *
from OSUT3Analysis.Configuration.histogramDefinitions import *

histSets = cms.VPSet (
    TrackHistograms,
    TrackBeamspotHistograms,
    TrackExtraHistograms,
    MetHistograms,
    JetHistograms
)

histSetsMetJet = cms.VPSet (
    MetHistograms,
    JetHistograms
)

histSetsMuon = copy.deepcopy(histSets)
histSetsMuon.append(MuonHistograms)
histSetsMuon.append(DiMuonHistograms)
histSetsMuon.append(TrackMuonHistograms)

histSetsElectron = copy.deepcopy(histSets)
histSetsElectron.append(ElectronHistograms)
histSetsElectron.append(DiElectronHistograms)
histSetsElectron.append(TrackElectronHistograms)

test = cms.PSet(
    name = cms.string("test"),
    triggers = cms.vstring(),
    cuts = cms.VPSet (
        cms.PSet (
            inputCollection = cms.vstring("primaryvertexs"),
            cutString = cms.string("isValid > 0 && ndof >= 4"),
            numberRequired = cms.string(">= 1")
        )
    )
)


################################################################################
##### Attach the channels and histograms to the process ########################
################################################################################

#  add_channels  (process,  [test],  cms.VPSet(),  weights,  collectionMap,  variableProducers,  True)

#  add_channels  (process,  [metMinimalSkim],         histSetsMetJet,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [nonIsoTrkSelection],     histSets,        weights,  collectionMap,  variableProducers,  False)
#  add_channels  (process,  [isoTrkSelection],        histSets,        weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [candTrkSelection],       histSets,        weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [disTrkSelection],        histSets,        weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [elecCtrlSelection],      histSets,        weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [muonCtrlSelection],      histSets,        weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [tauCtrlSelection],       histSets,        weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [candTrkEcaloSdband],     histSets,        weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [candTrkNMissOutSdband],  histSets,        weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [disTrkSelectionNHits3],  histSets,        weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [disTrkSelectionNHits4],  histSets,        weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [disTrkSelectionNHits5],  histSets,        weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [disTrkSelectionNHits6],  histSets,        weights,  collectionMap,  variableProducers,  True)

#  add_channels  (process,  [ZtoMuMu],                       histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuMu],                       histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuMuCandTrk],                histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuMuDisTrk],                 cms.VPSet      (),       weights,        collectionMap,      variableProducers,  True)
#  add_channels  (process,  [ZtoMuMuCandTrkEcaloSdband],     histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuMuCandTrkNMissOutSdband],  histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuProbeTrkNoMassCut],        cms.VPSet      (),       weights,        collectionMap,      variableProducers,  True)
#  add_channels  (process,  [ZtoMuMuDisTrkNHits3],           histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuMuDisTrkNHits4],           histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuMuDisTrkNHits5],           histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuMuDisTrkNHits6],           histSetsMuon,  weights,  collectionMap,  variableProducers,  True)

#  add_channels  (process,  [ZtoMuProbeTrk],                                histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuProbeTrkWithZCuts],                       histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuDisTrk],                                  histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuProbeTrkNoMissingOuterHitsCut],           histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuProbeTrkWithZCutsNoMissingOuterHitsCut],  histSetsMuon,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoMuDisTrkNoMissingOuterHitsCut],             histSetsMuon,  weights,  collectionMap,  variableProducers,  True)

#  add_channels  (process,  [ZtoEleProbeTrk],                                histSetsElectron,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoEleProbeTrkWithZCuts],                       histSetsElectron,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoEleDisTrk],                                  histSetsElectron,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoEleProbeTrkNoMissingOuterHitsCut],           histSetsElectron,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoEleProbeTrkWithZCutsNoMissingOuterHitsCut],  histSetsElectron,  weights,  collectionMap,  variableProducers,  True)
#  add_channels  (process,  [ZtoEleDisTrkNoMissingOuterHitsCut],             histSetsElectron,  weights,  collectionMap,  variableProducers,  True)

#add_channels  (process,  [test],  cms.VPSet  (),  weights,  collectionMap,  variableProducers,  False)

# uncomment to produce a full python configuration log file
#outfile = open('dumpedConfig.py','w'); print >> outfile,process.dumpPython(); outfile.close()

#process.Tracer = cms.Service("Tracer")

