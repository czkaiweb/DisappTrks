import FWCore.ParameterSet.Config as cms
import glob, sys, os

# The following are needed for the calculation of associated calorimeter energy
from Configuration.StandardSequences.GeometryRecoDB_cff import *
from Configuration.StandardSequences.MagneticField_38T_cff import *

dirName = sys.argv[2]
fileName = sys.argv[3]
use2018D = False
probeType = track # Options: track/electron/muon
print "processing \"" + dirName + "/*.root\"..."
print "writing histograms to \"" + fileName + "\""

###########################################################
##### Set up process #####
###########################################################

process = cms.Process ('PROBEANALYZER')
process.load ('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

from DisappTrks.StandardAnalysis.protoConfig_cfg import *
if use2018D:
  data_global_tag = '102X_dataRun2_Prompt_v13'
process.GlobalTag = GlobalTag(process.GlobalTag, data_global_tag, '')

process.maxEvents = cms.untracked.PSet (
    input = cms.untracked.int32 (100)
)
process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring (
        map (lambda a : "file:" + a, [f for f in glob.glob (dirName + "/*.root") if os.path.getsize(f) != 0])
    ),
)
process.TFileService = cms.Service ('TFileService',
    fileName = cms.string (fileName)
)

if probeType == track:
  probeCollection = "generalTracks"
  probeAnalyzerName   = "TrackProbeAnalyzer"
if probeType == electron:
  probeCollection = "gedGsfElectrons"
  probeAnalyzerName   = "ElectronProbeAnlyzer"
if probeType == muon:
  probeCollection = "muons"
  probeAnalyzerName   = "MuonProbeAnlyzer"

###########################################################
##### Set up the producer and the end path            #####
###########################################################

process.probeAnalyzer = cms.EDAnalyzer (probeAnalyzerName,
    probeTracks        =  cms.InputTag  ( probeCollection, ""),
    EBRecHits          =  cms.InputTag  ("reducedEcalRecHitsEB"),
    EERecHits          =  cms.InputTag  ("reducedEcalRecHitsEE"),
    HBHERecHits        =  cms.InputTag  ("reducedHcalRecHits", "hbhereco"),
    x_lo = cms.double (-0.5),
    x_hi = cms.double (0.5),
    n_x = cms.int32 (50),
    y_lo = cms.double (-0.5),
    y_hi = cms.double (0.5),
    n_y = cms.int32 (50),
)

process.myPath = cms.Path (process.ProbeAnalyzer)
