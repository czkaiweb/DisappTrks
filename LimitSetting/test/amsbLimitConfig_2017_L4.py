#!/usr/bin/env python

# Local options file to be used with makeDataCards.py
# Usage:
# > makeDatacards.py -l amsbLimitConfig.py -c test
#
# Copied from https://raw.github.com/DisplacedSUSY/DisplacedSUSY/master/LimitsCalculation/test/sampleLimitConfig.py

from DisappTrks.SignalMC.signalCrossSecs import *
from DisappTrks.StandardAnalysis.plotUtilities import *
from DisappTrks.StandardAnalysis.IntegratedLuminosity_cff import *

dirs = getUser ()

nLayersWord = 'NLayers4'

exec('from amsbLimitConfigBkgds_2017_' + nLayersWord + ' import *')

##################################
### Event Selection Parameters ###
##################################

#name of histogram to integrate to get yields
integrateHistogramName = "metPt"

#########################
### Signal Parameters ###
#########################

# a separate datacard will be produced with each value of MASS,TAU
# named "datacard_AMSB_mGravMASSK_TAUns.txt"

samplesByGravitinoMass = False

#NOTE: These are the chargino masses in GeV
masses = ['100', '200', '300', '400', '500', '600', '700', '800', '900']

#chargino tau values in cm
lifetimes = ['2', '3', '4', '5', '6', '7', '8', '9', '10',
             '20', '30', '40', '50', '60', '70', '80', '90', '100',
             '200', '300', '400', '500', '600', '700', '800', '900', '1000',
             '2000', '3000', '4000', '5000', '6000', '7000', '8000', '9000', '10000']

lumi = lumi["MET_2017"]

#condor directory in which to find signal root files
signal_condor_dir = dirs["Brian"] + '/2017/signalAcceptance_full/'

signal_suffix = signal_suffix_in_datacard = '94X'

#name of event selection from which to take signal yields
signal_channel = 'disTrkSelectionSmearedJets' + nLayersWord + 'Plotter/Met Plots'

#######################
### Data Parameters ###
#######################

#this just sets the observed number of events equal to the total background expectation
run_blind_limits = False

data_dataset = "MET_2017"

#condor directory in which to find data root file
data_condor_dir = dirs["Andrew"] + '/2016_final_prompt/disappearingTracks/'

#name of event selection from which to take observed events
data_channel = 'DisTrkSelectionPlotter/Met Plots'

actual_bin_name = 'Bin2017' + nLayersWord

#############################
### Systematic Uncertainties ###
#############################

external_systematic_uncertainties = [
    # Use order of AN
    "isr_Bin2017" + nLayersWord,
    "jec_Bin2017" + nLayersWord,
    "jer_Bin2017" + nLayersWord,
    "metVaryElectronEn_Bin2017" + nLayersWord,
    "metVaryJetEn_Bin2017" + nLayersWord,
    "metVaryJetRes_Bin2017" + nLayersWord,
    "metVaryPhotonEn_Bin2017" + nLayersWord,
    "metVaryTauEn_Bin2017" + nLayersWord,
    "metVaryUnclusteredEn_Bin2017" + nLayersWord,
    "nMissOut_Bin2017" + nLayersWord,
    "pileup_Bin2017" + nLayersWord,
    "trigger_grandOrWeightData_Bin2017" + nLayersWord,
    "trigger_grandOrWeightMC_Bin2017" + nLayersWord,
]

#uncertainties on signal only (we can alter this if we need to)
signal_systematic_uncertainties = {
    'lumi_Bin2017' + nLayersWord :  {
        'value' : '1.023',
    },
    'trkReco_Bin2017' + nLayersWord :  {
        'value' : '1.015',
    },
    'Ecalo_Bin2017' + nLayersWord : {
        'value' : str (1.0 + 0.626555065492 / 100.0),
    },
    'Nmissin_Bin2017' + nLayersWord :  {
        'value' : str (1.0 + 1.16347975146 / 100.0),
    },
    'Nmissmid_Bin2017' + nLayersWord :  {
        'value' : str (1.0 + 0.26902546109 / 100.0),
    },
}
