#!/usr/bin/env python

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.workArea = 'crab'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = ''

config.Data.outputPrimaryDataset = ''
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 100
NJOBS = 100  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/group/lpclonglived/DisappTrks/'
config.Data.publication = True
config.Data.outputDatasetTag = 'RunIIFall17DRPremix-93X_mc2017_realistic_v3-v1'

config.Site.storageSite = 'T3_US_FNALLPC'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    def forkAndSubmit(config):
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    reallySubmitMass = { x : False for x in range(100, 1000, 100)}
    reallySubmitLifetime = { x : False for x in [1, 10, 100, 1000, 10000]}
    numJobsPerLifetime = { x : (2500 if x == 10000 else 500) for x in [1, 10, 100, 1000, 10000]}

    for mass in range(100, 1000, 100):
        for ctau in [1, 10, 100, 1000, 10000]:
            config.General.requestName = 'AMSB_chargino%dGeV_ctau%dcm_step1' % (mass, ctau)
            config.JobType.psetName = 'step1/pythia8Decay/AMSB_chargino%dGeV_ctau%dcm_step1.py' % (mass, ctau)
            config.Data.outputPrimaryDataset = 'AMSB_chargino_M-%d_CTau-%d_TuneCP5_13TeV_pythia8' % (mass, ctau)
            config.Data.totalUnits = config.Data.unitsPerJob * numJobsPerLifetime[ctau]
            if reallySubmitMass[mass] and reallySubmitLifetime[ctau]:
                forkAndSubmit(config)
            else:
                print 'Skipping submission of request:', config.General.requestName
