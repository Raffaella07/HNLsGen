'''
Job option for the B-initiated HNL generation
'''

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')
# define the defaults here, changed from command line
options.maxEvents = -1 # -1 means all events, maxEvents considers the total over files considered
options.outputFile = 'BPH-test.root'
# add costum parameters
options.register ('severityLevel',
                  'ERROR', # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.string,          # string, int, or float
                  'severity level for log messages, DEBUG, INFO, WARNING, ERROR')
options.register('nThr',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 'Number of threads')
options.register('seedOffset',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 'Seed offset')
options.register('mass',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                 'mass of the HNL')
options.register('ctau',
                 100,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                 'ctau of the HNL [mm]')

#options.register ("doDirac",
#                  1, # default value
#                  VarParsing.multiplicity.singleton, # singleton or list
#                  VarParsing.varType.int,          # string, int, or float
#                  "do Dirac HNL? otherwise Majorana")
options.parseArguments()
print options

import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeVEarly2018Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('test'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(1),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Message logger
#process.MessageLogger = cms.Service("MessageLogger",
#    cout = cms.untracked.PSet(
#         threshold  = cms.untracked.string(options.severityLevel) 
#    ),
#    destinations = cms.untracked.vstring('cout')
#)

# Other statements
process.XMLFromDBSource.label = cms.string("Extended")
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v11', '')

process.MuFilter = cms.EDFilter("MCParticlePairFilter",
    MaxEta = cms.untracked.vdouble(2.45, 2.45),
    MinEta = cms.untracked.vdouble(-2.45, -2.45),
    MinPt = cms.untracked.vdouble(2.7, 2.7),
    ParticleID1 = cms.untracked.vint32(13),
    ParticleID2 = cms.untracked.vint32(13)
)


### Operates on all particles in the HepMC::GenEvent
### accpects events if:
###  - there is at least one particle with specified pdgID in the entire HepMC::GenEvent
###  - any status (but can be specified)
process.BpFilter = cms.EDFilter("PythiaFilter",
    #ParticleID = cms.untracked.int32(14) # nu_mu anti_numu filter
    ParticleID = cms.untracked.int32(521) # B+ B- filter 
)

#process.SingleMuFilter = cms.EDFilter("PythiaFilter", # using PythiaFilter instead of MCParticleFilter because the particleID is taken in abs value
#    MaxEta = cms.untracked.double(1.6),
#    MinEta = cms.untracked.double(-1.6),
#    MinPt = cms.untracked.double(5), # <=== keep it a bit lower than the pt cut at reco level... 
#    ParticleID = cms.untracked.int32(13), # abs value is taken
#    #Status = cms.untracked.int32(1), # new
#    MotherID = cms.untracked.int32(521), # require muon to come from B+/B- decay
#)
process.SingleMuFilter = cms.EDFilter("PythiaFilterMotherSister", 
    #MaxEta = cms.untracked.double(3),
    #MinEta = cms.untracked.double(-3),
    #MinPt = cms.untracked.double(0.5), # <=== keep it a bit lower than the pt cut at reco level... 
    MaxEta = cms.untracked.double(1.6),
    MinEta = cms.untracked.double(-1.6),
    MinPt = cms.untracked.double(5), # <=== keep it a bit lower than the pt cut at reco level... 
    ParticleID = cms.untracked.int32(13), # abs value is taken
    #Status = cms.untracked.int32(1),
    MotherID = cms.untracked.int32(521), # require muon to come from B+/B- decay
    SisterID = cms.untracked.int32(9900015), # require HNL sister
)


process.generator = cms.EDFilter("Pythia8GeneratorFilter",
    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            ### for info, see https://twiki.cern.ch/twiki/bin/view/CMS/EvtGenInterface
            convertPythiaCodes = cms.untracked.bool(False),
            
            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
            
            ### the list of particles that are aliased and forced to be decayed by EvtGen
            list_forced_decays = cms.vstring(       
                'myB+', 
                'myB-'
            ),
            
            ### the list of particles that remain undecayed by Pythia for EvtGen to operate on. 
            ### If the vector has a size 0 or size of 1 with a value of 0, the default list is used. 
            ### These are are hard-coded in: GeneratorInterface/EvtGenInterface/plugins/EvtGen/EvtGenInterface.cc., in the function SetDefault_m_PDGs().            
            operates_on_particles = cms.vint32(521, -521), 

            ### The file with properties of all particles
            particle_property_file = cms.FileInPath('HNLsGen/evtGenData/evt_2014_mass{m}_ctau{ctau}.pdl'.format(m=options.mass,ctau=options.ctau)), 
            #particle_property_file = cms.FileInPath('HNLsGen/evtGenData/evt_2014_mod.pdl'),
            #particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_2014_mod.pdl'),  
            #https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEdmFileInPath


            ### Decay chain
            #user_decay_embedded = cms.vstring("\nAlias myB+ B+\nAlias myB- B-\nAlias mytau+ tau+\nAlias mytau- tau-\nChargeConj myB+ myB-\nChargeConj mytau+ mytau-\n\nDecay myB-\n0.259     anti-D0       mytau-     nu_tau    ISGW2;\n0.592     anti-D*0      mytau-     nu_tau    ISGW2;\n0.074     anti-D_2*0    mytau-     nu_tau    ISGW2;\n0.074     anti-D\'_10    mytau-     nu_tau    ISGW2;\nEnddecay\nCDecay myB+\n\nDecay mytau-\n1.0 mu-    mu+    mu-             PHOTOS PHSP;\nEnddecay\nCDecay mytau+\n\nEnd\n")
            
            # decay to neutrino
            #user_decay_embedded = cms.vstring("\nAlias myB+ B+\nAlias myB- B-\nAlias myD0 D0\nAlias myAntiD0 anti-D0\nChargeConj myB+ myB-\nChargeConj myD0 myAntiD0\nDecay myB-\n1.0     myD0    mu-    anti-nu_mu    ISGW2;\nEnddecay\nCDecay myB+\nDecay myD0\n1.0    K-    pi+    PHSP;\nEnddecay\nCDecay myAntiD0\nEnd\n")
            

            # decay to HNL, majorana, no need to define charge conjugation
            #user_decay_embedded = cms.vstring("\nAlias myB+ B+\nAlias myB- B-\nAlias myD0 D0\nAlias myAntiD0 anti-D0\nChargeConj myB+ myB-\nChargeConj myD0 myAntiD0\nDecay myB-\n1.0     myD0    mu-    hnl    PHSP;\nEnddecay\nCDecay myB+\n\nDecay myD0\n1.0    K-    pi+    PHSP;\nEnddecay\nCDecay myAntiD0\n\nEnd\n")
            
      
            # decay to HNL, dirac: charge conjugate of hnl -> anti_hnl
            #user_decay_embedded = cms.vstring("\nAlias myB+ B+\nAlias myB- B-\nAlias myD0 D0\nAlias myAntiD0 anti-D0\nChargeConj myB+ myB-\nChargeConj myD0 myAntiD0\nChargeConj hnl anti_hnl\nDecay myB-\n1.0     myD0    mu-    anti_hnl    PHSP;\nEnddecay\nCDecay myB+\n\nDecay myD0\n1.0    K-    pi+    PHSP;\nEnddecay\nCDecay myAntiD0\n\nEnd\n")
            # decay to HNL, dirac, and decay also HNL 
            user_decay_embedded = cms.vstring("\nAlias myB+ B+\nAlias myB- B-\nAlias myD0 D0\nAlias myAntiD0 anti-D0\nChargeConj myB+ myB-\nChargeConj myD0 myAntiD0\nChargeConj hnl anti_hnl\nDecay myB-\n1.0     myD0    mu-    anti_hnl    PHSP;\nEnddecay\nCDecay myB+\n\nDecay myD0\n1.0    K-    pi+    PHSP;\nEnddecay\nCDecay myAntiD0\nDecay anti_hnl\n1.0     mu+    pi-    PHSP;\nEnddecay\nCDecay hnl\n\nEnd\n")

        ),
        parameterSets = cms.vstring('EvtGen130')
    ),
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring(
            'pythia8CommonSettings', 
            'pythia8CUEP8M1Settings',  # pythia8CP5Settings ?
                                       # pythia8PSweightsSettings ?
            'processParameters'
        ),
        processParameters = cms.vstring(
            ## 'SoftQCD' vs 'HardQCD' 
            ##     you want SoftQCD if you don#'t want to put any pT cut on the hard scatter process 
            ##     http://home.thep.lu.se/~torbjorn/pythia81html/QCDProcesses.html
            ##     eventually use SoftQCD if you#'re interested in the full bottom production at high energies

            ### original settings for tau->3mu   
            #'SoftQCD:nonDiffractive = on',            # default is off     
            #'SoftQCD:singleDiffractive = on',         # default is off
            #'SoftQCD:doubleDiffractive = on',         # default is off
            #'PTFilter:filter = on',                   # default is off  # could not find **ANYWHERE** in the Pythia code PTFilter 
            #'PTFilter:quarkToFilter = 5',                               # it's something that exists in CMSSW only, see Py8InterfaceBase.cc
            #'PTFilter:scaleToFilter = 1.0'            # default is 0.4 
           
            ### settings to generate bbar only as per tip https://twiki.cern.ch/twiki/bin/view/CMS/EvtGenInterface#Tips_for_Pythia8   
            'SoftQCD:nonDiffractive = off',            # 
            'SoftQCD:singleDiffractive = off',         #
            'SoftQCD:doubleDiffractive = off',         #
            'PTFilter:filter = off',                   #
            'HardQCD:gg2bbbar = on ',                  # default is off 
            'HardQCD:qqbar2bbbar = on ',               # default is off  
            'HardQCD:hardbbbar = off',                 # default is off  # should be set to off if gg2bbbar and hardbbbar on, otherwise double-counting
            'PhaseSpace:pTHatMin = 5.',               # default is 0    # minimum invariant pT
            ## 'PhaseSpace' to constrain the kinematics of a 2->2 process, 
            ##              for hard physics only, 
            ##              in the rest frame of the hard process, 
            ##              cross-section is adjusted to correspond for the allowed phase-space
        ),
        pythia8CUEP8M1Settings = cms.vstring( # these probably remain the same
            'Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:pT0Ref=2.4024',    # default is 2.28000 
            'MultipartonInteractions:ecmPow=0.25208',   # default is 0.21500
            'MultipartonInteractions:expPow=1.6'        # default is 1.85000
        ),
        pythia8CommonSettings = cms.vstring( 
            'Tune:preferLHAPDF = 2',                    # default is 1 
            'Main:timesAllowErrors = 10000',            # default is 10
            'Check:epTolErr = 0.01',                    # default is 1.0000e-04
            'Beams:setProductionScalesFromLHEF = off', 
            'SLHA:keepSM = on',                         # default is 100
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on',            # default is false
            'ParticleDecays:tau0Max = 10',   
            'ParticleDecays:allowPhotonRadiation = on'  # default is false
        ) 
        # do we want pythia8PSweightsSettings ?
    ),
    comEnergy = cms.double(13000.0),
    filterEfficiency = cms.untracked.double(0.0013),  # this will not be used by Pythia, only saved in GenInfo
    maxEventsToPrint = cms.untracked.int32(0),        # max events to print the complete event list information
    pythiaHepMCVerbosity = cms.untracked.bool(False), # to display HepMC information: vertices and particles (not interesting)
    pythiaPylistVerbosity = cms.untracked.int32(1)    # 1 for "normal" verbosity, 11 to display all Pythia Settings
)


process.ProductionFilterSequence = cms.Sequence(process.generator+process.BpFilter+process.SingleMuFilter)
#process.ProductionFilterSequence = cms.Sequence(process.generator+process.BpFilter) 

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(options.nThr)
process.options.numberOfStreams=cms.untracked.uint32(0)

# set a different offset seed, if you run multiple jobs 
process.RandomNumberGeneratorService.eventSeedOffset=cms.untracked.uint32(options.seedOffset)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions

# Customisation from command line
process.MessageLogger.cerr.FwkReport.reportEvery=100

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
