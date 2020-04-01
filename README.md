## Instructions to set up environment for HNL generation

First Installation
```
cmsrel CMSSW_10_2_3
cd CMSSW_10_2_3/src
cmsenv
git cms-init

git cms-addpkg GeneratorInterface/EvtGenInterface

# git cms-addpkg GeneratorInterface/ExternalDecays
# git cms-addpkg GeneratorInterface/Pythia8Interface

git clone git@github.com:BParkHNLs/HNLsGen.git

cp HNLsGen/evtGenData/evt_2014_mod.pdl GeneratorInterface/EvtGenInterface/data/.

export CMSSW_SEARCH_PATH=$CMSSW_BASE/src/GeneratorInterface/EvtGenInterface/data/:$CMSSW_SEARCH_PATH  # needed to use local evt.pdf file

scram b

git checkout -b mybranch

```

After first installation:
```
cd CMSSW_10_2_3/src
cmsenv
```
If you modify ```evtGenData/evt_2014_mod.pdl```, repeat the copy step above

## Instructions to set up a different version of Pythia within CMSSW
(WORK IN PROGRESS)

This needs to be started from clean CMSSW directory, before cmsenv

* On my laptop, cloned sonia's version of pythia
* compiled there (make)
* copied to t3
* then followed https://twiki.cern.ch/twiki/bin/viewauth/CMS/Pythia8Interface#How_to_setup_the_SCRAM_tool_with 
* however does not compile due to imcompatibility between versions 8.230 which is what is used by CMSSSW_10_2_0 and 8.240_sonia 

Therefore adopt different strategy:
* start from 8230 version (as obtained from pythia8 website)
* compile on laptop => but then removed in /lib the file
* copy it somewhere on t3
* copy there from sonia's pythia the files she added / changed
* the follow again above mentioend twiki page...



## Drivers 
```
BPH_start_cfg.py                  => mod tau->3mu  with Fall18 
BPH_mod_cfg.py                    => tentative HNL with Fall18
```

## Produce GEN-SIM
```
cd HNLsGen 
cmsRun cmsDrivers/BPH_mod_cfg.py maxEvents=100 outputFile=BPH-test.root
```

## Analyze
To visualize the decay chain in a tree (printout to screen), using ```vector<reco::genParticles>```
```
cd genLevelAnalysis
cmsRun test_ParticleTreeDrawer.py maxEvents=1 inputFiles=file:/work/mratti/GEN_HNL/CMSSW_10_2_3/src/HNLsGen/genSimFiles/BPH-test_HardQCDon.root
```

Proto-analyzer of ```edm::HepMCProduct```
```
cd genLevelAnalysis
cmsRun test_EvtGenTestAnalyzer.py
```

