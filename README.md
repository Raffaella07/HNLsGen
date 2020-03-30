## Instructions to set up environment for HNL generation

First Installation
```
cmsrel CMSSW_10_2_3
cd CMSSW_10_2_3/src
cmsenv
git cms-init

git cms-addpkg GeneratorInterface/EvtGenInterface

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


## Drivers 
```
BPH_start_cfg.py                  => mod tau->3mu  with Fall18 
BPH_mod_cfg.py                    => tentative HNL with Fall18
```

## Run
```
cmsRun BPH_mod_cfg.py maxEvents=100 
```

test 
```
GeneratorInterface/ExternalDecays/test
cmsRun evtgentest_cfg.py 
```
