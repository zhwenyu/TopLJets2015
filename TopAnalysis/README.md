# TopLJets2015

## Installation instructions

These installation instructions correspond to the 2017 data/MC production.
To install execute the following in your work area.
Notice: if you are not creating the ntuples, you can skip the part of the instructions 
marked with the `##OPTIONAL/##END OPTIONAL` markers

```
cmsrel CMSSW_9_4_2
cd CMSSW_9_4_2/src
cmsenv

##OPTIONAL

#photon/electron id+scale and smearing
git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP    
git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3
git cms-merge-topic cms-egamma:MiniAOD2017V2_940
scram b -j 8
cd $CMSSW_BASE/external
cd ${SCRAM_ARCH}/
git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
cd data/RecoEgamma/PhotonIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/external
cd ${SCRAM_ARCH}/
git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/external
cd ${SCRAM_ARCH}/
git clone git@github.com:Sam-Harper/EgammaAnalysis-ElectronTools.git data/EgammaAnalysis/ElectronTools/data
cd data/EgammaAnalysis/ElectronTools/data
git checkout ReReco17NovScaleAndSmearing 

##END OPTIONAL


cd $CMSSW_BASE/src
git clone git@github.com:pfs/TopLJets2015.git
cd TopLJets2015
git submodule init
git submodule update
scram b -j 8
```

## Running ntuple creation and checking the selection

The ntuplizer is steered with test/runMiniAnalyzer_cfg.py.
It takes several options from command line (see cfg for details).
To run locally the ntuplizer, for testing purposes do something like:

```
cmsRun test/runMiniAnalyzer_cfg.py runOnData=False era=era2017 outFilename=MC13TeV_TTJets.root
cmsRun test/runMiniAnalyzer_cfg.py runOnData=True  era=era2017 outFilename=Data13TeV_SingleMuon.root
```

To submit the ntuplizer to the grid start by setting the environment for crab3.
More details can be found in [CRAB3CheatSheet](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3CheatSheet#Environment_setup)

```
source /cvmfs/cms.cern.ch/crab3/crab.sh
```
The following script helps submitting a list of files described in a json file.
Partial submission can be made adding "-o csv_list" as an option.
Adding "-s" will trigger the submission to the grid (otherwise the script only writes down the crab cfg files)

```
python scripts/submitToGrid.py -j data/era2017/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py 
```

As soon as ntuple production starts to finish, to move from crab output directories to a simpler directory structure which can be easily parsed by the local analysis runThe merging can be run locally if needed by using the checkProductionIntegrity.py script

```
python scripts/submitCheckProductionIntegrity.py -i /store/cmst3/group/top/psilva/c29f431 -o /store/cmst3/group/top/RunIIFall17/c29f431
```

## Luminosity

After ntuples are processed, you can create the list of runs/lumi sections processed using crab as:
```
a=(`find grid/ -maxdepth 1 | grep crab_Data `)
for i in ${a[@]}; do
    crab kill ${i};
    crab status ${i};
    crab report ${i}; 
done
``` 
In case of failed jobs the missing lumis can be processed with the following script to wrap the tedious process of 
updating the cfg with a finer grain luminosity per job and the missing lumis json
```
for i in ${a[@]}; do
    python scripts/runMissingLumiSecs.py ${i}
done
```
You can then run the brilcalc tool to get the integrated luminosity in total and per run 
(see http://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html for more details).
The following script runs brilcalc inclusively and per trigger path, and stores the results in a ROOT file with the total integrated lumi per run.
It takes a bit to run, depending on the number of triggers configured to use in the analysis
```
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda/bin:$PATH
python scripts/convertLumiTable.py -o data/era2017/
```

## Preparing the analysis 

Correction and uncertainty files are stored under data by era directories (e.g data/era2017) in order no to mix different periods.

* Pileup weighting. To update the pileup distributions run the script below. It will store the data pileup distributions for different min.bias cross section in data/pileupWgts.root
```
python scripts/runPileupEstimation.py --json /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt --out data/era2017/pileupWgts.root
```
* B-tagging. To apply corrections to the simulation one needs the expected efficiencies stored somwewhere. The script below will project the jet pT spectrum from the TTbar sample before and after applying b-tagging, to compute the expecte efficiencies. The result will be stored in data/expTageff.root
```
python scripts/saveExpectedBtagEff.py -i /store/cmst3/group/top/RunIIFall17/c29f431/MC13TeV_TTJets -o data/era2017/expTageff.root;
```
* MC normalization. This will loop over all the samples available in EOS and produce a normalization cache (weights to normalize MC). The file will be available in data/genweights.pck
```
python scripts/produceNormalizationCache.py -i /store/cmst3/group/top/RunIIFall17/c29f431 -o data/era2017/genweights.root
```
The lepton/photon trigger/id/iso efficiencies should also be placed under data/era2017. 
The src/EfficiencyScaleFactorsWrapper.cc  should then be updated to handle the reading of the ROOT files and the application of the scale factors
event by event.

## Updating the code

Commit your changes regularly with
```
git commit -a -m'comment on the changes made'
```
Push to your forked repository
```
git push git@github.com:MYGITHUBLOGIN/TopLJets2015.git
```
From the github area of the repository cleak on the green button "Compare,review and create a pull request" to create the PR to merge with your colleagues.

# Local analyses

The ROOT trees created by MiniAnalyzer.cc can be analyzed with a simple executable.
See some examples under `src`.
The new executable should be included in `bin/analysisWrapper.cc` so that it can be used with the runLocalAnalysis.py script
which allows to run over single files or full directories.
See examples under test/ in ```steer*Analysis.sh```.
To plot the output of the local analysis you can run the following:
```
python scripts/plotter.py -i analysis/   -j data/era2017/samples.json  -l 12870
```

# Preparation of the data cards and workspaces

This part currently works only for the VBF analysis. It assumes that there are root files in the working directory which includes the plots created in the previous section. For the example below, the root file is called plotter_LowMJJ.root. Running the script below will make two data cards and a single corresponding workspace for signal region (here gamma+jets) and control region (here Z+jets) of the "LowMJJ" the category:
```
python scripts/makeWorkspace.py --Chan LowMJJ --nBin 5 --doSignalPH
```
The input histogram will be rebinned to have five bins. If you remove "--doSignalPH", the signal process in the signal and control region will NOT be connected via the transfer function.
