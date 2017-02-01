# TopLJets2015

## Installation instructions

These installation instructions correspond to the 2016 data/MC Moriond17 re-reco.
To install execute the following in your work area.

```
cmsrel CMSSW_8_0_25
cd CMSSW_8_0_25/src 
cmsenv

#EGM electron MVA
git cms-merge-topic ikrav:egm_id_80X_v2
scram b -j 10
cd $CMSSW_BASE/external/$SCRAM_ARCH
git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout egm_id_80X_v1
cd $CMSSW_BASE/src

#EGM regression + smearer
git cms-merge-topic rafaellopesdesa:Regression80XEgammaAnalysis_v2
git cms-merge-topic shervin86:Moriond2017_JEC_energyScales
cd EgammaAnalysis/ElectronTools/data
git clone git@github.com:ECALELFS/ScalesSmearings.git
cd - 

#MET
git cms-merge-topic cms-met:METRecipe_8020
git cms-merge-topic ahinzmann:METRecipe_8020_Moriond17
git cms-merge-topic cms-met:fromCMSSW_8_0_20_postICHEPfilter

#BTV DeepCSV
git cms-merge-topic -u mverzett:DeepFlavourCMVA-from-CMSSW_8_0_21
mkdir RecoBTag/DeepFlavour/data/
cd RecoBTag/DeepFlavour/data/
wget http://home.fnal.gov/~verzetti//DeepFlavour/training/DeepFlavourNoSL.json
wget http://mon.iihe.ac.be/~smoortga/DeepFlavour/CMSSW_implementation_DeepCMVA/Model_DeepCMVA.json
cd -

#ntuplizer
git clone git@github.com:pfs/TopLJets2015.git
cd TopLJets2015/TopAnalysis
git checkout 80x_rereco_rev

#delete, not used
cd -
rm -rf DQMOffline

#compile
scram b -j 8
```

## Running ntuple creation and checking the selection

The ntuplizer is steered with test/runMiniAnalyzer_cfg.py.
It takes several options from command line (see cfg for details).
To run locally the ntuplizer, for testing purposes do something like:

```
cmsRun test/runMiniAnalyzer_cfg.py runOnData=False era=era2016 outFilename=MC13TeV_TTJets.root
cmsRun test/runMiniAnalyzer_cfg.py runOnData=True  era=era2016 outFilename=Data13TeV_MuonEG.root
```

The default files point to the ones used in the TOP synchronization exercise
(see details in )
The output files can be analyzed with a simple executable (also a skeleton for the implementation of new analysis) as detailed below.
The output file will contain cutflow histograms which can be used to check the selection against the synchronization twiki.
Check the implementation of this analysis in src/TOPSynchExercise.cc. 
Other analysis should also be implemented there.

```
analysisWrapper --in MC13TeV_TTJets.root   --out mc_synch.root   --method TOPSynchExercise::RunTOPSynchExercise
analysisWrapper --in Data13TeV_MuonEG.root --out data_synch.root --method TOPSynchExercise::RunTOPSynchExercise
```

## Submitting ntuple creation through the grid

To submit the ntuplizer to the grid start by setting the environment for crab3.
More details can be found in [CRAB3CheatSheet](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3CheatSheet#Environment_setup)

```
source /cvmfs/cms.cern.ch/crab3/crab.sh
```
The following script helps submitting a list of files described in a json file.
Partial submission can be made adding "-o csv_list" as an option.
Adding "-s" will trigger the submission to the grid (otherwise the script only writes down the crab cfg files)

```
python scripts/submitToGrid.py -j data/era2016/samples.json -c ${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/test/runMiniAnalyzer_cfg.py 
```

As soon as ntuple production starts to finish, to move from crab output directories to a simpler directory structure which can be easily parsed by the local analysis runThe merging can be run locally if needed by using the checkProductionIntegrity.py script

```
python scripts/submitCheckProductionIntegrity.py -i /store/group/phys_top/psilva/8274336 -o /store/cmst3/group/top/ReReco2016/8274336
```

After ntuples are processed, you can create the list of runs/lumi sections processed using crab as:
```
a=(`find grid/ -maxdepth 1 | grep crab_Data `)
for i in ${a[@]}; do
    crab report ${i}; 
done
``` 
Then you can merge the json files for the same dataset to get the full list of run/lumi sections to analyse
```
mergeJSON.py grid/crab_Data13TeV_DoubleMuon_2016B/results/processedLumis.json grid/crab_Data13TeV_DoubleMuon_2015C/results/processedLumis.json grid/crab_Data13TeV_DoubleMuon_2015D/results/processedLumis.json --output data/era2016/Data13TeV_DoubleMuon_lumis.json
```


## Preparing the analysis 

Correction and uncertainty files are stored under data by era directories (e.g. data/era2015, data/era2016) in order no to mix different periods.

You can then run the brilcalc tool to get the integrated luminosity in total and per run (see https://twiki.cern.ch/twiki/bin/view/CMS/2015LumiNormtag for more details).
```
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i data/era2016/Data13TeV_DoubleMuon_lumis.json
```
Use the table which is printed out to update the "lumiPerRun" method in ReadTree.cc.
That will be used to monitor the event yields per run in order to identify outlier runs.
* Pileup weighting. To update the pileup distributions run the script below. It will store the data pileup distributions for different min.bias cross section in data/pileupWgts.root
```
python scripts/runPileupEstimation.py --json /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --out data/era2016/pileupWgts.root
```
* B-tagging. To apply corrections to the simulation one needs the expected efficiencies stored somwewhere. The script below will project the jet pT spectrum from the TTbar sample before and after applying b-tagging, to compute the expecte efficiencies. The result will be stored in data/expTageff.root
```
python scripts/saveExpectedBtagEff.py -i /store/cmst3/group/top/ReReco2016/8274336/MC13TeV_TTJets -o data/era2016/expTageff.root;
```
* MC normalization. This will loop over all the samples available in EOS and produce a normalization cache (weights to normalize MC). The file will be available in data/genweights.pck
```
python scripts/produceNormalizationCache.py -i /store/cmst3/group/top/ReReco2016/8274336 -o data/era2016/genweights.root
```
The lepton trigger/id/iso efficiencies should also be placed under data/era2016. 
The src/LeptonEfficiencyWrapper.cc  should then be updated to handle the reading of the ROOT files and the application of the scale factors
event by event.


## Running locally the analysis for testing

The analysis (histogram filling, final selection) is in src/ReadTree.cc.
Recompile (scram b) everytime you change it so that you can test the new features.
To test the code on a single file to produce plots.
```
python scripts/runLocalAnalysis.py -i MiniEvents.root
```
To run the code on a set of samples stored in EOS you can run it as shown below.
If "-q queue_name" is appended the jobs are submitted to the batch system instead of running locally. 
To check the status of your jobs run "bjobs" and then "bpeek job_number" if you want to inspect how the job is running in the cluster.
If "-n n_jobs" is passed the script runs locally using "n_jobs" parallel threads.
```
python scripts/runLocalAnalysis.py -i /store/cmst3/user/psilva/LJets2015/7e62835 -n 8 --runSysts -o analysis_muplus   --ch 13   --charge 1
```
If you want to suppress the mails sent automatically after job completion please do
```
export LSB_JOB_REPORT_MAIL=N
```
before submitting the jobs to the batch. After the jobs have run you can merge the outputs with
```
./scripts/mergeOutputs.py analysis_muplus
```
To plot the output of the local analysis you can run the following:
```
python scripts/plotter.py -i analysis_muplus/   -j data/era2016/samples.json  -l 12870
```
After the plotters are created one can run the QCD estimation normalization, by fitting the MET distribution.
The script will also produce the QCD templates using the data from the sideband region. It runs as
```
python scripts/runQCDEstimation.py --iso analysis_muplus/plots/plotter.root --noniso analysis_munoniso/plots/plotter.root --out analysis_muplus/
```
The output is a ROOT file called Data_QCDMultijets.root which can now be used in addition to the predictions of all the other backgrounds.
To include it in the final plots you can run the plotter script again (see instructions above).

## Submitting the full analysis to the batch system

A script wraps up the above procedure for all the signal and control regions used in the analyis.
To use it you can use the following script
```
sh scripts/steerAnalysis.sh <DISTS/MERGE/PLOT/BKG>
```

## Cross section fitting

We use the Higgs combination tool to perform the fit of the production cross section.
(cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit for details of the release to use).
It currently has to be run from a CMSSW_7_1_5 release. To create the datacard you can run the following script
```
python scripts/createDataCard.py -i analysis_muplus/plots/plotter.root -o  analysis_muplus/datacard  -q analysis_muplus/.qcdscalefactors.pck -d nbtags
```
The script can be used to create the datacard from any histogram stored in plotter.root.
For the systematic variations it expects a 2D histogram named as HISTONAMEshapes_{exp,gen} filled with alternative variations of the shape,
being exp/gen used for experimental/generator-level systematics.
Additional systematics from alternative samples can also be used to build the datacards using the --systInput option.
Other options are available to choose the categories to use.
The datacards can be further combined using the standard combineCards.py script provided by the Higgs Combination package.

To run the fits and show the results you can use the following script.
```
python scripts/fitCrossSection.py "#mu^{+}"=analysis_muplus/datacard/datacard.dat -o analysis_muplus/datacard &
```
If --noFit is passed it displays the results of the last fit. The script is a wrapper used to run combine 
to perform the fit with and without systematics, produce the post-fit nuisance parameters summary
and the likelihood scans.
For the standard analysis one can re-use the steerAnalysis.sh script with two options CinC/SHAPE
will run the Cut-in-Categories/Shape analyses.
```
sh scripts/steerAnalysis.sh CinC/SHAPE
```

## Updating the code

Commit your changes regularly with
```
git commit -a -m'comment on the changes made'
```
Push to your forked repository
```
git push git@github.com:MYGITHUBLOGIN/TopLJets2015.git
```
From the github area of the repository cleak on the green button "Compare,review and create a pull request"
to create the PR to merge with your colleagues.
