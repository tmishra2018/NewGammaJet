# Instructions for producing Gamma + jet step 3

Step 1 is [here](https://github.com/lucastorterotot/DijetRootTreeMaker/blob/master/instructions/GammaJetTree_Instruction.md).

Step 2 is [here](https://github.com/lucastorterotot/DijetRootTreeAnalyzer).

**Target** CMS Gamma + Jet  - JEC calculation

**Authors** Hugues Lattaud and Lucas Torterotot

**Last update** 13 Feb 2019

## Introduction

This is a documentation for the gamma + jets framework. This Step is use to perform HLT selection, PUreweighting and plots for the tuple created at the [previous steps](https://github.com/lucastorterotot/DijetRootTreeAnalyzer).

## Installation recipe
```
export SCRAM_ARCH=slc6_amd64_gcc530
cmsrel CMSSW_8_0_25
cd CMSSW_8_0_25/src
git clone https://github.com/lucastorterotot/NewGammaJet.git JetMETCorrections/GammaJetFilter/
scram b -j 8
```

## Instructions (13 TeV)
### STEP 1: merging
First get the different files to merge. It's the files named `*reduced_skim.root` that came as the output of [DijetTreeAnalyser](https://github.com/lucastorterotot/DijetRootTreeAnalyzer). You can construct the list, for example, with
```
find /eos/user/l/ltortero/JEC-task/HT_Condor_output/DijetRootTreeAnalyzer/lists_2018/Run2018C-17Sep2018-v1/ -iname \*_reduced_skim\* > list_to_merge.txt
```
#### MC

In order to produce the plots for the MC outputs from the previous steps, one has to merge and add the weight of the different files.
The weight that to be applied at MC is defined as `evtWeightTot = xsec / sum_of_generatorWeights`.
This has to be done  in a separate step because it's necessary to run once over the full dataset in order to calculate the sum of generator weights.
In the output of [DijetTreeAnalyser](https://github.com/lucastorterotot/DijetRootTreeAnalyzer) one stored an histogram and a TTree filled using generator weights, in order to extract the sum of weights at the end with `Integral()`.
The cross section is given from the outside, you can get it using [DAS](https://cmsweb.cern.ch/das/) or MCM.
```
python mergeAndAddWeights.py -i [list_to_merge.txt] -o [output_directory] --xsec [number_from_DAS]
```
For example:
```
python mergeAndAddWeights.py -i GJet_Pt-15To6000_RunIIAutumn18MiniAOD-102X.txt -o /eos/user/l/ltortero/JEC-task/Step3_outputs/2018/ --xsec 283000.0
```
The merging will update the tree `analysis` and  the TTree `puvariable`  with a new branch called `evtWeightTot`.
This number is used in the following steps to fill histograms, to draw plots and perform PUreweighting. 

#### Data
For data the outputs are merged and the luminosity from `BrilCalc` is uploaded.
In order to calculate the integrated luminosity, the official recipe is followed.

To calculate the integrated luminosity, follow the `BrilCalc` [recipe](http://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html).

1) Get `lumiSummary.json` from crab
```
crab report -d crab_folder
```
2) Execute `brilcalc`
```
brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json -u /pb -i lumi_summary.json
```

In the end you can merge the output for the data with the command:
```
python mergeData.py -i [list_to_merge.txt] -o [output_directory] --lumi_tot [integrated_luminosity] --run [suffix for the run you process]
```

#### Step 1 conclusion
You should now have a root file for each MC dataset and one for each data dataset, with a prefix `PhotonJet_2ndLevel_`.
Copy those files somewhere else. A good place could be the folder in your working area, these are massive files.


### Step 2: PileUp

The MC is reweighting according to data based on the number of vertices in the event, this is made for each different HLT path in use for the analysis in order to take into account differences between simulation and data scenario from PU.
All the utilities to do that are available in the folder `analysis/PUReweighting`.
The relevant scripts are `ComputePU_perHLT.py` and `generate_mc_pileup.c`.

#### MC
1) Create a list on MC sample for which you want to calculate the PU reweighting.
This list contains all the MC files produced at the merging step.
For example you can create a list as `files_GJet_plus_QCD.list` which contains the files
```
[path]/PhotonJet_2ndLevel_GJet_Pythia_25ns_ReReco_2016-02-15.root
[path]/PhotonJet_2ndLevel_QCD_Pt-20toInf_2016-02-26.root  
```
2) Execute `generate_mc_pileup.c`. You have to compile with `Makefile`, and then
type the command followed by the list name (only central name).
```
./generate_mc_pileup.exe files_GJet_plus_QCD
```

#### Data
The pile up in data is calculated following the official recipe for different HLT path, written in `ComputePU_perHLT.py` that use `pileupCalc.py`.
At this script must be passed the json file for which you want to calculate the pu reweighting.
```
python ComputePU_perHLT.py --Cert_json [your processed Json file from crab report] --whichrun [suffix of the data you running on] 
```
This script will download the `pileuplatest.txt` file and then create one PU root file per HLT with a name begining by `pu_truth_data2018_100bins_HLTphoton`.  For example:
```
python ComputePU_perHLT.py --Cert_json /afs/cern.ch/work/l/ltortero/JEC-task/CMSSW_9_4_10/src/CMSDIJET/DijetRootTreeMaker/prod/crab/CERN_crab/crab_Run2017E-17Nov2017-v1_V3/results/processedLumis.json --whichrun E
```

### Step 3: Trigger selection
#### Data
To avoid any bias in the selection one explicitely require that for each bin in pt_gamma, only one trigger was active. For that, we use an XML description of the trigger of the analysis, as you can find in the`'bin/` folder. The description is file named `triggers.xml`.

The format should be straightforward: you have a separation in run ranges, as well as in triggers.
The weight of each HLT is used to reweight various distribution for the prescale.
The prescale is saved in the miniAOD and saved in the ntuples from step 1.

#### MC
You have a similar file for MC, named `triggers_mc.xml`. On this file, you have no run range, only a list of HLT path.
This list is used in order to know with HLT the event should have fired if it was data.

2012 note:
You can also specify multiple HLT path for one pt bin if there were multiple active triggers during the data taking period.
In this case, you'll need to provide a weight for each trigger (of course, the sum of the weight must be 1). Each trigger will be choose randolmy in order to respect the probabilities.

### Step 4: Plotting
### For a given alpha value
To produce the histograms with the response you can use the macro
`gammaJetFinalizer.cpp` in the `bin/` directory.

* MC:
```
gammaJetFinalizer -i output_singleFile_GJet.root -d Test_GJet --algo ak4 --alpha 0.30 --mc
```
* Data
```
gammaJetFinalizer -i output_singleFile_Data.root -d Test_Data --algo ak4 --alpha 0.30
```

These commands will produce two files with all necessary histograms. You should now have at least two files (three if you have run on QCD) and they will call
```
PhotonJet_Test_MC_PFlowAK4chs.root
PhotonJet_Test_Data_PFlowAK4chs.root
```

These final are used to produce the final plots. The code to do this is in `analysis/draw`.
The first time that you used this code you need to compile it. 
Go in that directory and compile with `make all`.

The drawers work only if you are in the same directory of the files produce in step 2 (named PhotonJet_*).
One may move the files in the directory `analysis/draw/Plot`.

The first drawer uses 1D histograms of the response and creates graphics with the response in function of pT for each bin of eta.
To run it:
```
../draw/drawPhotonJet_2bkg [Data_file] [GJet_file] [QCD_file] [jet type] [algorithm] [Normalization]
```
For example:
```
../draw/drawPhotonJet_2bkg Test_Data Test_MC Test_MC pf ak4 LUMI
```

At the moment the option `pf` and `ak4` are mandatory, we have only those jets.
The option `LUMI` is to normalize the plots at the data luminosity, while if you want to
normalize the plots at the same area you can use the option `SHAPE`.
The GJet files is passed twice because there is the opportunity to pass also a QCD file.

This drawer produces a directory called `PhotonJetPlots_<Data>_vs_<MC>_PFlowAK4_LUMI` that
contains the plot (png format) and it produces also a root file.

The second drawer is to perform the extrapolation at `alpha = 0` (no secondary activity)
To run it:
```
../draw/drawPhotonJetExtrap --type pf --algo ak4 [Data_file] [GJet_file] [QCD_file]
```
For example:
```
../draw/drawPhotonJetExtrap --type pf --algo ak4 Test_Data Test_MC Test_MC
```

The third drawer gets the previous results to build the plots for the extrapolated responses 
in function of pT in each bin of eta.
To run it:
```
../draw/draw_ratios_vs_pt Test_Data Test_MC Test_MC pf ak4
```

The plots are saved in the directory `PhotonJetPlots_<Data>_vs_<MC>_PFlowAK4_LUMI/vs_pt`.

The last drawer produces plots with some comparison between the different responses (MPF and Balancing) before and after the extrapolation.
To run it:
```
../draw/draw_all_methods_vs_pt Test_Data Test_MC Test_MC pf ak4
```

If everything went fine, you should now have a *lot* of plots in the folder `PhotonJetPlots_Data_file_vs_GJet_file_plus_QCD_file_PFlowAK4_LUMI`, and some more useful in the folder `PhotonJetPlots_Data_file_vs_GJet_file_plus_QCD_file_PFlowAK4_LUMI/vs_pt`.
In this last directory a root file named `plots.root` will be also saved.
This root file is very important because is used by Mikko for the global fit.

#### Cuts on alpha

You have to run all analysis (from Finalizer to this last drawer) for different alpha cut: 0.10, 0.15, 0.20, 0.25, 0.30.
The last drawer produces in the directory "PhotonJetPlots...../vs_pt/" a root file named plots.root.
So you will have a plots.root for each alpha cut, these for files have to be added (simple hadd) 
and send to Mikko in order to perform the global fit.

##### Macro
You may use `Produce_Combination_File_and_plots.py` which do the whole plotting step (step 4). For example:
```
python Produce_Combination_File_and_plots.py --output=TMP_TEST_3 --IOV=0 --Pu_profile=B --Input_data=/eos/user/l/ltortero/JEC-task/Step3_outputs/2017/PhotonJet_2ndLevel_DATA_RUN_B_2019-02-13.root --Input_mc=/eos/user/l/ltortero/JEC-task/Step3_outputs/2017/PhotonJet_2ndLevel_MC_1_2019-02-13.root
```

#### Any other business

Others drawers could be found in the `draw` directory.
For example `draw_vs_run` which draw the time dependence study ie response vs run number (only for Data).
```
../../draw/draw_vs_run Data_file pf ak4
```