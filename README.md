# kaonAna
Repository for software connected to kaon flow analysis in HADES data. Tailored for Ag+Ag @ 1.58A GeV.

### Analysis structure

The analysis is divided into four main steps:
 - making TTrees with useful data from DST files which provides significant filesize reduction and an easy-to-use data structure;
 - filling necessary histograms from the TTrees, mainly mass distributions;
 - fitting mass distribution peaks, background subtraction, preparing flow (Dphi) histograms and fitting the Fourier series;
 - extracting results from fits and preparing plots
 
 Directories corresponding to these stages are (respectively): `makeTrees`, `fillHistFromTree`, `flowAna` and `results`.
 
 ### makeTrees
 
##### Launching

To launch this stage of the analysis, run the bash script `sendJobArrayScript.sh`. It will create or load jobarray files and then run the `jobArrayScript.sh` which controls the actual program - `analysis`. The command should look like this:

```. sendJobArrayScript.sh DAY DATASET PARTICLE PID TRACKCUTS OCPID```

DAY is the selected day of data (for Ag+Ag@1.58A data this should be 058 or any value selected from 062 - 088). DATASET is the particular type of files that you wish to use. Currently only `accepted` is supported for all days and `test` (small dataset to debug loopDST.C locally) for day 074. These two parameteres together select DST files which will be used for analysis. 

PARTICLE is the prefix for file name and so it can be anything. For consistency `Kcharged` is used for trees with K+ and K- while `Kzero` is used for TTrees for K0s analysis. 

PID, TRACKCUTS and OCPID are important parameters which are passed on to the `analysis` program. They select different pid/track cuts/occupancy correction scenarios.

From which machine do you do this???

##### Compilation



#####

