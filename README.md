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
 
