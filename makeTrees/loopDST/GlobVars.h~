#ifndef __GLOBVARS__
#define __GLOBVARS__

#include "mylibs.h"
#include "mystruc.h"

#define isRealData 1 //use 0 for simulation
#define isExtendedOutput 0 //use 0 for normal analysis, extended only for some use!!

TFile* out ;  // pointer to outputfile


//----------- BECAUSE WE NEED TO USE THESE IN loopDST_task.C and in *reco.h -------------
#if isEmbedding
HParticleEvtCharaBK evtChara;
#else
HParticleEvtChara evtChara;
#endif

HEnergyLossCorrPar dEdxCorr;

TF1 *fEnergyLossMdc;

#if isRealData
HParticleEventMixer eventmixer;
#else
HGenericEventMixer<HParticleCandSim> eventmixer;
#endif

//---------- NUMBER OF HISTS IN ARRAYS ------------
const int nCentralityClasses = 21;
const int nCutValues = 7;

//------------------------------ HISTOGRAMS -----------------------------

enum EventCounters_gen9
{
    cNumInputEv            =0,
    cisTriggerPT3          =1,
    cselectStart           =2,
    cselectStartPile       =3,
    cisGoodClusterVertex   =4,
    cisGoodVertexCand      =5,
    cisVeto                =6,
    cisStartVeto           =7,
    cisStartMeta           =8,
    cisGoodCentrality      =9,
    cNumEventCounters      =10

};

enum TrackCounters_gen9
{
    cAllTracks         =0,
    cisGoodSector      =1,
    cisUsed            =2,
    cisAtAnyEdge       =3,
    cisGoodSystem      =4,
    cChisqRK400        =5,
    cChisqRK400MM3     =6,
    cChisqRK100        =7,
    cChisqRK100MM2     =8,
    cNumTrackCounters  =9
};

TH1F *hEventCounter, *hTrackCounter;

TH2D *hTofElossMom, *hBetaMom_TOF, *hBetaMom_RPC, *hMdcEloss_TOF, *hMdcEloss_RPC, *hMassMom_TOF, *hMassMom_RPC, *hTofElossBeta_neg, *hZMdcEloss,
    *hZVerMult_noSelect, *hZVerMult_triggerSelect, *hZVerMult_allSelect, *hZVerMult_someSelect,
    *hTofElossBeta_OneRod_all[8], *hTofElossBeta_OneRod_massCut[8], *hTofElossMom_OneRod_all[8], *hTofElossMom_OneRod_massCut[8],
    *hMdcEloss_sim_TOF, *hMdcEloss_sim_RPC,*hMdcEloss_dataKcut_TOF, *hMdcEloss_dataKcut_RPC;
TH1D *hMass_TOF, *hMass_RPC, *hKaonMass_RPC[4], *hKaonMass_TOF[5];

TH2D *hBetaMom_TOF_massCut, *hBetaMom_RPC_massCut, *hMassMom_TOF_betamomCut, *hMassMom_RPC_betamomCut;
TH1D *hPionMass_TOF[4], *hPionMass_RPC[4];

TH1D *hChisqRK_proton, *hChisqMM_proton, *hChisqRK_pim, *hChisqMM_pim, *hChisqRK_pip, *hChisqMM_pip,
    *hChisqRK_km, *hChisqMM_km, *hChisqRK_kp, *hChisqMM_kp;

TH1D *hPiPairInvMass_all, *hPiPairInvMass_minTC,
    *hPiPairInvMass_mva50, *hPiPairInvMass_mva85, *hPiPairInvMass_mva95,
    *hPiPairInvMass_broadTC, *hPiPairInvMass_StoBratio, *hPiPairInvMass_Signif;


//control hist for centrality occurance in data set
TH1F *hCentralityCounter;

// 5% centrality class edges
#if isRealData || isEmbedding
Float_t CentralityEdges[]    = { 312., 180., 157., 136., 117.,  99., 82., 68., 55.};
#else
Float_t CentralityEdges[] = { 252., 183., 162., 142., 124., 107., 92., 78., 65.};
#endif

//event plane hists
TH1F *heventplane[nCentralityClasses];
TH1F *hphiAB[nCentralityClasses];

//occupancy efficiency hists
TH2 *effPiP_cc0 = NULL, *effPiP_cc1 = NULL, *effPiP_cc2 = NULL, *effPiP_cc3 = NULL,
    *effPiM_cc0 = NULL, *effPiM_cc1 = NULL, *effPiM_cc2 = NULL, *effPiM_cc3 = NULL,
    *effKp_cc0 = NULL, *effKp_cc1 = NULL, *effKp_cc2 = NULL, *effKp_cc3 = NULL,
    *effKm_cc0 = NULL, *effKm_cc1 = NULL, *effKm_cc2 = NULL, *effKm_cc3 = NULL,
    *effP_cc0 = NULL, *effP_cc1 = NULL, *effP_cc2 = NULL, *effP_cc3 = NULL;

//acceptance and reconstruction efficiency
TH2 *hAccCorr_cc0 = NULL, *hAccCorr_cc1 = NULL, *hAccCorr_cc2 = NULL, *hAccCorr_cc3 = NULL;
TH2 *hRecoEff_cc0 = NULL, *hRecoEff_cc1 = NULL, *hRecoEff_cc2 = NULL, *hRecoEff_cc3 = NULL;


//------- only for K0s ---------------
const Float_t preCuts[4][4] = { { 6., 12., 17., 13.}, { 7., 12., 20., 10.}, { 5., 12., 15., 20.}, {0.,100.,0.,100.} };
// 2 - my original precuts - not very well working
// 0 - my precuts 2 - woring well..
// 1 - Simon's precuts - ethalon!!!
// 3 -> effectively no cuts!! only for testing!!


//for gen8 -> might be needed to back-check
//const Int_t cutVals[7][4] = { {10,8,15,10},   //very broad cuts
//                                 {7,6,22,11},
//                                 {11,8,30,9},
//                                 {14,7,29,8}, //best by the s/sqrt(b) done in Feb18
//                                 {15,7,30,12},
//                                 {15,7,29,8},
//                                 {15,6,30,8} };

//for gen9 -> some research has been done..  (from ?? 2018 -> not that well controled bckg)
//const Int_t cutVals[7][4] = { {10,8,15,10},   //very broad cuts
//                                 {15,6,30,8}, //for gen8 comparison
//                                 {14,7,29,8}, //best by the s/sqrt(b) done in Feb18, tuned for gen8
//                                 {17,5,33,5},  //exponential bckg, best s/sqrt(b)
//				 {16,5,34,6},  //exponential bckg, almost as good s/sqrt(b) + more signal (s)
//                                 {20,4,16,3},  //linear bckg, best s/sqrt(b)
//                                 {18,4,31,3} }; //linear bckg, almost as good s/sqrt(b) + more signal (s)

// NEW Jan19 set of CV -> better control over the background (several types of bckg was used => best is Exp+Gauss for bckg and Gauss for signal)
//for gen9 -> all have similar significance (s/sqrt(s+b)) and I look for the most s as possible (reasonable well fitted..)
const Int_t cutVals[nCutValues][4] = { {5,12,21,13},
                                 {8,6,17,16},
                                 {5,9,24,10},
                                 {10,15,23,16},
				 {6,7,24,12},
                                 {8,9,21,15},
                                 {9,7,18,12} };

//------ only for K+ -----------------
const Float_t trackQualityCuts[6][2] = { {100., 2.}, {400., 3.}, {60., 1.5}, {100., 2.}, {100., 2.}, {100., 2.} };
const Bool_t trackElossCuts[6][2] = { {1,1}, {1,1}, {1,1}, {0,1}, {1,0}, {0,0} };
//chi2_RK, chi2_MM, dEdx_MDC, dEdx_TOF

//------- cut search ----------
HHistMap hm;


#endif
