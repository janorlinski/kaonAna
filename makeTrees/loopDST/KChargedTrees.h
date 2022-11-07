#ifndef __KCHARGEDTREES___     // good style to avaoid multiple includes
#define __KCHARGEDTREES___

#include "mylibs.h"
#include "mystruc.h"

#include "GlobVars.h"

#include "UsefulFunctions.h"

using namespace std;

class KChargedTrees : public HReconstructor
{
protected:
    // put all vars here which are not
    // local to a function

    Int_t PID_option=-1, iTrackCuts=-1, OccupCorr_pid=-1;

    HParticleTrackSorter sorter;

    CAND_KCHARGED pluskaon;
    TTree *pluskaontree;

    CAND_KCHARGED minuskaon;
    TTree *minuskaontree;

#if !isRealData
    GKINE_KCHARGED gkplus;
    TTree *gkplustree;

    GKINE_KCHARGED gkminus;
    TTree *gkminustree;
#endif
    
    //------------------------------------
    // file handling
    TFile *effcorFile;
    TFile *mdcdEdx_cutfile_new_rpc;
    TFile *mdcdEdx_cutfile_new_tof;
    TFile *mdcdEdx_cutfile_new_tof2;
    TFile *tofdEdx_cutfile;
    TFile *fileAccCorrRecEff;


    TCutG *mdcdEdx_tof_Kp_data_new, *mdcdEdx_rpc_Kp_data_new, *mdcdEdx_tof_Km_data_new, *mdcdEdx_rpc_Km_data_new, *tofdEdx_beta_noK;

    Float_t getEff(Int_t pid, Int_t CentralityClass, Float_t dphi, Float_t the)
    {
	Float_t eff = -1.;
	switch(pid){

	case 8:
	    if (effPiP_cc0 && CentralityClass==0) eff = effPiP_cc0->GetBinContent(effPiP_cc0->FindBin(dphi,the));
	    if (effPiP_cc1 && CentralityClass==1) eff = effPiP_cc1->GetBinContent(effPiP_cc1->FindBin(dphi,the));
	    if (effPiP_cc2 && CentralityClass==2) eff = effPiP_cc2->GetBinContent(effPiP_cc2->FindBin(dphi,the));
	    if (effPiP_cc3 && CentralityClass==3) eff = effPiP_cc3->GetBinContent(effPiP_cc3->FindBin(dphi,the));
	    break;

	case 9:
	    if (effPiM_cc0 && CentralityClass==0) eff = effPiM_cc0->GetBinContent(effPiM_cc0->FindBin(dphi,the));
	    if (effPiM_cc1 && CentralityClass==1) eff = effPiM_cc1->GetBinContent(effPiM_cc1->FindBin(dphi,the));
	    if (effPiM_cc2 && CentralityClass==2) eff = effPiM_cc2->GetBinContent(effPiM_cc2->FindBin(dphi,the));
	    if (effPiM_cc3 && CentralityClass==3) eff = effPiM_cc3->GetBinContent(effPiM_cc3->FindBin(dphi,the));
	    break;

	case 11:
	    if (effKp_cc0 && CentralityClass==0) eff = effKp_cc0->GetBinContent(effKp_cc0->FindBin(dphi,the));
	    if (effKp_cc1 && CentralityClass==1) eff = effKp_cc1->GetBinContent(effKp_cc1->FindBin(dphi,the));
	    if (effKp_cc2 && CentralityClass==2) eff = effKp_cc2->GetBinContent(effKp_cc2->FindBin(dphi,the));
	    if (effKp_cc3 && CentralityClass==3) eff = effKp_cc3->GetBinContent(effKp_cc3->FindBin(dphi,the));
	    break;
        //TODO -> Km occupancy corrections
//	case 12:
//	    if (effKm_cc0 && CentralityClass==0) eff = effKm_cc0->GetBinContent(effKm_cc0->FindBin(dphi,the));
//	    if (effKm_cc1 && CentralityClass==1) eff = effKm_cc1->GetBinContent(effKm_cc1->FindBin(dphi,the));
//	    if (effKm_cc2 && CentralityClass==2) eff = effKm_cc2->GetBinContent(effKm_cc2->FindBin(dphi,the));
//	    if (effKm_cc3 && CentralityClass==3) eff = effKm_cc3->GetBinContent(effKm_cc3->FindBin(dphi,the));
//	    break;

	case 14:
	    if (effP_cc0 && CentralityClass==0) eff = effP_cc0->GetBinContent(effP_cc0->FindBin(dphi,the));
	    if (effP_cc1 && CentralityClass==1) eff = effP_cc1->GetBinContent(effP_cc1->FindBin(dphi,the));
	    if (effP_cc2 && CentralityClass==2) eff = effP_cc2->GetBinContent(effP_cc2->FindBin(dphi,the));
	    if (effP_cc3 && CentralityClass==3) eff = effP_cc3->GetBinContent(effP_cc3->FindBin(dphi,the));
	    break;


	default:

	    cerr<< " PID "<<pid<<" not known -> functions only for pi+ (8) and pi- (9) and proton (14)"<<endl;
	    break;
	}

	return eff;
    }

    Float_t getAccCorr(Int_t CentralityClass, Float_t pt, Float_t y0)
    {
	Float_t eff = -1;
	switch(CentralityClass){

	case 0:
	    if(hAccCorr_cc0) eff = hAccCorr_cc0->GetBinContent(hAccCorr_cc0->FindBin(y0,pt));
            break;

	case 1:
	    if(hAccCorr_cc1) eff = hAccCorr_cc1->GetBinContent(hAccCorr_cc1->FindBin(y0,pt));
            break;

	case 2:
	    if(hAccCorr_cc2) eff = hAccCorr_cc2->GetBinContent(hAccCorr_cc2->FindBin(y0,pt));
            break;

	case 3:
	    if(hAccCorr_cc3) eff = hAccCorr_cc3->GetBinContent(hAccCorr_cc3->FindBin(y0,pt));
            break;


	default:
	    cerr<< "Centrality Class "<< CentralityClass <<" not known -> functions only for 0 (0-10%) up to 3 (30-40%)."<<endl;
	    break;

	}

        return eff;
    }

    Float_t getRecoEff(Int_t CentralityClass, Float_t pt, Float_t y0)
    {
	Float_t eff = -1;
	switch(CentralityClass){

	case 0:
	    if(hRecoEff_cc0) eff = hRecoEff_cc0->GetBinContent(hRecoEff_cc0->FindBin(y0,pt));
            break;

	case 1:
	    if(hRecoEff_cc1) eff = hRecoEff_cc1->GetBinContent(hRecoEff_cc1->FindBin(y0,pt));
            break;

	case 2:
	    if(hRecoEff_cc2) eff = hRecoEff_cc2->GetBinContent(hRecoEff_cc2->FindBin(y0,pt));
            break;

	case 3:
	    if(hRecoEff_cc3) eff = hRecoEff_cc3->GetBinContent(hRecoEff_cc3->FindBin(y0,pt));
            break;


	default:
	    cerr<< "Centrality Class "<< CentralityClass <<" not known -> functions only for 0 (0-10%) up to 3 (30-40%)."<<endl;
	    break;

	}

        return eff;
    }



public:

    KChargedTrees(const Text_t *name = "",const Text_t *title ="", Int_t PID_OPTION=-1, Int_t TRACK_CUTS=-1, Int_t OCCUP_CORR_PID=-1)
	: HReconstructor(name,title)
    {  // init your vars
	PID_option = PID_OPTION;
	iTrackCuts = TRACK_CUTS;
	OccupCorr_pid = OCCUP_CORR_PID;

    }

    virtual ~KChargedTrees()
    {
	// clean up all dynamically created objects

    }

    Bool_t init()
    {
	// this function is called once in the beginning
	// create histograms or get pointer to param containers
	// or data containers here. Tip: Create first your
	// output file and and after that your histograms/ ntuples.
	// In this way you make shure the object will end up in your
        // root file and not any other one.

	if(out) {
	    out->cd();
	    out->mkdir("kaonPlus");
	    out->cd("kaonPlus");

	    pluskaontree = new TTree("pluskaontree","K+ Tree");
	    pluskaontree->Branch("pluskaon",&pluskaon,"eventPlane:phiAB:"
			       "eventVerX:eventVerY:eventVerZ:evntVerChi2:"
			       "system:theta:phi:momX:momY:momZ:mom:"
			       "multiplicity:centralityclass:MDCtracks:METAclst:"
				 "invMass:mass:mt:mtm0:pt:y:ycm:y0:"
			       "flow:occupancyCorr:accCorr:recoEff:eventWeight:phiMass:downflag");

	    out->cd();
	    out->mkdir("kaonMinus");
	    out->cd("kaonMinus");

	    minuskaontree = new TTree("minuskaontree","K- Tree");
	    minuskaontree->Branch("minuskaon",&minuskaon,"eventPlane:phiAB:"
			       "eventVerX:eventVerY:eventVerZ:evntVerChi2:"
			       "system:theta:phi:momX:momY:momZ:mom:"
			       "multiplicity:centralityclass:MDCtracks:METAclst:"
			       "invMass:mass:mt:mtm0:pt:y:ycm:y0:"
			       "flow:occupancyCorr:accCorr:recoEff:eventWeight:phiMass:downflag");

#if !isRealData
	    out->cd();
	    out->mkdir("simulation");
	    out->cd("simulation");

	    gkplustree = new TTree("gkplustree","GeantKine Kaon Tree");
	    gkplustree->Branch("gkplus",&gkplus,"eventPlane:eventPlane_gkine:phiAB:"
			      "system:theta:phi:momX:momY:momZ:mom:"
			      "multiplicity:centralityclass:impactparam:MDCtracks:METAclst:"
			       "invMass:mt:mtm0:pt:y:ycm:y0:"
#if isEmbedding
			       "pcand_chi2:pcand_innerchi2:pcand_mm:pcand_mass:pcand_system:pcand_tofeloss:pcand_mdceloss:pcand_beta:pcand_mom:pcand_charge:"
#endif
			       "flow:flow_gkine:weight:"
			       "kAcc/B:kRec:"
#if isEmbedding
			       "pcand_mdcedge:pcand_goodflags:"
#endif
			       "hasDecayed");

	    gkminustree = new TTree("gkminustree","GeantKine Kaon Tree");
	    gkminustree->Branch("gkminus",&gkminus,"eventPlane:eventPlane_gkine:phiAB:"
				"system:theta:phi:momX:momY:momZ:mom:"
				"multiplicity:centralityclass:impactparam:MDCtracks:METAclst:"
				"invMass:mt:mtm0:pt:y:ycm:y0:"
#if isEmbedding
				"pcand_chi2:pcand_innerchi2:pcand_mm:pcand_mass:pcand_system:pcand_tofeloss:pcand_mdceloss:pcand_beta:pcand_mom:pcand_charge:"
#endif
				"flow:flow_gkine:weight:"
				"kAcc/B:kRec:"
#if isEmbedding
				"pcand_mdcedge:pcand_goodflags:"
#endif
				"hasDecayed");

#endif

	    const char *effcorFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/ALLMY_OC_matrix_apr12_day108_gen9_4mbins.root";
//	    const char *effcorFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/MY_KpIII_effcor_flow_matrix_apr12_day108_gen8_4mbins.root";
	    effcorFile = TFile::Open( effcorFileName );
	    effPiM_cc3 = (TH2F *)effcorFile->Get("heff4_id9");
	    effPiP_cc3 = (TH2F *)effcorFile->Get("heff4_id8");
	    effKp_cc3  = (TH2F *)effcorFile->Get("heff4_id11");
//	    effKm_cc3  = (TH2F *)effcorFile->Get("heff4_id12");
	    effP_cc3   = (TH2F *)effcorFile->Get("heff4_id14");
	    effPiM_cc2 = (TH2F *)effcorFile->Get("heff3_id9");
	    effPiP_cc2 = (TH2F *)effcorFile->Get("heff3_id8");
	    effKp_cc2  = (TH2F *)effcorFile->Get("heff3_id11");
//	    effKm_cc2  = (TH2F *)effcorFile->Get("heff3_id12");
	    effP_cc2   = (TH2F *)effcorFile->Get("heff3_id14");
	    effPiM_cc1 = (TH2F *)effcorFile->Get("heff2_id9");
	    effPiP_cc1 = (TH2F *)effcorFile->Get("heff2_id8");
	    effKp_cc1  = (TH2F *)effcorFile->Get("heff2_id11");
//	    effKm_cc1  = (TH2F *)effcorFile->Get("heff2_id12");
	    effP_cc1   = (TH2F *)effcorFile->Get("heff2_id14");
	    effPiM_cc0 = (TH2F *)effcorFile->Get("heff1_id9");
	    effPiP_cc0 = (TH2F *)effcorFile->Get("heff1_id8");
	    effKp_cc0  = (TH2F *)effcorFile->Get("heff1_id11");
//	    effKm_cc0  = (TH2F *)effcorFile->Get("heff1_id12");
	    effP_cc0   = (TH2F *)effcorFile->Get("heff1_id14");

	    const char *filenameAccCorrRecEff = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/AccCorr_and_RecoEff_fineBinning.root";
	    fileAccCorrRecEff = TFile::Open( filenameAccCorrRecEff );
            hAccCorr_cc0 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusAccCorr_centrality_0_10");
            hAccCorr_cc1 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusAccCorr_centrality_10_20");
            hAccCorr_cc2 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusAccCorr_centrality_20_30");
            hAccCorr_cc3 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusAccCorr_centrality_30_40");
            hRecoEff_cc0 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusRecoEff_centrality_0_10");
            hRecoEff_cc1 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusRecoEff_centrality_10_20");
            hRecoEff_cc2 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusRecoEff_centrality_20_30");
            hRecoEff_cc3 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusRecoEff_centrality_30_40");

	    const char *tofcutFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/Kp_cutfile_tof_gen8a.root";
	    mdcdEdx_cutfile_new_tof = TFile::Open( tofcutFileName );
	    mdcdEdx_tof_Kp_data_new =(TCutG*)mdcdEdx_cutfile_new_tof->Get("KP_mdcdEdx_cut_scaled");

	    const char *tofcutFileName2 = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/Kp_cutfile_tof_gen8a_15sigma.root";
	    mdcdEdx_cutfile_new_tof2 = TFile::Open( tofcutFileName2 );
	    mdcdEdx_tof_Km_data_new =(TCutG*)mdcdEdx_cutfile_new_tof2->Get("KP_mdcdEdx_cut_scaled");

            const char *rpccutFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/Kp_cutfile_rpc_gen8a_rpc25.root";
	    mdcdEdx_cutfile_new_rpc = TFile::Open( rpccutFileName );
	    mdcdEdx_rpc_Kp_data_new =(TCutG*)mdcdEdx_cutfile_new_rpc->Get("KP_mdcdEdx_cut_scaled");
	    mdcdEdx_rpc_Km_data_new =(TCutG*)mdcdEdx_cutfile_new_rpc->Get("KP_mdcdEdx_cut_scaled");

	    const char *kaonVetocutFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/betaTofdEdx_noKm.root";
	    tofdEdx_cutfile = TFile::Open( kaonVetocutFileName );
	    tofdEdx_beta_noK =(TCutG*)tofdEdx_cutfile->Get("betaTofdEdx_noKm");


	    out->cd();
	    out->mkdir("eventPlanes");
	    out->cd("eventPlanes");

	    for(int icc=0; icc<nCentralityClasses; icc++){
		heventplane[icc] = new TH1F(Form("heventplane_icc%i",icc),"event plane distribution",200,-200,200);
		hphiAB[icc] = new TH1F(Form("hphiAB_icc%i",icc),"difference of event plane for two subsets",180,0,180);
	    }


	    out->cd();
	    out->mkdir("histos");
	    out->cd("histos");

	    hCentralityCounter = new TH1F("hCentralityCounter","number of events in each centrality class",10,0,10);

	    hEventCounter = new TH1F( "hEventCounter", "Counter of events", cNumEventCounters, 0, cNumEventCounters);
//		hEventCounter->GetXaxis()->SetBinLabel( 1, "InputEv" );
//		hEventCounter->GetXaxis()->SetBinLabel( 2, "PT3" );
//		hEventCounter->GetXaxis()->SetBinLabel( 3, "selectStart" );
//		hEventCounter->GetXaxis()->SetBinLabel( 4, "StartPileUp" );
//		hEventCounter->GetXaxis()->SetBinLabel( 5, "GoodClusterVertex" );
//		hEventCounter->GetXaxis()->SetBinLabel( 6, "GoodCandVertex" );
//		hEventCounter->GetXaxis()->SetBinLabel( 7, "NoVeto" );
//		hEventCounter->GetXaxis()->SetBinLabel( 8, "VetoStart" );
//		hEventCounter->GetXaxis()->SetBinLabel( 9, "StartMeta" );
//		hEventCounter->GetXaxis()->SetBinLabel( 10, "GoodCentrality" );

	    hTrackCounter = new TH1F( "hTrackCounter", "Counter of tracks", cNumTrackCounters, 0, cNumTrackCounters);
//		hTrackCounter->GetXaxis()->SetBinLabel( 1, "AllTracks" );
//		hTrackCounter->GetXaxis()->SetBinLabel( 2, "GoodSector" );
//		hTrackCounter->GetXaxis()->SetBinLabel( 3, "GoodTrackSorter" );
//		hTrackCounter->GetXaxis()->SetBinLabel( 4, "notAtAnyEdge" );
//		hTrackCounter->GetXaxis()->SetBinLabel( 5, "GoodSystem" );
//		hTrackCounter->GetXaxis()->SetBinLabel( 6, "Chi2RK400" );
//		hTrackCounter->GetXaxis()->SetBinLabel( 7, "Chi2RK400MM3" );
//		hTrackCounter->GetXaxis()->SetBinLabel( 8, "Chi2RK100" );
//		hTrackCounter->GetXaxis()->SetBinLabel( 9, "Chi2RK100MM2" );


	    hTofElossMom  = new TH2D("hTofElossMom","TOF dE/dx vs mom",150,-1000,2000,100,0,20);
	    hBetaMom_TOF  = new TH2D("hBetaMom_TOF","beta vs mom -> TOF",300,-3000,3000,120,0,1.2);
	    hBetaMom_RPC  = new TH2D("hBetaMom_RPC","beta vs mom -> RPC",300,-3000,3000,120,0,1.2);
	    hMdcEloss_TOF = new TH2D("hMdcEloss_TOF","MDC dE/dx vs mom ->TOF",150,-1000,2000,100,0,20);
	    hMdcEloss_RPC = new TH2D("hMdcEloss_RPC","MDC dE/dx vs mom ->RPC",150,-1000,2000,100,0,20);
	    hMassMom_TOF  = new TH2D("hMassMom_TOF","mass vs mom -> TOF",500,-1000,4000,150,0,3000);
	    hMassMom_RPC  = new TH2D("hMassMom_RPC","mass vs mom -> RPC",500,-1000,4000,150,0,3000);
	    hMdcEloss_sim_TOF = new TH2D("hMdcEloss_sim_TOF","MDC dE/dx vs mom ->TOF",100,-1000,1000,100,0,20);
	    hMdcEloss_sim_RPC = new TH2D("hMdcEloss_sim_RPC","MDC dE/dx vs mom ->RPC",100,-1000,1000,100,0,20);
	    hMdcEloss_dataKcut_TOF = new TH2D("hMdcEloss_dataKcut_TOF","MDC dE/dx vs mom ->TOF",100,-1000,1000,100,0,20);
	    hMdcEloss_dataKcut_RPC = new TH2D("hMdcEloss_dataKcut_RPC","MDC dE/dx vs mom ->RPC",100,-1000,1000,100,0,20);

	    hTofElossBeta_neg = new TH2D("hTofElossBeta_neg"," TOF dE/dx vs beta -> negative charge",240,0,1.2,250,0,25);
	    for(int ih=0; ih<8; ih++){
		hTofElossBeta_OneRod_all[ih] = new TH2D(Form("hTofElossBeta_OneRod%i_all",ih)," TOF dE/dx vs beta -> all tracks but only one rod",240,0,1.2,250,0,25);
		hTofElossBeta_OneRod_massCut[ih] = new TH2D(Form("hTofElossBeta_OneRod%i_massCut",ih)," TOF dE/dx vs beta -> tracks with mass>600 but only one rod",240,0,1.2,250,0,25);
		hTofElossMom_OneRod_all[ih] = new TH2D(Form("hTofElossMom_OneRod%i_all",ih)," TOF dE/dx vs mom*q -> all tracks but only one rod",400,-1000,3000,250,0,25);
		hTofElossMom_OneRod_massCut[ih] = new TH2D(Form("hTofElossMom_OneRod%i_massCut",ih)," TOF dE/dx vs mom*q -> tracks with mass>600 but only one rod",400,-1000,3000,250,0,25);
	    }

	    hZMdcEloss = new TH2D("hZMdcEloss","Log_e( MDC dE/dx exp/proton) vs mass",400,0.,4000,300,-2,4);

            hZVerMult_noSelect = new TH2D("hZVerMult_noSelect","Vertex_z vs N_Tof+Rpc -> no event selection",240,-110.,10.,200.,5.,205.);
            hZVerMult_triggerSelect = new TH2D("hZVerMult_triggerSelect","Vertex_z vs N_Tof+Rpc -> no event selection",240,-110.,10.,200.,5.,205.);
            hZVerMult_allSelect = new TH2D("hZVerMult_allSelect","Vertex_z vs N_Tof+Rpc -> all event selection",240,-110.,10.,200.,5.,205.);
            hZVerMult_someSelect = new TH2D("hZVerMult_someSelect","Vertex_z vs N_Tof+Rpc -> some event selection",240,-110.,10.,200.,5.,205.);

            hMass_TOF = new TH1D("hMass_TOF","mass -> TOF",500,-1000,4000);
            hMass_RPC = new TH1D("hMass_RPC","mass -> RPC",500,-1000,4000);

	    for(int ih=0; ih<4; ih++)
                hKaonMass_RPC[ih] = new TH1D( Form("hKaonMass_RPC_%i",ih), "mass zoom -> RPC with conclusive cuts", 220,-1000,1200);

	    for(int ih=0; ih<5; ih++)
                hKaonMass_TOF[ih] = new TH1D( Form("hKaonMass_TOF_%i",ih), "mass zoom -> TOF with conclusive cuts", 220,-1000,1200);

            hChisqRK_proton = new TH1D("hChisqRK_proton", "chi2 of RK method for one particle", 80,0.,400.);
            hChisqRK_pim    = new TH1D("hChisqRK_pim", "chi2 of RK method for one particle", 80,0.,400.);
            hChisqRK_pip    = new TH1D("hChisqRK_pip", "chi2 of RK method for one particle", 80,0.,400.);
            hChisqRK_km     = new TH1D("hChisqRK_km", "chi2 of RK method for one particle", 80,0.,400.);
            hChisqRK_kp     = new TH1D("hChisqRK_kp", "chi2 of RK method for one particle", 80,0.,400.);

	    hChisqMM_proton = new TH1D("hChisqMM_proton", "chi2 of MM method for one particle", 80,0.,3.);
            hChisqMM_pim    = new TH1D("hChisqMM_pim", "chi2 of MM method for one particle", 80,0.,3.);
            hChisqMM_pip    = new TH1D("hChisqMM_pip", "chi2 of MM method for one particle", 80,0.,3.);
            hChisqMM_km     = new TH1D("hChisqMM_km", "chi2 of MM method for one particle", 80,0.,3.);
            hChisqMM_kp     = new TH1D("hChisqMM_kp", "chi2 of MM method for one particle", 80,0.,3.);

	}

	fEnergyLossMdc = new TF1("fEnergyLossMdc",EnergyLossMdc,0.,3000.,2);
        fEnergyLossMdc->SetParameter(1,0.6); //He fraction

	return kTRUE;
    }
    Bool_t reinit()
    {   // this function is called for each file
        // after init()
	// use this place if you need already initialized
	// containers


        //WHERE TO PUT THIS??
	//--------------------------CONFIGURATION---------------------------------------------------
	//At begin of the program (outside the event loop)
	//sorter.setDebug();                                            // for debug
	//sorter.setPrintLevel(3);                                      // max prints
	//sorter.setRICHMatching(HParticleTrackSorter::kUseRKRICHWindow,4.); // select matching RICH-MDC for selectLeptons() function
	//sorter.setUseYMatching(kTRUE);                                // use meta cell build in select functions
	//sorter.setMetaBoundary(4.0);                                   // match in 4 mm to meta cell
	//sorter.setUseBeta(kTRUE);                                     // require beta >0 in build in select functions
	//sorter.setBetaLeptonCut(0.9);                                 // lower beta cut for lepton select
	//sorter.setUseFakeRejection(kTRUE);                            // reject fakes in  build in select functions (default kTRUE)
	//sorter.setIgnoreInnerMDC();                                   // do not reject Double_t inner MDC hits
	//sorter.setIgnoreOuterMDC();                                   // do not reject Double_t outer MDC hits
	//sorter.setIgnoreMETA();                                       // do not reject Double_t META hits
	//sorter.setIgnorePreviousIndex();                              // do not reject indices from previous selctions
	sorter.init();                                                  // get catgegory pointers etc...

	//t0Reco.setUseFlaggedCandidates(kTRUE);
	//t0Reco.init();

        evtChara.reinit();
	//--------------------------------------------------------------------------------------------

	const char *effcorFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/ALLMY_OC_matrix_apr12_day108_gen9_4mbins.root";
//	const char *effcorFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/MY_KpIII_effcor_flow_matrix_apr12_day108_gen8_4mbins.root";
	effcorFile = TFile::Open( effcorFileName );
	effPiM_cc3 = (TH2F *)effcorFile->Get("heff4_id9");
	effPiP_cc3 = (TH2F *)effcorFile->Get("heff4_id8");
	effKp_cc3  = (TH2F *)effcorFile->Get("heff4_id11");
//	effKm_cc3  = (TH2F *)effcorFile->Get("heff4_id12");
	effP_cc3   = (TH2F *)effcorFile->Get("heff4_id14");
	effPiM_cc2 = (TH2F *)effcorFile->Get("heff3_id9");
	effPiP_cc2 = (TH2F *)effcorFile->Get("heff3_id8");
	effKp_cc2  = (TH2F *)effcorFile->Get("heff3_id11");
//	effKm_cc2  = (TH2F *)effcorFile->Get("heff3_id12");
	effP_cc2   = (TH2F *)effcorFile->Get("heff3_id14");
	effPiM_cc1 = (TH2F *)effcorFile->Get("heff2_id9");
	effPiP_cc1 = (TH2F *)effcorFile->Get("heff2_id8");
	effKp_cc1  = (TH2F *)effcorFile->Get("heff2_id11");
//	effKm_cc1  = (TH2F *)effcorFile->Get("heff2_id12");
	effP_cc1   = (TH2F *)effcorFile->Get("heff2_id14");
	effPiM_cc0 = (TH2F *)effcorFile->Get("heff1_id9");
	effPiP_cc0 = (TH2F *)effcorFile->Get("heff1_id8");
	effKp_cc0  = (TH2F *)effcorFile->Get("heff1_id11");
//	effKm_cc0  = (TH2F *)effcorFile->Get("heff1_id12");
	effP_cc0   = (TH2F *)effcorFile->Get("heff1_id14");

	const char *filenameAccCorrRecEff = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/AccCorr_and_RecoEff_fineBinning.root";
	fileAccCorrRecEff = TFile::Open( filenameAccCorrRecEff );
	hAccCorr_cc0 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusAccCorr_centrality_0_10");
	hAccCorr_cc1 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusAccCorr_centrality_10_20");
	hAccCorr_cc2 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusAccCorr_centrality_20_30");
	hAccCorr_cc3 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusAccCorr_centrality_30_40");
	hRecoEff_cc0 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusRecoEff_centrality_0_10");
	hRecoEff_cc1 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusRecoEff_centrality_10_20");
	hRecoEff_cc2 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusRecoEff_centrality_20_30");
	hRecoEff_cc3 = (TH2F*) fileAccCorrRecEff->Get("kaonPlus/hKplusRecoEff_centrality_30_40");

	const char *tofcutFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/Kp_cutfile_tof_gen8a.root";
	mdcdEdx_cutfile_new_tof = TFile::Open( tofcutFileName );
	mdcdEdx_tof_Kp_data_new =(TCutG*)mdcdEdx_cutfile_new_tof->Get("KP_mdcdEdx_cut_scaled");

	const char *tofcutFileName2 = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/Kp_cutfile_tof_gen8a_15sigma.root";
	mdcdEdx_cutfile_new_tof2 = TFile::Open( tofcutFileName2 );
	mdcdEdx_tof_Km_data_new =(TCutG*)mdcdEdx_cutfile_new_tof2->Get("KP_mdcdEdx_cut_scaled");

	const char *rpccutFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/Kp_cutfile_rpc_gen8a_rpc25.root";
	mdcdEdx_cutfile_new_rpc = TFile::Open( rpccutFileName );
	mdcdEdx_rpc_Kp_data_new =(TCutG*)mdcdEdx_cutfile_new_rpc->Get("KP_mdcdEdx_cut_scaled");
	mdcdEdx_rpc_Km_data_new =(TCutG*)mdcdEdx_cutfile_new_rpc->Get("KP_mdcdEdx_cut_scaled");

	const char *kaonVetocutFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/betaTofdEdx_noKm.root";
	tofdEdx_cutfile = TFile::Open( kaonVetocutFileName );
	tofdEdx_beta_noK =(TCutG*)tofdEdx_cutfile->Get("betaTofdEdx_noKm");


	return kTRUE;
    }

    Int_t execute()
    {   // this function is called once per event.
	// if the function returns kSkipEvent the event
	// will be skipped a. for all following tasks
	// and b. not appear in output (hades eventLoop())

        if(!gHades) cout << "NO gHades !!!" << endl;
	if(gHades){
	    HCategory* evtInfoCat = (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo);
	    HCategory* candCat    = (HCategory*)HCategoryManager::getCategory(catParticleCand);

	    HEventHeader* header  = gHades->getCurrentEvent()->getHeader();
	    HVertex vertex      = header->getVertexReco();

	    HCategory* evtPlaneCat = (HCategory*)HCategoryManager::getCategory(catWallEventPlane);

#if isRealData && isExtendedOutput
	    HCategory* tofclstCat = (HCategory*)HCategoryManager::getCategory(catTofCluster);
#elif !isRealData
	    HGeantKine*       kine;
	    HLinearCategory* kineCat = (HLinearCategory*)HCategoryManager::getCategory(catGeantKine);
	    HGeantHeader* fGeantHeader = NULL;
#endif

            pluskaon.eventVerX = vertex.getX();
            pluskaon.eventVerY = vertex.getY();
	    pluskaon.eventVerZ = vertex.getZ();
	    pluskaon.evntVerChi2 = vertex.getChi2();

	    minuskaon.eventVerX = vertex.getX();
            minuskaon.eventVerY = vertex.getY();
	    minuskaon.eventVerZ = vertex.getZ();
	    minuskaon.evntVerChi2 = vertex.getChi2();

	    Float_t downflag = 1.;

            hEventCounter->Fill(cNumInputEv);

            if(!evtInfoCat) cout << "NO evtInfoCat !!!" << endl;
	    if(evtInfoCat)
	    {
		HParticleEvtInfo* evtInfo=0;
		evtInfo = HCategoryManager::getObject(evtInfo,evtInfoCat,0 );

		//HWallEventPlane* event_plane;
		//event_plane = HCategoryManager::getObject(event_plane,evtPlaneCat,0);

		Bool_t isGoodEvent = kFALSE;
#if isRealData || isEmbedding
		//nice plots of z-ver vs N_TofRpc
/*		if(evtInfo){
		    hZVerMult_noSelect->Fill( vertex.getZ(), evtInfo->getSumTofMultCut()+evtInfo->getSumRpcMultHitCut() );
                    if(evtInfo->isGoodEvent(Particle::kGoodTRIGGER))
			hZVerMult_triggerSelect->Fill( vertex.getZ(), evtInfo->getSumTofMultCut()+evtInfo->getSumRpcMultHitCut() );

		    if(evtInfo->isGoodEvent(Particle::kGoodTRIGGER|
					    Particle::kGoodVertexClust|
					    Particle::kGoodVertexCand|
					    Particle::kGoodSTART|
					    Particle::kNoPileUpSTART|
					    Particle::kNoVETO|
					    Particle::kGoodSTARTVETO|
					    Particle::kGoodSTARTMETA)
		      ) hZVerMult_allSelect->Fill( vertex.getZ(), evtInfo->getSumTofMultCut()+evtInfo->getSumRpcMultHitCut() );

		    if(evtInfo->isGoodEvent(Particle::kGoodTRIGGER|
					    Particle::kGoodSTART|
					    Particle::kNoPileUpSTART|
					    Particle::kNoVETO|
					    Particle::kGoodSTARTVETO|
					    Particle::kGoodSTARTMETA)
		      ) hZVerMult_someSelect->Fill( vertex.getZ(), evtInfo->getSumTofMultCut()+evtInfo->getSumRpcMultHitCut() );
		}
*/
		//--------------------- Event selection --------------------------------
		/*if(!evtInfo->isGoodEvent(Particle::kGoodTRIGGER) ) return kSkipEvent;
		hEventCounter->Fill(cisTriggerPT3);
		if(!evtInfo->isGoodEvent(Particle::kGoodSTART)) return kSkipEvent;
		hEventCounter->Fill(cselectStart);
		if(!evtInfo->isGoodEvent(Particle::kNoPileUpSTART)) return kSkipEvent;
		hEventCounter->Fill(cselectStartPile);
		if(!evtInfo->isGoodEvent(Particle::kGoodVertexClust)) return kSkipEvent;
		hEventCounter->Fill(cisGoodClusterVertex);
		if(!evtInfo->isGoodEvent(Particle::kGoodVertexCand)) return kSkipEvent;
		hEventCounter->Fill(cisGoodVertexCand);
		if(!evtInfo->isGoodEvent(Particle::kNoVETO)) return kSkipEvent;
		hEventCounter->Fill(cisVeto);
		if(!evtInfo->isGoodEvent(Particle::kGoodSTARTVETO)) return kSkipEvent;
		hEventCounter->Fill(cisStartVeto);
		if(!evtInfo->isGoodEvent(Particle::kGoodSTARTMETA)) return kSkipEvent;
		hEventCounter->Fill(cisStartMeta);
*/
                isGoodEvent = kTRUE;
//		if(evtInfo && evtInfo->isGoodEvent(Particle::kGoodTRIGGER|
//						   Particle::kGoodVertexClust|
//						   Particle::kGoodVertexCand|
//						   Particle::kGoodSTART|
//						   Particle::kNoPileUpSTART|
//						   Particle::kNoVETO|
//						   Particle::kGoodSTARTVETO|
//						   Particle::kGoodSTARTMETA
//						  )
//		  ) isGoodEvent = kTRUE;
#else
		if(evtInfo) isGoodEvent = kTRUE;
#endif
                if(!isGoodEvent) cout << "NO GoodEvent = NO evtInfo !!!" << endl;
                if(isGoodEvent)
		{
		    Double_t mult = evtInfo->getSumTofMultCut()+evtInfo->getSumRpcMultHitCut();
//cout << "Event mult: " << mult << endl;
#if !isRealData
		    if(!(fGeantHeader = gLoop->getGeantHeader(kFALSE))){ cout << "NO HGEANT!!! " << endl; return 0; }

		    Double_t geantEventPlaneAngle = fGeantHeader->getEventPlane();
		    Double_t geantImpactParameter = fGeantHeader->getImpactParameter();

//		    cout << "EP: " << geantEventPlaneAngle  << "deg.    IP:  " <<  geantImpactParameter << " fm" <<endl;
#endif


		    Int_t behruzCentralityClass = 0;
		    Int_t multClassNo = 0;
                    Float_t phiRP=-10000., phiRP_A=-10000., phiRP_B=-10000.;

#if isEmbedding
		    behruzCentralityClass = evtChara.getCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k5);  //is saved in tree!!

		    multClassNo = evtChara.getCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k10);  //needed for occupancy corr.

		    phiRP = evtChara.getEventPlane(HParticleEvtCharaBK::kDefault)*TMath::RadToDeg();
		    phiRP_A = evtChara.getEventPlane(HParticleEvtCharaBK::kDefault,1)*TMath::RadToDeg();
		    phiRP_B = evtChara.getEventPlane(HParticleEvtCharaBK::kDefault,2)*TMath::RadToDeg();
#else
		    behruzCentralityClass = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k5);  //is saved in tree!!
//cout << "Centrality: " << behruzCentralityClass << endl;
		    multClassNo = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10);  //needed for occupancy corr.

		    if(evtChara.getEventPlane(HParticleEvtChara::kDefault) == -1){
#if isRealData
			return 1;
#else
			phiRP = geantEventPlaneAngle;
			phiRP_A = geantEventPlaneAngle;
			phiRP_B = geantEventPlaneAngle;
#endif
		    }
		    else{
			phiRP = evtChara.getEventPlane(HParticleEvtChara::kDefault)*TMath::RadToDeg();
			phiRP_A = evtChara.getEventPlane(HParticleEvtChara::kDefault,1)*TMath::RadToDeg();
			phiRP_B = evtChara.getEventPlane(HParticleEvtChara::kDefault,2)*TMath::RadToDeg();
		    }
//		    cout << "Behruz: RP= " << phiRP << " deg.    centrclass= " << behruzCentralityClass << endl;
#endif

		    Float_t  event_weight = evtChara.getEventWeight();

                    if(multClassNo>0 && multClassNo<5) hEventCounter->Fill(cisGoodCentrality);

                    Float_t phiRP_tree = -10000., phiRPAB_tree = -10000.;

		    if( phiRP >= 0. )
		    {
			if(phiRP > 180.) phiRP -= 360.;
                        phiRP_tree = phiRP;
		    }
		    if( phiRP_A >= 0. && phiRP_B >= 0. ){
			phiRPAB_tree = TMath::Abs( phiRP_A - phiRP_B );
			if(phiRPAB_tree > 180.) phiRPAB_tree = 360. - phiRPAB_tree;    //this makes sure that phiAB is in (0,180)
		    }

#if isRealData
		    //my NEW centrality class -> see UsefulFunctions.h for more info
		    vector<int> myVecCC = getMyCentralityClass( behruzCentralityClass );
		    if( myVecCC.size() == 0 ) return 1; //skipping event outside of region of interest in centrality
		    for(int j=0; j<myVecCC.size(); j++){
			heventplane[myVecCC[j]]->Fill(phiRP_tree);
			hphiAB[myVecCC[j]]->Fill(phiRPAB_tree);
		    }
#endif

                    //fill counter with appropriate weight
		    hCentralityCounter->Fill(behruzCentralityClass,event_weight);


		    pluskaon.eventPlane      = phiRP_tree;
		    pluskaon.phiAB           = phiRPAB_tree;
		    pluskaon.multiplicity    = mult;
		    pluskaon.centralityclass = behruzCentralityClass;
                    pluskaon.eventWeight     = event_weight;

		    minuskaon.eventPlane      = phiRP_tree;
		    minuskaon.phiAB           = phiRPAB_tree;
		    minuskaon.multiplicity    = mult;
		    minuskaon.centralityclass = behruzCentralityClass;
                    minuskaon.eventWeight     = event_weight;

#if !isRealData
		    gkplus.eventPlane = phiRP_tree;
		    gkplus.eventPlane_gkine = geantEventPlaneAngle;
		    gkplus.phiAB = phiRPAB_tree;
		    gkplus.multiplicity = mult;
		    gkplus.centralityclass = behruzCentralityClass;
		    gkplus.impactparam = geantImpactParameter;

		    gkminus.eventPlane = phiRP_tree;
		    gkminus.eventPlane_gkine = geantEventPlaneAngle;
		    gkminus.phiAB = phiRPAB_tree;
		    gkminus.multiplicity = mult;
		    gkminus.centralityclass = behruzCentralityClass;
		    gkminus.impactparam = geantImpactParameter;
#endif



#if !isRealData
       //---------------------------------- Kine --------------------------------------------------------
       Int_t particleID = -1;

       //---------------------------- KPlus -------------------------------
       particleID = 11;
       Bool_t goodKPlusFound = kFALSE;
       for(Int_t j = 0; j < kineCat->getEntries(); j ++) {
	   kine = HCategoryManager::getObject(kine,kineCat,j);
	   if( kine->isPrimary() && kine->getID() == particleID && kine->getGeneratorInfo() == (500+particleID) ){
               goodKPlusFound = kTRUE;

               gkplus.invMass = HPhysicsConstants::mass(particleID);
	       gkplus.pt   = kine->getTransverseMomentum();
	       gkplus.mt   = sqrt(gkplus.pt*gkplus.pt+gkplus.invMass*gkplus.invMass);
               gkplus.mtm0 = sqrt(gkplus.pt*gkplus.pt+gkplus.invMass*gkplus.invMass) - HPhysicsConstants::mass(particleID);
	       gkplus.y    = kine->getRapidity();
               gkplus.ycm  = kine->getRapidity() - 0.74;
	       gkplus.y0  = kine->getRapidity()/0.74 - 1.0;
               gkplus.system =  kine->getSystem();

               kine->getMomentum(gkplus.momX,gkplus.momY,gkplus.momZ);
	       gkplus.mom = kine->getTotalMomentum();

	       Int_t sec_kine;
	       Float_t phi_kine, theta_kine;
	       sec_kine = kine->getPhiThetaDeg(theta_kine,phi_kine);

	       gkplus.theta = theta_kine;
               if( phi_kine > 180. ) phi_kine -= 360.;
	       gkplus.phi = phi_kine;

	       gkplus.flow = phi_kine - phiRP_tree;
               gkplus.flow_gkine = phi_kine - geantEventPlaneAngle;
               if(gkplus.flow < -180.) gkplus.flow += 360.;
               if(gkplus.flow > +180.) gkplus.flow -= 360.;
               if(gkplus.flow_gkine < -180.) gkplus.flow_gkine += 360.;
               if(gkplus.flow_gkine > +180.) gkplus.flow_gkine -= 360.;

	       gkplus.weight = 1.;

	       vector<HGeantKine *> kaonDaughters;
	       Int_t numberOfDaughters = kine->getDaughters(kine,kaonDaughters);
	       if(numberOfDaughters>0) gkplus.hasDecayed = kTRUE;
               else gkplus.hasDecayed = kFALSE;


	       Bool_t test_kp_acc = kFALSE;
               if( kine->isInAcceptanceBit(4,4,4,4,1,1) && !kine->isAtAnyMdcEdge(2) ) test_kp_acc = kTRUE;

	       gkplus.kAcc = test_kp_acc;


	       Bool_t test_kp_reco = kFALSE;
	       HParticleCandSim* cand=0;
	       if(candCat){
		   Int_t size = candCat->getEntries();

		   for(Int_t j = 0; j < size; j++)
		   {
		       cand = HCategoryManager::getObject(cand,candCat,j);

		       if( cand->getGeantTrack() == kine->getTrack() ){
#if isEmbedding
		           gkplus.pcand_chi2 = cand->getChi2();
		           gkplus.pcand_innerchi2 = cand->getInnerSegmentChi2();
		           gkplus.pcand_mm = cand->getMetaMatchQuality();
		           gkplus.pcand_mass = cand->getMass();
		           gkplus.pcand_system = cand->getSystemUsed();
                           gkplus.pcand_tofeloss = cand->getTofdEdx();
                           gkplus.pcand_mdceloss = cand->getMdcdEdx();
                           gkplus.pcand_beta = cand->getBeta();
                           gkplus.pcand_mom = cand->getMomentum();
			   gkplus.pcand_charge = cand->getCharge();
			   gkplus.pcand_mdcedge = cand->isAtAnyMdcEdge();
			   gkplus.pcand_goodflags = cand->isFlagAND(4,Particle::kIsAcceptedHitInnerMDC,Particle::kIsAcceptedHitOuterMDC,Particle::kIsAcceptedHitMETA,Particle::kIsAcceptedRK);
#endif
		       }

		       Bool_t test=kFALSE;
		       if( cand->isFlagAND(4,
					   Particle::kIsAcceptedHitInnerMDC,
					   Particle::kIsAcceptedHitOuterMDC,
					   Particle::kIsAcceptedHitMETA,
					   Particle::kIsAcceptedRK
					  )
			  &&
			  cand->getInnerSegmentChi2() > 0
			  &&
			  cand->getChi2()             < 1000
			 ) test=kTRUE;
		       if(!test) continue;

		       if(cand->isAtAnyMdcEdge()) continue;
		       if(cand->getSystemUsed() == -1) continue;

		       if( cand->getGeantTrack() == kine->getTrack() )
			   test_kp_reco = kTRUE;
		   }
	       }
               gkplus.kRec = test_kp_reco;


	       gkplus.MDCtracks = -10000.;
	       gkplus.METAclst = -10000.;

               break; //if one good is found take the first one (TODO)
	   }
       }
       if(goodKPlusFound) gkplustree->Fill();
       else{
	   gkplus.eventPlane = -10000.;
	   gkplus.eventPlane_gkine = -10000.;
	   gkplus.phiAB = -10000.;
	   gkplus.system = -10000.;
	   gkplus.theta = -10000.;
	   gkplus.phi = -10000.;
	   gkplus.momX = -10000.;
	   gkplus.momY = -10000.;
	   gkplus.momZ = -10000.;
	   gkplus.mom = -10000.;
	   gkplus.multiplicity = -10000.;
	   gkplus.centralityclass = -10000.;
	   gkplus.impactparam = -10000.;
	   gkplus.invMass = -10000.;
	   gkplus.mt = -10000.;
	   gkplus.mtm0 = -10000.;
	   gkplus.pt = -10000.;
	   gkplus.y = -10000.;
	   gkplus.ycm = -10000.;
	   gkplus.y0 = -10000.;
	   gkplus.flow = -10000.;
	   gkplus.flow_gkine = -10000.;
	   gkplus.weight = -10000.;
	   gkplus.kAcc = kFALSE;
	   gkplus.kRec = kFALSE;
	   gkplus.hasDecayed = kFALSE;
	   gkplus.MDCtracks = -10000.;
	   gkplus.METAclst = -10000.;
#if isEmbedding
           gkplus.pcand_chi2 = -10000.;
	   gkplus.pcand_innerchi2 = -10000.;
           gkplus.pcand_mm = -10000.;
           gkplus.pcand_mass = -10000.;
           gkplus.pcand_system = -10000.;
           gkplus.pcand_tofeloss = -10000.;
           gkplus.pcand_mdceloss = -10000.;
           gkplus.pcand_beta = -10000.;
           gkplus.pcand_mom = -10000.;
	   gkplus.pcand_charge = -10000.;
	   gkplus.pcand_mdcedge = kFALSE;
           gkplus.pcand_goodflags = kFALSE;
#endif

           gkplustree->Fill();
       }

       //---------------------------- KMinus -------------------------------
       particleID = 12;
       Bool_t goodKMinusFound = kFALSE;
       for(Int_t j = 0; j < kineCat->getEntries(); j ++) {
	   kine = HCategoryManager::getObject(kine,kineCat,j);
	   if( kine->isPrimary() && kine->getID() == particleID && kine->getGeneratorInfo() == (500+particleID) ){
               goodKMinusFound = kTRUE;

               gkminus.invMass = HPhysicsConstants::mass(particleID);
	       gkminus.pt   = kine->getTransverseMomentum();
	       gkminus.mt   = sqrt(gkminus.pt*gkminus.pt+gkminus.invMass*gkminus.invMass);
               gkminus.mtm0 = sqrt(gkminus.pt*gkminus.pt+gkminus.invMass*gkminus.invMass) - HPhysicsConstants::mass(particleID);
	       gkminus.y    = kine->getRapidity();
               gkminus.ycm  = kine->getRapidity() - 0.74;
	       gkminus.y0  = kine->getRapidity()/0.74 - 1.0;
               gkminus.system =  kine->getSystem();

               kine->getMomentum(gkminus.momX,gkminus.momY,gkminus.momZ);
	       gkminus.mom = kine->getTotalMomentum();

	       Int_t sec_kine;
	       Float_t phi_kine, theta_kine;
	       sec_kine = kine->getPhiThetaDeg(theta_kine,phi_kine);

	       gkminus.theta = theta_kine;
               if( phi_kine > 180. ) phi_kine -= 360.;
	       gkminus.phi = phi_kine;

	       gkminus.flow = phi_kine - phiRP_tree;
               gkminus.flow_gkine = phi_kine - geantEventPlaneAngle;
               if(gkminus.flow < -180.) gkminus.flow += 360.;
               if(gkminus.flow > +180.) gkminus.flow -= 360.;
               if(gkminus.flow_gkine < -180.) gkminus.flow_gkine += 360.;
               if(gkminus.flow_gkine > +180.) gkminus.flow_gkine -= 360.;

	       gkminus.weight = 1.;

	       vector<HGeantKine *> kaonDaughters;
	       Int_t numberOfDaughters = kine->getDaughters(kine,kaonDaughters);
	       if(numberOfDaughters>0) gkminus.hasDecayed = kTRUE;
               else gkminus.hasDecayed = kFALSE;

               Bool_t test_km_acc = kTRUE;
	       if( kine->isAtAnyMdcEdge(2) ) test_km_acc = kFALSE;
	       for(int imm=0; imm<4; imm++)
		   if( kine->getNLayerMod(imm) < 4 ) test_km_acc = kFALSE;

	       if( !( kine->getSys(0) || kine->getSys(1) ) ) test_km_acc = kFALSE;
	       gkminus.kAcc = test_km_acc;


	       Bool_t test_km_reco = kFALSE;
	       HParticleCandSim* cand=0;
	       if(candCat){
		   Int_t size = candCat->getEntries();

		   for(Int_t j = 0; j < size; j++)
		   {
		       cand = HCategoryManager::getObject(cand,candCat,j);

		       Bool_t test=kFALSE;
		       if( cand->isFlagAND(4,
					   Particle::kIsAcceptedHitInnerMDC,
					   Particle::kIsAcceptedHitOuterMDC,
					   Particle::kIsAcceptedHitMETA,
					   Particle::kIsAcceptedRK
					  )
			  &&
			  cand->getInnerSegmentChi2() > 0
			  &&
			  cand->getChi2()             < 1000
			 ) test=kTRUE;
		       if(!test) continue;

		       if(cand->isAtAnyMdcEdge()) continue;
		       if(cand->getSystemUsed() == -1) continue;

		       if( cand->getGeantTrack() == kine->getTrack() ){
			   test_km_reco = kTRUE;
		       }
		   }
	       }
               gkminus.kRec = test_km_reco;


	       gkminus.MDCtracks = -10000.;
	       gkminus.METAclst = -10000.;
	   }
       }
       if(goodKMinusFound) gkminustree->Fill();
       else{
	   gkminus.eventPlane = -10000.;
	   gkminus.eventPlane_gkine = -10000.;
	   gkminus.phiAB = -10000.;
	   gkminus.system = -10000.;
	   gkminus.theta = -10000.;
	   gkminus.phi = -10000.;
	   gkminus.momX = -10000.;
	   gkminus.momY = -10000.;
	   gkminus.momZ = -10000.;
	   gkminus.mom = -10000.;
	   gkminus.multiplicity = -10000.;
	   gkminus.centralityclass = -10000.;
	   gkminus.impactparam = -10000.;
	   gkminus.invMass = -10000.;
	   gkminus.mt = -10000.;
	   gkminus.mtm0 = -10000.;
	   gkminus.pt = -10000.;
	   gkminus.y = -10000.;
	   gkminus.ycm = -10000.;
	   gkminus.y0 = -10000.;
	   gkminus.flow = -10000.;
	   gkminus.flow_gkine = -10000.;
	   gkminus.weight = -10000.;
	   gkminus.kAcc = kFALSE;
	   gkminus.kRec = kFALSE;
	   gkminus.hasDecayed = kFALSE;
	   gkminus.MDCtracks = -10000.;
	   gkminus.METAclst = -10000.;

           gkminustree->Fill();
       }

#endif

//cout << "number of pcands: " << candCat->getEntries() << endl;
      //--------------------------------------------------------------------------------------------------
		    // loop over particle candidates in event
		    if(candCat){
			TString filename;
			//-------------------------------------------------
			// get the sector list infos for HLoop
			Int_t sectors [6];
			gLoop->getSectors(sectors); // fill sector array
//cout << "test sector 1: " << sectors[1] << endl;

			//------------------------------------------------------------------------
			// prepare track sorting
			// clean vectors and index arrays
			sorter.cleanUp();
			sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);
			Int_t nCandHad     = sorter.fill(HParticleTrackSorter::selectHadrons);
			Int_t nCandHadBest = sorter.selectBest(Particle::kIsBestRKSorter,Particle::kIsHadronSorter);
//			Int_t nCandHadBest = sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsHadron);
//cout << "selected hadrons: " << nCandHadBest << endl;
			if( nCandHad<0 || nCandHadBest<0 ) return 1;

			vector<Int_t> IC_kp;
			vector<Int_t> IC_km;
#if isRealData
			vector<HParticleCand *> vkp;
			vector<HParticleCand *> vkm;
			HParticleCand* cand=0;
  #if isExtendedOutput
			HTofCluster* tofClst=0;
  #endif
#else
			vector<HParticleCandSim *> vkp;
			vector<HParticleCandSim *> vkm;
			HParticleCandSim* cand=0;
                        HGeantKine *kineRec=0;
#endif

			Int_t size = candCat->getEntries();
//                        cout << "PCand size= " << size << endl;
			for(Int_t j = 0; j < size; j++)
			{
			    cand = HCategoryManager::getObject(cand,candCat,j);
//#if !isRealData
//                          if(cand->getGeantPID()==11) cout << "K+ found" << endl;
//			    kineRec = HCategoryManager::getObject(kineRec,kineCat,cand->getGeantTrack()-1);
//			    if( !( (kine->isPrimary() && kine->getID()==11 && kine->getGeneratorInfo()==511) ||
//				   (kine->isPrimary() && kine->getID()==12 && kine->getGeneratorInfo()==512) ) )
//				continue;
//#endif

			    if(cand) {
                                hTrackCounter->Fill(cAllTracks);

				if(!gLoop->goodSector(cand->getSector())) continue;  // skipp inactive sectors
                                hTrackCounter->Fill(cisGoodSector);

				if(!cand->isFlagBit(kIsUsed)) continue;
                                hTrackCounter->Fill(cisUsed);

				if(cand->isAtAnyMdcEdge()) continue;
                                hTrackCounter->Fill(cisAtAnyEdge);

				if(cand->getSystemUsed() == -1) continue;
				hTrackCounter->Fill(cisGoodSystem);

				if(cand->getChi2()<=400){
				    hTrackCounter->Fill(cChisqRK400);

                                    if(cand->getMetaMatchQuality()<=3) hTrackCounter->Fill(cChisqRK400MM3);
				}
				if(cand->getChi2()<=100){
				    hTrackCounter->Fill(cChisqRK100);

                                    if(cand->getMetaMatchQuality()<=2) hTrackCounter->Fill(cChisqRK100MM2);
				}

				if(cand->getCharge()>0 && cand->getMass()>50. && cand->getMass()<300.){
				    hChisqRK_pip->Fill( cand->getChi2() );
				    hChisqMM_pip->Fill( cand->getMetaMatchQuality() );
				}
				if(cand->getCharge()>0 && cand->getMass()>300. && cand->getMass()<700.){
				    hChisqRK_kp->Fill( cand->getChi2() );
				    hChisqMM_kp->Fill( cand->getMetaMatchQuality() );
				}
				if(cand->getCharge()>0 && cand->getMass()>700. && cand->getMass()<1250.){
				    hChisqRK_proton->Fill( cand->getChi2() );
				    hChisqMM_proton->Fill( cand->getMetaMatchQuality() );
				}
				if(cand->getCharge()<0 && cand->getMass()>50. && cand->getMass()<300.){
				    hChisqRK_pim->Fill( cand->getChi2() );
				    hChisqMM_pim->Fill( cand->getMetaMatchQuality() );
				}
				if(cand->getCharge()<0 && cand->getMass()>300. && cand->getMass()<700.){
				    hChisqRK_km->Fill( cand->getChi2() );
				    hChisqMM_km->Fill( cand->getMetaMatchQuality() );
				}


				Int_t   pid = cand->getPID();
				Bool_t goodPBeta=false, goodPEloss=false;
				Bool_t goodMdcEloss=false, goodTofEloss=false, goodChi2RK=false, goodChi2MM=false, goodMass=false, goodMom=false, goodNLayers=false, goodDistVer=false;

				if(cand->getSystemUsed() == 1){
#if isRealData && isExtendedOutput
				    Int_t tofClstInd = cand->getTofClstInd();
                                    Int_t tofRod = -1;
				    if(tofClstInd >= 0){
					tofClst = HCategoryManager::getObject(tofClst,tofclstCat,tofClstInd);
					tofRod = 64*cand->getSector() + 8*tofClst->getModule() + tofClst->getCell();
				    }
				    if( (tofRod%8)==3 && cand->getSector()==0 &&
				       cand->getChi2()<=400 && cand->getMetaMatchQuality()<=3){
                                        hTofElossBeta_OneRod_all[tofRod/8]->Fill( cand->getBeta(), cand->getTofdEdx() );
                                        hTofElossMom_OneRod_all[tofRod/8]->Fill( cand->getCharge()*cand->getMomentum(), cand->getTofdEdx() );

					if(cand->getMass()>600.)
					    hTofElossBeta_OneRod_massCut[tofRod/8]->Fill( cand->getBeta(), cand->getTofdEdx() );
					    hTofElossMom_OneRod_massCut[tofRod/8]->Fill( cand->getCharge()*cand->getMomentum(), cand->getTofdEdx() );
				    }
#endif

#if !isRealData
				    if(cand->getGeantPID()==11 && cand->getGeantParentTrackNum()==0)
                                        hMdcEloss_sim_TOF->Fill( cand->getCharge()*cand->getMomentum(), cand->getMdcdEdx() );
#endif
				    if(cand->getChi2()<=400 && cand->getMetaMatchQuality()<=3){
					hTofElossMom->Fill( cand->getCharge()*cand->getMomentum(), cand->getTofdEdx() );
					hBetaMom_TOF->Fill( cand->getCharge()*cand->getMomentum(), cand->getBeta() );
					hMdcEloss_TOF->Fill( cand->getCharge()*cand->getMomentum(), cand->getMdcdEdx() );
					hMassMom_TOF->Fill( cand->getCharge()*cand->getMass(), cand->getMomentum() );
					if(cand->getCharge()<0) hTofElossBeta_neg->Fill( cand->getBeta(), cand->getTofdEdx() );

//                                      Int_t id_eloss = -1;
//					if( cand->getMass()>50. && cand->getMass()<300.) id_eloss = 8;
//                                        else if( cand->getMass()>700. && cand->getMass()<1200.) id_eloss = 14;
//                                        else if( cand->getMass()>1250. && cand->getMass()<1550.) id_eloss = 49;
//					else if( cand->getMass()>1600. && cand->getMass()<2150.) id_eloss = 45;
//
//					if( id_eloss!=-1 ){
//					    fEnergyLossMdc->SetParameter(0,id_eloss);
//					    Double_t eloss_theory = fEnergyLossMdc->Eval( cand->getMomentum() );
//					    if( eloss_theory!= 0){
//						Double_t Z_MdcEloss = TMath::Log( cand->getMdcdEdx() / eloss_theory );
//                                                hZMdcEloss->Fill( cand->getMass(), Z_MdcEloss );
//					    }
//					}

					fEnergyLossMdc->SetParameter(0,14);
					Double_t Z_MdcEloss = TMath::Log( cand->getMdcdEdx() / fEnergyLossMdc->Eval( cand->getMomentum() ) );
                                        hZMdcEloss->Fill( cand->getMass(), Z_MdcEloss );
				    }

                                    hMass_TOF->Fill( cand->getCharge()*cand->getMass() );
				    hKaonMass_TOF[0]->Fill( cand->getCharge()*cand->getMass() );
				    if(cand->getChi2()<=100){
					hKaonMass_TOF[1]->Fill( cand->getCharge()*cand->getMass() );

					if(cand->getMetaMatchQuality()<=2){
					    hKaonMass_TOF[2]->Fill( cand->getCharge()*cand->getMass() );

					    if(cand->getMass()>300. && cand->getMass()<700.)
						hMdcEloss_dataKcut_TOF->Fill(cand->getCharge()*cand->getMomentum(),cand->getMdcdEdx());


					    if( mdcdEdx_tof_Kp_data_new->IsInside(cand->getMomentum(),cand->getMdcdEdx())
					       || mdcdEdx_tof_Km_data_new->IsInside(cand->getMomentum(),cand->getMdcdEdx()) ){
						hKaonMass_TOF[3]->Fill( cand->getCharge()*cand->getMass() );

						if( !tofdEdx_beta_noK->IsInside(cand->getBeta(),cand->getTofdEdx()) )
						    hKaonMass_TOF[4]->Fill( cand->getCharge()*cand->getMass() );
					    }
					}
				    }
				}
				if(cand->getSystemUsed() == 0){
#if !isRealData
				    if(cand->getGeantPID()==11 && cand->getGeantParentTrackNum()==0)
                                        hMdcEloss_sim_RPC->Fill( cand->getCharge()*cand->getMomentum(), cand->getMdcdEdx() );
#endif

				    if(cand->getChi2()<=400 && cand->getMetaMatchQuality()<=3){
					hBetaMom_RPC->Fill( cand->getCharge()*cand->getMomentum(), cand->getBeta() );
					hMdcEloss_RPC->Fill( cand->getCharge()*cand->getMomentum(), cand->getMdcdEdx() );
					hMassMom_RPC->Fill( cand->getCharge()*cand->getMass(), cand->getMomentum() );

//					Int_t id_eloss = -1;
//					if( cand->getMass()>50. && cand->getMass()<300.) id_eloss = 8;
//                                        else if( cand->getMass()>700. && cand->getMass()<1200.) id_eloss = 14;
//                                        else if( cand->getMass()>1250. && cand->getMass()<1550.) id_eloss = 49;
//					else if( cand->getMass()>1600. && cand->getMass()<2150.) id_eloss = 45;
//
//					if( id_eloss!=-1 ){
//					    fEnergyLossMdc->SetParameter(0,id_eloss);
//					    Double_t eloss_theory = fEnergyLossMdc->Eval( cand->getMomentum() );
//					    if( eloss_theory!= 0){
//						Double_t Z_MdcEloss = TMath::Log( cand->getMdcdEdx() / eloss_theory );
//                                                hZMdcEloss->Fill( cand->getMass(), Z_MdcEloss );
//					    }
//					}

					fEnergyLossMdc->SetParameter(0,14);
					Double_t Z_MdcEloss = TMath::Log( cand->getMdcdEdx() / fEnergyLossMdc->Eval( cand->getMomentum() ) );
                                        hZMdcEloss->Fill( cand->getMass(), Z_MdcEloss );

				    }

				    hMass_RPC->Fill( cand->getCharge()*cand->getMass() );
				    hKaonMass_RPC[0]->Fill( cand->getCharge()*cand->getMass() );
				    if(cand->getChi2()<=100){
					hKaonMass_RPC[1]->Fill( cand->getCharge()*cand->getMass() );

					if(cand->getMetaMatchQuality()<=2){
					    hKaonMass_RPC[2]->Fill( cand->getCharge()*cand->getMass() );

					    if(cand->getMass()>300. && cand->getMass()<700.)
						hMdcEloss_dataKcut_RPC->Fill(cand->getCharge()*cand->getMomentum(),cand->getMdcdEdx());


					    if( mdcdEdx_rpc_Kp_data_new->IsInside(cand->getMomentum(),cand->getMdcdEdx())
					       || mdcdEdx_rpc_Km_data_new->IsInside(cand->getMomentum(),cand->getMdcdEdx()) )
						hKaonMass_RPC[3]->Fill( cand->getCharge()*cand->getMass() );
					}
				    }
				}

				switch(PID_option){
				case 1:
				    {
					//First method of particle identification - Georgy's default one in dst-files
					if(pid==11){
					    IC_kp.push_back(j);
#if isRealData
					    vkp.push_back(new HParticleCand(*cand));
#else
					    vkp.push_back(new HParticleCandSim(*cand));
#endif
					}
					if(pid==12){
					    IC_km.push_back(j);
#if isRealData
					    vkm.push_back(new HParticleCand(*cand));
#else
					    vkm.push_back(new HParticleCandSim(*cand));
#endif
					}
				    }
				    break;

				case 2:
				    {
					//track quality cuts
					if( cand->getChi2() <= trackQualityCuts[iTrackCuts][0] )
					    goodChi2RK = true;

					if( cand->getMetaMatchQuality() <= trackQualityCuts[iTrackCuts][1] )
					    goodChi2MM = true;

					//plus charge and mass check
					if( cand->getMass()>300 && cand->getMass()<700 )
					    goodMass = true;

					//TOF dEdx veto
					if( ( cand->getSystemUsed()==1 && !tofdEdx_beta_noK->IsInside(cand->getBeta(),cand->getTofdEdx()) )
					   || ( cand->getSystemUsed()==0 )
					  )
					    goodTofEloss = true;


					if( cand->getCharge()==1 ){
					    //MDC dEdx cut
					    if( ( cand->getSystemUsed()==1 && mdcdEdx_tof_Kp_data_new->IsInside(cand->getMomentum(),cand->getMdcdEdx()) )
					       || ( cand->getSystemUsed()==0 && mdcdEdx_rpc_Kp_data_new->IsInside(cand->getMomentum(),cand->getMdcdEdx()) )
					      )
						goodMdcEloss = true;

					    //plus momentum cut
					    if( cand->getMomentum()>0 && cand->getMomentum()<1200 )
						goodMom = true;
					}

					if( cand->getCharge()==-1 ){
					    //MDC dEdx cut
					    if( ( cand->getSystemUsed()==1 && mdcdEdx_tof_Km_data_new->IsInside(cand->getMomentum(),cand->getMdcdEdx()) )
					       || ( cand->getSystemUsed()==0 && mdcdEdx_rpc_Km_data_new->IsInside(cand->getMomentum(),cand->getMdcdEdx()) )
					      )
						goodMdcEloss = true;

					    //plus momentum cut
					    if( cand->getMomentum()>100 && cand->getMomentum()<800 )
						goodMom = true;
					}

					//Heidi's cuts together
					if( ( (trackElossCuts[iTrackCuts][0] && goodMdcEloss) || (!trackElossCuts[iTrackCuts][0]) ) &&
					    ( (trackElossCuts[iTrackCuts][1] && goodTofEloss) || (!trackElossCuts[iTrackCuts][1]) ) &&
					   goodChi2RK && goodChi2MM && goodMass && goodMom ){

					    if(cand->getCharge()==1){
						IC_kp.push_back(j);
#if isRealData
						vkp.push_back(new HParticleCand(*cand));
#else
						vkp.push_back(new HParticleCandSim(*cand));
#endif
					    }

					    if(cand->getCharge()==-1){
						IC_km.push_back(j);
#if isRealData
						vkm.push_back(new HParticleCand(*cand));
#else
						vkm.push_back(new HParticleCandSim(*cand));
#endif
					    }

					}
				    }
				    break;

				case 3:
				    {
					//Third method of particle identification - original one (just some mass cut)
					if( (cand->getCharge() > 0) && (cand->getChi2()<=400 && cand->getMetaMatchQuality()<=3) ){
					    Float_t mass = TMath::Sqrt( cand->getMass2());
					    if(mass > 460 && mass < 520){
						IC_kp.push_back(j);
#if isRealData
					    vkp.push_back(new HParticleCand(*cand));
#else
					    vkp.push_back(new HParticleCandSim(*cand));
#endif
					    }
					}
				    }
				    break;

				case 4:
				    {
					//Second method of particle identification - check with HParticleTool functions
					Float_t deltaEloss = 0., sigmaEloss = -1.;
					if( HParticleTool::isParticledEdx(11,cand,deltaEloss,sigmaEloss)) {
					    if( TMath::Abs(deltaEloss/sigmaEloss) < 2. )
						goodPEloss = true;
					}

					Float_t deltaTime = 0., sigmaTime = -1.;
					if( HParticleTool::isParticleBeta(11,cand,2.,0.,1500.,deltaTime,sigmaTime,"apr12"))
					    goodPBeta = true; //2. is sigma allowed, 0. minMom, 1500. maxMom

					//ParticleTools cuts together
					if( goodPEloss && goodPBeta ){
					    IC_kp.push_back(j);
#if isRealData
					    vkp.push_back(new HParticleCand(*cand));
#else
					    vkp.push_back(new HParticleCandSim(*cand));
#endif
					}

					goodPEloss=false;
                                        goodPBeta=false;
					deltaEloss = 0., sigmaEloss = -1.;
					if( HParticleTool::isParticledEdx(12,cand,deltaEloss,sigmaEloss)) {
					    if( TMath::Abs(deltaEloss/sigmaEloss) < 2. )
						goodPEloss = true;
					}

					deltaTime = 0., sigmaTime = -1.;
					if( HParticleTool::isParticleBeta(12,cand,2.,50.,900.,deltaTime,sigmaTime,"apr12"))
					    goodPBeta = true; //2. is sigma allowed, 0. minMom, 1500. maxMom

					//ParticleTools cuts together
					if( goodPEloss && goodPBeta ){
					    IC_km.push_back(j);
#if isRealData
					    vkm.push_back(new HParticleCand(*cand));
#else
					    vkm.push_back(new HParticleCandSim(*cand));
#endif
					}

				    }
				    break;
#if !isRealData
				case 5:
				    {
					if( cand->getGeantPID() == 11 ){
					    IC_kp.push_back(j);
					    vkp.push_back(new HParticleCandSim(*cand));
					}
					if( cand->getGeantPID() == 12 ){
					    IC_km.push_back(j);
					    vkm.push_back(new HParticleCandSim(*cand));
					}
				    }
                                    break;
#endif

				case 6:  //K. Piasecki
				    {
					//track quality cuts
					if( cand->getChi2() <= 100 )
					    goodChi2RK = true;

					if( cand->getMetaMatchQuality() <= 2 )
					    goodChi2MM = true;

					if( (cand->getNLayerMod(0)+cand->getNLayerMod(1)+cand->getNLayerMod(2)+cand->getNLayerMod(3))>19 )
					    goodNLayers = true;

					HGeomVector base, dir, eventVertex;
					eventVertex.setXYZ(vertex.getX(),vertex.getY(),vertex.getZ());
					CalcSegVector(cand->getZ(), cand->getR(), cand->getPhi()*TMath::DegToRad(), cand->getTheta()*TMath::DegToRad(), base, dir);
					Float_t distVer = calculateMinimumDistanceStraightToPoint(base, dir, eventVertex);
					if( distVer<= 20 )
                                            goodDistVer = true;

					//plus charge and mass check
					if( cand->getMass()>340 && cand->getMass()<660 )
					    goodMass = true;


					Float_t beta = (cand->getBeta()>1 ? 0. : cand->getBeta());
					Float_t gamma = 1./TMath::Sqrt(1-beta*beta);

					if( cand->getCharge()==1 ){
					    if( cand->getSystemUsed()==1 ){
						if( (cand->getMdcdEdx()>1.1) && (cand->getMdcdEdx()<30.) && (cand->getMdcdEdx()*cand->getMomentum()>500) && (cand->getMdcdEdx()*cand->getMomentum()<15000)
//						   && mdcdEdx_tof_Kp_data_new->IsInside(cand->getMomentum(),cand->getMdcdEdx())
						  )
						    goodMdcEloss = true;

						if( (cand->getTofdEdx()*beta*gamma >= 1.0) && (cand->getTofdEdx()*beta*gamma <= 2.2)
//						   && !tofdEdx_beta_noK->IsInside(cand->getBeta(),cand->getTofdEdx())
						  )
						    goodTofEloss = true;

						if( cand->getMomentum()>150. && cand->getMomentum()<900. )
                                                    goodMom = true;
					    }

					    if( cand->getSystemUsed()==0 ){
						goodTofEloss = true;

						if( (cand->getMdcdEdx()*cand->getMomentum()>550) && (cand->getMdcdEdx()*cand->getMomentum()<10000) && (cand->getMdcdEdx()<40.)
//						   && mdcdEdx_rpc_Kp_data_new->IsInside(cand->getMomentum(),cand->getMdcdEdx())
						  )
						    goodMdcEloss = true;

						if( cand->getMomentum()>200. && cand->getMomentum()<1200. )
                                                    goodMom = true;
					    }
					}

					if( cand->getCharge()==-1 ){
					    if( cand->getSystemUsed()==1 ){
						if( (cand->getMdcdEdx()>2.5) && (cand->getMdcdEdx()*cand->getMomentum()>500)
//						   && mdcdEdx_tof_Km_data_new->IsInside(cand->getMomentum(),cand->getMdcdEdx())
						  )
						    goodMdcEloss = true;

						if( (cand->getTofdEdx()*beta*gamma >= 1.0) && (cand->getTofdEdx()*beta*gamma <= 2.2)
//						   && !tofdEdx_beta_noK->IsInside(cand->getBeta(),cand->getTofdEdx())
						  )
						    goodTofEloss = true;

						if( cand->getMomentum()>250. && cand->getMomentum()<650. )
                                                    goodMom = true;
					    }

					    if( cand->getSystemUsed()==0 ){
						goodTofEloss = true;

						if( cand->getMdcdEdx()*cand->getMomentum()>800
//						   && mdcdEdx_rpc_Km_data_new->IsInside(cand->getMomentum(),cand->getMdcdEdx())
						  )
						    goodMdcEloss = true;

						if( cand->getMomentum()>250. && cand->getMomentum()<700. )
                                                    goodMom = true;
					    }
					}

					if( goodChi2RK && goodChi2MM && goodNLayers && goodDistVer && goodMass && goodMom && goodMdcEloss && goodTofEloss ){

					    if(cand->getCharge()==1){
						IC_kp.push_back(j);
#if isRealData
						vkp.push_back(new HParticleCand(*cand));
#else
						vkp.push_back(new HParticleCandSim(*cand));
#endif
					    }

					    if(cand->getCharge()==-1){
						IC_km.push_back(j);
#if isRealData
						vkm.push_back(new HParticleCand(*cand));
#else
						vkm.push_back(new HParticleCandSim(*cand));
#endif
					    }

					}
				    }
				    break;

				default:
				    cout << "This PID_option is not supported. Please select from 1 to 6." << endl;
                                    break;
				}

			    }
			} // end cand loop

                        //cout << "N(K+) = " << IC_kp.size() << endl;

			//K+ loop
			if( IC_kp.size() > 0 ){
#if isRealData
			    HParticleCand* kpCand;
#else
			    HParticleCandSim* kpCand;
#endif

			    TLorentzVector kplusLV(1.,1.,1.,1.);

			    for(unsigned int ikp=0; ikp<IC_kp.size(); ikp++){
				kpCand = HCategoryManager::getObject(kpCand,candCat,IC_kp[ikp]);
				HParticleTool::fillTLorentzVector(kplusLV,kpCand,11,true);

				Double_t flow     = kplusLV.Phi()*TMath::RadToDeg() - phiRP_tree;  // angle w.r.t. FW event plane
				if (flow<-180.) flow += 360.;
				if (flow>+180.) flow -= 360.;

				Double_t kpGammaInv = ( kpCand->getBeta()>=1. ? 0 : TMath::Sqrt(1-TMath::Power(kpCand->getBeta(),2)) );
				Double_t kpMass = ( kpCand->getCorrectedMomentumPID(11)*kpGammaInv )/kpCand->getBeta();

				pluskaon.invMass = kpMass;
                                pluskaon.mass = (kpCand->getMass2()<0 ? -1 : TMath::Sqrt(kpCand->getMass2()));
                                pluskaon.system = kpCand->getSystemUsed();
				pluskaon.theta = kplusLV.Theta()*TMath::RadToDeg();
				pluskaon.phi = ( kplusLV.Phi()*TMath::RadToDeg()>0 ? kplusLV.Phi()*TMath::RadToDeg() : kplusLV.Phi()*TMath::RadToDeg()+360.) ;
				pluskaon.mom = kplusLV.P();
				pluskaon.momX = kplusLV.Px();
				pluskaon.momY = kplusLV.Py();
				pluskaon.momZ = kplusLV.Pz();
				pluskaon.pt = kplusLV.Pt();
				pluskaon.mt = kplusLV.Mt();
				pluskaon.mtm0 = kplusLV.Mt() - HPhysicsConstants::mass(11);
				pluskaon.y = kplusLV.Rapidity();
				pluskaon.ycm = kplusLV.Rapidity() - 0.74;  //y_cm = y - 0.5* y_proj (y_proj = 1.48 from definition)
                                pluskaon.y0 = kplusLV.Rapidity()/0.74 - 1.; //y_0 = (y/y_proj)_cm = y_cm / y_proj,cm = (y - 0.5*y_proj)/(y_proj-0.5*y_proj) = (y - 0.5*y_proj)/(0.5*y_proj) = y/(0.5*y_proj) - 1
				pluskaon.flow = flow;

				// eff.  corrections - dependence on track multiplicity vs phi-phi_RP
				Int_t id = OccupCorr_pid;
				Float_t eff=0., weight=0.;

				eff=getEff(id,4-multClassNo,flow,kplusLV.Theta()*TMath::RadToDeg());
				if(eff!=0) weight=1./eff;
				if(weight>10||weight<0) weight=0.;
				pluskaon.occupancyCorr = weight;

				weight=0.;
				eff=getAccCorr(multClassNo-1,kplusLV.Pt(),kplusLV.Rapidity()/0.74 - 1.);
				if(eff!=0.) weight=1./eff;
				if(weight>10000. || weight<1.) weight=0.;
				pluskaon.accCorr = weight;

				weight=0.;
				eff=getRecoEff(multClassNo-1,kplusLV.Pt(),kplusLV.Rapidity()/0.74 - 1.);
				if(eff!=0.) weight=1./eff;
				if(weight>10000. || weight<1.) weight=0.;
				pluskaon.recoEff = weight;

				pluskaon.MDCtracks = -10000.;
				pluskaon.METAclst = -10000.;

                                pluskaon.downflag = downflag;

				if( pluskaon.invMass > 300 && pluskaon.invMass < 700 ){
				    pluskaontree->Fill();

                                    downflag = -1.;
				}

			    }
			}//end of K+ loop





			//K- loop
			if( IC_km.size() > 0 ){
#if isRealData
			    HParticleCand* kmCand;
#else
			    HParticleCandSim* kmCand;
#endif

			    TLorentzVector kminusLV(1.,1.,1.,1.);

			    for(unsigned int ikm=0; ikm<IC_km.size(); ikm++){
				kmCand = HCategoryManager::getObject(kmCand,candCat,IC_km[ikm]);
				HParticleTool::fillTLorentzVector(kminusLV,kmCand,12,true);

				Double_t flow     = kminusLV.Phi()*TMath::RadToDeg() - phiRP_tree;  // angle w.r.t. FW event plane
				if (flow<-180.) flow += 360.;
				if (flow>+180.) flow -= 360.;

				Double_t kmGammaInv = ( kmCand->getBeta()>=1. ? 0 : TMath::Sqrt(1-TMath::Power(kmCand->getBeta(),2)) );
				Double_t kmMass = ( kmCand->getCorrectedMomentumPID(12)*kmGammaInv )/kmCand->getBeta();

				minuskaon.invMass = kmMass;
                                minuskaon.mass = (kmCand->getMass2()<0 ? -1 : TMath::Sqrt(kmCand->getMass2()));
                                minuskaon.system = kmCand->getSystemUsed();
				minuskaon.theta = kminusLV.Theta()*TMath::RadToDeg();
				minuskaon.phi = ( kminusLV.Phi()*TMath::RadToDeg()>0 ? kminusLV.Phi()*TMath::RadToDeg() : kminusLV.Phi()*TMath::RadToDeg()+360.) ;
				minuskaon.mom = kminusLV.P();
				minuskaon.momX = kminusLV.Px();
				minuskaon.momY = kminusLV.Py();
				minuskaon.momZ = kminusLV.Pz();
				minuskaon.pt = kminusLV.Pt();
				minuskaon.mt = kminusLV.Mt();
				minuskaon.mtm0 = kminusLV.Mt() - HPhysicsConstants::mass(12);
				minuskaon.y = kminusLV.Rapidity();
				minuskaon.ycm = kminusLV.Rapidity() - 0.74;
                                minuskaon.y0 = kminusLV.Rapidity()/0.74 - 1.;
				minuskaon.flow = flow;

				// eff.  corrections - dependence on track multiplicity vs phi-phi_RP
//				Int_t id = 12;   //TODO k+ corrections!! , 8 -> pi+ , 14 -> proton
				Int_t id = 11;   //TODO k+ corrections!! , 8 -> pi+ , 14 -> proton
				Float_t eff=0., weight=0.;

				eff=getEff(id,4-multClassNo,flow,kminusLV.Theta()*TMath::RadToDeg());
				if(eff!=0) weight=1./eff;
				if(weight>10||weight<0) weight=0.;
				minuskaon.occupancyCorr = weight;

//				eff=getAccRecEff(12,kminusLV.Pt(),kminusLV.Rapidity());
//				if(eff!=0) weight=1./eff;
//				if(weight>10||weight<0) weight=0.;
//				minuskaon.accreceffCorr = weight;

				minuskaon.MDCtracks = -10000.;
				minuskaon.METAclst = -10000.;

				minuskaon.downflag = downflag;



				//check that it does not come from phi meson decay (find a K+ companion)
                                minuskaon.phiMass = -1; //means that there was no phi meson mother particle
				if( (minuskaon.invMass>440 && minuskaon.invMass<550)
				   && (IC_kp.size() > 0) ){
#if isRealData
				    HParticleCand* kpCand;
#else
				    HParticleCandSim* kpCand;
#endif

				    TLorentzVector kplusLV(1.,1.,1.,1.);

				    for(unsigned int ikp=0; ikp<IC_kp.size(); ikp++){
					kpCand = HCategoryManager::getObject(kpCand,candCat,IC_kp[ikp]);
					HParticleTool::fillTLorentzVector(kplusLV,kpCand,11,true);

					Double_t kpGammaInv = ( kpCand->getBeta()>=1. ? 0 : TMath::Sqrt(1-TMath::Power(kpCand->getBeta(),2)) );
					Double_t kpMass = ( kpCand->getCorrectedMomentumPID(11)*kpGammaInv )/kpCand->getBeta();

					if(kpMass<440. || kpMass>550) continue;
					if(kplusLV.P()>1000) continue;

					TLorentzVector sumLV(1.,1.,1.,1.);
					sumLV = kminusLV + kplusLV;

                                        if( sumLV.M()>1005. && sumLV.M()<1030. ) minuskaon.phiMass = 1.;

				    }
				}

				if( minuskaon.invMass>300 && minuskaon.invMass<700
//				   && minuskaon.phiMass!=1.
				  ){
				    minuskaontree->Fill();

                                    downflag = -1.;
				}

			    }
			}//end of K- loop

IC_kp.clear();
IC_km.clear();

		    } // end cand cat
		    //-------------------------------------------------
		}
		//-------------------------------------------------
	    }
	}

	return 0;
    }

    Bool_t finalize()
    {   // this function is called once in the end of the
	// runtime of the program. Save histograms to file
	// etc. Tip: Change to your ROOT output file before
	// writing our histograms. You might not be the only
	// user of root files and the current directory could
        // the one from another user.

        delete effcorFile;

//	if(!out) out = new TFile(outfile,"UPDATE");
	if(out) {
	    out->cd();
	    out->mkdir("macros");
	    out->cd("macros");
	    TMacro *mAna = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/analysis.cc");
	    TMacro *mLoop = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/loopDST_task.C");
	    TMacro *mGlobVars = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/GlobVars.h");
	    TMacro *mPartReco = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/KChargedTrees.h");

	    mAna->Write();
	    mLoop->Write();
	    mGlobVars->Write();
            mPartReco->Write();

            out->cd();
            out->cd("kaonPlus");
	    pluskaontree->Write();

	    out->cd();
            out->cd("kaonMinus");
	    minuskaontree->Write();

            out->cd("eventPlanes");
	    for(int icc=0; icc<nCentralityClasses; icc++){
		heventplane[icc]->Write();
                hphiAB[icc]->Write();
	    }

            out->cd("histos");
	    hCentralityCounter->Write();
	    hEventCounter->Write();
	    hTrackCounter->Write();
            hZVerMult_noSelect->Write();
            hZVerMult_triggerSelect->Write();
            hZVerMult_allSelect->Write();
            hZVerMult_someSelect->Write();
	    hTofElossMom->Write();
	    hBetaMom_TOF->Write();
            hBetaMom_RPC->Write();
	    hMdcEloss_TOF->Write();
	    hMdcEloss_RPC->Write();
            hZMdcEloss->Write();
	    hMassMom_TOF->Write();
	    hMassMom_RPC->Write();
	    hTofElossBeta_neg->Write();

	    for(int ih=0; ih<8; ih++){
		hTofElossBeta_OneRod_all[ih]->Write();
		hTofElossBeta_OneRod_massCut[ih]->Write();
		hTofElossMom_OneRod_all[ih]->Write();
		hTofElossMom_OneRod_massCut[ih]->Write();
	    }

	    hMass_TOF->Write();
	    hMass_RPC->Write();
	    for(int ih=0; ih<4; ih++) hKaonMass_RPC[ih]->Write();
	    for(int ih=0; ih<5; ih++) hKaonMass_TOF[ih]->Write();

            hChisqRK_proton->Write();
            hChisqRK_pim->Write();
	    hChisqRK_pip->Write();
            hChisqRK_km->Write();
	    hChisqRK_kp->Write();
            hChisqMM_proton->Write();
            hChisqMM_pim->Write();
	    hChisqMM_pip->Write();
            hChisqMM_km->Write();
	    hChisqMM_kp->Write();

	    hMdcEloss_sim_TOF->Write();
	    hMdcEloss_sim_RPC->Write();
	    hMdcEloss_dataKcut_TOF->Write();
	    hMdcEloss_dataKcut_RPC->Write();


#if !isRealData
            out->cd("simulation");
	    gkplustree->Write();
	    gkminustree->Write();
#endif

//	    out->Save();
//	    out->Close();
	}

	return kTRUE;
    }
   //ClassDef(KChargedTrees,0) // no streamer, do not use if you include blank header in ACLIC
};

#endif
