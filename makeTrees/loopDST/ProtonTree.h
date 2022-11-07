#ifndef __PROTONTREE___     // good style to avaoid multiple includes
#define __PROTONTREE___

#include "mylibs.h"
#include "mystruc.h"

#include "GlobVars.h"

#include "UsefulFunctions.h"

using namespace std;

class ProtonTree : public HReconstructor
{
protected:
    // put all vars here which are not
    // local to a function

    Int_t PID_option=-1, iTrackCuts=-1, OccupCorr_pid=-1;

    HParticleTrackSorter sorter;

    CAND_KCHARGED proton;
    TTree *protontree;

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

    ProtonTree(const Text_t *name = "",const Text_t *title ="", Int_t PID_OPTION=-1, Int_t TRACK_CUTS=-1, Int_t OCCUP_CORR_PID=-1)
	: HReconstructor(name,title)
    {  // init your vars
	PID_option = PID_OPTION;
	iTrackCuts = TRACK_CUTS;
	OccupCorr_pid = OCCUP_CORR_PID;

    }

    virtual ~ProtonTree()
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
	    out->mkdir("proton");
	    out->cd("proton");

	    protontree = new TTree("protontree","K+ Tree");
	    protontree->Branch("proton",&proton,"eventPlane:phiAB:"
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

	    const char *effcorFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/MY_KpIII_effcor_flow_matrix_apr12_day108_gen8_4mbins.root";
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

/*	    const char *filenameAccCorrRecEff = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/AccCorr_and_RecoEff_fineBinning.root";
	    fileAccCorrRecEff = TFile::Open( filenameAccCorrRecEff );
            hAccCorr_cc0 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusAccCorr_centrality_0_10");
            hAccCorr_cc1 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusAccCorr_centrality_10_20");
            hAccCorr_cc2 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusAccCorr_centrality_20_30");
            hAccCorr_cc3 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusAccCorr_centrality_30_40");
            hRecoEff_cc0 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusRecoEff_centrality_0_10");
            hRecoEff_cc1 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusRecoEff_centrality_10_20");
            hRecoEff_cc2 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusRecoEff_centrality_20_30");
            hRecoEff_cc3 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusRecoEff_centrality_30_40");
*/
 
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

	const char *effcorFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/MY_KpIII_effcor_flow_matrix_apr12_day108_gen8_4mbins.root";
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

/*	const char *filenameAccCorrRecEff = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/AccCorr_and_RecoEff_fineBinning.root";
	fileAccCorrRecEff = TFile::Open( filenameAccCorrRecEff );
	hAccCorr_cc0 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusAccCorr_centrality_0_10");
	hAccCorr_cc1 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusAccCorr_centrality_10_20");
	hAccCorr_cc2 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusAccCorr_centrality_20_30");
	hAccCorr_cc3 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusAccCorr_centrality_30_40");
	hRecoEff_cc0 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusRecoEff_centrality_0_10");
	hRecoEff_cc1 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusRecoEff_centrality_10_20");
	hRecoEff_cc2 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusRecoEff_centrality_20_30");
	hRecoEff_cc3 = (TH2F*) fileAccCorrRecEff->Get("proton/hKplusRecoEff_centrality_30_40");
*/

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

            proton.eventVerX = vertex.getX();
            proton.eventVerY = vertex.getY();
	    proton.eventVerZ = vertex.getZ();
	    proton.evntVerChi2 = vertex.getChi2();

	    Float_t downflag = 1.;

            hEventCounter->Fill(cNumInputEv);

            if(!evtInfoCat) cout << "NO evtInfoCat !!!" << endl;
	    if(evtInfoCat)
	    {
		HParticleEvtInfo* evtInfo=0;
		evtInfo = HCategoryManager::getObject(evtInfo,evtInfoCat,0 );

		HWallEventPlane* event_plane;
		event_plane = HCategoryManager::getObject(event_plane,evtPlaneCat,0);

		Bool_t isGoodEvent = kFALSE;
#if isRealData || isEmbedding

		//--------------------- Event selection --------------------------------
		if(!evtInfo->isGoodEvent(Particle::kGoodTRIGGER) ) return kSkipEvent;
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


		    proton.eventPlane      = phiRP_tree;
		    proton.phiAB           = phiRPAB_tree;
		    proton.multiplicity    = mult;
		    proton.centralityclass = behruzCentralityClass;
                    proton.eventWeight     = event_weight;

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
       particleID = 14;
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

 
#endif

      //--------------------------------------------------------------------------------------------------
		    // loop over particle candidates in event
		    if(candCat){
			TString filename;
			//-------------------------------------------------
			// get the sector list infos for HLoop
			Int_t sectors [6];
			gLoop->getSectors(sectors); // fill sector array


			//------------------------------------------------------------------------
			// prepare track sorting
			// clean vectors and index arrays
			sorter.cleanUp();
			sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);
			Int_t nCandHad     = sorter.fill(HParticleTrackSorter::selectHadrons);
			Int_t nCandHadBest = sorter.selectBest(Particle::kIsBestRKSorter,Particle::kIsHadronSorter);

			if( nCandHad<0 || nCandHadBest<0 ) return 1;

			vector<Int_t> IC_prot;
#if isRealData
			vector<HParticleCand *> vprot;
			HParticleCand* cand=0;
  #if isExtendedOutput
			HTofCluster* tofClst=0;
  #endif
#else
			vector<HParticleCandSim *> vprot;
			HParticleCandSim* cand=0;
                        HGeantKine *kineRec=0;
#endif

			Int_t size = candCat->getEntries();
			for(Int_t j = 0; j < size; j++)
			{
			    cand = HCategoryManager::getObject(cand,candCat,j);

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


				Int_t   pid = cand->getPID();
				Bool_t goodPBeta=false, goodPEloss=false;

				switch(PID_option){
				case 4:
				    {
					//Second method of particle identification - check with HParticleTool functions
					Float_t deltaEloss = 0., sigmaEloss = -1.;
					if( HParticleTool::isParticledEdx(14,cand,deltaEloss,sigmaEloss)) {
					    if( TMath::Abs(deltaEloss/sigmaEloss) < 2. )
						goodPEloss = true;
					}

					Float_t deltaTime = 0., sigmaTime = -1.;
					if( HParticleTool::isParticleBeta(14,cand,2.,0.,2500.,deltaTime,sigmaTime,"apr12"))
					    goodPBeta = true; //2. is sigma allowed, 0. minMom, 1500. maxMom

					//ParticleTools cuts together
					if( goodPEloss && goodPBeta ){
					    IC_prot.push_back(j);
#if isRealData
					    vprot.push_back(new HParticleCand(*cand));
#else
					    vprot.push_back(new HParticleCandSim(*cand));
#endif
					}

				    }
				    break;
#if !isRealData
				case 5:
				    {
					if( cand->getGeantPID() == 11 ){
					    IC_prot.push_back(j);
					    vprot.push_back(new HParticleCandSim(*cand));
					}
					if( cand->getGeantPID() == 12 ){
					    IC_km.push_back(j);
					    vkm.push_back(new HParticleCandSim(*cand));
					}
				    }
                                    break;
#endif

				default:
				    cout << "This PID_option is not supported. Please select from 4." << endl;
                                    break;
				}

			    }
			} // end cand loop

                        cout << "N(proton) = " << IC_prot.size() << endl;

			//proton loop
			if( IC_prot.size() > 0 ){
#if isRealData
			    HParticleCand* pCand;
#else
			    HParticleCandSim* pCand;
#endif

			    TLorentzVector protonLV(1.,1.,1.,1.);

			    for(unsigned int ikp=0; ikp<IC_prot.size(); ikp++){
				pCand = HCategoryManager::getObject(pCand,candCat,IC_prot[ikp]);
				HParticleTool::fillTLorentzVector(protonLV,pCand,14,true);

				Double_t flow     = protonLV.Phi()*TMath::RadToDeg() - phiRP_tree;  // angle w.r.t. FW event plane
				if (flow<-180.) flow += 360.;
				if (flow>+180.) flow -= 360.;

				Double_t kpGammaInv = ( pCand->getBeta()>=1. ? 0 : TMath::Sqrt(1-TMath::Power(pCand->getBeta(),2)) );
				Double_t kpMass = ( pCand->getCorrectedMomentumPID(14)*kpGammaInv )/pCand->getBeta();

				proton.invMass = kpMass;
                                proton.mass = (pCand->getMass2()<0 ? -1 : TMath::Sqrt(pCand->getMass2()));
                                proton.system = pCand->getSystemUsed();
				proton.theta = protonLV.Theta()*TMath::RadToDeg();
				proton.phi = ( protonLV.Phi()*TMath::RadToDeg()>0 ? protonLV.Phi()*TMath::RadToDeg() : protonLV.Phi()*TMath::RadToDeg()+360.) ;
				proton.mom = protonLV.P();
				proton.momX = protonLV.Px();
				proton.momY = protonLV.Py();
				proton.momZ = protonLV.Pz();
				proton.pt = protonLV.Pt();
				proton.mt = protonLV.Mt();
				proton.mtm0 = protonLV.Mt() - HPhysicsConstants::mass(14);
				proton.y = protonLV.Rapidity();
				proton.ycm = protonLV.Rapidity() - 0.74;  //y_cm = y - 0.5* y_proj (y_proj = 1.48 from definition)
                                proton.y0 = protonLV.Rapidity()/0.74 - 1.; //y_0 = (y/y_proj)_cm = y_cm / y_proj,cm = (y - 0.5*y_proj)/(y_proj-0.5*y_proj) = (y - 0.5*y_proj)/(0.5*y_proj) = y/(0.5*y_proj) - 1
				proton.flow = flow;

				// eff.  corrections - dependence on track multiplicity vs phi-phi_RP
				Int_t id = OccupCorr_pid;
				Float_t eff=0., weight=0.;

				eff=getEff(id,4-multClassNo,flow,protonLV.Theta()*TMath::RadToDeg());
				if(eff!=0) weight=1./eff;
				if(weight>10||weight<0) weight=0.;
				proton.occupancyCorr = weight;

/*				weight=0.;
				eff=getAccCorr(multClassNo-1,protonLV.Pt(),protonLV.Rapidity()/0.74 - 1.);
				if(eff!=0.) weight=1./eff;
				if(weight>10000. || weight<1.) weight=0.;
				proton.accCorr = weight;

				weight=0.;
				eff=getRecoEff(multClassNo-1,protonLV.Pt(),protonLV.Rapidity()/0.74 - 1.);
				if(eff!=0.) weight=1./eff;
				if(weight>10000. || weight<1.) weight=0.;
				proton.recoEff = weight;
*/
				proton.MDCtracks = -10000.;
				proton.METAclst = -10000.;

                                proton.downflag = downflag;

				if( proton.invMass > 700 && proton.invMass < 1200 ){
				    protontree->Fill();

                                    downflag = -1.;
				}

			    }
			}//end of proton loop



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
	    TMacro *mPartReco = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/ProtonTree.h");

	    mAna->Write();
	    mLoop->Write();
	    mGlobVars->Write();
            mPartReco->Write();

            out->cd();
            out->cd("proton");
	    protontree->Write();


            out->cd("eventPlanes");
	    for(int icc=0; icc<nCentralityClasses; icc++){
		heventplane[icc]->Write();
                hphiAB[icc]->Write();
	    }


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
   //ClassDef(ProtonTree,0) // no streamer, do not use if you include blank header in ACLIC
};

#endif
