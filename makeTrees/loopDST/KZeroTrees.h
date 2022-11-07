#ifndef __KZEROTREES___     // good style to avaoid multiple includes
#define __KZEROTREES___

#include "mylibs.h"
#include "mystruc.h"

#include "GlobVars.h"

#include "UsefulFunctions.h"

#include "MEclassifiers.h"
#include "MVAK0s.h"

using namespace std;

class KZeroTrees : public HReconstructor
{
protected:
    // put all vars here which are not
    // local to a function

    static const Int_t particleID = 16;
//    static const Int_t PID_option = 5; // 4 - geantPID (only sim!!)
//    static const Int_t iPrecuts = 2; //0 - my original, 1 - my working, 2 - Simon's
    Int_t PID_option=-1, iPrecuts=-1;
    bool useMomCorr = false;

    HParticleTrackSorter sorter;
    //HParticleT0Reco t0Reco("apr12");

    CAND_K0S sekaon, mekaon;  //se = same event, me = mix event
    TTree *sekaontree, *mekaontree;
#if !isRealData
    GKINE_K0S gkzero;
    TTree *gkzerotree;
#endif
    
    //------------------------------------
    // file handling
    TFile* effcorFile;
    TFile* bananacutsFile;
    TFile *fileAccCorrRecEff;

    TCutG *betamom_2sig_pim_tof_pionCmom, *betamom_2sig_pim_rpc_pionCmom,
	  *betamom_2sig_pip_tof_pionCmom, *betamom_2sig_pip_rpc_pionCmom;

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


    Float_t TMVAd1, TMVAd2, TMVAd3, TMVAdVer, TMVAdMin, TMVAkaonMom;


public:

    KZeroTrees(const Text_t *name = "",const Text_t *title ="", Int_t PID_OPTION=-1, Int_t PRECUTS_INDEX=-1, Int_t MOM_CORR=-1)
	: HReconstructor(name,title)
    {  // init your vars
	PID_option = PID_OPTION;
	iPrecuts = PRECUTS_INDEX;
	if(MOM_CORR==1) useMomCorr = true;
        else useMomCorr = false;
    }

    virtual ~KZeroTrees()
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
	    out->mkdir("kZero");
	    out->cd("kZero");

	    sekaontree = new TTree("sekaontree","Same Event Kaon Tree");
	    sekaontree->Branch("sekaon",&sekaon,"eventPlane:phiAB:"
			       "eventVerX:eventVerY:eventVerZ:evntVerChi2:decayVerX:decayVerY:decayVerZ:"
			       "dist1:dist2:dist3:distVer:distMin:oAngle:"
			       "theta:phi:momX:momY:momZ:mom:"
			       "multiplicity:centralityclass:"
			       "invMass:mt:mtm0:pt:y:ycm:y0:"
			       "flow:"
			       "occupancyCorr:accCorr:recoEff:eventWeight:"
			       "mvaResult:"
			       "downflag");

	    mekaontree = new TTree("mekaontree","Mix Event Kaon Tree");
	    mekaontree->Branch("mekaon",&mekaon,"eventPlane:phiAB:"
			       "eventVerX:eventVerY:eventVerZ:evntVerChi2:decayVerX:decayVerY:decayVerZ:"
			       "dist1:dist2:dist3:distVer:distMin:oAngle:"
			       "theta:phi:momX:momY:momZ:mom:"
			       "multiplicity:centralityclass:"
			       "invMass:mt:mtm0:pt:y:ycm:y0:"
			       "flow:"
			       "occupancyCorr:accCorr:recoEff:eventWeight:"
			       "mvaResult:"
			       "downflag");

#if !isRealData
	    out->cd();
	    out->mkdir("simulation");
	    out->cd("simulation");
	    gkzerotree = new TTree("gkzerotree","GeantKine Kaon Tree");
	    gkzerotree->Branch("gkzero",&gkzero,"eventPlane:eventPlane_gkine:phiAB:"
			       "theta:phi:momX:momY:momZ:mom:"
			       "multiplicity:centralityclass:impactparam:invMass:mt:pt:y:ycm:y0:flow:flow_gkine:weight:"
			       "dist1:dist2:dist3:distVer:distMin:dist1Org:dist2Org:dist3Org:distVerOrg:distMinOrg:mvaResult:"
			       "pipAcc/B:pimAcc:pipDecay:pimDecay:pipRec:pimRec:pichargedDecay:pizeroDecay");
#endif

	    const char *effcorFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/MY_KpIII_effcor_flow_matrix_apr12_day108_gen8_4mbins.root";
	    effcorFile = TFile::Open( effcorFileName );
	    effPiM_cc3 = (TH2F *)effcorFile->Get("heff4_id9");
	    effPiP_cc3 = (TH2F *)effcorFile->Get("heff4_id8");
	    effPiM_cc2 = (TH2F *)effcorFile->Get("heff3_id9");
	    effPiP_cc2 = (TH2F *)effcorFile->Get("heff3_id8");
	    effPiM_cc1 = (TH2F *)effcorFile->Get("heff2_id9");
	    effPiP_cc1 = (TH2F *)effcorFile->Get("heff2_id8");
	    effPiM_cc0 = (TH2F *)effcorFile->Get("heff1_id9");
	    effPiP_cc0 = (TH2F *)effcorFile->Get("heff1_id8");

	    const char *filenameAccCorrRecEff = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/AccCorr_and_RecoEff_fineBinning.root";
	    fileAccCorrRecEff = TFile::Open( filenameAccCorrRecEff );
            hAccCorr_cc0 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroAccCorr_centrality_0_10");
            hAccCorr_cc1 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroAccCorr_centrality_10_20");
            hAccCorr_cc2 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroAccCorr_centrality_20_30");
            hAccCorr_cc3 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroAccCorr_centrality_30_40");
            hRecoEff_cc0 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroRecoEff_centrality_0_10");
            hRecoEff_cc1 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroRecoEff_centrality_10_20");
            hRecoEff_cc2 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroRecoEff_centrality_20_30");
            hRecoEff_cc3 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroRecoEff_centrality_30_40");

	    const char *bananacutsFileName = "/lustre/hebe/hades/user/prozorov/mar19/Mar19AgAg1580_Gen5.root";
	    bananacutsFile = TFile::Open( bananacutsFileName );
		betamom_2sig_pim_tof_pionCmom = (TCutG *) bananacutsFile->Get("tcgPBetaPiMTOF2Sig");
		betamom_2sig_pim_rpc_pionCmom = (TCutG *) bananacutsFile->Get("tcgPBetaPiMRPC2Sig");
		betamom_2sig_pip_tof_pionCmom = (TCutG *) bananacutsFile->Get("tcgPBetaPiPTOF2Sig");
		betamom_2sig_pip_rpc_pionCmom = (TCutG *) bananacutsFile->Get("tcgPBetaPiPRPC2Sig");

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

	    hEventCounter = new TH1F( "hEventCounter", "Counter of events", cNumEventCounters+1, 0, cNumEventCounters+1);
//	        hEventCounter->GetXaxis()->SetBinLabel( 1, "InputEv" );
//		hEventCounter->GetXaxis()->SetBinLabel( 2, "PT3" );
//		hEventCounter->GetXaxis()->SetBinLabel( 3, "selectStart" );
//		hEventCounter->GetXaxis()->SetBinLabel( 4, "StartPileUp" );
//		hEventCounter->GetXaxis()->SetBinLabel( 5, "GoodClusterVertex" );
//		hEventCounter->GetXaxis()->SetBinLabel( 6, "GoodCandVertex" );
//		hEventCounter->GetXaxis()->SetBinLabel( 7, "NoVeto" );
//		hEventCounter->GetXaxis()->SetBinLabel( 8, "VetoStart" );
//		hEventCounter->GetXaxis()->SetBinLabel( 9, "StartMeta" );
//		hEventCounter->GetXaxis()->SetBinLabel( 10, "GoodCentrality" );
//		hEventCounter->GetXaxis()->SetBinLabel( 11,"MixedEventSelector" );

            hBetaMom_TOF_massCut = new TH2D("hBetaMom_TOF_massCut","beta vs mom -> TOF with mass<300",150,-1000,2000,120,0,1.2);
            hBetaMom_RPC_massCut = new TH2D("hBetaMom_RPC_massCut","beta vs mom -> RPC with mass<300",150,-1000,2000,120,0,1.2);
            hMassMom_TOF_betamomCut = new TH2D("hMassMom_TOF_betamomCut","mass vs mom -> TOF with 2 sigma banana cut",300,-1000,2000,100,0,2000);
            hMassMom_RPC_betamomCut = new TH2D("hMassMom_RPC_betamomCut","mass vs mom -> RPC with 2 sigma banana cut",300,-1000,2000,100,0,2000);

	    for(int ih=0; ih<4; ih++){
		hPionMass_RPC[ih] = new TH1D( Form("hPionMass_RPC_%i",ih), "mass zoom -> RPC with conclusive cuts", 160,-400,400);
		hPionMass_TOF[ih] = new TH1D( Form("hPionMass_TOF_%i",ih), "mass zoom -> TOF with conclusive cuts", 160,-400,400);
	    }

            hPiPairInvMass_all       = new TH1D("hPiPairInvMass_all","k0s cand inv mass",80,300,700);
            hPiPairInvMass_minTC     = new TH1D("hPiPairInvMass_minTC","k0s cand inv mass",80,300,700);
            hPiPairInvMass_mva50     = new TH1D("hPiPairInvMass_mva50","k0s cand inv mass",80,300,700);
            hPiPairInvMass_mva85     = new TH1D("hPiPairInvMass_mva85","k0s cand inv mass",80,300,700);
            hPiPairInvMass_mva95     = new TH1D("hPiPairInvMass_mva95","k0s cand inv mass",80,300,700);
            hPiPairInvMass_broadTC   = new TH1D("hPiPairInvMass_broadTC","k0s cand inv mass",80,300,700);
            hPiPairInvMass_StoBratio = new TH1D("hPiPairInvMass_StoBratio","k0s cand inv mass",80,300,700);
            hPiPairInvMass_Signif    = new TH1D("hPiPairInvMass_Signif","k0s cand inv mass",80,300,700);
	}

	InitMVAK0s(&TMVAd1, &TMVAd2, &TMVAd3, &TMVAdVer, &TMVAdMin, &TMVAkaonMom, iPrecuts);

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
	effPiM_cc2 = (TH2F *)effcorFile->Get("heff3_id9");
	effPiP_cc2 = (TH2F *)effcorFile->Get("heff3_id8");
	effPiM_cc1 = (TH2F *)effcorFile->Get("heff2_id9");
	effPiP_cc1 = (TH2F *)effcorFile->Get("heff2_id8");
	effPiM_cc0 = (TH2F *)effcorFile->Get("heff1_id9");
	effPiP_cc0 = (TH2F *)effcorFile->Get("heff1_id8");

	const char *filenameAccCorrRecEff = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/AccCorr_and_RecoEff_fineBinning.root";
	fileAccCorrRecEff = TFile::Open( filenameAccCorrRecEff );
	hAccCorr_cc0 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroAccCorr_centrality_0_10");
	hAccCorr_cc1 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroAccCorr_centrality_10_20");
	hAccCorr_cc2 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroAccCorr_centrality_20_30");
	hAccCorr_cc3 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroAccCorr_centrality_30_40");
	hRecoEff_cc0 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroRecoEff_centrality_0_10");
	hRecoEff_cc1 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroRecoEff_centrality_10_20");
	hRecoEff_cc2 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroRecoEff_centrality_20_30");
	hRecoEff_cc3 = (TH2F*) fileAccCorrRecEff->Get("kaonZero/hKzeroRecoEff_centrality_30_40");

	const char *bananacutsFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root";
	bananacutsFile = TFile::Open( bananacutsFileName );
	betamom_2sig_pim_tof_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutPiM_TOF_2.0");
	betamom_2sig_pim_rpc_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutPiM_RPC_2.0");
	betamom_2sig_pip_tof_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutPiP_TOF_2.0");
	betamom_2sig_pip_rpc_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutPiP_RPC_2.0");


	InitMVAK0s(&TMVAd1, &TMVAd2, &TMVAd3, &TMVAdVer, &TMVAdMin, &TMVAkaonMom, iPrecuts);

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

//	    cout << "VERTEX: x=" << vertex.getX() << "   y=" << vertex.getY() << "   z=" << vertex.getZ() << "    chi2=" << vertex.getChi2() << endl;

            sekaon.eventVerX = vertex.getX();
            sekaon.eventVerY = vertex.getY();
	    sekaon.eventVerZ = vertex.getZ();
            sekaon.evntVerChi2 = vertex.getChi2();

	    mekaon.evntVerChi2 = vertex.getChi2(); //mix event tree - has special event vertex, so i take only Chi2

	    Float_t downflag_se = 1.;
	    Float_t downflag_me = 1.;

	    hEventCounter->Fill(cNumInputEv);

#if !isRealData
	    HGeantKine*       kine;
	    HLinearCategory* kineCat = (HLinearCategory*)HCategoryManager::getCategory(catGeantKine);
	    HGeantHeader* fGeantHeader = NULL;
#endif
            if(!evtInfoCat) cout << "NO evtInfoCat !!!" << endl;
	    if(evtInfoCat)
	    {
		HParticleEvtInfo* evtInfo=0;
		evtInfo = HCategoryManager::getObject(evtInfo,evtInfoCat,0 );

//                cout << "N_TOFRPC= " << evtInfo->getSumTofMultCut()+evtInfo->getSumRpcMultHitCut() << endl;

		eventmixer.setUseLeptons(kFALSE);
		eventmixer.setPIDs(8,9,particleID);  // pi-,pi+,k0s  which PIDs and MotherID are stored
		eventmixer.setBuffSize(20);   // size of buffer has to chosen according to stat   //TODO how big is needed??
                eventmixer.setWaitForBuffer(kFALSE);
		eventmixer.setEventClassifier(eventClassifierMultTargetTime);
//		eventmixer.setEventClassifier(eventClassifierMultRPhi);

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

		if(! MixedEventSelectionCuts(header,evtInfo) ) return kSkipEvent;
		hEventCounter->Fill(cNumEventCounters); //fill one more than standard

		isGoodEvent = kTRUE;

//		if(evtInfo && evtInfo->isGoodEvent(Particle::kGoodTRIGGER|
//						   Particle::kGoodVertexClust|
//						   Particle::kGoodVertexCand|
//						   Particle::kGoodSTART|
//						   Particle::kNoPileUpSTART|
//						   Particle::kNoVETO|
//						   Particle::kGoodSTARTVETO|
//       					   Particle::kGoodSTARTMETA
//						  )
//                   && MixedEventSelectionCuts(header,evtInfo) //new condition - from Simon!!
//		  ) isGoodEvent = kTRUE;
#else
		if(evtInfo && MixedEventSelectionCuts(header,evtInfo)) isGoodEvent = kTRUE;
//		if(evtInfo) isGoodEvent = kTRUE;
#endif
                if(!isGoodEvent) cout << "NO GoodEvent = NO evtInfo !!!" << endl;
                if(isGoodEvent)
		{
		    Double_t mult = evtInfo->getSumTofMultCut()+evtInfo->getSumRpcMultHitCut();

#if !isRealData
		    if(!(fGeantHeader = gLoop->getGeantHeader(kFALSE))){ cout << "NO HGEANT!!! " << endl; return 0; }

		    Double_t geantEventPlaneAngle = fGeantHeader->getEventPlane();
                    Double_t geantImpactParameter = fGeantHeader->getImpactParameter();
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

		    phiRP = evtChara.getEventPlane(HParticleEvtChara::kDefault)*TMath::RadToDeg();
		    phiRP_A = evtChara.getEventPlane(HParticleEvtChara::kDefault,1)*TMath::RadToDeg();
		    phiRP_B = evtChara.getEventPlane(HParticleEvtChara::kDefault,2)*TMath::RadToDeg();
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

		    if( (phiRP == -1*TMath::RadToDeg()) || (phiRP_A == -1*TMath::RadToDeg()) || (phiRP_B == -1*TMath::RadToDeg()) ) return 1;
#endif

		    //fill counter with appropriate weight
		    hCentralityCounter->Fill(behruzCentralityClass,event_weight);


		    sekaon.eventPlane      = phiRP_tree;
		    sekaon.phiAB           = phiRPAB_tree;
		    sekaon.multiplicity    = mult;
		    sekaon.centralityclass = behruzCentralityClass;
                    sekaon.eventWeight     = event_weight;

		    mekaon.eventPlane      = phiRP_tree;
		    mekaon.phiAB           = phiRPAB_tree;
		    mekaon.multiplicity    = mult;
		    mekaon.centralityclass = behruzCentralityClass;
                    mekaon.eventWeight     = event_weight;

#if !isRealData
		    gkzero.eventPlane = phiRP_tree;
		    gkzero.eventPlane_gkine = geantEventPlaneAngle;
		    gkzero.phiAB = phiRPAB_tree;
		    gkzero.multiplicity = mult;
		    gkzero.centralityclass = behruzCentralityClass;
		    gkzero.impactparam = geantImpactParameter;
#endif



#if !isRealData
       //---------------------------------- Kine --------------------------------------------------------
       Bool_t goodKZeroFound = kFALSE;
       for(Int_t j = 0; j < kineCat->getEntries(); j ++) {
	   kine = HCategoryManager::getObject(kine,kineCat,j);
	   if( kine->isPrimary() && kine->getID() == particleID && kine->getGeneratorInfo() == (500+particleID) ){
               goodKZeroFound = kTRUE;

               gkzero.invMass = HPhysicsConstants::mass(particleID);
	       gkzero.pt   = kine->getTransverseMomentum();
	       gkzero.mt   = sqrt(gkzero.pt*gkzero.pt+gkzero.invMass*gkzero.invMass);
	       gkzero.y    = kine->getRapidity();
               gkzero.ycm  = kine->getRapidity() - 0.822;
               gkzero.y0   = kine->getRapidity()/0.822 - 1.;

               kine->getMomentum(gkzero.momX,gkzero.momY,gkzero.momZ);
	       gkzero.mom = kine->getTotalMomentum();

	       Int_t sec_kine;
	       Float_t phi_kine, theta_kine;
	       sec_kine = kine->getPhiThetaDeg(theta_kine,phi_kine);

	       gkzero.theta = theta_kine;
               if( phi_kine > 180. ) phi_kine -= 360.;
	       gkzero.phi = phi_kine;

	       gkzero.flow = phi_kine - phiRP_tree;
               gkzero.flow_gkine = phi_kine - geantEventPlaneAngle;
               if(gkzero.flow < -180.) gkzero.flow += 360.;
               if(gkzero.flow > +180.) gkzero.flow -= 360.;
               if(gkzero.flow_gkine < -180.) gkzero.flow_gkine += 360.;
               if(gkzero.flow_gkine > +180.) gkzero.flow_gkine -= 360.;

	       gkzero.weight = 1.;

	       vector<HGeantKine *> kaonDaughters;
	       Int_t numberOfDaughters = kine->getDaughters(kine,kaonDaughters);

	       HGeantKine *pipGK = 0;
	       HGeantKine *pimGK = 0;
	       Int_t counterPiZero = 0;

	       Bool_t pipAcc = kFALSE;
	       Bool_t pimAcc = kFALSE;
	       Bool_t pipRec = kFALSE;
	       Bool_t pimRec = kFALSE;

	       for(int id=0; id<numberOfDaughters; id++){
                   if( kaonDaughters[id]->getID() == 8 ) pipGK = kaonDaughters[id];
                   else if( kaonDaughters[id]->getID() == 9 ) pimGK = kaonDaughters[id];
		   else  counterPiZero++;
	       }


	       //if( counterPiZero == 2 ){}
	       //TODO: probehl rozpad K0s -> pi0 + pi0 (nemusi byt vzdy 2 z nejakeho divneho duvodu.. proto zkontrolovat ze mame 2 dcerine castice!!!)
	       if( counterPiZero > 0 ) gkzero.pizeroDecay = kTRUE;
               else gkzero.pizeroDecay = kFALSE;

	       if( (pipGK == 0) || (pimGK == 0) ) gkzero.pichargedDecay = kFALSE;
	       else{
                   gkzero.pichargedDecay = kTRUE;

                   HGeantKine *dPip = 0; dPip = pipGK->getChargedDecayDaughter(pipGK);
		   HGeantKine *dPim = 0; dPim = pimGK->getChargedDecayDaughter(pimGK);


		   //v134
//		   if( dPip != 0 ) gkzero.pipDecay = kTRUE;
//		   else gkzero.pipDecay = kFALSE;
//		   if( dPim != 0 ) gkzero.pimDecay = kTRUE;
//                   else gkzero.pimDecay = kFALSE;
//                   if( pipGK->isInAcceptanceBit(4,4,4,4,1,1) && !pipGK->isAtAnyMdcEdge(2) ) pipAcc = kTRUE;
//                   if( pimGK->isInAcceptanceBit(4,4,4,4,1,1) && !pimGK->isAtAnyMdcEdge(2) ) pimAcc = kTRUE;
//		   if( dPip != 0 ) gkzero.pipDecay = kTRUE;
//		   else gkzero.pipDecay = kFALSE;
//		   if( dPim != 0 ) gkzero.pimDecay = kTRUE;
//                   else gkzero.pimDecay = kFALSE;

		   //v135
		   if( dPip != 0 ){
		       gkzero.pipDecay = kTRUE;
                       if( (pipGK->isInAcceptanceBit(4,4,4,4,1,1) && !pipGK->isAtAnyMdcEdge(2))
			  || (dPip->isInAcceptanceBit(4,4,4,4,1,1) && !dPip->isAtAnyMdcEdge(2)) )
			   pipAcc = kTRUE;
		   }
		   else{
		       gkzero.pipDecay = kFALSE;
		       if( pipGK->isInAcceptanceBit(4,4,4,4,1,1) && !pipGK->isAtAnyMdcEdge(2) )
			   pipAcc = kTRUE;
		   }

		   if( dPim != 0 ){
		       gkzero.pimDecay = kTRUE;
		       if( (pimGK->isInAcceptanceBit(4,4,4,4,1,1) && !pimGK->isAtAnyMdcEdge(2))
			  || (dPim->isInAcceptanceBit(4,4,4,4,1,1) && !dPim->isAtAnyMdcEdge(2)) )
			   pimAcc = kTRUE;
		   }
		   else{
		       gkzero.pimDecay = kFALSE;
		       if( pimGK->isInAcceptanceBit(4,4,4,4,1,1) && !pimGK->isAtAnyMdcEdge(2) )
			   pimAcc = kTRUE;
		   }

                   //calculate d_x from event generator
		   HGeomVector pipBaseOrg, pipDirOrg, pimBaseOrg, pimDirOrg, decayVertexOrg, eventVertexOrg, motherDirOrg, pipDVOrg, pimDVOrg;
		   Double_t distVertOrg, dist1Org, dist2Org, dist3Org, distMinOrg;
/*
                   //firts approach
                   kine->getVertex(eventVertexOrg);
		   pipGK->getVertex(pipBaseOrg);
		   pimGK->getVertex(pimBaseOrg);

		   kine->getMomentum(motherDirOrg);
                   pipGK->getMomentum(pipDirOrg);
		   pimGK->getMomentum(pimDirOrg);

		   if(pipBaseOrg==pimBaseOrg) decayVertexOrg = pipBaseOrg;
                   else decayVertexOrg = calcVertexAnalytical(pipBaseOrg, pipDirOrg, pimBaseOrg, pimDirOrg);


		   distVertOrg = (decayVertexOrg-eventVertexOrg).length();
		   dist1Org = calculateMinimumDistanceStraightToPoint(pimBaseOrg,     pimDirOrg,    eventVertexOrg);
		   dist2Org = calculateMinimumDistanceStraightToPoint(pipBaseOrg,     pipDirOrg,    eventVertexOrg);
		   dist3Org = calculateMinimumDistanceStraightToPoint(decayVertexOrg, motherDirOrg, eventVertexOrg);
		   distMinOrg = calculateMinimumDistance(pimBaseOrg, pimDirOrg, pipBaseOrg, pipDirOrg);



                   //second approach
		   kine->getVertex(eventVertexOrg);
		   pipGK->getVertex(pipDVOrg);
		   pimGK->getVertex(pimDVOrg);

		   CalcSegVector(pipDVOrg.Z(), TMath::Sqrt(pipDVOrg.X()*pipDVOrg.X()+pipDVOrg.Y()*pipDVOrg.Y()), pipGK->getPhiDeg(kTRUE)*TMath::DegToRad(), pipGK->getThetaDeg()*TMath::DegToRad(), pipBaseOrg, pipDirOrg);
		   CalcSegVector(pimDVOrg.Z(), TMath::Sqrt(pimDVOrg.X()*pimDVOrg.X()+pimDVOrg.Y()*pimDVOrg.Y()), pimGK->getPhiDeg(kTRUE)*TMath::DegToRad(), pimGK->getThetaDeg()*TMath::DegToRad(), pimBaseOrg, pimDirOrg);

		   kine->getMomentum(motherDirOrg);

		   if(pipDVOrg==pimDVOrg) decayVertexOrg = pipDVOrg;
                   else decayVertexOrg = calcVertexAnalytical(pipBaseOrg, pipDirOrg, pimBaseOrg, pimDirOrg);

		   distVertOrg = (decayVertexOrg-eventVertexOrg).length();
		   dist1Org = calculateMinimumDistanceStraightToPoint(pimBaseOrg,     pimDirOrg,    eventVertexOrg);
		   dist2Org = calculateMinimumDistanceStraightToPoint(pipBaseOrg,     pipDirOrg,    eventVertexOrg);
		   dist3Org = calculateMinimumDistanceStraightToPoint(decayVertexOrg, motherDirOrg, eventVertexOrg);
		   distMinOrg = calculateMinimumDistance(pimBaseOrg, pimDirOrg, pipBaseOrg, pipDirOrg);

*/

		   //third approach
		   kine->getVertex(eventVertexOrg);
		   pipGK->getVertex(pipDVOrg);
		   pimGK->getVertex(pimDVOrg);

		   HGeomVector axisOrigin(0.,0.,0.);
		   HGeomVector zAxisDir(0.,0.,1.);
		   pipGK->getVertex(pipBaseOrg);
		   pipGK->getMomentum(pipDirOrg);
                   HGeomVector pipClosestPointToZ = calculatePointOfClosestApproach(axisOrigin,zAxisDir,pipBaseOrg,pipDirOrg);
		   pimGK->getVertex(pimBaseOrg);
		   pimGK->getMomentum(pimDirOrg);
                   HGeomVector pimClosestPointToZ = calculatePointOfClosestApproach(axisOrigin,zAxisDir,pimBaseOrg,pimDirOrg);

		   CalcSegVector(pipClosestPointToZ.Z(), TMath::Sqrt(pipClosestPointToZ.X()*pipClosestPointToZ.X()+pipClosestPointToZ.Y()*pipClosestPointToZ.Y()), pipGK->getPhiDeg(kTRUE)*TMath::DegToRad(), pipGK->getThetaDeg()*TMath::DegToRad(), pipBaseOrg, pipDirOrg);
		   CalcSegVector(pimClosestPointToZ.Z(), TMath::Sqrt(pimClosestPointToZ.X()*pimClosestPointToZ.X()+pimClosestPointToZ.Y()*pimClosestPointToZ.Y()), pimGK->getPhiDeg(kTRUE)*TMath::DegToRad(), pimGK->getThetaDeg()*TMath::DegToRad(), pimBaseOrg, pimDirOrg);

		   kine->getMomentum(motherDirOrg);

		   if(pipDVOrg==pimDVOrg) decayVertexOrg = pipDVOrg;
                   else decayVertexOrg = calcVertexAnalytical(pipBaseOrg, pipDirOrg, pimBaseOrg, pimDirOrg);

		   distVertOrg = (decayVertexOrg-eventVertexOrg).length();
		   dist1Org = calculateMinimumDistanceStraightToPoint(pimBaseOrg,     pimDirOrg,    eventVertexOrg);
		   dist2Org = calculateMinimumDistanceStraightToPoint(pipBaseOrg,     pipDirOrg,    eventVertexOrg);
		   dist3Org = calculateMinimumDistanceStraightToPoint(decayVertexOrg, motherDirOrg, eventVertexOrg);
		   distMinOrg = calculateMinimumDistance(pimBaseOrg, pimDirOrg, pipBaseOrg, pipDirOrg);


                   //save the data
		   gkzero.dist1Org = dist1Org;
		   gkzero.dist2Org = dist2Org;
		   gkzero.dist3Org = dist3Org;
		   gkzero.distMinOrg = distMinOrg;
		   gkzero.distVerOrg = distVertOrg;



		   HParticleCandSim* cand=0;
		   HParticleCandSim* pipCand=0;
		   HParticleCandSim* pimCand=0;

		   if(candCat){
		       Int_t size = candCat->getEntries();

		       for(Int_t j = 0; j < size; j++)
		       {
			   cand = HCategoryManager::getObject(cand,candCat,j);

			   Bool_t test=kFALSE;
			   if(
//			      cand->isFlagAND(4,Particle::kIsAcceptedHitInnerMDC,Particle::kIsAcceptedHitOuterMDC,Particle::kIsAcceptedHitMETA,Particle::kIsAcceptedRK) &&
			      cand->getInnerSegmentChi2() > 0 &&
			      cand->getChi2()             < 1000
			     ) test=kTRUE;
			   if(!test) continue;

			   if(cand->isAtAnyMdcEdge()) continue;
//			   if(cand->getSystemUsed() == -1) continue;

			   if(cand->getGeantParentTrackNum()==kine->getTrack() || cand->getGeantGrandParentTrackNum()==kine->getTrack()){
			       if(cand->getCharge()>0){
                                   pipRec = kTRUE;
				   pipCand = HCategoryManager::getObject(pipCand,candCat,j);
			       }
			       if(cand->getCharge()<0){
                                   pimRec = kTRUE;
				   pimCand = HCategoryManager::getObject(pimCand,candCat,j);
			       }
			    }

//			   if( cand->getGeantTrack() == pipGK->getTrack() ){
//			       pipRec = kTRUE;
//			       pipCand = HCategoryManager::getObject(pipCand,candCat,j);
//			   }
//			   if( cand->getGeantTrack() == pimGK->getTrack() ){
//			       pimRec = kTRUE;
//			       pimCand = HCategoryManager::getObject(pimCand,candCat,j);
//			   }
		       }

		   }

		   if(pipRec && pimRec){
		       HGeomVector pipBase, pipDir, pimBase, pimDir, decayVertex, eventVertex, motherDir;

                       kine->getVertex(eventVertex);
//		       eventVertex.setXYZ(vertex.getX(),vertex.getY(),vertex.getZ());

		       TLorentzVector pipLV(1.,1.,1.,1.), pimLV(1.,1.,1.,1.), sumLV(1.,1.,1.,1.);
		       Double_t distVert, dist1, dist2, dist3, distMin;

		       pipLV.SetXYZM(TMath::Abs(pipCand->getMomentum()) * TMath::Sin( TMath::DegToRad() * pipCand->getTheta() ) * TMath::Cos( TMath::DegToRad() * pipCand->getPhi() ),
				     TMath::Abs(pipCand->getMomentum()) * TMath::Sin( TMath::DegToRad() * pipCand->getTheta() ) * TMath::Sin( TMath::DegToRad() * pipCand->getPhi() ),
				     TMath::Abs(pipCand->getMomentum()) * TMath::Cos( TMath::DegToRad() * pipCand->getTheta() ),
				     HPhysicsConstants::mass(8) );
		       CalcSegVector(pipCand->getZ(), pipCand->getR(), pipCand->getPhi()*TMath::DegToRad(), pipCand->getTheta()*TMath::DegToRad(), pipBase, pipDir);

		       pimLV.SetXYZM(TMath::Abs(pimCand->getMomentum()) * TMath::Sin( TMath::DegToRad() * pimCand->getTheta() ) * TMath::Cos( TMath::DegToRad() * pimCand->getPhi() ),
				     TMath::Abs(pimCand->getMomentum()) * TMath::Sin( TMath::DegToRad() * pimCand->getTheta() ) * TMath::Sin( TMath::DegToRad() * pimCand->getPhi() ),
				     TMath::Abs(pimCand->getMomentum()) * TMath::Cos( TMath::DegToRad() * pimCand->getTheta() ),
				     HPhysicsConstants::mass(9) );
		       CalcSegVector(pimCand->getZ(), pimCand->getR(), pimCand->getPhi()*TMath::DegToRad(), pimCand->getTheta()*TMath::DegToRad(), pimBase, pimDir);

		       //reconstructed LV from one pair of pions
		       sumLV = pipLV + pimLV;

		       //topology variables from particle cands
		       motherDir.setXYZ(sumLV.X(),sumLV.Y(),sumLV.Z());
#if isEmbedding
		       decayVertex = calcVertexAnalytical(pipBase, pipDir, pimBase, pimDir);
#else
                       HGeomVector pipDV, pimDV;
		       pipGK->getVertex(pipDV);
		       pimGK->getVertex(pimDV);
		       if(pipDV==pimDV) decayVertex = pipDV;
                       else decayVertex = calcVertexAnalytical(pipBase, pipDir, pimBase, pimDir);
#endif
		       distVert = TMath::Sqrt( (decayVertex.getX()-eventVertex.getX())*(decayVertex.getX()-eventVertex.getX()) +
					      (decayVertex.getY()-eventVertex.getY())*(decayVertex.getY()-eventVertex.getY()) +
					      (decayVertex.getZ()-eventVertex.getZ())*(decayVertex.getZ()-eventVertex.getZ()) );
		       dist1 = calculateMinimumDistanceStraightToPoint(pimBase,     pimDir,    eventVertex);
		       dist2 = calculateMinimumDistanceStraightToPoint(pipBase,     pipDir,    eventVertex);
		       dist3 = calculateMinimumDistanceStraightToPoint(decayVertex, motherDir, eventVertex);
		       distMin = calculateMinimumDistance(pimBase, pimDir, pipBase, pipDir);

		       TMVAd1 = dist1;
		       TMVAd2 = dist2;
		       TMVAd3 = dist3;
		       TMVAdVer = distVert;
		       TMVAdMin = distMin;
		       TMVAkaonMom = sumLV.P();

		       gkzero.dist1 = dist1;
		       gkzero.dist2 = dist2;
		       gkzero.dist3 = dist3;
		       gkzero.distMin = distMin;
		       gkzero.distVer = distVert;
		       gkzero.mvaResult = EvalMVAK0s();
		   }
		   else{
		       gkzero.dist1 = -10000.;
		       gkzero.dist2 = -10000.;
		       gkzero.dist3 = -10000.;
		       gkzero.distMin = -10000.;
		       gkzero.distVer = -10000.;
		       gkzero.mvaResult = -10000.;
		   }

	       }

	       gkzero.pipAcc = pipAcc;
	       gkzero.pimAcc = pimAcc;
	       gkzero.pipRec = pipRec;
	       gkzero.pimRec = pimRec;
	   }
       }

       if(goodKZeroFound) gkzerotree->Fill();
       else{
	   gkzero.eventPlane = -10000.;
	   gkzero.eventPlane_gkine = -10000.;
	   gkzero.phiAB = -10000.;
	   gkzero.theta = -10000.;
	   gkzero.phi = -10000.;
	   gkzero.momX = -10000.;
	   gkzero.momY = -10000.;
	   gkzero.momZ = -10000.;
	   gkzero.mom = -10000.;
	   gkzero.multiplicity = -10000.;
	   gkzero.centralityclass = -10000.;
	   gkzero.impactparam = -10000.;
	   gkzero.invMass = -10000.;
	   gkzero.mt = -10000.;
	   gkzero.pt = -10000.;
	   gkzero.y = -10000.;
	   gkzero.ycm = -10000.;
	   gkzero.y0 = -10000.;
	   gkzero.flow = -10000.;
	   gkzero.flow_gkine = -10000.;
	   gkzero.weight = -10000.;
	   gkzero.dist1 = -10000.;
	   gkzero.dist2 = -10000.;
	   gkzero.dist3 = -10000.;
	   gkzero.distVer = -10000.;
	   gkzero.distMin = -10000.;
	   gkzero.dist1Org = -10000.;
	   gkzero.dist2Org = -10000.;
	   gkzero.dist3Org = -10000.;
	   gkzero.distVerOrg = -10000.;
	   gkzero.distMinOrg = -10000.;
	   gkzero.mvaResult = -10000.;
	   gkzero.pipAcc = kFALSE;
	   gkzero.pimAcc = kFALSE;
	   gkzero.pipRec = kFALSE;
	   gkzero.pimRec = kFALSE;
	   gkzero.pipDecay = kFALSE;
	   gkzero.pimDecay = kFALSE;
	   gkzero.pichargedDecay = kFALSE;
	   gkzero.pizeroDecay = kFALSE;

           gkzerotree->Fill();
       }
#endif

      //--------------------------------------------------------------------------------------------------
		    // loop over particle candidates in event
		    if(candCat){
			TString filename;
			Int_t sectors [6];
			gLoop->getSectors(sectors); // fill sector array


			//------------------------------------------------------------------------
			// prepare track sorting
			// clean vectors and index arrays
			sorter.cleanUp();
			sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);
			Int_t nCandHad     = sorter.fill(HParticleTrackSorter::selectHadrons);
			Int_t nCandHadBest = sorter.selectBest(Particle::kIsBestRKSorter,Particle::kIsHadronSorter);
//			Int_t nCandHadBest = sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsHadron);

			if( nCandHad<0 || nCandHadBest<0 ) return 1;

			//t0Reco.execute();

			vector<Int_t> IC_pim, IC_pip;

#if isRealData
			vector<HParticleCand *> vpim;
			vector<HParticleCand *> vpip;
			HParticleCand* cand=0;
#else
			vector<HParticleCandSim *> vpim;
			vector<HParticleCandSim *> vpip;
			HParticleCandSim* cand=0;
                        HGeantKine *kineRec=0;
#endif

			Int_t size = candCat->getEntries();
			for(Int_t j = 0; j < size; j++)
			{
			    cand = HCategoryManager::getObject(cand,candCat,j);
//#if !isRealData
//			    kineRec = HCategoryManager::getObject(kineRec,kineCat,cand->getGeantTrack()-1);
//#endif

			    if(cand) {

				if(!gLoop->goodSector(cand->getSector())) { continue;}  // skipp inactive sectors
				if(!cand->isFlagBit(kIsUsed)) continue;
				if(cand->isAtAnyMdcEdge()) continue;
                                if(cand->getSystemUsed() == -1) continue;

				Int_t   pid = cand->getPID();

				if(cand->getSystemUsed() == 1){
				    hPionMass_TOF[0]->Fill( cand->getCharge()*cand->getMass() );
				    if(cand->getChi2()<=400){
					hPionMass_TOF[1]->Fill( cand->getCharge()*cand->getMass() );

					if(cand->getMetaMatchQuality()<=3){
					    hPionMass_TOF[2]->Fill( cand->getCharge()*cand->getMass() );

					    if(cand->getMass()>0 && cand->getMass()<300 && cand->getMomentum()>0 && cand->getMomentum()<1000)
						hBetaMom_TOF_massCut->Fill( cand->getCharge()*cand->getMomentum(), cand->getBeta() );

					    if( betamom_2sig_pip_tof_pionCmom->IsInside(cand->getCharge()*cand->getMomentum(),cand->getBeta())
					       || betamom_2sig_pim_tof_pionCmom->IsInside(cand->getCharge()*cand->getMomentum(),cand->getBeta()) ){
						hPionMass_TOF[3]->Fill( cand->getCharge()*cand->getMass() );

						hMassMom_TOF_betamomCut->Fill( cand->getCharge()*cand->getMass(), cand->getMomentum() );
					    }
					}
				    }
				}
				if(cand->getSystemUsed() == 0){
				    hPionMass_RPC[0]->Fill( cand->getCharge()*cand->getMass() );
				    if(cand->getChi2()<=400){
					hPionMass_RPC[1]->Fill( cand->getCharge()*cand->getMass() );

					if(cand->getMetaMatchQuality()<=3){
					    hPionMass_RPC[2]->Fill( cand->getCharge()*cand->getMass() );

					    if(cand->getMass()>0 && cand->getMass()<300 && cand->getMomentum()>0 && cand->getMomentum()<1000)
						hBetaMom_RPC_massCut->Fill( cand->getCharge()*cand->getMomentum(), cand->getBeta() );

					    if( betamom_2sig_pip_rpc_pionCmom->IsInside(cand->getCharge()*cand->getMomentum(),cand->getBeta())
					       || betamom_2sig_pim_rpc_pionCmom->IsInside(cand->getCharge()*cand->getMomentum(),cand->getBeta()) ){
						hPionMass_RPC[3]->Fill( cand->getCharge()*cand->getMass() );

						hMassMom_RPC_betamomCut->Fill( cand->getCharge()*cand->getMass(), cand->getMomentum() );
					    }
					}
				    }
				}


				switch(PID_option){
				case 1:
				    {
					//First method of particle identification - Georgy's default one in dst-files
					if(pid==8){
					    IC_pip.push_back(j);
#if isRealData
					    vpip.push_back(new HParticleCand(*cand));
#else
					    vpip.push_back(new HParticleCandSim(*cand));
#endif
        				}
					if(pid==9){
					    IC_pim.push_back(j);
#if isRealData
					    vpim.push_back(new HParticleCand(*cand));
#else
					    vpim.push_back(new HParticleCandSim(*cand));
#endif
					}
				    }
				    break;

				case 2:
				    {
					//Second method of particle identification - check with HParticleTool functions
					if(cand->getCharge() > 0){
					    Float_t deltaEloss = 0., sigmaEloss = -1.;
					    if(! HParticleTool::isParticledEdx(8,cand,deltaEloss,sigmaEloss)) continue;
					    if( TMath::Abs(deltaEloss/sigmaEloss) > 2. ) continue;   //if more than 2 sigma difference not take it

					    Float_t deltaTime = 0., sigmaTime = -1.;
					    if(! HParticleTool::isParticleBeta(8,cand,2.,0.,5000.,deltaTime,sigmaTime,"apr12")) continue; //2. is sigma allowed, 0. minMom, 5000. maxMom

					    IC_pip.push_back(j);
#if isRealData
					    vpip.push_back(new HParticleCand(*cand));
#else
					    vpip.push_back(new HParticleCandSim(*cand));
#endif
					}
					if(cand->getCharge() < 0){
					    Float_t deltaEloss = 0., sigmaEloss = -1.;
					    if(! HParticleTool::isParticledEdx(9,cand,deltaEloss,sigmaEloss)) continue;
					    if( TMath::Abs(deltaEloss/sigmaEloss) > 2. ) continue;   //if more than 2 sigma difference not take it

					    Float_t deltaTime = 0., sigmaTime = -1.;
					    if(! HParticleTool::isParticleBeta(9,cand,2.,0.,5000.,deltaTime,sigmaTime,"apr12")) continue; //2. is sigma allowed, 0. minMom, 5000. maxMom

					    IC_pim.push_back(j);
#if isRealData
					    vpim.push_back(new HParticleCand(*cand));
#else
					    vpim.push_back(new HParticleCandSim(*cand));
#endif
					}
				    }
				    break;

				case 3:
				    {
					//Third method of particle identification - Timo's banana cuts
					if( (cand->getSystemUsed()==1 && betamom_2sig_pip_tof_pionCmom->IsInside(cand->getCharge()*cand->getMomentum(),cand->getBeta()))
					   || (cand->getSystemUsed()==0 && betamom_2sig_pip_rpc_pionCmom->IsInside(cand->getCharge()*cand->getMomentum(),cand->getBeta())) ){
					    IC_pip.push_back(j);
#if isRealData
					    vpip.push_back(new HParticleCand(*cand));
#else
					    vpip.push_back(new HParticleCandSim(*cand));
#endif

					}

					if( (cand->getSystemUsed()==1 && betamom_2sig_pim_tof_pionCmom->IsInside(cand->getCharge()*cand->getMomentum(),cand->getBeta()))
					   || (cand->getSystemUsed()==0 && betamom_2sig_pim_rpc_pionCmom->IsInside(cand->getCharge()*cand->getMomentum(),cand->getBeta())) ){
					    IC_pim.push_back(j);
#if isRealData
					    vpim.push_back(new HParticleCand(*cand));
#else
					    vpim.push_back(new HParticleCandSim(*cand));
#endif
        				}
				    }
				    break;

#if !isRealData
				case 4:
				    {
					if( cand->getGeantPID() == 8 ){ IC_pip.push_back(j); vpip.push_back(new HParticleCandSim(*cand)); }

					if( cand->getGeantPID() == 9 ){ IC_pim.push_back(j); vpim.push_back(new HParticleCandSim(*cand)); }
				    }
                                    break;
#endif

				case 5:
                                    //finaly used in Timo's disertation
				    {
					if( cand->getChi2()>0 && cand->getChi2()<400 && cand->getMetaMatchQuality()<3 ){
					    if( cand->getMass()>0 && cand->getMass()<300 && cand->getMomentum()>0 && cand->getMomentum()<1000){
						if( cand->getCharge()>0
//						   && cand->isOffVertexClust()
						  ){   
						    IC_pip.push_back(j);
#if isRealData
						    vpip.push_back(new HParticleCand(*cand));
#else
						    vpip.push_back(new HParticleCandSim(*cand));
#endif
						}

						if( cand->getCharge()<0
//						   && cand->isOffVertexClust()
						  ){
						    IC_pim.push_back(j);
#if isRealData
						    vpim.push_back(new HParticleCand(*cand));
#else
						    vpim.push_back(new HParticleCandSim(*cand));
#endif
						}
					    }
					}
				    }
				    break;

				default:
				    cout << "This PID_option is not supported. Please select from 1 to 5." << endl;
                                    break;
				}

			    }
			} // end cand loop

//                        cout << "Npip=" << IC_pip.size() << "   Npim=" << IC_pim.size() << "    ";

			eventmixer.nextEvent();
			eventmixer.addVector(vpim,9);
			eventmixer.addVector(vpip,8);

			//K0s loop
			if( (IC_pip.size() > 0) && (IC_pim.size() > 0) ){
#if isRealData
			    HParticleCand* pipCand;
			    HParticleCand* pimCand;
			    HVirtualCand* pipVCand;
			    HVirtualCand* pimVCand;
#else
			    HParticleCandSim* pipCand;
			    HParticleCandSim* pimCand;
			    HVirtualCandSim* pipVCand;
			    HVirtualCandSim* pimVCand;
#endif

			    TLorentzVector pipLV(1.,1.,1.,1.), pimLV(1.,1.,1.,1.), sumLV(1.,1.,1.,1.);
			    //Float_t metahit=-1, tof=-1000;

			    HGeomVector pipBase, pipDir, pimBase, pimDir, decayVertex, eventVertex, motherDir;
			    eventVertex.setXYZ(vertex.getX(),vertex.getY(),vertex.getZ());

			    Double_t distVert, dist1, dist2, dist3, distMin;

			    for(unsigned int ipip=0; ipip<IC_pip.size(); ipip++){
				pipCand = HCategoryManager::getObject(pipCand,candCat,IC_pip[ipip]);
				HParticleTool::fillTLorentzVector(pipLV,pipCand,8,useMomCorr);

				CalcSegVector(pipCand->getZ(), pipCand->getR(), pipCand->getPhi()*TMath::DegToRad(), pipCand->getTheta()*TMath::DegToRad(), pipBase, pipDir);


				for(unsigned int ipim=0; ipim<IC_pim.size(); ipim++){
				    pimCand = HCategoryManager::getObject(pimCand,candCat,IC_pim[ipim]);
				    HParticleTool::fillTLorentzVector(pimLV,pimCand,9,useMomCorr);

				    CalcSegVector(pimCand->getZ(), pimCand->getR(), pimCand->getPhi()*TMath::DegToRad(), pimCand->getTheta()*TMath::DegToRad(), pimBase, pimDir);

                                    //reconstructed LV from one pair of pions
				    sumLV = pipLV + pimLV;

                                    //getting the flow variable in range (-180,+180)
				    Double_t flow     = sumLV.Phi()*TMath::RadToDeg() - phiRP_tree;  // angle w.r.t. FW event plane
				    if (flow<-180.) flow += 360.;
				    if (flow>+180.) flow -= 360.;

				    // eff.  corrections - dependence on track multiplicity vs phi-phi_RP
				    Float_t eff=0., weight=0., weight_pim=0., weight_pip=0.;
				    eff=getEff(8,4-multClassNo,flow,pipLV.Theta()*TMath::RadToDeg());
				    if(eff!=0) weight_pip=1./eff;
				    eff=getEff(9,4-multClassNo,flow,pimLV.Theta()*TMath::RadToDeg());
				    if(eff!=0) weight_pim=1./eff;
				    if(weight_pip>10 || weight_pip<0 || weight_pim>10 || weight_pim<0) weight=0.;
				    else weight = weight_pip * weight_pim;
				    sekaon.occupancyCorr = weight;

                                    weight=0.;
				    eff=getAccCorr(multClassNo-1,sumLV.Pt(),sumLV.Rapidity()/0.822 - 1.);
				    if(eff!=0.) weight=1./eff;
				    if(weight>10000. || weight<1.) weight=0.;
				    sekaon.accCorr = weight;

				    weight=0.;
				    eff=getRecoEff(multClassNo-1,sumLV.Pt(),sumLV.Rapidity()/0.822 - 1.);
				    if(eff!=0.) weight=1./eff;
				    if(weight>10000. || weight<1.) weight=0.;
				    sekaon.recoEff = weight;


				    //topology variables
				    motherDir.setXYZ(sumLV.X(),sumLV.Y(),sumLV.Z());
				    decayVertex = calcVertexAnalytical(pipBase, pipDir, pimBase, pimDir);
        			    distVert = TMath::Sqrt( (decayVertex.getX()-eventVertex.getX())*(decayVertex.getX()-eventVertex.getX()) +
								    (decayVertex.getY()-eventVertex.getY())*(decayVertex.getY()-eventVertex.getY()) +
								    (decayVertex.getZ()-eventVertex.getZ())*(decayVertex.getZ()-eventVertex.getZ()) );
				    dist1 = calculateMinimumDistanceStraightToPoint(pimBase,     pimDir,    eventVertex);
				    dist2 = calculateMinimumDistanceStraightToPoint(pipBase,     pipDir,    eventVertex);
				    dist3 = calculateMinimumDistanceStraightToPoint(decayVertex, motherDir, eventVertex);
                                    distMin = calculateMinimumDistance(pimBase, pimDir, pipBase, pipDir);

				    TMVAd1 = dist1;
				    TMVAd2 = dist2;
				    TMVAd3 = dist3;
				    TMVAdVer = distVert;
				    TMVAdMin = distMin;
				    TMVAkaonMom = sumLV.P();

                                    sekaon.mvaResult = EvalMVAK0s();

                                    //storing the variables into tree
				    sekaon.invMass = sumLV.M();
				    sekaon.theta = sumLV.Theta()*TMath::RadToDeg();
				    sekaon.phi = ( sumLV.Phi()*TMath::RadToDeg()>0 ? sumLV.Phi()*TMath::RadToDeg() : sumLV.Phi()*TMath::RadToDeg()+360.) ;
				    sekaon.mom = sumLV.P();
                                    sekaon.momX = sumLV.Px();
                                    sekaon.momY = sumLV.Py();
				    sekaon.momZ = sumLV.Pz();
				    sekaon.pt = sumLV.Pt();
				    sekaon.mt = sumLV.Mt();
				    sekaon.y = sumLV.Rapidity();
				    sekaon.ycm = sumLV.Rapidity() - 0.822;
				    sekaon.y0 = sumLV.Rapidity()/0.822 - 1.;
				    sekaon.flow = flow;
                                    sekaon.decayVerX = decayVertex.getX();
                                    sekaon.decayVerY = decayVertex.getY();
                                    sekaon.decayVerZ = decayVertex.getZ();
				    sekaon.dist1 = dist1;
				    sekaon.dist2 = dist2;
				    sekaon.dist3 = dist3;
				    sekaon.distVer = distVert;
				    sekaon.distMin = distMin;
				    sekaon.oAngle = pimLV.Angle(pipLV.Vect());

                                    sekaon.downflag = downflag_se;


				    bool minimalTopologySet = false;
                                    //topology preCut set
				    if( dist1 > preCuts[iPrecuts][0] && dist2 > preCuts[iPrecuts][0] &&
				       dist3 < preCuts[iPrecuts][1] &&
				       distVert > preCuts[iPrecuts][2] &&
				       distMin < preCuts[iPrecuts][3] &&
				       decayVertex.getZ() > eventVertex.getZ()  )
					minimalTopologySet = true;

				    if( minimalTopologySet && sumLV.M() > 350. && sumLV.M() < 640. ){
					sekaontree->Fill();

                                        downflag_se = -1.;
				    }

				    hPiPairInvMass_all->Fill(sumLV.M());
				    if(minimalTopologySet) hPiPairInvMass_minTC->Fill(sumLV.M());
				    if(sekaon.mvaResult>0.5) hPiPairInvMass_mva50->Fill(sumLV.M());
				    if(sekaon.mvaResult>0.85) hPiPairInvMass_mva85->Fill(sumLV.M());
				    if(sekaon.mvaResult>0.95) hPiPairInvMass_mva95->Fill(sumLV.M());
				    if( dist1 > 10 && dist2 > 10 &&
				       dist3 < 8 &&
				       distVert > 15 &&
				       distMin < 10 &&
				       decayVertex.getZ() > eventVertex.getZ()  ) hPiPairInvMass_broadTC->Fill(sumLV.M());
				    if( dist1 > 14 && dist2 > 14 &&
				       dist3 < 7 &&
				       distVert > 29 &&
				       distMin < 8 &&
				       decayVertex.getZ() > eventVertex.getZ()  ) hPiPairInvMass_StoBratio->Fill(sumLV.M());
				    if( dist1 > 9 && dist2 > 9 &&
				       dist3 < 7 &&
				       distVert > 18 &&
				       distMin < 12 &&
				       decayVertex.getZ() > eventVertex.getZ()  ) hPiPairInvMass_Signif->Fill(sumLV.M());

				}
			    }

			    IC_pip.clear();
			    IC_pim.clear();
			    vpim.clear();
                            vpip.clear();



			    //Mixed Event background-------------------------------------------------
#if isRealData
//			    vector<pair<HParticleCand *, HParticleCand * > >& pairsVec = eventmixer.getMixedVector();
			    vector<HParticlePair>& pairsVec = eventmixer.getMixedVector();
#else
			    vector<pair<HParticleCandSim *, HParticleCandSim * > >& pairsVec = eventmixer.getMixedVector();
#endif
//                            cout << "nMEpairs=" << pairsVec.size() << "   ";

			    Int_t eventClass = eventmixer.currentEventClass();
//                            cout << "eventClass=" << eventClass << "   ";

			    HGeomVector eventVertexME(TMath::Floor((eventClass % 1500) / 150.) / 2. - 2.25,
						       TMath::Floor((eventClass %  150) /  15.) / 2. - 2.25,
						      ((eventClass % 15) * 3.577) - 54.7865);

                            mekaon.eventVerX = eventVertexME.getX();
                            mekaon.eventVerY = eventVertexME.getY();
			    mekaon.eventVerZ = eventVertexME.getZ();

//			    cout << "vertex: [ " << eventVertexME.getX() << ", " << eventVertexME.getY() << ", " << eventVertexME.getZ() << " ]" << endl;

			    for(UInt_t j = 0;j < pairsVec.size(); j++){

#if isRealData
				HParticlePair& pair = pairsVec[j];
                                //pair.setDoMomentumCorrection(useMomCorr);

                                //get candidates
				if( pair.getCandPID(0) == 8 && pair.getCandPID(1) == 9 ){
				    pipVCand = pair.getCand(0);
				    pimVCand = pair.getCand(1);
				}
				else if( pair.getCandPID(0) == 9 && pair.getCandPID(1) == 8 ){
				    pipVCand = pair.getCand(1);
				    pimVCand = pair.getCand(0);
				}
				else{
				    //cout << "Not correctly assigned particle cands!!" << endl;
                                    continue;
				}

				HParticleTool::fillTLorentzVector(pipLV,pipVCand,8,useMomCorr);
				HParticleTool::fillTLorentzVector(pimLV,pimVCand,9,useMomCorr);

				CalcSegVector(pipVCand->getZ(), pipVCand->getR(), pipVCand->getPhi()*TMath::DegToRad(), pipVCand->getTheta()*TMath::DegToRad(), pipBase, pipDir);
				CalcSegVector(pimVCand->getZ(), pimVCand->getR(), pimVCand->getPhi()*TMath::DegToRad(), pimVCand->getTheta()*TMath::DegToRad(), pimBase, pimDir);

#else
				pair<HParticleCandSim *, HParticleCandSim * >& pair = pairsVec[j];

				//get sim candidates
				HParticleCandSim *firstPCandSim = pair.first;
				HParticleCandSim *secondPCandSim = pair.second;

#if isEmbedding
				if( firstPCandSim->getCharge()>0 && secondPCandSim->getCharge()<0 ){
				    pipCand = pair.first;
				    pimCand = pair.second;
				}
				else if( firstPCandSim->getCharge()<0 && secondPCandSim->getCharge()>0 ){
				    pipCand = pair.second;
				    pimCand = pair.first;
				}
				else{
				    //cout << "Not correctly assigned particle cands!!" << endl;
                                    continue;
				}
#else
				if( firstPCandSim->getGeantPID() == 8 && secondPCandSim->getGeantPID() == 9 ){
				    pipCand = pair.first;
				    pimCand = pair.second;
				}
				else if( firstPCandSim->getGeantPID() == 9 && secondPCandSim->getGeantPID() == 8 ){
				    pipCand = pair.second;
				    pimCand = pair.first;
				}
				else{
				    //cout << "Not correctly assigned particle cands!!" << endl;
                                    continue;
				}
#endif

				HParticleTool::fillTLorentzVector(pipLV,pipCand,8,useMomCorr);
				HParticleTool::fillTLorentzVector(pimLV,pimCand,9,useMomCorr);

				CalcSegVector(pipCand->getZ(), pipCand->getR(), pipCand->getPhi()*TMath::DegToRad(), pipCand->getTheta()*TMath::DegToRad(), pipBase, pipDir);
				CalcSegVector(pimCand->getZ(), pimCand->getR(), pimCand->getPhi()*TMath::DegToRad(), pimCand->getTheta()*TMath::DegToRad(), pimBase, pimDir);

#endif

                                sumLV = pipLV + pimLV;

				motherDir.setXYZ(sumLV.X(),sumLV.Y(),sumLV.Z());
				decayVertex = calcVertexAnalytical(pipBase, pipDir, pimBase, pimDir);

				distVert = TMath::Sqrt( (decayVertex.getX()-eventVertexME.getX())*(decayVertex.getX()-eventVertexME.getX()) +
								(decayVertex.getY()-eventVertexME.getY())*(decayVertex.getY()-eventVertexME.getY()) +
								(decayVertex.getZ()-eventVertexME.getZ())*(decayVertex.getZ()-eventVertexME.getZ()) );
				dist1 = calculateMinimumDistanceStraightToPoint(pimBase,     pimDir,    eventVertexME);
				dist2 = calculateMinimumDistanceStraightToPoint(pipBase,     pipDir,    eventVertexME);
				dist3 = calculateMinimumDistanceStraightToPoint(decayVertex, motherDir, eventVertexME);
				distMin = calculateMinimumDistance(pimBase, pimDir, pipBase, pipDir);

                                //getting the flow value in range (-180,+180)
				Double_t flow_ME     = sumLV.Phi()*TMath::RadToDeg() - phiRP_tree;  // angle w.r.t. FW event plane
				if (flow_ME<-180.) flow_ME += 360.;
				if (flow_ME>+180.) flow_ME -= 360.;

				// eff.  corrections - dependence on track multiplicity vs phi-phi_RP
				Float_t eff=0., weight=0., weight_pim=0., weight_pip=0.;
				eff=getEff(8,4-multClassNo,flow_ME,pipLV.Theta()*TMath::RadToDeg());
				if(eff!=0) weight_pip=1./eff;
				eff=getEff(9,4-multClassNo,flow_ME,pimLV.Theta()*TMath::RadToDeg());
				if(eff!=0) weight_pim=1./eff;
				if(weight_pip>10 || weight_pip<0 || weight_pim>10 || weight_pim<0) weight=0.;
				else weight = weight_pip * weight_pim;
				mekaon.occupancyCorr = weight;

				weight=0.;
				eff=getAccCorr(multClassNo-1,sumLV.Pt(),sumLV.Rapidity()/0.822 - 1.);
				if(eff!=0.) weight=1./eff;
				if(weight>10000. || weight<1.) weight=0.;
				mekaon.accCorr = weight;

				weight=0.;
				eff=getRecoEff(multClassNo-1,sumLV.Pt(),sumLV.Rapidity()/0.822 - 1.);
				if(eff!=0.) weight=1./eff;
				if(weight>10000. || weight<1.) weight=0.;
				mekaon.recoEff = weight;


				TMVAd1 = dist1;
				TMVAd2 = dist2;
				TMVAd3 = dist3;
				TMVAdVer = distVert;
				TMVAdMin = distMin;
				TMVAkaonMom = sumLV.P();

				mekaon.mvaResult = EvalMVAK0s();


				//tree variables filling
				mekaon.invMass = sumLV.M();
				mekaon.theta = sumLV.Theta()*TMath::RadToDeg();
				mekaon.phi = ( sumLV.Phi()*TMath::RadToDeg()>0 ? sumLV.Phi()*TMath::RadToDeg() : sumLV.Phi()*TMath::RadToDeg()+360.) ;
				mekaon.mom = sumLV.P();
				mekaon.momX = sumLV.Px();
				mekaon.momY = sumLV.Py();
				mekaon.momZ = sumLV.Pz();
				mekaon.pt = sumLV.Pt();
				mekaon.mt = sumLV.Mt();
				mekaon.y = sumLV.Rapidity();
				mekaon.ycm = sumLV.Rapidity() - 0.822;
				mekaon.y0 = sumLV.Rapidity()/0.822 - 1.;
				mekaon.flow = flow_ME;
				mekaon.decayVerX = decayVertex.getX();
				mekaon.decayVerY = decayVertex.getY();
				mekaon.decayVerZ = decayVertex.getZ();
				mekaon.dist1 = dist1;
				mekaon.dist2 = dist2;
				mekaon.dist3 = dist3;
				mekaon.distVer = distVert;
				mekaon.distMin = distMin;
				mekaon.oAngle = pimLV.Angle(pipLV.Vect());

                                mekaon.downflag = downflag_me;


				bool minimalTopologySet = false;
				//topology preCut set
				if( dist1 > preCuts[iPrecuts][0] && dist2 > preCuts[iPrecuts][0] &&
				   dist3 < preCuts[iPrecuts][1] &&
				   distVert > preCuts[iPrecuts][2] &&
				   distMin < preCuts[iPrecuts][3] &&
				   decayVertex.getZ() > eventVertexME.getZ()  )
				    minimalTopologySet = true;

				if( minimalTopologySet && sumLV.M() > 350 && sumLV.M() < 640 ){
				    //fill the mix-event tree
				    mekaontree->Fill();

                                    downflag_me = -1.;
				}

			    }

			}

			eventmixer.RemoveObjectsToDelete();

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

//        delete effcorFile;

//	if(!out) out = new TFile(outfile,"UPDATE");
	if(out) {
	    out->cd();
	    out->mkdir("macros");
	    out->cd("macros");
	    TMacro *mAna = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/analysis.cc");
	    TMacro *mLoop = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/loopDST_task.C");
	    TMacro *mPartReco = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/KZeroTrees.h");
	    TMacro *mMVA = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/MVAK0s.h");
	    TMacro *mMEclass = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/MEclassifiers.h");
	    TMacro *mGlobVars = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/GlobVars.h");

            mAna->Write();
            mLoop->Write();
            mPartReco->Write();
	    mMVA->Write();
	    mMEclass->Write();
            mGlobVars->Write();

	    out->cd("kZero");
	    sekaontree->Write();
	    mekaontree->Write();


            out->cd("eventPlanes");
	    for(int icc=0; icc<nCentralityClasses; icc++){
		heventplane[icc]->Write();
                hphiAB[icc]->Write();
	    }

            out->cd("histos");
	    hCentralityCounter->Write();
            hEventCounter->Write();
	    hMassMom_TOF_betamomCut->Write();
            hMassMom_RPC_betamomCut->Write();
            hBetaMom_TOF_massCut->Write();
            hBetaMom_RPC_massCut->Write();
	    for(int ih=0; ih<4; ih++){
		hPionMass_TOF[ih]->Write();
                hPionMass_RPC[ih]->Write();
	    }
	    hPiPairInvMass_all->Write();
            hPiPairInvMass_minTC->Write();
            hPiPairInvMass_mva50->Write();
            hPiPairInvMass_mva85->Write();
            hPiPairInvMass_mva95->Write();
            hPiPairInvMass_broadTC->Write();
            hPiPairInvMass_StoBratio->Write();
            hPiPairInvMass_Signif->Write();



#if !isRealData
            out->cd("simulation");
	    gkzerotree->Write();
#endif

//	    out->Save();
//	    out->Close();

	}

	return kTRUE;
    }
   //ClassDef(KZeroTrees,0) // no streamer, do not use if you include blank header in ACLIC
};

#endif
