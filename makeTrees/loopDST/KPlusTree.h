#ifndef __KPLUSTREE___     // good style to avaoid multiple includes
#define __KPLUSTREE___

#include "mylibs.h"
#include "mystruc.h"

#include "GlobVars.h"

#include "UsefulFunctions.h"

using namespace std;

class KPlusTree : public HReconstructor
{
protected:
    // put all vars here which are not
    // local to a function

    static const Int_t particleID = 11;
    static const Int_t PID_option = 2; // 5 - geantPID (only sim!!)
    static const Int_t iTrackCuts = 5;
//    static const Int_t EP_method = 3; //1 and 2 - Pavel (debug mode), 3 - my, 4 - Gosia

    HParticleTrackSorter sorter;
    //HParticleT0Reco t0Reco("apr12");

    CAND_KPLUS sekaon;  //se = same event
    TTree *sekaontree;
#ifndef isRealData
    GKINE_KPLUS gkplus;
    TTree *gkplustree;
#endif
    
    //------------------------------------
    // file handling
    TFile* effcorFile;
    TFile *mdcdEdx_cutfile_new_rpc;
    TFile *mdcdEdx_cutfile_new_tof;
    TFile *tofdEdx_cutfile;
    TFile *kplusAccRecEffcorrFile;


    TCutG *mdcdEdx_tof_Kp_data_new, *mdcdEdx_rpc_Kp_data_new, *tofdEdx_beta_noK;

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

    Float_t getAccRecEff(Int_t pid, Float_t pt, Float_t y)
    {
	Float_t mt_m0 =  TMath::Sqrt(HPhysicsConstants::mass(pid)*HPhysicsConstants::mass(pid) + pt*pt) - HPhysicsConstants::mass(pid);

	Float_t eff = -1;
	switch(pid){

	case 11:
                  if(kplusAccRecEff) eff = kplusAccRecEff->GetBinContent(kplusAccRecEff->FindBin(y,mt_m0));
            break;

	case 16:
                  if(kzeroAccRecEff) eff = kzeroAccRecEff->GetBinContent(kzeroAccRecEff->FindBin(y,mt_m0));
            break;

	default:
	    cerr<< " PID "<<pid<<" not known -> functions only for K+ (11) and K0s (16)"<<endl;
	    break;

	}

        return eff;
    }


public:

    KPlusTree(const Text_t *name = "",const Text_t *title ="",TString outfile="KPlusTree.root")
	: HReconstructor(name,title)
    {  // init your vars

//	out = TFile::Open(TFile::AsyncOpen(outfile,"UPDATE"));

    }

    virtual ~KPlusTree()
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

//	if(!out) out = new TFile(outfile,"UPDATE");
	if(out) {
	    out->cd();
	    out->mkdir("kaonPlus");
	    out->cd("kaonPlus");

	    sekaontree = new TTree("sekaontree","Same Event Kaon Tree");
	    sekaontree->Branch("sekaon",&sekaon,"eventPlane:phiAB:"
			       "eventVerX:eventVerY:eventVerZ:evntVerChi2:"
			       "system:theta:phi:momX:momY:momZ:mom:"
			       "multiplicity:centralityclass:MDCtracks:METAclst:"
			       "invMass:mt:mtm0:pt:y:ycm:y0:"
			       "flow:occupancyCorr:accreceffCorr:eventWeight:downflag");

#ifndef isRealData
	    gkplustree = new TTree("gkplustree","GeantKine Kaon Tree");
	    gkplustree->Branch("gkplus",&gkplus,"eventPlane_dst:eventPlane_my:eventPlane_gkine:phiAB_dst:phiAB_my:"
			      "system:theta:phi:momX:momY:momZ:mom:"
			      "multiplicity:centralityclass:impactparam:MDCtracks:METAclst:"
			      "invMass:mt:mtm0:pt:y:ycm:flow:flowDST:flowGKINE:weight:"
			      "kpAcc/B:kpRec");
#endif

	    const char *effcorFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/MY_KpIII_effcor_flow_matrix_apr12_day108_gen8_4mbins.root";
	    effcorFile = TFile::Open( effcorFileName );
	    effPiM_cc3 = (TH2F *)effcorFile->Get("heff4_id9");
	    effPiP_cc3 = (TH2F *)effcorFile->Get("heff4_id8");
	    effKp_cc3  = (TH2F *)effcorFile->Get("heff4_id11");
	    effP_cc3   = (TH2F *)effcorFile->Get("heff4_id14");
	    effPiM_cc2 = (TH2F *)effcorFile->Get("heff3_id9");
	    effPiP_cc2 = (TH2F *)effcorFile->Get("heff3_id8");
	    effKp_cc2  = (TH2F *)effcorFile->Get("heff3_id11");
	    effP_cc2   = (TH2F *)effcorFile->Get("heff3_id14");
	    effPiM_cc1 = (TH2F *)effcorFile->Get("heff2_id9");
	    effPiP_cc1 = (TH2F *)effcorFile->Get("heff2_id8");
	    effKp_cc1  = (TH2F *)effcorFile->Get("heff2_id11");
	    effP_cc1   = (TH2F *)effcorFile->Get("heff2_id14");
	    effPiM_cc0 = (TH2F *)effcorFile->Get("heff1_id9");
	    effPiP_cc0 = (TH2F *)effcorFile->Get("heff1_id8");
	    effKp_cc0  = (TH2F *)effcorFile->Get("heff1_id11");
	    effP_cc0   = (TH2F *)effcorFile->Get("heff1_id14");

            const char *kplusAccRecEffcorrFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/Happy_Kp_fit_Rk100_MM2_40centr_draw.root";
	    kplusAccRecEffcorrFile = TFile::Open( kplusAccRecEffcorrFileName );
            kplusAccRecEff = (TH2F*) kplusAccRecEffcorrFile->Get("hmty_eff_sys0");

	    const char *tofcutFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/Kp_cutfile_tof_gen8a.root";
	    mdcdEdx_cutfile_new_tof = TFile::Open( tofcutFileName );
	    mdcdEdx_tof_Kp_data_new =(TCutG*)mdcdEdx_cutfile_new_tof->Get("KP_mdcdEdx_cut_scaled");

            const char *rpccutFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/Kp_cutfile_rpc_gen8a_rpc25.root";
	    mdcdEdx_cutfile_new_rpc = TFile::Open( rpccutFileName );
	    mdcdEdx_rpc_Kp_data_new =(TCutG*)mdcdEdx_cutfile_new_rpc->Get("KP_mdcdEdx_cut_scaled");

	    const char *kaonVetocutFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/betaTofdEdx_noKm.root";
	    tofdEdx_cutfile = TFile::Open( kaonVetocutFileName );
	    tofdEdx_beta_noK =(TCutG*)tofdEdx_cutfile->Get("betaTofdEdx_noKm");


	    for(int icc=0; icc<nCentralityClasses; icc++){
		heventplane[icc] = new TH1F(Form("heventplane_icc%i",icc),"event plane distribution",200,-200,200);
		hphiAB[icc] = new TH1F(Form("hphiAB_icc%i",icc),"difference of event plane for two subsets",180,0,180);
	    }

	    hCentralityCounter = new TH1I("hCentralityCounter","number of events in each centrality class",10,0,10);

	}

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
	effP_cc3   = (TH2F *)effcorFile->Get("heff4_id14");
	effPiM_cc2 = (TH2F *)effcorFile->Get("heff3_id9");
	effPiP_cc2 = (TH2F *)effcorFile->Get("heff3_id8");
	effKp_cc2  = (TH2F *)effcorFile->Get("heff3_id11");
	effP_cc2   = (TH2F *)effcorFile->Get("heff3_id14");
	effPiM_cc1 = (TH2F *)effcorFile->Get("heff2_id9");
	effPiP_cc1 = (TH2F *)effcorFile->Get("heff2_id8");
	effKp_cc1  = (TH2F *)effcorFile->Get("heff2_id11");
	effP_cc1   = (TH2F *)effcorFile->Get("heff2_id14");
	effPiM_cc0 = (TH2F *)effcorFile->Get("heff1_id9");
	effPiP_cc0 = (TH2F *)effcorFile->Get("heff1_id8");
	effKp_cc0  = (TH2F *)effcorFile->Get("heff1_id11");
	effP_cc0   = (TH2F *)effcorFile->Get("heff1_id14");

	const char *kplusAccRecEffcorrFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/Happy_Kp_fit_Rk100_MM2_40centr_draw.root";
	kplusAccRecEffcorrFile = TFile::Open( kplusAccRecEffcorrFileName );
	kplusAccRecEff = (TH2F*) kplusAccRecEffcorrFile->Get("hmty_eff_sys0");

	const char *tofcutFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/Kp_cutfile_tof_gen8a.root";
	mdcdEdx_cutfile_new_tof = TFile::Open( tofcutFileName );
	mdcdEdx_tof_Kp_data_new =(TCutG*)mdcdEdx_cutfile_new_tof->Get("KP_mdcdEdx_cut_scaled");

	const char *rpccutFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/Kp_cutfile_rpc_gen8a_rpc25.root";
	mdcdEdx_cutfile_new_rpc = TFile::Open( rpccutFileName );
	mdcdEdx_rpc_Kp_data_new =(TCutG*)mdcdEdx_cutfile_new_rpc->Get("KP_mdcdEdx_cut_scaled");

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

	if(gHades){
	    HCategory* evtInfoCat = (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo);
	    HCategory* candCat    = (HCategory*)HCategoryManager::getCategory(catParticleCand);

	    HEventHeader* header  = gHades->getCurrentEvent()->getHeader();
	    HVertex vertex      = header->getVertexReco();

	    HCategory* evtPlaneCat = (HCategory*)HCategoryManager::getCategory(catWallEventPlane);

            sekaon.eventVerX = vertex.getX();
            sekaon.eventVerY = vertex.getY();
	    sekaon.eventVerZ = vertex.getZ();
	    sekaon.evntVerChi2 = vertex.getChi2();

	    Float_t downflag = 1.;

#ifndef isRealData
	    HGeantKine*       kine;
	    HLinearCategory* kineCat = (HLinearCategory*)HCategoryManager::getCategory(catGeantKine);
//	    HGeantHeader* fSubHeader = NULL;
//	    HPartialEvent* fSimul = NULL;
#endif

	    if(evtInfoCat)
	    {
		HParticleEvtInfo* evtInfo=0;
		evtInfo = HCategoryManager::getObject(evtInfo,evtInfoCat,0 );

		HWallEventPlane* event_plane;
		event_plane = HCategoryManager::getObject(event_plane,evtPlaneCat,0);

		Bool_t isGoodEvent = kFALSE;
#ifdef isRealData
		if(evtInfo && evtInfo->isGoodEvent(Particle::kGoodTRIGGER|
						   Particle::kGoodVertexClust|
						   Particle::kGoodVertexCand|
						   Particle::kGoodSTART|
						   Particle::kNoPileUpSTART|
						   Particle::kNoVETO|
						   Particle::kGoodSTARTVETO|
						   Particle::kGoodSTARTMETA
						  )
		  ) isGoodEvent = kTRUE;
#endif
#ifndef isRealData
		if(evtInfo) isGoodEvent = kTRUE;
#endif

                if(isGoodEvent)
		{
		    //-------------------------------------------------
		    // loop over wall hits to get the Reaction Plane angle and resolution
		    Double_t mult = evtInfo->getSumTofMultCut()+evtInfo->getSumRpcMultHitCut();

/*#ifndef isRealData
		    //TODO: from Heidi how she get the evet plane from simulations
		    fSimul = ( (HRecEvent*)gHades->getCurrentEvent())->getPartialEvent(catSimul);
		    Double_t eventPlaneAngle = -1.;
                    Double_t impactParameter = -1.;

		    fSubHeader = (HGeantHeader*)(fSimul->getSubHeader());

		    if (fSubHeader) {
			impactParameter = fSubHeader->getImpactParameter();
			eventPlaneAngle = fSubHeader->getEventPlane();
		    }
#endif*/

		    Int_t behruzCentralityClass = 0;
		    behruzCentralityClass = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k5);  //is saved in tree!!

		    Int_t multClassNo = 0;
		    multClassNo = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10);  //needed for occupancy corr.

                    Float_t  event_weight = evtChara.getEventWeight();

		    Float_t phiRP = evtChara.getEventPlane(HParticleEvtChara::kDefault)*TMath::RadToDeg();
		    Float_t phiRP_A = evtChara.getEventPlane(HParticleEvtChara::kDefault,1)*TMath::RadToDeg();
		    Float_t phiRP_B = evtChara.getEventPlane(HParticleEvtChara::kDefault,2)*TMath::RadToDeg();

		    if( (phiRP == -1*TMath::RadToDeg()) || (phiRP_A == -1*TMath::RadToDeg()) || (phiRP_B == -1*TMath::RadToDeg()) ) return 1;  

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

		    //my NEW centrality class -> see UsefulFunctions.h for more info
		    vector<int> myVecCC = getMyCentralityClass( behruzCentralityClass );
		    if( myVecCC.size() == 0 ) return 1; //skipping event outside of region of interest in centrality
		    for(int j=0; j<myVecCC.size(); j++){
			heventplane[myVecCC[j]]->Fill(phiRP_tree);
			hphiAB[myVecCC[j]]->Fill(phiRPAB_tree);
		    }

                    //fill counter with appropriate weight
		    hCentralityCounter->Fill(behruzCentralityClass,event_weight);


		    sekaon.eventPlane      = phiRP_tree;
		    sekaon.phiAB           = phiRPAB_tree;
		    sekaon.multiplicity    = mult;
		    sekaon.centralityclass = behruzCentralityClass;
                    sekaon.eventWeight     = event_weight;

#ifndef isRealData
		    gkplus.eventPlane_my = phiRP_my;
		    gkplus.eventPlane_dst = phiRP_dst;
		    //gkplus.eventPlane_gkine = eventPlaneAngle;
                    gkplus.eventPlane_gkine = -10000.;
		    gkplus.phiAB_my = phiRPAB_my;
		    gkplus.phiAB_dst = phiRPAB_dst;
		    gkplus.multiplicity = mult;
		    gkplus.centralityclass = behruzCentralityClass;
		    //gkplus.impactparam = impactParameter;
                    gkplus.impactparam = -10000.;
#endif



#ifndef isRealData
       //---------------------------------- Kine --------------------------------------------------------
       for(Int_t j = 0; j < kineCat->getEntries(); j ++) {
	   kine = HCategoryManager::getObject(kine,kineCat,j);
	   if( kine->isPrimary() && kine->getID() == particleID && kine->getGeneratorInfo() == (500+particleID) ){

               gkplus.invMass = HPhysicsConstants::mass(particleID);
	       gkplus.pt   = kine->getTransverseMomentum();
	       gkplus.mt   = sqrt(gkplus.pt*gkplus.pt+gkplus.invMass*gkplus.invMass);
               gkplus.mtm0 = sqrt(gkplus.pt*gkplus.pt+gkplus.invMass*gkplus.invMass) - HPhysicsConstants::mass(particleID);
	       gkplus.y    = kine->getRapidity();
               gkplus.ycm  = kine->getRapidity() - 0.74;

               kine->getMomentum(gkplus.momX,gkplus.momY,gkplus.momZ);
	       gkplus.mom = kine->getTotalMomentum();

	       Int_t sec_kine;
	       Float_t phi_kine, theta_kine;
	       sec_kine = kine->getPhiThetaDeg(theta_kine,phi_kine);

	       gkplus.theta = theta_kine;
               if( phi_kine > 180. ) phi_kine -= 360.;
	       gkplus.phi = phi_kine;

	       gkplus.flow = phi_kine - phiRP_my;
	       gkplus.flowDST = phi_kine - phiRP_dst;
               //gkplus.flowGKINE = phi_kine - eventPlaneAngle;
	       gkplus.flowGKINE = -10000.;
               if(gkplus.flow < -180.) gkplus.flow += 360.;
               if(gkplus.flow > +180.) gkplus.flow -= 360.;
               if(gkplus.flowDST < -180.) gkplus.flowDST += 360.;
               if(gkplus.flowDST > +180.) gkplus.flowDST -= 360.;

	       gkplus.weight = 1.;

	       gkplus.MDCtracks = -10000.;
	       gkplus.METAclst = -10000.;

	       gkplus.kpAcc = false;
               gkplus.kpRec = false;

               gkplustree->Fill();
	   }
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
			Int_t nCandHadBest = sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsHadron);

			if( nCandHad<0 || nCandHadBest<0 ) return 1;

			//t0Reco.execute();



			vector<Int_t> IC_kp;
#ifdef isRealData
			vector<HParticleCand *> vkp;
			HParticleCand* cand=0;
#endif
#ifndef isRealData
			vector<HParticleCandSim *> vkp;
			HParticleCandSim* cand=0;
                        HGeantKine *kineRec=0;
#endif

			Int_t size = candCat->getEntries();
			for(Int_t j = 0; j < size; j++)
			{
			    cand = HCategoryManager::getObject(cand,candCat,j);
#ifndef isRealData
			    kineRec = HCategoryManager::getObject(kineRec,kineCat,cand->getGeantTrack()-1);
#endif

			    if(cand) {

				if(!gLoop->goodSector(cand->getSector())) { continue;}  // skipp inactive sectors
				if(!cand->isFlagBit(kIsUsed)) continue;
				if(cand->isAtAnyMdcEdge()) continue;
                                if(cand->getSystemUsed() == -1) continue;

				//if(cand->getPID()<0)          continue;   //THIS IS VERY BAD BAD LINE!!!!!!!!!

				//additional cuts..
				//if(cand->getChi2()<0 || cand->getChi2()>500) continue;
				//if(cand->getMetaMatchQuality()<0 || cand->getMetaMatchQuality()>5) continue;
                                //if(cand->getInnerSegmentChi2()<0 || cand->getOuterSegmentChi2()<0) continue;

				Int_t   pid = cand->getPID();
				Bool_t goodPBeta=false, goodPEloss=false;
				Bool_t goodMdcEloss=false, goodTofEloss=false, goodChi2RK=false, goodChi2MM=false, goodChargeMass=false, goodMom=false;

				switch(PID_option){
				case 1:
				    {
					//First method of particle identification - Georgy's default one in dst-files
					if(pid==11){
					    IC_kp.push_back(j);
#ifdef isRealData
					    vkp.push_back(new HParticleCand(*cand));
#endif
#ifndef isRealData
					    vkp.push_back(new HParticleCandSim(*cand));
#endif
					}
				    }
				    break;

/*				case 2:
				    {
					//Heidi's cuts on MDC dEdx and TOF dEdx veto
					if( ( cand->getSystemUsed()==1 && mdcdEdx_tof_Kp_data_new->IsInside(cand->getCharge()*cand->getMomentum(),cand->getMdcdEdx())
					                               && !tofdEdx_beta_noK->IsInside(cand->getBeta(),cand->getTofdEdx()) )
					   || ( cand->getSystemUsed()==0 && mdcdEdx_rpc_Kp_data_new->IsInside(cand->getCharge()*cand->getMomentum(),cand->getMdcdEdx()) )
					  ){
					    goodMdcEloss = true;
					    goodTofEloss = true;
					}

					//plus the track quality cuts
					if( cand->getChi2() <= 100. )
					    goodChi2RK = true;
					if( cand->getMetaMatchQuality() <= 2. )
					    goodChi2MM = true;

					//plus charge and mass check
//					if( cand->getCharge()==1 && cand->getMass()>60 && cand->getMass()<1200 )    //test of background subtraction for day108
					if( cand->getCharge()==1 && cand->getMass()>300 && cand->getMass()<700 )
					    goodChargeMass = true;

					//plus momentum cut
					if( cand->getMomentum()>0 && cand->getMomentum()<1200 )
                                            goodMom = true;

					//Heidi's cuts together
					if( goodMdcEloss &&
					    goodTofEloss &&
					    goodChi2RK && goodChi2MM && goodChargeMass && goodMom ){
					    IC_kp.push_back(j);
#ifdef isRealData
					    vkp.push_back(new HParticleCand(*cand));
#endif
#ifndef isRealData
					    vkp.push_back(new HParticleCandSim(*cand));
#endif
					}

				    }
				    break;
*/
/*				case 2:
				    {
					//Heidi's cuts on MDC dEdx and TOF dEdx veto
					if( iTrackCuts>=0 && iTrackCuts <=2 ){  //both MDC and TOF are used
					    if( ( cand->getSystemUsed()==1 && mdcdEdx_tof_Kp_data_new->IsInside(cand->getCharge()*cand->getMomentum(),cand->getMdcdEdx())
						                           && !tofdEdx_beta_noK->IsInside(cand->getBeta(),cand->getTofdEdx()) )
					       || ( cand->getSystemUsed()==0 && mdcdEdx_rpc_Kp_data_new->IsInside(cand->getCharge()*cand->getMomentum(),cand->getMdcdEdx()) )
					      )
						goodEloss = true;
					}
					if( iTrackCuts==3 ){ //only TOF is used
					    if( ( cand->getSystemUsed()==1 && !tofdEdx_beta_noK->IsInside(cand->getBeta(),cand->getTofdEdx()) )
                                               || ( cand->getSystemUsed()==0 )
					      )
						goodEloss = true;

					}
					if( iTrackCuts==4 ){ //only MDC is used
					    if( ( cand->getSystemUsed()==1 && mdcdEdx_tof_Kp_data_new->IsInside(cand->getCharge()*cand->getMomentum(),cand->getMdcdEdx()) )
					       || ( cand->getSystemUsed()==0 && mdcdEdx_rpc_Kp_data_new->IsInside(cand->getCharge()*cand->getMomentum(),cand->getMdcdEdx()) )
					      )
						goodEloss = true;
					}
					if( iTrackCuts==5 ){ //no Eloss cuts are used
					    goodEloss = true;
					}

					//plus the track quality cuts
					if( cand->getChi2() <= trackQualityCuts[iTrackCuts][0] )
					    goodChi2RK = true;

					if( cand->getMetaMatchQuality() <= trackQualityCuts[iTrackCuts][1] )
					    goodChi2MM = true;

					//plus charge and mass check
					if( cand->getCharge()==1 && cand->getMass()>300 && cand->getMass()<700 )
					    goodChargeMass = true;

					//plus momentum cut
					if( cand->getMomentum()>0 && cand->getMomentum()<1200 )
                                            goodMom = true;

					//Heidi's cuts together
					if( goodEloss && goodChi2RK && goodChi2MM && goodChargeMass && goodMom ){
					    IC_kp.push_back(j);
#ifdef isRealData
					    vkp.push_back(new HParticleCand(*cand));
#endif
#ifndef isRealData
					    vkp.push_back(new HParticleCandSim(*cand));
#endif
					}

				    }
				    break;
*/
				case 2:
				    {
					//Heidi's cuts on MDC dEdx and TOF dEdx veto
					if( ( cand->getSystemUsed()==1 && mdcdEdx_tof_Kp_data_new->IsInside(cand->getCharge()*cand->getMomentum(),cand->getMdcdEdx()) )
					   || ( cand->getSystemUsed()==0 && mdcdEdx_rpc_Kp_data_new->IsInside(cand->getCharge()*cand->getMomentum(),cand->getMdcdEdx()) )
					  )
					    goodMdcEloss = true;

					if( ( cand->getSystemUsed()==1 && !tofdEdx_beta_noK->IsInside(cand->getBeta(),cand->getTofdEdx()) )
					   || ( cand->getSystemUsed()==0 )
					  )
					    goodTofEloss = true;

					//plus the track quality cuts
					if( cand->getChi2() <= trackQualityCuts[iTrackCuts][0] )
					    goodChi2RK = true;

					if( cand->getMetaMatchQuality() <= trackQualityCuts[iTrackCuts][1] )
					    goodChi2MM = true;

					//plus charge and mass check
					if( cand->getCharge()==1 && cand->getMass()>300 && cand->getMass()<700 )
					    goodChargeMass = true;

					//plus momentum cut
					if( cand->getMomentum()>0 && cand->getMomentum()<1200 )
                                            goodMom = true;

					//Heidi's cuts together
					if( ( (trackElossCuts[iTrackCuts][0] && goodMdcEloss) || (!trackElossCuts[iTrackCuts][0]) ) &&
					    ( (trackElossCuts[iTrackCuts][1] && goodTofEloss) || (!trackElossCuts[iTrackCuts][1]) ) &&
					   goodChi2RK && goodChi2MM && goodChargeMass && goodMom ){

					    IC_kp.push_back(j);
#ifdef isRealData
					    vkp.push_back(new HParticleCand(*cand));
#endif
#ifndef isRealData
					    vkp.push_back(new HParticleCandSim(*cand));
#endif
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
#ifdef isRealData
					    vkp.push_back(new HParticleCand(*cand));
#endif
#ifndef isRealData
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
#ifdef isRealData
					    vkp.push_back(new HParticleCand(*cand));
#endif
#ifndef isRealData
					    vkp.push_back(new HParticleCandSim(*cand));
#endif
					}
				    }
				    break;
#ifndef isRealData
				case 5:
				    {
					if( cand->getGeantPID() == 11 ) vkp.push_back(new HParticleCandSim(*cand));
				    }
                                    break;
#endif

				default:
				    cout << "This PID_option is not supported. Please select 1/2/3." << endl;
                                    break;
				}
			    }
			} // end cand loop


			//K+ loop
			if( IC_kp.size() > 0 ){
#ifdef isRealData
			    HParticleCand* kpCand;
#endif
#ifndef isRealData
			    HParticleCandSim* kpCand;
#endif

			    TLorentzVector kplusLV(1.,1.,1.,1.);
			    //Float_t metahit=-1, tof=-1000;

			    for(unsigned int ikp=0; ikp<IC_kp.size(); ikp++){
				kpCand = HCategoryManager::getObject(kpCand,candCat,IC_kp[ikp]);
				HParticleTool::fillTLorentzVector(kplusLV,kpCand,11,true);

				Double_t flow     = kplusLV.Phi()*TMath::RadToDeg() - phiRP_tree;  // angle w.r.t. FW event plane
				if (flow<-180.) flow += 360.;
				if (flow>+180.) flow -= 360.;

				Double_t kpGammaInv = ( kpCand->getBeta()>=1. ? 0 : TMath::Sqrt(1-TMath::Power(kpCand->getBeta(),2)) );
				Double_t kpMass = ( kpCand->getCorrectedMomentumPID(11)*kpGammaInv )/kpCand->getBeta();

				sekaon.invMass = kpMass;
                                sekaon.system = kpCand->getSystemUsed();
				sekaon.theta = kplusLV.Theta()*TMath::RadToDeg();
				sekaon.phi = ( kplusLV.Phi()*TMath::RadToDeg()>0 ? kplusLV.Phi()*TMath::RadToDeg() : kplusLV.Phi()*TMath::RadToDeg()+360.) ;
				sekaon.mom = kplusLV.P();
				sekaon.momX = kplusLV.Px();
				sekaon.momY = kplusLV.Py();
				sekaon.momZ = kplusLV.Pz();
				sekaon.pt = kplusLV.Pt();
				sekaon.mt = kplusLV.Mt();
				sekaon.mtm0 = kplusLV.Mt() - HPhysicsConstants::mass(11);
				sekaon.y = kplusLV.Rapidity();
				sekaon.ycm = kplusLV.Rapidity() - 0.74;
                                sekaon.y0 = kplusLV.Rapidity()/0.74 - 1.;
				sekaon.flow = flow;

				// eff.  corrections - dependence on track multiplicity vs phi-phi_RP
//				Int_t id = 8;   //TODO k+ corrections!! , 8 -> pi+ , 14 -> proton
//				Int_t id = 14;   //TODO k+ corrections!! , 8 -> pi+ , 14 -> proton
				Int_t id = 11;   //TODO k+ corrections!! , 8 -> pi+ , 14 -> proton
				Float_t eff=0., weight=0.;

				eff=getEff(id,4-multClassNo,flow,kplusLV.Theta()*TMath::RadToDeg());
				if(eff!=0) weight=1./eff;
				if(weight>10||weight<0) weight=0.;
				sekaon.occupancyCorr = weight;

				eff=getAccRecEff(11,kplusLV.Pt(),kplusLV.Rapidity());
				if(eff!=0) weight=1./eff;
				if(weight>10||weight<0) weight=0.;
				sekaon.accreceffCorr = weight;

				sekaon.MDCtracks = -10000.;
				sekaon.METAclst = -10000.;

                                sekaon.downflag = downflag;

				if( sekaon.invMass > 300 && sekaon.invMass < 700 ){
				    sekaontree->Fill();

                                    downflag = -1.;
				}

			    }
			}//end of K+ loop


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
            out->cd("kaonPlus");

	    TMacro *mAna = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/analysis.cc");
	    TMacro *mLoop = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/loopDST_task.C");
	    TMacro *mPartReco = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/KPlusTree.h");
            sekaontree->GetUserInfo()->Add(mAna);
            sekaontree->GetUserInfo()->Add(mLoop);
            sekaontree->GetUserInfo()->Add(mPartReco);

	    sekaontree->Write();

	    for(int icc=0; icc<nCentralityClasses; icc++){
		heventplane[icc]->Write();
                hphiAB[icc]->Write();
	    }
	    hCentralityCounter->Write();

#ifndef isRealData
	    gkplustree->Write();
#endif

//	    out->Save();
//	    out->Close();
	}

	return kTRUE;
    }
   //ClassDef(KPlusTree,0) // no streamer, do not use if you include blank header in ACLIC
};

#endif
