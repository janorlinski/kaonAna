#ifndef __FILLNTUPLE___     // good style to avaoid multiple includes
#define __FILLNTUPLE___

#include "mylibs.h"
#include "mystruc.h"

#include "GlobVars.h"

#include "UsefulFunctions.h"


using namespace std;

class FillNtuple : public HReconstructor
{
protected:
    HParticleTrackSorter sorter;

    BRANCH pid;
    TTree *nt;

    static const int nMultClass = 4;
    TH1F *hMC_counter;
    TH2F *hThetaVsDphi_all_MC[nMultClass];
    
public:

    FillNtuple(const Text_t *name = "",const Text_t *title ="",TString outfile="FillNtuple.root")
	: HReconstructor(name,title)
    {  // init your vars

    }

    virtual ~FillNtuple()
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

	    nt = new TTree("nt","pid");
	    nt->Branch("pid",&pid,"eventnumber:downflag:trigbit:"
		                  "mult_tofrpc:mult_class:start_strip:system:sector:"
		                  "tof:theta:phi:mom:beta:chi2:charge:mass:isused:"
		                  "metamatch_quality:ishadron:zprime:id:"
		                  "mult_wall:RP:RPAB");


	    hMC_counter = new TH1F("hMC_counter","counting events in each multiplicity class",nMultClass,0,nMultClass);

	    for(int imc=0; imc<nMultClass; imc++){
		hThetaVsDphi_all_MC[imc] = new TH2F( Form("hThetaVsDphi_all_MC%i",imc), "theta Delta phi for one multiplicity class",180,-180,180,100,0,100);
                hThetaVsDphi_all_MC[imc]->Sumw2();
	    }
	}

	return kTRUE;
    }
    Bool_t reinit()
    {   // this function is called for each file
        // after init()
	// use this place if you need already initialized
	// containers


	//--------------------------CONFIGURATION---------------------------------------------------
	//At begin of the program (outside the event loop)
	//sorter.setDebug();                                            // for debug
	//sorter.setPrintLevel(3);                                      // max prints
	//sorter.setRICHMatching(HParticleTrackSorter::kUseRKRICHWindow,4.); // select matching RICH-MDC for selectLeptons() function
	sorter.setUseYMatching(kTRUE);                                // use meta cell build in select functions
	//sorter.setMetaBoundary(4.0);                                   // match in 4 mm to meta cell
	//sorter.setUseBeta(kTRUE);                                     // require beta >0 in build in select functions
	//sorter.setBetaLeptonCut(0.9);                                 // lower beta cut for lepton select
	//sorter.setUseFakeRejection(kTRUE);                            // reject fakes in  build in select functions (default kTRUE)
	//sorter.setIgnoreInnerMDC();                                   // do not reject Double_t inner MDC hits
	//sorter.setIgnoreOuterMDC();                                   // do not reject Double_t outer MDC hits
	//sorter.setIgnoreMETA();                                       // do not reject Double_t META hits
	//sorter.setIgnorePreviousIndex();                              // do not reject indices from previous selctions
	sorter.init();                                                  // get catgegory pointers etc...

	evtChara.reinit();
	//--------------------------------------------------------------------------------------------

	return kTRUE;
    }

    Int_t execute()
    {   // this function is called once per event.
	// if the function returns kSkipEvent the event
	// will be skipped a. for all following tasks
	// and b. not appear in output (hades eventLoop())

	if(gHades){

	    HEventHeader* header  = gHades->getCurrentEvent()->getHeader();

	    Int_t TBit   = (Int_t) header->getTBit();
	    Int_t DF     = (Int_t) header->getDownscalingFlag();
	    Int_t SeqNum = (Int_t) header->getEventSeqNumber();
	    Int_t TDec   = (Int_t) header->getTriggerDecision();
	    Float_t Vertex_z = header->getVertexZ();

	    HCategory* evtInfoCat   = (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo);
	    HCategory* candCat      = (HCategory*)HCategoryManager::getCategory(catParticleCand);
	    HCategory* evtPlaneCat  = (HCategory*)HCategoryManager::getCategory(catWallEventPlane);
	    HCategory* fCatWallHit  = (HCategory*)HCategoryManager::getCategory(catWallHit);
	    HCategory* fCatStartHit = (HCategory*)HCategoryManager::getCategory(catStart2Hit);

	    if(evtInfoCat)
	    {

		HParticleEvtInfo* evtInfo=0;
		evtInfo = HCategoryManager::getObject(evtInfo,evtInfoCat,0 );

		Double_t mult = evtInfo->getSumTofMultCut()+evtInfo->getSumRpcMultHitCut();



		HWallEventPlane* event_plane;
		event_plane = HCategoryManager::getObject(event_plane,evtPlaneCat,0);

		//-------------------------------------------------
		// loop over start hits to get start_strip
		HStart2Hit *start = 0;
		Int_t start_strip=-1;
		for(int istart=0; istart < fCatStartHit->getEntries(); istart++){
		    start = HCategoryManager::getObject(start,fCatStartHit,istart);
		    start_strip=start->getStrip();
		    if(start->getCorrFlag()==1) start_strip=16+start->getStrip();
		}
		//-------------------------------------------------

		//-------------------------------------------------
		// loop over wall hits to get multiplicity in ForwardWall
		HWallHit *wall = 0;
		Float_t mult_wall=0;
		for(int iwall=0; iwall < fCatWallHit->getEntries(); iwall++){
		    wall = HCategoryManager::getObject(wall,fCatWallHit,iwall);
		    if(wall->getTime()>22&&wall->getTime()<30) {
			if(wall->getCell()<144&&wall->getCharge()>83)
			    mult_wall++;

			else if(wall->getCell()>=144&&wall->getCell()<208&&wall->getCharge()>84)
			    mult_wall++;

			else if(wall->getCell()>=208&&wall->getCharge()>88)
			    mult_wall++;

		    }

		}
		//-------------------------------------------------

		//-------------------------------------------------
		// using Behruz event characteristics to get the Reaction Plane angle and resolution
		Int_t behruzCentralityClass = 0;
		behruzCentralityClass = evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10);
		Float_t phiRP = evtChara.getEventPlane(HParticleEvtChara::kDefault)*TMath::RadToDeg();

		if(phiRP > 180.) phiRP -= 360.;
		//-------------------------------------------------

		Bool_t isGoodEvent = kFALSE;
		if( (start_strip >= 0) &&
		   ((Vertex_z > -60) && (Vertex_z < 0)) &&
		   ((TBit&8192)==8192) &&
//		   (mult_wall >= 4) ) isGoodEvent = kTRUE;
		   kTRUE ) isGoodEvent = kTRUE;



                /*
		Bool_t isGoodEvent = kFALSE;
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
		*/



		if(isGoodEvent)
		{
		    //-------------------------------------------------
		    // loop over particle candidates in event
		    if(candCat){
			TString filename;

			//------------------------------------------------------------------------
			// prepare track sorting
			// clean vectors and index arrays
			sorter.cleanUp();
			sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);
			Int_t nCandHad     = sorter.fill(selectNegativeHadrons);
//			Int_t nCandHad     = sorter.fill(HParticleTrackSorter::selectHadrons);
			Int_t nCandHadBest = sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsHadron);

			HParticleCand* cand1=0;

			Int_t size = candCat->getEntries();

			Int_t N_PC = 0;

//                        Bool_t isGoodTrack = kTRUE;
                        Bool_t isGoodTrack = kFALSE;

			for(Int_t j = 0; j < size; j++)
			{
			    cand1 = HCategoryManager::getObject(cand1,candCat,j);

                            clearPID(pid);

			    if(cand1) {
				Short_t  sector      = cand1->getSector();
				Short_t  system      = cand1->getSystem();
				Short_t  charge      = cand1->getCharge();
				Float_t  beta        = cand1->getBeta() ;
				Float_t  mom         = cand1->getMomentum();
				Float_t  phi         = cand1->getPhi() ;
				Float_t  theta       = cand1->getTheta();
				Float_t  chi2        = cand1->getChi2();
				Float_t  MetaMatchQA = cand1->getMetaMatchQuality();
				Float_t  mass2       = cand1->getMass2();
                                Float_t  Tof         = cand1->getTof();

				Float_t flagHadron = 0.;
				if( cand1->isFlagAND(4,
						     Particle::kIsAcceptedHitInnerMDC,
						     Particle::kIsAcceptedHitOuterMDC,
						     Particle::kIsAcceptedHitMETA,
						     Particle::kIsAcceptedRK
						    )
				   &&
				   cand1->getInnerSegmentChi2() > 0
				  ) flagHadron = 1.;

				Float_t flagIsUsed = 0.;
				if (cand1->isFlagBit(Particle::kIsUsed)) flagIsUsed = 1.;

				//PID
				// hard cut: beta vs mom
				Int_t id=0;
				Float_t beta_calc=mom/TMath::Sqrt(mom*mom+139.*139.);
				Float_t beta_calc_p=mom/TMath::Sqrt(mom*mom+938.*938.);

				if(beta>0) {
				    if(system==0) {
					if(charge>0&&TMath::Abs(7./beta-7./beta_calc)<0.4&&mom<1200) id=8;
					if(charge<0&&TMath::Abs(7./beta-7./beta_calc)<0.4) id=9;
				    }
				    if(system>=1) {
					if(charge>0&&TMath::Abs(7./beta-7./beta_calc)<0.6&&mom<1000) id=8;
					if(charge<0&&TMath::Abs(7./beta-7./beta_calc)<0.6) id=9;
				    }

				    if(charge>0&&TMath::Abs(beta-beta_calc_p)<0.06) id=14;
				}

				N_PC++;
				pid.downflag =0;
				if(N_PC > 1) pid.downflag = -1; //only one pid per event has downlag > -1, important for next step

				pid.eventnumber       = SeqNum;
				pid.trigbit           = TBit;
				pid.zprime            = Vertex_z;
				pid.mult_tofrpc       = mult;
				pid.RP                = phiRP;
				pid.mult_class        = behruzCentralityClass;
				pid.start_strip       = start_strip;
				pid.mult_wall         = mult_wall;
				pid.system            = system;
				pid.sector            = sector;
				pid.tof               = Tof;
				pid.theta             = theta;
				pid.phi               = phi;
				pid.mom               = mom;
				pid.beta              = beta;
				pid.chi2              = chi2;
				pid.charge            = charge;
				pid.mass              = mass2;
				pid.metamatch_quality = MetaMatchQA;
				pid.ishadron          = flagHadron;
				pid.isused            = flagIsUsed;
				pid.id                = id;

			    }

			    //track selection
                            /*
			    Int_t flag_acc_mq=0;
			    if(pid.metamatch_quality<3) flag_acc_mq=1;
			    if(pid.id==9&&pid.metamatch_quality<4) flag_acc_mq=1;
			    if(flag_acc_mq==1&&pid.chi2<400&&pid.ishadron==1&&pid.system>-1&&pid.isused==1) isGoodTrack = kTRUE;
                            */
			    //track selection
			    Int_t flag_acc_mq=0;
			    if(pid.metamatch_quality<1.5) flag_acc_mq=1;
			    if(pid.id==9&&pid.metamatch_quality<2) flag_acc_mq=1;
			    if(flag_acc_mq==1&&pid.chi2<100&&pid.ishadron==1&&pid.system>-1&&pid.system<2&&pid.isused==1) isGoodTrack = kTRUE;

//			    if(isGoodTrack) nt->Fill();

//			    if(N_PC>=1 && behruzCentralityClass>=1 && behruzCentralityClass<=4) hMC_counter->Fill(behruzCentralityClass-0.5);

			    if(mult_wall>4 && isGoodTrack){
				Float_t dphi=pid.phi-pid.RP;
				if(dphi >= 180.) dphi -= 360.;
				if(dphi <= -180.) dphi += 360.;

				if(behruzCentralityClass>=1 && behruzCentralityClass<=4){
				    hThetaVsDphi_all_MC[behruzCentralityClass-1]->Fill(dphi,pid.theta);
				    if(N_PC>1) hMC_counter->Fill(behruzCentralityClass-0.5);
				}
			    }

			} // end cand loop

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

	if(out) {
	    out->cd();

	    //nt->Write();
	    hMC_counter->Write();
            for(int imc=0; imc<nMultClass; imc++) hThetaVsDphi_all_MC[imc]->Write();

	    out->Save();
	    out->Close();
	}

	return kTRUE;
    }
   //ClassDef(FillNtuple,0) // no streamer, do not use if you include blank header in ACLIC
};

#endif
