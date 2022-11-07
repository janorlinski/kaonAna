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

    static const int nMultClass = 4;
    TH1F *hMC_counter;
    TH2F *hCentrClassVsNTracks;
    TH2F *hThetaVsDphi_all_MC[nMultClass];
    
public:

    FillNtuple(const Text_t *name = "",const Text_t *title ="",TString outfile="FillNtuple.root", Int_t dummy1=0, Int_t dummy2=0, Int_t dummy3=0)
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

	    hMC_counter = new TH1F("hMC_counter","counting events in each multiplicity class",nMultClass,0,nMultClass);

            hCentrClassVsNTracks = new TH2F("hCentrClassVsNTracks","centrality and number of taken tracks",300,0,300,nMultClass,0,nMultClass);

	    for(int imc=0; imc<nMultClass; imc++){
		hThetaVsDphi_all_MC[imc] = new TH2F( Form("hThetaVsDphi_all_MC%i",imc), "theta Delta phi for one multiplicity class",180,-180,180,90,0,90);
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



//		HWallEventPlane* event_plane;
//		event_plane = HCategoryManager::getObject(event_plane,evtPlaneCat,0);

		//-------------------------------------------------
		// loop over start hits to get start_strip
//		HStart2Hit *start = 0;
//		Int_t start_strip=-1;
//		for(int istart=0; istart < fCatStartHit->getEntries(); istart++){
//		    start = HCategoryManager::getObject(start,fCatStartHit,istart);
//		    start_strip=start->getStrip();
//		    if(start->getCorrFlag()==1) start_strip=16+start->getStrip();
//		}
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

//		Bool_t isGoodEvent = kFALSE;
//		if( (start_strip >= 0) &&
//		   ((Vertex_z > -60) && (Vertex_z < 0)) &&
//       	   ((TBit&8192)==8192) &&
//		   (mult_wall >= 4) ) isGoodEvent = kTRUE;
//		   kTRUE ) isGoodEvent = kTRUE;



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
                   && ((Vertex_z > -60) && (Vertex_z < 0))
		  ) isGoodEvent = kTRUE;


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
			Int_t nCandHad     = sorter.fill(HParticleTrackSorter::selectHadrons);
//			Int_t nCandHadBest = sorter.selectBest(Particle::kIsBestRKSorter,Particle::kIsHadronSorter);
			Int_t nCandHadBest = sorter.selectBest(Particle::kIsBestHitMETASorter,Particle::kIsHadronSorter);

			HParticleCand* cand=0;

			Double_t nGoodTracks = 0.;

			Int_t size = candCat->getEntries();
			for(Int_t j = 0; j < size; j++)
			{
			    Bool_t isGoodTrack = kFALSE;
			    Float_t candPhi = 0.;
			    Float_t candTheta = 0.;

			    cand = HCategoryManager::getObject(cand,candCat,j);

			    if(cand) {

//				if(!gLoop->goodSector(cand->getSector())) continue;  // skipp inactive sectors
				if(!cand->isFlagBit(kIsUsed)) continue;
				if(cand->isAtAnyMdcEdge()) continue;
				if(cand->getSystemUsed() == -1) continue;
				if(cand->getChi2()>400) continue;

				//if(cand->getSystem() == 2) continue;

				Short_t  sector      = cand->getSector();
				Short_t  system      = cand->getSystem();
				Short_t  charge      = cand->getCharge();
				Float_t  beta        = cand->getBeta() ;
				Float_t  mom         = cand->getMomentum();
				Float_t  phi         = cand->getPhi() ;
				Float_t  theta       = cand->getTheta();
				Float_t  chi2        = cand->getChi2();
				Float_t  MetaMatchQA = cand->getMetaMatchQuality();
				Float_t  mass2       = cand->getMass2();
                                Float_t  Tof         = cand->getTof();

				Int_t flagHadron = 0.;
				if( cand->isFlagAND(4,
						     Particle::kIsAcceptedHitInnerMDC,
						     Particle::kIsAcceptedHitOuterMDC,
						     Particle::kIsAcceptedHitMETA,
						     Particle::kIsAcceptedRK
						    )
				   &&
				   cand->getInnerSegmentChi2() > 0
				  ) flagHadron = 1;

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

				//track selection
				Int_t flag_acc_mq=0;
				if(MetaMatchQA<3) flag_acc_mq=1;
				if(id==9&&MetaMatchQA<4) flag_acc_mq=1;
				if(flag_acc_mq==1 && flagHadron==1){
				    isGoodTrack = kTRUE;
				    candPhi = phi;
                                    candTheta = theta;
				}

			    }

			    if(mult_wall>4 && isGoodTrack){
				Float_t dphi= candPhi - phiRP;
				if(dphi >= 180.) dphi -= 360.;
				if(dphi <= -180.) dphi += 360.;

				if(behruzCentralityClass>=1 && behruzCentralityClass<=4){
				    hThetaVsDphi_all_MC[behruzCentralityClass-1]->Fill(dphi,candTheta);
                                    nGoodTracks += 1.;
				}
			    }

			} // end cand loop

                        hCentrClassVsNTracks->Fill(nGoodTracks-0.5,behruzCentralityClass-0.5);

			if(nGoodTracks>0) hMC_counter->Fill(behruzCentralityClass-0.5);
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

	    hMC_counter->Write();
            hCentrClassVsNTracks->Write();

            for(int imc=0; imc<nMultClass; imc++) hThetaVsDphi_all_MC[imc]->Write();

	    out->Save();
	    out->Close();
	}

	return kTRUE;
    }
   //ClassDef(FillNtuple,0) // no streamer, do not use if you include blank header in ACLIC
};

#endif
