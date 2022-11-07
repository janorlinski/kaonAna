#ifndef __MVATREES___     // good style to avaoid multiple includes
#define __MVATREES___

#include "mylibs.h"
#include "mystruc.h"

#include "GlobVars.h"
#include "UsefulFunctions.h"
#include "MEclassifiers.h"

#define Embedding 0
#define MixEvent 1

#define useBananaCuts 0

using namespace std;

class MVATrees : public HReconstructor
{
protected:
    // put all vars here which are not
    // local to a function

    HParticleTrackSorter sorter;
    //HParticleT0Reco t0Reco("apr12");

    static const Int_t particleID = 16;

    Float_t d1=0., d2=0., d3=1000., dVer=0., dMin=1000., kaonMom=-1.;

#if MixEvent>0
    TTree *mixeventTree;
#elif Embedding>0
    TTree *embTree;
#endif
    
    //------------------------------------
    // file handling
#if useBananaCuts>0
    TFile* bananacutsFile;
    TCutG *betamom_2sig_pim_tof_pionCmom, *betamom_2sig_pim_rpc_pionCmom,
	  *betamom_2sig_pip_tof_pionCmom, *betamom_2sig_pip_rpc_pionCmom,
          *betamom_2sig_p_tof_pionCmom, *betamom_2sig_p_rpc_pionCmom;
#endif

public:

    MVATrees(const Text_t *name = "",const Text_t *title ="",TString outfile="MVATrees.root")
	: HReconstructor(name,title)
    {  // init your vars

    }

    virtual ~MVATrees()
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

#if MixEvent>0
	    mixeventTree = new TTree("tMix","Mix Event Kaon Tree");
	    mixeventTree->Branch("d1",      &d1,      "d1/F",     32000);
	    mixeventTree->Branch("d2",      &d2,      "d2/F",     32000);
	    mixeventTree->Branch("d3",      &d3,      "d3/F",     32000);
	    mixeventTree->Branch("dVer",    &dVer,    "dVer/F",   32000);
	    mixeventTree->Branch("dMin",    &dMin,    "dMin/F",   32000);
	    mixeventTree->Branch("kaonMom", &kaonMom, "kaonMom/F",32000);

#elif Embedding>0
	    embTree = new TTree("tEmb","Embedded Kaon Tree");
	    embTree->Branch("d1",      &d1,      "d1/F",     32000);
	    embTree->Branch("d2",      &d2,      "d2/F",     32000);
	    embTree->Branch("d3",      &d3,      "d3/F",     32000);
	    embTree->Branch("dVer",    &dVer,    "dVer/F",   32000);
	    embTree->Branch("dMin",    &dMin,    "dMin/F",   32000);
	    embTree->Branch("kaonMom", &kaonMom, "kaonMom/F",32000);
#endif

#if useBananaCuts>0
	    const char *bananacutsFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root";
	    bananacutsFile = TFile::Open( bananacutsFileName );
	    betamom_2sig_pim_tof_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutPiM_TOF_2.0");
	    betamom_2sig_pim_rpc_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutPiM_RPC_2.0");
	    betamom_2sig_pip_tof_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutPiP_TOF_2.0");
	    betamom_2sig_pip_rpc_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutPiP_RPC_2.0");
	    betamom_2sig_p_tof_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutProton_TOF_2.0");
	    betamom_2sig_p_rpc_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutProton_RPC_2.0");
#endif
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
#if useBananaCuts>0
	const char *bananacutsFileName = "/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/rootFiles/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root";
	bananacutsFile = TFile::Open( bananacutsFileName );
	betamom_2sig_pim_tof_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutPiM_TOF_2.0");
	betamom_2sig_pim_rpc_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutPiM_RPC_2.0");
	betamom_2sig_pip_tof_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutPiP_TOF_2.0");
	betamom_2sig_pip_rpc_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutPiP_RPC_2.0");
	betamom_2sig_p_tof_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutProton_TOF_2.0");
	betamom_2sig_p_rpc_pionCmom=(TCutG*)bananacutsFile->Get("BetaCutProton_RPC_2.0");
#endif

	return kTRUE;
    }

    Int_t execute()
    {   // this function is called once per event.
	// if the function returns kSkipEvent the event
	// will be skipped a. for all following tasks
	// and b. not appear in output (hades eventLoop())

	if(gHades){
	    HCategory* evtInfoCat = (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo);
	    HEventHeader* header  = gHades->getCurrentEvent()->getHeader();
	    HVertex vertex        = header->getVertexReco();

            HCategory* candCat = (HCategory*)HCategoryManager::getCategory(catParticleCand);

#if Embedding>0
	    HGeantKine*       kine;
	    HLinearCategory* kineCat = (HLinearCategory*)HCategoryManager::getCategory(catGeantKine);
#endif

	    if(evtInfoCat)
	    {
		HParticleEvtInfo* evtInfo=0;
		evtInfo = HCategoryManager::getObject(evtInfo,evtInfoCat,0 );

		Bool_t isGoodEvent = kFALSE;
#if MixEvent>0
		eventmixer.setUseLeptons(kFALSE);
		eventmixer.setPIDs(8,9,particleID);  // pi-,pi+,k0s  which PIDs and MotherID are stored
		eventmixer.setBuffSize(100);   // size of buffer has to chosen according to stat   //TODO how big is needed??
                eventmixer.setWaitForBuffer(kFALSE);
		eventmixer.setEventClassifier(eventClassifierMultTargetTime);

		if(evtInfo && evtInfo->isGoodEvent(Particle::kGoodTRIGGER|
						   Particle::kGoodVertexClust|
						   Particle::kGoodVertexCand|
						   Particle::kGoodSTART|
						   Particle::kNoPileUpSTART|
						   Particle::kNoVETO|
						   Particle::kGoodSTARTVETO|
						   Particle::kGoodSTARTMETA
						  )
                   && MixedEventSelectionCuts(header,evtInfo) //new condition - from Simon!!
		  ) isGoodEvent = kTRUE;
#elif Embedding>0
		if(evtInfo) isGoodEvent = kTRUE;
#endif

                if(isGoodEvent)
		{


#if Embedding>0
       //---------------------------------- Kine --------------------------------------------------------
       for(Int_t j = 0; j < kineCat->getEntries(); j ++) {
	   kine = HCategoryManager::getObject(kine,kineCat,j);
	   if( kine->isPrimary() && kine->getID() == particleID && kine->getGeneratorInfo() == (500+particleID) ){

	       kaonMom = kine->getTotalMomentum();

               Float_t momX, momY, momZ;
	       kine->getMomentum(momX,momY,momZ);

	       vector<HGeantKine *> kaonDaughters;
	       Int_t numberOfDaughters = kine->getDaughters(kine,kaonDaughters);

	       HGeantKine *pipGK = 0;
	       HGeantKine *pimGK = 0;

	       Bool_t pipAcc = kFALSE;
	       Bool_t pimAcc = kFALSE;
	       Bool_t pipRec = kFALSE;
	       Bool_t pimRec = kFALSE;

	       HGeomVector pipBase, pipDir, pimBase, pimDir, decayVertex, eventVertex, motherDir, decayVertexGeantPIP, decayVertexGeantPIM;
	       eventVertex.setXYZ(vertex.getX(),vertex.getY(),vertex.getZ());
               motherDir.setXYZ(momX,momY,momZ);

	       for(int id=0; id<numberOfDaughters; id++){
                   if( kaonDaughters[id]->getID() == 8 ) pipGK = kaonDaughters[id];
                   if( kaonDaughters[id]->getID() == 9 ) pimGK = kaonDaughters[id];
	       }

	       if( (pipGK != 0) && (pimGK != 0) ){

//                 HGeantKine *dPip = 0; dPip = pipGK->getChargedDecayDaughter(pipGK);
//		   HGeantKine *dPim = 0; dPim = pimGK->getChargedDecayDaughter(pimGK);

		   //acceptance test: isInAcceptance(numbers)!!

                   pipGK->getVertex(decayVertexGeantPIP);
		   pimGK->getVertex(decayVertexGeantPIM);
                   if( decayVertexGeantPIP != decayVertexGeantPIM ) cout << "WARNING!!! Decay vertices from Geant are not the same!!!" << endl;

		   Bool_t test_pip_acc = kTRUE;
		   Bool_t test_pim_acc = kTRUE;

		   if( pipGK->isAtAnyMdcEdge(2) ) test_pip_acc = kFALSE;
		   if( pimGK->isAtAnyMdcEdge(2) ) test_pim_acc = kFALSE;
		   for(int imm=0; imm<4; imm++){
		       if( pipGK->getNLayerMod(imm) < 4 ) test_pip_acc = kFALSE;
		       if( pimGK->getNLayerMod(imm) < 4 ) test_pim_acc = kFALSE;
		   }

                   if( !( pipGK->getSys(0) || pipGK->getSys(1) ) ) test_pip_acc = kFALSE;
		   if( !( pimGK->getSys(0) || pimGK->getSys(1) ) ) test_pim_acc = kFALSE;

		   pipAcc = test_pip_acc;
		   pimAcc = test_pim_acc;


		   HParticleCandSim* cand=0;
                   if(candCat){
		       Int_t size = candCat->getEntries();

                       cout << size << endl;

		       for(Int_t j = 0; j < size; j++)
		       {
			   cand = HCategoryManager::getObject(cand,candCat,j);

			   Bool_t test=kFALSE;
			   if( cand->isFlagAND(4,
						Particle::kIsAcceptedHitInnerMDC, //inner chi2 >0
						Particle::kIsAcceptedHitOuterMDC, //outer chi2 >0
						Particle::kIsAcceptedHitMETA,     //has meta hit
						Particle::kIsAcceptedRK           //rk chi2>0
					       )
			      &&
			      cand->getMetaMatchQuality() < 3
			      &&
			      cand->getChi2()             < 400
			     ) test=kTRUE;
			   if(!test) continue;

			   if(cand->isAtAnyMdcEdge()) continue;
			   if(cand->getSystemUsed() == -1) continue;

			   if( cand->getGeantTrack() == pipGK->getTrack() ){
			       pipRec = kTRUE;
			       CalcSegVector(cand->getZ(), cand->getR(), cand->getPhi()*TMath::DegToRad(), cand->getTheta()*TMath::DegToRad(), pipBase, pipDir);
			   }

			   if( cand->getGeantTrack() == pimGK->getTrack() ){
			       pimRec = kTRUE;
			       CalcSegVector(cand->getZ(), cand->getR(), cand->getPhi()*TMath::DegToRad(), cand->getTheta()*TMath::DegToRad(), pimBase, pimDir);
			   }
		       }
		   }
	       }

	       if( pipAcc && pipRec && pimAcc && pimRec ){
		   decayVertex = calcVertexAnalytical(pipBase, pipDir, pimBase, pimDir);
		   if(   TMath::Abs( decayVertex.getX() - decayVertexGeantPIM.getX() ) > 10
		      || TMath::Abs( decayVertex.getY() - decayVertexGeantPIM.getY() ) > 10
		      || TMath::Abs( decayVertex.getZ() - decayVertexGeantPIM.getZ() ) > 10 ) cout << "WARNING!! Calculated decay vertex is far from Geant decay vertex!!" << endl;


		   dVer = TMath::Sqrt( (decayVertex.getX()-eventVertex.getX())*(decayVertex.getX()-eventVertex.getX()) +
					  (decayVertex.getY()-eventVertex.getY())*(decayVertex.getY()-eventVertex.getY()) +
					  (decayVertex.getZ()-eventVertex.getZ())*(decayVertex.getZ()-eventVertex.getZ()) );
		   d1 = calculateMinimumDistanceStraightToPoint(pimBase,     pimDir,    eventVertex);
		   d2 = calculateMinimumDistanceStraightToPoint(pipBase,     pipDir,    eventVertex);
		   d3 = calculateMinimumDistanceStraightToPoint(decayVertex, motherDir, eventVertex);
		   dMin = calculateMinimumDistance(pimBase, pimDir, pipBase, pipDir);
	       }

               embTree->Fill();
	   }
       }
#endif

#if MixEvent>0
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
	   Int_t nCandHadBest = sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsHadron);

	   if( nCandHad<0 || nCandHadBest<0 ) return 1;

	   //t0Reco.execute();

	   vector<HParticleCand *> vpim;
	   vector<HParticleCand *> vpip;
	   HParticleCand* cand=0;

	   Int_t size = candCat->getEntries();
	   for(Int_t j = 0; j < size; j++)
	   {
	       cand = HCategoryManager::getObject(cand,candCat,j);

	       if(cand) {

		   if(!gLoop->goodSector(cand->getSector())) { continue;}  // skipp inactive sectors
		   if(!cand->isFlagBit(kIsUsed)) continue;
		   if(cand->isAtAnyMdcEdge()) continue;
		   if(cand->getSystemUsed() == -1) continue;
		   if(cand->getChi2() > 400) continue;
		   if(cand->getMetaMatchQuality() > 3) continue;

#if useBananaCuts>0
		   if( (cand->getSystemUsed()==1 && betamom_2sig_pip_tof_pionCmom->IsInside(cand->getCharge()*cand->getMomentum(),cand->getBeta()))
		      || (cand->getSystemUsed()==0 && betamom_2sig_pip_rpc_pionCmom->IsInside(cand->getCharge()*cand->getMomentum(),cand->getBeta())) )
		       vpip.push_back(new HParticleCand(*cand));


		   if( (cand->getSystemUsed()==1 && betamom_2sig_pim_tof_pionCmom->IsInside(cand->getCharge()*cand->getMomentum(),cand->getBeta()))
		      || (cand->getSystemUsed()==0 && betamom_2sig_pim_rpc_pionCmom->IsInside(cand->getCharge()*cand->getMomentum(),cand->getBeta())) )
		       vpim.push_back(new HParticleCand(*cand));
#else
		   if( cand->getCharge()>0 && cand->getMass()<300 && cand->getMass()>0 && cand->getMomentum()>0 && cand->getMomentum()<1000 ) vpip.push_back(new HParticleCand(*cand));
		   if( cand->getCharge()<0 && cand->getMass()<300 && cand->getMass()>0 && cand->getMomentum()>0 && cand->getMomentum()<1000 ) vpim.push_back(new HParticleCand(*cand));
//		   if( cand->getCharge()>0 && cand->getCorrectedMass2PID(8)<300 && cand->getCorrectedMass2PID(8)>0 ) vpip.push_back(new HParticleCand(*cand));
//		   if( cand->getCharge()<0 && cand->getCorrectedMass2PID(9)<300 && cand->getCorrectedMass2PID(9)>0 ) vpim.push_back(new HParticleCand(*cand));
#endif

	       }
	   } // end cand loop

	   eventmixer.nextEvent();
	   eventmixer.addVector(vpim,9);
	   eventmixer.addVector(vpip,8);

	   HVirtualCand* pipVCand;
	   HVirtualCand* pimVCand;
	   TLorentzVector pipLV(1.,1.,1.,1.), pimLV(1.,1.,1.,1.), sumLV(1.,1.,1.,1.);

	   HGeomVector pipBase, pipDir, pimBase, pimDir, decayVertex, eventVertex, motherDir;
	   eventVertex.setXYZ(vertex.getX(),vertex.getY(),vertex.getZ());


	   vector<HParticlePair>& pairsVec = eventmixer.getMixedVector();

	   Int_t eventClass = eventmixer.currentEventClass();
	   HGeomVector eventVertexME(TMath::Floor((eventClass % 1500) / 150.) / 2. - 2.25,
				     TMath::Floor((eventClass %  150) /  15.) / 2. - 2.25,
				     ((eventClass % 15) * 3.577) - 54.7865);

	   for(UInt_t j = 0;j < pairsVec.size(); j++){
	       HParticlePair& pair = pairsVec[j];
	       //pair.setDoMomentumCorrection(kFALSE);

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

	       HParticleTool::fillTLorentzVector(pipLV,pipVCand,8,true);
	       HParticleTool::fillTLorentzVector(pimLV,pimVCand,9,true);
	       sumLV = pipLV + pimLV;

	       //calc topology variables
	       CalcSegVector(pipVCand->getZ(), pipVCand->getR(), pipVCand->getPhi()*TMath::DegToRad(), pipVCand->getTheta()*TMath::DegToRad(), pipBase, pipDir);
	       CalcSegVector(pimVCand->getZ(), pimVCand->getR(), pimVCand->getPhi()*TMath::DegToRad(), pimVCand->getTheta()*TMath::DegToRad(), pimBase, pimDir);

	       motherDir.setXYZ(sumLV.X(),sumLV.Y(),sumLV.Z());
	       decayVertex = calcVertexAnalytical(pipBase, pipDir, pimBase, pimDir);

	       dVer = TMath::Sqrt( (decayVertex.getX()-eventVertexME.getX())*(decayVertex.getX()-eventVertexME.getX()) +
				      (decayVertex.getY()-eventVertexME.getY())*(decayVertex.getY()-eventVertexME.getY()) +
				      (decayVertex.getZ()-eventVertexME.getZ())*(decayVertex.getZ()-eventVertexME.getZ()) );
	       d1 = calculateMinimumDistanceStraightToPoint(pimBase,     pimDir,    eventVertexME);
	       d2 = calculateMinimumDistanceStraightToPoint(pipBase,     pipDir,    eventVertexME);
	       d3 = calculateMinimumDistanceStraightToPoint(decayVertex, motherDir, eventVertexME);
	       dMin = calculateMinimumDistance(pimBase, pimDir, pipBase, pipDir);

	       kaonMom = sumLV.P();

	       mixeventTree->Fill();

	   }

	   eventmixer.RemoveObjectsToDelete();


       } // end cand cat
       //-------------------------------------------------
#endif

		} // end of GoodEvent selection
		//-------------------------------------------------

	    } // end evtInfoCat
	} // end gHades

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

	    TMacro *mAna = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/analysis.cc");
	    TMacro *mLoop = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/loopDST_task.C");
	    TMacro *mPartReco = new TMacro("/lustre/nyx/hades/user/lchlad/flow/scripts/apr12/loopDST/MVATrees.h");

#if MixEvent>0
	    mixeventTree->GetUserInfo()->Add(mAna);
	    mixeventTree->GetUserInfo()->Add(mLoop);
	    mixeventTree->GetUserInfo()->Add(mPartReco);

	    mixeventTree->Write();
#elif Embedding>0
	    embTree->GetUserInfo()->Add(mAna);
	    embTree->GetUserInfo()->Add(mLoop);
	    embTree->GetUserInfo()->Add(mPartReco);

	    embTree->Write();
#endif

	}

	return kTRUE;
    }
   //ClassDef(MVATrees,0) // no streamer, do not use if you include blank header in ACLIC
};

#endif
