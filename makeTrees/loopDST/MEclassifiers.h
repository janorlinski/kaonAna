#ifndef __MECLASSIFIERS__
#define __MECLASSIFIRES__

#include "hparticleevtinfo.h"
#include "hcategorymanager.h"


/*
//original version
Int_t eventClassifierMultTargetTime() {
    HParticleEvtInfo* hpei = NULL;
    hpei = HCategoryManager::getObject(hpei, catParticleEvtInfo, 0);

    // Multiplicity Index from 0 to 9
    Int_t IndexMult = (hpei->getSumSelectedParticleCandMult() / 10) - 1;

    // Vertex Index for X and Y from 0 to 9 and for Z from 0 to 14
    Int_t IndexVertexX = TMath::Max(TMath::Min(TMath::FloorNint(gHades->getCurrentEvent()->getHeader()->getVertexReco().getX() * 2 + 5.), 9), 0);
    Int_t IndexVertexY = TMath::Max(TMath::Min(TMath::FloorNint(gHades->getCurrentEvent()->getHeader()->getVertexReco().getY() * 2 + 5.), 9), 0);
    Int_t IndexVertexZ = TMath::FloorNint((gHades->getCurrentEvent()->getHeader()->getVertexReco().getZ() + 56.575) / 3.577);

    // Time Index representing the day of the collision
    Int_t IndexTime = TMath::FloorNint((gHades->getCurrentEvent()->getHeader()->getEventRunNumber() - 125368800) / 86400.) + 1;
    return 15000*IndexTime + 1500*IndexMult + 150*IndexVertexX + 15*IndexVertexY + IndexVertexZ;
}



inline Bool_t MixedEventSelectionCuts(HEventHeader* event_header, HParticleEvtInfo* event_info) {
    return (   event_header->getVertexReco().getX() >= -2.5
	    && event_header->getVertexReco().getX() <   2.5
	    && event_header->getVertexReco().getY() >= -2.5
            && event_header->getVertexReco().getY() <   2.5
	    && event_header->getVertexReco().getZ() >= -56.575
	    && event_header->getVertexReco().getZ() <   -2.92
	    && event_info->getSumSelectedParticleCandMult() >= 10
	    && event_info->getSumSelectedParticleCandMult() < 110);
}
*/

//new version
Int_t eventClassifierMultTargetTime() {
    const UInt_t   NTargetLayer                 =  15;
#if isRealData || isEmbedding
    const Double_t TargetStartPos               = -56.575;
    const Double_t TargetEndPos                 =  -2.920;
#else
    const Double_t TargetStartPos               = -10.;
    const Double_t TargetEndPos                 =  10.;
#endif

    const Int_t    NCentClasses                 = 8;
#if isRealData || isEmbedding
    const Int_t    CentClassesNToFRPCHits[]     = {173   , 101   , 76   ,  55   ,  38   , 25   ,  15  ,  9  ,   0   }; //10% classes!!
#else
    const Int_t    CentClassesNToFRPCHits[]     = {173   , 101   , 76   ,  55   ,  38   , 25   ,  15  ,  9  ,   0   }; //10% classes!!
#endif

    HParticleEvtInfo* hpei = NULL;
    hpei = HCategoryManager::getObject(hpei, catParticleEvtInfo, 0);

    // Multiplicity Index from 0 to 19
    Double_t IndexMultWidth = TMath::Ceil((CentClassesNToFRPCHits[0] - CentClassesNToFRPCHits[NCentClasses]) / 20.);
    Int_t IndexMult = TMath::FloorNint((hpei->getSumTofMultCut() + hpei->getSumRpcMultHitCut() - CentClassesNToFRPCHits[NCentClasses]) / IndexMultWidth);

    // Vertex Index for X and Y from 0 to 9 and for Z from 0 to NTargetLayer-1
    Int_t IndexVertexX = TMath::Max(TMath::Min(TMath::FloorNint(gHades->getCurrentEvent()->getHeader()->getVertexReco().getX() * 2 + 5.), 9), 0);
    Int_t IndexVertexY = TMath::Max(TMath::Min(TMath::FloorNint(gHades->getCurrentEvent()->getHeader()->getVertexReco().getY() * 2 + 5.), 9), 0);
    Int_t IndexVertexZ = TMath::FloorNint((gHades->getCurrentEvent()->getHeader()->getVertexReco().getZ() - TargetStartPos) / ((TargetEndPos - TargetStartPos) / NTargetLayer));

    // Time Index representing the day of the collision
    Int_t IndexTime = TMath::FloorNint((gHades->getCurrentEvent()->getHeader()->getEventRunNumber()) / 86400.);
    return 2000*NTargetLayer*IndexTime + 100*NTargetLayer*IndexMult + 10*NTargetLayer*IndexVertexX + NTargetLayer*IndexVertexY + IndexVertexZ;
}

inline Bool_t MixedEventSelectionCuts(HEventHeader* event_header, HParticleEvtInfo* event_info) {
#if isRealData || isEmbedding
    return (   event_header->getVertexReco().getX() >= -2.5
	    && event_header->getVertexReco().getX() <   2.5
	    && event_header->getVertexReco().getY() >= -2.5
            && event_header->getVertexReco().getY() <   2.5
	    && event_header->getVertexReco().getZ() >= -56.575
	    && event_header->getVertexReco().getZ() <   -2.92
	    && (event_info->getSumTofMultCut() + event_info->getSumRpcMultHitCut()) >= 10
	    && (event_info->getSumTofMultCut() + event_info->getSumRpcMultHitCut()) < 312);

#else
    return (   event_header->getVertexReco().getX() >= -5.
	    && event_header->getVertexReco().getX() <   5.
	    && event_header->getVertexReco().getY() >= -5.
            && event_header->getVertexReco().getY() <   5.
	    && event_header->getVertexReco().getZ() >= -100.
	    && event_header->getVertexReco().getZ() <   10.
	    && (event_info->getSumTofMultCut() + event_info->getSumRpcMultHitCut()) >= 15
	    && (event_info->getSumTofMultCut() + event_info->getSumRpcMultHitCut()) < 252);

#endif
}

#endif
