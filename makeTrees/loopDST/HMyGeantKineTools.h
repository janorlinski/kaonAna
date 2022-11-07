#ifndef __HMYGEANTKINETOOLS__
#define __HMYGEANTKINETOOLS__

#include "mylibs.h"

using namespace std;

Int_t mygetNMdcHits(HGeantKine* gk, Int_t module) {
    // Return number of MDC hits made by present track
    // for MDC module (range: [0,3]).

    Int_t n = 0;
    if (module >= 0 && module <4)
    {
	cout << "First MDC hit: " << gk->getFirstMdcHit() << endl;
	if(gk->getFirstMdcHit() > -1)
	{
	    HGeantMdc* hit = NULL;
	    gk->resetMdcIter();
	    while((hit = (HGeantMdc*) gk->nextMdcHit()) != NULL)
	    {
		cout << "hit module: " << hit->getModule() << endl;
		if (hit->getModule() == module)
		{
		    n++;
		}
	    }
	    gk->resetMdcIter();
	    return n;
	}
    }
    else
    {
	return -1;
    }
    return n;
}


Int_t mygetNTofHits(HGeantKine* gk) {
    // Return number of TOF hits made by present track.
    //
    Int_t n = 0;
    if(gk->getFirstTofHit() > -1) {
	gk->resetTofIter();
	while(gk->nextTofHit() != NULL) n++;
	gk->resetTofIter();
    }
    return n;
}

Int_t mygetNRpcHits(HGeantKine* gk) {
    // Return number of RPC hits made by present track.
    //
    Int_t n = 0;
    if(gk->getFirstRpcHit() > -1) {
	gk->resetRpcIter();
	while(gk->nextRpcHit() != NULL) n++;
	gk->resetRpcIter();
    }
    return n;
}

#endif
