#include "mylibs.h"
#include "mystruc.h"

#include "GlobVars.h"

//#include "KChargedTrees.h"
//#include "KPlusTree.h"
//#include "PionTrees.h"
//#include "ProtonTree.h"
//#include "FillNtuple.h"

#include "KZeroTrees.h"
//#include "MVATrees.h"

//#include "CUTSEARCHreco.h"

using namespace std;

Int_t loopDST_task(TString infileList="/lustre/nyx/hades/dst/mar19/gen5/ag158ag/3200A/058/root/be1905800144901.hld_dst_mar19_accepted.root",
		   TString outfile="test.root",
//		   Int_t nEvents=1000,
		   Int_t myInt1=-1,
		   Int_t myInt2=-1,
		   Int_t myInt3=-1
		  )
{
    //  infileList : comma seprated file list "file1.root,file2.root" or "something*.root"
    //  outfile    : optional (not used here) , used to store hists in root file
    //  nEvents    : number of events to processed. if  nEvents < entries or < 0 the chain will be processed

    //-------------------------------------------------
    // create loop obejct and hades
    HLoop loop(kTRUE);
    //-------------------------------------------------

    //-------------------------------------------------
    // setup the spectrometer for filter tasks
    // this part is needed if your task needs
    // detectors or parameter container
    if(0){
	TString beamtime="mar19";

	Int_t mdcMods[6][4]=
	{ {1,1,1,1},
	{1,1,1,1},
	{1,1,1,1},
	{1,1,1,1},
	{1,1,1,1},
	{1,1,1,1} };
	TString asciiParFile     = "";
	TString rootParFile      = "/cvmfs/hades.gsi.de/param/real/mar19/allParam_Mar19_gen5_13042021.root";
	TString paramSource      = "root"; // root, ascii, oracle
	TString paramrelease     = "Mar19_dst_gen5";  // now for oracle
	HDst::setupSpectrometer(beamtime,mdcMods,"rich,mdc,tof,rpc,emc,wall,start,tbox");
	HDst::setupParameterSources(paramSource,asciiParFile,rootParFile,paramrelease);
    }
    //-------------------------------------------------

    //-------------------------------------------------
    // list of all files with working sectors
#if isRealData || isEmbedding
    loop.readSectorFileList("/lustre/nyx/hades/dst/mar19/gen5/sector_selection/Mar19AgAg1580_Gen5.list",kFALSE,kFALSE);  //TODO for gen9 there is no such a list!!
#endif
    //-------------------------------------------------


    //-------------------------------------------------
    // reading input files and decalring containers
    Bool_t ret =kFALSE;
    if(infileList.Contains(",")){
	ret = loop.addMultFiles(infileList);      // file1,file2,file3
    } else if(infileList.Contains(".txt") || infileList.Contains(".list")){
        ret = loop.addFilesList(infileList); //in case of file-list (1 file per line with full path)
    } else{
	ret = loop.addFiles(infileList); // myroot*.root
    }

    if(ret == 0) {
	cout<<"READBACK: ERROR : cannot find inputfiles : "<<infileList.Data()<<endl;
	return 1;
    }

#if isRealData && !isExtendedOutput
    if(!loop.setInput("-*,+HParticleCand,+HParticleEvtInfo,+HWallHit")) { //select some cotegories
	cout<<"READBACK: ERROR : cannot read input !"<<endl;
	exit(1);
    }
#elif isRealData && isExtendedOutput
    if(!loop.setInput("-*,+HParticleCand,+HParticleEvtInfo,+HWallHit,+HTofCluster")) { //select some cotegories
	cout<<"READBACK: ERROR : cannot read input !"<<endl;
	exit(1);
    }
#elif !isRealData
    if(!loop.setInput("-*,+HParticleCandSim,+HParticleEvtInfo,+HWallHit,+HGeantKine")) { //select some cotegories
	cout<<"READBACK: ERROR : cannot read input !"<<endl;
	exit(1);
    }
#endif

    //-------------------------------------------------
    // Fix lustre many seek/read issue by reading larger
    // file junks into cache. Cache all branches and
    // switch of the learning phase
    loop.getChain()->SetCacheSize(-1); // 8Mb - 8000000
    //loop.getChain()->AddBranchToCache("-*,+HParticleCand,+HParticleEvtInfo,+HRpcCluster,+HTofCluster,+HStart2Hit,+HStart2Cal,+HWallHit,+HWallEventPlane",kTRUE);
    loop.getChain()->AddBranchToCache("*",kTRUE);
    loop.getChain()->StopCacheLearningPhase();
    //-------------------------------------------------

    loop.printCategories();
    loop.printChain();
    //-------------------------------------------------

    //-------------------------------------------------
    //parameters
    dEdxCorr.setDefaultPar("mar19");     //TODO use it for K0s??


    TString EventCharaParFile;
    //due to the overlap of day126  there is dedidacate param-file for the reverse-field runs
    //if(isRevFieldDATA)   EventCharaParFile = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_apr12_gen8_revfield_2019_02_pass29.root";
    //if(isSimulation)     EventCharaParFile = "/lustre/nyx/hades/user/bkardan/param/centrality_sim_au1230au_gen8a_UrQMD_minbias_2018_06.root";
#if isRealData || isEmbedding
    EventCharaParFile = "/lustre/nyx/hades/user/jorlinsk/MSc/centrality_epcorr_mar19_ag158ag_3200A_glauber_gen5_pass3_2021_08.root";  //
#else
    EventCharaParFile = "/lustre/nyx/hades/user/jorlinsk/MSc/centrality_epcorr_mar19_ag158ag_3200A_glauber_gen5_pass3_2021_08.root";
#endif

    if(!evtChara.setParameterFile(EventCharaParFile)){
	cout << "EventCharaParFile not found !!! " << endl;
	return kFALSE;
    }
    if(!evtChara.init()) {
	cout << "HParticleEvtChara not init!!! " << endl;
	return kFALSE;
    }
#if isEmbedding
    evtChara.printCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k10);
#else
    evtChara.printCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10);
#endif
    //-------------------------------------------------

    out = TFile::Open(outfile,"RECREATE");
//    out = TFile::Open(TFile::AsyncOpen(outfile,"UPDATE"));  //TODO: to be used if more than 1 task set is used!!!

    //-------------------------------------------------
    // added for the task
    HTaskSet *masterTaskSet = gHades->getTaskSet("all");

      masterTaskSet->add(new KZeroTrees("KZeroTrees","KZeroTrees",myInt1,myInt2,myInt3));
//    masterTaskSet->add(new MVATrees("MVATrees","MVATrees",outfile));

//    masterTaskSet->add(new KChargedTrees("KChargedTrees","KChargedTrees",myInt1,myInt2,myInt3));
//    masterTaskSet->add(new KPlusTree("KPlusTree","KPlusTree",outfile));
//    masterTaskSet->add(new PionTrees("PionTrees","PionTrees",myInt1,myInt2,myInt3));
//    masterTaskSet->add(new ProtonTree("ProtonTree","ProtonTree",myInt1,myInt2,myInt3));
//    masterTaskSet->add(new FillNtuple("FillNtuple","FillNtuple",outfile,myInt1,myInt2,myInt3));

//    masterTaskSet->add(new CUTSEARCHReco("CUTSEARCHReco","CUTSEARCHReco",outfile));

    //-------------------------------------------------

    //-------------------------------------------------
    // event loop starts here
    Int_t entries = loop.getEntries();
//    if(nEvents < entries && nEvents >= 0 ) entries = nEvents;

    for (Int_t i = 0; i < entries; i++) {
	Int_t nbytes =  loop.nextEvent(i);             // get next event. categories will be cleared before

//	if(nbytes == 0) { cout<<nbytes<<endl; break; } // last event reached -> end
//	if(nbytes < 0) { cout<<nbytes<<endl; continue; } // IoError -> skip (causing crash anyway)
	if(nbytes <= 0) { cout<<nbytes<<endl; break; } // last or broken event reached -> end


	if(i%1000 == 0) cout<<"event "<<i<<endl;

//	cout << "Event " << i << ":  ";

    } // end eventloop


    //-------------------------------------------------
    // added for the task
    if(gHades) gHades->finalizeTasks();
    //-------------------------------------------------

    out->Save();
    out->Close();

    delete gHades;
    return 0;
}
