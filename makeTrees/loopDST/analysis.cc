#include "mylibs.h"
#include "mystruc.h"
#include "loopDST_task.C"

using namespace std;

int main(int argc, char **argv)
{
    TROOT Analysis("Analysis","compiled analysis macro");

    // argc is the number of arguments in char* array argv
    // CAUTION: argv[0] contains the progname
    // argc has to be nargs+1
    cout<<argc<<" arguments "<<endl;
    if(argc>1) cout<<"arg1 ="<<argv[1]<<endl;
    if(argc>2) cout<<"arg2 ="<<argv[2]<<endl;
    if(argc>3) cout<<"arg3 ="<<argv[3]<<endl;
    if(argc>4) cout<<"arg4 ="<<argv[4]<<endl;
    if(argc>5) cout<<"arg5 ="<<argv[5]<<endl;

    TString nevts, myIntString1, myIntString2, myIntString3;
    switch (argc)
    {
    case 6:       // just inputfile name + nEvents
//	nevts  = argv[3];
        myIntString1 = argv[3];
        myIntString2 = argv[4];
        myIntString3 = argv[5];
	return loopDST_task(TString(argv[1]),TString(argv[2]),myIntString1.Atoi(),myIntString2.Atoi(),myIntString3.Atoi());
	break;
    default:
	cerr<<"ERROR: analysis() : WRONG NUMBER OF ARGUMENTS! TString infile="",TString outfile="", nevents=1000, myIntString1=-1, myIntString2=-1, myIntString3=-1"<<endl;

	return 1; // fail
    }
}
