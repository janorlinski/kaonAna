#ifndef __USEFULFUNCTIONS___     // good style to avaoid multiple includes
#define __USEFULFUNCTIONS___

#include "mylibs.h"
#include "mystruc.h"

#include "GlobVars.h"

using namespace std;

static Bool_t selectNegativeHadrons(HParticleCand* pcand )
{
    // build in selection function for hadron candidates.
    // Requires besides RK + META and fitted inner+outer MDC.

    if(//HParticleTrackSorter::kUseFakeRejection&&
       pcand->isFakeRejected()) return kFALSE;
    Bool_t test=kFALSE;
    if( pcand->isFlagAND(4,
                         Particle::kIsAcceptedHitInnerMDC,
                         Particle::kIsAcceptedHitOuterMDC,
                         Particle::kIsAcceptedHitMETA,
                         Particle::kIsAcceptedRK
                        )
       &&
       pcand->getInnerSegmentChi2() > 0
       &&
       pcand->getChi2()             < 1000
       &&
       pcand->getMetaMatchQuality() < 6
       &&
	//       HParticleTool::isGoodMetaCell(pcand,4,kTRUE)
       HParticleTool::isGoodMetaCell(pcand,10,kTRUE)
       &&
       pcand->getTof() < 60
      &&
       pcand->getBeta() > 0
	) test=kTRUE;
    if(!test) return kFALSE;

    return test;
}


vector<int> getMyCentralityClass(int originalCC){
    vector<int> myCC;
    myCC.clear();

    switch(originalCC){

    case 1:  // 0 - 5 %
	{
	    myCC.push_back( 0 ); // 0 - 5 %
	    myCC.push_back( 8 ); // 0 - 10 %
	    myCC.push_back( 12 ); // 0 - 20 %
            myCC.push_back( 13 ); // 0 - 40 %
	}
	break;

    case 2:  // 5 - 10 %
	{
	    myCC.push_back( 1 ); // 5 - 10 %
	    myCC.push_back( 8 ); // 0 - 10 %
	    myCC.push_back( 12 ); // 0 - 20 %
	    myCC.push_back( 13 ); // 0 - 40 %
	    myCC.push_back( 14 ); // 5 - 15 %
	    myCC.push_back( 15 ); // 5 - 20 %
	    myCC.push_back( 16 ); // 5 - 30 %
	    myCC.push_back( 17 ); // 5 - 40 %
	}
	break;

    case 3:  // 10 - 15 %
	{
	    myCC.push_back( 2 ); // 10 - 15 %
	    myCC.push_back( 9 ); // 10 - 20 %
	    myCC.push_back( 12 ); // 0 - 20 %
	    myCC.push_back( 13 ); // 0 - 40 %
	    myCC.push_back( 14 ); // 5 - 15 %
	    myCC.push_back( 15 ); // 5 - 20 %
	    myCC.push_back( 16 ); // 5 - 30 %
	    myCC.push_back( 17 ); // 5 - 40 %
	    myCC.push_back( 18 ); // 10 - 30 %
	    myCC.push_back( 19 ); // 10 - 40 %
	}
	break;

    case 4:  // 15 - 20 %
	{
	    myCC.push_back( 3 ); // 15 - 20 %
	    myCC.push_back( 9 ); // 10 - 20 %
	    myCC.push_back( 12 ); // 0 - 20 %
	    myCC.push_back( 13 ); // 0 - 40 %
	    myCC.push_back( 15 ); // 5 - 20 %
	    myCC.push_back( 16 ); // 5 - 30 %
	    myCC.push_back( 17 ); // 5 - 40 %
	    myCC.push_back( 18 ); // 10 - 30 %
	    myCC.push_back( 19 ); // 10 - 40 %
	}
	break;

    case 5:  // 20 - 25 %
	{
	    myCC.push_back( 4 ); // 20 - 25 %
	    myCC.push_back( 10 ); // 20 - 30 %
	    myCC.push_back( 13 ); // 0 - 40 %
	    myCC.push_back( 16 ); // 5 - 30 %
	    myCC.push_back( 17 ); // 5 - 40 %
	    myCC.push_back( 18 ); // 10 - 30 %
	    myCC.push_back( 19 ); // 10 - 40 %
	    myCC.push_back( 20 ); // 20 - 40 %
	}
	break;

    case 6:  // 25 - 30 %
	{
	    myCC.push_back( 5 ); // 25 - 30 %
	    myCC.push_back( 10 ); // 20 - 30 %
	    myCC.push_back( 13 ); // 0 - 40 %
	    myCC.push_back( 16 ); // 5 - 30 %
	    myCC.push_back( 17 ); // 5 - 40 %
	    myCC.push_back( 18 ); // 10 - 30 %
	    myCC.push_back( 19 ); // 10 - 40 %
	    myCC.push_back( 20 ); // 20 - 40 %
	}
	break;

    case 7:  // 30 - 35 %
	{
	    myCC.push_back( 6 ); // 30 - 35 %
	    myCC.push_back( 11 ); // 30 - 40 %
	    myCC.push_back( 13 ); // 0 - 40 %
	    myCC.push_back( 17 ); // 5 - 40 %
	    myCC.push_back( 19 ); // 10 - 40 %
	    myCC.push_back( 20 ); // 20 - 40 %
	}
	break;

    case 8:  // 35 - 40 %
	{
	    myCC.push_back( 7 ); // 35 - 40 %
	    myCC.push_back( 11 ); // 30 - 40 %
	    myCC.push_back( 13 ); // 0 - 40 %
	    myCC.push_back( 17 ); // 5 - 40 %
	    myCC.push_back( 19 ); // 10 - 40 %
	    myCC.push_back( 20 ); // 20 - 40 %
	}
	break;

    default:
	{
            //cout << "centrality is out of my region of interest" << endl;
	}
        break;
    }

    return myCC;
}

Double_t EnergyLossMdc(Double_t* x, Double_t* par)
{
    // Calculates the dEdx (MeV cm2/g ) of an particle with GEANT ID id
    // and momentum p (MeV) for He/i-butan gas mixture with He fraction hefr
    // (He (hefr) / i-butan (1-hefr)). The fomular is taken from PDG and doesn't
    // include the density correction term. The values for the mean excitation
    // energy are taken from Sauli.
    Double_t p    = x[0];
    Int_t id      = (Int_t)par[0];
    Double_t hefr = par[1];
    if(p==0)             return -1;
    if(hefr<0.||hefr>1.) return -1;
    Double_t mass    = HPhysicsConstants::mass(id);
    if(mass==0) return -1;
    Double_t Z_gas   = 2.*hefr+(1.-hefr)*34.;
    Double_t A_gas   = 4.*hefr+(1.-hefr)*58.;
    Double_t I_0_gas = 24.6*hefr+(1.-hefr)*10.8;
    Double_t I2      = pow(I_0_gas*Z_gas*(1.e-6),2); // sauli
    //Double_t I2     =pow(16.*pow(Z_gas,0.9),2); //C.Lippmann thesis

    Double_t K       = 0.307075; // MeV cm2 PDG, 4*pi*N_{A}*r_{e}^2*m_{e}^2*c^2
    Double_t mass2   = pow(mass,2);
    Double_t m_ec2   = HPhysicsConstants::mass(3);
    Double_t z2      = pow((Float_t)HPhysicsConstants::charge(id),2);
    Double_t p2      = pow(p,2);
    Double_t beta2   = 1./((mass2/p2)+1.);
    Double_t gamma2  = 1./(1-beta2);
    Double_t gamma   = sqrt(gamma2);
    Double_t Tmax    = (2.*m_ec2*beta2*gamma2)/(1.+ 2.*gamma*(m_ec2/mass)+pow((m_ec2/mass),2));
    Double_t term1   = K*z2*(Z_gas/A_gas)*(1./beta2);
    Double_t term2   = ((2.*m_ec2*beta2*gamma2*Tmax)/I2);
    Double_t dedx    = term1 * (0.5*log(term2)-beta2);

    return dedx;

}



/*
void clearPID(BRANCH &pid)
{
    pid.eventnumber = -1.;
    pid.downflag = -1.;
    pid.trigbit = -1.;
    pid.mult_tofrpc = -1.;
    pid.mult_class = -1.;
    pid.start_strip = -1.;
    pid.system = -1.;
    pid.sector = -1.;
    pid.tof = -1.;
    pid.theta = -1.;
    pid.phi = -1000.;
    pid.mom = -1.;
    pid.beta = -1.;
    pid.chi2 = 1000000000.;
    pid.charge = -1000.;
    pid.mass = -1.;
    pid.isused = -1.;
    pid.metamatch_quality = 1000.;
    pid.ishadron = -1.;
    pid.zprime = -1000.;
    pid.id = -1.;
    pid.mult_wall = -1.;
    pid.RP = -1000.;
    pid.RPAB = -1000.;
}
*/
#endif
