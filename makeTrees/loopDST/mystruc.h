#ifndef MYSTRUCT_H
#define MYSTRUCT_H

#include "mylibs.h"

typedef struct {
    Float_t p1_id, p1_geantId,
    p1_z, p1_r, p1_baseX, p1_baseY, p1_baseZ, p1_dirX, p1_dirY, p1_dirZ,
    p1_innerChi2, p1_outerChi2, p1_RKChi2, p1_metaMatchQa, p1_distMetaHit,
    p1_system,p1_sector,p1_theta,p1_phi,
    p1_mom,p1_mom_corr,p1_mass,p1_beta,p1_charge,p1_tof,p1_metahit,
    p1_dEdxTOF,p1_dEdxMDC,p1_isAtEdge,
    p2_id, p2_geantId,
    p2_z, p2_r, p2_baseX, p2_baseY, p2_baseZ, p2_dirX, p2_dirY, p2_dirZ,
    p2_innerChi2, p2_outerChi2, p2_RKChi2, p2_metaMatchQa, p2_distMetaHit,
    p2_system,p2_sector,p2_theta,p2_phi,
    p2_mom,p2_mom_corr,p2_mass,p2_beta,p2_charge,p2_tof,p2_metahit,
    p2_dEdxTOF,p2_dEdxMDC,p2_isAtEdge,
    eventVerX,eventVerY,eventVerZ,evntVerChi2,
    decayVerX,decayVerY,decayVerZ,
    dist1,dist2,dist3,distVer,distMin,
    multiplicity,invMass,oAngle;
} CAND_DOUBLET;

typedef struct {
    Float_t eventPlane,phiAB,
    eventVerX,eventVerY,eventVerZ,evntVerChi2,
    decayVerX,decayVerY,decayVerZ,
    dist1,dist2,dist3,distVer,distMin,oAngle,
    theta,phi,momX,momY,momZ,mom,
    multiplicity,centralityclass,
    invMass,mt,mtm0,pt,y,ycm,y0,
    flow,
    occupancyCorr,accCorr,recoEff,eventWeight,
    mvaResult,
    downflag;
} CAND_K0S;

typedef struct {
    Float_t eventPlane,phiAB,
    eventVerX,eventVerY,eventVerZ,evntVerChi2,
    system,theta,phi,momX,momY,momZ,mom,
    multiplicity,centralityclass,MDCtracks,METAclst,
    invMass,mass,mt,mtm0,pt,y,ycm,y0,
    flow,
    occupancyCorr,accCorr,recoEff,eventWeight,
    phiMass,downflag;
} CAND_KCHARGED;

typedef struct {
    Float_t eventPlane,phiAB,
    eventVerX,eventVerY,eventVerZ,evntVerChi2,
    system,theta,phi,momX,momY,momZ,mom,
    multiplicity,centralityclass,MDCtracks,METAclst,
    invMass,mt,mtm0,pt,y,ycm,y0,
    flow,
    occupancyCorr,accreceffCorr,eventWeight,
    downflag;
} CAND_KPLUS;

typedef struct {
    Float_t eventPlane,eventPlane_gkine,phiAB,
    //eventVerX,eventVerY,eventVerZ,evntVerChi2,
    //decayVerX,decayVerY,decayVerZ,
    //dist1,dist2,dist3,distVer,distMin,oAngle,
    theta,phi,momX,momY,momZ,mom,
    multiplicity,centralityclass,impactparam,
    invMass,mt,pt,y,ycm,y0,
    flow,flow_gkine,weight,
    dist1,dist2,dist3,distVer,distMin,dist1Org,dist2Org,dist3Org,distVerOrg,distMinOrg,mvaResult;
    Bool_t pipAcc,pimAcc,pipDecay,pimDecay,pipRec,pimRec,pichargedDecay,pizeroDecay;
} GKINE_K0S;

typedef struct {
    Float_t eventPlane,eventPlane_gkine,phiAB,
    //eventVerX,eventVerY,eventVerZ,evntVerChi2,
    system,theta,phi,momX,momY,momZ,mom,
    multiplicity,centralityclass,impactparam,MDCtracks,METAclst,
    invMass,mt,mtm0,pt,y,ycm,y0,
    flow,flow_gkine,weight;
    Bool_t kpAcc,kpRec,hasDecayed;
} GKINE_KPLUS;

typedef struct {
    Float_t eventPlane,eventPlane_gkine,phiAB,
    //eventVerX,eventVerY,eventVerZ,evntVerChi2,
    system,theta,phi,momX,momY,momZ,mom,
    multiplicity,centralityclass,impactparam,MDCtracks,METAclst,
    invMass,mt,mtm0,pt,y,ycm,y0,
#if isEmbedding
    pcand_chi2, pcand_innerchi2, pcand_mm, pcand_mass, pcand_system, pcand_tofeloss, pcand_mdceloss, pcand_beta, pcand_mom, pcand_charge,
#endif
    flow,flow_gkine,weight;
    Bool_t kAcc,kRec,
#if isEmbedding
    pcand_mdcedge, pcand_goodflags,
#endif
    hasDecayed;
} GKINE_KCHARGED;

//typedef struct {
//    Float_t eventnumber,downflag,trigbit,
//    mult_tofrpc,mult_class,start_strip,system,sector,
//    tof,theta,phi,mom,beta,chi2,charge,mass,isused,
//    metamatch_quality,ishadron,zprime,id,
//    mult_wall,RP,RPAB;
//} BRANCH;

typedef struct {
    Float_t chi2, metamatch, beta, betaOrg, tof, mom, momOrg, momcorr,
    charge, mass, system, theta, phi, tofdedx, mdcdedx,
    my_kaonmass, kaonmass,
    goodTofEloss, goodMdcEloss,
    mtm0,pt,y;
} BRANCH;

#endif
