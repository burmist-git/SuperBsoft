//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description:
//	Class BtaPidDstarDFinder: D(+-)* -> D0 pi(+-)_s, D0 ->  k pi finder.
//     Adapted from BetaPidCalib/BtaPidDstarSample.cc
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Nicolas ARNAUD (SuperB)
//      Giampiero Mancinelli                    Original author
//
// Copyright Information:
//	 (C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
//-----------------------
// This Class's Header --
//-----------------------
#include "PacPidCalib/PacPidDstarSample.hh"

//-------------
// C Headers --
//-------------
#include <assert.h>

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include <string>
using std::string;
#include "AbsEventTag/AbsEventTag.hh"
#include "AssocTools/AstNamedMapVector.hh"

#include "CLHEP/Alist/ConstAList.h"
#include "CLHEP/Alist/AList.h"
#include "CLHEP/Alist/ConstAIterator.h"
#include "CLHEP/HepPoint.h"

#include "AbsEventTag/AbsEventTag.hh"

#include "AbsEnv/AbsEnv.hh"
#include "AbsEvent/AbsEvent.hh"
#include "GenEnv/GenEnv.hh"

#include "Beta/EventInfo.hh"
#include "Beta/BtaCandidate.hh"
#include "Beta/BtaAbsVertex.hh"
#include "BetaCoreTools/BtaMcAssoc.hh"
#include "BetaEvent/BtaParam.hh"

#include "BetaMicroAdapter/BtaMicroAdapter.hh"
#include "BetaMicroAdapter/BtaPidInfo.hh"
#include "BetaMicroAdapter/BtaTrkQual.hh"

#include "RecVtx/BtaOpVtxKs.hh"
#include "FastVtx/BtaOpFastVtx.hh"
#include "FastVtx/BtaOpFastVtxV0.hh"
#include "BetaCoreTools/BtaOpAdd4.hh"

#include "PDT/PdtEntry.hh"
#include "PDT/Pdt.hh"

#include "HepTuple/TupleManager.h"
#include "HepTuple/Tuple.h"

#include "ProbTools/probab.hh"
#include "ProbTools/Consistency.hh"

#include "BbrGeom/BbrPointErr.hh"  


#include "TrkBase/TrkExchangePar.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkFit.hh"
#include "TrajGeom/TrkLineTraj.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "AbsPid/PidInfoSummary.hh"

#include "Framework/AppModule.hh"

#include "ProxyDict/Ifd.hh"
#include "ProxyDict/IfdStrKey.hh"
using std::cerr;
using std::endl;
using std::fstream;

using namespace std;

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//----------------
// Constructors --
//----------------


PacPidDstarSample::PacPidDstarSample( AppModule *callingModule ) :
  PacPidCalibSample(callingModule, "PacPidDstarSample"),
  _outlisttmp(new HepAList<BtaCandidate>),
  _outlistdstartmp(new HepAList<BtaCandidate>)
{
  HepString listName("ChargedTracks");
  _pidDstarInputList1 = addNewParmString("pidDstarInputList1", listName, callingModule);

  listName = "Dstar";
  _PidDstar = addNewParmString("PidDstar", listName, callingModule);

  listName = "D0";
  _PidDstarD0 = addNewParmString("PidDstarD0", listName, callingModule);

  listName = "SoftPion";
  _mapPidDstarSoftPion = addNewParmString("mapPidDstarSoftPion", listName, callingModule);

  listName = "Pion";
  _mapPidDstarPion = addNewParmString("mapPidDstarPion", listName, callingModule);

  listName = "Kaon";
  _mapPidDstarKaon = addNewParmString("mapPidDstarKaon", listName, callingModule);

  listName = "PidD0";
  _PidD0 = addNewParmString("PidD0", listName, callingModule);

  _pidDstarMassCutD0 = addNewParmDouble("pidDstarMassCutD0",1.845,callingModule);
  _pidDstarMassCutD0Upper = addNewParmDouble("pidDstarMassCutD0Upper",2.12,callingModule);
  _pidDstarLooseMassCutD0 = addNewParmDouble("pidDstarLooseMassCutD0",1.4,callingModule);
  _pidDstarLooseMassCutD0Upper = addNewParmDouble("pidDstarLooseMassCutD0Upper",2.32,callingModule);
  _pidDstarCosThetaD0Cut = addNewParmDouble("pidDstarCosThetaD0Cut",0.90,callingModule);
  _pidDstarMomD0Cut = addNewParmDouble("pidDstarMomD0Cut",1.4,callingModule);
  _pidDstarMassCutDD = addNewParmDouble("pidDstarMassCutDD",.14375,callingModule);
  _pidDstarMassCutDDUpper = addNewParmDouble("pidDstarMassCutDDUpper",.14715,callingModule);
  _pidDstarMassCutDstar = addNewParmDouble("pidDstarMassCutDstar",1.99,callingModule);
  _pidDstarMassCutDstarUpper = addNewParmDouble("pidDstarMassCutDstarUpper",2.03,callingModule);
  _pidDstarD0PionCut = addNewParmDouble("pidDstarD0PionCut",1000.,callingModule);
  _pidDstarD0KaonCut = addNewParmDouble("pidDstarD0KaonCut",1000.,callingModule);
  _pidDstarDocaPiKCut = addNewParmDouble("pidDstarDocaPiKCut",0.05,callingModule);
  _pidDstarDocaPisCut = addNewParmDouble("pidDstarDocaPisCut",0.3,callingModule);
  _pidDstarAper_dpsCut = addNewParmDouble("pidDstarAper_dpsCut",0.40,callingModule);
  _pidDstarAper_kpiCut = addNewParmDouble("pidDstarAper_kpiCut",2.00,callingModule);
  _pidDstarHelicCut = addNewParmDouble("pidDstarHelicCut",-0.90,callingModule);
  _pidDstarMaxPiSoftMom = addNewParmDouble("pidDstarMaxPiSoftMom",.5,callingModule);
  _pidDstarconstrDauMassKaon = addNewParmBool("pidDstarconstrDauMassKaon",true,callingModule);
  _pidDstarconstrDauMassPion = addNewParmBool("pidDstarconstrDauMassPion",true,callingModule);
  _pidDstarconstrDauMassD0 = addNewParmBool("pidDstarconstrDauMassD0",false,callingModule);
  _pidDstarconstrDauMassSoftPi = addNewParmBool("pidDstarconstrDauMassSoftPi",true,callingModule);
  _pidDstarFillNtp = addNewParmBool("pidDstarFillNtp",false,callingModule);
  _pidDstarCheckTruth = addNewParmBool("pidDstarCheckTruth",false,callingModule);
  _pidDstarUseCompOutput = addNewParmBool("pidDstarUseCompOutput",false,callingModule);
}

PacPidDstarSample::~PacPidDstarSample() {
  delete _outlisttmp;
  delete _outlistdstartmp;
}
//-----------------------------------------------------------------------
// The following function contains the way the pi and K are combined to
// get a D0 or an anti-D0
//-----------------------------------------------------------------------

void PacPidDstarSample::combine(BtaCandidate* cndT1, BtaCandidate* cndT2) 
{
  assert(cndT1 != 0);
  assert(cndT2 != 0);
  BtaCandidate* cnd1 = new BtaCandidate(*cndT1);
  BtaCandidate* cnd2 = new BtaCandidate(*cndT2);
      const PdtEntry* kaon;
      const PdtEntry* pion;
      
      kaon = (cnd1->charge()) < 0 ? Pdt::lookup(PdtLund::K_minus) : 
	Pdt::lookup(PdtLund::K_plus) ;
      cnd1->setType(kaon);
  
      pion = (cnd2->charge()) < 0 ? Pdt::lookup(PdtLund::pi_minus) : 
	Pdt::lookup(PdtLund::pi_plus) ;
      cnd2->setType(pion);

  const string d0Name( _PidD0->value() );

  //  BtaOpVtx makeVtxD;             // vertex operator for doing combinations
  BtaOpFastVtx makeVtxD; 
  const PdtEntry* dZero;
  
  BtaCandidate* newTrkptr = 0;
  
  newTrkptr = makeVtxD.create(*cnd1, *cnd2);
  // use the result only if a valid vertex was found
  if ( newTrkptr != 0  ) {
    if( newTrkptr->decayVtx() == 0){
      delete newTrkptr;  // not keeping, so must delete 
      delete cnd1;
      delete cnd2;
      newTrkptr = 0;
    }
  }
  else {
      delete cnd1;
      delete cnd2;
  }    
  if( 0 != newTrkptr )
    {
      double motherInvMass = constrainMass(cnd1, cnd2, newTrkptr, 1);
      newTrkptr->setMass(motherInvMass);
     
      if (cnd1->charge() < 0)  dZero=Pdt::lookup(PdtLund::D0);
      else dZero=Pdt::lookup(PdtLund::anti_D0) ; 
      newTrkptr->setType(dZero);
      if ( acceptD0(newTrkptr) ){
      //string d0Name( _PidD0->value() );
        bool useForFilter(false);
        selectThis( newTrkptr , d0Name, useForFilter);  
        _outlisttmp->append(newTrkptr);	 
	delete cnd1;
	delete cnd2;
      } else {
        delete newTrkptr; // not kept, delete it
	delete cnd1;
	delete cnd2;
      }
    }
}

bool PacPidDstarSample::acceptDstar(BtaCandidate* newTrk, const
				    BtaCandidate* cnd1, const
				    BtaCandidate* cnd2, const  
				    BtaCandidate* dauKaon, const
				    BtaCandidate* dauPion) const  
{
  
  bool isAccepted(true);
  double motherInvMass = newTrk->mass();
  motherInvMass = constrainMass(cnd1, cnd2, newTrk, 3);
  const PdtEntry* dstar;
  dstar = (newTrk->charge()) < 0 ? Pdt::lookup(PdtLund::D_star_minus) : Pdt::lookup(PdtLund::D_star_plus) ;
  newTrk->setType(dstar);
  newTrk->setMass(motherInvMass);
  double deltaM=motherInvMass - cnd1->mass();
  if (( deltaM < _cut1DD || 
	deltaM > _cut2DD ) || (motherInvMass < _cut1Dstar || 
			       motherInvMass > _cut2Dstar )){ 
    return false;
  }  
  // define extra cuts
  double d0pion(10000); 
  double d0kaon(10000); 
  // double d0pis(10000); 
  double docaPion(10000); 
  double docaKaon(10000); 
  double docaPis(10000); 
  double s=0, t=0, u=0;
  HepPoint vtxp(0.,0.,0.);
  vtxp = cnd1->decayVtx()->point();
  
  if (dauPion->recoTrk()) {
    if (dauPion->recoTrk()->fitResult()) {		  
      d0pion=fabs(dauPion->recoTrk()->fitResult()->helix(0.0).d0());
      if ( d0pion > _cutD0Pion ) return false;
      TrkPoca poca1(dauPion->recoTrk()->fitResult()->traj(), s, vtxp);  
      if (poca1.status().success()) {
	docaPion =poca1.doca();
      } else {
	return false;
      }
    }
  }

  if (d0pion > _cutD0Pion) return false;

  if (dauKaon->recoTrk()) {
    if (dauKaon->recoTrk()->fitResult()) {
      d0kaon=fabs(dauKaon->recoTrk()->fitResult()->helix(0.0).d0());
      if (d0kaon > _cutD0Kaon) return false;
      TrkPoca poca2(dauKaon->recoTrk()->fitResult()->traj(), t, vtxp);  
      if (poca2.status().success()) {
	docaKaon = poca2.doca();
      } else {
	return false;
      }
    }
  }
  if (d0kaon > _cutD0Kaon) return false;
  if ((docaPion+docaKaon) > _cutDocaPiK ) return false;

  if (cnd2->recoTrk()) {
    if (cnd2->recoTrk()->fitResult()) {
      // d0pis=abs(cnd2->recoTrk()->fitResult()->helix(0.0).d0());
      TrkPoca poca3(cnd2->recoTrk()->fitResult()->traj(), u, vtxp);  
      if (poca3.status().success()) {
	docaPis = poca3.doca();
      } else {
	return false;
      }
    }
  }

  
  Hep3Vector unit3(newTrk->p3().unit());
  
  HepLorentzVector pPip4(0);
  if (dauKaon->recoTrk()) {
    if (dauKaon->recoTrk()->fitResult()) {		  
      pPip4 = dauKaon->p4(vtxp);
    }
    else {
      pPip4 = dauKaon->p4();
    }
  }
  else {
    pPip4 = dauKaon->p4();
  }
  
  pPip4.boost(-(newTrk->p4().boostVector()));     
  
  
  double cosHelic;
  Hep3Vector pPip3(pPip4.x(),pPip4.y(),pPip4.z());
  pPip3 = pPip3.unit();
  cosHelic = pPip3.dot(unit3); 
  if  ( d0pion < _cutD0Pion &&  d0kaon < _cutD0Kaon  &&
	(docaPion+docaKaon) < _cutDocaPiK &&
	cosHelic > _cutHelic ) {
    if (docaPis > _cutDocaPis) {
      return false;
    } 
  } else {
    return false;
  }
  
  Hep3Vector unit1(dauPion->p3(vtxp).unit());
  Hep3Vector unit2(dauKaon->p3(vtxp).unit()); 
  Hep3Vector unit6(cnd2->p3(vtxp).unit()); 
  Hep3Vector unit7(cnd1->p3().unit()); 
  double apertureKPi(acos(unit1.dot(unit2))); 
  double apertureD0PiSoft(acos(unit6.dot(unit7))); 
  
  if  ( apertureKPi > _cutAper_kpi ||  apertureD0PiSoft > _cutAper_dps ) {
    isAccepted = false;
  }
  
  return isAccepted;
}

bool PacPidDstarSample::acceptD0(const BtaCandidate* newCand) const 
{
  bool isAccepted(true);
  
  // apply mass cut (tight) on D0
  if( newCand->mass() < _cut1D0 || newCand->mass() > _cut2D0 ){
    return false;
  }
  if( newCand->p3().cosTheta() > _cutCosThetaD0 || newCand->p() < _cutMomD0 ){
    return false;
  }
  return isAccepted;
}

void PacPidDstarSample::fillNtuple(const BtaCandidate* Dstar, const
				   BtaCandidate* D0,  const
				   BtaCandidate* piSoft,   const
				   BtaCandidate* pion1,   const
				   BtaCandidate* kaon1, const HepPoint
				   primVtx) const  
{
  
  _ntuple->column("mass_Dst", Dstar->mass());
  _ntuple->column("ener_Dst", Dstar->energy());
  _ntuple->column("q_Dst", Dstar->charge());
  _ntuple->column("p_Dstar",Dstar->p(), -10.);
  _ntuple->column("px_Dstar",Dstar->p3().x(), -10.); 
  _ntuple->column("py_Dstar",Dstar->p3().y(), -10.); 
  _ntuple->column("pz_Dstar",Dstar->p3().z(), -10.); 
  _ntuple->column("cost_Dst",Dstar->p3().cosTheta(), -10.); 
  _ntuple->column("phi_Dst",Dstar->p3().phi(), -10.); 
  if ( _truthMap ) { 
    BtaCandidate* candMC(0);
    candMC = _truthMap->mcFromReco(Dstar);
    if ( candMC ){
      _ntuple->column("DstMCmas",candMC->mass());
      _ntuple->column("DstMCtyp",candMC->pdtEntry()->lundId());
      _ntuple->column("DstMCpx",candMC->p3().x());
      _ntuple->column("DstMCpy",candMC->p3().y());
      _ntuple->column("DstMCpz",candMC->p3().z());
      
      BtaCandidate* Mother = candMC->theMother();
      
      _ntuple->column("GmDstMCm",Mother->mass());
      _ntuple->column("GmDstMCt",Mother->pdtEntry()->lundId());
    }
    else{
      _ntuple->column("DstMCmas",-10.);
      _ntuple->column("DstMCtyp",-15000000);
      _ntuple->column("DstMCpx",-100.);
      _ntuple->column("DstMCpy",-100.);
      _ntuple->column("DstMCpz",-100.);
      _ntuple->column("GmDstMCm",-10.);
      _ntuple->column("GmDstMCt",-15000000);    }
  }
  _ntuple->column("mass_D0", D0->mass());
  _ntuple->column("ener_D0", D0->energy());
  _ntuple->column("q_D0", D0->charge());
  _ntuple->column("p_D0",D0->p(), -10.);
  _ntuple->column("px_D0",D0->p3().x(), -10.); 
  _ntuple->column("py_D0",D0->p3().y(), -10.); 
  _ntuple->column("pz_D0",D0->p3().z(), -10.); 
  _ntuple->column("cost_D0",D0->p3().cosTheta(), -10.); 
  _ntuple->column("phi_D0",D0->p3().phi(), -10.); 
  if ( _truthMap ) { 
    BtaCandidate* candMC(0);
    candMC = _truthMap->mcFromReco(D0);
    if ( candMC ){
      _ntuple->column("D0MCmas",candMC->mass());
      _ntuple->column("D0MCtyp",candMC->pdtEntry()->lundId());
      _ntuple->column("D0MCpx",candMC->p3().x());
      _ntuple->column("D0MCpy",candMC->p3().y());
      _ntuple->column("D0MCpz",candMC->p3().z());
      
      BtaCandidate* Mother = candMC->theMother();
      
      _ntuple->column("GmD0MCm",Mother->mass());
      _ntuple->column("GmD0MCt",Mother->pdtEntry()->lundId());
    }
    else{
      _ntuple->column("D0MCmas",-10.);
      _ntuple->column("D0MCtyp",-15000000);
      _ntuple->column("D0MCpx",-100.);
      _ntuple->column("D0MCpy",-100.);
      _ntuple->column("D0MCpz",-100.);
      _ntuple->column("GmD0MCm",-10.);
      _ntuple->column("GmD0MCt",-15000000);
    }
  }
  _ntuple->column("mass_PiS", piSoft->mass());
  _ntuple->column("q_PiS", piSoft->charge());
  _ntuple->column("p_PiS",piSoft->p(), -10.);
  _ntuple->column("px_PiS",piSoft->p3().x(), -10.); 
  _ntuple->column("py_PiS",piSoft->p3().y(), -10.); 
  _ntuple->column("pz_PiS",piSoft->p3().z(), -10.); 
  _ntuple->column("cost_PiS",piSoft->p3().cosTheta(), -10.); 
  _ntuple->column("phi_PiS",piSoft->p3().phi(), -10.); 
  if ( _truthMap ) { 
    BtaCandidate* candMC(0);
    candMC = _truthMap->mcFromReco(piSoft);
    if ( candMC ){
      _ntuple->column("PiSMCmas",candMC->mass());
      _ntuple->column("PiSMCtyp",candMC->pdtEntry()->lundId());
      _ntuple->column("PiSMCpx",candMC->p3().x());
      _ntuple->column("PiSMCpy",candMC->p3().y());
      _ntuple->column("PiSMCpz",candMC->p3().z());
      
      BtaCandidate* Mother = candMC->theMother();
      if (Mother) {      
	_ntuple->column("GmPiSMCm",Mother->mass());
	_ntuple->column("GmPiSMCt",Mother->pdtEntry()->lundId());
      }
      else {
	_ntuple->column("GmPiSMCm", -10.);
	_ntuple->column("GmPiSMCt", -15000000);
      }
    }
    else{
      _ntuple->column("PiSMCmas",-10.);
      _ntuple->column("PiSMCtyp",-15000000);
      _ntuple->column("PiSMCpx",-100.);
      _ntuple->column("PiSMCpy",-100.);
      _ntuple->column("PiSMCpz",-100.);
      _ntuple->column("GmPiSMCm",-10.);
      _ntuple->column("GmPiSMCt",-15000000);
    }
  }
  _ntuple->column("mass_pi", pion1->mass());
  _ntuple->column("q_pi", pion1->charge());
  _ntuple->column("p_pi",pion1->p(), -10.);
  _ntuple->column("px_pi",pion1->p3().x(), -10.); 
  _ntuple->column("py_pi",pion1->p3().y(), -10.); 
  _ntuple->column("pz_pi",pion1->p3().z(), -10.); 
  _ntuple->column("cost_pi",pion1->p3().cosTheta(), -10.); 
  _ntuple->column("phi_pi",pion1->p3().phi(), -10.); 
  if ( _truthMap ) { 
    BtaCandidate* candMC(0);
    candMC = _truthMap->mcFromReco(pion1);
    if ( candMC ){
      _ntuple->column("piMCmas",candMC->mass());
      _ntuple->column("piMCtyp",candMC->pdtEntry()->lundId());
      _ntuple->column("piMCpx",candMC->p3().x());
      _ntuple->column("piMCpy",candMC->p3().y());
      _ntuple->column("piMCpz",candMC->p3().z());
      
      BtaCandidate* Mother = candMC->theMother();
      if (Mother) {      
	_ntuple->column("GmPiMCm",Mother->mass());
	_ntuple->column("GmPiMCt",Mother->pdtEntry()->lundId());
      }
      else {
	_ntuple->column("GmPiMCm", -10.);
	_ntuple->column("GmPiMCt", -15000000);
      }
    }
    else{
      _ntuple->column("piMCmas",-10.);
      _ntuple->column("piMCtyp",-15000000);
      _ntuple->column("piMCpx",-100.);
      _ntuple->column("piMCpy",-100.);
      _ntuple->column("piMCpz",-100.);
      _ntuple->column("GmPiMCm",-10.);
      _ntuple->column("GmPiMCt",-15000000);
    }
  }
  _ntuple->column("mass_Kn", kaon1->mass());
  _ntuple->column("q_Kn", kaon1->charge());
  _ntuple->column("p_Kn",kaon1->p(), -10.);
  _ntuple->column("px_Kn",kaon1->p3().x(), -10.); 
  _ntuple->column("py_Kn",kaon1->p3().y(), -10.); 
  _ntuple->column("pz_Kn",kaon1->p3().z(), -10.); 
  _ntuple->column("cost_Kn",kaon1->p3().cosTheta(), -10.); 
  _ntuple->column("phi_Kn",kaon1->p3().phi(), -10.); 
  HepPoint vtxpK(0.,0.,0.);
  HepPoint primPointK(primVtx.x(),primVtx.y(),primVtx.z());
  Hep3Vector vtxpIPK(0.,0.,0.);
  if ( _truthMap ) { 
    BtaCandidate* candMC(0);
    candMC = _truthMap->mcFromReco(kaon1);
    if ( candMC ){
      _ntuple->column("KnMCmas",candMC->mass());
      _ntuple->column("KnMCtyp",candMC->pdtEntry()->lundId());
      _ntuple->column("KnMCpx",candMC->p3().x());
      _ntuple->column("KnMCpy",candMC->p3().y());
      _ntuple->column("KnMCpz",candMC->p3().z());
      if ( candMC->decayVtx() ) {
	//         const BtaAbsVertex* vtxK = candMC->decayVtx();
         vtxpK = candMC->decayVtx()->point();
         vtxpIPK = vtxpK - primPointK ;
      }
      _ntuple->column("KnMCdec",vtxpIPK.mag());
      
      BtaCandidate* Mother = candMC->theMother();
      if (Mother) { 
	_ntuple->column("GmKnMCm",Mother->mass());
	_ntuple->column("GmKnMCt",Mother->pdtEntry()->lundId());
	BtaCandidate* GMother = Mother->theMother();
	if (GMother) { 
	  _ntuple->column("GGKnMCt",GMother->pdtEntry()->lundId());
	}
	else {
	  _ntuple->column("GGKnMCt", -15000000);
	}
      }
      else {
	_ntuple->column("GmKnMCm", -10.);
	_ntuple->column("GmKnMCt", -15000000);
	_ntuple->column("GGKnMCt", -15000000);
      }
    }
    else{
      _ntuple->column("KnMCmas",-10.);
      _ntuple->column("KnMCtyp",-15000000);
      _ntuple->column("KnMCpx",-100.);
      _ntuple->column("KnMCpy",-100.);
      _ntuple->column("KnMCpz",-100.);
      _ntuple->column("GmKnMCm",-10.);
      _ntuple->column("GmKnMCt",-15000000);
      _ntuple->column("GGKnMCt",-15000000);
    }
  }
  
  // to be implemented
  bool missed(false);
  
  if(missed)
    _ntuple->column("missed",1);
  else 
    _ntuple->column("missed",0);
  
  //
  Hep3Vector unit3(D0->p3().unit());
  Hep3Vector unit5(Dstar->p3().unit());
  HepPoint vtxp(0.,0.,0.);
  HepPoint primPoint(primVtx.x(),primVtx.y(),primVtx.z());
  Hep3Vector vtxpIP(0.,0.,0.);
  
  
  
  
  if ( D0->decayVtx() ) {
    const BtaAbsVertex* vtx = D0->decayVtx();
    vtxp = D0->decayVtx()->point();
    vtxpIP = vtxp - primPoint ;
    //	cout <<primPoint.x()<<primPoint.y()<<primPoint.z()<<endl;  
    int ndof = vtx->nDof();
    double chi2= vtx->chiSquared();
    Hep3Vector unit4(vtxp.x(), vtxp.y(), vtxp.z()); 
    Hep3Vector unit4IP(vtxpIP.x(), vtxpIP.y(), vtxpIP.z()); 
    unit4 = unit4.unit();
    unit4IP = unit4IP.unit();
    //
    // We will have to sort this out. What a dumb CLHEP!
    //
//              BbrError covMatrix(primVtx.covMatrix());
//	      HepSymMatrixD covMat(covMatrix);
//	      HepVector vec(3);
//	      vec = unit4;        
//	      HepDouble vtxErr = vec.dot(vec, covMat*vec);
//	      double vtxErrIP = unit4IP.dot(covMatrix*unit4IP);
   
    _ntuple->column("vdis",vtx->point().mag());
    //	_ntuple->column("vtxErr",vtxErr);
    _ntuple->column("vx",vtx->point().x());
    _ntuple->column("vy",vtx->point().y());
    _ntuple->column("vz",vtx->point().z());        
    _ntuple->column("vdisIP",vtxpIP.mag());
    //	_ntuple->column("vtxErrIP",vtxErrIP);
    _ntuple->column("vxIP",vtxpIP.x());
    _ntuple->column("vyIP",vtxpIP.y());
    _ntuple->column("vzIP",vtxpIP.z());        
    _ntuple->column("chi2",chi2/ndof, -10.);
    _ntuple->column("pchi2",(double) probab(ndof,chi2), -10.);    
  }
  else{
    // can't use default args, if the columns aren't filled the first
    // time an error occurs
    _ntuple->column("vdis",-10.);
    //	_ntuple->column("vtxErr",-10.);
    _ntuple->column("vx",0.);
    _ntuple->column("vy",0.);
    _ntuple->column("vz",0.);        
    _ntuple->column("vdisIP",-10.);
    //	_ntuple->column("vtxErrIP",-10.);
    _ntuple->column("vxIP",0.);
    _ntuple->column("vyIP",0.);
    _ntuple->column("vzIP",0.);        
    _ntuple->column("chi2", -10.);
    _ntuple->column("pchi2", -10.);    
  }
  
  
  if( piSoft->recoTrk() ) {
    if( piSoft->recoTrk()->fitResult() ) {
      _ntuple->column("d0PiSoft", fabs(piSoft->recoTrk()->fitResult()->helix(0.0).d0()));
    }
    else _ntuple->column("d0PiSoft",10.);
  }
  else _ntuple->column("d0PiSoft",10.);
  if( pion1->recoTrk() ) {
    if( pion1->recoTrk()->fitResult() ) {
      _ntuple->column("d0Pion", fabs(pion1->recoTrk()->fitResult()->helix(0.0).d0()));
    }
    else _ntuple->column("d0Pion",10.);
  }
  else _ntuple->column("d0Pion",10.);
  if( kaon1->recoTrk() ) {
    if( kaon1->recoTrk()->fitResult() ) {
      _ntuple->column("d0Kaon", fabs(kaon1->recoTrk()->fitResult()->helix(0.0).d0()));
    }
    else _ntuple->column("d0Kaon",10.);
  }
  else _ntuple->column("d0Kaon",10.);
  
  double s=0, t=0, u=0;
  if( piSoft->recoTrk() ) {
    if( piSoft->recoTrk()->fitResult() ) {
      TrkPoca poca1(piSoft->recoTrk()->fitResult()->traj(), s, vtxp );  
      if (poca1.status().success()) {
	_ntuple->column("doca_PiS", poca1.doca(), -10.);
      }
      else _ntuple->column("doca_PiS",-10.);
    }
    else _ntuple->column("doca_Pis",-10.);
  }
  else _ntuple->column("doca_Pis",-10.);
  if( pion1->recoTrk() ) {
    if( pion1->recoTrk()->fitResult() ) {
      TrkPoca poca2(pion1->recoTrk()->fitResult()->traj(), t, vtxp);  
      if (poca2.status().success()) {
	_ntuple->column("doca_Pi", poca2.doca(), -10.);
      }
      else _ntuple->column("doca_Pi",-10.);
    }
    else _ntuple->column("doca_Pi",-10.);
  }
  else _ntuple->column("doca_Pi",-10.);
  if( kaon1->recoTrk()) {
    if( kaon1->recoTrk()->fitResult() ) {
      TrkPoca poca3(kaon1->recoTrk()->fitResult()->traj(), u, vtxp);  
      if (poca3.status().success()) {
	_ntuple->column("doca_Kn", poca3.doca(), -10.);
      }
      else _ntuple->column("doca_Kn",-10.);
    }
    else _ntuple->column("doca_Kn",-10.);
  }
  else _ntuple->column("doca_Kn",-10.);
  
  Hep3Vector unit1(pion1->p3(vtxp).unit());
  Hep3Vector unit2(kaon1->p3(vtxp).unit()); 
  Hep3Vector unit6(piSoft->p3(vtxp).unit()); 
  Hep3Vector unit7(D0->p3().unit()); 
  Hep3Vector unit4(vtxp.x(), vtxp.y(), vtxp.z()); 
  Hep3Vector unit4IP(vtxpIP.x(), vtxpIP.y(), vtxpIP.z()); 
  unit4 = unit4.unit();
  unit4IP = unit4IP.unit();
  
  double apertureKPi(acos(unit1.dot(unit2))); 
  double apertureD0PiSoft(acos(unit6.dot(unit7))); 
  
  _ntuple->column("Aper_Kpi", apertureKPi, -10.); 
  _ntuple->column("Aper_Dps", apertureD0PiSoft, -10.); 
  

  double deltaVtxD0(acos(unit3.dot(unit4)));
  double deltaVtxDstar(acos(unit5.dot(unit4)));
  Hep3Vector unit4d2(vtxp.x(), vtxp.y(), 0.); 
  Hep3Vector unit3d2(D0->p3().x(), D0->p3().y(), 0.); 
  Hep3Vector unit5d2(Dstar->p3().x(), Dstar->p3().y(), 0.); 
  unit4d2 = unit4d2.unit();
  unit3d2 = unit3d2.unit();
  unit5d2 = unit5d2.unit();
  double deltaVtx2dD0(acos(unit3d2.dot(unit4d2)));
  double deltaVtx2dDstar(acos(unit5d2.dot(unit4d2)));
  double deltaVtxIPD0(acos(unit3.dot(unit4IP)));
  double deltaVtxIPDstar(acos(unit5.dot(unit4IP)));
  
  _ntuple->column("DVtxD0", deltaVtxD0, -10.); 
  _ntuple->column("DVtx2dD0", deltaVtx2dD0, -10.); 
  _ntuple->column("DVtxIPD0", deltaVtxIPD0, -10.); 
  _ntuple->column("DtaVtxDs", deltaVtxDstar, -10.); 
  _ntuple->column("DVtx2dDs", deltaVtx2dDstar, -10.); 
  _ntuple->column("DVtxIPDs", deltaVtxIPDstar, -10.); 
  
  HepLorentzVector pPip4(0);
  pPip4 = kaon1->p4(vtxp);
  pPip4.boost(-(D0->p4().boostVector()));     
  double cosHelic;
  Hep3Vector pPip3(pPip4.x(),pPip4.y(),pPip4.z());
  pPip3 = pPip3.unit();
  cosHelic = pPip3.dot(unit3); 
  
  _ntuple->column("cosHelic", cosHelic, -10.); 

  _ntuple->dumpData();  
}    

void PacPidDstarSample::fillPidNtuple(const BtaCandidate* Dstar, const
				   BtaCandidate* D0,  const
				   BtaCandidate* piSoft,   const
				   BtaCandidate* pion1,   const
				   BtaCandidate* kaon1, const HepPoint
				   primVtx) const  
{
  // Removed to eliminate XxxData dependencies
  // See PacPidDstarSampleCCOLD for contents
}

double PacPidDstarSample::mass(int mode) const 
{ 
  double amass(0);
  switch(mode){
  case 1: 
    amass = _mass1;
    break;
  case 2:
    amass = _mass2;
    break;
  case 3: 
    amass =  _mass3;
    break;
  case 4: 
    amass =  _mass4;
    break;
  case 0: default:
    break;  
  }
  return amass;
}

bool PacPidDstarSample::constrDauMass(int mode) const
{ 
  bool constr(false);
  switch(mode){
  case 1: 
    constr =  _constrDauMassKaon;
    break;
  case 2:
    constr =  _constrDauMassPion;
    break;
  case 3: 
    constr =  _constrDauMassD0;
    break;
  case 4: 
    constr =  _constrDauMassSoftPi;
    break;
  case 0: default:
    break;  
  }
  return constr;
}

double PacPidDstarSample::constrainMass(const BtaCandidate* cand1,
					const BtaCandidate* cand2, 
					BtaCandidate* candMoth, int
					mode) const  
{
  
  double invmass=0;
  double e1=0;
  double e2=0;
  double p1=0;
  double p2=0;
  double mom = candMoth->p();
//   if( candMoth->decayVtx()!=NULL ){
//     HepPoint vtx(candMoth->decayVtx()->point());
//     //This should be at the vertex, but let's try and see if we can shut
//     //up the warnings! 
//     //    p1 = cand1->p(vtx);
//     //    p2 = cand2->p(vtx);
//     p1 = cand1->p();
//     p2 = cand2->p();
//   }
//   else {
//     p1 = cand1->p();
//     p2 = cand2->p();
//   }

  // PDS 28 Feb 2000
  // If we do not make a calculation of p at the vertex, then we can 
  // avoid making unecessary p calculations
  
  double massMode = mass(mode);
  
  double massModePlusOne = mass(mode+1);

  if ( massMode>=0 && constrDauMass(mode) ) { 
    p1 = cand1->p();
    e1 = sqrt(p1*p1+massMode*massMode); 
  }else { 
    e1 = cand1->energy(); 
  }
  
  if ( massModePlusOne>=0 && constrDauMass(mode+1) ) { 
    p2 = cand2->p();
    e2 = sqrt(p2*p2+massModePlusOne*massModePlusOne); 
  } else { 
    e2 = cand2->energy(); 
  }
  
  invmass = sqrt ( fabs(e1*e1+e2*e2-mom*mom +2*e1*e2));
  
  return invmass;
}

void PacPidDstarSample::setUp()
{
  
  _mass1 = Pdt::mass(PdtLund::K_plus);
  _mass2 = Pdt::mass(PdtLund::pi_plus);
  _mass3 = Pdt::mass(PdtLund::D0);
  _mass4 = Pdt::mass(PdtLund::pi_plus);
  
  _constrDauMassKaon = _pidDstarconstrDauMassKaon->value();
  _constrDauMassPion = _pidDstarconstrDauMassPion->value();
  _constrDauMassSoftPi = _pidDstarconstrDauMassSoftPi->value();
  _constrDauMassD0 = _pidDstarconstrDauMassD0->value();
  
  // mass cut:
  // if only one number is provided by the user, it is taken as a deltaM
  // w.r.t. the motherType mass
  
  if( _pidDstarMassCutD0Upper->value() > 0 ) {
    _cut1D0 = _pidDstarMassCutD0->value();
    _cut2D0 = _pidDstarMassCutD0Upper->value();
  }
  else
    if( _pidDstarMassCutD0->value() > 0 ) {
      _cut1D0 =  Pdt::mass(PdtLund::D0)-_pidDstarMassCutD0->value();
      _cut2D0 =  Pdt::mass(PdtLund::D0)+_pidDstarMassCutD0->value();
    }
  // loose mass cut
  if( _pidDstarLooseMassCutD0Upper->value() > 0 ) {
    _looseCut1D0 = _pidDstarLooseMassCutD0->value();
    _looseCut2D0 = _pidDstarLooseMassCutD0Upper->value();
  }
  else
    if( _pidDstarLooseMassCutD0->value() > 0 ) {
      _looseCut1D0 =  Pdt::mass(PdtLund::D0)-_pidDstarLooseMassCutD0->value();
      _looseCut2D0 =  Pdt::mass(PdtLund::D0)+_pidDstarLooseMassCutD0->value();
    }
  
  if( _pidDstarMassCutDstarUpper->value() > 0 ) {
    _cut1Dstar = _pidDstarMassCutDstar->value();
    _cut2Dstar = _pidDstarMassCutDstarUpper->value();
  }
  else
    if( _pidDstarMassCutDstar->value() > 0 ) {
      _cut1Dstar =  Pdt::mass(PdtLund::D_star_plus)-_pidDstarMassCutDstar->value();
      _cut2Dstar =  Pdt::mass(PdtLund::D_star_plus)+_pidDstarMassCutDstar->value();
    }
  
  if( _pidDstarMassCutDDUpper->value() > 0 ) {
    _cut1DD = _pidDstarMassCutDD->value();
    _cut2DD = _pidDstarMassCutDDUpper->value();
  }
  else
    if( _pidDstarMassCutDD->value() > 0 ) {
      _cut1DD =  .1455-_pidDstarMassCutDD->value();
      _cut2DD =  .1455+_pidDstarMassCutDD->value();
    }
  
  _cutCosThetaD0 = _pidDstarCosThetaD0Cut->value();
  _cutMomD0 = _pidDstarMomD0Cut->value();
  _cutD0Pion = _pidDstarD0PionCut->value();
  _cutD0Kaon = _pidDstarD0KaonCut->value();
  _cutDocaPiK = _pidDstarDocaPiKCut->value();
  _cutDocaPis = _pidDstarDocaPisCut->value();
  _cutAper_kpi = _pidDstarAper_kpiCut->value();
  _cutAper_dps = _pidDstarAper_dpsCut->value();
  _cutHelic = _pidDstarHelicCut->value();
  _maxPiSoftMom = _pidDstarMaxPiSoftMom->value();
  
  
  HepTupleManager* manager = gblEnv->getGen()->ntupleManager();
  
  // books an n-tuple
  //  if(  _pidDstarFillNtp->value() )
    _ntuple = manager->ntuple(factoryName());
}


bool 
PacPidDstarSample::passesTagSelection( const AbsEventTag *theTag, AbsEvent*){
  // NA
  // No tagbit is set in FastSim (so far)
  //bool aBool=false;
  ////  theTag->getBool( aBool, "isMultiHadron" );
  //bool anotherBool=false;
  //theTag->getBool( anotherBool, "BGFMultiHadron" ); 
  //return aBool || anotherBool;
  return( true );  
}


void PacPidDstarSample::createCandidateLists( AbsEvent* event )
{
  // Get the tag info for the event
  AbsEventTag* tag = Ifd<AbsEventTag>::get( event );
  assert(tag!=0);
  
  if (!passesTagSelection(tag, event)) return;

  BtaMcAssoc* truthMap(0);
  
  if ( _pidDstarCheckTruth->value()) {
    // get list of MC truth particles
    static const IfdStrKey keyDflt("Default");
    truthMap = Ifd< BtaMcAssoc >::get( event, keyDflt );
    assert (truthMap != 0);
  }     
  
  _truthMap = truthMap;
  
  // get the primary vtx from the event (only on Reco, so we'll get it
  // from the Tag
  //  HepAList< EventInfo >* infoList=NULL;
  //  getTmpAList(event, infoList, IfdStrKey("Default"));
  //  EventInfo* eventInfo = infoList->first();
  //  assert(eventInfo != 0);
  //
  //  const BbrPointErr primVtx = eventInfo->primaryVertex();  

   float x,y,z;
   if ( tag==0 ) {
   // maybe reco?
     static const IfdStrKey defaultKey("Default");
     HepAList< EventInfo >* infoList= Ifd<HepAList<EventInfo> >::get(event, defaultKey);
     assert(infoList != 0);
     EventInfo* eventInfo = infoList->first();
     assert(eventInfo != 0);
     
     const BbrPointErr primaryVtx = eventInfo->primaryVertex();  
     x = primaryVtx.x(); 
     y = primaryVtx.y(); 
     z = primaryVtx.z(); 
   }
   else {
     bool status(false);
     status = tag->getFloat( x, "xPrimaryVtx" );
     status = tag->getFloat( y, "yPrimaryVtx" );
     status = tag->getFloat( z, "zPrimaryVtx" );
   }
   HepPoint primVtx(x,y,z);
   
 //  BtaOpVtx makeVtxD;             // vertex operator for doing combinations
   BtaOpFastVtx makeVtxD; 
   
   const PdtEntry* thePdtD1;         
   const PdtEntry* thePdtD2;         
   const PdtEntry* thePdtD01;         
   const PdtEntry* thePdtD02;   

   if (_pidDstarUseCompOutput->value()) {  
     IfdStrKey pidDstarInputList1Key( _pidDstarInputList1->value() );
     HepAList<BtaCandidate>* candTmpList1 =
       Ifd< HepAList<BtaCandidate> >::get( event, pidDstarInputList1Key);
     
     BtaCandidate* cnd1;
     
     HepAList<BtaCandidate> chargedCandidateList1(*candTmpList1);
     
     
     assert(candTmpList1 != 0);
     HepAListIterator<BtaCandidate> iterCnd1( chargedCandidateList1);
     while ( cnd1 = iterCnd1() ) {
       
       HepAListIterator<BtaCandidate> dauIter(cnd1->daughterIterator());
       BtaCandidate* dau1 = dauIter.next();
       BtaCandidate* dau2 = dauIter.next();
       BtaCandidate* dauD0(0);
       BtaCandidate* dauSoftPi(0);
       thePdtD1 = dau1->pdtEntry();
      thePdtD2 = dau2->pdtEntry();
      if  (thePdtD1 == Pdt::lookup(PdtLund::anti_D0) ||
	   thePdtD1 == Pdt::lookup(PdtLund::D0)) {
	dauD0 = dau1;
	dauSoftPi = dau2;
      }
      else if  (thePdtD2 == Pdt::lookup(PdtLund::anti_D0) ||
		thePdtD2 == Pdt::lookup(PdtLund::D0)) {
	dauD0 = dau2;
	dauSoftPi = dau1;
      }
      if (dauD0->nDaughters() == 2) {  
	HepAListIterator<BtaCandidate> dauIter1(dauD0->daughterIterator());
	BtaCandidate* dauD01 = dauIter1.next();
	BtaCandidate* dauD02 = dauIter1.next(); 
	
	BtaCandidate* dauKaon(0);
	BtaCandidate* dauPion(0);
	thePdtD01 = dauD01->pdtEntry();
	thePdtD02 = dauD02->pdtEntry();
	if  (thePdtD01 == Pdt::lookup(PdtLund::K_plus) ||
	     thePdtD01 == Pdt::lookup(PdtLund::K_minus)) {
	  dauKaon = dauD01;
	  dauPion = dauD02;
	}
	else if  (thePdtD02 == Pdt::lookup(PdtLund::K_plus) ||
		  thePdtD02 == Pdt::lookup(PdtLund::K_minus)) {
	  dauKaon = dauD02;
	  dauPion = dauD01;
	}
       
	BtaCandidate* dauD0b =  makeVtxD.create(*dauKaon, *dauPion);;
	if ( acceptDstar(cnd1, dauD0b, dauSoftPi, dauKaon, dauPion) ){
	  if( _pidDstarFillNtp->value()) {
	    fillNtuple(cnd1, dauD0b, dauSoftPi, dauPion, dauKaon, primVtx);
	    // fillPidNtuple(cnd1, dauD0b, dauSoftPi, dauPion, dauKaon, primVtx);
	  }
	  const string dstarName( _PidDstar->value() );
	  selectThis( cnd1 , dstarName);
	  const string d0Name( _PidDstarD0->value() );
	  selectThis( dauD0 , d0Name);
	  const string softpionName( _mapPidDstarSoftPion->value() );
	  selectThis( dauSoftPi , softpionName);
	  const string kaonName( _mapPidDstarKaon->value() );
	  selectThis( dauKaon , kaonName);
	  const string pionName( _mapPidDstarPion->value() );
	  selectThis( dauPion , pionName);
	}
	delete dauD0b;
      }
    }
  }
  else {

  
    // get list(s) of input candidates
    IfdStrKey pidDstarInputList1Key( _pidDstarInputList1->value() );
    HepAList<BtaCandidate>* candList1 =
      Ifd< HepAList<BtaCandidate> >::get( event, pidDstarInputList1Key);
    
    assert(candList1 != 0);
    HepAList<BtaCandidate>* candList2;
    HepAList<BtaCandidate>* candList3;
    candList2 = candList1;
    candList3 = candList1;    

    // in this way there might be double counting.... need to fix
    // this. FIXED, I think...
    
    if( candList1->length() != 0 ){
      
      double vCharge = 0.;
      
      BtaCandidate* cnd1;
      BtaCandidate* cnd2;
      HepAListIterator<BtaCandidate> iterCnd1(*candList1);
      
      // loop over candidates in first list
      while ( cnd1 = iterCnd1() ) {
	// start another loop over the candidates
	HepAListIterator<BtaCandidate> iterCnd2(*candList2);
	while ( cnd2 = iterCnd2() ) {
	  
	  if ( (cnd1->charge())+(cnd2->charge()) != vCharge ) continue;
	  if ( cnd1->overlaps(*cnd2) ) continue;
	  
	  // do loose mass cut (note this is done _before_ vertexing)
	  BtaOpAdd4 d;
	  BtaCandidate* newTrk = 0;
	  newTrk = d.create( *cnd1, *cnd2 );
	  // cout<< "cnd1->mass() in event" <<cnd1->mass()  << endl;
	  double motherInvMass = newTrk->mass();
	  motherInvMass = constrainMass(cnd1, cnd2, newTrk, 1);
	  delete newTrk;
	  if ( motherInvMass < _looseCut1D0 || 
	       motherInvMass > _looseCut2D0 )  continue;
	  
	  combine(cnd1, cnd2 );
	  
	} // end of loop over candidates in list 2
      }   // end of loop over candidates in list 1  
      
    } //if( candList1->length ...

    HepAList<BtaCandidate>* candOutList1 = _outlisttmp;
    if( candList3->length() != 0 && candOutList1->length() != 0 ){
      
      const PdtEntry* thePdtD;         
      double chargePi;
      
      BtaCandidate* cnd1;
      BtaCandidate* cnd2; 
      HepAListIterator<BtaCandidate> iterCnd1(*candOutList1);
      
      // loop over candidates in first list
      
      while ( cnd1 = iterCnd1() ) {
	thePdtD = cnd1->pdtEntry();
	//check there is a D0 candidate
	if  (thePdtD != Pdt::lookup(PdtLund::anti_D0) &&
	     thePdtD != Pdt::lookup(PdtLund::D0)){
	  cerr << "Warning from D* Finder: no D0 send to this class" << endl;
	  break;
	} else {
	  //sign combination : a Pi+ and a D0, or a Pi- and a D0bar
	  if (thePdtD == Pdt::lookup(PdtLund::D0))
	    chargePi=1.;
	  else
	    chargePi=-1.;
	}
	
	// start another loop over the candidates (pions)
	HepAListIterator<BtaCandidate> iterCnd2(*candList3);	
	while ( cnd2 = iterCnd2() ) {
	  
	  //cut on piSoft momentum
	  if (cnd2->p() < _maxPiSoftMom){
	    //check if this track does not overlap with D daughters
	    if ( ! cnd2->overlaps(*cnd1) ){
	      if (cnd2->charge()==chargePi){
		// do loose mass cut (note this is done _before_ vertexing)
		BtaOpAdd4 d;
		BtaCandidate* newTrk = 0;
		newTrk = d.create( *cnd1, *cnd2 );
		HepAListIterator<BtaCandidate>
		  dauIter(cnd1->daughterIterator()); 
		BtaCandidate* dau1 = dauIter.next();
		BtaCandidate* dau2 = dauIter.next();
		if (dau1->charge() == cnd2->charge()) {			
		  if ( acceptDstar(newTrk, cnd1, cnd2, dau2, dau1) ){
		    if(  _pidDstarFillNtp->value()) {
		      fillNtuple(newTrk, cnd1, cnd2, dau1, dau2, primVtx);
		      //fillPidNtuple(newTrk, cnd1, cnd2, dau1, dau2, primVtx);
		    }
		    const string dstarName( _PidDstar->value() );
		    selectThis( newTrk , dstarName);
		    const string d0Name( _PidDstarD0->value() );
		    selectThis( cnd1, d0Name);
		    const string softpionName( _mapPidDstarSoftPion->value() );
		    selectThis( cnd2 , softpionName);
		    const string kaonName( _mapPidDstarKaon->value() );
		    selectThis( dau2 , kaonName);
		    const string pionName( _mapPidDstarPion->value() );
		    selectThis( dau1 , pionName);
		    _outlistdstartmp->append(newTrk);
		    //		    delete newTrk;
		  }
		  else {
		    delete newTrk;
		  } 
		} else {
			
		  if ( acceptDstar(newTrk, cnd1, cnd2, dau1, dau2) ){
		    if(  _pidDstarFillNtp->value()) {
		      fillNtuple(newTrk, cnd1, cnd2, dau2, dau1, primVtx);
		      //fillPidNtuple(newTrk, cnd1, cnd2, dau2, dau1, primVtx);
		    }
		    const string dstarName( _PidDstar->value() );
		    selectThis( newTrk , dstarName);
		    const string d0Name( _PidDstarD0->value() );
		    selectThis( cnd1, d0Name);
		    const string softpionName( _mapPidDstarSoftPion->value() );
		    selectThis( cnd2 , softpionName);
		    const string kaonName( _mapPidDstarKaon->value() );
		    selectThis( dau1 , kaonName);
		    const string pionName( _mapPidDstarPion->value() );
		    selectThis( dau2 , pionName);
		    _outlistdstartmp->append(newTrk);
		    //		    delete newTrk;
		  }
		  else {
		    delete newTrk;
		  }
		}
	      }
	    }
	  }
	}
      }
      HepAListDeleteAll(*_outlisttmp);
      HepAListDeleteAll(*_outlistdstartmp);
    }
  }
}

//--------------// Operations --
//--------------




