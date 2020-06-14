//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description: Select very pure samples of Pions and Kaons from
// D*->D0+Pi_soft, D0->K Pi decays
// Adapted from BetaPidCalib/BtaPidDstarSample.hh
//       
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Nicolas ARNAUD (SuperB)
//      Giampiero Mancinelli                    Original author
// Copyright Information:
//	 (C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------

#ifndef PACPIDDSTARSAMPLE_HH
#define PACPIDDSTARSAMPLE_HH

#include "PacPidCalib/PacPidCalibSample.hh"
#include <string>


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class BtaCandidate;
class PdtEntry;
class EventInfo;
class HepPoint;
class HepTuple;
class BtaMcAssoc;
class AbsParmIfdStrKey;
template <class T> class HepAList;

//		---------------------
// 		-- Class Interface --
//		---------------------

class PacPidDstarSample : public PacPidCalibSample {
  
  //--------------------
  // Instance Members --
  //--------------------
  
public:
  
  // Constructors
  PacPidDstarSample( AppModule *) ;
  
  // Default need for Beta Factory stuff.  Does nothing 
  // except flag the error (in PacPidCalibSample)
  PacPidDstarSample(){};
  
  // Destructor
  virtual ~PacPidDstarSample( );
  
  virtual bool passesTagSelection( const AbsEventTag *, AbsEvent*);
  virtual void createCandidateLists( AbsEvent* );
  
  virtual void setUp();
protected:
  
  void combine(BtaCandidate*, BtaCandidate*);
  
  double constrainMass(const BtaCandidate*,const BtaCandidate*, 
		       BtaCandidate* candMoth, int mode) const;
  
  bool acceptDstar(BtaCandidate*, const BtaCandidate*, const
		   BtaCandidate*, const BtaCandidate*, const
		   BtaCandidate*) const; 
  bool acceptD0(const BtaCandidate*) const;
  
  bool constrDauMass(int mode) const;
  
  double mass(int mode) const;
  
  void fillNtuple(const BtaCandidate* Dstar, const BtaCandidate*
		  D0, const BtaCandidate* piSoft, const
		  BtaCandidate* pion, const BtaCandidate* kaon, 
		  const HepPoint primVtx) const; 
   void fillPidNtuple(const BtaCandidate* Dstar, const BtaCandidate*
		  D0, const BtaCandidate* piSoft, const
		  BtaCandidate* pion, const BtaCandidate* kaon, 
		  const HepPoint primVtx) const; 
 
  
private:
  
  HepTuple* _ntuple;
  BtaMcAssoc* _truthMap;
  
  double _looseCut1D0;
  double _looseCut2D0;
  
  
  double _cut1D0;  
  double _cut2D0;
  double _cutMomD0;
  double _cutCosThetaD0;
  double _cut1DD; 
  double _cut2DD;
  double _cut1Dstar;  
  double _cut2Dstar;
  
  double _cutD0Pion;
  double _cutD0Kaon;
  double _cutDocaPiK;
  double _cutDocaPis;
  double _cutAper_dps;
  double _cutAper_kpi;
  double _cutHelic;
  
  double _maxPiSoftMom;

  double _mass1; // mass of candidates to combine
  double _mass2; //  "
  double _mass3; //  " 
  double _mass4; //  "

  bool _constrDauMassKaon;
  bool _constrDauMassPion;
  bool _constrDauMassD0;
  bool _constrDauMassSoftPi;
  
  HepAList<BtaCandidate> *_outlisttmp;
  HepAList<BtaCandidate> *_outlistdstartmp;
  
  AbsParmIfdStrKey* _pidDstarEventInfoList;

  BtaParam<HepString>* _pidDstarInputList1 ;
  BtaParam<HepString>* _PidDstar ;
  BtaParam<HepString>* _PidDstarD0 ;
  BtaParam<HepString>* _mapPidDstarSoftPion ;
  BtaParam<HepString>* _mapPidDstarPion ;
  BtaParam<HepString>* _mapPidDstarKaon ;
  BtaParam<HepString>* _PidD0 ;
  BtaParam<double>* _pidDstarMassCutD0 ;
  BtaParam<double>* _pidDstarMassCutD0Upper ;
  BtaParam<double>* _pidDstarLooseMassCutD0 ;
  BtaParam<double>* _pidDstarLooseMassCutD0Upper ;
  BtaParam<double>* _pidDstarCosThetaD0Cut ;
  BtaParam<double>* _pidDstarMomD0Cut ;
  BtaParam<double>* _pidDstarMassCutDD ;
  BtaParam<double>* _pidDstarMassCutDDUpper ;
  BtaParam<double>* _pidDstarMassCutDstar ;
  BtaParam<double>* _pidDstarMassCutDstarUpper ;
  BtaParam<double>* _pidDstarD0PionCut ;
  BtaParam<double>* _pidDstarD0KaonCut ;
  BtaParam<double>* _pidDstarDocaPiKCut ;
  BtaParam<double>* _pidDstarDocaPisCut ;
  BtaParam<double>* _pidDstarAper_dpsCut ;
  BtaParam<double>* _pidDstarAper_kpiCut ;
  BtaParam<double>* _pidDstarHelicCut ;
  BtaParam<double>* _pidDstarMaxPiSoftMom ;
  BtaParam<bool>* _pidDstarconstrDauMassKaon ;
  BtaParam<bool>* _pidDstarconstrDauMassPion ;
  BtaParam<bool>* _pidDstarconstrDauMassD0 ;
  BtaParam<bool>* _pidDstarconstrDauMassSoftPi ;
  BtaParam<bool>* _pidDstarFillNtp ;
  BtaParam<bool>* _pidDstarCheckTruth ;
  BtaParam<bool>* _pidDstarUseCompOutput ;
 
public:
  
  static std::string factoryName() { return "PacPidDstarSample"; }
};

#endif

