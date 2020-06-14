//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description:
//	Class PacPidPCNtupleWriterMC - the barest outline of a Beta
//      Analysis, suitable for simple filling-in
//     Adapted from BetaPidCalibNtuple/BtaPCNtupleWriterMC.hh
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Nicolas ARNAUD (SuperB)
//      Bob Jacobsen                    Original author
//
// Copyright Information:
//	(C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------

#ifndef PACPIDPCNTUPLEWRITERMC_HH
#define PACPIDPCNTUPLEWRITERMC_HH

#include <string>
//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"
#include "PacPidCalib/PacPidCand2Ntuple.hh"
#include "AbsParm/AbsParmVector.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class HepTupleManager;
class HepHistogram;
class HepTuple;
class AbsParmIfdStrKey;
class AbsParmBool;

//		---------------------
// 		-- Class Interface --
//		---------------------
 
class PacPidPCNtupleWriterMC : public AppModule {

  //--------------------
  // Instance Members --
  //--------------------

public:

  // Constructors
  PacPidPCNtupleWriterMC( const char* const theName, const char* const theDescription );

  // Destructor
  virtual ~PacPidPCNtupleWriterMC( );

  // Operations

  virtual AppResult           beginJob( AbsEvent* anEvent );
  virtual AppResult           event   ( AbsEvent* anEvent );
  virtual AppResult           endJob  ( AbsEvent* anEvent );
    
protected:

  // The following are sample parameters

  PacPidCand2Ntuple* filler;
  AbsParmIfdStrKey* _trkList;
  HepTuple *_eTuple, *_piTuple, *_muTuple, *_KTuple, *_pTuple;
  AbsParmVector<std::string> *_pidLists;
  AbsParmVector<std::string> *_tagBits;
  AbsParmBool *_writeElectrons;
  AbsParmBool *_writePions;
  AbsParmBool *_writeKaons;
  AbsParmBool *_writeProtons;
  AbsParmBool *_writeMuons;
  
  int fElGlobalScale;
  int fPiGlobalScale;
  int fMuGlobalScale;
  int fKaGlobalScale;
  int fPrGlobalScale;

  int fElLowMomScale;
  int fPiLowMomScale;
  int fMuLowMomScale;
  int fKaLowMomScale;
  int fPrLowMomScale;

  int fElGlobalCounter;
  int fPiGlobalCounter;
  int fMuGlobalCounter;
  int fKaGlobalCounter;
  int fPrGlobalCounter;

  int fElLowMomCounter;
  int fPiLowMomCounter;
  int fMuLowMomCounter;
  int fKaLowMomCounter;
  int fPrLowMomCounter;

  //AbsParmIfdStrKey* _btaTrackList;
  //AbsParmIfdStrKey* _btaTruthMap;  
  //AbsParmBool *_loadMC;
  // define "local" variables to store from event to event
  // see the "begin" method for a discussion of these
  //HepHistogram*    _aHisto;
  
};

#endif

