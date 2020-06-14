//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description:
//	Class BtaPCNtupleWriterMC - based on BetaUser examples.
//     Adapted from BetaPidCalibNtuple/BtaPCNtupleWriterMC.cc
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
#include "BaBar/BaBar.hh"
//-----------------------
// This Class's Header --
//-----------------------

#include "AbsParm/AbsParmIfdStrKey.hh"
#include "TrgTools/TrgFctTimePointInspector.hh"
#include "PacPidCalib/PacPidPCNtupleWriterMC.hh"
// also defines the class variables

#include <string>
using std::string;
#include <vector>
using std::vector;

//-------------
// C Headers --
//-------------
#include <assert.h>

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <math.h>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "AbsEnv/AbsEnv.hh"
#include "AbsEvent/AbsEvent.hh"
#include "Beta/BtaCandidate.hh"
#include "Beta/EventInfo.hh"
#include "BetaCoreTools/BtaMcAssoc.hh"
#include "CLHEP/Alist/AList.h"
#include "CLHEP/Alist/AIterator.h"
#include "ErrLogger/ErrLog.hh"
#include "AbsParm/AbsParmIfdStrKey.hh"
#include "Framework/AbsParmBool.hh"
#include "GenEnv/GenEnv.hh"
#include "HepTuple/Tuple.h"
#include "HepTuple/TupleManager.h"
#include "HepTuple/Histogram.h"
#include "ProxyDict/IfdStrKey.hh"
#include "ProxyDict/IfdKey.hh"
#include "ProxyDict/Ifd.hh"
#include "PDT/Pdt.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//----------------
// Constructors --
//----------------

// in general, a module constructor should not do much.  The beginJob
// member function is a better places to put initialization
//
// This ctor initializes the three sample list-name parameters
PacPidPCNtupleWriterMC::PacPidPCNtupleWriterMC( const char* const theName, 
			const char* const theDescription )
  : AppModule( theName, theDescription )
    ,_trkList( new AbsParmIfdStrKey("TracksList",this,"ChargedTracks"))
    ,_writeElectrons( new AbsParmBool("writeElectrons",this,true))
    ,_writePions    ( new AbsParmBool("writePions",    this,true))
    ,_writeKaons    ( new AbsParmBool("writeKaons",    this,true))
    ,_writeProtons  ( new AbsParmBool("writeProtons",  this,true))
    ,_writeMuons    ( new AbsParmBool("writeMuons",    this,true))
    , fElGlobalScale(1)
    , fPiGlobalScale(10)
    , fMuGlobalScale(1)
    , fKaGlobalScale(3)
    , fPrGlobalScale(1)
    , fElLowMomScale(1)  //was 2
    , fPiLowMomScale(1)  //was 4
    , fMuLowMomScale(1)  //was 2
    , fKaLowMomScale(1)  //was 2
    , fPrLowMomScale(1)  //was 2
    , fElGlobalCounter(0)
    , fPiGlobalCounter(0)
    , fMuGlobalCounter(0)
    , fKaGlobalCounter(0)
    , fPrGlobalCounter(0)
    , fElLowMomCounter(0)
    , fPiLowMomCounter(0)
    , fMuLowMomCounter(0)
    , fKaLowMomCounter(0)
    , fPrLowMomCounter(0)
{
  
  filler = new PacPidCand2Ntuple(this);
  _pidLists = new AbsParmVector<string>("PidLists",this);
  _tagBits = new AbsParmVector<string>("TagBits",this);
  
  //_pidLists = new AbsParmVector<string>("PidLists",this);
  //_tagBits = new AbsParmVector<string>("TagBits",this);
  
  //commands()->append( _pidLists );
  //commands()->append( _tagBits );
  
  commands()->append(_trkList );
  commands()->append(_writeElectrons);
  commands()->append(_writePions);
  commands()->append(_writeKaons);
  commands()->append(_writeProtons);
  commands()->append(_writeMuons);
  commands()->append(_pidLists);
  commands()->append(_tagBits);
 
}

//--------------
// Destructor --
//--------------

// The destructor should be limited to undoing the work of the constructor
PacPidPCNtupleWriterMC::~PacPidPCNtupleWriterMC( ) {
  delete _trkList;
  delete filler;
  delete _pidLists;
  delete _tagBits;
  delete _writeElectrons;
  delete _writePions;
  delete _writeKaons;
  delete _writeMuons;
  delete _writeProtons;
}

//--------------
// Operations --
//--------------

// The begin(AppJob*) member function is run before any events are
// processed.  In this analysis, it opens the output histogram file
// and then books a number of histograms and a ntuple.

AppResult
PacPidPCNtupleWriterMC::beginJob( AbsEvent* anEvent ) {
  ErrMsg(routine)<<"begin Job"<<endmsg; 

  HepTupleManager* manager = gblEnv->getGen()->ntupleManager();
  assert(manager != 0);

  // Names starting with underscores are defined in the
  // PacPidPCNtupleWriterMC.hh file as class variables. They live for the
  // life of the analysis object (this one), so can be used to communicate
  // between the begin, event and end functions. If you want to book another
  // here, make another entry in the obvious way in SampleBetaAnalysis.hh
    
  // book histograms  (arguments are nBins, low edge, high edge,
  //                   comment is the PAW histogram number)

  if (_writeElectrons->value()) {
    _eTuple  = manager->ntuple("MCElectrons", HepHistID(100));}
  if (_writePions->value()) {
    _piTuple = manager->ntuple("MCPions", HepHistID(101));}
  if (_writeMuons->value()) {
    _muTuple = manager->ntuple("MCMuons", HepHistID(102));}
  if (_writeKaons->value()) {
    _KTuple  = manager->ntuple("MCKaons", HepHistID(103));}
  if (_writeProtons->value()) {
    _pTuple  = manager->ntuple("MCProtons", HepHistID(104));}

  //_aHisto      =manager->histogram("MC reco abs mtm difference", 50, -0.1,  0.1 );

  // ntuple filler has tcl parameters too, so use them here
  filler->beginJob(anEvent);

  return AppResult::OK;
}

AppResult
PacPidPCNtupleWriterMC::event( AbsEvent* anEvent ) {
  // For convenience later, get a pointer to the event summary info object
  
  filler->setEvent(anEvent,_pidLists,_tagBits);
  
  HepAList<BtaCandidate> *tracks =   
    Ifd< HepAList<BtaCandidate> >::get( anEvent, _trkList->value());
  if (tracks !=0) {
    
    BtaCandidate *cand = 0;
    HepAListIterator<BtaCandidate> iter(*tracks);
    while (cand=iter()) {
      float plab = cand->p();
      BtaMcAssoc* _theTruth = Ifd< BtaMcAssoc >::get( anEvent, "GHit" );
      assert(_theTruth);
      BtaCandidate *mc;
      mc = _theTruth->mcFromReco(cand);
      
      if (mc!=0) {
	
	const PdtEntry* themcPdt = mc->pdtEntry();
	assert(themcPdt !=0);
	
	if (_writeElectrons->value()) {
	  if( themcPdt->lundId()==PdtLund::e_minus || themcPdt->lundId()==PdtLund::e_plus) {
	    if( fElGlobalCounter++%fElGlobalScale == 0 ) {
	      if( plab<.5 ) {
		if( fElLowMomCounter++%fElLowMomScale == 0 ) {
		  filler->convert(cand,anEvent,_eTuple);
		  _eTuple->dumpData();
		}
	      } else {
		filler->convert(cand,anEvent,_eTuple);
		_eTuple->dumpData();
	      }
	    }
	  }
	}
	if (_writePions->value()) {
	  if( themcPdt->lundId()==PdtLund::pi_minus || themcPdt->lundId()==PdtLund::pi_plus) {
	    if( fPiGlobalCounter++%fPiGlobalScale == 0 ) {
	      if( plab<.5 ) {
		if( fPiLowMomCounter++%fPiLowMomScale == 0 ) {
		  filler->convert(cand,anEvent,_piTuple);
		  _piTuple->dumpData();
		}
	      } else {
		filler->convert(cand,anEvent,_piTuple);
		_piTuple->dumpData();
	      }
	    }
	  }
	}
	if (_writeMuons->value()) {
	  if( themcPdt->lundId()==PdtLund::mu_minus || themcPdt->lundId()==PdtLund::mu_plus) {

	    if( fMuGlobalCounter++%fMuGlobalScale == 0 ) {
	      if( plab<.5 ) {
		if( fMuLowMomCounter++%fMuLowMomScale == 0 ) {
		  filler->convert(cand,anEvent,_muTuple);
		  _muTuple->dumpData();
		}
	      } else {
		filler->convert(cand,anEvent,_muTuple);
		_muTuple->dumpData();
	      }
	    }
	  }
	}
	if (_writeKaons->value()) {
	  if( themcPdt->lundId()==PdtLund::K_minus || themcPdt->lundId()==PdtLund::K_plus) {

	    if( fKaGlobalCounter++%fKaGlobalScale == 0 ) {
	      if( plab<.5 ) {
		if( fKaLowMomCounter++%fKaLowMomScale == 0 ) {
		  filler->convert(cand,anEvent,_KTuple);
		  _KTuple->dumpData();
		}
	      } else {
		filler->convert(cand,anEvent,_KTuple);
		_KTuple->dumpData();
	      }
	    }
	  }
	}
	if (_writeProtons->value()) {
	  if( themcPdt->lundId()==PdtLund::anti_p_minus || themcPdt->lundId()==PdtLund::p_plus) {
	    if( fPrGlobalCounter++%fPrGlobalScale == 0 ) {
	      if( plab<.5 ) {
		if( fPrLowMomCounter++%fPrLowMomScale == 0 ) {
		  filler->convert(cand,anEvent,_pTuple);
		  _pTuple->dumpData();
		}
	      } else {
		filler->convert(cand,anEvent,_pTuple);
		_pTuple->dumpData();
	      }
	    }
	  }
	}	
      }
    }
  }
  
  
  return AppResult::OK;
}

AppResult
PacPidPCNtupleWriterMC::endJob( AbsEvent* anEvent )
{
  ErrMsg(routine)<<" end Job" << endmsg;
  return AppResult::OK;
}

 
