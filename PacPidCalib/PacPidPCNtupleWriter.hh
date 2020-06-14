//--------------------------------------------------------------------------
//
// File and Version Information:
//  $Id: $
//
// Description:
//	Class BtaPCNtupleWriter - analysis class to produce
//	ntuples from _selected_ PacPidCalibSamples
//     Adapted from BetaPidCalibNtuple/BtaPCNtupleWriter.hh
//
//      More information can be found in PacPidPCNtupleWriter.cc
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Nicolas ARNAUD (SuperB)
//      Bob Jacobsen                    Original author
//      Thorsten Brandt
//
// Copyright Information:
//	(C) 2009 CNRS-IN2P3
//                                 
//------------------------------------------------------------------------

#ifndef PACPIDPCNTUPLEWRITER_HH
#define PACPIDPCNTUPLEWRITER_HH

#include <string>

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"
#include "Framework/AbsParmBool.hh"
#include "AbsParm/AbsParmVector.hh"
#include "AbsParm/AbsParmIfdStrKey.hh"
#include "PacPidCalib/PacPidCand2Ntuple.hh"
//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class HepTupleManager;
class HepHistogram;
class HepTuple;
class AbsParmString;
class PidElectronMicroSelector;
class ItgBinnedPdf;
class BtaPidElectronIdAlgo;
class PidLHElectronSelector;
class DtgDtoKpipiLikelihood;
class BtaCandidate;
#include <iosfwd>
//		---------------------
// 		-- Class Interface --
//		---------------------
 
class PacPidPCNtupleWriter : public AppModule {

//--------------------
// Instance Members --
//--------------------

public:

    // Constructors
    PacPidPCNtupleWriter( const char* const theName, const char* const theDescription );

    // Destructor
    virtual ~PacPidPCNtupleWriter( );

    // Operations

    virtual AppResult           beginJob( AbsEvent* anEvent );
    virtual AppResult event( AbsEvent* anEvent );
    virtual AppResult           endJob  ( AbsEvent* anEvent );

  void extractLambdas(AbsEvent *anEvent,const IfdKey &listName,HepTuple *dest);
  void extractDstar(AbsEvent *anEvent);
  void extractKs(AbsEvent  *anEvent);
  void extractFastKs(AbsEvent  *anEvent);
  void extract4e(AbsEvent *anEvent);
  void extractConversion(AbsEvent *anEvent);
  void extractTau31(AbsEvent *anEvent);
  void extractVcs(AbsEvent *anEvent);

  void extractJpsiKee(AbsEvent *anEvent);
  void extractJpsiKmm(AbsEvent *anEvent);
  void extractBtoDpi(AbsEvent *anEvent);

  void extractDtoKpipi(AbsEvent *anEvent);

  void extractLambdaCtopKpi(AbsEvent *anEvent);

  void processBhabhas(const HepAList<BtaCandidate> *elList, AbsEvent *anEvent, HepTuple *ntuple,
		      const HepAList<BtaCandidate> *otherList=NULL, bool doMuons=false);

  void addVtxQuantities(const BtaCandidate *cand,AbsEvent *anEvent,HepTuple *ntuple);
  bool isYloose(const BtaCandidate *c);
  bool isYtight(const BtaCandidate *c);

  void doPhi(AbsEvent*);

protected:
  PacPidCand2Ntuple* filler;

  AbsParmIfdStrKey _trkBhabhaList;
  AbsParmIfdStrKey _radBhabhaList;
  AbsParmIfdStrKey _ifrBhabhaList;
  AbsParmIfdStrKey _pidBhabhaList;
  AbsParmIfdStrKey _trkBhabhaListOther;
  AbsParmIfdStrKey _radBhabhaListOther;
  AbsParmIfdStrKey _ifrBhabhaListOther;
  AbsParmIfdStrKey _pidBhabhaListOther;	
  AbsParmIfdStrKey _4eList,_4eListO;
  AbsParmIfdStrKey _vcsListName,_vcsGammaListName;
  AbsParmIfdStrKey _lambdaList;
  AbsParmIfdStrKey _trklambdaList;
  AbsParmIfdStrKey _dstarList;
  AbsParmIfdStrKey _ksList;
  AbsParmIfdStrKey _fastKsList;
  AbsParmIfdStrKey _convElectronList,_convGammaList;
  AbsParmIfdStrKey _mumuGammaList,_mumuGammaListOther;
  AbsParmIfdStrKey _mumuGamma2List,_mumuGamma2ListOther;	
  AbsParmIfdStrKey _eemumuList,_eemumuListOther;
  AbsParmIfdStrKey _tau31piList,_tau31OneProngList;
  AbsParmBool _ifrNNKernelDBMode;
  AbsParmIfdStrKey _ifrNNKernelFile;

  AbsParmString *_tupleName;
  AbsParmVector<std::string> *_pidLists;
  AbsParmVector<std::string> *_tagBits;

  	
  AbsParmBool *_doTwoTrk;
  AbsParmBool *_doBhabha,*_do4e,*_doMumugamma, *_doEemumu,*_doVcs;
 
  AbsParmBool *_doMuha;
  AbsParmBool *_doKs,*_doFastKs,*_doDstar,*_doLambda,*_doTrkLambda,*_doConv,*_doPhi;

  AbsParmBool *_doJpsiK;
  AbsParmBool *_doBDpi;

  AbsParmBool *_doDKpipi;

  AbsParmBool *_doTau31;

  AbsParmBool *_doLambdaCpKpi;

  AbsEvent *_currentEvent;

  void printConfig();

  // define "local" variables to store from event to event
  // see the "begin" method for a discussion of these
  HepHistogram*    _numSampleHisto;
  
  HepTuple*       _lamTuple;
  HepTuple*       _trklamTuple;
  HepTuple*       _piDstarTuple;
  HepTuple*       _kDstarTuple;
  HepTuple*       _pisoftDstarTuple;
  HepTuple*       _piksNtuple;
  HepTuple*       _piFastKsNtuple;
  HepTuple*       _trkBhaNtuple;
  HepTuple*       _ifrBhaNtuple;
  HepTuple*       _pidBhaNtuple;
  HepTuple*       _EmcRadBhaNtuple;
  HepTuple*       _convNtuple;
  HepTuple*       _eeeeNtuple; 
  HepTuple*       _vcsTuple; 
  HepTuple*       _tau31Ntuple;
  HepTuple*       _mumugammaTuple;
  HepTuple*       _mumugamma2Tuple;	
  HepTuple*       _eemumuTuple;
  HepTuple*       _kPhiTuple;
  HepTuple*       _jpsikeeTuple;
  HepTuple*       _jpsikmmTuple;
  HepTuple*       _btodkaonTuple;
  HepTuple*       _btodpionTuple;
  HepTuple*       _dtokaonTuple;
  HepTuple*       _dtopion1Tuple;
  HepTuple*       _dtopion2Tuple;
  HepTuple*       _lambdactopkpiTuple;

  DtgDtoKpipiLikelihood* _dtagfitter;

  PidElectronMicroSelector *_selector;
  int _nPrintout;
  std::ofstream *_theFile;

// Sasha Telnov, 2007/08/11: cache for ::addVtxQuantities()
  const BtaCandidate *_vtxCandCache;
  double _vtxCandEnergyCache;
  float _vtxipangle2dCache;
  float _vtxalphaminCache;
  float _vtxaperCache;
  float _vtxipangleCache;
  float _vtxprobCache;
  float _vtxmassCache;
  float _vtxflightCache;
  float _vtxXYCache;
//float _vtxflightBsCache;
//float _vtxflightBsErrCache;
  float _vtxflightBsXYCache;
  float _vtxflightBsXYErrCache;
//float _vtxflightPrimCache;
//float _vtxflightPrimErrCache;
//float _vtxflightPrimXYCache;
//float _vtxflightPrimXYErrCache;

};

#endif

