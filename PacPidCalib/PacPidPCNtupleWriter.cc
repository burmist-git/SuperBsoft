//--------------------------------------------------------------------------
// File and Version Information:
//  $Id: $
//
// Description:
//  Class PacPidPCNtupleWriter -analysis class to produce
//  ntuples from _selected_ PacPidPidCalibSamples
//     Adapted from BetaPidCalibNtuple/BtaPCNtupleWriter.cc
//
//      With the help of "PacPidCand2Ntuple" it produces ntuples from
//      the following samples :
//          - Radiative Bhabhas , Conversions, eeee
//          - mu mu gamma , e e mu mu
//          - Dstar (K , pi) in LOOSE mode
//          - K0s (pi) in LOOSE mode
//          - tau 3-1
//          - Lambda (p) in LOOSE mode
//
//          + higher-momentum K0s from KsLoose - see "fastKs" (A. Telnov, Aug 2007)
//
// Environment:
//  Software developed for the Super B project
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
#include "BaBar/BaBar.hh"
//-----------------------
// This Class's Header --
//-----------------------
#include "TrgTools/TrgFctTimePointInspector.hh"
#include "AbsParm/AbsParmIfdStrKey.hh"
#include "PacPidCalib/PacPidPCNtupleWriter.hh"

// NA
//#include "BetaPid/PidElectronMicroSelector.hh"
//#include "BetaPid/PidLHRatios.hh"

// also defines the class variables

//-------------
// C Headers --
//-------------
#include <assert.h>

#include "ProbTools/probab.hh"
#include "HepTuple/HepHistID.hh"
//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>
using std::string;

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Alist/AList.h"
#include "CLHEP/Alist/AIterator.h"
#include "ErrLogger/ErrLog.hh"
#include "HepTuple/TupleManager.h"
#include "HepTuple/Tuple.h"
#include "HepTuple/Histogram.h"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEvent/getTmpAList.hh"
#include "BaBar/Constants.hh"
#include "Framework/AbsParmString.hh"
#include "ProbTools/probab.hh"
#include "CLHEP/Alist/AList.h"
#include "ProxyDict/IfdStrKey.hh"
#include "ProxyDict/IfdKey.hh"
#include "ProxyDict/Ifd.hh"
#include "AbsEnv/AbsEnv.hh"
#include "GenEnv/GenEnv.hh"
#include "PDT/Pdt.hh"
#include "PDT/PdtEntry.hh"
#include "ErrLogger/ErrLog.hh"

#include "Beta/EventInfo.hh"
#include "Beta/BtaCandidate.hh"
#include "Beta/BtaAbsVertex.hh"
#include "BetaMicroAdapter/BtaMicroAdapter.hh"
#include "BetaMicroAdapter/BtaTrkQual.hh"
#include "BetaMicroAdapter/BtaCalQual.hh"
#include "BetaMicroAdapter/BtaIfrQual.hh"
#include "BetaMicroAdapter/BtaPidQual.hh"
#include "BetaMicroAdapter/BtaPidInfo.hh"
#include "BetaCoreTools/BtaOpMakeTree.hh"
#include "BetaCoreTools/BtaBooster.hh"

#include "BetaCoreTools/BtaOpAdd4.hh"
#include "BetaCoreTools/BtaBVariables.hh"
#include "BetaCoreTools/BtaHelicity.hh"

// NA
//#include "BetaPidCalibNtuple/BtaPCJPsiHeli.hh"

#include "FastVtx/BtaOpFastVtx.hh"
#include "FastVtx/BtaOpFastVtxV0.hh"

#include "PepEnv/PepEnv.hh"
#include "PepEnvData/PepBeams.hh"

#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkFit.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkFitter/TrkPocaXY.hh"

#include "VtxFitter/VtxLeastChiVertexer.hh"
#include "VtxFitter/VtxGammaConv.hh"

#include "VtxFitter/VtxFitterOper.hh"
#include "VtxFitter/VtxLeastChiAlgorithm.hh"
#include "VtxFitter/VtxGeoKinAlgorithm.hh"
#include "VtxFitter/VtxAdd4Algorithm.hh"
#include "VtxTreeFitter/VtxTreeFitterAlgorithm.hh"
#include "VtxTreeFitter/VtkFitter.hh"
#include "VtxBase/BtaAbsVertexer.hh"

#include "AbsEventTag/AbsEventTag.hh"
#include "AbsEventTag/AbsEventTagBoolIter.hh"

#include "IfrPidUtils/IfrNNKernelSetup.hh"

#include "ProbTools/probab.hh"
#include "DTaggingTools/DtgDtoKpipiLikelihood.hh"

#include "BbrGeom/BbrPointErr.hh"
#include "BbrGeom/BbrVectorErr.hh"

#include "PacPid/PacPidUtils.hh"

using std::cout;
using std::cerr;
using std::endl;
using std::fstream;
using std::ios;
using std::setw;


//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//----------------
// Constructors --
//----------------

// in general, a module constructor should not do much.  The begin(job) or
// begin(run) members are better places to put initialization
PacPidPCNtupleWriter::PacPidPCNtupleWriter(
    const char* const theName,
    const char* const theDescription )
    : AppModule( theName, theDescription )
  ,_trkBhabhaList("TrkBhabhaList",this,"BtaTrkBhabhaSampleElectrons")
  ,_radBhabhaList("RadBhabhaList",this,"BtaEmcRadBhabhaSampleElectrons")
  ,_ifrBhabhaList("IfrBhabhaList",this,"BtaKinIFRBhabhaSampleElectrons")
  ,_pidBhabhaList("PidBhabhaList",this,"BtaPidBhabhaSampleElectrons")
  ,_trkBhabhaListOther("TrkBhabhaListOther",this,"BtaTrkBhabhaSampleElectronsOther")
  ,_radBhabhaListOther("RadBhabhaListOther",this,"BtaEmcRadBhabhaSampleElectronsOther")
  ,_ifrBhabhaListOther("IfrBhabhaListOther",this,"BtaKinIFRBhabhaSampleElectronsOther")
  ,_pidBhabhaListOther("PidBhabhaListOther",this,"BtaPidBhabhaSampleElectronsOther")
  ,_4eList("4eList",this,"BtaPideeeeCalibSampleElectrons")
  ,_4eListO("4eListOther",this,"BtaPideeeeCalibSampleOtherTrack")
  ,_vcsListName("vcsListName",this,"BtaVcsSampleElectrons")
  ,_vcsGammaListName("vcsGammaListName",this,"BtaVcsSamplePhotons")
  ,_lambdaList("LambdaList",this,"BtaPidLambdaProtonSamplePiConstrainedLamda")
  ,_trklambdaList("trkLambdaList",this,"BtaLambdaTrkOnlySampleLambda0")

  // NA
  ,_dstarList("DstarList",this,"PacPidDstarSampleDstar")

  ,_ksList("ksList",this,"BtaPidKsSamplePidKs")
  ,_fastKsList("fastKsList",this,"KsLoose")
  ,_convElectronList("convElectronList",this,"BtaPidGammaeeSampleElectrons")
  ,_convGammaList("convGammaList",this,"GammaConvCalib")
  ,_mumuGammaList("mumuGammaList",this,"BtamumugammaSampleMuons")
  ,_mumuGammaListOther("mumuGammaListOther",this,"BtamumugammaSampleMuonsOther")
  ,_mumuGamma2List("mumuGamma2List",this,"Btamumugamma2SampleMuons")
  ,_mumuGamma2ListOther("mumuGamma2ListOther",this,"Btamumugamma2SampleMuonsOther")
  ,_eemumuList("eemumuList",this,"BtaeemumuSampleMuons")
  ,_eemumuListOther("eemumuListOther",this,"BtaeemumuSampleMuonsOther")
  ,_tau31piList("tau31piList",this,"BtaTau31SampleTightPions")
  ,_tau31OneProngList("tau31OneProngList",this,"BtaTau31SampleSinglePr")
  ,_ifrNNKernelDBMode("ifrNNKernelDBMode",this,true)
  ,_ifrNNKernelFile("ifrNNKernelFile",this,"noFileNameYet")
{
  
  filler = new PacPidCand2Ntuple(this);
  _doTwoTrk = new AbsParmBool("doTwoTrk",this,false);
  _doBhabha = new AbsParmBool("doBhabha",this,false);
  _do4e = new AbsParmBool("do4e",this,false);
  _doVcs = new AbsParmBool("doVcs",this,false);
  _doMumugamma = new AbsParmBool("doMumugamma",this,false);
  _doEemumu = new AbsParmBool("doEemumu",this,false);

  _doMuha = new AbsParmBool("doMuha",this,false);
  _doTau31 = new AbsParmBool("doTau31",this,false);
  _doKs = new AbsParmBool("doKs",this,false);
  _doFastKs = new AbsParmBool("doFastKs",this,false);
  _doDstar = new AbsParmBool("doDstar",this,false);
  _doLambda = new AbsParmBool("doLambda",this,false);
  _doTrkLambda = new AbsParmBool("doTrkLambda",this,false);
  _doPhi = new AbsParmBool("doPhi",this,false);
  _doConv = new AbsParmBool("doConversions",this,false);
  _doJpsiK = new AbsParmBool("doJpsiK",this,false);
  _doBDpi = new AbsParmBool("doBDpi",this,false);
  
  _doDKpipi = new AbsParmBool("doDKpipi",this,false);
  _doLambdaCpKpi = new AbsParmBool("doLambdaCpKpi",this,false);
  _pidLists = new AbsParmVector<string>("PidLists",this);
  _tagBits = new AbsParmVector<string>("TagBits",this);
  
  commands()->append(&_ifrNNKernelDBMode);
  commands()->append(&_ifrNNKernelFile);
  
  commands()->append(_doTwoTrk);
  commands()->append(_doBhabha);
  commands()->append(_do4e);
  commands()->append(_doVcs);
  commands()->append(_doMumugamma);
  commands()->append(_doEemumu);
  
  commands()->append(_doMuha);
  commands()->append(_doKs);
  commands()->append(_doFastKs);
  commands()->append(_doDstar);
  commands()->append(_doLambda);
  commands()->append(_doTrkLambda);
  commands()->append(_doPhi);
  commands()->append(_doConv);
  
  commands()->append(_doJpsiK);
  commands()->append(_doBDpi);
  
  commands()->append(_doDKpipi);
  commands()->append(_doLambdaCpKpi);
  commands()->append(_doTau31);

  commands()->append(_pidLists);
  commands()->append(_tagBits);

  commands()->append(&_trkBhabhaList); commands()->append(&_trkBhabhaListOther);
  commands()->append(&_radBhabhaList); commands()->append(&_radBhabhaListOther);
  commands()->append(&_ifrBhabhaList); commands()->append(&_ifrBhabhaListOther);
  commands()->append(&_pidBhabhaList); commands()->append(&_pidBhabhaListOther);
  commands()->append(&_4eList);
  commands()->append(&_vcsListName); commands()->append(&_vcsGammaListName);
  commands()->append(&_lambdaList);
  commands()->append(&_trklambdaList);
  commands()->append(&_dstarList);
  commands()->append(&_ksList);
  commands()->append(&_fastKsList);
  commands()->append(&_convElectronList); commands()->append(&_convGammaList);
  commands()->append(&_mumuGammaList); commands()->append(&_mumuGammaListOther);
  commands()->append(&_mumuGamma2List); commands()->append(&_mumuGamma2ListOther);
  commands()->append(&_eemumuList); commands()->append(&_eemumuListOther);
  commands()->append(&_tau31piList); commands()->append(&_tau31OneProngList);
  _nPrintout=0;

  // NA
  //_selector = new PidElectronMicroSelector();
  
  _dtagfitter=0;

// initialization required only to avoid getting warnings from valgrind
  _vtxCandCache=0;
  _vtxCandEnergyCache=-8888.;
  _vtxipangle2dCache=-8888.;
  _vtxalphaminCache=-8888.;
  _vtxaperCache=-8888.;
  _vtxipangleCache=-8888.;
  _vtxprobCache=-8888.;
  _vtxmassCache=-8888.;
  _vtxflightCache=-8888.;
  _vtxXYCache=-8888.;
  _vtxflightBsXYCache=-8888.;
  _vtxflightBsXYErrCache=-8888.;
  
}

//--------------
// Destructor --
//--------------

// The destructor should be limited to undoing the work of the constructor
PacPidPCNtupleWriter::~PacPidPCNtupleWriter( )
{
  delete filler;
  delete _doTwoTrk;
  delete _doBhabha;
  delete _do4e;
  delete _doVcs;
  delete _doMumugamma;
  delete _doEemumu;
  delete _doMuha;
  delete _doTau31;
  delete _doJpsiK;
  delete _doBDpi;
  delete _doDKpipi;
  delete _doKs;
  delete _doFastKs;
  delete _doDstar;
  delete _doLambda;
  delete _doTrkLambda;
  delete _doLambdaCpKpi;
  delete _doPhi;
  delete _doConv;
  delete _pidLists;
  delete _tagBits;
  
  // NA
  //delete _selector;

  if (_dtagfitter != 0) delete _dtagfitter;
}

//--------------
// Operations --
//--------------

// The begin(AppJob*) member function is run before any events are
// processed.  In this analysis, it opens the output histogram file
// and then books a number of histograms and a ntuple.

AppResult PacPidPCNtupleWriter::beginJob( AbsEvent* anEvent )
{
  ErrMsg(routine)<<"begin Job"<<endmsg; 
  HepTupleManager* manager = gblEnv->getGen()->ntupleManager();
  assert(manager != 0);
  bool had = _doMuha->value();
  if (had || _doLambda->value()) { _lamTuple      = manager->ntuple("Protons from Lambda ",HepHistID(100)); }
  if (had || _doTrkLambda->value()) { _trklamTuple      = manager->ntuple("Protons from Lambda, kinematic cuts only ",HepHistID(105)); }
  if (had || _doDstar->value() ) { _piDstarTuple =  manager->ntuple("Pions from Dstar",HepHistID(101)); }
  if (had || _doDstar->value() ) { _kDstarTuple =  manager->ntuple("Kaons from Dstar",HepHistID(102)); }
  if (had || _doPhi->value() )   { _kPhiTuple =  manager->ntuple("Kaons from Phi",HepHistID(103)); }
  if (had || _doKs->value() )    { _piksNtuple =  manager->ntuple("Pions from K0s",HepHistID(104)); }
  if (had || _doFastKs->value() )    { _piFastKsNtuple =  manager->ntuple("Pions from fast K0s",HepHistID(114)); }
  if (had || _doConv->value())   { _convNtuple = manager->ntuple("conversion electrons",HepHistID(209)); }
  if (had || _doDstar->value() ) { _pisoftDstarTuple =  manager->ntuple("Soft pions from Dstar",HepHistID(106)); }
  
  bool tt = _doTwoTrk->value();
  if (tt || _doBhabha->value()) { _trkBhaNtuple =  manager->ntuple("Trk Bhabhas",HepHistID(200)); }
  if (tt || _doBhabha->value()) { _ifrBhaNtuple =  manager->ntuple("Ifr Bhabhas",HepHistID(201)); }
  if (tt || _doBhabha->value()) { _EmcRadBhaNtuple =  manager->ntuple("Emc RadBhabhas",HepHistID(202)); }
  if (tt || _doBhabha->value()) { _pidBhaNtuple =  manager->ntuple("Control Sample PidBhabhas",HepHistID(208)); }
  if (tt || _do4e->value())     { _eeeeNtuple = manager->ntuple("eeee electrons",HepHistID(203)); }
  if (tt || _doMumugamma->value()) {  _mumugammaTuple = manager->ntuple("mumugamma muons",HepHistID(204));
  _mumugamma2Tuple = manager->ntuple("mumugamma2 muons",HepHistID(207));}
  if (tt || _doEemumu->value()) {   _eemumuTuple = manager->ntuple("eemumu muons",HepHistID(205)); }
  if (_doVcs->value()) {   _vcsTuple = manager->ntuple("Vcs electrons",HepHistID(206)); }

  if (_doTau31->value()) { _tau31Ntuple = manager->ntuple("tau31 pions",HepHistID(300)); }

  if (_doJpsiK->value()) {
    _jpsikeeTuple =  manager->ntuple("electrons from JpsiK",HepHistID(901));
    _jpsikmmTuple =  manager->ntuple("muons from JpsiK",HepHistID(902));
  }
  if (_doBDpi->value()) {
    _btodkaonTuple =  manager->ntuple("kaons from BtoDpi",HepHistID(903));
    _btodpionTuple =  manager->ntuple("pions from BtoDpi",HepHistID(904));
  }
  if(_doDKpipi->value()) {
    _dtokaonTuple = manager->ntuple("kaons from D",HepHistID(401));
    _dtopion1Tuple = manager->ntuple("pion1 from D",HepHistID(402));
    _dtopion2Tuple = manager->ntuple("pion2 from D",HepHistID(403));
    _dtagfitter = new DtgDtoKpipiLikelihood();
  } else _dtagfitter=0;

  if(_doLambdaCpKpi->value()) {
    _lambdactopkpiTuple = manager->ntuple("protons from LambdaC", HepHistID(501));
  }
  
  printConfig();
  
// ntuple filler has tcl parameters too, so use them here
  filler->beginJob(anEvent);
  
  if ( !_ifrNNKernelDBMode.value() ) {
    IfrNNKernelSetup *nnKernel=new IfrNNKernelSetup(_ifrNNKernelFile.valueString().c_str(),0,_ifrNNKernelFile.valueString().c_str(),0,0);
    filler->setIfrNNKernel(nnKernel);
  }

// debugging
// cout << "List used for fast Ks: " << _fastKsList.value() << endl;
  
  return AppResult::OK;
}


AppResult PacPidPCNtupleWriter::event( AbsEvent* anEvent )
{
  _currentEvent = anEvent;
  filler->setEvent(anEvent,_pidLists,_tagBits);
  _nPrintout++;
  
  bool had = _doMuha->value();
  if (had || _doLambda->value())  { extractLambdas(anEvent,_lambdaList.value(),_lamTuple); }
  if (had || _doTrkLambda->value())  { extractLambdas(anEvent,_trklambdaList.value(),_trklamTuple); }
  if (had || _doDstar->value())   { extractDstar(anEvent); }
  if (had || _doKs->value())      { extractKs(anEvent); }
  if (had || _doFastKs->value())  { extractFastKs(anEvent); }
  if (had || _doPhi->value())     { doPhi(anEvent); }
  if (had || _doConv->value())    { extractConversion(anEvent); }
  
  bool tt = _doTwoTrk->value();
  if ( tt || _doBhabha->value()) {
  // TRK Bhabhas
    HepAList<BtaCandidate> *trkBha =
      Ifd< HepAList<BtaCandidate> >::get( anEvent, _trkBhabhaList.value());
    if (trkBha!=0) { processBhabhas(trkBha,anEvent,_trkBhaNtuple); }
    
  // IFR Bhabhas
    HepAList<BtaCandidate> *ifrBha =
      Ifd< HepAList<BtaCandidate> >::get( anEvent, _ifrBhabhaList.value());
    if (ifrBha!=0) { processBhabhas(ifrBha,anEvent,_ifrBhaNtuple); }
    
  // EMC Rad Bhabhas
    HepAList<BtaCandidate> *EmcRadBha =
      Ifd< HepAList<BtaCandidate> >::get( anEvent, _radBhabhaList.value());
    if (EmcRadBha!=0) { processBhabhas(EmcRadBha,anEvent,_EmcRadBhaNtuple); }
    
  // PID Bhabhas
    HepAList<BtaCandidate> *pidBha =
      Ifd< HepAList<BtaCandidate> >::get( anEvent, _pidBhabhaList.value());
    if (pidBha!=0) { processBhabhas(pidBha,anEvent,_pidBhaNtuple); }
  }
  if (tt || _do4e->value()) {  extract4e(anEvent); }
  if (_doVcs->value()) {  extractVcs(anEvent); }
  
// Mu Mu Gamma Muons
  if (tt || _doMumugamma->value()) {
    
    HepAList<BtaCandidate> *mumuGamma =
      Ifd< HepAList<BtaCandidate> >::get( anEvent,_mumuGammaList.value() );
    HepAList<BtaCandidate> *mumuGammaO =
      Ifd< HepAList<BtaCandidate> >::get( anEvent, _mumuGammaListOther.value() );
    if (mumuGamma!=0) { processBhabhas(mumuGamma,anEvent,_mumugammaTuple,mumuGammaO,true); }
    
    HepAList<BtaCandidate> *mumuGamma2 =
      Ifd< HepAList<BtaCandidate> >::get( anEvent,_mumuGamma2List.value() );
    HepAList<BtaCandidate> *mumuGamma2O =
      Ifd< HepAList<BtaCandidate> >::get( anEvent, _mumuGamma2ListOther.value() );
    if (mumuGamma2 !=0) {
      processBhabhas(mumuGamma2,anEvent,_mumugamma2Tuple,mumuGamma2O,true);
    }
  }
  
// e e mu mu Muons
  if (tt || _doEemumu->value()) {
    HepAList<BtaCandidate> *eemumu =
      Ifd< HepAList<BtaCandidate> >::get( anEvent, _eemumuList.value());
    HepAList<BtaCandidate> *eemumuO =
      Ifd< HepAList<BtaCandidate> >::get( anEvent, _eemumuListOther.value());
    if (eemumu!=0) { processBhabhas(eemumu,anEvent,_eemumuTuple,eemumuO,true); }
  }
  
  if (_doTau31->value()) {
    extractTau31(anEvent);
  }
  
  if (_doJpsiK->value())   {
    extractJpsiKee(anEvent);
    extractJpsiKmm(anEvent);
  }
  
  if (_doBDpi->value())   {
    extractBtoDpi(anEvent);
  }
  
  if (_doDKpipi->value()) {
    extractDtoKpipi(anEvent);
  }
  if(_doLambdaCpKpi->value()) {
    extractLambdaCtopKpi(anEvent);
  }
  return AppResult::OK;
}


void PacPidPCNtupleWriter::extractLambdas(AbsEvent  *anEvent, const IfdKey &listName,HepTuple *dest)
{
// pick up protons
  HepAList<BtaCandidate>* LambdaList;
  LambdaList =
  //Ifd< HepAList<BtaCandidate> >::get( anEvent, _lambdaList.value());
    Ifd< HepAList<BtaCandidate> >::get( anEvent, listName);
  if (LambdaList != 0){
    HepAListIterator<BtaCandidate> iter(*LambdaList);
    BtaCandidate *theLambda;
    while (theLambda = iter()) {
    // sort out  proton and pion
      HepAListIterator<BtaCandidate> dauIter = theLambda->daughterIterator();
      BtaCandidate *pi,*proton;
      pi =  dauIter();
      if (pi->pdtEntry()->pidId() != PdtPid::pion) {
	proton = pi; pi = dauIter();
	if (pi->pdtEntry()->pidId() != PdtPid::pion) {
	  ErrMsg(error) << "ERROR ! Lambda has no pi-daughter !! " << endmsg;
	  return ;
	}
      } else {
	proton = dauIter();
      }
      if (proton->pdtEntry()->pidId() != PdtPid::proton) {
	ErrMsg(error) << "ERROR ! Lambda has no proton daughter !! " << endmsg;
	return ;
      }
      
      addVtxQuantities(theLambda,anEvent,dest);
      
    // pion quantities
      const BtaPidInfo *pidInfo = pi->getMicroAdapter()->getPidInfo();
      const Consistency& theCons = pidInfo->consistency(Pdt::lookup(PdtPid::pion),PidSystem::dch);
      dest->column("picons",PacPidCand2Ntuple::packFloat((float) theCons.consistency(),14)); // flat from 0 to 1
      if (pi->getMicroAdapter()->getTrkQual()) {
	int nSvt = pi->getMicroAdapter()->getTrkQual()->nSvtHits();
	dest->column("pinsvt",nSvt);
      } else {
	dest->column("pinsvt",0);
      }
      dest->column("pmtm",PacPidCand2Ntuple::packFloat((float) pi->p(),10));
      filler->convert(proton,anEvent,dest);
      
    // CMS-Momentum of Lambda
      HepLorentzVector p4Lam(theLambda->p4());
      Hep3Vector  UpsiBoost;
      const PepBeams *theBeams = gblEnv->getPep()->pepBeams();
      if (theBeams != 0) {
	UpsiBoost = theBeams->boostCMtoLab();
      } else  {
	UpsiBoost = Hep3Vector(0.,0.,0.486976);
      }
      p4Lam.boost(-UpsiBoost);
      dest->column("lambdapcms",PacPidCand2Ntuple::packFloat((float) p4Lam.rho(),11));
      dest->column("lambdap",PacPidCand2Ntuple::packFloat((float) theLambda->p(),11));
    // K0s - mass
      HepLorentzVector vec_p = proton->p4(theLambda->decayVtx()->point());
      HepLorentzVector vec_pi = pi->p4(theLambda->decayVtx()->point());
      vec_p.setVectM(vec_p.vect(),0.139570);
      vec_pi.setVectM(vec_pi.vect(),0.139570);
      HepLorentzVector sum = vec_p+vec_pi;
      dest->column("vtxmassks",PacPidCand2Ntuple::packFloat((float) sum.mag(),11));
    // Dump Data
      dest->dumpData();
    }
  }
}


void PacPidPCNtupleWriter::extractDstar(AbsEvent  *anEvent)
{
  HepAList<BtaCandidate>* DstarList;
  DstarList =
    Ifd< HepAList<BtaCandidate> >::get( anEvent, _dstarList.value());
  if (DstarList==0) { 
    return; 
  }
  HepAListIterator<BtaCandidate>   candidateListIter(*DstarList);
// Sort out K,Pi and Pi (soft)
// (Gampiero's Code )
  BtaCandidate *Dstar;
  while (Dstar = candidateListIter()) {
    HepAListIterator<BtaCandidate> dauIter(Dstar->daughterIterator());
    BtaCandidate* dau1 = dauIter.next();
    BtaCandidate* dau2 = dauIter.next();
    BtaCandidate* D0(0);
    BtaCandidate* piSoft(0);
    BtaCandidate* dauD01(0);
    BtaCandidate* dauD02(0);
    BtaCandidate* dauPion(0);
    BtaCandidate* dauKaon(0);
    const PdtEntry*       thePdtD1 = dau1->pdtEntry();
    const PdtEntry*       thePdtD2 = dau2->pdtEntry();
    if  (thePdtD1 == Pdt::lookup(PdtLund::anti_D0) ||
	 thePdtD1 == Pdt::lookup(PdtLund::D0)) {
      D0 = dau1;
      HepAListIterator<BtaCandidate> dauIter1(D0->daughterIterator());
      dauD01 = dauIter1.next();
      dauD02 = dauIter1.next();
      piSoft = dau2;
      if (dau2->charge() == dauD01->charge()) {
	dauPion = dauD01;
	dauKaon = dauD02;
      }
      else {
	dauPion = dauD02;
	dauKaon = dauD01;
      }
    }
    else if  (thePdtD2 == Pdt::lookup(PdtLund::anti_D0) ||
	      thePdtD2 == Pdt::lookup(PdtLund::D0)) {
      D0 = dau2;
      HepAListIterator<BtaCandidate> dauIter1(D0->daughterIterator());
      dauD01 = dauIter1.next();
      dauD02 = dauIter1.next();
      piSoft = dau1;
      if (dau1->charge() == dauD01->charge()) {
	dauPion = dauD01;
	dauKaon = dauD02;
      }
      else {
	dauPion = dauD02;
	dauKaon = dauD01;
      }
    }
  // Now, we have pointers to the following BtaCands  :
  // DStar, D0, piSoft, dauPion, dauKaon;
  // First of all, put Kaon and Pion  tracks into ntuple
  //Hep3Vector vD0  = D0->p3();
    Hep3Vector vDstar  = Dstar->p3();
    Hep3Vector vdaupisoft  = piSoft->p3();
    Hep3Vector vdauKaon  = dauKaon->p3();
    Hep3Vector vdauPion  = dauPion->p3();
    
    
    filler->convert(dauKaon,anEvent,_kDstarTuple);
    filler->convert(dauPion,anEvent,_piDstarTuple);
    filler->convert(piSoft,anEvent,_pisoftDstarTuple);
    
  // Mass difference M(D0) - M(D*)
    double deltaM = Dstar->mass() - D0->mass();
    _kDstarTuple->column("vtxmassdif",PacPidCand2Ntuple::packFloat((float) deltaM,10)); // 0.0007 MeV resolution with a mean of 0.1455
    _piDstarTuple->column("vtxmassdif",PacPidCand2Ntuple::packFloat((float) deltaM,10));
    _pisoftDstarTuple->column("vtxmassdif",PacPidCand2Ntuple::packFloat((float) deltaM,10));
  // Mass  of Dstar
    _kDstarTuple->column("vtxmassDstar",PacPidCand2Ntuple::packFloat((float) Dstar->mass(),10));
    _piDstarTuple->column("vtxmassDstar",PacPidCand2Ntuple::packFloat((float) Dstar->mass(),10));
    _pisoftDstarTuple->column("vtxmassDstar",PacPidCand2Ntuple::packFloat((float) Dstar->mass(),10));
  // CMS - Momentum of Dstar
    HepLorentzVector p4Dstar(Dstar->p4());
    HepLorentzVector p4D0(D0->p4());
    Hep3Vector  UpsiBoost;
    const PepBeams *theBeams = gblEnv->getPep()->pepBeams();
    if (theBeams != 0) {
      UpsiBoost = theBeams->boostCMtoLab();
    } else  {
      UpsiBoost = Hep3Vector(0.,0.,0.486976);
    }
    p4Dstar.boost(-UpsiBoost);
    p4D0.boost(-UpsiBoost);
    _kDstarTuple->column("dstarpcms",PacPidCand2Ntuple::packFloat((float)  p4Dstar.rho(),14));
    _piDstarTuple->column("dstarpcms",PacPidCand2Ntuple::packFloat((float)  p4Dstar.rho(),14));
    _pisoftDstarTuple->column("dstarpcms",PacPidCand2Ntuple::packFloat((float)  p4Dstar.rho(),14));
    _kDstarTuple->column("dstarthetacms",PacPidCand2Ntuple::packFloat((float)  p4Dstar.theta(),14));
    _piDstarTuple->column("dstarthetacms",PacPidCand2Ntuple::packFloat((float)  p4Dstar.theta(),14));
    _pisoftDstarTuple->column("dstarthetacms",PacPidCand2Ntuple::packFloat((float)  p4Dstar.theta(),14));
    _kDstarTuple->column("d0pcms",PacPidCand2Ntuple::packFloat((float)  p4D0.rho(),14));
    _piDstarTuple->column("d0pcms",PacPidCand2Ntuple::packFloat((float)  p4D0.rho(),14));
    _pisoftDstarTuple->column("d0pcms",PacPidCand2Ntuple::packFloat((float)  p4D0.rho(),14));
    
  // Mass  of D0
    _kDstarTuple->column("vtxmassD0",PacPidCand2Ntuple::packFloat((float) D0->mass(),9));
    _piDstarTuple->column("vtxmassD0",PacPidCand2Ntuple::packFloat((float) D0->mass(),9));
    _pisoftDstarTuple->column("vtxmassD0",PacPidCand2Ntuple::packFloat((float) D0->mass(),9));
  // Momentum , cos(theta) of D0
    _kDstarTuple->column("d0mtm",PacPidCand2Ntuple::packFloat((float) D0->p(),14)); 
    _piDstarTuple->column("d0mtm",PacPidCand2Ntuple::packFloat((float) D0->p(),14));
    _pisoftDstarTuple->column("d0mtm",PacPidCand2Ntuple::packFloat((float) D0->p(),14));
    _kDstarTuple->column("d0cth",PacPidCand2Ntuple::packFloat((float) cos(D0->p3().theta()),14)); 
    _piDstarTuple->column("d0cth",PacPidCand2Ntuple::packFloat((float) cos(D0->p3().theta()),14));
    _pisoftDstarTuple->column("d0cth",PacPidCand2Ntuple::packFloat((float) cos(D0->p3().theta()),14));
  // aperture
    _kDstarTuple->column("vtxaperD0",PacPidCand2Ntuple::packFloat((float) vdauKaon.angle(vdauPion),14));
    _kDstarTuple->column("vtxaperDstar",PacPidCand2Ntuple::packFloat((float) vDstar.angle(vdaupisoft),14));
    _piDstarTuple->column("vtxaperD0",PacPidCand2Ntuple::packFloat((float) vdauKaon.angle(vdauPion),14));
    _piDstarTuple->column("vtxaperDstar",PacPidCand2Ntuple::packFloat((float) vDstar.angle(vdaupisoft),14));
    _pisoftDstarTuple->column("vtxaperD0",PacPidCand2Ntuple::packFloat((float) vdauKaon.angle(vdauPion),14));
    _pisoftDstarTuple->column("vtxaperDstar",PacPidCand2Ntuple::packFloat((float) vDstar.angle(vdaupisoft),14));
    
  // probaility of D0 - Vertex
    int ndof = D0->decayVtx()->nDof();
    double chi2=  D0->decayVtx()->chiSquared();
    double pchi2(probab(ndof,chi2));
    _kDstarTuple->column("vtxD0prob",PacPidCand2Ntuple::packFloat((float) pchi2,14));
    _piDstarTuple->column("vtxD0prob",PacPidCand2Ntuple::packFloat((float) pchi2,14));
    _pisoftDstarTuple->column("vtxD0prob",PacPidCand2Ntuple::packFloat((float) pchi2,14));
  // Soft pion, momentum
    _kDstarTuple->column("pisoftmtm",PacPidCand2Ntuple::packFloat((float) piSoft->p(),14));
    _piDstarTuple->column("pisoftmtm",PacPidCand2Ntuple::packFloat((float) piSoft->p(),14));
  // Soft pion, d0
    HepPoint vtxp  = D0->decayVtx()->point();
    if (piSoft->recoTrk()!=0) {
      TrkPocaXY pocaxy(piSoft->trkAbsFit()->traj(),0,vtxp);
      HepPoint ca  = piSoft->trkAbsFit()->position(pocaxy.fltl1());
      _kDstarTuple->column("pisoftd0",PacPidCand2Ntuple::packFloat((float) sqrt(ca.x()*ca.x()+ca.y()*ca.y()),14));
      _piDstarTuple->column("pisoftd0",PacPidCand2Ntuple::packFloat((float) sqrt(ca.x()*ca.x()+ca.y()*ca.y()),14));
      _pisoftDstarTuple->column("pisoftd0",PacPidCand2Ntuple::packFloat((float) sqrt(ca.x()*ca.x()+ca.y()*ca.y()),14));
    } else {
      _kDstarTuple->column("pisoftd0",(float) -1.0);
      _piDstarTuple->column("pisoftd0",(float) -1.0);
      _pisoftDstarTuple->column("pisoftd0",(float) -1.0);
    }
  // Soft pion, #svt/dch - hits
    const BtaTrkQual *tQual = piSoft->getMicroAdapter()->getTrkQual();
    int nsvt = 0;
    int ndch = 0;
    if (tQual) nsvt = tQual->nSvtHits();
    if (tQual) ndch = tQual->nDchHits();
    _kDstarTuple->column("pisoftnsvt",nsvt);
    _kDstarTuple->column("pisoftndch",ndch);
    _piDstarTuple->column("pisoftnsvt",nsvt);
    _piDstarTuple->column("pisoftndch",ndch);
  // Soft pion, consistencies
    const BtaPidInfo *pidInfo = piSoft->getMicroAdapter()->getPidInfo();
    const Consistency& theDchCons = pidInfo->consistency(Pdt::lookup(PdtPid::pion),PidSystem::dch);
    _kDstarTuple->column("pisoftdchcons",PacPidCand2Ntuple::packFloat((float) theDchCons.consistency(),14));
    _piDstarTuple->column("pisoftdchcons",PacPidCand2Ntuple::packFloat((float) theDchCons.consistency(),14));
    const Consistency& theSvtCons = pidInfo->consistency(Pdt::lookup(PdtPid::pion),PidSystem::svt);
    _kDstarTuple->column("pisoftsvtcons",PacPidCand2Ntuple::packFloat((float) theSvtCons.consistency(),14));
    _piDstarTuple->column("pisoftsvtcons",PacPidCand2Ntuple::packFloat((float) theSvtCons.consistency(),14));
  // Kaon - data for pi
    const BtaPidInfo *pidInfo2 = dauKaon->getMicroAdapter()->getPidInfo();
    const Consistency& theDchCons2 = pidInfo2->consistency(Pdt::lookup(PdtPid::kaon),PidSystem::dch);
    _piDstarTuple->column("kaondchcons",PacPidCand2Ntuple::packFloat((float) theDchCons2.consistency(),14));
    _piDstarTuple->column("kaonmtm",PacPidCand2Ntuple::packFloat((float) dauKaon->p(),14));
    _pisoftDstarTuple->column("kaondchcons",PacPidCand2Ntuple::packFloat((float) theDchCons2.consistency(),14));
    _pisoftDstarTuple->column("kaonmtm",PacPidCand2Ntuple::packFloat((float) dauKaon->p(),14));
    
//
    
    static const IfdStrKey keyGTVL("GoodTracksVeryLoose");
    HepAList< BtaCandidate >* gtvllist = Ifd<HepAList<BtaCandidate> >::get(anEvent, keyGTVL);
    HepAListIterator<BtaCandidate> gtvlIter(*gtvllist);
    BtaCandidate *thegtvl;
    bool found=false;
    while (thegtvl=gtvlIter()) {
      if (thegtvl->uid() == dauKaon->uid()) { found=true; }
    }
    _piDstarTuple->column("kaongtvl",found);
    _pisoftDstarTuple->column("kaongtvla",found);
    
    static const IfdStrKey keyGTL("GoodTracksLoose");
    HepAList< BtaCandidate >* gtllist = Ifd<HepAList<BtaCandidate> >::get(anEvent, keyGTL);
    HepAListIterator<BtaCandidate> gtlIter(*gtllist);
    BtaCandidate *thegtl;
    found=false;
    while (thegtl=gtlIter()) {
      if (thegtl->uid() == dauKaon->uid()) { found=true; }
    }
    _piDstarTuple->column("kaongtl",found);
    _pisoftDstarTuple->column("kaongtla",found);
    
    // NA
    //PidLHRatios* likeRatios=new PidLHRatios();
    double pslikeKvsPi,pslikeKvsPro,pslikeKvsEle,pslikeProvsPi,pslikeProvsEle,pslikePivsEle;
    //likeRatios->ratios(dauKaon,pslikeKvsPi,pslikeKvsPro,pslikeKvsEle,pslikeProvsPi,pslikeProvsEle,pslikePivsEle);
    bool isGoodDrc;
    computeDircLikelihoods( dauKaon, pslikeKvsPi, pslikeKvsPro, pslikeKvsEle,
			    pslikeProvsPi, pslikeProvsEle, pslikePivsEle,
			    isGoodDrc, false );
        
    static const IfdStrKey keyKLHNP("KLHNotPion");
    static const IfdStrKey keyKLHVL("KLHVeryLoose");
    static const IfdStrKey keyKLHL("KLHLoose");
    static const IfdStrKey keyKLHT("KLHTight");
    static const IfdStrKey keyKLHVT("KLHVeryTight");
    const vector<string> &keyList1 = _pidLists->value();
    for (int i=0; i<keyList1.size(); i++) {
      const string name = keyList1[i].c_str();
      HepAList< BtaCandidate >* pidlist = 0;
      if (name == "KLHNotPion") pidlist = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyKLHNP);
      if (name == "KLHVeryLoose") pidlist = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyKLHVL);
      if (name == "KLHLoose") pidlist = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyKLHL);
      if (name == "KLHTight") pidlist = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyKLHT);
      if (name == "KLHVeryTight") pidlist = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyKLHVT);
      if ( pidlist ) {
	HepAListIterator<BtaCandidate> pidIter(*pidlist);
	BtaCandidate *thepid;
	found = false;
	while (thepid=pidIter()) {
	  if (thepid->uid()==dauKaon->uid()) { found=true; }
	}
	const string kstring = string("kaon") + name;
	_piDstarTuple->column(kstring,found);
	const string pstring = string("kaon") + name + string("a");
	_pisoftDstarTuple->column(pstring,found);
      }
    }
    
    _piDstarTuple->column("kaonlikeKvsPi",PacPidCand2Ntuple::packFloat((float)pslikeKvsPi,13));
    _piDstarTuple->column("kaonlikeKvsPro",PacPidCand2Ntuple::packFloat((float)pslikeKvsPro,13));
    _piDstarTuple->column("kaonlikeKvsEle",PacPidCand2Ntuple::packFloat((float)pslikeKvsEle,13));
    _piDstarTuple->column("kaonlikeProvsPi",PacPidCand2Ntuple::packFloat((float)pslikeProvsPi,13));
    _piDstarTuple->column("kaonlikeProvsEle",PacPidCand2Ntuple::packFloat((float)pslikeProvsEle,13));
    _piDstarTuple->column("kaonlikePivsEle",PacPidCand2Ntuple::packFloat((float)pslikePivsEle,13));
    _pisoftDstarTuple->column("kaonlikeKvsPi",PacPidCand2Ntuple::packFloat((float)pslikeKvsPi,13));
    _pisoftDstarTuple->column("kaonlikeKvsPro",PacPidCand2Ntuple::packFloat((float)pslikeKvsPro,13));
    _pisoftDstarTuple->column("kaonlikeKvsEle",PacPidCand2Ntuple::packFloat((float)pslikeKvsEle,13));
    _pisoftDstarTuple->column("kaonlikeProvsPi",PacPidCand2Ntuple::packFloat((float)pslikeProvsPi,13));
    _pisoftDstarTuple->column("kaonlikeProvsEle",PacPidCand2Ntuple::packFloat((float)pslikeProvsEle,13));
    _pisoftDstarTuple->column("kaonlikePivsEle",PacPidCand2Ntuple::packFloat((float)pslikePivsEle,13));

  // Pion - data for k
    const BtaPidInfo *pidInfo3 = dauPion->getMicroAdapter()->getPidInfo();
    const Consistency& theDchCons3 = pidInfo3->consistency(Pdt::lookup(PdtPid::pion),PidSystem::dch);
    _kDstarTuple->column("piondchcons",PacPidCand2Ntuple::packFloat((float) theDchCons3.consistency(),13));
    _kDstarTuple->column("pionmtm",PacPidCand2Ntuple::packFloat((float) dauPion->p(),13));
    _pisoftDstarTuple->column("piondchcons",PacPidCand2Ntuple::packFloat((float) theDchCons3.consistency(),13));
    _pisoftDstarTuple->column("pionmtm",PacPidCand2Ntuple::packFloat((float) dauPion->p(),13));
    const BtaCalQual* cQual = dauPion->getMicroAdapter()->getCalQual();
    double erawpi=0;
    if (cQual) { erawpi =  cQual->rawEnergy(); }
    ndch=0;
    nsvt=0;
    tQual = dauPion->getMicroAdapter()->getTrkQual();
    if (tQual) ndch = tQual->nDchHits();
    if (tQual) nsvt = tQual->nSvtHits();
    _kDstarTuple->column("pioneop",PacPidCand2Ntuple::packFloat((float)(erawpi/dauPion->p()),13));
    _kDstarTuple->column("pionndch",ndch);
    _kDstarTuple->column("pionnsvt",nsvt);
    _pisoftDstarTuple->column("pioneop",PacPidCand2Ntuple::packFloat((float)(erawpi/dauPion->p()),13));
    _pisoftDstarTuple->column("pionndch",ndch);
    _pisoftDstarTuple->column("pionnsvt",nsvt);
    
    // NA
    //likeRatios->ratios(dauPion,pslikeKvsPi,pslikeKvsPro,pslikeKvsEle,pslikeProvsPi,pslikeProvsEle,pslikePivsEle);
    computeDircLikelihoods( dauKaon, pslikeKvsPi, pslikeKvsPro, pslikeKvsEle,
			    pslikeProvsPi, pslikeProvsEle, pslikePivsEle,
			    isGoodDrc, false );
    
    static const IfdStrKey keyPiLHVL("piLHVeryLoose");
    static const IfdStrKey keyPiLHL("piLHLoose");
    static const IfdStrKey keyPiLHT("piLHTight");
    static const IfdStrKey keyPiLHVT("piLHVeryTight");
    for (int i=0; i<keyList1.size(); i++) {
      const string name = keyList1[i].c_str();
      HepAList< BtaCandidate >* pidlist = 0;
      if (name == "piLHVeryLoose") pidlist = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyPiLHVL);
      if (name == "piLHLoose") pidlist = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyPiLHL);
      if (name == "piLHTight") pidlist = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyPiLHT);
      if (name == "piLHVeryTight") pidlist = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyPiLHVT);
      if ( pidlist ) {
	HepAListIterator<BtaCandidate> pidIter(*pidlist);
	BtaCandidate *thepid;
	found = false;
	while (thepid=pidIter()) {
	  if (thepid->uid()==dauPion->uid()) { found=true; }
	}
	const string kstring = string("pi") + name;
	_kDstarTuple->column(kstring,found);
      }
    }
    
    _kDstarTuple->column("pilikeKvsPi",PacPidCand2Ntuple::packFloat((float)pslikeKvsPi,13));
    _kDstarTuple->column("pilikeKvsPro",PacPidCand2Ntuple::packFloat((float)pslikeKvsPro,13));
    _kDstarTuple->column("pilikeKvsEle",PacPidCand2Ntuple::packFloat((float)pslikeKvsEle,13));
    _kDstarTuple->column("pilikeProvsPi",PacPidCand2Ntuple::packFloat((float)pslikeProvsPi,13));
    _kDstarTuple->column("pilikeProvsEle",PacPidCand2Ntuple::packFloat((float)pslikeProvsEle,13));
    _kDstarTuple->column("pilikePivsEle",PacPidCand2Ntuple::packFloat((float)pslikePivsEle,13));
    
    HepAList< BtaCandidate >* gtvl1list = Ifd<HepAList<BtaCandidate> >::get(anEvent, keyGTL);
    HepAListIterator<BtaCandidate> gtvl1Iter(*gtvl1list);
    BtaCandidate *thegtvl1;
    found=false;
    while (thegtvl1=gtvl1Iter()) {
      if (thegtvl1->uid() == dauPion->uid()) { found=true; }
    }
    _kDstarTuple->column("piongtvl",found);
    _pisoftDstarTuple->column("piongtvla",found);
    
    HepAList< BtaCandidate >* gtl1list = Ifd<HepAList<BtaCandidate> >::get(anEvent, keyGTL);
    HepAListIterator<BtaCandidate> gtl1Iter(*gtl1list);
    BtaCandidate *thegtl1;
    found=false;
    while (thegtl1=gtl1Iter()) {
      if (thegtl1->uid() == dauPion->uid()) { found=true; }
    }
    _kDstarTuple->column("piongtl",found);
    _pisoftDstarTuple->column("piongtla",found);
    
    _piDstarTuple->dumpData();
    _kDstarTuple->dumpData();
    _pisoftDstarTuple->dumpData();
    
    // NA
    //delete likeRatios;      
  }
} // end the D* sample

 
  
///////////////////////////////////////////////////////////////////////////////// 
void PacPidPCNtupleWriter::extractKs(AbsEvent  *anEvent)
{
  HepAList<BtaCandidate>* K0sList;
  K0sList =
    Ifd< HepAList<BtaCandidate> >::get( anEvent, _ksList.value());
  if (K0sList==0) { return; }
//cout << "n(k0s) =  " << K0sList->length() << endl;
  HepAListIterator<BtaCandidate>   candidateListIter(*K0sList);
  BtaCandidate *K0s;

  while (K0s = candidateListIter()) {
    HepAListIterator<BtaCandidate> dauIter(K0s->daughterIterator());
    BtaCandidate *pi1 = dauIter();
    BtaCandidate *pi2 = dauIter();
    if (pi1->p()<pi2->p()) {
      BtaCandidate *tmp = pi1; pi1=pi2;  pi2=tmp;
    }
// masses as Lambdas
    const HepPoint decayVtxKs = K0s->decayVtx()->point(); 
    HepLorentzVector vec_pi1 = pi1->p4(decayVtxKs);
    HepLorentzVector vec_pi2 = pi2->p4(decayVtxKs);
// could the faster track be a proton?
    vec_pi1.setVectM(vec_pi1.vect(),0.938272);
    vec_pi2.setVectM(vec_pi2.vect(),0.139570);
    HepLorentzVector sum = vec_pi1+vec_pi2;
    float vtxmasslam1 = PacPidCand2Ntuple::packFloat((float) sum.mag(),10);
// could the slower track be a proton?
    vec_pi1.setVectM(vec_pi1.vect(),0.139570);
    vec_pi2.setVectM(vec_pi2.vect(),0.938272);
    sum = vec_pi1+vec_pi2;
    float vtxmasslam2 = PacPidCand2Ntuple::packFloat((float) sum.mag(),10);

    int ii=1;
    filler->convert(pi1,anEvent,_piksNtuple);
    
    addVtxQuantities(K0s,anEvent,_piksNtuple);
    _piksNtuple->column("othermtm",PacPidCand2Ntuple::packFloat((float) pi2->p(),10));
    _piksNtuple->column("ishigh",ii);
    const BtaCalQual* cQual = pi2->getMicroAdapter()->getCalQual();
    if (cQual) { _piksNtuple->column("othereraw",PacPidCand2Ntuple::packFloat((float) cQual->rawEnergy(),10)); } 
    else { _piksNtuple->column("othereraw",(float) 0.0); }
    _piksNtuple->column("vtxmasslam",vtxmasslam1);
    _piksNtuple->column("vtxmasslamOther",vtxmasslam2);
    _piksNtuple->dumpData();
    
  // second pi is the slower one, regardless of charge
    ii = 0;
    filler->convert(pi2,anEvent,_piksNtuple);
    
    addVtxQuantities(K0s,anEvent,_piksNtuple);
    _piksNtuple->column("othermtm",PacPidCand2Ntuple::packFloat((float) pi1->p(),10));
    _piksNtuple->column("ishigh",ii);
    cQual = pi1->getMicroAdapter()->getCalQual();
    if (cQual) { _piksNtuple->column("othereraw",PacPidCand2Ntuple::packFloat((float) cQual->rawEnergy(),10)); } 
    else { _piksNtuple->column("othereraw",(float) 0.0); }
    _piksNtuple->column("vtxmasslam",vtxmasslam2);
    _piksNtuple->column("vtxmasslamOther",vtxmasslam1);
    _piksNtuple->dumpData();
  }
}


///////////////////////////////////////////////////////////////////////////////// 
// Sasha Telnov, 2007/08/10
void PacPidPCNtupleWriter::extractFastKs(AbsEvent  *anEvent)
{
  HepAList<BtaCandidate>* K0sList;
  K0sList = Ifd< HepAList<BtaCandidate> >::get( anEvent, _fastKsList.value());
  if (K0sList==0) { 
  //  cout << "Empty Ks list " << _fastKsList.value() << endl;
    return; }
//cout << "n(k0s) =  " << K0sList->length() << endl;
  HepAListIterator<BtaCandidate> candidateListIter(*K0sList);
  BtaCandidate *K0s;
  while (K0s = candidateListIter()) {
  //if (K0s->p() < 1.5) continue;
    HepAListIterator<BtaCandidate> dauIter(K0s->daughterIterator());
    BtaCandidate *pi1 = dauIter();
    BtaCandidate *pi2 = dauIter();
    if (pi1->p()<pi2->p()) {
      BtaCandidate *tmp = pi1; pi1=pi2;  pi2=tmp;
    }
  // veto on cosmics and junk
    if (pi1->p() > 9.0) continue; 

  // first pi is the faster one, regardless of charge
    int ii=1;
  // Remove the most eggregiously bad Ks candidates
    addVtxQuantities(K0s,anEvent,_piFastKsNtuple);

  // remove poorly reconstructed vertices
    if (_vtxprobCache<0.1) continue;
  // 10 sigma significance cut on the XY flight length w.r.t. the beam spot
    if (_vtxflightBsXYCache/_vtxflightBsXYErrCache < 10.0) continue;
  // eliminate Ks and junk that come from the beampipe
    if (_vtxXYCache>2.35 && _vtxXYCache<2.9) continue; 

  // flight direction veto, 3-D. angle<0.1 = cos>0.995. For fast Ks, angle<0.03 would be OK
    if (_vtxipangleCache>0.03) continue;
  // flight direction veto, 2-D. angle<0.1 = cos>0.995. For fast Ks, angle<0.02 would be OK
    if (_vtxipangle2dCache>0.02) continue; // that's right, apply both

  // do keep a loose aperture cut - confirmed harmless to real Ks
    if (_vtxaperCache < 0.11) continue;
  // a sanity check
    if (_vtxXYCache>60.0 || _vtxflightCache>120.0) continue;

// masses as Lambdas
    const HepPoint decayVtxKs = K0s->decayVtx()->point(); 
    HepLorentzVector vec_pi1 = pi1->p4(decayVtxKs);
    HepLorentzVector vec_pi2 = pi2->p4(decayVtxKs);
// could the faster track be a proton?
    vec_pi1.setVectM(vec_pi1.vect(),0.938272);
    vec_pi2.setVectM(vec_pi2.vect(),0.139570);
    HepLorentzVector sum = vec_pi1+vec_pi2;
    float vtxmasslam1 = PacPidCand2Ntuple::packFloat((float) sum.mag(),10);
// could the slower track be a proton?
    vec_pi1.setVectM(vec_pi1.vect(),0.139570);
    vec_pi2.setVectM(vec_pi2.vect(),0.938272);
    sum = vec_pi1+vec_pi2;
    float vtxmasslam2 = PacPidCand2Ntuple::packFloat((float) sum.mag(),10);

// Lambda mass veto - also, turns out to be an effective aperture cut that removes most vtxaper < ~0.1
// Don't bother preserving the few good Ks that are below the Lambda mass peak
    if (vtxmasslam1<1.122 || vtxmasslam2<1.122) continue;

    _piFastKsNtuple->column("othermtm",PacPidCand2Ntuple::packFloat((float) pi2->p(),10));
    _piFastKsNtuple->column("ishigh",ii);
    const BtaCalQual* cQual = pi2->getMicroAdapter()->getCalQual();
    if (cQual) { _piFastKsNtuple->column("othereraw",PacPidCand2Ntuple::packFloat((float) cQual->rawEnergy(),10)); } 
    else { _piFastKsNtuple->column("othereraw",(float) 0.0); }
    _piFastKsNtuple->column("vtxmasslam",vtxmasslam1);
    _piFastKsNtuple->column("vtxmasslamOther",vtxmasslam2);
    filler->convert(pi1,anEvent,_piFastKsNtuple);
    _piFastKsNtuple->dumpData();
    
  // second pi is the slower one, regardless of charge
    ii = 0;    
    addVtxQuantities(K0s,anEvent,_piFastKsNtuple);
    _piFastKsNtuple->column("othermtm",PacPidCand2Ntuple::packFloat((float) pi1->p(),10));
    _piFastKsNtuple->column("ishigh",ii);
    cQual = pi1->getMicroAdapter()->getCalQual();
    if (cQual) { _piFastKsNtuple->column("othereraw",PacPidCand2Ntuple::packFloat((float) cQual->rawEnergy(),10)); } 
    else { _piFastKsNtuple->column("othereraw",(float) 0.0); }
    _piFastKsNtuple->column("vtxmasslam",vtxmasslam2);
    _piFastKsNtuple->column("vtxmasslamOther",vtxmasslam1);
    filler->convert(pi2,anEvent,_piFastKsNtuple);    
    _piFastKsNtuple->dumpData();
  }
}


//////////////////////////////////////////////////////////////////////////
// Sasha, 2007/08/11: some of the PID samples call addVtxQuantities
// 2 or 3 times for the same exact mother candidate. To avoid replication
// of identical calculation, add some primitive caching. Cached values
// are used if addVtxQuantities is called with the same cand pointer
// and the same cand energy as in the previous call (comparing energies is easier than momenta).
//   Also, add flight-length significance in 3-D and 2-D
void PacPidPCNtupleWriter::addVtxQuantities(const BtaCandidate *cand,AbsEvent *anEvent,HepTuple *ntuple)
{
  if (!(_vtxCandCache == cand && _vtxCandEnergyCache == cand->energy())) { // otherwise, used cached quantities
    _vtxCandCache = cand; // cache the pointer to the composite candidate
    _vtxCandEnergyCache = cand->energy(); // cache its energy
    HepAListIterator<BtaCandidate> dauIter(cand->daughterIterator());
    BtaCandidate *dau1 = dauIter();
    BtaCandidate *dau2 = dauIter(); 

// Get Primary Vertex from Tag
    float p1=0, p2=0, p2_bs=0, p3=0,aFloat=0;
    AbsEventTag* tag = Ifd<AbsEventTag>::get( anEvent ) ;    
    if ( tag->getFloat( aFloat, "xPrimaryVtx" ) ) { p1 = aFloat; }
    if ( tag->getFloat( aFloat, "yPrimaryVtx" ) ) { p2 = aFloat; }
    if ( tag->getFloat( aFloat, "zPrimaryVtx" ) ) { p3 = aFloat; }

  // To do flight significance calculation w.r.t. the beamspot, we need to get the real beam spot - with errors
    static const IfdStrKey keyDflt("Default");
  //getTmpAList(_currentEvent, infoList, keyDflt);
    HepAList< EventInfo >* infoList = Ifd<HepAList<EventInfo> >::get(_currentEvent,keyDflt);
    if (infoList == NULL) ErrMsg(fatal) << "No Default EventInfo in the event" << endmsg;

    EventInfo* eventInfo = infoList->first();
    const BbrPointErr beamSpotErr = eventInfo->beamSpot(); // beamSpotErr is beamspot with error matrix
  //BbrPointErr primVtxErr = beamSpotErr; // primVtxErr is primary vertex with error matrix
  //if (eventInfo->primaryVtx() != NULL) // fall back to beamSpotErr if primary vertex does not exist
  //  primVtxErr = eventInfo->primaryVertex(); 
   
    p2_bs = beamSpotErr.y();    

    const HepPoint decayVtxPoint = cand->decayVtx()->point(); // V0 decay point
    BbrPointErr decVtxErr(decayVtxPoint,cand->decayVtx()->xxCov()); // V0 decay point with error matrix
    BbrVectorErr flightVecErrBs(decVtxErr - beamSpotErr); // flight vector with error matrix w.r.t. the beamspot
  //BbrVectorErr flightVecErrPrim(decVtxErr - primVtxErr); // flight vector with error matrix w.r.t. the primary vertex
  // 3-D flight distance w.r.t. the beamspot
  //_vtxflightBsCache = PacPidCand2Ntuple::packFloat( flightVecErrBs.mag(), 14); 
  // 2-D flight distance w.r.t. the beamspot
    _vtxflightBsXYCache = PacPidCand2Ntuple::packFloat( flightVecErrBs.perp(), 12); 
  // error on the 3-D flight distance w.r.t. the beamspot
  //_vtxflightBsErrCache = PacPidCand2Ntuple::packFloat(sqrt( flightVecErrBs.covRTPMatrix()[BbrVectorErr::Rho][BbrVectorErr::Rho]), 16);
  // error on the 2-D flight distance w.r.t. the beamspot
    _vtxflightBsXYErrCache = PacPidCand2Ntuple::packFloat(sqrt( flightVecErrBs.covRZPMatrix()[BbrVectorErr::C_Rho][BbrVectorErr::C_Rho]), 14);
  // 3-D flight distance w.r.t. the primary vertex
  //_vtxflightPrimCache = PacPidCand2Ntuple::packFloat( flightVecErrPrim.mag(), 14); 
  // 2-D flight distance w.r.t. the primary vertex
  //_vtxflightPrimXYCache = PacPidCand2Ntuple::packFloat( flightVecErrPrim.perp(), 12); 
  // error on the 3-D flight distance w.r.t. the primary vertex 
  //_vtxflightPrimErrCache = PacPidCand2Ntuple::packFloat( sqrt(flightVecErrPrim.covRTPMatrix()[BbrVectorErr::Rho][BbrVectorErr::Rho]), 16);
  // error on the 2-D flight distance w.r.t. the primary vertex
  //_vtxflightPrimXYErrCache = PacPidCand2Ntuple::packFloat( sqrt(flightVecErrPrim.covRZPMatrix()[BbrVectorErr::C_Rho][BbrVectorErr::C_Rho]), 14);

  //BbrVectorErr momentumErr(cand->p3WCov()); // V0 lab momentum with error matrix


  // Now, back to the original calculation
  //Hep3Vector primVtx(p1,p2,p3);
  // y dimension of the beam spot is known much better than the y position of the primary vertex
    Hep3Vector primVtxBs(p1,p2_bs,p3);

// calculate vtx-quantities
    Hep3Vector decvtx( decayVtxPoint.x(),
		       decayVtxPoint.y(),
		       decayVtxPoint.z());
//Hep3Vector dif = decvtx  -  primVtx;
    Hep3Vector difBs = decvtx  -  primVtxBs;
    
    int ndof = cand->decayVtx()->nDof();
    double chi2=  cand->decayVtx()->chiSquared();
    double pchi2(probab(ndof,chi2));
    Hep3Vector vec1 = dau1->p3(decayVtxPoint);
    Hep3Vector vec2 = dau2->p3(decayVtxPoint);
// 2d  angle ---------------------------
    Hep3Vector u( difBs.unit() );
    Hep3Vector v( cand->p3().unit() );
    HepVector u2d(2);
    u2d[0]=u.x();u2d[1]=u.y();
    HepVector v2d(2);
    v2d[0]=v.x();v2d[1]=v.y();
    double ptot2 = u2d.normsq()*v2d.normsq();
    
    if(ptot2 <= 0) {
      _vtxipangle2dCache = 0;
    } else {
      HepDouble arg = dot(u2d,v2d)/sqrt(ptot2);
      if(arg >  1.0) arg =  1.0;
      if(arg < -1.0) arg = -1.0;
      _vtxipangle2dCache = PacPidCand2Ntuple::packFloat((float) acos(arg),14);
    }

// Min. Angle daugther <-> mother
    double a1 = v.angle(vec1);
    double a2 = v.angle(vec2);
    if (a1<a2) {
      _vtxalphaminCache = PacPidCand2Ntuple::packFloat((float) a1,14);
    } else {
      _vtxalphaminCache = PacPidCand2Ntuple::packFloat((float) a2,14);
    }

    _vtxaperCache = PacPidCand2Ntuple::packFloat((float) vec1.angle(vec2),14);
    _vtxipangleCache = PacPidCand2Ntuple::packFloat((float) difBs.angle(cand->p3()),14);
    _vtxprobCache = PacPidCand2Ntuple::packFloat((float) pchi2,14);
    _vtxmassCache = PacPidCand2Ntuple::packFloat((float) cand->p4().mag(),8);
    _vtxflightCache = PacPidCand2Ntuple::packFloat((float) difBs.mag(),14);
    _vtxXYCache = PacPidCand2Ntuple::packFloat((float) sqrt(decvtx.x()*decvtx.x() + decvtx.y()*decvtx.y()) ,12);
    
  } // end recalculating vertex quantitites

  ntuple->column("vtxipangle2d",_vtxipangle2dCache);
  ntuple->column("vtxalphamin",_vtxalphaminCache);
  ntuple->column("vtxaper",_vtxaperCache);
  ntuple->column("vtxipangle",_vtxipangleCache);
  ntuple->column("vtxchi2prob",_vtxprobCache);
  ntuple->column("vtxmass",_vtxmassCache);
  ntuple->column("vtxflight",_vtxflightCache);
  ntuple->column("vtxXY",_vtxXYCache);

//ntuple->column("vtxflightBs",_vtxflightBsCache);
//ntuple->column("vtxflightBsErr",_vtxflightBsErrCache);
  ntuple->column("vtxflightBsXY",_vtxflightBsXYCache);
  ntuple->column("vtxflightBsXYErr",_vtxflightBsXYErrCache);
//ntuple->column("vtxflightPrim",_vtxflightPrimCache);
//ntuple->column("vtxflightPrimErr",_vtxflightPrimErrCache);
//ntuple->column("vtxflightPrimXY",_vtxflightPrimXYCache);
//ntuple->column("vtxflightPrimXYErr",_vtxflightPrimXYErrCache);

} // end addVtxQuantities



////////////////////////////////////////////////////////////////
void PacPidPCNtupleWriter::extract4e(AbsEvent *anEvent)
{
  HepAList<BtaCandidate>* eeeeList =
    Ifd< HepAList<BtaCandidate> >::get( anEvent, _4eList.value());
  HepAList<BtaCandidate>* otherList =
    Ifd< HepAList<BtaCandidate> >::get( anEvent, _4eListO.value());;
  
  if ((eeeeList==0) || (otherList==NULL)) { return ; }
  if (eeeeList->length() != otherList->length()) {
    ErrMsg(warning) << "WARNING !! eeeeList->length() != otherList->length()" << endmsg;
    return;
  }
  
  for (int i=0; i<eeeeList->length(); i++) {
// Sasha Telnov, 2007/08/10:
// clone the electron candidate so that its type can be set to electron, 
// so that its momentum at IP can be taken from the electron-hypothesis fit
    BtaCandidate *elec = new BtaCandidate(*((*eeeeList)[i]));
    elec->setType(Pdt::lookup(PdtPid::electron, (int)elec->charge()));
// Cloning BtaCandidates is expensive, so don't bother doing so with the "other" cand
    BtaCandidate *other = (*otherList)[i];

    filler->convert(elec,anEvent,_eeeeNtuple);
    
    _eeeeNtuple->column("otherP",PacPidCand2Ntuple::packFloat((float) other->p(),10));
    _eeeeNtuple->column("otherTheta",PacPidCand2Ntuple::packFloat((float) other->p3().theta(),12));
    _eeeeNtuple->column("otherPhi",PacPidCand2Ntuple::packFloat((float) other->p3().phi(),12));
    const BtaCalQual* cQual = other->getMicroAdapter()->getCalQual();
    if (cQual) {
      _eeeeNtuple->column("otherE",PacPidCand2Ntuple::packFloat((float) cQual->rawEnergy(),12));
    } else {
      _eeeeNtuple->column("otherE",(float) -1.0);
    }
    const BtaPidInfo *pidInfo = other->getMicroAdapter()->getPidInfo();
    const Consistency& theCons = pidInfo->consistency(Pdt::lookup(PdtPid::electron),PidSystem::dch);
    _eeeeNtuple->column("otherDchecons",PacPidCand2Ntuple::packFloat((float) theCons.consistency(),12));
    _eeeeNtuple->dumpData();

    delete elec;
  }
}

///////////////////////////////////////////////////////////////////////////////
void PacPidPCNtupleWriter::extractConversion(AbsEvent *anEvent)
{
  HepAList<BtaCandidate>* elList =
    Ifd< HepAList<BtaCandidate> >::get( anEvent, _convElectronList.value());
  if (elList==NULL) { return ; }
  
  HepAListIterator<BtaCandidate> iter(*elList);
  BtaCandidate *theElectron;
  
  while (theElectron = iter()) {
    
  // find mother
    HepAList<BtaCandidate> *originalList =
      Ifd< HepAList<BtaCandidate> >::get(anEvent,_convGammaList.value());
    if (originalList==0) {
      ErrMsg(warning) << "Error ! No GammaConvCalib list " << endmsg;
      return;
    }
    bool found=false;
    HepAListIterator<BtaCandidate> iter2(*originalList);
    BtaCandidate *bc(0),*mother(0);
    while ((bc=iter2()) && !found) {
      mother=bc;
      HepAListIterator<BtaCandidate> dauIter = mother->daughterIterator();
      
      BtaCandidate *dau;
      while ((dau = dauIter()) && !found) {
	if (dau->uid() == theElectron->uid()) { found = true; }
      }
    }
    if (!found) {
      ErrMsg(warning)  << "Error ! electron not found as daughter of gamma " << endmsg;
      return;
    }
    BtaAbsVertex *vtx = mother->decayVtx();
    if (vtx==0) {
      ErrMsg(warning)  << "Error !  electron from conversion has no vertex ! " << endmsg;
      return;
    }

    BtaCandidate tempCand(*theElectron);
    tempCand.setType(Pdt::lookup(PdtPid::electron, (int)tempCand.charge()));
    filler->convert(&tempCand,anEvent,_convNtuple);
    
    double R,x,y, deltaxy, mass;
    x = mother->decayVtx()->point().x();
    y = mother->decayVtx()->point().y();
    R = sqrt(x*x+y*y);
    deltaxy = ((VtxGammaConv *) mother->decayVtx())->dxy();
    mass = mother->mass();
    _convNtuple->column("vtxR",PacPidCand2Ntuple::packFloat((float) R,10));
    _convNtuple->column("vtxMass",PacPidCand2Ntuple::packFloat((float) mass,10));
    _convNtuple->column("vtxDxy",PacPidCand2Ntuple::packFloat((float) deltaxy,10));
    _convNtuple->dumpData();
  } 
} // end extractConversions



// -----------------------------------------------------------------
// Bhabahs
// -----------------------------------------------------------------
void
PacPidPCNtupleWriter::processBhabhas(const HepAList<BtaCandidate> *elList, AbsEvent *anEvent, HepTuple *ntuple,
				  const HepAList<BtaCandidate> *otherList, bool doMuons)
{
  
// Pick up some Tag Quantities
  bool doOther=false;
  if (otherList!=0) {
    if (otherList->length()==elList->length()) {
      doOther=true;
    }
  }
  AbsEventTag* tag = Ifd<AbsEventTag>::get( anEvent ) ;
  bool radBhabha(false); tag->getBool(radBhabha,"DigiFRadiativeBhabha");
  float p1Mag(0.0); tag->getFloat(p1Mag,"p1Mag");
  float p2Mag(0.0); tag->getFloat(p2Mag,"p2Mag");
  int ntrks; tag->getInt(ntrks,"nTracks");
  static const IfdStrKey keyDflt("Default");
  HepAList< EventInfo >* infoList= Ifd<HepAList<EventInfo> >::get(_currentEvent, keyDflt);
  
  
  const EventInfo* eventInfo(infoList->first());
  
  
// Find the neutral cand wth highest E
  static const IfdStrKey keyCN("CalorNeutral");
  HepAList<BtaCandidate>* neutrList =
    Ifd<HepAList<BtaCandidate> >::get(anEvent, keyCN);
  
  if (neutrList==0) {
    ErrMsg(warning) << " ERROR : Neutral List == NULL " << endmsg; 
    return;
  }
  
  
  HepAListIterator<BtaCandidate> iter(*neutrList);
  BtaCandidate *maxCand = 0 ;
  BtaCandidate *theCand = 0 ;
  double Emax=-1;
  double eTot=0;
  
  while (theCand=iter()) {
    const BtaCalQual* calQual=theCand->getMicroAdapter()->getCalQual();
    if (calQual!=NULL) {
      double ecal = calQual->rawEnergy();
      eTot+=ecal;
      if (ecal>Emax) { maxCand = theCand; Emax = ecal; }
    }
  }
  if (maxCand==0) {
  //cout << " *** WARNING ! No neutrals ! " << endl;
    return;
  }
  
// ----------------------------------
// compute total (raw) energy
// ----------------------------------
  static const IfdStrKey keyGTL("GoodTracksLoose");
  HepAList<BtaCandidate>* chgdList =  Ifd<HepAList<BtaCandidate> >::get(anEvent, keyGTL);
  assert(chgdList!=0);
  if (chgdList->length()<2) { return ; }
  eTot += ((*chgdList)[0])->p() + ((*chgdList)[1])->p();
// ----------------------------------
// loop over list
// ----------------------------------
  int ctr=0;
  BtaCandidate *otherCand=NULL;
  HepAListIterator<BtaCandidate> iter2(*elList);
  while (theCand=iter2()) {
    
    if (doOther) { otherCand=(*otherList)[ctr]; }
    ctr++;
    const BtaCalQual* calQual =theCand->getMicroAdapter()->getCalQual();
    const BtaCalQual* calQualN = maxCand->getMicroAdapter()->getCalQual();
    double alpha;
  // compute angle between candidate's Cluster and Highest Neutral Cand
    if (calQual!=0) {
      HepPoint clus = calQual->centroid(); Hep3Vector clusV(clus.x(),clus.y(),clus.z());
      HepPoint clusN = calQualN->centroid(); Hep3Vector clusVN(clusN.x(),clusN.y(),clusN.z());
      alpha = clusV.angle(clusVN);
    } else {
      alpha=-1;
    }
  // compute angle between candidate's Cluster and closest neutral cand with E>0 / 0.1 / 0.2 / 0.5
    double alpha0=100.0,alpha1=100.0,alpha2=100.0,alpha3=100.0;
    HepAListIterator<BtaCandidate> iter3(*neutrList);
    BtaCandidate *theNeutral;
    if (calQual!=0) {
      HepPoint clus = calQual->centroid(); Hep3Vector clusV(clus.x(),clus.y(),clus.z());
      while (theNeutral=iter3()) {
	const BtaCalQual *neutralCal = theNeutral->getMicroAdapter()->getCalQual();
	if (neutralCal!=0) {
	  double neutralE = neutralCal->rawEnergy();
	  HepPoint p1 = neutralCal->centroid(); Hep3Vector p1v(p1.x(),p1.y(),p1.z());
	  double angle = p1v.angle(clusV);
	  if (angle<alpha0) { alpha0=angle; }
	  if ((neutralE>0.1) && (angle<alpha1)) { alpha1=angle; }
	  if ((neutralE>0.2) && (angle<alpha2)) { alpha2=angle; }
	  if ((neutralE>0.5) && (angle<alpha3)) { alpha3=angle; }
	}
      }
    } else {
      alpha0=-1; alpha1=-1; alpha2=-1; alpha3=-1;
    }
  // Compute critical values of 'other' can
    double otherE=-1,otherP=-1,otherPcms=-1;
    if (otherCand!=0) {
      if (otherCand->getMicroAdapter()->getCalQual()) {
	otherE = otherCand->getMicroAdapter()->getCalQual()->rawEnergy();
      }
      otherP =  otherCand->p();
      BtaCandidate *c2 = new BtaCandidate(*otherCand);
      c2->setType(Pdt::lookup(PdtPid::electron,(int) c2->charge()));
      BtaCandidate Y4S(eventInfo->cmFrame());
      BtaBooster theBooster(&Y4S);
      BtaCandidate* theBoostedCand = new BtaCandidate(theBooster.boostTo(*c2));
      otherPcms = theBoostedCand->p();
      delete c2;
      delete theBoostedCand;
    }
    ntuple->column("alpha",PacPidCand2Ntuple::packFloat((float) alpha,10));
    ntuple->column("alpha0",PacPidCand2Ntuple::packFloat((float) alpha0,10));
    ntuple->column("alpha1",PacPidCand2Ntuple::packFloat((float) alpha1,10));
    ntuple->column("alpha2",PacPidCand2Ntuple::packFloat((float) alpha2,10));
    ntuple->column("alpha3",PacPidCand2Ntuple::packFloat((float) alpha3,10));
    ntuple->column("ntrk",ntrks);
    ntuple->column("etot",PacPidCand2Ntuple::packFloat((float) eTot,10));
    ntuple->column("radbhabha",radBhabha);
    ntuple->column("p1mag",PacPidCand2Ntuple::packFloat((float) p1Mag,8));
    ntuple->column("p2mag",PacPidCand2Ntuple::packFloat((float) p2Mag,8));
    ntuple->column("othere",PacPidCand2Ntuple::packFloat((float) otherE,8));
    ntuple->column("otherp",PacPidCand2Ntuple::packFloat((float) otherP,8));
    ntuple->column("otherpcms",PacPidCand2Ntuple::packFloat((float) otherPcms,8));

// avtelnov, 20080327: enforce proper lundId for ntp201, ntp202, ntp204, ntp205, ntp207
// in an R24-compatible way
    BtaCandidate *tempCand = new BtaCandidate(*theCand);
    if (doMuons) 
      tempCand->setType(Pdt::lookup(PdtPid::muon, (int)tempCand->charge()));
    else
      tempCand->setType(Pdt::lookup(PdtPid::electron, (int)tempCand->charge()));
            
    filler->convert(tempCand,anEvent,ntuple);
    
    ntuple->dumpData();
    delete tempCand;
  }
} // end extractBhabhas



/////////////////////////////////////////////////////////////////////////
void PacPidPCNtupleWriter::extractTau31(AbsEvent *anEvent)
{
  
//HepAList<BtaCandidate> *pions =
//  Ifd< HepAList<BtaCandidate> >::get( anEvent, "BtaTau31SamplePions");
  
  HepAList<BtaCandidate> *pions =
    Ifd< HepAList<BtaCandidate> >::get( anEvent, _tau31piList.value());
  
  HepAList<BtaCandidate> *pr1L =
    Ifd< HepAList<BtaCandidate> >::get( anEvent, _tau31OneProngList.value());
  
  if ((pions==0) || (pr1L==0)) { return; }
  if (pions->length()!=3) {
    ErrMsg(warning) << " WARNING : != 3 Pions in Tau31 Sample " << endmsg;
    return;
  }
  if (pr1L->length()!=1) {
    ErrMsg(warning) << " WARNING : != 1 One-Prong in Tau31 Sample " << endmsg;
    return;
  }
  
  BtaCandidate *pr1 = (*pr1L)[0];
// ----------------------------------------------
// Get information about number of neutral cands
// ----------------------------------------------
  
  static const IfdStrKey keyCN("CalorNeutral");
  HepAList<BtaCandidate> *neut =
    Ifd< HepAList<BtaCandidate> >::get( anEvent, keyCN);
  
  int ngam1=0, ngam2=0, ngam3=0;
  HepAListIterator<BtaCandidate> iter(*neut);
  BtaCandidate *c;
  while (c=iter()) {
    if (c->energy()>0.2) { ngam1++; }
    if (c->energy()>0.5) { ngam2++; }
    if (c->energy()>1.0) { ngam3++; }
  }
  
// ----------------------------------------------
// pick out the "best" pi0
// ----------------------------------------------
  int ngam = neut->length();
  float mindif = 1000;
  float Pi0mass = 0;
  for (int i1=0; i1<ngam-1; i1++) {
    for (int i2=i1+1; i2<ngam; i2++) {
      BtaCandidate *g1 = (*neut)[i1];
      BtaCandidate *g2 = (*neut)[i2];
      if ((g1->energy()>0.2) && (g2->energy()>0.2)) {
	HepLorentzVector p4 = g1->p4() + g2->p4();
	float mass = p4.mag();
	if (fabs(mass-0.135)<mindif) { Pi0mass = mass; mindif = fabs(mass-0.135);  }
      }
    }
  }
  
  
// ----------------------------------------------
// Get Information about the 1 prong
// -----------------------------------------------
  float pr1eraw,pr1ecal,pr1lmom,pr1zmom42,pr1dedxdch,pr1p,pr1theta,pr1phi,pr1ifrlayhits,pr1ifrmil,pr1nsvt;
  pr1p = pr1->p(); pr1theta = pr1->p3().theta(); pr1phi = pr1->p3().phi();
  const BtaCalQual* calQual = pr1->getMicroAdapter()->getCalQual();
  if (calQual==0) {
    pr1eraw=-1; pr1ecal=-1; pr1lmom=-1; pr1zmom42=-1; pr1dedxdch=-1;
  } else {
    pr1eraw = calQual->rawEnergy();
    pr1ecal  = calQual->ecalEnergy();
    pr1lmom = calQual->lateralMoment();
    pr1zmom42 = calQual->absZernike42();
  }
  const BtaPidQual* pQual = pr1->getMicroAdapter()->getPidQual();
  if (pQual==0) {
    pr1dedxdch=-1;
  } else {
    pr1dedxdch = pQual->dEdXDch();
  }
  
  const BtaIfrQual* iQual = pr1->getMicroAdapter()->getIfrQual();
  if (iQual==0) {
    pr1ifrlayhits=-1;
    pr1ifrmil=-1;
  } else {
    pr1ifrlayhits=iQual->IfrLayHits();
    pr1ifrmil=iQual->measuredInteractionLengths();
  }
  
  const BtaTrkQual *pr1trk = pr1->getMicroAdapter()->getTrkQual();
  if (pr1trk==0) {
    pr1nsvt=0;
  } else {
    pr1nsvt=pr1trk->nSvtHits();
  }
  
// Examine vertices
  const VtxLeastChiVertexer theVertexer;
  BtaOpMakeTree theCombiner(theVertexer);
  BtaCandidate *pi1 = (*pions)[0];   BtaCandidate *pi2 = (*pions)[1];   BtaCandidate *pi3 = (*pions)[2];
  BtaCandidate *comb12 = theCombiner.create(*pi1,*pi2);
  BtaCandidate *comb13 = theCombiner.create(*pi1,*pi3);
  BtaCandidate *comb23 = theCombiner.create(*pi2,*pi3);
  float prob12 = probab(comb12->decayVtx()->nDof(),comb12->decayVtx()->chiSquared());
  float prob13 = probab(comb13->decayVtx()->nDof(),comb13->decayVtx()->chiSquared());
  float prob23 = probab(comb23->decayVtx()->nDof(),comb23->decayVtx()->chiSquared());
  BtaCandidate *pair; BtaCandidate *other;
  if      ((prob12 >= prob13) && (prob12 >=prob23)) { pair = comb12; other = pi3; }
  else if ((prob13 >= prob12) && (prob13 >=prob23)) { pair = comb13; other = pi2; }
  else if ((prob23 >= prob12) && (prob23 >=prob13)) { pair = comb23; other = pi1; }
  else { 
    ErrMsg(warning) << "WARINIG PacPidPCNtupleWriter::extractTau31 : Something is very wrong .... " << endmsg;  
    return;  
  }
  float pairprob =  probab(pair->decayVtx()->nDof(),pair->decayVtx()->chiSquared());
  float pairchi2 =  pair->decayVtx()->chiSquared();
  TrkPocaXY pocaxy(other->trkAbsFit()->traj(),0,pair->decayVtx()->point());
  HepPoint pairP = pair->decayVtx()->point();
  HepPoint ca = other->trkAbsFit()->position(pocaxy.fltl1());
  Hep3Vector caV(ca.x(),ca.y(),ca.z());
  Hep3Vector pairV(pairP.x(),pairP.y(),pairP.z());
  Hep3Vector dif = caV - pairV;
  float dxy = sqrt(dif.x()*dif.x()+dif.y()*dif.y());
  BtaCandidate *three = theCombiner.create(*pi1,*pi2,*pi3);
  float vtxprob =   probab(three->decayVtx()->nDof(),three->decayVtx()->chiSquared());
  float vtxchi2 =   three->decayVtx()->chiSquared();
  BtaCandidate *four = theCombiner.create(*pi1,*pi2,*pi3,*pr1);
  float vtx4prob =   probab(four->decayVtx()->nDof(),four->decayVtx()->chiSquared());
  float vtx4chi2 =   four->decayVtx()->chiSquared();
// distance pr1 <-> 3-vertex in XY plane
  TrkPocaXY pocaxy2(pr1->trkAbsFit()->traj(),0,three->decayVtx()->point());
  HepPoint pr1Poca = pr1->trkAbsFit()->position(pocaxy2.fltl1());
  Hep3Vector pr1PocaV(pr1Poca.x(),pr1Poca.y(),pr1Poca.z());
  Hep3Vector threeV(three->decayVtx()->point().x(),three->decayVtx()->point().y(),three->decayVtx()->point().z());
  Hep3Vector pr1PocaDif = pr1PocaV - threeV;
  float pr1dxy = sqrt(pr1PocaDif.x()*pr1PocaDif.x() + pr1PocaDif.y()*pr1PocaDif.y());
// distance pr1 <-> 3-vertex
  TrkPoca poca(pr1->trkAbsFit()->traj(),0,three->decayVtx()->point());
  float pr1dist = poca.doca();
  
// ----------------------------------------------
// How many pi's have SVT - info ?
// ----------------------------------------------
  
  int trksWithSvt =0;
  if (pi1->getMicroAdapter()->getTrkQual()->nSvtHits() > 0 ) { trksWithSvt++; }
  if (pi2->getMicroAdapter()->getTrkQual()->nSvtHits() > 0 ) { trksWithSvt++; }
  if (pi3->getMicroAdapter()->getTrkQual()->nSvtHits() > 0 ) { trksWithSvt++; }
  
// ---------------------------------------------------------------------
// Count Number of pi's passing loose electron selection and check 1pr
// ---------------------------------------------------------------------
  // NA
  /*
  _selector->setParmValue("criteria","loose");
  
  int nLooseElectrons=0;
  int pr1IsLooseElec=0;
  if (_selector->accept(pi1)) { nLooseElectrons++; }
  if (_selector->accept(pi2)) { nLooseElectrons++; }
  if (_selector->accept(pi3)) { nLooseElectrons++; }
  if (_selector->accept(pr1)) { pr1IsLooseElec=-1; }
  
  int nyloose=0, nytight=0;
  if (isYloose(pi1)) nyloose++;   if (isYloose(pi2)) nyloose++;
  if (isYloose(pi3)) nyloose++;   if (isYloose(pr1)) nyloose++;
  if (isYtight(pi1)) nytight++;   if (isYtight(pi2)) nytight++;
  if (isYtight(pi3)) nytight++;   if (isYtight(pr1)) nytight++;
  
  
// -------------------------------------------------------
// Count Number of pi's passing tight electron selection
// -------------------------------------------------------
  _selector->setParmValue("criteria","tight");
  int nTightElectrons=0;
  int pr1IsTightElec=0;
  if (_selector->accept(pi1)) { nTightElectrons++; }
  if (_selector->accept(pi2)) { nTightElectrons++; }
  if (_selector->accept(pi3)) { nTightElectrons++; }
  if (_selector->accept(pr1)) { pr1IsTightElec=-1; }
  */
  
// -------------------------------------------------------
// dump everything into ntuple
// -------------------------------------------------------
  for (int i=0; i<3; i++) {
    BtaCandidate *thePi = (*pions)[i];
    filler->convert(thePi,anEvent,_tau31Ntuple);
    
    _tau31Ntuple->column("ngam1",ngam1);
    _tau31Ntuple->column("ngam2",ngam2);
    _tau31Ntuple->column("ngam3",ngam3);
    _tau31Ntuple->column("pi0mass",PacPidCand2Ntuple::packFloat((float) Pi0mass,12)); // 7% resolution, packing adds 1e-4 in quadrature
    _tau31Ntuple->column("pr1p",PacPidCand2Ntuple::packFloat((float) pr1p,12));
    _tau31Ntuple->column("pr1theta",PacPidCand2Ntuple::packFloat((float) pr1theta,12));
    _tau31Ntuple->column("pr1phi",PacPidCand2Ntuple::packFloat((float) pr1phi,12));
    _tau31Ntuple->column("pr1dxy",PacPidCand2Ntuple::packFloat((float) pr1dxy,12));
    _tau31Ntuple->column("pr1dist",PacPidCand2Ntuple::packFloat((float) pr1dist,12));
    _tau31Ntuple->column("pr1eraw",PacPidCand2Ntuple::packFloat((float) pr1eraw,12)); // 20% resolution
    _tau31Ntuple->column("pr1ecal",PacPidCand2Ntuple::packFloat((float) pr1ecal,12));
    _tau31Ntuple->column("pr1lmom",PacPidCand2Ntuple::packFloat((float) pr1lmom,12));
    _tau31Ntuple->column("pr1zmom42",PacPidCand2Ntuple::packFloat((float) pr1zmom42,12));
    _tau31Ntuple->column("pr1dedxdch",PacPidCand2Ntuple::packFloat((float) pr1dedxdch,12));
    _tau31Ntuple->column("pr1nsvt",pr1nsvt);
    _tau31Ntuple->column("pr1ifrlayhits",pr1ifrlayhits);
    _tau31Ntuple->column("pr1ifrmil",PacPidCand2Ntuple::packFloat((float) pr1ifrmil,12));
    _tau31Ntuple->column("pairprob",PacPidCand2Ntuple::packFloat((float) pairprob,12));
    _tau31Ntuple->column("pairchi2",PacPidCand2Ntuple::packFloat((float) pairchi2,12));
    _tau31Ntuple->column("dxy",PacPidCand2Ntuple::packFloat((float) dxy,12));
    _tau31Ntuple->column("vtxprob",PacPidCand2Ntuple::packFloat((float) vtxprob,12));
    _tau31Ntuple->column("vtxchi2",PacPidCand2Ntuple::packFloat((float) vtxchi2,12));
    _tau31Ntuple->column("vtx4prob",PacPidCand2Ntuple::packFloat((float) vtx4prob,12));
    _tau31Ntuple->column("vtx4chi2",PacPidCand2Ntuple::packFloat((float) vtx4chi2,12));
    _tau31Ntuple->column("ntrksvt",trksWithSvt);

    // NA
    //_tau31Ntuple->column("nloose",nLooseElectrons);
    //_tau31Ntuple->column("ntight",nTightElectrons);
    //_tau31Ntuple->column("nyloose",nyloose);
    //_tau31Ntuple->column("nytight",nytight);
    //_tau31Ntuple->column("pr1loose",pr1IsLooseElec);
    //_tau31Ntuple->column("pr1tight",pr1IsTightElec);

    _tau31Ntuple->dumpData();
  }
  delete comb12;
  delete comb13;
  delete comb23;
  delete three;
  delete four;
  
   
} // end tau31 pions



/////////////////////////////////////////////////////////////
bool PacPidPCNtupleWriter::isYloose(const BtaCandidate *cand)
{
  bool answer = false;
  const BtaPidQual *PidQual;
  const BtaCalQual *CalQual;
  double momentum = cand->p();
  double ctheta = cand->p3().cosTheta();
  CalQual = cand->getMicroAdapter()->getCalQual();
  PidQual = cand->getMicroAdapter()->getPidQual();
//  Consider also calorimeter acceptance   19.2 -> 138 degrees
  if (CalQual  && PidQual && ctheta > -0.74314483 && ctheta < 0.94437637) {
    float Ecal = CalQual->ecalEnergy();
    int Ncrys = CalQual->nCrystals();
    float dedxdch = PidQual->dEdXDch();
    if ((Ecal/momentum>0.65) && (Ecal/momentum) <5.
	&& Ncrys>3 && 500.<dedxdch && dedxdch<1000.){
      answer=true;
    }
  }
  return answer;
}

////////////////////////////////////////
bool PacPidPCNtupleWriter::isYtight(const BtaCandidate *cand)
{
  bool answer = false;
  const BtaPidQual *PidQual;
  const BtaCalQual *CalQual;
  double momentum = cand->p();
  double ctheta = cand->p3().cosTheta();
  CalQual = cand->getMicroAdapter()->getCalQual();
  PidQual = cand->getMicroAdapter()->getPidQual();
//  Consider also calorimeter acceptance   19.2 -> 138 degrees
  if (CalQual  && PidQual && ctheta > -0.74314483 && ctheta < 0.94437637) {
    float Ecal = CalQual->ecalEnergy();
    float lat = CalQual->lateralMoment();
    int Ncrys = CalQual->nCrystals();
    float dedxdch = PidQual->dEdXDch();
    if ((Ecal/momentum>0.65) && (Ecal/momentum) <5.
	&& Ncrys>3 && 500.<dedxdch && dedxdch<1000.){
      if ((Ecal/momentum)>0.75 && (Ecal/momentum) <1.3 && lat<0.6) answer=true;
    }
  }
  return answer;
}

////////////////////////////////////////////////////////
void PacPidPCNtupleWriter::doPhi(AbsEvent* anEvent)
{
  BtaCandidate *theKaon,*theCand;

  static const IfdStrKey keyKMicroLoose("KMicroLoose");
  HepAList<BtaCandidate> *kaonList = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyKMicroLoose); // KaonSMS loose per default
  if (kaonList == NULL) ErrMsg(fatal) << "No list KMicroLoose in the event" << endmsg;

  static const IfdStrKey keyKMicroTight("KMicroTight");
  HepAList<BtaCandidate> *kaonListT = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyKMicroTight);
  if (kaonListT == NULL) ErrMsg(fatal) << "No list KMicroTight in the event" << endmsg;

  static const IfdStrKey keyKMicroVeryTight("KMicroVeryTight");
  HepAList<BtaCandidate> *kaonListVT = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyKMicroVeryTight);
  if (kaonListVT == NULL) ErrMsg(fatal) << "No list KMicroVeryTight in the event" << endmsg;

  static const IfdStrKey keyGoodTracksLoose("GoodTracksLoose");
  HepAList<BtaCandidate> *chargedList = Ifd<HepAList<BtaCandidate> >::get(anEvent,keyGoodTracksLoose); // Good tracks loose per default
  if (chargedList == NULL) ErrMsg(fatal) << "No list GoodTracksLoose in the event" << endmsg;

  HepAListIterator<BtaCandidate> kIter(*kaonList);
  HepAListIterator<BtaCandidate> cIter(*chargedList);
  VtxLeastChiVertexer lcV;
  BtaOpMakeTree phiMaker(lcV);
  while (theKaon=kIter()) {
    while (theCand=cIter()) {
      if (theKaon->charge()==-theCand->charge()) {
	BtaCandidate copy(*theCand);
	copy.setType(Pdt::lookup(PdtPid::kaon,(int) theCand->charge()));
	BtaCandidate *aPhi = phiMaker.create(*theKaon,copy);
	double mass = aPhi->mass();
	double prob = probab((aPhi->decayVtx()->nDof()), (aPhi->decayVtx()->chiSquared()) );
	if ((mass>0.99) && (mass<1.05) && (prob>0)) {
	  filler->convert(theCand,anEvent,_kPhiTuple);
	  
	  addVtxQuantities(aPhi,anEvent,_kPhiTuple);
	  bool otherSmsTight=false,otherSmsVeryTight=false;
	  int i;
	  for (i=0; i<kaonListT->length(); i++) {
	    if ((*kaonListT)[i]->uid() == theKaon->uid()) {  otherSmsTight=true; }
	  }
	  for (i=0; i<kaonListVT->length(); i++) {
	    if ((*kaonListVT)[i]->uid() == theKaon->uid()) {  otherSmsVeryTight=true; }
	  }
	  _kPhiTuple->column("othersmstight",otherSmsTight);
	  _kPhiTuple->column("othersmsverytight",otherSmsVeryTight);
	  _kPhiTuple->column("otherp",PacPidCand2Ntuple::packFloat((float) theKaon->p(),12));
	  _kPhiTuple->column("othertheta",PacPidCand2Ntuple::packFloat((float) theKaon->p3().theta(),14));
	  _kPhiTuple->column("otherphi",PacPidCand2Ntuple::packFloatPhi((float) theKaon->p3().phi(),14));
	  _kPhiTuple->column("mPhi",PacPidCand2Ntuple::packFloat((float) mass,14));
	  const BtaPidQual *pQual = theKaon->getMicroAdapter()->getPidQual();
	  if (pQual) {
	    _kPhiTuple->column("othernphote",pQual->ringNExPhot(Pdt::lookup(PdtPid::electron)));
	  } else {
	    _kPhiTuple->column("othernphote",(int) 99);
	  }
	  _kPhiTuple->dumpData();
	}
	delete aPhi;
      }
    }
  }
} // end the phi sample


///////////////////////////////////////////////////////////////////////////////////
void PacPidPCNtupleWriter::extractVcs(AbsEvent *anEvent)
{
  static unsigned int callCounter = 0;
  callCounter++;
  static char msg[] = " in the first processed event! This is the only warning you get.";
  HepAList<BtaCandidate> *vcsList = Ifd<HepAList<BtaCandidate> >::get(anEvent, _vcsListName.value());
  if (vcsList == NULL && callCounter==1) ErrMsg(error) << "No list " << _vcsListName.value() << msg << endmsg;
  HepAList<BtaCandidate> *vcsGammaList = Ifd<HepAList<BtaCandidate> >::get(anEvent, _vcsGammaListName.value());
  if (vcsGammaList == NULL && callCounter==1) ErrMsg(error) << "No list " << _vcsGammaListName.value() << msg << endmsg;
//getTmpAList(anEvent,vcsList,_vcsListName.value());
//getTmpAList(anEvent,vcsGammaList,_vcsGammaListName.value());
  if (vcsList->length()==0) { return ; }
  if (vcsGammaList->length()==0) { return ; }
  BtaCandidate *elec = (*vcsList)[0];
  BtaCandidate *gam = (*vcsGammaList)[0];
  filler->convert(elec,anEvent,_vcsTuple);

  HepLorentzVector p4Gam = gam->p4();  p4Gam.setVectM(p4Gam.vect(),0.0);
  HepLorentzVector p4Elec = elec->p4();

  // 3 - coordinates of gamma
  Hep3Vector pGam = p4Gam.vect();
  _vcsTuple->column("egam",PacPidCand2Ntuple::packFloat((float) pGam.mag(),12));
  _vcsTuple->column("thegam",PacPidCand2Ntuple::packFloat((float) pGam.theta(),14));
  _vcsTuple->column("phigam",PacPidCand2Ntuple::packFloatPhi((float) pGam.phi(),14));
  
// CMS - coordinates of gamma
  static const IfdStrKey keyDflt("Default");
  HepAList< EventInfo >* infoList= Ifd<HepAList<EventInfo> >::get(anEvent, keyDflt);
  assert(infoList != 0);
  EventInfo* eventInfo = infoList->first();
  HepLorentzVector cmframe = eventInfo->cmFrame();
  p4Gam.boost(-(cmframe.boostVector()) );
  p4Elec.boost(-(cmframe.boostVector()) );
  pGam = p4Gam.vect();
  _vcsTuple->column("egamcms",PacPidCand2Ntuple::packFloat((float) pGam.mag(),12));
  _vcsTuple->column("thegamcms",PacPidCand2Ntuple::packFloat((float) pGam.theta(),14));
  _vcsTuple->column("phigamcms",PacPidCand2Ntuple::packFloatPhi((float) pGam.phi(),14));
  
// Accolinearity
  double PI = Constants::pi;
  double DeltaThetaCms = PI-fabs (p4Gam.vect().theta()+p4Elec.vect().theta());
  double DeltaPhiCms = fabs(p4Gam.vect().phi()-p4Elec.vect().phi())-PI;
  _vcsTuple->column("deltathetacms",PacPidCand2Ntuple::packFloat((float) DeltaThetaCms,12));
  _vcsTuple->column("deltaphicms",PacPidCand2Ntuple::packFloat((float) DeltaPhiCms,12));
  _vcsTuple->dumpData();
} // end Vcs electron sample



//////////////////////////////////////
void PacPidPCNtupleWriter::extractJpsiKmm(AbsEvent *anEvent)
{
  AbsEventTag* tag = Ifd<AbsEventTag>::get( anEvent );
  
  int nTracks(0);
  float r2(0);
  bool bgfmultihadron;
  bool status;
  
  bool gotjpsikstlplus = false;
  bool gotjpsikstlminus = false;
  
  if ( tag != 0 )
    {
      status  = tag->getBool( bgfmultihadron, "BGFMultiHadron" );
      status &= tag->getInt( nTracks, "nTracks" );
      status &= tag->getFloat( r2, "R2");
      
      if((!status) || (nTracks < 4) || (r2 > 0.7) || (bgfmultihadron == false))
	return;
    }
  
  HepAList< BtaCandidate >* muList;
  HepAList< BtaCandidate >* muTightList;
  HepAList< BtaCandidate >* lepNoIdList;
  HepAList< BtaCandidate >* kList;
  HepAList< BtaCandidate >* ksList;
  HepAList< BtaCandidate >* kstList;
 
  static const IfdStrKey keyGTL("GoodTracksLoose");
  getTmpAList(anEvent, lepNoIdList, keyGTL );

  static const IfdStrKey keyMu("muNNVeryLoose");
  getTmpAList(anEvent, muList, keyMu );

  static const IfdStrKey keyMuTight("muNNTight");
  getTmpAList(anEvent, muTightList, keyMuTight );
  
  static const IfdStrKey keyKLHLoose("KLHLoose");
  getTmpAList(anEvent, kList, keyKLHLoose );
  static const IfdStrKey keyKsDefault("KsDefault");
  getTmpAList(anEvent, ksList, keyKsDefault );
  static const IfdStrKey keyKstarKPiDefaultPID("KstarKPiDefaultPID");
  getTmpAList(anEvent, kstList, keyKstarKPiDefaultPID );
  
  HepAListIterator<BtaCandidate> itMuTag( *muList );
  HepAListIterator<BtaCandidate> itMuTightTag( *muTightList );
  HepAListIterator<BtaCandidate> itTrk( *lepNoIdList );
  HepAListIterator<BtaCandidate> itK( *kList );
  HepAListIterator<BtaCandidate> itKs( *ksList );
  HepAListIterator<BtaCandidate> itKst( *kstList );
  
  BtaCandidate *theMuTag(0), *theTrk(0), *theK(0), *theKs(0), *theKst(0),
    *theMuTightTag(0);
  
  
  BtaOpFastVtx theCombiner;
  BtaOpFastVtx theBCombiner;
  BtaOpFastVtx theBCombinerb;
  BtaOpFastVtx theBCombinerc;
  
  while(theMuTag = itMuTag()) {
    itTrk.rewind();
    while(theTrk = itTrk()) {
      if(theTrk->charge() + theMuTag->charge() != 0)
	continue;
      if(theTrk->overlaps(*theMuTag))
	continue;
      
      BtaCandidate theLocalTrk( *theTrk );
      
      if(theLocalTrk.charge() == 1)
	theLocalTrk.setType(Pdt::lookup("mu+"));
      else if(theLocalTrk.charge() == -1)
	theLocalTrk.setType(Pdt::lookup("mu-"));
      
      BtaCandidate *theCC =
	theCombiner.create(*theMuTag,theLocalTrk);
      
    // Check for a valid vertex
      if(theCC->decayVtx() == 0)
	{
	  delete theCC;
	  continue;
	}
      
      if(theCC->mass() < 2.9 || theCC->mass() > 3.2)
	{
	  delete theCC;
	  continue;
	}
      
      
      int otherismunnt = 0;
      itMuTightTag.rewind();
      while(theMuTightTag = itMuTightTag())
	{
	// "overlaps" is expensive, so check charges first
	  if(theMuTightTag->charge()==theMuTag->charge() && theMuTightTag->overlaps(*theMuTag))
	    otherismunnt = 1;
	}
      
    // Jpsi K+
      itK.rewind();
      while(theK = itK()) {
	if( (theK->charge()==theMuTag->charge() && theK->overlaps(*theMuTag)) 
	    || (theK->charge()==theTrk->charge() && theK->overlaps(*theTrk)) )
	  continue;
	
	BtaCandidate *theB =
	  theBCombiner.create(theLocalTrk,*theMuTag,*theK);
	
      // Check for a valid vertex
	if(theB->decayVtx() == 0)
	  {
	    delete theB;
	    continue;
	  }
	
	// NA
	//JPsiHeli *theHelicity = new JPsiHeli(theMuTag, &theLocalTrk, theK);
	//float Helicity = theHelicity->helicity();
	
	HepLorentzVector BP4 = theB->p4();
	BtaBVariables bVars(BP4, true);
	float m_ES = bVars.m_ES();
	float deltaE = bVars.deltaE();
	
	if(m_ES < 5.2 || deltaE < -0.15 || deltaE > 0.15 || m_ES > 5.3 )
	  {
	    delete theB;
	    continue;
	  }
	
	addVtxQuantities(theB,anEvent,_jpsikmmTuple);
	
	_jpsikmmTuple->column("otherp",PacPidCand2Ntuple::packFloat(theMuTag->p(),8));
	_jpsikmmTuple->column("othertheta",PacPidCand2Ntuple::packFloat(theMuTag->p4().theta(),10));
	_jpsikmmTuple->column("otherphi",PacPidCand2Ntuple::packFloatPhi(theMuTag->p4().phi(),10));
	_jpsikmmTuple->column("othercharge",theMuTag->charge());
	_jpsikmmTuple->column("otherismunnt",PacPidCand2Ntuple::packFloat(otherismunnt,8));
	_jpsikmmTuple->column("kmass",PacPidCand2Ntuple::packFloat(theK->mass(),8));
	_jpsikmmTuple->column("kmode",1);
	_jpsikmmTuple->column("kp",PacPidCand2Ntuple::packFloat(theK->p(),8));
	_jpsikmmTuple->column("mPsi",PacPidCand2Ntuple::packFloat(theCC->mass(),8));
	_jpsikmmTuple->column("mES",PacPidCand2Ntuple::packFloat(m_ES,6));
	_jpsikmmTuple->column("deltaE",PacPidCand2Ntuple::packFloat(deltaE,8));
      // _jpsikmmTuple->column("R2",r2);
	_jpsikmmTuple->column("nTrk",nTracks);
	
	// NA
	//_jpsikmmTuple->column("helicity",PacPidCand2Ntuple::packFloat((float) Helicity,10));
	
	filler->convert(&theLocalTrk,anEvent,_jpsikmmTuple);
	filler->convertJpsiKMC(&theLocalTrk,theMuTag,theK,anEvent,_jpsikmmTuple);
	_jpsikmmTuple->dumpData();
	
	delete theB;
	
	// NA
	//delete theHelicity;
      }
	
      // Jpsi Ks
      while(theKs = itKs()) {
	if(theKs->overlaps(*theMuTag) || theKs->overlaps(*theTrk))
	  continue;
	
	BtaCandidate *theBb =
	  theBCombinerb.create(theLocalTrk,*theMuTag,*theKs);
	
      // Check for a valid vertex
	if(theBb->decayVtx() == 0)
	  {
	    delete theBb;
	    continue;
	  }
	
	// NA
	//JPsiHeli *theHelicity = new JPsiHeli(theMuTag, &theLocalTrk, theKs);
	//float Helicity = theHelicity->helicity();
	
	HepLorentzVector BP4 = theBb->p4();
	BtaBVariables bVarsb(BP4, true);
	float m_ES = bVarsb.m_ES();
	float deltaE = bVarsb.deltaE();
	
	if(m_ES < 5.2 || deltaE < -0.15 || deltaE > 0.15 || m_ES > 5.3)
	  {
	    delete theBb;
	    continue;
	  }
	
	addVtxQuantities(theBb,anEvent,_jpsikmmTuple);
	
	_jpsikmmTuple->column("otherp",PacPidCand2Ntuple::packFloat(theMuTag->p(),8));
	_jpsikmmTuple->column("othertheta",PacPidCand2Ntuple::packFloat(theMuTag->p4().theta(),10));
	_jpsikmmTuple->column("otherphi",PacPidCand2Ntuple::packFloatPhi(theMuTag->p4().phi(),10));
	_jpsikmmTuple->column("othercharge",theMuTag->charge());
	_jpsikmmTuple->column("otherismunnt",PacPidCand2Ntuple::packFloat(otherismunnt,8));
	_jpsikmmTuple->column("kmass",PacPidCand2Ntuple::packFloat(theKs->mass(),8));
	_jpsikmmTuple->column("kmode",2);
	_jpsikmmTuple->column("kp",PacPidCand2Ntuple::packFloat(theKs->p(),8));
	_jpsikmmTuple->column("mPsi",PacPidCand2Ntuple::packFloat(theCC->mass(),8));
	_jpsikmmTuple->column("mES",PacPidCand2Ntuple::packFloat(m_ES,6));
	_jpsikmmTuple->column("deltaE",PacPidCand2Ntuple::packFloat(deltaE,8));
      // _jpsikmmTuple->column("R2",r2);
	_jpsikmmTuple->column("nTrk",nTracks);

	// NA
	//_jpsikmmTuple->column("helicity",PacPidCand2Ntuple::packFloat((float) Helicity,10));
	
	filler->convert(&theLocalTrk,anEvent,_jpsikmmTuple);
	filler->convertJpsiKMC(&theLocalTrk,theMuTag,theKs,anEvent,_jpsikmmTuple);
	_jpsikmmTuple->dumpData();
	
	delete theBb;
	
	// NA
	//delete theHelicity;
      }
	
    // Multiple cands. are an issue for the K* mode. Only keep
    // the first K* candidate in the event.
      if((theTrk->charge() == 1) && (gotjpsikstlplus == true))
	continue;
      if((theTrk->charge() == -1) && (gotjpsikstlminus == true))
	continue;
      
    // Jpsi K*
      itKst.rewind();
      while(theKst = itKst()) {
	if(theKst->overlaps(*theMuTag) ||
	   theKst->overlaps(theLocalTrk))
	  continue;
	
	BtaCandidate *theBc =
	  theBCombinerc.create(theLocalTrk,*theMuTag,*theKst);
	
      // Check for a valid vertex
	if(theBc->decayVtx() == 0)
	  {
	    delete theBc;
	    continue;
	  }

	// NA
	//JPsiHeli *theHelicity = new JPsiHeli(theMuTag, &theLocalTrk, theKst);
	//float Helicity = theHelicity->helicity();
	
	HepLorentzVector BP4 = theBc->p4();
	BtaBVariables bVarsc(BP4, true);
	float m_ES = bVarsc.m_ES();
	float deltaE = bVarsc.deltaE();
	
	if(m_ES < 5.2 || deltaE < -0.15 || deltaE > 0.15 || m_ES > 5.3)
	  {
	    delete theBc;
	    continue;
	  }
	
	if(theTrk->charge() == 1)
	  gotjpsikstlplus = true;
	if(theTrk->charge() == -1)
	  gotjpsikstlminus = true;
	
	addVtxQuantities(theBc,anEvent,_jpsikmmTuple);
	
	_jpsikmmTuple->column("otherp",PacPidCand2Ntuple::packFloat(theMuTag->p(),8));
	_jpsikmmTuple->column("othertheta",PacPidCand2Ntuple::packFloat(theMuTag->p4().theta(),10));
	_jpsikmmTuple->column("otherphi",PacPidCand2Ntuple::packFloatPhi(theMuTag->p4().phi(),10));
	_jpsikmmTuple->column("othercharge",theMuTag->charge());
	_jpsikmmTuple->column("otherismunnt",PacPidCand2Ntuple::packFloat(otherismunnt,8));
	_jpsikmmTuple->column("kmass",PacPidCand2Ntuple::packFloat(theKst->mass(),8));
	_jpsikmmTuple->column("kmode",3);
	_jpsikmmTuple->column("kp",PacPidCand2Ntuple::packFloat(theKst->p(),8));
	_jpsikmmTuple->column("mPsi",PacPidCand2Ntuple::packFloat(theCC->mass(),8));
	_jpsikmmTuple->column("mES",PacPidCand2Ntuple::packFloat(m_ES,6));
	_jpsikmmTuple->column("deltaE",PacPidCand2Ntuple::packFloat(deltaE,8));
      // _jpsikmmTuple->column("R2",r2);
	_jpsikmmTuple->column("nTrk",nTracks);

	// NA
	//_jpsikmmTuple->column("helicity",PacPidCand2Ntuple::packFloat((float) Helicity,10));

	filler->convert(&theLocalTrk,anEvent,_jpsikmmTuple);
	filler->convertJpsiKMC(&theLocalTrk,theMuTag,theKst,anEvent,_jpsikmmTuple);
	_jpsikmmTuple->dumpData();
	
	delete theBc;
	
	// NA
	//delete theHelicity;
      }
      delete theCC;
      
    }
  }
}  // end J/psi mu sample





//////////////////////////////////////////
void PacPidPCNtupleWriter::extractJpsiKee(AbsEvent *anEvent)
{
  AbsEventTag* tag = Ifd<AbsEventTag>::get( anEvent );
  
  int nTracks(0);
  float r2(0);
  bool bgfmultihadron;
  bool status;
  
  bool gotjpsikstlplus = false;
  bool gotjpsikstlminus = false;
  
  if ( tag != 0 )
    {
      status  = tag->getBool( bgfmultihadron, "BGFMultiHadron" );
      status &= tag->getInt( nTracks, "nTracks" );
      status &= tag->getFloat( r2, "R2");
      
      if((!status) || (nTracks < 4) || (r2 > 0.7) || (bgfmultihadron == false))
	return;
    }
  
  
  HepAList< BtaCandidate >* eList;
  HepAList< BtaCandidate >* eTightList;
  HepAList< BtaCandidate >* lepNoIdBremList;
  HepAList< BtaCandidate >* kList;
  HepAList< BtaCandidate >* ksList;
  HepAList< BtaCandidate >* kstList;
  
  static const IfdStrKey keyGTL("GoodTracksLoose");
  static const IfdStrKey key_eBremRecoELNC("eBremRecoELNC");
  static const IfdStrKey keyPidLHElectrons("PidLHElectrons");
  static const IfdStrKey keyKLHLoose("KLHLoose");
  static const IfdStrKey keyKsDefault("KsDefault");
  static const IfdStrKey keyKstarKPiDefaultPID("KstarKPiDefaultPID");
  
  getTmpAList(anEvent, lepNoIdBremList, keyGTL );
  getTmpAList(anEvent, eList, key_eBremRecoELNC );
  getTmpAList(anEvent, eTightList, keyPidLHElectrons );
  
  getTmpAList(anEvent, kList, keyKLHLoose );
  getTmpAList(anEvent, ksList, keyKsDefault );
  getTmpAList(anEvent, kstList, keyKstarKPiDefaultPID );
  
  HepAListIterator<BtaCandidate> itElTag( *eList );
  HepAListIterator<BtaCandidate> itElTightTag( *eTightList );
  HepAListIterator<BtaCandidate> itTrk2( *lepNoIdBremList );
  HepAListIterator<BtaCandidate> itK2( *kList );
  HepAListIterator<BtaCandidate> itKs2( *ksList );
  HepAListIterator<BtaCandidate> itKst2( *kstList );
  
  BtaCandidate *theElTag(0), *theTrk2(0), *theK2(0), *theKs2(0), *theKst2(0),
    *theElTightTag(0);
  
  HepLorentzVector BP4;
  
  BtaOpFastVtx theCombiner2;
  BtaOpFastVtx theBCombiner2;
  BtaOpFastVtx theBCombiner2b;
  BtaOpFastVtx theBCombiner2c;
  
  while(theElTag = itElTag()) {
    itTrk2.rewind();
    while(theTrk2 = itTrk2())
      {
	if(theTrk2->charge() + theElTag->charge() != 0)
	  continue;
	if(theTrk2->overlaps(*theElTag))
	  continue;
	
	BtaCandidate theLocalTrk2( *theTrk2 );
	
	if(theLocalTrk2.charge() == 1)
	  theLocalTrk2.setType(Pdt::lookup("e+"));
	else if(theLocalTrk2.charge() == -1)
	  theLocalTrk2.setType(Pdt::lookup("e-"));
	
	BtaCandidate *theCC2 =
	  theCombiner2.create(*theElTag,theLocalTrk2);
	
	
      // Check for a valid vertex
	if(theCC2->decayVtx() == 0)
	  {
	    delete theCC2;
	    continue;
	  }
	
	if(theCC2->mass() < 2.9 || theCC2->mass() > 3.2)
	  {
	    delete theCC2;
	    continue;
	  }
	
	
	int otheriselh = 0;
	itElTightTag.rewind();
	while(theElTightTag = itElTightTag())
	  {
	  //        if(theTightElTag->uid() == theElTag->uid())
	    if(theElTightTag->charge()==theElTag->charge() && theElTightTag->overlaps(*theElTag))
	      otheriselh = 1;
	  }
	
      // Jpsi K+
	itK2.rewind();
	while(theK2 = itK2()) {
	  if( (theK2->charge()==theElTag->charge() && theK2->overlaps(*theElTag)) 
	      || (theK2->charge()==theTrk2->charge() && theK2->overlaps(*theTrk2)) )
	    continue;
	  
	  BtaCandidate *theB2 =
	    theBCombiner2.create(theLocalTrk2,*theElTag,*theK2);
	  
        // Check for a valid vertex
	  if(theB2->decayVtx() == 0)
	    {
	      delete theB2;
	      continue;
	    }
	  
	  
	  // NA
	  //JPsiHeli *theHelicity = new JPsiHeli(theElTag, &theLocalTrk2, theK2);
	  //float Helicity = theHelicity->helicity();
	  
	  BP4 = theB2->p4();
	  BtaBVariables bVars2(BP4, true);
	  float m_ES = bVars2.m_ES();
	  float deltaE = bVars2.deltaE();
	  
	  if(m_ES < 5.2 || deltaE < -0.15 || deltaE > 0.15 || m_ES > 5.3)
	    {
	      delete theB2;
	      continue;
	    }
	  
	  addVtxQuantities(theB2,anEvent,_jpsikeeTuple);
	  
	  _jpsikeeTuple->column("otherp",PacPidCand2Ntuple::packFloat(theElTag->p(),8));
	  _jpsikeeTuple->column("othertheta",PacPidCand2Ntuple::packFloat(theElTag->p4().theta(),10));
	  _jpsikeeTuple->column("otherphi",PacPidCand2Ntuple::packFloatPhi(theElTag->p4().phi(),10));
	  _jpsikeeTuple->column("othercharge",theElTag->charge());
	  _jpsikeeTuple->column("otheriselh", PacPidCand2Ntuple::packFloat(otheriselh,8));
	  _jpsikeeTuple->column("kmass",PacPidCand2Ntuple::packFloat(theK2->mass(),8));
	  _jpsikeeTuple->column("kmode",1);
	  _jpsikeeTuple->column("kp",PacPidCand2Ntuple::packFloat(theK2->p(),8));
	  _jpsikeeTuple->column("mPsi",PacPidCand2Ntuple::packFloat(theCC2->mass(),8));
	  _jpsikeeTuple->column("mES",PacPidCand2Ntuple::packFloat(m_ES,6));
	  _jpsikeeTuple->column("deltaE",PacPidCand2Ntuple::packFloat(deltaE,8));
        // _jpsikeeTuple->column("R2",r2);
	  _jpsikeeTuple->column("nTrk",nTracks);

	  // NA
	  //_jpsikeeTuple->column("helicity",PacPidCand2Ntuple::packFloat((float)Helicity,10));

	  filler->convert(&theLocalTrk2,anEvent,_jpsikeeTuple);
	  filler->convertJpsiKMC(&theLocalTrk2,theElTag,theK2,anEvent,_jpsikeeTuple);
	  _jpsikeeTuple->dumpData();
	  
	  delete theB2;

	  // NA
	  //delete theHelicity;
	}
	
      // Jpsi Ks
	itKs2.rewind();
	while(theKs2 = itKs2()) {
	  if(theKs2->overlaps(*theElTag) || theKs2->overlaps(theLocalTrk2))
	    continue;
	  
	  BtaCandidate *theB2b =
	    theBCombiner2b.create(theLocalTrk2,*theElTag,*theKs2);
	  
	// Check for a valid vertex
	  if(theB2b->decayVtx() == 0)
	    {
	      delete theB2b;
	      continue;
	    }
	  
	  // NA
	  //JPsiHeli *theHelicity = new JPsiHeli(theElTag, &theLocalTrk2, theKs2);
	  //float Helicity = theHelicity->helicity();
	  
	  HepLorentzVector BP4a = theB2b->p4();
	  BtaBVariables bVars2b(BP4a, true);
	  float m_ES = bVars2b.m_ES();
	  float deltaE = bVars2b.deltaE();
	  
	  if(m_ES < 5.2 || deltaE < -0.15 || deltaE > 0.15 || m_ES > 5.3)
	    {
	      delete theB2b;
	      continue;
	    }
	  
	  addVtxQuantities(theB2b,anEvent,_jpsikeeTuple);
	  
	  _jpsikeeTuple->column("otherp",PacPidCand2Ntuple::packFloat(theElTag->p(),8));
	  _jpsikeeTuple->column("othertheta",PacPidCand2Ntuple::packFloat(theElTag->p4().theta(),10));
	  _jpsikeeTuple->column("otherphi",PacPidCand2Ntuple::packFloat(theElTag->p4().phi(),10));
	  _jpsikeeTuple->column("othercharge",theElTag->charge());
	  _jpsikeeTuple->column("otheriselh", PacPidCand2Ntuple::packFloat(otheriselh,8));
	  _jpsikeeTuple->column("kmass",PacPidCand2Ntuple::packFloat(theKs2->mass(),8));
	  _jpsikeeTuple->column("kmode",2);
	  _jpsikeeTuple->column("kp",PacPidCand2Ntuple::packFloat(theKs2->p(),8));
	  _jpsikeeTuple->column("mPsi",PacPidCand2Ntuple::packFloat(theCC2->mass(),8));
	  _jpsikeeTuple->column("mES",PacPidCand2Ntuple::packFloat(m_ES,6));
	  _jpsikeeTuple->column("deltaE",PacPidCand2Ntuple::packFloat(deltaE,8));
        // _jpsikeeTuple->column("R2",r2);
	  _jpsikeeTuple->column("nTrk",nTracks);
	  
	  // NA
	  //_jpsikeeTuple->column("helicity",PacPidCand2Ntuple::packFloat((float) Helicity,10));
	  
	  filler->convert(&theLocalTrk2,anEvent,_jpsikeeTuple);
	  filler->convertJpsiKMC(&theLocalTrk2,theElTag,theKs2,anEvent,_jpsikeeTuple);
	  _jpsikeeTuple->dumpData();
	  
	  delete theB2b;

	  // NA
	  // delete theHelicity;
	}
	
      // Multiple cands. are an issue for the K* mode. Only keep
      // the first K* candidate in the event.
	if((theTrk2->charge() == 1) && (gotjpsikstlplus == true))
	  continue;
	if((theTrk2->charge() == -1) && (gotjpsikstlminus == true))
	  continue;
	
      // Jpsi K*
	itKst2.rewind();
	while(theKst2 = itKst2()) {
	  if(theKst2->overlaps(*theElTag) ||
	     theKst2->overlaps(theLocalTrk2))
	    continue;
	  
	  BtaCandidate *theB2c =
	    theBCombiner2c.create(theLocalTrk2,*theElTag,*theKst2);
	  
	// Check for a valid vertex
	  if(theB2c->decayVtx() == 0)
	    {
	      delete theB2c;
	      continue;
	    }
	  
	  // NA
	  //JPsiHeli *theHelicity = new JPsiHeli(theElTag, &theLocalTrk2, theKst2);
	  //float Helicity = theHelicity->helicity();
	  
	  HepLorentzVector BP4a = theB2c->p4();
	  BtaBVariables bVars2c(BP4a, true);
	  float m_ES = bVars2c.m_ES();
	  float deltaE = bVars2c.deltaE();
	  
	  if(m_ES < 5.2 || deltaE < -0.15 || deltaE > 0.15 || m_ES > 5.3)
	    {
	      delete theB2c;
	      continue;
	    }
	  
	  if(theTrk2->charge() == 1)
	    gotjpsikstlplus = true;
	  if(theTrk2->charge() == -1)
	    gotjpsikstlminus = true;
	  
	  addVtxQuantities(theB2c,anEvent,_jpsikeeTuple);
	  
	  _jpsikeeTuple->column("otherp",PacPidCand2Ntuple::packFloat(theElTag->p(),8));
	  _jpsikeeTuple->column("othertheta",PacPidCand2Ntuple::packFloat(theElTag->p4().theta(),10));
	  _jpsikeeTuple->column("otherphi",PacPidCand2Ntuple::packFloatPhi(theElTag->p4().phi(),10));
	  _jpsikeeTuple->column("othercharge",theElTag->charge());
	  _jpsikeeTuple->column("otheriselh", PacPidCand2Ntuple::packFloat(otheriselh,8));
	  _jpsikeeTuple->column("kmass",PacPidCand2Ntuple::packFloat(theKst2->mass(),8));
	  _jpsikeeTuple->column("kmode",3);
	  _jpsikeeTuple->column("kp",PacPidCand2Ntuple::packFloat(theKst2->p(),8));
	  _jpsikeeTuple->column("mPsi",PacPidCand2Ntuple::packFloat(theCC2->mass(),8));
	  _jpsikeeTuple->column("mES",PacPidCand2Ntuple::packFloat(m_ES,6));
	  _jpsikeeTuple->column("deltaE",PacPidCand2Ntuple::packFloat(deltaE,8));
        // _jpsikeeTuple->column("R2",r2);
	  _jpsikeeTuple->column("nTrk",nTracks);

	  // NA
	  //_jpsikeeTuple->column("helicity",PacPidCand2Ntuple::packFloat((float) Helicity,10));
	  
	  filler->convert(&theLocalTrk2,anEvent,_jpsikeeTuple);
	  filler->convertJpsiKMC(&theLocalTrk2,theElTag,theKst2,anEvent,_jpsikeeTuple);
	  _jpsikeeTuple->dumpData();
	  
	  delete theB2c;
	  
	  // NA
	  //delete theHelicity;
	}
	delete theCC2;
      }
  }
} // end J/psi e sample





//////////////////////////////////////////// 
void PacPidPCNtupleWriter::extractBtoDpi(AbsEvent *anEvent)
{
  AbsEventTag* tag = Ifd<AbsEventTag>::get( anEvent );
  
  int nTracks(0);
  float r2All(0), thrustCosTh(0);
  bool bgfmultihadron;
  bool status;
  
  if ( tag != 0 )
    {
      status  = tag->getBool( bgfmultihadron, "BGFMultiHadron" );
      status &= tag->getInt( nTracks, "nTracks" );
      status &= tag->getFloat( r2All, "R2All");
      status &= tag->getFloat(thrustCosTh,"thrustCosTh");
      if((!status) || (nTracks < 4) || (r2All > 0.7) || (bgfmultihadron == false))
	return;
    }
  
  HepAList< BtaCandidate >* bachPiList;
  HepAList< BtaCandidate >* piList;
  HepAList< BtaCandidate >* kList;
  
  static const IfdStrKey key_piLHTight("piLHTight");
  static const IfdStrKey keyGTVL("GoodTracksVeryLoose");

  getTmpAList(anEvent, bachPiList, key_piLHTight );
  getTmpAList(anEvent, piList, keyGTVL );
  getTmpAList(anEvent, kList, keyGTVL );
  
  HepAListIterator<BtaCandidate> itPi( *piList );
  HepAListIterator<BtaCandidate> itK( *kList );
  HepAListIterator<BtaCandidate> itBachPi( *bachPiList );

  BtaCandidate *thePi(0), *theK(0), *theBachPi(0);
  HepLorentzVector BP4;
  
  BtaOpFastVtx theCombiner;
  BtaOpFastVtx theBCombiner;
  
  while(thePi = itPi()) {
    itK.rewind();
    while(theK = itK())
      {
	if(theK->charge() + thePi->charge() != 0)
	  continue;
	if(theK->overlaps(*thePi))
	  continue;
	
	BtaCandidate theLocalK( *theK );
	BtaCandidate theLocalPi( *thePi );
	
	if(theLocalK.charge() == 1)
	  {
	    theLocalK.setType(Pdt::lookup("K+"));
	    theLocalPi.setType(Pdt::lookup("pi-"));
	  }
	else if(theLocalK.charge() == -1)
	  {
	    theLocalK.setType(Pdt::lookup("K-"));
	    theLocalPi.setType(Pdt::lookup("pi+"));
	  }
	
	BtaCandidate *theD = theCombiner.create(theLocalPi,theLocalK);
	
      // Check for a valid vertex
	if(theD->decayVtx() == 0)
	  {
	    delete theD;
	    continue;
	  }
	
	if(theD->mass() < 1.82 || theD->mass() > 1.90)
	  {
	    delete theD;
	    continue;
	  }
	
	itBachPi.rewind();
	while(theBachPi = itBachPi()) {
	  if(theBachPi->charge() == thePi->charge())
	    continue;
	  
	  if(theBachPi->overlaps(*thePi) || theBachPi->overlaps(*theK))
	    continue;
	  
	  BtaCandidate *theB =
	    theBCombiner.create(theLocalK,theLocalPi,*theBachPi);
	  
        // Check for a valid vertex
	  if(theB->decayVtx() == 0)
	    {
	      delete theB;
	      continue;
	    }
	  
	  BP4 = theB->p4();
	  BtaBVariables bVars(BP4, true);
	  
	  float m_ES = bVars.m_ES();
	  float deltaE = bVars.deltaE();
	  
	  if(m_ES < 5.2 || deltaE < -0.15 || deltaE > 0.15 || m_ES > 5.3)
	    {
	      delete theB;
	      continue;
	    }
	  
	  addVtxQuantities(theB,anEvent,_btodkaonTuple);
	  addVtxQuantities(theB,anEvent,_btodpionTuple);
	  
	  _btodkaonTuple->column("mD0",PacPidCand2Ntuple::packFloat(theD->mass(),8));
	  _btodkaonTuple->column("mES",PacPidCand2Ntuple::packFloat(m_ES,6));
	  _btodkaonTuple->column("deltaE",PacPidCand2Ntuple::packFloat(deltaE,8));
	  _btodkaonTuple->column("R2All",PacPidCand2Ntuple::packFloat(r2All,10));
	  _btodkaonTuple->column("cosThrust",PacPidCand2Ntuple::packFloat(thrustCosTh,10));
	  _btodkaonTuple->column("nTrk",nTracks);
	  
	  filler->convert(&theLocalK,anEvent,_btodkaonTuple);
	  
	  _btodkaonTuple->dumpData();
	  
	  _btodpionTuple->column("mD0",PacPidCand2Ntuple::packFloat(theD->mass(),8));
	  _btodpionTuple->column("mES",PacPidCand2Ntuple::packFloat(m_ES,6));
	  _btodpionTuple->column("deltaE",PacPidCand2Ntuple::packFloat(deltaE,8));
	  _btodpionTuple->column("R2All",PacPidCand2Ntuple::packFloat(r2All,10));
	  _btodpionTuple->column("cosThrust",PacPidCand2Ntuple::packFloat(thrustCosTh,10));
	  _btodpionTuple->column("nTrk",nTracks);
	  
	  filler->convert(&theLocalPi,anEvent,_btodpionTuple);
	  
	  _btodpionTuple->dumpData();
	  
	  delete theB;
	}
	delete theD;
      }
  }
} // end B -> D pi sample


//////////////////////////////////////////////////////////////////////////////////
// extractDtoKpipi method for D+ to K-pi+pi+ control sample
// Written by Ryan Mackenzie White
void PacPidPCNtupleWriter::extractDtoKpipi(AbsEvent *anEvent)
{
  AbsEventTag* tag = Ifd<AbsEventTag>::get( anEvent );
  
  int nTracks(0);
  float R2CTrk(0), R2All(0), ThrustMagCTrk(0), ThrustMagAll(0), SphericityAll(0);
  bool bgfmultihadron;
  bool status;
  
//Require good events: multihadron events
//R2 > 0.2; Sphericity < 0.5; Thrust > 0.7; #tracks > 3
  if ( tag != 0 )
    {
      status  = tag->getBool( bgfmultihadron, "BGFMultiHadron" );
      status &= tag->getInt( nTracks, "nTracks" );
      status &= tag->getFloat( R2All, "R2All");
      status &= tag->getFloat( R2CTrk, "R2");
      status &= tag->getFloat(ThrustMagAll,"thrustMagAll");
      status &=tag->getFloat(ThrustMagCTrk, "thrustMag");
      status &=tag->getFloat(SphericityAll, "sphericityAll");
      if((!status) || (nTracks < 3) || ( (R2All < 0.2) && (R2CTrk < 0.2) &&
					 (SphericityAll > 0.5) && (ThrustMagAll < 0.7) &&
					 (ThrustMagCTrk < 0.7) ) || (bgfmultihadron == false))
        return;
    }
  
  
  
  HepAList< BtaCandidate >* DplusList;
  static const IfdStrKey keyKPiPi("BtaPidDplusToKPiPi");
  getTmpAList(anEvent, DplusList, keyKPiPi );
  
  if (DplusList==0) { return; }
  
  HepAListIterator<BtaCandidate> candidateListIter( *DplusList );
  
  BtaCandidate *theDplus;
  while(theDplus = candidateListIter()) {
    HepAListIterator<BtaCandidate> dauIter(theDplus->daughterIterator());
    BtaCandidate * theK = dauIter.next();
    BtaCandidate * thePi1 = dauIter.next();
    BtaCandidate * thePi2 = dauIter.next();
    
    double likelihoodBeam = _dtagfitter->doBeamFit(theDplus, anEvent);
    if(likelihoodBeam < 1) continue;
    double likelihoodDalitz = _dtagfitter->doDalitzFit(theDplus, anEvent);
    
  //Get the CM momentum of the D
    HepLorentzVector p4D(theDplus->p4());
    Hep3Vector  UpsiBoost;
    const PepBeams *theBeams = gblEnv->getPep()->pepBeams();
    if (theBeams != 0) {
      UpsiBoost = theBeams->boostCMtoLab();
    } else {
      UpsiBoost = Hep3Vector(0.,0.,0.486976);
    }
    p4D.boost(-UpsiBoost);
    
    double theDplusVtxProb = probab(theDplus->decayVtx()->nDof(),theDplus->decayVtx()->chiSquared());
    
    BtaOpAdd4 op4;
    BtaCandidate kpi1Cand = op4.combine(*theK,*thePi1);
    float kpi1Mass = kpi1Cand.mass();
    
    BtaCandidate kpi2Cand = op4.combine(*theK,*thePi2);
    float kpi2Mass = kpi2Cand.mass();
    
    BtaCandidate pi1pi2Cand = op4.combine(*thePi1,*thePi2);
    float pi1pi2Mass = pi1pi2Cand.mass();
    
    
    addVtxQuantities(theDplus,anEvent,_dtokaonTuple);
    _dtokaonTuple->column("dpcms",PacPidCand2Ntuple::packFloat((float) p4D.rho(),10));
    _dtokaonTuple->column("vtxprobD",PacPidCand2Ntuple::packFloat((float) theDplusVtxProb,12));
    _dtokaonTuple->column("likelihdDalitz",PacPidCand2Ntuple::packFloat((float) likelihoodDalitz,8));
    _dtokaonTuple->column("likelihdBeam",PacPidCand2Ntuple::packFloat((float) likelihoodBeam,8));
    _dtokaonTuple->column("nTrk",nTracks);
    _dtokaonTuple->column("kpi1Mass",PacPidCand2Ntuple::packFloat((float) kpi1Mass,8));
    _dtokaonTuple->column("kpi2Mass",PacPidCand2Ntuple::packFloat((float) kpi2Mass,8));
    _dtokaonTuple->column("pi1pi2Mass",PacPidCand2Ntuple::packFloat((float) pi1pi2Mass,8));
    filler->convert(theK, anEvent, _dtokaonTuple);
    _dtokaonTuple->dumpData();
    
    addVtxQuantities(theDplus,anEvent,_dtopion1Tuple);
    _dtopion1Tuple->column("dpcms",PacPidCand2Ntuple::packFloat((float) p4D.rho(),10));
    _dtopion1Tuple->column("vtxprobD",PacPidCand2Ntuple::packFloat((float) theDplusVtxProb,12));
    _dtopion1Tuple->column("likelihdDalitz",PacPidCand2Ntuple::packFloat((float) likelihoodDalitz,8));
    _dtopion1Tuple->column("likelihdBeam",PacPidCand2Ntuple::packFloat((float) likelihoodBeam,8));
    _dtopion1Tuple->column("nTrk",nTracks);
    _dtopion1Tuple->column("kpi1Mass",PacPidCand2Ntuple::packFloat((float) kpi1Mass,8));
    _dtopion1Tuple->column("kpi2Mass",PacPidCand2Ntuple::packFloat((float) kpi2Mass,8));
    _dtopion1Tuple->column("pi1pi2Mass",PacPidCand2Ntuple::packFloat((float) pi1pi2Mass,8));
    filler->convert(thePi1, anEvent, _dtopion1Tuple);
    _dtopion1Tuple->dumpData();
    
    addVtxQuantities(theDplus,anEvent,_dtopion2Tuple);
    _dtopion2Tuple->column("dpcms",PacPidCand2Ntuple::packFloat((float) p4D.rho(),10));
    _dtopion2Tuple->column("vtxprobD",PacPidCand2Ntuple::packFloat((float) theDplusVtxProb,12));
    _dtopion2Tuple->column("likelihdDalitz",PacPidCand2Ntuple::packFloat((float) likelihoodDalitz,8));
    _dtopion2Tuple->column("likelihdBeam",PacPidCand2Ntuple::packFloat((float) likelihoodBeam,8));
    _dtopion2Tuple->column("nTrk",nTracks);
    _dtopion2Tuple->column("kpi1Mass",PacPidCand2Ntuple::packFloat((float) kpi1Mass,8));
    _dtopion2Tuple->column("kpi2Mass",PacPidCand2Ntuple::packFloat((float) kpi2Mass,8));
    _dtopion2Tuple->column("pi1pi2Mass",PacPidCand2Ntuple::packFloat((float) pi1pi2Mass,8));
    filler->convert(thePi2, anEvent, _dtopion2Tuple);
    _dtopion2Tuple->dumpData();
    
    delete thePi1;
    delete thePi2;
    delete theK;
    
  } // Dplus iterator
  delete theDplus;
} //End method to D to Kpipi



void PacPidPCNtupleWriter::extractLambdaCtopKpi(AbsEvent *anEvent){
  AbsEventTag* tag = Ifd<AbsEventTag>::get( anEvent );
  
  int nTracks(0);
  float R2CTrk(0), R2All(0), ThrustMagCTrk(0), ThrustMagAll(0), SphericityAll(0);
  bool bgfmultihadron;
  bool status;
  
//Require good events: multihadron events
//R2 > 0.2; Sphericity < 0.5; Thrust > 0.7; #tracks > 3
  if ( tag != 0 )
    {
      status  = tag->getBool( bgfmultihadron, "BGFMultiHadron" );
      status &= tag->getInt( nTracks, "nTracks" );
      status &= tag->getFloat( R2All, "R2All");
      status &= tag->getFloat( R2CTrk, "R2");
      status &= tag->getFloat(ThrustMagAll,"thrustMagAll");
      status &=tag->getFloat(ThrustMagCTrk, "thrustMag");
      status &=tag->getFloat(SphericityAll, "sphericityAll");
      if((!status) || (nTracks < 3) || ( (R2All < 0.2) && (R2CTrk < 0.2) &&
					 (SphericityAll > 0.5) && (ThrustMagAll < 0.7) &&
					 (ThrustMagCTrk < 0.7) ) || (bgfmultihadron == false))
	return;
    }

  HepAList< BtaCandidate >* piList;
  HepAList< BtaCandidate >* LambdaCList;

  static const IfdStrKey key_piLHTight("piLHTight");
  static const IfdStrKey key_pKpi("BtaPCPidLambdaCTopKpi");
  getTmpAList(anEvent, piList, key_piLHTight );
  getTmpAList(anEvent, LambdaCList, key_pKpi );

  if (LambdaCList==0) { return; }

  HepAListIterator<BtaCandidate> itPi( *piList );

  HepAListIterator<BtaCandidate> candidateListIter( *LambdaCList );

  //Get the event boost
  Hep3Vector  UpsiBoost;
  const PepBeams *theBeams = gblEnv->getPep()->pepBeams();
  if (theBeams != 0) {
    UpsiBoost = theBeams->boostCMtoLab();
  } else {
    UpsiBoost = Hep3Vector(0.,0.,0.486976);
  }

  BtaCandidate *thePi(0), *theLambdaC(0);

  while(thePi = itPi()) {
    candidateListIter.rewind();
    while(theLambdaC = candidateListIter()) {
      BtaOpAdd4 theCombiner;
      BtaCandidate *theSigmaC = theCombiner.create(*thePi, *theLambdaC);
      
      //Get the CM momentum of the SigmaC
      HepLorentzVector p4SigmaC(theSigmaC->p4());
      p4SigmaC.boost(-UpsiBoost);
      
      if(p4SigmaC.rho() < 3.2) 
	{
	  delete theSigmaC;
	  continue;
	}

      double massDif = theSigmaC->mass() - theLambdaC->mass();

      
      if(massDif < 0.160 || massDif > 0.175)
	{
	  delete theSigmaC;
	  continue;
	}
      
      
      HepAListIterator<BtaCandidate> dauIter(theLambdaC->daughterIterator());
      BtaCandidate * theProton = dauIter.next();
      BtaCandidate * theK = dauIter.next();
      BtaCandidate * theLambdaPi = dauIter.next();

      if(theProton->charge()==thePi->charge() && theProton->overlaps(*thePi))
	{
	  delete theSigmaC;
	  continue;
	}
      if(theLambdaPi->charge()==thePi->charge() && theLambdaPi->overlaps(*thePi))
	{
	  delete theSigmaC;
	  continue;
	}

      if(theK->charge()==thePi->charge() && theK->overlaps(*thePi))
	{
	  delete theSigmaC;
	  continue;
	}

      //Need to check for reflections from D, Ds, D*
      HepLorentzVector theFakeProtonK(theProton->p4().x(), theProton->p4().y(), theProton->p4().z(), sqrt(0.49368*0.49368 + theProton->p4().rho()*theProton->p4().rho()));
      HepLorentzVector theFakeProtonPi(theProton->p4().x(), theProton->p4().y(), theProton->p4().z(), sqrt(0.139570*0.139570 + theProton->p4().rho()*theProton->p4().rho()));

      HepLorentzVector theRealK(theK->p4());
      HepLorentzVector theRealPi(theLambdaPi->p4());

      HepLorentzVector theFakeLambdaToKpipi = theFakeProtonPi + theRealK + theRealPi;
      double KpipiMass = theFakeLambdaToKpipi.mag();

      HepLorentzVector theKpi = theRealK + theFakeProtonPi;
      double theKpiMass = theKpi.mag();
      double theSlowPiMass = theRealPi.mag();
      double theSlowPiP = theRealPi.rho();
      double Qval = KpipiMass - theKpiMass - theSlowPiMass;

      HepLorentzVector theFakeLambdaToKKpi = theFakeProtonK + theRealK + theRealPi;
      double KKpiMass = theFakeLambdaToKKpi.mag();
      
      //LamdaC variables
      HepLorentzVector p4Lc2(theLambdaC->p4());
      p4Lc2.boost(-UpsiBoost);

      addVtxQuantities(theLambdaC,anEvent,_lambdactopkpiTuple);
      _lambdactopkpiTuple->column("LambdaCpcms",PacPidCand2Ntuple::packFloat((float) p4Lc2.rho(),12));
      _lambdactopkpiTuple->column("SigmaCpcms",PacPidCand2Ntuple::packFloat((float) p4SigmaC.rho(),12));
      _lambdactopkpiTuple->column("vtxmassSigmaC",PacPidCand2Ntuple::packFloat((float) theSigmaC->p4().mag(),8));
      _lambdactopkpiTuple->column("massDiff",PacPidCand2Ntuple::packFloat((float) massDif,8));
      _lambdactopkpiTuple->column("KKpiMass",PacPidCand2Ntuple::packFloat((float) KKpiMass,8));
      _lambdactopkpiTuple->column("KpipiMass",PacPidCand2Ntuple::packFloat((float) KpipiMass,8));
      _lambdactopkpiTuple->column("KpiMass",PacPidCand2Ntuple::packFloat((float) theKpiMass,8));
      _lambdactopkpiTuple->column("Qvalue",PacPidCand2Ntuple::packFloat((float) Qval,8));
      _lambdactopkpiTuple->column("pisoftmass",PacPidCand2Ntuple::packFloat((float) theSlowPiMass,8));
      _lambdactopkpiTuple->column("pisoftp",PacPidCand2Ntuple::packFloat((float) theSlowPiP,8));
      filler->convert(theProton, anEvent,_lambdactopkpiTuple);
      _lambdactopkpiTuple->dumpData();
      delete theSigmaC;
    } //End Lambda list
  }
} // end of Lambda_c -> p K pi sample



//////////////////////////////////////
void PacPidPCNtupleWriter::printConfig()
{
  ErrMsg(routine) << endmsg;
  ErrMsg(routine) << "  +============================================================+" << endmsg;
  ErrMsg(routine) << "  |           Configuration of PacPidPCNtupleWriter               |" << endmsg;
  ErrMsg(routine) << "  +============================================================+" << endmsg;
  ErrMsg(routine) << "  | Electrons :                                                |" << endmsg;
  bool tt = _doTwoTrk->value();
  bool muHa = _doMuha->value();
  cout.flags(ios::left);

  if (tt || _doBhabha->value()) {
    IfdStrKey& k = (IfdStrKey&) _trkBhabhaList.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k.asString() << "ntp200|" << endmsg;
    IfdStrKey& k2 = (IfdStrKey&) _ifrBhabhaList.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k2.asString() << "ntp201|" << endmsg;
    IfdStrKey& k3 = (IfdStrKey&) _radBhabhaList.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k3.asString() << "ntp202|" << endmsg;
  }
  if (tt || _do4e->value()) {
    IfdStrKey& k = (IfdStrKey&) _4eList.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k.asString() << "ntp203|" << endmsg;
    IfdStrKey& k2  = (IfdStrKey&) _4eListO.value();
    ErrMsg(routine) << "  |       " << setw(47) << k2.asString() << "      |" << endmsg;
  }
  if (tt || _doVcs->value()) {
    IfdStrKey& k = (IfdStrKey&) _vcsListName.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k.asString() << "ntp206|" << endmsg;
    IfdStrKey& k2  = (IfdStrKey&) _vcsGammaListName.value();
    ErrMsg(routine) << "  |       " << setw(47) << k2.asString() << "      |" << endmsg;
  }


  if (muHa || _doConv->value()) {
      IfdStrKey& k1 = (IfdStrKey&) _convElectronList.value();
      IfdStrKey& k2 = (IfdStrKey&) _convGammaList.value();
      ErrMsg(routine) << "  |    *  " << setw(47) << k1.asString() << "ntp206|" << endmsg;
      ErrMsg(routine) << "  |       " << setw(47) << k2.asString() << "      |" << endmsg;
  }
  ErrMsg(routine) << "  |                                                            |" << endmsg;
  ErrMsg(routine) << "  | Muons :                                                    |" << endmsg;
  if (tt || _doMumugamma->value()) {
    IfdStrKey& k1 = (IfdStrKey&) _mumuGammaList.value();
    IfdStrKey& k2 = (IfdStrKey&) _mumuGammaListOther.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k1.asString() << "ntp204|" << endmsg;
    ErrMsg(routine) << "  |       " << setw(47) << k2.asString() << "      |" << endmsg;
    IfdStrKey& k3 = (IfdStrKey&) _mumuGamma2List.value();
    IfdStrKey& k4 = (IfdStrKey&) _mumuGamma2ListOther.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k3.asString() << "ntp207|" << endmsg;
    ErrMsg(routine) << "  |       " << setw(47) << k4.asString() << "      |" << endmsg;
  }
  if (tt || _doEemumu->value()) {
    IfdStrKey& k1 = (IfdStrKey&) _eemumuList.value();
    IfdStrKey& k2 = (IfdStrKey&) _eemumuListOther.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k1.asString() << "ntp205|" << endmsg;
    ErrMsg(routine) << "  |       " << setw(47) << k2.asString() << "      |" << endmsg;

  }

  ErrMsg(routine) << "  |                                                            |" << endmsg;
  ErrMsg(routine) << "  | Pions :                                                    |" << endmsg;
  if (muHa || _doKs->value()) {
    IfdStrKey& k = (IfdStrKey&) _ksList.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k.asString() << "ntp104|" << endmsg;
  }
  if (muHa || _doDstar->value()) {
    IfdStrKey& k = (IfdStrKey&) _dstarList.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k.asString() << "ntp101|" << endmsg;
    ErrMsg(routine) << "  |    *  " << setw(47) << k.asString() << "ntp106|" << endmsg;
  }
  if (_doTau31->value()) {
    IfdStrKey& k1 = (IfdStrKey&) _tau31piList.value();
    IfdStrKey& k2 = (IfdStrKey&) _tau31OneProngList.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k1.asString() << "ntp300|" << endmsg;
    ErrMsg(routine) << "  |       " << setw(47) << k2.asString() << "      |" << endmsg;
  }
  ErrMsg(routine) << "  |                                                            |" << endmsg;
  ErrMsg(routine) << "  | Kaons :                                                    |" << endmsg;
  if (muHa || _doDstar->value()) {
    IfdStrKey& k = (IfdStrKey&) _dstarList.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k.asString() << "ntp102|" << endmsg;
  }
  ErrMsg(routine) << "  |                                                            |" << endmsg;
  ErrMsg(routine) << "  | Protons :                                                  |" << endmsg;
  if (muHa || _doLambda->value()) {
    IfdStrKey& k = (IfdStrKey&) _lambdaList.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k.asString() << "ntp100|" << endmsg;
  }
  if (muHa || _doTrkLambda->value()) {
    IfdStrKey& k = (IfdStrKey&) _trklambdaList.value();
    ErrMsg(routine) << "  |    *  " << setw(47) << k.asString() << "ntp105|" << endmsg;
  }

  ErrMsg(routine) << "  +------------------------------------------------------------+" << endmsg;
  ErrMsg(routine) << "  +  Common content of ntuples                                 +" << endmsg;
  ErrMsg(routine) << "  +------------------------------------------------------------+" << endmsg;
  ErrMsg(routine) << "  |      * Basic kinematic quantities (momentum,POCA,q)        |" << endmsg;
  if (filler->getParmValue(string("writeMcTruth"))) {
    ErrMsg(routine) << "  |      * MC Truth (machted track , mother )                  |" << endmsg;
  }
  if (filler->getParmValue(string("writeElectronCms"))) {
    ErrMsg(routine) << "  |      * CMS momentum under electron hypothesis              |" << endmsg;
  }
  if (filler->getParmValue(string("writeCalQual"))) {
    ErrMsg(routine) << "  |      * Quantities from BtaCalQual                          |" << endmsg;
  }
  if (filler->getParmValue(string("writeTrkQual"))) {
    ErrMsg(routine) << "  |      * Quantities from BtaTrkQual                          |" << endmsg;
  }
  if (filler->getParmValue(string("writePidQual"))) {
    ErrMsg(routine) << "  |      * Quantities from BtaPidQual                          |" << endmsg;
  }
  if (filler->getParmValue(string("writeIfrQual"))) {
    ErrMsg(routine) << "  |      * Quantities from BtaIfrQual                          |" << endmsg;
  }
  if (filler->getParmValue(string("writePidInfo"))) {
    ErrMsg(routine) << "  |      * Quantities from BtaPidInfo                          |" << endmsg;
  }
  if (filler->getParmValue(string("writeExpectedDedx"))) {
    ErrMsg(routine) << "  |      * expected DCH-de/dx from Cond-DB                     |" << endmsg;
  }
  if (filler->getParmValue(string("writeSMSprobs"))) {
    ErrMsg(routine) << "  |      * probabilities from SMS-selector                     |" << endmsg;
  }
  if (filler->getParmValue(string("writeSelectorBits"))) {
    ErrMsg(routine) << "  |      * Selector Bits :                                     |" << endmsg;
    const vector<string> &keyList = _pidLists->value();
    for (int i=0; i<keyList.size(); i++) {
      ErrMsg(routine) << "  |        - " << setw(44) << keyList[i] << "      |" << endmsg;
    }
  }
  if (filler->getParmValue(string("writeTagBits"))) {
    ErrMsg(routine) << "  |      * Tag Bits                                            |" << endmsg;
    const vector<string> &keyList = _tagBits->value();
    for (int i=0; i<keyList.size(); i++) {
      ErrMsg(routine) << "  |        - " << setw(44) << keyList[i] << "      |" << endmsg;
    }
  }
  ErrMsg(routine) << "  +============================================================+\n" << endmsg;
}

///////////////////////////////////////////////////////////////////////////////////
AppResult PacPidPCNtupleWriter::endJob( AbsEvent* anEvent ) //
{
  ErrMsg(routine) << name( ) << " end Job" << endmsg;
// Workaround for crashes at endJob. Saves our ntuples!
//  HepTupleManager* manager = gblEnv->getGen()->ntupleManager();
//  assert(manager != 0);
//  manager->write();
//  manager->resetAll();
//  ErrMsg(warning) << name( ) << " endJob closed hbook file!" << endmsg;
  return AppResult::OK;
}


