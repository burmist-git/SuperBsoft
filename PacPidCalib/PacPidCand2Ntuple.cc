//--------------------------------------------------------------------------
// File and Version Information:
//  $Id: $
//
// Description:
//  Dump Micro information to an ntuple.
//  Adapted from BetaPidCalibNtuple/BtaCand2Ntuple.hh
//
// Environment:
//  Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//   Nicolas ARNAUD (SuperB)
//   Dave Aston, Thorsten Brandt, Kevin Flood, Alexandre Telnov, ... 
//
// Copyright
//   (C) 2009 CNRS-IN2P3
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "PacPidCalib/PacPidCand2Ntuple.hh"

// moved here from class header file
#include "AbsParm/AbsParmIfdStrKey.hh"

// NA
//#include "TrgTools/TrgFctTimePointInspector.hh"

#include "AbsEventTag/AbsEventTag.hh"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEvent/AbsEventID.hh"
#include "AbsEvent/getTmpAList.hh"

#include "IfrPidUtils/IfrMicroPidInfo.hh"
#include "DchCalib/DchBetheBloch.hh"
#include "ErrLogger/ErrLog.hh"
#include "AbsParm/AbsParmVector.hh"
#include "Framework/AbsParmBool.hh"
#include "Framework/AbsParmGeneral.hh"
#include "Framework/AppModule.hh"

#include "HepTuple/Tuple.h"

#include "CLHEP/Alist/AList.h"
#include "CLHEP/Alist/AIterator.h"

#include "Beta/EventInfo.hh"
#include "Beta/BtaCandidate.hh"
#include "BetaRecoAdapter/BtaAbsRecoObject.hh"
#include "BetaMicroAdapter/BtaMicroAdapter.hh"
#include "BetaMicroAdapter/BtaTrkQual.hh"
#include "BetaMicroAdapter/BtaCalQual.hh"
#include "BetaMicroAdapter/BtaIfrQual.hh"
#include "BetaMicroAdapter/BtaPidQual.hh"
#include "BetaMicroAdapter/BtaPidInfo.hh"
#include "BetaTools/BtaMcAssocChiSq.hh"
#include "BetaCoreTools/BtaBooster.hh"

// NA
//#include "BetaPid/PidMuonSPR.hh"
//#include "BetaPid/PidLHRatios.hh"
//#include "BetaPid/PidKaonSMSSelector.hh"
//#include "BetaPid/PidKaonMicroSelector.hh"
//#include "BetaPid/PidKaonBDTSelector.hh"
//#include "BetaPid/PidKMSelector.hh"
//#include "BetaPid/PidELHRatios.hh"

#include "AbsEnv/AbsEnv.hh"
#include "GenEnv/GenEnv.hh"
#include "BtaEnv/BtaEnv.hh"
#include "PepEnvData/PepBeams.hh"
#include "OdfCommon/odfTime.hh"
#include "EidData/EidEventTriplet.hh"
#include "CommonUtils/ComPathNameSearch.hh"
#include "IfrPidUtils/IfrNNKernelSetup.hh"

#include "EidData/EidCondKeyTriplet.hh"

// NA
//#include "L1FctData/L1FctTimePoint.hh"
//#include "OprMonBase/OprTrickleRegions.hh"

#include "PDT/PdtEntry.hh"
#include "PDT/Pdt.hh"
#include "ProbTools/NumRecipes.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

// Detailed DCH and SVT dE/dx parameterizations, first from ASCII files
// NA
//#include "PidDchSvtDrcCalib/PidConvertMajorID.hh"
//#include "PidDchSvtDrcCalib/PidExpectedBetheBlochDch_AVT.hh"
//#include "PidDchSvtDrcCalib/PidExpectedBetheBlochSvt_AVT.hh"
//#include "PidDchSvtDrcCalib/PidBetheBlochErrorDch_AVT.hh"
//#include "PidDchSvtDrcCalib/PidBetheBlochErrorSvt_AVT.hh"
//#include "PidDchSvtDrcCalib/PidBetheBlochPhiCorrDch_AVT.hh"
//#include "PidDchSvtDrcCalib/PidBetheBlochPhiCorrSvt_AVT.hh"
// and now also from CDB! (new in October 2007)
// NA
//#include "BetaPid/PidDEDXCdbDch.hh"
//#include "BetaPid/PidDEDXCdbSvt.hh"

// avtelnov, 2008/02/10: a tool to tell if a track from this event was used 
// in the training of the given BDT classifier. Compute and store a bool
// for each of the six classifiers (KMe, KMpi, KMk, KMp, MUmu, MUpi);
// use it later in the making of PID performance plots and tables 
// (BetaPidUser to be modified accordingly at a later time)
// NA
//#include "BetaPidCalibNtuple/UsedForBdtTraining.hh" 

#include "BetaMini/BetaMiniTools.hh"

// Added to accomodate the 2007 DCH dE/dx developments
#include "AbsPid/PidInfoSummary.hh"
#include "DchPid/DchPidInfo.hh"
#include "SvtPid/SvtPidInfo.hh"

// Added at Al Eisner's request
//Merged pi0 stuff
// NA
//#include "EmcPid/EmcMergedPi0Identifier.hh"
//#include "BetaPid/BtaMergedPi0Algo.hh"

#include "EmcPid/EmcIDFactory.hh"
// Consistency Stuff
#include "ProbTools/Consistency.hh"

// NA
#include "PacPid/PacPidUtils.hh"

// NA
#include "Framework/AbsParmBool.hh"

using std::string;
using std::fstream;
using std::ostream;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// bypass packing if set to false
AbsParmBool* PacPidCand2Ntuple::_doFloatPacking = 0;
// remove no more bits than the value of this parameter
AbsParmGeneral<int>* PacPidCand2Ntuple::_maxFloatPacking = 0;

// NA
//AbsParmBool* PacPidCand2Ntuple::_useCDBforDEDX = 0;

#if defined(__SUNPRO_CC)

// Needed by buggy Sun compiler so that the template constructor will be
// instantiated.  Otherwise, the linker can't find the template.  This function
// does nothing and is never called; its only purpose is to get around the
// Sun compiler bug.

static void dummy()
{
  AppModule dummy2("dummy", "unused");
  PacPidCand2Ntuple dummy1(&dummy2);
  PidMuonSPR unused(&dummy1, &dummy2);
}
#endif

// BEWARE! The maximum size of the PidLists arrays is fixed here (currently at 140)
// dyaeb 20060605 -- we have just exceeded 60!
// anticipating more selectors in the near future, increased max from 80 -> 100, 17 nov 06 (kflood)
// avtelnov, Aug 2007: with all the new BDT selectors, the number of selectors is precisely 100! 
// avtelnov, Oct 2, 2007: increase to 140 to squeeze in the ChargedTracks flag, the combined lists, 
// and to leave room for the future

PacPidCand2Ntuple::PacPidCand2Ntuple(AppModule *aModule)
  :_listOfPidLists(140),
   _namesOfPidLists(140),
   _namesOfTagBits(),
   _mcTruthList("MCTruthList",aModule,"GHit")

   // NA
   //_muonTree(this, aModule), // for the muon BDT selector
				    //_muonLoPTree(this, aModule),  // for the muon BDT selector

   // NA
   //_getTrickleInfo("GetTrickleInfo",aModule,false),
   //_timePointLER("TimePointLER",aModule,"LERInjectionTimePoint"),
   //_timePointHER("TimePointHER",aModule,"HERInjectionTimePoint")
   
   // NA
   //_MergedPi0Workhorse(0),
   //_MergedPi0algo(0)
  
{
  
  _doFloatPacking = new AbsParmBool("doFloatPacking",aModule,false);
  _maxFloatPacking = new AbsParmGeneral<int>("maxFloatPacking",aModule,8);

// set to false to use detailed dE/dx from the ASCII files in PidDchSvtDrcCalib
  // NA
  //_useCDBforDEDX = new AbsParmBool("useCDBforDEDX",aModule,true); 
  
  aModule->commands()->append(&_mcTruthList);
  aModule->commands()->append(_doFloatPacking);
  aModule->commands()->append(_maxFloatPacking);
  
  // NA
  //aModule->commands()->append(_useCDBforDEDX);
  
  // NA
  //aModule->commands()->append(&_getTrickleInfo);
  //aModule->commands()->append(&_timePointLER);
  //aModule->commands()->append(&_timePointHER);
  
  _theAssoc=NULL;
  
  // NA
  //_theSms = new PidKaonSMSSelector();
  //_theNN = new PidKaonMicroSelector();
  //_theKaonBDT = new PidKaonBDTSelector();
  //_theKM      = new PidKMSelector();

  _ifrNNKernel=0;
  addNewParmBool(string("writeTagBits"),true,aModule);
  addNewParmBool(string("writeSelectorBits"),true,aModule);
  addNewParmBool(string("writeMcTruth"),false,aModule);
  addNewParmBool(string("writeElectronCms"),false,aModule);
  addNewParmBool(string("writeCalQual"),true,aModule);
  addNewParmBool(string("writeTrkQual"),true,aModule);
  addNewParmBool(string("writePidQual"),true,aModule);
  
  // NA: default value set to false
  addNewParmBool(string("writeIfrQual"),false,aModule);
  
  addNewParmBool(string("writePidInfo"),true,aModule);

  // NA: default value set to false
  addNewParmBool(string("writeExpectedDedx"),false,aModule);

  addNewParmBool(string("writeSMSprobs"),true,aModule);
  addNewParmBool(string("writeKLHratios"),true,aModule);
  addNewParmBool(string("writeELHratios"),true,aModule);

  // NA
  //_trickleRegion_LER = 0;
  //_trickleRegion_HER = 0;

  _warningCounter = 0;
}


//--------------
// Destructor --
//--------------
PacPidCand2Ntuple::~PacPidCand2Ntuple()
{
  // NA
  //delete _theSms;
  //delete _theNN;
  //delete _theKaonBDT;
  //delete _theKM;
  
  if ( _ifrNNKernel ) delete _ifrNNKernel;

  // NA
  //if (_trickleRegion_LER) delete _trickleRegion_LER;
  //if (_trickleRegion_HER) delete _trickleRegion_HER;
  
  // NA
  //if (_MergedPi0Workhorse) delete _MergedPi0Workhorse;
  //if (_MergedPi0algo) delete _MergedPi0algo;

}


//////////////////////////////////////////////////////////////////////////////
void
PacPidCand2Ntuple::beginJob(const AbsEvent *ev)
{
  // NA
  /*
  if (_getTrickleInfo.value())  {
  _timePointInspector_LER.use( _timePointLER.valueString() );
  _timePointInspector_HER.use( _timePointHER.valueString() );
  }
  */

  // NA
  /*
  if (_useCDBforDEDX->value() == false) {
    PidBetheBlochErrorDch_AVT::initialize(ComPathNameSearch("PidDchSvtDrcCalib/dEdx_errors_DCH_R22.txt"));
    PidBetheBlochErrorSvt_AVT::initialize(ComPathNameSearch("PidDchSvtDrcCalib/dEdx_errors_SVT_R22.txt"));
    PidBetheBlochPhiCorrDch_AVT::initialize(ComPathNameSearch("PidDchSvtDrcCalib/dEdx_phiCorrections_DCH_R22.txt"));
    PidBetheBlochPhiCorrSvt_AVT::initialize(ComPathNameSearch("PidDchSvtDrcCalib/dEdx_phiCorrections_SVT_R22.txt"));
    PidExpectedBetheBlochDch_AVT::initialize(ComPathNameSearch("PidDchSvtDrcCalib/dEdx_parameters_DCH_R22.txt"));
    PidExpectedBetheBlochSvt_AVT::initialize(ComPathNameSearch("PidDchSvtDrcCalib/dEdx_parameters_SVT_R22.txt"));
  }
  */

  // NA
  //UsedForBdtTraining::initialize();
}


//////////////////////////////////////////////////////////////////////////
void
PacPidCand2Ntuple::setEvent(AbsEvent *ev,AbsParmVector<string> *pidVector,
       AbsParmVector<string> *tagVector)
{

  _theAssoc = Ifd< BtaMcAssoc >::get( ev, _mcTruthList.value());
  // For MC Truth
  if (getParmValue("writeMcTruth")) {
    if (_theAssoc==NULL) {
      ErrMsg(fatal) << " *** ERROR in PacPidCand2Ntuple::setEvent : " << "\n"
        << "     MCTruth-Map \"" <<  _mcTruthList.value() << "\" not in event ! " << endmsg;
    }
    assert(_theAssoc);
  }

  _currentEvent = ev;
  HepAList<BtaCandidate> *list;
  const vector<string> &keyList1 = pidVector->value();
  const vector<string> &keyList2 = tagVector->value();
  _nSelectors = keyList1.size();
  for (int i=0; i<keyList1.size(); i++) {
    _namesOfPidLists[i] = keyList1[i].c_str();
    IfdStrKey theKey(_namesOfPidLists[i].c_str());
    list = Ifd<HepAList<BtaCandidate> >::get(ev,theKey);
    _listOfPidLists[i] = list;
    // cout << " List : " << _namesOfPidLists[i]
    //  << " at " << _listOfPidLists[i] << endl;
  }
  _namesOfTagBits.clear();
  for (int k=0; k<keyList2.size(); k++) {
    _namesOfTagBits.push_back(string(keyList2[k].c_str()));
  }
}


////////////////////////////////////////////////////////////////////////////
// Carl Vuosalo, August 2007:
//
// makeBDTCutCol -- creates the columns that hold the cut values used by muon BDT.
//
// These columns can be used to verify that the muon BDT boolean columns are
// set correctly.  They have names like ifrmuBDTLooseCut and ifrmuBDTLoPLooseCut.
// 
// momentum -- particle momentum (p).
// theta -- angle.
// runNum -- run number of candidate.
// charge -- charge of candidate.
// str -- the name of the list, muBDTLoose or muBDTLoPLoose, for example.  
//        The list names must have this format.
// nt -- the ntuple we're adding to.

void PacPidCand2Ntuple::makeBDTCutCol(double momentum, double theta,
				   double runNum, double charge, 
				   const HepString &str, HepTuple *nt)
{
  if (getParmValue("writeIfrQual") == false || str.index("BDT") < 0 || str.index("BDTKaon") >= 0)
    return;
  unsigned int critLen = str.length();
  unsigned int trimLen = 5;
  if (str(7) == 'P')  // Looking for the "P" in "LoP".
    trimLen = 8;
// Trim off leading "muBDT" or "muBDTLoP".
  HepString endCrit = str(trimLen, critLen - trimLen);
  
  // NA
  //double bdtCut;
  //if (str.index("LoP") >= 0)
  //bdtCut = _muonLoPTree.cuts.get(momentum, theta, runNum, charge, endCrit);
  //else
  //bdtCut = _muonTree.cuts.get(momentum, theta, runNum, charge, endCrit);
  //string fldNam = string("ifr") + (const char *) str + string("Cut");
  //nt->column(fldNam, bdtCut);
}


//////////////////////////////////////////////////
void
PacPidCand2Ntuple::convert(const BtaCandidate *c, AbsEvent *ev, HepTuple *nt)
{

  // avtelnov, 2007/12/06: printDebuggingInfo() is to be commented out only in private test releases
  // printDebuggingInfo(c);

  // Added at Al Eisner's request -- Sasha, June 2007
  // Placing this either in the constructor or beginJob does not work 
  // because Pdt initialization has to occur prior to calling _MergedPi0Workhorse->setUp()
  // NA
  /*
  if (!_MergedPi0Workhorse) {
  // get the factory from the environment the first time it's needed
    EmcIDFactory* factory=  EmcIDFactory::instance();
    if (!factory) {
      ErrMsg(fatal) << "PidMergedPi0MicroSelector"
                    << ": EmcIDFactory not created.  "
                    << "Check that EmcLoadPid has been run in EmcInitSequence."
                    << endmsg;
    }
    
    _MergedPi0Workhorse=
      dynamic_cast<EmcMergedPi0Identifier*>(factory->getSpecificIdentifier("MergedPi0Identifier"));
    
    if ( !_MergedPi0Workhorse ) {
      ErrMsg(fatal) << "Workhorse for MergedPi0Identifier couldn't be found."
                    << endl
                    << "Configuration of EmcIDFactory has changed."
                    << endmsg;
    }
    
    _MergedPi0Workhorse->setUp();
    _MergedPi0algo = new BtaMergedPi0Algo(PidSystem::emc);
  } // end setting up pi0 consistency stuff
  */

  static const IfdStrKey keyAbsEventID("AbsEventID");
  const AbsEventID* eventID = Ifd<AbsEventID>::get( ev , keyAbsEventID );
  static const IfdStrKey keyDflt("Default");
  HepAList< EventInfo >* infoList= Ifd<HepAList<EventInfo> >::get(_currentEvent, keyDflt);
  const EventInfo* eventInfo(infoList->first());

  int nrun(0);
  int date(0);
  int mday(0); // added for Run 7 - day of the month, 1..31
  int upperID(0), lowerID(0); // these are actually unsigned ints, but this does not matter

  if (eventID != 0) {
    nrun = eventID->run();
    nt->column("runNumber",(int)eventID->run(),0,"EventID");
    odfTime timeStamp = eventID->eventIdTriplet().timeStamp();
    //    nt->column("platform",(int)eventID->eventIdTriplet().platform(),0,"EventID");


    // NA
    /*
    // deal with LER & HER trickle
    const L1FctTimePoint* tp = 0;

    unsigned m1bin = 0;
    unsigned m2bin = 0;

    // LER first
    if (_getTrickleInfo.value())  {
      if ( _trickleRegion_LER == 0 ) {
  _trickleRegion_LER = new OprTrickleRegions(  _timePointLER.valueString().c_str(), gblEnv );
  ostream& ooo=  ErrMsg(warning);
  ooo << _timePointLER.valueString() << " Define Trickle Regions";
  _trickleRegion_LER->printLimits(ooo);
  _trickleRegion_LER->printPhysicsLimits(ooo);
  ooo << endmsg;
      }

      tp = _timePointInspector_LER.get( ev );
    }

    if  ( tp != 0 ) {
      //    int phase        = _trickleRegion_LER->phase(timeStamp - tp->timeStamp());
      //    int nrevolutions = _trickleRegion_LER->revolutions(timeStamp - tp->timeStamp());
      m1bin = _trickleRegion_LER->region(timeStamp - tp->timeStamp());
      m2bin = _trickleRegion_LER->regionPhysics(timeStamp - tp->timeStamp());
    }
    // platform is only 1 byte, TrickleRegion is a 4-bit enum, so stick them all together
    unsigned platform= eventID->eventIdTriplet().platform();
    platform= platform | ( m1bin << 8 ) | ( m2bin << 16 );
    
  // there is still space to deal with HER trickle
    if (_getTrickleInfo.value())  {
      if ( _trickleRegion_HER == 0 ) {
	_trickleRegion_HER = new OprTrickleRegions(  _timePointHER.valueString().c_str(), gblEnv );
	ostream& ooo=  ErrMsg(warning);
	ooo << _timePointHER.valueString() << " Define Trickle Regions";
	_trickleRegion_HER->printLimits(ooo);
	_trickleRegion_HER->printPhysicsLimits(ooo);
	ooo << endmsg;
      }
      tp = _timePointInspector_HER.get( ev );
    }
    
    if  ( tp != 0 ) {
      m1bin = _trickleRegion_HER->region(timeStamp - tp->timeStamp());
      m2bin = _trickleRegion_HER->regionPhysics(timeStamp - tp->timeStamp());
      platform= platform | ( m1bin << 12 ) | ( m2bin << 20 );
    }
    */

  // now find out year and month
    const EidCondKeyTriplet* triplet = gblEnv->getGen()->primaryCondKey();
    if (triplet==NULL) {
      ErrMsg(debugging) << " No Event time ! "<< endmsg;
    } else {
      BdbTime bTime = triplet->key();
      struct tm myTime;
      bTime.tm(&myTime,BdbTime::Local);
      int y = myTime.tm_year + 1900;
      if (_theAssoc) { y+=25; }
      int m = myTime.tm_mon + 1;
      mday = myTime.tm_mday;
      date = y*100 + m;
    }
    nt->column("date",date,0,"EventID"); // YYYYMM
    nt->column("mday",mday,0,"EventID"); // DD - day of the month, added for Run 7
    
    // NA
    //nt->column("platform",(int)platform,0,"EventID");
    
    nt->column("partition",(int)eventID->eventIdTriplet().partitionMask() ,0,"EventID");
    upperID = (int)timeStamp.binary().upper;
    lowerID = (int)timeStamp.binary().lower;
    nt->column("upperID",(int)upperID,0,"EventID");
    nt->column("lowerID",(int)lowerID,0,"EventID");
  } else {
    ErrMsg(fatal) << " Event information not available! " << endmsg;
    nt->column("runNumber",-1,0,"EventID");
    nt->column("date",0,0,"EventID");
    nt->column("mday",0,0,"EventID");
    nt->column("platform",-1,0,"EventID");
    nt->column("partition",-1,0,"EventID");
    nt->column("upperID",-1,0,"EventID");
    nt->column("lowerID",-1,0,"EventID");
  }

// background key

  AbsEventID *eventID2(0);
  static const IfdStrKey keyAbsEventIDList("AbsEventIDList");
  HepAList<AbsEventID> *idList = Ifd<HepAList<AbsEventID> >::get(ev, keyAbsEventIDList);
  if ( 0 != idList && idList->length() > 1 ) {
    eventID2 = (*idList)[1];
    ErrMsg(trace) << " switching to background event for timestamp: "
                  << eventID->eventIdTriplet() << endmsg;
//    cout << " switching to background event for timestamp: "
//                  << eventID->eventIdTriplet() << "\n";
  }

  if (eventID2 != 0) {
    odfTime timeStamp2 = eventID2->eventIdTriplet().timeStamp();
    nt->column("bkgdUpperID",(int)timeStamp2.binary().upper,0,"EventID");
    nt->column("bkgdLowerID",(int)timeStamp2.binary().lower,0,"EventID");
  } else {
    nt->column("bkgdUpperID",-1,0,"EventID");
    nt->column("bkgdLowerID",-1,0,"EventID");
  }



  //
  // Tag bits
  //
  //

  AbsEventTag* tag = Ifd<AbsEventTag>::get( ev ) ;
  if (getParmValue("writeTagBits")) {
    for (int i=0; i<_namesOfTagBits.size(); i++) {
      string tagBit = _namesOfTagBits[i];
      bool b(false);
      if ( !(tag->getBool( b,tagBit ))) {
	// NA
        //ErrMsg(error) << " *** ERROR : Tag Bit # " << i << " with name "  << tagBit << " not in Event ! " << endmsg;
        b = false;
      }
      nt->column(tagBit,b);
    }
  }


// write some particular tag quantities regardless of value of "writeTagBits"

// SASHA, 2008/01/20: TEMPORARILY TURN OFF SOME VARIABLES THAT WE DO NOT USE
/*

    float thrustMag(-100.0);
    float thrustMagAll(-100.0);
    float thrustCosTh(-100.0);
    float thrustCosThAll(-100.0);
    float thrustPhi(-100.0);
    float thrustPhiAll(-100.0);
    float R2(-100.0);
    float R2All(-100.0);
    float sphericityAll(-100.0);
    int nTracks(-1);
    int nGoodTrkLoose(-1);
    int nGoodTrkTight(-1);

    if (tag) {
      tag->getFloat(thrustMag,"thrustMag");
      tag->getFloat(thrustMagAll,"thrustMagAll");
      tag->getFloat(thrustCosTh,"thrustCosTh");
      tag->getFloat(thrustCosThAll,"thrustCosThAll");
      tag->getFloat(thrustPhi,"thrustPhi");
      tag->getFloat(thrustPhiAll,"thrustPhiAll");
      tag->getFloat(R2,"R2");
      tag->getFloat(R2All,"R2All");
      tag->getFloat(sphericityAll,"sphericityAll");
      tag->getInt(nTracks,"nTracks");
      tag->getInt(nGoodTrkLoose,"nGoodTrkLoose");
      tag->getInt(nGoodTrkTight,"nGoodTrkTight");
      
// packing for thrustMag and thrustMagAll is limited by very sharp peaking toward 1 for Bhabhas
      nt->column("thrustMag",packFloat(thrustMag,8)); 
      nt->column("thrustMagAll",packFloat(thrustMagAll,8));
// thrustCosTh RMS is 0.07 near the peak (at 0.7) for Bhabhas. 
// 8.0e-4 * 0.7 would be added to that in quadrature from 15-bit truncation (0.41/2^(24-bits_to_kill))
// which worsens the thrustCosTh resolution by less than 1 part in 10,000. 
// Similar logic applied to the other variables: RMS from packing should be less than 0.01
// of the smallest RMS characteristic of the distributions of the particular variable
// for all PID samples. Unless, of course, it may be used in a more complex way, e.g., 
// for comparison with MC Truth or to calculate invariant mass, etc. - then we keep more precision. 
      nt->column("thrustCosTh",packFloat(thrustCosTh,14)); 
      nt->column("thrustCosThAll",packFloat(thrustCosThAll,14));
      nt->column("thrustPhi",packFloatPhi(thrustPhi,16)); // who cares if it is wrong by at most 0.7 degrees
      nt->column("thrustPhiAll",packFloatPhi(thrustPhiAll,16));
      nt->column("R2",packFloat(R2,7)); // could have packed more if not for two-prongs
      nt->column("R2All",packFloat(R2All,7));
      nt->column("sphericityAll",packFloat(sphericityAll,14)); // broad; for Bhabhas, sharply peaks at 0.
 
   } else { // do not pack "failed" values
      nt->column("thrustMag",thrustMag);
      nt->column("thrustMagAll",thrustMagAll);
      nt->column("thrustCosTh",thrustCosTh);
      nt->column("thrustCosThAll",thrustCosThAll);
      nt->column("thrustPhi",thrustPhi);
      nt->column("thrustPhiAll",thrustPhiAll);
      nt->column("R2",R2); 
      nt->column("R2All",R2All);
      nt->column("sphericityAll",sphericityAll);
    }    

    nt->column("nTracks",nTracks);
    nt->column("nGoodTrkLoose",nGoodTrkLoose);
    nt->column("nGoodTrkTight",nGoodTrkTight);

*/ // END SASHA 2008/01/20

  //
  // MCTruth
  //
    if (_theAssoc!=NULL) {
      BtaCandidate *tc = _theAssoc->mcFromReco(c);
      if (tc!=0) {
	nt->column("lundIdTruth",tc->pdtEntry()->lundId());
	if (tc->theMother()!=NULL) {
	  nt->column("motherIdTruth",tc->theMother()->pdtEntry()->lundId());
	  if ((tc->theMother())->theMother()!=NULL) {
	    nt->column("grandmotherIdTruth",(tc->theMother())->theMother()->pdtEntry()->lundId());
	  } else {
	    nt->column("grandmotherIdTruth",0);
	  }
	} else {
	  nt->column("motherIdTruth",0);
	  nt->column("grandmotherIdTruth",0);
	}
      //  nt->column("energyTruth",(float) tc->energy());
      //  nt->column("pTruth",(float) tc->p());
      //  nt->column("thetaTruth",(float) tc->p3().theta());
      //  nt->column("phiTruth",(float) tc->p3().phi());
      //  nt->column("chargeTruth",(float) tc->charge());
	nt->column("energyTruth",packFloat((float)tc->energy(),6)); // for reco-truth comparisons, keep much of precision
	nt->column("pTruth",packFloat((float)tc->p(),6));
	nt->column("thetaTruth",packFloat((float)tc->p3().theta(),6));
	nt->column("phiTruth",packFloat((float)tc->p3().phi(),6));
	nt->column("chargeTruth",(float)tc->charge()); // no point in packing a float that stores an integer!
      } else {
	nt->column("lundIdTruth",0);
	nt->column("motherIdTruth",0);
	nt->column("grandmotherIdTruth",0);
	nt->column("energyTruth",(float) 0.0);
	nt->column("pTruth",(float) 0.0);
	nt->column("thetaTruth",(float) 0.0);
	nt->column("phiTruth",(float) 0.0);
	nt->column("chargeTruth",(float) 0.0);
      }
    }


  // kinematics (for the hypo this track comes with)

  const int ncharge = (int)c->charge();
  const float candMomentum = c->p();
  const float candTheta = c->p3().theta(); // in dE/dx calculations, we will be using the one from the pion hypothesis 
  const float candPhi = c->p3().phi();

// Momentum at IP dependent on particle type. We do not know what Pdt type c comes with, so we record all 5
// Charge, actually, also may depend on the particle-type hypothesis used in the Kalman fit. The charge 
// discrepancy occurs particularly often in protons. Record the hypo-dependent charge as well
  BtaCandidate tempCand(*c);

  tempCand.setType(Pdt::lookup(PdtPid::kaon, ncharge));
  const float pIPk = tempCand.p();
  const int charge_k = static_cast<int> (tempCand.charge());

  tempCand.setType(Pdt::lookup(PdtPid::proton, ncharge));
  const float pIPp = tempCand.p();
  const int charge_p = static_cast<int> (tempCand.charge());

  tempCand.setType(Pdt::lookup(PdtPid::electron, ncharge));
  const float pIPe = tempCand.p();
  const int charge_e = static_cast<int> (tempCand.charge());

  tempCand.setType(Pdt::lookup(PdtPid::muon, ncharge));
  const float pIPmu = tempCand.p();
  const int charge_mu = static_cast<int> (tempCand.charge());

  tempCand.setType(Pdt::lookup(PdtPid::pion, ncharge));
  const float pIPpi = tempCand.p();
  const float candThetaPi = tempCand.p3().theta();
  
  // NA: not used in uncommented code
  //const float candPhiPi = tempCand.p3().phi();
  
  const int charge_pi = static_cast<int> (tempCand.charge());
  float pdchPi = -1.0;
  const BtaPidQual* pQualPi = tempCand.getMicroAdapter()->getPidQual();
  if (pQualPi) pdchPi = pIPpi + pQualPi->deltaDchMomentum();
 
  nt->column("energy",packFloat((float) c->energy(),6));
  nt->column("p",packFloat(candMomentum,6)); // for the hypo this track comes with
// SASHA, 2008/01/20: TURN OFF SOME VARIABLES THAT WE DO NOT USE
  nt->column("pIPe",packFloat(pIPe,6));
  nt->column("pIPmu",packFloat(pIPmu,6));
  nt->column("pIPpi",packFloat(pIPpi,6));
  nt->column("pIPk",packFloat(pIPk,6));
  nt->column("pIPp",packFloat(pIPp,6));
  nt->column("theta",packFloat(candTheta,6));
  nt->column("thetaPi",packFloat(candThetaPi,6));
  nt->column("phi",packFloat(candPhi,6));
  nt->column("charge",(float) c->charge()); // for the hypo this track comes with
  nt->column("charge_e",charge_e);
  nt->column("charge_mu",charge_mu);
  nt->column("charge_pi",charge_pi);
  nt->column("charge_k",charge_k);
  nt->column("charge_p",charge_p);

// Sasha Telnov, 2007/08/10 -- 
// this is the lundId of the candidate that was passed to the "convert" method.
// Depending on the version of BetaPidCalibNtuple, it may or may not match 
// the intended particle type. In particular, in official PID ntuples up to and 
// including r22a (summer 2007), the following samples had lundId = "pion":
// electrons (ntp201 and ntp202), muons (204 and 205, but not 207),
// dE/dx slow kaons (ntp803) and protons (ntp804), Jpsi leptons (901, 902).
// This has now been fixed
  nt->column("lundId",c->pdtEntry()->lundId()); 

  if (c->getMicroAdapter()==NULL) { return; }
  const BtaTrkQual* tQual = c->getMicroAdapter()->getTrkQual();
  const BtaCalQual* cQual = c->getMicroAdapter()->getCalQual();
  
  // NA
  //const BtaIfrQual* iQual = c->getMicroAdapter()->getIfrQual();
  const BtaIfrQual* iQual = NULL;
  
  const BtaPidQual* pQual = c->getMicroAdapter()->getPidQual();
  const BtaPidInfo* pInfo = c->getMicroAdapter()->getPidInfo();
  if (tQual==0&&cQual==0&&iQual==0&&pQual==0&&pInfo==0) ErrMsg(error) << "BAD Micro cand " << endmsg;

  if (c->recoObject()!=NULL) {
    const BbrLorentzVectorErr & ip(gblEnv->getBta()->pepBeams()->interactionPoint());
    const HepPoint beamSpotPoint(ip.x(),ip.y(),ip.z());
    HepPoint closest(c->recoObject()->position(beamSpotPoint,BtaAbsRecoObject::XY));
    const double docax(closest.x()-beamSpotPoint.x());
    const double docay(closest.y()-beamSpotPoint.y());
    const double docaz(closest.z()-beamSpotPoint.z());
    nt->column("xpoca",packFloat((float) docax,16)); // peaks around zero, so effect of packing is negligible
    nt->column("ypoca",packFloat((float) docay,16));
    nt->column("zpoca",packFloat((float) docaz,16));
  } else {
    nt->column("xpoca",(float) -99.0);
    nt->column("ypoca",(float) -99.0);
    nt->column("zpoca",(float) -99.0);
  }

  // kinematics for CMS (boosted under electron hypo)
  BtaCandidate c2(*c);
  c2.setType(Pdt::lookup(PdtPid::electron,(int) c->charge()));
  BtaCandidate Y4S(eventInfo->cmFrame());
  BtaBooster theBooster(&Y4S);
  BtaCandidate* theBoostedCand = new BtaCandidate(theBooster.boostTo(c2));
  nt->column("electronpcms",packFloat((float) theBoostedCand->p(),14)); // this is already sufficiently imprecise
  nt->column("electronthetacms",packFloat((float) theBoostedCand->p3().theta(),14));
//nt->column("electronphicms",packFloat((float) theBoostedCand->p3().phi(),14));
  delete theBoostedCand;

  // EMC

// add EMC longitudinal shower depth

  double showerDepth(0);
  double showerDepthErr(0);
  bool OK = BetaMiniTools::showerDepth(c, showerDepth, showerDepthErr);
  if (OK) {
  // for electrons, mean 10, RMS 6, so adding 8e-4 in quadrature to RMS will not hurt a bit 
    nt->column("emcdepth",packFloat((float)showerDepth,13));
    nt->column("emcdeptherr",packFloat((float)showerDepthErr,13));
  } else {
    nt->column("emcdepth",(float)-10000.0);
    nt->column("emcdeptherr",(float)-10000.0);
  }

// add inputs for EMC NN for muon PID
  double ecalnn(0.0);
  double latcalnn(0.0);
  double z20calnn(0.0);
  double z42calnn(0.0);
  double s1s9calnn(0.0);
  double s9s25calnn(0.0);
  double smomcalnn(0.0);
  double ncrycalnn(0.0);
  double eopcalnn(0.0);
  
  // NA
  //bool doemcnn(false);

  if (getParmValue("writeCalQual")) {

// Added at Al Eisner's request
    float pi0Cons=-1.;
    float pi0LogLike=-1., gammaLogLike=-1., pi0LikeRatio=-999.;
    int pi0ConsStat=-1;
    
    if  (cQual){
      nt->column("eraw",packFloat((float)cQual->rawEnergy(),10));
      nt->column("ecal",packFloat((float)cQual->ecalEnergy(),10));
      nt->column("lmom",packFloat((float)cQual->lateralMoment(),10));
      nt->column("zmom20",packFloat((float)cQual->absZernike20(),8));
      nt->column("zmom42",packFloat((float)cQual->absZernike42(),10));
      nt->column("s1s9",packFloat((float)cQual->s1s9(),8));
      nt->column("s9s25",packFloat((float)cQual->s9s25(),8));
      nt->column("secmom",packFloat((float)cQual->secondMomentTP(),12));
      HepPoint p = cQual->centroid();
      nt->column("phicluster",packFloatPhi((float)p.phi(),12));
      nt->column("thetacluster",packFloat((float)p.theta(),12));
      nt->column("phicmat",packFloat((float)cQual->phiMatchConsistency().consistency(),10));
      nt->column("ncry",cQual->nCrystals());
      nt->column("nbump",cQual->nBumps());
      nt->column("emcstat",cQual->status());

      if ( candMomentum > 0.5 &&
            cQual->ecalEnergy() > 0.0 && cQual->ecalEnergy() < 2.0 &&
            cQual->secondMomentTP() > 0.0 &&
            cQual->absZernike42() > 0.0 &&
            cQual->nCrystals() < 13 ) {
        
	// NA
	//doemcnn    = true;
        
	ecalnn     = cQual->ecalEnergy()/2.000001;
        latcalnn   = cQual->lateralMoment()/1.000001;
        z20calnn   = cQual->absZernike20()*cQual->absZernike20()*cQual->absZernike20()/1.000001;
        z42calnn   = sqrt(cQual->absZernike42())/1.000001;
        s1s9calnn  = cQual->s1s9()/1.000001;
        s9s25calnn = cQual->s9s25()*cQual->s9s25()*cQual->s9s25()/1.000001;
        smomcalnn  = 10.0*sqrt(cQual->secondMomentTP())/1.000001;
        ncrycalnn  = cQual->nCrystals()/12.000001;
        eopcalnn   = (cQual->ecalEnergy()/candMomentum)/2.000001;
        if (candMomentum > 4.0) eopcalnn = (cQual->ecalEnergy()/4.0)/2.000001;
        if (smomcalnn > 1.0) smomcalnn = 0.9999999;
        if (eopcalnn > 1.0) eopcalnn = 0.9999999;
      }
      
    //pi0 Likelihood...
      // NA
      /*
      const PdtEntry* gammaEntry;
      gammaEntry = Pdt::lookup((PdtPid::PidNeutralType)PdtPid::gamma,0);
      const PdtEntry* pi0Entry;
      pi0Entry = Pdt::lookup((PdtPid::PidNeutralType)PdtPid::pi0,0);
      Consistency *gammaConsistency = _MergedPi0algo->consistency(*c,gammaEntry);
      Consistency *pi0Consistency = _MergedPi0algo->consistency(*c,pi0Entry);
      if (pi0Consistency!=NULL){
	pi0Cons = pi0Consistency->consistency();
	double pi0LL = pi0Consistency->likelihood()>0 ? log10(pi0Consistency->likelihood()) : 0;
	double gammaLL = gammaConsistency->likelihood()>0 ? log10(gammaConsistency->likelihood()) : 0;
	pi0LogLike = pi0LL;
	gammaLogLike = gammaLL;
	pi0LikeRatio = pi0LL - gammaLL;
	pi0ConsStat = pi0Consistency->status();
      } 
      */

    /////
    } else { // CalQual does not exist
      if (pQual==0) {
	fillEmptyCalQual(nt);
      } else {
	if (pQual->thetaAtEMC()==0.0) {
	  fillEmptyCalQual(nt);
	} else {
	  HepAList<BtaCandidate> *chargedList;
	  HepAList<BtaCandidate> *neutralList;
	  static const IfdStrKey keyCT("ChargedTracks");
	  static const IfdStrKey keyCN("CalorNeutral");
	  getTmpAList(ev,chargedList,keyCT);
	  getTmpAList(ev,neutralList,keyCN);
	  
	  HepAListIterator<BtaCandidate> iterC(*chargedList);
	  HepAListIterator<BtaCandidate> iterN(*neutralList);
	  BtaCandidate *cand;
	  bool charged(false);
	  double shortestDistance=0.0;
	  double mytheta = pQual->thetaAtEMC();
	  double myphi = pQual->phiAtEMC();
	  const BtaCalQual *bestCal = NULL;
	  while (cand=iterC()) {
	    const  BtaCalQual *cal = cand->getMicroAdapter()->getCalQual();
	    if (cal!=0) {
	      double th = cal->centroid().theta();
	      double ph = cal->centroid().phi();
	      double dist = sqrt((mytheta-th)*(mytheta-th)+(myphi-ph)*(myphi-ph));
	      if ((dist<shortestDistance)  || (bestCal==NULL)) {
		shortestDistance=dist;
		bestCal = cal;
		charged=true;
	      }
	    }
	  }
	  while (cand=iterN()) {
	    const BtaCalQual *cal = cand->getMicroAdapter()->getCalQual();
	    if (cal!=0) {
	      double th = cal->centroid().theta();
	      double ph = cal->centroid().phi();
	      double dist = sqrt((mytheta-th)*(mytheta-th)+(myphi-ph)*(myphi-ph));
	      if ((dist<shortestDistance)  || (bestCal==NULL)) {
		shortestDistance=dist;
		bestCal = cal;
		charged=false;
	      }
	    }
	  }
	  if (bestCal==NULL) {
	    fillEmptyCalQual(nt);
	  } else {
	    nt->column("eraw",packFloat((float)-bestCal->rawEnergy(),10));
	    nt->column("ecal",packFloat((float)-bestCal->ecalEnergy(),10));
	    nt->column("lmom",packFloat((float)-bestCal->lateralMoment(),10));
	    nt->column("zmom20",packFloat((float)-bestCal->absZernike20(),8));
	    nt->column("zmom42",packFloat((float)-bestCal->absZernike42(),10));
	    nt->column("s1s9",packFloat((float)-bestCal->s1s9(),8));
	    nt->column("s9s25",packFloat((float)-bestCal->s9s25(),8));
	    nt->column("secmom",packFloat((float)-bestCal->secondMomentTP(),12));
	    HepPoint p = bestCal->centroid();
	    nt->column("phicluster",packFloatPhi((float)p.phi(),12));
	    nt->column("thetacluster",packFloat((float)p.theta(),12));
	    if (charged) {  nt->column("phicmat",(float)-1.0); } else {  nt->column("phicmat",(float)-2.0); }
	    nt->column("ncry",-bestCal->nCrystals());
	    nt->column("nbump",-bestCal->nBumps());
	    nt->column("emcstat",bestCal->status());
	  }
	}
      }
    }

    nt->column("pi0Cons",packFloat(pi0Cons,10));
    nt->column("pi0LogLike",packFloat(pi0LogLike,10));
    nt->column("gammaLogLike",packFloat(gammaLogLike,10));
    nt->column("pi0LikeRatio",packFloat(pi0LikeRatio,8));
    nt->column("pi0ConsStat",pi0ConsStat);
  }
  
// The pi0 stuff...
  







////////////////////
  // TRKQUAL
// 2007/07/24, avtelnov: added SVT hit pattern
//   The 10 bits of the word correspond to whether one of the 10 views of the SVT registered a hit for this candidate:
//     Layer 5 z-view = Most Significant Bit 
//     Layer 4 z-view 
//     Layer 3 z-view 
//     Layer 2 z-view 
//     Layer 1 z-view 
//     Layer 5 r/phi-view 
//     Layer 4 r/phi-view 
//     Layer 3 r/phi-view 
//     Layer 2 r/phi-view 
//     Layer 1 r-phi-view = Least Significant Bit
// You can use a bit-wise '&' operator to check the bit values

  if (getParmValue("writeTrkQual")) {
    const BtaTrkQual* tQual2 = c->getMicroAdapter()->getTrkQual();
    if (tQual2) {
      nt->column("ndch",(int) tQual2->nDchHits());
      nt->column("nsvt",(int) tQual2->nSvtHits());
    //nt->column("isvt",(int) tQual2->nSvtHits()); why duplicate stuff?
      nt->column("svtpatt",(int) tQual2->SvtPattern()); 
      nt->column("fhit",(int) tQual2->firstDchHit());
      nt->column("lhit",(int) tQual2->lastDchHit());
    } else { // this never happens to charged tracks!
      nt->column("ndch",(int) -1);
      nt->column("nsvt",(int) -1);
    // nt->column("isvt",(int) -1);
      nt->column("svtpatt",(int) -1); 
      nt->column("fhit",(int) -1);
      nt->column("lhit",(int) -1);
    }
  }

  // PIDQUAL

  float dedxdch(0);
  float dedxsvt(0);
  int  nsampdch(0);
  int  nsampsvt(0);
  float pdch = candMomentum;
  float deltaDchMom = 0;

  if (getParmValue("writePidQual")) {
    if (pQual) {
      dedxdch = pQual->dEdXDch();
      dedxsvt = pQual->dEdXSvt();
      nsampdch = pQual->nSamplesDeDxDch();
      nsampsvt = pQual->nSamplesDeDxSvt();
      deltaDchMom = pQual->deltaDchMomentum();
      pdch = candMomentum + deltaDchMom;
      if (pdchPi == -1.0) pdchPi = pdch; 
      nt->column("dedxdch",packFloat(dedxdch,10)); // 10% resolution, so adding 2.5e-5 in quadrature does nothing
      nt->column("dedxsvt",packFloat(dedxsvt,12)); // 15% resolution
      nt->column("pdch",packFloat(pdch,8)); // hypothesis-independent - in R24, too!
      nt->column("deltaDchMom",packFloat(deltaDchMom,12));      
      nt->column("nsampdch",nsampdch);
      nt->column("nsampsvt",nsampsvt);
      nt->column("nphot",pQual->ringNPhot());
      nt->column("nbgphot",pQual->ringNBkgd());
      nt->column("cthe",packFloat((float)pQual->thetaC(),9)); // 0.8 with 0.3% resolution, add 1.2e-5 in quadrature
      nt->column("sthe",packFloat((float)pQual->thetaCErr(),13)); // error on the error ~ the error itself, hence crude packing
      nt->column("nphote",pQual->ringNExPhot(Pdt::lookup(PdtPid::electron)));
      nt->column("nphotmu",pQual->ringNExPhot(Pdt::lookup(PdtPid::muon)));
      nt->column("nphotpi",pQual->ringNExPhot(Pdt::lookup(PdtPid::pion)));
      nt->column("nphotk",pQual->ringNExPhot(Pdt::lookup(PdtPid::kaon)));
      nt->column("nphotp",pQual->ringNExPhot(Pdt::lookup(PdtPid::proton)));
      addNphotFromGlbLikelihood(c,nt);
      nt->column("pdrc",packFloat((float) (candMomentum + pQual->deltaDrcMomentum()),9));
      nt->column("drcinbar", (int) pQual->drcInBar());
      nt->column("drcoutbar", (int) pQual->drcExitBar());
      nt->column("drcxpos",(float)  pQual->drcXPos()); // looks like this is already very tightly packed in CM2
      nt->column("phiatemc",packFloatPhi((float)pQual->phiAtEMC(),13));
      nt->column("thetaatemc",packFloat((float)pQual->thetaAtEMC(),13));

      // NA, forward pid stuff
      nt->column( "fwdpidmeas",
		  packFloat( (float)pQual->forwardPidMeasurement(), 13 ) );
      nt->column( "fwdpiderr",
		  packFloat( (float)pQual->forwardPidError(), 13 ) );

    } else {
      nt->column("dedxdch",(float)-1.0);
      nt->column("dedxsvt",(float)-1.0);
      nt->column("pdch",(float)-1.0);
      nt->column("nsampdch",-1);
      nt->column("nsampsvt",-1);
      nt->column("fhit",-1);
      nt->column("lhit",-1);
      nt->column("nphot",-1);
      nt->column("nbgphot",-1);
      nt->column("cthe",(float)-1.0);
      nt->column("sthe",(float)-1.0);
      nt->column("nphote",-1);
      nt->column("nphotmu",-1);
      nt->column("nphotpi",-1);
      nt->column("nphotk",-1);
      nt->column("nphotp",-1);
      nt->column("pdrc",(float)-1.0);
      nt->column("phiatemc",(float)-1.0);
      nt->column("thetaatemc",(float)-1.0);
      nt->column("drcinbar", (int) 0);
      nt->column("drcoutbar", (int) 0);
      nt->column("drcxpos",(float) 0.0);
      nt->column("phiatemc",(float) 0.0);
      nt->column("thetaatemc",(float) 0.0);

      // NA, forward pid stuff
      nt->column( "fwdpidmeas", (float)0.0 );
      nt->column( "fwdpiderr",  (float)0.0 );
    }
  }

  // IFR
  tempCand.setType(Pdt::lookup(PdtPid::pion, ncharge)); // same as in the PID sequence
  if (getParmValue("writeIfrQual")) {
    
  // Write out the values of the BDT tree
  //PidMuonSPR::MuList muVars(c, _theSms, _muonTree._useDchLH);

    // NA
    //PidMuonSPR::MuList muVars(&tempCand, _theSms, _muonTree._useDchLH);
    //double ifrBDTVal = _muonTree.getClsfrOutput(pIPpi, candThetaPi, muVars);
    //double ifrBDTLoPVal = _muonLoPTree.getClsfrOutput(pIPpi, candThetaPi, muVars);
    //nt->column("ifrBDTVal", ifrBDTVal); 
    //nt->column("ifrBDTLoPVal", ifrBDTLoPVal);
    //nt->column("ifrcrackphi", PidMuonSPR::MuList::getIfrcrackphi(candMomentum, candTheta, candPhi, ncharge));
		
    float logMuLike(-1000.);
    float logPiLike(+1000.);
    if (iQual) {
      // IFR likelihood code
      IfrMicroPidInfo pid( c );

      // NA
      //const PdtEntry * mu =  Pdt::lookup( PdtPid::muon, (int) c->charge() );
      //const PdtEntry * pi =  Pdt::lookup( PdtPid::pion, (int) c->charge() );
      //const Consistency& muCons = pid.consistency( mu );
      //const Consistency& piCons = pid.consistency( pi );
      //if ( muCons.likelihood() > 0)
      //logMuLike = log10( muCons.likelihood() );
      //if ( piCons.likelihood() > 0)
      //logPiLike = log10( piCons.likelihood() );

      //
      double continuity = 0.;
      double ifrLay =  iQual->IfrLayHits();
      double ifrFirst = iQual->firstHit();
      double ifrLast = iQual->lastHit();
      if(ifrLay>1) {
  if(iQual->hasInner()){
    continuity = double(ifrLay) / double(ifrLast-ifrFirst);
  }
  else{
    continuity = double(ifrLay) / double(ifrLast-ifrFirst+1);
  }
      }

      //Lange - Sept 25, 2001 -add sigma multiplicity variable

      double sigmaMulti=50.;

      if( ifrLay > 1 ) {

  int sum = 0;
  int sum2 = 0;
  for(int j=0;j<20;j++)
    {
      const int numStrips = std::min(14,iQual->nStrips(j));
      sum  += numStrips;
      sum2 += (numStrips*numStrips);
    }

  const double normalization = (double)(ifrLay-1);
  const double sumMean = ((double)(sum*sum))/ ((double)ifrLay);
  sigmaMulti =
    sqrt( ( (double)sum2 - sumMean ) / normalization);
      }


      // nn code
      double ifrNNVal=-1.;

      IfrNNKernelSetup *ifrNNSetup;
      if ( _ifrNNKernel )
  ifrNNSetup = _ifrNNKernel;
      else
  ifrNNSetup = (IfrNNKernelSetup*)gblEnv->getBta()->ifrNNKernelSetup();

      if ( ifrNNSetup && cQual) {
  IfrNNKernel *selKernel = ifrNNSetup->SelectKernel(candMomentum, c->p3().theta());
  vector<double> ifrNNVars;
  ifrNNVars.push_back((iQual->expectedInteractionLengths()-iQual->measuredInteractionLengths())/5.0);
  ifrNNVars.push_back(cQual->ecalEnergy());
  ifrNNVars.push_back(iQual->IfrTrkMatchChi2()/7.0); // changed to add correct normalization
  ifrNNVars.push_back(sigmaMulti/5.0);               // changed to add correct normalization
  ifrNNVars.push_back(iQual->clusterFitChi2()/20.0); // changed to add correct normalization
  ifrNNVars.push_back(continuity);
  ifrNNVars.push_back(double(iQual->IfrNStrips())/double(iQual->IfrLayHits())/10.0);
  ifrNNVars.push_back(iQual->measuredInteractionLengths()/5.0);

  if ( selKernel ) {
    selKernel->readInVar(ifrNNVars);
    selKernel->ComputeWT(0);
    ifrNNVal=selKernel->EvtOutputVal();
  }
      }
      nt->column("logmulike",packFloat(logMuLike,11)); // RMS ~30%, +/-1000 are "magic values", so keep at least 9 bits in the significand
      nt->column("logpilike",packFloat(logPiLike,11));
      nt->column("ifrns",iQual->IfrNStrips());
      nt->column("ifrexpintlen",packFloat((float)iQual->expectedInteractionLengths(),11)); // these are pretty crude variables
      nt->column("ifrmeasintlen",packFloat((float)iQual->measuredInteractionLengths(),11)); // what matters is the difference
      nt->column("ifrintlenbeforeiron",packFloat((float)iQual->interactionLengthsBeforeIron(),13)); // 30% spread
      nt->column("ifrmatchchi2",packFloat((float)iQual->IfrTrkMatchChi2(),12)); // resolution ~ 100%
      nt->column("ifremcmatch",(float)iQual->IfrEmcMatch()); // a float that is, in fact, very discrete 
      nt->column("ifrfitchi2",packFloat((float)iQual->clusterFitChi2(),12)); // resolution ~ 100%
      nt->column("ifrlayhits",iQual->IfrLayHits());
      nt->column("ifrfirsthit",iQual->firstHit());
      nt->column("ifrlasthit",iQual->lastHit());
      nt->column("ifrhasinner",iQual->hasInner());
      nt->column("ifrhasbarrel",iQual->hasBarrel());
      nt->column("ifrhasfwd",iQual->hasFWD());
      nt->column("ifrhasbwd",iQual->hasBWD());
      nt->column("ifrcont",(float)continuity);  // a float that is, in fact, quite discrete
      nt->column("ifrsigmu",(float)sigmaMulti); // resolution ~ 100%
      nt->column("ifrNNVal",packFloat((float)ifrNNVal,10)); // rms ~1% for muons

      // NA
      /*
      if (doemcnn) {

        float* pattern; pattern = new float[9];
        pattern[0] = latcalnn;
        pattern[1] = ecalnn;
        pattern[2] = eopcalnn;
        pattern[3] = ncrycalnn;
        pattern[4] = z20calnn;
        pattern[5] = z42calnn;
        pattern[6] = smomcalnn;
        pattern[7] = s1s9calnn;
        pattern[8] = s9s25calnn;

	float* outpattern; outpattern = new float[1];
        outpattern[0] = 0.0;

        if ( candMomentum <= 0.75 ) {
          if ( candTheta < 0.80 )                   NNout_p01_t1(pattern, outpattern);
          if ( candTheta >= 0.80 && candTheta <= 1.20 ) NNout_p01_t2(pattern, outpattern);
          if ( candTheta > 1.20 )                   NNout_p01_t3(pattern, outpattern);
        } else if ( candMomentum > 0.75 && candMomentum <= 1.25 ) {
          if ( candTheta < 0.80 )                   NNout_p02_03_t1(pattern, outpattern);
          if ( candTheta >= 0.80 && candTheta <= 1.20 ) NNout_p02_03_t2(pattern, outpattern);
          if ( candTheta > 1.20 )                   NNout_p02_03_t3(pattern, outpattern);
        } else if ( candMomentum > 1.25 && candMomentum <= 2.00 ) {
          if ( candTheta < 0.80 )                   NNout_p04_06_t1(pattern, outpattern);
          if ( candTheta >= 0.80 && candTheta <= 1.20 ) NNout_p04_06_t2(pattern, outpattern);
          if ( candTheta > 1.20 )                   NNout_p04_06_t3(pattern, outpattern);
        } else if ( candMomentum > 2.00 && candMomentum <= 2.75 ) {
          if ( candTheta < 0.80 )                   NNout_p11_13_t1(pattern, outpattern);
          if ( candTheta >= 0.80 && candTheta <= 1.20 ) NNout_p11_13_t2(pattern, outpattern);
          if ( candTheta > 1.20 )                   NNout_p11_13_t3(pattern, outpattern);
        } else if ( candMomentum > 2.75 ) {
          if ( candTheta < 0.80 )                   NNout_p14_18_t1(pattern, outpattern);
          if ( candTheta >= 0.80 && candTheta <= 1.20 ) NNout_p14_18_t2(pattern, outpattern);
          if ( candTheta > 1.20 )                   NNout_p14_18_t3(pattern, outpattern);
        }

        nt->column("ifremcNNVal",packFloat((float)outpattern[0],12)); // resolution ~100%

        delete [] pattern;
        delete [] outpattern;
      } else {
        nt->column("ifremcNNVal",(float)-0.5);
      }
      */

    } else {
      nt->column("logmulike",(float) -1000.);
      nt->column("logpilike",(float) +1000.);
      nt->column("ifrns",-1);
      nt->column("ifrexpintlen",(float)-1.0);
      nt->column("ifrmeasintlen",(float)-1.0);
      nt->column("ifrintlenbeforeiron",(float)-1.0);
      nt->column("ifrmatchchi2",(float)-1.0);
      nt->column("ifremcmatch",(float)-1.0);
      nt->column("ifrfitchi2",(float)-1.0);
      nt->column("ifrlayhits",(int) 1);
      nt->column("ifrfirsthit",(int) -1);
      nt->column("ifrlasthit",(int) -1);
      nt->column("ifrhasinner",false);
      nt->column("ifrhasbarrel",false);
      nt->column("ifrhasfwd",false);
      nt->column("ifrhasbwd",false);
      nt->column("ifrcont",(float)-1.0);
      nt->column("ifrsigmu",(float)-1.0);
      nt->column("ifrNNVal",(float)-1.0);
      nt->column("ifremcNNVal",(float)-1.0);
    }
  }

  // PID
  if (getParmValue("writePidInfo")) {
    if (pInfo) {
      nt->column("svtecons",packFloat((float)consval(c,PdtPid::electron,PidSystem::svt),12)); // -1 or from 0 to 1, ~100% resolution
      nt->column("svtmucons",packFloat((float)consval(c,PdtPid::muon,PidSystem::svt),12));
      nt->column("svtpicons",packFloat((float)consval(c,PdtPid::pion,PidSystem::svt),12));
      nt->column("svtkcons",packFloat((float)consval(c,PdtPid::kaon,PidSystem::svt),12));
      nt->column("svtpcons",packFloat((float)consval(c,PdtPid::proton,PidSystem::svt),12));
      nt->column("svteprob",packFloat((float) lhood(c,PdtPid::electron,PidSystem::svt),12));
      nt->column("svtmuprob",packFloat((float) lhood(c,PdtPid::muon,PidSystem::svt),12));
      nt->column("svtpiprob",packFloat((float) lhood(c,PdtPid::pion,PidSystem::svt),12));
      nt->column("svtkprob",packFloat((float) lhood(c,PdtPid::kaon,PidSystem::svt),12));
      nt->column("svtpprob",packFloat((float) lhood(c,PdtPid::proton,PidSystem::svt),12));
      
      nt->column("dchecons",packFloat((float)consval(c,PdtPid::electron,PidSystem::dch),12));
      nt->column("dchmucons",packFloat((float)consval(c,PdtPid::muon,PidSystem::dch),12));
      nt->column("dchpicons",packFloat((float)consval(c,PdtPid::pion,PidSystem::dch),12));
      nt->column("dchkcons",packFloat((float)consval(c,PdtPid::kaon,PidSystem::dch),12));
      nt->column("dchpcons",packFloat((float)consval(c,PdtPid::proton,PidSystem::dch),12));
      nt->column("dcheprob",packFloat((float) lhood(c,PdtPid::electron,PidSystem::dch),12));
      nt->column("dchmuprob",packFloat((float) lhood(c,PdtPid::muon,PidSystem::dch),12));
      nt->column("dchpiprob",packFloat((float) lhood(c,PdtPid::pion,PidSystem::dch),12));
      nt->column("dchkprob",packFloat((float) lhood(c,PdtPid::kaon,PidSystem::dch),12));
      nt->column("dchpprob",packFloat((float) lhood(c,PdtPid::proton,PidSystem::dch),12));
      
      nt->column("drcecons", (float)consval(c,PdtPid::electron,PidSystem::drc)); // already packed in CM2 - pretty crudely, 8 bits or so
      nt->column("drcmucons",(float)consval(c,PdtPid::muon,PidSystem::drc));
      nt->column("drcpicons",(float)consval(c,PdtPid::pion,PidSystem::drc));
      nt->column("drckcons",(float)consval(c,PdtPid::kaon,PidSystem::drc));
      nt->column("drcpcons",(float)consval(c,PdtPid::proton,PidSystem::drc));
      nt->column("drceprob",packFloat((float) lhood(c,PdtPid::electron,PidSystem::drc),12)); // 
      nt->column("drcmuprob",packFloat((float) lhood(c,PdtPid::muon,PidSystem::drc),12));
      nt->column("drcpiprob",packFloat((float) lhood(c,PdtPid::pion,PidSystem::drc),12));
      nt->column("drckprob",packFloat((float) lhood(c,PdtPid::kaon,PidSystem::drc),12));
      nt->column("drcpprob",packFloat((float) lhood(c,PdtPid::proton,PidSystem::drc),12));
      
      nt->column("ifrecons",(float)consval(c,PdtPid::electron,PidSystem::ifr));  // already packed in CM2 - pretty crudely, 8 bits or so
      nt->column("ifrmucons",(float)consval(c,PdtPid::muon,PidSystem::ifr));
      nt->column("ifrpicons",(float)consval(c,PdtPid::pion,PidSystem::ifr));
      nt->column("ifrkcons",(float)consval(c,PdtPid::kaon,PidSystem::ifr));
      nt->column("ifrpcons",(float)consval(c,PdtPid::proton,PidSystem::ifr));
      nt->column("ifreprob",packFloat((float) lhood(c,PdtPid::electron,PidSystem::ifr),12)); // 
      nt->column("ifrmuprob",packFloat((float) lhood(c,PdtPid::muon,PidSystem::ifr),12));
      nt->column("ifrpiprob",packFloat((float) lhood(c,PdtPid::pion,PidSystem::ifr),12));
      nt->column("ifrkprob",packFloat((float) lhood(c,PdtPid::kaon,PidSystem::ifr),12));
      nt->column("ifrpprob",packFloat((float) lhood(c,PdtPid::proton,PidSystem::ifr),12));

    } else {
      nt->column("svtecons",(float)-99.0);
      nt->column("svtmucons",(float)-99.0);
      nt->column("svtpicons",(float)-99.0);
      nt->column("svtkcons",(float)-99.0);
      nt->column("svtpcons",(float)-99.0);
      nt->column("svteprob",(float) -99.0);
      nt->column("svtmuprob",(float) -99.0);
      nt->column("svtpiprob",(float) -99.0);
      nt->column("svtkprob",(float) -99.0);
      nt->column("svtpprob",(float) -99.0);
      nt->column("dchecons",(float)-99.0);
      nt->column("dchmucons",(float)-99.0);
      nt->column("dchpicons",(float)-99.0);
      nt->column("dchkcons",(float)-99.0);
      nt->column("dchpcons",(float)-99.0);
      nt->column("dcheprob",(float) -99.0);
      nt->column("dchmuprob",(float) -99.0);
      nt->column("dchpiprob",(float) -99.0);
      nt->column("dchkprob",(float) -99.0);
      nt->column("dchpprob",(float) -99.0);
      nt->column("drcecons",(float)-99.0);
      nt->column("drcmucons",(float)-99.0);
      nt->column("drcpicons",(float)-99.0);
      nt->column("drckcons",(float)-99.0);
      nt->column("drcpcons",(float)-99.0);
      nt->column("drceprob",(float) -99.0);
      nt->column("drcmuprob",(float) -99.0);
      nt->column("drcpiprob",(float) -99.0);
      nt->column("drckprob",(float) -99.0);
      nt->column("drcpprob",(float) -99.0);
      nt->column("ifrecons",(float)-99.0);
      nt->column("ifrmucons",(float)-99.0);
      nt->column("ifrpicons",(float)-99.0);
      nt->column("ifrkcons",(float)-99.0);
      nt->column("ifrpcons",(float)-99.0);
      nt->column("ifreprob",(float) -99.0);
      nt->column("ifrmuprob",(float) -99.0);
      nt->column("ifrpiprob",(float) -99.0);
      nt->column("ifrkprob",(float) -99.0);
      nt->column("ifrpprob",(float) -99.0);
    }
  }


  // expected dedx
  if (getParmValue("writeExpectedDedx")) {
  // avtelnov, 2007/05/13
  // Direct calls to DchBetheBloch::ionization(mom,type) should be limited to internal DCH code.  
  //  const DchBetheBloch*  BetheBloch = gblEnv->getBta()->dchBetheBloch();
  //  nt->column("dedxeOld",(float)BetheBloch->ionization(candMomentum+pQual->deltaDchMomentum(),PdtPid::electron));
  //  nt->column("dedxpiOld",(float)BetheBloch->ionization(candMomentum+pQual->deltaDchMomentum(),PdtPid::pion));
  //  nt->column("dedxkOld",(float)BetheBloch->ionization(candMomentum+pQual->deltaDchMomentum(),PdtPid::kaon));
  //  nt->column("dedxpOld",(float)BetheBloch->ionization(candMomentum+pQual->deltaDchMomentum(),PdtPid::proton));
  //  nt->column("dedxmuOld",(float)BetheBloch->ionization(candMomentum+pQual->deltaDchMomentum(),PdtPid::muon));

// Also, store momentum in the middle of SVT (presumably, at 7 cm for the pion hypothesis)
// Let the default be the momentum at the IP. 
    float psvt = pIPpi;
    
    float dedxe=0., dedxpi=0., dedxk=0., dedxp=0., dedxmu=0.; // expected DCH dE/dx from DchPidInfo
  // expected DCH dE/dx error from DchPidInfo
    float dedxErr=999.;
  //float dedxeErr=999., dedxpiErr=999., dedxkErr=999., dedxpErr=999., dedxmuErr=999.;
  // New with Korneliy Todushev's OptimalTM: hypothesis-dependent ***measured*** dE/dx, nsampdch, and pdch 
  // (dedxdch = dedxdchTMpi, nsampdch = nsampdchTMpi, pdch = pdchTMpi)
  // The "TM" variables are hypothesis-dependent, would be applicable if the OptimalTM 
  //float dedxdchTMe=0., dedxdchTMpi=0., dedxdchTMk=0., dedxdchTMp=0., dedxdchTMmu=0.;
  //int nsampdchTMe=0, nsampdchTMpi=0, nsampdchTMk=0, nsampdchTMp=0, nsampdchTMmu=0;
  //float pdchTMe=pdch, pdchTMpi=pdch, pdchTMk=pdch, pdchTMp=pdch, pdchTMmu=pdch;

    const PidInfoSummary *realPidInfo = c->pidInfoSummary();
    if (realPidInfo != 0) {
      const DchPidInfo* dchPid = realPidInfo->dchPidInfo();
      if (dchPid != 0) {
	dedxe = dchPid->getExpected(PdtPid::electron);
	dedxpi = dchPid->getExpected(PdtPid::pion);
	dedxk = dchPid->getExpected(PdtPid::kaon);
	dedxp = dchPid->getExpected(PdtPid::proton);
	dedxmu = dchPid->getExpected(PdtPid::muon);

      // as of 24.2.1, there is only one dE/dx error in DchPidInfo
      //dedxeErr = dchPid->getdEdxErr(PdtPid::electron);
      //dedxpiErr = dchPid->getdEdxErr(PdtPid::pion);
      //dedxkErr = dchPid->getdEdxErr(PdtPid::kaon);
      //dedxpErr = dchPid->getdEdxErr(PdtPid::proton);
      //dedxmuErr = dchPid->getdEdxErr(PdtPid::muon);
	dedxErr = dchPid->getdEdxErr(PdtPid::pion);

      // as of 24.2.1, there is only one measured dE/dx, it's dedxdch  
      //dedxdchTMe = dchPid->getdEdxVal(PdtPid::electron);
      //dedxdchTMpi = dchPid->getdEdxVal(PdtPid::pion); // == dedxdch
      //dedxdchTMk = dchPid->getdEdxVal(PdtPid::kaon);
      //dedxdchTMp = dchPid->getdEdxVal(PdtPid::proton);
      //dedxdchTMmu = dchPid->getdEdxVal(PdtPid::muon);

      // ... and only one number of dE/dx samples, it's nsampdch
      //nsampdchTMe = dchPid->getSamples(PdtPid::electron);
      //nsampdchTMpi = dchPid->getSamples(PdtPid::pion); // == nsampdch
      //nsampdchTMk = dchPid->getSamples(PdtPid::kaon);
      //nsampdchTMp = dchPid->getSamples(PdtPid::proton);
      //nsampdchTMmu = dchPid->getSamples(PdtPid::muon);

      // ... and only one momentum at DCH (which is fine), it's pdch 
      //pdchTMe = dchPid->getDchMomentum(PdtPid::electron);
      //pdchTMpi = dchPid->getDchMomentum(PdtPid::pion); // == pdch
      //pdchTMk = dchPid->getDchMomentum(PdtPid::kaon);
      //pdchTMp = dchPid->getDchMomentum(PdtPid::proton);
      //pdchTMmu = dchPid->getDchMomentum(PdtPid::muon);
      }
      const SvtPidInfo* svtPid = realPidInfo->svtPidInfo();
      if (svtPid != 0 && svtPid->momentum()<3.9 ) { // psvt is chopped off in CM2 at 4.0 GeV/c!!! 
	psvt = svtPid->momentum();
      }
    }

// This is what we get from the DCH:
    nt->column("dedxe",packFloat(dedxe,12));  // DCH dE/dx resolution is 10%, so adding 1e-4 in quadrature does nothing
    nt->column("dedxpi",packFloat(dedxpi,12));
    nt->column("dedxk",packFloat(dedxk,12));
    nt->column("dedxp",packFloat(dedxp,12));
    nt->column("dedxmu",packFloat(dedxmu,12));

  //nt->column("dedxeErr",packFloat(dedxeErr,14)); // Error on the error is ~50%, so adding 4e-4 in quadrature does nothing
  //nt->column("dedxpiErr",packFloat(dedxpiErr,14));
  //nt->column("dedxkErr",packFloat(dedxkErr,14));
  //nt->column("dedxpErr",packFloat(dedxpErr,14));
  //nt->column("dedxmuErr",packFloat(dedxmuErr,14));
    nt->column("dedxErr",packFloat(dedxErr,14));

  //nt->column("dedxdchTMe",packFloat(dedxdchTMe,10)); // 10% resolution, so adding 2.5e-5 in quadrature does nothing
  //nt->column("dedxdchTMpi",packFloat(dedxdchTMpi,10));
  //nt->column("dedxdchTMk",packFloat(dedxdchTMk,10));
  //nt->column("dedxdchTMp",packFloat(dedxdchTMp,10));
  //nt->column("dedxdchTMmu",packFloat(dedxdchTMmu,10));

  //nt->column("nsampdchTMe",nsampdchTMe);
  //nt->column("nsampdchTMpi",nsampdchTMpi);
  //nt->column("nsampdchTMk",nsampdchTMk);
  //nt->column("nsampdchTMp",nsampdchTMp);
  //nt->column("nsampdchTMmu",nsampdchTMmu);

  //nt->column("pdchTMe",packFloat(pdchTMe,8));
  //nt->column("pdchTMpi",packFloat(pdchTMpi,8));
  //nt->column("pdchTMk",packFloat(pdchTMk,8));
  //nt->column("pdchTMp",packFloat(pdchTMp,8));
  //nt->column("pdchTMmu",packFloat(pdchTMmu,8));

    nt->column("psvt",packFloat(psvt,8));


  // sasha telnov's dedx parameters, coded by kevin 

    //static const float me(0.000511); // particles masses 
    //static const float mmu(0.105658);
    //static const float mpi(0.13957);
    //static const float mk(0.493677);
    //static const float mp(0.938272);

    float dedxdche(0.); // expected values of DCH dE/dx
    float dedxdchmu(0.);
    float dedxdchpi(0.);
    float dedxdchk(0.);
    float dedxdchp(0.);
    float dedxsvte(0.); // same for SVT
    float dedxsvtmu(0.);
    float dedxsvtpi(0.);
    float dedxsvtk(0.);
    float dedxsvtp(0.);

  //float phidedxdche(0.); // phi corrections to the expected values of DCH dE/dx, have no business being in PID ntuples
  //float phidedxdchmu(0.);
  //float phidedxdchpi(0.);
  //float phidedxdchk(0.);
  //float phidedxdchp(0.);
  //float phidedxsvte(0.); // same for SVT
  //float phidedxsvtmu(0.);
  //float phidedxsvtpi(0.);
  //float phidedxsvtk(0.);
  //float phidedxsvtp(0.);

    float dEdxdchPulle(-999.); // pull = ((measured-expected)/measured - relShift)/relError
  //float dEdxdchErrore(0.); // error = measured * relError 
  //float dEdxdchRelativeErrore(0.); // relative error = rms of (measured-expected)/measured
  //float dEdxdchRelativeShifte(0.); // relative shift = mean of (measured-expected)/measured
    float dEdxsvtPulle(-999.); // same for SVT
  //float dEdxsvtErrore(0.);
  //float dEdxsvtRelativeErrore(0.);
  //float dEdxsvtRelativeShifte(0.);

    float dEdxdchPullmu(-999.); // same for muons
  //float dEdxdchErrormu(0.);
  //float dEdxdchRelativeErrormu(0.);
  //float dEdxdchRelativeShiftmu(0.);
    float dEdxsvtPullmu(-999.);
  //float dEdxsvtErrormu(0.);
  //float dEdxsvtRelativeErrormu(0.);
  //float dEdxsvtRelativeShiftmu(0.);

    float dEdxdchPullpi(-999.); // same for pions
  //float dEdxdchErrorpi(0.);
  //float dEdxdchRelativeErrorpi(0.);
  //float dEdxdchRelativeShiftpi(0.);
    float dEdxsvtPullpi(-999.);
  //float dEdxsvtErrorpi(0.);
  //float dEdxsvtRelativeErrorpi(0.);
  //float dEdxsvtRelativeShiftpi(0.);

    float dEdxdchPullk(-999.); // same for kaons
  //float dEdxdchErrork(0.);
  //float dEdxdchRelativeErrork(0.);
  //float dEdxdchRelativeShiftk(0.);
    float dEdxsvtPullk(-999.);
  //float dEdxsvtErrork(0.);
  //float dEdxsvtRelativeErrork(0.);
  //float dEdxsvtRelativeShiftk(0.);

    float dEdxdchPullp(-999.); // same for protons
  //float dEdxdchErrorp(0.);
  //float dEdxdchRelativeErrorp(0.);
  //float dEdxdchRelativeShiftp(0.);
    float dEdxsvtPullp(-999.);
  //float dEdxsvtErrorp(0.);
  //float dEdxsvtRelativeErrorp(0.);
  //float dEdxsvtRelativeShiftp(0.);

    // NA
    /*
    if (_useCDBforDEDX->value() == false) {
      dedxdche  = PidBetheBlochPhiCorrDch_AVT::dEdxExpected(pdchPi, me, candThetaPi, candPhiPi, charge_pi, nrun, date);
      dedxdchmu = PidBetheBlochPhiCorrDch_AVT::dEdxExpected(pdchPi, mmu, candThetaPi, candPhiPi, charge_pi, nrun, date);
      dedxdchpi = PidBetheBlochPhiCorrDch_AVT::dEdxExpected(pdchPi, mpi, candThetaPi, candPhiPi, charge_pi, nrun, date);
      dedxdchk  = PidBetheBlochPhiCorrDch_AVT::dEdxExpected(pdchPi, mk, candThetaPi, candPhiPi, charge_pi, nrun, date);
      dedxdchp  = PidBetheBlochPhiCorrDch_AVT::dEdxExpected(pdchPi, mp, candThetaPi, candPhiPi, charge_pi, nrun, date);
      dedxsvte  = PidBetheBlochPhiCorrSvt_AVT::dEdxExpected(psvt, me, candThetaPi, candPhiPi, charge_pi, nrun, date);
      dedxsvtmu = PidBetheBlochPhiCorrSvt_AVT::dEdxExpected(psvt, mmu, candThetaPi, candPhiPi, charge_pi, nrun, date);
      dedxsvtpi = PidBetheBlochPhiCorrSvt_AVT::dEdxExpected(psvt, mpi, candThetaPi, candPhiPi, charge_pi, nrun, date);
      dedxsvtk  = PidBetheBlochPhiCorrSvt_AVT::dEdxExpected(psvt, mk, candThetaPi, candPhiPi, charge_pi, nrun, date);
      dedxsvtp  = PidBetheBlochPhiCorrSvt_AVT::dEdxExpected(psvt, mp, candThetaPi, candPhiPi, charge_pi, nrun, date);
    //phidedxdche  = PidBetheBlochPhiCorrDch_AVT::dEdxPhiCorrection(pdchPi, me, candThetaPi, candPhiPi, charge_pi, nrun, date);
    //phidedxdchmu = PidBetheBlochPhiCorrDch_AVT::dEdxPhiCorrection(pdchPi, mmu, candThetaPi, candPhiPi, charge_pi, nrun, date);
    //phidedxdchpi = PidBetheBlochPhiCorrDch_AVT::dEdxPhiCorrection(pdchPi, mpi, candThetaPi, candPhiPi, charge_pi, nrun, date);
    //phidedxdchk  = PidBetheBlochPhiCorrDch_AVT::dEdxPhiCorrection(pdchPi, mk, candThetaPi, candPhiPi, charge_pi, nrun, date);
    //phidedxdchp  = PidBetheBlochPhiCorrDch_AVT::dEdxPhiCorrection(pdchPi, mp, candThetaPi, candPhiPi, charge_pi, nrun, date);
    //phidedxsvte  = PidBetheBlochPhiCorrSvt_AVT::dEdxPhiCorrection(psvt, me, candThetaPi, candPhiPi, charge_pi, nrun, date);
    //phidedxsvtmu = PidBetheBlochPhiCorrSvt_AVT::dEdxPhiCorrection(psvt, mmu, candThetaPi, candPhiPi, charge_pi, nrun, date);
    //phidedxsvtpi = PidBetheBlochPhiCorrSvt_AVT::dEdxPhiCorrection(psvt, mpi, candThetaPi, candPhiPi, charge_pi, nrun, date);
    //phidedxsvtk  = PidBetheBlochPhiCorrSvt_AVT::dEdxPhiCorrection(psvt, mk, candThetaPi, candPhiPi, charge_pi, nrun, date);
    //phidedxsvtp  = PidBetheBlochPhiCorrSvt_AVT::dEdxPhiCorrection(psvt, mp, candThetaPi, candPhiPi, charge_pi, nrun, date);
    }
    else {
      dedxdche  = PidDEDXCdbDch::dEdxExpected(pdchPi, me, candThetaPi, candPhiPi, charge_pi);
      dedxdchmu = PidDEDXCdbDch::dEdxExpected(pdchPi, mmu, candThetaPi, candPhiPi, charge_pi);
      dedxdchpi = PidDEDXCdbDch::dEdxExpected(pdchPi, mpi, candThetaPi, candPhiPi, charge_pi);
      dedxdchk  = PidDEDXCdbDch::dEdxExpected(pdchPi, mk, candThetaPi, candPhiPi, charge_pi);
      dedxdchp  = PidDEDXCdbDch::dEdxExpected(pdchPi, mp, candThetaPi, candPhiPi, charge_pi);
      dedxsvte  = PidDEDXCdbSvt::dEdxExpected(psvt, me, candThetaPi, candPhiPi, charge_pi);
      dedxsvtmu = PidDEDXCdbSvt::dEdxExpected(psvt, mmu, candThetaPi, candPhiPi, charge_pi);
      dedxsvtpi = PidDEDXCdbSvt::dEdxExpected(psvt, mpi, candThetaPi, candPhiPi, charge_pi);
      dedxsvtk  = PidDEDXCdbSvt::dEdxExpected(psvt, mk, candThetaPi, candPhiPi, charge_pi);
      dedxsvtp  = PidDEDXCdbSvt::dEdxExpected(psvt, mp, candThetaPi, candPhiPi, charge_pi);
    }
    */
    
    if (_warningCounter<5 && nsampdch>0 && dedxdch<10) {
      _warningCounter++;
      ErrMsg(warning) << "nsampdch=" << nsampdch << " and dedxdch=" << dedxdch << ", which is far too low and is a known, very rare problem." << endmsg;  
      ErrMsg(warning) << "I will record this information anyway. Disregard the \"should not be used\" warning that follows." << endmsg;  
    }

    // NA
    /*
    if (_useCDBforDEDX->value() == false) {
      if (nsampdch>0) {
	dEdxdchPulle = PidBetheBlochErrorDch_AVT::dEdxPull(dedxdch, nsampdch, pdchPi, me, candThetaPi, candPhiPi, charge_pi, nrun, date);
      //dEdxdchErrore = PidBetheBlochErrorDch_AVT::dEdxError(dedxdch,        nsampdch, pdchPi, me, candThetaPi,      charge_pi, nrun, date);
      //dEdxdchRelativeErrore = PidBetheBlochErrorDch_AVT::dEdxRelativeError(nsampdch, pdchPi, me, candThetaPi,      charge_pi, nrun, date);
      //dEdxdchRelativeShifte = PidBetheBlochErrorDch_AVT::dEdxRelativeShift(nsampdch, pdchPi, me, candThetaPi,      charge_pi, nrun, date);

      // use pi-hypo charge for muons. There is a one-in-a-million difference. 
	dEdxdchPullmu = PidBetheBlochErrorDch_AVT::dEdxPull(dedxdch, nsampdch, pdchPi, mmu, candThetaPi, candPhiPi, charge_pi, nrun, date);
      //dEdxdchErrormu = PidBetheBlochErrorDch_AVT::dEdxError(dedxdch,        nsampdch, pdchPi, mmu, candThetaPi,      charge_pi, nrun, date);
      //dEdxdchRelativeErrormu = PidBetheBlochErrorDch_AVT::dEdxRelativeError(nsampdch, pdchPi, mmu, candThetaPi,      charge_pi, nrun, date);
      //dEdxdchRelativeShiftmu = PidBetheBlochErrorDch_AVT::dEdxRelativeShift(nsampdch, pdchPi, mmu, candThetaPi,      charge_pi, nrun, date);

	dEdxdchPullpi = PidBetheBlochErrorDch_AVT::dEdxPull(dedxdch, nsampdch, pdchPi, mpi, candThetaPi, candPhiPi, charge_pi, nrun, date);
      //dEdxdchErrorpi = PidBetheBlochErrorDch_AVT::dEdxError(dedxdch,        nsampdch, pdchPi, mpi, candThetaPi,      charge_pi, nrun, date);
      //dEdxdchRelativeErrorpi = PidBetheBlochErrorDch_AVT::dEdxRelativeError(nsampdch, pdchPi, mpi, candThetaPi,      charge_pi, nrun, date);
      //dEdxdchRelativeShiftpi = PidBetheBlochErrorDch_AVT::dEdxRelativeShift(nsampdch, pdchPi, mpi, candThetaPi,      charge_pi, nrun, date);

	dEdxdchPullk = PidBetheBlochErrorDch_AVT::dEdxPull(dedxdch, nsampdch, pdchPi, mk, candThetaPi, candPhiPi, charge_pi, nrun, date);
      //dEdxdchErrork = PidBetheBlochErrorDch_AVT::dEdxError(dedxdch,        nsampdch, pdchPi, mk, candThetaPi,      charge_pi, nrun, date);
      //dEdxdchRelativeErrork = PidBetheBlochErrorDch_AVT::dEdxRelativeError(nsampdch, pdchPi, mk, candThetaPi,      charge_pi, nrun, date);
      //dEdxdchRelativeShiftk = PidBetheBlochErrorDch_AVT::dEdxRelativeShift(nsampdch, pdchPi, mk, candThetaPi,      charge_pi, nrun, date);
	
	dEdxdchPullp = PidBetheBlochErrorDch_AVT::dEdxPull(dedxdch, nsampdch, pdchPi, mp, candThetaPi, candPhiPi, charge_pi, nrun, date);
      //dEdxdchErrorp = PidBetheBlochErrorDch_AVT::dEdxError(dedxdch,        nsampdch, pdchPi, mp, candThetaPi,      charge_pi, nrun, date);
      //dEdxdchRelativeErrorp = PidBetheBlochErrorDch_AVT::dEdxRelativeError(nsampdch, pdchPi, mp, candThetaPi,      charge_pi, nrun, date);
      //dEdxdchRelativeShiftp = PidBetheBlochErrorDch_AVT::dEdxRelativeShift(nsampdch, pdchPi, mp, candThetaPi,      charge_pi, nrun, date);
      }
      if (nsampsvt>0) {
	dEdxsvtPulle  = PidBetheBlochErrorSvt_AVT::dEdxPull(dedxsvt, nsampsvt, psvt, me, candThetaPi, candPhiPi, charge_pi, nrun, date);
      //dEdxsvtErrore = PidBetheBlochErrorSvt_AVT::dEdxError(dedxsvt,        nsampsvt, psvt, me, candThetaPi,      charge_pi, nrun, date);
      //dEdxsvtRelativeErrore = PidBetheBlochErrorSvt_AVT::dEdxRelativeError(nsampsvt, psvt, me, candThetaPi,      charge_pi, nrun, date);
      //dEdxsvtRelativeShifte = PidBetheBlochErrorSvt_AVT::dEdxRelativeShift(nsampsvt, psvt, me, candThetaPi,      charge_pi, nrun, date);
	
	dEdxsvtPullmu  = PidBetheBlochErrorSvt_AVT::dEdxPull(dedxsvt, nsampsvt, psvt, mmu, candThetaPi, candPhiPi, charge_pi, nrun, date);
      //dEdxsvtErrormu = PidBetheBlochErrorSvt_AVT::dEdxError(dedxsvt,        nsampsvt, psvt, mmu, candThetaPi,      charge_pi, nrun, date);
      //dEdxsvtRelativeErrormu = PidBetheBlochErrorSvt_AVT::dEdxRelativeError(nsampsvt, psvt, mmu, candThetaPi,      charge_pi, nrun, date);
      //dEdxsvtRelativeShiftmu = PidBetheBlochErrorSvt_AVT::dEdxRelativeShift(nsampsvt, psvt, mmu, candThetaPi,      charge_pi, nrun, date);
	
	dEdxsvtPullpi  = PidBetheBlochErrorSvt_AVT::dEdxPull(dedxsvt, nsampsvt, psvt, mpi, candThetaPi, candPhiPi, charge_pi, nrun, date);
      //dEdxsvtErrorpi = PidBetheBlochErrorSvt_AVT::dEdxError(dedxsvt,        nsampsvt, psvt, mpi, candThetaPi,      charge_pi, nrun, date);
      //dEdxsvtRelativeErrorpi = PidBetheBlochErrorSvt_AVT::dEdxRelativeError(nsampsvt, psvt, mpi, candThetaPi,      charge_pi, nrun, date);
      //dEdxsvtRelativeShiftpi = PidBetheBlochErrorSvt_AVT::dEdxRelativeShift(nsampsvt, psvt, mpi, candThetaPi,      charge_pi, nrun, date);
	
	dEdxsvtPullk  = PidBetheBlochErrorSvt_AVT::dEdxPull(dedxsvt, nsampsvt, psvt, mk, candThetaPi, candPhiPi, charge_pi, nrun, date);
      //dEdxsvtErrork = PidBetheBlochErrorSvt_AVT::dEdxError(dedxsvt,        nsampsvt, psvt, mk, candThetaPi,      charge_pi, nrun, date);
      //dEdxsvtRelativeErrork = PidBetheBlochErrorSvt_AVT::dEdxRelativeError(nsampsvt, psvt, mk, candThetaPi,      charge_pi, nrun, date);
      //dEdxsvtRelativeShiftk = PidBetheBlochErrorSvt_AVT::dEdxRelativeShift(nsampsvt, psvt, mk, candThetaPi,      charge_pi, nrun, date);
	
	dEdxsvtPullp  = PidBetheBlochErrorSvt_AVT::dEdxPull(dedxsvt, nsampsvt, psvt, mp, candThetaPi, candPhiPi, charge_pi, nrun, date);
      //dEdxsvtErrorp = PidBetheBlochErrorSvt_AVT::dEdxError(dedxsvt,        nsampsvt, psvt, mp, candThetaPi,      charge_pi, nrun, date);
      //dEdxsvtRelativeErrorp = PidBetheBlochErrorSvt_AVT::dEdxRelativeError(nsampsvt, psvt, mp, candThetaPi,      charge_pi, nrun, date);
      //dEdxsvtRelativeShiftp = PidBetheBlochErrorSvt_AVT::dEdxRelativeShift(nsampsvt, psvt, mp, candThetaPi,      charge_pi, nrun, date);
      }
    } else { // get detailed dE/dx from CDB; comment the errors out once the R24 dust settles
      if (nsampdch>0) {
	dEdxdchPulle = PidDEDXCdbDch::dEdxPull(dedxdch, nsampdch, pdchPi, me, candThetaPi, candPhiPi, charge_pi);
	dEdxdchPullmu = PidDEDXCdbDch::dEdxPull(dedxdch, nsampdch, pdchPi, mmu, candThetaPi, candPhiPi, charge_pi);
	dEdxdchPullpi = PidDEDXCdbDch::dEdxPull(dedxdch, nsampdch, pdchPi, mpi, candThetaPi, candPhiPi, charge_pi);
	dEdxdchPullk = PidDEDXCdbDch::dEdxPull(dedxdch, nsampdch, pdchPi, mk, candThetaPi, candPhiPi, charge_pi);
	dEdxdchPullp = PidDEDXCdbDch::dEdxPull(dedxdch, nsampdch, pdchPi, mp, candThetaPi, candPhiPi, charge_pi);
	
      //dEdxdchErrore = PidDEDXCdbDch::dEdxError(dedxdch, nsampdch, pdchPi, me, candThetaPi, charge_pi, candPhiPi);
      //dEdxdchErrormu = PidDEDXCdbDch::dEdxError(dedxdch, nsampdch, pdchPi, mmu, candThetaPi, charge_pi, candPhiPi);
      //dEdxdchErrorpi = PidDEDXCdbDch::dEdxError(dedxdch, nsampdch, pdchPi, mpi, candThetaPi, charge_pi, candPhiPi);
      //dEdxdchErrork = PidDEDXCdbDch::dEdxError(dedxdch, nsampdch, pdchPi, mk, candThetaPi, charge_pi, candPhiPi);
      //dEdxdchErrorp = PidDEDXCdbDch::dEdxError(dedxdch, nsampdch, pdchPi, mp, candThetaPi, charge_pi, candPhiPi);
      }
      if (nsampsvt>0) {
	dEdxsvtPulle  = PidDEDXCdbSvt::dEdxPull(dedxsvt, nsampsvt, psvt, me, candThetaPi, candPhiPi, charge_pi);
	dEdxsvtPullmu  = PidDEDXCdbSvt::dEdxPull(dedxsvt, nsampsvt, psvt, mmu, candThetaPi, candPhiPi, charge_pi);
	dEdxsvtPullpi  = PidDEDXCdbSvt::dEdxPull(dedxsvt, nsampsvt, psvt, mpi, candThetaPi, candPhiPi, charge_pi);
	dEdxsvtPullk  = PidDEDXCdbSvt::dEdxPull(dedxsvt, nsampsvt, psvt, mk, candThetaPi, candPhiPi, charge_pi);
	dEdxsvtPullp  = PidDEDXCdbSvt::dEdxPull(dedxsvt, nsampsvt, psvt, mp, candThetaPi, candPhiPi, charge_pi);
	
      //dEdxsvtErrore = PidDEDXCdbSvt::dEdxError(dedxsvt, nsampsvt, psvt, me, candThetaPi, charge_pi, candPhiPi);
      //dEdxsvtErrormu = PidDEDXCdbSvt::dEdxError(dedxsvt, nsampsvt, psvt, mmu, candThetaPi, charge_pi, candPhiPi);
      //dEdxsvtErrorpi = PidDEDXCdbSvt::dEdxError(dedxsvt, nsampsvt, psvt, mpi, candThetaPi, charge_pi, candPhiPi);
      //dEdxsvtErrork = PidDEDXCdbSvt::dEdxError(dedxsvt, nsampsvt, psvt, mk, candThetaPi, charge_pi, candPhiPi);
      //dEdxsvtErrorp = PidDEDXCdbSvt::dEdxError(dedxsvt, nsampsvt, psvt, mp, candThetaPi, charge_pi, candPhiPi);
      }  
    }
    */
 
// These are values computed with standalone dE/dx code, not subsystem code   
    nt->column("dedxdche",packFloat((float)dedxdche,12));
    nt->column("dedxdchmu",packFloat((float)dedxdchmu,12));
    nt->column("dedxdchpi",packFloat((float)dedxdchpi,12));
    nt->column("dedxdchk",packFloat((float)dedxdchk,12));
    nt->column("dedxdchp",packFloat((float)dedxdchp,12));
    nt->column("dedxsvte",packFloat((float)dedxsvte,12));
    nt->column("dedxsvtmu",packFloat((float)dedxsvtmu,12));
    nt->column("dedxsvtpi",packFloat((float)dedxsvtpi,12));
    nt->column("dedxsvtk",packFloat((float)dedxsvtk,12));
    nt->column("dedxsvtp",packFloat((float)dedxsvtp,12));

  //nt->column("phidedxdche",(float)phidedxdche);
  //nt->column("phidedxdchmu",(float)phidedxdchmu);
  //nt->column("phidedxdchpi",(float)phidedxdchpi);
  //nt->column("phidedxdchk",(float)phidedxdchk);
  //nt->column("phidedxdchp",(float)phidedxdchp);
  //nt->column("phidedxsvte",(float)phidedxsvte);
  //nt->column("phidedxsvtmu",(float)phidedxsvtmu);
  //nt->column("phidedxsvtpi",(float)phidedxsvtpi);
  //nt->column("phidedxsvtk",(float)phidedxsvtk);
  //nt->column("phidedxsvtp",(float)phidedxsvtp);

// pull width is 1, so packing with 10 bits for the significand should be more than enough 
// (an rms of 2e-4 is added in quadrature)
    nt->column("dEdxdchPulle",packFloat((float)dEdxdchPulle,13)); 
  //nt->column("dEdxdchErrore",(float)dEdxdchErrore);
  //nt->column("dEdxdchRelativeErrore",(float)dEdxdchRelativeErrore);
  //nt->column("dEdxdchRelativeShifte",(float)dEdxdchRelativeShifte);
    nt->column("dEdxsvtPulle",packFloat((float)dEdxsvtPulle,13));
  //nt->column("dEdxsvtErrore",(float)dEdxsvtErrore);
  //nt->column("dEdxsvtRelativeErrore",(float)dEdxsvtRelativeErrore);
  //nt->column("dEdxsvtRelativeShifte",(float)dEdxsvtRelativeShifte);

    nt->column("dEdxdchPullmu",packFloat((float)dEdxdchPullmu,13));
  //nt->column("dEdxdchErrormu",(float)dEdxdchErrormu);
  //nt->column("dEdxdchRelativeErrormu",(float)dEdxdchRelativeErrormu);
  //nt->column("dEdxdchRelativeShiftmu",(float)dEdxdchRelativeShiftmu);
    nt->column("dEdxsvtPullmu",packFloat((float)dEdxsvtPullmu,13));
  //nt->column("dEdxsvtErrormu",(float)dEdxsvtErrormu);
  //nt->column("dEdxsvtRelativeErrormu",(float)dEdxsvtRelativeErrormu);
  //nt->column("dEdxsvtRelativeShiftmu",(float)dEdxsvtRelativeShiftmu);

    nt->column("dEdxdchPullk",packFloat((float)dEdxdchPullk,13));
  //nt->column("dEdxdchErrork",(float)dEdxdchErrork);
  //nt->column("dEdxdchRelativeErrork",(float)dEdxdchRelativeErrork);
  //nt->column("dEdxdchRelativeShiftk",(float)dEdxdchRelativeShiftk);
    nt->column("dEdxsvtPullk",packFloat((float)dEdxsvtPullk,13));
  //nt->column("dEdxsvtErrork",(float)dEdxsvtErrork);
  //nt->column("dEdxsvtRelativeErrork",(float)dEdxsvtRelativeErrork);
  //nt->column("dEdxsvtRelativeShiftk",(float)dEdxsvtRelativeShiftk);

    nt->column("dEdxdchPullpi",packFloat((float)dEdxdchPullpi,13));
  //nt->column("dEdxdchErrorpi",(float)dEdxdchErrorpi);
  //nt->column("dEdxdchRelativeErrorpi",(float)dEdxdchRelativeErrorpi);
  //nt->column("dEdxdchRelativeShiftpi",(float)dEdxdchRelativeShiftpi);
    nt->column("dEdxsvtPullpi",packFloat((float)dEdxsvtPullpi,13));
  //nt->column("dEdxsvtErrorpi",(float)dEdxsvtErrorpi);
  //nt->column("dEdxsvtRelativeErrorpi",(float)dEdxsvtRelativeErrorpi);
  //nt->column("dEdxsvtRelativeShiftpi",(float)dEdxsvtRelativeShiftpi);

    nt->column("dEdxdchPullp",packFloat((float)dEdxdchPullp,13));
  //nt->column("dEdxdchErrorp",(float)dEdxdchErrorp);
  //nt->column("dEdxdchRelativeErrorp",(float)dEdxdchRelativeErrorp);
  //nt->column("dEdxdchRelativeShiftp",(float)dEdxdchRelativeShiftp);
    nt->column("dEdxsvtPullp",packFloat((float)dEdxsvtPullp,13));
  //nt->column("dEdxsvtErrorp",(float)dEdxsvtErrorp);
  //nt->column("dEdxsvtRelativeErrorp",(float)dEdxsvtRelativeErrorp);
  //nt->column("dEdxsvtRelativeShiftp",(float)dEdxsvtRelativeShiftp);
    
  }
  
// NN - val for kaon
// tempCand is setType'ed to pion, as the chargedTracks candidates PID selectors 
// see in the production sequence
  // NA
  //nt->column("knnval",packFloat((float) _theNN->nn_result(&tempCand),12));

//
  // NA
  //nt->column("kaonBDToutput",packFloat((float) _theKaonBDT->result(&tempCand),12));
  //double classout = -1.0;
  //_theKM->bdt_result(&tempCand, 0, classout);
  //nt->column("KMoutputKaon",(float) classout);
  //_theKM->bdt_result(&tempCand, 1, classout);
  //nt->column("KMoutputPion",(float) classout);
  //_theKM->bdt_result(&tempCand, 2, classout);
  //nt->column("KMoutputProton",(float) classout);
  //_theKM->bdt_result(&tempCand, 3, classout);
  //nt->column("KMoutputElectron",(float) classout);
  
// booleans to tell if a track from this event was used 
// in the training of the given BDT classifier. Compute and store a bool
// for each of the six classifiers (KMe, KMpi, KMk, KMp, MUmu, MUpi);
// use it later in the making of PID performance plots and tables 
  // NA
  /*
  nt->column("usedKMe",  (bool) UsedForBdtTraining::isUsed(UsedForBdtTraining::KMe, upperID,lowerID));
  nt->column("usedKMpi", (bool) UsedForBdtTraining::isUsed(UsedForBdtTraining::KMpi,upperID,lowerID));
  nt->column("usedKMk",  (bool) UsedForBdtTraining::isUsed(UsedForBdtTraining::KMk, upperID,lowerID));
  nt->column("usedKMp",  (bool) UsedForBdtTraining::isUsed(UsedForBdtTraining::KMp, upperID,lowerID));
  nt->column("usedMUmu", (bool) UsedForBdtTraining::isUsed(UsedForBdtTraining::MUmu,upperID,lowerID));
  nt->column("usedMUpi", (bool) UsedForBdtTraining::isUsed(UsedForBdtTraining::MUpi,upperID,lowerID));
  nt->column("usedBDTpi", (bool) UsedForBdtTraining::isUsed(UsedForBdtTraining::BDTpi,upperID,lowerID));
  nt->column("usedBDTk",  (bool) UsedForBdtTraining::isUsed(UsedForBdtTraining::BDTk, upperID,lowerID));
  */

// sms probabilities for de/dx
  if (getParmValue("writeSMSprobs")) {
    /*
    nt->column("smssvteprob",packFloat((float) _theSms->SvtProba(c,PdtPid::electron),12));
    nt->column("smssvtpiprob",packFloat((float) _theSms->SvtProba(c,PdtPid::pion),12));
    nt->column("smssvtkprob",packFloat((float) _theSms->SvtProba(c,PdtPid::kaon),12));
    nt->column("smssvtpprob",packFloat((float) _theSms->SvtProba(c,PdtPid::proton),12));
    nt->column("smssvtmuprob",packFloat((float) _theSms->SvtProba(c,PdtPid::muon),12));

    nt->column("smsdcheprob",packFloat((float) _theSms->DchProba(c,PdtPid::electron),12));
    nt->column("smsdchpiprob",packFloat((float) _theSms->DchProba(c,PdtPid::pion),12));
    nt->column("smsdchkprob",packFloat((float) _theSms->DchProba(c,PdtPid::kaon),12));
    nt->column("smsdchpprob",packFloat((float) _theSms->DchProba(c,PdtPid::proton),12));
    nt->column("smsdchmuprob",packFloat((float) _theSms->DchProba(c,PdtPid::muon),12));
    
    nt->column("smsdrceprob",packFloat((float) _theSms->DrcProba(c,PdtPid::electron),12));
    nt->column("smsdrcpiprob",packFloat((float) _theSms->DrcProba(c,PdtPid::pion),12));
    nt->column("smsdrckprob",packFloat((float) _theSms->DrcProba(c,PdtPid::kaon),12));
    nt->column("smsdrcpprob",packFloat((float) _theSms->DrcProba(c,PdtPid::proton),12));
    nt->column("smsdrcmuprob",packFloat((float) _theSms->DrcProba(c,PdtPid::muon),12));
    */
  }
  // likelihoods ratios for Kaons
  if (getParmValue("writeKLHratios")) {
    
    // NA
    //PidLHRatios* likeRatios=new PidLHRatios();
    //double likeKvsPi,likeKvsPro,likeKvsEle,likeProvsPi,likeProvsEle,likePivsEle;
    //likeRatios->ratios(c,likeKvsPi,likeKvsPro,likeKvsEle,likeProvsPi,likeProvsEle,likePivsEle);
    double likeKvsPi;
    double likeKvsPro;
    double likeKvsEle; 
    double likeProvsPi; 
    double likeProvsEle;
    double likePivsEle;
    bool isGoodDrc;
    computeDircLikelihoods( c, likeKvsPi, likeKvsPro, likeKvsEle,
			    likeProvsPi, likeProvsEle, 
			    likePivsEle, isGoodDrc );
    
    //cout << "Right after computeDircLikelihoods:"
    //<< " likeKvsPi = " << likeKvsPi 
    //<< " likeKvsPro = " << likeKvsPro 
    //<< " likeKvsEle = " << likeKvsEle 
    //<< " likeProvsPi = " << likeProvsPi
    //<< " likeProvsEle = " << likeProvsEle
    //<< " likePiVsEle = " << likePivsEle
    //<< endl;
    
    //if( likeKvsPi < 0 ) {
    //cout << "Problem: likeKvsPi = " << likeKvsPi << endl;
    //}
    //if( likeKvsPro < 0 ) {
    //cout << "Problem: likeKvsPro = " << likeKvsPro << endl;
    //}
    //if( likeKvsEle < 0 ) {
    //cout << "Problem: likeKvsEle = " << likeKvsEle << endl;
    //}
    //if( likeProvsPi < 0 ) {
    //cout << "Problem: likeProvsPi = " << likeProvsPi << endl;
    //}
    //if( likeProvsEle < 0 ) {
    //cout << "Problem: likeProvsEle = " << likeProvsEle << endl;
    //}
    //if( likePivsEle < 0 ) {
    //cout << "Problem: likePiVsEle = " << likePivsEle << endl;
    //}
    
    nt->column("likeKvsPi",packFloat((float)likeKvsPi,8));
    nt->column("likeKvsPro",packFloat((float)likeKvsPro,8));
    nt->column("likeKvsEle",packFloat((float)likeKvsEle,8));
    nt->column("likeProvsPi",packFloat((float)likeProvsPi,8));
    nt->column("likeProvsEle",packFloat((float)likeProvsEle,8));
    nt->column("likePivsEle",packFloat((float)likePivsEle,8));
    
    // NA
    //delete likeRatios;
  }

  // check selectors
  if (getParmValue("writeSelectorBits")) {
    for (int k=0; k<_nSelectors; k++) {
      HepAList<BtaCandidate> *list = _listOfPidLists[k];
      string &str =  _namesOfPidLists[k];
      
      // NA
      //double runNum = nrun;
      //if (runNum >= 200000.0 || runNum < 9000.0) // MC
      //runNum = _muonTree.getRunNumFromDate(date);
      //makeBDTCutCol(candMomentum, candThetaPi, runNum, ncharge, HepString(str), nt); // muBDT diagnostic
      
      string ntupleName = string("is") + str;
      if (list!=NULL) {
	HepAListIterator<BtaCandidate> iter(*list);
	BtaCandidate *identifiedCand;
	bool found=false;
	while (identifiedCand=iter()) {
	  if (identifiedCand->uid()==c->uid()) { found=true; }
	}
	nt->column(ntupleName,found);
      } else {
	ErrMsg(warning) << " *** ERROR in PacPidCand2Ntuple : PidList \"" <<  str << " not in Event ! " << endmsg;
	nt->column(ntupleName,false);
      }
    }
  }
  
  
}

/////////////////////////////////////////////////////////////////////////////////
float
PacPidCand2Ntuple::consval(const BtaCandidate *cand,PdtPid::PidType hypo,PidSystem::System d)
{
  const BtaPidInfo* pInfo = cand->getMicroAdapter()->getPidInfo();
  const Consistency& theCons = pInfo->consistency(Pdt::lookup(hypo,(int) cand->charge()),d);
  if (theCons.status()==Consistency::OK) {
    return theCons.consistency();
  } else if (theCons.status()==Consistency::noMeasure) {
    return -1;
  } else if (theCons.status()==Consistency::unPhysical) {
    return -2;
  } else if (theCons.status()==Consistency::underFlow) {
    return 0;
  }
  return -4;
}

float
PacPidCand2Ntuple::lhood(const BtaCandidate *cand,PdtPid::PidType hypo,PidSystem::System d)
{
  const BtaPidInfo* pInfo = cand->getMicroAdapter()->getPidInfo();
  const Consistency& theCons = pInfo->consistency(Pdt::lookup(hypo,(int) cand->charge()),d);
  if (theCons.status()==Consistency::OK) {
    return theCons.likelihood();
  } else if (theCons.status()==Consistency::noMeasure) {
    return -1;
  } else if (theCons.status()==Consistency::unPhysical) {
    return -2;
  } else if (theCons.status()==Consistency::underFlow) {
    return 0;
  }
  return -4;
}

void
PacPidCand2Ntuple::fillEmptyCalQual(HepTuple *nt)
{
  nt->column("eraw",(float)-99.0);
  nt->column("ecal",(float)-99.0);
  nt->column("lmom",(float)-99.0);
  nt->column("zmom20",(float)-99.0);
  nt->column("zmom42",(float)-99.0);
  nt->column("s1s9",(float)-99.0);
  nt->column("s9s25",(float)-99.0);
  nt->column("secmom",(float)-99.0);
  nt->column("phicluster",(float)-99.0);
  nt->column("thetacluster",(float)-99.0);
  nt->column("phicmat",(float)-99.0);
  nt->column("ncry",-99);
  nt->column("nbump",-99);
  nt->column("emcstat",-99);
}



void
PacPidCand2Ntuple::addNphotFromGlbLikelihood(const BtaCandidate *cand, HepTuple *ntp)
{
  float cele, cmu, cpi, cka, cpro;
  const BtaPidInfo* pInfo = cand->getMicroAdapter()->getPidInfo();
  //if (pInfo!=NULL) {
  PidSystem::System j = PidSystem::drc;
  cele = pInfo->consistency(Pdt::lookup(PdtPdg::e_minus),(PidSystem::System)j).consistency()
    *pInfo->consistency(Pdt::lookup(PdtPdg::e_minus),(PidSystem::System)j).sign();
  cmu = pInfo->consistency(Pdt::lookup(PdtPdg::mu_minus),(PidSystem::System)j).consistency()
    *pInfo->consistency(Pdt::lookup(PdtPdg::mu_minus),(PidSystem::System)j).sign();
  cpi = pInfo->consistency(Pdt::lookup(PdtPdg::pi_minus),(PidSystem::System)j).consistency()
    *pInfo->consistency(Pdt::lookup(PdtPdg::pi_minus),(PidSystem::System)j).sign();
  cka = pInfo->consistency(Pdt::lookup(PdtPdg::K_minus),(PidSystem::System)j).consistency()
    *pInfo->consistency(Pdt::lookup(PdtPdg::K_minus),(PidSystem::System)j).sign();
  cpro = pInfo->consistency(Pdt::lookup(PdtPdg::p_plus),(PidSystem::System)j).consistency()
    *pInfo->consistency(Pdt::lookup(PdtPdg::p_plus),(PidSystem::System)j).sign();

  //To get the number of photons (and put them in you ntuple):

   int nexPhotons[5];
   int nexPhotFix[5];

   const BtaPidQual* pQual = cand->getMicroAdapter()->getPidQual();
   if (pQual) {
     int k;
     for (k=0;k<5;k++) {
       nexPhotons[k] = pQual->ringNExPhot(k);
       nexPhotFix[k] = nexPhotons[k];
     }
     for (k=PdtPid::muon; k<=PdtPid::proton; k++ ) {
       if ( nexPhotFix[k] > nexPhotFix[k-1]) {
   if (k<PdtPid::proton) {
     nexPhotFix[k]=min((nexPhotFix[k-1]+nexPhotFix[k+1])/2, nexPhotFix[k-1]);
   }
   else {
     nexPhotFix[k]=nexPhotFix[k-1];
   }
       }
     }
     ntp->column("phogele",calNPhot(cele,nexPhotFix[0]));
     ntp->column("phogmu",calNPhot(cmu,nexPhotFix[1]));
     ntp->column("phogpi",calNPhot(cpi,nexPhotFix[2]));
     ntp->column("phogka",calNPhot(cka,nexPhotFix[3]));
     ntp->column("phogpro",calNPhot(cpro,nexPhotFix[4]));
   } else {
     ntp->column("phogele",(int) -1);
     ntp->column("phogmu",(int) -1);
     ntp->column("phogpi",(int) -1);
     ntp->column("phogka",(int) -1);
     ntp->column("phogpro",(int) -1);
   }
}

int
PacPidCand2Ntuple::calNPhot(double sigLvl, double nExp)
{
  nExp=nExp+2;
  float sigLvlTmp = 0;
   int i;
   if ( sigLvl <= 0) {
     i=-1;
     while ( sigLvlTmp < fabs(sigLvl/2) ) {
       i++;
       sigLvlTmp +=exp(-nExp+i*log(nExp)-NumRecipes::gammln(i+1));
     }
   }
   else {
     i=101;
     while ( sigLvlTmp < fabs(sigLvl/2) ) {
       i--;
       sigLvlTmp +=exp(-nExp+i*log(nExp)-NumRecipes::gammln(i+1));
     }
   }

   if ((i-2)>0) i-=2;
   else i=0;

   return i;

}

// NA
/*
void PacPidCand2Ntuple::NNout_p01_t1(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 27;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {2.622376,-0.346356,-0.680478,0.344862,-8.312274
			,-3.556094,0.768033,-3.594284,3.137520,-0.411721
			,14.303416,-0.855723,1.224183,1.760244,-6.226032
			,-0.614359,-1.018419,-0.906303,1.505687,-2.959635
			,-0.505513,-0.323116,4.835679,-0.907459,-0.376395
			,-0.039271,-1.227050,1.983309,12.732339,-0.231886
			,6.217088,-35.419064,4.141907,-12.195803,3.103076
			,-9.047935,-0.050220,-12.009645,-7.831938,-1.830352
			,-14.788012,-9.167117,-0.677867,-0.106471,1.025831
			,0.953186,-51.610020,0.060380,-0.418406,3.455338
			,0.009907,-0.234000,-0.559804,-0.490202,6.101725
			,45.788918,-0.458547,-1.680955,6.423030,-1.587581
			,-5.098583,-0.850734,-10.011388,-0.413536,8.776006
			,-4.415711,2.263852,-15.665514,7.466830,-0.614345
			,-0.428016,0.123363,1.131708,-3.715494,-0.636355
			,-0.735451,12.245057,-0.142554,2.171990,-1.217659
			,21.875633,-1.803506,-0.546894,-1.324546,-0.671429
			,3.679008,1.481578,-0.221416,-4.776346,1.303885
			,-1.356779,4.810832,-4.465165,5.497930,0.922637
			,1.671210,-1.400846,-1.127323,-0.570077,0.280727
			,3.577192,-1.289760,-1.557461,4.212109,-0.499161
			,-0.348781,-0.863671,7.511251,-6.025905,15.805678
			,-0.766300,-3.557365,-5.948771,3.179984,0.629269
			,-1.778380,-11.627464,-0.585102,4.274675,12.166986
			,-8.578757,2.337592,-6.203811,-0.400937,-1.077998
			,-0.941256,8.173085,-1.249825,-0.502284,-0.748484
			,5.737576,-0.625782,-2.565711,-0.044529,-0.994554
			,-1.076960,4.766891,0.295108,1.253181,-8.755034
			,1.565516,12.245161,2.960113,0.989054,0.275806
			,-3.140359,-3.301782,-11.568694,12.191068,2.333565
			,0.537040,0.650621,1.305628,-0.835930,-13.388430
			,0.049336,0.283083,-0.839047,0.829192,0.145601
			,-0.139836,-16.181623,0.035279,24.588383,-0.566478
			,-2.786978,5.388798,3.489367,-20.754118,1.394884
			,-3.807115,-0.160088,4.776802,-9.282251,2.630679
			,0.693081,-7.471769,-0.907215,-0.591970,0.131233
			,-3.065881,14.490877,-0.799109,-0.609207,-8.461136
			,-0.034861,1.202081,-0.765454,-13.220679,1.948892
			,5.900075,-1.472271,-4.528150,-8.647439,-4.411249
			,-3.293719,2.737428,6.687663,-1.303168,-12.528218
			,-10.690815,-5.917727,-5.103814,-0.897648,-0.927213
			,-1.312165,-2.075957,-1.811835,0.596563,-0.788435
			,-1.507283,11.762023,-1.265563,-1.295226,-0.907817
			,-6.755562,-1.991946,-0.392121,-0.801991,6.283622
			,1.273399,-1.633806,-0.710667,-6.283702,5.081472
			,-0.755063,-4.616141,-9.483772,1.878681,-1.789076
			,2.151485,-0.599549,-0.817854,0.234303,-8.281994
			,-2.088239,-0.689707,-0.515878,-3.830216,-0.833296
			,-1.083447,-1.128398,-0.045520};
  float H_Thresh[] = {-2.427187,-31.703320,-0.903970,-4.250201,9.396487
		      ,-1.331846,2.044437,0.878752,-3.589209,-1.014655
		      ,-11.923672,6.107040,-4.104749,-2.562412,1.362763
		      ,-1.283800,-0.834009,-1.018144,-1.371281,7.605524
		      ,-1.252402,-1.235169,-13.046421,-0.780539,-1.345140
		      ,-1.073318,3.568023};
  float H_O_Weight[] = {0.786953,-1.859783,0.049580,2.062247,-5.444054
			,1.365082,-0.871603,0.793754,-2.223660,0.304216
			,1.908344,-4.797303,3.435134,-2.596937,-2.438969
			,-0.440117,0.188516,1.023708,1.846282,-0.909258
			,-0.213889,-0.097433,-1.154391,0.664340,0.597750
			,-0.555613,-0.945491};
  float O_Thresh[] = {0.445608};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p01_t2(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 27;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {0.060011,-3.100123,0.333306,-5.214972,3.635129
			,0.909349,0.616686,0.058561,1.292090,1.138687
			,-6.714002,7.265969,1.032141,2.432893,-0.762659
			,-0.551170,-1.460012,1.842007,-1.991080,3.576232
			,0.518394,1.014913,-5.730227,-6.956591,0.772043
			,11.202765,-2.349404,3.587160,64.272141,-3.915412
			,-5.892306,32.538086,-8.741849,-2.284773,8.905869
			,-0.198289,0.194589,29.211657,12.965218,0.288076
			,-1.819663,19.015871,-5.413591,-5.040876,37.378521
			,8.166199,2.860725,-1.536411,-0.890970,-6.151333
			,-11.455548,-2.174905,2.785183,-29.238401,2.513933
			,45.460541,2.088664,24.989204,-7.134973,-0.538130
			,1.856221,4.096461,0.471033,0.430749,-0.069526
			,-26.866199,15.296914,-2.966318,-35.362556,2.075291
			,1.737696,15.111458,4.393198,2.557533,1.750148
			,1.010232,1.204579,2.603149,1.623047,-13.961191
			,8.028332,0.540918,1.571058,1.001433,0.877003
			,1.712628,3.382681,0.733695,1.433007,0.243197
			,0.297677,3.183094,0.955564,-6.746020,-5.634700
			,4.898459,1.500908,1.911572,0.403808,4.181181
			,1.243830,0.736863,0.464652,2.056954,9.369237
			,0.660145,-1.414760,7.236490,3.213126,0.686895
			,1.327999,-2.310526,6.414876,0.434478,1.224854
			,5.025589,0.214244,0.325065,-10.551279,16.678713
			,-0.425788,-10.393348,2.383971,1.488801,1.137267
			,7.733640,6.657104,15.067707,1.294396,0.688268
			,5.896295,-2.966306,0.994084,5.778960,-6.644119
			,5.671271,-5.919317,4.189513,-3.475550,-1.188257
			,5.537776,3.851329,6.342812,1.876130,2.343299
			,-0.747324,16.239407,10.636902,-18.995253,6.804368
			,5.386840,6.120283,4.780146,7.508995,-13.979384
			,4.137677,2.884032,2.376631,12.588356,3.400898
			,14.514387,-8.855587,-3.871174,0.006869,-1.648399
			,8.630911,11.267335,2.286182,-2.585672,-5.250882
			,-4.428108,-4.427617,5.977787,7.056985,5.063258
			,-2.675153,-0.799260,-0.383606,-0.445533,-1.146737
			,-4.906273,18.135612,-2.839047,-3.696486,2.539525
			,-1.111963,-2.832636,-17.773046,-12.053376,-2.094540
			,1.868875,-0.588621,8.420815,7.583140,-4.117611
			,-0.471090,-3.199038,-0.056094,-0.123262,-16.511988
			,1.114470,5.296480,-0.530819,-8.867537,-1.077126
			,-1.526106,0.496914,-5.360281,-5.146626,-0.560180
			,-0.202389,-2.003912,-7.451159,-0.355801,-5.162149
			,2.811965,-0.363140,1.408603,-0.811548,4.054521
			,-2.047897,-1.130980,-0.881745,-0.016317,-1.080826
			,-1.056707,-1.746895,1.074928,1.273929,10.531051
			,-3.085929,-0.646026,-0.486579,2.030849,0.575972
			,1.594451,-0.873652,-0.986198,2.120945,-8.000321
			,-0.916091,-7.649036,-1.078164};
  float H_Thresh[] = {-4.582959,-17.312899,-3.651601,-12.701927,-16.737064
		      ,-4.174128,-3.235123,-7.195346,-1.264238,-1.495265
		      ,6.677439,-22.521162,-14.155692,3.344338,2.766194
		      ,-4.320561,-4.246005,-16.354303,-7.829980,-13.871190
		      ,-3.284778,-2.213552,-7.560011,-0.945704,-2.862385
		      ,-0.553573,7.932677};
  float H_O_Weight[] = {0.514112,-2.648056,1.533433,-0.966428,-2.062227
			,1.276829,1.517496,-0.765538,1.636671,1.632118
			,-1.156585,-0.616118,3.873174,0.909304,-0.514026
			,1.659836,1.782993,1.204791,0.546600,-1.119453
			,1.512775,1.581817,-0.534175,1.009816,1.543508
			,-1.033285,-1.393946};
  float O_Thresh[] = {0.099190};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}



void PacPidCand2Ntuple::NNout_p01_t3(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 27;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {2.059283,-2.279034,0.750525,4.409858,-12.252503
			,1.495085,1.282756,3.204771,2.555260,-0.869178
			,2.216016,-3.002443,1.988724,2.250317,-0.356577
			,1.176134,0.932896,4.764035,2.175550,18.647181
			,1.113759,0.530939,-0.901009,-2.508357,-14.071325
			,0.173587,0.781407,11.376888,68.834465,4.107863
			,-14.328069,3.571196,7.751875,5.922430,18.781977
			,-8.874579,-3.837852,28.070335,-10.446446,11.195434
			,39.797562,-14.628906,5.703568,5.358143,16.705145
			,11.567447,9.572852,5.477983,2.781979,-6.053731
			,-4.089122,-2.190906,-1.469680,10.944974,5.062540
			,16.188084,4.762970,-15.402204,30.360046,5.041898
			,5.058812,6.087170,3.112065,2.492149,-23.784800
			,18.476852,5.311727,13.999017,-6.404410,5.222691
			,5.156560,6.058814,5.692225,-0.786071,4.268360
			,3.995503,1.817430,4.946680,11.562528,2.143012
			,-9.488118,-0.228813,0.116974,-0.522331,2.048943
			,-7.880021,-0.430789,-0.479294,0.070384,-0.255351
			,-1.773653,0.400358,4.089321,-0.250948,2.660325
			,-12.786448,-0.489941,-0.534666,1.939042,-0.116179
			,11.189636,-0.588656,-0.449916,-3.409268,-0.529030
			,10.155430,-1.387192,0.197977,5.153930,4.988123
			,0.617345,25.650026,-13.075147,2.980844,2.127386
			,10.192122,-9.120069,-3.195542,3.611017,-1.122117
			,5.032712,11.396843,-7.901510,1.971679,1.439377
			,8.582200,5.772889,-1.706970,1.211461,-0.338927
			,-4.385277,-4.016715,1.401840,-1.987938,9.150970
			,2.337822,-13.322887,0.974028,7.543377,-0.906312
			,1.731203,1.419781,4.491830,-1.668704,-0.087095
			,11.033727,-7.703509,2.443171,10.680549,5.858517
			,1.495476,1.443168,3.552843,2.856063,-11.259571
			,0.911296,0.440007,-1.125861,2.663589,1.375558
			,-0.806755,-3.834949,-1.850593,18.922165,-2.413683
			,21.516386,7.073940,-2.437956,-2.498109,-0.468224
			,4.568557,-0.733938,3.363798,1.848563,-1.936189
			,-2.187085,3.985531,-2.540545,-2.604373,-0.216645
			,-1.772782,3.958527,-2.429211,-2.007285,0.166357
			,-3.367847,1.594462,-0.466228,20.740377,-1.620741
			,-1.491090,-0.883755,1.291846,7.810455,-1.363137
			,-1.206023,-1.616483,0.989913,0.289979,2.378990
			,2.130044,-1.560603,2.191647,-1.342966,-1.139491
			,-1.031868,-2.782471,-1.579897,0.876763,-1.074140
			,-0.650993,0.364675,0.386188,-2.147277,0.144406
			,-6.826795,-0.699891,0.285481,-0.689227,2.626888
			,6.067005,-0.810351,-0.820369,-0.449057,0.394306
			,0.323277,-0.624304,0.375692,-0.689208,4.458147
			,-0.603669,-0.804070,-0.761588,-1.746977,-0.662780
			,-5.570108,-0.737484,-0.536331,1.060023,-0.146580
			,-4.536347,0.092302,-0.790589};
  float H_Thresh[] = {-7.554464,-17.187012,-3.427471,-29.136574,-5.242956
		      ,-5.529793,-4.713884,-12.317471,-0.660140,0.568945
		      ,-11.300864,-2.865648,-7.606428,-24.632904,2.584152
		      ,-4.683124,-4.262850,-9.280909,-8.471108,-12.894922
		      ,-3.721127,-2.517323,1.690473,0.380008,-5.323619
		      ,-0.363383,-9.227159};
  float H_O_Weight[] = {0.478993,-15.881866,0.584324,-0.405483,-0.683075
			,0.504306,0.492154,0.476084,0.644623,1.781121
			,-4.031495,-2.414008,0.495229,0.743645,7.555398
			,0.494102,0.522801,0.417513,0.434370,-0.732962
			,0.597557,0.710281,1.718444,2.399987,5.261827
			,1.676175,-0.636640};
  float O_Thresh[] = {0.209552};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p02_03_t1(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 27;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {-1.051430,3.737711,-1.539590,5.570607,0.535346
			,-3.720299,-0.872152,5.202000,3.887331,0.381193
			,7.396149,0.039511,2.771572,5.006593,-3.796318
			,0.516316,-2.091779,6.167094,2.490716,1.578823
			,-2.722751,0.665716,12.986141,0.992560,-2.222446
			,3.228211,-6.289504,3.933572,57.474678,2.419127
			,49.144840,-22.885622,1.775856,-60.735020,1.465299
			,8.802070,-0.310008,16.545218,-13.303639,59.409443
			,3.824812,-23.884365,-9.007070,0.241278,-6.573823
			,3.113803,-35.105045,-6.474595,-9.150476,-13.515204
			,-11.934220,1.486088,-20.155405,-11.541180,9.368827
			,41.191963,2.635324,12.629782,17.421860,2.680463
			,33.922745,9.633759,16.628998,-0.197880,8.828441
			,1.670542,-38.948811,-33.342346,19.461128,-2.188923
			,0.080104,11.581970,18.641348,-6.172005,-0.261771
			,-2.871108,-33.278572,-1.294957,2.144067,-3.727947
			,35.413929,-1.071277,0.184073,-0.717022,1.098694
			,-1.878933,-0.229736,1.831933,-2.171383,-21.652441
			,-0.720281,7.230062,2.684494,1.435198,1.078786
			,9.498843,3.070882,-0.647620,-4.740373,2.975658
			,1.522000,0.809650,3.182245,2.960005,3.899916
			,-0.121231,7.024380,3.866066,7.591809,15.891614
			,1.198458,23.572935,-3.976592,1.991780,8.693024
			,-13.258075,-5.204395,-0.776575,2.994528,-1.854997
			,3.156298,-14.778959,8.391713,-0.974361,-1.576927
			,0.197812,4.946000,18.072361,1.775853,-0.885591
			,28.395164,0.690029,0.614931,2.510830,-7.280754
			,1.893029,-1.551301,0.430737,17.447201,2.907760
			,0.586728,-3.252469,1.311677,-9.649301,-0.293432
			,-14.087135,-1.246141,10.674542,11.652440,-8.072358
			,2.968821,1.608604,18.443407,5.139914,4.820967
			,5.124738,2.067262,33.902153,5.410214,0.905692
			,-4.026092,-11.234615,-2.938491,14.614432,1.847521
			,9.445656,-10.510033,0.190842,2.989273,1.375551
			,-5.384264,0.032666,-13.545283,0.887277,9.487100
			,-12.272650,2.175329,-9.545884,1.622378,-0.220949
			,-6.543411,25.718103,-11.016820,-7.253171,0.107923
			,-13.757201,1.787170,2.063877,1.947335,-3.709504
			,0.926954,-3.134851,0.302520,1.565294,-3.423791
			,-1.456142,13.715613,-12.256047,-0.657221,-4.343851
			,0.338478,3.545903,2.003224,7.951712,-0.408865
			,-0.013679,0.503497,-2.024894,-1.633951,-1.854866
			,-0.206684,-2.856582,-1.011764,-2.588732,2.217410
			,2.121784,-5.717931,-0.873934,-0.542062,10.732372
			,-0.868980,-2.864122,-2.777858,-11.106639,-4.146719
			,-1.162242,-2.410786,1.955043,-0.599553,3.576746
			,-15.244066,-0.175969,-1.461884,2.955413,-1.960287
			,0.440850,-1.137473,0.448436,-0.411907,0.190827
			,-1.222120,-1.994886,1.797129};
  float H_Thresh[] = {-3.338475,-29.736528,-1.536733,-44.087700,2.196026
		      ,-0.816679,-2.678292,-5.529889,15.272588,-0.757486
		      ,-1.808241,-2.297474,-16.645435,-0.282031,-4.271160
		      ,-0.455491,-0.321774,-15.795039,-7.304166,-17.766346
		      ,-0.468779,-1.388039,-30.362162,-1.858687,-1.245273
		      ,-3.956205,-1.396482};
  float H_O_Weight[] = {4.802175,-1.776136,1.303460,0.435496,-2.009027
			,1.775442,-0.898594,2.124969,-9.601557,-0.744082
			,-0.867924,-1.082995,-6.585219,-3.742794,-1.713243
			,-2.211260,1.079900,9.887923,1.036745,-0.447121
			,-1.861805,-1.812412,-0.973317,-2.049599,1.144779
			,-1.110770,-6.016970};
  float O_Thresh[] = {0.234051};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p02_03_t2(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 27;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {4.700779,3.788246,-1.958382,-4.868440,6.224218
			,0.069424,0.095150,-1.217474,4.461841,-2.300685
			,-2.415920,-3.293389,7.952304,-21.879404,5.858496
			,-0.479132,-3.393609,8.770898,7.658543,-3.843240
			,-2.310434,-3.654654,1.161171,-0.007903,16.711744
			,-4.032805,-8.003207,61.047516,72.734955,11.641995
			,0.445839,0.177128,-10.261880,-1.400116,-5.032188
			,-24.784567,-0.589127,42.834026,-6.409650,76.561676
			,-30.534512,-1.888550,-0.206666,17.171957,55.663593
			,13.895306,23.915438,2.640229,0.488883,-5.158698
			,-5.291561,-22.667364,7.976601,-19.071920,0.332458
			,52.513737,23.842226,-1.578716,-22.554731,3.333412
			,-12.432005,32.018181,-10.814531,10.318067,21.377411
			,-6.499785,-58.565414,-11.889613,-16.927559,10.913691
			,-4.969024,3.444570,4.902532,6.862912,12.296255
			,11.506325,-11.915471,6.232499,-9.530120,-20.061659
			,66.404785,-5.796897,-2.144321,4.352976,1.512231
			,1.790486,3.588929,8.647932,1.021412,-3.313155
			,0.218648,4.766817,2.024085,-1.000715,14.784370
			,2.199343,2.407416,-6.409677,5.343894,0.510323
			,3.057294,0.144298,0.058310,0.599568,-0.253769
			,-1.552768,6.099288,3.239326,13.421889,15.583735
			,-2.430351,-4.246175,6.640174,2.140415,3.397738
			,3.006432,5.433804,4.388855,-6.465994,19.853081
			,13.571531,7.902466,-2.486516,-7.846948,-4.358871
			,29.334768,21.028629,-8.987032,-10.990654,1.298282
			,16.382952,-0.530743,0.873128,4.333436,-14.038733
			,-1.059915,-1.848211,-0.004269,-5.969730,12.026539
			,1.005515,-4.396564,7.310098,-2.386615,2.858741
			,-8.520425,5.567416,21.047613,2.845961,-3.890684
			,2.929632,21.354269,30.624924,16.411491,-3.628963
			,-6.772887,4.366673,8.953755,2.943284,12.068361
			,3.533985,-1.673323,8.750462,8.443762,-3.308982
			,2.691503,-13.415664,-8.916556,1.932439,-0.534179
			,1.121951,-0.396024,16.283569,18.850006,10.764607
			,8.268440,-2.851423,-7.223288,-13.119864,-9.188894
			,11.116436,-5.599998,-0.308092,0.292662,2.143317
			,-2.728349,-12.010427,-0.677050,-6.078368,-2.384439
			,-2.386265,3.293910,-6.314860,-1.809324,-5.082350
			,0.243563,-1.176241,-0.512097,-4.118711,1.807500
			,-0.850108,-0.512415,-2.109417,3.144979,3.271348
			,10.834007,1.522722,4.281784,2.315349,2.839577
			,-1.354145,-0.326255,-0.993570,-21.680481,-2.126068
			,14.534749,-2.594288,-2.552207,0.004716,9.844038
			,-1.413753,-8.605643,-3.134916,-0.278213,-8.147614
			,-4.690478,9.867951,-5.820439,0.514922,-5.925007
			,-0.889716,5.762610,5.977902,11.151290,-0.010973
			,-6.779505,20.619097,-4.700175,-8.201453,-1.011511
			,-0.176248,-1.528455,4.661776};
  float H_Thresh[] = {-14.409073,-26.980417,-5.601151,-4.938764,-4.991288
		      ,5.466325,-2.284545,-8.203089,3.263127,-0.387593
		      ,-19.857214,-18.353365,-25.465406,-4.274438,-0.734866
		      ,-4.765926,-16.176027,-50.686485,-32.083771,5.651153
		      ,-14.810565,-0.927249,-9.516143,-0.660593,2.603766
		      ,-3.487276,-10.131201};
  float H_O_Weight[] = {0.290542,-1.815312,0.351662,2.842907,-3.980795
			,4.126702,-0.571612,0.516705,1.695250,0.782196
			,-0.982176,-6.900722,-1.359651,1.108354,-1.416938
			,0.693282,4.568102,0.359741,0.765549,-2.640126
			,0.616803,1.076582,-2.190164,1.921755,0.769589
			,-0.686399,-1.225401};
  float O_Thresh[] = {0.003098};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p02_03_t3(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 27;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {-3.976975,0.286208,1.414621,0.881895,8.624770
			,3.215213,-9.809915,4.208948,9.435684,-3.334914
			,-0.217901,-3.287714,2.495372,1.015826,9.705297
			,1.696160,4.952572,8.480878,0.235819,-2.453387
			,-12.757239,-0.537479,7.784166,3.760798,-4.295581
			,3.338349,-7.758401,33.332981,77.026527,10.184397
			,17.182602,23.916393,-1.439763,-0.534412,7.792762
			,-3.686719,0.909146,73.102661,-0.000484,13.075744
			,8.206031,-0.450077,8.038611,8.977380,60.468052
			,19.147758,13.793731,5.133171,-3.513934,4.205285
			,-9.639578,-1.768350,13.574569,-12.998301,3.304946
			,49.131756,17.196146,-6.965580,-30.092966,7.484136
			,9.317901,11.887822,1.405873,9.671651,-38.828182
			,-10.640529,-7.125463,13.189467,-6.099749,16.114943
			,21.371208,14.797401,13.907582,12.968040,15.593643
			,11.704529,-15.867204,5.563485,-6.406273,-29.430107
			,64.589958,7.221457,-0.145973,-0.359031,2.192799
			,2.270631,-0.333919,2.550459,4.830898,-11.321748
			,-1.382959,1.581873,9.382557,1.510066,0.479272
			,-5.680426,1.274066,-0.105599,0.686946,5.679493
			,4.413431,11.708328,-1.267854,3.012748,-2.721281
			,-3.502964,1.483041,3.055446,16.321171,6.135687
			,4.275986,-11.187180,28.287266,3.447702,5.471018
			,4.084300,0.099503,-13.619280,-4.926514,15.632079
			,11.612224,3.531461,-5.094695,1.638765,-0.891317
			,27.304804,7.034790,-5.463281,-7.587438,-7.594675
			,3.533863,7.037813,5.427835,11.207474,-12.321970
			,-4.728426,-2.375612,6.609157,14.914864,27.469299
			,4.459276,0.395938,17.043545,-3.516221,-8.732965
			,8.192495,3.369210,19.797285,12.817239,-6.610083
			,6.110439,8.396561,22.645443,12.098929,-7.344576
			,12.118646,-0.313750,-5.735699,3.130340,-0.260177
			,12.747986,-6.490775,8.966516,10.305584,-1.546946
			,-35.555725,10.725844,-10.447294,-1.345410,-6.846353
			,7.456450,0.369152,7.182516,10.465590,-15.949421
			,-4.598981,3.757485,-5.066532,-13.725561,-7.851324
			,-6.092204,13.877939,-3.286835,-3.831798,2.809387
			,-5.699144,11.286193,2.701882,2.483039,-3.371398
			,0.186870,0.827082,-1.306895,-5.295185,-2.331592
			,5.787676,7.814324,12.270384,3.036394,-1.508933
			,-3.846463,3.490126,2.692736,0.980178,0.882754
			,-0.688937,-1.255975,0.920401,0.881939,-7.620941
			,2.421618,3.959882,-5.742014,-27.181713,0.763470
			,17.099323,5.891332,-0.785301,-4.377535,-0.665729
			,-2.155199,-4.435039,-13.425159,3.293980,-2.394776
			,8.884869,-0.020040,-4.477934,-7.810154,1.712162
			,-6.566752,-1.212231,1.288499,-0.648772,1.724545
			,9.402183,3.133296,0.910410,-0.693182,-6.959942
			,-1.658680,0.999143,1.144052};
  float H_Thresh[] = {-25.553637,-23.932095,-8.552048,10.787523,-31.330830
		      ,-0.722712,-1.538590,-21.367289,-11.192607,-0.100069
		      ,-8.011582,-12.319598,-9.730624,-12.595702,2.139218
		      ,-5.986263,-3.494727,-33.542580,-18.106792,-14.237093
		      ,-3.658335,1.932811,-7.366811,3.299145,5.725353
		      ,-16.316902,-11.722410};
  float H_O_Weight[] = {0.509234,-1.897415,1.503419,-0.903647,-1.206583
			,7.656332,0.988745,0.934083,3.515934,1.069366
			,-1.366220,-0.437028,-2.473332,0.509856,-6.330029
			,0.414327,0.426290,0.453650,-3.503709,-1.023448
			,0.631957,1.245540,-0.652931,0.919505,1.540553
			,-1.662005,-2.022648};
  float O_Thresh[] = {-0.048117};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p04_06_t1(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 27;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {7.647419,2.776640,-3.313729,17.305798,-1.504545
			,3.055646,15.023169,-4.360277,6.117315,-1.103525
			,22.056192,-3.972637,10.381137,12.988437,4.275847
			,0.779253,-7.120701,-5.513398,-7.682601,5.527375
			,1.293314,3.831113,-3.723697,4.117398,2.860862
			,3.244184,-3.100491,14.229985,55.077579,3.465662
			,-14.453705,-38.863773,-9.124185,-31.839663,3.326692
			,-2.229332,6.738578,-4.125888,7.182136,6.981361
			,23.795895,-3.800690,5.652484,2.138645,2.718099
			,3.712069,1.239048,3.493181,3.609848,-13.842639
			,-7.600157,21.424902,12.323227,8.538328,2.679967
			,32.794052,9.403129,5.084379,-60.489761,-9.952276
			,-35.653881,9.373093,9.989541,6.456274,0.243581
			,6.027714,-0.460997,-38.680637,-10.357758,7.358878
			,9.532847,9.664682,1.592089,-21.741152,10.078996
			,11.568951,-8.693162,-13.235929,2.385986,-20.635052
			,41.786537,-0.230950,-7.235461,-0.170136,-14.306744
			,3.935179,3.720631,6.599820,0.545775,-2.073852
			,-2.075480,0.818647,-2.274340,3.946939,-1.893811
			,3.672233,-1.949818,3.199135,2.146569,6.239397
			,4.008390,-1.043654,1.568762,2.260434,5.168362
			,4.537064,0.231167,7.358438,-10.448877,-0.833808
			,-0.302939,-2.497015,26.621483,-6.384554,9.480236
			,-0.122456,0.628318,0.992186,-8.184154,1.842599
			,1.114510,19.329872,-3.469657,1.221171,-3.525817
			,-0.440826,2.764849,3.717029,2.257154,1.595915
			,3.622661,11.537406,-3.253274,-1.162666,0.714588
			,-1.921187,2.440172,-2.576428,-1.211743,12.899721
			,1.668848,-16.024376,-2.668327,4.611773,0.864581
			,-2.204620,1.469202,0.753472,30.207403,5.999987
			,0.459445,-7.156719,-4.871691,10.324296,0.616886
			,0.806922,1.688537,-0.626095,13.063306,-12.549866
			,-2.799500,-0.240057,-3.104431,4.005985,2.843044
			,-3.395843,26.414665,-3.371767,7.526113,3.413617
			,6.847472,0.663769,6.121067,1.283003,10.515012
			,4.442589,-16.016872,0.703757,3.552856,3.866471
			,-3.701737,-5.068143,3.242543,-1.833290,7.252023
			,-13.001966,7.174603,15.698257,-5.491320,-8.582149
			,3.178304,-1.089735,1.018400,0.401596,1.694960
			,9.087839,-0.958001,6.730281,-1.065762,-1.465438
			,-1.991649,1.709540,-0.070920,1.776012,-0.759012
			,-0.373478,-1.118664,5.925251,-0.569053,0.241866
			,-5.598356,-2.618413,-2.602058,-4.526628,2.145278
			,0.827435,-1.286487,-5.641528,2.648863,7.872199
			,-7.480210,-1.155804,0.538229,3.504055,4.071605
			,-0.091372,-16.570177,-1.451725,-8.198981,7.179627
			,-3.343099,1.139162,6.176603,4.019566,-2.062552
			,-0.695353,2.341938,0.496878,-7.667631,-3.564581
			,1.001744,3.570620,2.215865};
  float H_Thresh[] = {4.222156,-8.976188,-4.416965,-6.395539,-16.349157
		      ,0.473509,-14.650753,-5.820154,-14.111280,-2.781895
		      ,8.335790,-1.831182,-13.092547,-34.461483,3.468302
		      ,-4.345779,-5.633445,-5.898817,-12.095002,-3.594775
		      ,-7.778479,-4.017400,3.967710,-7.297278,-2.693477
		      ,-7.925068,-14.543273};
  float H_O_Weight[] = {-5.057774,-1.684997,0.943479,2.026345,-0.289386
			,-5.563109,-0.369860,1.110152,0.922012,0.545749
			,0.805836,0.660854,4.144252,-1.652393,-0.689539
			,0.629048,1.361484,1.186382,1.588225,-0.471942
			,0.791290,1.132771,-2.683344,-0.545952,-0.793751
			,-0.644132,-11.217312};
  float O_Thresh[] = {0.182899};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p04_06_t2(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 27;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {10.518368,-5.258147,9.976996,-26.675854,24.614170
			,-1.170944,2.360745,-1.100834,-0.856301,2.018574
			,14.316254,1.160844,3.803314,2.532972,4.363242
			,-9.296311,-7.394849,2.617290,-1.292509,-1.797855
			,0.174030,-0.261703,5.060271,-0.441469,-9.566301
			,-10.001276,2.836937,-19.519604,53.478390,-2.538461
			,-9.746003,-18.753983,-17.948971,-19.395172,24.969494
			,-13.152670,-25.286278,10.760536,-1.563255,55.523468
			,-37.833488,-35.248493,14.242386,31.379284,21.529114
			,29.205290,-21.586277,31.669809,-19.737297,-11.978035
			,-21.527130,-2.044139,23.998539,6.146924,5.726713
			,53.403355,-2.859573,-59.624687,24.907225,7.234847
			,11.183745,-16.663780,0.376528,14.707035,10.290067
			,64.406891,13.251337,6.450013,26.603191,-18.362843
			,-57.279465,-81.453308,-35.038132,32.211124,-53.950661
			,12.056500,-3.772716,9.865427,11.980204,-32.012215
			,18.373829,2.261981,-8.729837,-3.128893,13.586120
			,0.693352,1.542658,-0.191523,-5.997388,2.614023
			,-3.215875,-2.976595,1.860124,-4.643592,0.059332
			,-6.388258,0.378147,1.213115,-1.114085,-4.302296
			,1.146157,0.103193,0.738093,5.867854,1.785180
			,0.107406,2.125002,2.028430,18.308964,0.103942
			,6.439246,6.230838,4.055483,13.605909,15.530104
			,-3.075600,13.619011,15.719556,-0.486796,8.437402
			,-3.121372,22.808020,15.659815,-9.274564,7.603892
			,9.065927,-0.355179,17.519321,9.155167,16.776781
			,8.534307,18.065540,-11.837345,-10.187612,-2.368269
			,30.362061,-2.550386,3.546347,-3.168740,-8.898213
			,-0.020827,4.993957,5.340093,1.747989,5.574124
			,-0.068361,4.648885,9.923824,13.852451,-17.220385
			,-4.015054,2.974653,3.909542,-4.031349,-1.579038
			,0.141742,8.950109,6.187891,0.602257,0.532021
			,-9.320897,-5.206374,0.222939,6.852871,-4.961308
			,0.074560,4.760335,7.717213,5.054600,-2.371832
			,26.323582,10.906875,7.595564,-3.507409,1.975392
			,5.537010,14.980271,-7.553007,3.965679,7.020234
			,-1.641871,17.932623,5.729939,10.218737,-3.474675
			,9.678618,7.266812,1.817647,1.040022,0.867516
			,10.828560,-2.725621,-5.381927,14.082308,-2.392204
			,-2.822238,1.973090,0.654791,-1.276165,13.463037
			,0.128349,2.427580,-2.224599,-2.556541,-0.194240
			,-0.953110,0.068748,0.459519,24.609112,-2.432652
			,-1.142776,-1.500530,-3.055646,7.847236,0.301433
			,-0.054048,4.255130,-1.590971,5.047879,-1.812040
			,15.620035,-2.579449,-8.765755,1.569721,-2.684007
			,-18.038025,3.707329,-3.962686,5.367682,-2.548885
			,-4.640941,6.619658,-10.025311,-3.852302,-1.267747
			,-4.329815,-5.763411,-2.193319,-2.655575,-3.204348
			,-13.425458,4.549561,0.394046};
  float H_Thresh[] = {-30.409254,-22.163227,-10.842057,0.867717,-34.654316
		      ,-10.231858,-6.839666,-2.492036,-15.894042,1.488957
		      ,-21.180571,-9.379153,-16.810724,-17.681669,-6.423572
		      ,0.181871,1.059586,-5.355598,-0.608128,-38.610352
		      ,-3.715908,-17.900520,-5.687150,-13.400617,8.197948
		      ,2.589552,-2.986198};
  float H_O_Weight[] = {-2.238692,-7.145565,4.977997,0.452758,-0.818577
			,-0.945719,-5.178027,0.344166,-0.427852,-0.553701
			,1.017888,0.406181,-1.437476,-0.568988,-0.386425
			,6.481688,1.273310,1.044639,1.220211,-0.440926
			,2.529138,-0.838232,-0.559334,-0.696214,-2.060392
			,2.072581,-1.390766};
  float O_Thresh[] = {0.162147};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p04_06_t3(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 27;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {10.873427,-14.539496,-3.021569,-11.706572,20.297688
			,0.769345,4.277472,4.762082,-5.313585,-6.503864
			,-0.348007,-1.426157,4.583822,0.525045,-1.961784
			,-2.522666,-1.383923,-7.174457,1.239825,-4.436516
			,6.332713,4.175361,-4.738801,0.420251,6.955471
			,-4.072864,3.022919,-25.842121,33.139671,21.282499
			,-12.735579,-39.630104,14.922155,10.009824,50.810097
			,-0.125658,-15.873994,63.175980,18.861038,20.576038
			,10.115094,-53.634800,22.511377,46.183460,50.087009
			,49.545006,8.582643,33.159386,-14.511734,8.348323
			,8.238430,9.495115,36.425144,-35.255829,5.596075
			,62.990936,-31.042429,-79.458847,23.220779,-3.185385
			,9.680047,-76.804611,23.415352,1.495229,20.590797
			,-13.598745,3.440242,15.526725,-0.824236,-40.310730
			,-89.879433,-78.768127,-82.436127,-0.367826,-88.814926
			,-8.013309,-0.408537,13.314787,12.459519,-59.163269
			,1.577556,2.515153,-7.250366,2.002997,6.367617
			,0.508216,-0.557158,-1.395316,-0.338470,1.088978
			,-3.354341,-7.533962,0.307582,-1.514161,1.672532
			,0.929226,2.450845,5.416008,3.106148,-0.774043
			,-2.811515,-1.405180,4.539915,1.673968,-0.891398
			,-3.249414,2.643153,7.751611,12.273778,-2.018523
			,-6.539083,13.045151,33.602642,-1.898833,-0.622334
			,-7.929263,-6.945366,21.646221,-2.033526,-2.861197
			,-12.422793,0.751346,39.629555,-6.074426,12.083822
			,16.750566,1.816008,-2.469772,7.693617,25.201548
			,-4.233734,-1.204559,2.643348,-13.970142,9.694752
			,47.296246,0.431085,-3.782957,-1.879682,-6.082032
			,4.167597,3.638556,1.165090,-3.165350,11.475713
			,4.748649,5.637335,4.507379,1.253583,18.019537
			,-3.133969,64.934586,53.774323,16.201805,-2.293644
			,2.073631,2.133940,-3.730639,1.250197,4.974655
			,-3.595292,-4.786320,-32.273251,0.781525,-5.289290
			,19.101789,23.502684,-2.317063,-1.551976,-11.127206
			,-2.353763,22.047642,12.042022,-5.255700,-5.767258
			,-0.132190,16.157085,-5.272397,-55.524094,-17.512484
			,-13.398760,-2.504638,4.316283,34.616756,-3.553299
			,-2.310932,1.580272,-7.726980,11.276610,5.823737
			,10.412008,-1.356338,-1.006603,30.615469,1.742137
			,2.029044,-1.101512,0.331472,1.702606,2.485987
			,2.185036,0.517226,0.966195,-4.038646,-1.049021
			,7.262940,17.310589,-3.068621,-3.596426,-0.435644
			,-1.613549,-1.455858,-0.117862,5.182494,1.002608
			,-2.176620,5.863804,-6.544684,9.920363,10.269640
			,8.344255,-0.164650,3.240958,4.201968,10.095146
			,-3.916811,10.384247,-0.350149,12.091197,-6.502905
			,-1.725020,10.182952,-5.944910,-13.321486,-1.514495
			,2.736452,0.514489,2.498733,5.376270,2.826505
			,-0.822971,23.464025,3.737623};
  float H_Thresh[] = {-25.488752,-12.465183,-4.287776,-22.159544,-69.992706
		      ,-4.683199,-7.482095,3.820189,-6.016065,-24.752647
		      ,-24.537828,-2.801567,-6.720068,-0.789123,-33.117764
		      ,-4.725405,-15.161522,-25.924175,-2.203373,-1.735872
		      ,-7.605275,-31.161697,-2.653242,-5.296488,-10.628063
		      ,-12.097428,-14.321425};
  float H_O_Weight[] = {-1.053622,-5.699197,0.753377,0.817478,-0.564779
			,-1.462100,0.435461,0.463512,0.539111,-4.545678
			,-1.180139,1.986080,-1.306598,-2.733591,-0.415149
			,0.645910,0.570144,0.362571,3.127333,6.691851
			,0.643333,-0.629039,0.894001,0.060123,0.970591
			,0.527369,-1.466021};
  float O_Thresh[] = {0.170734};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p11_13_t1(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 27;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {4.251773,1.071867,3.246481,-16.150843,7.483552
			,-1.417884,8.903205,-0.889573,-2.893770,-0.743038
			,-13.363668,-10.115088,7.207873,12.416153,4.781842
			,-1.092394,-4.945232,-5.992517,1.079020,-5.335440
			,-1.186292,1.390261,1.953962,-1.974153,-2.181745
			,-5.625131,8.909069,-10.067166,34.261917,1.801860
			,-1.827337,-13.829367,-11.996469,-9.332087,-10.322402
			,10.631377,-9.446750,7.206497,9.006998,60.104080
			,-43.384491,-11.631485,10.948232,11.827361,10.483948
			,-10.663581,-13.253948,10.488242,-10.669033,-11.353561
			,-12.921470,-7.244703,9.337849,1.017721,3.620646
			,9.277870,-3.548563,-9.429460,6.980902,4.977860
			,5.656138,3.673776,-5.784054,2.042804,7.599953
			,-8.141398,5.785372,0.867857,5.428246,-1.820991
			,-6.092402,-7.286266,4.035685,5.860745,-4.441041
			,3.497503,4.363162,6.309521,2.261640,-9.162250
			,5.768519,2.213994,-5.889475,-6.319467,4.678794
			,3.150579,-6.726159,1.769121,-0.922458,-2.593227
			,-0.425327,7.844189,0.918247,-9.926564,3.969331
			,2.357299,0.558055,-8.428244,-2.797957,1.885683
			,0.231693,-1.753323,2.078979,2.371230,-1.198862
			,0.780696,-1.647091,18.485353,3.497899,-2.928108
			,-5.115322,8.552238,3.763454,4.606990,2.219959
			,2.873684,2.272832,3.562647,-5.927258,0.484959
			,-8.758559,4.959258,3.900704,3.590003,1.373602
			,2.932852,2.673467,2.861807,2.608151,2.960573
			,3.165215,2.786712,-0.336239,5.104011,-0.766858
			,6.332221,0.515002,-5.916963,-10.572214,-3.459361
			,4.121755,11.931316,0.932976,-2.920390,1.089118
			,4.873620,-2.874457,-0.221257,-5.888495,7.252131
			,-4.566091,-0.762500,-2.790528,1.259296,-8.203853
			,-3.043900,1.649829,2.180043,-1.241957,-0.897980
			,-5.481543,-12.589396,-10.544699,1.433923,5.910949
			,8.735081,7.434902,8.526247,-19.645988,4.941188
			,6.248499,6.313047,1.982747,14.084028,-0.317726
			,-6.069146,-11.719019,3.884647,6.341482,6.217166
			,-1.122978,14.964791,4.931365,-1.707659,-3.385193
			,9.434432,10.958682,7.720162,-5.578452,0.158920
			,-1.629397,-0.661962,-4.721695,1.543437,-1.149131
			,1.790373,-1.353044,-1.082860,-1.474457,3.373431
			,0.910716,1.378188,4.245530,0.650017,-2.684091
			,-3.003381,-1.137442,-1.179406,0.770415,-1.274914
			,-1.069016,-0.978711,-1.993336,-1.888240,-2.382761
			,0.493739,-0.075891,-0.742469,3.481589,-2.313084
			,8.553950,-6.317563,1.304884,-2.485517,-0.848336
			,-1.712222,-1.862944,-2.625910,2.410107,0.509088
			,0.296288,-0.969054,-0.111418,-1.215948,-1.203091
			,-7.708010,-1.081438,-0.957418,-0.939890,-2.971243
			,0.513466,-1.727337,-0.796082};
  float H_Thresh[] = {-5.487309,-4.436758,0.714064,-3.184883,-15.313585
		      ,1.184713,-6.370799,-1.256285,-2.640965,-2.979884
		      ,-5.436469,-3.260153,-4.770226,-5.778254,-6.617346
		      ,-3.165884,-1.270606,-2.950084,-2.523870,3.387232
		      ,-2.431225,-3.212064,-3.418547,0.162968,-2.989942
		      ,-3.055797,-10.679008};
  float H_O_Weight[] = {-1.458606,-10.476552,1.558367,3.014992,-0.427610
			,-1.463846,-0.571851,-1.741040,0.968767,-2.032391
			,-5.853713,5.184241,-1.035765,-0.381534,-1.311383
			,1.219275,2.105433,1.324164,-1.586549,-1.487149
			,0.681569,-1.637906,-1.488508,-1.249044,-3.293720
			,1.276096,-0.320367};
  float O_Thresh[] = {0.138727};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p11_13_t2(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 27;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {-3.498198,0.834567,-5.733747,-0.498988,2.028677
			,-1.470873,0.086060,3.533501,-7.388380,8.516199
			,11.044760,6.002933,4.966320,1.463009,8.017811
			,6.964903,-11.359413,-10.780630,-1.090845,-14.063545
			,5.540856,22.416897,5.736620,0.062569,3.013125
			,-21.595701,1.635702,-42.857769,21.196442,13.901727
			,24.857183,-18.654778,-37.964863,-39.762756,-25.461048
			,18.466232,26.267277,12.653717,11.814047,48.640461
			,-52.158314,-38.313282,34.048878,40.730221,12.764567
			,-36.518497,3.540216,31.757172,-25.463987,-5.645867
			,-36.469296,-12.071112,27.100845,-41.037399,42.310425
			,26.098379,-10.800766,-98.528709,20.883520,54.996132
			,57.274540,27.918028,-26.009169,14.125729,18.048010
			,-122.058578,37.618866,61.025661,47.859032,-61.204258
			,-51.309475,-42.329063,52.437496,12.792183,-81.200134
			,24.956324,29.868883,45.978306,10.453690,-59.370728
			,58.631359,0.251579,-3.913804,-6.748681,0.338159
			,2.470271,1.736284,1.397320,1.652779,-4.501653
			,-10.293432,-8.490906,4.876348,-11.258936,-3.961212
			,0.158670,-2.211490,-28.226070,3.638247,1.799931
			,2.649294,-2.637757,1.987990,3.646040,3.554461
			,5.931115,6.795036,1.083025,5.503755,-6.811599
			,6.182138,10.175610,1.831083,-0.880411,0.522061
			,1.848512,28.462446,-4.274661,-6.910533,-1.160608
			,-6.085409,3.117068,4.519095,19.791685,3.229485
			,37.609627,-1.109030,-18.966841,15.950263,6.595416
			,-1.652150,-1.407349,2.994419,10.644972,1.950137
			,7.714367,14.648829,-2.724665,1.258698,-19.678711
			,-0.746938,1.239566,4.723058,-6.910842,-10.426063
			,-5.641014,16.890230,6.902866,4.403218,12.259186
			,0.957829,-7.196922,-6.029761,-0.856604,0.480004
			,2.627630,-2.251314,-2.290157,-1.312859,8.836378
			,-8.684787,3.746673,0.691361,-2.242930,2.383703
			,23.215685,8.638739,-1.591493,-2.028348,-5.037926
			,18.819338,-0.092543,-2.464194,-23.736925,-2.690840
			,-0.273370,-6.733482,10.149232,18.187675,21.511713
			,-1.695574,13.799703,5.351448,-2.304530,-4.698022
			,-0.862197,2.046785,27.907106,-2.894871,-2.245411
			,7.804303,1.222790,0.921965,12.886399,1.049142
			,0.977953,1.301740,4.753991,-0.611819,6.879083
			,0.827556,1.138841,4.722095,1.498317,-2.671881
			,-2.268554,4.645917,1.311075,22.682962,-3.211045
			,7.201231,-2.291462,1.722142,0.819182,-2.926863
			,0.967445,-2.897535,11.256960,-3.956401,2.634122
			,-10.197317,-1.068609,-0.578130,2.384021,-12.988188
			,-5.295470,-2.860793,6.058325,-0.174712,-5.590160
			,1.588967,-0.489316,-3.297811,-15.364578,-0.265245
			,11.378469,-5.773838,10.227625,-3.395113,2.300832
			,3.023327,-0.044966,-0.032210};
  float H_Thresh[] = {-1.139159,-22.071535,-6.305600,-17.291628,-6.963919
		      ,-0.131478,-2.167703,-6.212231,-21.478693,5.797341
		      ,-0.124357,-4.212658,-4.164402,-0.305321,-9.332055
		      ,-19.193649,1.859761,-29.141768,-0.841450,-20.863886
		      ,-8.563315,-24.149193,1.320901,-2.017896,-10.968766
		      ,-13.036150,-4.577792};
  float H_O_Weight[] = {-0.505142,-11.173621,9.066170,0.735692,-5.741280
			,-0.149015,0.021006,-0.772682,0.908176,-1.964467
			,1.856981,0.969065,-0.644269,-0.441946,-3.008178
			,0.813867,0.434074,2.966810,-0.215067,-0.777315
			,1.596130,-0.654612,-0.996872,-0.322435,-0.493111
			,0.405931,-0.192212};
  float O_Thresh[] = {0.003774};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p11_13_t3(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 29;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {0.254255,-0.887242,-1.579081,2.012368,-0.723303
			,8.791495,8.812758,-4.127661,11.873775,11.905257
			,6.773123,-6.299494,6.027431,-1.506644,1.820717
			,-0.323544,-0.699962,-3.978797,-11.411961,1.844618
			,1.424098,-8.800743,-9.508690,-1.145149,0.386247
			,-1.117924,13.554929,-0.872496,14.832987,17.093241
			,28.841957,15.049623,24.710712,23.862129,28.834431
			,30.125513,-46.734810,-131.242874,-28.966066,-5.046007
			,-35.733128,-23.993155,15.216646,-26.223097,22.712505
			,13.020645,38.665573,17.103165,11.041135,24.197248
			,20.120499,-11.458268,17.599365,22.575502,3.968131
			,-42.172119,46.044949,23.475180,-25.347374,-79.076683
			,3.005533,-68.172760,-46.916573,-95.699478,15.307354
			,64.565292,72.481796,46.031441,-92.676750,39.370129
			,1.568504,22.636440,38.364830,-48.256660,-9.527985
			,-84.479095,-47.374386,-8.177492,-52.344959,-105.089554
			,-53.127419,-28.326437,-61.084137,10.215949,46.852234
			,18.518543,-68.150368,-1.475277,-0.397568,-1.917422
			,-0.458916,0.053713,-2.120355,-7.637156,-2.122148
			,5.662523,2.086807,-0.342418,-9.331912,6.105079
			,3.033706,1.538050,0.032651,-1.447235,-1.421295
			,14.073545,0.184878,-1.046800,0.561960,13.473214
			,-0.518217,0.865498,-2.585414,2.578852,-9.300818
			,2.652600,-7.063679,-5.390216,-7.913442,-0.632750
			,-4.190047,7.404703,4.153655,19.407093,21.649860
			,26.166649,16.481375,-6.574354,6.646088,1.374714
			,8.755839,-4.778144,-6.321439,16.052916,-31.781424
			,-5.798396,-3.788064,6.076385,-10.309029,-4.091912
			,-3.668854,-3.199612,22.079184,-9.538126,-3.182677
			,-0.192261,4.521028,7.595740,3.810792,-0.623895
			,1.076497,-5.100084,-1.339660,-11.902866,2.790999
			,14.265047,-1.397019,-1.100532,11.054142,14.156887
			,0.395109,-1.092444,49.914307,122.018433,-14.598628
			,2.460727,-6.314129,-13.549213,-1.016039,1.128922
			,-6.097106,41.270264,28.782763,23.914177,-8.091448
			,-14.166061,9.592519,-6.153887,-1.569943,5.039526
			,4.161227,20.448524,26.802233,23.454138,16.656370
			,5.348643,3.527432,7.032623,0.814405,-4.551698
			,-6.204314,-22.075813,-127.240952,-10.988944,-6.051487
			,12.442115,-21.289907,-2.050966,-5.785337,8.992185
			,-12.118058,-14.446738,-32.590393,-1.663913,-3.463872
			,4.412180,-0.918080,0.394848,-0.293913,-4.905710
			,-4.894585,42.957592,-2.202898,20.963943,8.705330
			,0.733694,-1.480834,1.397624,-0.065858,-2.880416
			,14.483887,27.288408,-7.024310,0.088213,-1.283816
			,-8.271889,-0.612837,-0.008541,-4.300315,3.301595
			,14.444401,6.744441,4.641940,8.476648,2.997131
			,-1.813961,15.865555,1.400457,-4.180232,5.594414
			,10.796182,-3.523435,-5.785693,10.690939,-1.860972
			,7.459911,-7.830221,11.005470,7.181901,6.200663
			,8.672671,5.513043,-0.226515,-2.980952,2.660800
			,8.378541,3.495450,9.732493,8.106985,6.137670
			,-11.620546};
  float H_Thresh[] = {1.964007,-1.876447,-7.507443,1.177902,-14.177758
		      ,-8.902144,-4.243049,-21.057533,-66.500862,-24.204742
		      ,-31.920439,-11.877408,-7.599308,-17.888878,-6.884020
		      ,-8.285959,0.202711,-40.799282,-7.576574,3.575759
		      ,1.620553,-2.443382,16.710852,-6.127536,-2.029591
		      ,-8.612695,-36.866417,-19.184275,6.989274};
  float H_O_Weight[] = {0.350949,4.692076,-0.930832,4.369088,0.795556
			,0.668950,-10.188153,-0.313093,-0.483956,-0.595064
			,0.769318,-2.063986,-1.732364,-0.547622,-1.495908
			,0.724697,0.403709,2.041215,0.506583,10.231455
			,1.872764,1.390064,0.716342,0.687884,3.047024
			,-6.927502,-0.727881,-0.647027,0.709552};
  float O_Thresh[] = {0.072621};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p14_18_t1(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 21;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {3.199570,-7.039454,-1.144994,-10.345678,-0.717025
			,7.320136,-5.030087,-0.471650,-0.850079,9.631785
			,-6.921426,-0.686601,1.873436,-2.154675,10.859210
			,1.925977,-26.649389,-2.043921,-7.008679,-3.261639
			,7.174717,15.790805,-37.939011,-29.025307,19.175768
			,-27.677109,-69.693565,34.212555,48.104027,-22.542387
			,-22.222198,4.097265,-37.730129,-36.242706,-35.782383
			,-27.873999,-31.774633,29.225889,-31.374199,46.495590
			,-33.124477,5.115947,50.890839,70.161224,59.155071
			,-52.732906,53.604446,65.739586,37.760151,-209.138458
			,60.978600,56.785110,-13.109457,82.421654,75.612717
			,74.518967,59.907005,64.342590,-101.193031,63.546146
			,-131.051651,62.431229,22.463268,-6.937952,1.266285
			,0.925315,3.835039,3.895917,0.374872,-5.612980
			,-1.008276,2.130567,3.523959,0.132295,-0.241762
			,0.972317,1.572409,5.147573,1.394641,-3.096695
			,1.570799,1.493759,3.492492,-1.505181,12.818562
			,9.732446,-5.013268,8.876111,-6.839103,12.830342
			,-7.847266,0.935070,0.086978,11.528108,2.840357
			,1.123855,3.100554,-0.611391,-1.764915,1.570104
			,31.608175,-2.977247,22.294550,-2.003173,-0.078075
			,19.301519,1.848029,-0.149269,-7.588900,-1.138545
			,3.588528,-0.069752,0.138873,-10.547857,13.781713
			,0.798318,-4.384416,-0.235496,-2.705637,14.953203
			,0.290100,-12.555643,-1.559338,-4.543153,-2.341845
			,3.051731,-15.376122,2.042122,-3.547815,14.749159
			,-9.056931,-3.147963,1.970811,7.221807,8.249779
			,-26.997108,-5.838334,7.214774,3.096680,5.608147
			,-16.491066,0.490846,30.075756,1.085220,25.627024
			,-0.801178,-2.166731,1.320557,-2.610417,1.593009
			,-4.559985,3.845048,-3.030705,0.546408,-0.815875
			,-1.302372,-1.580154,6.941426,1.019056,0.491379
			,0.578757,4.767103,0.803197,14.025827,0.693916
			,-2.722209,0.106347,-1.311133,-2.732947,-13.296060
			,-0.121512,12.196682,0.721255,-1.133802,-2.531398
			,-0.351705,6.552813,-0.245258,-8.958188,-1.458536
			,-1.616653,0.289676,15.897521,-1.029540,-20.167858
			,0.103730,-6.489370,-12.436177,2.003604};
  float H_Thresh[] = {-17.277275,4.114820,1.644992,-20.793680,0.362656
		      ,-3.801091,1.755117,0.305569,-7.647297,-11.678584
		      ,-3.514322,-1.368163,-3.802866,-1.722315,-20.890432
		      ,-3.250178,-25.110348,0.174834,-15.118277,13.674134
		      ,-5.911378};
  float H_O_Weight[] = {-0.906263,-5.391995,-1.157650,1.442591,-3.071737
			,-0.378128,-0.457259,0.721263,-2.270399,-1.523669
			,3.944212,-1.053538,-0.662563,-0.955973,-0.229678
			,-0.718471,0.632109,-0.734910,0.422043,-0.428271
			,-1.472066};
  float O_Thresh[] = {-0.069997};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p14_18_t2(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 24;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {-3.380136,-1.135163,-3.599914,-24.973217,2.491927
			,10.319602,0.700984,0.055328,-9.572863,6.023386
			,0.814572,-2.976593,-0.347965,17.731113,3.994063
			,-0.148170,11.340511,-5.731233,-12.400249,1.387344
			,0.419512,8.165179,3.715097,-11.615696,2.583365
			,6.549341,12.308177,27.596363,35.415401,-5.784149
			,-50.876057,11.432244,-81.865028,13.948427,-24.707508
			,63.277699,22.259415,-31.583187,-59.167519,-31.484743
			,-2.329501,9.675364,20.572268,-28.105270,5.707983
			,-25.801613,-44.580891,22.639652,30.235987,9.185473
			,-88.890343,-94.314583,-141.435776,22.793385,95.432610
			,-77.492928,103.264023,-186.747421,36.525829,-229.365311
			,43.732288,37.176689,104.642227,75.088799,-0.232394
			,-65.489143,-98.664024,53.365681,54.826149,56.072266
			,103.694855,-58.759464,0.906420,16.208824,-5.626152
			,0.067859,5.345044,2.475055,0.445470,-0.481771
			,-0.763114,-7.990931,2.744881,0.559805,-9.049590
			,16.053728,-8.645160,3.148981,1.545267,-12.826681
			,-9.806294,2.783127,-7.915681,-1.289790,0.448950
			,-10.768787,-32.532497,-9.188558,29.270319,113.358360
			,97.300461,-6.964613,28.064625,1.289544,9.819228
			,26.838778,-12.593746,38.378563,-7.249500,-4.088984
			,13.190929,-11.278858,-6.278062,11.217149,56.918072
			,-6.076216,-0.778584,2.470195,6.828117,30.078596
			,42.015171,-6.739448,-12.792689,-10.843225,1.157812
			,0.034288,8.902040,3.735509,-26.610142,-24.099178
			,11.444094,-5.181591,-0.017643,-7.669642,-8.942697
			,4.607644,10.875838,-15.986262,-19.115093,3.327044
			,-1.272684,33.104183,-3.577001,-17.319698,-45.432732
			,4.316665,13.289506,52.827980,61.903408,-10.064597
			,37.050766,10.640875,40.803890,25.276575,-13.318869
			,44.728725,9.553937,-30.010485,12.762658,0.582986
			,-24.202328,6.828737,32.670704,-8.499022,2.509053
			,-30.729357,8.695580,24.105362,22.771059,6.320170
			,17.622465,17.923342,22.693565,2.044944,-2.419853
			,1.465247,17.109268,-9.631621,3.707756,-6.453905
			,-0.786013,3.857229,-3.985719,7.177620,-1.822190
			,7.409497,14.235882,1.646543,-4.153118,20.812519
			,-1.145312,5.990725,4.650846,-6.507271,-32.010006
			,-48.384125,-16.624491,8.015402,0.597853,2.297649
			,1.088388,-20.028086,5.852522,-3.261894,-2.902347
			,6.806984,-7.914451,8.351144,5.858640,-7.267718
			,-24.557932,-0.863479,-5.425487,-6.777747,-1.406155
			,-12.362561};
  float H_Thresh[] = {3.028005,-12.574837,-15.112863,-85.192223,-107.008102
		      ,-5.883706,-33.496773,-8.240110,-24.168364,6.741688
		      ,-1.855648,-33.780293,2.697949,-7.089548,0.276746
		      ,-7.738922,2.236442,-9.714922,-46.229763,1.065543
		      ,3.629447,-19.733215,-5.017787,-25.786951};
  float H_O_Weight[] = {-0.266490,-7.469507,4.390029,6.224319,0.792311
			,-0.368458,-1.366321,2.612978,-0.397671,0.553956
			,-3.888159,1.054949,-0.774808,-0.282426,-1.299502
			,-0.988678,-0.382289,0.553578,1.122869,-3.642761
			,-5.019010,-2.826515,-0.764428,7.668901};
  float O_Thresh[] = {-0.007092};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}


void PacPidCand2Ntuple::NNout_p14_18_t3(float* pattern, float* output){
  
// NN node numbers: Input, hidden, and output
  int MaxIndim = 9;
// int Indim = 9;
  int Hidden = 21;
  int Outdim = 1;
  int ivalid = 0;
  float I_SUM = 0.0;
  float* H_O_SUM;
  H_O_SUM = new float[Outdim];
  
// Input scales
  float scale[] = {1.000000,1.000000,1.000000,1.000000,1.000000
		   ,1.000000,1.000000,1.000000,1.000000};
  
// veto input patterns (0 == DO NOT USE)
  int veto[] = {1,1,1,1,1,1,1,1,1};
  
// NN node weights and threshholds
  float I_H_Weight[] = {-0.912031,-6.343918,-3.809599,-40.121124,-7.029084
			,-17.260618,-2.641922,-1.701540,1.255093,-4.856995
			,-7.192831,2.624448,-0.467317,5.070120,12.715511
			,12.346029,-4.105708,5.999301,3.129243,-18.542389
			,7.237540,11.431018,88.404564,-37.317745,26.030373
			,13.768045,8.529660,30.672899,29.324869,-34.735893
			,10.190025,30.615320,-92.020607,11.509786,-42.178940
			,-36.746513,4.476262,15.084729,-5.514809,46.271587
			,-25.447002,22.278847,50.684597,-212.302307,96.633217
			,-91.922089,-49.656883,-27.164574,50.598442,-108.529037
			,63.915901,8.590190,-136.588013,149.516693,-48.940422
			,115.970032,64.714149,7.563162,-70.933502,37.584641
			,-222.493332,82.137367,-87.225197,1.818411,6.595863
			,1.032799,4.464560,4.539032,11.311253,-6.827978
			,0.561561,6.151027,-1.997022,2.003586,-7.844569
			,2.472588,6.269403,-0.420084,10.285390,-0.614373
			,0.702889,-1.235545,0.751984,7.513115,13.715720
			,-44.803493,35.243210,-4.305363,-82.847595,-2.937542
			,-1.457211,28.684668,12.988535,-10.065824,-4.442765
			,19.035326,-1.854351,58.778969,24.785460,3.061210
			,-6.376042,-3.793427,11.267786,58.532623,-29.492636
			,47.677044,20.803528,1.958240,-12.863443,144.854996
			,-20.432734,-3.630363,-0.582949,5.827174,-6.683339
			,-7.418994,-18.346725,-11.753176,-8.676314,42.126667
			,-3.728024,-0.200821,39.276146,-1.824340,22.042871
			,30.240257,-2.179069,-36.110802,63.128891,10.715034
			,-164.130875,-16.376284,9.447223,39.343739,16.373848
			,7.363682,6.003677,20.306566,-18.617861,55.921288
			,-8.381248,10.904469,-0.544363,-8.389043,9.313037
			,1.322868,-56.940548,11.784877,10.911110,-1.527281
			,-0.398703,49.823559,-13.994326,-0.133103,-3.059897
			,-1.657161,-1.349131,-3.035515,-2.297844,-18.236456
			,33.681427,5.653768,-9.519320,-0.835475,11.653018
			,-1.100756,26.835239,3.024655,1.282355,19.581619
			,0.409949,4.384077,32.826736,16.880127,0.644331
			,-1.085496,1.957085,0.559308,3.130708,-15.359979
			,8.398887,5.389077,5.448095,-6.619625,0.234655
			,19.898464,-3.953809,-24.841444,2.417015};
  float H_Thresh[] = {-39.276360,10.554577,-43.185463,-2.994719,-1.587278
		      ,-2.609708,-6.006000,-31.126774,-18.205200,4.923743
		      ,1.489551,3.308110,8.168584,-101.212311,-40.828812
		      ,-4.025290,3.163165,-36.132591,-6.410570,-60.864758
		      ,22.581148};
  float H_O_Weight[] = {-0.858317,0.466127,-0.924788,9.863556,1.571438
			,0.851746,-0.589514,1.736075,-1.205241,-2.529839
			,2.578126,-0.731793,1.684136,-1.104784,-8.942947
			,-0.570735,8.624597,-1.574590,3.274760,-0.863343
			,0.613412};
  float O_Thresh[] = {0.377187};
  
// Loop over output, hidden, and input nodes
// to produce the neural net output
//============================================
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] = 0.0;
    output[out] = 0.0;
  }
  for (int h = 0; h<Hidden; h++){
    I_SUM = 0.;
    ivalid = 0;
    for (int i = 0; i<MaxIndim; i++){
      if(veto[i]){
        I_SUM += I_H_Weight[Hidden*ivalid+h] * pattern[i]/scale[i];
        ivalid++;
      }
    }
    I_SUM += H_Thresh[h];
    for (int out = 0; out<Outdim; out++) H_O_SUM[out] += H_O_Weight[Outdim*h+out]*Afunction(I_SUM);
  }
  for (int out = 0; out<Outdim; out++){
    
    H_O_SUM[out] += O_Thresh[out];
    output[out] = Afunction(H_O_SUM[out]);
  }
  delete[] H_O_SUM;
  return;
}

// Activation function
float PacPidCand2Ntuple::Afunction(float INPUT){
  return( 1.0 / ( 1.0 + exp(-2.0*INPUT) ) );
}
*/

void PacPidCand2Ntuple::convertJpsiKMC(const BtaCandidate *c, const BtaCandidate *cTag, const BtaCandidate *mc, AbsEvent *ev, HepTuple *nt){
  // c points to the muon(electron) from charged track, cTag points to muon(electron) from muNNLoose (eBremReco) track, mc points to the kaon

  //From the muon from charged track get the parent lundid

  // Fill ntuple with id, mother id, and grandmother id
  if (_theAssoc!=NULL) {
    BtaCandidate *tc = _theAssoc->mcFromReco(c);
  //BtaCandidate *grandmotherLepton(0);
    BtaCandidate *motherLepton(0);
    if (tc!=0){
      if (tc->theMother()!=NULL) {
	motherLepton = tc->theMother();
      //if ((tc->theMother())->theMother()!=NULL) {
      //  grandmotherLepton = motherLepton->theMother();
      //}
      }
    }

    BtaCandidate *tcTag = _theAssoc->mcFromReco(cTag);
    BtaCandidate *grandmotherLeptonTag(0);
    BtaCandidate *motherLeptonTag(0);
    BtaCandidate *motherKaon(0);
    int otherLundIdTruth = 0;
    int otherMotherIdTruth = 0;
    int otherMotherNDaughters = 0;
    int otherGrandmotherIdTruth = 0;
    int otherGrandmotherNDaughters = 0;
    
    if (tcTag!=0) {
      otherLundIdTruth = tcTag->pdtEntry()->lundId();
      if (tcTag->theMother()!=NULL) {
	motherLeptonTag = tcTag->theMother();
        otherMotherIdTruth = motherLeptonTag->pdtEntry()->lundId();
	otherMotherNDaughters = motherLeptonTag->nDaughters();
        if ((tcTag->theMother())->theMother()!=NULL) {
          grandmotherLeptonTag = motherLeptonTag->theMother();
	  otherGrandmotherIdTruth = grandmotherLeptonTag->pdtEntry()->lundId();
	  otherGrandmotherNDaughters = grandmotherLeptonTag->nDaughters();
	}
      }
    }
    nt->column("otherLundIdTruth", otherLundIdTruth);
    nt->column("otherMotherIdTruth", otherMotherIdTruth);
    nt->column("otherMotherNDaughters", otherMotherNDaughters);
    nt->column("otherGrandmotherIdTruth", otherGrandmotherIdTruth);
    nt->column("otherGrandmotherNDaughters", otherGrandmotherNDaughters);

    // Fill ntuple with sigmode (0 if not signal, otherwise kaon id)
    BtaCandidate *tmc = _theAssoc->mcFromReco(mc);
    int kaonLundIdTruth = 0;
    int kaonMotherIdTruth = 0;
    int kaonGrandmotherIdTruth = 0;
    int kaonMotherNDaughters = 0;
    if(tmc!=0) {
      kaonLundIdTruth = tmc->pdtEntry()->lundId();
      if (tmc->theMother()!=NULL) {
        motherKaon = tmc->theMother();
        kaonMotherIdTruth = motherKaon->pdtEntry()->lundId();
        kaonMotherNDaughters = motherKaon->nDaughters();
        if (motherKaon->theMother()!=NULL) {
          kaonGrandmotherIdTruth = motherKaon->theMother()->pdtEntry()->lundId();
        }
      }
    }
    nt->column("kaonLundIdTruth", kaonLundIdTruth);
    nt->column("kaonMotherIdTruth", kaonMotherIdTruth);
    nt->column("kaonMotherNDaughters", kaonMotherNDaughters);
    nt->column("kaonGrandmotherIdTruth", kaonGrandmotherIdTruth);

    bool sameBcand = (grandmotherLeptonTag != 0 && motherKaon != 0 && grandmotherLeptonTag->energy() == motherKaon->energy());
    nt->column("isBCand", sameBcand);

    bool sameJpsiCand = (motherLepton != 0 && motherLeptonTag != 0 && motherLepton->energy() == motherLeptonTag->energy());
    nt->column("isJpsiCand", sameJpsiCand);

    int sigmode = 0;
    if(otherGrandmotherNDaughters == 2 && otherMotherNDaughters == 2) {
      if(sameBcand && sameJpsiCand) {
        if (abs(otherMotherIdTruth) == 443) {
          if (abs(kaonMotherIdTruth) == 521 || abs(kaonMotherIdTruth) == 511) {
      int absId = abs(kaonLundIdTruth);
            if( absId == 321 ||  absId == 313 || absId == 323 ||  absId == 311 || absId == 310 ){
              sigmode = kaonLundIdTruth;
            }
          }
        }
      }
    }
    nt->column("sigmode", sigmode);
  }
}


/////////////////////////////////////////////
// Sasha Telnov, 2007/05/18: add a routine to reduce precision of floats 
// (and thus to improve packing of ROOT trees) 
// Will change the variable at most by 1/2^(24-bits_to_kill) of its magnitude. 
// The average change is zero.
// The RMS of the change is approximately 0.41/2^(24-bits_to_kill):
//
//                bits_to_kill  relative |change|    r.m.s.   
//  32-bit float       0              0                0
// "24-bit float"      8         < 1.5e-5           6.3e-6
// "20-bit float"      12        < 2.5e-4           1.0e-4 
// "16-bit float"      16        < 3.9e-3           1.6e-3
//
// Allowed values: 0 <= bits_to_kill <= 18, otherwise no action taken (same as 0)
float PacPidCand2Ntuple::packFloat2(const float &fl, int bits_to_kill) {
  if (_doFloatPacking->value() == false) return fl; // do no packing 
  if (fl != fl || fl==0) return fl; // the first check is for NAN
  if (_maxFloatPacking->value() < bits_to_kill) bits_to_kill = _maxFloatPacking->value();
  int temp;
  switch (bits_to_kill) {
  case 1:  temp = *((int*)(&fl)) & 0xFFFFFFFE; break;
  case 2:  temp = *((int*)(&fl)) & 0xFFFFFFFC; break;
  case 3:  temp = *((int*)(&fl)) & 0xFFFFFFF8; break;
  case 4:  temp = *((int*)(&fl)) & 0xFFFFFFF0; break;
  case 5:  temp = *((int*)(&fl)) & 0xFFFFFFE0; break;
  case 6:  temp = *((int*)(&fl)) & 0xFFFFFFC0; break;
  case 7:  temp = *((int*)(&fl)) & 0xFFFFFF80; break;
  case 8:  temp = *((int*)(&fl)) & 0xFFFFFF00; break;
  case 9:  temp = *((int*)(&fl)) & 0xFFFFFE00; break;
  case 10: temp = *((int*)(&fl)) & 0xFFFFFC00; break;
  case 11: temp = *((int*)(&fl)) & 0xFFFFF800; break;
  case 12: temp = *((int*)(&fl)) & 0xFFFFF000; break;
  case 13: temp = *((int*)(&fl)) & 0xFFFFE000; break;
  case 14: temp = *((int*)(&fl)) & 0xFFFFC000; break;
  case 15: temp = *((int*)(&fl)) & 0xFFFF8000; break;
  case 16: temp = *((int*)(&fl)) & 0xFFFF0000; break;
  case 17: temp = *((int*)(&fl)) & 0xFFFE0000; break;
  case 18: temp = *((int*)(&fl)) & 0xFFFC0000; break;
  //default: temp = *((int*)(&fl)) & 0xFFFFFFFF; break;
  default: return fl;
  }
  float tempFloat = fl - *(float*)(&temp) + fl; // we want to round off, not just truncate!
  switch (bits_to_kill) {
  case 1:  temp = *((int*)(&tempFloat)) & 0xFFFFFFFE; break;
  case 2:  temp = *((int*)(&tempFloat)) & 0xFFFFFFFC; break;
  case 3:  temp = *((int*)(&tempFloat)) & 0xFFFFFFF8; break;
  case 4:  temp = *((int*)(&tempFloat)) & 0xFFFFFFF0; break;
  case 5:  temp = *((int*)(&tempFloat)) & 0xFFFFFFE0; break;
  case 6:  temp = *((int*)(&tempFloat)) & 0xFFFFFFC0; break;
  case 7:  temp = *((int*)(&tempFloat)) & 0xFFFFFF80; break;
  case 8:  temp = *((int*)(&tempFloat)) & 0xFFFFFF00; break;
  case 9:  temp = *((int*)(&tempFloat)) & 0xFFFFFE00; break;
  case 10: temp = *((int*)(&tempFloat)) & 0xFFFFFC00; break;
  case 11: temp = *((int*)(&tempFloat)) & 0xFFFFF800; break;
  case 12: temp = *((int*)(&tempFloat)) & 0xFFFFF000; break;
  case 13: temp = *((int*)(&tempFloat)) & 0xFFFFE000; break;
  case 14: temp = *((int*)(&tempFloat)) & 0xFFFFC000; break;
  case 15: temp = *((int*)(&tempFloat)) & 0xFFFF8000; break;
  case 16: temp = *((int*)(&tempFloat)) & 0xFFFF0000; break;
  case 17: temp = *((int*)(&tempFloat)) & 0xFFFE0000; break;
  case 18: temp = *((int*)(&tempFloat)) & 0xFFFC0000; break;
  default: temp = *((int*)(&tempFloat)) & 0xFFFFFFFF; break;
  }
  return *(float*)(&temp);
}

// The variable is the azimuthal angle
float PacPidCand2Ntuple::packFloatPhi2(const float &fl, int bits_to_kill) {
  static const float floatPI = 3.141592654; 
  float tempFloat = packFloat(fl,bits_to_kill);
  if (tempFloat < -floatPI) tempFloat = tempFloat + floatPI + floatPI;
  if (tempFloat > floatPI) tempFloat = tempFloat - floatPI - floatPI;
  return tempFloat; 
}


// avtelnov, ::printDebuggingInfo() and ::printDebuggingInfo2() added on 2007/12/06
void PacPidCand2Ntuple::printDebuggingInfo(const BtaCandidate *c) {
  ErrMsg(warning) << " =============================================================================" << endmsg;
  ErrMsg(warning) << " Got cand with LundId " << c->pdtEntry()->lundId() << " and the following info in DchPidInfo and BtaPidInfo:" << endmsg;
  printDebuggingInfo2(c);
  
  ErrMsg(warning) << " Now, I will clone this cand, call setType, and see what I get now..." << endmsg;
  int ncharge = static_cast<int>(c->charge());
  BtaCandidate tempCand(*c);

  ErrMsg(warning) << " ... ELECTRON:" << endmsg;
  tempCand.setType(Pdt::lookup(PdtPid::electron, ncharge));
  printDebuggingInfo2(&tempCand);

  ErrMsg(warning) << " ... MUON:" << endmsg;
  tempCand.setType(Pdt::lookup(PdtPid::muon, ncharge));
  printDebuggingInfo2(&tempCand);

  ErrMsg(warning) << " ... PION:" << endmsg;
  tempCand.setType(Pdt::lookup(PdtPid::pion, ncharge));
  printDebuggingInfo2(&tempCand);

  ErrMsg(warning) << " ... KAON:" << endmsg;
  tempCand.setType(Pdt::lookup(PdtPid::kaon, ncharge));
  printDebuggingInfo2(&tempCand);

  ErrMsg(warning) << " ... PROTON:" << endmsg;
  tempCand.setType(Pdt::lookup(PdtPid::proton, ncharge));
  printDebuggingInfo2(&tempCand);

}

// printDebuggingInfo2() is to be called from printDebuggingInfo()
void PacPidCand2Ntuple::printDebuggingInfo2(const BtaCandidate *c) {
  ErrMsg(warning) << " momentum at IP = " << c->p() << endmsg;
  const BtaPidQual* pQual = c->getMicroAdapter()->getPidQual();
  if (pQual) {
    ErrMsg(warning) << " pQual->thetaAtEMC() = " << pQual->thetaAtEMC() << endmsg;
    ErrMsg(warning) << " pQual->dEdXDch() = " << pQual->dEdXDch() << endmsg;
    ErrMsg(warning) << " pQual->dEdXSvt() = " << pQual->dEdXSvt() << endmsg;
    ErrMsg(warning) << " pQual->deltaDchMomentum() = " << pQual->deltaDchMomentum() 
		    << "; +p = " << c->p() + pQual->deltaDchMomentum() << endmsg;
    ErrMsg(warning) << " pQual->thetaC() = " << pQual->thetaC() << endmsg;
    ErrMsg(warning) << " pQual->thetaCErr() = " << pQual->thetaCErr() << endmsg;
    ErrMsg(warning) << " pQual->deltaDrcMomentum() = " << pQual->deltaDrcMomentum() << endmsg;
  //ErrMsg(warning) << "  = " <<  << endmsg;

  } else ErrMsg(warning) << " No BtaPidQual!" << endmsg;

  const BtaPidInfo* pInfo = c->getMicroAdapter()->getPidInfo();
  if (pInfo) {


  } else ErrMsg(warning) << " No BtaPidInfo!" << endmsg;
  
  const PidInfoSummary *realPidInfo = c->pidInfoSummary();
  if (realPidInfo != 0) {
    
    const DchPidInfo* dchPid = realPidInfo->dchPidInfo();
    if (dchPid != 0) {
      ErrMsg(warning) << " DCH dE/dx expected for e hypo = " << dchPid->getExpected(PdtPid::electron) << endmsg;
      ErrMsg(warning) << " DCH dE/dx expected for mu hypo = " << dchPid->getExpected(PdtPid::muon) << endmsg;
      ErrMsg(warning) << " DCH dE/dx expected for pi hypo = " << dchPid->getExpected(PdtPid::pion) << endmsg;
      ErrMsg(warning) << " DCH dE/dx expected for K hypo = " << dchPid->getExpected(PdtPid::kaon) << endmsg;
      ErrMsg(warning) << " DCH dE/dx expected for p hypo = " << dchPid->getExpected(PdtPid::proton) << endmsg;
      ErrMsg(warning) << " DCH dE/dx error for e hypo = " << dchPid->getdEdxErr(PdtPid::electron) << endmsg;
      ErrMsg(warning) << " DCH dE/dx error for mu hypo  = " << dchPid->getdEdxErr(PdtPid::muon) << endmsg;
      ErrMsg(warning) << " DCH dE/dx error for pi hypo  = " << dchPid->getdEdxErr(PdtPid::pion) << endmsg;
      ErrMsg(warning) << " DCH dE/dx error for K hypo  = " << dchPid->getdEdxErr(PdtPid::kaon) << endmsg;
      ErrMsg(warning) << " DCH dE/dx error for p hypo  = " << dchPid->getdEdxErr(PdtPid::proton) << endmsg;
      ErrMsg(warning) << " p at DCH for e hypo = " << dchPid->getDchMomentum(PdtPid::electron) << endmsg;
      ErrMsg(warning) << " p at DCH for mu hypo = " << dchPid->getDchMomentum(PdtPid::muon) << endmsg;
      ErrMsg(warning) << " p at DCH for pi hypo = " << dchPid->getDchMomentum(PdtPid::pion) << endmsg;
      ErrMsg(warning) << " p at DCH for K hypo = " << dchPid->getDchMomentum(PdtPid::kaon) << endmsg;
      ErrMsg(warning) << " p at DCH for p hypo = " << dchPid->getDchMomentum(PdtPid::proton) << endmsg;
    //ErrMsg(warning) << "  = " <<  << endmsg;
    } else ErrMsg(warning) << " No DchPidInfo!" << endmsg;
    
    const SvtPidInfo* svtPid = realPidInfo->svtPidInfo();
    if (svtPid != 0) {
      ErrMsg(warning) << " p at SVT = " << svtPid->momentum() << endmsg;
    //ErrMsg(warning) << "  = " <<  << endmsg;
      
    } else ErrMsg(warning) << " No SvtPidInfo!" << endmsg;
    
  } else ErrMsg(warning) << " No PidInfoSummary!" << endmsg;

}
