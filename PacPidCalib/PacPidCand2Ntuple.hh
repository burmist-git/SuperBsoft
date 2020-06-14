//--------------------------------------------------------------------------
// File and Version Information:
//  $Id: $
//
// Description:
//  Dump Micro information to an ntuple.
//  Adapted from BetaPidCalibNtuple/BtaCand2Ntuple.hh
//
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
//
//------------------------------------------------------------------------

#ifndef PACPIDCAND2NTUPLE_HH
#define PACPIDCAND2NTUPLE_HH

#include "PDT/PdtPid.hh"
#include "AbsPid/PidSystem.hh"

// All of these are declared as actual data objects below, so must have
// full headers.  They ought to be pointers (and "new"ed in the .cc file).
#include "AbsParm/AbsParmIfdStrKey.hh"
#include "BetaEvent/BtaParametrizable.hh"

// NA
//#include "BetaPid/PidMuonSPR.hh"
//#include "BetaPid/PidMuonBDTLoPSelector.hh"
//#include "TrgTools/TrgFctTimePointInspector.hh"

#include <fstream>
#include <vector>
#include <string>

class BtaCandidate;
class BtaMcAssocChiSq;
class BtaMcAssoc;
class AbsEvent;
class HepTuple;

// NA
//class PidLHElectronSelector;
//class PidKaonSMSSelector;

class AppModule;

// NA
//class PidKaonMicroSelector;
//class PidKaonBDTSelector;
//class PidKMSelector;

class IfrNNKernelSetup;

// NA
//class OprTrickleRegions;

class AbsParmBool;
template <class T> class AbsParmGeneral;

// Added at Al Eisner's request
// NA
//class EmcMergedPi0Identifier;
//class BtaMergedPi0Algo;

template <class T> class AbsParmVector;
template <class T> class HepAList;

class PacPidCand2Ntuple : public BtaParametrizable {

public :
  PacPidCand2Ntuple(AppModule *aModule) ;
  virtual ~PacPidCand2Ntuple();
  void beginJob(const AbsEvent* ev);
  void convert(const BtaCandidate *c, AbsEvent *ev, HepTuple *nt);
  void convertJpsiKMC(const BtaCandidate *c1, const BtaCandidate *c2, const BtaCandidate *mc, AbsEvent *ev, HepTuple *nt);
  void setEvent(AbsEvent *ev,AbsParmVector<std::string> *pidVector=0,AbsParmVector<std::string> *tagVector=0);
  float lhood(const BtaCandidate *cand,PdtPid::PidType hypo,PidSystem::System d);
  float consval(const BtaCandidate *cand,PdtPid::PidType hypo,PidSystem::System d);

  void fillEmptyCalQual(HepTuple*);

  // From Giampi ----------
  void addNphotFromGlbLikelihood(const BtaCandidate *c, HepTuple *ntp);
  int calNPhot(double sigLvl, double nExp);
  int min(int x1,int x2) {
    if (x1<x2) { return x1; } else { return x2; }
  };
  void setIfrNNKernel(IfrNNKernelSetup *nn) { _ifrNNKernel=nn;}
  // ----------------------


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
// Packing is controlled by two parameters: 
//   _doFloatPacking : bypass packing if set to false
//   _maxFloatPacking : remove no more bits than the value of this parameter
  static float packFloat2(const float &fl, int bits_to_kill);
  static float packFloatPhi2(const float &fl, int bits_to_kill); // the variable is the azimuthal angle
// and fakes that do nothing, to speed things up without changing a lot of code
  static inline float packFloat(const float &fl, int bits_to_kill) {return fl;}
  static inline float packFloatPhi(const float &fl, int bits_to_kill) {return fl;}

// add EMC NN used in muon PID
  // NA
  /*
  static void NNout_p01_t1(float* , float* );
  static void NNout_p01_t2(float* , float* );
  static void NNout_p01_t3(float* , float* );
  static void NNout_p02_03_t1(float* , float* );
  static void NNout_p02_03_t2(float* , float* );
  static void NNout_p02_03_t3(float* , float* );
  static void NNout_p04_06_t1(float* , float* );
  static void NNout_p04_06_t2(float* , float* );
  static void NNout_p04_06_t3(float* , float* );
  static void NNout_p11_13_t1(float* , float* );
  static void NNout_p11_13_t2(float* , float* );
  static void NNout_p11_13_t3(float* , float* );
  static void NNout_p14_18_t1(float* , float* );
  static void NNout_p14_18_t2(float* , float* );
  static void NNout_p14_18_t3(float* , float* );
  static float Afunction(float);
  */

private :
  BtaMcAssoc* _theAssoc;
  AbsEvent *_currentEvent;
  std::vector< HepAList<BtaCandidate>* > _listOfPidLists;
  std::vector<std::string> _namesOfPidLists;
  std::vector<std::string> _namesOfTagBits;
  AbsParmIfdStrKey _mcTruthList;
  int _nSelectors;
  
  // NA
  //PidKaonSMSSelector *_theSms;
  //PidKaonMicroSelector *_theNN;
  //PidKaonBDTSelector*  _theKaonBDT;
  //PidKMSelector*  _theKM;
  
  IfrNNKernelSetup *_ifrNNKernel;

// muBDT stuff
  // NA
  //PidMuonBDTSelector _muonTree;
  //PidMuonBDTLoPSelector _muonLoPTree;
  //friend class PidMuonSPR;

// avtelnov, 2007/08/28: control packing
  static AbsParmBool* _doFloatPacking; // bypass packing if set to false
  static AbsParmGeneral<int>* _maxFloatPacking; // remove no more bits than the value of this parameter

// if true, take detailed dE/dx parameterizations from CDB, 
// not the ASCII files in PidDchSvtDrcCalib
  static AbsParmBool* _useCDBforDEDX; 

  // NA
  /*
  AbsParmBool _getTrickleInfo;

  AbsParmIfdStrKey _timePointLER;
  TrgFctTimePointInspector _timePointInspector_LER;
  OprTrickleRegions* _trickleRegion_LER;
  AbsParmIfdStrKey _timePointHER;
  TrgFctTimePointInspector _timePointInspector_HER;
  OprTrickleRegions* _trickleRegion_HER;
  */

// Added at Al Eisner's request
  // NA
  //EmcMergedPi0Identifier *_MergedPi0Workhorse;
  //BtaMergedPi0Algo *_MergedPi0algo;

  void makeBDTCutCol(double momentum, double theta,
		     double runNum, double charge, const HepString &str, HepTuple *nt);

  int _warningCounter;

// avtelnov, 2007/12/06: print to ErrMsg(warning) all kinds of information about the BtaCandidate.
// This method is not intended ever to be called in production, just in private test releases   
  void printDebuggingInfo(const BtaCandidate *c); 
  void printDebuggingInfo2(const BtaCandidate *c); 

};

#endif
