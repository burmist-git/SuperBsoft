//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: $
//
// Description:
//     Adapted from BetaPidCalib/BtaPidCSConcreteFactory.hh
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Nicolas ARNAUD (SuperB)
//      Gautier Hamel de Monchenault
//
// Copyright Information:
//	 (C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------

#ifndef PACPIDCSCONCRETEFACTORY_HH
#define PACPIDCSCONCRETEFACTORY_HH

#include "BaBar/BaBar.hh"

#include "PacPidCalib/PacPidCalibSample.hh"
#include "BetaEvent/BtaConcreteFactory.hh"

template <class Q>
class PacPidCSConcreteFactory : public BtaConcreteFactory<Q, PacPidCalibSample>
{
public:

  // constructor
  PacPidCSConcreteFactory() {}

  // virtual destructor
  virtual ~PacPidCSConcreteFactory() {}

};

// declare & instantiate

#include "PacPidCalib/PacPidDstarSample.hh"
static  PacPidCSConcreteFactory<PacPidDstarSample> thePacPidDstarSampleFactory;

//#include "BetaPidCalib/BtaPidMCCalibSample.hh"
//static PacPidCSConcreteFactory<BtaPidMCCalibSample>  theBtaPidMCCalibSampleFactory;

//#include "BetaPidCalib/BtaPidAllCalibSample.hh"
//static PacPidCSConcreteFactory<BtaPidAllCalibSample>  theBtaPidAllCalibSampleFactory;

//#include "BetaPidCalib/BtaEmcBhabhaSample.hh"
//static PacPidCSConcreteFactory<BtaEmcBhabhaSample>  theBtaEmcBhabhaSampleFactory;

//#include "BetaPidCalib/BtaEmcRadBhabhaSample.hh"
//static PacPidCSConcreteFactory<BtaEmcRadBhabhaSample>  theBtaEmcRadBhabhaSampleFactory;

//#include "BetaPidCalib/BtaTrkBhabhaSample.hh"
//static PacPidCSConcreteFactory<BtaTrkBhabhaSample>  theBtaTrkBhabhaSampleFactory;

//#include "BetaPidCalib/BtaPidBhabhaSample.hh"
//static PacPidCSConcreteFactory<BtaPidBhabhaSample>  theBtaPidBhabhaSampleFactory;

//#include "BetaPidCalib/BtaVcsSample.hh"
//static PacPidCSConcreteFactory<BtaVcsSample>  theBtaVcsSampleFactory;

//#include "BetaPidCalib/BtaKinIFRBhabhaSample.hh"
//static PacPidCSConcreteFactory<BtaKinIFRBhabhaSample>  theBtaKinIFRBhabhaSampleFactory;

//#include "BetaPidCalib/BtaPidKsSample.hh"
//static PacPidCSConcreteFactory<BtaPidKsSample>  theBtaPidKsSampleFactory;

//#include "BetaPidCalib/BtaPidDstarSample.hh"
//static PacPidCSConcreteFactory<BtaPidDstarSample>  theBtaPidDstarSampleFactory;

//#include "BetaPidCalib/BtaPideeeeCalibSample.hh"
//static PacPidCSConcreteFactory<BtaPideeeeCalibSample>  theBtaPideeeeCalibSampleFactory;

//#include "BetaPidCalib/BtaeemumuSample.hh"
//static PacPidCSConcreteFactory<BtaeemumuSample>  theBtaeemumuSampleFactory;

//#include "BetaPidCalib/BtamumuSample.hh"
//static PacPidCSConcreteFactory<BtamumuSample>  theBtamumuSampleFactory;

//#include "BetaPidCalib/BtamumuEmcSample.hh"
//static PacPidCSConcreteFactory<BtamumuEmcSample>  theBtamumuEmcSampleFactory;

//#include "BetaPidCalib/BtamumugammaSample.hh"
//static PacPidCSConcreteFactory<BtamumugammaSample>  theBtamumugammaSampleFactory;

//#include "BetaPidCalib/Btamumugamma2Sample.hh"
//static PacPidCSConcreteFactory<Btamumugamma2Sample>  theBtamumugamma2SampleFactory;

//#include "BetaPidCalib/BtaTau31Sample.hh"
//static PacPidCSConcreteFactory<BtaTau31Sample> theBtaTau31SampleFactory;

//#include "BetaPidCalib/BtaPidProtonSample.hh"
//static PacPidCSConcreteFactory<BtaPidProtonSample> theBtaPidProtonSampleFactory;

//#include "BetaPidCalib/BtaPidLambdaProtonSample.hh"
//static PacPidCSConcreteFactory<BtaPidLambdaProtonSample> theBtaPidLambdaProtonSampleFactory;

//#include "BetaPidCalib/BtaLambdaTrkOnlySample.hh"
//static PacPidCSConcreteFactory<BtaLambdaTrkOnlySample> theBtaLambdaTrkOnlySampleFactory;

//#include "BetaPidCalib/BtaPidV0SampleModule.hh"
//static PacPidCSConcreteFactory<BtaPidV0SampleModule> theBtaPidV0SampleModuleFactory;

//#include "BetaPidCalib/BtadEdxProtonSample.hh"
//static PacPidCSConcreteFactory<BtadEdxProtonSample> theBtadEdxProtonSampleFactory;

//#include "BetaPidCalib/BtadEdxKaonSample.hh"
//static PacPidCSConcreteFactory<BtadEdxKaonSample> theBtadEdxKaonSampleFactory;

//#include "BetaPidCalib/BtaDircKaonSample.hh"
//static PacPidCSConcreteFactory<BtaDircKaonSample> theBtaDircKaonSampleFactory;

//#include "BetaPidCalib/BtaPidgg4piCalibSample.hh"
//static PacPidCSConcreteFactory<BtaPidgg4piCalibSample> theBtaPidgg4piCalibSampleFactory;

//#include "BetaPidCalib/BtaPidD0AsymSample.hh"
//static PacPidCSConcreteFactory<BtaPidD0AsymSample> theBtaPidD0AsymSampleFactory;

//#include "BetaPidCalib/BtaPidGammaeeSample.hh"
//static PacPidCSConcreteFactory<BtaPidGammaeeSample> theBtaPidGammaeeSampleFactory;

//#include "BetaPidCalib/BtaPhiGammaSample.hh"
//static PacPidCSConcreteFactory<BtaPhiGammaSample> theBtaPhiGammaSampleFactory;

//#include "BetaPidCalib/BtaPidTau31Pi0MergedSample.hh"
//static BtaPidCSConcreteFactory<BtaPidTau31Pi0MergedSample> theBtaPidTau31Pi0MergedSampleFactory;

#endif




