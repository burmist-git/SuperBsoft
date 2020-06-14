//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description:
//	Class PacPidExampleCalibProcessor. This is a simple example of a user module. It
//	just prints out each entrypoint (operation) as it is accessed.
//     Adapted from BetaPidCalib/BtaExampleCalibProcessor.cc
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Nicolas ARNAUD (SuperB)
//      David R. Quarrie		Original Author
//
// Copyright Information:
//	 (C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "Experiment/Experiment.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "PacPidCalib/PacPidExampleCalibProcessor.hh"

//-------------
// C Headers --
//-------------
#include <assert.h>

//---------------
// C++ Headers --
//---------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "AbsEvent/AbsEvent.hh"
#include "PacPidCalib/PacPidAbsPidCalibModule.hh"
#include "PacPidCalib/PacPidExampleCMFactoryManager.hh"
#include "PacPidCalib/PacPidCSFactoryManager.hh"
#include "ErrLogger/ErrLog.hh"
#include <string>
using std::string;

struct PacPidExampleCalibProcessorData{
};

//----------------
// Constructors --
//----------------

PacPidExampleCalibProcessor::PacPidExampleCalibProcessor( const char* const theName, 
						    const char* const theDescription )
  : PacPidAbsCalibProcessor( theName, theDescription ),
    _privateData(new PacPidExampleCalibProcessorData())
{
  // We register all possible samples and modules this Processor could 
  // want to use.   These are switched on and off using the moduleNames
  // and sampleNames tcl lists in the base class.  

  // NB these are _not_ tcl parameters. 

  // NA
  string sampleName;
  sampleName = "PacPidDstarSample";
  registerSample(sampleName);

  //string moduleName("BtaExampleCalibModule");
  //registerModule(moduleName);

  // MC 
  //string sampleName("BtaPidMCCalibSample");
  //registerSample(sampleName);
  
  // Various electrons
  //sampleName = "BtaEmcBhabhaSample";
  //registerSample(sampleName);
 
  //sampleName = "BtaEmcRadBhabhaSample";
  //registerSample(sampleName);

  //sampleName = "BtaTrkBhabhaSample";
  //registerSample(sampleName);

  //sampleName = "BtaPidBhabhaSample";
  //registerSample(sampleName);

  //sampleName = "BtaKinIFRBhabhaSample";
  //registerSample(sampleName);

  //sampleName = "BtaPideeeeCalibSample";
  //registerSample(sampleName);

  //sampleName = "BtaPidV0SampleModule";
  //registerSample(sampleName);

  //sampleName = "BtaPideeeeCalibSample";
  //registerSample(sampleName);

  //sampleName = "BtaPidLambdaProtonSample";
  //registerSample(sampleName);

  //sampleName = "BtaPidDstarSample";
  //registerSample(sampleName);

  //sampleName = "BtaPidGammaeeSample";
  //registerSample(sampleName);

  //sampleName = "BtaPidKsSample";
  //registerSample(sampleName);

  //sampleName = "BtaTau31Sample";
  //registerSample(sampleName);

  //sampleName = "BtaLambdaTrkOnlySample";
  //registerSample(sampleName);

  //sampleName = "BtaeemumuSample";
  //registerSample(sampleName);

  //sampleName = "BtamumugammaSample";
  //registerSample(sampleName);

  //sampleName = "Btamumugamma2Sample";
  //registerSample(sampleName);

  //sampleName = "BtaVcsSample";
  //registerSample(sampleName);
}


//--------------
// Destructor --
//--------------

PacPidExampleCalibProcessor::~PacPidExampleCalibProcessor( )
{
  delete _privateData;
}

//--------------
// Operations --
//--------------


AppModule*
PacPidExampleCalibProcessor::clone(const char* cloneName){
  return new PacPidExampleCalibProcessor(cloneName,"this module is a clone of user"); 
}

string 
PacPidExampleCalibProcessor::processorName() const{
  static string theName("PacPidExampleCalibProcessor");
  return theName;
}
