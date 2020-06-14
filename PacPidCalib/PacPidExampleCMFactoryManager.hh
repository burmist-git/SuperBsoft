//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description : 
//     Adapted from BetaPidCalib/BtaExampleCMFactoryManager.hh
//      The Selector Factory Manager maintains a list a concrete
//      Selector Factories.  
//      The Factory Manager is a singleton class. 
//      A global function is provided to 'hide' the factory manager from
//      the user.  Here is an example:
//      Suppose that you have registered your Selector, say PidKillerSelector
//      (see the README file for details on how to do that).  
//      To use this selector in your analysis code :
//		#include "BtaExample/BtaAbsPidCalibModule.hh"
//		#include "BtaExample/BtaExampleCMFactoryManager.hh"
//		MyAnalysisModule::event(..)
//		{
//		  ...
//		  BtaAbsPidCalibModule* myModule = 
//		     newPacPidExampleCalibModule("PidKillerModule");
//		  mySelector->setParmValue("myCut",3.1416);
//			etc.
//              }
//
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Nicolas ARNAUD (SuperB)
//      Gautier Hamel de Monchenault,   LBNL & Saclay    Originator
//
// Copyright Information:
//	 (C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------

#ifndef PACPIDEXAMPLECMFACTORYMANAGER_HH
#define PACPIDEXAMPLECMFACTORYMANAGER_HH

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------
#include "PacPidCalib/PacPidAbsPidCalibModule.hh"
#include "BetaEvent/BtaFactoryManager.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------

class PacPidExampleCMFactoryManager : public BtaFactoryManager<PacPidAbsPidCalibModule>
{
  
public:

  // constructor
  PacPidExampleCMFactoryManager();

  // Destructor
  virtual ~PacPidExampleCMFactoryManager( );

  static PacPidAbsPidCalibModule* get(const char* const moduleName, 
				   const char* const givenName,
				   AppModule* aModule=0 );
};

//
// global function to hide the Factory manager to the user
// 
PacPidAbsPidCalibModule* newPacPidExampleCalibModule( const char* const str, 
					       AppModule* aModule=0 );

PacPidAbsPidCalibModule* newPacPidExampleCalibModule( const char* const moduleName,
						const char* const givenName,
					       AppModule* aModule=0 );

#endif // PACPIDEXAMPLECMFACTORYMANAGER_HH
