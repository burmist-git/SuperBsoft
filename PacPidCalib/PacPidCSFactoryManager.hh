//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description : 
//     Adapted from BetaPidCalib/BtaPidCSFactoryManager.hh
//      The Selector Factory Manager maintains a list a concrete
//      Selector Factories.  
//      The Factory Manager is a singleton class. 
//      A global function is provided to 'hide' the factory manager from
//      the user.  Here is an example:
//      Suppose that you have registered your Selector, say PidKillerSelector
//      (see the README file for details on how to do that).  
//      To use this selector in your analysis code :
//		#include "PidTools/BtaPidCSector.hh"
//		#include "PidTools/BtaPidCSFactoryManager.hh"
//		MyAnalysisModule::event(..)
//		{
//		  ...
//		  BtaPidCSector* mySelector = 
//		     newSelector("PidKillerSelector");
//		  mySelector->setParmValue("myCut",3.1416);
//			etc.
//              }
//
// Question :
//      How to make sure that the concrete factories are instanciated 
//      when the Factory Manager is invoked ?
// Present solution :
//      Every compiler insures that all the statics will be instanciated
//      before the first line of code where they are to be used.
//      The trick is to include the BtaPidCSConcreteFactory.hh file
//      in the BtaPidCSFactoryManager.cc file.  This is not really a dependance,
//      the Manager never really uses the factories themselves.  
//      It is not an inclusion in the header file, just in the .cc file.
//      It is just there to make sure that the static Concrete factories 
//      will be instanciated BEFORE the Manager.  
//      In fact what happens is the following :
//      The AppUserBuild is executed first.  In it there is one or several 
//      BtaPidCSectCand modules to be instanciated.  The constructor of 
//      BtaPidCSectCand invokes BtaPidCSFactoryManager::instance() which
//      will creates the unique instance of the Manager.  Thus the statics
//      in  BtaPidCSFactoryManager.cc MUST be instanciated.  The first concrete
//      factory to be instanciated registers itself, thereby trigging the
//      creation of the Manager.
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Nicolas ARNAUD (SuperB)
//      Gautier Hamel de Monchenault,   LBNL & Saclay
//
// Copyright Information:
//	 (C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------

#ifndef PACPIDCSFACTORYMANAGER_HH
#define PACPIDCSFACTORYMANAGER_HH

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------
#include "PacPidCalib/PacPidCalibSample.hh"
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

class PacPidCSFactoryManager : public BtaFactoryManager<PacPidCalibSample>
{
  
public:

  // constructor
  PacPidCSFactoryManager();

  // Destructor
  virtual ~PacPidCSFactoryManager( );
  
  // Additional member function to allow renaming of Samples
  // so that more than one can be run within a given Processor
  static PacPidCalibSample* get(const char* const sampleName, const char* const givenName,
				AppModule* aModule=0 );
};

//
// global function to hide the Factory manager to the user
// 
PacPidCalibSample* newPacPidCalibSample( const char* const str, AppModule* aModule=0 );

PacPidCalibSample* newPacPidCalibSample( const char* const sampleName,
					 const char* const givenName, AppModule* aModule=0 );

#endif // PACPIDCSFACTORYMANAGER_HH
