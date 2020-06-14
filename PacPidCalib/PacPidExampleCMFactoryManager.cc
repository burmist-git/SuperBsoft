//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description:
//       See the .hh file for details
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
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "PacPidCalib/PacPidExampleCMFactoryManager.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

// The following include is necessary in order to make sure that
// all the concrete factories are instantiated when the first Selection 
// module is created- Do NOT remove it (eventhough it is not a dependancy)-  
// Gautier-
// Instantiation of static concrete factories
#include "PacPidCalib/PacPidExampleCMConcreteFactory.hh" 

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//----------------
// Constructors --
//----------------
PacPidExampleCMFactoryManager::PacPidExampleCMFactoryManager()
{
}

//--------------
// Destructor --
//--------------

PacPidExampleCMFactoryManager::~PacPidExampleCMFactoryManager()
{
}

PacPidAbsPidCalibModule*
newPacPidExampleCalibModule( const char* const str,  AppModule* aModule )
{
  return BtaFactoryManager<PacPidAbsPidCalibModule>::get( str, aModule );
}

PacPidAbsPidCalibModule*
newPacPidExampleCalibModule( const char* const moduleName,  
			  const char* const givenName,
			  AppModule* aModule )
{
  return PacPidExampleCMFactoryManager::get( moduleName, givenName, aModule );
}

PacPidAbsPidCalibModule*
PacPidExampleCMFactoryManager::get(const char* const sampleName,  
				const char* const givenName, AppModule* aModule){
  PacPidAbsPidCalibModule* t;
  for( int i=0; i<list().size(); i++ ){
    if( ( t=list()[i]->create( sampleName, givenName, aModule ) ) ) return t;
  }
  return (PacPidAbsPidCalibModule*)0;

}
