//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: $
//
// Description:
//      AppModule filling the PID 'monitoring' ntuple(s)
//
// Environment:
//	Software developed for the Super B project
//
// Author List:
//      Nicolas ARNAUD
//
// Copyright Information:
//	 (C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"

//-------------
// C++ Headers:
//-------------
#include <signal.h>

//-----------------------
// This Class's Header --
//-----------------------
#include "Framework/AppFramework.hh"
#include "Framework/AppUserBuild.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "ErrLogger/ErrLog.hh"

// actions
#include "PacMC/PacActions.hh"

// simulation
#include "PacMC/PmcSimulationSequence.hh"

// physics sequence
#include "PacMC/PmcPhysicsSequence.hh"

// modules
#include "PacPidCalib/PacPidExampleCalibProcessor.hh"
#include "PacPidCalib/PacPidPCNtupleWriter.hh"
#include "PacPidCalib/PacPidPCNtupleWriterMC.hh"
#include "EmcPid/EmcLoadPid.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//----------------
// Constructors --
//----------------

AppUserBuild::AppUserBuild( AppFramework* theFramework )
    : AppBuild( theFramework )
{

  //  main simulation sequence
  PmcSimulationSequence(this);

  // minimal physics sequence
  PmcPhysicsSequence(this);

  // actions
  PacActions(theFramework);
 
  add( new EmcLoadPid( "EmcLoadPid",
		       "Loads environment with Emc PID factory" ) );

  add( new PacPidExampleCalibProcessor( "PacPidExampleCalibProcessor",
					"Example calibration sample processor" ) );

  add( new PacPidPCNtupleWriter( "PacPidPCNtupleWriter",
				 "Writing ntuples from selected sample" ) );

  add( new PacPidPCNtupleWriterMC( "PacPidPCNtupleWriterMC",
				   "Write MC ntuples for e,mu,pi,K,p using Truth information" ) );
}


//--------------
// Destructor --
//--------------

AppUserBuild::~AppUserBuild()
{
}



  
