//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description:
//	Class PacPidExampleCalibProcessor. This is a simple example of a user module. It
//	just prints out each entrypoint (operation) as it is accessed.
//     Adapted from BetaPidCalib/BtaExampleCalibProcessor.hh
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Nicolas ARNAUD (SuperB)
//      David R. Quarrie		Original Author
//
// Copyright Information:
//	 (C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------

#ifndef PACPIDEXAMPLECALIBPROCESSOR_HH
#define PACPIDEXAMPLECALIBPROCESSOR_HH

//----------------------
// Base Class Headers --
//----------------------
#include "PacPidCalib/PacPidAbsCalibProcessor.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

#include <string>
class PacPidAbsPidCalibModule;
class PacPidCalibSample;

//		---------------------
// 		-- Class Interface --
//		---------------------
class PacPidExampleCalibProcessorData;
 
class PacPidExampleCalibProcessor : public PacPidAbsCalibProcessor {

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  PacPidExampleCalibProcessor( const char* const theName, 
			const char* const theDescription );

  // Destructor
  virtual ~PacPidExampleCalibProcessor( );

  // Operations


  virtual AppModule* clone(const char* cloneName);

protected:
  virtual std::string processorName() const;

private:
  
  PacPidExampleCalibProcessorData *_privateData;

};

#endif
