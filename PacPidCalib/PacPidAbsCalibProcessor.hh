//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description:
//     Adapted from BetaPidCalib/BtaAbsCalibProcessor.hh
//	Class PacPidAbsCalibProcessor. This is a simple example of a user module. It
//	just prints out each entrypoint (operation) as it is accessed.
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Nicolas ARNAUD (SuperB)
//	David R. Quarrie		Original Author
//
// Copyright Information:
//	 (C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------

#ifndef PACPIDABSCALIBPROCESSOR_HH
#define PACPIDABSCALIBPROCESSOR_HH

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppFilterModule.hh"
#include "TagData/TagAccessor.hh"
#include <string>
#include <vector>
//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class PacPidAbsPidCalibModule;
class PacPidCalibSample;

//		---------------------
// 		-- Class Interface --
//		---------------------
class PacPidAbsCalibProcessorData;
 
class PacPidAbsCalibProcessor : public AppFilterModule,
			     public TagAccessor
{

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  PacPidAbsCalibProcessor( const char* const theName, 
			const char* const theDescription );

  // Destructor
  virtual ~PacPidAbsCalibProcessor( );

  // Operations

  virtual AppResult beginJob ( AbsEvent* anEvent );
  virtual AppResult event    ( AbsEvent* anEvent );
  virtual AppResult endJob   ( AbsEvent* anEvent );    
  virtual AppResult abortJob ( AbsEvent* anEvent );
    
  virtual AppModule* clone(const char* cloneName);

protected:

  virtual std::string processorName() const=0;


  // Modules will generally be very XxxData dependent.  No default
  /// implementation since a factory is expected in BetaXxxPidCalib

  virtual PacPidAbsPidCalibModule *generateCalibModule(const std::string &);
  
  // This is not PV to ensure backwards compatibility
  virtual PacPidAbsPidCalibModule *generateCalibModule(const std::string &,
						    const std::string &);


  // Samples will generally not be very XxxData dependent.  Default 
  // implementation is to use the BtaPidCSConcreteFactory.

  virtual PacPidCalibSample *generateCalibSample(const std::string &) ;
  virtual PacPidCalibSample *generateCalibSample(const std::string &,
						 const std::string &) ;

  const std::vector<PacPidCalibSample*> &executedSamples() const;
  const std::vector<PacPidAbsPidCalibModule*> &executedModules() const;

  const std::vector<PacPidCalibSample*> &knownSamples() const;
  const std::vector<PacPidAbsPidCalibModule*> &knownModules() const;
  
  const std::vector<std::string*> &eventListNames() const;

  // Register samples and modules with default names
  bool registerSample(const std::string &sampleName) ;
  bool registerModule(const std::string &moduleName);
  
  // Register samples and modules with given names 
  bool registerSample(const std::string &sampleName, const std::string &newSampleName) ;
  bool registerModule(const std::string &moduleName, const std::string &newModuleName);

private:
  

  PacPidAbsCalibProcessorData *_privateData;

};

#endif
