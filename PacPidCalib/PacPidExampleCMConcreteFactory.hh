//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: $
//
// Description:
//     Adapted from BetaPidCalib/BtaExampleCMConcreteFactory.hh
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

#ifndef PACPIDEXAMPLECMCONCRETEFACTORY_HH
#define PACPIDEXAMPLECMCONCRETEFACTORY_HH

#include "BaBar/BaBar.hh"

#include "PacPidCalib/PacPidAbsPidCalibModule.hh"
#include "BetaEvent/BtaConcreteFactory.hh"

#include <string>

template <class Q>
class PacPidExampleCMConcreteFactory : 
  public BtaConcreteFactory<Q, PacPidAbsPidCalibModule>
{
public:

  // constructor
  PacPidExampleCMConcreteFactory() {}

  // virtual destructor
  virtual ~PacPidExampleCMConcreteFactory() {}

protected:
  virtual PacPidAbsPidCalibModule* newInstance(AppModule *aModule=0) const;
  virtual PacPidAbsPidCalibModule* newInstance(const char* const newName,
					    AppModule *aModule=0) const;
};

template <class Q>
PacPidAbsPidCalibModule*
PacPidExampleCMConcreteFactory<Q>::newInstance(AppModule *aModule) const{
  return BtaConcreteFactory<Q, PacPidAbsPidCalibModule>::newInstance(aModule);
}

template <class Q>
PacPidAbsPidCalibModule*
PacPidExampleCMConcreteFactory<Q>::newInstance(const char* const newName,
					 AppModule *aModule) const{
  std::string name(newName);
  if( aModule ) 
    return new Q( aModule, newName);
  else
    return new Q( ); 

}

// declare & instantiate

// NA: commented out
//#include "BetaPidCalib/BtaExampleCalibModule.hh"
//static PacPidExampleCMConcreteFactory<BtaExampleCalibModule> theBtaExampleCalibModuleFactory;
//#include "BetaPidCalib/BtaOnePartMicroModule.hh"
//static PacPidExampleCMConcreteFactory<BtaOnePartMicroModule> theBtaOnePartMicroModuleFactory;

#endif
