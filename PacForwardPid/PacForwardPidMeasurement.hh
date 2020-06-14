//-----------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// File and Version Information:
//     PacForwardPidMeasurement
//
//  Copyright Information:
//      Copyright (C) 2008      INFN Padova
//                (C) 2009      CNRS-IN2P3
//
//  Authors:  A.Berdyugin    2009/02/02
//            Nicolas ARNAUD 2009
//-----------------------------------------------------------------------

#ifndef PacForwardPidMeasurement_HH
#define PacForwardPidMeasurement_HH

#include "PacGeom/PacMeasurement.hh"
#include "PacForwardPid/PacForwardPidResponse.hh"
#include "PacSim/PacSimHit.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/HepPoint.h"

class PacForwardPidMeasurement : public PacMeasurement {
public:
  PacForwardPidMeasurement(); 
  virtual ~PacForwardPidMeasurement();

  virtual void createForwardPidHit( const PacSimHit& hit,
				    PacForwardPidResponse* res ) const;
  
private:
};

#endif
