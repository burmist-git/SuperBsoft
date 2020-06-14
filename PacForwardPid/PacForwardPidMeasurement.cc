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

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "PacForwardPid/PacForwardPidMeasurement.hh"

// NA Not used for now
//#include "PacForwardPid/PacForwardPidHit.hh"

using namespace std;

PacForwardPidMeasurement::PacForwardPidMeasurement() :
  PacMeasurement( PacMeasurement::ForwardPID ) 
{
    cout << "Creating PacForwardPidMeasurement\n"; 
}

PacForwardPidMeasurement::~PacForwardPidMeasurement(){}

void PacForwardPidMeasurement::createForwardPidHit( const PacSimHit& hit, 
						    PacForwardPidResponse* res) const {
  
  /*
  //   std::cout << "I'm in PacForwardPidMeasurement::createPidHits\n";
 
  // Create the Pid hit from the PacSimHit
  PacForwardPidHit* forwardpidhit = new PacForwardPidHit( hit );
  */

  // NA
  res->setMeasurement( hit.time() ); // time measurement for now
  res->setError( 20.0E-12 ); // assuming timing error of 20 ps

  res->setX( hit.position().x() );
  res->setY( hit.position().y() );
  res->setZ( hit.position().z() );

  res->setFlight( hit.detIntersection().pathlen );
}

