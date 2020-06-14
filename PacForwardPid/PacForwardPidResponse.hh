//-------------------------------------------------------------------------- 
// File and Version Information:
// 	$Id: $
//
// Name: 
//   PacForwardPidResponse 
// Description: 
//   Simulate the reconstruction of a PID response. 
// Environment: 
//   Software developed for PACRAT / SuperB 
// 
//  Copyright Information: 
//      Copyright (C) 2008    
//                (C) 2009      CNRS-IN2P3   
// 
//  Authors:   A.Berdyugin  2009/02/02
//             Nicolas ARNAUD 2009
//             Leonid BURMISTROV 2009
//----------------------------------------------------------------------- 

#ifndef PacForwardPidResponse_HH 
#define PacForwardPidResponse_HH

#include <vector> 


using std::vector;

class PacDetector; 
class GTrack; 
class PacSimTrack; 
class PacSimHit; 
class TrkRecoTrk;

class PacForwardPidResponse {

public:

  PacForwardPidResponse();
  virtual ~PacForwardPidResponse();

  double measurement() const { return( _measurement ); }
  double error() const { return( _error ); }
  double x() const { return( _x ); }
  double y() const { return( _y ); }
  double z() const { return( _z ); }
  double flight() const { return( _flight ); }

  void setMeasurement( double measurement ) {
    _measurement = measurement;
  }
  void setError( double error ) {
    _error = error;
  }

  void setX( double x ) { _x = x; }
  void setY( double y ) { _y = y; }
  void setZ( double z ) { _z = z; }
  void setFlight( double flight ) { _flight = flight; }
  
  // NA: not used so far
  //std::vector<PacForwardPidResponse*> makePidResponse(const PacSimTrack *simtrk) ;  

private:
  double _measurement;
  double _error;

  // position of the hit
  double _x;
  double _y;
  double _z;

  double _flight; 
};
#endif

