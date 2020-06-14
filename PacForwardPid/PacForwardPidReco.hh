//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description:
//     Simulates the physics response of the Forward PID to
//     a charged track
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

#ifndef PACFORWARDPIDRECO_HH
#define PACFORWARDPIDRECO_HH

class PacDetector;
class GTrack;
class PacSimTrack;
class PacSimHit;
class TrkRecoTrk;
class PacForwardPidResponse;

class PacForwardPidReco {
public:
  // Constructor
  PacForwardPidReco();
  // Destructor
  ~PacForwardPidReco();
  PacForwardPidResponse* findForwardPidHits( const PacSimTrack* simtrk ) const;
private:
};

#endif // PACFORWARDPIDRECO_HH
