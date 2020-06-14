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

#include "BaBar/BaBar.hh"

#include "PacForwardPid/PacForwardPidReco.hh"
#include "PacForwardPid/PacForwardPidResponse.hh"

#include "PacSim/PacSimTrack.hh"
#include "PacGeom/PacDetElem.hh"
#include "PacGeom/PacMeasurement.hh"
#include "TrkBase/TrkRecoTrk.hh"

// Remove when PacSimTrack has fuller functionality. 
//#include "G3Data/GTrack.hh"
#include "PDT/PdtEntry.hh"

// Constructor
PacForwardPidReco::PacForwardPidReco() {}

// Destructor
PacForwardPidReco::~PacForwardPidReco() {}


PacForwardPidResponse*
PacForwardPidReco::findForwardPidHits( const PacSimTrack* simtrk ) const
{
  PacForwardPidResponse* res( 0 );

  // Safety checks copied from PacDircHitFinder
  if ( !simtrk ) return( res ); 
  if ( !simtrk->getGTrack() ) return( res );
  if ( 0 == simtrk->getGTrack()->pdt()->charge() ) return( res );

  // Loop over the simulated hits
  const std::vector<PacSimHit>& simhits = simtrk->getHitList();
  for ( std::vector<PacSimHit>::const_iterator ihit = simhits.begin(); 
	ihit != simhits.end(); ++ihit ) {
    const DetIntersection& dinter = ihit->detIntersection();
    const PacDetElem* pelem = dynamic_cast<const PacDetElem*>( dinter.delem );
    if ( !pelem ) continue;
    if ( !pelem->measurement() ) continue;
    if ( pelem->measurement()->measurementType() != 
	 PacMeasurement::ForwardPID) continue;
    if ( !pelem->activeRegion( ihit->position() ) ) continue; 

    res = new PacForwardPidResponse();
    pelem->measurement()->createForwardPidHit( *ihit, res );
    //cout << "Hit in the forward PID region" << endl;
    //ihit->print( cout );
  }

  return( res );
}
