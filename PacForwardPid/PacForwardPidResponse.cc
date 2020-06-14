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
//  Authors:  A.Berdyugin   2009/02/02
//            Nicolas ARNAUD 2009
//-----------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "PacForwardPid/PacForwardPidResponse.hh"
#include "PacGeom/PacMeasurement.hh"
#include "PacSim/PacSimTrack.hh"
#include "PacGeom/PacDetElem.hh"
#include "TrkBase/TrkRecoTrk.hh"


#include <iostream>
using namespace std;

PacForwardPidResponse::PacForwardPidResponse() :
    _measurement( 0.0 )
  , _error( 1E10 )
  , _x( 0.0 )
  , _y( 0.0 )
  , _z( -9999 )
  , _flight( 0.0 )
{ }

PacForwardPidResponse::~PacForwardPidResponse() { }

// NA: not used so far
/*
std::vector<PacForwardPidResponse*> PacForwardPidResponse::makePidResponse(const PacSimTrack *simtrk)  {
    
    bool foundPid(false);
    bool firstHit(true);

    // Basic sanity checks

    std::vector<PacForwardPidResponse*> list;

    if(simtrk != 0 && simtrk->getGTrack() != 0 && simtrk->getTraj() != 0 ) {
        
        const std::vector<PacSimHit> hitlist= simtrk->getHitList();
        for ( std::vector<PacSimHit>::const_iterator ihit= hitlist.begin();
              ihit!= hitlist.end(); ihit++) {
            
            const DetIntersection& dinter = ihit->detIntersection();
            const PacDetElem* pelem = dynamic_cast<const PacDetElem *>(dinter.delem);
            if( pelem != 0 && pelem->measurement()!= 0 &&
                pelem->measurement()->measurementType()==PacMeasurement::ForwardPID &&
                firstHit ) {
                
                foundPid = true;
                firstHit = false;
               
                // create a PID response and append it to the list of PacForwardPidResponse objects

                PacForwardPidResponse *pidres = new PacForwardPidResponse(); 
                pelem->measurement()->createForwardPidHit(*ihit,pidres);
                list.push_back(pidres);
            }
            
        }
    }
    
    return list;
    
}
*/
