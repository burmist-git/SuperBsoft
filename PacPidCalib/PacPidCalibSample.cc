//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: $
//
// Description:
//     Adapted from BetaPidCalib/BtaPidCalibSample.cc
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Nicolas ARNAUD (SuperB)
//      Luca Lista          2 May 1999
//
// Copyright Information:
//	 (C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "PacPidCalib/PacPidCalibSample.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

#include <string>
using std::string;
#include <set>


//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "AssocTools/AstNamedMapVector.hh"
#include "Beta/BtaCandidate.hh"
#include "CLHEP/String/Strings.h"
#include "ErrLogger/ErrLog.hh"
#include "BbrStdUtils/String.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//              ----------------------------------------
//              -- Public Function Member Definitions --
//              ----------------------------------------


//----------------
// Constructors --
//----------------

PacPidCalibSample::PacPidCalibSample(AppModule *theModule, 
				     const string &name) :
  //  _map ( new Map ( &babar::String::rwHash, 5 ) ),
  _map ( new Map () ),
  _mapStringsGarbageCollection(new std::vector<string*>()),
  _callingModule(theModule),
  _myName(new string(name)),
  _filterMap(new std::set<string>())
{
}

PacPidCalibSample::PacPidCalibSample()
{
  ErrMsg(fatal) << "Error!  Should not be calling default constructor of PacPidCalibSample"
		<< endmsg;
}

//--------------
// Destructor --
//--------------

PacPidCalibSample::~PacPidCalibSample()
{
  delete _mapStringsGarbageCollection;

  delete _filterMap;

  delete _map;

  delete _myName;
}

//-------------
// Methods   --
//-------------

bool 
PacPidCalibSample::sampleUsedForFilter(const string &theString)const {

  if (_filterMap->find(theString) != _filterMap->end() )
    return true;
  else 
    return false;
  //  return _filterMap->contains(theString);
}

void PacPidCalibSample::clearMap()
{
  //  _mapStringsGarbageCollection->clearAndDestroy();  
  for_each(_mapStringsGarbageCollection->begin(),_mapStringsGarbageCollection->end(),
	   DeleteObject());

  _mapStringsGarbageCollection->clear();

  _map->clear();

  // We do not need to clear the filterMap.
  // All this does is to say whether or not a 
  // sample is being used for filtering.
  // This logic does not change from event to event, 
  // so we can keep the status of the map 
  // untill the end of processing.
  // _filterMap->clear();

}

void 
PacPidCalibSample::setUp(){}

    
//-------------
// Operators --
//-------------

bool 
PacPidCalibSample::operator==(const PacPidCalibSample& other)const {
  return (*_myName == *other._myName);
}

//-------------
// Selectors --
//-------------
PacPidCalibSample::MapIterator 
PacPidCalibSample::iterator() const
{
  //  return MapIterator( * _map );
  return _map->begin();
}

PacPidCalibSample::MapIterator 
PacPidCalibSample::beginIterator() const
{
  return _map->begin();
}

PacPidCalibSample::MapIterator 
PacPidCalibSample::endIterator() const
{
  return _map->end();
}
    
//-------------
// Modifiers --
//-------------

//              -----------------------------------------------
//              -- Static Data & Function Member Definitions --
//              -----------------------------------------------


//              -------------------------------------------
//              -- Protected Function Member Definitions --
//              -------------------------------------------

void 
PacPidCalibSample::selectThis(BtaCandidate *theCand, const string &theMapName,
			      bool useForFilter){

  // Check to see if BtaCandidate pointer is not null.
  if (theCand != 0) {

    Map::iterator iter=_map->find((std::string*)&theMapName);


    AstNamedMapVector<BtaCandidate, BtaCandidate> *theNV;
	//theNV =  _map->findValue(&theMapName);
    //    if (theNV == 0){
    if (iter == _map->end()){
      
      string *stringPtr = new string(theMapName);
      _mapStringsGarbageCollection->push_back(stringPtr);
      theNV = new AstNamedMapVector<BtaCandidate, BtaCandidate>(theMapName.c_str());
      (*_map)[stringPtr] = theNV;
    } else {
      theNV = iter->second;
    }
    
    BtaCandidate *foundCand = new BtaCandidate(*theCand);
    theNV->append(theCand, foundCand);
    
    if (useForFilter) _filterMap->insert(theMapName);
  }
}

void 
PacPidCalibSample::addNewParm( const string& parmName, 
			       double parmValue, AppModule* aModule )
{
  this->addNewParmDouble(parmName, parmValue, _callingModule);
}

void 
PacPidCalibSample::addNewParm( const string& parmName, 
			       int    parmValue, AppModule* aModule )
{
  this->addNewParmInt(parmName, parmValue, _callingModule);
}


void 
PacPidCalibSample::addNewParm( const string& parmName, 
			       bool   parmValue, AppModule* aModule )
{
  this->addNewParmBool(parmName, parmValue, _callingModule);
}

void 
PacPidCalibSample::addNewParm( const string& parmName, 
			       const HepString& parmValue, AppModule* aModule )
{
  this->addNewParmString(parmName, parmValue, _callingModule);
}

BtaParam<double>* 
PacPidCalibSample::addNewParmDouble( const string& parmName, 
				     double parmValue, AppModule* aModule )
{
  checkModule( parmName, aModule ) ;
  return BtaParametrizable::addNewParmDouble(parmName, parmValue, _callingModule);
}

BtaParam<int>* 
PacPidCalibSample::addNewParmInt( const string& parmName, 
				  int    parmValue, AppModule* aModule )
{
  checkModule( parmName, aModule ) ;
  return BtaParametrizable::addNewParmInt(parmName, parmValue, _callingModule);
}

BtaParam<bool>* 
PacPidCalibSample::addNewParmBool( const string& parmName, 
				   bool   parmValue, AppModule* aModule ) 
{
  checkModule( parmName, aModule ) ;
  return BtaParametrizable::addNewParmBool(parmName, parmValue, _callingModule);
}

BtaParam<HepString>* 
PacPidCalibSample::addNewParmString( const string& parmName, 
				     const HepString&   parmValue, AppModule* aModule )
{
  checkModule( parmName, aModule ) ;
  return BtaParametrizable::addNewParmString(parmName, parmValue, _callingModule);
}


//              -----------------------------------------
//              -- Private Function Member Definitions --
//              -----------------------------------------
void
PacPidCalibSample::checkModule( const string& parmName, 
				AppModule* aModule ) 
{
  if (aModule !=0  && aModule != _callingModule){
    ErrMsg(error) << "Error.  Parameter " 
		  << parmName << " added to PacPidCalibSample of name "
		  << sampleName() << " with different Framework module to original"
		  << " constructing module." << endmsg;
  }
}
//              -----------------------------------
//              -- Internal Function Definitions --
//              -----------------------------------

/* static */ unsigned 
PacPidCalibSample::hashFunction(const PacPidCalibSample &theSample){
  return babar::String::rwHash(theSample.sampleName());
}




