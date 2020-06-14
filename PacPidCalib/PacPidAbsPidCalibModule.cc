//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description:
//     Adapted from BetaPidCalib/BtaAbsPidCalibModule.cc
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Nicolas ARNAUD (SuperB)
//      Luca Lista         21 Apr 1999 
//
// Copyright Information:
//	 (C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
//-----------------------
// This Class's Header --
//-----------------------
#include "PacPidCalib/PacPidAbsPidCalibModule.hh"

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
#include <vector>
using std::vector;

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include <string>
using std::string;
#include "AbsEvent/AbsEvent.hh"
#include "AssocTools/AstNamedMapVector.hh"
#include "AbsParm/AbsParmVector.hh"
#include "Beta/BtaCandidate.hh"
#include "CLHEP/String/Strings.h"
#include "ErrLogger/ErrLog.hh"
#include "Framework/AbsParmString.hh"
#include "PacPidCalib/PacPidCalibSample.hh"
#include "ProxyDict/Ifd.hh"
#include "ProxyDict/IfdStrKey.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//              ----------------------------------------
//              -- Public Function Member Definitions --
//              ----------------------------------------


//----------------
// Constructors --
//----------------



PacPidAbsPidCalibModule::PacPidAbsPidCalibModule( AppModule *callingModule,
					    const string &name ):
  _callingModule(callingModule),
  _myName(new string(name))
{
  string listName(*_myName);
  listName+="InputLists";
  _inputLists = new AbsParmVector<string>(listName.c_str(), callingModule);
  callingModule->commands()->append(_inputLists);
}

PacPidAbsPidCalibModule::PacPidAbsPidCalibModule()
{
  ErrMsg(fatal) << "Error!  Should not be calling default constructor of "
		<< " PacPidAbsPidCalibModule." << endmsg;
}

PacPidAbsPidCalibModule::~PacPidAbsPidCalibModule( )
{
  delete _myName;
  delete _inputLists;
}


void 
PacPidAbsPidCalibModule::setUp(){}

//-------------
// Modifiers --
//-------------

void 
PacPidAbsPidCalibModule::process(AbsEvent *anEvent){

  vector<string> samples(_inputLists->value());
  size_t i;
  for (i=0;i<samples.size();i++){
    // For now we assume the lists are independent!
    string inputListName(samples[i].c_str());
    IfdStrKey aKey(inputListName.c_str());
    
    AstNamedMapVector<BtaCandidate, BtaCandidate> *theMap =
      Ifd<AstNamedMapVector<BtaCandidate, BtaCandidate> >::get(anEvent, aKey);
  
    if (theMap != 0 ) {
      
      
      HepAList<BtaCandidate> *theList = Ifd<HepAList<BtaCandidate> >::get(anEvent, aKey);
      
      if (theList == 0 ) {
	
	if(ErrLogging(warning)) ErrMsg(warning) << "Failed to locate list of name " << inputListName
						<< " in the event. " << endmsg;
      } else {
	
	// We have a list and a Map.  Process it!
	
	HepAListIterator<BtaCandidate> iterator(*theList);
	
	this->processList(iterator, *theMap, anEvent);
      }
    }
  }
}

bool 
PacPidAbsPidCalibModule::operator==(const PacPidAbsPidCalibModule& other) const {
  return (moduleName() == other.moduleName());
}
//              -----------------------------------------------
//              -- Static Data & Function Member Definitions --
//              -----------------------------------------------

//              -------------------------------------------
//              -- Protected Function Member Definitions --
//              -------------------------------------------

void 
PacPidAbsPidCalibModule::addNewParm( const string& parmName, 
				  double parmValue, AppModule* aModule )
{
  this->addNewParmDouble(parmName, parmValue, _callingModule);
}

void 
PacPidAbsPidCalibModule::addNewParm( const string& parmName, 
			       int    parmValue, AppModule* aModule )
{
  this->addNewParmInt(parmName, parmValue, _callingModule);
}


void 
PacPidAbsPidCalibModule::addNewParm( const string& parmName, 
			       bool   parmValue, AppModule* aModule )
{
  this->addNewParmBool(parmName, parmValue, _callingModule);
}

void 
PacPidAbsPidCalibModule::addNewParm( const string& parmName, 
			       const HepString& parmValue, AppModule* aModule )
{
  this->addNewParmString(parmName, parmValue, _callingModule);
}

BtaParam<double>* 
PacPidAbsPidCalibModule::addNewParmDouble( const string& parmName, 
				     double parmValue, AppModule* aModule )
{
  checkModule( parmName, aModule ) ;
  return BtaParametrizable::addNewParmDouble(parmName, parmValue, _callingModule);
}

BtaParam<int>* 
PacPidAbsPidCalibModule::addNewParmInt( const string& parmName, 
				  int    parmValue, AppModule* aModule )
{
  checkModule( parmName, aModule ) ;
  return BtaParametrizable::addNewParmInt(parmName, parmValue, _callingModule);
}

BtaParam<bool>* 
PacPidAbsPidCalibModule::addNewParmBool( const string& parmName, 
				   bool   parmValue, AppModule* aModule ) 
{
  checkModule( parmName, aModule ) ;
  return BtaParametrizable::addNewParmBool(parmName, parmValue, _callingModule);
}

BtaParam<HepString>* 
PacPidAbsPidCalibModule::addNewParmString( const string& parmName, 
				     const HepString&   parmValue, AppModule* aModule )
{
  checkModule( parmName, aModule ) ;
  return BtaParametrizable::addNewParmString(parmName, parmValue, _callingModule);
}


void
PacPidAbsPidCalibModule::checkModule( const string& parmName, 
				AppModule* aModule ) 
{
  if (aModule !=0  && aModule != _callingModule){
    ErrMsg(error) << "Error.  Parameter " 
                  << parmName << " added to PacPidAbsPidCalibModule of name "
		  << moduleName() << " with different Framework module to original"
                  << " constructing module." << endmsg;
  }
}








