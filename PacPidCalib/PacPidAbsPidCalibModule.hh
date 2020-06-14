//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description:
//     Adapted from BetaPidCalib/BtaAbsPidCalibModule.hh
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
#ifndef PACPIDABSPIDCALIBMODULE_HH
#define PACPIDABSPIDCALIBMODULE_HH

//---------------
// C++ Headers --
//---------------
#include <string>
//----------------------
// Base Class Headers --
//----------------------

#include "BetaEvent/BtaParametrizable.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class PacPidCalibSample;
class AbsParmString;
class BtaCandidate;
#include <string>
template<class T> class HepAListIterator;
template<class T1, class T2> class AstNamedMapVector;
template<class T> class AbsParmVector;
template <class T> class BtaParam;


//              ---------------------
//              -- Class Interface --
//              ---------------------

class PacPidAbsPidCalibModule:  public BtaParametrizable
{

public:

  PacPidAbsPidCalibModule( AppModule *, const std::string &);

  // Default constructor needed by BetaEvent factory 
  // pattern, but is not used (error flagged on execution).
  PacPidAbsPidCalibModule();

  virtual ~PacPidAbsPidCalibModule( );
  
  virtual void setUp();  // Empty default implementation

  void process( AbsEvent *);

  const std::string &moduleName() const {return *_myName;}
  

  bool operator==(const PacPidAbsPidCalibModule &) const;

protected:

  virtual void processList( HepAListIterator<BtaCandidate>&,
			    const AstNamedMapVector<BtaCandidate, BtaCandidate> &,
			    AbsEvent *) = 0;

  // Hmm....this is beginning to have so much in common with
  // BtaPidCalibSample that I wonder 
  // Override the addNewParm functions to insist that the 
  // parameters belong to a module, sepcifically the calling module.

  // These four should stay here only during migration 
  void addNewParm( const std::string& parmName, 
                   double parmValue, AppModule* aModule=0 );
  void addNewParm( const std::string& parmName, 
                   int    parmValue, AppModule* aModule=0 );
  void addNewParm( const std::string& parmName, 
                   bool   parmValue, AppModule* aModule=0 );
  void addNewParm( const std::string& parmName, 
                   const HepString& parmValue, AppModule* aModule=0 );


  BtaParam<double>* addNewParmDouble( const std::string& parmName, 
				      double parmValue, AppModule* aModule=0 );
  BtaParam<int>* addNewParmInt( const std::string& parmName, 
				int    parmValue, AppModule* aModule=0 );
  BtaParam<bool>* addNewParmBool( const std::string& parmName, 
				  bool   parmValue, AppModule* aModule=0 );
  BtaParam<HepString>* addNewParmString( const std::string& parmName, 
					 const HepString&   parmValue, AppModule* aModule=0 );

  AppModule *myModule() const { return _callingModule;}

private:
  
  AppModule *_callingModule;

  std::string *_myName;

  AbsParmVector<std::string> *_inputLists;


  // check that the module is the same as _callingModule
  void checkModule( const std::string& parmName, AppModule* aModule ) ; 

};

#endif









