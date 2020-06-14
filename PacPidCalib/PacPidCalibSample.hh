//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description:
//     Adapted from BetaPidCalib/BtaPidCalibSample.hh
//      Pid Calibration sample. Essentially a container of 
//      lists of BtaCandidates stored by key.
//      The concrete way of filling the lists is left to
//      the conc#ifndef BTAPIDCALIBSAMPLE_HHrete class implementation.
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
#ifndef PACPIDCALIBSAMPLE_HH
#define PACPIDCALIBSAMPLE_HH

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------

#include "BetaEvent/BtaParametrizable.hh"

#include "BbrStdUtils/String.hh"
#include "BbrStdUtils/BbrCollectionUtils.hh"
#include "BbrStdUtils/CollectionUtils.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class AbsEvent;
class AbsEventTag;
#include <string>
#include <set>
#include <map>
#include <vector>

template <class T1, class T2> class AstNamedMapVector;
class BtaCandidate;

//              ---------------------
//              -- Class Interface --
//              ---------------------


class PacPidCalibSample: public BtaParametrizable{

public:
  
  // Constructor

  // Default constructor does nothing except print a warning message.
  // Should only exist within the context of a module (specifically 
  // an instance of BtaCalibProcessor

  PacPidCalibSample();
  
  PacPidCalibSample(AppModule *, const std::string &);

  // Destructor 

  virtual ~PacPidCalibSample();
 
  // setUp acts like beginJob.  Not made pure virtual for
  // backward compatibility during development phase. 
  // Base class implementation does nothing.  *DON'T* 
  // add new params here.  They should go in the constructor.

  virtual void setUp();
  
  // passesTagSelection is designed for making a first pass
  // using tag bits.   Default implementation is to set true
  virtual bool passesTagSelection( const AbsEventTag *,
				   AbsEvent *)=0;

  // The following pure virtual function is responsible for
  // creating the lists of candidates.  For each candidate
  // selected, call the selectThis() function.
  // *All* book keeping (maps, lists, MC, memory management)
  // will be done inside this function.  Use the original
  // candidate - cloning will be taken care of for you. 

  virtual void createCandidateLists( AbsEvent* ) = 0;

  const std::string & sampleName() const {return *_myName;}

  bool sampleUsedForFilter(const std::string &) const;

  void clearMap();

  bool operator==(const PacPidCalibSample &) const;

  typedef  std::map<std::string *,
    AstNamedMapVector<BtaCandidate, BtaCandidate> *, babar::Collection::PtrLess >::iterator
  MapIterator;

  MapIterator iterator() const;
  MapIterator beginIterator() const;
  MapIterator endIterator()   const;

  static unsigned hashFunction(const PacPidCalibSample &);

protected:

  // The first argument is the candidate you have selected
  // The second is the list/map name that it should be appended to.
  // The last argument can be set to false if you are selecting 
  // candidates which are _not_ useful for the control sample
  // but are useful for cross checks. etc. 
  // WARNING for micro filtering modes:  if the selectThis function
  // is called then the Sample is taken as having selected the event
  // unless the useForFilter is set to false. 

  void selectThis(BtaCandidate *, const std::string &, bool useForFilter=true);

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

  typedef std::map<std::string *,
    AstNamedMapVector<BtaCandidate, BtaCandidate> *, babar::Collection::PtrLess > Map;

  Map * _map;

  std::vector<std::string *> *_mapStringsGarbageCollection;

private:
  
  AppModule *_callingModule;

  std::string *_myName;

  std::set<std::string> * _filterMap;

  // check that the module is the same as _callingModule
  void checkModule( const std::string& parmName, AppModule* aModule ) ; 

};

#endif
