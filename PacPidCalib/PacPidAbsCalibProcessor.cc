//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: $
//
// Description:
//     Adapted from BetaPidCalib/BtaAbsCalibProcessor.cc
//	Class PacPidAbsCalibProcessor. This is a simple example of a user module. It
//	just prints out each entrypoint (operation) as it is accessed.
//
// Environment:
//	Software developed for the Super B project
//  adapted from software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Nicolas ARNAUD (SuperB)
//	David R. Quarrie		Original Author
//
// Copyright Information:
//	 (C) 2009 CNRS-IN2P3
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "Experiment/Experiment.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "PacPidCalib/PacPidAbsCalibProcessor.hh"

//-------------
// C Headers --
//-------------
#include <assert.h>

//---------------
// C++ Headers --
//---------------
#include <iostream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "AbsEvent/AbsEvent.hh"
#include "AbsEvent/AbsEventID.hh"
#include "AbsEvent/getTmpAList.hh"
#include "AbsEventTag/AbsEventTag.hh"
#include "AbsParm/AbsParmVector.hh"
#include "Beta/EventInfo.hh"
#include "AssocTools/AstNamedMapVector.hh"
#include "PacPidCalib/PacPidCalibSample.hh"
#include "PacPidCalib/PacPidAbsPidCalibModule.hh"
#include "PacPidCalib/PacPidCSFactoryManager.hh"
#include "PacPidCalib/PacPidExampleCMFactoryManager.hh"
#include "Beta/BtaCandidate.hh"
#include "CommonUtils/ComTemplatedBitMap.hh"
#include "ErrLogger/ErrLog.hh"
#include "Framework/AbsParmBool.hh"
#include "AbsParm/AbsParmIfdStrKey.hh"
#include "Framework/AbsParmString.hh"

#include "BbrStdUtils/String.hh"
#include "BbrStdUtils/BbrCollectionUtils.hh"

#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
#include <set>
using std::cout;
using std::endl;
using std::ostream;

struct PacPidAbsCalibProcessorData{
  AbsParmBool *tagFilterMode;
  AbsParmBool *tagBitFilterMode;
  AbsParmBool *microFilterMode;
  AbsParmBool *makeTagBits;
  AbsParmBool *putListsInEvent;
  AbsParmString *listNameSuffix;
  AbsParmVector<string> *sampleNames;
  std::map<string*, PacPidCalibSample*, babar::Collection::PtrLess> *sampleDict;
  std::map<string, int> *sampleStats;
  vector<PacPidCalibSample*> *executedSamples;
  vector<PacPidCalibSample*> *knownSamples;
  AbsParmVector<string> *moduleNames;
  std::map<string*,  PacPidAbsPidCalibModule*, babar::Collection::PtrLess> *moduleDict;
  vector<PacPidAbsPidCalibModule*> *executedModules;
  vector<PacPidAbsPidCalibModule*> *knownModules;
  vector<string*> *eventListNames;
  std::set<string*, BbrPtrLess> *allKnownListNames;
  // Looks like this is not being used anyway... don't migrate it to STL (VB - 05/14/02)
  //  RWTValHashDictionary<PacPidCalibSample*, int> *sampleIndex;

  //  ComTemplatedBitMap<PacPidCalibSample> *sampleMap;
  unsigned long compositeTagWord;

  int numEventsSeen;
  int numEventsPassed;
};

unsigned int absCalibProcessorHashFun(PacPidCalibSample * const & samplePtr){
  return PacPidCalibSample::hashFunction(*samplePtr);
}

unsigned int absCalibProcessorHashFunConst(const PacPidCalibSample * const & samplePtr){
  return PacPidCalibSample::hashFunction(*samplePtr);
}

//----------------
// Constructors --
//----------------

PacPidAbsCalibProcessor::PacPidAbsCalibProcessor( const char* const theName, 
					    const char* const theDescription )
  : AppFilterModule( theName, theDescription ),
    _privateData(new PacPidAbsCalibProcessorData())
{
  _privateData->tagFilterMode = new AbsParmBool("tagFilterMode", this, false);
  commands()->append(_privateData->tagFilterMode);
  
  _privateData->tagBitFilterMode = new AbsParmBool("tagBitFilterMode", this, false);
  commands()->append(_privateData->tagBitFilterMode);

  _privateData->microFilterMode = new AbsParmBool("microFilterMode", this, false);
  commands()->append(_privateData->microFilterMode);


  _privateData->sampleNames = new AbsParmVector<string>("sampleNames", this);
  commands()->append(_privateData->sampleNames);

  _privateData->sampleDict = 
    new std::map<string*, PacPidCalibSample*, babar::Collection::PtrLess>();

  _privateData->sampleStats = 
    new std::map<string, int>();

  _privateData->executedSamples = new vector<PacPidCalibSample*>;
  _privateData->knownSamples = new  vector<PacPidCalibSample*>;

  _privateData->moduleNames = new AbsParmVector<string>("moduleNames", this);
  commands()->append(_privateData->moduleNames);
  
  _privateData->moduleDict = 
    new std::map<string*, PacPidAbsPidCalibModule*, babar::Collection::PtrLess>();

  _privateData->executedModules = new vector<PacPidAbsPidCalibModule*> ;
  _privateData->knownModules = new vector<PacPidAbsPidCalibModule*> ;

  _privateData->allKnownListNames = new std::set<string*,BbrPtrLess>();


  _privateData->eventListNames = new vector<string*>;

  _privateData->numEventsSeen = 0;
  _privateData->numEventsPassed = 0;

  //  _privateData->sampleIndex = 
  //    new RWTValHashDictionary<PacPidCalibSample*, int>(&absCalibProcessorHashFun,25);

  _privateData->makeTagBits = new AbsParmBool("makeTagBits", this, false);

  commands()->append(_privateData->makeTagBits);

  
  _privateData->putListsInEvent = new AbsParmBool("putListsInEvent", this, true);

  commands()->append(_privateData->putListsInEvent);

  _privateData->listNameSuffix = new AbsParmString("listNameSuffix", this ,"");
  
  commands()->append(_privateData->listNameSuffix);

  //  _privateData->sampleMap = 0;

}


//--------------
// Destructor --
//--------------

PacPidAbsCalibProcessor::~PacPidAbsCalibProcessor( )
{
  //  _privateData->allKnownListNames->clearAndDestroy();
  for_each(_privateData->allKnownListNames->begin(), 
	   _privateData->allKnownListNames->end(),
	   DeleteObject());
  //  for (std::set<string*, BbrPtrLess>::iterator aln= _privateData->allKnownListNames->begin();
  //       aln != _privateData->allKnownListNames->end(); aln++) {
  //    delete *aln;
  //  }
  _privateData->allKnownListNames->clear();
  //
  delete _privateData->executedSamples;
  int i;
  for (i=0;i<_privateData->knownSamples->size();i++){
    delete (*_privateData->knownSamples)[i];
  }

  delete _privateData->knownSamples;
  delete _privateData->executedModules;
  for (i=0;i<_privateData->knownModules->size();i++){
    delete (*_privateData->knownModules)[i];
  }

  for ( std::map<string*, PacPidCalibSample*, babar::Collection::PtrLess>::iterator 
        it = _privateData->sampleDict->begin(); it != _privateData->sampleDict->end(); 
        it++) {
           delete it->first;
  }

  for ( std::map<string*, PacPidAbsPidCalibModule*, babar::Collection::PtrLess>::iterator
        it = _privateData->moduleDict->begin(); it != _privateData->moduleDict->end();
        it++) {
           delete it->first;
  }

  for ( i = 0; i< _privateData->eventListNames->size(); i++){
    delete (*_privateData->eventListNames)[i];
  }
   
  //  if (_privateData->sampleMap != 0) delete _privateData->sampleMap;

  delete _privateData->knownModules;
  delete _privateData->tagFilterMode;
  delete _privateData->tagBitFilterMode;
  delete _privateData->sampleNames;
  delete _privateData->sampleStats;
  delete _privateData->sampleDict;
  delete _privateData->moduleNames;
  delete _privateData->moduleDict;
  delete _privateData->eventListNames;
  delete _privateData->microFilterMode;
  //  delete _privateData->sampleIndex;
  delete _privateData->makeTagBits;
  delete _privateData->putListsInEvent;
  delete _privateData->listNameSuffix;
  delete _privateData->allKnownListNames;
  delete _privateData;
}

//--------------
// Operations --
//--------------

AppResult
PacPidAbsCalibProcessor::beginJob( AbsEvent* anEvent )
{
  
  if (_verbose.value() )  {
    ErrMsg(routine)  << name() << " begin Job"<< endmsg;
  }
  
  vector<string> samples(_privateData->sampleNames->value());
  size_t i;

  _privateData->executedSamples->reserve(samples.size());
  for (i=0;i<samples.size();i++){
    PacPidCalibSample *theCalibSample(0);

    string tmpStr=samples[i].c_str();
    
    std::map<string*, PacPidCalibSample*, 
      babar::Collection::PtrLess>::iterator iter;
    iter = _privateData->sampleDict->find(&tmpStr);
    if (iter != _privateData->sampleDict->end()) { // Is this check needed ??? (VB - 05/13/02)
    //      theCalibSample = _privateData->sampleDict->findValue(&tmpStr);
      theCalibSample = iter->second;
    }    
    if (theCalibSample==0){
      ErrMsg(fatal) << "Could not find PacPidCalibSample of name " 
		    << samples[i] << " in PacPidCalibSampleFactory." << endmsg;
    }
    
    (* _privateData->sampleStats)[theCalibSample->sampleName()] = 0;
    
  //    _privateData->executedSamples->append(theCalibSample);
    _privateData->executedSamples->push_back(theCalibSample);
    theCalibSample->setUp();
    
  }
  
  vector<string> modules(_privateData->moduleNames->value());
  _privateData->executedModules->reserve(modules.size());
  for (i=0;i<modules.size();i++){
    PacPidAbsPidCalibModule *theCalibModule(0);
    
    string tmpStr=modules[i].c_str();
    
    std::map<string*, PacPidAbsPidCalibModule*,
      babar::Collection::PtrLess>::iterator iter;
    iter = _privateData->moduleDict->find(&tmpStr);
    if (iter != _privateData->moduleDict->end()) { // Is this check needed ??? (VB - 05/13/02)
      theCalibModule = iter->second;
    }    
    
    if (theCalibModule==0){
      ErrMsg(fatal) << "Could not find PacPidCalibModule of name " 
		    << modules[i] << " in PacPidCalibModuleFactory." << endmsg;
    }
    _privateData->executedModules->push_back(theCalibModule);
    theCalibModule->setUp();
    
  }
  
  
  if (_privateData->tagFilterMode->value() &&
      _privateData->microFilterMode->value() ){
  // Make the error message big so that people will notice it.  
  // The beginJob output is very large, so this has to stand out
    ErrMsg(error) << "*************************************************************"
		  << endl
		  << "Warning.  You have configured your BtaCalibProcessor subclass"
		  << endl
		  << "to be in both tag and micro filter mode.  Resetting to micro"
		  << "filter mode. " << endl
		  << "*************************************************************"
		  << endmsg;
    _privateData->tagFilterMode->set(false);
  }
  
  const vector<PacPidCalibSample*> &knownSamples = 
    this->knownSamples();
  
  
// Finally, check that if we're making or decoding tag bits, 
// that the processor name is the one which knows about all the 
// samples.  If not, the tag word will make no sense. 
  if( _privateData->makeTagBits->value()) {
    
  //  declare one bit per sample with the sample name as key
    int ii;
    for (ii=0; ii<knownSamples.size(); ii++)
      {
	const char  * name = knownSamples[ii]->sampleName().c_str();
	declareBool(name);
      }
  // this was the old way: 2 packed integers for everything
  //    declareInt("BtaPidCalibLeft16");
  //    declareInt("BtaPidCalibRight16");
    
    
  //      _privateData->sampleMap = 
  //        new ComTemplatedBitMap<PacPidCalibSample>(knownSamples,
  //  						&absCalibProcessorHashFunConst);
    
    if (this->processorName() != "BtaAllCalibSamplesProcessor"){
      ErrMsg(fatal) << "Error.  Cannot use processor name " 
		    << this->processorName() << " to make or decode tag bits for "
		    << "PacPidCalib samples.  Can only use BtaAllCalibSamplesProcessor"
		    << endl
		    << "If you have set makeTagBits by mistake, please unset it" 
		    << endmsg;
      ::abort();
    }
  }
  
  return AppResult::OK;
}

AppResult
PacPidAbsCalibProcessor::event( AbsEvent* anEvent )
{
  
  if (_verbose.value() )  {
    cout  << name() << " event"<< endl;
  }

  _privateData->numEventsSeen++;

  if ( ! TagAccessor::setup(anEvent) ) {
    ErrMsg(error) << "TagAccessor::setup() failed" << endmsg;
  }

  const AbsEventTag *theTag = tag();

  bool passed;
  if (_privateData->tagFilterMode->value()){
    passed = false;
    
    if (theTag == 0) {
      ErrMsg(error) << "Error in " << name() << ". Cannot locate Tag in event.  Bypassing "
		    << "tag selection for this event." << endmsg;
      passed = true;
    } else {

      // We take the or of any of the samples. 
      size_t i;
      for ( i = 0; i< _privateData->executedSamples->size(); i++){
	PacPidCalibSample *theCalibSample = (*_privateData->executedSamples)[i];
	if (_privateData->tagBitFilterMode->value()) {
	  bool aBool;
	  if (theTag->getBool(aBool,theCalibSample->sampleName())) {
	    if (aBool) { passed=true; };
	  }
	} else {
	  if (theCalibSample->passesTagSelection(theTag, anEvent)){
	    passed = true;
	  }
	}
	if (passed) break;
      }  
    }
  } else {

    // Unless the micro filtering parameter is set, we should
    // always pass the event.  If it is set, then we set the 
    // passed variable to be false only if no selected candidates 
    // at all are found. 
    if (_privateData->microFilterMode->value()){
      passed = false;
    } else {
      passed = true;
    }

    if( _privateData->makeTagBits->value()){
      _privateData->compositeTagWord = 0;
      //      _privateData->sampleMap->resetAllBits();
    }

    size_t i;
    for ( i = 0; i< _privateData->eventListNames->size(); i++){
      delete (*_privateData->eventListNames)[i];
    }
    _privateData->eventListNames->clear();

    int totalNumberOfSelectedGroobers(0);

    for ( i = 0; i< _privateData->executedSamples->size(); i++){

      PacPidCalibSample *theCalibSample = (*_privateData->executedSamples)[i];


      if (_verbose.value() )  {
	cout  << name() << " Checking sample " << theCalibSample->sampleName()
	      << endl;
      }

      // If we are making tag bits, we have not employed the two stage
      // filter.  But we do not want to process every event!  Check 
      // to see if the tag selection is passed first.
      if( _privateData->makeTagBits->value()){
	if (!theCalibSample->passesTagSelection(theTag, anEvent)){
	  continue;
	}
      }	

      theCalibSample->createCandidateLists(anEvent);

      bool checkedThisSample = false;
      int totalNumberOfSelectedGroobersForThisSample(0);
    
      for (PacPidCalibSample::MapIterator theIter = theCalibSample->beginIterator();
	   theIter != theCalibSample->endIterator(); 
	   theIter++) {
      
	string listname(theCalibSample->sampleName());
	listname+=*(theIter->first);
	
	string suffix(_privateData->listNameSuffix->value());
	listname+= suffix;

	if (_privateData->allKnownListNames->find(&listname)  == _privateData->allKnownListNames->end() )
	  _privateData->allKnownListNames->insert(new string(listname));

	AstNamedMapVector<BtaCandidate, BtaCandidate> * theMap 
	  = theIter->second;

	if (theCalibSample->sampleUsedForFilter(*theIter->first)) {
	  totalNumberOfSelectedGroobers+=theMap->members();
	  totalNumberOfSelectedGroobersForThisSample+=theMap->members();
	}


	if (_privateData->microFilterMode->value()){

	  if (totalNumberOfSelectedGroobers > 0 && !checkedThisSample) {
	    // We do not break out of the loop, as the lists need to 
	    // be added to the event for QA/DB purposes
	    checkedThisSample = true;
	    std::map<string, int>::iterator iter = 
	      _privateData->sampleStats->find(theCalibSample->sampleName());
	    iter->second++;  // old RW code did not check for validity 
	    passed = true;
	  }
	}

	//  VB - 8/22/02: MOVED THIS CODE TO THE END OF THE LOOP:
	//       we need to set the tagbit only once / sample !
	//       and don't need to demand that all maps have something, one is enough!
	//       (that may have caused some trouble)

	//	// If we are in the business of making tag bits, 
	//	// find out if the calib sample found anything, 
	//	// set the necessary bit in the tag word.  
	//	
	//  	if (_privateData->makeTagBits->value()) {
	//  	  if(theMap->members() != 0
	//  	     && totalNumberOfSelectedGroobersForThisSample > 0)  {
	//  	    const char * name = theCalibSample->sampleName().c_str();
	//  	    bool result = setBool(name);
	//  	    if (! result) {
	//  	      ErrMsg(error) << "Failed to set tag bit for sample name " 
	//  			    << theCalibSample->sampleName() << endmsg;
	//  	    }
	//  	  }
	//  	}
	//  Remove this whole comment when sure it is correct!
	//
	//   END OF VB - 8/22/02 COMMENT

	
	if (_privateData->putListsInEvent->value()){
	  _privateData->eventListNames->push_back(new string(listname));

	  IfdStrKey key(listname.c_str());
	  HepAList<BtaCandidate> *sampleList;
	  getTmpAList(anEvent, sampleList, key);
	  if (sampleList->length() != 0) {
	    ErrMsg(error) << "Error in PacPidAbsCalibProcessor.  List of name " 
			  <<  *(theIter->first) << " already has " 
			  << sampleList->length() << " entries" 
			  << " my name : " << name()
			  << "\t sample name : " << theCalibSample->sampleName() << endmsg;
	  }
	  
	  // We happend to know that the cloned list is in list 2. YUK. 
	  // However, it'll do for now.
	  
	  const vector<BtaCandidate*> &newList(theMap->list2());
	  for (size_t j=0;j<newList.size();j++) sampleList->append(newList[j]);
      
	  IfdDataProxy<AstNamedMapVector<BtaCandidate, BtaCandidate> > *theDP = 
	    new IfdDataProxy<AstNamedMapVector<BtaCandidate, BtaCandidate> >(theMap);
      
	  bool putIt=Ifd<AstNamedMapVector<BtaCandidate, BtaCandidate> >::put(anEvent, 
									      theDP,key);
	  if (!putIt) {
	    delete theDP;
	    ErrMsg(fatal) << "Failed to place Map of name " << *(theIter->first) 
			  << " into event. " << endmsg;
	  }
	}

      
	// Process MC....
      }  //     for (PacPidCalibSample::MapIterator theIter = theCalibSample->beginIterator(); ...)

      // If we are in the business of making tag bits, 
      // find out if the calib sample found anything, 
      // set the necessary bit in the tag word.  

      // Moved this out of the upper loop: do it only once per sample
      // if at least one of the lists used for the filter got filled...
	
      if (_privateData->makeTagBits->value()) {
	
	const char * name = theCalibSample->sampleName().c_str();
	
	bool initresult = setBool(name,false);
	if (! initresult) {
	  ErrMsg(error) << "Failed to reset tag bit for sample name " 
			<< theCalibSample->sampleName() << endmsg;
	}
	
	if (totalNumberOfSelectedGroobersForThisSample > 0)  {
	  bool result = setBool(name);
	  if (! result) {
	    ErrMsg(error) << "Failed to set tag bit for sample name " 
			  << theCalibSample->sampleName() << endmsg;
	  }
	  if (_verbose.value() ) {
	    static const IfdStrKey theKey("AbsEventID");
	    AbsEventID* eventID = Ifd<AbsEventID>::get( anEvent ,theKey );
	    //	    HepAList< EventInfo >* infoList= Ifd<HepAList<EventInfo> >
	    //  ::get(anEvent, IfdStrKey("Default"));
	    //      const EventInfo* eventInfo(infoList->first());
	    cout << "Event accepted by sample: " << theCalibSample->sampleName()
		 << *eventID << endl;
	  }
	}
      }

      if (_privateData->microFilterMode->value() && passed){
	// Don't do anything else if it's already passed!
	break;
      }
    }
    
    if (_privateData->makeTagBits->value()){
      // this was the old way of doing it: remove it in December 2000 if new scheme working...
      //      theTag->putInt(_privateData->sampleMap->leftHandRep(), "BtaPidCalibLeft16");
      //      theTag->putInt(_privateData->sampleMap->rightHandRep(), "BtaPidCalibRight16");

    }
      

    if (_privateData->microFilterMode->value()){
      setPassed(passed);
    }
      
    // Now process all the attached modules.

    for ( i = 0; i< _privateData->executedModules->size(); i++){
      PacPidAbsPidCalibModule *theModule = (*_privateData->executedModules)[i];
      theModule->process(anEvent);
    }

    // Cleanup

    if (!_privateData->putListsInEvent->value()){
      // Have to clean up the work of the Samples
      // if they don't end up in the event
      for ( i = 0; i< _privateData->executedSamples->size(); i++){

	PacPidCalibSample *theCalibSample = (*_privateData->executedSamples)[i];
	for (PacPidCalibSample::MapIterator theIter = theCalibSample->beginIterator();
	     theIter != theCalibSample->endIterator(); 
	     theIter++) {
	  
	  AstNamedMapVector<BtaCandidate, BtaCandidate> * theMap 
	    = theIter->second;
	  const vector<BtaCandidate*> &newList(theMap->list2());
	  for (size_t j=0;j<newList.size();j++) delete newList[j];
	  
	  delete theMap;
	}
      }
    }
    
    // I'm not _sure_ this is necessary.  Leave it in for now.

    for (i = 0; i< _privateData->executedSamples->size(); i++){

      PacPidCalibSample *theCalibSample = (*_privateData->executedSamples)[i];

      theCalibSample->clearMap();
    }
  }

  if (passed) _privateData->numEventsPassed++;
  setPassed(passed);

  return AppResult::OK;
}

AppResult
PacPidAbsCalibProcessor::endJob( AbsEvent* anEvent )
{

  if ( _privateData->allKnownListNames->size() > 0 ) {
    ostream &myOS = ErrMsg(error);
    myOS << "Calib processor " << name() << " put the following samples in the event:" 
	 << endl;
    for (std::set<string*, BbrPtrLess>::iterator myHSIter = _privateData->allKnownListNames->begin();
	 myHSIter != _privateData->allKnownListNames->end(); myHSIter++ ) {
      myOS << **myHSIter << endl;
    }
    myOS << endmsg;
  }
  if (_privateData->tagFilterMode->value() ||
      _privateData->microFilterMode->value()) {
    ErrMsg(error) << name() << " processed " << _privateData->numEventsSeen
		  << " events and passed " << _privateData->numEventsPassed << endmsg;
  }

  if (_privateData->microFilterMode->value()){
    std::map<string, int>::iterator iter;
    ErrMsg(error) << "Stats for " << name() << "'s Samples" << endmsg;
    for (iter = _privateData->sampleStats->begin(); 
	 iter != _privateData->sampleStats->end(); iter++) {
      ErrMsg(error) << "Sample name: " << iter->first << " passed " 
		    << iter->second << " events." << endmsg;
    }
  }
	
  return AppResult::OK;
}

AppResult
PacPidAbsCalibProcessor::abortJob( AbsEvent* anEvent )
{
  return AppResult::OK;
}

AppModule*
PacPidAbsCalibProcessor::clone(const char* cloneName){
  ErrMsg(error) << "Error PacPidAbsCalibProcessor::clone function now overriden in module name " 
		<< name() << endmsg;
  return 0;
}

PacPidAbsPidCalibModule *
PacPidAbsCalibProcessor::generateCalibModule(const string &name) {
  // Use global function from PacPidCSFactoryManager
  return newPacPidExampleCalibModule(name.c_str(), this);
}

PacPidAbsPidCalibModule *
PacPidAbsCalibProcessor::generateCalibModule(const string &name,
					  const string &newName){
  return newPacPidExampleCalibModule(name.c_str(), newName.c_str(), this);
}

PacPidCalibSample *
PacPidAbsCalibProcessor::generateCalibSample(const string &name) {
  // Use global function from PacPidCSFactoryManager
  return newPacPidCalibSample(name.c_str(), this);
}

PacPidCalibSample *
PacPidAbsCalibProcessor::generateCalibSample(const string &name, const string &newName)  {
  // Use global function from PacPidCSFactoryManager
  return newPacPidCalibSample(name.c_str(), newName.c_str(), this);
}


bool
PacPidAbsCalibProcessor::registerSample(const string &sampleName){
  string *name = new string(sampleName);
  std::map<string*, PacPidCalibSample*, 
    babar::Collection::PtrLess>::iterator iter = 
    _privateData->sampleDict->find(name);
  if (iter != _privateData->sampleDict->end()) {
     delete name;
     return false;
  }
  PacPidCalibSample *theSample = this->generateCalibSample(sampleName);
  assert(theSample != 0);
  (* _privateData->sampleDict)[name] = theSample;
  _privateData->knownSamples->push_back(theSample);
  return true;
}

bool
PacPidAbsCalibProcessor::registerSample(const string &sampleName,
				     const string &newSampleName){
  string *name = new string(newSampleName);
  std::map<string*, PacPidCalibSample*, 
    babar::Collection::PtrLess>::iterator iter = 
    _privateData->sampleDict->find(name);
  if (iter != _privateData->sampleDict->end()) {
     delete name;
     return false;
  }
  PacPidCalibSample *theSample = this->generateCalibSample(sampleName, newSampleName);
  assert(theSample != 0);
  (* _privateData->sampleDict)[name] = theSample;
  _privateData->knownSamples->push_back(theSample);
  return true;
}

bool
PacPidAbsCalibProcessor::registerModule(const string &moduleName){
  string *name = new string(moduleName);
  std::map<string*, PacPidAbsPidCalibModule*,
    babar::Collection::PtrLess>::iterator iter = 
    _privateData->moduleDict->find(name);
  if (iter != _privateData->moduleDict->end()) return false;

  PacPidAbsPidCalibModule *theModule = this->generateCalibModule(moduleName);
  assert(theModule != 0);
  (* _privateData->moduleDict)[name] = theModule;
  _privateData->knownModules->push_back(theModule);
  return true;
}

bool
PacPidAbsCalibProcessor::registerModule(const string &moduleName,
				     const string &newModuleName){
  string *name = new string(newModuleName);
  std::map<string*, PacPidAbsPidCalibModule*,
    babar::Collection::PtrLess>::iterator iter = 
    _privateData->moduleDict->find(name);
  if (iter != _privateData->moduleDict->end()) return false;
  PacPidAbsPidCalibModule *theModule = this->generateCalibModule(moduleName, newModuleName);
  assert(theModule != 0);
  (* _privateData->moduleDict)[name] = theModule;
  _privateData->knownModules->push_back(theModule);
  return true;
}


const std::vector<PacPidCalibSample*> &
PacPidAbsCalibProcessor::executedSamples() const{
  return *_privateData->executedSamples;
}

const std::vector<PacPidAbsPidCalibModule*> &
PacPidAbsCalibProcessor::executedModules() const{
  return *_privateData->executedModules;
}

const vector<PacPidCalibSample*> &
PacPidAbsCalibProcessor::knownSamples() const{
  return *_privateData->knownSamples;
}

const std::vector<PacPidAbsPidCalibModule*> &
PacPidAbsCalibProcessor::knownModules() const{
  return *_privateData->knownModules;
}

const std::vector<std::string*> &
PacPidAbsCalibProcessor::eventListNames() const{
  return *_privateData->eventListNames;
}


