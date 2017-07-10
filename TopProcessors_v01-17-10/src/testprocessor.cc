#include "testprocessor.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <sstream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/VertexImpl.h>

#include <marlin/Exceptions.h>


//Add in type of lcio data you want to use here AW

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA


using namespace lcio ;
using namespace marlin ;
using namespace EVENT;

///replace scratchprocessor with name of your processor (which should match that given for the class in your header file)


testprocessor atestprocessor ;


testprocessor::testprocessor() : Processor("testprocessor") 
{
  // modify processor description
  _description = "Checks there are enough particles to complete lepton isolation and jet finding" ;


  // register steering parameters: name, description, class-variable, default value

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "recParticleIn",
			  "a list of all reconstructed particles we are searching for jets in.",
			  _lcParticleInName,
			  "PandoraPFANewPFOs");
}



void testprocessor::init() ///define all histograms and clear vectors in here
{ 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  _nskipped=0;

}


void testprocessor::processRunHeader( LCRunHeader* run) //won't need changed
{
  _nRun++ ;
} 


void testprocessor::processEvent( LCEvent * evt ) 
{ 
  // this gets called for every event 
  // usually the working horse ...

LCCollection* particlelist = evt->getCollection(_lcParticleInName);
//std::cout<<"Event contains "<<particlelist->getNumberOfElements()<<" events \n";
if(particlelist->getNumberOfElements() < 5)
    {
      _nskipped++;
      streamlog_out( MESSAGE ) << "skipping event as there are too few particles to contain 4 jets"<< std::endl;
      throw SkipEventException( this ) ;
    }

  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		       << "   in run:  " << evt->getRunNumber() << std::endl ;

 

  _nEvt ++ ; ///last line of actual processing
  //  std::cout<<"Nevts processed= "<<_nEvt<<std::endl;

}



void testprocessor::check( LCEvent * evt ) 
{ 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void testprocessor::end()
{ 

  //   std::cout << "ProcessorTemplate::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
  std::cout<<"Nskipped = "<<_nskipped<<std::endl;
  std::cout<<"Nevts processed= "<<_nEvt<<std::endl;
  ////add in all data to histograms here as final step
}






////////// any extra functions needed can be added here for convenience 

