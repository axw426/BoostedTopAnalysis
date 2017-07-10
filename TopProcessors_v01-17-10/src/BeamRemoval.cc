#include "BeamRemoval.h"
#include <iostream>
#include <math.h>
#include <vector>

#include <EVENT/LCCollection.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/LCCollectionVec.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include <marlin/Exceptions.h>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA


using namespace lcio ;
using namespace marlin ;


BeamRemoval aBeamRemoval ;


BeamRemoval::BeamRemoval() : Processor("BeamRemoval") 
{

  // modify processor description
  _description = "BeamRemoval looks at top decays and removes event in which tops produce a W which decays into a Tau: choice=1 removes Taus, choice=0 removes everything else";


  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "CollectionName" , 
			   "Name of the Jet collection"  ,
			   _colName ,
			   std::string("kt_6Jets")
			   );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "OutputJets",
			    "Jets without ",
			    _outputCol,
			    std::string("FakeJets") );
}

void BeamRemoval::init()
{ 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

}

void BeamRemoval::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 



void BeamRemoval::processEvent( LCEvent * evt ) 
{ 

  LCCollection* col = NULL;
  
  try
    {
      col = evt->getCollection( _colName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _colName << " collection not available" << std::endl;
      col = NULL;
    }

  std::vector<ReconstructedParticle*> jets;
  jets.resize(col->getNumberOfElements());
  if( col != NULL )
    {
      int nMCP = col->getNumberOfElements()  ;
      for(int i=0; i< nMCP ; i++)
	{
	  ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( col->getElementAt( i ) ) ;
	  jets[i]=jet;
	}
    }

  
  LCCollectionVec* CrudeJets = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;
  CrudeJets->setSubset(true) ;
  for(unsigned int i=0; i<jets.size(); i++)
    {
      for(unsigned int j=0; j<jets[i]->getParticles().size(); j++)
	{
	  CrudeJets->addElement(jets[i]->getParticles()[j]);
	}
    }
  evt->addCollection( CrudeJets, _outputCol.c_str());


 
  
  _nEvt ++ ;

  std::cout<<"processed event "<<_nEvt<<std::endl ;
   
}

void BeamRemoval::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void BeamRemoval::end(){ 

  // std::cout << "BeamRemoval::end()  " << name() 
  //   	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  //  	    << std::endl ;

}




