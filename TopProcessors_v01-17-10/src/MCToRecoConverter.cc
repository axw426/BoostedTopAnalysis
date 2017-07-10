#include "MCToRecoConverter.h"
#include <iostream>
#include <math.h>
#include <vector>

#include <EVENT/LCCollection.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCRelation.h>
#include <UTIL/LCRelationNavigator.h>
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCRelationImpl.h"

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


MCToRecoConverter aMCToRecoConverter ;


MCToRecoConverter::MCToRecoConverter() : Processor("MCToRecoConverter") 
{

  // modify processor description
  _description = "MCToRecoConverter looks at top decays and removes event in which tops produce a W which decays into a Tau: choice=1 removes Taus, choice=0 removes everything else";


  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::MCPARTICLE,
			   "CollectionName" , 
			   "Name of the MCParticle collection"  ,
			   _colName ,
			   std::string("MCParticlesSkimmed")
			   );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "OutputCollectionName",
			    "OutputCollectionName",
			    _OutColName,
			    std::string("FakeReco") );

  registerOutputCollection( LCIO::LCRELATION,
			    "FakeRecoRelation" , 
			    "FakeRecoRelation"  ,
			    _RelColName ,
			    std::string("FakeRecoRelation") ) ;

}

void MCToRecoConverter::init()
{ 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

}

void MCToRecoConverter::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 



void MCToRecoConverter::processEvent( LCEvent * evt ) 
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

  //create outputcollection
  LCCollectionVec* outCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;
  //outCol->setSubset(true) ;

  //create relation collection
  LCCollectionVec* relCol = new LCCollectionVec( LCIO::LCRELATION ) ;

  
  if( col != NULL )
    {
      int nMCP = col->getNumberOfElements()  ;
      for(int i=0; i< nMCP ; i++)
	{
	  mcp = dynamic_cast<MCParticle*>( col->getElementAt( i ) ) ;
	  if(mcp->getGeneratorStatus()!=1 || mcp->getPDG()==12 || mcp->getPDG()==14 || mcp->getPDG()==16) {continue;}
	  else
	    {
	      rp=ConvertToReconstructedParticle(mcp);
	      outCol->addElement( rp );
	      relCol->addElement( new LCRelationImpl( rp , mcp , 1 )    ) ;
	    }
	} //end of MC Loop
    }

  evt->addCollection( outCol, _OutColName.c_str() );
  evt->addCollection( relCol, _RelColName.c_str() );

  std::cout<<"processed event "<<_nEvt<<std::endl ;
  _nEvt ++ ;


}

void MCToRecoConverter::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void MCToRecoConverter::end(){ 



}



ReconstructedParticle* MCToRecoConverter::ConvertToReconstructedParticle(MCParticle *p1 )
{
  ReconstructedParticleImpl* newParticle = new ReconstructedParticleImpl();
  newParticle->setMass( p1->getMass() );
  newParticle->setEnergy( p1->getEnergy() );
  newParticle->setMomentum( p1->getMomentum() );
  newParticle->setType( p1->getPDG() );
  newParticle->setCharge( p1->getCharge() );
 
  return newParticle;
}
