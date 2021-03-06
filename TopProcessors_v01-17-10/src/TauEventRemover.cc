#include "TauEventRemover.h"
#include <iostream>
#include <math.h>
#include <vector>

#include <EVENT/LCCollection.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>

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


TauEventRemover aTauEventRemover ;


TauEventRemover::TauEventRemover() : Processor("TauEventRemover") 
{

  // modify processor description
  _description = "TauEventRemover looks at top decays and removes event in which tops produce a W which decays into a Tau: choice=1 removes Taus, choice=0 only accepts Taus";


  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::MCPARTICLE,
			   "CollectionName" , 
			   "Name of the MCParticle collection"  ,
			   _colName ,
			   std::string("MCParticlesSkimmed")
			   );

  registerProcessorParameter( "Choice",
			      "Choice of process to allow through",
			      _choice,
			      int(1));

}

void TauEventRemover::init()
{ 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  nskipped=0;
  npassed=0;

  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );
  
  if(  pHistogramFactory!=0 )
    {
      if (!(pTree->cd( "/" + name() + "/"))) 
	{
	  pTree->mkdir( "/" + name() + "/" );
	  pTree->cd( "/" + name() + "/");
	}
			

	_ndecaytype = new TH1F("decaytype", "Decay type label", 5, -0.5, 4.5);
     	
    }
}

void TauEventRemover::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 



void TauEventRemover::processEvent( LCEvent * evt ) 
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
    
  //assume success, look for failure
  int decaytype = 1; 

  if( col != NULL )
    {
      int nMCP = col->getNumberOfElements()  ;
      for(int i=6; i< 12 ; i++)
	{
	  MCParticle* mcp = dynamic_cast<MCParticle*>( col->getElementAt( i ) ) ;

	  if(abs(mcp->getPDG())== 15 ) //check if there is a Tau at the generator level
	    {
	      decaytype=0;
	      break;
	    } 
	} //end of MC Loop
    }
	   
  //fill histograms 
  _ndecaytype->Fill(decaytype);

  //  std::cout<<"processed event "<<_nEvt<<std::endl ;
  _nEvt ++ ;
  
  
  if(decaytype!=_choice)
    {
  	 
      nskipped++;
      throw SkipEventException( this );
    }

  else 
    {
      npassed++;
	  
    }
   
}

void TauEventRemover::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TauEventRemover::end(){ 

  // std::cout << "TauEventRemover::end()  " << name() 
  //   	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  //  	    << std::endl ;
  std::cout<< "Events passed: "<< npassed<< "\n"<<"Events rejected: "<<nskipped<<std::endl;

}




