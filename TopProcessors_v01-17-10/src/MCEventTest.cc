#include "MCEventTest.h"
#include "HelperClass.h"
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


MCEventTest aMCEventTest ;


MCEventTest::MCEventTest() : Processor("MCEventTest") 
{

  // modify processor description
  _description = "MCEventTest filter events to find those with two bjets, two qjets and a sensible top mass from the combination of 1 bjet and the 2qjets";


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

    registerProcessorParameter( "CutonEnergy",
			      "Choice of process to allow through",
			      _CutOnEnergy,
			      int(1));
}

void MCEventTest::init()
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
			

	_ndecaytype = new TH1F("decaytype_MC", "Decay type label MC", 5, -0.5, 4.5);
     	
    }
}

void MCEventTest::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 



void MCEventTest::processEvent( LCEvent * evt ) 
{ 
  decaytype=1;
  nBJets=0;
  nQJets=0;
  _vBJets.clear();
  _vQJets.clear();
  
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

  MCParticle *IsolatedLepton, *Neutrino;
  int nIsoLep=0;
  int nNeutrino=0;
  if( col != NULL )
    {
      int nMCP = col->getNumberOfElements()  ;
      for(int i=0; i< nMCP ; i++)
	{
	  MCParticle* mcp = dynamic_cast<MCParticle*>( col->getElementAt( i ) ) ;

	  //count number of b/q jets produced by intial ee interaction
	  if(abs(mcp->getPDG())== 5 && abs(mcp->getParents()[0]->getPDG())==11)
	    {
	      nBJets++;
	    }

	  if(abs(mcp->getPDG())<= 4 && abs(mcp->getPDG())>= 1 && abs(mcp->getParents()[0]->getPDG())==11)
	    {
	      nQJets++;
	    }

	  if(IsCorrectLepton(mcp)){nIsoLep++;}
	  //if(abs(mcp->getPDG())== 11 && mcp->getParents().size()>0 && abs(mcp->getParents()[0]->getPDG())==11){nIsoLep++;}
	  //if(IsCorrectNeutrino(mcp)){nNeutrino++;}
	  if(abs(mcp->getPDG())== 12 && abs(mcp->getParents()[0]->getPDG())==11){nNeutrino++;}
  
	}
    }

  std::cout<<"Neutrinos= "<<nNeutrino<<" Leptons= "<<nIsoLep<<" bjets= " <<nBJets<<" qjets= "<<nQJets<<std::endl;
  std::cout<<"processed event "<<_nEvt<<std::endl ;
  _nEvt ++ ;

    
}

void MCEventTest::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void MCEventTest::end(){ 

  // std::cout << "MCEventTest::end()  " << name() 
  //   	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  //  	    << std::endl ;

}


ReconstructedParticle* MCEventTest::CombineThreeParticle(MCParticle* a, MCParticle* b, MCParticle* c)
{

  ReconstructedParticle* dummy;
  dummy=HelperClass::CombineParticles(a,b);
  dummy=HelperClass::CombineParticles(dummy,c);

  return dummy;
}

bool MCEventTest::IsCorrectLepton(MCParticle *mcp)
{
  bool passedtest =false;

  if((abs(mcp->getPDG())==11 || abs(mcp->getPDG())==13) && mcp->getGeneratorStatus()==1) //if is really a lepton
    {
      //here on it becomes complicated to define the true lepton- MC not very consistant. Here have used prefix g to denote generator level particles (status 2 or 102)
      for (MCParticle* TestMCP = mcp; TestMCP->getParents().size()>0; TestMCP=TestMCP->getParents()[0]) 
	{
	  if(abs(TestMCP->getPDG())==5 || abs(TestMCP->getPDG())==22){break;} // Isolated lepton should never have come from a photon or b jet

	  //Case 1:  Initial Electron + Positron -> gPositron+gElectron -> final lepton
	  if((abs(TestMCP->getParents()[0]->getPDG())==11
	      || abs(TestMCP->getParents()[0]->getPDG())==13) // parent is a lepton
	     && TestMCP->getParents()[0]->getGeneratorStatus()==102 // parent is generator level particle
	     && TestMCP->getParents()[0]->getParents().size()==2 //generator particle produced directly from the e+e- collision
	     && abs(TestMCP->getParents()[0]->getParents()[0]->getPDG())==11
	     && abs(TestMCP->getParents()[0]->getParents()[1]->getPDG())==11
	     && TestMCP->getParents()[0]->getParents()[0]->getParents().size()==2
	     && (abs(TestMCP->getPDG())==11
		 || abs(TestMCP->getPDG())==13))
	    {
	      passedtest=true;
	      break;
	    }

	  //Case 2: gW -> final electron
	  else if(abs(TestMCP->getParents()[0]->getPDG())==24     //came from W decay
		  && (TestMCP->getParents()[0]->getGeneratorStatus()==102
		      || TestMCP->getParents()[0]->getGeneratorStatus()==2) // W present at generator level
		  && (abs(TestMCP->getPDG())==11
		      || abs(TestMCP->getPDG())==13)) // W decayed leptonically
	    {
	      passedtest=true;
	      break;
	    }

	  //case 3: gElectron->g94->final electron
	  else if(abs(TestMCP->getParents()[0]->getPDG())==94     //came from W like cluster
		  && (abs(TestMCP->getParents()[0]->getParents()[0]->getPDG())==11
		      || abs(TestMCP->getParents()[0]->getParents()[0]->getPDG())==13)// W appeared from nowhere....
		  && (TestMCP->getParents()[0]->getGeneratorStatus()==102
		      || TestMCP->getParents()[0]->getGeneratorStatus()==2) // ...but was produced by the generator
		  && (abs(TestMCP->getPDG())==11
		      || abs(TestMCP->getPDG())==13)) // W decayed leptonically
	    {
	      passedtest=true;
	      break;
	    }

		  		  
	}
	      
    }

    
  return passedtest;
}


bool MCEventTest::IsCorrectNeutrino(MCParticle *mcp)
{
  bool passedtest =false;


  if((abs(mcp->getPDG())==12 || abs(mcp->getPDG())==14) && mcp->getGeneratorStatus()==1) //if is really a lepton
    {
      //here on it becomes complicated to define the true lepton- MC not very consistant. Here have used prefix g to denote generator level particles (status 2 or 102)
      for (MCParticle* TestMCP = mcp; TestMCP->getParents().size()>0; TestMCP=TestMCP->getParents()[0]) 
	{
	  if(abs(TestMCP->getPDG())==5 || abs(TestMCP->getPDG())==22){break;} // Isolated lepton should never have come from a photon or b jet

	  //Case 1:  Initial Electron + Positron -> gPositron+gElectron -> final lepton
	  if((abs(TestMCP->getParents()[0]->getPDG())==12
	      || abs(TestMCP->getParents()[0]->getPDG())==14) // parent is a lepton
	     && TestMCP->getParents()[0]->getGeneratorStatus()==102 // parent is generator level particle
	     && TestMCP->getParents()[0]->getParents().size()==2 //generator particle produced directly from the e+e- collision
	     && abs(TestMCP->getParents()[0]->getParents()[0]->getPDG())==11
	     && abs(TestMCP->getParents()[0]->getParents()[1]->getPDG())==11
	     && TestMCP->getParents()[0]->getParents()[0]->getParents().size()==2
	     && (abs(TestMCP->getPDG())==12
		 || abs(TestMCP->getPDG())==14))
	    {
	      passedtest=true;
	      break;
	    }

	  //Case 2: gW -> final electron
	  else if(abs(TestMCP->getParents()[0]->getPDG())==24     //came from W decay
		  && (TestMCP->getParents()[0]->getGeneratorStatus()==102
		      || TestMCP->getParents()[0]->getGeneratorStatus()==2) // W present at generator level
		  && (abs(TestMCP->getPDG())==12
		      || abs(TestMCP->getPDG())==14)) // W decayed leptonically
	    {
	      passedtest=true;
	      break;
	    }

	  //case 3: gElectron->g94->final electron
	  else if(abs(TestMCP->getParents()[0]->getPDG())==94     //came from W like cluster
		  && (abs(TestMCP->getParents()[0]->getParents()[0]->getPDG())==12
		      || abs(TestMCP->getParents()[0]->getParents()[0]->getPDG())==14)// W appeared from nowhere....
		  && (TestMCP->getParents()[0]->getGeneratorStatus()==102
		      || TestMCP->getParents()[0]->getGeneratorStatus()==2) // ...but was produced by the generator
		  && (abs(TestMCP->getPDG())==12
		      || abs(TestMCP->getPDG())==14)) // W decayed leptonically
	    {
	      passedtest=true;
	      break;
	    }

		  		  
	}
	      
    }

    
  return passedtest;
}
