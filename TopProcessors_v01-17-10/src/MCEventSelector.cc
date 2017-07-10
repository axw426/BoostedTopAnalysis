#include "MCEventSelector.h"
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


MCEventSelector aMCEventSelector ;


MCEventSelector::MCEventSelector() : Processor("MCEventSelector") 
{

  // modify processor description
  _description = "MCEventSelector filters events to find various top decays, 1=single hadronic top, 2= hadronic + leptonic top, 3=single leptonic top, 4=tau 0=not top";


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
			      int(2));
}

void MCEventSelector::init()
{ 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  nskipped=0;
  npassed=0;

  //define accepted mass range for top
  LowerTopLimit=166.16;
  UpperTopLimit=180.26;


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

void MCEventSelector::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 



void MCEventSelector::processEvent( LCEvent * evt ) 
{ 
  //reset variables to zero for every event
  decaytype=0;
  nJets=0;
  _vJets.clear();
  IsolatedLepton=NULL; 
  Neutrino=NULL;
  nTaus=0;
  nLeptons=0;
  nNeutrinos=0;
  hadronictopjet=false;
  leptonictopjet=false;
  col = NULL;
  
  try
    {
      col = evt->getCollection( _colName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _colName << " collection not available" << std::endl;
      col = NULL;
    }


  if( col != NULL )
    {
      for(int i=6; i< 12 ; i++) //look only at primary 6 fermions
	{
	  MCParticle* mcp = dynamic_cast<MCParticle*>( col->getElementAt( i ) ) ;

	  //sort into jets/leptons/neutrinos 
	  if(abs(mcp->getPDG())<= 5 && abs(mcp->getPDG())>= 1)
	    {
	      nJets++;
	      _vJets.push_back(mcp);
	    }

  	  else if(abs(mcp->getPDG())== 12 || abs(mcp->getPDG())== 14  || abs(mcp->getPDG())== 16)
	    {
	      Neutrino=mcp;
	      nNeutrinos++;
	    }

	  else if(abs(mcp->getPDG())==11 || abs(mcp->getPDG())==13 || abs(mcp->getPDG())== 15 )
	    {
	      IsolatedLepton=mcp;
	      nLeptons++;
	      if(abs(mcp->getPDG())== 15){nTaus++;}
	    }
	  

	}
    }

  //sanity check!!
  if(nLeptons!=1 || nNeutrinos!=1 || nJets!=4){std::cout<<"Incorrect number of initial fermions! "<<nLeptons<<" Lep "<<nNeutrinos<<" Neut "<<nJets<<" Jets "<<std::endl;}
  else
    { 
      float test0, test1, test2, test3;
  
      //assume lepton and neutrino come from W
      WBoson=CombineParticles(Neutrino,IsolatedLepton);
      
      //check for lepttonic top decay
      for(int i=0;i<4;i++)
	{
	  topmass=CombineParticles(WBoson,_vJets[i])->getMass();
	  if(topmass>LowerTopLimit && topmass<UpperTopLimit)
	    {
	      leptonictopjet=true;
	      _vJets.erase(_vJets.begin()+i); //avoid double counting jet
	      break;
	    }
	}

      //check for second top- only one possible combination left for hadronic top
      if(leptonictopjet==true)
	{
	  topmass=CombineThreeParticles(_vJets[0],_vJets[1],_vJets[2])->getMass();
	  if(topmass>LowerTopLimit && topmass<UpperTopLimit)
	    {
	      hadronictopjet=true;
	    }
	}

      //no leptonic top, need to try all jet combinations for hadronic top
      else
	{
	  test0=CombineThreeParticles(_vJets[1],_vJets[2],_vJets[3])->getMass();
	  test1=CombineThreeParticles(_vJets[0],_vJets[2],_vJets[3])->getMass();
	  test2=CombineThreeParticles(_vJets[0],_vJets[1],_vJets[3])->getMass();
	  test3=CombineThreeParticles(_vJets[0],_vJets[1],_vJets[2])->getMass();

	  if(test0>LowerTopLimit && test0<UpperTopLimit){hadronictopjet=true;}
	  else if(test1>LowerTopLimit && test1<UpperTopLimit){hadronictopjet=true;}
	  else if(test2>LowerTopLimit && test2<UpperTopLimit){hadronictopjet=true;}
	  else if(test3>LowerTopLimit && test3<UpperTopLimit){hadronictopjet=true;}
	}
    }
  if(leptonictopjet==true && hadronictopjet==true){decaytype=2;}
  else if(leptonictopjet==true){decaytype=3;}
  else if(hadronictopjet==true){decaytype=1;}
  else {decaytype=0;}

  if(decaytype>0 && nTaus>0){decaytype=4;}
  
  _ndecaytype->Fill(decaytype);


  //  std::cout<<"decaytype= "<<decaytype<<std::endl; //debugging
  if(_nEvt%1000==0){  std::cout<<"processed event "<<_nEvt<<std::endl ;}
  _nEvt ++ ;

  //skip unwanted events
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

void MCEventSelector::check( LCEvent * evt ) { 
}


void MCEventSelector::end(){ 

  std::cout<< "Events passed: "<< npassed<< "\n"<<"Events rejected: "<<nskipped<<std::endl;

}


ReconstructedParticle* MCEventSelector::CombineThreeParticles(MCParticle* a, MCParticle* b, MCParticle* c)
{

  ReconstructedParticle* dummy;
  dummy=CombineParticles(a,b);
  dummy=CombineParticles(dummy,c);

  return dummy;
}

double MCEventSelector::ReconstructedMass( MCParticle* p1, MCParticle* p2 )
{
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  double mass_squared = pow( p1->getEnergy()+p2->getEnergy(), 2 ) 
    - pow( mom1[0]+mom2[0],2 ) 
    - pow( mom1[1]+mom2[1],2 )
    - pow( mom1[2]+mom2[2],2 );
											
  return pow( mass_squared, 0.5 );
}

ReconstructedParticle* MCEventSelector::CombineParticles( MCParticle *p1, MCParticle *p2 )
{
  // create a ReconstructedParticle that saves the jet
  ReconstructedParticleImpl* newParticle = new ReconstructedParticleImpl();
	
  double Energy = p1->getEnergy() + p2->getEnergy();

  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  //sum the momentums together
  double *Momentum = new double[3];
  Momentum[0] = mom1[0]+mom2[0];
  Momentum[1] = mom1[1]+mom2[1];
  Momentum[2] = mom1[2]+mom2[2];
	
  double mass = ReconstructedMass( p1, p2 );
	
  //set the parameters
  newParticle->setMass( mass );
  newParticle->setEnergy( Energy );
  newParticle->setMomentum( Momentum );
	

  //now return the jet
  delete Momentum;
  return newParticle ;
}


double MCEventSelector::ReconstructedMass( ReconstructedParticle* p1, MCParticle* p2 )
{
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  double mass_squared = pow( p1->getEnergy()+p2->getEnergy(), 2 ) 
    - pow( mom1[0]+mom2[0],2 ) 
    - pow( mom1[1]+mom2[1],2 )
    - pow( mom1[2]+mom2[2],2 );
											
  return pow( mass_squared, 0.5 );
}

ReconstructedParticle* MCEventSelector::CombineParticles( ReconstructedParticle *p1, MCParticle *p2 )
{
  // create a ReconstructedParticle that saves the jet
  ReconstructedParticleImpl* newParticle = new ReconstructedParticleImpl();
	
  double Energy = p1->getEnergy() + p2->getEnergy();

  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  //sum the momentums together
  double *Momentum = new double[3];
  Momentum[0] = mom1[0]+mom2[0];
  Momentum[1] = mom1[1]+mom2[1];
  Momentum[2] = mom1[2]+mom2[2];
	
  double mass = ReconstructedMass( p1, p2 );
	
  //set the parameters
  newParticle->setMass( mass );
  newParticle->setEnergy( Energy );
  newParticle->setMomentum( Momentum );
	

  //now return the jet
  delete Momentum;
  return newParticle ;
}
