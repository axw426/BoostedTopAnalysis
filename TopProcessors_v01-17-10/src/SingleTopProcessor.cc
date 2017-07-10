#include "SingleTopProcessor.h"

#include <iostream>
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <cstring>
//#include <initializer_list>
#include <sstream>
#include <fstream>
#include <utility>

#include <EVENT/LCCollection.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <marlin/Exceptions.h>
#include <UTIL/PIDHandler.h>
#include "UTIL/LCRelationNavigator.h"

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA

using namespace lcio ;
using namespace marlin ;


SingleTopProcessor aSingleTopProcessor ;


SingleTopProcessor::SingleTopProcessor() : Processor("SingleTopProcessor") 
{
  // modify processor description
  _description = "DummyProcessor for now, just looks at two fat jets and gets their masses " ;


  registerInputCollection( LCIO::MCPARTICLE,
			   "MCCollectionName" , 
			   "Name of the MCParticle collection"  ,
			   _MCColName ,
			   std::string("MCParticlesSkimmed"));

}


void SingleTopProcessor::init() 
{ 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  //aw Sect1/////////////////////////////////////// Set up Histograms //////////////////////////////////////////////
	
  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );
	
  if( pHistogramFactory!=0 )
    {
      if (!(pTree->cd( "/" + name() + "/"))) 
	{
	  pTree->mkdir( "/" + name() + "/" );
	  pTree->cd( "/" + name() + "/");
	}
		
      //create a tree to contain results
      _tRecoMasses = new TTree("RecoData", "Tree containing final reconstructed masses");
    }

}


void SingleTopProcessor::processRunHeader( LCRunHeader* run) 
{
  _nRun++ ;
} 


void SingleTopProcessor::processEvent( LCEvent * evt ) 
{

  
/////////////////////////////////////// Get Collections ////////////////////////////////////////////
  
  _AllCollectionsExist = true ;

  MCCol = NULL;
  try
    {
      MCCol = evt->getCollection( _MCColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _MCColName << " collection not available" << std::endl;
      MCCol = NULL;
      _AllCollectionsExist = false ;
    }

  
  if( !_AllCollectionsExist )
    {
      streamlog_out(WARNING) << "Very sad...Skipping event " << _nEvt 
			     << " - Cannot find all required collections" << std::endl;
      throw SkipEventException( this );

    }

  // Process Event
 
  nMCP = MCCol->getNumberOfElements();
  Top=NULL;
  for(int i=0; i<nMCP; i++)
    {
      MCParticle* mcp = dynamic_cast<MCParticle*>( MCCol->getElementAt(i) );

      if(mcp->getGeneratorStatus()==1)
	{
	  if( Top==NULL ) Top=ConvertToReconstructedParticle(mcp);
	  Top=CombineParticles(Top,mcp);
	}
    }
  _vTopMass.push_back(Top->getMass());
    
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		       << "   in run:  " << evt->getRunNumber() << std::endl ;
  
  _nEvt ++ ;
      
  std::cout << "Events processed: " << _nEvt << std::endl;
  
}

void SingleTopProcessor::check( LCEvent * evt ) 
{ 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void SingleTopProcessor::end()
{ 
  AddBranch(_tRecoMasses, "TopMass", _vTopMass);


}






void SingleTopProcessor::AddBranch( TTree* input_tree, std::string title, std::vector<float> data )
{
  if( (input_tree->GetEntries() != (int)data.size()) && (input_tree->GetEntries() != 0) )
    {
      std::cout << "ERROR: CANNOT ADD BRANCH " << title << " DATA SIZE DID NOT MATCH TREE." 
      	<< "Data size: " << data.size() << ". Tree size: "<< _vTopMass.size() << std::endl;
      return;
    }

  float thisValue=0.;
  TString name=title.c_str();
  name.Append("/F");
  TBranch* output_branch_obj = input_tree->Branch( title.c_str(), &thisValue, name );

  input_tree->SetEntries( (int)data.size() );

  for( unsigned int i=0; i< data.size(); ++i )
    {
      thisValue = data[i];
      output_branch_obj->Fill();
    }
}

double SingleTopProcessor::dira(ReconstructedParticle *p1, ReconstructedParticle *p2)
{
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();

  double mag1 = pow(pow(mom1[0],2) + pow(mom1[1],2) + pow(mom1[2],2),0.5);
  double mag2 = pow(pow(mom2[0],2) + pow(mom2[1],2) + pow(mom2[2],2),0.5);

  double theta = acos((mom1[0]*mom2[0]+mom1[1]*mom2[1]+mom1[2]*mom2[2])/(mag1*mag2));

  //if(theta<0) //Takes into account negative angles. Just in case
  //  {
  //    return abs(theta);
  //  }
  //else return theta;

  return theta;
}


ReconstructedParticle* SingleTopProcessor::CombineParticles( MCParticle *p1, MCParticle *p2 )
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

double SingleTopProcessor::ReconstructedMass( MCParticle* p1, MCParticle* p2 )
{
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  double mass_squared = pow( p1->getEnergy()+p2->getEnergy(), 2 ) 
    - pow( mom1[0]+mom2[0],2 ) 
    - pow( mom1[1]+mom2[1],2 )
    - pow( mom1[2]+mom2[2],2 );
											
  return pow( mass_squared, 0.5 );
}

ReconstructedParticle* SingleTopProcessor::CombineParticles( ReconstructedParticle *p1, MCParticle *p2 )
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

double SingleTopProcessor::ReconstructedMass( ReconstructedParticle* p1, MCParticle* p2 )
{
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  double mass_squared = pow( p1->getEnergy()+p2->getEnergy(), 2 ) 
    - pow( mom1[0]+mom2[0],2 ) 
    - pow( mom1[1]+mom2[1],2 )
    - pow( mom1[2]+mom2[2],2 );
											
  return pow( mass_squared, 0.5 );
}

ReconstructedParticle* SingleTopProcessor::ConvertToReconstructedParticle(MCParticle *p1 )
{
  ReconstructedParticleImpl* newParticle = new ReconstructedParticleImpl();
  newParticle->setMass( p1->getMass() );
  newParticle->setEnergy( p1->getEnergy() );
  newParticle->setMomentum( p1->getMomentum() );
  return newParticle;
}
