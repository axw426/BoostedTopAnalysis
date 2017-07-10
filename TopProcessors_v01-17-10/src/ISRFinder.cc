#include "ISRFinder.h"
#include "HelperClass.h"

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
using namespace EVENT;

///replace scratchprocessor with name of your processor (which should match that given for the class in your header file)


ISRFinder aISRFinder ;


ISRFinder::ISRFinder() : Processor("ISRFinder") 
{
  // modify processor description
  _description = "Checks there are enough particles to complete lepton isolation and jet finding" ;


  // register steering parameters: name, description, class-variable, default value

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "recParticleIn",
			  "a list of all reconstructed particles we are searching for jets in.",
			  _LoosePFOColName,
			  "PandoraPFANewPFOs");

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "recParticleIn",
			  "a list of all reconstructed particles we are searching for jets in.",
			  _TightPFOColName,
			  "TightSelectedPandoraPFANewPFOs");

    registerInputCollection( LCIO::LCRELATION,
			   "RecoMCTruthLinkCollectionName" , 
			   "Name of the LCRelation collection containing the linker information"  ,
			   _LinkerColName ,
			   std::string("RecoMCTruthLink"));

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCCollectionName" , 
			   "Name of the MCParticle collection"  ,
			   _MCColName ,
			   std::string("MCParticlesSkimmed"));
}



void ISRFinder::init() ///define all histograms and clear vectors in here
{ 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

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
      _tMCData = new TTree("MCPhotonData", "Tree containing info on MC photons");
      _tRecoData = new TTree("RecoPhotonData", "Tree containing info on MC photons");
      _tEventData = new TTree("EventData", "Tree containing info on MC photons");
    }
  

}


void ISRFinder::processRunHeader( LCRunHeader* run) //won't need changed
{
  _nRun++ ;
} 


void ISRFinder::processEvent( LCEvent * evt ) 
{ 
  // this gets called for every event 
  // usually the working horse ...

LCCollection* loosePFOCol = evt->getCollection(_LoosePFOColName);
LCCollection* tightPFOCol = evt->getCollection(_TightPFOColName);
LCCollection* MCCol = evt->getCollection(_MCColName);
LCCollection* LinkerCol = evt->getCollection(_LinkerColName);

 int nLoosePhotons, nTightPhotons;

 for(int i=0;i<loosePFOCol->getNumberOfElements();i++)
   {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( loosePFOCol->getElementAt(i) );
      if(pfo->getType()==22)
	{
	nLoosePhotons++;
	_vRecoEnergy.push_back(pfo->getEnergy());
	}
      
   }

  for(int i=0;i<tightPFOCol->getNumberOfElements();i++)
   {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( tightPFOCol->getElementAt(i) );
      if(pfo->getType()==22){nTightPhotons++;}
   }
  
  totalMCEnergy=0;
  for(int i=0;i<MCCol->getNumberOfElements();i++)
   {
      MCParticle* mcp = dynamic_cast<MCParticle*>( MCCol->getElementAt(i) );
      if(mcp->getGeneratorStatus()==1){totalMCEnergy+=mcp->getEnergy();}
      if(i==12){ISREnergy=mcp->getEnergy();}
      if(i==13){ISREnergy+=mcp->getEnergy();}
   }

  _vISREnergy.push_back(ISREnergy);
  _vTotalMCEnergy.push_back(totalMCEnergy);


  
  //MC linking
  LCRelationNavigator* relationNavigatorPFOMC = new LCRelationNavigator( LinkerCol );
  MCParticle* mcp = dynamic_cast<MCParticle*>( MCCol->getElementAt(12) );
  _vPhotonEnergy.push_back(mcp->getEnergy());
  _vCosTheta.push_back(HelperClass::GetCosTheta(mcp));
  
  EVENT::LCObjectVec relobjMC = relationNavigatorPFOMC->getRelatedFromObjects(mcp);
  ReconstructedParticle *rp=NULL;
  if(relobjMC.size()>0)
    {
      rp=dynamic_cast <ReconstructedParticle*>(relobjMC[0]);
      if (rp!=NULL){_vHasRecoMatch.push_back(1);}
      else {_vHasRecoMatch.push_back(0);}
    }
  else {_vHasRecoMatch.push_back(0);}

  mcp = dynamic_cast<MCParticle*>( MCCol->getElementAt(13) );
  _vPhotonEnergy.push_back(mcp->getEnergy());
  _vCosTheta.push_back(HelperClass::GetCosTheta(mcp));

  relobjMC = relationNavigatorPFOMC->getRelatedFromObjects(mcp);
  rp=NULL;
  if(relobjMC.size()>0)
    {
      rp=dynamic_cast <ReconstructedParticle*>(relobjMC[0]);
      if (rp!=NULL){_vHasRecoMatch.push_back(1);}
      else {_vHasRecoMatch.push_back(0);}
    }
  else {_vHasRecoMatch.push_back(0);}
  delete relationNavigatorPFOMC;

   

   
  _nEvt ++ ; ///last line of actual processing
   std::cout<<"Nevts processed= "<<_nEvt<<std::endl;

}



void ISRFinder::check( LCEvent * evt ) 
{ 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ISRFinder::end()
{ 

  HelperClass::AddBranch(_tRecoData, "RecoEnergy", _vRecoEnergy);

  HelperClass::AddBranch(_tEventData, "ISREnergy", _vISREnergy);
  HelperClass::AddBranch(_tEventData, "TotalMCEnergy", _vTotalMCEnergy);

  HelperClass::AddBranch(_tMCData, "PhotonEnergy", _vPhotonEnergy);
  HelperClass::AddBranch(_tMCData, "HasRecoMatch", _vHasRecoMatch);
  HelperClass::AddBranch(_tMCData, "CosTheta", _vCosTheta);
}






////////// any extra functions needed can be added here for convenience 

