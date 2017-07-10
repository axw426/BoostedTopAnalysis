#include "TopBoostProcessor.h"
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

///replace scratchprocessor with name of your processor (which should match that given for the class in your header file)


TopBoostProcessor aTopBoostProcessor ;


TopBoostProcessor::TopBoostProcessor() : Processor("TopBoostProcessor") 
{
  // modify processor description
  _description = "Boosts to rest frame of ttbar system to look at changes in AFB caused by ISR" ;


  // register steering parameters: name, description, class-variable, default value

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "MCParticlesSkimmed",
			  "a list of all reconstructed particles we are searching for jets in.",
			  _MCParticlesSkimmed,
			  "MCParticlesSkimmed");
}



void TopBoostProcessor::init() ///define all histograms and clear vectors in here
{ 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  _nskipped=0;
	
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


void TopBoostProcessor::processRunHeader( LCRunHeader* run) //won't need changed
{
  _nRun++ ;
} 


void TopBoostProcessor::processEvent( LCEvent * evt ) 
{ 

  _vMCJets.clear();
  TotalEnergy=0;
  
  LCCollection* col = evt->getCollection(_MCParticlesSkimmed);

  for(int i=0; i<col->getNumberOfElements(); i++)
    {
      MCParticle* mcp = dynamic_cast<MCParticle*>( col->getElementAt(i) );
      if(i==2){MCElectron=mcp;}
      if(i==3){MCPositron=mcp;}

      if(i>=6 && i<12) // These correspond to the 6 fermion final state particles before decaying
	{
	  if(abs(mcp->getPDG())<= 5 && abs(mcp->getPDG())>= 1)
	    {
	      _vMCJets.push_back(mcp);
	    }

	  else if(abs(mcp->getPDG())== 12 || abs(mcp->getPDG())== 14)
	    {
	      mcNeutrino=mcp;
	    }

	  else if((abs(mcp->getPDG())==11 || abs(mcp->getPDG())==13))
	    {
	      mcIsolatedLepton=mcp;
	    }
	  
	}
    }
  

  //group particles into hadronic and leptonic top
  FindMCTops(_vMCJets, mcNeutrino, mcIsolatedLepton);
  _vRawSeparation.push_back(HelperClass::dira(FinalMCTop,FinalMCLeptonicTop));
  if(FinalMCTop->getCharge()>0){_vRawTheta.push_back(HelperClass::GetCosTheta(FinalMCTop));}
  if(FinalMCLeptonicTop->getCharge()>0){_vRawTheta.push_back(HelperClass::GetCosTheta(FinalMCLeptonicTop));}
  
  //boost to top rest frame
  BoostTops(MCElectron, MCPositron, FinalMCTop, FinalMCLeptonicTop);
  _vBoostedSeparation.push_back(HelperClass::dira(BoostedTop,BoostedAntiTop));
  if(BoostedTop->getCharge()>0){_vBoostedTheta.push_back(HelperClass::GetCosTheta(BoostedTop));}
  if(BoostedAntiTop->getCharge()>0){_vBoostedTheta.push_back(HelperClass::GetCosTheta(BoostedAntiTop));}
  
  //get event energy after ISR
  TotalEnergy+=MCElectron->getEnergy();
  TotalEnergy+=MCPositron->getEnergy();
  _vTotalEnergy.push_back(TotalEnergy);
 

  _nEvt ++ ; 
  std::cout<<"Nevts processed= "<<_nEvt<<std::endl;

}



void TopBoostProcessor::check( LCEvent * evt ) 
{ 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TopBoostProcessor::end()
{ 


  HelperClass::AddBranch(_tRecoMasses, "TotalEnergy", _vTotalEnergy);
  HelperClass::AddBranch(_tRecoMasses, "RawSeparation", _vRawSeparation);
  HelperClass::AddBranch(_tRecoMasses, "BoostedSeparation", _vBoostedSeparation);
  HelperClass::AddBranch(_tRecoMasses, "RawTheta", _vRawTheta);
  HelperClass::AddBranch(_tRecoMasses, "BoostedTheta", _vBoostedTheta);
}


void TopBoostProcessor::FindMCTops(std::vector<MCParticle*> _vMCJets, MCParticle* mcNeutrino, MCParticle* mcIsolatedLepton)
{
  ReconstructedParticle* MCTop=NULL;
  ReconstructedParticle* MCLeptonicTop=NULL;
  ReconstructedParticle* MCLeptonicTopTest=NULL;
  ReconstructedParticle* TopTest=NULL;
  float bestTopDiff=10000;
      
  TopTest=HelperClass::CombineThreeParticles(_vMCJets[0],_vMCJets[1],_vMCJets[2]);
  MCLeptonicTopTest=HelperClass::CombineThreeParticles(mcNeutrino,mcIsolatedLepton,_vMCJets[3]);
  if(abs(TopTest->getMass()-173.21)+abs(MCLeptonicTopTest->getMass()-173.21)<bestTopDiff)
    {
      bestTopDiff=abs(TopTest->getMass()-173.21)+abs(MCLeptonicTopTest->getMass()-173.21);
      MCTop=TopTest;
      MCLeptonicTop=MCLeptonicTopTest;
    }

  TopTest=HelperClass::CombineThreeParticles(_vMCJets[0],_vMCJets[1],_vMCJets[3]);
  MCLeptonicTopTest=HelperClass::CombineThreeParticles(mcNeutrino,mcIsolatedLepton,_vMCJets[2]);
  if(abs(TopTest->getMass()-173.21)+abs(MCLeptonicTopTest->getMass()-173.21)<bestTopDiff)
    {
      bestTopDiff=abs(TopTest->getMass()-173.21)+abs(MCLeptonicTopTest->getMass()-173.21);
      MCTop=TopTest;
      MCLeptonicTop=MCLeptonicTopTest;

    }
      
  TopTest=HelperClass::CombineThreeParticles(_vMCJets[0],_vMCJets[3],_vMCJets[2]);
  MCLeptonicTopTest=HelperClass::CombineThreeParticles(mcNeutrino,mcIsolatedLepton,_vMCJets[1]);
  if(abs(TopTest->getMass()-173.21)+abs(MCLeptonicTopTest->getMass()-173.21)<bestTopDiff)
    {
      bestTopDiff=abs(TopTest->getMass()-173.21)+abs(MCLeptonicTopTest->getMass()-173.21);
      MCTop=TopTest;
      MCLeptonicTop=MCLeptonicTopTest;
	  
    }
      
  TopTest=HelperClass::CombineThreeParticles(_vMCJets[3],_vMCJets[1],_vMCJets[2]);
  MCLeptonicTopTest=HelperClass::CombineThreeParticles(mcNeutrino,mcIsolatedLepton,_vMCJets[0]);
  if(abs(TopTest->getMass()-173.21)+abs(MCLeptonicTopTest->getMass()-173.21)<bestTopDiff)
    {
      bestTopDiff=abs(TopTest->getMass()-173.21)+abs(MCLeptonicTopTest->getMass()-173.21);
      MCTop=TopTest;
      MCLeptonicTop=MCLeptonicTopTest;

    }

  FinalMCTop=MCTop;
  FinalMCLeptonicTop=MCLeptonicTop;

}


void TopBoostProcessor::BoostTops(MCParticle* electron, MCParticle* positron, ReconstructedParticle* top, ReconstructedParticle* antitop)
{

  ReconstructedParticle* CoM=HelperClass::CombineParticles(electron, positron);

  BoostedTop=Boost(top,CoM);
  BoostedAntiTop=Boost(antitop,CoM);
  
}

ReconstructedParticle* TopBoostProcessor::Boost(ReconstructedParticle* Top, ReconstructedParticle* CentreOfMass)
{
  //aw- use TLorentz vector for doing boost
  
  TLorentzVector Top4Vec(Top->getMomentum()[0],Top->getMomentum()[1],Top->getMomentum()[2],Top->getEnergy());
  TLorentzVector CentreOfMass4Vec(-CentreOfMass->getMomentum()[0],-CentreOfMass->getMomentum()[1],-CentreOfMass->getMomentum()[2],CentreOfMass->getEnergy());
  TVector3 b=CentreOfMass4Vec.BoostVector();
  Top4Vec.Boost(b);
		
  ReconstructedParticleImpl* newParticle = new ReconstructedParticleImpl();
	
  double Energy = Top4Vec.E();
  double Charge=Top->getCharge();
  
  double *Momentum = new double[3];
  Momentum[0] = Top4Vec.Px();
  Momentum[1] = Top4Vec.Py();
  Momentum[2] = Top4Vec.Pz();	
	
  newParticle->setEnergy( Energy );
  newParticle->setMomentum( Momentum );
  newParticle->setCharge( Charge );

  //now return the jet
  delete Momentum;
  return newParticle;
}

