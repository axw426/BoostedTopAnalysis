#include "LeptonOptimizationProcessorTest.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCRelation.h>
#include <EVENT/ReconstructedParticle.h>
#include <marlin/Exceptions.h>
#include <UTIL/PIDHandler.h>
#include "UTIL/LCRelationNavigator.h"

// ----- include for verbosity dependend logging ---------
#include <marlin/VerbosityLevels.h>

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA

using namespace lcio ;
using namespace marlin ;

LeptonOptimizationProcessorTest aLeptonOptimizationProcessorTest ;

LeptonOptimizationProcessorTest::LeptonOptimizationProcessorTest()
  : Processor("LeptonOptimizationProcessorTest") {

  // Processor description
  _description = "Isolated Lepton Finder ProcessorTest-> aw jet isolation performs box vuts on min z and xt" ;

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputCollection" ,
			   "Input collection of Leptons",
			   _inputLepCollection,
			   std::string("IsoLep"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "PFOCollection" ,
			   "Input collection of PFOs",
			   _PFOCollection,
			   std::string("PandoraPFANewPFOs"));

    registerInputCollection( LCIO::MCPARTICLE,
			   "MCCollectionName" , 
			   "Name of the MCParticle collection"  ,
			   _MCColName ,
			   std::string("MCParticlesSkimmed"));

    registerInputCollection( LCIO::LCRELATION,
			     "RecoMCTruthLinkCollectionName" , 
			     "Name of the LCRelation collection containing the linker information"  ,
			     _LinkerColName ,
			     std::string("RecoMCTruthLink"));


   
}


void LeptonOptimizationProcessorTest::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  nEvt=0;
  nReal=0;
  nFake=0;
  nDuplicate=0;
  // usually a good idea to
  printParameters() ;

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
      _tData = new TTree("Duplicates", "Tree containing final reconstructed masses");
      _tCorrectData = new TTree("NoLinks", "Tree containing final reconstructed masses");
    }
}

void LeptonOptimizationProcessorTest::processRunHeader( LCRunHeader* run) { 
} 

void LeptonOptimizationProcessorTest::processEvent( LCEvent * evt ) { 


  try{
    _LepCol = evt->getCollection( _inputLepCollection ) ;
  }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _inputLepCollection << " collection not available, exiting IsolatedLeptonProcessorTest..." << std::endl;
      _LepCol = NULL;
      return;
    }

  try{
    _PFOCol = evt->getCollection( _PFOCollection ) ;
  }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _PFOCollection << " collection not available, exiting IsolatedLeptonProcessorTest..." << std::endl;
      _PFOCol = NULL;
      return;
    }

  try{
    _MCCol = evt->getCollection( _MCColName ) ;
  }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _MCColName << " collection not available, exiting IsolatedLeptonProcessorTest..." << std::endl;
      _MCCol = NULL;
      return;
    }

  try{
    LinkerCol = evt->getCollection( _LinkerColName ) ;
  }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _LinkerColName << " collection not available, exiting IsolatedLeptonProcessorTest..." << std::endl;
      LinkerCol = NULL;
      return;
    }

  int nLeptons=_MCCol->getNumberOfElements();
  
  E_TrueLep=0;
  E_CurrentRP=0;
  bool foundlep=false;
  for (int i=0; i<nLeptons; i++)
    {
      rp= dynamic_cast<MCParticle*>( _MCCol->getElementAt(i) );
       
      if(IsCorrectLepton(rp)==true)
	{
	  if(foundlep==false)
	    {
	      theMCLepton=rp;
	      nReal++;
	    }
	  else{nDuplicate++;}
	      
	  foundlep=true;
	}

      else
	{
	  nFake++;
	}
    }

  relationNavigatorPFOMC = new LCRelationNavigator( LinkerCol );
  relobjMC = relationNavigatorPFOMC->getRelatedFromObjects(theMCLepton);

  if(relobjMC.size()==0)
    {
      //std::cout<<"No link to MC particle!!"<<std::endl;
      nReal--;
      _vNoLinksPDG.push_back(theMCLepton->getPDG());
      Pt=sqrt(theMCLepton->getMomentum()[0]*theMCLepton->getMomentum()[0]+theMCLepton->getMomentum()[1]*theMCLepton->getMomentum()[1]);
      Pz=theMCLepton->getMomentum()[2];
      _vNoLinksPt.push_back(Pt);
      _vNoLinksCosTheta.push_back(Pz/sqrt(Pt*Pt+Pz*Pz));
    }
  else if(relobjMC.size()>1)
    {
      //std::cout<<"MCParticle links to multiple reco!!"<<std::endl;
      nMultiReco+=relobjMC.size()-1;
      
      for (unsigned int i=0;i<relobjMC.size();i++)
	{
	  _vMCDuplicatesPDG.push_back( dynamic_cast<ReconstructedParticle*>(relobjMC[i])->getType() );
	}
    }
  
  delete relationNavigatorPFOMC;

  std::cout<<"Processing Event "<<nEvt<<std::endl;
  if(foundlep==false){std::cout<<"no lep found in event"<<nEvt<<std::endl;}

  nEvt++;

}

void LeptonOptimizationProcessorTest::check( LCEvent * evt ) { 
}

void LeptonOptimizationProcessorTest::end() {

  AddBranch(_tData,"PDG", _vMCDuplicatesPDG);
  AddBranch(_tCorrectData,"PDG", _vNoLinksPDG);
  AddBranch(_tCorrectData,"Pt", _vNoLinksPt);
  AddBranch(_tCorrectData,"CosTheta", _vNoLinksCosTheta);


  std::cout<<"Efficiency= "<<100*(float)nReal/nEvt<<std::endl;
  std::cout<<"Purity= "<<100*(float)nReal/(nReal+nFake)<<std::endl;
  std::cout<<"Duplicates= "<<nDuplicate<<std::endl;

  std::cout<<"Real Leptons Passed= "<<nReal<<std::endl;
  std::cout<<"Fake Leptons Passed= "<<nFake<<std::endl;

  std::cout<<"Number of duplicates expected in RECO= "<<nMultiReco<<std::endl;
}


bool LeptonOptimizationProcessorTest::IsCorrectLepton(ReconstructedParticle *lepton)
{
  passedtest =false;
  LCRelationNavigator* relationNavigatorPFOMC = new LCRelationNavigator( LinkerCol );
  EVENT::LCObjectVec relobjMC = relationNavigatorPFOMC->getRelatedToObjects(lepton);
  for(unsigned int i=0;i<relobjMC.size() ; i++) //always has size one
    {
      MCParticle *mcp=NULL;
      mcp=dynamic_cast <MCParticle*>(relobjMC[i]);
      if (mcp!=NULL) //sometimes no link exists-> tends not to happen for particles within the process definition 
	{
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

		   /*/case 3: gElectron->g94->final electron
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
		    }*/
		   
		  		  
		}
	      
	    }
	}
    }
  delete relationNavigatorPFOMC;
  return passedtest;
}
//used to test definition of IsoLep -> conditions found to allow axactly one particle per event when run on MC
bool LeptonOptimizationProcessorTest::IsCorrectLepton(MCParticle *mcp)
{
  passedtest =false;

  if((abs(mcp->getPDG())==11 || abs(mcp->getPDG())==13) && mcp->getGeneratorStatus()==1) //if is really a lepton
    {
      //here on it becomes dodgy to define the true lepton- MC not very consistant
      for (MCParticle* TestMCP = mcp; TestMCP->getParents().size()>0; TestMCP=TestMCP->getParents()[0]) 
	{
	  if(abs(TestMCP->getPDG())==5 || abs(TestMCP->getPDG())==22){break;}
	  
	  //Case 1: straightforward MC- lepton goes directly from generator to final state particle
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
	      //std::cout<<"LeptonEnergy= "<<mcp->getEnergy()<<std::endl;
	      break;
	      }
	    
	    //Case 2: W present at the generator level instead of Isolep- W not connected to original e+e- in MC
	   if(abs(TestMCP->getParents()[0]->getPDG())==24     //came from W
		  && (TestMCP->getParents()[0]->getGeneratorStatus()==102
		      || TestMCP->getParents()[0]->getGeneratorStatus()==2) // ...but W was produced by the generator
		  && (abs(TestMCP->getPDG())==11
		      || abs(TestMCP->getPDG())==13)) // W decayed leptonically
	    {
	      passedtest=true;
	      //std::cout<<"LeptonEnergy= "<<mcp->getEnergy()<<std::endl;
	      break;
	    }

	   /*/case 3: return of the mc
	     else if(abs(TestMCP->getParents()[0]->getPDG())==94     //came from Wish cluster
	     && (abs(TestMCP->getParents()[0]->getParents()[0]->getPDG())==11
	     || abs(TestMCP->getParents()[0]->getParents()[0]->getPDG())==13)// W appeared from nowhere....
	     && (TestMCP->getParents()[0]->getGeneratorStatus()==102
	     || TestMCP->getParents()[0]->getGeneratorStatus()==2) // ...but was produced by the generator
	     && (abs(TestMCP->getPDG())==11
	     || abs(TestMCP->getPDG())==13)) // W decayed leptonically
	     {
	     passedtest=true;
	     // std::cout<<"LeptonEnergy= "<<mcp->getEnergy()<<std::endl;
	     break;
	     }*/
	   	  
	}
    }	      
  return passedtest;
}

void LeptonOptimizationProcessorTest::AddBranch( TTree* input_tree, std::string title, std::vector<float> data )
{
  if( (input_tree->GetEntries() != (int)data.size()) && (input_tree->GetEntries() != 0) )
    {
      std::cout << "ERROR: CANNOT ADD BRANCH " << title << " DATA SIZE DID NOT MATCH TREE." 
      	<< "Data size: " << data.size() << ". Tree size: "<< _vCorrectLeptonYN.size() << std::endl;
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
