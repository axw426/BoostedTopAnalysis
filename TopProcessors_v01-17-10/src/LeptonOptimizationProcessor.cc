#include "LeptonOptimizationProcessor.h"

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

LeptonOptimizationProcessor aLeptonOptimizationProcessor ;

LeptonOptimizationProcessor::LeptonOptimizationProcessor()
  : Processor("LeptonOptimizationProcessor") {

  // Processor description
  _description = "Isolated Lepton Finder Processor-> aw jet isolation performs box vuts on min z and xt" ;

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

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "JetCollection" ,
			   "Input collection of jets for isolation",
			   _jetCollectionName,
			   std::string("JetsForIsolation"));
   
}


void LeptonOptimizationProcessor::init() { 
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
      _tData = new TTree("RecoData", "Tree containing final reconstructed masses");
      _tCorrectData = new TTree("CorrectRecoData", "Tree containing final reconstructed masses");
    }
}

void LeptonOptimizationProcessor::processRunHeader( LCRunHeader* run) { 
} 

void LeptonOptimizationProcessor::processEvent( LCEvent * evt ) { 


  try{
    _LepCol = evt->getCollection( _inputLepCollection ) ;
  }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _inputLepCollection << " collection not available, exiting IsolatedLeptonProcessor..." << std::endl;
      _LepCol = NULL;
      return;
    }

  try{
    _PFOCol = evt->getCollection( _PFOCollection ) ;
  }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _PFOCollection << " collection not available, exiting IsolatedLeptonProcessor..." << std::endl;
      _PFOCol = NULL;
      return;
    }

  try{
    _MCCol = evt->getCollection( _MCColName ) ;
  }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _MCColName << " collection not available, exiting IsolatedLeptonProcessor..." << std::endl;
      _MCCol = NULL;
      return;
    }

  try{
    LinkerCol = evt->getCollection( _LinkerColName ) ;
  }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _LinkerColName << " collection not available, exiting IsolatedLeptonProcessor..." << std::endl;
      LinkerCol = NULL;
      return;
    }

  nLeptons=_LepCol->getNumberOfElements();
  foundlep=false;
  bestenergy=0;
  bestenergyposition=-1;
  lowestimpactparameter=10000;
  bestimpactparameterposition=-1;
  _cosConeAngle=0.98;

  
  for (int i=0; i<nLeptons; i++)
    {
      rp= dynamic_cast<ReconstructedParticle*>( _LepCol->getElementAt(i) );

      //check if lepton is real or fake
      if(IsCorrectLepton(rp)==true)
	{
	  if(foundlep==false)//check if we had already found a lepton
	    {
	      nReal++;
	    }
	  else{nDuplicate++;}
	      
	  foundlep=true;
	}

      else
	{
	  nFake++;
	}

      //find highest energy candidate
       if(rp->getEnergy()>bestenergy)
	{
	  bestenergyposition=i;
	  bestenergy=rp->getEnergy();
	}

       //jet properties
      _rpJetMap.clear();
      colJet = evt->getCollection(_jetCollectionName);
      njet = colJet->getNumberOfElements();
      for (int j=0; j<njet; ++j) {
	jet = dynamic_cast<ReconstructedParticle*>( colJet->getElementAt(j) );
	for (ReconstructedParticleVec::const_iterator iter = jet->getParticles().begin();
	     iter != jet->getParticles().end(); ++iter) {
	  _rpJetMap.insert( std::make_pair( *iter, jet ) ); 
	}
      }
    
      jetxt=0;
      jetz=0;
      if ( _rpJetMap.find( rp ) != _rpJetMap.end() ) {
	jet = _rpJetMap[rp];
	TVector3 vec1( rp->getMomentum() );
	TVector3 jetmom( jet->getMomentum() );
	TLorentzVector jetmom4( jet->getMomentum(), jet->getEnergy() );
	jetxt = vec1.Pt( jetmom )/jetmom4.M();
	jetz = rp->getEnergy()/jet->getEnergy();
      }


      //calorimter deposits
      ecal = 0;
      hcal = 0;
      std::vector<lcio::Cluster*> clusters = rp->getClusters();
      for ( std::vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();
	    iCluster!=clusters.end();
	    ++iCluster) {
	ecal += (*iCluster)->getSubdetectorEnergies()[0];
	hcal += (*iCluster)->getSubdetectorEnergies()[1];
      }

      //impact parameters
      const EVENT::TrackVec & trkvec = rp->getTracks();
      if (trkvec.size()==0)
	{
	  // TODO: more sophisticated pfo/track matching
	  d0 = 0;
	  z0 = 0;
	  r0 = 0;
	  dEdX= 0;
	  omega= 0;
	  phi= 0;
	  trkType= 0;
	}
      else
	{
	  // TODO: more sophisticated pfo/track matching
	  d0 = fabs(trkvec[0]->getD0());
	  z0 = fabs(trkvec[0]->getZ0());
	  r0 = sqrt( d0*d0 + z0*z0 );
	  dEdX= trkvec[0]->getdEdx();
	  omega= trkvec[0]-> getOmega();
	  phi= trkvec[0]-> getPhi();
	  trkType= trkvec[0]->getType();
	}

      if(d0 < lowestimpactparameter)
	{
	  lowestimpactparameter=d0;
	  bestimpactparameterposition=i;
	}
      
      _vCorrectLeptonYN.push_back(IsCorrectLepton(rp));
      _vPDG.push_back(rp->getType());
      _vEnergy.push_back(rp->getEnergy());
      _vPx.push_back(rp->getMomentum()[0]);
      _vPy.push_back(rp->getMomentum()[1]);
      _vPz.push_back(rp->getMomentum()[2]);
      _vPt.push_back(sqrt(rp->getMomentum()[0]*rp->getMomentum()[0]+rp->getMomentum()[1]*rp->getMomentum()[1]));
      PTot=sqrt(rp->getMomentum()[0]*rp->getMomentum()[0]+rp->getMomentum()[1]*rp->getMomentum()[1]+rp->getMomentum()[2]*rp->getMomentum()[2]);
      _vPTot.push_back(PTot);
      _vCharge.push_back(rp->getCharge());
      _vJetXt.push_back(jetxt);
      _vJetZ.push_back(jetz);
      _vECALToHCALFraction.push_back(ecal/(ecal+hcal));
      _vCalByP.push_back((ecal + hcal)/PTot);
      _vConeEnergy.push_back(getConeEnergy(rp));
      _vD0.push_back(d0);
      _vZ0.push_back(z0);
      _vR0.push_back(r0);
      _vdEdX.push_back(dEdX);
      _vomega.push_back(omega);
      _vphi.push_back(phi);
      _vtrkType.push_back(trkType);
      _vHighestEnergyYN.push_back(0);
      _vLowestImpactParameterYN.push_back(0);

    }
  if(nLeptons>0)
    {
      _vHighestEnergyYN.at(_vHighestEnergyYN.size()-nLeptons+bestenergyposition)=1;
      _vLowestImpactParameterYN.at(_vHighestEnergyYN.size()-nLeptons+bestimpactparameterposition)=1;
    }
  
  std::cout<<"Processing Event "<<nEvt<<std::endl;
  if(foundlep==false){std::cout<<"no lep found in event"<<nEvt<<std::endl;}

  nEvt++;

}

void LeptonOptimizationProcessor::check( LCEvent * evt ) { 
}

void LeptonOptimizationProcessor::end() {

  AddBranch(_tData, "CorrectLeptonYN",_vCorrectLeptonYN );
  AddBranch(_tData, "PDG",_vPDG );
  AddBranch(_tData, "Energy",_vEnergy );
  AddBranch(_tData, "Px",_vPx );
  AddBranch(_tData, "Py",_vPy );
  AddBranch(_tData, "Pz",_vPz );
  AddBranch(_tData, "Pt",_vPt );
  AddBranch(_tData, "PTot",_vPTot );
  AddBranch(_tData, "Charge",_vCharge );
  AddBranch(_tData, "JetXt",_vJetXt );
  AddBranch(_tData, "JetZ",_vJetZ );
  AddBranch(_tData, "ECALToHCALFraction",_vECALToHCALFraction );
  AddBranch(_tData, "CalByP",_vCalByP );
  AddBranch(_tData, "ConeEnergy",_vConeEnergy );
  AddBranch(_tData, "D0",_vD0 );
  AddBranch(_tData, "Z0",_vZ0 );
  AddBranch(_tData, "R0",_vR0 );
  AddBranch(_tData, "dEdX",_vdEdX );
  AddBranch(_tData, "omega",_vomega );
  AddBranch(_tData, "phi",_vphi );
  AddBranch(_tData, "trkType",_vtrkType );
  AddBranch(_tData, "HighestEnergyYN",_vHighestEnergyYN );
  AddBranch(_tData, "LowestImpactParameterYN",_vLowestImpactParameterYN );

  //std::cout<<"Efficiency= "<<100*(float)nReal/191260<<std::endl; //191260 is the maximum findable (non-beamline) number of MC leptons
  std::cout<<"True Efficiency= "<<100*(float)nReal/nEvt<<std::endl; 
  std::cout<<"Purity= "<<100*(float)nReal/(nReal+nFake)<<std::endl;
  std::cout<<"Duplicates= "<<nDuplicate<<std::endl;

  std::cout<<"Real Leptons Passed= "<<nReal<<std::endl;
  std::cout<<"Fake Leptons Passed= "<<nFake<<std::endl;
}


bool LeptonOptimizationProcessor::IsCorrectLepton(ReconstructedParticle *lepton)
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
bool LeptonOptimizationProcessor::IsCorrectLepton(MCParticle *mcp)
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
	      }
	   */	  
	}
    }	      
  return passedtest;
}

void LeptonOptimizationProcessor::AddBranch( TTree* input_tree, std::string title, std::vector<float> data )
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

float LeptonOptimizationProcessor::getConeEnergy( ReconstructedParticle* pfo ) {
	float coneE = 0;

	TVector3 P( pfo->getMomentum() );
	int npfo = _PFOCol->getNumberOfElements();
	for ( int i = 0; i < npfo; i++ ) {
		ReconstructedParticle* pfo_i = dynamic_cast<ReconstructedParticle*>( _PFOCol->getElementAt(i) );

		// don't add itself to the cone energy
		if ( pfo == pfo_i ) continue; 

		TVector3 P_i( pfo_i->getMomentum() );
		float cosTheta = P.Dot( P_i )/(P.Mag()*P_i.Mag());
		if ( cosTheta >= _cosConeAngle )
			coneE += pfo_i->getEnergy(); 
	}

	return coneE;
}
