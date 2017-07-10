#include "TopPurityEfficiencyProcessor.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>

#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCRelation.h>
#include <marlin/Exceptions.h>
#include <UTIL/PIDHandler.h>
#include "UTIL/LCRelationNavigator.h"


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


TopPurityEfficiencyProcessor aTopPurityEfficiencyProcessor ;


TopPurityEfficiencyProcessor::TopPurityEfficiencyProcessor() : Processor("TopPurityEfficiencyProcessor") 
{

  // modify processor description
  _description = "Generates a purity vs efficiency plot for btagging";


  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
  			   "FlavourJetName" , 
  			   "Name of the flavour tag jets"  ,
  			   _FlavourJetsName ,
  			   std::string("RefinedJets"));

  registerInputCollection( LCIO::LCRELATION,
			   "RecoMCTruthLinkCollectionName" , 
			   "Name of the LCRelation collection containing the linker information"  ,
			   _LinkerColName ,
			     std::string("RecoMCTruthLink"));

  registerProcessorParameter( "nDataPoints" , 
			      "number of different btag cut values to examine"  ,
			      _npoints ,
			      float(100));

  registerProcessorParameter( "MinBJetness" , 
			      "fraction of jet particles associated with a b"  ,
			      _minBJetness ,
			      float(0.5));
}

void TopPurityEfficiencyProcessor::init()
{ 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;
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
      _tJetTree = new TTree("JetData", "Tree containing btag information for jets");
    }
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  

  npoints = _npoints;
  _vAnalysisSummary.clear(); // will contain the complete set of results from all events and cut
  _vBTag.clear();
  _vBJetness.clear();
}

void TopPurityEfficiencyProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 



void TopPurityEfficiencyProcessor::processEvent( LCEvent * evt ) 
{

  flavourtype=0;
  _vEventSummary.clear(); // will contain the results from all cuts in one event

  
  FlavourJets = NULL;
   try
    {
      FlavourJets = evt->getCollection( _FlavourJetsName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _FlavourJetsName << " collection not available" << std::endl;
      FlavourJets = NULL;
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
    
   testvector.clear();
   PIDHandler pidh( FlavourJets );
   algo = pidh.getAlgorithmID( "lcfiplus" );
   ibtag = pidh.getParameterIndex(algo, "BTag");
   for(int i=0; i < FlavourJets->getNumberOfElements(); i++) 
     {
       part = dynamic_cast<ReconstructedParticle*>( FlavourJets->getElementAt( i ) );
       const ParticleID &pid = pidh.getParticleID(part, algo);
       btag=pid.getParameters()[ibtag];
       bjetness= BJetness(part); //calculates the fraction of particles in the jet that came from a b quark
       _vBTag.push_back(btag);
       testvector.push_back(btag);
       _vBJetness.push_back(bjetness);
       
       for(int j=0; j<npoints; j++)
	 {
	   _vData.clear();
	   _vData.resize(4,0); // four spaces, entry 0: passes cut + should pass, 1: passes cut + should fail 2: fails + should pass 4: fails + should fail; 
  
	   if(btag>j/npoints) //if passes cut
	     {
	       if(bjetness>_minBJetness) //if should've passed
		 {
		   _vData[0]++;
		 }
	       else
		 {
		   _vData[1]++;
		 }
	     }

	   else
	     {
	       if(bjetness>_minBJetness)
		 {
		   _vData[2]++;
		 }
	       else
		 {
		   _vData[3]++;
		 }
	     }

	   _vEventSummary.push_back(_vData);
	   
	 }

     }

   for(int i=0; i < FlavourJets->getNumberOfElements(); i++)
     {
      
       if(testvector[i]=*std::max_element(testvector.begin(),testvector.end()))
	 {
	   _vHighestBTag.push_back(1);
	 }
       else
	 {
	   _vHighestBTag.push_back(0);
 
	 }
     }
   _vAnalysisSummary.push_back(_vEventSummary);
			       
   
  _nEvt ++ ;
  std::cout << "Events processed: " << _nEvt << std::endl;
}



void TopPurityEfficiencyProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TopPurityEfficiencyProcessor::end(){ 


  AddBranch(_tJetTree, "BTag", _vBTag);
  AddBranch(_tJetTree, "BJetness", _vBJetness);
  AddBranch(_tJetTree, "HighestBTag", _vHighestBTag);


  
  int nEvt= _vAnalysisSummary.size();
  std::cout<<"Size of analysis vector= "<<nEvt<<std::endl;
  std::cout<<"Size of event vector= "<<_vAnalysisSummary[0].size()<<std::endl;
  std::cout<<"First Element is "<<_vAnalysisSummary[0][0][0]<<std::endl;
  int size= (int)npoints;
  float _aPurity[size];
  float _aEfficiency[size];
  float _aBCut [size] ;
  
  for(int j=0; j<npoints; j++) //loop through cuts to get data points
    {
      totalPassShouldPass=0;
      totalPassShouldFail=0;
      totalFailShouldPass=0;
      totalFailShouldFail=0;

      for(int i=0; i<nEvt; i++) //event loop
	{
	  totalPassShouldPass+=_vAnalysisSummary[i][j][0];
	  totalPassShouldFail+=_vAnalysisSummary[i][j][1];
	  totalFailShouldPass+=_vAnalysisSummary[i][j][2];
	  totalFailShouldFail+=_vAnalysisSummary[i][j][3];
	}
      
      _aEfficiency[j]=totalPassShouldPass/(totalPassShouldPass+totalFailShouldPass);
      _aPurity[j]=totalPassShouldPass/(totalPassShouldPass+totalPassShouldFail);
      _aBCut[j]=j/npoints;
    }

  
  TGraph *P = new TGraph(size, _aBCut,_aPurity );
  P->Draw("");
  P->SetTitle("Purity vs Btag");
  P->GetXaxis()->SetTitle("BCut");
  P->GetYaxis()->SetTitle("Purity");
  P->SetName("PurityVsBCut");
  P->Write();

  TGraph *E = new TGraph(size, _aBCut, _aEfficiency);
  E->Draw("");
  E->SetTitle("Efficiency vs Btag");
  E->GetXaxis()->SetTitle("BCut");
  E->GetYaxis()->SetTitle("Efficiency");
  E->SetName("EfficiencyVsBCut");
  E->Write();

  TGraph *PvE = new TGraph(size, _aPurity, _aEfficiency);
  PvE->SetTitle("Purity vs Efficiency");
  PvE->GetXaxis()->SetTitle("Purity");
  PvE->GetYaxis()->SetTitle("Efficiency");
  PvE->SetName("PurityVsEfficiency");
  PvE->SetMarkerStyle(24);
  PvE->SetMarkerSize(1);
  PvE->Draw("ALP");
  PvE->Write();
  
  // std::cout << "TopPurityEfficiencyProcessor::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
  
}



void TopPurityEfficiencyProcessor::AddBranch( TTree* input_tree, std::string title, std::vector<float> data )
{
  if( (input_tree->GetEntries() != (int)data.size()) && (input_tree->GetEntries() != 0) )
    {
      std::cout << "ERROR: CANNOT ADD BRANCH " << title << " DATA SIZE DID NOT MATCH TREE." 
		<< "Data size: " << data.size() << std::endl;
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

float TopPurityEfficiencyProcessor::BJetness( ReconstructedParticle* jet)
{
  float totaljetsize= 0;
  float nparticlesfromb=0;

  for(unsigned int i=0; i<jet->getParticles().size(); i++)
    {
      jetparticle= jet->getParticles()[i];
      LCRelationNavigator* relationNavigatorPFOMC = new LCRelationNavigator( LinkerCol );
      relobjMC.clear();
      relobjMC = relationNavigatorPFOMC->getRelatedToObjects(jetparticle);

      for(unsigned int j=0; j<relobjMC.size(); j++)
	{
	  jetmcp=NULL;
	  jetmcp=dynamic_cast <MCParticle*>(relobjMC[j]);
	  if(jetmcp==NULL){continue;}
	  else
	    {
	      totaljetsize++;
	      for (MCParticle* TestMCP = jetmcp; TestMCP->getParents().size()>0; TestMCP=TestMCP->getParents()[0]) 
		{
		  if(abs(TestMCP->getPDG())==5)
		    {
		      nparticlesfromb++;
		      break;
		    }
		}
	    }
	}
      delete relationNavigatorPFOMC;
    }

  return nparticlesfromb/totaljetsize;
  
}
