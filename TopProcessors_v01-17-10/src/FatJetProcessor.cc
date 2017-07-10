#include "FatJetProcessor.h"

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


FatJetProcessor aFatJetProcessor ;


FatJetProcessor::FatJetProcessor() : Processor("FatJetProcessor") 
{
  // modify processor description
  _description = "Determine which fat jet is from the hadronic decay of a top " ;


  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "FatJetCollectionName" , 
			   "Name of the two jet collection"  ,
			   _FatJetColName ,
			   std::string("FatJets"));

  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "IsoLep",
			    "IsoLep",
			    _IsoLepColName,
			    std::string("IsoLep") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
  			   "FlavourJetName" , 
  			   "Name of the flavour tag jets"  ,
  			   _FlavourJetsName ,
  			   std::string("RefinedJets"));

  registerProcessorParameter( "SortingType",
			      "How to decide which jet is the hadronic jet",
			      _SortingType,
			      std::string("Energy"));

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "HadronicFatJet",
			    "HadronicFatJet",
			    _HadronicFatJetCol,
			    std::string("HadronicFatJet") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "HadronicFatJetConstituents",
			    "HadronicFatJetConstituents",
			    _HadronicFatJetConstituentsCol,
			    std::string("HadronicFatJetConstituents") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "LeptonicFatJet",
			    "LeptonicFatJet",
			    _LeptonicFatJetCol,
			    std::string("LeptonicFatJet") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "LeptonicFatJetConstituents",
			    "LeptonicFatJetConstituents",
			    _LeptonicFatJetConstituentsCol,
			    std::string("LeptonicFatJetConstituents") );
    

  

}


void FatJetProcessor::init() 
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


void FatJetProcessor::processRunHeader( LCRunHeader* run) 
{
  _nRun++ ;
} 


void FatJetProcessor::processEvent( LCEvent * evt ) 
{

  
/////////////////////////////////////// Get Collections ////////////////////////////////////////////
  
  _AllCollectionsExist = true ;

  IsoLepCol = NULL;
  try
    {
      IsoLepCol = evt->getCollection( _IsoLepColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _IsoLepColName << " collection not available" << std::endl;
      IsoLepCol=NULL;
      _AllCollectionsExist = false ;
    }

  
  FatJetCol = NULL;
  try
    {
      FatJetCol = evt->getCollection( _FatJetColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _FatJetColName << " collection not available" << std::endl;
      FatJetCol=NULL;
      _AllCollectionsExist = false ;
    }

  FlavourJets=NULL;
  try
    {
      FlavourJets = evt->getCollection( _FlavourJetsName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _FlavourJetsName << " collection not available" << std::endl;
      FlavourJets = NULL;
    }
  
  
  if( !_AllCollectionsExist )
    {
      streamlog_out(WARNING) << "Very sad...Skipping event " << _nEvt 
			     << " - Cannot find all required collections" << std::endl;
      throw SkipEventException( this );

    }

  //set up output collections
  HadronicFatJetCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
  HadronicFatJetCol->setSubset(true);

  HadronicFatJetConstituentsCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
  HadronicFatJetConstituentsCol->setSubset(true);

  LeptonicFatJetCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
  LeptonicFatJetCol->setSubset(true);

  LeptonicFatJetConstituentsCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
  LeptonicFatJetConstituentsCol->setSubset(true);
  
  // Process Event
  if(FatJetCol->getNumberOfElements()>1)
    {
      JetA = dynamic_cast<ReconstructedParticle*>( FatJetCol->getElementAt(0) );
      JetB = dynamic_cast<ReconstructedParticle*>( FatJetCol->getElementAt(1) );
    }
  else
    {
      std::cout<<"This event just isn't fat enough!!"<<std::endl;
      throw SkipEventException( this );
    }
  Lepton=dynamic_cast<ReconstructedParticle*>( IsoLepCol->getElementAt(0) );

  // get PIDHandler associated with the jet collection
  PIDHandler pidh( FlavourJets );
  // get algorithm ID associated with LCFIPlus
  int algo = pidh.getAlgorithmID( "lcfiplus" );
  // get index number for flavor tagging
  int ibtag = pidh.getParameterIndex(algo, "BTag");
  int ictag = pidh.getParameterIndex(algo, "CTag");
  BJet=dynamic_cast<ReconstructedParticle*>( FlavourJets->getElementAt( 0 ) );
  float bestBTag=0.0;
  for(int i=1; i < FlavourJets->getNumberOfElements(); i++) 
    {
      ReconstructedParticle *part = dynamic_cast<ReconstructedParticle*>( FlavourJets->getElementAt( i ) );
      const ParticleID &pid = pidh.getParticleID(part, algo);
      float btag=pid.getParameters()[ibtag];
      if(btag>bestBTag)
	{
	  BJet=part;
	  bestBTag=btag;
	}
    }
  
  _vJetDIRA.push_back(dira(JetA,JetB));


  //perform sorting
  if(_SortingType=="lepton")
    {

      if(abs(dira(JetA,Lepton))>abs(dira(JetB,Lepton)))
	{

	  HadronicJet=JetA;
	  LeptonicJet=JetB;
	}

     else
       {
	  HadronicJet=JetB;
	  LeptonicJet=JetA;
       }
	  
    }

  //perform sorting
  else if(_SortingType=="btag")
    {

      if(abs(dira(JetA,BJet))>abs(dira(JetB,BJet)))
	{

	  HadronicJet=JetA;
	  LeptonicJet=JetB;
	}

      else
	{
	  HadronicJet=JetB;
	  LeptonicJet=JetA;
	}
	  
    }

  else if(_SortingType=="energy")
    {

      if(JetA->getEnergy()>JetB->getEnergy())
	{

	  HadronicJet=JetA;
	  LeptonicJet=JetB;
	}

     else
       {
	  HadronicJet=JetB;
	  LeptonicJet=JetA;
       }
	  
    }

  else if(_SortingType=="particlenumber")
   {

     if(JetA->getParticles().size()>JetB->getParticles().size())
	{

	  HadronicJet=JetA;
	  LeptonicJet=JetB;
	}

     else
       {
	  HadronicJet=JetB;
	  LeptonicJet=JetA;
       }
	  
    }

  else if(_SortingType=="mass")
   {
     if(JetA->getMass()>JetB->getMass())
	{

	  HadronicJet=JetA;
	  LeptonicJet=JetB;
	}

     else
       {
	  HadronicJet=JetB;
	  LeptonicJet=JetA;

	  
       }
   }

 else if(_SortingType=="topmass")
   {
     if(abs(JetA->getMass()-174)<abs(JetB->getMass()-174))
	{

	  HadronicJet=JetA;
	  LeptonicJet=JetB;
	}

     else
       {
	  HadronicJet=JetB;
	  LeptonicJet=JetA;

	  
       }
   }

  else  if(_SortingType=="random")
   {
     if(_nEvt%2==1)
       {
	 HadronicJet=JetA;
	 LeptonicJet=JetB;
       }
     else
       {
	 HadronicJet=JetB;
	 LeptonicJet=JetA;
       }
   }
  else{std::cout<<"Sorting method not recognised"<<std::endl;}
    
  _vHadronicTopJetMass.push_back(HadronicJet->getMass());
  _vLeptonicTopJetMass.push_back(LeptonicJet->getMass());

  //fill collections
  HadronicFatJetCol->addElement(HadronicJet);
  LeptonicFatJetCol->addElement(LeptonicJet);

  for(int i=0; i<HadronicJet->getParticles().size();i++)
    {
      HadronicFatJetConstituentsCol->addElement(HadronicJet->getParticles()[i]);
    }
  for(int i=0; i<LeptonicJet->getParticles().size();i++)
    {
      LeptonicFatJetConstituentsCol->addElement(LeptonicJet->getParticles()[i]);
    }


 
  //add collections to events
  evt->addCollection( HadronicFatJetCol, _HadronicFatJetCol.c_str() );
  evt->addCollection( HadronicFatJetConstituentsCol, _HadronicFatJetConstituentsCol.c_str() );
  evt->addCollection( LeptonicFatJetCol, _LeptonicFatJetCol.c_str() );
  evt->addCollection( LeptonicFatJetConstituentsCol, _LeptonicFatJetConstituentsCol.c_str() );
  
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		       << "   in run:  " << evt->getRunNumber() << std::endl ;
  
  _nEvt ++ ;
      
  std::cout << "Events processed: " << _nEvt << std::endl;
  
}

void FatJetProcessor::check( LCEvent * evt ) 
{ 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void FatJetProcessor::end()
{ 
    
  // In here we fill the trees with the information stored within the vectors
  // Note: All of the vectors must be the same size else the branch will not be added

  //HadronicTop
  AddBranch(_tRecoMasses, "HadronicJetMass", _vHadronicTopJetMass);
  AddBranch(_tRecoMasses, "LeptonicJetMass", _vLeptonicTopJetMass);
  AddBranch(_tRecoMasses, "JetDira", _vJetDIRA);
}






void FatJetProcessor::AddBranch( TTree* input_tree, std::string title, std::vector<float> data )
{
  if( (input_tree->GetEntries() != (int)data.size()) && (input_tree->GetEntries() != 0) )
    {
      std::cout << "ERROR: CANNOT ADD BRANCH " << title << " DATA SIZE DID NOT MATCH TREE." 
      	<< "Data size: " << data.size() << ". Tree size: "<< _vHadronicTopJetMass.size() << std::endl;
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

double FatJetProcessor::dira(ReconstructedParticle *p1, ReconstructedParticle *p2)
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
