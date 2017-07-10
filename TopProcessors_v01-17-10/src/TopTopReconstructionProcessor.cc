#include "TopTopReconstructionProcessor.h"

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


TopTopReconstructionProcessor aTopTopReconstructionProcessor ;


TopTopReconstructionProcessor::TopTopReconstructionProcessor() : Processor("TopTopReconstructionProcessor") 
{
  // modify processor description
  _description = "TopTopReconstructionProcessor reconstructs Top pairs assuming a lvbqqb final state " ;


  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter( "EventEnergy",
			      "Energy of the event",
			      eventEnergy,
			      double(380));

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCCollectionName" , 
			   "Name of the MCParticle collection"  ,
			   _MCColName ,
			   std::string("MCParticlesSkimmed"));
	
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "TwoJetCollectionName" , 
			   "Name of the two jet collection"  ,
			   _TwoJetColName ,
			   std::string("Durham_2Jets"));

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "PFOCollectionName" , 
			   "Name of the PFO collection"  ,
			   _PFOColName ,
			   std::string("PandoraPFANewPFOs"));

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE	,
			   "PFOWithoutIsolepCollectionName",
			   "Name of the PFO collection without Isolated Leptons",
			   _PFONoIsoColName ,
			   std::string("PandoraPFOsWithoutIsolep"));

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE ,
			   "ISOLEPCollectionName" , 
			   "Name of the Isolated Lepton collection" ,
			   _IsoLepColName ,
			   std::string("Isolep"));

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "FourJetCollectionName" , 
			   "Name of the four jet collection"  ,
			   _FourJetColName ,
			   std::string("Durham_4Jets"));

  registerInputCollection( LCIO::LCRELATION,
			   "RecoMCTruthLinkCollectionName" , 
			   "Name of the LCRelation collection containing the linker information"  ,
			   _LinkerColName ,
			   std::string("RecoMCTruthLink"));

  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "FlavourJetName" , 
			   "Name of the flavour tag jets"  ,
			   _FlavourJetColName ,
			   std::string("RefinedJets"));

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "FakeJetColName" , 
			   "Name of the flavour tag jets"  ,
			   _FakeJetColName ,
			   std::string("RefinedJets"));

  registerProcessorParameter( "DoFlavourTagging",
			      "Use flavour information for analysis",
			      _DoFlavourTagging,
			      bool(false));

  registerProcessorParameter( "UseChiSquaredMatching",
			      "Use flavour information for analysis",
			      _UseChiSquaredMatching,
			      bool(false));
    

}


void TopTopReconstructionProcessor::init() 
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
      _tTestTree = new TTree("TestTree", "Tree containing final reconstructed masses");
    }

}


void TopTopReconstructionProcessor::processRunHeader( LCRunHeader* run) 
{
  _nRun++ ;
} 


void TopTopReconstructionProcessor::processEvent( LCEvent * evt ) 
{
  //clear variables 
  nPFOs=0;
  nIsolatedLeptons=0;
  nLinks=0;
  nMCP=0;
  nJets=0;
  nFlavourJets=0;
  
  pthval=0;
  majthval=0;
  minthval=0;
  SigPx=0;
  SigPy=0;
  SigPz=0;
  SigPt=0;
  SigE=0;
  MissPx=0;
  MissPy=0;
  MissPz=0;
  MissPt=0;
  MissE=0;
  Y12=0;
  Y23=0;
  Y34=0;
  Y45=0;

  Wqq=NULL;
  isolep=NULL;
  pseudoHiggs=NULL;

  jets.clear();
  _vBPrimaryTags.clear();
  _vCPrimaryTags.clear();
  
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

  TwoJetCol = NULL;
  try
    {
      TwoJetCol = evt->getCollection( _TwoJetColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _TwoJetColName << " collection not available" << std::endl;
      TwoJetCol=NULL;
      _AllCollectionsExist = false ;
    }

  PFOCol = NULL;
  try
    {
      PFOCol = evt->getCollection( _PFOColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _PFOColName << " collection not available" << std::endl;
      PFOCol = NULL;
      _AllCollectionsExist = false ;
    }

  PFONoIsoCol = NULL;
  try
    {
      PFONoIsoCol = evt->getCollection( _PFONoIsoColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _PFONoIsoColName << " collection not available" << std::endl;
      PFONoIsoCol = NULL;
      _AllCollectionsExist = false ;
    }

  IsoLepCol = NULL;
  try
    {
      IsoLepCol = evt->getCollection( _IsoLepColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _IsoLepColName << " collection not available" << std::endl;
      IsoLepCol = NULL;
      _AllCollectionsExist = false ;
    }

  FourJetCol = NULL;
  try
    {
      FourJetCol = evt->getCollection( _FourJetColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _FourJetColName << " collection not available" << std::endl;
      FourJetCol = NULL;
      _AllCollectionsExist = false ;
    }

  FakeJetCol = NULL;
  try
    {
      FakeJetCol = evt->getCollection( _FakeJetColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _FakeJetColName << " collection not available" << std::endl;
      FakeJetCol = NULL;
      _AllCollectionsExist = false ;
    }
  
  if(_DoFlavourTagging)
    {
      FlavourJetCol = NULL;
      try
	{
	  FlavourJetCol = evt->getCollection( _FlavourJetColName );
	}
      catch( lcio::DataNotAvailableException e )
	{
	  streamlog_out(WARNING) << _FlavourJetColName << " collection not available" << std::endl;
	  FlavourJetCol = NULL;
	  _AllCollectionsExist = false ;
	}
      nFlavourJets=FlavourJetCol->getNumberOfElements();
    }
  
  LinkerCol = NULL;
  try
    {
      LinkerCol = evt->getCollection( _LinkerColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _LinkerColName << " collection not available" << std::endl;
      LinkerCol = NULL;
      _AllCollectionsExist = false ;
    }

  nLinks = LinkerCol->getNumberOfElements() ;	


  if( !_AllCollectionsExist )
    {
      streamlog_out(WARNING) << "Very sad...Skipping event " << _nEvt 
			     << " - Cannot find all required collections" << std::endl;
      throw SkipEventException( this );

    }

  //get size of collections
   nPFOs=PFONoIsoCol->getNumberOfElements();
   nIsolatedLeptons = IsoLepCol->getNumberOfElements() ;
   nMCP = MCCol->getNumberOfElements();
  
  //////////////////////////////////// MCTruth Information ///////////////////////////////////////////////

  //////////////////////////////// Reconstructed Particle Information ///////////////////////////////////

   _vnPFOs.push_back(nPFOs);
   
  pthval=0., majthval=0., minthval=0.;

  //PFO Collection- thrust parameters
  LCParameters& pfocolParameters((LCParameters &)PFOCol->parameters() );
  pthval = pfocolParameters.getFloatVal( "principleThrustValue" ); 
  majthval = pfocolParameters.getFloatVal( "majorThrustValue" );
  minthval = pfocolParameters.getFloatVal( "minorThrustValue" );
  _vpthval.push_back( pthval );
  _vmajthval.push_back( majthval );
  _vminthval.push_back( minthval );
      
  pthval=0., majthval=0., minthval=0.; //reset values
  
  //PFO Collection without isolep- thrust parameters
  LCParameters& pfoisocolParameters((LCParameters &)PFONoIsoCol->parameters() );
  pthval 	 = pfoisocolParameters.getFloatVal( "principleThrustValue" ); 
  majthval = pfoisocolParameters.getFloatVal( "majorThrustValue" );
  minthval = pfoisocolParameters.getFloatVal( "minorThrustValue" );
  _vpthnisoval.push_back( pthval );
  _vmajthnisoval.push_back( majthval );
  _vminthnisoval.push_back( minthval );

  //Jet Parameters
  LCParameters& col2jParameters((LCParameters &)TwoJetCol->parameters() );
  LCParameters& col4jParameters((LCParameters &)FourJetCol->parameters() );
  Y12 = col2jParameters.getFloatVal( "y_{n-1,n}" );
  Y23 = col2jParameters.getFloatVal( "y_{n,n+1}" );
  Y34 = col4jParameters.getFloatVal( "y_{n-1,n}" );
  Y45 = col4jParameters.getFloatVal( "y_{n,n+1}" );
  _vY12.push_back(-(std::log(Y12)));
  _vY23.push_back(-(std::log(Y23)));
  _vY34.push_back(-(std::log(Y34)));
  _vY45.push_back(-(std::log(Y45)));


  ////////////////////////////////IsoLepInfo////////////////////////////////////////

  if (nIsolatedLeptons>0)
    {
      leptonfound=false;
      isolep = dynamic_cast<ReconstructedParticle*>( IsoLepCol->getElementAt(0) );
      if (IsCorrectLepton(isolep)==true) {_vCorrectLeptonYN.push_back(1);}
      else{_vCorrectLeptonYN.push_back(-1);}
      _vLEPMass.push_back(isolep->getMass());
      _vLEPEnergy.push_back(isolep->getEnergy());
      isolepmom = isolep->getMomentum();
      _vLEPMomx.push_back(isolepmom[0]);
      _vLEPMomy.push_back(isolepmom[1]);
      _vLEPMomz.push_back(isolepmom[2]);
      _vLEPMomt.push_back(pow((pow(isolepmom[0],2)+pow(isolepmom[1],2)),0.5));
      _vLEPdira.push_back(atan((pow((pow(isolepmom[0],2)+pow(isolepmom[1],2)),0.5))/isolepmom[2]));
      _vLEPCharge.push_back(isolep->getCharge());
      _vLEPID.push_back(isolep->getType());
    }

  else
    {
      _vCorrectLeptonYN.push_back(0);
      _vLEPMass.push_back(-0.1);
      _vLEPEnergy.push_back(-5.);
      _vLEPMomx.push_back(-500.);
      _vLEPMomy.push_back(-500.);
      _vLEPMomz.push_back(-500.);
      _vLEPMomt.push_back(-5.);
      _vLEPdira.push_back(-7.);
      _vLEPCharge.push_back(-5.);
      _vLEPID.push_back(0);
    }


  ///////////////////////// B and W Reconstruction ////////////////////////////////////////

  nJets = FlavourJetCol->getNumberOfElements();
  jets.resize(nJets);
  
  std::map<float,ReconstructedParticle*> BTagMap;
  
  for(int i=0; i< nJets ; i++)
    {
      ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( FlavourJetCol->getElementAt( i ) ) ;
      _vJetMass.push_back(jet->getMass());
      _vJetPt.push_back(sqrt(jet->getMomentum()[0]*jet->getMomentum()[0]+jet->getMomentum()[1]*jet->getMomentum()[1]));
      _vJetDIRA.push_back(atan((pow((pow(jet->getMomentum()[0],2)+pow(jet->getMomentum()[1],2)),0.5))/jet->getMomentum()[2]));
      jets[i]=jet;
      if(_UseChiSquaredMatching==false)
	{
	  _vBPrimaryTags.push_back(GetBTag(jet));
	  _vCPrimaryTags.push_back(GetCTag(jet));
	  _vAllBPrimaryTags.push_back(GetBTag(jet));
	  BTagMap.insert ( std::pair<float,ReconstructedParticle*>(_vBPrimaryTags.back(),jet) );
	}
    }

  if(_UseChiSquaredMatching==false)
    {
      jets.clear();
      //sort jets by flavour tag in ascending order. Assume lowest two are the W jets
      std::sort (_vBPrimaryTags.begin(), _vBPrimaryTags.end());
      for(int i=0; i<4; i++)
	{
	  jets[i]=BTagMap.find(_vBPrimaryTags[i])->second; 
	}
    }
  
  std::vector<int> jetindex;
  Wqq = NULL;
  ReconstructedParticle* HadronicTop=NULL;
  ReconstructedParticle* LoneBJet=NULL;
  float topMass=173.21;
  float WMass=80.385;

  if(_UseChiSquaredMatching==false)// if combining lowest btag jets produces a sensible W Mass
    {
      //assign lowest btagged jets to be from the W
      Wqq = CombineParticles(jets[0], jets[1]);

      //find which third jet produces a sensible top mass
      if(abs(CombineParticles(Wqq,jets[2])->getMass()-topMass) > abs(CombineParticles(Wqq,jets[3])->getMass()-topMass))
	{
	  LoneBJet= jets[2];
	  HadronicTop=CombineParticles(Wqq,jets[3]);
	}
      else
	{
	  LoneBJet= jets[3];
	  HadronicTop=CombineParticles(Wqq,jets[2]);
	}
    }
  
  else
    {
      //perform a chi squared fit of the W and Top Mass for all jet combinations
      jetindex= SortJetVector(jets); //first two indices are W jet positions, 3rd is hadronic b, 4th is leptonic b ("lone b")
      Wqq=CombineParticles(jets[jetindex[0]],jets[jetindex[1]]);
      HadronicTop=CombineParticles(Wqq,jets[jetindex[2]]);
      LoneBJet=jets[jetindex[3]];
    }

      
  //Get hadronic W, top and b jet properties
  _vWqqMass.push_back(Wqq->getMass());
  _vWqqEnergy.push_back(Wqq->getEnergy());
  const double *wqqmom = Wqq->getMomentum();
  _vWqqMomx.push_back(wqqmom[0]);
  _vWqqMomy.push_back(wqqmom[1]);
  _vWqqMomz.push_back(wqqmom[2]);
  _vWqqMomt.push_back(pow((pow(wqqmom[0],2)+pow(wqqmom[1],2)),0.5));
  _vWqqdira.push_back(atan((pow((pow(wqqmom[0],2)+pow(wqqmom[1],2)),0.5))/wqqmom[2]));

  //TopProperties
  _vHadronicTopMass.push_back(HadronicTop->getMass());
  _vHadronicTopEnergy.push_back(HadronicTop->getEnergy());
  _vHadronicTopPx.push_back(HadronicTop->getMomentum()[0]);
  _vHadronicTopPy.push_back(HadronicTop->getMomentum()[1]);
  _vHadronicTopPz.push_back(HadronicTop->getMomentum()[2]);

  _vBJetMass.push_back(LoneBJet->getMass());

  _vWqqTopSeperation.push_back(dira(Wqq,HadronicTop));
  _vLoneBJetTopSeperation.push_back(dira(LoneBJet,HadronicTop));
  
  //analyse jets in general
  nParticlesInJets=0;
  for(int i=0; i< nJets ; i++)
    {
      nParticlesInJets+=jets[i]->getParticles().size();
      SigE+=jets[i]->getEnergy();
      SigPx+=jets[i]->getMomentum()[0];
      SigPy+=jets[i]->getMomentum()[1];
      SigPz+=jets[i]->getMomentum()[2];
    }

  _vParticlesInJets.push_back(nParticlesInJets);
  SigE+=isolep->getEnergy();
  SigPx+=isolep->getMomentum()[0];
  SigPy+=isolep->getMomentum()[1];
  SigPz+=isolep->getMomentum()[2];
  MissE=eventEnergy-SigE;
  _vMissE.push_back(MissE);
  _vMissPx.push_back(-SigPx);
  _vMissPy.push_back(-SigPy);
  _vMissPz.push_back(-SigPz);

  //Reconstruct the neutrino based on missing energy and momentum

  ReconstructedParticleImpl* neutrino = new ReconstructedParticleImpl();
  neutrino->setEnergy( MissE );
  double Momentum_dummy [3] = {-SigPx,-SigPy,-SigPx};
  neutrino->setMomentum(Momentum_dummy );

  ReconstructedParticle* LeptonicW= CombineParticles(isolep,neutrino);
  _vLeptonicWMass.push_back(LeptonicW->getMass());
  _vLeptonicWEnergy.push_back(LeptonicW->getEnergy());
  _vLeptonicWPx.push_back(LeptonicW->getMomentum()[0]);
  _vLeptonicWPy.push_back(LeptonicW->getMomentum()[1]);
  _vLeptonicWPz.push_back(LeptonicW->getMomentum()[2]);
  
  ReconstructedParticle* LeptonicTop= CombineParticles(LeptonicW,LoneBJet);
  _vLeptonicTopMass.push_back(LeptonicTop->getMass());
  _vLeptonicTopEnergy.push_back(LeptonicTop->getEnergy());
  _vLeptonicTopPx.push_back(LeptonicTop->getMomentum()[0]);
  _vLeptonicTopPy.push_back(LeptonicTop->getMomentum()[1]);
  _vLeptonicTopPz.push_back(LeptonicTop->getMomentum()[2]);

  //info from mcjets can be added in here??
  
 
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		       << "   in run:  " << evt->getRunNumber() << std::endl ;
  
  _nEvt ++ ;
      
  std::cout << "Events processed: " << _nEvt << std::endl;
  
}

void TopTopReconstructionProcessor::check( LCEvent * evt ) 
{ 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TopTopReconstructionProcessor::end()
{ 
    
  // In here we fill the trees with the information stored within the vectors
  // Note: All of the vectors must be the same size else the branch will not be added

  //HadronicTop
  AddBranch(_tRecoMasses, "HadronicTopMass", _vHadronicTopMass);
  AddBranch(_tRecoMasses, "HadronicTopEnergy", _vHadronicTopEnergy);
  AddBranch(_tRecoMasses, "HadronicTopPx", _vHadronicTopPx);
  AddBranch(_tRecoMasses, "HadronicTopPy", _vHadronicTopPy);
  AddBranch(_tRecoMasses, "HadronicTopPz", _vHadronicTopPz);
  AddBranch(_tRecoMasses, "WqqTopSeperation", _vWqqTopSeperation);
  AddBranch(_tRecoMasses, "LoneBJetTopSeperation", _vLoneBJetTopSeperation);

  //Leptonic Top
  AddBranch(_tRecoMasses, "LeptonicTopMass", _vLeptonicTopMass);
  AddBranch(_tRecoMasses, "LeptonicTopEnergy", _vLeptonicTopEnergy);
  AddBranch(_tRecoMasses, "LeptonicTopPx", _vLeptonicTopPx);
  AddBranch(_tRecoMasses, "LeptonicTopPy", _vLeptonicTopPy);
  AddBranch(_tRecoMasses, "LeptonicTopPz", _vLeptonicTopPz);

  //IsoLep
  AddBranch(_tRecoMasses, "CorrectLeptonYN",_vCorrectLeptonYN );
  AddBranch(_tRecoMasses, "M_Lep", _vLEPMass);
  AddBranch(_tRecoMasses, "E_Lep", _vLEPEnergy);
  AddBranch(_tRecoMasses, "Px_Lep", _vLEPMomx);
  AddBranch(_tRecoMasses, "Py_Lep", _vLEPMomy);
  AddBranch(_tRecoMasses, "Pz_Lep", _vLEPMomz);
  AddBranch(_tRecoMasses, "Pt_Lep", _vLEPMomt);
  AddBranch(_tRecoMasses, "DIRA_Lep", _vLEPdira);
  AddBranch(_tRecoMasses, "Q_Lep", _vLEPCharge);
  AddBranch(_tRecoMasses, "ID_Lep", _vLEPID);

  //Hadronic W Boson
  AddBranch(_tRecoMasses, "HadronicWMass", _vWqqMass);
  AddBranch(_tRecoMasses, "HadronicWEnergy", _vWqqEnergy);
  AddBranch(_tRecoMasses, "HadronicWPx", _vWqqMomx);
  AddBranch(_tRecoMasses, "HadronicWPy", _vWqqMomy);
  AddBranch(_tRecoMasses, "HadronicWPz", _vWqqMomz);
  AddBranch(_tRecoMasses, "HadronicWPt", _vWqqMomt);
  AddBranch(_tRecoMasses, "HadronicWDIRA", _vWqqdira);

  //Leptonic W Mass
  AddBranch(_tRecoMasses, "LeptonicWMass", _vLeptonicWMass);
  AddBranch(_tRecoMasses, "LeptonicWEnergy", _vLeptonicWEnergy);
  AddBranch(_tRecoMasses, "LeptonicWPx", _vLeptonicWPx);
  AddBranch(_tRecoMasses, "LeptonicWPy", _vLeptonicWPy);
  AddBranch(_tRecoMasses, "LeptonicWPz", _vLeptonicWPz);

  //Jets
  AddBranch(_tRecoMasses, "log_Y12", _vY12);
  AddBranch(_tRecoMasses, "log_Y23", _vY23);
  AddBranch(_tRecoMasses, "log_Y34", _vY34);
  AddBranch(_tRecoMasses, "log_Y45", _vY45);
  AddBranch(_tRecoMasses, "ParticlesInJets",_vParticlesInJets );

  //Thrusts
  AddBranch(_tRecoMasses, "principleTh_Val", _vpthval);
  AddBranch(_tRecoMasses, "majorTh_Val", _vmajthval);	
  AddBranch(_tRecoMasses, "minorTh_Val", _vminthval);	
  AddBranch(_tRecoMasses, "principleThnISO_Val", _vpthnisoval);
  AddBranch(_tRecoMasses, "majorThnISO_Val", _vmajthnisoval);
  AddBranch(_tRecoMasses, "minorThnISO_Val", _vminthnisoval);

  //Missing
  AddBranch(_tRecoMasses, "E_Missing", _vMissE);
  AddBranch(_tRecoMasses, "Px_Missing", _vMissPx);
  AddBranch(_tRecoMasses, "Py_Missing", _vMissPy);
  AddBranch(_tRecoMasses, "Pz_Missing", _vMissPz);
  
  //Assorted
  AddBranch(_tRecoMasses, "n_PFOs", _vnPFOs);
  AddBranch(_tRecoMasses, "BJetMass", _vBJetMass);

  //flavour tags
  if(_DoFlavourTagging)
    {
      AddBranch(_tTestTree, "BTag",_vAllBPrimaryTags);
      AddBranch(_tTestTree, "CTag",_vCPrimaryTags);
    }
   //combined Leptonic and Hadronic Top Mass
  _vHadronicTopMass.reserve(_vHadronicTopMass.size()+_vLeptonicTopMass.size());
  _vHadronicTopMass.insert(_vHadronicTopMass.end(), _vLeptonicTopMass.begin(),_vLeptonicTopMass.end()); 
  // AddBranch(_tTestTree, "AllTopMass",_vHadronicTopMass);
}



/////////////////////////////////////// Functions ////////////////////////////////////


//calculates the mass of two particles
double TopTopReconstructionProcessor::ReconstructedMass( ReconstructedParticle* p1, ReconstructedParticle* p2 )
{
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  double mass_squared = pow( p1->getEnergy()+p2->getEnergy(), 2 ) 
    - pow( mom1[0]+mom2[0],2 ) 
    - pow( mom1[1]+mom2[1],2 )
    - pow( mom1[2]+mom2[2],2 );
											
  return pow( mass_squared, 0.5 );
}

//function that takes 2 particles and combines them calculating the new energy, mass and momentum
ReconstructedParticle* TopTopReconstructionProcessor::CombineParticles( ReconstructedParticle *p1, ReconstructedParticle *p2 )
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
	
  //add the particles in
  newParticle->addParticle( p1 );
  newParticle->addParticle( p2 );

  //now return the jet
  delete Momentum;
  return newParticle ;
}

//aw ->  Overload CombineParticle and ReconstructedMass functions to combine any combination of MCParticle and Reconstructed Particle (jet+neutrinos) 
double TopTopReconstructionProcessor::ReconstructedMass( ReconstructedParticle* p1, MCParticle* p2 )
{
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  double mass_squared = pow( p1->getEnergy()+p2->getEnergy(), 2 ) 
    - pow( mom1[0]+mom2[0],2 ) 
    - pow( mom1[1]+mom2[1],2 )
    - pow( mom1[2]+mom2[2],2 );
											
  return pow( mass_squared, 0.5 );
}

ReconstructedParticle* TopTopReconstructionProcessor::CombineParticles( ReconstructedParticle *p1, MCParticle *p2 )
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

double TopTopReconstructionProcessor::ReconstructedMass( MCParticle* p1, MCParticle* p2 )
{
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  double mass_squared = pow( p1->getEnergy()+p2->getEnergy(), 2 ) 
    - pow( mom1[0]+mom2[0],2 ) 
    - pow( mom1[1]+mom2[1],2 )
    - pow( mom1[2]+mom2[2],2 );
											
  return pow( mass_squared, 0.5 );
}

ReconstructedParticle* TopTopReconstructionProcessor::CombineParticles( MCParticle *p1, MCParticle *p2 )
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

ReconstructedParticle* TopTopReconstructionProcessor::ConvertToReconstructedParticle(MCParticle *p1 )
{
  ReconstructedParticleImpl* newParticle = new ReconstructedParticleImpl();
  newParticle->setMass( p1->getMass() );
  newParticle->setEnergy( p1->getEnergy() );
  newParticle->setMomentum( p1->getMomentum() );
  return newParticle;
}
//aw end

void TopTopReconstructionProcessor::AddBranch( TTree* input_tree, std::string title, std::vector<float> data )
{
  if( (input_tree->GetEntries() != (int)data.size()) && (input_tree->GetEntries() != 0) )
    {
      std::cout << "ERROR: CANNOT ADD BRANCH " << title << " DATA SIZE DID NOT MATCH TREE." 
      	<< "Data size: " << data.size() << ". Tree size: "<< _vnPFOs.size() << std::endl;
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

/*DXG: Angular separation of two particles in the detector*/
double TopTopReconstructionProcessor::dira(ReconstructedParticle *p1, ReconstructedParticle *p2)
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

float TopTopReconstructionProcessor::GetBTag(ReconstructedParticle* jet)
{
  PIDHandler pidh( FlavourJetCol );
  int algo = pidh.getAlgorithmID( "lcfiplus" );
  int ibtag = pidh.getParameterIndex(algo, "BTag");
  const ParticleID &pid = pidh.getParticleID(jet, algo);
  float btag = pid.getParameters()[ibtag];

  return btag;
}

float TopTopReconstructionProcessor::GetCTag(ReconstructedParticle* jet)
{
  PIDHandler pidh( FlavourJetCol );
  int algo = pidh.getAlgorithmID( "lcfiplus" );
  int ictag = pidh.getParameterIndex(algo, "CTag");
  const ParticleID &pid = pidh.getParticleID(jet, algo);
  float ctag = pid.getParameters()[ictag];

  return ctag;
}

bool TopTopReconstructionProcessor::IsCorrectLepton(ReconstructedParticle *lepton)
{
  passedtest =false;
  LCRelationNavigator* relationNavigatorPFOMC = new LCRelationNavigator( LinkerCol );
  EVENT::LCObjectVec relobjMC = relationNavigatorPFOMC->getRelatedToObjects(lepton);
  for(unsigned int i=0;i<relobjMC.size() ; i++) // relobjMC always has size one but just to be safe...
    {
      MCParticle *mcp=NULL;
      mcp=dynamic_cast <MCParticle*>(relobjMC[i]);
      if (mcp!=NULL) //sometimes no link exists-> tends not to happen for particles within the process definition 
	{
	  if((abs(mcp->getPDG())==11 || abs(mcp->getPDG())==13) && mcp->getGeneratorStatus()==1) //if is really a lepton
	    {
	      //here on it becomes complicated to define the true lepton- MC not very consistant. Here have used prefix g to denote generator level particles
	      for (MCParticle* TestMCP = mcp; TestMCP->getParents().size()>0; TestMCP=TestMCP->getParents()[0]) 
		{
		  if(abs(TestMCP->getPDG())==5 || abs(TestMCP->getPDG())==22){break;} // Isolated lepton should never have come from a photon or b jet
		  
		  //Case 1:  Initial Electron + Positron -> gPositron+gElectron -> final electron
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
	}
    }
  delete relationNavigatorPFOMC;
  return passedtest;
}

std::vector<int> TopTopReconstructionProcessor::SortJetVector(std::vector<ReconstructedParticle*> vector)
{
  std::vector<ReconstructedParticle*> dummyvector;
  dummyvector=vector;
  float bestchisquared=1000000000000;
  float chisquared= 1000000000000;
  float WMass =80.385;
  float TopMass =173.21;
  ReconstructedParticle *WBoson=NULL;
  ReconstructedParticle *Top=NULL;
  std::vector<int> resultvector;
  resultvector.resize(4);
  
  //i will be the first W jet index, j the second, and k is the b jet associated with the hadronic W
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	{
	  if(i==j){continue;}
	  WBoson=CombineParticles(dummyvector[i],dummyvector[j]);

	  for(int k=0; k<3; k++)
	    {
	      if(k==i || k==j){continue;}
	      Top=CombineParticles(WBoson,dummyvector[k]);
	      chisquared=(pow(WBoson->getMass()-WMass ,2)/WMass)+(pow(Top->getMass()-TopMass ,2)/TopMass);

	      if(chisquared<bestchisquared)
		{
		  resultvector[0]=i;
		  resultvector[1]=j;
		  resultvector[2]=k;
		  resultvector[3]=6-i-j-k;

		  bestchisquared=chisquared;
		}
	    }
	}
    }
      
  return resultvector;
  
}


