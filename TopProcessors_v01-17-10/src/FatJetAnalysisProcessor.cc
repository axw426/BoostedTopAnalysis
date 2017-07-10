#include "FatJetAnalysisProcessor.h"
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


FatJetAnalysisProcessor aFatJetAnalysisProcessor ;


FatJetAnalysisProcessor::FatJetAnalysisProcessor() : Processor("FatJetAnalysisProcessor") 
{
  // modify processor description
  _description = "DummyProcessor for now, just looks at two fat jets and gets their masses " ;



  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "IsoLep",
			    "IsoLep",
			    _IsoLepColName,
			    std::string("IsoLep") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "ee_kt5Jets",
			   "ee_kt5Jets",
			   _ee_kt5JetsColName,
			   std::string("JetsForIsolation") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "HadronicFatJet",
			    "HadronicFatJet",
			    _HadronicFatJetColName,
			    std::string("HadronicFatJet") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "NHadronicParticles",
			    "NHadronicParticles",
			    _NHadronicParticlesColName,
			    std::string("NHadronicParticles") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "LeptonicFatJet",
			    "LeptonicFatJet",
			    _LeptonicFatJetColName,
			    std::string("LeptonicFatJet") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "NLeptonicParticles",
			    "NLeptonicParticles",
			    _NLeptonicParticlesColName,
			    std::string("NLeptonicParticles") );
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "Hadronic_2Jet",
			   "Hadronic_2Jet",
			   _Hadronic_2JetColName,
			    std::string("Hadronic_2Jet") );
    
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "Hadronic_3Jet",
			    "Hadronic_3Jet",
			    _Hadronic_3JetColName,
			    std::string("Hadronic_3Jet") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "Hadronic_4Jet",
			    "Hadronic_4Jet",
			    _Hadronic_4JetColName,
			    std::string("Hadronic_4Jet") );

  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "Leptonic_2Jet",
			   "Leptonic_2Jet",
			   _Leptonic_2JetColName,
			    std::string("Leptonic_2Jet") );
    
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "Leptonic_3Jet",
			    "Leptonic_3Jet",
			    _Leptonic_3JetColName,
			    std::string("Leptonic_3Jet") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "Leptonic_4Jet",
			    "Leptonic_4Jet",
			    _Leptonic_4JetColName,
			   std::string("Leptonic_4Jet") );

  registerProcessorParameter( "SubJettyNormalization",
			      "Normalization factor R0 for the NSubjettiness calculation",
			      SubJettyNormalization,
			      float(0.3));

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

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "TightPFOs" , 
			   "TightPFOs"  ,
			   _TightPFOColName ,
			   std::string("TightSelectedPandoraPFANewPFOs"));
  

  registerProcessorParameter( "getMCInfo",
			      "getMCInfo",
			      _getMCInfo,
			      bool(false));

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
  			   "FlavourJetName" , 
  			   "Name of the flavour tag jets"  ,
  			   _FlavourJetsName ,
  			   std::string("RefinedJets"));
}

void FatJetAnalysisProcessor::init() 
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
      _tTest= new TTree("Test", "Tree containing final reconstructed masses");
    }

}


void FatJetAnalysisProcessor::processRunHeader( LCRunHeader* run) 
{
  _nRun++ ;
} 


void FatJetAnalysisProcessor::processEvent( LCEvent * evt ) 
{

  
/////////////////////////////////////// Get Collections ////////////////////////////////////////////
  
  _AllCollectionsExist = true ;
  _h2JetExist=true;
  _h3JetExist=true;
  _h4JetExist=true;
  _l2JetExist=true;
  _l3JetExist=true;
  _l4JetExist=true;

  _vHJet.clear();
  _vLJet.clear();
  _vh2SubJets.clear();
  _vh3SubJets.clear();
  _vh4SubJets.clear();
  _vl2SubJets.clear();
  _vl3SubJets.clear();
  _vl4SubJets.clear();

  TightPFOCol = NULL;
  try
    {
      TightPFOCol = evt->getCollection( _TightPFOColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _TightPFOColName << " collection not available" << std::endl;
      TightPFOCol = NULL;
      _AllCollectionsExist = false ;
    }

    ee_kt5JetsCol = NULL;
  try
    {
      ee_kt5JetsCol = evt->getCollection( _ee_kt5JetsColName);
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _ee_kt5JetsColName << " collection not available" << std::endl;
      ee_kt5JetsCol = NULL;
      _AllCollectionsExist = false ;
    }
  
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

  HadronicFatJetCol = NULL;
  try
    {
      HadronicFatJetCol = evt->getCollection( _HadronicFatJetColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _HadronicFatJetColName << " collection not available" << std::endl;
      HadronicFatJetCol=NULL;
      _AllCollectionsExist = false ;
    }


    LeptonicFatJetCol = NULL;
  try
    {
      LeptonicFatJetCol = evt->getCollection( _LeptonicFatJetColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _LeptonicFatJetColName << " collection not available" << std::endl;
      LeptonicFatJetCol=NULL;
      _AllCollectionsExist = false ;
    }

  Hadronic_2JetCol = NULL;
  try
    {
      Hadronic_2JetCol = evt->getCollection( _Hadronic_2JetColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _Hadronic_2JetColName << " collection not available" << std::endl;
      Hadronic_2JetCol=NULL;
      _h2JetExist=false;
    }

  Hadronic_3JetCol = NULL;
  try
    {
      Hadronic_3JetCol = evt->getCollection( _Hadronic_3JetColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _Hadronic_3JetColName << " collection not available" << std::endl;
      Hadronic_3JetCol=NULL;
       _h3JetExist=false;
   }

  Hadronic_4JetCol = NULL;
  try
    {
      Hadronic_4JetCol = evt->getCollection( _Hadronic_4JetColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _Hadronic_4JetColName << " collection not available" << std::endl;
      Hadronic_4JetCol=NULL;
      _h4JetExist=false;
    }
  
  Leptonic_2JetCol = NULL;
  try
    {
      Leptonic_2JetCol = evt->getCollection( _Leptonic_2JetColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _Leptonic_2JetColName << " collection not available" << std::endl;
      Leptonic_2JetCol=NULL;
        _l2JetExist=false;
  }

  Leptonic_3JetCol = NULL;
  try
    {
      Leptonic_3JetCol = evt->getCollection( _Leptonic_3JetColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _Leptonic_3JetColName << " collection not available" << std::endl;
      Leptonic_3JetCol=NULL;
        _l3JetExist=false;
  }

  Leptonic_4JetCol = NULL;
  try
    {
      Leptonic_4JetCol = evt->getCollection( _Leptonic_4JetColName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _Leptonic_4JetColName << " collection not available" << std::endl;
      Leptonic_4JetCol=NULL;
        _l4JetExist=false;
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

  LCCollection* FlavourJets = NULL;

  try
    {
      FlavourJets = evt->getCollection( _FlavourJetsName );
    }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _FlavourJetsName << " collection not available" << std::endl;
      FlavourJets = NULL;
    }


  //////////////////////// MC Info ////////////////////////
  _vMCJets.clear();

  int nMCP = MCCol->getNumberOfElements();
  MCEnergy=0;
  nTau=0;

  for(int i=0; i<nMCP; i++)
    {
      MCParticle* mcp = dynamic_cast<MCParticle*>( MCCol->getElementAt( i ) ) ;
      if (mcp->getGeneratorStatus()==1){MCEnergy+=mcp->getEnergy();}

      if(i>=6 && i<12) // These correspond to the 6 fermion final state
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
	  
	  else if(abs(mcp->getPDG())==15){nTau+=1;}

	}
    }

  _vMCEnergy.push_back(MCEnergy);
  _vNTau.push_back(nTau);
  if(_getMCInfo)
    {
      MCTop=NULL;
      MCLeptonicTop=NULL;
      MCLeptonicTopTest=NULL;
      TopTest=NULL;
      bestTopDiff=10000;
      
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
      
      _vMCTopMass.push_back(MCTop->getMass());
      _vMCTopEnergy.push_back(MCTop->getEnergy());
      _vMCTopCosTheta.push_back(HelperClass::GetCosTheta(MCTop));
      _vMCLeptonCharge.push_back(mcIsolatedLepton->getCharge());
      _vMCLeptonID.push_back(mcIsolatedLepton->getPDG());
      _vMCLeptonCosTheta.push_back(HelperClass::GetCosTheta(mcIsolatedLepton));
      _vMCLeptonMomentum.push_back(HelperClass::getTotalMomentum(mcIsolatedLepton));
      _vMCQByP.push_back(mcIsolatedLepton->getCharge()/HelperClass::getTotalMomentum(mcIsolatedLepton));
      _vMCLeptonicTopMass.push_back(MCLeptonicTop->getMass());
      _vMCLeptonicTopAngle.push_back(HelperClass::GetCosTheta(MCLeptonicTop));
      _vMCTopSeparation.push_back(HelperClass::dira(MCTop,MCLeptonicTop));
    }

  
  ///////////////////////////  Reco Info  ////////////////////////////////////////////////

  LCParameters& pfocolParameters((LCParameters &)TightPFOCol->parameters() );
  pthval = pfocolParameters.getFloatVal( "principleThrustValue" ); 
  majthval = pfocolParameters.getFloatVal( "majorThrustValue" );
  minthval = pfocolParameters.getFloatVal( "minorThrustValue" );
  _vpthval.push_back( pthval );
  _vmajthval.push_back( majthval );
  _vminthval.push_back( minthval );
 
  VisibleEnergy=0;
  VisiblePt=0;
  VisiblePx=0;
  VisiblePy=0;
  VisiblePz=0;
  nLeptonCandidates=0;
  
  _vTightPFOs.push_back(TightPFOCol->getNumberOfElements());
  for(int i=0; i < (TightPFOCol->getNumberOfElements()) ; i++)
    {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( TightPFOCol->getElementAt(i) );
      VisibleEnergy+=pfo->getEnergy();
      VisiblePt+=HelperClass::getTransverseMomentum(pfo);
      VisiblePz+=pfo->getMomentum()[2];
      if((abs(pfo->getType()) ==11 || abs(pfo->getType()) ==13) && pfo->getEnergy()>30){nLeptonCandidates++;}

      //collect info for ISR photon searches
      if(abs(pfo->getType()==22))
	{
	  _vPhotonEnergy.push_back(pfo->getEnergy());
	  _vPhotonPt.push_back(HelperClass::getTransverseMomentum(pfo));
	  _vPhotonTheta.push_back(HelperClass::GetCosTheta(pfo));
	}
    }
  
  _vVisibleEnergy.push_back(VisibleEnergy);
  _vVisiblePt.push_back(VisiblePt);
  _vVisiblePz.push_back(VisiblePz);
  _vNLeptonCandidates.push_back(nLeptonCandidates);


  //reconstruct the neutrino 
  neutrinoEnergy = 1400-VisibleEnergy;
  neutrinoPx =-VisiblePx;
  neutrinoPy =-VisiblePy;
  neutrinoPz =-VisiblePz;
  Neutrino=HelperClass::CreateNewParticle(neutrinoEnergy,neutrinoPx,neutrinoPy,neutrinoPz);

  
  //extract fat jet properties
  HadronicJet = dynamic_cast<ReconstructedParticle*>( HadronicFatJetCol->getElementAt(0) );
  LeptonicJet = dynamic_cast<ReconstructedParticle*>( LeptonicFatJetCol->getElementAt(0) );

  _vHJet.push_back(HadronicJet);
  _vLJet.push_back(LeptonicJet);

		      
  _vHadronic1SubJettiness.push_back(HelperClass::CalculateNSubJettiness(HadronicJet, _vHJet , SubJettyNormalization ));
  _vLeptonic1SubJettiness.push_back(HelperClass::CalculateNSubJettiness(LeptonicJet, _vLJet , SubJettyNormalization ));
  
  _vHadronicTopJetMass.push_back(HadronicJet->getMass());
  _vHadronicEnergy.push_back(HadronicJet->getEnergy());
  //_vHadronicNParticles.push_back(HadronicJet->getParticles().size());
  _vHadronicNParticles.push_back(evt->getCollection( _NHadronicParticlesColName )->getNumberOfElements());
  _vHadronicPt.push_back(HelperClass::getTransverseMomentum(HadronicJet));

  _vLeptonicTopJetMass.push_back(LeptonicJet->getMass());
  _vLeptonicEnergy.push_back(LeptonicJet->getEnergy());
  // _vLeptonicNParticles.push_back(LeptonicJet->getParticles().size());
  _vLeptonicNParticles.push_back(evt->getCollection( _NLeptonicParticlesColName )->getNumberOfElements());
  _vLeptonicPt.push_back(HelperClass::getTransverseMomentum(LeptonicJet));

  _vJetDIRA.push_back(HelperClass::dira(HadronicJet,LeptonicJet));
  _vTopCosTheta.push_back(HelperClass::GetCosTheta(HadronicJet));
  //calculate nsubjettiness

  // hadronic 2 jet
  if(_h2JetExist && Hadronic_2JetCol->getNumberOfElements()>0)
    {
      _vh2SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Hadronic_2JetCol->getElementAt(0) ));
      _vh2SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Hadronic_2JetCol->getElementAt(1) ));
      _vHadronic2SubJettiness.push_back(HelperClass::CalculateNSubJettiness(HadronicJet, _vh2SubJets , SubJettyNormalization ));
    }
  else{_vHadronic2SubJettiness.push_back(-2);}


  //hadronic 3 jets
  if( _h3JetExist && Hadronic_3JetCol->getNumberOfElements()>0 )
    {
      _vh3SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Hadronic_3JetCol->getElementAt(0) ));
      _vh3SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Hadronic_3JetCol->getElementAt(1) ));
      _vh3SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Hadronic_3JetCol->getElementAt(2) ));
      _vHadronic3SubJettiness.push_back(HelperClass::CalculateNSubJettiness(HadronicJet, _vh3SubJets , SubJettyNormalization ));

      //sort subjets by energy
      std::sort(_vh3SubJets.begin(),_vh3SubJets.end(), HelperClass::EnergySort);

      //get subjettiness for single subjets
      _vLowE1SubJet.push_back(HelperClass::CalculateNSubJettiness(_vh3SubJets[0], _vh3SubJets[0] , SubJettyNormalization ));
      _vHighE1SubJet.push_back(HelperClass::CalculateNSubJettiness(_vh3SubJets[2], _vh3SubJets[2] , SubJettyNormalization ));

      //get nParticles in the highest and lowest energy jets 
      _vLowENParticles.push_back(_vh3SubJets[0]->getParticles().size());
      _vHighENParticles.push_back(_vh3SubJets[2]->getParticles().size());
      
      _vHighEtoMidEDira.push_back(HelperClass::dira(_vh3SubJets[2],_vh3SubJets[1]));
      _vHighEtoLowEDira.push_back(HelperClass::dira(_vh3SubJets[2],_vh3SubJets[0]));
      _vMidEtoLowEDira.push_back(HelperClass::dira(_vh3SubJets[1],_vh3SubJets[0]));

      //std::vector<ReconstructedParticle*> _vSortedJets;
      std::vector<ReconstructedParticle*> _vSortedJets=HelperClass::SortFatJet(_vh3SubJets);
      _vChi2.push_back(HelperClass::TopChi2(_vSortedJets));
      WBoson=HelperClass::CombineParticles(_vSortedJets[0],_vSortedJets[1]);
      _vWMass.push_back(WBoson->getMass());
      _vWEnergy.push_back(WBoson->getEnergy());
      _vWPt.push_back(HelperClass::getTransverseMomentum(WBoson));
    }
  else
    {
    _vHadronic3SubJettiness.push_back(-2);
    _vLowE1SubJet.push_back(-2);
    _vHighE1SubJet.push_back(-2);
    _vLowENParticles.push_back(-2);
    _vHighENParticles.push_back(-2);
    _vHighEtoMidEDira.push_back(-2);
    _vHighEtoLowEDira.push_back(-2);
    _vMidEtoLowEDira.push_back(-2);
    _vChi2.push_back(-2);
    _vWMass.push_back(-2);
    _vWEnergy.push_back(-2);
    _vWPt.push_back(-2);
    }

  // hadronic 4 jets
  if( _h4JetExist && Hadronic_4JetCol->getNumberOfElements()>0 )
    {
      _vh4SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Hadronic_4JetCol->getElementAt(0) ));
      _vh4SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Hadronic_4JetCol->getElementAt(1) ));
      _vh4SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Hadronic_4JetCol->getElementAt(2) ));
      _vh4SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Hadronic_4JetCol->getElementAt(3) ));
      _vHadronic4SubJettiness.push_back(HelperClass::CalculateNSubJettiness(HadronicJet, _vh4SubJets , SubJettyNormalization ));
    }
  else{_vHadronic4SubJettiness.push_back(-2);}


  
  //leptonic 2 jets
  if(_l2JetExist && Leptonic_2JetCol->getNumberOfElements()>0 )
    {
      _vl2SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Leptonic_2JetCol->getElementAt(0) ));
      _vl2SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Leptonic_2JetCol->getElementAt(1) ));
      _vLeptonic2SubJettiness.push_back(HelperClass::CalculateNSubJettiness(LeptonicJet, _vl2SubJets , SubJettyNormalization ));
    }
  else{_vLeptonic2SubJettiness.push_back(-2);}

  //leptonic 3 jets
  if( _l3JetExist && Leptonic_3JetCol->getNumberOfElements()>0)
    {
      _vl3SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Leptonic_3JetCol->getElementAt(0) ));
      _vl3SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Leptonic_3JetCol->getElementAt(1) ));
      _vl3SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Leptonic_3JetCol->getElementAt(2) ));
      _vLeptonic3SubJettiness.push_back(HelperClass::CalculateNSubJettiness(LeptonicJet, _vl3SubJets , SubJettyNormalization ));
    }
  else{_vLeptonic3SubJettiness.push_back(-2);}

  //leptonic 4 jets
  if( _l4JetExist && Leptonic_4JetCol->getNumberOfElements()>0)
    {
      _vl4SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Leptonic_4JetCol->getElementAt(0) ));
      _vl4SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Leptonic_4JetCol->getElementAt(1) ));
      _vl4SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Leptonic_4JetCol->getElementAt(2) ));
      _vl4SubJets.push_back(dynamic_cast<ReconstructedParticle*>(Leptonic_4JetCol->getElementAt(3) ));
      _vLeptonic4SubJettiness.push_back(HelperClass::CalculateNSubJettiness(LeptonicJet, _vl4SubJets , SubJettyNormalization ));
      
    }
  else{_vLeptonic4SubJettiness.push_back(-2);}
  

  //calculate subjettiness ratios

  if(Hadronic_2JetCol->getNumberOfElements()>0){_vHadronic12Ratio.push_back(_vHadronic2SubJettiness.back()/_vHadronic1SubJettiness.back());}
  else{_vHadronic12Ratio.push_back(-2);}
    
  if(Hadronic_3JetCol->getNumberOfElements()>0 && Hadronic_2JetCol->getNumberOfElements()>0){ _vHadronic23Ratio.push_back(_vHadronic3SubJettiness.back()/_vHadronic2SubJettiness.back());}
  else{_vHadronic23Ratio.push_back(-2);}

  if(Hadronic_4JetCol->getNumberOfElements()>0 && Hadronic_3JetCol->getNumberOfElements()>0){ _vHadronic34Ratio.push_back(_vHadronic4SubJettiness.back()/_vHadronic3SubJettiness.back());}
  else{_vHadronic34Ratio.push_back(-2);}

  
  if(Leptonic_2JetCol->getNumberOfElements()>0){_vLeptonic12Ratio.push_back(_vLeptonic2SubJettiness.back()/_vLeptonic1SubJettiness.back());}
  else{_vLeptonic12Ratio.push_back(-2);}
    
  if(Leptonic_3JetCol->getNumberOfElements()>0 && Leptonic_2JetCol->getNumberOfElements()>0){ _vLeptonic23Ratio.push_back(_vLeptonic3SubJettiness.back()/_vLeptonic2SubJettiness.back());}
  else{_vLeptonic23Ratio.push_back(-2);}

  if(Leptonic_4JetCol->getNumberOfElements()>0 && Leptonic_3JetCol->getNumberOfElements()>0){ _vLeptonic34Ratio.push_back(_vLeptonic4SubJettiness.back()/_vLeptonic3SubJettiness.back());}
  else{_vLeptonic34Ratio.push_back(-2);}
  
  
  //isolep properties
  IsoLep = dynamic_cast<ReconstructedParticle*>( IsoLepCol->getElementAt(0) );
  _vIsoLepCosTheta.push_back(HelperClass::GetCosTheta(IsoLep));
  _vIsoLepCharge.push_back(IsoLep->getCharge());
  _vIsoLepEnergy.push_back(IsoLep->getEnergy());
  _vIsoLepMomx.push_back(IsoLep->getMomentum()[0]);
  _vIsoLepMomy.push_back(IsoLep->getMomentum()[1]);
  _vIsoLepMomz.push_back(IsoLep->getMomentum()[2]);
  _vIsoLepMomPt.push_back(HelperClass::getTransverseMomentum(IsoLep));
  _vIsoLepMomentum.push_back(HelperClass::getTotalMomentum(IsoLep));
  _vQByP.push_back(IsoLep->getCharge()/HelperClass::getTotalMomentum(IsoLep));
  _vIsoLepPID.push_back(IsoLep->getType());
  _vIsoLepPIDGoodness.push_back(IsoLep->getGoodnessOfPID());
  _vTopLepDira.push_back(HelperClass::dira(IsoLep,HadronicJet));
  _vCorrectLeptonYN.push_back(IsCorrectLepton(IsoLep));

  if(IsoLep->getCharge()==-1)
    {
      _vTopMass.push_back(HadronicJet->getMass());
      _vAntiTopMass.push_back(-2);
    }
  else if(IsoLep->getCharge()==1)
    {
      _vAntiTopMass.push_back(HadronicJet->getMass());
      _vTopMass.push_back(-2);
    }
  else
    {
      _vTopMass.push_back(-2);
      _vAntiTopMass.push_back(-2);
    }

  //reconstruct the leptonic top
  LeptonicTop=HelperClass::CombineThreeParticles(Neutrino,IsoLep,LeptonicJet);
  _vLeptonicTopMass.push_back(LeptonicTop->getMass());
  _vLeptonicTopEnergy.push_back(LeptonicTop->getEnergy());
  _vLeptonicTopPt.push_back(HelperClass::getTransverseMomentum(LeptonicTop));

  
  LCParameters& col5jParameters((LCParameters &)ee_kt5JetsCol->parameters() );
  float Y45 = col5jParameters.getFloatVal( "y_{n-1,n}" );
  float Y56 = col5jParameters.getFloatVal( "y_{n,n+1}" );
  _vY45.push_back(-(std::log(Y45)));
  _vY56.push_back(-(std::log(Y56)));


  //////Flavour Tagging/////////////////////////////////
  
  // get PIDHandler associated with the jet collection
  PIDHandler pidh( FlavourJets );
  // get algorithm ID associated with LCFIPlus
  int algo = pidh.getAlgorithmID( "lcfiplus" );
  // get index number for flavor tagging
  int ibtag = pidh.getParameterIndex(algo, "BTag");
  int ictag = pidh.getParameterIndex(algo, "CTag");
  _vBTags.clear();
  for(int i=0; i < FlavourJets->getNumberOfElements(); i++) 
    {
      ReconstructedParticle *part = dynamic_cast<ReconstructedParticle*>( FlavourJets->getElementAt( i ) );
      const ParticleID &pid = pidh.getParticleID(part, algo);
     _vBTags.push_back(pid.getParameters()[ibtag]);
     //_vCTags.push_back(pid.getParameters()[ictag]);
    }

  std::sort(_vBTags.begin(), _vBTags.end());
  _vHighestBTag.push_back(_vBTags[3]);
  _vNextHighestBTag.push_back(_vBTags[2]);
  _vSummedBTags.push_back(_vBTags[3]+_vBTags[2]);

  
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		       << "   in run:  " << evt->getRunNumber() << std::endl ;
  
  _nEvt ++ ;
      
  //std::cout << "Events processed: " << _nEvt << std::endl;
  
}

void FatJetAnalysisProcessor::check( LCEvent * evt ) 
{ 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void FatJetAnalysisProcessor::end()
{ 
    
  // In here we fill the trees with the information stored within the vectors
  // Note: All of the vectors must be the same size else the branch will not be added

  HelperClass::AddBranch(_tRecoMasses, "MCEnergy", _vMCEnergy);
  HelperClass::AddBranch(_tRecoMasses, "VisibleEnergy", _vVisibleEnergy);
  HelperClass::AddBranch(_tRecoMasses, "VisiblePt", _vVisiblePt);
  HelperClass::AddBranch(_tRecoMasses, "VisiblePz", _vVisiblePz);
  HelperClass::AddBranch(_tRecoMasses, "WMass", _vWMass);
  HelperClass::AddBranch(_tRecoMasses, "WEnergy", _vWEnergy);
  HelperClass::AddBranch(_tRecoMasses, "WPt", _vWPt);
  HelperClass::AddBranch(_tRecoMasses, "HadronicJetMass", _vHadronicTopJetMass);
  HelperClass::AddBranch(_tRecoMasses, "HadronicEnergy", _vHadronicEnergy);
  HelperClass::AddBranch(_tRecoMasses, "HadronicNParticles", _vHadronicNParticles);
  HelperClass::AddBranch(_tRecoMasses, "HadronicPt", _vHadronicPt);
  HelperClass::AddBranch(_tRecoMasses, "LeptonicJetMass", _vLeptonicTopJetMass);
  HelperClass::AddBranch(_tRecoMasses, "LeptonicEnergy", _vLeptonicEnergy);
  HelperClass::AddBranch(_tRecoMasses, "LeptonicNParticles", _vLeptonicNParticles);
  HelperClass::AddBranch(_tRecoMasses, "LeptonicPt", _vLeptonicPt);
  HelperClass::AddBranch(_tRecoMasses, "JetDira", _vJetDIRA);
  HelperClass::AddBranch(_tRecoMasses, "Hadronic1SubJettiness", _vHadronic1SubJettiness);
  HelperClass::AddBranch(_tRecoMasses, "Hadronic2SubJettiness", _vHadronic2SubJettiness);
  HelperClass::AddBranch(_tRecoMasses, "Hadronic3SubJettiness", _vHadronic3SubJettiness);
  HelperClass::AddBranch(_tRecoMasses, "Hadronic4SubJettiness", _vHadronic4SubJettiness);
  HelperClass::AddBranch(_tRecoMasses, "Leptonic1SubJettiness", _vLeptonic1SubJettiness);
  HelperClass::AddBranch(_tRecoMasses, "Leptonic2SubJettiness", _vLeptonic2SubJettiness);
  HelperClass::AddBranch(_tRecoMasses, "Leptonic3SubJettiness", _vLeptonic3SubJettiness);
  HelperClass::AddBranch(_tRecoMasses, "Leptonic4SubJettiness", _vLeptonic4SubJettiness);
  HelperClass::AddBranch(_tRecoMasses, "Hadronic12Ratio", _vHadronic12Ratio);
  HelperClass::AddBranch(_tRecoMasses, "Hadronic23Ratio", _vHadronic23Ratio);
  HelperClass::AddBranch(_tRecoMasses, "Hadronic34Ratio", _vHadronic34Ratio);
  HelperClass::AddBranch(_tRecoMasses, "Leptonic12Ratio", _vLeptonic12Ratio);
  HelperClass::AddBranch(_tRecoMasses, "Leptonic23Ratio", _vLeptonic23Ratio);
  HelperClass::AddBranch(_tRecoMasses, "Leptonic34Ratio", _vLeptonic34Ratio);
  HelperClass::AddBranch(_tRecoMasses, "HighE1SubJet", _vHighE1SubJet);
  HelperClass::AddBranch(_tRecoMasses, "LowE1SubJet", _vLowE1SubJet);
  HelperClass::AddBranch(_tRecoMasses, "LowENParticles", _vLowENParticles);
  HelperClass::AddBranch(_tRecoMasses, "HighENParticles", _vHighENParticles);
  HelperClass::AddBranch(_tRecoMasses, "HighEtoLowEDira", _vHighEtoLowEDira);
  HelperClass::AddBranch(_tRecoMasses, "HighEtoMidEDira", _vHighEtoMidEDira);
  HelperClass::AddBranch(_tRecoMasses, "MidEtoLowEDira", _vMidEtoLowEDira);
  HelperClass::AddBranch(_tRecoMasses, "IsoLepCosTheta", _vIsoLepCosTheta);
  HelperClass::AddBranch(_tRecoMasses, "TopCosTheta", _vTopCosTheta);
  HelperClass::AddBranch(_tRecoMasses, "IsoLepCharge", _vIsoLepCharge);
  HelperClass::AddBranch(_tRecoMasses, "IsoLepEnergy", _vIsoLepEnergy);
  HelperClass::AddBranch(_tRecoMasses, "IsoLepMomx", _vIsoLepMomx);
  HelperClass::AddBranch(_tRecoMasses, "IsoLepMomy", _vIsoLepMomy);
  HelperClass::AddBranch(_tRecoMasses, "IsoLepMomz", _vIsoLepMomz);
  HelperClass::AddBranch(_tRecoMasses, "IsoLepMomPt", _vIsoLepMomPt);
  HelperClass::AddBranch(_tRecoMasses, "IsoLepMomentum", _vIsoLepMomentum);
  HelperClass::AddBranch(_tRecoMasses, "QByP", _vQByP);
  HelperClass::AddBranch(_tRecoMasses, "IsoLepPID", _vIsoLepPID);
  HelperClass::AddBranch(_tRecoMasses, "IsoLepPIDGoodness", _vIsoLepPIDGoodness);
  HelperClass::AddBranch(_tRecoMasses, "NLeptonCandidates", _vNLeptonCandidates);
  HelperClass::AddBranch(_tRecoMasses, "CorrectLeptonYN", _vCorrectLeptonYN);
  HelperClass::AddBranch(_tRecoMasses, "TopLepDira", _vTopLepDira);
  HelperClass::AddBranch(_tRecoMasses, "Chi2", _vChi2);
  HelperClass::AddBranch(_tRecoMasses, "TightPFOs", _vTightPFOs);
  HelperClass::AddBranch(_tRecoMasses, "Y45", _vY45);
  HelperClass::AddBranch(_tRecoMasses, "Y56", _vY56);
  HelperClass::AddBranch(_tRecoMasses, "pthval", _vpthval);
  HelperClass::AddBranch(_tRecoMasses, "minthval", _vminthval);
  HelperClass::AddBranch(_tRecoMasses, "majthval", _vmajthval);
  HelperClass::AddBranch(_tRecoMasses, "LeptonicTopMass", _vLeptonicTopMass);
  HelperClass::AddBranch(_tRecoMasses, "LeptonicTopEnergy", _vLeptonicTopEnergy);
  HelperClass::AddBranch(_tRecoMasses, "LeptonicTopPt", _vLeptonicTopPt);
  HelperClass::AddBranch(_tRecoMasses, "HighestBTag", _vHighestBTag);
  HelperClass::AddBranch(_tRecoMasses, "NextHighestBTag", _vNextHighestBTag);
  HelperClass::AddBranch(_tRecoMasses, "SummedBTags", _vSummedBTags);



  //HelperClass::AddBranch(_tRecoMasses, "", _);

  if(_getMCInfo)
    {
      HelperClass::AddBranch(_tRecoMasses, "MCTopMass", _vMCTopMass);
      HelperClass::AddBranch(_tRecoMasses, "MCTopCosTheta", _vMCTopCosTheta);
      HelperClass::AddBranch(_tRecoMasses, "MCTopEnergy", _vMCTopEnergy);
      HelperClass::AddBranch(_tRecoMasses, "MCLeptonCharge", _vMCLeptonCharge);
      HelperClass::AddBranch(_tRecoMasses, "MCLeptonCosTheta", _vMCLeptonCosTheta);
      HelperClass::AddBranch(_tRecoMasses, "MCLeptonID", _vMCLeptonID);
      HelperClass::AddBranch(_tRecoMasses, "MCLeptonMomentum", _vMCLeptonMomentum);
      HelperClass::AddBranch(_tRecoMasses, "MCQByP", _vMCQByP);
      HelperClass::AddBranch(_tRecoMasses, "MCLeptonicTopAngle", _vMCLeptonicTopAngle);
      HelperClass::AddBranch(_tRecoMasses, "MCLeptonicTopMass", _vMCLeptonicTopMass);
      HelperClass::AddBranch(_tRecoMasses, "MCTopSeparation", _vMCTopSeparation);
      HelperClass::AddBranch(_tRecoMasses, "NTau", _vNTau);

    }

    HelperClass::AddBranch(_tTest, "PhotonEnergy", _vPhotonEnergy);
    HelperClass::AddBranch(_tTest, "PhotonPt", _vPhotonPt);
    HelperClass::AddBranch(_tTest, "PhotonTheta", _vPhotonTheta);
}

bool FatJetAnalysisProcessor::IsCorrectLepton(ReconstructedParticle *lepton)
{
  bool passedtest =false;
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

bool FatJetAnalysisProcessor::IsCorrectLepton(MCParticle *lepton)
{
  bool passedtest =false;

  MCParticle *mcp=lepton;

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
    
    
  return passedtest;
}





