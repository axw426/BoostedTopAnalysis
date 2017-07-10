#include "LeptonFinderProcessor.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>

// ----- include for verbosity dependend logging ---------
#include <marlin/VerbosityLevels.h>

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

using namespace lcio ;
using namespace marlin ;

LeptonFinderProcessor aLeptonFinderProcessor ;

LeptonFinderProcessor::LeptonFinderProcessor()
  : Processor("LeptonFinderProcessor") {

  // Processor description
  _description = "Isolated Lepton Finder Processor-> aw jet isolation performs box cuts on min z and xt" ;

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputCollection" ,
			   "Input collection of ReconstructedParticles",
			   _inputPFOsCollection,
			   std::string("PandoraPFANewPFOs"));


  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "OutputCollectionWithoutIsolatedLepton",
			    "Copy of input collection but without the isolated leptons",
			    _outputPFOsRemovedIsoLepCollection,
			    std::string("PandoraPFOsWithoutIsoLep") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "OutputCollectionIsolatedLeptons",
			    "Output collection of isolated leptons",
			    _outputIsoLepCollection,
			    std::string("IsoLep") );


  registerProcessorParameter( "CosConeAngle",
			      "Cosine of the half-angle of the cone used in isolation criteria",
			      _cosConeAngle,
			      float(0.98));

  registerProcessorParameter( "UseCharge",
			      "Use primitive particle ID based on calorimeter energy deposits",
			      _useCharge,
			      bool(true));
  
  registerProcessorParameter( "UsePID",
			      "Use primitive particle ID based on calorimeter energy deposits",
			      _usePID,
			      bool(true));

  registerProcessorParameter( "UsePandoraPID",
			      "Use PID based on Pandora",
			      _usePandoraPID,
			      bool(true));

  registerProcessorParameter( "SortByE",
			      "Use Energy For Sorting Leptons",
			      _sortByE,
			      bool(false));

  registerProcessorParameter( "SortByD0",
			      "Use Energy For Sorting Leptons",
			      _sortByD0,
			      bool(true));

  registerProcessorParameter( "SortByJetZ",
			      "Use Energy For Sorting Leptons",
			      _sortByJetZ,
			      bool(false));

  registerProcessorParameter( "ElectronMinEnergyDepositByMomentum",
			      "Electron ID: Minimum energy deposit divided by momentum",
			      _electronMinEnergyDepositByMomentum,
			      float(0.7));

  registerProcessorParameter( "ElectronMaxEnergyDepositByMomentum",
			      "Electron ID: Maximum energy deposit divided by momentum",
			      _electronMaxEnergyDepositByMomentum,
			      float(1.4));

  registerProcessorParameter( "ElectronMinEcalToHcalFraction",
			      "Electron ID: Minimum Ecal deposit divided by sum of Ecal and Hcal deposits",
			      _electronMinEcalToHcalFraction,
			      float(0.9));

  registerProcessorParameter( "ElectronMaxEcalToHcalFraction",
			      "Electron ID: Maximum Ecal deposit divided by sum of Ecal and Hcal deposits",
			      _electronMaxEcalToHcalFraction,
			      float(1.0));

  registerProcessorParameter( "MuonMinEnergyDepositByMomentum",
			      "Muon ID: Minimum energy deposit divided by momentum",
			      _muonMinEnergyDepositByMomentum,
			      float(0.0));

  registerProcessorParameter( "MuonMaxEnergyDepositByMomentum",
			      "Muon ID: Maximum energy deposit divided by momentum",
			      _muonMaxEnergyDepositByMomentum,
			      float(0.3));

  registerProcessorParameter( "MuonMinEcalToHcalFraction",
			      "Muon ID: Minimum Ecal deposit divided by sum of Ecal and Hcal deposits",
			      _muonMinEcalToHcalFraction,
			      float(0.0));

  registerProcessorParameter( "MuonMaxEcalToHcalFraction",
			      "Muon ID: Maximum Ecal deposit divided by sum of Ecal and Hcal deposits",
			      _muonMaxEcalToHcalFraction,
			      float(0.4));

  registerProcessorParameter( "UseImpactParameter",
			      "Use impact parameter cuts for consistency with primary/secondary track",
			      _useImpactParameter,
			      bool(true));

  registerProcessorParameter( "ImpactParameterMinD0",
			      "Minimum d0 impact parameter",
			      _minD0,
			      float(0.0));

  registerProcessorParameter( "ImpactParameterMaxD0",
			      "Maximum d0 impact parameter",
			      _maxD0,
			      float(1e20));

  registerProcessorParameter( "ImpactParameterMinZ0",
			      "Minimum z0 impact parameter",
			      _minZ0,
			      float(0.0));

  registerProcessorParameter( "ImpactParameterMaxZ0",
			      "Maximum z0 impact parameter",
			      _maxZ0,
			      float(1e20));

  registerProcessorParameter( "ImpactParameterMin3D",
			      "Minimum impact parameter in 3D",
			      _minR0,
			      float(0.0));

  registerProcessorParameter( "ImpactParameterMax3D",
			      "Maximum impact parameter in 3D",
			      _maxR0,
			      float(0.01));

  registerProcessorParameter( "UseImpactParameterSignificance",
			      "Use impact parameter significance cuts for consistency with primary/secondary track",
			      _useImpactParameterSignificance,
			      bool(true));

  registerProcessorParameter( "ImpactParameterMinD0Significance",
			      "Minimum d0 impact parameter significance",
			      _minD0Sig,
			      float(0.0));

  registerProcessorParameter( "ImpactParameterMaxD0Significance",
			      "Maximum d0 impact parameter significance",
			      _maxD0Sig,
			      float(1e20));

  registerProcessorParameter( "ImpactParameterMinZ0Significance",
			      "Minimum z0 impact parameter significance",
			      _minZ0Sig,
			      float(0.0));

  registerProcessorParameter( "ImpactParameterMaxZ0Significance",
			      "Maximum z0 impact parameter significance",
			      _maxZ0Sig,
			      float(1e20));

  registerProcessorParameter( "ImpactParameterMin3DSignificance",
			      "Minimum impact parameter significance in 3D",
			      _minR0Sig,
			      float(0.0));

  registerProcessorParameter( "ImpactParameterMax3DSignificance",
			      "Maximum impact parameter significance in 3D",
			      _maxR0Sig,
			      float(1e20));

  registerProcessorParameter( "UseRectangularIsolation",
			      "Use rectangular cuts on track and cone energy",
			      _useRectangularIsolation,
			      bool(true));

  registerProcessorParameter( "IsolationMinimumTrackEnergy",
			      "Minimum track energy for isolation requirement",
			      _isoMinTrackEnergy,
			      float(15));

  registerProcessorParameter( "IsolationMaximumTrackEnergy",
			      "Maximum track energy for isolation requirement",
			      _isoMaxTrackEnergy,
			      float(1e20));

  registerProcessorParameter( "IsolationMinimumConeEnergy",
			      "Minimum cone energy for isolation requirement",
			      _isoMinConeEnergy,
			      float(0));

  registerProcessorParameter( "IsolationMaximumConeEnergy",
			      "Maximum cone energy for isolation requirement",
			      _isoMaxConeEnergy,
			      float(1e20));

  registerProcessorParameter( "UsePolynomialIsolation",
			      "Use polynomial cuts on track and cone energy",
			      _usePolynomialIsolation,
			      bool(true));

  registerProcessorParameter( "IsolationPolynomialCutA",
			      "Polynomial cut (A) on track energy and cone energy: Econe^2 < A*Etrk^2+B*Etrk+C",
			      _isoPolynomialA,
			      float(0));

  registerProcessorParameter( "IsolationPolynomialCutB",
			      "Polynomial cut (B) on track energy and cone energy: Econe^2 < A*Etrk^2+B*Etrk+C",
			      _isoPolynomialB,
			      float(20));

  registerProcessorParameter( "IsolationPolynomialCutC",
			      "Polynomial cut (C) on track energy and cone energy: Econe^2 < A*Etrk^2+B*Etrk+C",
			      _isoPolynomialC,
			      float(-300));

  registerProcessorParameter( "UseJetIsolation",
			      "Use jet-based isolation",
			      _useJetIsolation,
			      bool(true));

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "JetCollection" ,
			   "Input collection of jets for isolation",
			   _jetCollectionName,
			   std::string("JetsForIsolation"));

  registerProcessorParameter( "JetIsolationVetoMinimumXt",
			      "Minimum Xt in jet-based isolation",
			      _jetIsoVetoMinXt,
			      float(0.));

  registerProcessorParameter( "JetIsolationVetoMaximumXt",
			      "Maximum Xt in jet-based isolation",
			      _jetIsoVetoMaxXt,
			      float(0.25));

  registerProcessorParameter( "JetIsolationVetoMinimumZ",
			      "Mininum Z in jet-based isolation",
			      _jetIsoVetoMinZ,
			      float(0.));

  registerProcessorParameter( "JetIsolationVetoMaximumZ",
			      "Maximum Z in jet-based isolation",
			      _jetIsoVetoMaxZ,
			      float(0.6));

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


void LeptonFinderProcessor::init() { 
  printParameters() ;
  nEvt=0;

  _tRecoData = new TTree("RecoData", "Tree containing final reconstructed masses");


  _vCanbefoundinfirstwave.clear();
  _vRejectedbyfirstwave.clear();
  _vCanbefoundinsecondwave.clear();
  _vRejectedbyfirstwave.clear();
  _vCanbefoundinthirdwave.clear();
  _vSelectedcorrectlepton.clear();
}

void LeptonFinderProcessor::processRunHeader( LCRunHeader* run) { 
} 

void LeptonFinderProcessor::processEvent( LCEvent * evt ) {

  try{
    _pfoCol = evt->getCollection( _inputPFOsCollection ) ;
  }
  catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << _inputPFOsCollection << " collection not available, exiting IsolatedLeptonProcessor..." << std::endl;
      _pfoCol = NULL;
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
 
  // Output PFOs with isolated leptons removed 
  otPFOsRemovedIsoLepCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;
  otPFOsRemovedIsoLepCol->setSubset(true) ;
  
  // Isolated leptons
  otIsoLepCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
  otIsoLepCol->setSubset(true);


  _vSelected.clear();
  _vRejected.clear();
  _vSortedParticles.clear();



  int nparticlesinjets=0;
  ///if (_useJetIsolation) {
  _rpJetMap.clear();
    LCCollection *colJet = evt->getCollection(_jetCollectionName);
    int njet = colJet->getNumberOfElements();
    for (int i=0; i<njet; ++i) {
      ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( colJet->getElementAt(i) );
      nparticlesinjets+=jet->getParticles().size();
      for (ReconstructedParticleVec::const_iterator iter = jet->getParticles().begin();
	   iter != jet->getParticles().end(); ++iter) {
	_rpJetMap.insert( std::make_pair( *iter, jet ) );
      }
    }
    //}

  // look through all PFOs and find lep candidates
  int npfo = _pfoCol->getNumberOfElements();
  
  for (int i = 0; i < npfo; i++ )
    {
      ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(i) );
      if ( IsPandoraLepton( pfo ) && IsIsolatedRectangular( pfo ))
	{
	  _vSelected.push_back( pfo );
	}
      else
	{
	  _vRejected.push_back( pfo );
	}
    }
  
  //pick candidate from all selected pfos
  if(_vSelected.size()>0)
    {
      _vSortedParticles=LeptonFinderProcessor::sortByJetZ(_vSelected);
    }


  else
    {
      _vSelected=_vRejected;
      _vRejected.clear();
      _vSortedParticles=LeptonFinderProcessor::sortByJetZ(_vSelected);

    }
    
      
  //add elements to collections
  
  otIsoLepCol->addElement( _vSortedParticles.back() );
  for(unsigned int i=0; i<(_vSortedParticles.size()-1); i++){otPFOsRemovedIsoLepCol->addElement( _vSortedParticles[i] );}
  for(unsigned int i=0; i<_vRejected.size(); i++){otPFOsRemovedIsoLepCol->addElement( _vRejected[i] );}
    
  
  //add collections to events
  evt->addCollection( otIsoLepCol, _outputIsoLepCollection.c_str() );
  evt->addCollection( otPFOsRemovedIsoLepCol, _outputPFOsRemovedIsoLepCollection.c_str() );

 streamlog_out(DEBUG) << "   nkw processing event: " << evt->getEventNumber() 
		       << "   in run:  " << evt->getRunNumber() 
		       << std::endl ;
  nEvt++;


    
}
void LeptonFinderProcessor::check( LCEvent * evt ) { 
}

void LeptonFinderProcessor::end() {


}






bool LeptonFinderProcessor::IsCharged( ReconstructedParticle* pfo ) {
  if ( pfo->getCharge() == 0 ) return false;
  return true;
}

bool LeptonFinderProcessor::IsPandoraLepton( ReconstructedParticle* pfo ) {
  if ( abs(pfo->getType()) ==11 || abs(pfo->getType()) ==13 ) return true;
  return false;
}

bool LeptonFinderProcessor::IsLepton( ReconstructedParticle* pfo ) {

  float CalE[2];
  getCalEnergy( pfo , CalE );
  double ecale  = CalE[0];
  double hcale  = CalE[1];
  double p      = TVector3( pfo->getMomentum() ).Mag();
  double calByP = p>0 ? (ecale + hcale)/p : 0;
  double calSum = ecale+hcale;
  double ecalFrac = calSum>0 ? ecale / calSum : 0;

  // electron
  if ( calByP >= _electronMinEnergyDepositByMomentum
       && calByP <= _electronMaxEnergyDepositByMomentum
       && ecalFrac >= _electronMinEcalToHcalFraction
       && ecalFrac <= _electronMaxEcalToHcalFraction )
    return true;

  // muon
    if ( calByP >= _muonMinEnergyDepositByMomentum
	 && calByP <= _muonMaxEnergyDepositByMomentum
	 && ecalFrac >= _muonMinEcalToHcalFraction
	 && ecalFrac <= _muonMaxEcalToHcalFraction )
      return true;

    return false;
  }

  bool LeptonFinderProcessor::IsIsolatedLepton( ReconstructedParticle* pfo ) {
    if ( _useCharge && !IsCharged(pfo) )
      return false;

    if ( _usePID && !IsLepton(pfo) )
      return false;
	
    if ( _usePandoraPID && !IsPandoraLepton(pfo) )
      return false;

    if ( _useImpactParameter && !PassesImpactParameterCuts(pfo) )
      return false ;

    if ( _useImpactParameterSignificance && !PassesImpactParameterSignificanceCuts(pfo) )
      return false ;

    if ( _useRectangularIsolation && !IsIsolatedRectangular(pfo) )
      return false;

    if ( _usePolynomialIsolation && !IsIsolatedPolynomial(pfo) )
      return false;

    if ( _useJetIsolation && !IsIsolatedJet(pfo) )
      return false;

    return true;
  }

  bool LeptonFinderProcessor::IsIsolatedRectangular( ReconstructedParticle* pfo ) {
    float E     = pfo->getEnergy() ;
    float coneE = getConeEnergy( pfo );

    if (E < _isoMinTrackEnergy) return false;
    if (E > _isoMaxTrackEnergy) return false;
    if (coneE < _isoMinConeEnergy) return false;
    if (coneE > _isoMaxConeEnergy) return false;

    return true;
  }

  bool LeptonFinderProcessor::IsIsolatedPolynomial( ReconstructedParticle* pfo ) {
    float E     = pfo->getEnergy() ;
    float coneE = getConeEnergy( pfo );

    if ( coneE*coneE <= _isoPolynomialA*E*E + _isoPolynomialB*E + _isoPolynomialC )
      return true ;
    return false;
  }

bool LeptonFinderProcessor::IsIsolatedJet( ReconstructedParticle* pfo ) {
  // jet-based isolated lepton (LAL algorithm)

  if ( _rpJetMap.find( pfo ) == _rpJetMap.end() ) {
    // this is often the case when jet finding fails e.g. due to too few particles in event
    return false;
  }

  ReconstructedParticle* jet = NULL;
  jet = _rpJetMap[pfo];
  
  TVector3 vec1( pfo->getMomentum() );
  TVector3 jetmom( jet->getMomentum() );
  TLorentzVector jetmom4( jet->getMomentum(), jet->getEnergy() );
    
    
  float jetxt = vec1.Pt( jetmom )/jetmom4.M();
  float jetz = pfo->getEnergy()/jet->getEnergy();

  if (jetxt >= _jetIsoVetoMinXt && jetxt < _jetIsoVetoMaxXt
      && jetz >= _jetIsoVetoMinZ && jetz < _jetIsoVetoMaxZ) {
    //printf("xt=%f z=%f (not pass)\n",jetxt,jetz);
    return false;
  }

  //printf("xt=%f z=%f (PASS)\n",jetxt,jetz);
  return true;


}

bool LeptonFinderProcessor::PassesImpactParameterCuts( ReconstructedParticle* pfo ) {
	const EVENT::TrackVec & trkvec = pfo->getTracks();

	if (trkvec.size()==0) return false;

	// TODO: more sophisticated pfo/track matching
	float d0 = fabs(trkvec[0]->getD0());
	float z0 = fabs(trkvec[0]->getZ0());
	float r0 = sqrt( d0*d0 + z0*z0 );

	if ( d0 < _minD0 ) return false;
	if ( d0 > _maxD0 ) return false;
	if ( z0 < _minZ0 ) return false;
	if ( z0 > _maxZ0 ) return false;
	if ( r0 < _minR0 ) return false;
	if ( r0 > _maxR0 ) return false;

	return true;
}

bool LeptonFinderProcessor::PassesImpactParameterSignificanceCuts( ReconstructedParticle* pfo ) {
	const EVENT::TrackVec & trkvec = pfo->getTracks();

	if (trkvec.size()==0) return false;

	// TODO: more sophisticated pfo/track matching
	float d0 = fabs(trkvec[0]->getD0());
	float z0 = fabs(trkvec[0]->getZ0());
	float d0err = sqrt(trkvec[0]->getCovMatrix()[0]);
	float z0err = sqrt(trkvec[0]->getCovMatrix()[9]);

	float d0sig = d0err != 0 ? d0/d0err : 0;
	float z0sig = z0err != 0 ? z0/z0err : 0;
	float r0sig = sqrt( d0sig*d0sig + z0sig*z0sig );

	if ( d0sig < _minD0Sig ) return false;
	if ( d0sig > _maxD0Sig ) return false;
	if ( z0sig < _minZ0Sig ) return false;
	if ( z0sig > _maxZ0Sig ) return false;
	if ( r0sig < _minR0Sig ) return false;
	if ( r0sig > _maxR0Sig ) return false;

	return true;
}

float LeptonFinderProcessor::getConeEnergy( ReconstructedParticle* pfo ) {
	float coneE = 0;

	TVector3 P( pfo->getMomentum() );
	int npfo = _pfoCol->getNumberOfElements();
	for ( int i = 0; i < npfo; i++ ) {
		ReconstructedParticle* pfo_i = dynamic_cast<ReconstructedParticle*>( _pfoCol->getElementAt(i) );

		// don't add itself to the cone energy
		if ( pfo == pfo_i ) continue; 

		TVector3 P_i( pfo_i->getMomentum() );
		float cosTheta = P.Dot( P_i )/(P.Mag()*P_i.Mag());
		if ( cosTheta >= _cosConeAngle )
			coneE += pfo_i->getEnergy(); 
	}

	return coneE;
}

void LeptonFinderProcessor::getCalEnergy( ReconstructedParticle* pfo , float* cale) {
	float ecal = 0;
	float hcal = 0;
	std::vector<lcio::Cluster*> clusters = pfo->getClusters();
	for ( std::vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();
			iCluster!=clusters.end();
			++iCluster) {
		ecal += (*iCluster)->getSubdetectorEnergies()[0];
		hcal += (*iCluster)->getSubdetectorEnergies()[1];
	}
	cale[0] = ecal;
	cale[1] = hcal;
}


float LeptonFinderProcessor::getJetZ(ReconstructedParticle* inputparticle)
{
  ReconstructedParticle* jet=NULL;
  float jetz;
  jet = _rpJetMap[inputparticle];
  if (jet!=NULL) {jetz = inputparticle->getEnergy()/jet->getEnergy();}
  else {jetz=0;}

  return jetz;

}

std::vector<ReconstructedParticle*> LeptonFinderProcessor::sortByJetZ(std::vector<ReconstructedParticle*> inputvector)
{
  int bestJetZPosition =0;
  float bestJetZ=0;
  float newJetZ=0;
  for(int i=0; i<inputvector.size(); i++)
    {
      newJetZ=getJetZ(inputvector[i]);
      if(newJetZ>bestJetZ)
      {
	bestJetZPosition=i;
	bestJetZ=newJetZ;
      }
    }
  ReconstructedParticle* dummy=inputvector[bestJetZPosition];
  inputvector.erase(inputvector.begin()+bestJetZPosition); 
  inputvector.push_back(dummy);

  return inputvector;
}

std::vector<ReconstructedParticle*> LeptonFinderProcessor::sortByE(std::vector<ReconstructedParticle*> inputvector)
{
  int bestEPosition =0;
  float bestE=0;
  float newE=0;
  for(int i=0; i<inputvector.size(); i++)
    {
      newE=inputvector[i]->getEnergy();
      if(newE>bestE)
	{
	  bestEPosition=i;
	  bestE=newE;
	}
    }
  ReconstructedParticle* dummy=inputvector[bestEPosition];
  inputvector.erase(inputvector.begin()+bestEPosition); 
  inputvector.push_back(dummy);

  return inputvector;
}

std::vector<ReconstructedParticle*> LeptonFinderProcessor::sortByD0(std::vector<ReconstructedParticle*> inputvector)
{
  int bestD0Position =0;
  float bestD0=1100000;
  for(int i=0; i<inputvector.size(); i++)
    {
      if(inputvector[i]->getTracks().size()>0 && fabs(inputvector[i]->getTracks()[0]->getD0())<bestD0)
	{
	  bestD0Position=i;
	  bestD0=fabs(inputvector[i]->getTracks()[0]->getD0());
	}
    }
  ReconstructedParticle* dummy=inputvector[bestD0Position];
  inputvector.erase(inputvector.begin()+bestD0Position); 
  inputvector.push_back(dummy);

  return inputvector;
}

std::vector<ReconstructedParticle*> LeptonFinderProcessor::Sort(std::vector<ReconstructedParticle*> inputvector)
{
  std::vector<ReconstructedParticle*> dummyvector;
  if(_sortByE){dummyvector=LeptonFinderProcessor::sortByE(inputvector);}
  else if(_sortByD0){dummyvector=LeptonFinderProcessor::sortByD0(inputvector);}
  else if(_sortByJetZ){dummyvector=LeptonFinderProcessor::sortByJetZ(inputvector);}

  return dummyvector;
}

bool LeptonFinderProcessor::IsCorrectLepton(ReconstructedParticle *lepton)
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

void LeptonFinderProcessor::AddBranch( TTree* input_tree, std::string title, std::vector<float> data )
{
  if( (input_tree->GetEntries() != (int)data.size()) && (input_tree->GetEntries() != 0) )
    {
      std::cout << "ERROR: CANNOT ADD BRANCH " << title << " DATA SIZE DID NOT MATCH TREE." 
      	<< "Data size: " << data.size() << ". Tree size: "<< _vCanbefoundinfirstwave.size() << std::endl;
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
