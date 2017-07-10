/**
 * @brief Marlin processor for finding isolated leptons.
 * @author Ryo Yonamine <yonamine@post.kek.jp>
 * @author Tomohiko Tanabe <tomohiko@icepp.s.u-tokyo.ac.jp>
 * @version $Id:$
 *
 * Given a list of ReconstructedParticle, identify isolated leptons
 * based on the track cone energy, lepton identification,
 * and the track impact parameters (optional).
 */
#ifndef LeptonFinderProcessor_h
#define LeptonFinderProcessor_h 1

#include <string>
#include <map>

#include <marlin/Processor.h>
#include <lcio.h>

#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <TTree.h>

using namespace lcio ;
using namespace marlin ;

class LeptonFinderProcessor : public Processor {

 public:

  virtual Processor*  newProcessor() { return new LeptonFinderProcessor ; }

  LeptonFinderProcessor() ;

  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  bool IsIsolatedLepton( ReconstructedParticle* pfo ) ;
  bool IsIsolatedRectangular( ReconstructedParticle* pfo ) ;
  bool IsIsolatedPolynomial( ReconstructedParticle* pfo ) ;
  bool IsIsolatedJet( ReconstructedParticle* pfo ) ;
  bool IsCharged( ReconstructedParticle* pfo ) ;
  bool IsPandoraLepton( ReconstructedParticle* pfo );
  bool IsLepton( ReconstructedParticle* pfo ) ;
  bool PassesImpactParameterCuts( ReconstructedParticle* pfo ) ; 
  bool PassesImpactParameterSignificanceCuts( ReconstructedParticle* pfo ) ; 
  float getConeEnergy( ReconstructedParticle* pfo ) ;
  void getCalEnergy( ReconstructedParticle* pfo , float* cale) ;
  float getJetZ(ReconstructedParticle* inputparticle);
  std::vector<ReconstructedParticle*> sortByJetZ(std::vector<ReconstructedParticle*> inputvector);
  std::vector<ReconstructedParticle*> sortByE(std::vector<ReconstructedParticle*> inputvector);
  std::vector<ReconstructedParticle*> sortByD0(std::vector<ReconstructedParticle*> inputvector);
  std::vector<ReconstructedParticle*> Sort(std::vector<ReconstructedParticle*> inputvector);
  bool IsCorrectLepton(ReconstructedParticle *lepton);
  void AddBranch( TTree* input_tree, std::string title, std::vector<float> data );

 protected:

	
  /** Input collection */
  std::string _inputPFOsCollection, _MCColName, _LinkerColName;

  /** Output collection (all input with isolated leptons removed) */
  std::string _outputPFOsRemovedIsoLepCollection;

  /** Output collection of isolated leptons */
  std::string _outputIsoLepCollection;

  TTree* _tRecoData;
  
  LCCollection* _pfoCol, *_InLepCol, *_InPFONoLepCol,*_MCCol, *LinkerCol;
  std::vector<ReconstructedParticle*> _vSelected, _vRejected, _vSortedParticles;

  int nEvt;
  float Pt;


  double BestImpact;
  int PositionOfBestMatch;
	
  float _cosConeAngle;
  bool _useCharge, passedtest;
  /** If set to true, uses PID cuts */
  bool _usePID;
  bool _usePandoraPID, _sortByE, _sortByD0, _sortByJetZ;
  float _electronMinEnergyDepositByMomentum;
  float _electronMaxEnergyDepositByMomentum;
  float _electronMinEcalToHcalFraction;
  float _electronMaxEcalToHcalFraction;
  float _muonMinEnergyDepositByMomentum;
  float _muonMaxEnergyDepositByMomentum;
  float _muonMinEcalToHcalFraction;
  float _muonMaxEcalToHcalFraction;

  /** If set to true, uses impact parameter cuts */
  bool _useImpactParameter;
  float _minD0;
  float _maxD0;
  float _minZ0;
  float _maxZ0;
  float _minR0;
  float _maxR0;

  /** If set to true, uses impact parameter significance cuts */
  bool _useImpactParameterSignificance;
  float _minD0Sig;
  float _maxD0Sig;
  float _minZ0Sig;
  float _maxZ0Sig;
  float _minR0Sig;
  float _maxR0Sig;

  /** If set to true, uses rectangular cuts for isolation */
  bool _useRectangularIsolation;
  float _isoMinTrackEnergy;
  float _isoMaxTrackEnergy;
  float _isoMinConeEnergy;
  float _isoMaxConeEnergy;

  /** If set to true, uses polynomial cuts for isolation */
  bool _usePolynomialIsolation;
  float _isoPolynomialA;
  float _isoPolynomialB;
  float _isoPolynomialC;

  /** If set to true, uses jet-based isolation (LAL algorithm) */
  bool _useJetIsolation;
  std::string _jetCollectionName;
  std::map<ReconstructedParticle*,ReconstructedParticle*> _rpJetMap;
  float _jetIsoVetoMinXt;
  float _jetIsoVetoMaxXt;
  float _jetIsoVetoMinZ;
  float _jetIsoVetoMaxZ;

  LCCollectionVec* otPFOsRemovedIsoLepCol;
  LCCollectionVec* otIsoLepCol;

  float canbefoundinfirstwave, canbefoundinsecondwave, canbefoundinthirdwave, selectedcorrectlepton,rejectedbyfirstwave,rejectedbysecondwave;
  std::vector<float> _vCanbefoundinfirstwave, _vRejectedbyfirstwave, _vCanbefoundinsecondwave,_vRejectedbysecondwave, _vCanbefoundinthirdwave, _vSelectedcorrectlepton;
} ;

#endif

