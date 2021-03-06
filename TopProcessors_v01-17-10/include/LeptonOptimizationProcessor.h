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
#ifndef LeptonOptimizationProcessor_h
#define LeptonOptimizationProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <map>
#include <math.h>
#include <vector>
#include <cstring>

#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/LCObject.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <UTIL/PIDHandler.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TBranch.h>
#include <TVector3.h>
#include <TLorentzVector.h>

using namespace lcio ;
using namespace marlin ;

class LeptonOptimizationProcessor : public Processor {

	public:

		virtual Processor*  newProcessor() { return new LeptonOptimizationProcessor ; }

		LeptonOptimizationProcessor() ;

		virtual void init() ;
		virtual void processRunHeader( LCRunHeader* run ) ;
		virtual void processEvent( LCEvent * evt ) ; 
		virtual void check( LCEvent * evt ) ; 
		virtual void end() ;
		bool IsCorrectLepton(ReconstructedParticle *lepton);
		bool IsCorrectLepton(MCParticle *lepton);
		void AddBranch( TTree* input_tree, std::string title, std::vector<float> data );
		float getConeEnergy(EVENT::ReconstructedParticle*);
	protected:

	
		/** Input collection */
		std::string _inputLepCollection, _PFOCollection;

		/** Output collection (all input with isolated leptons removed) */
		std::string _outputPFOsRemovedIsoLepCollection;

		/** Output collection of isolated leptons */
		std::string _outputIsoLepCollection, _MCColName, _LinkerColName, _jetCollectionName;

		LCCollection *_LepCol, *_PFOCol, *_MCCol, *LinkerCol, *colJet ;
		std::map<ReconstructedParticle*,ReconstructedParticle*> _rpJetMap;

		std::vector<Float_t> _vCorrectLeptonYN, _vPDG, _vEnergy, _vPx, _vPy, _vPz, _vPt, _vPTot, _vCharge, _vJetXt, _vJetZ;
		std::vector<Float_t> _vCorrectPDG, _vCorrectEnergy, _vCorrectPx, _vCorrectPy, _vCorrectPz, _vCorrectPt, _vCorrectPTot, _vCorrectCharge, _vCorrectJetXt, _vCorrectJetZ;
		std::vector<Float_t> _vECALToHCALFraction, _vConeEnergy, _vD0,_vZ0, _vR0, _vdEdX, _vomega, _vphi, _vtrkType, _vHighestEnergyYN, _vCalByP, _vLowestImpactParameterYN;
		  
		int nEvt, nReal, nFake, nDuplicate, nLeptons, njet;
		float Pt, E_TrueLep, E_CurrentRP;
		float ecal, hcal, PTot;
		float d0, z0, r0, dEdX, omega, phi, trkType;
		float _cosConeAngle;
		float jetxt, jetz;
		bool passedtest, foundlep;
		TTree* _tData;
		TTree* _tCorrectData;
		float	bestenergy, lowestimpactparameter;
		int bestenergyposition, bestimpactparameterposition;


		
		ReconstructedParticle* rp, *pfo, *jet;
} ;      

#endif

