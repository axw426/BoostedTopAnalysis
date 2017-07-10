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
#ifndef LeptonOptimizationProcessorTest_h
#define LeptonOptimizationProcessorTest_h 1

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
#include "UTIL/LCRelationNavigator.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TBranch.h>

using namespace lcio ;
using namespace marlin ;

class LeptonOptimizationProcessorTest : public Processor {

	public:

		virtual Processor*  newProcessor() { return new LeptonOptimizationProcessorTest ; }

		LeptonOptimizationProcessorTest() ;

		virtual void init() ;
		virtual void processRunHeader( LCRunHeader* run ) ;
		virtual void processEvent( LCEvent * evt ) ; 
		virtual void check( LCEvent * evt ) ; 
		virtual void end() ;
		bool IsCorrectLepton(ReconstructedParticle *lepton);
		bool IsCorrectLepton(MCParticle *lepton);
		void AddBranch( TTree* input_tree, std::string title, std::vector<float> data );

	protected:

	
		/** Input collection */
		std::string _inputLepCollection, _PFOCollection;

		/** Output collection (all input with isolated leptons removed) */
		std::string _outputPFOsRemovedIsoLepCollection;

		/** Output collection of isolated leptons */
		std::string _outputIsoLepCollection, _MCColName, _LinkerColName, _jetCollectionName;

		LCCollection *_LepCol, *_PFOCol, *_MCCol, *LinkerCol;
		std::map<ReconstructedParticle*,ReconstructedParticle*> _rpJetMap;

		std::vector<Float_t> _vCorrectLeptonYN, _vPDG, _vEnergy, _vPx, _vPy, _vPz, _vPt, _vPTot, _vCharge, _vJetXt, _vJetZ;
		std::vector<Float_t> _vMCDuplicatesPDG, _vNoLinksPDG, _vNoLinksPt, _vNoLinksCosTheta;

		int nEvt, nReal, nFake, nDuplicate, nMultiReco;
		float Pt, Pz, E_TrueLep, E_CurrentRP;
		bool passedtest;
		TTree* _tData;
		TTree* _tCorrectData;

		MCParticle* rp, *pfo, *theMCLepton;

		LCRelationNavigator* relationNavigatorPFOMC;
		EVENT::LCObjectVec relobjMC;
} ;      

#endif

