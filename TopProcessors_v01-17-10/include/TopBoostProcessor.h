#ifndef TopBoostProcessor_h
#define TopBoostProcessor_h 1

#include "lcio.h"
#include <string>
#include <math.h>
#include <vector>
#include <cstring>

#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/LCCollectionVec.h"
#include "LCIOSTLTypes.h"


#include <TH1F.h>
#include <TTree.h>
#include <TBranch.h>

using namespace lcio ;
using namespace marlin ;
using namespace std;

class TopBoostProcessor : public Processor {

	public:
	
		virtual Processor*  newProcessor() { return new TopBoostProcessor ; }
		
		TopBoostProcessor() ;
		
		/** Called at the begin of the job before anything is read.
		 * Use to initialize the processor, e.g. book histograms.
		 */
		virtual void init() ;
		
		/** Called for every run.
		 */
		virtual void processRunHeader( LCRunHeader* run ) ;
		
		/** Called for every event - the working horse.
		 */
		virtual void processEvent( LCEvent * evt ) ; 
		
		
		virtual void check( LCEvent * evt ) ; 
		
		
		/** Called after data processing for clean up.
		 */
		virtual void end() ;
		void FindMCTops(std::vector<MCParticle*> _vMCJets, MCParticle* mcNeutrino, MCParticle* mcIsolatedLepton);
		void BoostTops(MCParticle* electron, MCParticle* positron, ReconstructedParticle* top, ReconstructedParticle* antitop);
		ReconstructedParticle* Boost(ReconstructedParticle* Top, ReconstructedParticle* CentreOfMass);

		string _MCParticlesSkimmed;

	protected:
		bool _AllCollectionsExist ; 	

		TTree* _tRecoMasses;
		int _nRun ;
		int _nEvt ;
		int _nskipped;

		MCParticle* MCElectron, *MCPositron, *mcNeutrino, *mcIsolatedLepton;
		ReconstructedParticle* FinalMCLeptonicTop, *FinalMCTop, *BoostedTop, *BoostedAntiTop;
		std::vector<MCParticle*> _vMCJets;
		std::vector<Float_t> _vRawSeparation, _vBoostedSeparation, _vTotalEnergy, _vRawTheta, _vBoostedTheta;

		float TotalEnergy;

} ;

#endif
