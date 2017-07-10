#ifndef ISRFinder_h
#define ISRFinder_h 1

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

class ISRFinder : public Processor {

	public:
	
		virtual Processor*  newProcessor() { return new ISRFinder ; }
		
		ISRFinder() ;
		
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
		
		string _LoosePFOColName, _TightPFOColName, _LinkerColName, _MCColName;

	protected:
		bool _AllCollectionsExist ; 	

		
		int _nRun ;
		int _nEvt ;
		int _nskipped;

		TTree* _tMCData, *_tRecoData, *_tEventData;

		float ISREnergy, totalMCEnergy;

		std::vector<Float_t> _vPhotonEnergy, _vHasRecoMatch, _vCosTheta;
		std::vector<Float_t> _vISREnergy, _vTotalMCEnergy;
		std::vector<Float_t> _vRecoEnergy;
		
} ;

#endif
