#ifndef SingleTopProcessor_h
#define SingleTopProcessor_h 1

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

using namespace lcio ;
using namespace marlin ;

class SingleTopProcessor : public Processor {

	public:
	
		virtual Processor*  newProcessor() { return new SingleTopProcessor ; }
		
		SingleTopProcessor() ;
		
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

		virtual void end();

		void AddBranch( TTree* input_tree, std::string title, std::vector<float> data );

		double dira(ReconstructedParticle *p1, ReconstructedParticle *p2);
		ReconstructedParticle* CombineParticles( MCParticle *p1, MCParticle *p2 );
		ReconstructedParticle* CombineParticles( ReconstructedParticle *p1, MCParticle *p2 );
		double ReconstructedMass( MCParticle* p1, MCParticle* p2 );
		double ReconstructedMass( ReconstructedParticle* p1, MCParticle* p2 );
		ReconstructedParticle* ConvertToReconstructedParticle(MCParticle *p1 );

		
	protected:

		bool _AllCollectionsExist;

		std::string _MCColName;
		LCCollection *MCCol;

		int _nEvt, _nRun, nMCP ;
		
		TTree* _tRecoMasses;

		std::vector<float> _vTopMass;
		ReconstructedParticle *Top; 
};		
#endif
