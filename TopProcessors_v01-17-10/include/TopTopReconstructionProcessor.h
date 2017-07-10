#ifndef TopTopReconstructionProcessor_h
#define TopTopReconstructionProcessor_h 1

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

class TopTopReconstructionProcessor : public Processor {

	public:
	
		virtual Processor*  newProcessor() { return new TopTopReconstructionProcessor ; }
		
		TopTopReconstructionProcessor() ;
		
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
		
		double ReconstructedMass( ReconstructedParticle* p1, ReconstructedParticle* p2 );
		ReconstructedParticle* CombineParticles( ReconstructedParticle *p1, ReconstructedParticle *p2 );
		void AddBranch( TTree* input_tree, std::string title, std::vector<float> data );
		double dira(ReconstructedParticle *p1, ReconstructedParticle *p2);
		ReconstructedParticle* ConvertToReconstructedParticle(MCParticle *p1 );
		float GetBTag(ReconstructedParticle* jet);
		float GetCTag(ReconstructedParticle* jet);

		
		//operator overloads
		double ReconstructedMass( ReconstructedParticle* p1, MCParticle* p2 );
		ReconstructedParticle* CombineParticles( ReconstructedParticle *p1, MCParticle *p2 );
		double ReconstructedMass( MCParticle* p1, MCParticle* p2 );
		ReconstructedParticle* CombineParticles( MCParticle *p1, MCParticle *p2 );


		bool IsCorrectLepton(ReconstructedParticle *lepton);
		std::vector<int> SortJetVector(std::vector<ReconstructedParticle*> vector);

  
	protected:
		bool _AllCollectionsExist ; 	

		/** Input collection name.
		 */
		std::string _MCColName, _TwoJetColName, _PFOColName, _PFONoIsoColName, _IsoLepColName, _FourJetColName, _FlavourJetColName, _LinkerColName, _FakeJetColName;
		LCCollection *MCCol, *TwoJetCol, *PFOCol, *PFONoIsoCol, *IsoLepCol, *FourJetCol, *LinkerCol, *FlavourJetCol, *FakeJetCol;

		ReconstructedParticle* Wqq, *isolep, *pseudoHiggs, *JetA, *JetB, *MCWBoson;

		bool _DoFlavourTagging, _UseChiSquaredMatching;
		bool leptonfound, passedtest ;

		int nPFOs, nIsolatedLeptons, nLinks, nMCP, nJets, nFlavourJets ;
		int _nEvt, _nRun ;

		double eventEnergy, pthval, majthval, minthval;
		const double *isolepmom;
		double SigPx,SigPy,SigPz,SigPt,SigE;
		double MissPx,MissPy,MissPz,MissPt,MissE;
		double JetCosTheta, thetaAverage;
		float Y12, Y23, Y34, Y45;
		float nParticlesInJets;
		
		std::vector<ReconstructedParticle*> jets;
		std::vector<Float_t> _vJetMass, _vJetPt,  _vJetDIRA;

		TTree* _tRecoMasses;
		TTree* _tTestTree;

		//vectors to go into TTree
		std::vector<Float_t> _vHadronicTopMass, _vHadronicTopEnergy, _vHadronicTopPx, _vHadronicTopPy, _vHadronicTopPz;
		std::vector<Float_t> _vLeptonicTopMass, _vLeptonicTopEnergy, _vLeptonicTopPx, _vLeptonicTopPy, _vLeptonicTopPz;
		std::vector<Float_t> _vLEPMass, _vLEPEnergy, _vLEPMomx, _vLEPMomy, _vLEPMomz, _vLEPMomt, _vLEPdira, _vLEPCharge, _vLEPID, _vCorrectLeptonYN;
		std::vector<Float_t> _vWqqMass, _vWqqEnergy, _vWqqMomx, _vWqqMomy, _vWqqMomz, _vWqqMomt, _vWqqdira;
		std::vector<Float_t> _vLeptonicWMass, _vLeptonicWEnergy, _vLeptonicWPx, _vLeptonicWPy, _vLeptonicWPz;
		std::vector<Float_t>  _vY12, _vY23, _vY34, _vY45, _vParticlesInJets;
		std::vector<Float_t> _vpthval, _vmajthval, _vminthval, _vpthnisoval, _vmajthnisoval, _vminthnisoval;
		std::vector<Float_t> _vMissPx, _vMissPy, _vMissPz, _vMissPt, _vMissE;
		std::vector<Float_t> _vWqqTopSeperation, _vLoneBJetTopSeperation;
		std::vector<Float_t> _vBJetMass;
		std::vector<Float_t> _vnPFOs;
		std::vector<Float_t> _vBPrimaryTags, _vCPrimaryTags, _vAllBPrimaryTags;
		std::vector<Float_t> _vQJets, _vBJets;
		
		
		
} ;

#endif
