#ifndef FatJetAnalysisProcessor_h
#define FatJetAnalysisProcessor_h 1

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
#include <IMPL/LCCollectionVec.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TBranch.h>

using namespace lcio ;
using namespace marlin ;

class FatJetAnalysisProcessor : public Processor {

	public:
	
		virtual Processor*  newProcessor() { return new FatJetAnalysisProcessor ; }
		
		FatJetAnalysisProcessor() ;
		
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

		bool IsCorrectLepton(ReconstructedParticle *lepton);
		bool IsCorrectLepton(MCParticle *lepton);

		
	protected:

		bool _AllCollectionsExist, _h2JetExist, _h3JetExist,_h4JetExist, _l2JetExist, _l3JetExist,  _l4JetExist;
		bool _getMCInfo;

		std::string _IsoLepColName, _HadronicFatJetColName, _LeptonicFatJetColName, _Hadronic_2JetColName ,_Hadronic_3JetColName,_Hadronic_4JetColName, _Leptonic_2JetColName ,_Leptonic_3JetColName,_Leptonic_4JetColName, _NLeptonicParticlesColName, _NHadronicParticlesColName, _LinkerColName, _MCColName,  _TightPFOColName, _ee_kt5JetsColName, _FlavourJetsName;
		LCCollection *IsoLepCol, *HadronicFatJetCol, *LeptonicFatJetCol, *Hadronic_2JetCol, *Hadronic_3JetCol, *Hadronic_4JetCol, *Leptonic_2JetCol, *Leptonic_3JetCol, *Leptonic_4JetCol, *NHadronicParticlesJetCol, *NLeptonicParticlesJetCol, *LinkerCol, *MCCol, *TightPFOCol, *ee_kt5JetsCol, *FlavourJets ;
		  
		int _nEvt, _nRun, nLinks ;
		int nLeptonCandidates, nTau;
		
		TTree* _tRecoMasses;
		TTree* _tTest;

		ReconstructedParticle *HadronicJet, *LeptonicJet, *IsoLep, *MCTop, *WBoson, *Neutrino, *LeptonicTop, *TopTest, *MCLeptonicTop, *MCLeptonicTopTest;
		MCParticle *mcIsolatedLepton, *mcNeutrino;
		
		std::vector<float> _vHadronicTopJetMass, _vHadronicEnergy, _vHadronicNParticles, _vHadronicPt;
		std::vector<float> _vLeptonicTopJetMass, _vLeptonicEnergy, _vLeptonicNParticles, _vLeptonicPt;
		std::vector<float> _vJetDIRA;
		std::vector<float> _vTopMass, _vAntiTopMass;
		

		std::vector<ReconstructedParticle*> _vHJet, _vLJet, _vh2SubJets, _vh3SubJets,_vh4SubJets, _vl2SubJets, _vl3SubJets,_vl4SubJets;
		std::vector<float> _vHadronic1SubJettiness, _vHadronic2SubJettiness, _vHadronic3SubJettiness, _vHadronic4SubJettiness;
		std::vector<float> _vLeptonic1SubJettiness,  _vLeptonic2SubJettiness, _vLeptonic3SubJettiness, _vLeptonic4SubJettiness;
		std::vector<float> _vHadronic12Ratio, _vHadronic23Ratio, _vHadronic34Ratio, _vLeptonic12Ratio, _vLeptonic23Ratio,  _vLeptonic34Ratio;
		std::vector<float> _vHighE1SubJet, _vLowE1SubJet, _vLowENParticles, _vHighENParticles, _vHighEtoMidEDira, _vHighEtoLowEDira, _vMidEtoLowEDira;
		std::vector<float> _vIsoLepCosTheta, _vIsoLepCharge, _vIsoLepEnergy, _vIsoLepMomx, _vIsoLepMomy, _vIsoLepMomz, _vIsoLepMomPt, _vIsoLepMomentum, _vIsoLepPID, _vCorrectLeptonYN, _vIsoLepPIDGoodness, _vQByP;
		std::vector<float> _vTopCosTheta, _vTopLepDira, _vMCTopMass, _vMCTopCosTheta, _vMCLeptonCharge, _vMCLeptonID,_vMCLeptonCosTheta, _vMCLeptonMomentum, _vMCQByP, _vMCLeptonicTopAngle, _vMCLeptonicTopMass, _vMCTopSeparation, _vMCTopEnergy;
		std::vector<float> _vChi2, _vTightPFOs;
		std::vector<float> _vWMass, _vWEnergy, _vWPt;
		std::vector<MCParticle*> _vMCJets;
		std::vector<float> _vVisibleEnergy, _vVisiblePt, _vVisiblePz, _vMCEnergy, _vNLeptonCandidates;
		std::vector<float> _vY45, _vY56;
		std::vector<float> _vpthval,  _vmajthval,  _vminthval;
		std::vector<float> _vNTau;
		std::vector<float> _vLeptonicTopMass, _vLeptonicTopEnergy, _vLeptonicTopPt;
		std::vector<float> _vPhotonEnergy, _vPhotonPt, _vPhotonTheta;
		std::vector<float> _vBTags, _vCTags, _vHighestBTag, _vSummedBTags, _vNextHighestBTag;
		  
		float pthval,  majthval,  minthval;
		float SubJettyNormalization;
		float VisibleEnergy, VisiblePt, VisiblePx, VisiblePy, VisiblePz, MCEnergy, MCtopmass;
		float neutrinoEnergy, neutrinoPx, neutrinoPy, neutrinoPz;
		float bestTopDiff;
};		
#endif
