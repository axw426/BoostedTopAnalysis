#ifndef HelperClass_h
#define HelperClass_h 1

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

class HelperClass  {

 public:
	
  static float TopChi2(std::vector<ReconstructedParticle*> jets);
  static float CalculateNSubJettiness(ReconstructedParticle* FatJet,  std::vector<ReconstructedParticle*> subjets, float R0 );
  static float CalculateNSubJettiness(ReconstructedParticle* FatJet, ReconstructedParticle* subjets, float R0 );
  static float getTransverseMomentum(ReconstructedParticle* p);
  static float getTotalMomentum(MCParticle* p);
  static float getTotalMomentum(ReconstructedParticle* p);
  static float getRelativeTransverseMomentum(ReconstructedParticle* p1, ReconstructedParticle* p2);
  static float FindDeltaR(ReconstructedParticle* p1, ReconstructedParticle* p2);
  static float Rapidity(ReconstructedParticle* p1);
  static float Phi(ReconstructedParticle* p1);
  static void AddBranch( TTree* input_tree, std::string title, std::vector<float> data );
  static double dira(ReconstructedParticle *p1, ReconstructedParticle *p2);
  static double GetCosTheta(ReconstructedParticle *p1);
  static double GetCosTheta(MCParticle *p1);
  static double ReconstructedMass( MCParticle* p1, MCParticle* p2 );
  static ReconstructedParticle* CombineParticles( MCParticle *p1, MCParticle *p2 );
  static double ReconstructedMass( ReconstructedParticle* p1, MCParticle* p2 );
  static ReconstructedParticle* CombineParticles(  ReconstructedParticle*p1, MCParticle *p2 );
  static double ReconstructedMass( ReconstructedParticle* p1, ReconstructedParticle* p2 );
  static ReconstructedParticle* CombineParticles(  ReconstructedParticle*p1, ReconstructedParticle *p2 );
  static bool EnergySort( ReconstructedParticle *p1, ReconstructedParticle *p2 );
  static std::vector<ReconstructedParticle*> SortFatJet(std::vector<ReconstructedParticle*> jets);
  static ReconstructedParticle* CombineThreeParticles(MCParticle* a, MCParticle* b, MCParticle* c);
  static ReconstructedParticle* CreateNewParticle( float energy, float px, float py, float pz);
  static ReconstructedParticle* CombineThreeParticles(ReconstructedParticle* a, ReconstructedParticle* b, ReconstructedParticle* c);

 protected:
  //float d0, pT, minDeltaR, deltaR, SumPtDeltaR;
  //int nSubJets;
  //ReconstructedParticleVec _vFatJetPFOs;
};		
#endif
