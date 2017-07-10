#ifndef TopPurityEfficiencyProcessor_h
#define TopPurityEfficiencyProcessor_h 1

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
#include <TGraph.h>
#include <TTree.h>
#include <TBranch.h>


using namespace lcio ;
using namespace marlin ;


/**  Example processor for marlin.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author F. Gaede, DESY
 * @version $Id: TopPurityEfficiencyProcessor.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class TopPurityEfficiencyProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TopPurityEfficiencyProcessor ; }
  
  
  TopPurityEfficiencyProcessor() ;
  
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

  void AddBranch( TTree* input_tree, std::string title, std::vector<float> data );
  float BJetness( ReconstructedParticle* jet);

  
 protected:

  /** Input collection name.
   */
  std::string _LinkerColName; 
  LCCollection* LinkerCol;
  LCCollection* FlavourJets;

  TTree* _tJetTree;
  
  std::string _FlavourJetsName ;
  std::string _ProcessID ;

  std::vector<Float_t> _vEfficiency, _vPurity, _vBTag, _vBJetness, _vHighestBTag;
  std::vector<Float_t> testvector;

  std::vector <float> _vData;
  std::vector< std::vector <float> > _vEventSummary;
  std::vector< std::vector< std::vector <float> > > _vAnalysisSummary;

 
  int _nRun ;
  int _nEvt ;
  float npoints, _npoints;
  float totalPassShouldPass, totalPassShouldFail, totalFailShouldPass, totalFailShouldFail;
  int _whichJet;
  float flavourtype;
  double btag;
  float _minBJetness, bjetness;
  int algo;
  int ibtag;
  ReconstructedParticle* jetparticle;
  ReconstructedParticle *part;
  MCParticle* jetmcp;
  EVENT::LCObjectVec relobjMC;
} ;

#endif



