#ifndef BeamRemoval_h
#define BeamRemoval_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <map>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include <TH1F.h>

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
 * @version $Id: BeamRemoval.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class BeamRemoval : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new BeamRemoval ; }
  
  
  BeamRemoval() ;
  
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

  ReconstructedParticle* CombineParticles( MCParticle *p1, MCParticle *p2 );
  double ReconstructedMass( MCParticle* p1, MCParticle* p2 );

  
  int _choice;  

 protected:

  /** Input collection name.
   */
  std::string _colName, _outputCol ;

  int _nRun ;
  int _nEvt ;

} ;

#endif



