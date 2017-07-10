#include "HelperClass.h"

#include <iostream>
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <cstring>
//#include <initializer_list>
#include <sstream>
#include <fstream>
#include <utility>

#include <EVENT/LCCollection.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <marlin/Exceptions.h>
#include <UTIL/PIDHandler.h>
#include "UTIL/LCRelationNavigator.h"

#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA

using namespace lcio ;
using namespace marlin ;


HelperClass aHelperClass ;


std::vector<ReconstructedParticle*> HelperClass::SortFatJet(std::vector<ReconstructedParticle*> jets)
{
  //returns a vector of jets where the first two elements are the W jets and the third is the bjet 
  std::vector<ReconstructedParticle*> outjets(3);
  std::vector<ReconstructedParticle*> dummy(3);
  double bestchisquared=10000000000000000;
  float chisquared;
  ReconstructedParticle *WBoson=NULL;
  ReconstructedParticle *Top=NULL;

  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  if(j==i){continue;}
	  int k= 2-i-j;

	  dummy[0]=jets[i];
	  dummy[1]=jets[j];
	  dummy[2]=jets[k];

	  chisquared=TopChi2(dummy);
	  if(chisquared<bestchisquared)
	    {
	      bestchisquared=chisquared;
	      outjets[0]=jets[i];
	      outjets[1]=jets[j];
	      outjets[2]=jets[k];
	    }

	}
    }


  return outjets;

}

float HelperClass::TopChi2(std::vector<ReconstructedParticle*> jets)
{
  //only accepts vectors of three particles
  ReconstructedParticle* WReco= CombineParticles(jets[0],jets[1]);
  ReconstructedParticle* TopReco= CombineParticles(jets[2],WReco);

  double WMass=80.4;
  double WWidth=2.1;
  double TopMass=174;
  double TopWidth=1.4;

  double WChi2= pow((WReco->getMass()-WMass)/WWidth,2);
  double TopChi2= pow((TopReco->getMass()-TopMass)/TopWidth,2);

  return TopChi2+WChi2;
}

ReconstructedParticle* HelperClass::CombineThreeParticles(MCParticle* a, MCParticle* b, MCParticle* c)
{

  ReconstructedParticle* dummy;
  dummy=CombineParticles(a,b);
  dummy=CombineParticles(dummy,c);

  return dummy;
}

ReconstructedParticle* HelperClass::CombineThreeParticles(ReconstructedParticle* a, ReconstructedParticle* b, ReconstructedParticle* c)
{

  ReconstructedParticle* dummy;
  dummy=CombineParticles(a,b);
  dummy=CombineParticles(dummy,c);

  return dummy;
}

float HelperClass::CalculateNSubJettiness(ReconstructedParticle* FatJet,  std::vector<ReconstructedParticle*> subjets, float R0 )
{

  float d0=0; //normalization parameter
  int nSubJets=subjets.size();
  float SumPtDeltaR=0;
  float minDeltaR;
  float deltaR;
  float pT=0;
  // std::cout<<"nSubjets= "<<nSubJets<<std::endl;
  for(unsigned int i=0; i<FatJet->getParticles().size(); i++)
    {
 
      //find minimum seperation between pfo and any subjet
      minDeltaR=9999;
      deltaR=10000;
      for(int j=0; j<nSubJets; j++)
	{
	  deltaR=FindDeltaR(FatJet->getParticles()[i], subjets[j]);
	  //	  std::cout<<"deltaR= "<<deltaR<<std::endl;
	  if(deltaR<minDeltaR)
	    {
	      minDeltaR=deltaR;
	      pT=getTransverseMomentum(FatJet->getParticles()[i]);
	    }
	}
      //    std::cout<<"mindeltaR= "<<minDeltaR<<std::endl;

      

      SumPtDeltaR+=pT*minDeltaR;
      d0+=pT*R0;
    }

  return SumPtDeltaR/d0;
}

float HelperClass::CalculateNSubJettiness(ReconstructedParticle* FatJet, ReconstructedParticle* subjets, float R0 )
{

  float d0=0; //normalization parameter
  float SumPtDeltaR=0;
  float deltaR;
  float pT=0;
  
  for(unsigned int i=0; i<FatJet->getParticles().size(); i++)
    {
 
      deltaR=FindDeltaR(FatJet->getParticles()[i], subjets);
      //  std::cout<<"deltaR= "<<deltaR<<std::endl;
      pT=getTransverseMomentum(FatJet->getParticles()[i]);

      SumPtDeltaR+=pT*deltaR;
      d0+=pT*R0;
    }

  return SumPtDeltaR/d0;
}


float HelperClass::getTransverseMomentum(ReconstructedParticle* p)
{
  return sqrt(p->getMomentum()[0]*p->getMomentum()[0] + p->getMomentum()[1]*p->getMomentum()[1]);

}
float HelperClass::getTotalMomentum(MCParticle* p)
{
  return sqrt(p->getMomentum()[0]*p->getMomentum()[0] + p->getMomentum()[1]*p->getMomentum()[1]  + p->getMomentum()[2]*p->getMomentum()[2]);

}

float HelperClass::getTotalMomentum(ReconstructedParticle* p)
{
  return sqrt(p->getMomentum()[0]*p->getMomentum()[0] + p->getMomentum()[1]*p->getMomentum()[1]  + p->getMomentum()[2]*p->getMomentum()[2]);

}

float HelperClass::getRelativeTransverseMomentum(ReconstructedParticle* p1, ReconstructedParticle* p2)
{
  return sqrt( pow( (p1->getMomentum()[0] - p2->getMomentum()[0]),2) + pow( (p1->getMomentum()[1] - p2->getMomentum()[1]),2));

}

float HelperClass::FindDeltaR(ReconstructedParticle* p1, ReconstructedParticle* p2)
{
  return sqrt( pow(Rapidity(p1)-Rapidity(p2),2) + pow(Phi(p1)-Phi(p2),2)  );
}

float HelperClass::Rapidity(ReconstructedParticle* p1)
{
  //really pseudorapidity eta...
  float pt= sqrt(p1->getMomentum()[0]*p1->getMomentum()[0]+p1->getMomentum()[1]*p1->getMomentum()[1]);
  float pz=p1->getMomentum()[2];
  return asinh(pz/pt);
  
  //float theta =acos(GetCosTheta(p1));
  //return -log(tan(theta/2));

}

float HelperClass::Phi(ReconstructedParticle* p1)
{
  float px= p1->getMomentum()[0];
  float pt= sqrt(p1->getMomentum()[0]*p1->getMomentum()[0]+p1->getMomentum()[1]*p1->getMomentum()[1]);
  return acos(px/pt);

  //return atan(pY/pX);
}


void HelperClass::AddBranch( TTree* input_tree, std::string title, std::vector<float> data )
{
  if( (input_tree->GetEntries() != (int)data.size()) && (input_tree->GetEntries() != 0) )
    {
      std::cout << "ERROR: CANNOT ADD BRANCH " << title << " DATA SIZE DID NOT MATCH TREE." << std::endl;
      return;
    }

  float thisValue=0.;
  TString name=title.c_str();
  name.Append("/F");
  TBranch* output_branch_obj = input_tree->Branch( title.c_str(), &thisValue, name );

  input_tree->SetEntries( (int)data.size() );

  for( unsigned int i=0; i< data.size(); ++i )
    {
      thisValue = data[i];
      output_branch_obj->Fill();
    }
}

double HelperClass::dira(ReconstructedParticle *p1, ReconstructedParticle *p2)
{
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();

  double mag1 = pow(pow(mom1[0],2) + pow(mom1[1],2) + pow(mom1[2],2),0.5);
  double mag2 = pow(pow(mom2[0],2) + pow(mom2[1],2) + pow(mom2[2],2),0.5);

  double theta = acos((mom1[0]*mom2[0]+mom1[1]*mom2[1]+mom1[2]*mom2[2])/(mag1*mag2));

  //if(theta<0) //Takes into account negative angles. Just in case
  //  {
  //    return abs(theta);
  //  }
  //else return theta;

  return theta;
}

double HelperClass::GetCosTheta(ReconstructedParticle *p1)
{

  const double *p =p1->getMomentum();

  double pTot=sqrt( p[0]*p[0] +p[1]*p[1] + p[2]*p[2] );

  return ( p[2]/pTot );
 
}

double HelperClass::GetCosTheta(MCParticle *p1)
{

  const double *p =p1->getMomentum();

  double pTot=sqrt( p[0]*p[0] +p[1]*p[1] + p[2]*p[2] );

  return ( p[2]/pTot );
 
}

double HelperClass::ReconstructedMass( MCParticle* p1, MCParticle* p2 )
{
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  double mass_squared = pow( p1->getEnergy()+p2->getEnergy(), 2 ) 
    - pow( mom1[0]+mom2[0],2 ) 
    - pow( mom1[1]+mom2[1],2 )
    - pow( mom1[2]+mom2[2],2 );
											
  return pow( mass_squared, 0.5 );
}

ReconstructedParticle* HelperClass::CombineParticles( MCParticle *p1, MCParticle *p2 )
{
  // create a ReconstructedParticle that saves the jet
  ReconstructedParticleImpl* newParticle = new ReconstructedParticleImpl();
	
  double Energy = p1->getEnergy() + p2->getEnergy();
  double Charge=p1->getCharge() + p2->getCharge();
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  //sum the momentums together
  double *Momentum = new double[3];
  Momentum[0] = mom1[0]+mom2[0];
  Momentum[1] = mom1[1]+mom2[1];
  Momentum[2] = mom1[2]+mom2[2];
	
  double mass = ReconstructedMass( p1, p2 );
	
  //set the parameters
  newParticle->setMass( mass );
  newParticle->setEnergy( Energy );
  newParticle->setMomentum( Momentum );
  newParticle->setCharge(Charge);

  //now return the jet
  delete Momentum;
  return newParticle ;
}

double HelperClass::ReconstructedMass( ReconstructedParticle* p1, MCParticle* p2 )
{
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  double mass_squared = pow( p1->getEnergy()+p2->getEnergy(), 2 ) 
    - pow( mom1[0]+mom2[0],2 ) 
    - pow( mom1[1]+mom2[1],2 )
    - pow( mom1[2]+mom2[2],2 );
											
  return pow( mass_squared, 0.5 );
}

ReconstructedParticle* HelperClass::CombineParticles( ReconstructedParticle *p1, MCParticle *p2 )
{
  // create a ReconstructedParticle that saves the jet
  ReconstructedParticleImpl* newParticle = new ReconstructedParticleImpl();
	
  double Energy = p1->getEnergy() + p2->getEnergy();
  double Charge= p1->getCharge() + p2->getCharge();

  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  //sum the momentums together
  double *Momentum = new double[3];
  Momentum[0] = mom1[0]+mom2[0];
  Momentum[1] = mom1[1]+mom2[1];
  Momentum[2] = mom1[2]+mom2[2];
	
  double mass = ReconstructedMass( p1, p2 );
	
  //set the parameters
  newParticle->setMass( mass );
  newParticle->setEnergy( Energy );
  newParticle->setMomentum( Momentum );
  newParticle->setCharge(Charge);	

  //now return the jet
  delete Momentum;
  return newParticle ;
}

double HelperClass::ReconstructedMass( ReconstructedParticle* p1, ReconstructedParticle* p2 )
{
  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  double mass_squared = pow( p1->getEnergy()+p2->getEnergy(), 2 ) 
    - pow( mom1[0]+mom2[0],2 ) 
    - pow( mom1[1]+mom2[1],2 )
    - pow( mom1[2]+mom2[2],2 );
											
  return pow( mass_squared, 0.5 );
}

ReconstructedParticle* HelperClass::CombineParticles( ReconstructedParticle *p1, ReconstructedParticle *p2 )
{
  // create a ReconstructedParticle that saves the jet
  ReconstructedParticleImpl* newParticle = new ReconstructedParticleImpl();
	
  double Energy = p1->getEnergy() + p2->getEnergy();
  double Charge=p1->getCharge() + p2->getCharge();

  const double *mom1 = p1->getMomentum();
  const double *mom2 = p2->getMomentum();
	
  //sum the momentums together
  double *Momentum = new double[3];
  Momentum[0] = mom1[0]+mom2[0];
  Momentum[1] = mom1[1]+mom2[1];
  Momentum[2] = mom1[2]+mom2[2];
	
  double mass = ReconstructedMass( p1, p2 );
	
  //set the parameters
  newParticle->setMass( mass );
  newParticle->setEnergy( Energy );
  newParticle->setMomentum( Momentum );
  newParticle->setCharge(Charge);	
	

  //now return the jet
  delete Momentum;
  return newParticle ;
}


bool HelperClass::EnergySort( ReconstructedParticle *p1, ReconstructedParticle *p2 )
{
  return (p1->getEnergy()<p2->getEnergy());
}

ReconstructedParticle* HelperClass::CreateNewParticle( float energy, float px, float py, float pz)
{
  // create a ReconstructedParticle 
  ReconstructedParticleImpl* newParticle = new ReconstructedParticleImpl();

  double *Momentum = new double[3];
  Momentum[0] = px;
  Momentum[1] = py;
  Momentum[2] = pz;

  float mass = pow(energy*energy-px*px-py*py-pz*pz,0.5);
	  
  //set the parameters
  newParticle->setMass( mass );
  newParticle->setEnergy( energy );
  newParticle->setMomentum( Momentum );
	
  return newParticle ;
}


