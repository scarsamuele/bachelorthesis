#include "TComplex.h"
#include "THStack.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TMath.h"
#include <TRandom3.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>


// Total number of particles:
const Int_t nParticles = 500;
const Int_t nEvents = 1e6;

// Azimuthal angles:
Double_t angles[nParticles]; 

// Particle weights:
Double_t weights[nParticles]; 

// Use or not the non-unit particle weights above:
Bool_t bUseWeights = kFALSE; 

// Book Q-vector components: 
const Int_t maxCorrelator = 8; // We will not go beyond 6-p correlations
const Int_t maxHarmonic = 6; 
const Int_t maxPower = maxCorrelator+1; 
TComplex Qvector[maxHarmonic*maxCorrelator+1][maxCorrelator+1]; // All needed Q-vector components

// Store the final results here:
//  Remark: [2][maxCorrelator] => [Cos,Sin][<2>,<3>,<4>,<5>,<6>]
TProfile *uniform[2][maxCorrelator] = {{NULL}}; // Correlations calculated from Q-vector componets using standalone formulas
TProfile *nonuniform[2][maxCorrelator] = {{NULL}}; // Correlations calculated from Q-vector componets using standalone formulas
TProfile *nonuniformw[2][maxCorrelator] = {{NULL}}; // Correlations calculated from Q-vector componets using standalone formulas

//=======================================================================================================================

void Cosmetics()
{
 // Book everything here.
  
 for(Int_t cs=0;cs<2;cs++) 
 {
  for(Int_t c=0;c<maxCorrelator;c++)
  {
   nonuniformw[cs][c] = new TProfile("","",1,0.,1.);
   nonuniformw[cs][c]->Sumw2();
   nonuniform[cs][c] = new TProfile("","",1,0.,1.);
   nonuniform[cs][c]->Sumw2();
   uniform[cs][c] = new TProfile("","",1,0.,1.);
   uniform[cs][c]->Sumw2();
  } // end of for(Int_t c=0;c<maxCorrelator;c++)
 } // end of for(Int_t cs=0;cs<2;cs++) 

} // void Cosmetics()

//==========================================================================================================================

void CalculateQvectors()
{
 // Calculate Q-vectors.

 // a) Make sure all Q-vectors are initially zero;
 // b) Calculate Q-vectors for available angles and weights. 

 // a) Make sure all Q-vectors are initially zero:
 for(Int_t h=0;h<maxHarmonic*maxCorrelator+1;h++)
 {
  for(Int_t p=0;p<maxCorrelator+1;p++)
  {
   Qvector[h][p] = TComplex(0.,0.);
  } //  for(Int_t p=0;p<maxPower;p++)
 } // for(Int_t h=0;h<maxHarmonic*maxCorrelator+1;h++)

 // b) Calculate Q-vectors for available angles and weights: 
 Double_t dPhi = 0.; // particle angle
 Double_t wPhi = 1.; // particle weight
 Double_t wPhiToPowerP = 1.; // particle weight raised to power p
 for(Int_t i=0;i<nParticles;i++) // loop over particles
 {
  dPhi = angles[i];
  if(bUseWeights){wPhi = weights[i];}
  for(Int_t h=0;h<maxHarmonic*maxCorrelator+1;h++)
  {
   for(Int_t p=0;p<maxCorrelator+1;p++)
   {
    if(bUseWeights){wPhiToPowerP = pow(wPhi,p);}
    Qvector[h][p] += TComplex(wPhiToPowerP*TMath::Cos(h*dPhi),wPhiToPowerP*TMath::Sin(h*dPhi));
   } //  for(Int_t p=0;p<maxPower;p++)
  } // for(Int_t h=0;h<maxHarmonic*MaxCorrelator+1;h++)
 } //  for(Int_t i=0;i<nParticles;i++) // loop over particles

} // void CalculateQvectors()

//==========================================================================================================================

TComplex Q(Int_t n, Int_t p)
{
 // Using the fact that Q{-n,p} = Q{n,p}^*. 
 
 if(n>=0){return Qvector[n][p];} 
 return TComplex::Conjugate(Qvector[-n][p]);
 
} // TComplex Q(Int_t n, Int_t p)

TComplex Two(Int_t n1, Int_t n2)
{
 // Generic two-particle correlation <exp[i(n1*phi1+n2*phi2)]>.

 TComplex two = Q(n1,1)*Q(n2,1)-Q(n1+n2,2);

 return two;

} // TComplex Two(Int_t n1, Int_t n2)

//=======================================================================================================================

TComplex Three(Int_t n1, Int_t n2, Int_t n3)
{
 // Generic three-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3)]>.

TComplex three = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)-Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3); 

 return three;

} // TComplex Three(Int_t n1, Int_t n2, Int_t n3)

//=======================================================================================================================

TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)
{
 // Generic four-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4)]>.

 TComplex four = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
               - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
               + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
               + 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
               + 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);

 return four;

} // TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)

//=======================================================================================================================

TComplex Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5)
{
 // Generic five-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5)]>.

 TComplex five = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
               - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
               + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
               + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
               + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
               - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
               + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
               - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
               + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
               + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
               - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
               + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
               - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
               - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
               + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
               + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
               + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
               + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
               - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
               + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
               + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
               + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
               + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
               - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3) 
               - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
               - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5);
 
 return five;

} // TComplex Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5)

//=======================================================================================================================

TComplex Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6)
{
 // Generic six-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6)]>.

 TComplex six = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
              - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
              + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
              + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
              + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
              - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
              - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
              + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
              + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
              - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
              + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
              - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
              - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
              + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
              + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
              + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
              - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
              + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
              + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
              + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
              - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
              - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
              - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
              - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
              + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
              - 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
              - Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
              - Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
              + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
              - 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
              + 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
              + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
              + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
              + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
              + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
              - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
              - Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
              - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
              + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
              - 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
              + 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
              + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
              + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
              + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
              + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
              - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
              - Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
              - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
              + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
              - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
              + 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
              - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
              - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
              - 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
              - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
              + Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
              - 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
              - Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
              - Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
              + Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
              - 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
              + 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
              - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
              - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
              + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
              - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
              + 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
              + 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
              - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
              - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
              - 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
              - 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
              + 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
              + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
              + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
              + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
              + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
              - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
              - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
              - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
              + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
              - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
              + 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
              - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
              - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
              - 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
              + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
              - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
              + 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
              + 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
              + 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
              + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
              - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
              + 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
              + 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
              + 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
              - 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
              + 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
              - 120.*Q(n1+n2+n3+n4+n5+n6,6);

 return six;

} // TComplex Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6)
//=======================================================================================================================


//=======================================================================================================================
TComplex Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0) 
{
 // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by 
 // Kristjan Gulbrandsen (gulbrand@nbi.dk). 

  Int_t nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip) return c;

  Int_t multp1 = mult+1;
  Int_t nm2 = n-2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n-3;
  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-Double_t(mult)*c2;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0) 

//=======================================================================================================================


TComplex Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7)
{
 // Generic seven-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7)]>.

 Int_t harmonic[7] = {n1,n2,n3,n4,n5,n6,n7};

 TComplex seven = Recursion(7,harmonic); 

 return seven;

} // end of TComplex AliAnalysisTaskMuPa::Seven(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7)

//=======================================================================================================================

TComplex Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8)
{
 // Generic eight-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6+n7*phi7+n8*phi8)]>.

 Int_t harmonic[8] = {n1,n2,n3,n4,n5,n6,n7,n8};

 TComplex eight = Recursion(8,harmonic); 

 return eight;

} // end of TComplex AliAnalysisTaskMuPa::Eight(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6, Int_t n7, Int_t n8)

//=======================================================================================================================

void GF1()
{


delete gRandom;
gRandom = new TRandom3(0);


TF1 * fvarphi = new TF1("fvarphi","(1./TMath::TwoPi())*(1.+2*0.05*cos(x-[0])+2*0.06*cos(2*(x-[0]))+2*0.07*cos(3*(x-[0]))+2*0.08*cos(4*(x-[0]))+2*0.09*cos(5*(x-[0]))+2*0.10*cos(6*(x-[0])))",0.,TMath::TwoPi());

auto gr1= new TCanvas();
THStack *hs = new THStack("hs","");
// distribution non-uniform acceptance
TH1F *hh2= new TH1F("hh2","",100,0.,TMath::TwoPi());
// distribution uniform acceptance
TH1F *hh3= new TH1F("hh3","",100,0.,TMath::TwoPi());
//filling non uniform distribution
Double_t angle=0; //angle to sample
for (int l=0;l<nEvents;l++){ // #events
	fvarphi->SetParameter(0,gRandom->Uniform(TMath::TwoPi())); // sample randomly reaction plane
	for (int m=0;m<nParticles;m++){ // #particles
		angle=fvarphi->GetRandom();
		if (angle > TMath::TwoPi()/6 && angle < TMath::TwoPi()/3){
			if (gRandom -> Uniform(0,1)<0.5){
				hh2->Fill(angle);
			}		
		}	
		else {
			hh2->Fill(angle);
		}
	}
}
hh2->SetMinimum(0);
hh2->SetLineColor(kRed);
hs->Add(hh2);

//filling uniform distribution
for (int l=0;l<nEvents;l++){ // #events
	fvarphi->SetParameter(0,gRandom->Uniform(TMath::TwoPi())); // sample randomly reaction plane
	for (int m=0;m<nParticles;m++){ // #particles
		hh3->Fill(fvarphi->GetRandom());
	}
}
hh3->SetMinimum(0);
hh3->SetLineColor(kBlue);
hs->Add(hh3);
//Drawing and saving
hs->Draw("NOSTACK");
hs->GetXaxis()->SetTitle("#varphi");
hs->GetXaxis()->SetTitleSize(0.05);
hs->GetYaxis()->CenterTitle(true);
hs->GetYaxis()->SetTitle("Counts");
hs->GetYaxis()->SetTitleOffset(0.95);
hs->GetYaxis()->SetTitleSize(0.05);

auto legend = new TLegend();
legend->AddEntry(hh2,"Nonuniform accetance");
legend->AddEntry(hh3,"Uniform acceptance");
legend->Draw("AL*");
gr1->SaveAs("distributions.png");

// obtaining the weights from the histogram
TH1F *w= new TH1F("w","",100,0.,TMath::TwoPi());
// #set the content of the eachbin 1/hh2 
for (Int_t l=1;l<101;l++){
w->SetBinContent(l,(1./(hh2->GetBinContent(l))));
  }
auto gr2= new TCanvas();
w->SetStats(0);
w->GetXaxis()->SetTitle("#varphi");
w->GetXaxis()->SetTitleSize(0.05);
w->GetYaxis()->SetTitle("W_{#varphi}");
w->GetYaxis()->SetTitleOffset(0.8);
w->GetYaxis()->SetTitleSize(0.05);
w->Draw();
gr2->SaveAs("weights.png");

// b) Book all objects:
Cosmetics();

//=======================================================================================================================

//UNIFORM Filling

for (int t=0;t<nEvents;t++){ //events
//Filling angles 
	fvarphi->SetParameter(0,gRandom->Uniform(TMath::TwoPi())); // sample randomly reaction plane
	for(int r=0;r<nParticles;r++){ //particles 
	 // Azimuthal angles:
		angles[r] =fvarphi->GetRandom();
		weights[r]=1.;
	}
	
	bUseWeights = kFALSE; 
	
	 // c) Calculate Q-vectors for available angles and weights;
	CalculateQvectors();

	 // d) Calculate n-particle correlations from Q-vectors (using standalone formulas), and with the defined harmonics "h" from the paper 
	 //  2-p correlations:
  	 TComplex two = Two(-2,2)/Two(0,0).Re();
	 Double_t wTwo = Two(0,0).Re(); // Weight is 'number of combinations' by default
	 uniform[0][0]->Fill(0.5,two.Re(),wTwo); // <<cos(h1*phi1+h2*phi2)>>
	 //  3-p correlations:
	 TComplex three = Three(-5,-1,6)/Three(0,0,0).Re();
	 Double_t wThree = Three(0,0,0).Re(); // Weight is 'number of combinations' by default
	 uniform[0][1]->Fill(0.5,three.Re(),wThree); // <<cos(h1*phi1+h2*phi2+h3*phi3)>>
	 //  4-p correlations:
	 TComplex four = Four(-3,-2,2,3)/Four(0,0,0,0).Re();
	 Double_t wFour = Four(0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 uniform[0][2]->Fill(0.5,four.Re(),wFour); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
	 //  5-p correlations:
	 TComplex five = Five(-5,-4,3,3,3)/Five(0,0,0,0,0).Re();
	 Double_t wFive = Five(0,0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 uniform[0][3]->Fill(0.5,five.Re(),wFive); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>>
	 //  6-p correlations:
	 TComplex six = Six(-2,-2,-1,-1,3,3)/Six(0,0,0,0,0,0).Re();
	 Double_t wSix = Six(0,0,0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 uniform[0][4]->Fill(0.5,six.Re(),wSix); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
	//  7-p correlations:
	 TComplex seven = Seven(-6,-5,-1,1,2,3,6)/Seven(0,0,0,0,0,0,0).Re();
	 Double_t wSeven = Seven(0,0,0,0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 uniform[0][5]->Fill(0.5,seven.Re(),wSeven); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>>
 //  8-p correlations:
	 TComplex eight = Eight(-6,-6,-5,2,3,3,4,5)/Eight(0,0,0,0,0,0,0,0).Re();
	 Double_t wEight = Eight(0,0,0,0,0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 uniform[0][6]->Fill(0.5,eight.Re(),wEight); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
  
}

//=======================================================================================================================

//NON UNIFORM Filling

for (int t=0;t<nEvents;t++){ //events

	fvarphi->SetParameter(0,gRandom->Uniform(TMath::TwoPi())); // sample randomly reaction plane
//Filling angles and weights
	for(int r=0;r<nParticles;r++){ //particles 
		 // Azimuthal angles:
		angle=fvarphi->GetRandom();
		while (angle > TMath::TwoPi()/6 && angle < TMath::TwoPi()/3){
			if (gRandom -> Uniform(0,1)<0.5){
				break;
			}		
			else {
				angle=fvarphi->GetRandom();
			}
		}	
		angles[r]=angle;
		weights[r]=1.;

	}
	
	//To do for each event:

	bUseWeights = kFALSE; 
	
	 // c) Calculate Q-vectors for available angles and weights;
	CalculateQvectors();

	 // d) Calculate n-particle correlations from Q-vectors (using standalone formulas), and with the defined harmonics "h" from the paper 
	 //  2-p correlations:
  	 TComplex two = Two(-2,2)/Two(0,0).Re();
	 Double_t wTwo = Two(0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniform[0][0]->Fill(0.5,two.Re(),wTwo); // <<cos(h1*phi1+h2*phi2)>>
	 //  3-p correlations:
	 TComplex three = Three(-5,-1,6)/Three(0,0,0).Re();
	 Double_t wThree = Three(0,0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniform[0][1]->Fill(0.5,three.Re(),wThree); // <<cos(h1*phi1+h2*phi2+h3*phi3)>>
	 //  4-p correlations:
	 TComplex four = Four(-3,-2,2,3)/Four(0,0,0,0).Re();
	 Double_t wFour = Four(0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniform[0][2]->Fill(0.5,four.Re(),wFour); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
	 //  5-p correlations:
	 TComplex five = Five(-5,-4,3,3,3)/Five(0,0,0,0,0).Re();
	 Double_t wFive = Five(0,0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniform[0][3]->Fill(0.5,five.Re(),wFive); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>>
	 //  6-p correlations:
	 TComplex six = Six(-2,-2,-1,-1,3,3)/Six(0,0,0,0,0,0).Re();
	 Double_t wSix = Six(0,0,0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniform[0][4]->Fill(0.5,six.Re(),wSix); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
	//  7-p correlations:
	 TComplex seven = Seven(-6,-5,-1,1,2,3,6)/Seven(0,0,0,0,0,0,0).Re();
	 Double_t wSeven = Seven(0,0,0,0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniform[0][5]->Fill(0.5,seven.Re(),wSeven); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>>
 //  8-p correlations:
	 TComplex eight = Eight(-6,-6,-5,2,3,3,4,5)/Eight(0,0,0,0,0,0,0,0).Re();
	 Double_t wEight = Eight(0,0,0,0,0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniform[0][6]->Fill(0.5,eight.Re(),wEight); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
  }	

//=======================================================================================================================
//NONUNIFORM + Weights Filling
for (int t=0;t<nEvents;t++){ //events
	fvarphi->SetParameter(0,gRandom->Uniform(TMath::TwoPi())); // sample randomly reaction plane

//Filling angles and weights
	for(int r=0;r<nParticles;r++){ //particles 
		 // Azimuthal angles:
		angle=fvarphi->GetRandom();
		while (angle > TMath::TwoPi()/6 && angle < TMath::TwoPi()/3){
			if (gRandom -> Uniform(0,1)<0.5){
				break;
			}		
			else {
				angle=fvarphi->GetRandom();
			}
		}	
		angles[r]=angle;
		//Find the corresponding weights
		weights[r] = w->GetBinContent(w->FindBin(angles[r]));
	}
	
	//To do for each event:

	bUseWeights = kTRUE; 
	
	 // c) Calculate Q-vectors for available angles and weights;
	CalculateQvectors();

	 // d) Calculate n-particle correlations from Q-vectors (using standalone formulas), and with the defined harmonics "h" from the paper 
	 //  2-p correlations:
  	 TComplex two = Two(-2,2)/Two(0,0).Re();
	 Double_t wTwo = Two(0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniformw[0][0]->Fill(0.5,two.Re(),wTwo); // <<cos(h1*phi1+h2*phi2)>>
	 //  3-p correlations:
	 TComplex three = Three(-5,-1,6)/Three(0,0,0).Re();
	 Double_t wThree = Three(0,0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniformw[0][1]->Fill(0.5,three.Re(),wThree); // <<cos(h1*phi1+h2*phi2+h3*phi3)>>
	 //  4-p correlations:
	 TComplex four = Four(-3,-2,2,3)/Four(0,0,0,0).Re();
	 Double_t wFour = Four(0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniformw[0][2]->Fill(0.5,four.Re(),wFour); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
	 //  5-p correlations:
	 TComplex five = Five(-5,-4,3,3,3)/Five(0,0,0,0,0).Re();
	 Double_t wFive = Five(0,0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniformw[0][3]->Fill(0.5,five.Re(),wFive); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>>
	 //  6-p correlations:
	 TComplex six = Six(-2,-2,-1,-1,3,3)/Six(0,0,0,0,0,0).Re();
	 Double_t wSix = Six(0,0,0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniformw[0][4]->Fill(0.5,six.Re(),wSix); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
	//  7-p correlations:
	 TComplex seven = Seven(-6,-5,-1,1,2,3,6)/Seven(0,0,0,0,0,0,0).Re();
	 Double_t wSeven = Seven(0,0,0,0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniformw[0][5]->Fill(0.5,seven.Re(),wSeven); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>>
 //  8-p correlations:
	 TComplex eight = Eight(-6,-6,-5,2,3,3,4,5)/Eight(0,0,0,0,0,0,0,0).Re();
	 Double_t wEight = Eight(0,0,0,0,0,0,0,0).Re(); // Weight is 'number of combinations' by default
	 nonuniformw[0][6]->Fill(0.5,eight.Re(),wEight); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
  }	


		

//=======================================================================================================================


//Drawing the graph

auto c3 = new TCanvas();

//Adding input values
auto points = new TGraph();
points->SetPoint(0,2,3.6e-1);
points->SetPoint(1,3,4.5e-1);
points->SetPoint(2,4,1.764e-1);
points->SetPoint(3,5,2.4696e-1);
points->SetPoint(4,6,4.41e-2);
points->SetPoint(5,7,9.45e-2);
points->SetPoint(6,8,1.90512e-1);

//graph points uniform
   auto gruniform = new TGraphErrors(); 
//graph points nonuniform
   auto grnonuniform  = new TGraphErrors();
//graph points nonuniform+weights
   auto grnonuniformw = new TGraphErrors();

//Adding points
for (Int_t v=0;v<7;v++){
	//uniform points+errors
	gruniform->SetPoint(v,v+2,uniform[0][v]->GetBinContent(1)/pow(10,-(2+v)));
	gruniform->SetPointError(v,0.,uniform[0][v]->GetBinError(1)/pow(10,-(2+v)));

	//nonuniform points+errors
	grnonuniform->SetPoint(v,v+2,nonuniform[0][v]->GetBinContent(1)/pow(10,-(2+v)));
	grnonuniform->SetPointError(v,0.,nonuniform[0][v]->GetBinError(1)/pow(10,-(2+v)));
	
	//nonuniform+weights points+errors
	 grnonuniformw->SetPoint(v,v+2,nonuniformw[0][v]->GetBinContent(1)/pow(10,-(2+v)));
	grnonuniformw->SetPointError(v,0.,nonuniformw[0][v]->GetBinError(1)/pow(10,-(2+v)));
}

points->SetMarkerColor(kBlack);
points->SetMarkerStyle(kFullDotMedium);
gruniform->SetMarkerColor(kBlack);
gruniform->SetMarkerStyle(kOpenSquare);
grnonuniform->SetMarkerColor(kRed);
grnonuniform->SetLineColor(kRed);
grnonuniform->SetMarkerStyle(kFullSquare);
grnonuniformw->SetMarkerColor(kBlue);
grnonuniformw->SetLineColor(kBlue);
grnonuniformw->SetMarkerStyle(kCircle);

TH1F *h = new TH1F("h","",7,1.5,8.5);
h->GetYaxis()->SetTitle("#LT#LTk#GT#GT/10^{-k}");
h->GetYaxis()->SetTitleOffset(0.8);
h->GetYaxis()->SetTitleSize(0.05);
h->GetYaxis()->SetRangeUser(-0.8,1.1);
h->SetStats(0);
h->GetXaxis()->TAttAxis::SetLabelSize(0.07);
h->GetXaxis()->SetBinLabel(1,"#LT#LT2#GT#GT");
h->GetXaxis()->SetBinLabel(2,"#LT#LT3#GT#GT");
h->GetXaxis()->SetBinLabel(3,"#LT#LT4#GT#GT");
h->GetXaxis()->SetBinLabel(4,"#LT#LT5#GT#GT");
h->GetXaxis()->SetBinLabel(5,"#LT#LT6#GT#GT");
h->GetXaxis()->SetBinLabel(6,"#LT#LT7#GT#GT");
h->GetXaxis()->SetBinLabel(7,"#LT#LT8#GT#GT");
h->Draw();
points->Draw("P SAME");
gruniform->Draw("p  SAME");
grnonuniform->Draw("p  SAME");
grnonuniformw->Draw("p e1 SAME");

auto legend2 = new TLegend(0.115,0.15,0.3,0.32);
legend2->AddEntry(points,"input values","p");
legend2->AddEntry(gruniform,"uniform acceptance","p");
legend2->AddEntry(grnonuniform,"non-uniform acceptance","p");
legend2->AddEntry(grnonuniformw,"non-uniform acceptance + #varphi-weights","p");
legend2->SetFillStyle(0);
legend2->SetTextAlign(12);
legend2->SetTextSize(0.045);
//set cordinates of legend
legend2->SetBorderSize(0.0);
legend2->Draw("p");
c3->SaveAs("framework.png");

}
