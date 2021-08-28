#include "Riostream.h"

#include "Interference.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "Constants.h"

#include <cmath>

// ClassImp(Interference);

Interference::Interference(const char* name, const char* title,
                           RooAbsReal& _m,
                           RooAbsReal& _mphi,
                           RooAbsReal& _gphi,
                           RooAbsReal& _A,
                           RooAbsReal& _B,
                           RooAbsReal& _L,
                           RooAbsReal& _U,
                           RooAbsReal& _Ang) :
   RooAbsPdf(name,title),
   m("m","m",this,_m),
   mphi("mphi","mphi",this,_mphi),
   gphi("gphi","gphi",this,_gphi),
   A("A","A",this,_A),
   B("B","B",this,_B),
   L("L","L",this,_L),
   U("U","U",this,_U),
   Ang("Ang","Ang",this,_Ang) {
}

Interference::Interference(const Interference& other, const char* name) :
   RooAbsPdf(other,name),
   m("m",this,other.m),
   mphi("mphi",this,other.mphi),
   gphi("gphi",this,other.gphi),
   A("A",this,other.A),
   B("B",this,other.B),
   L("L",this,other.L),
   U("U",this,other.U),
   Ang("Ang",this,other.Ang) {
}

//-------------------------------------------------------------------------
Double_t Interference::evaluate() const {
//-------------------------------------------------------------------------
// and BreitWignerPhi.cxx and RevArgus.cxx
//
// Return value: simple sum of argus + resonant + interference
//-------------------------------------------------------------------------

   if ( m < dL ) { // phase space becomes zero
      return 0.;
   }

   double m2    = SQ(m);
   double mphi2 = SQ(mphi);
   double gphi2 = SQ(gphi);

   double gam = mphi*sqrt(mphi2+gphi2);
   double kap = 2*M_SQRT2*mphi*gphi*gam / (M_PI*sqrt(mphi2+gam));

   double p_phi0 = sqrt(SQ(Mjpsi2-mphi2-Meta2)-4*mphi2*Meta2) / (2*Mjpsi);
   double p_phi  = sqrt(SQ(Mjpsi2-m2-Meta2)-4*m2*Meta2) / (2*Mjpsi);
   double r_phi  = p_phi / p_phi0;

   double p_k0 = 0.5*sqrt(mphi2-4*Mk2);
   double p_k  = 0.5*sqrt(m2-4*Mk2);
   double r_k  = p_k / p_k0;

   // B(m)/B(m0)
   double BB_phi = sqrt( (1+SQ(R*p_phi0)) / (1+SQ(R*p_phi)) );
   double BB_k = sqrt( (1+SQ(R*p_k0)) / (1+SQ(R*p_k)) );

   double fm = r_phi*r_k*BB_phi*BB_k;

   double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k; // == m*G(m)

   // common multiplier
   double com_mult =  fm / ( SQ(m2-mphi2) + SQ(GM) );

   // Breit-Wigner:
   double BW = kap * fm * com_mult;

   // Argus function:
   double tmp = 1 - (m-L)/U;
   tmp = 1 - tmp*tmp; // = v
   double Ar = m/(L*L) * exp(A*log(tmp) - B*tmp);

   // interference:
   double Comb = 2 * sqrt(Ar*kap) * com_mult *
                 ( (m2-mphi2)*cos(Ang) + GM*sin(Ang) );

   double result = Ar + BW + Comb;

   return result;
}
