#include "Riostream.h"

#include "BreitWignerPhi.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "Constants.h"

#include <cmath>

// ClassImp(BreitWignerPhi);

BreitWignerPhi::BreitWignerPhi(const char* name, const char* title,
                               RooAbsReal& _m,
                               RooAbsReal& _mphi,
                               RooAbsReal& _gphi) :
   RooAbsPdf(name,title),
   m("m","m",this,_m),
   mphi("mphi","mphi",this,_mphi),
   gphi("gphi","gphi",this,_gphi) {
}

BreitWignerPhi::BreitWignerPhi(const BreitWignerPhi& other, const char* name) :
   RooAbsPdf(other,name),
   m("m",this,other.m),
   mphi("mphi",this,other.mphi),
   gphi("gphi",this,other.gphi) {
}

//-------------------------------------------------------------------------
Double_t BreitWignerPhi::evaluate() const {
//-------------------------------------------------------------------------
// see BAM-00117:  PhysRevD.91.112001.pdf (arXiv:1504.03194v2)
//
// we assume that Mphi and Gphi are the fit parameters

   if ( m < dL ) { // phase space becomes zero (see p_k)
      return 0;
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

   return (kap * SQ(fm)) / ( SQ(m2-mphi2) + SQ(GM) );
}



