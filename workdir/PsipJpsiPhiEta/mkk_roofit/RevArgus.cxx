#include "Riostream.h"

#include "RevArgus.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <cmath>

// ClassImp(RevArgus);

RevArgus::RevArgus(const char* name, const char* title,
                   RooAbsReal& _m,
                   RooAbsReal& _A,
                   RooAbsReal& _B,
                   RooAbsReal& _L,
                   RooAbsReal& _U) :
   RooAbsPdf(name,title),
   m("m","m",this,_m),
   A("A","A",this,_A),
   B("B","B",this,_B),
   L("L","L",this,_L),
   U("U","U",this,_U) {
}


RevArgus::RevArgus(const RevArgus& other, const char* name) :
   RooAbsPdf(other,name),
   m("m",this,other.m),
   A("A",this,other.A),
   B("B",this,other.B),
   L("L",this,other.L),
   U("U",this,other.U) {
}

//----------------------------------------------------------------------
Double_t RevArgus::evaluate() const {
//----------------------------------------------------------------------
//
//   AR(m) = m/L^2 * exp(A*log(v) − B*v), v = 1-(1-(m-L)/U)^2
//
   if ( m < L ) {
      return 0;           // by definition
   }
   double tmp = 1 - (m-L)/U;
   tmp = 1 - tmp*tmp; // = v
   return m/(L*L) * exp(A*log(tmp) - B*tmp);
}
