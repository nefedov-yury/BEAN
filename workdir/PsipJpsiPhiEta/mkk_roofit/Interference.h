#ifndef INTERFERENCE
#define INTERFERENCE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class Interference : public RooAbsPdf {
   public:
      Interference() = default;
      Interference(const char* name, const char* title,
                   RooAbsReal& _m,
                   RooAbsReal& _mphi,
                   RooAbsReal& _gphi,
                   RooAbsReal& _A,
                   RooAbsReal& _B,
                   RooAbsReal& _L,
                   RooAbsReal& _U,
                   RooAbsReal& _Ang);
      Interference(const Interference& other, const char* name=0);
      virtual TObject* clone(const char* newname) const {
         return new Interference(*this,newname);
      }
      virtual ~Interference() = default;

   protected:
      RooRealProxy m;
      RooRealProxy mphi;
      RooRealProxy gphi;
      RooRealProxy A;
      RooRealProxy B;
      RooRealProxy L;
      RooRealProxy U;
      RooRealProxy Ang;

      Double_t evaluate() const;

//   ClassDef(Interference,1)
};
#endif
