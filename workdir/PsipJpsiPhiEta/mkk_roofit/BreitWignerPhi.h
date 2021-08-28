#ifndef BREITWIGNERPHI
#define BREITWIGNERPHI

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class BreitWignerPhi : public RooAbsPdf {
   public:
      BreitWignerPhi() = default;
      BreitWignerPhi(const char* name, const char* title,
                     RooAbsReal& _m,
                     RooAbsReal& _mphi,
                     RooAbsReal& _gphi);
      BreitWignerPhi(const BreitWignerPhi& other, const char* name=0);
      virtual TObject* clone(const char* newname) const {
         return new BreitWignerPhi(*this,newname);
      }
      virtual ~BreitWignerPhi() = default;

   protected:
      RooRealProxy m;
      RooRealProxy mphi;
      RooRealProxy gphi;

      Double_t evaluate() const;

//       ClassDef(BreitWignerPhi,1)
};
#endif
