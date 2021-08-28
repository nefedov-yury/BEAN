#ifndef REVARGUS
#define REVARGUS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RevArgus : public RooAbsPdf {
   public:
      RevArgus() = default;
      RevArgus(const char* name, const char* title,
               RooAbsReal& _m,
               RooAbsReal& _A,
               RooAbsReal& _B,
               RooAbsReal& _L,
               RooAbsReal& _U);
      RevArgus(const RevArgus& other, const char* name=0);
      virtual TObject* clone(const char* newname) const {
         return new RevArgus(*this,newname);
      }
      virtual ~RevArgus() = default;

   protected:
      RooRealProxy m;
      RooRealProxy A;
      RooRealProxy B;
      RooRealProxy L;
      RooRealProxy U;

      Double_t evaluate() const;

//       ClassDef(RevArgus,1)
};
#endif
