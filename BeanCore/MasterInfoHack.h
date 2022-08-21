//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MasterInfoHack                                                       //
//                                                                      //
// The class to get workdir location on master when using PROOF         //
// This is done using the 'fake' Merge method which executes on master  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef _MasterInfoHack_h
#define _MasterInfoHack_h

#include "TNamed.h"
#include "TSystem.h"
#include "TCollection.h"

class MasterInfoHack : public TNamed
{
    public:
        MasterInfoHack();
        ~MasterInfoHack() {};

         Long64_t Merge(TCollection* list);
         std::string  MasterWorkdir() const        {return fDir;}

    private:
        std::string fDir;
        ClassDef(MasterInfoHack,1);
};

#endif
