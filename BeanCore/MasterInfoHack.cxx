//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MasterInfoHack                                                       //
//                                                                      //
// The class to get workdir location on master when using PROOF         //
// This is done using the 'fake' Merge method which executes on master  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "MasterInfoHack.h"

#include <iostream>
using namespace std;

ClassImp(MasterInfoHack);


MasterInfoHack::MasterInfoHack() 
{
    fDir = "";
    fName = "MasterInfoHack";

    //~ fDir = gSystem->   WorkingDirectory(); 
}

Long64_t MasterInfoHack::Merge(TCollection* list) 
{ 
    fDir = gSystem->   WorkingDirectory(); 
    return 0;
}
