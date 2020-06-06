//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMergeableMap                                                        //
//                                                                      //
// The map with merging support                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef _TMergeableMap_h
#define _TMergeableMap_h

#include "TMap.h"

class TMergeableMap : public TMap
{
   public:
      Long64_t Merge(TCollection* list);

      ClassDef(TMergeableMap,1); // The map with merging support
};
#endif
