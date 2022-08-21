//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMergeableMap                                                        //
//                                                                      //
// The map (TObjString,TObjString) with merging support                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMergeableMap.h"
#include "TNamed.h"
#include "TObjString.h"

#include <iostream>
using namespace std;
ClassImp(TMergeableMap);

Long64_t TMergeableMap::Merge(TCollection* li)
{
   if (!li) return 0;
   TIter next(li);
   TMap *map;
   while ((map = (TMap*)next())) {
      if (map==this) continue;
      TMapIter next_key(map);
      TObjString* key;

      while (( key = (TObjString*)next_key() )) {

         //check whether this key already exists in map
         TObjString* present_value =
            (TObjString*) GetValue( key->GetName() );
         TObjString* new_value = (TObjString*) map->GetValue( key );
         if (present_value) {
            if ( new_value->GetString() != present_value->GetString() ) {
               Error("TMergeableMaps::Merge",
                     "Trying to merge TMergeableMaps with different"
                     " values(%s,%s) of the same key %s!",
                     new_value->GetString().Data(),
                     present_value->GetString().Data(),
                     key->GetString().Data()                         );
            }
         } else {
            Add(key, new_value);
         }

      }
   }
   return 0;
}
