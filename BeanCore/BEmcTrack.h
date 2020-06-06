#ifndef _BEmcTrack_h
#define _bEmcTrack_h 1

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// BEmcTrack                                                            //
//                                                                      //
// This is adapter between RootEventData/TEmcTrack.h and                //
// DstEvent/DstEmcTrack.h                                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
   typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include "RootEventData/TEmcTrack.h"

class BEmcTrack : public TEmcTrack
{
 public:
        BEmcTrack(const TEmcTrack* trk) : TEmcTrack(*trk) {;}
       ~BEmcTrack() {;}

  // "CLHEP-functions" absent in TEmcTrack
  const HepPoint3D              position()      const;
  const CLHEP::HepSymMatrix     errorMatrix()   const;
};
#endif
