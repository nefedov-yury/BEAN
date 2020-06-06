#ifndef _BExtTrack_h
#define _BExtTrack_h

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// BExtTrack                                                            //
//                                                                      //
// This is adapter between RootEventData/TExtTrack.h and                //
// DstEvent/DstExtTrack.h                                               //
//                                                                      //
// Root file contains information only for extrapolation of pions,      //
// therefore all functions with (int parID) are refused                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <string>

#include <CLHEP/Matrix/SymMatrix.h>
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Vector/ThreeVector.h>
using CLHEP::HepSymMatrix;
using CLHEP::HepVector;
using CLHEP::Hep3Vector;

#include "RootEventData/TExtTrack.h"

class BExtTrack
{
 public:
                        BExtTrack(TExtTrack* trk) {m_trk = trk;}

   const TExtTrack*     getTExtTrack()  const {return m_trk;}

   const int            GetTrackId()    const {return m_trk->GetTrackId();}
   const int            trackId()       const {return m_trk->GetTrackId();}

   //Get track extrapolation data @ Tof layer1.
   const Hep3Vector     tof1Position()          const {
        return Hep3Vector(m_trk->GetTof1PositionX(),
               m_trk->GetTof1PositionY(),m_trk->GetTof1PositionZ());}
   const Hep3Vector     tof1Momentum()          const {
        return Hep3Vector(m_trk->GetTof1MomentumX(),
               m_trk->GetTof1MomentumY(),m_trk->GetTof1MomentumZ());}
   const std::string    tof1VolumeName()        const {
        return std::string((m_trk->GetTof1VolumeName()).Data());}
   const int            tof1VolumeNumber()      const {
        return m_trk->GetTof1VolumeNumber();}
   const double         tof1()                  const {
        return m_trk->GetTof1();}
   const double         tof1Path()              const {
        return m_trk->GetTof1Path();}
   const double         tof1PosSigmaAlongZ()    const {
        return m_trk->GetTof1PosSigmaAlongZ();}
   const double         tof1PosSigmaAlongT()    const {
        return m_trk->GetTof1PosSigmaAlongT();}
   const double         tof1PosSigmaAlongX()    const {
        return m_trk->GetTof1PosSigmaAlongX();}
   const double         tof1PosSigmaAlongY()    const {
        return m_trk->GetTof1PosSigmaAlongY();}
   const HepSymMatrix   tof1ErrorMatrix()       const;

   //Get track extrapolation data @ Tof layer2.
   const Hep3Vector     tof2Position()          const {
        return Hep3Vector(m_trk->GetTof2PositionX(),
               m_trk->GetTof2PositionY(),m_trk->GetTof2PositionZ());}
   const Hep3Vector     tof2Momentum()          const {
        return Hep3Vector(m_trk->GetTof2MomentumX(),
               m_trk->GetTof2MomentumY(),m_trk->GetTof2MomentumZ());}
   const std::string    tof2VolumeName()        const {
        return std::string((m_trk->GetTof2VolumeName()).Data());}
   const int            tof2VolumeNumber()      const {
        return m_trk->GetTof2VolumeNumber();}
   const double         tof2()                  const {
        return m_trk->GetTof2();}
   const double         tof2Path()              const {
        return m_trk->GetTof2Path();}
   const double         tof2PosSigmaAlongZ()    const {
        return m_trk->GetTof2PosSigmaAlongZ();}
   const double         tof2PosSigmaAlongT()    const {
        return m_trk->GetTof2PosSigmaAlongT();}
   const double         tof2PosSigmaAlongX()    const {
        return m_trk->GetTof2PosSigmaAlongX();}
   const double         tof2PosSigmaAlongY()    const {
        return m_trk->GetTof2PosSigmaAlongY();}
   const HepSymMatrix   tof2ErrorMatrix()       const;

   //Get track extrapolation data @ EMC.
   const Hep3Vector     emcPosition()           const {
        return Hep3Vector(m_trk->GetEmcPositionX(),
               m_trk->GetEmcPositionY(),m_trk->GetEmcPositionZ());}
   const Hep3Vector     emcMomentum()           const {
        return Hep3Vector(m_trk->GetEmcMomentumX(),
               m_trk->GetEmcMomentumY(),m_trk->GetEmcMomentumZ());}
   const std::string    emcVolumeName()         const {
        return std::string((m_trk->GetEmcVolumeName()).Data());}
   const int            emcVolumeNumber()       const {
        return m_trk->GetEmcVolumeNumber();}
   const double         emcPosSigmaAlongTheta() const {
        return m_trk->GetEmcPosSigmaAlongTheta();}
   const double         emcPosSigmaAlongPhi()   const {
        return m_trk->GetEmcPosSigmaAlongPhi();}
   const HepSymMatrix   emcErrorMatrix()        const;
   const double         emcPath()               const {
        return m_trk->emcPath();}

   //Get track extrapolation data @ MUC.
   const Hep3Vector     mucPosition()           const {
        return Hep3Vector(m_trk->GetMucPositionX(),
               m_trk->GetMucPositionY(),m_trk->GetMucPositionZ());}
   const Hep3Vector     mucMomentum()           const {
        return Hep3Vector(m_trk->GetMucMomentumX(),
               m_trk->GetMucMomentumY(),m_trk->GetMucMomentumZ());}
   const std::string    mucVolumeName()         const {
        return std::string((m_trk->GetMucVolumeName()).Data());}
   const int            mucVolumeNumber()       const {
        return m_trk->GetMucVolumeNumber();}
   const double         mucPosSigmaAlongZ()     const {
        return m_trk->GetMucPosSigmaAlongZ();}
   const double         mucPosSigmaAlongT()     const {
        return m_trk->GetMucPosSigmaAlongT();}
   const double         mucPosSigmaAlongX()     const {
        return m_trk->GetMucPosSigmaAlongX();}
   const double         mucPosSigmaAlongY()     const {
        return m_trk->GetMucPosSigmaAlongY();}
   const HepSymMatrix   mucErrorMatrix()        const;

 private:
   TExtTrack*           m_trk;
//    int          myParticleType;// it is always = 2

};
#endif
