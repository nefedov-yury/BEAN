#ifndef _BMdcKalTrack_h
#define _BMdcKalTrack_h

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// BMdcKalTrack                                                         //
//                                                                      //
// This is adapter between RootEventData/TMdcKalTrack.h and             //
// DstEvent/DstMdcKalTrack.h and partialy MdcRecEvent/RecMdcKalTrack.h  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <CLHEP/Matrix/SymMatrix.h>
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Geometry/Point3D.h>
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
   typedef HepGeom::Point3D<double> HepPoint3D;
#endif
using CLHEP::HepSymMatrix;
using CLHEP::HepVector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

#include "RootEventData/TMdcKalTrack.h"

class BMdcKalTrack
{
 public:
   enum PidType
   {
     null = -1,
     electron = 0,
     muon = 1,
     pion = 2,
     kaon = 3,
     proton = 4
   };

 public:
   BMdcKalTrack(TMdcKalTrack* root_trk) {
     m_trk = root_trk; m_pid = pion;
   }

   const TMdcKalTrack* getTMdcKalTrack() const {return m_trk;}

   void setPidType(PidType pidType = pion) {m_pid = pidType;} 
   PidType getPidType()                    {return m_pid;}
   
   const int trackId()  const {return m_trk->getTrackId();}
   // const double mass()  const { return m_mass[m_pid];  }
   const int charge()   const {return getCharge(m_pid);}
   const double pxy()   const;
   const double px()    const;
   const double py()    const;
   const double pz()    const;
   const double p()     const;
   const double theta() const;
   const double phi()   const;
     
   const double x()     const   {return x(m_pid);}   
   const double y()     const   {return y(m_pid);}
   const double z()     const   {return z(m_pid);}     

   const double x(int pid)  const;
   const double y(int pid)  const;
   const double z(int pid)  const;
   const double r()         const;

   const int    stat()        const {return m_trk->getStat(m_pid);      }
   const double chi2()        const {return m_trk->getChisq(m_pid);     }
   const int    ndof()        const {return m_trk->getNdf(m_pid);       }
#if BOSS_VER < 663
   const int    nster()       const {return m_trk->getNster(m_pid);     }
   const int    firstLayer()  const {return m_trk->getFirstLayer(m_pid);}
   const int    lastLayer()   const {return m_trk->getLastLayer(m_pid); }
#endif
   
   const double dr()    const {return getZHelix_pid(0,m_pid);}
   const double fi0()   const {return getZHelix_pid(1,m_pid);}
   const double kappa() const {return getZHelix_pid(2,m_pid);}
   const double dz()    const {return getZHelix_pid(3,m_pid);}
   const double tanl()  const {return getZHelix_pid(4,m_pid);}

   const HepVector helix()    const {return getZHelix(m_pid);}
   const HepSymMatrix err()   const {return getZError(m_pid);}
   const HepVector fhelix()   const {return getFHelix(m_pid);}
   const HepSymMatrix ferr()  const {return getFError(m_pid);}
//    const HepPoint3D poca()    const {return getPoca(m_pid);}
   const Hep3Vector p3()      const;
   const HepPoint3D x3()      const;

//    const HepLorentzVector p4() const;
   const HepLorentzVector p4(double mass) const;
   
   const int getTrackId()                 const {return m_trk->getTrackId();}
   const int getCharge(const int pid)     const;
   const int getStat(const int pid)       const {return m_trk->getStat(pid);}
   const double getChisq(const int pid)   const {return m_trk->getChisq(pid);}
   const int getNdf(const int pid)        const {return m_trk->getNdf(pid);}
#if BOSS_VER < 663
   const int getNster(const int pid)      const {return m_trk->getNster(pid);}
   const int getFirstLayer(const int pid) const 
                {return m_trk->getFirstLayer(pid);}
   const int getLastLayer(const int pid)  const 
                {return m_trk->getLastLayer(pid); }
#endif

   const HepVector getZHelix(const int pid)     const;
   const HepSymMatrix getZError(const int pid)  const;
   const HepVector getFHelix(const int pid)     const;
   const HepSymMatrix getFError(const int pid)  const;
//    const HepPoint3D getPoca(const int pid)      const; 
   
   // this is functions from MdcRecEvent/RecMdcKalTrack.h
   const HepVector    getZHelix() const { return getZHelix(pion); }
   const HepSymMatrix getZError() const { return getZError(pion); }
   const HepVector    getFHelix() const { return getFHelix(pion); }
   const HepSymMatrix getFError() const { return getFError(pion); }

   const HepVector    getZHelixE() const { return getZHelix(electron); }
   const HepSymMatrix getZErrorE() const { return getZError(electron); }
   const HepVector    getFHelixE() const { return getFHelix(electron); }
   const HepSymMatrix getFErrorE() const { return getFError(electron); }

   const HepVector    getZHelixMu() const { return getZHelix(muon); }
   const HepSymMatrix getZErrorMu() const { return getZError(muon); }
   const HepVector    getFHelixMu() const { return getFHelix(muon); }
   const HepSymMatrix getFErrorMu() const { return getFError(muon); }

   const HepVector    getZHelixK() const { return getZHelix(kaon); }
   const HepSymMatrix getZErrorK() const { return getZError(kaon); }
   const HepVector    getFHelixK() const { return getFHelix(kaon); }
   const HepSymMatrix getFErrorK() const { return getFError(kaon); }

   const HepVector    getZHelixP() const { return getZHelix(proton); }
   const HepSymMatrix getZErrorP() const { return getZError(proton); }
   const HepVector    getFHelixP() const { return getFHelix(proton); }
   const HepSymMatrix getFErrorP() const { return getFError(proton); }

//    const HepPoint3D getPocaE()  const { return getPoca(electron); }
//    const HepPoint3D getPocaMu() const { return getPoca(muon); }
//    const HepPoint3D getPoca()   const { return getPoca(pion); }
//    const HepPoint3D getPocaK()  const { return getPoca(kaon); }
//    const HepPoint3D getPocaP()  const { return getPoca(proton); }

 private:
   TMdcKalTrack*        m_trk;
   PidType              m_pid;

   double getZHelix_pid(int i,int pid = 2) const {
      switch( pid ) {
         case 0: return m_trk->getZHelixE(i);
         case 1: return m_trk->getZHelixMu(i);
         case 2: return m_trk->getZHelix(i);
         case 3: return m_trk->getZHelixK(i);
         case 4: return m_trk->getZHelixP(i);
      }
      return m_trk->getZHelix(i);
   }

   double getZError_pid(int i,int j,int pid = 2) const {
      switch( pid ) {
         case 0: return m_trk->getZErrorE(i,j);
         case 1: return m_trk->getZErrorMu(i,j);
         case 2: return m_trk->getZError(i,j);
         case 3: return m_trk->getZErrorK(i,j);
         case 4: return m_trk->getZErrorP(i,j);
      }
      return m_trk->getZError(i,j);
   }

// ther are no getPoca.. functions for RootEventData >= 6.5.5
//    double getPoca_pid(int i,int pid = 2) const {
//       switch( pid ) {
//          case 0: return m_trk->getPocaE(i);
//          case 1: return m_trk->getPocaMu(i);
//          case 2: return m_trk->getPoca(i);
//          case 3: return m_trk->getPocaK(i);
//          case 4: return m_trk->getPocaP(i);
//       }
//       return m_trk->getPoca(i);
//    }
         
   double getFHelix_pid(int i,int pid = 2) const {
      switch( pid ) {
         case 0: return m_trk->getFHelixE(i);
         case 1: return m_trk->getFHelixMu(i);
         case 2: return m_trk->getFHelix(i);
         case 3: return m_trk->getFHelixK(i);
         case 4: return m_trk->getFHelixP(i);
      }
      return m_trk->getFHelix(i);
   }

   double getFError_pid(int i,int j,int pid = 2) const {
      switch( pid ) {
         case 0: return m_trk->getFErrorE(i,j);
         case 1: return m_trk->getFErrorMu(i,j);
         case 2: return m_trk->getFError(i,j);
         case 3: return m_trk->getFErrorK(i,j);
         case 4: return m_trk->getFErrorP(i,j);
      }
      return m_trk->getFError(i,j);
   }
};
#endif
