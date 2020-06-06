// #include <iostream>
#include <cmath>

#include "BMdcKalTrack.h"

using namespace std;
using namespace CLHEP;

//-----------------------------------------------------------------------------
const int BMdcKalTrack::getCharge(const int pid) const 
//-----------------------------------------------------------------------------
{ 
  int charge = 0;
  double kappa = getZHelix_pid(2,pid);
  if (kappa > 0.0000000001)
    charge = 1 ;
  else if (kappa < -0.0000000001)
    charge = -1;
  else
    charge = 0;
  return charge; 
}

//-----------------------------------------------------------------------------
const double BMdcKalTrack::pxy() const 
//-----------------------------------------------------------------------------
{ 
  double kappa = getZHelix_pid(2,m_pid);
  if(kappa != 0) 
    return 1./fabs(kappa);
  else return 0.;
}
 
//-----------------------------------------------------------------------------
const double BMdcKalTrack::px() const 
//-----------------------------------------------------------------------------
{ 
  double phi0 = getZHelix_pid(1,m_pid);
  return pxy()*(-sin(phi0));
}
 
//-----------------------------------------------------------------------------
const double BMdcKalTrack::py() const 
//-----------------------------------------------------------------------------
{ 
  double phi0 = getZHelix_pid(1,m_pid);
  return pxy()*cos(phi0);
}
 
//-----------------------------------------------------------------------------
const double BMdcKalTrack::pz() const 
//-----------------------------------------------------------------------------
{ 
  double tanl = getZHelix_pid(4,m_pid);
  return  pxy()*tanl;
}

//-----------------------------------------------------------------------------
const double BMdcKalTrack::p() const 
//-----------------------------------------------------------------------------
{ 
  double tanl = getZHelix_pid(4,m_pid);
  return  pxy()*sqrt(1. + tanl*tanl);
}

//-----------------------------------------------------------------------------
const double BMdcKalTrack::theta() const 
//-----------------------------------------------------------------------------
{
  double tanl = getZHelix_pid(4,m_pid);
  return acos(tanl/sqrt(1. + tanl*tanl));
}

//-----------------------------------------------------------------------------
const double BMdcKalTrack::phi() const 
//-----------------------------------------------------------------------------
{
  double phi0 = getZHelix_pid(1,m_pid);
//   return atan2(py(),px());
  return atan2( cos(phi0),-sin(phi0) );
}

//-----------------------------------------------------------------------------
const double BMdcKalTrack::x(int pid) const 
//-----------------------------------------------------------------------------
{
  double dr   = getZHelix_pid(0,pid);
  double phi0 = getZHelix_pid(1,pid);
  return  dr*cos(phi0);
}

//-----------------------------------------------------------------------------
const double BMdcKalTrack::y(int pid) const 
//-----------------------------------------------------------------------------
{
  double dr   = getZHelix_pid(0,pid);
  double phi0 = getZHelix_pid(1,pid);
  return  dr*sin(phi0);
}

//-----------------------------------------------------------------------------
const double BMdcKalTrack::z(int pid) const
//-----------------------------------------------------------------------------
{
  double dz   = getZHelix_pid(3,pid);
  return  dz;
}

//-----------------------------------------------------------------------------
const double BMdcKalTrack::r() const
//-----------------------------------------------------------------------------
{ 
  double dr   = getZHelix_pid(0,m_pid);
  return  fabs(dr);
}

//-----------------------------------------------------------------------------
const HepVector BMdcKalTrack::getZHelix(const int pid) const 
//-----------------------------------------------------------------------------
{
  HepVector tmp(5);
  for(int i = 0; i < 5; i++) {
    tmp[i] = getZHelix_pid(i,pid);
  }
  return tmp;
}

//-----------------------------------------------------------------------------
const HepSymMatrix BMdcKalTrack::getZError(const int pid) const 
//-----------------------------------------------------------------------------
{
  HepSymMatrix tmp(5);
  for(int i = 0; i < 5; i++) {
    for(int j = 0; j <= i; j++) {
      tmp[i][j] = getZError_pid(i,j,pid);
    }
  }
  return tmp;
}

//-----------------------------------------------------------------------------
const HepVector BMdcKalTrack::getFHelix(const int pid) const 
//-----------------------------------------------------------------------------
{
  HepVector tmp(5);
  for(int i = 0; i < 5; i++) {
    tmp[i] = getFHelix_pid(i,pid);
  }
  return tmp;
}

//-----------------------------------------------------------------------------
const HepSymMatrix BMdcKalTrack::getFError(const int pid) const 
//-----------------------------------------------------------------------------
{
  HepSymMatrix tmp(5);
  for(int i = 0; i < 5; i++) {
    for(int j = 0; j <= i; j++) {
      tmp[i][j] = getFError_pid(i,j,pid);
    }
  }
  return tmp;
}

// ther are no getPoca.. functions for RootEventData >= 6.5.5
// //-----------------------------------------------------------------------------
// const HepPoint3D BMdcKalTrack::getPoca(const int pid) const 
// //-----------------------------------------------------------------------------
// {
//   return HepPoint3D(getPoca_pid(0,pid),getPoca_pid(1,pid),getPoca_pid(2,pid));
// }

//-----------------------------------------------------------------------------
const Hep3Vector BMdcKalTrack::p3() const 
//-----------------------------------------------------------------------------
{
  double pt = pxy();
  double phi0 = getZHelix_pid(1,m_pid);
  double tanl = getZHelix_pid(4,m_pid);
  return Hep3Vector(-pt*sin(phi0),pt*cos(phi0),pt*tanl);
}

//-----------------------------------------------------------------------------
const HepPoint3D BMdcKalTrack::x3() const 
//-----------------------------------------------------------------------------
{
  return HepPoint3D(x(),y(),z());
}

//-----------------------------------------------------------------------------
const HepLorentzVector BMdcKalTrack::p4(double mass) const 
//-----------------------------------------------------------------------------
{
  double pt = pxy();
  double phi0 = getZHelix_pid(1,m_pid);
  double tanl = getZHelix_pid(4,m_pid);
  double E = sqrt(pt*pt*(1. + tanl*tanl) + mass*mass);
  return
    HepLorentzVector(-pt*sin(phi0),pt*cos(phi0),pt*tanl,E);
}
