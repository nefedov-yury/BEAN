// #include <iostream>
#include <cmath>

#include "BMdcTrack.h"

using namespace std;
using namespace CLHEP;

//-----------------------------------------------------------------------------
const HepVector BMdcTrack::helix() const 
//-----------------------------------------------------------------------------
{
  HepVector tmp(5);
  for(int i = 0; i < 5; i++) {
    tmp[i] = helix(i);
  }
  return tmp;
}

//-----------------------------------------------------------------------------
const HepSymMatrix BMdcTrack::err() const 
//-----------------------------------------------------------------------------
{
  HepSymMatrix tmp(5);
   // m_err[] is stored as a lower triangular matrix (see DstMdcTrack.cxx)
  int k = 0;
  for(int i = 0; i < 5; i++) {
    for(int j = 0; j <= i; j++) {
      tmp[i][j] = err(k++);
    }
  }
  return tmp;
}

//-----------------------------------------------------------------------------
const Hep3Vector BMdcTrack::p3() const 
//-----------------------------------------------------------------------------
{
  double pt = pxy();
  double fi0 = helix(1);
  double tanl = helix(4);
  return Hep3Vector(-pt*sin(fi0),pt*cos(fi0),pt*tanl);
}

//-----------------------------------------------------------------------------
const HepLorentzVector BMdcTrack::p4(const double mass) const 
//-----------------------------------------------------------------------------
{
  double pt = pxy();
  double fi0 = helix(1);
  double tanl = helix(4);
  double E = sqrt(pt*pt*(1. + tanl*tanl) + mass*mass);

  return HepLorentzVector(-pt*sin(fi0),pt*cos(fi0),pt*tanl,E);
}

//-----------------------------------------------------------------------------
const HepPoint3D BMdcTrack::poca() const 
//-----------------------------------------------------------------------------
{
  // see RootCnvSvc: Dst/MdcTrackCnv.cxx poca == x3 ???
  return HepPoint3D(x(),y(),z());
}

//-----------------------------------------------------------------------------
const HepPoint3D BMdcTrack::x3() const 
//-----------------------------------------------------------------------------
{
  return HepPoint3D(x(),y(),z());
}
