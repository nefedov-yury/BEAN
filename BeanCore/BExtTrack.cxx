// #include <iostream>

#include "BExtTrack.h"

using namespace std;
using namespace CLHEP;

//-----------------------------------------------------------------------------
const HepSymMatrix BExtTrack::tof1ErrorMatrix() const 
//-----------------------------------------------------------------------------
{
  // see ExtTrackCnv.cxx
  HepSymMatrix tmp(6);
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j <= i; j++) {
      tmp[i][j] = m_trk->GetTof1ErrorMatrix(i,j);
    }
  }
  return tmp;
}

//-----------------------------------------------------------------------------
const HepSymMatrix BExtTrack::tof2ErrorMatrix() const 
//-----------------------------------------------------------------------------
{
  // see ExtTrackCnv.cxx
  HepSymMatrix tmp(6);
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j <= i; j++) {
      tmp[i][j] = m_trk->GetTof2ErrorMatrix(i,j);
    }
  }
  return tmp;
}

//-----------------------------------------------------------------------------
const HepSymMatrix BExtTrack::emcErrorMatrix() const 
//-----------------------------------------------------------------------------
{
  // see ExtTrackCnv.cxx
  HepSymMatrix tmp(6);
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j <= i; j++) {
      tmp[i][j] = m_trk->GetEmcErrorMatrix(i,j);
    }
  }
  return tmp;
}

//-----------------------------------------------------------------------------
const HepSymMatrix BExtTrack::mucErrorMatrix() const 
//-----------------------------------------------------------------------------
{
  // see ExtTrackCnv.cxx
  HepSymMatrix tmp(6);
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j <= i; j++) {
      tmp[i][j] = m_trk->GetMucErrorMatrix(i,j);
    }
  }
  return tmp;
}
