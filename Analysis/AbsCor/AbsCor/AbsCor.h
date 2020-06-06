#ifndef Analysis_AbsCor_H
#define Analysis_AbsCor_H

#if   (BOSS_VER >= 661 && BOSS_VER <= 665)
#define AbsCor_00_00_28
#elif (BOSS_VER >= 702 && BOSS_VER <= 704)
#define AbsCor_00_00_36
#else
#error "unknown BOSS version"
#endif

#include <string>
#include "ReadDst.h"

#ifdef  AbsCor_00_00_36
#include "TGraph2DErrors.h"
#endif

//-----------------------------------------------------------------------------
// This class is a replacement for class AbsCor from the Boss
// AbsCor: <=> boss/Analysis/PhotonCor/AbsCor
//-----------------------------------------------------------------------------

class AbsCor {

public:
                // for default values see declareProperty() in BOSS
                AbsCor( const std::string& path,
                        bool _usetof = true,
                        bool _dodatacor = true,
                        bool _edgecor = false
#ifdef  AbsCor_00_00_36
                        ,
                        bool _hotcellmask = false,
                        bool _dopi0Cor = true,
                        bool _MCuseTof = true
#endif
                      );

  void          AbsorptionCorrection(ReadDst* selector);
  void          SuppressHotCrystals(ReadDst* selector);

private:
  bool          usetof;
  bool          dodatacor;
  bool          edgecor;

  double        ai[4];

  int hrunstart[10];
  int hrunend[10];
  int hcell[10];

#ifdef  AbsCor_00_00_36
  bool          hotcellmask;
  bool          dopi0Cor;
  bool          MCuseTof;

  double e25min[28];
  double e25max[28];
  // Shower energy correction
  TGraph2DErrors *dt;
  // Energy error
  TGraph2DErrors *dtErr;

  double ECorrMC(double eg, double theid) const;
  double ErrMC(double eg, double theid) const;
  double E25min(int n) const { return e25min[n]; }
  double E25max(int n) const { return e25max[n]; }
#endif

};
#endif
