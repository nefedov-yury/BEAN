#ifndef CONSTANTSH
#define CONSTANTSH

//----------------------------------------------------------------------
static constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

static const double Mjpsi = 3.096916; // 3096.916 +/- 0.011 MeV
static const double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
static const double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV
static const double Meta  = 0.547862; // 547.862  +/- 0.017 MeV
static const double Mk    = 0.493677; // 493.677  +/- 0.016 MeV
static const double Mjpsi2 = SQ(Mjpsi);
static const double Meta2  = SQ(Meta);
static const double Mk2    = SQ(Mk);

static const double R      = 3.; // GeV^{-1} for Blatt-Weisskopf ff

static const double dL = 2*Mk;
static const double dU = 1.10;
#endif
