#include <iostream>
#include <vector>

#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "DstEvtRecTracks.h"   // for RecMdcKalTrack

#include "TrackCorrection/TrackCorrection.h"

using CLHEP::HepVector;
using CLHEP::HepMatrix;
using CLHEP::HepSymMatrix;
using namespace std;

// #define DEBUGPRT
#undef DEBUGPRT

// Helper functions
//-------------------------------------------------------------------------
inline double SQ(double x) {
//-------------------------------------------------------------------------
   return x*x;
}

// *****************************************************************
// ** A macro to create correlated Gaussian-distributed variables **
// *****************************************************************

// corset(): sets up the generation by calculating C from V
//-------------------------------------------------------------------------
static HepMatrix corset(const HepSymMatrix& V) {
//-------------------------------------------------------------------------
   int n = V.num_row();
   HepMatrix C(n,n,0);

   // The Cholesky–Crout algorithm V = C*C^{transpose}, where
   // C is a lower triangular matrix with real and positive diagonal
   // entries (https://en.wikipedia.org/wiki/Cholesky_decomposition)
   for ( int j = 0; j < n; j++ ) {
      double sum = 0;
      for ( int k = 0; k < j; k++ ) {
         sum += SQ(C[j][k]);
      }
      C[j][j] = sqrt(fabs(V[j][j] - sum));
      // Off Diagonal terms
      for ( int i = j+1; i < n; i++ ) {
         double sum = 0;
         for ( int k = 0; k < j; k++ ) {
            sum += C[i][k]*C[j][k];
         }
         C[i][j] = (V[i][j] - sum)/C[j][j];
      }
   }

   return C;
}

// corgen(): generates a set of n random numbers Gaussian-distributed
// with covariance matrix V (V = C*C') and mean values zero
//-------------------------------------------------------------------------
static HepVector corgen(const HepMatrix& C) {
//-------------------------------------------------------------------------
   int n = C.num_row();
   HepVector x(n,0);

   // CLHEP::RandGauss::shootArray is static method depends on
   // the random-engine inside the 'HepRandom' defined globaly
   vector<double> tmp(n);
   CLHEP::RandGauss::shootArray(n,tmp.data());

   // x = C * tmp, where tmp is a vector whose components are
   // independent standard normal variates
   for ( int i = 0; i < n; i++ ) {
      x[i] = 0.0;
      for (int j = 0; j <= i; j++ ) {
         x[i] += C[i][j] * tmp[j];
      }
   }

   return x;
}

// class functions
//-------------------------------------------------------------------------
TrackCorrection::TrackCorrection(int period) {
//-------------------------------------------------------------------------
// see
// docbes3.ihep.ac.cn/~charmoniumgroup/index.php/File:Helix_zhangjielei.pdf
// for Psi' -> K+ K− pi+ pi− in 2009 (TABLE VI) and 2012(TABLE V)

   if ( period == 0 ) {
      cout << "INFO: TrackCorrection for Psi' -> K+K-pi+pi- (2009)"
           << endl;
      pip_mu  = {0.04,  0.12, 0.49}; // pi+
      pip_sig = {1.17,  1.20, 1.14};
      pim_mu  = {0.03, -0.10, 0.49}; // pi-
      pim_sig = {1.17,  1.20, 1.14};
      Kp_mu   = {0.01,  0.12, 0.53}; // K+
      Kp_sig  = {1.17,  1.18, 1.14};
      Km_mu   = {0.06, -0.07, 0.53}; // K-
      Km_sig  = {1.18,  1.18, 1.12};
   } else if ( period == 1 ) {
      cout << "INFO: TrackCorrection for Psi' -> K+K-pi+pi- (2012)"
           << endl;
      pip_mu  = {-0.02,  0.30, -0.06}; // pi+
      pip_sig = { 1.18,  1.21,  1.13};
      pim_mu  = { 0.01, -0.31, -0.07}; // pi-
      pim_sig = { 1.18,  1.21,  1.12};
      Kp_mu   = {-0.05,  0.31, -0.07}; // K+
      Kp_sig  = { 1.20,  1.21,  1.14};
      Km_mu   = { 0.05, -0.30, -0.06}; // K-
      Km_sig  = { 1.20,  1.21,  1.14};
   } else {
      cout << "FATAL ERROR in TrackCorrection:: period= "
           << period << " is not defined yet" << endl;
      exit(EXIT_FAILURE);
   }

   // print table
   auto prtline = [](string name,const vector<double>& mean,
                                 const vector<double>& sig) -> void {
      string fmt = "%5s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n"; 
      printf( fmt.c_str(), name.c_str(), 
              mean[0],sig[0], mean[1],sig[1], mean[2],sig[2] );
   };
   prtline("pi+",pip_mu,pip_sig);
   prtline("pi-",pim_mu,pim_sig);
   prtline(" K+",Kp_mu, Kp_sig);
   prtline(" K-",Km_mu, Km_sig);
}

//-------------------------------------------------------------------------
TrackCorrection::~TrackCorrection() {}
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
void TrackCorrection::calibration(RecMdcKalTrack* trk) const {
//-------------------------------------------------------------------------
   int pid = trk -> getPidType();

   const vector<double> *mean = nullptr, *sigma = nullptr;
   if ( pid == 2 ) { // pion
      if ( trk -> charge() > 0 ) {
         mean = &pip_mu;
         sigma = &pip_sig;
      } else {
         mean = &pim_mu;
         sigma = &pim_sig;
      }
   } else if (pid == 3 ) { // kaon
      if ( trk -> charge() > 0 ) {
         mean = &Kp_mu;
         sigma = &Kp_sig;
      } else {
         mean = &Km_mu;
         sigma = &Km_sig;
      }
   }

   if ( !mean || !sigma ) {
      return; // no corections
   }

   HepSymMatrix h_zerr = trk -> getZError(pid);

   // order of helix parameters: d0, phi0, kappa = 1/Pt, dz, tan(lambda)
   // i=0,1,2 -> j=1,2,4 : j = i+1+i/2;
   HepSymMatrix h_zcal(3,0);
   h_zcal[0][0] = (SQ((*sigma)[0])-1) * h_zerr[1][1]; // phi0
   h_zcal[1][1] = (SQ((*sigma)[1])-1) * h_zerr[2][2]; // kappa
   h_zcal[2][2] = (SQ((*sigma)[2])-1) * h_zerr[4][4]; // tan(lambda)
#ifdef DEBUGPRT
   cout << " h_zcal= " << h_zcal << endl;
#endif

   HepMatrix h_zerr_corr = corset(h_zcal);
   HepVector h_zgen = corgen(h_zerr_corr);
#ifdef DEBUGPRT
   cout << " h_zerr_corr= " << h_zerr_corr << endl;
   cout << " h_zgen= " << h_zgen << endl;
#endif

   HepVector zhelix = trk -> getZHelix(pid);
#ifdef DEBUGPRT
   cout << " befor corrections: zhelix= " << zhelix << endl;
#endif
   zhelix[1] += (*mean)[0] * sqrt(h_zerr[1][1]) + h_zgen[0];
   zhelix[2] += (*mean)[1] * sqrt(h_zerr[2][2]) + h_zgen[1];
   zhelix[4] += (*mean)[2] * sqrt(h_zerr[4][4]) + h_zgen[2];
   trk -> setZHelix(zhelix,pid);

#ifdef DEBUGPRT
   cout << " after corrections: zhelix= " << zhelix << endl;
#endif
}

