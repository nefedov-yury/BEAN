#ifndef TrackCorrection_H
#define TrackCorrection_H

#include <vector>

class BMdcKalTrack; // == RecMdcKalTrack in bean

class TrackCorrection {
   public:
      TrackCorrection(int period);
     ~TrackCorrection();
      void calibration(RecMdcKalTrack* trk) const;

   private:
      // vectors of correction for mean (_mu) and sigmas (_sig)
      // of helix parameters: [0]->phi0, [1]->kappa, [2]->tan(lambda)
      std::vector<double> pip_mu, pip_sig, pim_mu, pim_sig; //pi+/-
      std::vector<double> Kp_mu,  Kp_sig,  Km_mu,  Km_sig;  //K+/-
};

#endif
