// RewTrkPiK.hpp
// RewTrkPi() and RewTrk_K() functions for correction of
// reconstruction efficiency, including PID, kaons and pions.
// This is for Monte Carlo events.

// {{{1 Common function for calculating B-splines
//--------------------------------------------------------------------
double deBoor(int l, double x, int k,
      const std::vector<double>& t, const std::vector<double>& c) {
//--------------------------------------------------------------------
   // Evaluates S(x) = Sum c_i * B_i,k (x)

   // Arguments
   // ---------
   // l: Index of knot interval that contains x
   // x: Position
   // t: Array of knot positions
   // c: Array of control points
   // k: Degree of B-spline

   std::vector<double> d(k+1);
   for( int j = 0; j < k+1; ++j ) {
      d[j] = c[j+l-k];
   }

   for( int r = 1; r < k+1; ++r ) {
      for( int j = k; j > r-1; --j ) {
         double alpha = (x-t[j+l-k]) / (t[j+1+l-r] - t[j+l-k]);
         d[j] = (1.0 - alpha) * d[j-1] + alpha * d[j];
      }
   }

   return d[k];
}

//--------------------------------------------------------------------
double bspline(double x, int k,
      const std::vector<double>& t, const std::vector<double>& c) {
//--------------------------------------------------------------------
   int n = t.size() - k - 1;
   // find index l in range [k,n) what t[l] <= x < t[l+1]
   int l = k+1;
   while( x >= t[l] && l != n ) {
      l += 1;
   }
   return deBoor(l-1, x, k, t, c);
}

// {{{1 RewTrkPi()
double Bsp_Pim_2009(const double* xx, const double* p);
double Bsp_Pim_2012(const double* xx, const double* p);
double Bsp_Pim_2021(const double* xx, const double* p);
double Bsp_Pip_2009(const double* xx, const double* p);
double Bsp_Pip_2012(const double* xx, const double* p);
double Bsp_Pip_2021(const double* xx, const double* p);

//--------------------------------------------------------------------
double RewTrkPi(int DataPeriod, double Pt, double Z) {
//--------------------------------------------------------------------
   // Corrections for the efficiency of reconstruction a pion having
   // transverse momentum Pt and sign Z.
   // The return value is the weight for the MC event.
   // v709, DelClonedTrk, helix corrections
   // Note: This function works well for MC without helix corrections

   const double Ptmin = 0.025, Ptmax = 0.425;
   Pt = max( Ptmin, Pt );
   Pt = min( Ptmax, Pt );
   const double xx[] {Pt};
   const double* p = nullptr;

   double W = 1.;
   if ( DataPeriod == 2021 ) {
      if ( Z > 0 ) {
         W = Bsp_Pip_2021(xx,p);
      } else {
         W = Bsp_Pim_2021(xx,p);
      }
   } else if ( DataPeriod == 2012 ) {
      if ( Z > 0 ) {
         W = Bsp_Pip_2012(xx,p);
      } else {
         W = Bsp_Pim_2012(xx,p);
      }
   } else if ( DataPeriod == 2009 ) {
      if ( Z > 0 ) {
         W = Bsp_Pip_2009(xx,p);
      } else {
         W = Bsp_Pim_2009(xx,p);
      }
   }
   return W;
}

// {{{1 RewTrk_K()
double Bsp_Km_2009(const double* xx, const double* p);
double Bsp_Km_2012(const double* xx, const double* p);
double Bsp_Km_2021(const double* xx, const double* p);
double Bsp_Kp_2009(const double* xx, const double* p);
double Bsp_Kp_2012(const double* xx, const double* p);
double Bsp_Kp_2021(const double* xx, const double* p);

//--------------------------------------------------------------------
double RewTrk_K(int DataPeriod, double Pt, double Z) {
//--------------------------------------------------------------------
   // Corrections for the efficiency of reconstruction a kaon having
   // transverse momentum Pt and sign Z.
   // The return value is the weight for the MC event.
   // v709, DelClonedTrk, helix corrections

   const double Ptmin = 0.05, Ptmax = 1.45;
   Pt = max( Ptmin, Pt );
   Pt = min( Ptmax, Pt );
   const double xx[] {Pt};
   const double* p = nullptr;

   double W = 1.;
   if ( DataPeriod == 2021 ) {
      if ( Z > 0 ) {
         W = Bsp_Kp_2021(xx,p);
      } else {
         W = Bsp_Km_2021(xx,p);
      }
   } else if ( DataPeriod == 2012 ) {
      if ( Z > 0 ) {
         W = Bsp_Kp_2012(xx,p);
      } else {
         W = Bsp_Km_2012(xx,p);
      }
   } else if ( DataPeriod == 2009 ) {
      if ( Z > 0 ) {
         W = Bsp_Kp_2009(xx,p);
      } else {
         W = Bsp_Km_2009(xx,p);
      }
   }

   return W;
}

// {{{1 Bsp_K... functions v709n4_eff

double Bsp_Kp_2021(const double* xx, const double* p){
   const std::vector<double> t {
   0.13467, 0.13467, 0.13467, 0.13467, 0.304449, 0.4018, 0.500252, 0.798304, 1.38142, 1.38142, 1.38142, 1.38142
   };
   const std::vector<double> c {
   0.8273154162878921, 1.0402039207840328, 0.9290896739664779, 0.9776059512908271, 0.9897608209891136, 0.9579184055637083, 0.9809913906409862, 0.9757542667346214, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return 0.9995*bspline(xx[0],k,t,c);
};


double Bsp_Km_2021(const double* xx, const double* p){
   const std::vector<double> t {
   0.134633, 0.134633, 0.134633, 0.134633, 0.30467, 0.500292, 0.798233, 1.380349, 1.380349, 1.380349, 1.380349
   };
   const std::vector<double> c {
   0.9930696133617953, 1.0299200191880586, 0.9448655141046931, 1.0036705201166607, 0.9651979902160243, 0.9613484313621373, 0.9729196513714798, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};


double Bsp_Kp_2012(const double* xx, const double* p){
   const std::vector<double> t {
   0.135089, 0.135089, 0.135089, 0.135089, 1.381341, 1.381341, 1.381341, 1.381341
   };
   const std::vector<double> c {
   0.9542646773905216, 1.0282897771104642, 0.9397031992269845, 0.9882640524187752, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return 0.999*bspline(xx[0],k,t,c);
};


double Bsp_Km_2012(const double* xx, const double* p){
   const std::vector<double> t {
   0.13456, 0.13456, 0.13456, 0.13456, 0.304414, 0.500039, 0.798589, 1.380842, 1.380842, 1.380842, 1.380842
   };
   const std::vector<double> c {
   1.1148395728875884, 1.0082820203852463, 0.9629666287457717, 1.0042826460926535, 0.9703461349766389, 0.9621819851406338, 0.9955085707248393, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return 0.9985*bspline(xx[0],k,t,c);
};


double Bsp_Kp_2009(const double* xx, const double* p){
   const std::vector<double> t {
   0.135264, 0.135264, 0.135264, 0.135264, 1.382833, 1.382833, 1.382833, 1.382833
   };
   const std::vector<double> c {
   0.9594053434164469, 1.0452079399874556, 0.9501202365052208, 0.9847378231958747, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return 0.9975*bspline(xx[0],k,t,c);
};


double Bsp_Km_2009(const double* xx, const double* p){
   const std::vector<double> t {
   0.133898, 0.133898, 0.133898, 0.133898, 1.374697, 1.374697, 1.374697, 1.374697
   };
   const std::vector<double> c {
   0.9742589670446634, 1.0727083889774363, 0.8889879804687788, 1.037858558733724, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};

// {{{1 Bsp_Pi... functions v709n4_eff

double Bsp_Pip_2021(const double* xx, const double* p){
   const std::vector<double> t {
   0.064741, 0.064741, 0.064741, 0.064741, 0.382702, 0.382702, 0.382702, 0.382702
   };
   const std::vector<double> c {
   0.9746388546604619, 1.0025200838915187, 0.9654166110874602, 0.9838484961597891, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return 0.9995*bspline(xx[0],k,t,c);
};


double Bsp_Pim_2021(const double* xx, const double* p){
   const std::vector<double> t {
   0.064765, 0.064765, 0.064765, 0.064765, 0.382994, 0.382994, 0.382994, 0.382994
   };
   const std::vector<double> c {
   0.9899979904584183, 0.9792846062293079, 0.9993223171938992, 0.978371895930028, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return 0.9998*bspline(xx[0],k,t,c);
};


double Bsp_Pip_2012(const double* xx, const double* p){
   const std::vector<double> t {
   0.065247, 0.065247, 0.065247, 0.065247, 0.382331, 0.382331, 0.382331, 0.382331
   };
   const std::vector<double> c {
   1.0219546209256012, 0.9698895724616781, 0.978592112743844, 0.9883104984875716, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return 0.999*bspline(xx[0],k,t,c);
};


double Bsp_Pim_2012(const double* xx, const double* p){
   const std::vector<double> t {
   0.065473, 0.065473, 0.065473, 0.065473, 0.249208, 0.382422, 0.382422, 0.382422, 0.382422
   };
   const std::vector<double> c {
   1.0472791996405002, 0.9783651107097368, 0.9989427729181226, 0.9880863977272629, 0.9895290715065959, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return 0.999*bspline(xx[0],k,t,c);
};


double Bsp_Pip_2009(const double* xx, const double* p){
   const std::vector<double> t {
   0.064557, 0.064557, 0.064557, 0.064557, 0.383181, 0.383181, 0.383181, 0.383181
   };
   const std::vector<double> c {
   0.977090824618079, 1.0251174636345728, 0.9482690797103407, 0.9982213536371004, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return 0.999*bspline(xx[0],k,t,c);
};


double Bsp_Pim_2009(const double* xx, const double* p){
   const std::vector<double> t {
   0.064522, 0.064522, 0.064522, 0.064522, 0.249287, 0.383226, 0.383226, 0.383226, 0.383226
   };
   const std::vector<double> c {
   1.018848895922333, 0.9602742621692294, 1.0084634280073208, 0.9923457896714548, 0.9904332666330907, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return 0.9985*bspline(xx[0],k,t,c);
};

