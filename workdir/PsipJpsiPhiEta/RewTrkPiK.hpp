// RewTrkPiK.hpp - functions with corrections for the efficiency of
// Kaon/Pion reconstruction,including PID, for the MC events

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

// {{{1 Bsp_Pi.. functions

double Bsp_Pim_2009(const double* xx, const double* p){
   const std::vector<double> t {
   0.058022, 0.058022, 0.058022, 0.058022, 0.151757, 0.200466, 0.249164, 0.382652, 0.382652, 0.382652, 0.382652
   };
   const std::vector<double> c {
   1.0316555558631755, 0.9865314701326018, 0.9747768800627127, 0.9891757114964591, 1.0006523836776509, 0.9904527953279562, 0.9927324003081914, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};


double Bsp_Pim_2012(const double* xx, const double* p){
   const std::vector<double> t {
   0.057535, 0.057535, 0.057535, 0.057535, 0.15189, 0.249109, 0.383025, 0.383025, 0.383025, 0.383025
   };
   const std::vector<double> c {
   1.0436166312664574, 1.0065266232048038, 0.9873677021707664, 0.9967643150045482, 0.985904029951522, 0.9911131487590454, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};


double Bsp_Pim_2021(const double* xx, const double* p){
   const std::vector<double> t {
   0.057602, 0.057602, 0.057602, 0.057602, 0.151986, 0.249138, 0.382836, 0.382836, 0.382836, 0.382836
   };
   const std::vector<double> c {
   0.9776852365159849, 0.9933964171034397, 0.9856783072562179, 0.9873150959289043, 0.9908367733596619, 0.9785150689442943, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};


double Bsp_Pip_2009(const double* xx, const double* p){
   const std::vector<double> t {
   0.056894, 0.056894, 0.056894, 0.056894, 0.382392, 0.382392, 0.382392, 0.382392
   };
   const std::vector<double> c {
   0.9655878123381549, 1.0321984624332223, 0.9492114950396141, 0.9959800358884756, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};


double Bsp_Pip_2012(const double* xx, const double* p){
   const std::vector<double> t {
   0.057352, 0.057352, 0.057352, 0.057352, 0.383106, 0.383106, 0.383106, 0.383106
   };
   const std::vector<double> c {
   1.0221329021709886, 0.9670032969363239, 0.9843540593903949, 0.9855822686515928, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};


double Bsp_Pip_2021(const double* xx, const double* p){
   const std::vector<double> t {
   0.057348, 0.057348, 0.057348, 0.057348, 0.152039, 0.249067, 0.382645, 0.382645, 0.382645, 0.382645
   };
   const std::vector<double> c {
   0.9689421551447139, 0.9903289648878736, 0.9890672521943833, 0.9736522384127888, 0.9871843950297, 0.9839459214685715, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};

// {{{1 Bsp_K.. functions

double Bsp_Km_2009(const double* xx, const double* p){
   const std::vector<double> t {
   0.118622, 0.118622, 0.118622, 0.118622, 1.382442, 1.382442, 1.382442, 1.382442
   };
   const std::vector<double> c {
   0.9784535992948784, 1.0693050402458204, 0.893493674479572, 1.0448066546414638, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};


double Bsp_Km_2012(const double* xx, const double* p){
   const std::vector<double> t {
   0.117178, 0.117178, 0.117178, 0.117178, 0.500067, 0.798367, 1.379634, 1.379634, 1.379634, 1.379634
   };
   const std::vector<double> c {
   1.0799010496446086, 0.9473427969836942, 1.0077263482219598, 0.9774186583594245, 0.9627662306828378, 1.0027369075089614, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};


double Bsp_Km_2021(const double* xx, const double* p){
   const std::vector<double> t {
   0.117332, 0.117332, 0.117332, 0.117332, 0.303751, 0.499794, 0.79848, 1.380862, 1.380862, 1.380862, 1.380862
   };
   const std::vector<double> c {
   0.990537714100703, 1.0387779401627641, 0.9434804097989746, 1.005341088573763, 0.9666672319662621, 0.9642185942297706, 0.9735501992240055, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};


double Bsp_Kp_2009(const double* xx, const double* p){
   const std::vector<double> t {
   0.113101, 0.113101, 0.113101, 0.113101, 1.376656, 1.376656, 1.376656, 1.376656
   };
   const std::vector<double> c {
   0.9604929064711097, 1.0467127931608082, 0.9516843852269559, 0.9891691275626162, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};


double Bsp_Kp_2012(const double* xx, const double* p){
   const std::vector<double> t {
   0.113044, 0.113044, 0.113044, 0.113044, 0.2082, 0.303237, 0.40132, 0.500595, 0.798721, 1.380111, 1.380111, 1.380111, 1.380111
   };
   const std::vector<double> c {
   0.959938402702726, 0.993963804376123, 0.9915589986059254, 0.9504956406850832, 0.9857602985764935, 0.9940152435966706, 0.9771501193927454, 0.965161710580633, 0.98691548872554, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};


double Bsp_Kp_2021(const double* xx, const double* p){
   const std::vector<double> t {
   0.114426, 0.114426, 0.114426, 0.114426, 0.303828, 0.401088, 0.499803, 0.798711, 1.381018, 1.381018, 1.381018, 1.381018
   };
   const std::vector<double> c {
   0.8225745729636734, 1.038263692572686, 0.9356949133683078, 0.9771068098547581, 0.9929286762684737, 0.9607878255322894, 0.9807011197315845, 0.977736586412017, 0.0, 0.0, 0.0, 0.0
   };
   int k = 3;
   return bspline(xx[0],k,t,c);
};

