//----------------------------------------------------------------------
double ReWeightTrkPid(int DataPeriod, int Kp, double Pt) {
//----------------------------------------------------------------------
// The corrections are based on production-11(helix corrections for MC).
// They do not dependent of the sign of the particle and cos(Theta) of
// the track. Input parameters are following.
// Kp is the type of the particle: 1 for kaon and 0 for pion.
// Pt is the transverse momentum if the particle.
// The return value is the weight of MC event with such a particle. 

   // parameters
   static const vector<double> K09 {0.991,-0.018};
   static const double K09_first = 0.922;
   static const vector<double> pi09 {0.9872,0.0243};

   static const vector<double> K12 {0.9891,-0.0005,0.0074,0.0111,0.0102};
   static const vector<double> pi12 {0.9851,0.032};

   double W = 1.;

   if ( Kp == 1 ) {             // kaons
      const double Ptmin = 0.1, Ptmax = 1.4;
      Pt = max( Ptmin, Pt );
      Pt = min( Ptmax, Pt );

      int nch = (DataPeriod == 2009) ? 1 : 4;
      double xmin = Ptmin, xmax = Ptmax;
      auto Lchb = [nch,xmin,xmax](double xx, const double* p) {
         if (nch == 0) { return p[0]; }
         // [xmin,xmax] -> [-1,+1]
         double x = (2*xx-xmin-xmax)/(xmax-xmin);
         double sum = p[0] + x*p[1];
         if (nch == 1) { return sum; }
         double T0 = 1, T1 = x;
         for ( int i = 2; i <= nch; ++i ) {
            double Tmp = 2*x*T1 - T0;
            sum += p[i]*Tmp;
            T0 = T1;
            T1 = Tmp;
         }
         return sum;
      };

      if ( DataPeriod == 2009 ) {
         W = Lchb(Pt,K09.data());
         if ( Pt < 0.2 ) {
            W = K09_first;
         }
      } else if ( DataPeriod == 2012 ) {
         W = Lchb( Pt, K12.data() );
      }
   } else if ( Kp == 0 ) {      // pions
      const double Ptmin = 0.05, Ptmax = 0.4;
      Pt = max( Ptmin, Pt );
      Pt = min( Ptmax, Pt );

      auto CUBE = [](double x)-> double{return x*x*x;};
      if ( DataPeriod == 2009 ) {
         W = pi09[0] + CUBE(pi09[1]/Pt);
      } else if ( DataPeriod == 2012 ) {
         W = pi12[0] + CUBE(pi12[1]/Pt);
      }
   }
   return W;
}
