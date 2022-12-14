// unbinned LH fit of M(K+K-) distributions
// -> Mkk_fit.pdf

#include "gsl/gsl_errno.h"   // GSL error handler

// {{{1 helper function
//----------------------------------------------------------------------
// adapter class to use lambda functions with closures in ROOT::MATH
// Using:  rmath_fun< decltype(Lambda) > Functor(Lambda);
//----------------------------------------------------------------------
template< typename F >
class rmath_fun: public ROOT::Math::IBaseFunctionOneDim {
//----------------------------------------------------------------------
   public:
      rmath_fun(const F& lam) : _f(lam) {}
      rmath_fun* Clone() const {
         return new rmath_fun(*this);
      }
   private:
      const F& _f;
      double DoEval(double x) const {
         return static_cast<double>(_f(x));
      }
};

//----------------------------------------------------------------------
constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

//--------------------------------------------------------------------
void SetHstFace(TH1* hst) {
//--------------------------------------------------------------------
   TAxis* X = hst->GetXaxis();
   if ( X ) {
      X->SetLabelFont(62);
      X->SetLabelSize(0.04);
      X->SetTitleFont(62);
      X->SetTitleSize(0.04);
   }
   TAxis* Y = hst->GetYaxis();
   if ( Y ) {
      Y->SetLabelFont(62);
      Y->SetLabelSize(0.04);
      Y->SetTitleFont(62);
      Y->SetTitleSize(0.04);
   }
   TAxis* Z = hst->GetZaxis();
   if ( Z ) {
      Z->SetLabelFont(62);
      Z->SetLabelSize(0.04);
      Z->SetTitleFont(62);
      Z->SetTitleSize(0.04);
   }
}

//----------------------------------------------------------------------
// adapter class for handling parameters and errors
//----------------------------------------------------------------------
struct ParAndErr {
   vector<double> Fpar; // parameters
   vector<double> Perr; // symmetric errors
   vector<double> Uerr; // upper minos errors
   vector<double> Lerr; // lower minos errors

   // ctor
   //-------------------------------------------------------------------
   ParAndErr(const ROOT::Fit::FitResult& res) {
   //-------------------------------------------------------------------
      Fpar = res.Parameters();
      Perr = res.Errors();

      int Npar = res.NPar();
      Lerr.resize(Npar,0.);
      Uerr.resize(Npar,0.);

      // if the upper and lower errors differ no more than the
      // 'threshold', the maximum of them is used as a symmetric error
      constexpr double threshold = 0.1;
      for ( int np = 0; np < Npar; ++np ) {
         if ( res.HasMinosError(np) ) {
            double uerr = fabs(res.UpperError(np));
            double lerr = fabs(res.LowerError(np));
            if ( uerr*lerr == 0 ) { // one of the errors is not defined
               continue;
            }
            Uerr[np] = uerr;
            Lerr[np] = lerr;
            if ( fabs(1-lerr/uerr) < threshold ) {
               Perr[np] = max(lerr,uerr);
            } else {
               Perr[np] = -1; // use asymmetric errors
               if ( res.IsParameterFixed(np) ) {
                  Perr[np] = 0;
               }
            }
         }
//          cout << " set Perr[" << np << "] = " << Perr[np] << endl;
      }
   }

   // pretty print parameters with errors
   //-------------------------------------------------------------------
   const char* Eform (int i,string fm, double X = 1) {
   //-------------------------------------------------------------------
      if ( Perr[i] > 0 ) {
         string fmt = "%" +fm+ " #pm %" +fm;
         return Form(fmt.c_str(),Fpar[i]*X,Perr[i]*X);
      } else if ( Perr[i] < 0 ) {
         string fmt = "%" + fm
            + "^{#lower[0.1]{#kern[0.1]{%+" + fm + "}}" + "}"
            + "_{#lower[-0.1]{#kern[0.2]{%+" + fm + "}}" + "}";
         return Form(fmt.c_str(),Fpar[i]*X,Uerr[i]*X,-Lerr[i]*X);
      }
      string fmt = "%" +fm+ " (fixed)";
      return Form(fmt.c_str(),Fpar[i]*X);
   }

   //-------------------------------------------------------------------
   const char* Tform (int i,string fm, double X = 1) {
   //-------------------------------------------------------------------
      if ( Perr[i] > 0 ) {
         string fmt = "%" +fm+ " \\pm %" +fm;
         return Form(fmt.c_str(),Fpar[i]*X,Perr[i]*X);
      } else if ( Perr[i] < 0 ) {
         string fmt = "%" + fm + "^{%+" + fm + "}" + "_{%+" + fm + "}";
         return Form(fmt.c_str(),Fpar[i]*X,Uerr[i]*X,-Lerr[i]*X);
      }
      string fmt = "%" +fm;
      return Form(fmt.c_str(),Fpar[i]*X);
   }
};

//----------------------------------------------------------------------
struct ValErr {
   double val;
   double err;
   const char* prt(string fm) {
      string fmt = "%" + fm + " +/- %" + fm;
      return Form(fmt.c_str(),val,err);
   }
};

// {{{1 Breigt Wigner for phi -> KK
//----------------------------------------------------------------------
// Breigt Wigner for e+ e- -> phi eta -> K+K- eta
//----------------------------------------------------------------------

//----------------------------------------------------------------------
double BreitWigner( double m, double Eee, double mphi, double gphi ) {
//----------------------------------------------------------------------
//    constexpr double Mjpsi = 3.096916; // 3096.916 +/- 0.011 MeV
   constexpr double Meta  = 0.547862; // 547.862  +/- 0.017 MeV
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV

//    constexpr double Mjpsi2 = SQ(Mjpsi);
   constexpr double Meta2  = SQ(Meta);
   constexpr double Mk2    = SQ(Mk);

   constexpr double dL = 2*Mk; // the left cutoff = 0.987354
   if ( m < dL ) { // phase space becomes zero (see p_k)
      return 0;
   }

   double m2    = SQ(m);
   double Eee2  = SQ(Eee);
   double mphi2 = SQ(mphi);
   double gphi2 = SQ(gphi);

   double gam = mphi*sqrt(mphi2+gphi2);
   double kap = 2*M_SQRT2*mphi*gphi*gam / (M_PI*sqrt(mphi2+gam));

   double p_phi0 = sqrt(SQ(Eee2-mphi2-Meta2)-4*mphi2*Meta2) / (2*Eee);
   double p_phi  = sqrt(SQ(Eee2-m2-Meta2)-4*m2*Meta2) / (2*Eee);
   double r_phi  = p_phi / p_phi0;

   double p_k0 = 0.5*sqrt(mphi2-4*Mk2);
   double p_k  = 0.5*sqrt(m2-4*Mk2);
   double r_k  = p_k / p_k0;

   // B(m)/B(m0)
   constexpr double R  = 3.; // GeV^{-1} for Blatt-Weisskopf ff
   double BB_phi = sqrt( (1+SQ(R*p_phi0)) / (1+SQ(R*p_phi)) );
   double BB_k2  = (1+SQ(R*p_k0)) / (1+SQ(R*p_k));
   double BB_k = sqrt( BB_k2 );

   double fm = r_phi*r_k*BB_phi*BB_k;

//?    double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k; // == m*G(m)
   double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k2; // == m*G(m)

   return (kap * SQ(fm)) / ( SQ(m2-mphi2) + SQ(GM) );
}

//----------------------------------------------------------------------
double BreitWignerGauss( double m, double Eee,
                         double mphi, double gphi, // BW parameters
                         double sigma              // Gauss
                       ) {
//----------------------------------------------------------------------
// Function Fun(x) folded with a normal distribution:
//
//        1                        [                  (X'-X)^2    ]
// ---------------- *   Integral   [ Fun(X') * exp( - --------- ) ] dX'
// sqrt(2*pi)*sigma   (-oo<X'<+oo) [                  2*sigma^2   ]
//
//----------------------------------------------------------------------
// sigma ~ 1/3 gphi => +/-5 sigma integration
//----------------------------------------------------------------------
   constexpr double one_over_sqrt2pi = 0.5*M_2_SQRTPI*M_SQRT1_2;

   // integrand lambda function
   double pp[] = { m, Eee, mphi, gphi, sigma };
   auto Lint = [](double t, void* pp) -> double{ // t == (X'-X)/sigma
      double* p = static_cast<double*>(pp);
      return exp(-0.5*t*t) * BreitWigner( (p[0]+t*p[4]),p[1],p[2],p[3]);
   };

   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);
   double result = gsl_int.Integral(Lint,pp,-5.,+5.);

   return one_over_sqrt2pi * result;
}

//----------------------------------------------------------------------
double BreitWignerGaussN( double m, double Eee,
                          double mphi, double gphi, double sigma
                        ) {
//----------------------------------------------------------------------
// Numeric normalisation on one for range [dL,dU]
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV
   constexpr double dL = 2*Mk; // the left cutoff = 0.987354
   constexpr double dU = 1.08; // upper limit

   double norm = 0;
   // cash parameters: Eee is the constant
   static double cacheN = 0;
   static double cacheM = 0;
   static double cacheG = 0;
   static double cacheS = 0;
   if ( cacheN > 0 &&
         mphi == cacheM && gphi == cacheG && sigma == cacheS ) {
      norm = cacheN;
   } else {
      // integrand lambda function
      double p[] = {Eee,mphi,gphi,sigma};
      auto Lint = [](double x, void* pp) -> double{
         double* p = static_cast<double*>(pp);
         return BreitWignerGauss(x,p[0],p[1],p[2],p[3]);
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-7, 1.e-6, 1000);
      norm = gsl_int.Integral(Lint,p,dL,dU);
      cacheN = norm;
      cacheM = mphi;
      cacheG = gphi;
      cacheS = sigma;
   }

   return BreitWignerGauss(m,Eee,mphi,gphi,sigma) / norm;
}

//-------------------------------------------------------------------------
void test_BreitWigner() {
//-------------------------------------------------------------------------
   constexpr double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
   constexpr double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV

   constexpr double dL = 0.98;
   constexpr double dU = 1.08;

   // 1) BreitWigner (internally normalized to 1)
   constexpr double Mjpsi = 3.096916; // 3096.916 +/- 0.011 MeV
//    double Eee = Mjpsi;
   double Eee = 3.08;
   auto Lbw = [Eee](const double* x,const double* p) -> double {
      return p[0]*BreitWigner(x[0],Eee,p[1],p[2]);
   };
   TF1* bw = new TF1("bw", Lbw, dL, dU, 3);
   bw -> SetParNames("Norm","Mphi","Gphi");
   bw -> SetParameters(1., Mphi, Gphi);
   bw -> SetLineWidth(1);
   bw -> SetLineColor(kBlue);
   bw -> SetNpx(500);

   double norm = 1./bw -> Integral(dL,dU,1e-8);
   printf("norm = %.7f\n",norm);
   bw -> SetParameter(0, norm );

   // 2) BreitWignerGaussN
   auto Lbwgn = [Eee](const double* x,const double* p) -> double {
      return p[0]*BreitWignerGaussN(x[0],Eee,p[1],p[2],p[3]);
   };
   TF1* bwgn = new TF1("bwgn", Lbwgn, dL, dU, 4);
   bwgn -> SetParNames("Norm","Mphi","Gphi","Sigma");
   bwgn -> SetParameters(1., Mphi, Gphi, 1.2e-3);
   bwgn -> SetLineWidth(2);
   bwgn -> SetLineColor(kRed);
   bwgn -> SetLineStyle(kDashed);
   bwgn -> SetNpx(500);

   double normN = 1./bwgn->Integral(dL,dU,1e-8);
   printf("normN = %.6f\n",normN);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   bw -> Draw();
   bwgn -> Draw("SAME");
}

// {{{1 Argus functions
//----------------------------------------------------------------------
// Argus functions
//----------------------------------------------------------------------

//----------------------------------------------------------------------
double Argus(double x, double a) {
//----------------------------------------------------------------------
// ARGUS distribution: https://en.wikipedia.org/wiki/ARGUS_distribution
// This is a regular ARGUS distribution with one parameter 'a'.

   if ( x < 0 || x > 1 ) {
      return 0;           // by definition
   }

   if ( fabs(a)<1e-5 ) {
      return 3*x*sqrt(1-x*x);
   }

   // normalization constant
   constexpr double one_over_sqrt2pi = 0.5*M_2_SQRTPI*M_SQRT1_2;
   double a2 = a*a;
   double Psi = 0.5*(1+erf(a*M_SQRT1_2))
      - a*one_over_sqrt2pi*exp(-0.5*a2)
      -0.5;
   double norm = a2*a*one_over_sqrt2pi/Psi;

   // regular Argus distribution:
   double u = 1 - x*x;
   double Ar = norm * x * sqrt(u) * exp(-0.5*a2*u);
   return Ar;
}

//----------------------------------------------------------------------
double RevArgus(double m, double A) {
//----------------------------------------------------------------------
// Argus function with left cutoff (L)
//      U is the upper limit of fitting
//      L,U - must be const for fit

   constexpr double Mjpsi = 3.096916; // 3096.916 +/- 0.011 MeV
   constexpr double Meta  = 0.547862; // 547.862  +/- 0.017 MeV
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV
   constexpr double L = 2*Mk;
   constexpr double U = Mjpsi - Meta;

   double x = (U-m)/(U-L); // transformation: (L,U) -> (1,0)
   return Argus(x,A)/(U-L);
}

//----------------------------------------------------------------------
double RevArgusN(double m, double A) {
//----------------------------------------------------------------------
// Numeric normalisation on one for range [dL,dU]
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV
   constexpr double dL = 2*Mk; // the left cutoff = 0.987354
   constexpr double dU = 1.08; // upper limit

   double norm = 0;
   // cash parameters
   static double cacheN = 0;
   static double cacheA = 0;
   if ( cacheN > 0 && A == cacheA ) {
      norm = cacheN;
   } else {
      // integrand lambda function
      double p[] = {A};
      auto Lint = [](double x, void* pp) -> double{
         double* p = static_cast<double*>(pp);
         return RevArgus(x,p[0]);
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-9, 1.e-6, 1000);
      norm = gsl_int.Integral(Lint,p,dL,dU);
      cacheN = norm;
      cacheA = A;
   }

   return RevArgus(m,A) / norm;
}

//-------------------------------------------------------------------------
void test_RevAgrus() {
//-------------------------------------------------------------------------
   constexpr double dL = 0.98;
   constexpr double dU = 1.08;

   // ---------------------------------------------------------------------
   // create Unuran 1D distribution object
   auto Lar = [](const double* x,const double* p) -> double {
      return RevArgus(x[0],p[0]);
   };
   TF1* far = new TF1("far", Lar, dL,dU, 1);
//    far -> SetParameter(0,0.98);
   far -> SetParameter(0,0.);
   TUnuranContDist dist(far);
   dist.SetDomain(dL,dU);
   // to use a different random number engine do:
   TRandom2* random = new TRandom2();
   int logLevel = 2;
   TUnuran unr(random,logLevel);
   // select unuran method for generating the random numbers
   string method = "tdr";
   if ( !unr.Init(dist,method) ) {
      cout << "Error initializing unuran" << endl;
      exit(0);
   }

   // ---------------------------------------------------------------------
   // generate
   TH1D* h1 = new TH1D("h1Ar","", 100,dL,dU);
   h1 -> SetLineWidth(2);
   h1 -> SetLineColor(kBlack);
   h1 -> SetMarkerStyle(20);

   int Ngen = 2000;
   for (int i = 0; i < Ngen; ++i) {
      double x = unr.Sample();
      h1 -> Fill(x);
   }

   // ---------------------------------------------------------------------
   // to draw
   double Wb = Ngen * (dU-dL) / h1 -> GetNbinsX();
   auto Larg = [Wb](const double* x,const double* p) -> double {
      return Wb*RevArgusN(x[0],p[0]);
   };

   TF1* fbkg = new TF1("fbkg", Larg, dL, dU, 1);
   fbkg -> SetParameter(0, far -> GetParameter(0) );
   fbkg -> SetLineWidth(2);
   fbkg -> SetLineColor(kBlue);


   TLegend* leg = new TLegend(0.13,0.74,0.43,0.89);
   leg -> SetHeader( Form("RevArgus a=%.2f",fbkg -> GetParameter(0)),"C" );

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   h1 -> Draw("EP");
   leg -> AddEntry( h1,"Unuran generated","PLE");

   fbkg -> DrawCopy("SAME");
   leg -> AddEntry( fbkg -> Clone(), "Function line","L");

   leg -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();

   // test normalization:
   printf(" RevArgusN norm= %.6f\n",fbkg -> Integral(dL,dU)/Wb );

}

// {{{1 Interference BW with Argus bkg.
//----------------------------------------------------------------------
// Interference BW with Argus bkg.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
vector<double> IntfrBWAR( double m, double Eee,
                  double mphi, double gphi,      // B-W
                  double A, double F, double Ang // Argus & interference
               ) {
//----------------------------------------------------------------------
// see BreitWigner() function
// The return vector contains:
//                   0    1     2        3
//        ret     { Sum, B-W, Argus, Interference }
//----------------------------------------------------------------------
//    constexpr double Mjpsi = 3.096916; // 3096.916 +/- 0.011 MeV
   constexpr double Meta  = 0.547862; // 547.862  +/- 0.017 MeV
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV

//    constexpr double Mjpsi2 = SQ(Mjpsi);
   constexpr double Meta2  = SQ(Meta);
   constexpr double Mk2    = SQ(Mk);

   constexpr double dL = 2*Mk; // the left cutoff = 0.987354
   if ( m < dL ) { // phase space becomes zero
      vector<double> ret(4,0.);
      return ret;
   }

   double m2    = SQ(m);
   double Eee2  = SQ(Eee);
   double mphi2 = SQ(mphi);
   double gphi2 = SQ(gphi);

   double gam = mphi*sqrt(mphi2+gphi2);
   double kap = 2*M_SQRT2*mphi*gphi*gam / (M_PI*sqrt(mphi2+gam));

   double p_phi0 = sqrt(SQ(Eee2-mphi2-Meta2)-4*mphi2*Meta2) / (2*Eee);
   double p_phi  = sqrt(SQ(Eee2-m2-Meta2)-4*m2*Meta2) / (2*Eee);
   double r_phi  = p_phi / p_phi0;

   double p_k0 = 0.5*sqrt(mphi2-4*Mk2);
   double p_k  = 0.5*sqrt(m2-4*Mk2);
   double r_k  = p_k / p_k0;

   // B(m)/B(m0)
   constexpr double R  = 3.; // GeV^{-1} for Blatt-Weisskopf ff
   double BB_phi = sqrt( (1+SQ(R*p_phi0)) / (1+SQ(R*p_phi)) );
   double BB_k2  = (1+SQ(R*p_k0)) / (1+SQ(R*p_k));
   double BB_k   = sqrt( BB_k2 );

   double fm = r_phi*r_k*BB_phi*BB_k;
//?    double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k; // == m*G(m)
   double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k2; // == m*G(m)

   // common multipliers
   double tmp =  fm / ( SQ(m2-mphi2) + SQ(GM) );

   // signal
   double BW = tmp*fm*kap;

   // background:
   double Ar = RevArgus(m,A)*SQ(F);

   // interference:
   double Intfr = 2 * tmp*sqrt(kap*Ar) *
                  ( (m2-mphi2)*cos(Ang) - GM*sin(Ang) );

   double Sum = BW + Ar + Intfr;
   vector<double> ret { Sum, BW, Ar, Intfr };

   return ret;
}

//----------------------------------------------------------------------
double IntfrBWARG( double m, double Eee,
                   const double p[], // {mphi,gphi,sigma,A,F,Ang}
                   unsigned int idx = 0 // what to return
                 ) {
//----------------------------------------------------------------------
// Function Fun(x) folded with a normal distribution:
//
//        1                        [                  (X'-X)^2    ]
// ---------------- *   Integral   [ Fun(X') * exp( - --------- ) ] dX'
// sqrt(2*pi)*sigma   (-oo<X'<+oo) [                  2*sigma^2   ]
//
//----------------------------------------------------------------------
// integration in +/-5 sigma interval
//----------------------------------------------------------------------
//    static const double Ngauss = 1./sqrt(2*M_PI);
   constexpr double one_over_sqrt2pi = 0.5*M_2_SQRTPI*M_SQRT1_2;

   // integrand lambda function
   auto Lint = [m,Eee,p,idx](double t) -> double{
      // t == (X'-X)/sigma
      double x = m+t*p[2];
      vector<double> Fun = IntfrBWAR(x, Eee, p[0],p[1],p[3],p[4],p[5]);
      return exp(-0.5*t*t) * Fun[idx];
   };
   rmath_fun< decltype(Lint) > Fint(Lint);

   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);
   double result = gsl_int.Integral(Fint,-5.,+5.);

   return one_over_sqrt2pi * result;
}

//----------------------------------------------------------------------
double IntfrBWARGN( double m, double Eee,
                    const double p[], //{mphi,gphi,sigma,A,F,Ang}
                    int idx = 0       // what to return
                  ) {
//----------------------------------------------------------------------
// Numerical normalisation to one for range [dL,dU]
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV
   constexpr double dL = 2*Mk; // the left cutoff = 0.987354
   constexpr double dU = 1.08; // upper limit

   // cash parameters: Eee is the constant

   static double cacheN = 0;
   constexpr int nP = 6;
   static double cP[nP]; // cache p[]

   double norm = cacheN;
   bool need_to_calc = !(cacheN > 0);
   if ( !need_to_calc ) {
      for (int i = 0; i < nP; ++i) {
         if ( p[i] != cP[i] ) { // are parameters the same?
            need_to_calc = true;
            break;
         }
      }
   }
   if ( need_to_calc ) {
      // integrand lambda function
      auto Lint = [Eee,p](double x) -> double {
         return IntfrBWARG(x,Eee,p,0);  // MUST be zero!
      };
      rmath_fun< decltype(Lint) > Fint(Lint);

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-5, 1.e-4, 1000);
      norm = gsl_int.Integral(Fint,dL,dU);

      // save to cache
      cacheN = norm;
      for (int i = 0; i < nP; ++i) {
         cP[i] = p[i];
      }
   }

   return IntfrBWARG(m,Eee,p,idx) / norm;
}

//-------------------------------------------------------------------------
void test_Intfr() {
//-------------------------------------------------------------------------
   constexpr double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
   constexpr double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV

   constexpr double dL = 0.98;
   constexpr double dU = 1.08;

   double Eee = 3.097; // ~ M(J/Psi)
   // ---------------------------------------------------------------------
   // create Unuran 1D distribution object
   auto Lgen = [Eee](const double* x,const double* p) -> double {
      vector<double> res = IntfrBWAR( x[0], Eee, p[0],p[1],p[2],p[3],p[4] );
      return res[0];
   };
   TF1* fgen = new TF1("fgen", Lgen, dL,dU, 5);
   double ag = 0, Fg = 1.0, Ang = 0.8;
   fgen -> SetParameters(Mphi,Gphi,ag,Fg,Ang);

   double Xmax = fgen -> GetMaximumX(Mphi-0.4e-3,Mphi+0.4e-3);
   printf(" Xmax= %.7f (delta= %.1e)\n",Xmax,Xmax-Mphi);

   TUnuranContDist dist(fgen);
   dist.SetDomain(dL,dU);
   dist.SetMode(Xmax);
   // to use a different random number engine do:
   TRandom2* random = new TRandom2();
   int logLevel = 2;
   TUnuran unr(random,logLevel);
   // select unuran method for generating the random numbers
   string method = "method=nrou";
   cout << " start Unuran initializing with method: " << method << endl;
   if ( !unr.Init(dist,method) ) {
      cout << "Error initializing unuran" << endl;
      exit(0);
   }

   // ---------------------------------------------------------------------
   // generate
   TH1D* h1 = new TH1D("h1Intfr","", 100,dL,dU);
   h1 -> SetLineWidth(2);
   h1 -> SetLineColor(kBlack);
   h1 -> SetMarkerStyle(20);

   int Ngen = 3000;
   for (int i = 0; i < Ngen; ++i) {
      double x = unr.Sample();
      h1 -> Fill(x);
   }

   // ---------------------------------------------------------------------
   // to draw
   double Wb = Ngen * (dU-dL) / h1 -> GetNbinsX();
   double mphi = Mphi, gphi = Gphi;
   auto Lintfr = [Eee,Wb,mphi,gphi,ag,Fg,Ang](double* x, double* p)
      -> double {
      vector<double> res = IntfrBWAR( x[0], Eee, mphi,gphi,ag,Fg,Ang );
      return Wb*res[int(p[0])];
   };

   TF1* fun = new TF1("fun", Lintfr, dL, dU, 1);
   fun -> SetLineWidth(2);
   fun -> SetNpx(500);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   TLegend* leg = new TLegend(0.59,0.59,0.89,0.89);
   leg -> SetHeader("Interference BW #otimes Argus","C");

   h1 -> SetMinimum(-50);
   h1 -> Draw("EP");
   leg -> AddEntry( h1,"Unuran generated","PLE");

   fun -> SetParameter(0, 0.); // draw all
   fun -> SetLineColor(kRed);
   fun -> DrawCopy("SAME");
   leg -> AddEntry( fun -> Clone(), "Function line", "L");

   fun -> SetParameter(0, 1); // draw BW
   fun -> SetLineColor(kGreen+2);
   fun -> DrawCopy("SAME");
   leg -> AddEntry( fun -> Clone(), "Breit-Wigner", "L");

   fun -> SetParameter(0, 2); // draw Argus
   fun -> SetLineColor(kBlue);
   fun -> SetLineStyle(kDashed);
   fun -> DrawCopy("SAME");
   leg -> AddEntry( fun -> Clone(), "Argus", "L");

   fun -> SetParameter(0, 3); // draw interference
   fun -> SetLineColor(kMagenta);
   fun -> SetLineStyle(kSolid);
   fun -> DrawCopy("SAME");
   leg -> AddEntry( fun -> Clone(), "Interference", "L");

   leg -> Draw();
   c1->Update();
}

// {{{1  Data
//-------------------------------------------------------------------------
 vector<double> get_mkk_hist( string fname, string hname,
                              TH1D* hst[], int type = 0 ) {
//-------------------------------------------------------------------------
#include "cuts.h"

   constexpr double dL = 0.98;
   constexpr double dU = 1.08;
   constexpr int Nbins = 50;
   constexpr double bW = (dU-dL)/Nbins; // bin width

   // name of folder with root files
   static string dir("Ntpls/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot -> cd("SelectKKgg");
   TTree* a4c = (TTree*)gDirectory -> Get("a4c");

   TCut c_here = c_chi2;
   c_here += TCut(Form("%f<=Mkk&&Mkk<%f",dL,dU));
   if ( type == 0 ) {        // central part
      c_here += c_cpgg;
   } else if ( type == 1 ) { // side-band
      c_here += c_sbgg;
   }
   int n = a4c -> Draw("Mkk", c_here, "goff");
//    cout << " n= " << n << endl;
   double* buffer = a4c -> GetVal(0);
   vector<double> mkk(buffer, buffer+n);

   string title(";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2}");
   int ibW = int(bW*1e5);
   title += string(";Entries / ") +
      ( (ibW%10 == 0) ? string(Form("%.0f MeV/c^{2}",bW*1e3))
                      : string(Form("%.2f MeV/c^{2}",bW*1e3)) );
   hst[0] = new TH1D(hname.c_str(),title.c_str(),Nbins,dL,dU);

   hst[0] -> Sumw2(true);
   for ( const auto& mk : mkk ) {
      hst[0] -> Fill(mk);
   }
   return mkk;
}

//-------------------------------------------------------------------------
TH1D* get_mkk_hist( string fname, string hname, int type = 0 ) {
//-------------------------------------------------------------------------
   TH1D* hist[1];
   get_mkk_hist( fname, hname, hist, type );
   return hist[0];
}

// {{{1  Fit BW + flat bkg
//-------------------------------------------------------------------------
void do_fit(string fname, double Eee, string title, string pdf="") {
//-------------------------------------------------------------------------
   constexpr double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
   constexpr double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV

   constexpr double dL = 0.98;
   constexpr double dU = 1.08;
   const double Wul = 1./(dU-dL);

   // Get un-binned data
   TH1D* hist[1];
   vector<double> mkk = get_mkk_hist( fname, "mkk", hist );
   TH1D* hst = hist[0];
   int Nb = hst -> GetNbinsX();
   double bW = (dU-dL)/Nb; // bin width

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Sig(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Sig.Add(mkk[i]);
   }
//    cout << " n= " << n << " Integral= " << hist[0]->Integral() << endl;

   //-----------------------------------------------------------------------
   // Fit signal MC(phi,eta) by Breit-Wigner convoluted with Gauss
   // function represents the total number of events (extended LH)
   auto Lfit = [Eee,Wul](const double* x,const double* p) -> double {
      return p[3] * BreitWignerGaussN(x[0],Eee,p[0],p[1],p[2])
             + p[4] * Wul;
   };

   vector<string> par_name { "Mphi", "Gphi", "Sigma", "Nphi", "Nbg" };
   vector<double> par_ini  {  Mphi,   Gphi,   1.0e-3, 0.99*n, 0.01*n };

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 Wfit( *Ffit );
   fitter.SetFunction( Wfit, false); // false == no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
//    fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   fitter.Config().ParSettings(0).SetValue(1.01952); // like MC-sig
   fitter.Config().ParSettings(0).Fix();

//    fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3, Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();

//    fitter.Config().ParSettings(2).SetValue(1e-3);
   fitter.Config().ParSettings(2).SetLimits(0.1e-3, 3e-3); // Sigma

   fitter.Config().ParSettings(3).SetLimits(0.8*n,1.2*n); // Nphi
   fitter.Config().ParSettings(4).SetLimits(0.,0.2*n);    // Nbg

   fitter.LikelihoodFit( Sig, true ); // true == extended likelihood fit
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof = [Lfit,Fpar](double x) -> double {
      return Lfit(&x,Fpar.data());
   };
   rmath_fun< decltype(Lgof) > ftest(Lgof);
   ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(
         mkk.size(),mkk.data(),ftest,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS = goftest -> KolmogorovSmirnovTest();
   cout << " pvalueKS= " << pvalueKS << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   auto Ldr = [Lfit,Fpar,bW](double* x,double* p) -> double {
      return bW * Lfit(x,Fpar.data());
   };
   TF1* fdr = new TF1("fdr", Ldr, dL, dU, 0);
   fdr -> SetLineWidth(2);
   fdr -> SetLineColor(kRed);
   fdr -> SetNpx(500);

   auto Lbkg = [Fpar,Wul,bW](double* x,double* p) -> double {
      return bW * Fpar[4] * Wul;
   };
   TF1* fbkg = new TF1("fbkg", Lbkg, dL, dU, 0);
   fbkg -> SetLineWidth(2);
   fbkg -> SetLineColor(kBlue);
   fbkg -> SetLineStyle(kDashed);

   //-----------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();
//    gPad -> SetLogy();

   SetHstFace(hst);
   hst -> GetYaxis() -> SetMaxDigits(3);
   hst -> GetYaxis() -> SetTitleOffset(1.2);
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
   hst -> SetMarkerStyle(20);

   hst -> Draw("EP");
   fdr -> Draw("SAME");

   TLegend* leg = new TLegend(0.60,0.79,0.89,0.89);
   leg -> AddEntry(hst,title.c_str(),"LEP");
   leg -> AddEntry(fdr,"BW #otimes Gauss","L");
   if ( Fpar[4] > 1 ) { // Nbg > 1
      fbkg -> Draw("SAME");
      leg -> AddEntry(fbkg,"flat background","L");
   }
   leg -> Draw();

   TPaveText* pt = new TPaveText(0.60,0.57,0.89,0.79,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
   pt -> AddText( Form("#it{p-value(K-S) = %.3f}",pvalueKS) );
   pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt -> AddText( Form("#sigma= %s MeV", PE.Eform(2,".2f",1e3)) );
   pt -> AddText( Form("N_{#phi}= %s",PE.Eform(3,".1f")) );
   pt -> AddText( Form("N_{bkg}= %s",PE.Eform(4,".1f")) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      pdf = string("mkk_pict/") + pdf + ".pdf";
      c1 -> Print(pdf.c_str());
   }

   // print for table
   printf("%s\n no interference: p(KS)= %.3f\n"
          " sigma= %s MeV  Nphi= %s  Nbg= %s\n",
         title.c_str(), pvalueKS,
         PE.Tform(2,".2f",1e3),
         PE.Tform(3,".1f"), PE.Tform(4,".1f") );
}

// {{{1  Fit Interference(BW + RevArgus)
//----------------------------------------------------------------------
ValErr CalcIntEr( const ROOT::Fit::FitResult& res,
                  double Eee, unsigned int idx ) {
//----------------------------------------------------------------------
// Numerical calculation of integrals and their errors
   // desired errors:
   constexpr double dL = 0.98;
   constexpr double dU = 1.08;
   constexpr double eps_abs = 1.e-6;
   constexpr double eps_rel = 1.e-6;
   constexpr double epsilon = 1.e-3; // max numerical error

   // get parameters and covariation matrix of error
   int Npar = res.NPar();
   vector<double> Fpar = res.Parameters();

   int covMatrStatus = res.CovMatrixStatus();
   if ( covMatrStatus != 3 ) {
      cout << " WARNING: " << __func__ << " covariance matrix"
         " status code is " << covMatrStatus << endl;
   }
   vector<double> cov_m(Npar*Npar,0);
   for ( int i = 0; i < Npar; ++i ) {
      for ( int j = 0; j < Npar; ++j ) {
         cov_m[i*Npar+j] = res.CovMatrix(i,j);
      }
   }

   // 1) calculate integral of normalized component
   auto Lfc = [Eee,idx](const double* x,const double* p) -> double {
      return (1-p[6])*IntfrBWARGN( x[0], Eee, p, idx );
   };
   TF1 FC("FC", Lfc, dL, dU, Npar);
   FC.SetParameters(Fpar.data());

   double err_num = 0;
   double Integ = FC.IntegralOneDim(dL,dU,eps_rel,eps_abs, err_num);
//    cout << " err_num(" << idx << ") = " << err_num << endl;
   if ( err_num > epsilon ) {
      cout << " WARNING: " << __func__ << " numerical error of"
         " integration of F_C(" << idx << ") is too big "
         << err_num << endl;
   }

   // 2) calculate error of this integral
   TMatrixDSym covMatrix(Npar);
   covMatrix.Use(Npar,cov_m.data());
//    printf("\ncovMatrix-%d : ",idat);
//    covMatrix.Print();

   TVectorD IntGrad(Npar);
   double err_num2 = 0;
   for ( int i = 0; i < Npar; ++i ) {
      // skip parameters with zero error
      if ( covMatrix(i,i) == 0 ) {
         continue;
      }

      auto LdF = [FC,i](const double* x, const double* p) -> double {
         TF1 tmp(FC); // may modify 'tmp' but not FC
         return tmp.GradientPar(i,x);
      };
      TF1 dF("dF",LdF,dL,dU,0);
      double err_num = 0;
      IntGrad[i] = dF.IntegralOneDim(dL,dU,eps_rel,eps_abs, err_num);
      err_num = covMatrix(i,i) * IntGrad[i] * err_num;
//       cout << " abs err_num(" << i << ") = " << err_num << endl;
      err_num2 += SQ(err_num);
//       cout << " IntGrad[" << i << "] = " << IntGrad[i] << endl;
   }

   double err_Int = sqrt( covMatrix.Similarity(IntGrad) );
   err_num2 = sqrt( err_num2 / err_Int ); // abs numerical error
//    cout << " err_num(" << idx << ") = " << err_num2 << endl;
   if ( err_num2 > epsilon ) {
      cout << " WARNING: " << __func__ << " numerical error "
         "of integration of dF is too big " << err_num2 << endl;
   }

   return (ValErr {Integ,err_Int});
}

//-------------------------------------------------------------------------
void do_fitI(string fname, double Eee, string title, string pdf="") {
//-------------------------------------------------------------------------
   constexpr double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
   constexpr double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV

   constexpr double dL = 0.98;
   constexpr double dU = 1.08;
   const double Wul = 1./(dU-dL);

   // Get un-binned data
   TH1D* hist[1];
   vector<double> mkk = get_mkk_hist( fname, "mkk", hist );
   TH1D* hst = hist[0];
   int Nb = hst -> GetNbinsX();
   double bW = (dU-dL)/Nb; // bin width

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Sig(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Sig.Add(mkk[i]);
   }
   cout << " n= " << n << " Integral= " << hist[0]->Integral() << endl;

   //-----------------------------------------------------------------------
   // Fit signal MC(phi,eta) by Breit-Wigner convoluted with Gauss
   // function MUST be normalized to 1 on the fit range
   auto Lfit = [Eee,Wul](const double* x,const double* p) -> double {
      return (1-p[6])*IntfrBWARGN( x[0], Eee, p, 0 ) + p[6]*Wul;
   };

   vector<string> par_name { "Mphi", "Gphi", "Sigma",
                             "A", "F", "vartheta", "Bkg" };
//    vector<double> par_ini { Mphi,  Gphi, 1.9e-3,
//                              0.0, 0.4, -0.81, 0.01 }; //p
   vector<double> par_ini {  Mphi,  Gphi, 2.0e-3,
                             0.0, 0.2, 0.76, 0.01 }; //n

   bool neg_pos = par_ini[5] > 0; // true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 Wfit( *Ffit );
   fitter.SetFunction( Wfit, false); // false == no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
//    fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   fitter.Config().ParSettings(0).SetValue(1.01952); // like MC-sig
   fitter.Config().ParSettings(0).Fix();
//    fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3, Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();
   fitter.Config().ParSettings(2).SetLimits(0.1e-3, 3e-3);  // Sigma
   fitter.Config().ParSettings(3).Fix();                    // A
   fitter.Config().ParSettings(4).SetLimits(0., 10.);       // F
   fitter.Config().ParSettings(5).SetLimits(-M_PI, M_PI);   // vartheta
   if ( neg_pos ) {
      fitter.Config().ParSettings(5).SetValue(0.76);        // n-MEMO
   } else {
      fitter.Config().ParSettings(5).SetValue(-0.81);       // p-MEMO
   }
   fitter.Config().ParSettings(5).Fix();
   fitter.Config().ParSettings(6).SetLimits(0., 0.5);       // Bkg
//    fitter.Config().ParSettings(6).SetValue(0.);
//    fitter.Config().ParSettings(6).Fix();

   fitter.LikelihoodFit( Sig, false ); // true == extended likelihood fit
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   ValErr Integ = CalcIntEr( res, Eee, 1 ); // BW
   double Nphi = n * Integ.val;
   double err_Nphi = fabs(Nphi) * sqrt( 1./n + SQ(Integ.err/Integ.val) );
   printf("Nphi = %.1f +/- %.1f\n",Nphi,err_Nphi);
   double Nbg = Fpar[6]*n;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof = [Lfit,Fpar](double x) -> double {
      return Lfit(&x,Fpar.data());
   };
   rmath_fun< decltype(Lgof) > ftest(Lgof);
   ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(
         mkk.size(),mkk.data(),ftest,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS = goftest -> KolmogorovSmirnovTest();
   cout << " pvalueKS= " << pvalueKS << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   double bWn = bW*n;
   auto Ldr = [bWn,Eee,Wul,Fpar](double* x,double* p) -> double {
      int idx = int(p[0]);
      return ( idx == 0 )
         ? bWn*( (1-Fpar[6]) * IntfrBWARGN(x[0],Eee,Fpar.data(),0)
               + Fpar[6] * Wul)
         : bWn*(1-Fpar[6])*IntfrBWARGN(x[0],Eee,Fpar.data(),idx);
   };
   TF1* fdr = new TF1("fdr", Ldr, dL, dU, 1);
   fdr -> SetLineWidth(2);
   fdr -> SetLineColor(kRed);
   fdr -> SetNpx(500);

   // find min to draw
   double hmax = 0, hmin = 0;
   if ( neg_pos ) {
      double x  = Mphi, inx=1;
      hmax = int(1.15*Ldr(&x,&inx))+1;
      inx=3;
      hmin = int(1.1*Ldr(&x,&inx))-2;
   } else {
      hmin = min(-1.,-hst -> GetMaximum() / 25.);
   }

   auto Lbkg = [Fpar,bWn,Wul](double* x,double* p) -> double {
      return bWn * Fpar[3] * Wul;
   };
   TF1* fbkg = new TF1("fbkg", Lbkg, dL, dU, 0);
   fbkg -> SetLineWidth(1);
   fbkg -> SetLineColor(kBlue+3);
   fbkg -> SetLineStyle(kDashed);

   //-----------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   SetHstFace(hst);
   if ( hmax > 0 ) {
      hst -> SetMaximum(hmax);
   }
   if ( hmin < 0 ) {
      hst -> SetMinimum(hmin);
   }
   hst -> GetYaxis() -> SetMaxDigits(3);
   hst -> GetYaxis() -> SetTitleOffset(1.2);
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
   hst -> SetMarkerStyle(20);

   hst -> Draw("EP");

   TLegend* leg = new TLegend(0.60,0.74,0.89,0.89);
   leg -> AddEntry(hst,title.c_str(),"LEP");

   fdr -> SetParameter(0, 0); // SUM
   fdr -> DrawCopy("SAME");
   leg -> AddEntry( fdr -> Clone(), "Result of the fit", "L");

   // draw components
   fdr -> SetLineWidth(1);
   fdr -> SetLineStyle(kDashed);

   fdr -> SetParameter(0, 1);   // BW
   fdr -> SetLineColor(kGreen+2);
   fdr -> DrawCopy("SAME");
   leg -> AddEntry( fdr -> Clone(), "Breit-Wigner", "L");

   fdr -> SetParameter(0, 2);   // Argus
   fdr -> SetLineColor(kBlue);
   fdr -> DrawCopy("SAME");
   leg -> AddEntry( fdr -> Clone(), "Argus", "L");

   fdr -> SetParameter(0, 3);   // interference
   fdr -> SetLineColor(kMagenta+1);
   fdr -> DrawCopy("SAME");
   leg -> AddEntry( fdr -> Clone(), "Interference", "L");

   if ( Nbg > 1 ) { // Bkg
      fbkg -> Draw("SAME");
      leg -> AddEntry(fbkg,"flat background","L");
   }
   leg -> Draw();

   TPaveText* pt = new TPaveText(0.60,0.50,0.89,0.74,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
   pt -> AddText( Form("#it{p-value(K-S) = %.3f}",pvalueKS) );
   pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt -> AddText( Form("#sigma= %s MeV", PE.Eform(2,".2f",1e3)) );
   pt -> AddText( Form("a= %s",PE.Eform(3,".1f")) );
   pt -> AddText( Form("F= %s",PE.Eform(4,".2f")) );
   pt -> AddText( Form("#vartheta= %s",PE.Eform(5,".2f")) );
//    if ( fabs(PE.Perr[6]) > 1e-3 ) {
//       pt -> AddText( Form("Bkg= %s", PE.Eform(6,".3f")) );
//    }
   pt -> AddText( Form("N_{#phi}= %.1f #pm %.1f",Nphi,err_Nphi) );
   pt -> AddText( Form("N_{bg}= %s", PE.Eform(6,".1f",n)) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      pdf = string("mkk_pict/") + pdf +
         (neg_pos ? "_n.pdf" : "_p.pdf");
      c1 -> Print(pdf.c_str());
   }

   // print for table
   printf("%s\n Interference:  vartheta= %.2f  p(KS)= %.3f\n"
          " sigma= %s MeV  F= %s  Nphi= %.1f \\pm %.1f  Nbg= %s\n",
         title.c_str(), Fpar[5], pvalueKS,
         PE.Tform(2,".2f",1e3), PE.Tform(4,".2f"),
         Nphi,err_Nphi, PE.Tform(6,".1f",n));
}

// {{{1 MAIN:
//-------------------------------------------------------------------------
void mass_KK_fit() {
//-------------------------------------------------------------------------
   gROOT -> Reset();
   gStyle -> SetOptStat(0);
   gStyle -> SetOptFit(111); // do not print fixed parameters!
//    gStyle -> SetOptFit(112); // print all parameters (fixed)
   gStyle -> SetStatFont(62);
   gStyle -> SetLegendFont(42);

//    gStyle->SetStatH(0.25);
//    gStyle->SetStatW(0.23);

   //-------------------------------------------------------------------
   // = define GSL error handler which does nothing =
   gsl_set_error_handler_off();

   //--------------------------------------------
//    test_BreitWigner();
//    test_RevAgrus();
//    test_Intfr();

   //--------------------------------------------
   // Fit: BW + flat background
   //--------------------------------------------
//    do_fit( "ntpl_2900_rs.root", 2.90,
//            "R-scan: 2900MeV","mkk_2900_rs");

//    do_fit( "ntpl_3080_rs.root", 3.08,
//            "R-scan: 3080MeV","mkk_3080_rs");

//    do_fit( "ntpl_3080_2019.root", 3.08,
//            "2019: 3080MeV","mkk_3080_2019");

//    do_fit( "ntpl_J1.root", 3.087659,
//            "2018: 3087.659MeV","mkk_3087.659" );

//    do_fit( "ntpl_3090.root", 3.088868,
//            "2012: 3088.868MeV","mkk_3088.868" );

//    do_fit( "ntpl_3093.root", 3.091774,
//            "2012: 3091.774MeV","mkk_3091.774" );

//    do_fit( "ntpl_3094.root", 3.094711,
//            "2012: 3094.711MeV","mkk_3094.711" );

//    do_fit( "ntpl_3095.root", 3.095444,
//            "2012: 3095.444MeV","mkk_3095.444" );

//    do_fit( "ntpl_J2.root", 3.095726,
//            "2018: 3095.726MeV","mkk_3095.726" );

//    do_fit( "ntpl_3096.root", 3.095840,
//            "2012: 3095.840MeV","mkk_3095.840" );

//    do_fit( "ntpl_J3.root", 3.096203,
//            "2018: 3096.203MeV","mkk_3096.203" );

//    do_fit( "ntpl_J4.root", 3.096986,
//            "2018: 3096.986MeV","mkk_3096.986" );

//    do_fit( "ntpl_J5.root", 3.097226,
//            "2018: 3097.226MeV","mkk_3097.226");

//    do_fit( "ntpl_3097.root", 3.097227,
//            "2012: 3097.227MeV","mkk_3097.227" );

//    do_fit( "ntpl_J6.root", 3.097654,
//            "2018: 3097.654MeV","mkk_3097.654");

//    do_fit( "ntpl_3098.root", 3.098354,
//            "2012: 3098.354MeV","mkk_3098.354" );

//    do_fit( "ntpl_J7.root", 3.098728,
//            "2018: 3098.728MeV","mkk_3098.728" );

//    do_fit( "ntpl_3099.root", 3.099056,
//            "2012: 3099.056MeV","mkk_3099.056" );

//    do_fit( "ntpl_J8.root", 3.104,
//            "2018: 3104.000MeV","mkk_3104.000" );

//    do_fit( "ntpl_3112.root", 3.112065,
//            "2012: 3112.065MeV","mkk_3112.065" );

   //--------------------------------------------
   // Fit: Interference(BW + RevArgus)
   //--------------------------------------------
//    do_fitI( "ntpl_2900_rs.root", 2.90,
//             "R-scan: 2900MeV","mkkI_2900_rsF");

//    do_fitI( "ntpl_3080_rs.root", 3.08,
//             "R-scan: 3080MeV","mkkI_3080_rsF");

//    do_fitI( "ntpl_3080_2019.root", 3.08,
//             "2019: 3080MeV","mkkI_3080_2019F");

   do_fitI( "ntpl_3090.root", 3.088868,
            "2012: 3088.868MeV","mkkI_3088.868" );

//    do_fitI( "ntpl_3093.root", 3.091774,
//             "2012: 3091.774MeV","mkkI_3091.774" );

//    do_fitI( "ntpl_3094.root", 3.094711,
//             "2012: 3094.711MeV","mkkI_3094.711" );

//    do_fitI( "ntpl_3095.root", 3.095444,
//            "2012: 3095.444MeV","mkkI_3095.444" );

//    do_fitI( "ntpl_J2.root", 3.095726,
//            "2018: 3095.726MeV","mkkI_3095.726" );

//    do_fitI( "ntpl_3096.root", 3.095840,
//            "2012: 3095.840MeV","mkkI_3095.840" );

//    do_fitI( "ntpl_J3.root", 3.096203,
//             "2018: 3096.203MeV","mkkI_3096.203" );

//    do_fitI( "ntpl_J4.root", 3.096986,
//             "2018: 3096.986MeV","mkkI_3096.986" );

//    do_fitI( "ntpl_J5.root", 3.097226,
//             "2018: 3097.226MeV","mkkI_3097.226");

//    do_fitI( "ntpl_3097.root", 3.097227,
//             "2012: 3097.227MeV","mkkI_3097.227" );

//    do_fitI( "ntpl_J6.root", 3.097654,
//             "2018: 3097.654MeV","mkkI_3097.654");

//    do_fitI( "ntpl_3098.root", 3.098354,
//             "2012: 3098.354MeV","mkkI_3098.354" );

//    do_fitI( "ntpl_J7.root", 3.098728,
//             "2018: 3098.728MeV","mkkI_3098.728" );

//    do_fitI( "ntpl_3099.root", 3.099056,
//             "2012: 3099.056MeV","mkkI_3099.056" );

//    do_fitI( "ntpl_J8.root", 3.104,
//             "2018: 3104.000MeV","mkkI_3104.000" );

//    do_fitI( "ntpl_3112.root", 3.112065,
//             "2012: 3112.065MeV","mkkI_3112.065" );
}
