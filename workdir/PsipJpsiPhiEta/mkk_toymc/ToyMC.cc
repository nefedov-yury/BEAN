// toy MC for J/Psi -> K+ K- gamma gamma
// This is ROOT version
//
#include "gsl/gsl_errno.h"   // GSL error handler

// {{{1 helper functions
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

//----------------------------------------------------------------------
// adapter class for handling parameters and errors
//----------------------------------------------------------------------
struct ParAndErr {
//----------------------------------------------------------------------
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
         string fmt = "%" +fm +"^{%+" +fm+ "}_{%+" +fm+ "}";
         return Form(fmt.c_str(),Fpar[i]*X,Uerr[i]*X,-Lerr[i]*X);
      }
      string fmt = "%" +fm+ " (fixed)";
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

//----------------------------------------------------------------------
void SetHstFace(TH1* hst) {
//----------------------------------------------------------------------
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

// {{{1 Breigt Wigner for phi -> KK
//----------------------------------------------------------------------
// Breigt Wigner for phi -> KK
//----------------------------------------------------------------------

//----------------------------------------------------------------------
double BreitWigner(double m, double mphi, double gphi) {
//----------------------------------------------------------------------
   constexpr double Mjpsi = 3.096916; // 3096.916 +/- 0.011 MeV
   constexpr double Meta  = 0.547862; // 547.862  +/- 0.017 MeV
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV

   constexpr double Mjpsi2 = SQ(Mjpsi);
   constexpr double Meta2  = SQ(Meta);
   constexpr double Mk2    = SQ(Mk);

   constexpr double dL = 2*Mk; // the left cutoff = 0.987354
   if ( m < dL ) { // phase space becomes zero (see p_k)
      return 0;
   }

   double m2    = SQ(m);
   double mphi2 = SQ(mphi);
   double gphi2 = SQ(gphi);

   double gam = mphi*sqrt(mphi2+gphi2);
   double kap = 2*M_SQRT2*mphi*gphi*gam / (M_PI*sqrt(mphi2+gam));

   double p_phi0 = sqrt(SQ(Mjpsi2-mphi2-Meta2)-4*mphi2*Meta2) / (2*Mjpsi);
   double p_phi  = sqrt(SQ(Mjpsi2-m2-Meta2)-4*m2*Meta2) / (2*Mjpsi);
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

//    double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k; // == m*G(m)
   double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k2; // == m*G(m)

   return (kap * SQ(fm)) / ( SQ(m2-mphi2) + SQ(GM) );
}

//----------------------------------------------------------------------
double BreitWignerGauss( double m,
                         double mphi, double gphi, // BW parameters
                         double sigma,             // Gauss
                         double slope              // eff(m)
                       ) {
//----------------------------------------------------------------------
// Function Fun(x) folded with a normal distribution:
//
//        1                        [                  (X'-X)^2    ]
// ---------------- *   Integral   [ Fun(X') * exp( - --------- ) ] dX'
// sqrt(2*pi)*sigma   (-oo<X'<+oo) [                  2*sigma^2   ]
//
//----------------------------------------------------------------------
// sigma ~ 1/3 gphi (bes3-exp)
// integration in +/-5 sigma interval
//----------------------------------------------------------------------
   static const double Ngauss = 1./sqrt(2*M_PI);

   // integrand lambda function
   double pp[] { m, mphi, gphi, sigma, slope };
   auto Lint = [](double t, void* pp) -> double{ // t == (X'-X)/sigma
      double* p = static_cast<double*>(pp);
      double x = p[0] + t*p[3];
      return exp(-0.5*t*t) *
         BreitWigner(x,p[1],p[2]) * (1+p[4]*(x-1.02));
   };

   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-7, 1.e-6, 1000);
   double result = gsl_int.Integral(Lint,pp,-5.,+5.);

   return Ngauss * result;
}

//----------------------------------------------------------------------
double BreitWignerGaussN( double m,
                          double mphi, double gphi, double sigma,
                          double slope
                        ) {
//----------------------------------------------------------------------
// Numeric normalisation on one for range [dL,dU]
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV
   constexpr double dL = 2*Mk; // the left cutoff = 0.987354
   constexpr double dU = 1.08; // upper limit

   double norm = 0;
   // cash parameters
   static double cacheN = 0;
   static double cacheM = 0;
   static double cacheG = 0;
   static double cacheS = 0;
   static double cacheSl = 0;
   if ( cacheN > 0 && mphi == cacheM && gphi == cacheG
         && sigma == cacheS && slope == cacheSl ) {
      norm = cacheN;
   } else {
      // integrand lambda function
      double p[] = {mphi,gphi,sigma,slope};
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
      cacheSl= slope;
   }

   return BreitWignerGauss(m,mphi,gphi,sigma,slope) / norm;
}

//-------------------------------------------------------------------------
void test_BreitWigner() {
//-------------------------------------------------------------------------
   constexpr double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
   constexpr double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV

   constexpr double dL = 0.98;
   constexpr double dU = 1.08;

   // 1) BreitWigner (internally normalized to 1)
   auto Lbw = [](const double* x,const double* p) -> double {
      return p[0]*BreitWigner(x[0],p[1],p[2]);
   };
   TF1* bw = new TF1("bw", Lbw, dL, dU, 3);
   bw -> SetParNames("Norm","Mphi","Gphi");
   bw -> SetParameters(1., Mphi, Gphi);
   bw -> SetLineWidth(1);
   bw -> SetLineColor(kBlue);
   bw -> SetNpx(500);

   double norm = 1./bw -> Integral(dL,dU,1e-7);
   printf("norm = %.7f\n",norm);
   bw -> SetParameter(0, norm );

   // 2) BreitWignerGaussN
   auto Lbwgn = [](const double* x,const double* p) -> double {
      return p[0]*BreitWignerGaussN(x[0],p[1],p[2],p[3],p[4]);
   };
   TF1* bwgn = new TF1("bwgn", Lbwgn, dL, dU, 5);
   bwgn -> SetParNames("Norm","Mphi","Gphi","Sigma","Slope");
   bwgn -> SetParameters(1., Mphi, Gphi, 1.2e-3, 0.);
   bwgn -> SetLineWidth(2);
   bwgn -> SetLineColor(kRed);
   bwgn -> SetLineStyle(kDashed);
   bwgn -> SetNpx(500);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg -> SetHeader("Breit-Wigner (M_{#phi}, #Gamma_{#phi})","C");
   leg -> AddEntry(bw, "BW", "L");

   bw -> Draw();
   bwgn -> DrawCopy("SAME");
   leg -> AddEntry(bwgn -> Clone(), Form("#sigma=%.2f sl=%.2f",
            bwgn->GetParameter(3)*1e3, bwgn->GetParameter(4) ), "L");

   bwgn -> SetParameter(4, -1.9 );  // slope
   bwgn -> SetLineColor(kGreen);
   bwgn -> DrawCopy("SAME");
   leg -> AddEntry(bwgn -> Clone(), Form("#sigma=%.2f sl=%.2f",
            bwgn->GetParameter(3)*1e3, bwgn->GetParameter(4) ), "L");

   leg -> Draw();
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
   far -> SetParameter(0,0.98);
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
   fbkg -> SetParameter(0, 0.98);
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
vector<double> IntfrBWAR( double m,
                  double mphi, double gphi,      // B-W
                  double A, double F, double Ang // Argus & interference
               ) {
//----------------------------------------------------------------------
// see memo and BreitWigner() function
// The return vector contains:
//                   0    1     2        3
//        ret     { Sum, B-W, Argus, Interference }
//----------------------------------------------------------------------
   constexpr double Mjpsi = 3.096916; // 3096.916 +/- 0.011 MeV
   constexpr double Meta  = 0.547862; // 547.862  +/- 0.017 MeV
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV

   constexpr double Mjpsi2 = SQ(Mjpsi);
   constexpr double Meta2  = SQ(Meta);
   constexpr double Mk2    = SQ(Mk);

   constexpr double dL = 2*Mk; // the left cutoff = 0.987354
   if ( m < dL ) { // phase space becomes zero
      vector<double> ret(4,0.);
      return ret;
   }

   double m2    = SQ(m);
   double mphi2 = SQ(mphi);
   double gphi2 = SQ(gphi);

   double gam = mphi*sqrt(mphi2+gphi2);
   double kap = 2*M_SQRT2*mphi*gphi*gam / (M_PI*sqrt(mphi2+gam));

   double p_phi0 = sqrt(SQ(Mjpsi2-mphi2-Meta2)-4*mphi2*Meta2) / (2*Mjpsi);
   double p_phi  = sqrt(SQ(Mjpsi2-m2-Meta2)-4*m2*Meta2) / (2*Mjpsi);
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
//    double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k; // == m*G(m)
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

// remove slope parameter for ToyMC
//----------------------------------------------------------------------
double IntfrBWARG( double m,
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
   auto Lint = [m,p,idx](double t) -> double{
      // t == (X'-X)/sigma
      double x = m+t*p[2];
      vector<double> Fun = IntfrBWAR(x, p[0],p[1],p[3],p[4],p[5]);
      return exp(-0.5*t*t) * Fun[idx];
   };
   rmath_fun< decltype(Lint) > Fint(Lint);

   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);
   double result = gsl_int.Integral(Fint,-5.,+5.);

   return one_over_sqrt2pi * result;
}

//----------------------------------------------------------------------
double IntfrBWARGN( double m,
                    const double p[], //{mphi,gphi,sigma,A,F,Ang}
                    int idx = 0       // what to return
                  ) {
//----------------------------------------------------------------------
// Numerical normalisation to one for range [dL,dU]
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV
   constexpr double dL = 2*Mk; // the left cutoff = 0.987354
   constexpr double dU = 1.08; // upper limit

   // cash parameters
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
      auto Lint = [](double x, void* pp) -> double{
         const double* p = static_cast<const double*>(pp);
         return IntfrBWARG(x,p,0);  // MUST be zero!
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-5, 1.e-4, 1000);
      norm = gsl_int.Integral(Lint,(void *)p,dL,dU);

      // save to cache
      cacheN = norm;
      for (int i = 0; i < nP; ++i) {
         cP[i] = p[i];
      }
   }

   return IntfrBWARG(m,p,idx) / norm;
}

//-------------------------------------------------------------------------
void test_Intfr() {
//-------------------------------------------------------------------------
   constexpr double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
   constexpr double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV

   constexpr double dL = 0.98;
   constexpr double dU = 1.08;

   // ---------------------------------------------------------------------
   // create Unuran 1D distribution object
   auto Lgen = [](const double* x,const double* p) -> double {
      vector<double> res = IntfrBWAR( x[0], p[0],p[1],p[2],p[3],p[4] );
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
   auto Lintfr = [Wb,mphi,gphi,ag,Fg,Ang](double* x, double* p) -> double {
      vector<double> res = IntfrBWAR( x[0], mphi,gphi,ag,Fg,Ang );
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

// {{{1 Monte Carlo
//-------------------------------------------------------------------------
struct GPar{
   double mphi;
   double gphi;
   double sigma_mkk;
   double a;
   double fb;
   double ang;
   double bg_rate; // flat background ???
};
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
double IntegNphi_1( const GPar* gp ) {
//-------------------------------------------------------------------------
   constexpr double Lphi = 0.98, Uphi = 1.08; // boundaries
   constexpr double eps_abs = 1.e-6;
   constexpr double eps_rel = 1.e-6;
   auto Lfc = [](const double* x,const double* p) -> double {
      return IntfrBWARGN(x[0],p,1); // 1 == BW
   };
   TF1 FC("FC", Lfc, Lphi, Uphi, 6);
   FC.SetParameters( gp -> mphi, gp -> gphi, gp -> sigma_mkk,
                     gp -> a, gp -> fb, gp -> ang );

   double err_num = 0;
   double Integ = FC.IntegralOneDim(Lphi,Uphi,eps_rel,eps_abs, err_num);
   if ( err_num > 1e-3 ) {
      cout << " WARNING: " << __func__
         << " numerical error of integration is too big "
         << err_num << endl;
   }
   return Integ;
}

//-------------------------------------------------------------------------
double IntegNphi( const GPar* gp ) {
//-------------------------------------------------------------------------
// Numerical integration to calculate B-W part
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV
   constexpr double dL = 2*Mk; // the left cutoff = 0.987354
   constexpr double dU = 1.08; // upper limit

   auto Lall = [](double x, void* pp) -> double {
      const double* p = static_cast<const double*>(pp);
      vector<double> ret = IntfrBWAR(x, p[0],p[1],p[2],p[3],p[4]);
      return ret[0];
   };
   const double p[] { gp -> mphi, gp -> gphi,
                      gp -> a, gp -> fb, gp -> ang };
   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-6, 1.e-6, 1000);
   double norm = gsl_int.Integral(Lall,(void *)p,dL,dU);

   auto Lbw = [](double x, void* pp) -> double {
      const double* p = static_cast<const double*>(pp);
      vector<double> ret = IntfrBWAR(x, p[0],p[1],p[2],p[3],p[4]);
      return ret[1];
   };
   double BW = gsl_int.Integral(Lbw,(void *)p,dL,dU);

   return BW/norm;
}

//-------------------------------------------------------------------------
vector<double> get_ToyMC(int Nevents, const GPar* gp) {
//-------------------------------------------------------------------------
   constexpr double Lphi = 0.98, Uphi = 1.08; // boundaries

   // ---------------------------------------------------------------------
   // set UNU.RAN ( http://statistik.wu-wien.ac.at/unuran/ )
   int logLevel = 2;
   TUnuran unr(gRandom,logLevel);

   // create Unuran 1D distribution object
   auto Lgen = [](const double* x,const double* p) -> double {
      vector<double> res = IntfrBWAR( x[0], p[0],p[1],p[2],p[3],p[4] );
      return res[0];
   };
   TF1* fgen = new TF1("fgen", Lgen, Lphi,Uphi, 5);
   fgen -> SetParameters( gp -> mphi, gp -> gphi,
                          gp -> a, gp -> fb, gp -> ang );

   double Xmax = fgen -> GetMaximumX(gp->mphi-0.4e-3,gp->mphi+0.4e-3);
//    printf(" Xmax= %.7f (delta= %.1e)\n",Xmax,Xmax-mphi);

   TUnuranContDist dist(fgen);
   dist.SetDomain(Lphi,Uphi);
   dist.SetMode(Xmax);

   // select unuran method for generating the random numbers
   string method = "method=nrou";
   cout << " start Unuran initializing with method: " << method << endl;
   if ( !unr.Init(dist,method) ) {
      cout << "Error initializing unuran" << endl;
      exit(0);
   }

   // ---------------------------------------------------------------------
   // generate
   vector<double> Vm;
   Vm.reserve( Nevents );

   for(int i = 0; i < Nevents; i++) {
      // 1) flat background
//       if ( random -> Uniform() < bg_rate ) {
//          double mphi = random -> Uniform(Lphi,Uphi);
//          Vmphi.push_back(mphi);
//          hst[0] -> Fill(mphi);
//          continue;
//       }

      // 2) interference + Gauss
      double mphi = unr.Sample();
      double mkk = gRandom -> Gaus(mphi,gp -> sigma_mkk);
      if ( mkk >= Lphi && mkk < Uphi ) {
         Vm.push_back(mkk);
      }
   }

   return Vm;
}

//-------------------------------------------------------------------------
TH1D* get_hst( const vector<double>& mkk, string hname) {
//-------------------------------------------------------------------------
   constexpr double Lphi = 0.98, Uphi = 1.08; // boundaries
   string title(";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2}");
   title += string(";Entries / 1 MeV/c^{2}");
   TH1D* hst = new TH1D(hname.c_str(), title.c_str(), 100,Lphi,Uphi);
   for ( auto m : mkk ) {
         hst -> Fill(m);
   }
   return hst;
}

// {{{1 Fit
//----------------------------------------------------------------------
ValErr CalcIntEr( const ROOT::Fit::FitResult& res,
                  unsigned int idx, double dL, double dU ) {
//----------------------------------------------------------------------
// Numerical calculation of integrals and their errors
   // desired errors:
   constexpr double eps_abs = 1.e-6;
   constexpr double eps_rel = 1.e-6;
   constexpr double epsilon  = 1.e-3; // max numerical error

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
   auto Lfc = [idx](const double* x,const double* p) -> double {
      return IntfrBWARGN(x[0],p,idx);
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

//----------------------------------------------------------------------
void do_fit(const GPar* gp, const vector<double>& mkk,
            TH1D* hst, double NphiMC, int np, string pdf) {
//----------------------------------------------------------------------
   constexpr double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
   constexpr double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV
   constexpr double Lphi = 0.98, Uphi = 1.08; // boundaries

   // prepare data
   double norm = hst -> Integral();
   int Nb = hst -> GetNbinsX();
   double bW = (Uphi-Lphi)/Nb; // bin width

   int n = mkk.size();
   ROOT::Fit::DataRange dr(Lphi,Uphi);
   ROOT::Fit::UnBinData Dat(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Dat.Add(mkk[i]);
   }
   cout << " norm= " << norm << " n= " << n << endl;
//    return;

   //--------------------------------------------------------------------
   // function MUST be normalized to 1 on the fit range
   auto Lintfr = [](const double* x,const double* p) -> double {
      const double pp[] { p[0],p[1],p[2],p[3],p[4],p[5] };
      return IntfrBWARGN( x[0], pp, 0 );
   };

   vector<string> par_name { "M#phi", "G#phi", "#sigma", "A",
                             "F", "#vartheta" };
   vector<double> par_ini {Mphi,Gphi,1.2e-3, 0.,1.0,0.8};
   if ( np > 0 ) {
      par_ini.back() *= -1;
   }

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* fintfr = new TF1("fintfr", Lintfr, Lphi, Uphi, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 WFun( *fintfr );
   fitter.SetFunction( WFun,false); // false == no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   fitter.Config().ParSettings(0).SetValue(1.01953);        // like MC
   fitter.Config().ParSettings(0).Fix();                    // Mphi
   fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3, Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();                    // Gphi
   fitter.Config().ParSettings(2).SetLimits(0.5e-3, 2.e-3); // sigma
   fitter.Config().ParSettings(3).Fix();                    // A
   fitter.Config().ParSettings(4).SetLimits(0.01, 10.);     // F
   fitter.Config().ParSettings(5).SetLimits(-M_PI, M_PI);   // vartheta

   fitter.LikelihoodFit( Dat, false ); // true == extended likelihood fit
   fitter.CalculateMinosErrors();

   const ROOT::Fit::FitResult& res = fitter.Result();
   res.Print(cout);
//    res.PrintCovMatrix(cout); // print error matrix and correlations
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   vector<string> names { "N(KK)", "Nphi", "Nnonphi", "Nifr"  };
   vector<ValErr> Nos;
   Nos.reserve(4);
   for ( int idx = 0; idx <= 3; ++idx ) {
      ValErr Integ = CalcIntEr( res, idx, Lphi, Uphi );
      double Num = norm * Integ.val;
      double Err = fabs(Num) * sqrt( 1./norm + SQ(Integ.err/Integ.val) );
      Nos.push_back( ValErr {Num,Err} );
      printf("%s = %s\n",names[idx].c_str(),Nos[idx].prt(".1f"));
   }
   double Nphi = Nos[1].val, err_Nphi = Nos[1].err;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof = [Lintfr,Fpar](double x) -> double {
      return Lintfr(&x,Fpar.data());
   };
   rmath_fun< decltype(Lgof) > ftest(Lgof);
   ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(
         mkk.size(),mkk.data(),ftest,ROOT::Math::GoFTest::kPDF,Lphi,Uphi);
   double pvalueKS = goftest -> KolmogorovSmirnovTest();
   cout << " pvalueKS= " << pvalueKS << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   double Wn = bW*norm;
   auto Lfit = [Wn,Fpar](double* x,double* p) -> double {
      const double pp[]
         { Fpar[0],Fpar[1],Fpar[2],Fpar[3],Fpar[4],Fpar[5] };
      return Wn*IntfrBWARGN( x[0], pp, int(p[0]) );
   };
   TF1* ffit = new TF1("ffit", Lfit, Lphi, Uphi, 1);
   ffit -> SetLineWidth(2);
   ffit -> SetLineColor(kRed);
   ffit -> SetNpx(500);

   //--------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd(1);
   gPad -> SetGrid();

   SetHstFace(hst);
   hst -> SetMinimum(-50.);
   hst -> SetMaximum(400. ); // the same scale!
//    if ( np < 0 ) {
//       hst -> SetMinimum(-40.);
//       hst -> SetMaximum( 380. ); // TODO ??
//    } else {
//       hst -> SetMinimum(-10.);
//    }
   hst -> GetYaxis() -> SetTitleOffset(1.3);
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
   hst -> SetMarkerStyle(20);

   hst -> Draw("EP");

   TLegend* leg = new TLegend(0.52,0.72,0.89,0.89);
   leg -> AddEntry(hst,"Toy MC","LEP");
//    leg -> AddEntry( (TObject*)0, "M_{#phi}, #Gamma_{#phi} - PDG", "");
//    leg -> AddEntry( (TObject*)0,
//          Form("#sigma= %.1f MeV",gp -> sigma_mkk*1e3), "");
//    leg -> AddEntry( (TObject*)0,
//          Form("a= %.1f, F= %.1f, #vartheta= %.1f",gp->a,gp->fb,gp->ang),"");
//    leg -> AddEntry( (TObject*)0, Form("N_{#phi}= %.1f",NphiMC), "");

   ffit -> SetParameter(0, 0); // SUM
   ffit -> SetLineWidth(2);
   ffit -> SetLineColor(kRed);
   ffit -> DrawCopy("SAME");
   leg -> AddEntry( ffit -> Clone(), "Result of fit", "L");

   // draw components
   ffit -> SetLineWidth(1);
   ffit -> SetLineStyle(kDashed);

   ffit -> SetParameter(0, 1); // BW
   ffit -> SetLineColor(kGreen+2);
   ffit -> DrawCopy("SAME");
   leg -> AddEntry( ffit -> Clone(), "Breit-Wigner #phi#eta", "L");

   ffit -> SetParameter(0, 2); // Argus
   ffit -> SetLineColor(kBlue);
   ffit -> DrawCopy("SAME");
//    leg -> AddEntry( ffit -> Clone(), "Argus", "L");
   leg -> AddEntry( ffit -> Clone(), "Non-#phi KK#eta", "L");

   ffit -> SetParameter(0, 3); // interference
   ffit -> SetLineColor(kMagenta+1);
   ffit -> DrawCopy("SAME");
   leg -> AddEntry( ffit -> Clone(), "Interference", "L");

   leg -> Draw();

   TPaveText* pt = new TPaveText(0.52,0.38,0.89,0.71,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
//    pt -> AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt -> AddText( Form("#it{p-value(K-S) = %.3f}",pvalueKS) );
   pt -> AddText( Form("N_{#phi} = %.1f #pm %.1f",Nphi,err_Nphi) );
   pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt -> AddText( Form("#sigma= %s MeV", PE.Eform(2,".2f",1e3)) );
   pt -> AddText( Form("a= %s",PE.Eform(3,".1f")) );
   pt -> AddText( Form("F= %s",PE.Eform(4,".2f")) );
   pt -> AddText( Form("#vartheta= %s",PE.Eform(5,".2f")) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//----------------------------------------------------------------------
void ToyMC_fit(int np, string pdf) { // np = // +/- 1 pos/neg solutions
//----------------------------------------------------------------------
   constexpr double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
   constexpr double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV
   constexpr double Lphi = 0.98, Uphi = 1.08; // boundaries

   GPar gp1 {Mphi,Gphi,1.2e-3, 0.,1.,0.8}; // ~ 2012
   int Nevents = 2810;
//    int NphiMC_1 = Nevents * IntegNphi_1(&gp1);
//    cout << " NphiMC_1= " << NphiMC_1 << endl;
   int NphiMC = Nevents * IntegNphi(&gp1);
   cout << " NphiMC= " << NphiMC << endl;

   vector<double> mkk = get_ToyMC(Nevents, &gp1);
   TH1D* hst = get_hst( mkk, "toy_gp1" );

   do_fit(&gp1, mkk,hst,NphiMC, np, pdf);
}

// {{{1 MAIN:
//----------------------------------------------------------------------
void ToyMC() {
//----------------------------------------------------------------------
   gROOT -> Reset();
   gStyle -> SetOptStat(0);
   gStyle -> SetStatFont(62);
   gStyle -> SetLegendFont(42);

   // = define GSL error handler which does nothing =
   gsl_set_error_handler_off();

   // set integrator: ROOT::Math::GSLIntegrator adaptive method (QAG)
   ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Adaptive");

   // set random number generator: TRandom3 by default?
//    gRandom -> SetSeed(0); // seed is set to a random value

   ULong_t Iseed = 10;
   gRandom -> SetSeed(Iseed);

   // ------------- tests ---------------
//    test_BreitWigner();
//    test_RevAgrus();
//    test_Intfr();

   // ------------- ToyMC ---------------
   int np = +1; // +/- 1 => pos/neg solutions
   string pdf("toyMC_");
   pdf = pdf + to_string(Iseed) + ((np<0) ? "_n" : "_p") + ".pdf";

   ToyMC_fit(np,pdf);
}
