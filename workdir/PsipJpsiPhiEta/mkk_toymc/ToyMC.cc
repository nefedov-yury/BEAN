// toy MC for J/Psi -> K+ K- gamma gamma
// This is ROOT version
//
#include <time.h>
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
struct Gpar {
   Gpar(double brKKeta, double brphieta); // ctor

   void Print();
   vector<double> signal_ToyMC() const;
   vector<double> bg_ToyMC() const;

   double BrKKeta;
   double Brphieta;

   // signal
   double mphi;
   double gphi;
   double sig;
   double ar;
   double F;
   double ang;
   int NetaKK;

   // background
   double ar_bg;
   int Nbg;

   // boundaries
   double Lphi;
   double Uphi;
   double TwoMk;

   private:
   double CalcF();
};

//-------------------------------------------------------------------------
Gpar::Gpar(double brKKeta, double brphieta) {
//-------------------------------------------------------------------------
   constexpr double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
   constexpr double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV

   // 2012 with bg: fig.23-24
//    constexpr double Breta = 0.3941; // Br(eta->2gamma) = 39.41%
   constexpr double Nini = 6.04e6; // ~ N(Jpsi)*eff*Breta

   // save branchings
   BrKKeta = brKKeta;
   Brphieta = brphieta;

   // boundaries
   Lphi = 0.98;
   Uphi = 1.08;
   TwoMk = 2*Mk;

   // signal
   mphi = Mphi;
   gphi = Gphi;
   sig = 1.2e-3;
   ar = 0.;
   F = 0.;
   ang = 0;

   NetaKK = int(Nini * BrKKeta);  // ~ 2718
   F = this->CalcF();

   // background
   ar_bg = 0;
   Nbg = 35;
};

//-------------------------------------------------------------------------
void Gpar::Print() {
//-------------------------------------------------------------------------
   cout << " -------- Toy MC parameters --------- " << endl;
   cout << " Br(KKeta) = " << BrKKeta*1e4 << " *10^{-4}" << endl;
   cout << " Br(phi eta) = " << Brphieta *1e4 << " *10^{-4}" << endl;
   cout << " Mphi = " << mphi << " GeV" << endl;
   cout << " Gphi = " << gphi << " GeV" << endl;
   cout << " sigma(mkk) = " << sig*1e3 << " MeV" << endl;
   cout << " a (Argus) = " << ar << endl;
   cout << " F = " << F << endl;
   cout << " ang = " << ang << endl;
   cout << " NetaKK  = " << NetaKK << endl;
   cout << " Nbg     = " << Nbg << endl;
   cout << " a (Bgr Argus) = " << ar_bg << endl;
   cout << " ---- End of Toy MC parameters ------ " << endl << endl;
};

//-------------------------------------------------------------------------
double Gpar::CalcF() {
//-------------------------------------------------------------------------
   constexpr double Brphi = 0.492;  // Br(phi->K+K-) = 49.2%

   // integrand lambda function
   double pp[] { mphi, gphi, sig, ar, 1.,0.,1.}; // F=1,Ang=0,idx=1
   auto Lint = [](double x, void* pp) -> double {
      const double* p = static_cast<const double*>(pp);
      int idx = int(p[6]);
      return IntfrBWARG(x,p,idx);
   };

   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-9, 1.e-9, 1000);

   // calculate integrals
   auto checkInt = [](string t, double Int, double pp[]) -> void {
      if ( !isfinite(Int) ) {
         printf("%s: pp= %g,%g,%g,%g,%g,%g,%g\n", t.data(),
               pp[0],pp[1],pp[2],pp[3],pp[4],pp[5],pp[6]);
         exit(1);
      }
   };
   pp[6] = 1; // idx=1 BW-part
   double IBW = gsl_int.Integral(Lint,(void *)pp,TwoMk,Uphi);
   checkInt("IBW",IBW,pp);
   pp[6] = 2; // idx=2 Argus-part
   double IAr = gsl_int.Integral(Lint,(void *)pp,TwoMk,Uphi);
   checkInt("IAr",IAr,pp);
   pp[6] = 3; // idx=3 - Interference
   pp[5] = 0.; // Ang=0 -> cos* part
   double IIc = gsl_int.Integral(Lint,(void *)pp,TwoMk,Uphi);
   checkInt("IIc",IIc,pp);
   pp[5] = M_PI/2; // Ang -> sin* part
   double IIs = gsl_int.Integral(Lint,(void *)pp,TwoMk,Uphi);
   checkInt("IIs",IIs,pp);

   // solve quadratic equation for F
   double A = IAr;
   double B = IIc * cos(ang) + IIs * sin(ang);
   double C = IBW*(1. - BrKKeta/(Brphieta*Brphi));
   double Dis = B*B-4*A*C;
   double Ff = 0;
   if ( Dis >= 0 ) {
      Ff = (-B + sqrt(Dis))/(2*A); // only positive root
   }
   if ( Dis < 0 || Ff < 0 ) {
      printf(" FATAL ERROR in Gpar::CalcF:\n BrKKeta= %e Brphieta=%e"
            " A=%.2g B=%.2g C=%.2g Dis=%.2g Ff=%g\n",
            BrKKeta,Brphieta,A,B,C,Dis,Ff);
      exit(1);
   }

   return Ff;
}

// interference BW with Argus
//-------------------------------------------------------------------------
vector<double> Gpar::signal_ToyMC() const {
//-------------------------------------------------------------------------
   // ---------------------------------------------------------------------
   // set UNU.RAN ( http://statistik.wu-wien.ac.at/unuran/ )
   // ---------------------------------------------------------------------
   static TUnuran* unr = nullptr;
   if ( !unr ) {
      int logLevel = 2;
      unr = new TUnuran(gRandom,logLevel);
      // create Unuran 1D distribution object
      auto Lgen = [](const double* x,const double* p) -> double {
         vector<double> res = IntfrBWAR( x[0],p[0],p[1],p[2],p[3],p[4] );
         return res[0];
      };
      TF1* fgen = new TF1("fgen", Lgen, Lphi, Uphi, 5);
      fgen -> SetParameters( mphi, gphi, ar, F, ang );

      double Xmax = fgen -> GetMaximumX( mphi-0.4e-3, mphi+0.4e-3 );
//    printf(" Xmax= %.7f (delta= %.1e)\n",Xmax,Xmax-mphi);

      TUnuranContDist dist(fgen);
      dist.SetDomain( Lphi, Uphi );
      dist.SetMode(Xmax);

      // select unuran method for generating the random numbers
      string method = "method=nrou";
      cout << "start unuran initializing with method: " << method
           << endl;
      if ( !unr -> Init(dist,method) ) {
         cout << "Error initializing unr" << endl;
         exit(0);
      }
   }

   // ---------------------------------------------------------------------
   // generate
   unsigned int n_etaKK = gRandom -> Poisson(NetaKK);
   vector<double> Vm;
   Vm.reserve( n_etaKK + 100 ); // additional place for background
   for( ; Vm.size() < n_etaKK; ) {
      double m_phi = unr -> Sample();
      double mkk = gRandom -> Gaus(m_phi,sig);
      if ( mkk >= Lphi && mkk < Uphi ) {
         Vm.push_back(mkk);
      }
   }
   return Vm;
}

// Argus for background
//-------------------------------------------------------------------------
vector<double> Gpar::bg_ToyMC() const {
//-------------------------------------------------------------------------
   // ---------------------------------------------------------------------
   // set UNU.RAN ( http://statistik.wu-wien.ac.at/unuran/ )
   // ---------------------------------------------------------------------
   static TUnuran* unrA = nullptr;
   if ( !unrA ) {
      int logLevel = 2;
      unrA = new TUnuran(gRandom,logLevel);
      // create Unuran 1D distribution object
      auto LgenA = [](const double* x,const double* p) -> double {
         return RevArgus( x[0], p[0] );
      };
      TF1* fgenA = new TF1("fgenA", LgenA, Lphi, Uphi, 1);
      fgenA -> SetParameter( 0, ar_bg );

      TUnuranContDist distA(fgenA);
      distA.SetDomain( Lphi, Uphi );

      // select unuran method for generating the random numbers
      string method = "tdr";
      cout << "start unuran initializing with method: " << method
           << endl;
      if ( !unrA -> Init(distA,method) ) {
         cout << "Error initializing unrA" << endl;
         exit(0);
      }
   }

   // ---------------------------------------------------------------------
   // generate
   unsigned int n_bg = gRandom -> Poisson(Nbg);
   vector<double> Vm;
   Vm.reserve( n_bg );
   for( ; Vm.size() < n_bg; ) {
      double m = unrA -> Sample();
      Vm.push_back(m);
   }
   return Vm;
}

// {{{1 do_fit_toy()
//--------------------------------------------------------------------
struct myFCN_toy {
   const double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV
   const double Mk    = 0.493677; // 493.677  +/- 0.016 MeV
   const double dL = 2*Mk, dU = 1.08; // boundaries

   vector<double> mkk;  // data central part
   vector<double> sb;   // data side band

   const double mphi = 1.01953; // Mphi like in MC
   const double gphi = Gphi;
   const double ar = 0.; // Argus parameter
   double Br2Nkk = 0.;   // convertion Br(J/Psi -> KK eta) -> Nkk
   double Br2Nphi = 0.;  // convertion Br(J/Psi -> phi eta) -> Nphi

   // integrals
   double IBW = 0, IAr = 0, IIc = 0, IIs = 0;

   // Argus for side-band:
   const double arsb = 0.; // Argus parameter for SB
   double normArsb = 0;    // normalization for Argus

   //-----------------------------------------------------------------
   myFCN_toy() { // set parameters
   //-----------------------------------------------------------------
      // Br(eta->2gamma) = 39.41%
//       const double breta = 0.3941;
      // Br(phi->K+K-) = 49.2%
      const double brphi = 0.492;

      constexpr double Nini = 6.04e6; // ~ N(Jpsi)*eff*Breta
      Br2Nkk = Nini;
      Br2Nphi = Br2Nkk * brphi;

      // normalization for Argus SB
      double par[] = {arsb};
      auto Lintar = [](double x, void* pp) -> double{
         double* p = static_cast<double*>(pp);
         return RevArgus(x,p[0]);
      };
      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);
      normArsb = gsl_int.Integral(Lintar,(void *)par,dL,dU);
   }

   //-----------------------------------------------------------------
   void calcIntegrals(double sig) {
   //-----------------------------------------------------------------
      // cach old values
      static double sig_save = 0.;
      if ( fabs(sig_save-sig) < 1e-6*fabs(sig_save) ) {
         return;
      }
      sig_save = sig;

      // integrand lambda function
      const double F = 1.;
      const double ang = 0.; //0(PI/2) for cos(sin) members of Infr.
      double pp[] { mphi,gphi,sig,ar,F,ang, 1 };

      auto Lint = [](double x, void* pp) -> double{
         const double* p = static_cast<const double*>(pp);
         int idx = int(p[6]);
         return IntfrBWARG(x,p,idx);
      };
      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);

      // calculate integrals
      auto checkInt = [](string t, double Int, double pp[]) -> void {
         if ( !isfinite(Int) ) {
            printf("%s: pp= %g,%g,%g,%g,%g,%g,%g\n", t.data(),
                  pp[0],pp[1],pp[2],pp[3],pp[4],pp[5],pp[6]);
            exit(1);
         }
      };
      pp[6] = 1; // idx
      IBW = gsl_int.Integral(Lint,(void *)pp,dL,dU);
      checkInt("IBW",IBW,pp);
      pp[6] = 2;
      IAr = gsl_int.Integral(Lint,(void *)pp,dL,dU);
      checkInt("IAr",IAr,pp);
      pp[6] = 3;
      IIc = gsl_int.Integral(Lint,(void *)pp,dL,dU);
      checkInt("IIc",IIc,pp);
      pp[5] = M_PI/2; // ang
      IIs = gsl_int.Integral(Lint,(void *)pp,dL,dU);
      checkInt("IIs",IIs,pp);

      // debug print
//       printf(" IeBW= %.3g IeAr= %.3g IeIc= %.3g IeIs= %.3g\n",
//             IeBW, IeAr, IeIc, IeIs);
//       printf(" IBW= %.3g IAr= %.3g IIc= %.3g IIs= %.3g\n",
//             IBW, IAr, IIc, IIs);
   }

   // intergrals MUST be calculeted before calling this function
   //-----------------------------------------------------------------
   tuple<double,double> calcF(const double* p) const {
   //-----------------------------------------------------------------
      double Brkk = p[0];
      double Nkk = Brkk * Br2Nkk;

      double Brphi = p[1];
      double Nphi = Brphi * Br2Nphi;

      double ang = p[2];

      // calculate F:
      double A = IAr;
      double B = IIc * cos(ang) + IIs * sin(ang);
      double C = IBW*(1. - Nkk/Nphi);
      double Dis = B*B-4*A*C;
      double Ff = 0., penalty = 0.;
      if ( Dis < 0 ) { // return "minimum"
         Ff = -B/(2*A);
         penalty += 1e5*fabs(Dis);
      } else {
         Ff = (-B + sqrt(Dis))/(2*A);
      }
      if ( Ff < 0. ) {
         penalty += 10.;
      }

      // debug print
      if ( Dis < 0 || Ff < 0 ) {
         printf(" calcF: Nkk=%.1f Nphi=%.1f A=%.2g B=%.2g C=%.2g"
               " Dis=%.2g Ff=%g %g\n",Nkk,Nphi,A,B,C,Dis,Ff,penalty);
      }
      return make_tuple(Ff,penalty);
   }

   // minimization function
   // the signature of this operator() MUST be exactly this:
   // THE EXTENDED MAXIMUM LIKELIHOOD ESTIMATOR
   // Kai-Feng Chen "FittingMinuitRooFit" p.35
   //-----------------------------------------------------------------
   double operator() (const double* p) {
   //-----------------------------------------------------------------
      if ( !isfinite(p[0]) || !isfinite(p[1]) ||
           !isfinite(p[2]) || !isfinite(p[3]) || !isfinite(p[4]) ) {
         printf("NAN: p[]= %g,%g,%g,%g,%g\n",
               p[0],p[1],p[2],p[3],p[4]);
         return DBL_MAX;
      }

      double Brkk = p[0];
      double Brphi = p[1];
      double Nkk = Brkk * Br2Nkk;
      double Nphi = Brphi * Br2Nphi;

      double ang = p[2];
      double sigma = p[3];
      double Nsb = p[4];

      calcIntegrals(sigma);
      double Ff = 0., penalty = 0.;
      tie(Ff,penalty) = calcF(p);
//       if ( Ff < 0. ) {
//          printf("F NEG: p[]= %.2g,%.2g,%.2g,%.2g,%.2g\n",
//                p[0],p[1],p[2],p[3],p[4]);
//       }

//+ this long calculation depends only on Br(eta->2gamma)...
//+       double normI = IBW+Ff*((IIc*cos(ang)+IIs*sin(ang))+Ff*IAr);
//+       double NFit = Nkk/normI;
//+ short computation is mathematically completely equivalent to long
      double NFit = Nphi/IBW; // === Nkk/normI;

      // for toyMC normE == normI => NkkFit == Nkk
      double normE = IBW+Ff*((IIc*cos(ang)+IIs*sin(ang))+Ff*IAr); 
      double NkkFit = NFit*normE;

      double res = 2*(NkkFit+Nsb) + penalty;

      const double pp[] { mphi,gphi,sigma,ar,Ff,ang };
      const double NsbNorm = Nsb / normArsb;
      int n_mkk = mkk.size();
      for ( int i = 0; i < n_mkk; ++i ) {
         double m = mkk[i];
         double L = NFit * IntfrBWARG( m,pp,0 ) +
                    NsbNorm * RevArgus(m,arsb);
         if (L > 0.) {
            res -= 2*log(L);
         } else {
            res += FLT_MAX; // ~3e38
         }
      }

      // fit SB
      int n_sb = sb.size();
      res += 2*Nsb;
      for ( int i = 0; i < n_sb; ++i ) {
         double m = sb[i];
         double L = NsbNorm * RevArgus(m,arsb);
         if (L > 0.) {
            res -= 2*log(L);
         } else {
            res += FLT_MAX; // ~3e38
         }
      }

//       printf("RES: p[]= %g,%g,%g,%g,%g -> %g\n",
//               p[0],p[1],p[2],p[3],p[4],res);
      return res;
   }
};

//--------------------------------------------------------------------
void do_fit_toy(myFCN_toy& my_fcn, TH1D* hist[], string pdf="" ) {
//--------------------------------------------------------------------
   TH1D* hst = hist[0];
   const vector<double>& mkk = my_fcn.mkk;
   int n = mkk.size();

   TH1D* hsb = hist[1];
   const vector<double>& sb = my_fcn.sb;
   int nsb = sb.size();

   const double& dL = my_fcn.dL; 
   const double& dU = my_fcn.dU; 
   double bW = hst -> GetBinWidth(1); // bin width
   double bL = hst -> GetBinLowEdge(1); // left boundary of hst

   //-----------------------------------------------------------------
   // Fit data
   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   vector<string> par_name { "Brkk", "Brphi", "angle", "sigma", "Nbg" };

   vector<double> par_ini {4.5e-4, 8.5e-4, 0., 1.2e-3, 35.};

   const unsigned int Npar = par_name.size(); // number of parameters

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(3e-4, 6e-4);  // Brkk
   fitter.Config().ParSettings(1).SetLimits(5e-4, 12e-4); // Brphi
   fitter.Config().ParSettings(2).SetLimits(-M_PI, M_PI); // angle
   fitter.Config().ParSettings(3).SetLimits(0.5e-3,2.e-3);// sigma
//    fitter.Config().ParSettings(3).Fix();
//    fitter.Config().ParSettings(4).SetLimits(0., 1.5*nsb); // Nbg

   // == Fit
   int Ndat = n + nsb;
   fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // false=likelihood
//    fitter.CalculateMinosErrors();
   fitter.CalculateHessErrors(); //in case of Minos find a new minimum

   const ROOT::Fit::FitResult& res = fitter.Result();
   res.Print(cout);

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test (see goftest from ROOT-tutorial)
   auto Lcr = [Fpar,my_fcn](double x) -> double {
      double Nphi = Fpar[1] * my_fcn.Br2Nphi;
      double NFit = Nphi/my_fcn.IBW;
      double Ff,penalty;
      tie(Ff,penalty) = my_fcn.calcF(Fpar.data());

      // {mphi,gphi,sigma,A,F,Ang,sl}
      double ang = Fpar[2];
      double sig = Fpar[3];
      const double pp[] {
         my_fcn.mphi, my_fcn.gphi, sig, my_fcn.ar, Ff, ang
      };
      double NsbNorm = Fpar[4] / my_fcn.normArsb;
      return NFit*IntfrBWARG( x,pp,0 ) +
         NsbNorm * RevArgus(x,my_fcn.arsb);
   };
   rmath_fun< decltype(Lcr) > fcr(Lcr);
   ROOT::Math::GoFTest* gofcr =
      new ROOT::Math::GoFTest( mkk.size(),mkk.data(),fcr,
            ROOT::Math::GoFTest::kPDF, dL,dU );
   double pvalueKS = gofcr -> KolmogorovSmirnovTest();
   cout << " pvalueKS(cr)= " << pvalueKS << endl;

   auto Lsb = [Fpar,my_fcn](double x) -> double {
      double NsbNorm = Fpar[4] / my_fcn.normArsb;
      return NsbNorm * RevArgus(x,my_fcn.arsb);
   };
   rmath_fun< decltype(Lsb) > fsb(Lsb);
   ROOT::Math::GoFTest* gofsb =
      new ROOT::Math::GoFTest( sb.size(),sb.data(),fsb,
         ROOT::Math::GoFTest::kPDF, dL,dU );
   double pvalueKSsb = gofsb -> KolmogorovSmirnovTest();
   cout << " pvalueKS(sb)= " << pvalueKSsb << endl;

   //-----------------------------------------------------------------
   // Functions to draw
   auto Ldr = [Fpar,my_fcn,bW](double* x,double* p) -> double {
      int isw = int(p[0]+0.5);

      double Nphi = Fpar[1] * my_fcn.Br2Nphi;
      double NFit = Nphi/my_fcn.IBW;
      double Ff,penalty;
      tie(Ff,penalty) = my_fcn.calcF(Fpar.data());
      double ang = Fpar[2];
      double sig = Fpar[3];
      const double pp[] {
         my_fcn.mphi, my_fcn.gphi, sig, my_fcn.ar, Ff, ang
      };

      double BWARG = NFit * IntfrBWARG( x[0],pp,isw );
      if ( isw > 0 && isw < 4 ) { return bW * BWARG; }

      double NsbNorm = Fpar[4] / my_fcn.normArsb;
      double Bg = NsbNorm * RevArgus(x[0],my_fcn.arsb);
      if ( isw == 4 ) { return bW * Bg; }

      return bW * (BWARG + Bg);
   };
   TF1* fdr = new TF1("fdr", Ldr, bL, dU, 1);
   fdr -> SetNpx(500);

   // find max/min to draw
   fdr -> SetParameter(0, 1); // BW
   int maxbin = hst -> GetMaximumBin();
   double hmax = hst -> GetBinContent(maxbin) +
      hst -> GetBinError(maxbin);
   double fmax = fdr -> GetMaximum( 1.01, 1.03, 1e-6);
   hmax = floor(1.1 * max( hmax, fmax ));
   if ( hmax > 0 ) {
      hst -> SetMaximum(hmax);
   }

   fdr -> SetParameter(0, 3); // interference
   double hmin = fdr -> GetMinimum( 1.01, 1.03, 1e-6);
   hmin = floor(1.15 * min( 0., hmin ));
   if ( hmin < 0 ) {
      hst -> SetMinimum(hmin);
   }

   //-----------------------------------------------------------------
   // Draw results
   TLegend* leg = new TLegend(0.58,0.70,0.89,0.89);
   TPaveText* pt = new TPaveText(0.58,0.41,0.89,0.69,"NDC");
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();

   gPad -> SetGrid();

   SetHstFace(hst);
   hst -> GetXaxis() -> SetTitleOffset(1.1);
   hst -> GetYaxis() -> SetTitleOffset(1.3);
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
   hst -> SetMarkerStyle(20);

   hst -> Draw("EP");
   leg -> AddEntry(hst,"Toy MC","LEP");

   fdr -> SetParameter(0, 0); // SUM
   fdr -> SetLineWidth(2);
   fdr -> SetLineColor(kRed+1);
   fdr -> DrawCopy("SAME");
   leg -> AddEntry( fdr -> Clone(), "Result of fit", "L");

   fdr -> SetParameter(0, 1); // BW
   fdr -> SetLineWidth(2);
   fdr -> SetLineStyle(kDashed);
   fdr -> SetLineColor(kGreen+2);
   fdr -> DrawCopy("SAME");
   leg -> AddEntry( fdr -> Clone(), "Breit-Wigner #phi#eta", "L");

   fdr -> SetParameter(0, 2); // Argus
   fdr -> SetLineWidth(2);
   fdr -> SetLineStyle(kDashed);
   fdr -> SetLineColor(kBlue);
   fdr -> DrawCopy("SAME");
   leg -> AddEntry( fdr -> Clone(), "Non-#phi KK#eta", "L");

   fdr -> SetParameter(0, 3); // interference
   fdr -> SetLineWidth(2);
   fdr -> SetLineStyle(kDashed);
   fdr -> SetLineColor(kMagenta+1);
   fdr -> DrawCopy("SAME");
   leg -> AddEntry( fdr -> Clone(), "Interference", "L");

   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
   pt -> AddText( Form("#it{p-value = %.3f}",pvalueKS) );
   pt -> AddText( Form("Br(KK#eta)= %s #times10^{-4}",
            PE.Eform(0,".2f",1e4)) );
   pt -> AddText( Form("Br(#phi#eta)= %s #times10^{-4}",
            PE.Eform(1,".2f",1e4)) );
   pt -> AddText( Form("#vartheta = %s",PE.Eform(2,".2f")) );
   pt -> AddText( Form("#sigma = %s MeV", PE.Eform(3,".2f",1e3)) );
//    pt -> AddText( Form("NKK(fit)= %s", PE.Eform(3,".1f")) );
   pt -> AddText( Form("#lower[-0.1]{Nbg = %s}",PE.Eform(4,".1f")) );
   pt -> Draw();
   leg -> Draw();
   gPad -> RedrawAxis();

   c1 -> Update();
   if ( !pdf.empty() ) {
      pdf += ".pdf";
      c1 -> Print(pdf.c_str());
   }
}

// {{{1  ToyMC_fit

//-------------------------------------------------------------------------
TH1D* get_hst( const vector<double>& mkk, string hname) {
//-------------------------------------------------------------------------
   constexpr double Lphi = 0.98, Uphi = 1.08; // boundaries
   string title(";M^{ inv}_{ K^{#plus}K^{#minus }}, GeV/c^{2}");
   title += string(";Entries / 1 MeV/c^{2}");
   TH1D* hst = new TH1D(hname.c_str(), title.c_str(), 100,Lphi,Uphi);
   for ( auto m : mkk ) {
         hst -> Fill(m);
   }
   return hst;
}

//-------------------------------------------------------------------------
void plot_hst(TH1D* hst[], string pdf) {
//-------------------------------------------------------------------------
   gStyle->SetOptStat(1001110);//-name,+entrs,+mean,+rms,-under,-over,+int

   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,500);
   c1 -> Divide(2,1);

   c1 -> cd(1);
   gPad -> SetGrid();

   SetHstFace(hst[0]);
   hst[0] -> GetXaxis() -> SetTitleOffset(1.1);
   hst[0] -> GetYaxis() -> SetTitleOffset(1.3);
   hst[0] -> Draw("E");

   c1 -> cd(2);
   gPad -> SetGrid();

   SetHstFace(hst[1]);
   hst[1] -> GetXaxis() -> SetTitleOffset(1.1);

   TF1* pl0 = (TF1*)gROOT -> GetFunction("pol0") -> Clone();
   pl0 -> SetLineColor(kRed);
   pl0 -> SetLineWidth(2);
   pl0 -> SetLineStyle(kDashed);
   hst[1] -> Fit(pl0,"L","E");

   c1 -> Update();
   pdf += ".pdf";
   c1 -> Print( pdf.c_str() );
}

//----------------------------------------------------------------------
void ToyMC_fit(string pdf) {
//----------------------------------------------------------------------

   double BrKKeta  = 4.5e-4;
   double Brphieta = 8.5e-4;
   Gpar* gp = new Gpar( BrKKeta, Brphieta );
   gp -> Print();

   TH1D* hst[10];
   myFCN_toy my_fcn;  // class for 'FitFCN'
   my_fcn.mkk = gp -> signal_ToyMC();
   vector<double> mkk_bg = gp -> bg_ToyMC(); 
   my_fcn.mkk.insert( my_fcn.mkk.end(), mkk_bg.begin(),mkk_bg.end() );
   hst[0] = get_hst( my_fcn.mkk, "toy_sig" );

   my_fcn.sb = gp -> bg_ToyMC(); // sideband
   hst[1] = get_hst( my_fcn.sb, "toy_sb" );

//    plot_hst(hst,pdf);

   do_fit_toy(my_fcn,hst,pdf);
}

// {{{1 MAIN:
//----------------------------------------------------------------------
void ToyMC() {
//----------------------------------------------------------------------
   gROOT -> Reset();
   gStyle -> SetOptStat(0);
   gStyle -> SetLegendFont(42);

   // = define GSL error handler which does nothing =
   gsl_set_error_handler_off();

   // set integrator: ROOT::Math::GSLIntegrator adaptive method (QAG)
   ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Adaptive");

   // set random number generator: TRandom3 by default?
//    gRandom -> SetSeed(0); // seed is set to a random value

   ULong_t Iseed = (ULong_t)time(NULL);
//    ULong_t Iseed= 1638656110L;
   printf("   ULong_t Iseed= %luL;\n",Iseed);
   gRandom -> SetSeed(Iseed);

   // ------------- tests ---------------
//    test_BreitWigner();
//    test_RevAgrus();
//    test_Intfr();

   // ------------- ToyMC ---------------
   string pdf("toyMC");
   ToyMC_fit(pdf);
}
