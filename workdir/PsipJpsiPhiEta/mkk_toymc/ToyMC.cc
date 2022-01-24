// Toy MC for checking: fitting of M(KK); error estimation accuracy
// this is a case of the combined fit of central region and side-band
// try to reproduce combined fit of 2012 and 2009

#include <iostream>
#include <cmath>
#include <vector>
#include <getopt.h>
#include <time.h>

#include "gsl/gsl_errno.h"   // GSL error handler

#include <Fit/Fitter.h>
#include <Math/MinimizerOptions.h>
#include <Math/GSLIntegrator.h>
#include <Math/WrappedTF1.h>
#include <Math/GoFTest.h>

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TNtupleD.h>
#include <TRandom.h>
#include <TUnuran.h>
#include <TUnuranContDist.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>

#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TROOT.h>

using namespace std;

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

//--------------------------------------------------------------------
// adapter class for handling parameters and errors
//--------------------------------------------------------------------
struct ParAndErr {
   vector<double> Fpar; // parameters
   vector<double> Perr; // symmetric errors
   vector<double> Uerr; // upper minos errors
   vector<double> Lerr; // lower minos errors

   // ctor
   //-----------------------------------------------------------------
   ParAndErr(const ROOT::Fit::FitResult& res, double threshold = 0.1) {
   //-----------------------------------------------------------------
      Fpar = res.Parameters();
      Perr = res.Errors();

      int Npar = res.NPar();
      Lerr.resize(Npar,0.);
      Uerr.resize(Npar,0.);

      // if the upper and lower errors differ no more than the
      // 'threshold', the maximum of them is used as a symmetric error
//       constexpr double threshold = 0.1;
      for ( int np = 0; np < Npar; ++np ) {
         if ( res.HasMinosError(np) ) {
            double uerr = fabs(res.UpperError(np));
            double lerr = fabs(res.LowerError(np));
            if ( uerr*lerr == 0 ) {// one of the errors is not defined
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
   //-----------------------------------------------------------------
   const char* Eform (int i,string fm, double X = 1) {
   //-----------------------------------------------------------------
      if ( Perr[i] > 0 ) {
         string fmt = "%" +fm+ " #pm %" +fm;
         return Form(fmt.c_str(),Fpar[i]*X,Perr[i]*X);
      } else if ( Perr[i] < 0 ) {
         string fmt = "%" + fm
            + "^{#lower[0.15]{#kern[0.15]{#plus}#kern[0.05]{%"
            + fm + "}}" + "}"
            + "_{#lower[-0.15]{#kern[0.15]{#minus}#kern[0.05]{%"
            + fm + "}}" + "}";
         return Form(fmt.c_str(),Fpar[i]*X,Uerr[i]*X,Lerr[i]*X);
      }
      string fmt = "%" +fm+ " (fixed)";
      return Form(fmt.c_str(),Fpar[i]*X);
   }

   // return middle point = par(i) + 0.5(UpperError+LowerError)
   //-----------------------------------------------------------------
   double Middle (int i) {
   //-----------------------------------------------------------------
      return Fpar[i] + 0.5*(Uerr[i]-Lerr[i]);
   }
};

//----------------------------------------------------------------------
constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

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
   TRandom* random = new TRandom();
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
   TRandom* random = new TRandom();
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
   Gpar(double brKKeta, double brphieta, double sigma,
         double nj, double nbg, double arsb, double setF); // ctor

   void Print(string info);
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
   double Arsb;
   int Nbg;

   // boundaries
   double Lphi;
   double Uphi;
   double TwoMk;

   private:
   double CalcF();
};

// nj = N(Jpsi)*eff*Breta
//-------------------------------------------------------------------------
Gpar::Gpar(double brKKeta, double brphieta, double sigma,
      double nj, double nbg, double arsb, double setF=-1) {
//-------------------------------------------------------------------------
   constexpr double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
   constexpr double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV

//    constexpr double Breta = 0.3941; // Br(eta->2gamma) = 39.41%

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
   sig = sigma;
   ar = 0.;
   F = 0.;
   ang = 0;

   NetaKK = int(nj * BrKKeta);
   if ( setF > 0 ) {
      F = setF;
   } else {
      F = this->CalcF();
   }

   // background
   Nbg = nbg;
   Arsb = arsb;
};

//-------------------------------------------------------------------------
void Gpar::Print(string info) {
//-------------------------------------------------------------------------
   cout << " -------- Toy MC parameters: " << info << " --------- " << endl;
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
   cout << " a (Argus side-band) = " << Arsb << endl;
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
      fgenA -> SetParameter( 0, Arsb );

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

   vector<double> mkk09;  // data central part
   vector<double> mkk12;
   vector<double> sb09;   // data side band
   vector<double> sb12;

   const double mphi = 1.01953; // Mphi like in MC
   const double gphi = Gphi;
   const double ar = 0.; // Argus parameter

   double Br2Nkk09 = 0.;   // conversion Br(J/Psi -> KK eta) -> Nkk
   double Br2Nkk12 = 0.;
   double Br2Nphi09 = 0.;  // conversion Br(J/Psi -> phi eta) -> Nphi
   double Br2Nphi12 = 0.;

   // integrals 0 -> 09, 1 -> 12
   vector<double> IBW, IAr, IIc, IIs;

   // normalizations for Argus in side-band
   vector<double> normArsb;

   //-----------------------------------------------------------------
   myFCN_toy(double Nj09, double Nj12) { // set parameters
   //-----------------------------------------------------------------
      // Br(eta->2gamma) = 39.41%
//       const double breta = 0.3941;
      // Br(phi->K+K-) = 49.2%
      const double brphi = 0.492;

      Br2Nkk09 = Nj09;
      Br2Nkk12 = Nj12;
      Br2Nphi09 = Br2Nkk09 * brphi;
      Br2Nphi12 = Br2Nkk12 * brphi;

      IBW.resize(2,0.);
      IAr.resize(2,0.);
      IIc.resize(2,0.);
      IIs.resize(2,0.);

      // normalization for Argus SB
      normArsb.resize(2,1.);
   }

   // normalization for Argus SB
   //-----------------------------------------------------------------
   void calcNormAr(const double arsb[]) { // 0 -> 09, 1 -> 12
   //-----------------------------------------------------------------
      // cache for calculated values
      static double arsb_save[2] = {-1.,-1.};

      // integrand lambda function
      double par[] { arsb[0] };
      auto Lintar = [](double x, void* pp) -> double{
         double* p = static_cast<double*>(pp);
         return RevArgus(x,p[0]);
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);

      // integrals 0 -> 09, 1 -> 12
      for ( int i = 0; i < 2; ++i ) {
         if ( arsb_save[i] == arsb[i] ) {
            continue;
         }
         arsb_save[i] = arsb[i];
         par[0] = arsb[i]; // set argus
         normArsb[i] = gsl_int.Integral(Lintar,(void *)par,dL,dU);

         // debug print
//          printf(" %i: arsb= %.3f, normArsb= %.3g\n",
//                i, arsb[i], normArsb[i]);
      }
   }

   //-----------------------------------------------------------------
   void calcIntegrals(const double sig[]) { // 0 -> 09, 1 -> 12
   //-----------------------------------------------------------------
      // cache for calculated values
      static double sig_save[2] = {0.,0.};

      // integrand lambda function: mphi,gphi,sig,ar,F,ang,sl,idx
      const double F = 1.;
      const double ang = 0.; //0(PI/2) for cos(sin) members of Infr.
      double pp[] { mphi,gphi,0.,ar,F,ang, 1. };

      auto Lint = [](double x, void* pp) -> double{
         const double* p = static_cast<const double*>(pp);
         int idx = int(p[6]);
         return IntfrBWARG(x,p,idx);
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);

      // function for check results
      [[maybe_unused]]
      auto checkInt = [](string t, double Int, double pp[]) -> void {
         if ( !isfinite(Int) ) {
            printf("%s: pp= %g,%g,%g,%g,%g,%g,%g\n", t.data(),
                  pp[0],pp[1],pp[2],pp[3],pp[4],pp[5],pp[6]);
            exit(1);
         }
      };

      // integrals 0 -> 09, 1 -> 12
      for ( int i = 0; i < 2; ++i ) {
         if ( sig_save[i] == sig[i] ) {
            continue;
         }
         sig_save[i] = sig[i];
         pp[2] = sig[i]; // sigma

         pp[6] = 1; // idx
         IBW[i] = gsl_int.Integral(Lint,(void *)pp,dL,dU);
#ifndef BATCH
         checkInt(string("IBW_")+to_string(i),IBW[i],pp);
#endif
         pp[6] = 2;
         IAr[i] = gsl_int.Integral(Lint,(void *)pp,dL,dU);
#ifndef BATCH
         checkInt(string("IAr_")+to_string(i),IAr[i],pp);
#endif
         pp[6] = 3;
         pp[5] = 0; // ang=0 for cos
         IIc[i] = gsl_int.Integral(Lint,(void *)pp,dL,dU);
#ifndef BATCH
         checkInt(string("IIc_")+to_string(i),IIc[i],pp);
#endif
         pp[5] = M_PI/2; // ang = pi/2 for sin
         IIs[i] = gsl_int.Integral(Lint,(void *)pp,dL,dU);
#ifndef BATCH
         checkInt(string("IIs_")+to_string(i),IIs[i],pp);
#endif

         // debug print
//          printf(" %i: IBW= %.3g IAr= %.3g IIc= %.3g IIs= %.3g\n",
//                i, IBW[i], IAr[i], IIc[i], IIs[i]);
      }
   }

   // integrals MUST be calculated before calling this function
   //-----------------------------------------------------------------
   tuple<double,double> calcF12(const double* p) const {
   //-----------------------------------------------------------------
   // NOTE: the ratio Nkk/Nphi = Brkk / (Brphi*Br(phi->KK)) does not
   // depend on the year, but the integrals do, however very slightly
   // therefore we take integrals for 2012 year => sys: F09

      double Brkk = p[0];
      double Nkk = Brkk * Br2Nkk12;

      double Brphi = p[1];
      double Nphi = Brphi * Br2Nphi12;

      double ang = p[2];

      // calculate F:
      const int iy = 1; // 1 for 2012 and sys: 0 for 2009 (F09)
      double A = IAr[iy];
      double B = IIc[iy] * cos(ang) + IIs[iy] * sin(ang);
      double C = IBW[iy]*(1. - Nkk/Nphi);
      double Dis = B*B-4*A*C;
      double Ff = 0., penalty = 0.;
      if ( Dis < 0 ) { // return "minimum"
         Ff = max(0.,-B/(2*A));
         penalty += 1e5*fabs(Dis);
      } else {
         Ff = (-B + sqrt(Dis))/(2*A);
      }
      if ( Ff < 0. ) {
         penalty = -1e5*Ff;
         Ff = 0.;
      }

      // debug print
      if ( penalty > 0 ) {
         printf(" calcF penalty= %g: Nkk=%.1f Nphi=%.1f"
               " A=%.2g B=%.2g C=%.2g Dis=%.2g Ff=%g\n",
               penalty,Nkk,Nphi,A,B,C,Dis,Ff);
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
      if ( !isfinite(p[0]) || !isfinite(p[1]) || !isfinite(p[2]) ||
           !isfinite(p[3]) || !isfinite(p[4]) || !isfinite(p[5]) ||
           !isfinite(p[6]) || !isfinite(p[7]) || !isfinite(p[8]) ) {
         printf("NAN: p[]= %g,%g,%g,%g,%g,%g,%g\n",
               p[0],p[1],p[2],p[3],p[4],p[5],p[6]);
         return DBL_MAX;
      }

      double Brkk = p[0];
     [[maybe_unused]]
      double Brphi = p[1]; // Ff depends of it !
      double Nkk09  = Brkk  * Br2Nkk09;
     [[maybe_unused]]
      double Nphi09 = Brphi * Br2Nphi09;
      double Nkk12  = Brkk  * Br2Nkk12;
     [[maybe_unused]]
      double Nphi12 = Brphi * Br2Nphi12;

      double ang = p[2];
      double sig09 = p[3];
      double sig12 = p[4];
      calcIntegrals(&p[3]); // p[3],p[4]

      double Ff = 0., penalty = 0.;
      tie(Ff,penalty) = calcF12(p);

      double Nsb09 = p[5];
      double Nsb12 = p[6];
      double arsb09 = p[7];
      double arsb12 = p[8];
      calcNormAr(&p[7]); // p[7],p[8]
      const double NsbNorm09 = Nsb09 / normArsb[0];
      const double NsbNorm12 = Nsb12 / normArsb[1];

      // 2009
//+ see myFCN_sbbr::calcF for explanation of short computation
//+       double NFit09 = Nphi09/IBW[0]; // =? Nkk/normI;
//+ the calculation of F is ambiguous, it depends on the sigma
//+ different for different years, so I used long way
      double normI09 =
         IBW[0]+Ff*((IIc[0]*cos(ang)+IIs[0]*sin(ang))+Ff*IAr[0]);
      double NFit09 = Nkk09/normI09;

      // for toyMC normE == normI => NkkFit == Nkk
      double NkkFit09 = Nkk09;

      double res = 2*(NkkFit09+Nsb09) + penalty;

      const double p09[] { mphi,gphi,sig09,ar,Ff,ang };
      int n_mkk09 = mkk09.size();
      for ( int i = 0; i < n_mkk09; ++i ) {
         double m = mkk09[i];
         double L = NFit09 * IntfrBWARG( m,p09,0 ) +
                    NsbNorm09 * RevArgus(m,arsb09);
         if (L > 0.) {
            res -= 2*log(L);
         } else {
            res += FLT_MAX; // ~3e38
         }
      }

      // fit SB
      int n_sb09 = sb09.size();
      res += 2*Nsb09;
      for ( int i = 0; i < n_sb09; ++i ) {
         double m = sb09[i];
         double L = NsbNorm09 * RevArgus(m,arsb09);
         if (L > 0.) {
            res -= 2*log(L);
         } else {
            res += FLT_MAX; // ~3e38
         }
      }

      // 2012
      double normI12 =
         IBW[1]+Ff*((IIc[1]*cos(ang)+IIs[1]*sin(ang))+Ff*IAr[1]);
      double NFit12 = Nkk12/normI12;

      // for toyMC normE == normI => NkkFit == Nkk
      double NkkFit12 = Nkk12;

      res += 2*(NkkFit12+Nsb12);

      const double p12[] { mphi,gphi,sig12,ar,Ff,ang };
      int n_mkk12 = mkk12.size();
      for ( int i = 0; i < n_mkk12; ++i ) {
         double m = mkk12[i];
         double L = NFit12 * IntfrBWARG( m,p12,0 ) +
                    NsbNorm12 * RevArgus(m,arsb12);
         if (L > 0.) {
            res -= 2*log(L);
         } else {
            res += FLT_MAX; // ~3e38
         }
      }

      // fit SB
      int n_sb12 = sb12.size();
      res += 2*Nsb12;
      for ( int i = 0; i < n_sb12; ++i ) {
         double m = sb12[i];
         double L = NsbNorm12 * RevArgus(m,arsb12);
         if (L > 0.) {
            res -= 2*log(L);
         } else {
            res += FLT_MAX; // ~3e38
         }
      }

      return res;
   }
};

//--------------------------------------------------------------------
vector<double> do_fit_toy( myFCN_toy& my_fcn, TH1D* hist[],
      string pdf="" ) {
//--------------------------------------------------------------------
   TH1D* h09 = hist[0];
   const vector<double>& mkk09 = my_fcn.mkk09;
   int n09 = mkk09.size();

   [[maybe_unused]]
   TH1D* hsb09 = hist[1];
   const vector<double>& sb09 = my_fcn.sb09;
   int nsb09 = sb09.size();

   [[maybe_unused]]
   TH1D* h12 = hist[2];
   const vector<double>& mkk12 = my_fcn.mkk12;
   int n12 = mkk12.size();

   [[maybe_unused]]
   TH1D* hsb12 = hist[3];
   const vector<double>& sb12 = my_fcn.sb12;
   int nsb12 = sb12.size();

   const double& dL = my_fcn.dL;
   const double& dU = my_fcn.dU;

   bool isBatch = pdf.empty();

   //-----------------------------------------------------------------
   // Fit data
   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   vector<string> par_name { "Brkk", "Brphi", "angle",
      "sig09", "sig12", "Nbg09", "Nbg12", "Arsb09", "Arsb12" };

   vector<double> par_ini {4.5e-4, 8.5e-4, 0.,
      1.4e-3, 1.1e-3, 1.*nsb09, 1.*nsb12, 8., 5. };

   const unsigned int Npar = par_name.size(); // number of parameters

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(3e-4, 6e-4);   // Brkk
   fitter.Config().ParSettings(1).SetLimits(5e-4, 12e-4);  // Brphi
   fitter.Config().ParSettings(2).SetLimits(-M_PI, M_PI);  // angle
   fitter.Config().ParSettings(3).SetLimits(0.3e-3,3.e-3); // sig09
   fitter.Config().ParSettings(4).SetLimits(0.3e-3,3.e-3); // sig12
   fitter.Config().ParSettings(5).SetLimits(0.,2.*nsb09);  // Nbg09
   fitter.Config().ParSettings(6).SetLimits(0.,2.*nsb12);  // Nbg12
   fitter.Config().ParSettings(7).SetLimits(0.,15.);       // Arsb09
   fitter.Config().ParSettings(8).SetLimits(0.,15.);       // Arsb12

   // == Fit
   int Ndat = n09 + nsb09 + n12 + nsb12;
   fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // false=likelihood
   fitter.CalculateHessErrors();
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   vector<double> Fpar = res.Parameters();

   bool do_refit = false;
   double Ff,penalty;
   tie(Ff,penalty) = my_fcn.calcF12(Fpar.data());
   if ( !isfinite(penalty) || penalty > 0. ) {
      do_refit = true;
      printf("\n=> INVALID RESULT: penalty= %f <=\n",penalty);
   } else {
      printf("\n=> FINAL RESULT <=\n");
   }

#ifdef BATCH
   if ( !do_refit ) {
      // additional check upper and lower bounds after MINOS
      for ( unsigned int ip = 0; ip < Npar; ++ip ) {
         if ( !res.HasMinosError(ip) ||
               res.UpperError(ip)*res.LowerError(ip) == 0 ) {
            do_refit = true;
            printf("\n=> INVALID RESULT: minos errors <=\n");
            break;
         }
      }
   }
   if ( do_refit ) {
      res.Print(cout);
      par_ini = res.Parameters();
      if ( fabs(par_ini[2]) > 0.1 ) {
         par_ini[2] *= -1.;
      } else {
         par_ini[2] = copysign(0.5, par_ini[2]);
      }
      fitter.Config().SetParamsSettings(Npar,par_ini.data());
      fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false);
      fitter.CalculateHessErrors();
      fitter.CalculateMinosErrors();

      res = fitter.Result();
      Fpar = res.Parameters();
      tie(Ff,penalty) = my_fcn.calcF12(Fpar.data());
   }
   printf("\n=> FINAL RESULT AFTER REFIT <=\n");
#endif

   res.Print(cout);

   double arsb09 = Fpar[7];
   double arsb12 = Fpar[8];
   my_fcn.calcNormAr(&Fpar[7]);

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res,0.03); // ignore 3% upper/lower errors

   vector<double> F_ret;  // vector to return
   F_ret.reserve(4*PE.Fpar.size()+20);
   F_ret.insert(F_ret.end(),PE.Fpar.begin(),PE.Fpar.end()); // values
   F_ret.insert(F_ret.end(),PE.Uerr.begin(),PE.Uerr.end()); // upper err
   F_ret.insert(F_ret.end(),PE.Lerr.begin(),PE.Lerr.end()); // lower err
   const vector<double> & Perr = res.Errors();
   F_ret.insert(F_ret.end(),Perr.begin(),Perr.end()); // symmetric err
   F_ret.push_back(Lmin);

   // status of minimization
   double st = ((res.IsValid()) ? 1. : 0.);
   if ( st == 1 && do_refit ) {
      st = 0.5;
   }
   // something wrong with that solution
   if ( !isfinite(penalty) || penalty > 0. ) {
      st = -2.;
   }
   F_ret.push_back(st);

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test (see goftest from ROOT-tutorial)
   bool calc_pval = true;
   double pvalueKS09 = 0, pvalueKS09sb = 0;
   double pvalueKS12 = 0, pvalueKS12sb = 0;
   if ( calc_pval ) {
      // central part
      double Nkk09 = Fpar[0] * my_fcn.Br2Nkk09;
      double Nphi09 = Fpar[1] * my_fcn.Br2Nphi09;
      double NFit09 = Nphi09/my_fcn.IBW[0];
      double Nkk12 = Fpar[0] * my_fcn.Br2Nkk12;
      double Nphi12 = Fpar[1] * my_fcn.Br2Nphi12;
      double NFit12 = Nphi12/my_fcn.IBW[1];

      double ang = Fpar[2];
      double sig09 = Fpar[3];
      double sig12 = Fpar[4];

      // {mphi,gphi,sigma,A,F,Ang}
      const double pp09[] {
         my_fcn.mphi, my_fcn.gphi, sig09, my_fcn.ar, Ff, ang
      };
      const double pp12[] {
         my_fcn.mphi, my_fcn.gphi, sig12, my_fcn.ar, Ff, ang
      };

      double Nsb09 = Fpar[5];
      double Nsb12 = Fpar[6];

      // normalization on 1
      double Ntot09 = Nkk09 + Nsb09;
      NFit09 /= Ntot09;
      double Ntot12 = Nkk12 + Nsb12;
      NFit12 /= Ntot12;

      double NsbNorm09 = Nsb09 / my_fcn.normArsb[0];
      NsbNorm09 /= Ntot09;
      double NsbNorm12 = Nsb12 / my_fcn.normArsb[1];
      NsbNorm12 /= Ntot12;

      auto Lcr09 =
         [NFit09,pp09,NsbNorm09,arsb09](double x) -> double {
         return NFit09 * IntfrBWARG( x,pp09,0 ) +
            NsbNorm09 * RevArgus(x,arsb09);
      };
      auto Lcr12 =
         [NFit12,pp12,NsbNorm12,arsb12](double x) -> double {
         return NFit12 * IntfrBWARG( x,pp12,0 ) +
            NsbNorm12 * RevArgus(x,arsb12);
      };

      // linear extrapolation in [dL,dU] divided into  n points
      int n = 10000;
      vector<double> vLcr(n);
      double dm = (dU-dL)/(n-1);
      for ( int i = 0; i < n; ++i ) {
         double m = dL + i * dm;
         vLcr[i] = Lcr09(m);
      }
      auto Lcr = [&vLcr,dm,dL,dU](double x) -> double {
         if ( x < dL || x >= dU ) {
            return 0.;
         }
         int i = (x-dL) / dm;
         double mi = dL + i * dm;
         double f = vLcr[i] + (vLcr[i+1]-vLcr[i])*((x-mi)/dm);
         return f;
      };
      rmath_fun< decltype(Lcr) > fcr(Lcr);

      ROOT::Math::GoFTest* gofcr09 =
         new ROOT::Math::GoFTest( mkk09.size(),mkk09.data(),fcr,
               ROOT::Math::GoFTest::kPDF, dL,dU );
      pvalueKS09 = gofcr09 -> KolmogorovSmirnovTest();
      if ( !isBatch ) {
         cout << " pvalueKS09(cr)= " << pvalueKS09 << endl;
      }

      for ( int i = 0; i < n; ++i ) {
         double m = dL + i * dm;
         vLcr[i] = Lcr12(m);
      }
      ROOT::Math::GoFTest* gofcr12 =
         new ROOT::Math::GoFTest( mkk12.size(),mkk12.data(),fcr,
               ROOT::Math::GoFTest::kPDF, dL,dU );
      pvalueKS12 = gofcr12 -> KolmogorovSmirnovTest();
      if ( !isBatch ) {
         cout << " pvalueKS12(cr)= " << pvalueKS12 << endl;
      }

      // side-band
      double NsbNorm = 1. / my_fcn.normArsb[0]; // norm to 1 for Lsb
      double arsb = arsb09;
      auto Lsb = [&NsbNorm,&arsb](double x) -> double {
         return NsbNorm * RevArgus(x,arsb);
      };
      rmath_fun< decltype(Lsb) > fsb(Lsb);
      ROOT::Math::GoFTest* gofsb09 =
         new ROOT::Math::GoFTest( sb09.size(),sb09.data(),fsb,
            ROOT::Math::GoFTest::kPDF, dL,dU );
      pvalueKS09sb = gofsb09 -> KolmogorovSmirnovTest();
      if ( !isBatch ) {
         cout << " pvalueKS09(sb)= " << pvalueKS09sb << endl;
      }

      NsbNorm = 1. / my_fcn.normArsb[1];
      arsb = arsb12;
      ROOT::Math::GoFTest* gofsb12 =
         new ROOT::Math::GoFTest( sb12.size(),sb12.data(),fsb,
            ROOT::Math::GoFTest::kPDF, dL,dU );
      pvalueKS12sb = gofsb12 -> KolmogorovSmirnovTest();
      if ( !isBatch ) {
         cout << " pvalueKS12(sb)= " << pvalueKS12sb << endl;
      }
   }

   F_ret.push_back(pvalueKS09);
   F_ret.push_back(pvalueKS09sb);
   F_ret.push_back(pvalueKS12);
   F_ret.push_back(pvalueKS12sb);
   if ( isBatch ) {
      return F_ret;
   }

   //-----------------------------------------------------------------
   // Functions to draw
   double bW = h09 -> GetBinWidth(1); // bin width
   double bL = h09 -> GetBinLowEdge(1); // left boundary of hst

   auto Ldr09 =
      [Fpar,my_fcn,Ff,arsb09,bW](double* x,double* p) -> double {
      int isw = int(p[0]+0.5);

      double Nphi09 = Fpar[1] * my_fcn.Br2Nphi09;
      double NFit09 = Nphi09/my_fcn.IBW[0];
      double ang = Fpar[2];
      double sig09 = Fpar[3];
      const double pp[] {
         my_fcn.mphi,my_fcn.gphi,sig09,my_fcn.ar,Ff,ang
      };

      double BWARG = NFit09 * IntfrBWARG( x[0],pp,isw );
      if ( isw > 0 && isw < 4 ) { return bW * BWARG; }

      double NsbNorm09 = Fpar[5] / my_fcn.normArsb[0];
      double Bg = NsbNorm09 * RevArgus(x[0],arsb09);
      if ( isw == 4 ) { return bW * Bg; }

      return bW * (BWARG + Bg);
   };
   TF1* f09 = new TF1("f09", Ldr09, bL, dU, 1);
   f09 -> SetNpx(500);

   auto Ldr12 = [Fpar,my_fcn,Ff,arsb12,bW](double* x,double* p) -> double {
      int isw = int(p[0]+0.5);

      double Nphi12 = Fpar[1] * my_fcn.Br2Nphi12;
      double NFit12 = Nphi12/my_fcn.IBW[1];
      double ang = Fpar[2];
      double sig12 = Fpar[4];
      const double pp[] {
         my_fcn.mphi,my_fcn.gphi,sig12,my_fcn.ar,Ff,ang
      };

      double BWARG = NFit12 * IntfrBWARG( x[0],pp,isw );
      if ( isw > 0 && isw < 4 ) { return bW * BWARG; }

      double NsbNorm12 = Fpar[6] / my_fcn.normArsb[1];
      double Bg = NsbNorm12 * RevArgus(x[0],arsb12);
      if ( isw == 4 ) { return bW * Bg; }

      return bW * (BWARG + Bg);
   };
   TF1* f12 = new TF1("f12", Ldr12, bL, dU, 1);
   f12 -> SetNpx(500);

   // find max/min to draw
   f09 -> SetParameter(0, 1); // BW
   int mb09 = h09 -> GetMaximumBin();
   double hmax09 = h09 -> GetBinContent(mb09) +
      h09 -> GetBinError(mb09);
   double fmax09 = f09 -> GetMaximum( 1.01, 1.03, 1e-6);
   hmax09 = floor(1.1 * max( hmax09, fmax09 ));
   if ( hmax09 > 0 ) {
      h09 -> SetMaximum(hmax09);
   }

   f09 -> SetParameter(0, 3); // interference
   double hmin09 = f09 -> GetMinimum( 1.01, 1.03, 1e-6);
   hmin09 = floor(1.2 * min( 0., hmin09 ));
   if ( hmin09 < 0 ) {
      h09 -> SetMinimum(hmin09);
   }

   f12 -> SetParameter(0, 1); // BW
   int mb12 = h12 -> GetMaximumBin();
   double hmax12 = h12 -> GetBinContent(mb12) +
      h12 -> GetBinError(mb12);
   double fmax12 = f12 -> GetMaximum( 1.01, 1.03, 1e-6);
   hmax12 = floor(1.1 * max( hmax12, fmax12 ));
   if ( hmax12 > 0 ) {
      h12 -> SetMaximum(hmax12);
   }

   f12 -> SetParameter(0, 3); // interference
   double hmin12 = f12 -> GetMinimum( 1.01, 1.03, 1e-6);
   hmin12 = floor(1.2 * min( 0., hmin12 ));
   if ( hmin12 < 0 ) {
      h12 -> SetMinimum(hmin12);
   }

   //-----------------------------------------------------------------
   // Draw results
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000);
   c1 -> Divide(1,2);

   c1 -> cd(1);
   gPad -> SetGrid();

   SetHstFace(h09);
   h09 -> GetXaxis() -> SetTitleOffset(1.1);
   h09 -> GetYaxis() -> SetTitleOffset(1.);
   h09 -> SetLineWidth(2);
   h09 -> SetLineColor(kBlack);
   h09 -> SetMarkerStyle(20);

   h09 -> Draw("EP");

   TLegend* leg = new TLegend(0.65,0.53,0.99,0.89);
   leg -> SetTextSize(0.045);
   leg -> SetHeader("Toy MC:2009(top),2012(bot)","C");
   leg -> AddEntry( h09,"Toy MC","LEP" );

   f09 -> SetParameter(0, 0); // Sum
   f09 -> SetLineWidth(2);
   f09 -> SetLineStyle(kSolid);
   f09 -> SetLineColor(kRed+1);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Combined fit", "L" );

   f09 -> SetParameter(0, 1); // BW
   f09 -> SetLineWidth(2);
   f09 -> SetLineStyle(kDashed);
   f09 -> SetLineColor(kGreen+2);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Breit-Wigner #phi#eta", "L");

   f09 -> SetParameter(0, 2); // Argus
   f09 -> SetLineWidth(2);
   f09 -> SetLineStyle(kDashed);
   f09 -> SetLineColor(kBlue);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Non-#phi KK#eta", "L");

   f09 -> SetParameter(0, 3); // interference
   f09 -> SetLineWidth(2);
   f09 -> SetLineStyle(kDashed);
   f09 -> SetLineColor(kMagenta+1);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Interference", "L");
   leg -> Draw();

   TPaveText* pt = new TPaveText(0.65,0.32,0.99,0.52,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
   pt -> AddText( Form("Br(KK#eta)= %s #times10^{-4} ",
            PE.Eform(0,".2f",1e4)) );
   pt -> AddText( Form("Br(#phi#eta)= %s #times10^{-4}",
            PE.Eform(1,".2f",1e4)) );
   pt -> AddText( Form("#vartheta = %s",PE.Eform(2,".2f")) );
   pt -> Draw();
   gPad -> RedrawAxis();

   c1 -> cd(2);
   gPad -> SetGrid();

   SetHstFace(h12);
   h12 -> GetXaxis() -> SetTitleOffset(1.1);
   h12 -> GetYaxis() -> SetTitleOffset(1.);
   h12 -> SetLineWidth(2);
   h12 -> SetLineColor(kBlack);
   h12 -> SetMarkerStyle(20);

   h12 -> Draw("EP");

   f12 -> SetParameter(0, 0); // Sum
   f12 -> SetLineWidth(2);
   f12 -> SetLineColor(kRed+1);
   f12 -> DrawCopy("SAME");

   f12 -> SetParameter(0, 1); // BW
   f12 -> SetLineWidth(2);
   f12 -> SetLineStyle(kDashed);
   f12 -> SetLineColor(kGreen+2);
   f12 -> DrawCopy("SAME");

   f12 -> SetParameter(0, 2); // Argus
   f12 -> SetLineWidth(2);
   f12 -> SetLineStyle(kDashed);
   f12 -> SetLineColor(kBlue);
   f12 -> DrawCopy("SAME");

   f12 -> SetParameter(0, 3); // interference
   f12 -> SetLineWidth(2);
   f12 -> SetLineStyle(kDashed);
   f12 -> SetLineColor(kMagenta+1);
   f12 -> DrawCopy("SAME");

   TPaveText* pt09 = new TPaveText(0.65,0.69,0.99,0.99,"NDC");
   pt09 -> SetTextAlign(12);
   pt09 -> SetTextFont(42);
   pt09 -> AddText( Form("#it{p-value(2009) = %.3f}",pvalueKS09) );
   pt09 -> AddText( Form("#sigma(2009) = %s MeV",
            PE.Eform(3,".2f",1e3)) );
   pt09 -> AddText( Form("#lower[-0.1]{Nbg(2009) = %s}",
            PE.Eform(5,".1f")) );
   pt09 -> AddText( Form("a(sb09)= %s", PE.Eform(7,".1f")) );
   pt09 -> Draw();

   TPaveText* pt12 = new TPaveText(0.65,0.39,0.99,0.69,"NDC");
   pt12 -> SetTextAlign(12);
   pt12 -> SetTextFont(42);
   pt12 -> AddText( Form("#it{p-value(2012) = %.3f}", pvalueKS12) );
   pt12 -> AddText( Form("#sigma(2012) = %s MeV",
            PE.Eform(4,".2f",1e3)) );
   pt12 -> AddText( Form("#lower[-0.1]{Nbg(2012) = %s}",
            PE.Eform(6,".1f")) );
   pt12 -> AddText( Form("a(sb12)= %s", PE.Eform(8,".1f")) );
   pt12 -> Draw();

   gPad -> RedrawAxis();

   c1 -> Update();
   c1 -> Print(pdf.c_str());

   return F_ret;
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

   TCanvas* c1 = new TCanvas("c1","...",0,0,1300,1000);
   c1 -> Divide(2,2);

   TF1* pl0 = (TF1*)gROOT -> GetFunction("pol0") -> Clone();
   pl0 -> SetLineColor(kRed);
   pl0 -> SetLineWidth(2);
   pl0 -> SetLineStyle(kDashed);

   for ( int i = 0; i < 3; i+=2 ) {
      c1 -> cd(1+i);
      gPad -> SetGrid();

      SetHstFace(hst[i]);
      hst[i] -> GetXaxis() -> SetTitleOffset(1.1);
      hst[i] -> GetYaxis() -> SetTitleOffset(1.3);
      hst[i] -> Draw("E");

      c1 -> cd(2+i);
      gPad -> SetGrid();

      SetHstFace(hst[1+i]);
      hst[1+i] -> GetXaxis() -> SetTitleOffset(1.1);
      hst[1+i] -> Fit(pl0,"L","E");
   }

   c1 -> Update();
   c1 -> Print( pdf.c_str() );
}

//----------------------------------------------------------------------
void ToyMC_fit(int Ntoys, string file_name) {
//----------------------------------------------------------------------
   bool isBatch = Ntoys > 1;

   string pdf = ( (isBatch) ? "" : file_name );
   TFile* c_out = nullptr;
   TNtupleD* tuple = nullptr;

   if ( isBatch ) {
      // Histograms
      c_out = new TFile(file_name.c_str(),"recreate");
      if ( !c_out -> IsOpen() ) {
         printf(" can not open %s file\n",file_name.c_str());
         exit(1);
      }
      c_out -> cd();
      tuple = new TNtupleD("tmc","Toy MC",
            "bkk:bphi:ang:sig09:sig12:nbg09:nbg12:ar09:ar12:"
            "ubkk:ubphi:uang:usig09:usig12:unbg09:unbg12:uar09:uar12:"
            "lbkk:lbphi:lang:lsig09:lsig12:lnbg09:lnbg12:lar09:lar12:"
            "pbkk:pbphi:pang:psig09:psig12:pnbg09:pnbg12:par09:par12:"
            "Lmin:st:pv09:pv09sb:pv12:pv12sb:"
            "nkk09g:nbg09g:nsb09g:"
            "nkk12g:nbg12g:nsb12g"
            );
   }

   // Set parameters
   double BrKKeta  = 4.5e-4;
   double Brphieta = 8.5e-4;

   // Nj = N(Jpsi)*eff*Breta:
   double Nj09 = 2.02e6; // 2009 : 2.02e6*4.5e-4 = 909
   double Nj12 = 6.14e6; // 2012 : 6.14e6*4.5e-4 = 2763

   double sig09 = 1.4e-3;
   double sig12 = 1.1e-3;

   double nbg09 = 17;
   double nbg12 = 54;
   double arsb09 = 8.;
   double arsb12 = 5.;

   Gpar* gp12 = new Gpar(BrKKeta, Brphieta,
         sig12, Nj12, nbg12, arsb12 );
   gp12 -> Print("2012");
   Gpar* gp09 = new Gpar(BrKKeta, Brphieta,
         sig09, Nj09, nbg09, arsb09, gp12 -> F);
   gp09 -> Print("2009");

   TStopwatch t;
   TH1D* hst[4] {nullptr, nullptr, nullptr, nullptr};
   for ( int itoy = 0; itoy < Ntoys; ++itoy ) {
      t.Start();

      myFCN_toy my_fcn(Nj09,Nj12);  // class for 'FitFCN'

      my_fcn.mkk09 = gp09 -> signal_ToyMC();
      vector<double> mkk_bg09 = gp09 -> bg_ToyMC();
      my_fcn.mkk09.insert( my_fcn.mkk09.end(),
            mkk_bg09.begin(),mkk_bg09.end() );
      delete hst[0];
      hst[0] = get_hst( my_fcn.mkk09, "toy_sig09" );
      my_fcn.sb09 = gp09 -> bg_ToyMC(); // side-band 09
      delete hst[1];
      hst[1] = get_hst( my_fcn.sb09, "toy_sb09" );

      my_fcn.mkk12 = gp12 -> signal_ToyMC();
      vector<double> mkk_bg12 = gp12 -> bg_ToyMC();
      my_fcn.mkk12.insert( my_fcn.mkk12.end(),
            mkk_bg12.begin(),mkk_bg12.end() );
      delete hst[2];
      hst[2] = get_hst( my_fcn.mkk12, "toy_sig12" );
      my_fcn.sb12 = gp12 -> bg_ToyMC(); // side-band 12
      delete hst[3];
      hst[3] = get_hst( my_fcn.sb12, "toy_sb12" );

      if ( !isBatch ) {
         printf(" nkk09= %zu, nbg09= %zu, nsb09= %zu\n",
           my_fcn.mkk09.size(), mkk_bg09.size(), my_fcn.sb09.size());
         printf(" nkk12= %zu, nbg12= %zu, nsb12= %zu\n",
           my_fcn.mkk12.size(), mkk_bg12.size(), my_fcn.sb12.size());
      }
//       plot_hst(hst,pdf); // for debug

      vector<double> ret;
      ret = do_fit_toy(my_fcn,hst,pdf);

      if ( isBatch ) {
         // add generated numbers
         ret.push_back( double(my_fcn.mkk09.size()) );
         ret.push_back( double(mkk_bg09.size()) );
         ret.push_back( double(my_fcn.sb09.size()) );
         ret.push_back( double(my_fcn.mkk12.size()) );
         ret.push_back( double(mkk_bg12.size()) );
         ret.push_back( double(my_fcn.sb12.size()) );
         tuple -> Fill( ret.data() );
      }

      t.Stop();
      printf(" end of MCtoy# %i: ",itoy+1);
      t.Print();
   }

   if ( isBatch ) { // save in root file
      tuple -> Write();
      c_out -> Close();
   }
}

// {{{1 MAIN for interpreter
#ifndef BATCH
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

   // set random number generator:
   ULong_t Iseed = (ULong_t)time(NULL);
   printf("   ULong_t Iseed= %luL;\n",Iseed);
   gRandom -> SetSeed(Iseed);

   // ------------- ToyMC ---------------
   ToyMC_fit(1, "ToyMC_cf_1.pdf" ); // must be Ntoys=1
}
#endif

// {{{1 MAIN for batch mode:
#ifdef BATCH
//----------------------------------------------------------------------
int main(int argc, char* argv[]) {
//----------------------------------------------------------------------
   ULong_t Iseed = 0;
   int Ntoys = 10;
   string hst_file("ToyMC_cf.root");

   //-------------------------------------------------------------------
   // => getopt()
   bool is_error = false;
   int oc; // option
   while( (oc = getopt(argc,argv,":n:h:s:")) != -1 ) {
      switch( oc ) {

      case 'h':  // change default file name for histograms
         hst_file = string(optarg);
         break;

      case 's':  // change seed for random generator
         Iseed = ULong_t(atoi(optarg));
         break;

      case 'n':  // number of toy MC generated
         Ntoys = atoi(optarg);
         break;

      case ':':  // no argument (first character of optstring MUST be ':')
         is_error = true;
         printf(" option `-%c' requires an argument\n",optopt);
         break;

      case '?':  // errors
      default:
         is_error = true;
         printf(" invalid option `-%c'\n",optopt);
         break;
      }
   }

   if( is_error ) {
      cout << endl;
      cout << " Usage: " << argv[0] << " -h hst_file.root "
           << " -s SEED(rndm)" << endl;
      exit(0);
   }

   //-------------------------------------------------------------------
   // set random number generator: TRandom3 by default
   if ( !Iseed ) {
      printf(" invalid SEED `-%c'\n",optopt);
      exit(0);
   }

   printf(" -- getopt parameters --\n");
   printf(" Ntoys= %i\n", Ntoys);
   printf(" Set gRrandom with ULong_t Iseed= %luL;\n",Iseed);
   printf(" output root file: %s\n", hst_file.c_str());
   printf(" -- end getopt parameters --\n");

//    return 0; // test getopt

   gRandom -> SetSeed(Iseed);

   //-------------------------------------------------------------------
   // = define GSL error handler which does nothing =
   gsl_set_error_handler_off();

   // set integrator: ROOT::Math::GSLIntegrator adaptive method (QAG)
   ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Adaptive");

   //-------------------------------------------------------------------
   TTimeStamp ts;
   cout << " Start time: " << ts.AsString() << endl;
   ToyMC_fit(Ntoys,hst_file);
   ts.Set(); // update
   cout << " Stop time: " << ts.AsString() << endl;
}
#endif
