// mass_kk_fit.cc
// unbinned LH fit of M(K+K-) distributions
// Formulas for:
// *) Breit-Wigner for psi -> K+K- (psi from J/Psi -> phi eta)
// *) convolution of BW with the Gauss distribution
// *) ARGUS function (and "reverse" Argus)
// *) interference
//    -> mkkYY_fit??.pdf
//    -> mkk_cf_???.pdf (combined fit)

#include <algorithm>
#include <iterator>
#include <functional>        // for std::function
#include "gsl/gsl_errno.h"   // GSL error handler

#include "masses.h"

// {{{1 constants and helper functions
//--------------------------------------------------------------------
// {{{2 GLOBAL constants
// name of folder with root files
string Dir;

// the meson radial parameter in Blatt-Weisskopf form factor
static const double R_BW = 3.0; // GeV^{-1}; +/-1 uncertainties study

// the left cutoff = 0.987354
static const double dL = 2*Mk;

// binning for Mkk: bin width 1.0 MeV
static const int Nbins = 100;
static const double bL = 0.98; // first bin < dL
static const double dU = 1.08; // MUST BE < 1.0835 !!!
static const double bW = (dU-bL)/Nbins; // bin width

// {{{2 SQ() and SetHstFace
//--------------------------------------------------------------------
constexpr double SQ(double x)
//--------------------------------------------------------------------
{
   return x*x;
}

//--------------------------------------------------------------------
void SetHstFace(TH1* hst)
//--------------------------------------------------------------------
{
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

// {{{2 Handler for GSL exceptions
//--------------------------------------------------------------------
void my_gslhandler(const char* reason, const char* file, int line,
      int gsl_errno)
//--------------------------------------------------------------------
{
   static size_t Nerr = 0;
   const size_t maxNerr = 100;
   // const size_t maxNerr = 10'000;
   // printf("-- GSL: %s:%d, error %s, reason %s\n",
         // file,line,gsl_strerror(gsl_errno), reason);
   if ( reason ) {
      printf("-- GSL error in %s, reason: %s\n", file,reason);
   } else {
      printf("-- GSL total number of errors is %zu\n",Nerr);
   }

   if ( ++Nerr > maxNerr ) {
      printf("-- GSL FATAL: number of errors is %zu\n",Nerr);
      exit(EXIT_FAILURE);
   }
   return;
}

// {{{2 Adapter: rmath_fun< decltype(Lambda) > Func(Lambda);
//--------------------------------------------------------------------
// adapter class to use lambda functions with closures in ROOT::MATH
//--------------------------------------------------------------------
template< typename F >
class rmath_fun: public ROOT::Math::IBaseFunctionOneDim {
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

// {{{2 Adapter class for handling parameters and errors
//--------------------------------------------------------------------
struct ParAndErr {
   vector<double> Fpar; // parameters
   vector<double> Perr; // symmetric errors
   vector<double> Uerr; // upper minos errors
   vector<double> Lerr; // lower minos errors

   // if the upper and lower errors differ no more than the
   // 'threshold', the maximum of them is used as a symmetric error
   double threshold = 0.1;

   // ctor
   //-----------------------------------------------------------------
   ParAndErr(const ROOT::Fit::FitResult& res,double thr = 0.1)
   //-----------------------------------------------------------------
   {
      Fpar = res.Parameters();
      Perr = res.Errors();

      int Npar = res.NPar();
      Lerr.resize(Npar,0.);
      Uerr.resize(Npar,0.);

      threshold = thr;

      for ( int np = 0; np < Npar; ++np ) {
         if ( res.IsParameterFixed(np) ) {
            Perr[np] = 0;
            continue;
         }
         if ( res.HasMinosError(np) ) {
            Uerr[np] = fabs(res.UpperError(np));
            Lerr[np] = fabs(res.LowerError(np));
         }
      }
   }

   // pretty print parameters with errors
   //-----------------------------------------------------------------
   const char* Eform (int i, string fm, double X = 1)
   //-----------------------------------------------------------------
   {
      auto SymF = [fm,X](double par, double err) {
         string fmt = "%" +fm+ " #pm %" +fm;
         return Form(fmt.c_str(),par*X,err*X);
      };

      auto FixF = [fm,X](double par) {
         string fmt = "%" +fm+ " (fixed)";
         return Form(fmt.c_str(),par*X);
      };

      auto AsymF = [fm,X](double par, double uerr, double lerr) {
         string fmt = "%" + fm
            + "^{#lower[0.15]{#kern[0.20]{#plus}#kern[0.05]{%"
            + fm + "}}" + "}"
            + "_{#lower[-0.15]{#kern[0.20]{#minus}#kern[0.05]{%"
            + fm + "}}" + "}";
         return Form(fmt.c_str(),par*X,uerr*X,lerr*X);
      };

      if ( Perr[i] == 0 ) {
         return FixF(Fpar[i]);
      }

      if ( Uerr[i]*Lerr[i] != 0 ) {
         // there are both upper and lower errors
         if ( fabs(1-Lerr[i]/Uerr[i]) < threshold ) {
            double err = max(Lerr[i],Uerr[i]);
            return SymF(Fpar[i],err);
         } else {
            return AsymF(Fpar[i],Uerr[i],Lerr[i]);
         }
      }

      return SymF(Fpar[i],Perr[i]);
   }

   // return middle point = par(i) + 0.5(UpperError+LowerError)
   //-----------------------------------------------------------------
   double Middle (int i)
   //-----------------------------------------------------------------
   {
      return Fpar[i] + 0.5*(Uerr[i]-Lerr[i]);
   }
};

//--------------------------------------------------------------------
struct ValErr {
   double val = 0.;
   double err = 0.;
   const char* prt(string fm) const
   {
      string fmt = "%" + fm + " +/- %" + fm;
      return Form(fmt.c_str(),val,err);
   }
};
//--------------------------------------------------------------------

// {{{2 Adapter function to spead up the Kolmogorov-Smirnov test
//--------------------------------------------------------------------
double FastKSTest(const function<double (double)>& Fun,
      vector<double>& Mkk)
//--------------------------------------------------------------------
{
   // NOTE: Mkk vector is sorted here!
   // Fun is kind of PDF function for Mkk-data;
   size_t n = 10000;
   double dx = (dU-dL)/(n-1);
   vector<double> vec(n,0.);
   // CDF: Integral_dL^{x_i} Fun(x) dx, trapezoidal rule
   double fo = 0.;
   for ( size_t i = 1; i < n; ++i ) {
      double x = dL + i*dx;
      double fx = Fun(x);
      vec[i] = vec[i-1] + (fo+fx); // 0.5*dx is common for all terms
      fo = fx;
   }
   // normalise CDF on one
   double norm = 1. / vec[n-1];
   for ( size_t i = 1; i < n; ++i ) {
      vec[i] *= norm;
   }
   // linear extrapolation of CDF on a uniform grid into [dL,dU]
   auto Lcdf = [&vec,dx](double x) {
      if ( x <= dL ) {
         return 0.;
      }
      if ( x >= dU ) {
         return 1.;
      }
      int i = (x-dL) / dx;
      double xi = dL + i * dx;
      double f = vec[i] + (vec[i+1]-vec[i])*((x-xi)/dx);
      return f;
   };

   // Kolmogorov-Smirnov test: https://numerical.recipes/book.html
   sort(begin(Mkk),end(Mkk)); // sort in ascending order
   size_t nmkk = Mkk.size();
   double en = nmkk;
   fo = 0.;
   double d = 0.;
   for ( size_t j = 0; j < nmkk; ++j ) {
      double fn = (j+1)/en;
      double ff = Lcdf(Mkk[j]);
      double dt = max( fabs(fo-ff), fabs(fn-ff) );
      if ( dt > d ) {
         d = dt;
         // debug print
         //  printf("dt= %g, m=%g, ff=%g, ff=%g, fo=%g\n",
               //  dt,Mkk[j],ff,fn,fo);
      }
      fo = fn;
   }
   en = sqrt(en);
   double pval = TMath::KolmogorovProb( (en+0.12+0.11/en)*d );
   //  printf("K-S test debug print: pval= %g, d= %g\n", pval,d);

   return pval;
}

//--------------------------------------------------------------------
double FastKSTestOLD(const function<double (double)>& Fun,
      const vector<double>& Mkk, size_t NMkk = 0)
//--------------------------------------------------------------------
{
   // linear extrapolation of Fun on a uniform grid of n points in the
   // region [dL,dU] is used
   size_t n = 10000;
   vector<double> vec(n);
   double dx = (dU-dL)/(n-1);
   for ( size_t i = 0; i < n; ++i ) {
      double x = dL + i * dx;
      vec[i] = Fun(x);
   }
   auto Lin = [&vec,dx](double x) -> double {
      if ( x < dL || x >= dU ) {
         return 0.;
      }
      int i = (x-dL) / dx;
      double xi = dL + i * dx;
      double f = vec[i] + (vec[i+1]-vec[i])*((x-xi)/dx);
      return f;
   };
   rmath_fun< decltype(Lin) > ftest(Lin); // lambda -> Root1Dfunc
   if ( NMkk == 0 ) {
      NMkk = Mkk.size();
   } else {
      NMkk = min(NMkk, Mkk.size());
   }
   ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(
         NMkk,Mkk.data(),ftest,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pval = 0, Dn = 0;
   goftest->KolmogorovSmirnovTest(pval, Dn);
   printf("K-S test debug print: pval= %g, Dn= %g\n", pval,Dn);

   return pval;
   // double pvalueAD = gofcr->AndersonDarlingTest();
   // cout << " pvalueAD= " << pvalueAD << endl;
}

// {{{2 Function to calculate chi squared and ndof
//--------------------------------------------------------------------
tuple<double,int> GetCh2Ndof(TH1D* hst, TF1* fun, double norm=1.)
      // const function<double (double)>& Fun)
//--------------------------------------------------------------------
{
   double ch2 = 0;
   int ndof = 0;
   int nbins = hst->GetNbinsX();
   for ( int ib = 1; ib <= nbins; ++ib ) {
      double x = hst->GetBinCenter(ib);
      if ( x < dL || x > dU ) { // mkk range
         continue;
      }
      double y = hst->GetBinContent(ib);
      if ( y < 1e-6 ) { // skip empty bins
         continue;
      }
      double e = hst->GetBinError(ib);
      double f = norm*fun->Integral(x-bW/2,x+bW/2);
      // double f = norm*bW*fun->Eval(x);
      // double f = norm*bW*Fun(x);
      double d =   SQ((y-f)/e);
      // debug print
      //  printf("x= %.4f, ch2= %.1f, y= %.1f f=%.1f e= %.1f\n",
            //  x,d,y,f,e);
      ch2 += d;
      ndof += 1;
   }
   return make_tuple(ch2,ndof);
}

// {{{1 data processing
//--------------------------------------------------------------------
vector<double> get_mkk_vec(string fname, int type)
//--------------------------------------------------------------------
{
#include "cuts.h"
   bool isMC = (fname.find("mcsig") != string::npos);

   fname = Dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   TCut c_here = c_Mrec+c_chi2;
   // c_here += TCut("chsq3g>ch2");// with the search for third gamma
   c_here += TCut( Form("%f<=Mkk&&Mkk<%f",dL,dU) );
   if ( type == 0 ) {        // central part
      c_here += c_cpgg;
   } else if ( type == 1 ) { // side-band
      c_here += c_sbgg;
   }
   // if ( isMC ) {
      // cout << " THIS IS MONTE CARLO EVENTS" << endl;
      // c_here += c_MCmkk;
   // }

   int n = a4c->Draw("Mkk", c_here, "goff");
   // cout << " n= " << n << endl;
   double* buffer = a4c->GetVal(0);
   vector<double> mkk(buffer, buffer+n);

   return mkk;
}

//--------------------------------------------------------------------
TH1D* get_hst( const vector<double>& mkk, string hname )
//--------------------------------------------------------------------
{
   string title(";M^{ inv}_{ K^{#plus}K^{#minus }}, GeV/c^{2}");
   int ibW = int(bW*1e5+0.1);
   title += string(";Entries / ") +
      ( (ibW%10 == 0)
        ? string(Form("%.0f MeV/c^{2}",bW*1e3))
        : string(Form("%.2f MeV/c^{2}",bW*1e3)) );
   TH1D* hst = new TH1D(hname.c_str(),title.c_str(),Nbins,bL,dU);

   hst->Sumw2(true);
   for ( const auto& mk : mkk ) {
      hst->Fill(mk);
   }

   return hst;
}

//--------------------------------------------------------------------
void ReoderMkk(vector<double>& mkk)
//--------------------------------------------------------------------
{
   // reorder the mkk vector so that the likelihood summation goes
   // from smaller values to larger ones
   auto it = partition( begin(mkk), end(mkk),
         [](double m){return (m < 1.016 || 1.024 < m);} );
   // sort in ascending order
   sort( begin(mkk), it );
   // sort by distance from 1.02 in descending order
   sort( it, end(mkk),
         [](double x,double y) {return abs(1.02-x)>abs(1.02-y);} );
   // debug print
   // int nit = distance(begin(mkk),it);
   // cout << "nit=" << nit << endl;
   // int n = mkk.size();
   // for (int i = 0; i < n; ++i ) {
      // if ( i < 5 || abs(i-nit) < 5 || i > n-5 ) {
         // printf("i= %d  mkk= %.9f\n", i,mkk[i]);
      // }
   // }
}

//--------------------------------------------------------------------
TH1D* get_mcMkk(int date)
//--------------------------------------------------------------------
{
   string fname( Form("mcsig_kkmc_%02i.root",date%100) );
   fname = Dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( !froot ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TH1D* hmcmkk = (TH1D*) gDirectory->Get("mc_Mkk");
   if( !hmcmkk ) {
      cerr << "can not find mc_Mkk" << endl;
      exit(EXIT_FAILURE);
   }

   string title(";M^{ inv}_{ K^{#plus}K^{#minus }}, GeV/c^{2}");
   title += string(";Entries / 1 MeV/c^{2}");
   hmcmkk->SetTitle(title.c_str());
   hmcmkk->Sumw2(true);

   return hmcmkk;
}

// {{{1 1.Breigt Wigner for phi -> KK
//--------------------------------------------------------------------
// {{{2 Functions
//--------------------------------------------------------------------
double BreitWigner(double m, double mphi, double gphi, double R=R_BW)
//--------------------------------------------------------------------
{
   // we assume that mphi, gphi and R are the fit parameters
   if ( m < dL ) { // phase space becomes zero (see p_k)
      return 0;
   }

   constexpr double Mjpsi2 = SQ(Mjpsi);
   constexpr double Meta2  = SQ(Meta);
   constexpr double Mk2    = SQ(Mk);

   double m2    = SQ(m);
   double mphi2 = SQ(mphi);
   double gphi2 = SQ(gphi);

   double gam = mphi*sqrt(mphi2+gphi2);
   double kap = 2*M_SQRT2*mphi*gphi*gam / (M_PI*sqrt(mphi2+gam));

   double p_phi0 = sqrt(SQ(Mjpsi2-mphi2-Meta2)-4*mphi2*Meta2)
      / (2*Mjpsi);
   double p_phi  = sqrt(SQ(Mjpsi2-m2-Meta2)-4*m2*Meta2) / (2*Mjpsi);
   double r_phi  = p_phi / p_phi0;

   double p_k0 = 0.5*sqrt(mphi2-4*Mk2);
   double p_k  = 0.5*sqrt(m2-4*Mk2);
   double r_k  = p_k / p_k0;

   // B(m)/B(m0)
   double BB_phi = sqrt( (1+SQ(R*p_phi0)) / (1+SQ(R*p_phi)) );
   double BB_k2  = (1+SQ(R*p_k0)) / (1+SQ(R*p_k));
   double BB_k = sqrt( BB_k2 );

   double fm = r_phi*r_k*BB_phi*BB_k;
   double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k2; // == m*G(m)

   return (kap * SQ(fm)) / ( SQ(m2-mphi2) + SQ(GM) );
}

//--------------------------------------------------------------------
double BreitWignerGauss( double m,
      double mphi, double gphi, // BW parameters
      double sigma,             // Gauss
      double slope              // eff(m)
      )
//--------------------------------------------------------------------
// Function Fun(x) folded with a normal distribution:
//
//        1                       [                  (X'-X)^2    ]
// ---------------- *  Integral   [ Fun(X') * exp( - --------- ) ] dX'
// sqrt(2*pi)*sigma   (-oo<X'<+oo)[                  2*sigma^2   ]
//
//--------------------------------------------------------------------
// sigma ~ 1/3 gphi (bes3-exp)
// integration in +/-5 sigma interval
//--------------------------------------------------------------------
{
   static const double Ngauss = 1./sqrt(2*M_PI);

   // integrand lambda function
   double pp[] { m, mphi, gphi, sigma, slope };
   auto Lint = [](double t, void* pp) { // t is (X'-X)/sigma
      double* p = static_cast<double*>(pp);
      double x = p[0] + t*p[3];
      return exp(-0.5*t*t) *
         BreitWigner(x,p[1],p[2]) * (1+p[4]*(x-1.02));
   };

   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-12, 1.e-9, 1000);
   double result = gsl_int.Integral(Lint,pp,-5.,+5.);

   if ( gsl_int.Status() ) { // true for error
      double abs_err = gsl_int.Error();
      printf("%s error for m=%.3f: result=%g, error=%g\n",
            __func__,m,result,abs_err);
      printf("           mphi= %g, gphi=%g, sig=%g, sl= %g\n",
            mphi,gphi,sigma,slope);
   }

   return Ngauss * result;
}

//--------------------------------------------------------------------
double BreitWignerGaussN( double m,
      double mphi, double gphi, double sigma,
      double slope
      )
//--------------------------------------------------------------------
{
   // Numeric normalisation on one for range [dL,dU]

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
      ROOT::Math::GSLIntegrator gsl_int(1.e-12, 1.e-9, 1000);
      norm = gsl_int.Integral(Lint,p,dL,dU);

      if ( gsl_int.Status() ) { // true for error
         double abs_err = gsl_int.Error();
         printf("%s error: norm=%g, error=%g\n",
               __func__,norm,abs_err);
      }

      cacheN = norm;
      cacheM = mphi;
      cacheG = gphi;
      cacheS = sigma;
      cacheSl= slope;
   }

   return BreitWignerGauss(m,mphi,gphi,sigma,slope) / norm;
}

//--------------------------------------------------------------------
tuple<double,double> BreitWignerGaussN( double m,
      double mphi, double gphi, double sig1, double sig2,
      double slope
      )
//--------------------------------------------------------------------
{
   // Normalisation on one: case of two sigma

   double norm1 = 0, norm2 = 0;
   // cash parameters
   static double cacheN1 = 0;
   static double cacheN2 = 0;
   static double cacheM = 0;
   static double cacheG = 0;
   static double cacheS1 = 0;
   static double cacheS2 = 0;
   static double cacheSl = 0;
   if ( cacheN1 > 0 && mphi == cacheM && gphi == cacheG
         && sig1 == cacheS1 && sig2 == cacheS2
         && slope == cacheSl ) {
      norm1 = cacheN1;
      norm2 = cacheN2;
   } else {
      // integrand lambda function
      double p[] = {mphi,gphi,sig1,slope};
      auto Lint = [](double x, void* pp) -> double{
         double* p = static_cast<double*>(pp);
         return BreitWignerGauss(x,p[0],p[1],p[2],p[3]);
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-8, 1000);
      norm1 = gsl_int.Integral(Lint,p,dL,dU);
      p[2] = sig2;
      norm2 = gsl_int.Integral(Lint,p,dL,dU);

      cacheN1 = norm1;
      cacheN2 = norm2;
      cacheM = mphi;
      cacheG = gphi;
      cacheS1 = sig1;
      cacheS2 = sig2;
      cacheSl= slope;
   }

   return make_tuple(
         BreitWignerGauss(m,mphi,gphi,sig1,slope)/norm1,
         BreitWignerGauss(m,mphi,gphi,sig2,slope)/norm2
         );
}

// {{{2 Tests
//--------------------------------------------------------------------
void test_BW(size_t Cx, size_t Cy)
//--------------------------------------------------------------------
{
   // BreitWigner (internally normalized to 1)
   auto Lbw = [](const double* x,const double* p) -> double {
      return p[0]*BreitWigner(x[0],p[1],p[2],p[3]);
   };
   TF1* bw = new TF1("bw", Lbw, dL, dU, 4);
   bw->SetParNames("Norm","Mphi","Gphi","Rff");
   bw->SetParameters(1., Mphi, Gphi,R_BW);
   bw->SetLineWidth(1);
   bw->SetLineColor(kBlack);
   bw->SetNpx(500);

   double norm = 1. / bw->Integral(dL,dU,1e-8);
   printf("norm = %.7f\n",norm);
   bw->SetParameter(0, norm );

   TCanvas* c1 = new TCanvas("c1_BW","...",0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(true);

   TLegend* leg = new TLegend(0.54,0.69,0.892,0.89);
   leg->SetHeader("Breit-Wigner (M#phi, #Gamma#phi)","C");

   bw->DrawCopy("L");
   c1->Update();
   leg->AddEntry(bw->Clone(), "Rff=3GeV^{-1}", "L");

   bw->Update();
   bw->SetParameters(1., Mphi, Gphi, 6.); // Blatt-Weisskopf x2
   bw->SetLineColor(kBlue);
   bw->DrawCopy("LSAME");
   leg->AddEntry(bw->Clone(), "Rff=6GeV^{-1}", "L");

   bw->Update();
   bw->SetParameters(1., Mphi, Gphi, 1.5); // Blatt-Weisskopf /2
   bw->SetLineColor(kRed);
   bw->DrawCopy("LSAME");
   leg->AddEntry(bw->Clone(), "Rff=1.5GeV^{-1}", "L");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
}

//--------------------------------------------------------------------
void test_BWGN(size_t Cx, size_t Cy)
//--------------------------------------------------------------------
{
   // BreitWignerGaussN

   auto Lbwgn = [](const double* x,const double* p) -> double {
      return p[0]*BreitWignerGaussN(x[0],p[1],p[2],p[3],p[4]);
   };
   TF1* bwgn = new TF1("bwgn", Lbwgn, dL, dU, 5);
   bwgn->SetParNames("Norm","Mphi","Gphi","Sigma","Slope");
   bwgn->SetParameters(1., Mphi, Gphi, 1.e-3, 0.);
   bwgn->SetLineWidth(2);
   bwgn->SetLineColor(kRed);
   bwgn->SetLineStyle(kDashed);
   bwgn->SetNpx(500);

   TCanvas* c1 = new TCanvas("c1_BWGN","...",Cx/2,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(true);

   TLegend* leg = new TLegend(0.49,0.69,0.892,0.89);
   leg->SetHeader("BW(M#phi,#Gamma#phi) (+) Gauss(#sigma)","C");

   bwgn->DrawCopy("L");
   leg->AddEntry( bwgn->Clone(), Form("#sigma=%.2f sl=%.2f",
            bwgn->GetParameter(3)*1e3,
            bwgn->GetParameter(4) ), "L" );

   bwgn->Update();
   bwgn->SetParameter(4, -2.0 );  // slope
   bwgn->SetLineColor(kBlue);
   bwgn->DrawCopy("LSAME");
   leg->AddEntry(bwgn->Clone(), Form("#sigma=%.2f sl=%.2f",
            bwgn->GetParameter(3)*1e3,
            bwgn->GetParameter(4) ), "L");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
}

// {{{2 Fit Mkk(MC-truth) histogram by Breit-Wigner
//--------------------------------------------------------------------
void sig_fit_mc(int date, string pdf, size_t Cx, size_t Cy)
//--------------------------------------------------------------------
{
   TH1D* hst = nullptr;
   if (date > 0) {
      hst = get_mcMkk(date);
   } else {
      // auto h = get_mcMkk(2009);
      // hst = (TH1D*)h->Clone("mc_mkk_sum");
      hst = get_mcMkk(2009);
      hst->Add( get_mcMkk(2012) );
      hst->Add( get_mcMkk(2021) );
   }

   double bin_width = hst->GetBinWidth(1);
   auto Lbwg = [bin_width](const double* x,const double* p) {
      return bin_width * 1e6*p[0] *
         BreitWigner(x[0], 1e-3*p[1], 1e-3*p[2], p[3]);
      // return bin_width*p[0]*BreitWigner(x[0],p[1],p[2],p[3]);
   };
   TF1* bwg = new TF1(Form("bwg_%i",date), Lbwg, 0.98, 2., 4);
   // bwg->SetParNames("Norm","Mphi","Gphi","r");
   // bwg->SetParameters(hst->Integral(), Mphi, Gphi, R_BW);
   bwg->SetParNames("Norm, 10^{6}",
         "M_{#phi}, MeV",
         "#Gamma_{#phi}, MeV",
         "#scale[1.2]{r}, GeV^{-1}");
   bwg->SetParameters(hst->Integral()*1e-6, Mphi*1e3, Gphi*1e3, R_BW);

   // bwg->FixParameter(1, 1.01946e3); // BesEventGen v709
   // bwg->FixParameter(2, 0.00426e3); // BesEventGen v709
   // bwg->FixParameter(3, R_BW);
   // bwg->SetParLimits(3, 0.5, 5.);

   bwg->SetLineWidth(2);
   bwg->SetLineColor(kRed);
   bwg->SetNpx(500);

   // gStyle->SetOptFit(111); // do not print fixed parameters
   gStyle->SetOptFit(112); // print all parameters
   gStyle->SetStatFont(42);
   // gStyle->SetFitFormat(".7g");
   gStyle->SetFitFormat(".3f");
   gStyle->SetStatX(0.892);
   gStyle->SetStatY(0.89);

   auto name = Form("c1_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   SetHstFace(hst);
   hst->SetMinimum(1.);
   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetTitleOffset(1.1);
   hst->SetLineWidth(3);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);
   hst->SetMarkerSize(0.8);

   hst->Draw("EP");
   double maxMkk=1.08336;
   int bin = hst->FindBin(maxMkk);
   double maxMkk_fit = hst->GetXaxis()->GetBinCenter(bin)-bin_width/2;
   hst->Fit(bwg,"IE","",dL,maxMkk_fit);

   TLine* lB = new TLine;
   lB->SetLineColor(kBlue+1);
   lB->SetLineWidth(2);
   lB->SetLineStyle(kDashed);
   lB->DrawLine(maxMkk,1,maxMkk,bwg->Eval(maxMkk));

   TF1* bwg1 = (TF1*) bwg->Clone();
   bwg1->SetLineStyle(kDashed);
   bwg1->DrawF1(maxMkk_fit,1.125,"SAME");

   double genEv = hst->GetEntries();
   double corr = 1e6*bwg->GetParameter(0) / genEv;
   double dcorr = 1e6*bwg->GetParError(0) / genEv;
   cout << " genEv= " << genEv << " corr= " << corr << endl;

   double intg108 = bwg->Integral(dL,1.08,1e-8) / bin_width;
   double norm108 = bwg->GetParameter(0) / intg108;
   cout << " norm108= " << norm108 << endl;

   TLegend* leg = new TLegend(0.24,0.17,0.60,0.39);
   leg->SetTextSize(0.032);
   leg->SetHeader(
         Form("Cut-off correction: %.3f #pm %.3f",corr,dcorr),"C" );
   if ( date > 0 ) {
      leg->AddEntry(hst, Form("MC truth %i",date),"PE");
   } else {
      leg->AddEntry(hst, "MC truth", "PE");
   }
   leg->AddEntry(bwg,"Breit-Wigner","L");
   leg->AddEntry(lB,"Cut-off","L");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{2 Fit Mkk(MC-reconstructed) by BW x Gauss
// {{{3 * 'myFCN_fitSig' class
//--------------------------------------------------------------------
struct myFCN_fitSig {
   const double sl = -1.97; // v709n4
   vector<double> mkk;   // data in central part

   //-----------------------------------------------------------------
   myFCN_fitSig(int date)
   //-----------------------------------------------------------------
   {
      // Get non-binned data in central part
      if ( date > 0 ) {
         string fname( Form("mcsig_kkmc_%02i.root",date%100) );
         mkk = get_mkk_vec(fname,0);
      } else {
         for ( int date : {2021,2012,2009} ) {
            auto file = Form("mcsig_kkmc_%02i.root",date%100);
            auto&& mk= get_mkk_vec(file,0);
            mkk.insert(end(mkk),begin(mk),end(mk));
         }
      }
      ReoderMkk(this->mkk); // the same result
   }

   // minimization function
   // the signature of this operator() MUST be exactly this:
   // * Kahan-Babushka-Neumaier summation
   //   https://en.wikipedia.org/wiki/Kahan_summation_algorithm
   //-----------------------------------------------------------------
   double operator() (const double* p)
   //-----------------------------------------------------------------
   {
      const double& mphi = p[0];
      const double& gphi = p[1];
      const double& sigma = p[2];

      double res = 0;
      double cres = 0; // compensation for lost low-order bits of res
      for ( const auto& m : mkk ) {
         double L = BreitWignerGaussN(m,mphi,gphi,sigma,sl);

         // -> two sigma Mkk resolution
         // auto [bw1,bw2] = BreitWignerGaussN(m,mphi,gphi,p[2],p[3],sl);
         // double L = p[4] * bw1 + (1-p[4]) * bw2;

         double lh = (L>0) ? -2*log(L) : FLT_MAX;
         double t = res + lh;
         if ( fabs(res) >= fabs(lh) ) {
            cres += (res - t) + lh;
         } else {
            cres += (lh - t) + res;
         }
         res = t;
      }
      res += cres; // Correction only applied once in the very end

      return res;
   }
};

// {{{3 * sig_fit()
//--------------------------------------------------------------------
void sig_fit(int date, string pdf, size_t Cx, size_t Cy)
//--------------------------------------------------------------------
{
   myFCN_fitSig my_fcn(date);  // class for 'FitFCN'
   const auto& sl = my_fcn.sl; // efficiency parameters
   auto& mkk = my_fcn.mkk;

   TH1D* hst = get_hst(mkk,Form("mkksig_%d",date));

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Sig(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Sig.Add(mkk[i]);
   }

   double norm = hst->Integral();
   cout << " norm= " << norm << " n= " << n << endl;

   //-----------------------------------------------------------------
   // == fit configuration
   ROOT::Fit::Fitter fitter;
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   // Possible printing levels are from 0 (minimal) to 3 (maximum)
   fitter.Config().MinimizerOptions().SetPrintLevel(2);

   vector<string> par_name { "Mphi",  "Gphi", "sigma" };
   vector<double> par_ini  {   Mphi, 0.00426, 1.2e-3  };

   // -> two sigma Mkk resolution
   // vector<string> par_name { "Mphi", "Gphi", "sig1", "sig2",  "f" };
   // vector<double> par_ini  {  Mphi,  0.0044, 0.8e-3, 1.7e-3, 0.75 };

   const unsigned int Npar = par_name.size(); // number of parameters

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(Mphi-1e-3,Mphi+1e-3);
   // fitter.Config().ParSettings(0).Fix();
   fitter.Config().ParSettings(1).SetLimits(Gphi-1e-3,Gphi+1e-3);
   //  fitter.Config().ParSettings(1).Fix();
   fitter.Config().ParSettings(2).SetLimits(0.2e-3, 2e-3);
   // -> two sigma Mkk resolution
   // fitter.Config().ParSettings(3).SetLimits(1.e-3, 10e-3);
   // fitter.Config().ParSettings(4).SetLimits(0., 1.);

   // == Fit
   int Ndat = mkk.size();
   fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // false=likelihood
   fitter.CalculateHessErrors(); // without MINOS here
   // fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   // double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   //  vector<double> Fpar { 1.019524029, 0.00426, 0.001207400594 };
   printf("vector<double> Fpar { ");
   for ( const auto& fp : Fpar ) {
      printf("%.10g, ",fp);
   }
   printf("}; // to debug\n");

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   auto& mkkt = mkk;
   auto& hstt = hst;
   // Test for limited statistic
   //  size_t Nmkk = 30'000;
   //  vector<double> mkkt(begin(mkk),begin(mkk)+Nmkk);
   //  TH1D* hstt = get_hst(mkkt,Form("mkkt"));

   auto Fun = [Fpar,sl](double m) {
      return BreitWignerGaussN(m,Fpar[0],Fpar[1],Fpar[2],sl);
      // -> two sigma Mkk resolution
      // auto [bw1,bw2] =
         // BreitWignerGaussN(m,Fpar[0],Fpar[1],Fpar[2],Fpar[3],sl);
      // return Fpar[4] * bw1 + (1-Fpar[4]) * bw2;
   };
   double pvalueKS = FastKSTest(Fun,mkkt);
   printf("pvalueKS= %.3g\n", pvalueKS);

   // Get chi-square of histogram
   auto Fch = new TF1("Fch",
         [&Fun](double* x,double* p){return Fun(x[0]);}, dL,dU,0);
   auto [ch2_hst,ndof_hst] = GetCh2Ndof(hstt,Fch,hstt->Integral());
   ndof_hst -= res.NFreeParameters();
   printf("chi^2/Ndof= %.1f / %d\n",ch2_hst,ndof_hst);

   //-----------------------------------------------------------------
   // Functions to draw
   auto Lfit = [norm,Fun](double* x,double* p) -> double {
      return  bW * norm * Fun(x[0]);
   };
   TF1* ffit = new TF1("ffit", Lfit, dL, dU, 0);
   ffit->SetLineWidth(2);
   ffit->SetLineColor(kRed);
   ffit->SetNpx(500);

   //-----------------------------------------------------------------
   auto name = Form("c1_mkksig_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   SetHstFace(hst);
   if ( date == 0 ) {
      hst->SetMinimum(4.);
      hst->SetMaximum(2.e5);
   } else {
      hst->SetMinimum(1.);
   }
   hst->GetYaxis()->SetTitleOffset(1.1);
   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
   ffit->Draw("SAME");

   TLegend* leg = new TLegend(0.57,0.815,0.892,0.89);
   if ( date > 0 ) {
      leg->AddEntry(hst,Form("MC signal #phi#eta %i",date),"EP");
   } else {
      leg->AddEntry(hst,"MC signal #phi#eta","EP");
   }
   leg->AddEntry(ffit,"(#it{BW}#times#it{E}) #otimes #it{Gauss}","L");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.57,0.61,0.892,0.81,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#it{p-value(KS) = %.3f}",pvalueKS) );
   pt->AddText( Form("#chi^{2} / ndof = %.1f / %d",
            ch2_hst,ndof_hst) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".3f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt->AddText( Form("#sigma = %s MeV", PE.Eform(2,".2f",1e3)) );
   // -> two sigma Mkk resolution
   // pt->AddText( Form("#sigma_{1} = %s MeV",PE.Eform(2,".2f",1e3)) );
   // pt->AddText( Form("#sigma_{2} = %s MeV",PE.Eform(3,".2f",1e3)) );
   // pt->AddText( Form("f = %s %%", PE.Eform(4,".1f",1e2)) );
   pt->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 2.Argus functions
//--------------------------------------------------------------------
// {{{2 Functions
//--------------------------------------------------------------------
double Argus(double x, double a)
//--------------------------------------------------------------------
{
   // ARGUS distribution:
   // https://en.wikipedia.org/wiki/ARGUS_distribution
   // This is a regular ARGUS distribution with one parameter 'a'

   if ( x < 0 || x > 1 ) {
      return 0;           // by definition
   }

   double a2 = a*a;
   double y2 = 1-x*x;
   double tmp = x * sqrt(y2) * exp(-0.5*a2*y2);

   // normalization to one
   // 1) small 'a'
   constexpr double a2_eps = 1e-4;
   if ( a2 < a2_eps ) {
      return tmp * 3/(1-0.3*a2); // a2_eps = 1e-4
      // constexpr double c4 = 3./56.;
      // return tmp * 3/(1+a2*(-0.3+a2*c4)); // a2_eps = 1e-3
   }

   // 2) general formula
   constexpr double one_over_sqrt2pi = 0.5*M_2_SQRTPI*M_SQRT1_2;
   double Psi = 0.5*(1+erf(a*M_SQRT1_2))
      - a*one_over_sqrt2pi*exp(-0.5*a2)
      -0.5;
   double norm = a2*a*one_over_sqrt2pi/Psi;
   return tmp*norm;
}

//--------------------------------------------------------------------
double RevArgus(double m, double A)
//--------------------------------------------------------------------
{
   // Argus function with left cutoff (L)
   //      U is the upper limit of fitting
   //      L,U - must be const for fit

   static const double L = 2*Mk;
   static const double U = Mjpsi - Meta;

   double x = (U-m)/(U-L); // transformation: (L,U) -> (1,0)
   return Argus(x,A)/(U-L);
}

//--------------------------------------------------------------------
double RevArgusN(double m, double A, double slope)
//--------------------------------------------------------------------
{
   // Numeric normalisation on one for range [dL,dU]
   // * multiply Argus by function of efficiency(m) (parameter slope)

   double norm = 0;
   // cash parameters
   static double cacheN = 0;
   static double cacheA = 0;
   static double cacheS = 0;
   if ( cacheN > 0 && A == cacheA && slope == cacheS ) {
      norm = cacheN;
   } else {
      // integrand lambda function
      double p[] = {A,slope};
      auto Lint = [](double x, void* pp) -> double{
         double* p = static_cast<double*>(pp);
         return RevArgus(x,p[0]) * (1+p[1]*(x-1.02));
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-12, 1.e-9, 1000);
      norm = gsl_int.Integral(Lint,p,dL,dU);

      if ( gsl_int.Status() ) { // true for error
         double abs_err = gsl_int.Error();
         printf("%s error for A=%.3f: norm=%g, error=%g\n",
               __func__,A,norm,abs_err);
      }

      cacheN = norm;
      cacheA = A;
      cacheS = slope;
   }

   return RevArgus(m,A) * (1+slope*(m-1.02)) / norm;
}

// {{{2 Tests
//--------------------------------------------------------------------
void test_RevArgus(size_t Cx, size_t Cy)
//--------------------------------------------------------------------
{
   auto Largus = [](const double* x,const double* p) -> double {
      // return p[0]*Argus(x[0],p[1]);
      return p[0]*RevArgus(x[0],p[1]);
   };

   // TF1* fbkg = new TF1("fbkg", Largus, 0., 1., 2);
   TF1* fbkg = new TF1("fbkg", Largus, 0.98, 1.08, 2);
   fbkg->SetParNames("N","a");
   fbkg->SetParameters(1., 5.);
   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);

   TCanvas* c1 = new TCanvas("c1_RevArg","...",0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   TLegend* leg = new TLegend(0.54,0.69,0.892,0.89);
   leg->SetHeader("RevArgus(mkk,a)","C");

   fbkg->DrawCopy();
   leg->AddEntry(fbkg->Clone(),
         Form("a=%.1f",fbkg->GetParameter(1)), "L");

   // test normalization:
   printf(" norm= %.15f\n",fbkg->Integral(2*Mk,Mjpsi-Meta));

   fbkg->SetParameters(10., 0.);
   fbkg->SetLineColor(kRed);
   fbkg->DrawCopy("SAME");
   leg->AddEntry(fbkg->Clone(),
         Form("a=%.1f",fbkg->GetParameter(1)), "L");
   leg->Draw();

   // test normalization:
   printf(" norm= %.15f\n",fbkg->Integral(2*Mk,Mjpsi-Meta));

   gPad->RedrawAxis();
   c1->Update();
}

//--------------------------------------------------------------------
void test_RevArgusN(size_t Cx, size_t Cy) {
//--------------------------------------------------------------------
   auto Lrargus = [](const double* x,const double* p) -> double {
      return p[0]*RevArgusN(x[0],p[1],p[2]);
   };

   TF1* fbkg = new TF1("fbkg", Lrargus, 0.98, dU, 3);
   fbkg->SetParNames("N","A","Sl");
   fbkg->SetParameters(1., 0., -2.0);
   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);

   TCanvas* c1 = new TCanvas("c1_RevArgN","...",0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   TLegend* leg = new TLegend(0.35,0.20,0.75,0.40);
   leg->SetHeader("RevArgusN(mkk,a)*Eff(sl)","C");

   fbkg->DrawCopy();
   leg->AddEntry(fbkg->Clone(), Form("a=%.2f sl=%.2f",
            fbkg->GetParameter(1),
            fbkg->GetParameter(2)), "L");

   // test normalization:
   printf(" RevArgusN norm= %.15f\n",fbkg->Integral( dL,dU ) );

   fbkg->SetParameters(1., 5., -2.0);
   fbkg->SetLineColor(kRed);
   fbkg->SetLineStyle(kDashed);
   fbkg->DrawCopy("SAME");
   leg->AddEntry(fbkg->Clone(), Form("a=%.2f sl=%.2f",
            fbkg->GetParameter(1),
            fbkg->GetParameter(2)), "L");

   printf(" RevArgusN norm2= %.15f\n",fbkg->Integral( dL,dU ) );

   fbkg->SetParameters(1., 0., 0);
   fbkg->SetLineColor(kGreen);
   fbkg->DrawCopy("SAME");
   leg->AddEntry(fbkg->Clone(), Form("a=%.2f sl=%.2f",
            fbkg->GetParameter(1),
            fbkg->GetParameter(2)), "L");

   leg->Draw();
}

//--------------------------------------------------------------------
void test_Argus_smallA(size_t Cx, size_t Cy)
//--------------------------------------------------------------------
{
   auto Largus = [](const double* x,const double* p) -> double {
      return p[0]*Argus(x[0],p[1]);
   };

   TF1* far = new TF1("far", Largus, 0., 1., 2);
   far->SetParNames("N","a");
   far->SetParameters(1., 0.5);
   far->SetLineWidth(2);
   far->SetLineColor(kBlue);

   TCanvas* c1 = new TCanvas("c1_Arg","...",0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   TLegend* leg = new TLegend(0.11,0.75,0.41,0.89);
   leg->SetHeader("Argus(x,a)","C");

   far->DrawCopy();
   leg->AddEntry(far->Clone(),
         Form("a=%.2g",far->GetParameter(1)), "L");

   // test normalization:
   printf(" norm= %.15f\n",far->Integral(0.,1.));

   far->SetParameters(1., 0.009);
   far->SetLineColor(kRed);
   far->DrawCopy("SAME");
   leg->AddEntry(far->Clone(),
         Form("a=%.2g",far->GetParameter(1)), "L");
   leg->Draw();

   // test normalization:
   printf(" norm= %.15f\n",far->Integral(0.,1.));

   gPad->RedrawAxis();
   c1->Update();
}

// {{{2 Fit data mkk side-band by Argus + mcsig(sb)
//--------------------------------------------------------------------
void sideband_fit(string pdf, size_t Cx, size_t Cy)
//--------------------------------------------------------------------
{
   // Get un-binned side-band data
   vector<double> mkk;
   for ( int d : {2021,2012,2009} ) {
      auto&& mk= get_mkk_vec(Form("data_%02ipsip_all.root",d%100),1);
      mkk.insert(end(mkk),begin(mk),end(mk));
   }
   TH1D* hst = get_hst(mkk,"data_sb");

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Sb(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Sb.Add(mkk[i]);
   }

   double norm = hst->Integral();
   cout << "norm= " << norm << " n= " << n << endl;

   // 2. Get side-band MC-signal as histogram
   vector<double> mcs;
   for ( int d : {2021,2012,2009} ) {
      auto&& mk= get_mkk_vec(Form("mcsig_kkmc_%02i.root",d%100),1);
      mcs.insert(end(mcs),begin(mk),end(mk));
   }
   TH1D* hmcS = get_hst(mcs,"mcs_sb");
   double nmc = hmcS->Integral();
   hmcS->Scale(1./(nmc*bW));
   cout << "MC-sig: norm= " << nmc << " n= " << mcs.size() << endl;

   //-----------------------------------------------------------------
   // Fit side-band by reversed Argus function
   // using extended likelihood here
   double slope = -1.97; // v709n4
   auto Lan = [slope,hmcS](const double* x,const double* p) {
      int bin = hmcS->GetXaxis()->FindBin(x[0]);
      return p[1]*RevArgusN(x[0],p[0],slope)
         + p[2]*hmcS->GetBinContent(bin);

   };
   vector<string> par_name { "a", "Nsb", "Ns" };
   // vector<double> par_ini  {  5., 550, norm-500 };
   vector<double> par_ini  {  5., norm*0.8, norm*0.2 };

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* sb = new TF1("sb", Lan, dL, dU, Npar);

   // == fit configuration
   ROOT::Fit::Fitter fitter;
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   // Possible printing levels are from 0 (minimal) to 3 (maximum)
   fitter.Config().MinimizerOptions().SetPrintLevel(1);

   ROOT::Math::WrappedTF1 fsb( *sb );
   fitter.SetFunction(fsb,false); // false=no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(0., 10.); // a
   // fitter.Config().ParSettings(0).Fix();

   fitter.LikelihoodFit(Sb,true); // true=extended likelihood fit
   fitter.CalculateHessErrors();
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   auto Fun = [Lan,Fpar](double x) -> double {
      return Lan(&x,Fpar.data());
   };
   double pvalueKS = FastKSTest(Fun,mkk);
   printf("pvalueKS= %.3g\n", pvalueKS);

   //-----------------------------------------------------------------
   // Functions to draw
   auto Lfit = [Lan,Fpar](double* x,double* p) {
      return bW * Lan(x,Fpar.data());
   };
   TF1* ffit = new TF1("ffit", Lfit, dL, dU, 0);
   ffit->SetLineWidth(2);
   ffit->SetLineColor(kRed);
   ffit->SetNpx(200);

   auto Lsig = [Fpar,hmcS](double* x,double* p) {
      int bin = hmcS->GetXaxis()->FindBin(x[0]);
      return  bW*Fpar[2]*hmcS->GetBinContent(bin);
   };
   TF1* fsig = new TF1("fsig", Lsig, dL, dU, 0);
   fsig->SetLineWidth(2);
   fsig->SetLineColor(kGreen+3);
   fsig->SetNpx(200);

   auto Lar = [Fpar,slope](double* x,double* p) {
      return bW * Fpar[1]*RevArgusN(x[0],Fpar[0],slope);
   };
   TF1* far = new TF1("far", Lar, dL, dU, 0);
   far->SetLineWidth(2);
   far->SetLineColor(kBlue);
   far->SetNpx(200);

   //-----------------------------------------------------------------
   auto name = Form("c1_sbfit");
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst);
   double hst_max = hst->GetMaximum();
   hst->SetMaximum( 5*floor(1.4*hst_max/5) );
   hst->GetYaxis()->SetMaxDigits(3);
   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetTitleOffset(1.1);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
   ffit->Draw("SAME");
   fsig->Draw("SAME");
   far->Draw("SAME");

   TLegend* leg = new TLegend(0.11,0.70,0.35,0.89);
   leg->AddEntry(hst,"Side-band data","LEP");
   leg->AddEntry(ffit,"Fitting result","L");
   leg->AddEntry(far,"Argus function","L");
   leg->AddEntry(fsig,"MC signal #phi#eta","L");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.57,0.70,0.892,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#it{p-value(KS) = %.3f}",pvalueKS) );
   pt->AddText( Form("a_{bg} = %s",PE.Eform(0,".2f")) );
   pt->AddText( Form("N_{bg} = %s",PE.Eform(1,".0f")) );
   pt->AddText( Form("N_{sl} = %s",PE.Eform(2,".0f")) );
   pt->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 3.No Interference: BW(x)Gauss + Argus
//--------------------------------------------------------------------
void data_fit(size_t Cx, size_t Cy)
//--------------------------------------------------------------------
{
   // Get un-binned data in central part of Mgg
   vector<double> mkk;
   for ( int d : {2021,2012,2009} ) {
      auto&& mk= get_mkk_vec(Form("data_%02ipsip_all.root",d%100),0);
      mkk.insert(end(mkk),begin(mk),end(mk));
   }
   TH1D* hst = get_hst(mkk,"data_cp");
   ReoderMkk(mkk);

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Dat(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Dat.Add(mkk[i]);
   }

   double norm = hst->Integral();
   cout << " norm= " << norm << " n= " << n << endl;

   //-----------------------------------------------------------------
   // Fit data
   double slope = -1.97; // v709n4
   // use extended LH
   auto Lsum = [slope](const double* x,const double* p) -> double {
      return p[4] * BreitWignerGaussN(x[0],p[0],p[1],p[2],slope) +
         p[5] * RevArgusN(x[0],fabs(p[3]),slope);
   };
   vector<string> par_name {"Mphi","Gphi","sigma","a","Nphi","Nbg"};
   vector<double> par_ini  { Mphi, Gphi, 1e-3, 5.2, norm-554, 554.};

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* fsum = new TF1("fsum", Lsum, dL, dU, Npar);

   // == fit configuration
   ROOT::Fit::Fitter fitter;
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   // Possible printing levels are from 0 (minimal) to 3 (maximum)
   fitter.Config().MinimizerOptions().SetPrintLevel(2);

   ROOT::Math::WrappedTF1 WFsum( *fsum );
   fitter.SetFunction(WFsum,false); // false=no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(Mphi-1e-3,Mphi+1e-3);
   // fitter.Config().ParSettings(0).Fix(); // Mphi
   fitter.Config().ParSettings(1).SetLimits(Gphi-1e-3,Gphi+1e-3);
   // fitter.Config().ParSettings(1).Fix(); // Gphi
   fitter.Config().ParSettings(2).SetLimits(0.2e-3, 2.e-3);// sigma
   fitter.Config().ParSettings(3).SetLimits(0.,10.);       // Argus
   // fitter.Config().ParSettings(3).Fix(); // Argus
   // fitter.Config().ParSettings(4).SetLimits(0.,1.5*norm);  // Nphi
   // fitter.Config().ParSettings(5).SetLimits(0.,norm);      // Nbg

   fitter.LikelihoodFit(Dat,true); // true=extended likelihood fit
   fitter.CalculateHessErrors();
   // fitter.CalculateMinosErrors(); // does not work for many par-s

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   auto Fun = [Lsum,Fpar](double x) -> double {
      return Lsum(&x,Fpar.data());
   };
   double pvalueKS = FastKSTest(Fun,mkk);
   printf("pvalueKS= %.3g\n", pvalueKS);

   // Get chisquare of histogram
   auto [ch2_hst,ndof_hst] = GetCh2Ndof(hst,fsum);
   ndof_hst -= res.NFreeParameters();
   printf("hst: chi^2/Ndof= %.1f / %d\n",ch2_hst,ndof_hst);

   //-----------------------------------------------------------------
   // Functions to draw
   auto Lfit = [Lsum,Fpar](double* x,double* p) -> double {
      return bW * Lsum(x,Fpar.data());
   };
   TF1* ffit = new TF1("ffit", Lfit, dL, dU, 0);
   ffit->SetLineWidth(2);
   ffit->SetLineColor(kRed);
   ffit->SetNpx(500);

   auto Lbw = [Fpar,slope](double* x,double* p) -> double {
      return bW * Fpar[4] *
         BreitWignerGaussN(x[0],Fpar[0],Fpar[1],Fpar[2],slope);
   };
   TF1* fbw = new TF1("fbw", Lbw, dL, dU, 0);
   fbw->SetLineWidth(2);
   fbw->SetLineStyle(kDashed);
   fbw->SetLineColor(kGreen+2);

   auto Lar = [Fpar,slope](double* x,double* p) -> double {
      return bW * Fpar[5] * RevArgusN(x[0],fabs(Fpar[3]),slope);
   };
   TF1* far = new TF1("far", Lar, dL, dU, 0);
   far->SetLineWidth(2);
   far->SetLineStyle(kDashed);
   far->SetLineColor(kCyan+2);

   //-----------------------------------------------------------------
   auto name = Form("c1_data_fit");
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst);
   hst->SetMaximum(3e3);
   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetMaxDigits(3);
   hst->GetYaxis()->SetTitleOffset(1.1);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
   fbw->Draw("SAME");
   far->Draw("SAME");
   ffit->Draw("SAME");
   hst->Draw("EP SAME");

   TLegend* leg = new TLegend(0.12,0.72,0.34,0.89);
   leg->AddEntry(hst, "Data","LEP");
   leg->AddEntry(ffit,"Fitting result","L");
   leg->AddEntry(fbw, "Breit-Wigner","L");
   leg->AddEntry(far, "Background","L");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.57,0.54,0.892,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#it{p-value(KS) = %.3f}",pvalueKS) );
   pt->AddText( Form("#chi^{2} / ndof = %.0f / %d",
            ch2_hst,ndof_hst) );
   pt->AddText( Form("N_{#phi}= %s",PE.Eform(4,".0f")) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt->AddText( Form("#sigma = %s MeV", PE.Eform(2,".2f",1e3)) );
   pt->AddText( Form("N_{bg}= %s",PE.Eform(5,".0f")) );
   pt->AddText( Form("a_{bg} = %s",PE.Eform(3,".1f")) );
   pt->Draw();

   gPad->RedrawAxis();
   c1->Update();
   const char* pdf = "mkk_fit.pdf";
   c1->Print( pdf );
}

// {{{1 4.Interference BW with Argus
//--------------------------------------------------------------------
// {{{2 Functions
//--------------------------------------------------------------------
vector<double> IntfrBWAR( double m,
      double mphi, double gphi,      // B-W
      double A, double F, double Ang // Argus & interference
      )
//--------------------------------------------------------------------
// The return vector contains:
//                   0    1     2        3
//        ret     { Sum, B-W, Argus, Interference }
//--------------------------------------------------------------------
{
   if ( m < dL ) { // phase space becomes zero
      vector<double> ret(4,0.);
      return ret;
   }

   constexpr double Mjpsi2 = SQ(Mjpsi);
   constexpr double Meta2  = SQ(Meta);
   constexpr double Mk2    = SQ(Mk);

   double m2    = SQ(m);
   double mphi2 = SQ(mphi);
   double gphi2 = SQ(gphi);

   double R     = R_BW; // Blatt-Weisskopf ff

   double gam = mphi*sqrt(mphi2+gphi2);
   double kap = 2*M_SQRT2*mphi*gphi*gam / (M_PI*sqrt(mphi2+gam));

   double p_phi0 = sqrt(SQ(Mjpsi2-mphi2-Meta2)-4*mphi2*Meta2) /
      (2*Mjpsi);
   double p_phi  = sqrt(SQ(Mjpsi2-m2-Meta2)-4*m2*Meta2) / (2*Mjpsi);
   double r_phi  = p_phi / p_phi0;

   double p_k0 = 0.5*sqrt(mphi2-4*Mk2);
   double p_k  = 0.5*sqrt(m2-4*Mk2);
   double r_k  = p_k / p_k0;

   // B(m)/B(m0)
   double BB_phi = sqrt( (1+SQ(R*p_phi0)) / (1+SQ(R*p_phi)) );
   double BB_k2  = (1+SQ(R*p_k0)) / (1+SQ(R*p_k));
   double BB_k   = sqrt( BB_k2 );

   double fm = r_phi*r_k*BB_phi*BB_k;
   double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k2; // == m*G(m)

   // common multipliers
   double tmp =  fm / ( SQ(m2-mphi2) + SQ(GM) );

   // signal
   double BW = tmp*fm*kap;

   // background:
   double Ar_fun = RevArgus(m,A);
   double Ar = Ar_fun*SQ(F);

   // interference:
   double Intfr = 2 * tmp*F*sqrt(kap*Ar_fun) *
                  ( (m2-mphi2)*cos(Ang) - GM*sin(Ang) );

   double Sum = BW + Ar + Intfr;
   vector<double> ret { Sum, BW, Ar, Intfr };
   return ret;
}

//--------------------------------------------------------------------
double IntfrBWARG( double m,
      const double p[],    // {mphi,gphi,sigma,A,F,Ang,sl}
      unsigned int idx = 0 // what to return
      )
//--------------------------------------------------------------------
// Function Fun(x) folded with a normal distribution:
//
//        1                       [                  (X'-X)^2    ]
// ---------------- *  Integral   [ Fun(X') * exp( - --------- ) ] dX'
// sqrt(2*pi)*sigma   (-oo<X'<+oo)[                  2*sigma^2   ]
//
//--------------------------------------------------------------------
// integration in +/-5 sigma interval
//--------------------------------------------------------------------
{
   static const double Ngauss = 1./sqrt(2*M_PI);

   // integrand lambda function
   auto Lint = [m,p,idx](double t) -> double{
      // t == (X'-X)/sigma
      double x = m+t*p[2];
      vector<double> Fun = IntfrBWAR(x, p[0],p[1],p[3],p[4],p[5]);
      // efficiency correction: slope
      double effcor = 1 + p[6]*(x-1.02);
      return exp(-0.5*t*t) * Fun[idx] * effcor;
   };

   rmath_fun< decltype(Lint) > Fint(Lint);

   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-12, 1.e-9, 1000);
   double result = gsl_int.Integral(Fint,-5.,+5.);

   if ( gsl_int.Status() ) { // true for error
      double abs_err = gsl_int.Error();
      printf("%s error for m=%.3f: result=%g, error=%g\n",
            __func__,m,result,abs_err);
      printf("           mphi= %g, gphi=%g, sig=%g\n",p[0],p[1],p[2]);
      printf("           a= %g, F=%g, Ang=%g, sl= %g\n",
            p[3],p[4],p[5],p[6]);
      // exit(0);
   }

   return Ngauss * result;
}

// Numerical normalisation to one for range [dL,dU]
//--------------------------------------------------------------------
double IntfrBWARGN( double m,
      const double p[], // {mphi,gphi,sigma,A,F,Ang,sl}
      int idx = 0       // what to return
      )
//--------------------------------------------------------------------
{
   // cash parameters
   static double cacheN = 0;
   constexpr int npar = 7;
   static double cP[npar]; // cache p[]

   double norm = cacheN;
   bool need_to_calc = !(cacheN > 0);
   if ( !need_to_calc ) {
      for (int i = 0; i < npar; ++i) {
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
         return IntfrBWARG(x,p,0);  // MUST be idx=0
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-12, 1.e-9, 1000);
      norm = gsl_int.Integral(Lint,(void *)p,dL,dU);

      if ( gsl_int.Status() ) { // true for error
         double abs_err = gsl_int.Error();
         printf("%s error: norm=%g, error=%g\n",
               __func__,norm,abs_err);
      }

      // save to cache
      cacheN = norm;
      for (int i = 0; i < npar; ++i) {
         cP[i] = p[i];
      }
   }

   // 'idx == 100' means to return normalization only
   if ( idx == 100 ) {
      return norm;
   }

   return IntfrBWARG(m,p,idx) / norm;
}

// {{{2 Tests
//--------------------------------------------------------------------
void test_Intfr(size_t Cx, size_t Cy)
//--------------------------------------------------------------------
{
   double mphi = Mphi, gphi = Gphi;
   auto Lintfr = [mphi,gphi](const double* x, const double* p) {
      double pp[] { mphi, gphi, p[0],p[1],p[2],p[3],p[4] };
      return IntfrBWARGN( x[0], pp, int(p[5]) );
   };

   int npar = 6;
   TF1* fun = new TF1("fun", Lintfr, dL, dU, npar);
   fun->SetParNames("Sigma", "A", "F", "Ang",  "Sl", "T");
   fun->SetParameters(1.2e-3, 0., 0.8,   0.7,  -1.98,  0.);
   fun->SetNpx(500);

   TCanvas* c1 = new TCanvas("c1_Intfr","...",0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   TH1D* tmp = new TH1D("tmp","",500,bL,dU);
   tmp->SetAxisRange(-20.,130.,"Y");
   tmp->Draw();

   TLegend* leg = new TLegend(0.54,0.69,0.892,0.89);
   leg->SetHeader("Breit-Wigner #odot Argus","C");

   fun->SetParameter(npar-1, 0); // draw all
   fun->SetLineWidth(2);
   fun->SetLineStyle(kSolid);
   fun->SetLineColor(kRed);
   fun->DrawCopy("SAME");
   printf(" Norm= %.6f\n", fun->Integral(dL,dU,1e-4) );
   leg->AddEntry(fun->Clone(), "Sum", "L");

   // draw components
   fun->SetLineWidth(1);
   fun->SetLineStyle(kDashed);

   fun->SetParameter(npar-1, 1); // BW
   fun->SetLineColor(kGreen+2);
   fun->DrawCopy("SAME");
   leg->AddEntry(fun->Clone(), "Breit-Wigner", "L");

   fun->SetParameter(npar-1, 2); // Argus
   fun->SetLineColor(kBlue);
   fun->DrawCopy("SAME");
   leg->AddEntry(fun->Clone(), "Argus", "L");

   fun->SetParameter(npar-1, 3); // interference
   fun->SetLineColor(kMagenta+1);
   fun->DrawCopy("SAME");
   leg->AddEntry(fun->Clone(), "Interference", "L");

   leg->Draw();
   gPad->RedrawAxis();
   c1->Update();
}

// {{{2 -- data_Intfr_fit(): BW(x)Gauss (*) Argus
// {{{3 * Calculation of integrals and their errors for data_Intfr_fit
//--------------------------------------------------------------------
ValErr CalcIntEr( const ROOT::Fit::FitResult& res,
      double sl, unsigned int idx )
//--------------------------------------------------------------------
{
   // desired errors:
   constexpr double eps_abs = 1.e-6;
   constexpr double eps_rel = 1.e-6;
   constexpr double epsilon = 1.e-3; // max numerical error

   // get parameters and covariation matrix of error
   int Npar = res.NPar();
   vector<double> Fpar = res.Parameters();

   if ( Npar != 6 ) {
      cerr << " FATAL ERROR in " << __func__
         << " Npar= " << Npar << endl;
      exit(EXIT_FAILURE);
   }

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
   auto Lfc = [sl,idx](const double* x,const double* p) {
      double m = x[0];
      double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl }; // real slope
      double normI = IntfrBWARGN( m,pp,100 );
      pp[6] = 0.; // set slope to zero
      double intfr = IntfrBWARG( m,pp,idx ) / normI;
      return intfr;
   };
   TF1 FC("FC", Lfc, dL, dU, Npar);
   FC.SetParameters( Fpar.data() );

   double err_num = 0;
   double Integ = FC.IntegralOneDim(dL,dU,eps_rel,eps_abs, err_num);
   // cout << " err_num(" << idx << ") = " << err_num << endl;
   if ( err_num > epsilon ) {
      cout << " WARNING: " << __func__ << " numerical error of"
         " integration of F_C(" << idx << ") is too big "
         << err_num << endl;
   }

   // 2) calculate error of this integral
   TMatrixDSym covMatrix(Npar);
   covMatrix.Use(Npar,cov_m.data());
   // printf("\ncovMatrix : ");
   // covMatrix.Print();

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
      // cout << " abs err_num(" << i << ") = " << err_num << endl;
      err_num2 += SQ(err_num);
      // cout << " IntGrad[" << i << "] = " << IntGrad[i] << endl;
   }

   double err_Int = sqrt( covMatrix.Similarity(IntGrad) );
   err_num2 = sqrt( err_num2 / err_Int ); // abs numerical error
   // cout << " err_num(" << idx << ") = " << err_num2 << endl;
   if ( err_num2 > epsilon ) {
      cout << " WARNING: " << __func__ << " numerical error "
         "of integration of dF is too big " << err_num2 << endl;
   }

   return (ValErr {Integ,err_Int});
}

// {{{3 * Fit data mkk by BW(x)Gauss (*) Argus: data_Intfr_fit
//--------------------------------------------------------------------
void data_Intfr_fit(size_t Cx, size_t Cy)
//--------------------------------------------------------------------
{
   // Get un-binned data in central part of Mgg
   vector<double> mkk;
   for ( int d : {2009,2012,2021} ) {
      auto&& mk= get_mkk_vec(Form("data_%02ipsip_all.root",d%100),0);
      mkk.insert(end(mkk),begin(mk),end(mk));
   }
   TH1D* hst = get_hst(mkk,"data_cp");

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Dat(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Dat.Add(mkk[i]);
   }

   double norm = hst->Integral();
   cout << " norm= " << norm << " n= " << n << endl;

   //-----------------------------------------------------------------
   // Fit data
   double sl = -1.97; // v709n4
   // function MUST be normalized to 1 on the fit range
   auto Lfit = [sl](const double* x,const double* p) -> double {
      const double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl };
      return IntfrBWARGN( x[0], pp, 0 );
   };
   vector<string> par_name { "Mphi","Gphi","sig","a","F","angle" };

   vector<double> par_ini;
   par_ini = {Mphi, Gphi, 1.2e-3, 3., 0.6, 1.}; //n
   // par_ini = {Mphi, Gphi, 1.2e-3, 3., 0.6, -1.1}; //p
   bool neg_pos = par_ini.back() > 0; //true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   // == fit configuration
   ROOT::Fit::Fitter fitter;
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   // Possible printing levels are from 0 (minimal) to 3 (maximum)
   fitter.Config().MinimizerOptions().SetPrintLevel(2);

   ROOT::Math::WrappedTF1 WFit( *Ffit );
   fitter.SetFunction(WFit,false); // false=no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(Mphi-1e-3, Mphi+1e3);
   // fitter.Config().ParSettings(0).Fix(); // Mphi
   fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3, Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix(); // Gphi
   fitter.Config().ParSettings(2).SetLimits(0.5e-3, 2.e-3); // sig
   fitter.Config().ParSettings(3).SetLimits(0.,10.); // Argus
   // fitter.Config().ParSettings(3).SetValue(5.2);
   // fitter.Config().ParSettings(3).Fix(); // Argus
   fitter.Config().ParSettings(4).SetLimits(0.01, 10.); // F
   fitter.Config().ParSettings(5).SetLimits(-M_PI, M_PI); // angle

   fitter.LikelihoodFit(Dat,false); // true=extended likelihood fit
   // fitter.CalculateMinosErrors();
   fitter.CalculateHessErrors(); // for correct integral errors

   const ROOT::Fit::FitResult& res = fitter.Result();
   res.Print(cout);
   // res.PrintCovMatrix(cout);

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res,0.05);
   vector<double>& Fpar = PE.Fpar;

   vector<string> names { "NKK", "Nphi", "Nnonphi", "Nifr"  };
   vector<ValErr> Nos;
   Nos.reserve(4);
   cout << "n= " << n << endl;
   for ( int idx = 0; idx <= 3; ++idx ) {
      ValErr Integ = CalcIntEr( res, sl, idx );
      // printf("Int[%i] = %s\n",idx,Integ.prt(".4f")); // debug print
      double Num = n * Integ.val;
      double Err = fabs(Num) * sqrt( 1./n + SQ(Integ.err/Integ.val) );
      Nos.push_back( ValErr {Num,Err} );
      printf("%7s = %s\n",names[idx].c_str(),Nos[idx].prt("7.1f"));
   }
   double Nphi = Nos[1].val, err_Nphi = Nos[1].err;

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   auto Fun = [Lfit,Fpar](double x) -> double {
      return Lfit( &x,Fpar.data() );
   };
   double pvalueKS = FastKSTest(Fun,mkk);
   printf("pvalueKS= %.3g\n", pvalueKS);

   // Get chisquare of histogram
   auto [ch2_hst,ndof_hst] = GetCh2Ndof(hst,Ffit,norm);
   ndof_hst -= res.NFreeParameters();
   printf("hst: chi^2/Ndof= %.1f / %d\n",ch2_hst,ndof_hst);

   //-----------------------------------------------------------------
   // Functions to draw
   auto Ldr = [n,Fpar,sl](double* x,double* p) -> double {
      const double pp[]
        { Fpar[0],Fpar[1],Fpar[2],Fpar[3],Fpar[4],Fpar[5],sl };
      return bW*n*IntfrBWARGN( x[0], pp, int(p[0]) );
   };
   TF1* ffit = new TF1("ffit", Ldr, dL, dU, 1);
   ffit->SetNpx(500);

   // find max/min to draw
   ffit->SetParameter(0, 1); // BW
   int maxbin = hst->GetMaximumBin();
   double hmax = hst->GetBinContent(maxbin)+hst->GetBinError(maxbin);
   double fmax = ffit->GetMaximum(1.01, 1.03, 1e-6);
   if ( fmax > hmax ) { // negative interference
      hmax = max(3e3,1.1*fmax);
      hst->SetMaximum(hmax);
   }

   ffit->SetParameter(0, 3); // interference
   double fmin = ffit->GetMinimum(1.01, 1.03, 1e-6);
   if ( fmin < 0 ) {
      fmin = min( -0.2e3, 1.2*fmin );
      hst->SetMinimum(fmin);
   }

   //-----------------------------------------------------------------
   auto name = Form("c1_data_intfr");
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst);

   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetMaxDigits(3);
   hst->GetYaxis()->SetTitleOffset(1.1);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");

   TLegend* leg = new TLegend(0.125,0.66,0.355,0.89);
   leg->AddEntry(hst,"Data","LEP");

   ffit->SetParameter(0, 0); // SUM
   ffit->SetLineWidth(2);
   ffit->SetLineStyle(kSolid);
   ffit->SetLineColor(kRed);
   ffit->DrawCopy("SAME");
   leg->AddEntry( ffit->Clone(), "Fitting result", "L");

   // draw components
   ffit->SetLineWidth(1);
   ffit->SetLineStyle(kDashed);

   ffit->SetParameter(0, 1); // BW
   ffit->SetLineColor(kGreen+2);
   ffit->DrawCopy("SAME");
   leg->AddEntry( ffit->Clone(), "Breit-Wigner", "L");

   ffit->SetParameter(0, 2); // Argus
   ffit->SetLineColor(kBlue);
   ffit->DrawCopy("SAME");
   leg->AddEntry( ffit->Clone(), "Argus", "L");

   ffit->SetParameter(0, 3); // interference
   ffit->SetLineColor(kMagenta+1);
   ffit->DrawCopy("SAME");
   leg->AddEntry( ffit->Clone(), "Interference", "L");

   leg->Draw();

   TPaveText* pt = new TPaveText(0.57,0.54,0.892,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#chi^{2} / ndof = %.1f / %d",
            ch2_hst,ndof_hst) );
   pt->AddText( Form("#it{p-value(KS) = %.3f}",pvalueKS) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt->AddText( Form("#sigma = %s MeV", PE.Eform(2,".2f",1e3)) );
   pt->AddText( Form("a = %s",PE.Eform(3,".1f")) );
   pt->AddText( Form("F = %s",PE.Eform(4,".2f")) );
   pt->AddText( Form("#vartheta = %s",PE.Eform(5,".2f")) );
   pt->AddText( Form("N_{#phi} = %.0f #pm %.0f",Nphi,err_Nphi) );
   pt->Draw();
   // TPaveText* ptr = new TPaveText(0.58,0.545,0.892,0.59,"NDC");
   // ptr->SetTextAlign(12);
   // ptr->SetTextFont(42);
   // ptr->Draw();

   gPad->RedrawAxis();
   c1->Update();
   auto pdf = Form("mkk_ifr_%s.pdf", (neg_pos ? "n" : "p") );
   c1->Print(pdf);
}

// {{{2 BW(x)Gauss (*) Argus + Argus(SB)
// {{{3 * 'myFCN_IntSB' class
//--------------------------------------------------------------------
struct myFCN_IntSB {
   const double sl = -1.97; // v709n4

   vector<double> mkk;   // data in central part
   vector<double> sb;    // side-band data
   TH1D* hmcS = nullptr; // side-band MC-signal

   //-----------------------------------------------------------------
   myFCN_IntSB()
   //-----------------------------------------------------------------
   {
      // 1. Get un-binned data in central part and in side-band
      for ( int date : {2021,2012,2009} ) {
         auto file = Form("data_%02ipsip_all.root",date%100);
         auto&& mk= get_mkk_vec(file,0);
         mkk.insert(end(mkk),begin(mk),end(mk));
         auto&& sbb= get_mkk_vec(file,1);
         sb.insert(end(sb),begin(sbb),end(sbb));
      }
      ReoderMkk(this->mkk);

      // 2. Get side-band MC-signal as histogram
      vector<double> mcs;
      for ( int date : {2021,2012,2009} ) {
         auto file = Form("mcsig_kkmc_%02i.root",date%100);
         auto&& mk= get_mkk_vec(file,1);
         mcs.insert(end(mcs),begin(mk),end(mk));
      }
      hmcS = get_hst(mcs,"mcs_sb");
      double nmc = hmcS->Integral();
      hmcS->Scale(1./(nmc*bW));
      printf("MC-signal SB: norm= %.0f n=%zu\n",nmc,mcs.size());
   }

   // minimization function
   // the signature of this operator() MUST be exactly this:
   // * THE EXTENDED MAXIMUM LIKELIHOOD ESTIMATOR
   //   Kai-Feng Chen "FittingMinuitRooFit" p.35
   // * Kahan-Babushka-Neumaier summation
   //   https://en.wikipedia.org/wiki/Kahan_summation_algorithm
   //-----------------------------------------------------------------
   double operator() (const double* p)
   //-----------------------------------------------------------------
   {
      // debug
      // for (int i = 0; i < 10; ++i) {
         // if ( isnan(p[i]) ) {
            // printf("---> WARNING myFCN_IntSB(): p[%i]=NaN\n",i);
         // }
      // }

      // const double& mphi = p[0];
      // const double& gphi = p[1];
      // const double& sigma = p[2];
      const double a_int = fabs(p[3]);
      // const double& F = p[4];
      // const double& angle = p[5];

      const double& Nint = p[6];
      const double& a_bg = p[7];
      const double& Nbg = p[8];
      // const double& Nsl = p[9];

      double res = 2*(Nint+Nbg);
      double cres = 0; // compensation for lost low-order bits of res

      const double pp[] { p[0],p[1],p[2],a_int,p[4],p[5],sl };
      for ( const auto& m : mkk ) {
         double L = Nint * IntfrBWARGN( m,pp,0 );
         L += Nbg * RevArgusN(m,a_bg,sl);
         double lh = (L>0) ? -2*log(L) : FLT_MAX;
         double t = res + lh;
         if ( fabs(res) >= fabs(lh) ) {
            cres += (res - t) + lh;
         } else {
            cres += (lh - t) + res;
         }
         res = t;
      }

      // fit sideband by Argus + SignalMC leaking in sideband
      /*
      double lhsb = 2*(Nbg+Nsl);
      {
         double t = res + lhsb;
         if ( fabs(res) >= fabs(lhsb) ) {
            cres += (res - t) + lhsb;
         } else {
            cres += (lhsb - t) + res;
         }
         res = t;
      }

      for ( const auto& m : sb ) {
         double L = Nbg * RevArgusN(m,a_bg,sl);
         int bin = hmcS->GetXaxis()->FindBin(m);
         L += Nsl * hmcS->GetBinContent(bin);
         double lh = (L>0) ? -2*log(L) : FLT_MAX;
         double t = res + lh;
         if ( fabs(res) >= fabs(lh) ) {
            cres += (res - t) + lh;
         } else {
            cres += (lh - t) + res;
         }
         res = t;
      }
      */

      res += cres; // Correction only applied once in the very end
      return res;
   }
};

// {{{3 * Integrals and their errors for data_IntSB_fit
//--------------------------------------------------------------------
ValErr CalcSBIntEr( const ROOT::Fit::FitResult& res,
      double sl, unsigned int idx )
//--------------------------------------------------------------------
{
   // desired errors:
   constexpr double eps_abs = 1.e-2;
   constexpr double eps_rel = 1.e-6;

   // get parameters and covariation matrix of error
   int Npar = res.NPar();
   vector<double> Fpar = res.Parameters();
   vector<double> Ferr = res.Errors();

   if ( Npar != 9 ) {
      printf("FATAL ERROR in %s: Npar= %i\n",__func__,Npar);
      exit(EXIT_FAILURE);
   }
   Npar = 7; // skip a_bg & Nbg

   int covMatrStatus = res.CovMatrixStatus();
   if ( covMatrStatus != 3 ) {
      printf("WARNING %s: Status of CovMatrix is %i\n",
            __func__,covMatrStatus);
   }
   vector<double> cov_m(Npar*Npar,0);
   for ( int i = 0; i < Npar; ++i ) {
      for ( int j = 0; j < Npar; ++j ) {
         cov_m[i*Npar+j] = res.CovMatrix(i,j);
      }
   }

   // 1) calculate integral of normalized component
   auto Lfc = [sl,idx](const double* x,const double* p) {
      double m = x[0];
      double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl }; // real slope
      double normI = IntfrBWARGN( m,pp,100 );
      pp[6] = 0.; // set slope to zero
      double intfr = p[6]*IntfrBWARG( m,pp,idx ) / normI;
      return intfr;
   };
   TF1 FC("FC", Lfc, dL, dU, Npar);
   FC.SetParameters( Fpar.data() );
   FC.SetParErrors( Ferr.data() ); // for GradientPar()

   double err_num = 0;
   double Integ = FC.IntegralOneDim(dL,dU,eps_rel,eps_abs, err_num);
   // cout << " Integ(" << idx << ") = " << Integ << endl;
   // cout << " err_num(" << idx << ") = " << err_num << endl;
   if ( err_num > eps_abs ) {
      printf("WARNING %s: numerical integration error of F_C(%i)"
            " is too large: %.3g\n",__func__,idx,err_num);
   }

   // 2) calculate error of this integral
   TMatrixDSym covMatrix(Npar);
   covMatrix.Use(Npar,cov_m.data());
   // printf("\ncovMatrix : ");
   // covMatrix.Print();

   TVectorD IntGrad(Npar);
   double err_num2 = 0;
   for ( int i = 0; i < Npar; ++i ) {
      // skip parameters with zero error
      if ( covMatrix(i,i) == 0 ) {
         continue;
      }

      auto LdF = [FC,i](const double* x, const double* p) -> double {
         TF1 tmp(FC); // may modify 'tmp' but not FC
         const double eps = 0.1; // the step used in numerical
                                 // differentiation is eps*par_error
         return tmp.GradientPar(i,x,eps);
      };
      TF1 dF("dF",LdF,dL,dU,0);
      double err_num = 0;
      IntGrad[i] = dF.IntegralOneDim(dL,dU,eps_rel,eps_abs, err_num);
      // cout << " ntGrad[" << i << "] = " << IntGrad[i] << endl;
      err_num = covMatrix(i,i) * IntGrad[i] * err_num;
      err_num2 += SQ(err_num);
      // cout << " abs err_num(" << i << ") = " << err_num << endl;
   }

   double err_Int = sqrt( covMatrix.Similarity(IntGrad) );
   err_num2 = sqrt( err_num2 / err_Int ); // abs numerical error
   // cout << " err_num(" << idx << ") = " << err_num2 << endl;
   if ( err_num2 > 0.1*err_Int ) {
      printf("WARNING %s: numerical integration error of dF(%i)"
            " is too large: %.3g\n",__func__,idx,err_num2);
   }

   return (ValErr {Integ,err_Int});
}

// {{{3 * data_IntSB_fit()
//--------------------------------------------------------------------
void data_IntSB_fit(size_t Cx, size_t Cy)
//--------------------------------------------------------------------
{
   myFCN_IntSB my_fcn;  // class for 'FitFCN'
   const auto& sl = my_fcn.sl; // efficiency parameters
   auto& mkk = my_fcn.mkk;
   auto& sb = my_fcn.sb;

   TH1D* hst = get_hst(mkk,"data_cp");
   TH1D* hsb = get_hst(sb,"data_sb");
   const auto& hmcS = my_fcn.hmcS;

   int n = mkk.size();
   double norm = hst->Integral();
   cout << " norm= " << norm << " n= " << n << endl;
   int nsb = sb.size();
   double norm_sb = hsb->Integral();
   cout << " norm_sb= " << norm_sb << " nsb= " << nsb << endl;

   //-----------------------------------------------------------------
   // == fit configuration
   ROOT::Fit::Fitter fitter;
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   // Possible printing levels are from 0 (minimal) to 3 (maximum)
   fitter.Config().MinimizerOptions().SetPrintLevel(2);
   // fitter.Config().MinimizerOptions().Print(cout);

   vector<string> par_name { "Mphi", "Gphi", "sigma",
      "a_int", "F", "angle", "Nint", "a_bg", "Nbg" };

   vector<double> par_ini =
   {Mphi,Gphi,1e-3,0.,0.9, 0.6,norm-554,5.2,554.}; //n
   // {Mphi,Gphi,1e-3,0.,0.9,-0.6,norm-554,5.2,554.}; //p

   bool neg_pos = par_ini[5] > 0; // true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(Mphi-1e-3, Mphi+1e3);
   // fitter.Config().ParSettings(0).Fix(); // Mphi
   fitter.Config().ParSettings(1).SetLimits(Gphi-1e-3, Gphi+1e-3);
   // fitter.Config().ParSettings(1).Fix(); // Gphi
   fitter.Config().ParSettings(2).SetLimits(0.2e-3, 2.e-3); // sig
   // fitter.Config().ParSettings(3).SetLimits(0.,5); // Argus intfr
   fitter.Config().ParSettings(3).Fix(); // Argus intfr
   fitter.Config().ParSettings(4).SetLimits(0.01, 10.); // F
   fitter.Config().ParSettings(5).SetLimits(-M_PI, M_PI); // angle
   // fitter.Config().ParSettings(6).SetLimits(0,norm); // Nint
   // fitter.Config().ParSettings(7).SetLimits(0.,10.); // Argus bg
   fitter.Config().ParSettings(7).Fix(); // Argus bg
   fitter.Config().ParSettings(8).Fix(); // Nbg

   // == Fit
   // fitter.Config().MinimizerOptions().SetTolerance(0.02);
   int Ndat = mkk.size(); // + sb.size();
   fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // false=likelihood
   if ( fitter.Result().Status() != 0 ) {
      printf("ERROR %s: Status of result is %i\n",
            __func__,fitter.Result().Status());
      exit(EXIT_FAILURE);
   }
   fitter.CalculateHessErrors(); // improving errors
   // MINOS: parameter indices for running Minos
   // fitter.Config().SetMinosErrors({4,5}); // F, angle
   fitter.Config().MinimizerOptions().SetTolerance(0.1); //Minos only
   fitter.CalculateMinosErrors();
   fitter.CalculateHessErrors(); // for correct integral errors

   const ROOT::Fit::FitResult& res = fitter.Result();
   res.Print(cout);
   // res.PrintCovMatrix(cout);

   // double Lmin = res.MinFcnValue();
   ParAndErr PE(res,0.05);
   vector<double>& Fpar = PE.Fpar;
   Fpar[3] = fabs(Fpar[3]); // |a_int|

   vector<string> names { "NKK", "Nphi", "Nnonphi", "Nifr"  };
   vector<ValErr> Nos;
   Nos.reserve(4);
   for ( int idx = 0; idx <= 3; ++idx ) {
      ValErr Integ = CalcSBIntEr( res, sl, idx );
      Nos.push_back( Integ );
      printf("%7s = %s\n",names[idx].c_str(),Nos[idx].prt("7.1f"));
   }

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   auto Fun = [Fpar,sl](double x) -> double {
      const double pp[]
      { Fpar[0],Fpar[1],Fpar[2],Fpar[3],Fpar[4],Fpar[5],sl };
      const double& Nint = Fpar[6];
      const double& a_bg = Fpar[7];
      const double& Nbg = Fpar[8];
      return Nint * IntfrBWARGN( x,pp,0 ) +
         Nbg * RevArgusN(x,a_bg,sl);
   };
   double pvalueKS = FastKSTest(Fun,mkk);
   printf("pvalueKS= %.3g\n", pvalueKS);

   // Chisquare of histogram
   auto Fch = new TF1("Fch",
         [&Fun](double* x,double* p){return Fun(x[0]);}, dL,dU,0);
   auto [ch2_hst,ndof_hst] = GetCh2Ndof(hst,Fch);
   ndof_hst -= res.NFreeParameters();
   printf("hst: chi^2/Ndof= %.1f / %d\n",ch2_hst,ndof_hst);

   //-----------------------------------------------------------------
   // Functions to draw
   auto Ldr = [Fpar,sl](double* x,double* p) -> double {
      int isw = int(p[0]+0.5);
      if ( isw < 0 || isw > 4) return 0;

      const double pp[]
      { Fpar[0],Fpar[1],Fpar[2],Fpar[3],Fpar[4],Fpar[5],sl };
      const double& Nint = Fpar[6];
      const double& a_bg = Fpar[7];
      const double& Nbg = Fpar[8];

      double Argus_bg = bW * Nbg*RevArgusN(x[0],a_bg,sl);
      if ( isw == 4 ) return Argus_bg;

      double BWARGN = bW * Nint*IntfrBWARGN( x[0],pp,isw);
      if ( isw > 0 ) return BWARGN;
      return BWARGN + Argus_bg;
   };
   TF1* ffit = new TF1("ffit", Ldr, dL, dU, 1);
   ffit->SetNpx(500);

   // find max/min to draw
   ffit->SetParameter(0, 1); // BW
   int maxbin = hst->GetMaximumBin();
   double hmax = hst->GetBinContent(maxbin)+hst->GetBinError(maxbin);
   double fmax = ffit->GetMaximum(1.01, 1.03, 1e-6);
   if ( fmax > hmax ) { // negative interference
      hmax = max(3e3,1.1*fmax);
      hst->SetMaximum(hmax);
   }

   ffit->SetParameter(0, 3); // interference
   double fmin = ffit->GetMinimum(1.01, 1.03, 1e-6);
   if ( fmin < 0 ) {
      fmin = min( -0.2e3, 1.2*fmin );
      hst->SetMinimum(fmin);
   }

   //-----------------------------------------------------------------
   auto name = Form("c1_data_intSB");
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst);

   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetMaxDigits(3);
   hst->GetYaxis()->SetTitleOffset(1.1);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");

   TLegend* leg = new TLegend(0.125,0.62,0.355,0.89);
   leg->AddEntry(hst,"Data","LEP");

   ffit->SetParameter(0, 0); // SUM
   ffit->SetLineWidth(2);
   ffit->SetLineStyle(kSolid);
   ffit->SetLineColor(kRed);
   ffit->DrawCopy("SAME");
   leg->AddEntry( ffit->Clone(), "Fitting result", "L");

   // draw components
   ffit->SetLineWidth(1);
   ffit->SetLineStyle(kDashed);

   ffit->SetParameter(0, 1); // BW
   ffit->SetLineColor(kGreen+2);
   ffit->DrawCopy("SAME");
   leg->AddEntry( ffit->Clone(), "Breit-Wigner", "L");
   // leg->AddEntry( ffit->Clone(), "Breit-Wigner #phi#eta", "L");

   ffit->SetParameter(0, 2); // Argus
   ffit->SetLineColor(kBlue);
   ffit->DrawCopy("SAME");
   leg->AddEntry( ffit->Clone(), "Argus", "L");

   ffit->SetParameter(0, 3); // interference
   ffit->SetLineColor(kMagenta+1);
   ffit->DrawCopy("SAME");
   leg->AddEntry( ffit->Clone(), "Interference", "L");

   ffit->SetParameter(0, 4); // Argus(SB)
   ffit->SetLineStyle(kSolid);
   ffit->SetLineColor(kCyan+2);
   ffit->DrawCopy("SAME");
   leg->AddEntry( ffit->Clone(), "Background", "L");

   leg->Draw();

   // TPaveText* pt = new TPaveText(0.56,0.45,0.892,0.89,"NDC");
   TPaveText* pt = new TPaveText(0.58,0.45,0.892,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#it{p-value(KS) = %.3f}",pvalueKS) );
   pt->AddText( Form("#chi^{2} / ndof = %.1f / %d",
            ch2_hst,ndof_hst) );
   pt->AddText( Form("N_{int} = %s",PE.Eform(6,".0f")) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".2f",1e3)) );
   pt->AddText( Form("#sigma = %s MeV", PE.Eform(2,".2f",1e3)) );
   pt->AddText( Form("a_{int} = %s",PE.Eform(3,".1f")) );
   pt->AddText( Form("F = %s",PE.Eform(4,".2f")) );
   pt->AddText( Form("#vartheta = %s",PE.Eform(5,".2f")) );
   pt->AddText( Form("N_{bg} = %s",PE.Eform(8,".0f")) );
   pt->AddText( Form("a_{bg} = %s",PE.Eform(7,".1f")) );

   // const double& Nphi = Nos[1].val;
   // const double& err_Nphi = Nos[1].err;
   // pt->AddText( Form("N_{#phi} = %.0f #pm %.0f",Nphi,err_Nphi) );
   // const double& Nkk = Nos[0].val;
   // const double& err_Nkk = Nos[0].err;
   // pt->AddText( Form("N_{KK} = %.0f #pm %.0f",Nkk,err_Nkk) );
   pt->Draw();

   gPad->RedrawAxis();
   c1->Update();
   auto pdf = Form("mkk_intSB_%s.pdf", (neg_pos ? "n" : "p") );
   c1->Print(pdf);
}

// {{{2 BR as a fitting parameter: BW(x)Gauss (*) Argus + Argus(SB)
// {{{3 * 'myFCN_fitBR' class
//--------------------------------------------------------------------
struct myFCN_fitBR {
   const double sl = -1.97; // v709n4 +/- 0.13
   // const double sl = -1.84; // sys_kp
   // const double sl = -2.10; // sys_km

   vector<double> mkk;   // data in central part
   vector<double> sb;    // side-band data
   TH1D* hmcS = nullptr; // side-band MC-signal

   // Br(eta->2gamma) = 39.36 +/- 0.18% PDG
   const double BReta = 0.3936;
   // Br(phi->K+K-) = 49.1 +/- 0.5% PDG
   const double BRphi = 0.491;

   // v709n4
   const double NJpsi = 390.848e6; // +/- 0.848
   const double e_phi = 0.3380; // eff chi2<100, +/- (4)
   // const double e_phi = 0.3385; // NoHC, chi2<100
   // const double e_phi = 0.3182; // weta=2*seta, chi2<100
   // const double e_phi = 0.3315; // weta=2.5*seta, chi2<100
   // const double e_phi = 0.3416; // weta=3.5*seta, chi2<100
   // const double e_phi = 0.3439; // weta=4*seta, chi2<100

   // conversion Br(J/Psi -> KK eta) -> Nkk
   double Br2Nkk = NJpsi * e_phi * BReta;
   // conversion Br(J/Psi -> phi eta) -> Nphi
   double Br2Nphi = Br2Nkk * BRphi;

   // integrals
   double IeBW = 0, IeAr = 0, IeIc = 0, IeIs = 0;
   double IBW = 0, IAr = 0, IIc = 0, IIs = 0;

   // normalizations for Argus in side-band
   double normArsb = 1;

   //-----------------------------------------------------------------
   myFCN_fitBR()
   //-----------------------------------------------------------------
   {
      // 1. Get un-binned data in central part and in side-band
      for ( int date : {2021,2012,2009} ) {
         auto file = Form("data_%02ipsip_all.root",date%100);
         auto&& mk= get_mkk_vec(file,0);
         mkk.insert(end(mkk),begin(mk),end(mk));
         auto&& sbb= get_mkk_vec(file,1);
         sb.insert(end(sb),begin(sbb),end(sbb));
      }
      ReoderMkk(this->mkk);

      // 2. Get side-band MC-signal as histogram
      vector<double> mcs;
      for ( int date : {2021,2012,2009} ) {
         auto file = Form("mcsig_kkmc_%02i.root",date%100);
         auto&& mk= get_mkk_vec(file,1);
         mcs.insert(end(mcs),begin(mk),end(mk));
      }
      hmcS = get_hst(mcs,"mcs_sb");
      double nmc = hmcS->Integral();
      hmcS->Scale(1./(nmc*bW));
      printf("MC-signal SB: norm= %.0f n=%zu\n",nmc,mcs.size());
   }

   // normalization for Argus SB (see RevArgusN)
   //-----------------------------------------------------------------
   void calcNormAr(const double A)
   //-----------------------------------------------------------------
   {
      // cache for calculated values
      static double A_save = -1.;
      if ( A_save == A ) {
         return;
      }
      A_save = A;

      double p[] = {A,sl};
      auto Lint = [](double x, void* pp) -> double{
         double* p = static_cast<double*>(pp);
         return RevArgus(x,p[0]) * (1+p[1]*(x-1.02));
      };
      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1e-12, 1e-9, 1000);
      normArsb = gsl_int.Integral(Lint,(void *)p,dL,dU);

      if ( gsl_int.Status() ) { // true for error
         double abs_err = gsl_int.Error();
         printf("%s error for A=%.3f: norm=%g, error=%g\n",
               __func__,A,normArsb,abs_err);
      }
   }

   //-----------------------------------------------------------------
   void calcIntegrals(double mphi,double gphi,double sig,double a_int)
   //-----------------------------------------------------------------
   {
      // cach old values
      static double mphi_save = 0., gphi_save = 0., sig_save = 0.;
      static double aint_save = -1.;
      if ( mphi_save == mphi && gphi_save == gphi
            && sig_save == sig && aint_save == a_int ) {
         return;
      }
      mphi_save = mphi;
      gphi_save = gphi;
      sig_save = sig;
      aint_save = a_int;

      const double F = 1.;
      const double ang = 0.; //0(PI/2) for cos(sin) members of Infr.
      // integrand lambda function: mphi,gphi,sig,ar,F,ang,sl,idx
      double pp[] { mphi,gphi,sig,a_int,F,ang,sl, 1 };

      auto Lint = [](double x, void* pp) -> double{
         const double* p = static_cast<const double*>(pp);
         int idx = int(p[7]);
         return IntfrBWARG(x,p,idx);
      };
      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1e-12, 1e-9, 1000);

      // function for check results
      auto checkInt = [&gsl_int](string t, double res, double pp[]) {
         if ( gsl_int.Status() ) { // true for error
            printf("%s: result= %g\n", t.data(),res);
            printf("pp= %g,%g,%g,%g,%g,%g,%g,%g\n",
                  pp[0],pp[1],pp[2],pp[3],pp[4],pp[5],pp[6],pp[7]);
            exit(EXIT_FAILURE);
         }
      };

      // calculate integrals
      pp[6] = sl; // for Ie..

      pp[7] = 1; // idx
      IeBW = gsl_int.Integral(Lint,(void *)pp,dL,dU);
      checkInt("IeBW",IeBW,pp);
      pp[7] = 2;
      IeAr = gsl_int.Integral(Lint,(void *)pp,dL,dU);
      checkInt("IeAr",IeAr,pp);
      pp[7] = 3;
      pp[5] = 0; // ang=0 for cos
      IeIc = gsl_int.Integral(Lint,(void *)pp,dL,dU);
      checkInt("IeIc",IeIc,pp);
      pp[5] = M_PI/2; // ang = pi/2 for sin
      IeIs = gsl_int.Integral(Lint,(void *)pp,dL,dU);
      checkInt("IeIs",IeIs,pp);

      pp[6] = 0; // for I..

      pp[7] = 1; // idx
      IBW = gsl_int.Integral(Lint,(void *)pp,dL,dU);
      checkInt("IBW",IBW,pp);
      pp[7] = 2;
      IAr = gsl_int.Integral(Lint,(void *)pp,dL,dU);
      checkInt("IAr",IAr,pp);
      pp[7] = 3;
      pp[5] = 0; // ang=0 for cos
      IIc = gsl_int.Integral(Lint,(void *)pp,dL,dU);
      checkInt("IIc",IIc,pp);
      pp[5] = M_PI/2; // ang = pi/2 for sin
      IIs = gsl_int.Integral(Lint,(void *)pp,dL,dU);
      checkInt("IIs",IIs,pp);

      // debug print
      // printf(" IeBW= %.3g IeAr= %.3g IeIc= %.3g IeIs= %.3g\n",
      // IeBW, IeAr, IeIc, IeIs);
      // printf(" IBW= %.3g IAr= %.3g IIc= %.3g IIs= %.3g\n",
      // IBW, IAr, IIc, IIs);
   }

   // intergrals MUST be calculeted before calling this function
   //-----------------------------------------------------------------
   tuple<double,double> calcF(const double* p) const
   //-----------------------------------------------------------------
   {
      const double& Brkk = p[0];
      double Nkk = Brkk * Br2Nkk;

      const double& Brphi = p[1];
      double Nphi = Brphi * Br2Nphi;

      const double& ang = p[6];

      // calculate F:
      double A = IAr;
      double B = IIc * cos(ang) + IIs * sin(ang);
      double C = IBW*(1. - Nkk/Nphi);
      double Dis = B*B-4*A*C;
      double Ff = 0., penalty = 0.;
      if ( Dis < 0 ) { // return "minimum"
         Ff = max(0.,-B/(2*A));
         penalty = 1e5*fabs(Dis);
      } else {
         Ff = (-B + sqrt(Dis))/(2*A);
      }
      if ( Ff < 0. ) {
         penalty = -10*Ff;
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
   // * THE EXTENDED MAXIMUM LIKELIHOOD ESTIMATOR
   //   Kai-Feng Chen "FittingMinuitRooFit" p.35
   // * Kahan-Babushka-Neumaier summation
   //   https://en.wikipedia.org/wiki/Kahan_summation_algorithm
   //-----------------------------------------------------------------
   double operator() (const double* p)
   //-----------------------------------------------------------------
   {
      const double& Brkk = p[0];
      const double& Brphi = p[1];
      double Nkk = Brkk * Br2Nkk;
      double Nphi = Brphi * Br2Nphi;

      const double& mphi = p[2];
      const double& gphi = p[3];
      const double& sigma = p[4];
      const double a_int = fabs(p[5]);
      calcIntegrals(mphi,gphi,sigma,a_int);

      const double& ang = p[6];
      auto [Ff,penalty] = calcF(p);

      const double& a_bg = p[7];
      const double& Nbg = p[8];
      const double& Nsl = p[9];
      calcNormAr(a_bg);
      const double NbgNorm = Nbg / normArsb;
      double NFit = Nphi/IBW; // == Nint / normE

      double normE = IeBW+Ff*((IeIc*cos(ang)+IeIs*sin(ang))+Ff*IeAr);
      double Nint = NFit*normE;

      double res = 2*(Nint+Nbg) + penalty;
      double cres = 0; // compensation for lost low-order bits of res

      const double pp[] { mphi,gphi,sigma,a_int,Ff,ang,sl };
      for ( const auto& m : mkk ) {
         double L = NFit * IntfrBWARG( m,pp,0 );
         L += NbgNorm * RevArgus(m,a_bg)*(1+sl*(m-1.02));
         double lh = (L>0) ? -2*log(L) : FLT_MAX;
         double t = res + lh;
         if ( fabs(res) >= fabs(lh) ) {
            cres += (res - t) + lh;
         } else {
            cres += (lh - t) + res;
         }
         res = t;
      }

      // fit SB by Argus
      double lhsb = 2*(Nbg+Nsl);
      {
         double t = res + lhsb;
         if ( fabs(res) >= fabs(lhsb) ) {
            cres += (res - t) + lhsb;
         } else {
            cres += (lhsb - t) + res;
         }
         res = t;
      }

      for ( const auto& m : sb ) {
         double L = NbgNorm * RevArgus(m,a_bg)*(1+sl*(m-1.02));
         int bin = hmcS->GetXaxis()->FindBin(m);
         L += Nsl * hmcS->GetBinContent(bin);
         double lh = (L>0) ? -2*log(L) : FLT_MAX;
         double t = res + lh;
         if ( fabs(res) >= fabs(lh) ) {
            cres += (res - t) + lh;
         } else {
            cres += (lh - t) + res;
         }
         res = t;
      }

      res += cres; // Correction only applied once in the very end
      return res;
   }
};

// {{{3 * data_IntSB_fitBR()
//--------------------------------------------------------------------
void data_IntSB_fitBR(size_t Cx, size_t Cy)
//--------------------------------------------------------------------
{
   myFCN_fitBR my_fcn;  // class for 'FitFCN'
   const auto& sl = my_fcn.sl; // efficiency parameters
   auto& mkk = my_fcn.mkk;
   auto& sb = my_fcn.sb;

   TH1D* hst = get_hst(mkk,"data_cp");
   TH1D* hsb = get_hst(sb,"data_sb");
   const auto& hmcS = my_fcn.hmcS;

   int n = mkk.size();
   double norm = hst->Integral();
   cout << " norm= " << norm << " n= " << n << endl;
   int nsb = sb.size();
   double norm_sb = hsb->Integral();
   cout << " norm_sb= " << norm_sb << " nsb= " << nsb << endl;
   // cout << " sl= " << sl << endl;
   // cout << " R_BW= " << R_BW << endl;

   //-----------------------------------------------------------------
   // == fit configuration
   ROOT::Fit::Fitter fitter;
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   // Possible printing levels are from 0 (minimal) to 3 (maximum)
   fitter.Config().MinimizerOptions().SetPrintLevel(2);

   vector<string> par_name { "Brkk", "Brphi",
      "Mphi", "Gphi", "sig", "a_int", "angle",
      "a_bg", "Nbg", "Nsl" };

   // double MphiMC = 1.01953;
   vector<double> par_ini =
   {4.4e-4,8.9e-4, Mphi,Gphi, 1e-3, 0., 0.7,    // n
   // {4.4e-4,7.9e-4, Mphi,Gphi, 1e-3, 0.,-0.7,    // p
      5.1, 0.77*norm_sb, 0.23*norm_sb};

   bool neg_pos = par_ini[6] > 0; // true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(3e-4, 6e-4);  // Brkk
   fitter.Config().ParSettings(1).SetLimits(5e-4, 12e-4); // Brphi
   fitter.Config().ParSettings(2).SetLimits(Mphi-1e-3, Mphi+1e3);
   fitter.Config().ParSettings(3).SetLimits(Gphi-1e-3, Gphi+1e-3);
   // fitter.Config().ParSettings(3).Fix();
   fitter.Config().ParSettings(4).SetLimits(0.5e-3,2.e-3);// sigma
   // fitter.Config().ParSettings(5).SetValue(1.7); // sys a_int
   fitter.Config().ParSettings(5).Fix(); // a_int
   fitter.Config().ParSettings(6).SetLimits(-M_PI, M_PI); // angle
   fitter.Config().ParSettings(7).SetLimits(0., 10);// a_bg
   // fitter.Config().ParSettings(8).SetLimits(416.,692.); // Nbg
   // fitter.Config().ParSettings(9).Fix(); // Nsl

   // == Fit
   int Ndat = mkk.size() + sb.size();
   fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // false=likelihood
   fitter.CalculateHessErrors(); // fast
   // MINOS: parameter indices for running Minos
   // fitter.Config().SetMinosErrors({0,1,6}); // BrKK, BrPhi, angle
   fitter.Config().SetMinosErrors({0,1}); // BrKK, BrPhi
   fitter.Config().MinimizerOptions().SetTolerance(0.1); //Minos only
   fitter.CalculateMinosErrors();

   const ROOT::Fit::FitResult& res = fitter.Result();
   ParAndErr PE(res,0.05); // ignore 5% upper/lower errors
   vector<double>& Fpar = PE.Fpar;
   Fpar[5] = fabs(Fpar[5]); // |a_int|

   // variables for KS-test and draw-function
   auto [Ff,penalty] = my_fcn.calcF(Fpar.data());
   my_fcn.calcNormAr(Fpar[7]); // a_bg
   const auto& normArsb = my_fcn.normArsb;

   // print result
   if ( !isfinite(penalty) || penalty > 0. ) {
      printf("\n=> INVALID RESULT: penalty= %f <=\n",penalty);
   } else {
      printf("\n=> FINAL RESULT <=\n");
   }
   res.Print(cout);
   // res.PrintCovMatrix(cout);
   // printf("debug print: Ff=%g\n",Ff);
   // printf("debug print: Mphi=%.4f MeV\n",Fpar[2]*1e3);
   // double Lmin = res.MinFcnValue();

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   // central part
   double Nkk = Fpar[0] * my_fcn.Br2Nkk;
   double Nphi = Fpar[1] * my_fcn.Br2Nphi;
   // printf("debug print: Nkk= %g, Nphi=%g\n",Nkk,Nphi);
   double NFit = Nphi/my_fcn.IBW; // == Nint/normE

   const double& mphi = Fpar[2];
   const double& gphi = Fpar[3];
   const double& sig = Fpar[4];
   const double& a_int = Fpar[5];
   const double& ang = Fpar[6];

   // {mphi,gphi,sigma,A,F,Ang,sl}
   const double pp[] { mphi, gphi, sig, a_int, Ff, ang, sl };

   // normalization on 1
   const double& a_bg = Fpar[7];
   const double& Nbg = Fpar[8];
   double Ntot = Nkk + Nbg;
   NFit /= Ntot;

   double NbgNorm = Nbg/normArsb;
   NbgNorm /= Ntot;
   auto Fcr = [NFit,pp,NbgNorm,a_bg,sl](double x) {
      double ret = NFit*IntfrBWARG( x,pp,0 );
      ret += NbgNorm * RevArgus(x,a_bg)*(1+sl*(x-1.02));
      return ret;
   };
   double pvalueKS = FastKSTest(Fcr,mkk);
   printf("pvalueKS(cr)= %.3g\n", pvalueKS);

   // Side-band
   double Nsl = Fpar[9]; // signal leak to sideband
   Ntot = Nbg + Nsl;
   NbgNorm = Nbg / normArsb;
   NbgNorm /= Ntot;
   Nsl /= Ntot;
   auto Fsb = [NbgNorm,a_bg,sl,Nsl,hmcS](double x) {
      double ret = NbgNorm * RevArgus(x,a_bg)*(1+sl*(x-1.02));
      int bin = hmcS->GetXaxis()->FindBin(x);
      ret += Nsl * hmcS->GetBinContent(bin);
      return ret;
   };
   double pvalueKSsb = FastKSTest(Fsb,sb);
   printf("pvalueKS(sb)= %.3g\n", pvalueKSsb);

   // Chisquare of histograms
   auto FcrT = new TF1("FcrT",
         [&Fcr](double* x,double* p){return Fcr(x[0]);}, dL,dU,0);
   auto [ch2_hst,ndof_hst] = GetCh2Ndof(hst,FcrT,norm);
   ndof_hst -= res.NFreeParameters() - 1; // Nsl for SB only
   printf("hst: chi^2/Ndof= %.1f / %d\n",ch2_hst,ndof_hst);

   auto FsbT = new TF1("FsbT",
         [&Fsb](double* x,double* p){return Fsb(x[0]);}, dL,dU,0);
   auto [ch2_hsb,ndof_hsb] = GetCh2Ndof(hsb,FsbT,norm_sb);
   ndof_hsb -= 3; // N_bg, a_bg and Nsl
   printf("hsb: chi^2/Ndof= %.1f / %d\n",ch2_hsb,ndof_hsb);

   //-----------------------------------------------------------------
   // Functions to draw
   // Ff=Ff to avoid warning about 'captured structured bindings'
   auto Ldr = [Fpar,my_fcn,Ff=Ff](double* x,double* p) -> double {
      int isw = int(p[0]+0.5);

      const double& sl = my_fcn.sl;
      const double& normArsb = my_fcn.normArsb;

      const double& a_bg = Fpar[7];
      const double& Nbg = Fpar[8];
      double NbgNorm = Nbg / normArsb;
      double Argus_bg = bW *
         NbgNorm * RevArgus(x[0],a_bg)*(1+sl*(x[0]-1.02));
      if ( isw == 4  ) return Argus_bg;

      if ( isw < 4 ) {
         double Nphi = Fpar[1] * my_fcn.Br2Nphi;
         double NFit = Nphi/my_fcn.IBW;
         const double& mphi = Fpar[2];
         const double& gphi = Fpar[3];
         const double& sig = Fpar[4];
         const double& a_int = Fpar[5];
         const double& ang = Fpar[6];
         const double pp[] { mphi, gphi, sig, a_int, Ff, ang, sl };
         double BWARG = bW * NFit * IntfrBWARG( x[0],pp,isw );
         if ( isw == 0 ) return BWARG + Argus_bg;
         return BWARG;
      }

      const double& Nsl = Fpar[9];
      const auto& hmcS = my_fcn.hmcS;
      int bin = hmcS->GetXaxis()->FindBin(x[0]);
      double Nmc_sig = bW * Nsl*hmcS->GetBinContent(bin);
      if ( isw == 5 ) return Nmc_sig;
      return Nmc_sig + Argus_bg;
   };
   TF1* fdr = new TF1("fdr", Ldr, bL, dU, 1);
   fdr->SetNpx(500);

   // find max/min to draw
   fdr->SetParameter(0, 1); // BW
   int maxbin = hst->GetMaximumBin();
   double hmax = hst->GetBinContent(maxbin)+hst->GetBinError(maxbin);
   double fmax = fdr->GetMaximum( 1.01, 1.03, 1e-6);
   if ( fmax > hmax ) { // negative interference
      hmax = max(3e3,1.1*fmax);
      hst->SetMaximum(hmax);
   }

   fdr->SetParameter(0, 3); // interference
   double fmin = fdr->GetMinimum( 1.01, 1.03, 1e-6);
   if ( fmin < 0 ) {
      fmin = min( -0.2e3, 1.2*fmin );
      hst->SetMinimum(fmin);
   }

   //-----------------------------------------------------------------
   auto name = Form("c1_data_fitBR");
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst);
   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetMaxDigits(3);
   hst->GetYaxis()->SetTitleOffset(1.1);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");

   TLegend* leg = new TLegend(0.125,0.62,0.355,0.89);
   leg->AddEntry(hst,"Data","LEP");

   fdr->SetParameter(0, 0); // SUM
   fdr->SetLineWidth(2);
   fdr->SetLineStyle(kSolid);
   fdr->SetLineColor(kRed);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Fitting result", "L");

   // draw components
   fdr->SetLineWidth(1);
   fdr->SetLineStyle(kDashed);

   fdr->SetParameter(0, 1); // BW
   fdr->SetLineColor(kGreen+2);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Breit-Wigner", "L");

   fdr->SetParameter(0, 2); // Argus
   fdr->SetLineColor(kBlue);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Argus", "L");

   fdr->SetParameter(0, 3); // interference
   fdr->SetLineColor(kMagenta+1);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Interference", "L");
   leg->Draw();

   fdr->SetParameter(0, 4); // Argus(SB)
   fdr->SetLineStyle(kSolid);
   fdr->SetLineColor(kCyan+2);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Background", "L");

   // TPaveText* pt = new TPaveText(0.56,0.45,0.892,0.89,"NDC");
   TPaveText* pt = new TPaveText(0.57,0.40,0.892,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#it{p-value(KS) = %.3f}",pvalueKS) );
   pt->AddText( Form("#chi^{2} / ndof = %.1f / %d",
            ch2_hst,ndof_hst) );
   // the space after -4 is intentional
   pt->AddText( Form("Br(KK#eta)= %s #times10^{-4 }",
            PE.Eform(0,".2f",1e4)) );
   pt->AddText( Form("Br(#phi#eta)= %s #times10^{-4}",
            PE.Eform(1,".2f",1e4)) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(2,".2f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(3,".2f",1e3)) );
   pt->AddText( Form("#sigma = %s MeV", PE.Eform(4,".2f",1e3)) );
   pt->AddText( Form("a_{int} = %s",PE.Eform(5,".1f")) );
   pt->AddText( Form("#vartheta = %s",PE.Eform(6,".2f")) );
   pt->AddText( Form("N_{bg} = %s",PE.Eform(8,".0f")) );
   pt->AddText( Form("a_{bg} = %s",PE.Eform(7,".1f")) );
   pt->Draw();

   gPad->RedrawAxis();
   c1->Update();
   auto pdf = Form("mkk_fitBR_%s.pdf", (neg_pos ? "n" : "p"));
   // auto pdf = Form("mkk_fitBR_%s_km.pdf", (neg_pos ? "n" : "p"));
   // auto pdf = Form("mkk_fitBR_%s_nohc.pdf", (neg_pos ? "n" : "p"));
   // auto pdf = Form("mkk_fitBR_%s_r4.pdf", (neg_pos ? "n" : "p"));
   // auto pdf = Form("mkk_fitBR_%s_a17.pdf", (neg_pos ? "n" : "p"));
   // auto pdf = Form("mkk_fitBR_%s_w35.pdf", (neg_pos ? "n" : "p"));
   c1->Print(pdf);

   // Side-band
   //-----------------------------------------------------------------
   auto name2 = Form("c2_fitBR_SB");
   TCanvas* c2 = new TCanvas(name2,name2,Cx/2,0,Cx,Cy);
   c2->cd();
   gPad->SetGrid();

   SetHstFace(hsb);
   double hsb_max = hsb->GetMaximum();
   hsb->SetMaximum( 5*floor(1.4*hsb_max/5) );
   hsb->GetXaxis()->SetTitleOffset(1.1);
   hsb->GetYaxis()->SetTitleOffset(1.1);
   hsb->SetLineWidth(2);
   hsb->SetLineColor(kBlack);
   hsb->SetMarkerStyle(20);

   hsb->Draw("EP");

   TLegend* leg2 = new TLegend(0.57,0.70,0.892,0.89);
   leg2->AddEntry(hsb,"Side-band data","LEP");

   fdr->SetLineWidth(2);
   fdr->SetLineStyle(kSolid);

   fdr->SetParameter(0, 6); // Fitting result
   fdr->SetLineColor(kRed);
   fdr->DrawCopy("SAME");
   leg2->AddEntry(fdr->Clone(),"Fitting result","L");

   fdr->SetParameter(0, 4); // Argus SB
   fdr->SetLineColor(kCyan+2);
   fdr->DrawCopy("SAME");
   leg2->AddEntry(fdr->Clone(),"Argus function","L");

   fdr->SetParameter(0, 5); // MC signal in side-band
   fdr->SetLineWidth(1);
   fdr->SetLineStyle(kDashed);
   fdr->SetLineColor(kGreen+3);
   fdr->DrawCopy("SAME");
   leg2->AddEntry(fdr->Clone(),"MC signal #phi#eta","L");
   leg2->Draw();

   TPaveText* ptsb = new TPaveText(0.57,0.505,0.892,0.695,"NDC");
   ptsb->SetTextAlign(12);
   ptsb->SetTextFont(42);
   ptsb->AddText( Form("#it{p-value(KS) = %.3f}", pvalueKSsb) );
   ptsb->AddText( Form("a_{bg} = %s",PE.Eform(7,".1f")) );
   ptsb->AddText( Form("N_{bg} = %s",PE.Eform(8,".0f")) );
   ptsb->AddText( Form("N_{sl} = %s",PE.Eform(9,".0f")) );
   ptsb->Draw();

   gPad->RedrawAxis();
   c2->Update();
   auto pdf2 = Form("mkk_fitBR_%s_SB.pdf", (neg_pos ? "n" : "p"));
   // auto pdf2 = Form("mkk_fitBR_%s_SB_km.pdf", (neg_pos ? "n" : "p"));
   // auto pdf2 = Form("mkk_fitBR_%s_SB_nohc.pdf", (neg_pos ? "n" : "p"));
   // auto pdf2 = Form("mkk_fitBR_%s_SB_r4.pdf", (neg_pos ? "n" : "p"));
   // auto pdf2 = Form("mkk_fitBR_%s_SB_a17.pdf", (neg_pos ? "n" : "p"));
   // auto pdf2 = Form("mkk_fitBR_%s_SB_w35.pdf", (neg_pos ? "n" : "p"));
   c2->Print(pdf2);
}

// {{{2  ? dataSB_Intfr_scan()
#ifdef OLD_IFRSCAN
//--------------------------------------------------------------------
void dataSB_Intfr_scan(int date, string sname) {
//--------------------------------------------------------------------
   string fname( Form("data_%02ipsip_all.root",date%100) );

   myFCN_inter my_fcn;             // class for 'FitFCN'
   const int FunSB = my_fcn.funSB; // 0 - constant, 1 - Argus
   const double sl = my_fcn.sl;    // efficiency parameters

   // Get un-binned data
   TH1D* hist[2];
   my_fcn.mkk = get_mkk_hist( fname, "mkk_cp", &hist[0] );
   TH1D* hst = hist[0];
   const vector<double>& mkk = my_fcn.mkk;
   int n = mkk.size();

   my_fcn.sb = get_mkk_hist( fname, "mkk_sb", &hist[1], 1 );
   TH1D* hsb = hist[1];
   const vector<double>& sb = my_fcn.sb;
   int nsb = sb.size();

   //-----------------------------------------------------------------
   // Fit data
   ROOT::Fit::Fitter fitter;
   // fitter.Config().MinimizerOptions().SetPrintLevel(3);

   vector<string> par_name
      { "Mphi","Gphi","sig","a","F","angle","NKK","Nbg","asb" };

   // parameters of loop for angles (degrees)
   int angmin = -90, angmax = 90, angstep = 181; // test

   // initial values for scan
   vector<double> par_ini;
   if ( date == 2009 ) {
      angmin = -90; angmax = 90; angstep = 1;
      par_ini = {Mphi,Gphi,1.4e-3,0., 1.6, 0., 890., 10., 7. };
   } else if ( date == 2012 ) {
      angmin = -90; angmax = 90; angstep = 1;
      par_ini = {Mphi,Gphi,1.1e-3,0., 1.5, 0., 2750.,35., 4. };
   } else if ( date == 2021 ) {
      angmin = -90; angmax = 90; angstep = 1;
      par_ini = {Mphi,Gphi,1.1e-3,0., 1.5, 0., 17500.,195.,6.};
   }

   const unsigned int Npar = par_name.size(); // number of parameters

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetValue(1.01952);        // like MC
   fitter.Config().ParSettings(0).Fix();                    // Mphi
   fitter.Config().ParSettings(1).Fix();                    // Gphi
   fitter.Config().ParSettings(2).SetLimits(0.5e-3, 2.e-3); // sig
   fitter.Config().ParSettings(3).Fix();                    // a
   fitter.Config().ParSettings(4).SetLimits(0.01, 10.);     // F
   // fitter.Config().ParSettings(6).SetLimits(0., 1.5*n);     // NKK
   // fitter.Config().ParSettings(7).SetLimits(0., 2.*nsb);    // Nbg
   if ( FunSB == 0 ) {                                      // asb
      fitter.Config().ParSettings(8).SetValue(0.);
      fitter.Config().ParSettings(8).Fix();
   } else {
      fitter.Config().ParSettings(8).SetLimits(0.,15.);
   }

   //-----------------------------------------------------------------
   // Functions to draw
   vector<double> Fpar;
   auto Ldr = [&Fpar,sl,FunSB](double* x,double* p) -> double {
      int isw = int(p[0]+0.5);

      double pp[]
         {Fpar[0],Fpar[1],Fpar[2],Fpar[3],Fpar[4],Fpar[5],sl};
      double BWARG = Fpar[6] * IntfrBWARGN(x[0],pp,isw);
      if ( isw > 0 && isw < 4 ) { return bW * BWARG; }

      double Bg = 0;
      if ( FunSB == 0 ) {
         Bg = Fpar[7]/(dU-dL);
      } else {
         Bg = Fpar[7] * RevArgusN( x[0],Fpar[8],sl );
      }
      if ( isw == 4 ) { return bW * Bg; }

      return bW * (BWARG + Bg);
   };
   TF1* fdr = new TF1("fdr", Ldr, dL, dU, 1);
   fdr->SetNpx(500);

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   sname += (FunSB == 0 ) ? "_l" : "_ar";
   string hfile = sname + ".root";
   TFile c_out(hfile.c_str(),"recreate"); // file for scan results
   c_out.cd();

   string pdf = sname + ".pdf";
   c1->Print((pdf+"[").c_str()); // open pdf-file

   SetHstFace(hst);
   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetTitleOffset(1.3);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   //-----------------------------------------------------------------
   // start the cycle
   vector<double> angl, Lmin, nphi;
   for( int ang = angmin; ang <= angmax; ang+=angstep ) {

      fitter.Config().ParSettings(5).SetValue(M_PI*ang/180.); // angle
      fitter.Config().ParSettings(5).Fix();

      // == Fit
      int Ndat = n + nsb;
      fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // likelihood
      // fitter.CalculateMinosErrors();
      // fitter.CalculateHessErrors();

      const ROOT::Fit::FitResult& res = fitter.Result();
      Fpar = res.Parameters();
      // res.Print(cout);
      // res.PrintCovMatrix(cout);
      ParAndErr PE(res);

      double lmin = res.MinFcnValue();
      angl.push_back(ang);
      Lmin.push_back(lmin);

      // calculate Nphi
      double Nfit = Fpar[6];
      ValErr Integ = CalcIntErSb( res, sl, 1,true );
      double Nphi = Nfit * Integ.val;
      nphi.push_back(Nphi);

      //--------------------------------------------------------------
      // Draw intermidiate results
      TLegend* leg = new TLegend(0.58,0.725,0.89,0.89);

      hst->Draw("EP");
      leg->AddEntry(hst,Form("Data %i",date),"LEP");

      fdr->SetParameter(0, 0); // SUM
      fdr->SetLineWidth(2);
      fdr->SetLineColor(kRed+1);
      fdr->SetLineStyle(kSolid);
      fdr->DrawCopy("SAME");
      leg->AddEntry( fdr->Clone(), "Result of fit", "L");

      fdr->SetParameter(0, 1); // BW
      fdr->SetLineWidth(2);
      fdr->SetLineColor(kGreen+2);
      fdr->SetLineStyle(kDashed);
      fdr->DrawCopy("SAME");
      leg->AddEntry( fdr->Clone(), "Breit-Wigner #phi#eta", "L");

      fdr->SetParameter(0, 2); // Argus
      fdr->SetLineWidth(2);
      fdr->SetLineColor(kBlue);
      fdr->SetLineStyle(kDashed);
      fdr->DrawCopy("SAME");
      leg->AddEntry( fdr->Clone(), "Non-#phi KK#eta", "L");

      fdr->SetParameter(0, 3); // interference
      fdr->SetLineWidth(2);
      fdr->SetLineColor(kMagenta+1);
      fdr->SetLineStyle(kDashed);
      fdr->DrawCopy("SAME");
      leg->AddEntry( fdr->Clone(), "Interference", "L");
      leg->Draw();

      TPaveText* pt = new TPaveText(0.58,0.48,0.89,0.72,"NDC");
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->AddText( Form("#vartheta = %s deg",
               PE.Eform(5,".0f",180./M_PI)) );
      pt->AddText( Form("#it{L_{min} = %.1f}",lmin) );
      pt->AddText( Form("M_{#phi}= %s MeV", PE.Eform(0,".2f",1e3)) );
      pt->AddText( Form("#Gamma_{#phi}= %s MeV",
               PE.Eform(1,".3f",1e3)) );
      pt->AddText( Form("#sigma = %s MeV", PE.Eform(2,".2f",1e3)) );
      // pt->AddText( Form("a = %s",PE.Eform(3,".1f")) );
      pt->AddText( Form("F = %s",PE.Eform(4,".2f")) );
      pt->AddText( Form("NKK(fit) = %s", PE.Eform(6,".1f")) );
      pt->AddText(Form("N_{#phi } = %.1f",Nphi));
      pt->AddText( Form("#lower[-0.1]{Nbg = %s}",PE.Eform(7,".1f")) );
      if ( FunSB != 0 ) {
         pt->AddText( Form("#lower[-0.1]{a(sb) = %s}",
                  PE.Eform(8,".1f")) );
      }
      pt->Draw();

      gPad->RedrawAxis();
      c1->Update();
      c1->Print(pdf.c_str()); // add to pdf-file
   }
   //-----------------------------------------------------------------
   // final draw
   int nch = angl.size();
   if ( nch < 3 ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
      return;
   }

   // 1 - Nphi vs angl
   c1->cd();
   auto grN = new TGraph( nch, angl.data(), nphi.data() );
   grN->SetName( Form("Nphi_scan_%02i",date%100) );
   grN->SetTitle(";#vartheta, degrees;N_{#phi}");
   grN->SetMarkerColor(kBlue);
   grN->SetMarkerStyle(21);
   grN->SetLineWidth(2);
   grN->GetYaxis()->SetTitleOffset(1.);
   grN->Draw("APL");
   grN->Write();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   // 2 - Lmin vs angle
   // subtract minimal value
   double minL = *(min_element(Lmin.begin(),Lmin.end()));
   for( auto &l : Lmin ) {
      l -= minL;
   }

   auto grL = new TGraph( nch, angl.data(), Lmin.data() );
   grL->SetTitle(";#vartheta, degrees;#it{-2log(L/L_{max})}");
   grL->GetYaxis()->SetMaxDigits(3);
   grL->GetYaxis()->SetTitleOffset(1.);
   grL->SetMarkerColor(kBlue);
   grL->SetMarkerStyle(20);
   grL->SetLineWidth(2);
   grL->Draw("APL");
   grL->SetName( Form("L_scan_%02i",date%100) );
   grL->Write();
   gPad->RedrawAxis();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
   c_out.Close(); // close root-file with scan results
}
#endif

// {{{1 MAIN:
//--------------------------------------------------------------------
void mass_kk_fit()
//--------------------------------------------------------------------
{
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(42);
   gStyle->SetLegendFont(42);

   //-----------------------------------------------------------------
   // = define GSL error handler which does nothing =
   // gsl_set_error_handler_off();

   // = set handler =
   gsl_set_error_handler( my_gslhandler );

   //-----------------------------------------------------------------
   // set integrator: ROOT::Math::GSLIntegrator adaptive method (QAG)
   ROOT::Math::IntegratorOneDimOptions::
      SetDefaultIntegrator("Adaptive");
   ROOT::Math::IntegratorOneDimOptions::
      SetDefaultRelTolerance(1e-9);
   ROOT::Math::IntegratorOneDimOptions::
      SetDefaultWKSize(2000);

   //-----------------------------------------------------------------
   // default Minimizer
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

   //========================================================
   // set the name of the folder with the root files
   Dir = "prod_v709n4/";
   // Dir = "prod_v709n4/NoHC/";
   //========================================================

   size_t Cx = 800, Cy = 640; // canvas sizes, X/Y = 1.25

   // ------------- signal sec.6.1 -----------------------------------
   // test_BW(Cx,Cy);
   // test_BWGN(Cx,Cy);

   // * Fit Mkk(MC-truth) by Breit-Wigner, Fig 17
   // for ( int date : {2009,2012,2021} ) {
      // string pdf_mcsig; //( Form("mcmkk%02i_sig.pdf",date%100) );
      // sig_fit_mc(date, pdf_mcsig, Cx, Cy);
   // }
   // sig_fit_mc(0, "mcmkk_all.pdf", Cx, Cy);

   // * Fit Mkk(MC-reconstructed) by BW x Gauss, Fig 18
   // -> or by BW x (Gauss1+Gauss2), just a test...
   // for ( int date : {2009,2012,2021} ) {
      // string pdf_mc( Form("mkk%02i_sig_log.pdf",date%100) );
      // sig_fit(date,pdf_mc,Cx,Cy);
   // }
   // sig_fit(0,"mkk_sig_all_log.pdf",Cx,Cy);


   // ------------- background sec.6.2 -------------------------------
   // test_RevArgus(Cx,Cy);
   // test_RevArgusN(Cx,Cy);
   // test_Argus_smallA(Cx,Cy);

   // * Fit data mkk side-band by Argus + mcsig(sb)
   // sideband_fit("mkk_sb_fit.pdf",Cx,Cy);

   // ------------- data: no interference sec.6.3.1 ------------------
   // * Fit data mkk by BW(x)Gauss + Argus
   // data_fit(Cx,Cy);

   // ------------- data: interference sec 6.3.2 ---------------------
   // test_Intfr(Cx,Cy);

   // --- Fit data mkk by BW(x)Gauss (*) Argus
   // --- data_Intfr_fit(Cx,Cy);

   // * Fit data mkk by BW(x)Gauss (*) Argus(Int) + Argus(SB)
   // data_IntSB_fit(Cx,Cy);

   // * BR as a fitting parameter: BW(x)Gauss (*) Argus + Argus(SB)
   data_IntSB_fitBR(Cx,Cy);

   // ++ ? scan
   // string pdf_sc( Form("mkk%02i_ifrSB_scan",date%100) );
   // dataSB_Intfr_scan(date,pdf_sc);

   // ------------- completion of the program  ---------------
   my_gslhandler(nullptr,__func__,0,0);
}
