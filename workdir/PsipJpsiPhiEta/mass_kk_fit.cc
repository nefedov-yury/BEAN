// unbinned LH fit of M(K+K-) distributions
// after the cuts (see cuts.h): Mrec + chi^2(4C) + Mgg
// Formulas for:
// *) Breit-Wigner for psi -> K+K- (psi from J/Psi -> phi eta)
// *) convolution of BW with the Gauss distribution
// *) ARGUS function (and "reverse" Argus)
// *) interference B-W and Argus amplitudes
//    -> mkkYY_fit??.pdf
//    -> mkk_cf_???.pdf (combined fit)

#include "gsl/gsl_errno.h"   // GSL error handler

// {{{1 ===== Handler for GSL exceptions =====
static char gsl_integral_func[91] {"not used"};
// histograms for debug
static bool DEBUG = false;
static TH1D* htt[10];

//----------------------------------------------------------------------
void my_handler(const char* reason, const char* file, int line,
                int gsl_errno) {
//----------------------------------------------------------------------
   printf("-- GSL EXCEPTION: %s in func %s\n",reason,gsl_integral_func);
   return;
}

// {{{1 helper functions and constants
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
};

//----------------------------------------------------------------------
struct ValErr {
   double val;
   double err;
   const char* prt(string fm) const {
      string fmt = "%" + fm + " +/- %" + fm;
      return Form(fmt.c_str(),val,err);
   }
};

//----------------------------------------------------------------------
constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

static const double Mjpsi = 3.096916; // 3096.916 +/- 0.011 MeV
static const double Mphi  = 1.019461; //1019.461  +/- 0.016 MeV
static const double Gphi  = 4.247e-3; //   4.249  +/- 0.013 MeV
static const double Meta  = 0.547862; // 547.862  +/- 0.017 MeV
static const double Mk    = 0.493677; // 493.677  +/- 0.016 MeV

static const double Mjpsi2 = SQ(Mjpsi);
static const double Meta2  = SQ(Meta);
static const double Mk2    = SQ(Mk);
// the meson radial parameter in Blatt-Weisskopf form factor:
static const double R_BW   = 3.; // GeV^{-1} for
// static const double R_BW   = 4.; // GeV^{-1} for

static const double dL = 2*Mk; // the left cutoff = 0.987354

// binning
static const int Nbins = 100;
static const double bL = 0.98; // first bin < dL
// ------ Breit-Wigner: bin width 1.0 MeV ------
static const double dU = 1.08; // MUST BE < 1.0835 !!!
// ---------------------------------------------
static const double bW = (dU-bL)/Nbins; // bin width

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

// {{{1 DATA
//-------------------------------------------------------------------------
// DATA
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
vector<double> get_mkk_hist( string fname, string hname,
                             TH1D* hst[], int type = 0 ) {
//-------------------------------------------------------------------------
#include "cuts.h"
   bool isMC = (fname.find("mcsig") != string::npos);

   // name of folder with root files
//    static string dir("prod-9/");
   static string dir("prod-10/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory -> Get("a4c");

   TCut c_here = c_Mrec+c_chi2;
   c_here += TCut(Form("%f<=Mkk&&Mkk<%f",bL,dU));
   if ( type == 0 ) {        // central part
      c_here += c_cpgg;
   } else if ( type == 1 ) { // side-band
      c_here += c_sbgg;
   }
   if ( isMC ) {
      cout << " THIS IS MONTE CARLO EVENTS" << endl;
      c_here += c_MCmkk;
   }

   int n = a4c -> Draw("Mkk", c_here, "goff");
//    cout << " n= " << n << endl;
   double* buffer = a4c -> GetVal(0);
   vector<double> mkk(buffer, buffer+n);

   string title(";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2}");
   int ibW = int(bW*1e5+0.1);
   title += string(";Entries / ") +
      ( (ibW%10 == 0) ? string(Form("%.0f MeV/c^{2}",bW*1e3))
                      : string(Form("%.2f MeV/c^{2}",bW*1e3)) );
   hst[0] = new TH1D(hname.c_str(),title.c_str(),Nbins,bL,dU);

   hst[0] -> Sumw2(true);
   for ( const auto& mk : mkk ) {
      hst[0] -> Fill(mk);
   }
   return mkk;
}

//-------------------------------------------------------------------------
TH1D* plot_Mkk( string fname, string hname, int type = 0 ) {
//-------------------------------------------------------------------------
   TH1D* hist[1];
   get_mkk_hist( fname, hname, hist, type );
   return hist[0];
}

//-------------------------------------------------------------------------
TH1D* get_mcMkk(string fname) {
//-------------------------------------------------------------------------
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( !froot ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TH1D* mcmkk = (TH1D*) gDirectory->Get("mc_Mkk");
   if( !mcmkk ) {
      cerr << "can not find mc_Mkk" << endl;
      exit(0);
   }

   string title(";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2}");
   title += string(";Entries/1.0 MeV/c^{2}");
   mcmkk->SetTitle(title.c_str());
//    mcmkk->Sumw2(true);

   return mcmkk;
}

// {{{1 Breigt Wigner for phi -> KK
//----------------------------------------------------------------------
// Breigt Wigner for phi -> KK
//----------------------------------------------------------------------

//----------------------------------------------------------------------
double BreitWigner(double m, double mphi, double gphi, double R = R_BW) {
//----------------------------------------------------------------------
// see BAM-00117:  PhysRevD.91.112001.pdf (arXiv:1504.03194v2)
//
// we assume that mphi and gphi are the fit parameters

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
   double BB_phi = sqrt( (1+SQ(R*p_phi0)) / (1+SQ(R*p_phi)) );
   double BB_k2  = (1+SQ(R*p_k0)) / (1+SQ(R*p_k));
   double BB_k = sqrt( BB_k2 );

   double fm = r_phi*r_k*BB_phi*BB_k;

//--    double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k; // == m*G(m)
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

   if ( DEBUG ) {
      // int status = gsl_int.Status();
      // the estimate of the absolute Error
      if ( gsl_int.Error() > 1e-6 ) {
         htt[0]->Fill(log10(gsl_int.Error()));
      } else {
         htt[0]->Fill(-7.);
      }
   }

   return Ngauss * result;
}

//----------------------------------------------------------------------
double BreitWignerGaussN( double m,
                          double mphi, double gphi, double sigma,
                          double slope
                        ) {
//----------------------------------------------------------------------
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
      ROOT::Math::GSLIntegrator gsl_int(1.e-7, 1.e-6, 1000);
      norm = gsl_int.Integral(Lint,p,dL,dU);
      if ( DEBUG ) {
         htt[1]->Fill(norm);
         if ( gsl_int.Error() > 1e-6 ) {
            htt[2]->Fill(log10(gsl_int.Error()));
         } else {
            htt[2]->Fill(-7.);
         }
      }
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

   // 1) BreitWigner (internally normalized to 1)
   auto Lbw = [](const double* x,const double* p) -> double {
      return p[0]*BreitWigner(x[0],p[1],p[2],p[3]);
   };
   TF1* bw = new TF1("bw", Lbw, dL, dU, 4);
   bw -> SetParNames("Norm","Mphi","Gphi","Rff");
   bw -> SetParameters(1., Mphi, Gphi,R_BW);
   bw -> SetLineWidth(1);
   bw -> SetLineColor(kBlue);
   bw -> SetNpx(500);

   double norm = 1./bw->Integral(dL,dU,1e-8);
   printf("norm = %.7f\n",norm);
   bw -> SetParameter(0, norm );

   // test of stability
//    for ( int i = 0; i < 100; i++ ) { // mphi = Mphi +/- 1MeV
//       double mphi = Mphi + (-1. + 0.02*i)*1e-3;
//       for ( int j = 0; j < 100; j++ ) { // sigma from 1.1 to 1.3 MeV
//          double sigma = (1.1 + 0.02*j)*1e-3;
//          BreitWignerGaussN(mphi-4*sigma,mphi,Gphi,sigma);
//       }
//    }

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
   gPad -> SetLogy(true);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg -> SetHeader("Breit-Wigner (M_{#phi}, #Gamma_{#phi})","C");
   leg -> AddEntry(bw -> Clone(), "BW, Rff=3GeV^{-1}", "L");

   bw -> DrawCopy("L");
   c1 -> Update();

   bw -> Update();
   bw -> SetParameters(1., Mphi, Gphi, 6.); // Blatt-Weisskopf x2
   bw -> SetLineColor(kGreen);
   bw -> DrawCopy("LSAME");
   leg -> AddEntry(bw -> Clone(), "BW, Rff=6GeV^{-1}", "L");

   bw -> Update();
   bw -> SetParameters(1., Mphi, Gphi, 1.); // Blatt-Weisskopf ff
   bw -> SetLineColor(kRed);
   bw -> DrawCopy("LSAME");
   leg -> AddEntry(bw -> Clone(), "BW, Rff=1GeV^{-1}", "L");


//    bwgn -> DrawCopy("SAME");
//    leg -> AddEntry(bwgn -> Clone(), Form("#sigma=%.2f sl=%.2f",
//             bwgn->GetParameter(3)*1e3, bwgn->GetParameter(4) ), "L");

//    bwgn -> SetParameter(4, -1.9 );  // slope
//    bwgn -> SetLineColor(kGreen);
//    bwgn -> DrawCopy("SAME");
//    leg -> AddEntry(bwgn -> Clone(), Form("#sigma=%.2f sl=%.2f",
//             bwgn->GetParameter(3)*1e3, bwgn->GetParameter(4) ), "L");

   leg -> Draw();
   gPad -> RedrawAxis();
   c1 -> Update();
}

//-------------------------------------------------------------------------
void sig_fit(string fname, string title, string pdf="") {
//-------------------------------------------------------------------------
   bool is2009 = (fname.find("_09") != string::npos);

   // Get un-binned data
   TH1D* hist[1];
   vector<double> mkk = get_mkk_hist( fname, "mkk_phi", hist );
   TH1D* hst = hist[0];
   double norm = hst -> Integral();

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Sig(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Sig.Add(mkk[i]);
   }
   cout << " norm= " << norm << " n= " << n << endl;

   //-----------------------------------------------------------------------
   // Fit signal MC(phi,eta) by Breit-Wigner convoluted with Gauss
   // function MUST be normalized to 1 on the fit range
   double slope = -1.87;
   if ( is2009 ) {
      slope = -1.78;
   }
   auto Lbwg = [slope](const double* x,const double* p) -> double {
      return BreitWignerGaussN(x[0],p[0],p[1],p[2],slope);
   };
   vector<string> par_name { "M#phi", "G#phi", "#sigma" };
   vector<double> par_ini  {   Mphi,   Gphi,    1.e-3   };

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* bwg = new TF1("bwg", Lbwg, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 fbwg( *bwg );
   fitter.SetFunction( fbwg, false); // false == no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(Mphi-0.01,Mphi+0.01); // Mphi
//    fitter.Config().ParSettings(0).Fix();
//    fitter.Config().ParSettings(1).SetLimits(Gphi-1e-4,Gphi+1e-4); // Gphi
//    fitter.Config().ParSettings(1).SetValue(4.3e-3); // MC
   fitter.Config().ParSettings(1).Fix();
   fitter.Config().ParSettings(2).SetLimits(0.2e-3, 2.e-3);       // sigma

   fitter.LikelihoodFit( Sig, false ); // true == extended likelihood fit
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof = [Lbwg,Fpar](double x) -> double {
      return Lbwg(&x,Fpar.data());
   };
   rmath_fun< decltype(Lgof) > ftest(Lgof); // lambda -> Root1Dfunc
   // we use statistics five times larger than the data
//    size_t max_mkk = 14000; // 5*2800
   size_t max_mkk = 10000;
   if ( is2009 ) {
      max_mkk = 4700;
//       max_mkk = 1000;
   }
   mkk.resize(max_mkk);
   ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(
         mkk.size(),mkk.data(),ftest,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS = goftest -> KolmogorovSmirnovTest();
   cout << " pvalueKS= " << pvalueKS << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   double Wn = bW*n;
   auto Lfit = [Wn,Lbwg,Fpar](double* x,double* p) -> double {
      return  Wn * Lbwg(x,Fpar.data());
   };
   TF1* ffit = new TF1("ffit", Lfit, dL, dU, 0);
   ffit -> SetLineWidth(2);
   ffit -> SetLineColor(kRed);
   ffit -> SetNpx(500);

   //-----------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();
   gPad -> SetLogy();

   SetHstFace(hst);
   hst -> SetMinimum(1.);
   hst -> GetYaxis() -> SetTitleOffset(1.25);
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
   hst -> SetMarkerStyle(20);

   hst -> Draw("EP");
   ffit -> Draw("SAME");

   TLegend* leg = new TLegend(0.55,0.81,0.89,0.89);
   leg -> AddEntry(hst,title.c_str(),"LEP");
   leg -> AddEntry(bwg,"(BW#timesE) #otimes Gauss","L");
   leg -> Draw();

   TPaveText* pt = new TPaveText(0.55,0.63,0.89,0.80,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
//    pt -> AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt -> AddText( Form("#it{p-value(K-S) = %.3f}",pvalueKS) );
   pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt -> AddText( Form("#sigma= %s MeV", PE.Eform(2,".2f",1e3)) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void sig_fit_mc(string fname, string pdf="") {
//-------------------------------------------------------------------------
   TH1D* hst = get_mcMkk(fname);
   double bin_width = hst->GetBinWidth(1);

   // Fit signal MC(phi,eta) by Breit-Wigner
   auto Lbwg = [bin_width](const double* x,const double* p) -> double {
//       return bin_width*p[0]*BreitWigner(x[0],p[1],p[2]);
      return bin_width*p[0]*BreitWigner(x[0],p[1],p[2],p[3]);
   };
   TF1* bwg = new TF1("bwg", Lbwg, 0.98, 2., 4);
   bwg->SetParNames("Norm","Mphi","Gphi","r");
   bwg->SetParameters(1.e5, Mphi, Gphi, R_BW);

//    bwg->FixParameter(1, Mphi);
//    bwg->FixParameter(2, Gphi);
//    bwg->FixParameter(3, R_BW);
   bwg->SetParLimits(3, 0.5, 5.);

   bwg->SetLineWidth(2);
   bwg->SetLineColor(kRed);
   bwg->SetNpx(500);

//    gStyle->SetOptFit(111); // do not print fixed parameters!
   gStyle->SetOptFit(112); // do not print fixed parameters!

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();
   gStyle->SetFitFormat(".7g");
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);

   SetHstFace(hst);
   hst->SetMinimum(1.);
   hst->GetYaxis()->SetTitleOffset(1.25);
   hst->SetLineWidth(3);
   hst->SetLineColor(kBlack);

   hst->Draw("EP");
   double maxMkk=1.0835;
   hst->Fit("bwg","E","",dL,maxMkk);

   TLine* lB = new TLine;
   lB -> SetLineColor(kBlue+1);
   lB -> SetLineWidth(2);
   lB -> SetLineStyle(kDashed);
   lB -> DrawLine(maxMkk,1,maxMkk,hst->GetBinContent(103));

   TF1* bwg1 = (TF1*) bwg -> Clone();
   bwg1 -> SetLineStyle(kDashed);
   bwg1 -> DrawF1(maxMkk,1.125,"SAME");

   double genEv = hst -> GetEntries();
   double corr = bwg -> GetParameter(0) / genEv;
   double dcorr = bwg -> GetParError(0) / genEv;
   cout << " genEv= " << genEv << " corr= " << corr << endl;

   double intg108 = bwg -> Integral(dL,1.08,1e-8) / bin_width;
   double intg12 = bwg -> Integral(dL,1.2,1e-8) / bin_width;
   double norm108 = bwg -> GetParameter(0) / intg108;
   double norm12 = bwg -> GetParameter(0) / intg12;
   cout << " norm108= " << norm108 << endl;
   cout << " norm12= " << norm12 << endl;

   gStyle->SetLegendTextSize(0.03);
   TLegend* leg = new TLegend(0.22,0.20,0.62,0.40);
   if ( genEv > 8e5 ) {
      leg->AddEntry(hst,"MC truth 2012","LEP");
   } else {
      leg->AddEntry(hst,"MC truth 2009","LEP");
   }
   leg->AddEntry(bwg,"Breit-Wigner","L");
   leg->AddEntry(lB,"Cut-off","L");
   string txt(Form(" Cut-off corr.: %.3f #pm %.3f",corr,dcorr));
   leg -> SetHeader( txt.c_str(), "C" );
   leg->Draw();

   c1->Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
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

   static const double L = 2*Mk;
   static const double U = Mjpsi - Meta;

   double x = (U-m)/(U-L); // transformation: (L,U) -> (1,0)
   return Argus(x,A)/(U-L);
}

//----------------------------------------------------------------------
double RevArgusN(double m, double A, double slope) {
//----------------------------------------------------------------------
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
      ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);
      norm = gsl_int.Integral(Lint,p,dL,dU);
      if ( DEBUG ) {
         htt[3]->Fill(norm);
         if ( gsl_int.Error() > 1e-6 ) {
            htt[4]->Fill(log10(gsl_int.Error()));
         } else {
            htt[4]->Fill(-7.);
         }
      }
      cacheN = norm;
      cacheA = A;
      cacheS = slope;
   }

   return RevArgus(m,A) * (1+slope*(m-1.02)) / norm;
}

//-------------------------------------------------------------------------
void test_Agrus() {
//-------------------------------------------------------------------------
   auto Largus = [](const double* x,const double* p) -> double {
//       return p[0]*Argus(x[0],p[1]);
      return p[0]*RevArgus(x[0],p[1]);
   };

//    TF1* fbkg = new TF1("fbkg", Largus, 0., 1., 2);
   TF1* fbkg = new TF1("fbkg", Largus, 0.98, 2., 2);
   fbkg -> SetParNames("N","A");
   fbkg -> SetParameters(1., 2.);
   fbkg -> SetLineWidth(2);
   fbkg -> SetLineColor(kBlue);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1 -> cd();
   gPad -> SetGrid();
   fbkg -> DrawCopy();

   // test normalization:
//    printf(" norm= %.15f\n",fbkg->Integral(0.,1.));
   printf(" norm= %.15f\n",fbkg -> Integral(2*Mk,Mjpsi-Meta));

   fbkg -> SetParameters(1., 1.);
   fbkg -> SetLineColor(kRed);
   fbkg -> DrawCopy("SAME");
}

//-------------------------------------------------------------------------
void test_RevAgrus() {
//-------------------------------------------------------------------------
   auto Lrargus = [](const double* x,const double* p) -> double {
      return p[0]*RevArgusN(x[0],p[1],p[2]);
   };

   TF1* fbkg = new TF1("fbkg", Lrargus, 0.98, dU, 3);
   fbkg -> SetParNames("N","A","Sl");
   fbkg -> SetParameters(1., 0.98, 0.);
   fbkg -> SetLineWidth(2);
   fbkg -> SetLineColor(kBlue);

   TLegend* leg = new TLegend(0.35,0.20,0.75,0.40);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad -> SetGrid();
   fbkg -> DrawCopy();
   leg -> AddEntry(fbkg -> Clone(), Form("RevArgus A=%.2f Sl=%.2f",
            fbkg -> GetParameter(1),fbkg -> GetParameter(2)),"L");

   // test normalization:
   printf(" RevArgusN norm= %.15f\n",fbkg -> Integral( dL,dU ) );

   fbkg -> SetParameter(1, 0. );  // A
   fbkg -> SetParameter(2, -1.9 ); // Slope
   fbkg -> SetLineColor(kRed);
   fbkg -> SetLineStyle(kDashed);
   fbkg -> DrawCopy("SAME");
   leg -> AddEntry(fbkg -> Clone(), Form("RevArgus A=%.2f Sl=%.2f",
            fbkg -> GetParameter(1),fbkg -> GetParameter(2)),"L");
   printf(" RevArgusN norm2= %.15f\n",fbkg -> Integral( dL,dU ) );

   fbkg -> SetParameter(1, 0.98 );  // A
   fbkg -> SetLineColor(kGreen);
   fbkg -> DrawCopy("SAME");
   leg -> AddEntry(fbkg -> Clone(), Form("RevArgus A=%.2f Sl=%.2f",
            fbkg -> GetParameter(1),fbkg -> GetParameter(2)),"L");

   leg -> Draw();
}

//-------------------------------------------------------------------------
void bkg_fit(string fname, string title, string pdf="") {
//-------------------------------------------------------------------------
   bool is2009 = (fname.find("_09") != string::npos);

   // Get un-binned data
   TH1D* hist[1];
   vector<double> mkk = get_mkk_hist( fname, "mkk_kketa", hist );
   TH1D* hst = hist[0];
   double norm = hst -> Integral();

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Bkg(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Bkg.Add(mkk[i]);
   }
   cout << " norm= " << norm << " n= " << n << endl;

   //-----------------------------------------------------------------------
   // Fit KKeta background by reversed Argus function
   // function MUST be normalized to 1 on the fit range
   double slope = -1.87;
   if ( is2009 ) {
      slope = -1.78;
   }
   auto Lan = [slope](const double* x,const double* p) -> double {
      return RevArgusN(x[0],p[0],slope);
   };
   vector<string> par_name { "A" };
   vector<double> par_ini  {  0. };

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* bkg = new TF1("bkg", Lan, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 fbkg( *bkg );
   fitter.SetFunction( fbkg, false); // false == no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(-2.5, 2.5);

   fitter.LikelihoodFit( Bkg, false ); // true == extended likelihood fit
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof = [Lan,Fpar](double x) -> double {
      return Lan(&x,Fpar.data());
   };
   if ( mkk.size() > 4000 ) {
      mkk.resize(4000);
   }
   rmath_fun< decltype(Lgof) > ftest(Lgof); // lambda -> Root1Dfunc
   ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(
         mkk.size(),mkk.data(),ftest,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS = goftest -> KolmogorovSmirnovTest();
   cout << " pvalueKS= " << pvalueKS << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   double Wn = bW*n;
   auto Lfit = [Wn,Lan,Fpar](double* x,double* p) -> double {
      return Wn * Lan(x,Fpar.data());
   };
   TF1* ffit = new TF1("ffit", Lfit, dL-1e-2, dU, 0);
   ffit -> SetLineWidth(2);
   ffit -> SetLineColor(kBlue);
   ffit -> SetNpx(200);

   //-----------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   SetHstFace(hst);
   hst -> SetMaximum( 1.3*hst->GetMaximum() );
   hst -> GetYaxis() -> SetMaxDigits(3);
   hst -> GetYaxis() -> SetTitleOffset(1.3);
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
   hst -> SetMarkerStyle(20);

   hst -> Draw("EP");
   ffit -> Draw("SAME");

   TLegend* leg = new TLegend(0.11,0.77,0.50,0.89);
   leg -> AddEntry(hst,title.c_str(),"LEP");
   leg -> AddEntry(ffit,"Argus function","L");
   leg -> Draw();

   TPaveText* pt = new TPaveText(0.55,0.78,0.89,0.89,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
//    pt -> AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt -> AddText( Form("#it{p-value(K-S) = %.3f}",pvalueKS) );
   pt -> AddText( Form("a= %s",PE.Eform(0,".2f")) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void bkg_fit2(string fname09, string fname12, string pdf="") {
//-------------------------------------------------------------------------
   double norm09 = 107.0/1.2;
   double norm12 = 341.1/3.6;

   // Get un-binned data
   TH1D* hist[2];
   vector<double> mkk09 = get_mkk_hist( fname09, "mkk_kketa09", &hist[0] );
   vector<double> mkk12 = get_mkk_hist( fname12, "mkk_kketa12", &hist[1] );

   TH1D* hst = hist[0];
   hst -> Add(hist[0], hist[1], norm09, norm12);
   double norm = hst -> Integral();

   int n09 = mkk09.size();
   int n12 = mkk12.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Bkg(dr, n09+n12, 1, true); // weighted !
   for ( int i = 0; i < n09; ++i ) {
      Bkg.Add(mkk09[i], norm09);
   }
   for ( int i = 0; i < n12; ++i ) {
      Bkg.Add(mkk12[i], norm12);
   }

   //-----------------------------------------------------------------------
   // Fit KKeta background by reversed Argus function
   // function MUST be normalized to 1 on the fit range
   auto Lan = [](const double* x,const double* p) -> double {
      double slope = -1.9;
      return RevArgusN(x[0],p[0],slope);
   };
   vector<string> par_name { "A" };
   vector<double> par_ini  {  0. };

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* bkg = new TF1("bkg", Lan, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 fbkg( *bkg );
   fitter.SetFunction( fbkg, false); // false == no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(-2.5, 2.5);

   fitter.LikelihoodFit( Bkg, false ); // true == extended likelihood fit
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof = [Lan,Fpar](double x) -> double {
      return Lan(&x,Fpar.data());
   };
   rmath_fun< decltype(Lgof) > ftest(Lgof); // lambda -> Root1Dfunc
   vector<double> mkk;
   mkk.reserve(n09+n12);
   mkk.insert( mkk.end(), mkk09.begin(), mkk09.end() );
   mkk.insert( mkk.end(), mkk12.begin(), mkk12.end() );
   ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(
         mkk.size(),mkk.data(),ftest,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS = goftest -> KolmogorovSmirnovTest();
   cout << " pvalueKS= " << pvalueKS << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   double Wn = bW*norm;
   auto Lfit = [Wn,Lan,Fpar](double* x,double* p) -> double {
      return Wn * Lan(x,Fpar.data());
   };
   TF1* ffit = new TF1("ffit", Lfit, dL-1e-2, dU, 0);
   ffit -> SetLineWidth(2);
   ffit -> SetLineColor(kBlue);
   ffit -> SetNpx(200);

   //-----------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   SetHstFace(hst);
   hst -> SetMaximum( 1.3*hst->GetMaximum() );
   hst -> GetYaxis() -> SetMaxDigits(3);
   hst -> GetYaxis() -> SetTitleOffset(1.2);
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
   hst -> SetMarkerStyle(20);

   hst -> Draw("EP");
   ffit -> Draw("SAME");

   TLegend* leg = new TLegend(0.11,0.71,0.47,0.89);
   leg -> SetHeader("2009 & 2012", "C");
   leg -> AddEntry(hst,"MC KK#eta","LEP");
   leg -> AddEntry(ffit,"Argus function","L");
   leg -> Draw();

   TPaveText* pt = new TPaveText(0.53,0.77,0.89,0.89,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
//    pt -> AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt -> AddText( Form("#it{p-value(K-S) = %.3f}",pvalueKS) );
   pt -> AddText( Form("a= %s", PE.Eform(0,".2f")) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

// {{{1 data_fit() no interference
//-------------------------------------------------------------------------
void data_fit(string fname,string title,int date,int Mphix,string pdf="") {
//-------------------------------------------------------------------------
   bool is2009 = (fname.find("_09") != string::npos);

   // Get un-binned data
   TH1D* hist[1];
   vector<double> mkk = get_mkk_hist( fname, "mkk_data_cp", hist );
   TH1D* hst = hist[0];
   double norm = hst -> Integral();

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Dat(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Dat.Add(mkk[i]);
   }
   cout << " norm= " << norm << " n= " << n << endl;

   //-----------------------------------------------------------------------
   // Fit data
   double slope = -1.87;
   if ( is2009 ) {
      slope = -1.78;
   }
   // function represents the total number of events (extended LH)
   auto Lsum = [slope](const double* x,const double* p) -> double {
      return p[4] * BreitWignerGaussN(x[0],p[0],p[1],p[2],slope) +
             p[5] * RevArgusN(x[0],fabs(p[3]),slope);
   };
   vector<string> par_name { "M#phi", "#Gamma#phi", "#sigma",
                             "Ar", "N#phi", "Nbkg"                     };
   vector<double> par_ini  { Mphi, Gphi, 1.16e-3,
                             0., 0.94*norm, 0.06*norm };

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* fsum = new TF1("fsum", Lsum, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 WFsum( *fsum );
   fitter.SetFunction( WFsum, false); // false == no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   if ( Mphix == 1 ) {
      fitter.Config().ParSettings(0).SetValue(1.01953); // like MC
      fitter.Config().ParSettings(0).Fix();
   } else {
      fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   }
   fitter.Config().ParSettings(1).Fix(); // Gphi
   if ( is2009 ) {
      fitter.Config().ParSettings(2).SetValue(1.4e-3); // sigma09
   }
   fitter.Config().ParSettings(2).SetLimits(0.5e-3, 2.e-3); // sigma
   fitter.Config().ParSettings(3).SetLimits(-5.,5.);        // Argus
   fitter.Config().ParSettings(3).Fix();                    // Ar
//    fitter.Config().ParSettings(4).SetLimits(0.,1.5*norm);   // Nphi
//    fitter.Config().ParSettings(5).SetLimits(0.,norm);       // Nbkg

   fitter.LikelihoodFit( Dat, true ); // true == extended likelihood fit
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof = [Lsum,Fpar](double x) -> double {
      return Lsum(&x,Fpar.data());
   };
   rmath_fun< decltype(Lgof) > ftest(Lgof); // adapter lambda -> Root1Dfunc
   ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(
         mkk.size(),mkk.data(),ftest,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS = goftest -> KolmogorovSmirnovTest();
   cout << " pvalueKS= " << pvalueKS << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   auto Lfit = [Lsum,Fpar](double* x,double* p) -> double {
      return bW * Lsum(x,Fpar.data());
   };
   TF1* ffit = new TF1("ffit", Lfit, dL, dU, 0);
   ffit -> SetLineWidth(2);
   ffit -> SetLineColor(kRed);
   ffit -> SetNpx(500);

   auto Lbkg = [Fpar,slope](double* x,double* p) -> double {
      return bW * Fpar[5] * RevArgusN(x[0],fabs(Fpar[3]),slope);
   };
   TF1* fbkg = new TF1("fbkg", Lbkg, dL, dU, 0);
   fbkg -> SetLineWidth(2);
   fbkg -> SetLineColor(kBlue);
   fbkg -> SetLineStyle(kDashed);

   //-----------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   SetHstFace(hst);
   hst -> GetYaxis() -> SetTitleOffset(1.3);
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
   hst -> SetMarkerStyle(20);

   hst -> Draw("EP");
   ffit -> Draw("SAME");
   fbkg -> Draw("SAME");

   TLegend* leg = new TLegend(0.55,0.76,0.89,0.89);
   leg -> AddEntry(hst,title.c_str(),"LEP");
   leg -> AddEntry(ffit,"Fit result","L");
   leg -> AddEntry(fbkg,"non-#phi KK#eta","L");
   leg -> Draw();

   TPaveText* pt = new TPaveText(0.55,0.45,0.89,0.75,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
//    pt -> AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt -> AddText( Form("#it{p-value(K-S) = %.3f}",pvalueKS) );
   pt -> AddText( Form("N_{#phi}= %s",PE.Eform(4,".1f")) );
   pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt -> AddText( Form("#sigma= %s MeV", PE.Eform(2,".2f",1e3)) );
   pt -> AddText( Form("N_{non-#phi}= %s",PE.Eform(5,".1f")) );
   pt -> AddText( Form("a= %s",PE.Eform(3,".2f")) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
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
   if ( m < dL ) { // phase space becomes zero
      vector<double> ret(4,0.);
      return ret;
   }

   double m2    = SQ(m);
   double mphi2 = SQ(mphi);
   double gphi2 = SQ(gphi);

   double R     = R_BW; // Blatt-Weisskopf ff

   double gam = mphi*sqrt(mphi2+gphi2);
   double kap = 2*M_SQRT2*mphi*gphi*gam / (M_PI*sqrt(mphi2+gam));

   double p_phi0 = sqrt(SQ(Mjpsi2-mphi2-Meta2)-4*mphi2*Meta2) / (2*Mjpsi);
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
//    double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k; // == m*G(m)
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

//----------------------------------------------------------------------
double IntfrBWARG( double m,
                   const double p[],    // {mphi,gphi,sigma,A,F,Ang,sl}
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
   static const double Ngauss = 1./sqrt(2*M_PI);

   // paranoic check
//    if ( idx > 3 ) {
//       cout << " FATAL ERROR: " << __func__ << " idx= " << idx << endl;
//       exit;
//    }

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
   ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);
   double result = gsl_int.Integral(Fint,-5.,+5.);

   return Ngauss * result;
}

//----------------------------------------------------------------------
double IntfrBWARGN( double m,
                    const double p[], // {mphi,gphi,sigma,A,F,Ang,sl}
                    int idx = 0       // what to return
                  ) {
//----------------------------------------------------------------------
// Numerical normalisation to one for range [dL,dU]

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
         return IntfrBWARG(x,p,0);  // MUST be zero!
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-6, 1.e-5, 1000);
      norm = gsl_int.Integral(Lint,(void *)p,dL,dU);

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

//-------------------------------------------------------------------------
void test_Intfr() {
//-------------------------------------------------------------------------
   double mphi = Mphi, gphi = Gphi;
   auto Lintfr = [mphi,gphi](const double* x, const double* p) -> double {
      double pp[] { mphi, gphi, p[0],p[1],p[2],p[3],p[4] };
      return IntfrBWARGN( x[0], pp, int(p[5]) );
   };

   int npar = 6;
   TF1* fun = new TF1("fun", Lintfr, dL, dU, npar);
   fun -> SetParNames("Sigma", "A", "F", "Ang",  "Sl", "T");
   fun -> SetParameters(1.2e-3, 0., 0.8,   0.7,  -1.99,  0.);
   fun -> SetLineWidth(2);
   fun -> SetNpx(500);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   TH1D* tmp = new TH1D("tmp","",500,bL,dU);
   tmp -> SetAxisRange(-20.,130.,"Y");
   tmp -> Draw();

   fun -> SetParameter(npar-1, 0); // draw all
   fun -> SetLineColor(kRed);
   fun -> DrawCopy("SAME");
   printf(" Norm= %.6f\n", fun -> Integral(dL,dU,1e-4) );

   fun -> SetParameter(npar-1, 2); // draw Argus
   fun -> SetLineColor(kBlue);
   fun -> SetLineStyle(kDashed);
   fun -> DrawCopy("SAME");

   fun -> SetParameter(npar-1, 1); // draw BW
   fun -> SetLineColor(kGreen+2);
   fun -> DrawCopy("SAME");

   fun -> SetParameter(npar-1, 3); // draw interference
   fun -> SetLineColor(kMagenta);
   fun -> SetLineStyle(kSolid);
   fun -> DrawCopy("SAME");

   c1->Update();
}


//----------------------------------------------------------------------
ValErr CalcIntEr( const ROOT::Fit::FitResult& res,
                  double sl, unsigned int idx ) {
//----------------------------------------------------------------------
// Numerical calculation of integrals and their errors
   // desired errors:
   constexpr double eps_abs = 1.e-6;
   constexpr double eps_rel = 1.e-6;
   constexpr double epsilon  = 1.e-3; // max numerical error

   // get parameters and covariation matrix of error
   int Npar = res.NPar();
   vector<double> Fpar = res.Parameters();

   if ( Npar != 6 && Npar != 7 ) {
      cerr << " FATAL ERROR in " << __func__ << " Npar= " << Npar << endl;
      exit(1);
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
   auto Lfc = [sl,idx,Npar](const double* x,const double* p) -> double {
      double m = x[0];
      double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl }; // use real slope
      double normI = IntfrBWARGN( m,pp,100 );
      pp[6] = 0.; // set slope to zero
      double intfr = IntfrBWARG( m,pp,idx ) / normI;
      if ( Npar == 6 ) {
         return intfr;
      }
      if ( idx != 0 ) {
         return (1-p[6])*intfr;
      }
      double bg = RevArgusN(m,0.,sl) / (1+sl*(m-1.02)); // forget sl
      return (1-p[6])*intfr + p[6]*bg;
   };
   TF1 FC("FC", Lfc, dL, dU, Npar);
   FC.SetParameters( Fpar.data() );

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
//    printf("\ncovMatrix : ");
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
void data_Intfr_fit(string fname, string title, string pdf="") {
//-------------------------------------------------------------------------
   bool is2009 = (fname.find("_09") != string::npos);

   // Get un-binned data
   TH1D* hist[1];
   vector<double> mkk = get_mkk_hist( fname, "mkk_data_cp", hist );
   TH1D* hst = hist[0];

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Dat(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Dat.Add(mkk[i]);
   }

   //-----------------------------------------------------------------------
   // Fit data
   // efficiency parameters
   double sl = -1.87;
   if ( is2009 ) {
      sl = -1.78;
   }
   // function MUST be normalized to 1 on the fit range
   auto Lfit = [sl](const double* x,const double* p) -> double {
      const double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl };
      return IntfrBWARGN( x[0], pp, 0 );
   };
   vector<string> par_name { "M#phi", "G#phi", "#sigma", "A",
                             "F", "#vartheta" };

   vector<double> par_ini {Mphi,Gphi,1.5e-3,0.,1.0,1.0}; // 09 neg
//    vector<double> par_ini {Mphi,Gphi,1.5e-3,0.,1.0,-1.1}; // 09 pos
//    vector<double> par_ini {Mphi,Gphi,1.2e-3,0.,1.0,0.7}; // 12 neg
//    vector<double> par_ini {Mphi,Gphi,1.2e-3,0.,1.0,-0.7}; // 12 pos

   bool neg_pos = par_ini.back() > 0; // true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 WFit( *Ffit );
   fitter.SetFunction( WFit, false); // false == no parameter derivatives

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
//    fitter.Config().ParSettings(5).SetValue(-95*M_PI/180);
//    fitter.Config().ParSettings(5).Fix();

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
   cout << " n= " << n << endl;
   for ( int idx = 0; idx <= 3; ++idx ) {
      ValErr Integ = CalcIntEr( res, sl, idx );
//       printf("Int[%i] = %s\n",idx,Integ.prt(".4f")); // DEBUG
      double Num = n * Integ.val;
      double Err = fabs(Num) * sqrt( 1./n + SQ(Integ.err/Integ.val) );
      Nos.push_back( ValErr {Num,Err} );
      printf("%s = %s\n",names[idx].c_str(),Nos[idx].prt(".1f"));
   }
   double Nphi = Nos[1].val, err_Nphi = Nos[1].err;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof = [Lfit,Fpar](double x) -> double {
      return Lfit( &x,Fpar.data() );
   };
   rmath_fun< decltype(Lgof) > ftest(Lgof); // adapter lambda -> Root1Dfunc
   ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(
         mkk.size(),mkk.data(),ftest,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS = goftest -> KolmogorovSmirnovTest();
   cout << " pvalueKS= " << pvalueKS << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   auto Ldr = [n,Fpar,sl](double* x,double* p) -> double {
      const double pp[]
        { Fpar[0],Fpar[1],Fpar[2],Fpar[3],Fpar[4],Fpar[5],sl };
      return bW*n*IntfrBWARGN( x[0], pp, int(p[0]) );
   };
   TF1* ffit = new TF1("ffit", Ldr, dL, dU, 1);
   ffit -> SetNpx(500);

   // find max/min to draw
   ffit -> SetParameter(0, 1); // BW
   int maxbin = hst -> GetMaximumBin();
   double hmax = hst -> GetBinContent(maxbin) + hst -> GetBinError(maxbin);
   double fmax = ffit -> GetMaximum( 1.01, 1.03, 1e-6);
   hmax = floor(1.1 * max( hmax, fmax ));
   if ( hmax > 0 ) {
      hst -> SetMaximum(hmax);
   }

   ffit -> SetParameter(0, 3); // interference
   double hmin = ffit -> GetMinimum( 1.01, 1.03, 1e-6);
   hmin = floor(1.15 * min( 0., hmin ));
   if ( hmin < 0 ) {
      hst -> SetMinimum(hmin);
   }

   //-----------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   SetHstFace(hst);

   hst -> GetYaxis() -> SetTitleOffset(1.3);
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
   hst -> SetMarkerStyle(20);

   hst -> Draw("EP");

   TLegend* leg = new TLegend(0.52,0.67,0.89,0.89);
   leg -> AddEntry(hst,title.c_str(),"LEP");

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
   leg -> AddEntry( ffit -> Clone(), "Breit-Wigner", "L");

   ffit -> SetParameter(0, 2); // Argus
   ffit -> SetLineColor(kBlue);
   ffit -> DrawCopy("SAME");
   leg -> AddEntry( ffit -> Clone(), "Argus", "L");

   ffit -> SetParameter(0, 3); // interference
   ffit -> SetLineColor(kMagenta+1);
   ffit -> DrawCopy("SAME");
   leg -> AddEntry( ffit -> Clone(), "Interference", "L");

   leg -> Draw();

   TPaveText* pt = new TPaveText(0.52,0.33,0.89,0.66,"NDC");
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
      pdf += (neg_pos ? "_n.pdf" : "_p.pdf");
      c1 -> Print(pdf.c_str());
   }
}

// combined with side-band
//----------------------------------------------------------------------
struct myFCN {
   double sl;           // efficiency parameter
   vector<double> mkk;  // data central part
   vector<double> sb;   // data side band

   // minimization function
   // the signature of this operator() MUST be exactly this:
   // THE EXTENDED MAXIMUM LIKELIHOOD ESTIMATOR
   // Kai-Feng Chen "FittingMinuitRooFit" p.35
   //------------------------------------------------------------------
   double operator() (const double* p) {
   //------------------------------------------------------------------
      double Nkk = p[6];
      double Nsb = p[7];
      double res = 2*(Nkk+Nsb);

      const double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl };
      int n_mkk = mkk.size();
      for ( int i = 0; i < n_mkk; ++i ) {
         double m = mkk[i];
         double L = Nkk*IntfrBWARGN( m,pp,0 ) + Nsb/(dU-bL);
         if (L > 0.) {
            res -= 2*log(L);
         } else {
            return DBL_MAX;
         }
      }

//       int n_sb = sb.size();
//       res += 2*Nsb - n_sb*2*log(Nsb/(dU-bL)); // fit SB by constant

      return res;
   }
};
//----------------------------------------------------------------------

//----------------------------------------------------------------------
ValErr CalcIntErSb( const ROOT::Fit::FitResult& res,
                    double sl, unsigned int idx ) {
//----------------------------------------------------------------------
// Numerical calculation of integrals and their errors
   // desired errors:
   constexpr double eps_abs = 1.e-6;
   constexpr double eps_rel = 1.e-6;
   constexpr double epsilon  = 1.e-3; // max numerical error

   // get parameters and covariation matrix of error
   int Npar = res.NPar();
   vector<double> Fpar = res.Parameters();

   if ( Npar != 7 ) {
      cerr << " FATAL ERROR in " << __func__ << " Npar= " << Npar << endl;
      exit(1);
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
   auto Lfc = [sl,idx](const double* x,const double* p) -> double {
      double m = x[0];
      double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl }; // use real slope
      double normI = IntfrBWARGN( m,pp,100 );
      pp[6] = 0.; // set slope to zero
      double intfr = IntfrBWARG( m,pp,idx ) / normI;
      // remove background from consideration
      return (1-p[6])*intfr;
   };
   TF1 FC("FC", Lfc, dL, dU, Npar);
   FC.SetParameters( Fpar.data() );

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
//    printf("\ncovMatrix : ");
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
void data_Intfr_sb(string fname, string title, string pdf="") {
//-------------------------------------------------------------------------
   bool is2009 = (fname.find("_09") != string::npos);

   myFCN my_fcn; // class for 'FitFCN'

   // efficiency parameters
   double sl = -1.87;
   if ( is2009 ) {
      sl = -1.78;
   }
   my_fcn.sl = sl;

   // Get un-binned data
   TH1D* hist[2];
   my_fcn.mkk = get_mkk_hist( fname, "mkk_data_cp", hist );
   const vector<double>& mkk = my_fcn.mkk;
   TH1D* hst = hist[0];

   // side-band:
   my_fcn.sb = get_mkk_hist( fname, "mkk_sb", &hist[1], 1 );
   const vector<double>& sb = my_fcn.sb;
   TH1D* hsb = hist[1];

   int n = my_fcn.mkk.size();
   int nsb = my_fcn.sb.size();

   //-----------------------------------------------------------------------
   // Fit data
   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   vector<string> par_name { "M#phi", "G#phi", "#sigma", "A",
                             "F", "#vartheta","NKK","Nbg" };

   vector<double> par_ini {Mphi,Gphi,1.5e-3,0.,0.75,0.75,
                           double(n-nsb),double(nsb)};  // 09 neg
//    vector<double> par_ini {Mphi,Gphi,1.5e-3,0.,0.75,-0.8,
//                            double(n-nsb),double(nsb)};  // 09 pos
//    vector<double> par_ini {Mphi,Gphi,1.2e-3,0.,0.8,0.8,
//                            double(n-nsb),double(nsb)};  // 12 neg
//    vector<double> par_ini {Mphi,Gphi,1.2e-3,0.,1.0,-0.7,
//                            double(n-nsb),double(nsb)};  // 12 pos

   bool neg_pos = par_ini[5] > 0; // true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters

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
//    fitter.Config().ParSettings(5).SetValue(-95*M_PI/180);
//    fitter.Config().ParSettings(5).Fix();
//    fitter.Config().ParSettings(6).SetLimits(0., 1.5*n);        // NKK
   fitter.Config().ParSettings(7).SetLimits(0., 2.*nsb);       // Nbg

   // == Fit
   int Ndat = n+nsb;
   fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // false == likelihood fit
   fitter.CalculateMinosErrors();

   const ROOT::Fit::FitResult& res = fitter.Result();
   res.Print(cout);
//    res.PrintCovMatrix(cout); // print error matrix and correlations

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   cout << " n= " << n << " nsb= " << nsb << endl;
//    vector<string> names { "N(KK)", "Nphi", "Nnonphi", "Nifr"  };
//    vector<ValErr> Nos;
//    Nos.reserve(4);
//    for ( int idx = 0; idx <= 3; ++idx ) {
//       ValErr Integ = CalcIntErSb( res, sl, idx );
//       printf("Int[%i] = %s\n",idx,Integ.prt(".4f")); // DEBUG
//       double Num = n * Integ.val;
//       double Err = fabs(Num) * sqrt( 1./n + SQ(Integ.err/Integ.val) );
//       Nos.push_back( ValErr {Num,Err} );
//       printf("%s = %s\n",names[idx].c_str(),Nos[idx].prt(".1f"));
//    }
//    double Nphi = Nos[1].val, err_Nphi = Nos[1].err;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto LgofCR = [Fpar,sl](double x) -> double {
      double pp[] {Fpar[0],Fpar[1],Fpar[2],Fpar[3],Fpar[4],Fpar[5],sl};
      return Fpar[6]*IntfrBWARGN( x,pp,0 ) + Fpar[7]/(dU-bL);
   };
   rmath_fun< decltype(LgofCR) > ftestCR(LgofCR);
   ROOT::Math::GoFTest* goftestCR = new ROOT::Math::GoFTest(
         mkk.size(),mkk.data(),ftestCR,ROOT::Math::GoFTest::kPDF,bL,dU);
   double pvalueKScr = goftestCR -> KolmogorovSmirnovTest();
   cout << " pvalueKScr= " << pvalueKScr << endl;

   auto LgofSB = [Fpar](double x) -> double {
      return Fpar[7]/(dU-bL);
   };
   rmath_fun< decltype(LgofSB) > ftestSB(LgofSB);
   ROOT::Math::GoFTest* goftestSB = new ROOT::Math::GoFTest(
         sb.size(),sb.data(),ftestSB,ROOT::Math::GoFTest::kPDF,bL,dU);
   double pvalueKSsb = goftestSB -> KolmogorovSmirnovTest();
   cout << " pvalueKSsb= " << pvalueKSsb << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   auto Ldr = [Fpar,sl](double* x,double* p) -> double {
      double pp[] {Fpar[0],Fpar[1],Fpar[2],Fpar[3],Fpar[4],Fpar[5],sl};
      return bW * Fpar[6]*IntfrBWARGN(x[0],pp,int(p[0]));
   };
   TF1* ffit = new TF1("ffit", Ldr, bL, dU, 1);
   ffit -> SetNpx(500);

   auto Ldr_bg = [Fpar](double* x,double* p) -> double {
      return bW * Fpar[7]/(dU-dL);
   };
   TF1* f_bg = new TF1("f_bg", Ldr_bg, bL, dU, 0);

   // find max/min to draw
   ffit -> SetParameter(0, 1); // BW
   int maxbin = hst -> GetMaximumBin();
   double hmax = hst -> GetBinContent(maxbin) + hst -> GetBinError(maxbin);
   double fmax = ffit -> GetMaximum( 1.01, 1.03, 1e-6);
   hmax = floor(1.1 * max( hmax, fmax ));
   if ( hmax > 0 ) {
      hst -> SetMaximum(hmax);
   }

   ffit -> SetParameter(0, 3); // interference
   double hmin = ffit -> GetMinimum( 1.01, 1.03, 1e-6);
   hmin = floor(1.15 * min( 0., hmin ));
   if ( hmin < 0 ) {
      hst -> SetMinimum(hmin);
   }

   //-----------------------------------------------------------------------
   // Draw results
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000);
   c1 -> Divide(1,2);

   c1 -> cd(1);
   gPad -> SetGrid();
   SetHstFace(hsb);

   hsb -> SetLineWidth(2);
   hsb -> SetLineColor(kBlack);
   hsb -> SetMarkerStyle(20);

   hsb -> Draw("EP");

   f_bg -> SetParameter(0, double(n));
   f_bg -> SetLineWidth(2);
   f_bg -> SetLineColor(kBlue+2);
   f_bg -> DrawCopy("SAME");

   c1 -> cd(2);
   gPad -> SetGrid();
   SetHstFace(hst);

   hst -> GetYaxis() -> SetTitleOffset(1.3);
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
   hst -> SetMarkerStyle(20);

   hst -> Draw("EP");

   TLegend* leg = new TLegend(0.51,0.45,0.89,0.89);
   leg -> SetHeader("side-band(top), central(bottom)}","C");
   leg -> AddEntry(hst,title.c_str(),"LEP");

   ffit -> SetParameter(0, 0); // SUM
   ffit -> SetLineWidth(2);
   ffit -> SetLineColor(kRed);
   ffit -> DrawCopy("SAME");
   leg -> AddEntry( ffit -> Clone(), "Result of fit", "L");

   ffit -> SetParameter(0, 1); // BW
   ffit -> SetLineWidth(1);
   ffit -> SetLineStyle(kDashed);
   ffit -> SetLineColor(kGreen+2);
   ffit -> DrawCopy("SAME");
   leg -> AddEntry( ffit -> Clone(), "Breit-Wigner #phi#eta", "L");

   ffit -> SetParameter(0, 2); // Argus
   ffit -> SetLineColor(kBlue);
   ffit -> DrawCopy("SAME");
   leg -> AddEntry( ffit -> Clone(), "Non-#phi KK#eta", "L");

   ffit -> SetParameter(0, 3); // interference
   ffit -> SetLineColor(kMagenta+1);
   ffit -> DrawCopy("SAME");
   leg -> AddEntry( ffit -> Clone(), "Interference", "L");

   leg -> AddEntry( f_bg -> Clone(), "Background", "L" );

   TPaveText* pt = new TPaveText(0.51,0.34,0.89,0.99,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
//    pt -> AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt -> AddText( Form("#it{p-value = %.3f}",pvalueKScr) );
   pt -> AddText( Form("#it{p-val(SB) = %.3f}",pvalueKSsb) );
//    pt -> AddText( Form("N_{#phi} = %.1f #pm %.1f",Nphi,err_Nphi) );
   pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt -> AddText( Form("#sigma= %s MeV", PE.Eform(2,".2f",1e3)) );
   pt -> AddText( Form("a= %s",PE.Eform(3,".1f")) );
   pt -> AddText( Form("F= %s",PE.Eform(4,".2f")) );
   pt -> AddText( Form("#vartheta= %s",PE.Eform(5,".2f")) );
   pt -> AddText( Form("NKK= %s", PE.Eform(6,".1f")) );
   pt -> AddText( Form("Nbg= %s", PE.Eform(7,".1f")) );
   pt -> Draw();

   gPad -> RedrawAxis();

   c1 -> cd(1);
   leg -> Draw();
   gPad -> RedrawAxis();

   c1 -> Update();
   if ( !pdf.empty() ) {
      pdf += (neg_pos ? "_n.pdf" : "_p.pdf");
      c1 -> Print(pdf.c_str());
   }
}

// {{{1 Combined fit
//-------------------------------------------------------------------------
// Combined fit
//-------------------------------------------------------------------------
//----------------------------------------------------------------------
vector<ValErr> CombIntEr( const ROOT::Fit::FitResult& res,
                          double sl09, double sl12, unsigned int idx ) {
//----------------------------------------------------------------------
// Numerical calculation of integrals and their errors
   // desired errors:
   constexpr double eps_abs = 1.e-6;
   constexpr double eps_rel = 1.e-6;
   constexpr double epsilon = 1.e-3; // max numerical error

   const double Mkk_max = dU; //dU -> std, test: 1.04 & 1.03

   // get parameters and covariation matrix of error
   int Npar = res.NPar();
   vector<double> Fp = res.Parameters();

   if ( Npar != 7 && Npar != 9 ) {
      cerr << " FATAL ERROR in " << __func__ << " Npar= " << Npar << endl;
      exit(1);
   }

   // parameters for 2009 and 2012 are considered separately
   int npar = Npar - 1; // 6
   if ( Npar == 9 ) {
      npar = Npar - 2;  // 7
   }

   vector<double> Fpar[2];
   if ( npar == 6 ) {
      Fpar[0] = { Fp[0],Fp[1],Fp[5],Fp[2],Fp[3],Fp[4],sl09 };
      Fpar[1] = { Fp[0],Fp[1],Fp[6],Fp[2],Fp[3],Fp[4],sl12 };
   } else {
      Fpar[0] = { Fp[0],Fp[1],Fp[5],Fp[2],Fp[3],Fp[4],Fp[7],sl09 };
      Fpar[1] = { Fp[0],Fp[1],Fp[6],Fp[2],Fp[3],Fp[4],Fp[8],sl12 };
   }

   int covMatrStatus = res.CovMatrixStatus();
   if ( covMatrStatus != 3 ) {
      cout << " WARNING: " << __func__ << " covariance matrix"
         " status code is " << covMatrStatus << endl;
   }
   vector<double> cov_m[2];
   cov_m[0].resize(npar*npar,0);
   cov_m[1].resize(npar*npar,0);
   vector<int> I09, I12;
   if ( npar == 6 ) {
      I09 = {0,1,3,4,5, 2,-1};
      I12 = {0,1,3,4,5,-1, 2};
   } else {
      I09 = {0,1,3,4,5, 2,-1, 6,-1};
      I12 = {0,1,3,4,5,-1, 2,-1, 6};
   }
   for ( int i = 0; i < Npar; ++i ) {
      int i09 = I09[i];
      int i12 = I12[i];
      for ( int j = 0; j < Npar; ++j ) {
         int j09 = I09[j];
         if ( i09 != -1 && j09 != -1 ) {
            cov_m[0][i09*npar+j09] = res.CovMatrix(i,j);
         }
         int j12 = I12[j];
         if ( i12 != -1 && j12 != -1 ) {
            cov_m[1][i12*npar+j12] = res.CovMatrix(i,j);
         }
      }
   }

   // 1) calculate integrals of normalized component
   auto Lfc = [idx,npar](const double* x,const double* p) -> double {
      double m = x[0];
      double sl = p[npar]; // slope!
      double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl };
      double normI = IntfrBWARGN( m,pp,100 );
      pp[6] = 0; // set slope to zero
      double intfr = IntfrBWARG( m,pp,idx ) / normI;
      if ( npar == 6 ) {
         return intfr;
      }
      return (1-p[6])*intfr;
   };
   TF1 FC = TF1("FC",Lfc, dL, dU, npar+1); // add sl in parameters

   double Integ[2] {0.,0.}, err_Int[2] {0.,0.};
   for ( int idat = 0; idat < 2; ++idat ) {
      FC.SetParameters( Fpar[idat].data() );
      double err_num = 0;
      Integ[idat] =
         FC.IntegralOneDim(dL,Mkk_max,eps_rel,eps_abs, err_num);
//       cout << " err_num(" << idat << ") = " << err_num << endl;
      if ( err_num > epsilon ) {
         cout << " WARNING: " << __func__ << " numerical error"
            " of integration of FC-" << idat << " is too big "
            << err_num << endl;
      }

      // 2) calculate error of this integral
      TMatrixDSym covMatrix(npar);
      covMatrix.Use(npar,cov_m[idat].data());
//       printf("\ncovMatrix-%d : ",idat);
//       covMatrix.Print();

      TVectorD IntGrad(npar);
      double err_num2 = 0;
      for ( int i = 0; i < npar; ++i ) {
         // skip parameters with zero error
         if ( covMatrix(i,i) == 0 ) {
            continue;
         }

         auto L_dF = [FC,i](const double* x, const double* p) -> double {
            TF1 tmp(FC); // may modify 'tmp' but not FC
            return tmp.GradientPar(i,x);
         };
         TF1 dF("dF",L_dF,dL,dU,0);
         double err_num = 0;
         IntGrad[i] =
            dF.IntegralOneDim(dL,Mkk_max,eps_rel,eps_abs, err_num);
         err_num = covMatrix(i,i) * IntGrad[i] * err_num;
//          cout << " abs err_num(dF_" << i << ") = " << err_num << endl;
         err_num2 += SQ(err_num);
//          cout << " IntGrad-"<<idat<<"[" << i << "] = "<<IntGrad[i]<<endl;
      }

      err_Int[idat] = sqrt( covMatrix.Similarity(IntGrad) );
      err_num2 = sqrt( err_num2 / err_Int[idat] ); // abs numerical error
//       cout << " err_num(dBW-" << idat << ") = " << err_num2 << endl;
      if ( err_num2 > epsilon ) {
         cout << " WARNING: " << __func__ << " numerical error"
            " of integration of dF-" << idat << " is too big "
            << err_num2 << endl;
      }
   }

   vector<ValErr> ret { ValErr {Integ[0], err_Int[0]},
                        ValErr {Integ[1], err_Int[1]}  };
   return ret;
}

//-------------------------------------------------------------------------
void combine_Intfr(string fname09, string fname12, string pdf="") {
//-------------------------------------------------------------------------
   // Get un-binned data
   TH1D* hist[2];
   vector<double> mkk09 = get_mkk_hist( fname09, "mkk_09_cp", &hist[0] );
   vector<double> mkk12 = get_mkk_hist( fname12, "mkk_12_cp", &hist[1] );
   TH1D* h09 = hist[0];
   TH1D* h12 = hist[1];

   int n09 = mkk09.size();
   int n12 = mkk12.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Dat(dr, n09+n12);
   for ( int i = 0; i < n09; ++i ) {
      Dat.Add(-mkk09[i]);
   }
   for ( int i = 0; i < n12; ++i ) {
      Dat.Add(mkk12[i]);
   }

   //-----------------------------------------------------------------------
   // Fit data
   // efficiency parameters
   double sl09 = -1.78;
   double sl12 = -1.87;
   // systematic study: (sl09 = -1.53,-2.03; sl12 = -1.72,-2.02)
//    double sl09 = -2.03;
//    double sl12 = -2.02;
//    printf(" sl09= %.2f, sl12= %.2f\n",sl09,sl12);
   //
   // function MUST be normalized to 1 on the fit range
   auto Lfit = [sl09,sl12](double* x,double* p) -> double {
      double m = x[0];
      if ( m < 0 ) { // 2009
         const double pp[] = { p[0],p[1],p[5],p[2],p[3],p[4],sl09 };
         return IntfrBWARGN( -m,pp,0 );
      } else { // 2012
         const double pp[] = { p[0],p[1],p[6],p[2],p[3],p[4],sl12 };
         return IntfrBWARGN( m,pp,0 );
      }
   };
   vector<string> par_name { "M#phi", "G#phi", "A", "F", "#vartheta",
                             "#sigma{09}", "#sigma{12}"  };
   vector<double> par_ini {Mphi,Gphi,0.,1.0,0.8,1.5e-3,1.2e-3};  // neg
//    vector<double> par_ini {Mphi,Gphi,0.,1.0,-0.8,1.5e-3,1.2e-3}; // pos

   bool neg_pos = par_ini[4] > 0; // true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters
//    cerr<< " Npar= " << Npar << endl;
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 WFit( *Ffit );
   fitter.SetFunction( WFit, false); // false == no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
//    fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   fitter.Config().ParSettings(0).SetValue(1.01953); // like MC-sig
   fitter.Config().ParSettings(0).Fix();                    // Mphi
   fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3, Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();                    // Gphi
   fitter.Config().ParSettings(2).SetLimits(-5., 5.);       // A
   fitter.Config().ParSettings(2).Fix();                    // A
   fitter.Config().ParSettings(3).SetLimits(0.01, 10.);     // F
   fitter.Config().ParSettings(4).SetLimits(-M_PI, M_PI);   // vartheta
   fitter.Config().ParSettings(5).SetLimits(0.2e-3, 2.e-3); // sigma09
   fitter.Config().ParSettings(6).SetLimits(0.2e-3, 2.e-3); // sigma12

   // == Fit
   fitter.LikelihoodFit( Dat, false ); // true == extended likelihood fit
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
//    res.PrintCovMatrix(cout); // print error matrix and correlations

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   vector<string> names { "N(KK)", "Nphi", "Nnonphi", "Nifr"  };
   vector<ValErr> Nos;
   cout << " n09= " << n09 << "  n12= " << n12 << endl;
   vector<double> norm { double(n09), double(n12) };
   Nos.reserve(8);
   for ( int idx = 0; idx <= 3; ++idx ) {
      vector<ValErr> Vinteg = CombIntEr( res,sl09,sl12,idx );
      for ( int j = 0; j < 2 ; ++j ) {
         const auto& Integ = Vinteg[j];
         double Num = norm[j] * Integ.val;
         double Err = fabs(Num) *
            sqrt( 1./norm[j] + SQ(Integ.err/Integ.val) );
         Nos.push_back( ValErr {Num,Err} );
      }
      printf("%s09 = %s %s12 = %s\n",
            names[idx].c_str(),Nos[2*idx].prt(".1f"),
            names[idx].c_str(),Nos[2*idx+1].prt(".1f") );
   }
   double Nphi09 = Nos[2].val, err_Nphi09 = Nos[2].err;
   double Nphi12 = Nos[3].val, err_Nphi12 = Nos[3].err;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof09 = [Fpar,sl09](double x) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[5],Fpar[2],Fpar[3],Fpar[4], sl09};
      return IntfrBWARGN( x,pp,0 );
   };
   rmath_fun< decltype(Lgof09) > ftest09(Lgof09);
   ROOT::Math::GoFTest* goftest09 = new ROOT::Math::GoFTest(
         mkk09.size(),mkk09.data(),ftest09,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS09 = goftest09 -> KolmogorovSmirnovTest();
   cout << " pvalueKS09= " << pvalueKS09 << endl;

   auto Lgof12 = [Fpar,sl12](double x) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[6],Fpar[2],Fpar[3],Fpar[4], sl12};
      return IntfrBWARGN( x,pp,0 );
   };
   rmath_fun< decltype(Lgof12) > ftest12(Lgof12);
   ROOT::Math::GoFTest* goftest12 = new ROOT::Math::GoFTest(
         mkk12.size(),mkk12.data(),ftest12,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS12 = goftest12 -> KolmogorovSmirnovTest();
   cout << " pvalueKS12= " << pvalueKS12 << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   auto Ldr09 = [n09,Fpar,sl09](double* x,double* p) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[5],Fpar[2],Fpar[3],Fpar[4],sl09};
      return bW*n09 * IntfrBWARGN( x[0],pp, int(p[0]) );
   };
   TF1* ff09 = new TF1("ff09", Ldr09, dL, dU, 1);
   ff09 -> SetNpx(500);

   auto Ldr12 = [n12,Fpar,sl12](double* x,double* p) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[6],Fpar[2],Fpar[3],Fpar[4],sl12};
      return bW*n12 * IntfrBWARGN( x[0],pp, int(p[0]) );
   };
   TF1* ff12 = new TF1("ff12", Ldr12, dL, dU, 1);
   ff12 -> SetNpx(500);

   // find max/min to draw
   bool presentation_plots = false;
   if ( presentation_plots ) {
      // must be the same parameters for positive and negative!
      h09 -> SetMinimum( -15.);
      h09 -> SetMaximum(+135.);
      h12 -> SetMinimum( -40.);
      h12 -> SetMaximum(+390.);
   } else {
      ff09 -> SetParameter(0, 1); // BW
      int mb09 = h09 -> GetMaximumBin();
      double hmax09 = h09 -> GetBinContent(mb09) + h09 -> GetBinError(mb09);
      double fmax09 = ff09 -> GetMaximum( 1.01, 1.03, 1e-6);
      hmax09 = floor(1.1 * max( hmax09, fmax09 ));
      if ( hmax09 > 0 ) {
         h09 -> SetMaximum(hmax09);
      }

      ff09 -> SetParameter(0, 3); // interference
      double hmin09 = ff09 -> GetMinimum( 1.01, 1.03, 1e-6);
      hmin09 = floor(1.15 * min( 0., hmin09 ));
      if ( hmin09 < 0 ) {
         h09 -> SetMinimum(hmin09);
      }

      ff12 -> SetParameter(0, 1); // BW
      int mb12 = h12 -> GetMaximumBin();
      double hmax12 = h12 -> GetBinContent(mb12) + h12 -> GetBinError(mb12);
      double fmax12 = ff12 -> GetMaximum( 1.01, 1.03, 1e-6);
      hmax12 = floor(1.1 * max( hmax12, fmax12 ));
      if ( hmax12 > 0 ) {
         h12 -> SetMaximum(hmax12);
      }

      ff12 -> SetParameter(0, 3); // interference
      double hmin12 = ff12 -> GetMinimum( 1.01, 1.03, 1e-6);
      hmin12 = floor(1.15 * min( 0., hmin12 ));
      if ( hmin12 < 0 ) {
         h12 -> SetMinimum(hmin12);
      }
   }

   //-----------------------------------------------------------------------
   // Draw results
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000);
   c1 -> Divide(1,2);

   c1 -> cd(1);
   gPad -> SetGrid();

   SetHstFace(h09);
   h09 -> GetYaxis() -> SetTitleOffset(1.);
   h09->SetLineWidth(2);
   h09->SetLineColor(kBlack);
   h09->SetMarkerStyle(20);

   h09->Draw("EP");

   TLegend* leg = new TLegend(0.51,0.45,0.89,0.89);
   leg -> SetHeader("#bf{2009(top)   2012(bottom)}","C");
   leg -> AddEntry( h09,"Data","LEP" );

   ff09 -> SetParameter(0, 0); // Sum
   ff09 -> SetLineWidth(2);
   ff09 -> SetLineColor(kRed+1);
   ff09 -> DrawCopy("SAME");
   leg -> AddEntry( ff09 -> Clone(), "Combined fit", "L" );

   ff09 -> SetParameter(0, 1); // BW
   ff09 -> SetLineWidth(2);
   ff09 -> SetLineStyle(kDashed);
   ff09 -> SetLineColor(kGreen+3);
   ff09 -> DrawCopy("SAME");
//    leg -> AddEntry( ff09 -> Clone(), "Breit-Wigner", "L");
   leg -> AddEntry( ff09 -> Clone(), "Breit-Wigner #phi#eta", "L");

   ff09 -> SetParameter(0, 2); // Argus
   ff09 -> SetLineWidth(1);
   ff09 -> SetLineColor(kBlue);
   ff09 -> DrawCopy("SAME");
//    leg -> AddEntry( ff09 -> Clone(), "Argus", "L");
   leg -> AddEntry( ff09 -> Clone(), "Non-#phi KK#eta", "L");

   ff09 -> SetParameter(0, 3); // interference
   ff09 -> SetLineWidth(1);
   ff09 -> SetLineColor(kMagenta+1);
   ff09 -> DrawCopy("SAME");
   leg -> AddEntry( ff09 -> Clone(), "Interference", "L");
   leg -> Draw();
   gPad -> RedrawAxis();

   c1 -> cd(2);
   gPad -> SetGrid();
   SetHstFace(h12);
   h12 -> GetYaxis() -> SetTitleOffset(1.);
   h12 -> SetLineWidth(2);
   h12 -> SetLineColor(kBlack);
   h12 -> SetMarkerStyle(20);

   h12 -> Draw("EP");

   ff12 -> SetParameter(0, 0); // Sum
   ff12 -> SetLineWidth(2);
   ff12 -> SetLineColor(kRed+1);
   ff12 -> DrawCopy("SAME");

   ff12 -> SetParameter(0, 1); // BW
   ff12 -> SetLineWidth(2);
   ff12 -> SetLineStyle(kDashed);
   ff12 -> SetLineColor(kGreen+3);
   ff12 -> DrawCopy("SAME");

   ff12 -> SetParameter(0, 2); // Argus
   ff12 -> SetLineWidth(1);
   ff12 -> SetLineColor(kBlue);
   ff12 -> DrawCopy("SAME");

   ff12 -> SetParameter(0, 3); // interference
   ff12 -> SetLineWidth(1);
   ff12 -> SetLineColor(kMagenta+1);
   ff12 -> DrawCopy("SAME");

   TPaveText* pt09 = new TPaveText(0.51,0.82,0.89,0.99,"NDC");
   pt09 -> SetTextAlign(12);
   pt09 -> SetTextFont(42);
   pt09 -> AddText( Form("#it{p-value(2009) = %.3f}", pvalueKS09) );
   pt09 -> AddText( Form("N_{#phi}(2009) = %.1f #pm %.1f",
                    Nphi09,err_Nphi09) );
   pt09 -> AddText( Form("#sigma(2009)= %s MeV", PE.Eform(5,".2f",1e3)) );
   pt09 -> Draw();

   TPaveText* pt12 = new TPaveText(0.51,0.65,0.89,0.82,"NDC");
   pt12 -> SetTextAlign(12);
   pt12 -> SetTextFont(42);
   pt12 -> AddText( Form("#it{p-value(2012) = %.3f}", pvalueKS12) );
   pt12 -> AddText( Form("N_{#phi}(2012) = %.1f #pm %.1f",
                    Nphi12,err_Nphi12) );
   pt12 -> AddText( Form("#sigma(2012)= %s MeV", PE.Eform(6,".2f",1e3)) );
   pt12 -> Draw();

   TPaveText* pt = new TPaveText(0.51,0.35,0.89,0.65,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
//    pt -> AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt -> AddText( Form("a= %s",PE.Eform(2,".2f")) );
   pt -> AddText( Form("F= %s",PE.Eform(3,".2f")) );
   pt -> AddText( Form("#vartheta= %s",PE.Eform(4,".2f")) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      pdf += (neg_pos ? "_n.pdf" : "_p.pdf");
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void combine_Intfr_bg(string fname09, string fname12, string pdf="") {
//-------------------------------------------------------------------------
   // Get un-binned data
   TH1D* hist[2];
   vector<double> mkk09 = get_mkk_hist( fname09, "mkk_09_cp", &hist[0] );
   vector<double> mkk12 = get_mkk_hist( fname12, "mkk_12_cp", &hist[1] );
   TH1D* h09 = hist[0];
   TH1D* h12 = hist[1];

   int n09 = mkk09.size();
   int n12 = mkk12.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Dat(dr, n09+n12);
   for ( int i = 0; i < n09; ++i ) {
      Dat.Add(-mkk09[i]);
   }
   for ( int i = 0; i < n12; ++i ) {
      Dat.Add(mkk12[i]);
   }

   //-----------------------------------------------------------------------
   // Fit data
   // efficiency parameters
   double sl09 = -1.78;
   double sl12 = -1.87;
   // function MUST be normalized to 1 on the fit range
   auto Lfit = [sl09,sl12](double* x,double* p) -> double {
      double m = x[0];
      if ( m < 0 ) { // 2009
         double pp[] = { p[0],p[1],p[5],p[2],p[3],p[4],sl09 };
         return (1-p[7])*IntfrBWARGN( -m,pp,0 ) + p[7]/(dU-dL);
      } else { // 2012
         double pp[] = { p[0],p[1],p[6],p[2],p[3],p[4],sl12 };
         return (1-p[8])*IntfrBWARGN( m,pp,0 ) + p[8]/(dU-dL);
      }
   };
   vector<string> par_name { "M#phi", "G#phi", "A", "F", "#vartheta",
                             "#sigma{09}", "#sigma{12}", "Bg09", "Bg12" };

   vector<double> par_ini { Mphi, Gphi, 0., 1.0, 0.8,
                            1.5e-3,1.2e-3, 1e-3,1e-3 };         // neg
//    vector<double> par_ini { Mphi, Gphi, 0., 1.0, -0.8,
//                             1.5e-3,1.2e-3, 3.e-3,1e-4 };        // pos

   bool neg_pos = par_ini[4] > 0; // true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 WFit( *Ffit );
   fitter.SetFunction( WFit, false); // false == no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   fitter.Config().ParSettings(0).SetValue(1.01953); // like MC-sig
   fitter.Config().ParSettings(0).Fix();                    // Mphi
   fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3, Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();                    // Gphi
   fitter.Config().ParSettings(2).Fix();                    // A
   fitter.Config().ParSettings(3).SetLimits(0.01, 10.);     // F
   fitter.Config().ParSettings(4).SetLimits(-M_PI, M_PI);   // vartheta
   fitter.Config().ParSettings(5).SetLimits(0.2e-3, 2.e-3); // sigma09
   fitter.Config().ParSettings(6).SetLimits(0.2e-3, 2.e-3); // sigma12
   fitter.Config().ParSettings(7).SetLimits(0., 1.);        // Bg09
//    fitter.Config().ParSettings(7).SetValue(0);
//    fitter.Config().ParSettings(7).Fix();
   fitter.Config().ParSettings(8).SetLimits(0., 1.);        // Bg12
//    fitter.Config().ParSettings(8).SetValue(0);
//    fitter.Config().ParSettings(8).Fix();

   // == Fit
   fitter.LikelihoodFit( Dat, false ); // true == extended likelihood fit
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
//    res.PrintCovMatrix(cout); // print error matrix and correlations

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   vector<string> names { "N(KK)", "Nphi", "Nnonphi", "Nifr"  };
   vector<ValErr> Nos;
   cout << " n09= " << n09 << "  n12= " << n12 << endl;
   vector<double> norm { double(n09), double(n12) };
   Nos.reserve(8);
   for ( int idx = 0; idx <= 3; ++idx ) {
      vector<ValErr> Vinteg = CombIntEr( res,sl09,sl12,idx );
      for ( int j = 0; j < 2 ; ++j ) {
         const auto& Integ = Vinteg[j];
//          printf("Int[%i,%i] = %s\n",idx,j,Integ.prt(".4f")); // DEBUG
         double Num = norm[j] * Integ.val;
         double Err = fabs(Num) *
            sqrt( 1./norm[j] + SQ(Integ.err/Integ.val) );
         Nos.push_back( ValErr {Num,Err} );
      }
      printf("%s09 = %s %s12 = %s\n",
            names[idx].c_str(),Nos[2*idx].prt(".1f"),
            names[idx].c_str(),Nos[2*idx+1].prt(".1f") );
   }
   double Nphi09 = Nos[2].val, err_Nphi09 = Nos[2].err;
   double Nphi12 = Nos[3].val, err_Nphi12 = Nos[3].err;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof09 = [Fpar,sl09](double x) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[5],Fpar[2],Fpar[3],Fpar[4], sl09};
      return (1-Fpar[7])*IntfrBWARGN( x,pp,0 ) + Fpar[7]/(dU-dL);
   };
   rmath_fun< decltype(Lgof09) > ftest09(Lgof09);
   ROOT::Math::GoFTest* goftest09 = new ROOT::Math::GoFTest(
         mkk09.size(),mkk09.data(),ftest09,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS09 = goftest09 -> KolmogorovSmirnovTest();
   cout << " pvalueKS09= " << pvalueKS09 << endl;

   auto Lgof12 = [Fpar,sl12](double x) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[6],Fpar[2],Fpar[3],Fpar[4], sl12};
      return (1-Fpar[8])*IntfrBWARGN( x,pp,0 ) + Fpar[8]/(dU-dL);
   };
   rmath_fun< decltype(Lgof12) > ftest12(Lgof12);
   ROOT::Math::GoFTest* goftest12 = new ROOT::Math::GoFTest(
         mkk12.size(),mkk12.data(),ftest12,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS12 = goftest12 -> KolmogorovSmirnovTest();
   cout << " pvalueKS12= " << pvalueKS12 << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   auto Ldr09 = [n09,Fpar,sl09](double* x,double* p) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[5],Fpar[2],Fpar[3],Fpar[4],sl09};
      return bW*n09 * (1-Fpar[7])*IntfrBWARGN( x[0],pp, int(p[0]) );
   };
   TF1* f09 = new TF1("f09", Ldr09, dL, dU, 1);
   f09 -> SetNpx(500);

   auto Ldr09_bg = [n09,Fpar](double* x,double* p) -> double {
      return bW*n09 * Fpar[7]/(dU-dL);
   };
   TF1* f09_bg = new TF1("f09_bg", Ldr09_bg, dL, dU, 0);
   cout << " N(f09_bg)= " << f09_bg -> Integral(dL,dU,1e-2)/bW << endl;

   auto Ldr12 = [n12,Fpar,sl12](double* x,double* p) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[6],Fpar[2],Fpar[3],Fpar[4],sl12};
      return bW*n12 * (1-Fpar[8])*IntfrBWARGN( x[0],pp, int(p[0]) );
   };
   TF1* f12 = new TF1("f12", Ldr12, dL, dU, 1);
   f12 -> SetNpx(500);

   auto Ldr12_bg = [n12,Fpar](double* x,double* p) -> double {
      return bW*n12 * Fpar[8]/(dU-dL);
   };
   TF1* f12_bg = new TF1("f12_bg", Ldr12_bg, bL, dU, 0);
   cout << " N(f12_bg)= " << f12_bg -> Integral(dL,dU,1e-2)/bW << endl;

   // find max/min to draw
   f09 -> SetParameter(0, 1); // BW
   int mb09 = h09 -> GetMaximumBin();
   double hmax09 = h09 -> GetBinContent(mb09) + h09 -> GetBinError(mb09);
   double fmax09 = f09 -> GetMaximum( 1.01, 1.03, 1e-6);
   hmax09 = floor(1.1 * max( hmax09, fmax09 ));
   if ( hmax09 > 0 ) {
      h09 -> SetMaximum(hmax09);
   }

   f09 -> SetParameter(0, 3); // interference
   double hmin09 = f09 -> GetMinimum( 1.01, 1.03, 1e-6);
   hmin09 = floor(1.15 * min( 0., hmin09 ));
   if ( hmin09 < 0 ) {
      h09 -> SetMinimum(hmin09);
   }

   f12 -> SetParameter(0, 1); // BW
   int mb12 = h12 -> GetMaximumBin();
   double hmax12 = h12 -> GetBinContent(mb12) + h12 -> GetBinError(mb12);
   double fmax12 = f12 -> GetMaximum( 1.01, 1.03, 1e-6);
   hmax12 = floor(1.1 * max( hmax12, fmax12 ));
   if ( hmax12 > 0 ) {
      h12 -> SetMaximum(hmax12);
   }

   f12 -> SetParameter(0, 3); // interference
   double hmin12 = f12 -> GetMinimum( 1.01, 1.03, 1e-6);
   hmin12 = floor(1.15 * min( 0., hmin12 ));
   if ( hmin12 < 0 ) {
      h12 -> SetMinimum(hmin12);
   }

   //-----------------------------------------------------------------------
   // Draw results
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000);
   c1 -> Divide(1,2);

   c1 -> cd(1);
   gPad -> SetGrid();

   SetHstFace(h09);
   h09 -> GetYaxis() -> SetTitleOffset(1.);
   h09 -> SetLineWidth(2);
   h09 -> SetLineColor(kBlack);
   h09 -> SetMarkerStyle(20);

   h09->Draw("EP");

   TLegend* leg = new TLegend(0.51,0.45,0.89,0.89);
   leg -> SetHeader("#bf{2009(top)   2012(bottom)}","C");
   leg -> AddEntry( h09,"Data","LEP" );

   f09 -> SetParameter(0, 0); // Sum
   f09 -> SetLineWidth(2);
   f09 -> SetLineColor(kRed+1);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Combined fit", "L" );

   f09 -> SetParameter(0, 1); // BW
   f09 -> SetLineWidth(1);
   f09 -> SetLineStyle(kDashed);
   f09 -> SetLineColor(kGreen+2);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Breit-Wigner #phi#eta", "L");

   f09 -> SetParameter(0, 2); // Argus
   f09 -> SetLineColor(kBlue);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Non-#phi KK#eta", "L");

   f09 -> SetParameter(0, 3); // interference
   f09 -> SetLineColor(kMagenta+1);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Interference", "L");

//    f09_bg -> SetLineWidth(2);
//    f09_bg -> SetLineColor(kBlue+2);
//    f09_bg -> DrawCopy("SAME");
//    leg -> AddEntry( f09_bg, "Flat background", "L" );

   leg -> Draw();
   gPad -> RedrawAxis();

   c1 -> cd(2);
   gPad -> SetGrid();
   SetHstFace(h12);
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
   f12 -> SetLineWidth(1);
   f12 -> SetLineStyle(kDashed);
   f12 -> SetLineColor(kGreen+2);
   f12 -> DrawCopy("SAME");

   f12 -> SetParameter(0, 2); // Argus
   f12 -> SetLineColor(kBlue);
   f12 -> DrawCopy("SAME");

   f12 -> SetParameter(0, 3); // interference
   f12 -> SetLineColor(kMagenta+1);
   f12 -> DrawCopy("SAME");

//    f12_bg -> SetLineWidth(2);
//    f12_bg -> SetLineColor(kBlue+2);
//    f12_bg -> DrawCopy("SAME");

   TPaveText* pt09 = new TPaveText(0.51,0.79,0.89,0.99,"NDC");
   pt09 -> SetTextAlign(12);
   pt09 -> SetTextFont(42);
   pt09 -> AddText( Form("#it{p-value(2009) = %.3f}", pvalueKS09) );
   pt09 -> AddText( Form("N_{#phi}(2009) = %.1f #pm %.1f",
                    Nphi09,err_Nphi09) );
   pt09 -> AddText( Form("#sigma(2009)= %s MeV", PE.Eform(5,".2f",1e3)) );
   pt09 -> AddText( Form("Bkg(2009)= %s", PE.Eform(7,".4f")) );
   pt09 -> Draw();

   TPaveText* pt12 = new TPaveText(0.51,0.59,0.89,0.79,"NDC");
   pt12 -> SetTextAlign(12);
   pt12 -> SetTextFont(42);
   pt12 -> AddText( Form("#it{p-value(2012) = %.3f}", pvalueKS12) );
   pt12 -> AddText( Form("N_{#phi}(2012) = %.1f #pm %.1f",
                    Nphi12,err_Nphi12) );
   pt12 -> AddText( Form("#sigma(2012)= %s MeV", PE.Eform(6,".2f",1e3)) );
   pt12 -> AddText( Form("Bkg(2012)= %s", PE.Eform(8,".4f")) );
   pt12 -> Draw();

   TPaveText* pt = new TPaveText(0.51,0.34,0.89,0.59,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
   pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt -> AddText( Form("a= %s",PE.Eform(2,".2f")) );
   pt -> AddText( Form("F= %s",PE.Eform(3,".2f")) );
   pt -> AddText( Form("#vartheta= %s",PE.Eform(4,".2f")) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      pdf += (neg_pos ? "_n.pdf" : "_p.pdf");
      c1 -> Print(pdf.c_str());
   }
}

// combined plus side-band condition
//----------------------------------------------------------------------
struct myFCNsb {
   // efficiency parameters
   double sl09 = -1.78;
   double sl12 = -1.87;
   // data central part
   vector<double> mkk09;
   vector<double> mkk12;
   // data side band
   vector<double> sb09;
   vector<double> sb12;

   // minimization function
   // the signature of this operator() MUST be exactly this:
   // THE EXTENDED MAXIMUM LIKELIHOOD ESTIMATOR
   // Kai-Feng Chen "FittingMinuitRooFit" p.35
   //------------------------------------------------------------------
   double operator() (const double* p) {
   //------------------------------------------------------------------
      double Nkk09 = p[7];
      double Nsb09 = p[8];
      double res = 2*(Nkk09+Nsb09);

      const double p09[] = { p[0],p[1],p[5],p[2],p[3],p[4],sl09 };
      int n_mkk09 = mkk09.size();
      for ( int i = 0; i < n_mkk09; ++i ) {
         double m = mkk09[i];
         double L = Nkk09*IntfrBWARGN( m,p09,0 ) + Nsb09/(dU-dL);
         if (L > 0.) {
            res -= 2*log(L);
         } else {
            return DBL_MAX;
         }
      }

      int n_sb09 = sb09.size();
      res += 2*Nsb09 - n_sb09*2*log(Nsb09/(dU-dL)); // fit SB by constant

      double Nkk12 = p[9];
      double Nsb12 = p[10];
      res += 2*(Nkk12+Nsb12);

      const double p12[] = { p[0],p[1],p[6],p[2],p[3],p[4],sl12 };
      int n_mkk12 = mkk12.size();
      for ( int i = 0; i < n_mkk12; ++i ) {
         double m = mkk12[i];
         double L = Nkk12*IntfrBWARGN( m,p12,0 ) + Nsb12/(dU-dL);
         if (L > 0.) {
            res -= 2*log(L);
         } else {
            return DBL_MAX;
         }
      }

      int n_sb12 = sb12.size();
      res += 2*Nsb12 - n_sb12*2*log(Nsb12/(dU-dL)); // fit SB by constant

      return res;
   }
};
//----------------------------------------------------------------------

//-------------------------------------------------------------------------
void combine_Intfr_sb(string fname09, string fname12, string pdf="") {
//-------------------------------------------------------------------------

   myFCNsb my_fcn; // class for 'FitFCN'

   // efficiency parameters
   double sl09 = -1.78;
   double sl12 = -1.87;
   my_fcn.sl09 = sl09;
   my_fcn.sl12 = sl12;

   // Get un-binned data
   TH1D* hist[4];
   my_fcn.mkk09 = get_mkk_hist( fname09, "mkk_09_cp", &hist[0] );
   const vector<double>& mkk09 = my_fcn.mkk09;
   TH1D* h09 = hist[0];
   int n09 = mkk09.size();

   my_fcn.mkk12 = get_mkk_hist( fname12, "mkk_12_cp", &hist[1] );
   const vector<double>& mkk12 = my_fcn.mkk12;
   TH1D* h12 = hist[1];
   int n12 = mkk12.size();

   my_fcn.sb09 = get_mkk_hist( fname09, "mkk_sb09", &hist[2], 1 );
   const vector<double>& sb09 = my_fcn.sb09;
   TH1D* hsb09 = hist[2];
   int nsb09 = sb09.size();

   my_fcn.sb12 = get_mkk_hist( fname12, "mkk_sb12", &hist[3], 1 );
   const vector<double>& sb12 = my_fcn.sb12;
   TH1D* hsb12 = hist[3];
   int nsb12 = sb12.size();

   //-----------------------------------------------------------------------
   // Fit data
   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   vector<string> par_name { "M#phi", "G#phi", "A", "F", "#vartheta",
        "#sigma{09}", "#sigma{12}", "Nkk09", "Nbg09", "Nkk12", "Nbg12" };

   vector<double> par_ini { Mphi, Gphi, 0., 1.0, 0.8, 1.5e-3,1.2e-3,
                            double(n09-nsb09),double(nsb09),
                            double(n12-nsb12),double(nsb12) };     // neg
//    vector<double> par_ini { Mphi, Gphi, 0., 1.0, -0.8,
//                             1.5e-3,1.2e-3, 3.e-3,1e-4 };        // pos

   bool neg_pos = par_ini[4] > 0; // true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   fitter.Config().ParSettings(0).SetValue(1.01953); // like MC-sig
   fitter.Config().ParSettings(0).Fix();                    // Mphi
   fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3, Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();                    // Gphi
   fitter.Config().ParSettings(2).Fix();                    // A
   fitter.Config().ParSettings(3).SetLimits(0.01, 10.);     // F
   fitter.Config().ParSettings(4).SetLimits(-M_PI, M_PI);   // vartheta
   fitter.Config().ParSettings(5).SetLimits(0.2e-3, 2.e-3); // sigma09
   fitter.Config().ParSettings(6).SetLimits(0.2e-3, 2.e-3); // sigma12

   // == Fit
   int Ndat = n09+nsb09+n12+nsb12;
   fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // false == likelihood fit
//    fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
//    res.PrintCovMatrix(cout); // print error matrix and correlations

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   cout << " n09= " << n09 << "  n12= " << n12;
   cout << " nsb09= " << nsb09 << "  nsb12= " << nsb12 << endl;
//    vector<string> names { "N(KK)", "Nphi", "Nnonphi", "Nifr"  };
//    vector<ValErr> Nos;
//    vector<double> norm { double(n09), double(n12) };
//    Nos.reserve(8);
//    for ( int idx = 0; idx <= 3; ++idx ) {
//       vector<ValErr> Vinteg = CombIntEr( res,sl09,sl12,idx );
//       for ( int j = 0; j < 2 ; ++j ) {
//          const auto& Integ = Vinteg[j];
//          printf("Int[%i,%i] = %s\n",idx,j,Integ.prt(".4f")); // DEBUG
//          double Num = norm[j] * Integ.val;
//          double Err = fabs(Num) *
//             sqrt( 1./norm[j] + SQ(Integ.err/Integ.val) );
//          Nos.push_back( ValErr {Num,Err} );
//       }
//       printf("%s09 = %s %s12 = %s\n",
//             names[idx].c_str(),Nos[2*idx].prt(".1f"),
//             names[idx].c_str(),Nos[2*idx+1].prt(".1f") );
//    }
//    double Nphi09 = Nos[2].val, err_Nphi09 = Nos[2].err;
//    double Nphi12 = Nos[3].val, err_Nphi12 = Nos[3].err;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lcr09 = [Fpar,sl09](double x) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[5],Fpar[2],Fpar[3],Fpar[4], sl09};
      return Fpar[7]*IntfrBWARGN( x,pp,0 ) + Fpar[8]/(dU-dL);
   };
   rmath_fun< decltype(Lcr09) > fcr09(Lcr09);
   ROOT::Math::GoFTest* goftest09 = new ROOT::Math::GoFTest(
         mkk09.size(),mkk09.data(),fcr09,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS09 = goftest09 -> KolmogorovSmirnovTest();
   cout << " pvalueKS09= " << pvalueKS09 << endl;

   auto Lsb09 = [Fpar](double x) -> double {
      return Fpar[8]/(dU-dL);
   };
   rmath_fun< decltype(Lsb09) > fsb09(Lsb09);
   ROOT::Math::GoFTest* goftest09sb = new ROOT::Math::GoFTest(
         sb09.size(),sb09.data(),fsb09,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS09sb = goftest09sb -> KolmogorovSmirnovTest();
   cout << " pvalueKS09sb= " << pvalueKS09sb << endl;

   auto Lcr12 = [Fpar,sl12](double x) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[6],Fpar[2],Fpar[3],Fpar[4], sl12};
      return Fpar[9]*IntfrBWARGN( x,pp,0 ) + Fpar[10]/(dU-dL);
   };
   rmath_fun< decltype(Lcr12) > fcr12(Lcr12);
   ROOT::Math::GoFTest* goftest12 = new ROOT::Math::GoFTest(
         mkk12.size(),mkk12.data(),fcr12,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS12 = goftest12 -> KolmogorovSmirnovTest();
   cout << " pvalueKS12= " << pvalueKS12 << endl;

   auto Lsb12 = [Fpar](double x) -> double {
      return Fpar[10]/(dU-dL);
   };
   rmath_fun< decltype(Lsb12) > fsb12(Lsb12);
   ROOT::Math::GoFTest* goftest12sb = new ROOT::Math::GoFTest(
         sb12.size(),sb12.data(),fsb12,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS12sb = goftest12sb -> KolmogorovSmirnovTest();
   cout << " pvalueKS12sb= " << pvalueKS12sb << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   auto Ldr09 = [Fpar,sl09](double* x,double* p) -> double {
      double pp[] {Fpar[0],Fpar[1],Fpar[5],Fpar[2],Fpar[3],Fpar[4],sl09};
      return bW * Fpar[7]*IntfrBWARGN( x[0],pp,int(p[0]) );
   };
   TF1* f09 = new TF1("f09", Ldr09, dL, dU, 1);
   f09 -> SetNpx(500);

   auto Ldr09_bg = [Fpar](double* x,double* p) -> double {
      return bW * Fpar[8]/(dU-dL);
   };
   TF1* f09_bg = new TF1("f09_bg", Ldr09_bg, dL, dU, 0);

   auto Ldr12 = [Fpar,sl12](double* x,double* p) -> double {
      double pp[] {Fpar[0],Fpar[1],Fpar[6],Fpar[2],Fpar[3],Fpar[4],sl12};
      return bW * Fpar[9]*IntfrBWARGN( x[0],pp,int(p[0]) );
   };
   TF1* f12 = new TF1("f12", Ldr12, dL, dU, 1);
   f12 -> SetNpx(500);

   auto Ldr12_bg = [Fpar](double* x,double* p) -> double {
      return bW * Fpar[10]/(dU-dL);
   };
   TF1* f12_bg = new TF1("f12_bg", Ldr12_bg, bL, dU, 0);

   // find max/min to draw
   f09 -> SetParameter(0, 1); // BW
   int mb09 = h09 -> GetMaximumBin();
   double hmax09 = h09 -> GetBinContent(mb09) + h09 -> GetBinError(mb09);
   double fmax09 = f09 -> GetMaximum( 1.01, 1.03, 1e-6);
   hmax09 = floor(1.1 * max( hmax09, fmax09 ));
   if ( hmax09 > 0 ) {
      h09 -> SetMaximum(hmax09);
   }

   f09 -> SetParameter(0, 3); // interference
   double hmin09 = f09 -> GetMinimum( 1.01, 1.03, 1e-6);
   hmin09 = floor(1.15 * min( 0., hmin09 ));
   if ( hmin09 < 0 ) {
      h09 -> SetMinimum(hmin09);
   }

   f12 -> SetParameter(0, 1); // BW
   int mb12 = h12 -> GetMaximumBin();
   double hmax12 = h12 -> GetBinContent(mb12) + h12 -> GetBinError(mb12);
   double fmax12 = f12 -> GetMaximum( 1.01, 1.03, 1e-6);
   hmax12 = floor(1.1 * max( hmax12, fmax12 ));
   if ( hmax12 > 0 ) {
      h12 -> SetMaximum(hmax12);
   }

   f12 -> SetParameter(0, 3); // interference
   double hmin12 = f12 -> GetMinimum( 1.01, 1.03, 1e-6);
   hmin12 = floor(1.15 * min( 0., hmin12 ));
   if ( hmin12 < 0 ) {
      h12 -> SetMinimum(hmin12);
   }

   //-----------------------------------------------------------------------
   // Draw results
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000);
   c1 -> Divide(1,2);

   c1 -> cd(1);
   gPad -> SetGrid();

   SetHstFace(h09);
   h09 -> GetYaxis() -> SetTitleOffset(1.);
   h09 -> SetLineWidth(2);
   h09 -> SetLineColor(kBlack);
   h09 -> SetMarkerStyle(20);

   h09->Draw("EP");

   TLegend* leg = new TLegend(0.51,0.45,0.89,0.89);
   leg -> SetHeader("#bf{2009(top)   2012(bottom)}","C");
   leg -> AddEntry( h09,"Data","LEP" );

   f09 -> SetParameter(0, 0); // Sum
   f09 -> SetLineWidth(2);
   f09 -> SetLineColor(kRed+1);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Combined fit", "L" );

   f09 -> SetParameter(0, 1); // BW
   f09 -> SetLineWidth(1);
   f09 -> SetLineStyle(kDashed);
   f09 -> SetLineColor(kGreen+2);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Breit-Wigner #phi#eta", "L");

   f09 -> SetParameter(0, 2); // Argus
   f09 -> SetLineColor(kBlue);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Non-#phi KK#eta", "L");

   f09 -> SetParameter(0, 3); // interference
   f09 -> SetLineColor(kMagenta+1);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Interference", "L");

//    f09_bg -> SetLineWidth(2);
//    f09_bg -> SetLineColor(kBlue+2);
//    f09_bg -> DrawCopy("SAME");
//    leg -> AddEntry( f09_bg, "Background", "L" );

   leg -> Draw();
   gPad -> RedrawAxis();

   c1 -> cd(2);
   gPad -> SetGrid();
   SetHstFace(h12);
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
   f12 -> SetLineWidth(1);
   f12 -> SetLineStyle(kDashed);
   f12 -> SetLineColor(kGreen+2);
   f12 -> DrawCopy("SAME");

   f12 -> SetParameter(0, 2); // Argus
   f12 -> SetLineColor(kBlue);
   f12 -> DrawCopy("SAME");

   f12 -> SetParameter(0, 3); // interference
   f12 -> SetLineColor(kMagenta+1);
   f12 -> DrawCopy("SAME");

//    f12_bg -> SetLineWidth(2);
//    f12_bg -> SetLineColor(kBlue+2);
//    f12_bg -> DrawCopy("SAME");

   TPaveText* pt09 = new TPaveText(0.51,0.79,0.89,0.99,"NDC");
   pt09 -> SetTextAlign(12);
   pt09 -> SetTextFont(42);
   pt09 -> AddText( Form("#it{p-value(2009) = %.3f}", pvalueKS09) );
   pt09 -> AddText( Form("#it{p-v(sb2009) = %.3f}", pvalueKS09sb) );
//    pt09 -> AddText( Form("N_{#phi}(2009) = %.1f #pm %.1f",
//                     Nphi09,err_Nphi09) );
   pt09 -> AddText( Form("#sigma(2009)= %s MeV", PE.Eform(5,".2f",1e3)) );
   pt09 -> AddText( Form("NKK(2009)= %s", PE.Eform(7,".1f")) );
   pt09 -> AddText( Form("Nbg(2009)= %s", PE.Eform(8,".1f")) );
   pt09 -> Draw();

   TPaveText* pt12 = new TPaveText(0.51,0.59,0.89,0.79,"NDC");
   pt12 -> SetTextAlign(12);
   pt12 -> SetTextFont(42);
   pt12 -> AddText( Form("#it{p-value(2012) = %.3f}", pvalueKS12) );
   pt12 -> AddText( Form("#it{p-v(sb2012) = %.3f}", pvalueKS12sb) );
//    pt12 -> AddText( Form("N_{#phi}(2012) = %.1f #pm %.1f",
//                     Nphi12,err_Nphi12) );
   pt12 -> AddText( Form("#sigma(2012)= %s MeV", PE.Eform(6,".2f",1e3)) );
   pt12 -> AddText( Form("NKK(2012)= %s", PE.Eform(9,".1f")) );
   pt12 -> AddText( Form("Nbg(2012)= %s", PE.Eform(10,".1f")) );
   pt12 -> Draw();

   TPaveText* pt = new TPaveText(0.51,0.34,0.89,0.59,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
   pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt -> AddText( Form("a= %s",PE.Eform(2,".2f")) );
   pt -> AddText( Form("F= %s",PE.Eform(3,".2f")) );
   pt -> AddText( Form("#vartheta= %s",PE.Eform(4,".2f")) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      pdf += (neg_pos ? "_n.pdf" : "_p.pdf");
      c1 -> Print(pdf.c_str());
   }
}

// {{{1 Interference + incoherent sum
//-------------------------------------------------------------------------
// Interference + incoherent sum
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void data_Intfr_Ar(string fname, string title, string pdf="") {
//-------------------------------------------------------------------------
   bool is2009 = (fname.find("_09") != string::npos);

   // Get un-binned data
   TH1D* hist[1];
   vector<double> mkk = get_mkk_hist( fname, "mkk_data_cp", hist );
   TH1D* hst = hist[0];

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Dat(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Dat.Add(mkk[i]);
   }

   //-----------------------------------------------------------------------
   // Fit data
   // efficiency parameters
   double sl = -1.87;
   if ( is2009 ) {
      sl = -1.78;
   }
   // function MUST be normalized to 1 on the fit range
   auto Lfit = [sl](const double* x,const double* p) -> double {
      const double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl };
      double intfr = IntfrBWARGN(x[0],pp,0);
      double bg = RevArgusN(x[0],0.,sl); // only slope!
      return (1-p[6])*intfr + p[6]*bg;
   };
   vector<string> par_name { "M#phi", "G#phi", "#sigma", "A",
                             "F", "#vartheta", "Abg" };

   vector<double> par_ini {Mphi,Gphi,1.5e-3,0.,1.0,1.0,0.}; // 09neg 1
//    vector<double> par_ini {Mphi,Gphi,1.5e-3,0.,1.0,-1.1,0.}; // 09 pos 1
//    vector<double> par_ini {Mphi,Gphi,1.2e-3,0.,1.0,0.7,0.}; // 12 neg 1
//    vector<double> par_ini {Mphi,Gphi,1.2e-3,0.,1.0,-0.7,0.}; // 12 pos 1

   bool neg_pos = par_ini[5] > 0; // true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 WFit( *Ffit );
   fitter.SetFunction( WFit, false); // false == no parameter derivatives

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
   fitter.Config().ParSettings(6).SetLimits(0., 1.);        // Abg
//    fitter.Config().ParSettings(6).Fix();

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
   cout << " n= " << n << endl;
   for ( int idx = 0; idx <= 3; ++idx ) {
      ValErr Integ = CalcIntEr( res, sl, idx );
      double Num = n * Integ.val;
      double Err = fabs(Num) * sqrt( 1./n + SQ(Integ.err/Integ.val) );
      Nos.push_back( ValErr {Num,Err} );
      printf("%s = %s\n",names[idx].c_str(),Nos[idx].prt(".1f"));
   }
   double Nphi = Nos[1].val, err_Nphi = Nos[1].err;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof = [Lfit,Fpar](double x) -> double {
      return Lfit( &x,Fpar.data() );
   };
   rmath_fun< decltype(Lgof) > ftest(Lgof); // adapter lambda -> Root1Dfunc
   ROOT::Math::GoFTest* goftest = new ROOT::Math::GoFTest(
         mkk.size(),mkk.data(),ftest,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS = goftest -> KolmogorovSmirnovTest();
   cout << " pvalueKS= " << pvalueKS << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   auto Ldr = [n,Fpar,sl](double* x,double* p) -> double {
      const double pp[]
        { Fpar[0],Fpar[1],Fpar[2],Fpar[3],Fpar[4],Fpar[5],sl };
      return bW*n*(1-Fpar[6])*IntfrBWARGN( x[0],pp,int(p[0]) );
   };
   TF1* ffit = new TF1("ffit", Ldr, dL, dU, 1);
   ffit -> SetNpx(500);

   auto Lfbg = [n,Fpar,sl](double* x,double* p) -> double {
      return bW*n*Fpar[6]*RevArgusN( x[0],0.,sl );
   };
   TF1* fbg = new TF1("fbg", Lfbg, bL, dU, 0);

   // find max/min to draw
   ffit -> SetParameter(0, 1); // BW
   int maxbin = hst -> GetMaximumBin();
   double hmax = hst -> GetBinContent(maxbin) + hst -> GetBinError(maxbin);
   double fmax = ffit -> GetMaximum( 1.01, 1.03, 1e-6);
   hmax = floor(1.1 * max( hmax, fmax ));
   if ( hmax > 0 ) {
      hst -> SetMaximum(hmax);
   }

   ffit -> SetParameter(0, 3); // interference
   double hmin = ffit -> GetMinimum( 1.01, 1.03, 1e-6);
   hmin = floor(1.15 * min( 0., hmin ));
   if ( hmin < 0 ) {
      hst -> SetMinimum(hmin);
   }

   //-----------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   SetHstFace(hst);

   hst -> GetYaxis() -> SetTitleOffset(1.3);
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
   hst -> SetMarkerStyle(20);

   hst -> Draw("EP");

   TLegend* leg = new TLegend(0.52,0.67,0.89,0.89);
   leg -> AddEntry(hst,title.c_str(),"LEP");

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
   leg -> AddEntry( ffit -> Clone(), "Breit-Wigner", "L");

   ffit -> SetParameter(0, 2); // interfering Argus
   ffit -> SetLineColor(kBlue);
   ffit -> DrawCopy("SAME");
   leg -> AddEntry( ffit -> Clone(), "Interfering Argus", "L");

   ffit -> SetParameter(0, 3); // interference
   ffit -> SetLineColor(kMagenta+1);
   ffit -> DrawCopy("SAME");
   leg -> AddEntry( ffit -> Clone(), "Interference", "L");

   fbg -> SetLineWidth(2); // non-interfering bkg
   fbg -> SetLineColor(kBlue+2);
   fbg -> DrawCopy("SAME");
   leg -> AddEntry( fbg, "Non-interfering bkg", "L" );

   leg -> Draw();

   TPaveText* pt = new TPaveText(0.52,0.33,0.89,0.66,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
//    pt -> AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt -> AddText( Form("#it{p-value(K-S) = %.3f}",pvalueKS) );
   pt -> AddText( Form("N_{#phi} = %.1f #pm %.1f",Nphi,err_Nphi) );
   pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt -> AddText( Form("#sigma= %s MeV", PE.Eform(2,".2f",1e3)) );
//    pt -> AddText( Form("a= %s",PE.Eform(3,".1f")) );
   pt -> AddText( Form("F= %s",PE.Eform(4,".2f")) );
   pt -> AddText( Form("#vartheta= %s",PE.Eform(5,".2f")) );
   pt -> AddText( Form("Bkg= %s", PE.Eform(6,".4f")) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      pdf += (neg_pos ? "_n.pdf" : "_p.pdf");
      c1 -> Print(pdf.c_str());
   }
}

//----------------------------------------------------------------------
vector<ValErr> CombIntErAr( const ROOT::Fit::FitResult& res,
                          double sl09, double sl12, unsigned int idx ) {
//----------------------------------------------------------------------
// Numerical calculation of integrals and their errors
   // desired errors:
   constexpr double eps_abs = 1.e-6;
   constexpr double eps_rel = 1.e-6;
   constexpr double epsilon = 1.e-3; // max numerical error

   const double Mkk_max = dU; //dU -> std, test: 1.04 & 1.03

   // get parameters and covariation matrix of error
   int Npar = res.NPar();
   vector<double> Fp = res.Parameters();

//    if ( Npar != 7 && Npar != 8 ) {
   if ( Npar != 8 ) { // only for combine_Intfr_Ar
      cerr << " FATAL ERROR in " << __func__ << " Npar= " << Npar << endl;
      exit(1);
   }

   // parameters for 2009 and 2012 are considered separately
   int npar = Npar - 1;
   vector<double> Fpar[2];
   if ( npar == 6 ) {
      Fpar[0] = { Fp[0],Fp[1],Fp[5],Fp[2],Fp[3],Fp[4],sl09 };
      Fpar[1] = { Fp[0],Fp[1],Fp[6],Fp[2],Fp[3],Fp[4],sl12 };
   } else {
      Fpar[0] = { Fp[0],Fp[1],Fp[5],Fp[2],Fp[3],Fp[4],Fp[7],sl09 };
      Fpar[1] = { Fp[0],Fp[1],Fp[6],Fp[2],Fp[3],Fp[4],Fp[7],sl12 };
   }

   int covMatrStatus = res.CovMatrixStatus();
   if ( covMatrStatus != 3 ) {
      cout << " WARNING: " << __func__ << " covariance matrix"
         " status code is " << covMatrStatus << endl;
   }
   vector<double> cov_m[2];
   cov_m[0].resize(npar*npar,0);
   cov_m[1].resize(npar*npar,0);
   vector<int> I09, I12;
   if ( npar == 6 ) {
      I09 = {0,1,3,4,5, 2,-1};
      I12 = {0,1,3,4,5,-1, 2};
   } else {
      I09 = {0,1,3,4,5, 2,-1, 6};
      I12 = {0,1,3,4,5,-1, 2, 6};
   }
   for ( int i = 0; i < Npar; ++i ) {
      int i09 = I09[i];
      int i12 = I12[i];
      for ( int j = 0; j < Npar; ++j ) {
         int j09 = I09[j];
         if ( i09 != -1 && j09 != -1 ) {
            cov_m[0][i09*npar+j09] = res.CovMatrix(i,j);
         }
         int j12 = I12[j];
         if ( i12 != -1 && j12 != -1 ) {
            cov_m[1][i12*npar+j12] = res.CovMatrix(i,j);
         }
      }
   }

   // 1) calculate integrals of normalized component
   auto Lfc = [idx,npar](const double* x,const double* p) -> double {
      double m = x[0];
      double sl = p[npar]; // slope!
      double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl };
      double normI = IntfrBWARGN( m,pp,100 );
      pp[6] = 0; // set slope to zero
      double intfr = IntfrBWARG( m,pp,idx ) / normI;
      if ( npar == 6 ) {
         return intfr;
      }
      if ( idx != 0 ) {
         return (1-p[6])*intfr;
      }
      double bg = RevArgusN(m,0.,sl) / (1+sl*(m-1.02)); // forget sl
      return (1-p[6])*intfr + p[6]*bg;
   };
   TF1 FC = TF1("FC",Lfc, dL, dU, npar+1); // add sl in parameters

   double Integ[2] {0.,0.}, err_Int[2] {0.,0.};
   for ( int idat = 0; idat < 2; ++idat ) {
      FC.SetParameters( Fpar[idat].data() );
      double err_num = 0;
      Integ[idat] =
         FC.IntegralOneDim(dL,Mkk_max,eps_rel,eps_abs, err_num);
//       cout << " err_num(" << idat << ") = " << err_num << endl;
      if ( err_num > epsilon ) {
         cout << " WARNING: " << __func__ << " numerical error"
            " of integration of FC-" << idat << " is too big "
            << err_num << endl;
      }

      // 2) calculate error of this integral
      TMatrixDSym covMatrix(npar);
      covMatrix.Use(npar,cov_m[idat].data());
//       printf("\ncovMatrix-%d : ",idat);
//       covMatrix.Print();

      TVectorD IntGrad(npar);
      double err_num2 = 0;
      for ( int i = 0; i < npar; ++i ) {
         // skip parameters with zero error
         if ( covMatrix(i,i) == 0 ) {
            continue;
         }

         auto L_dF = [FC,i](const double* x, const double* p) -> double {
            TF1 tmp(FC); // may modify 'tmp' but not FC
            return tmp.GradientPar(i,x);
         };
         TF1 dF("dF",L_dF,dL,dU,0);
         double err_num = 0;
         IntGrad[i] =
            dF.IntegralOneDim(dL,Mkk_max,eps_rel,eps_abs, err_num);
         err_num = covMatrix(i,i) * IntGrad[i] * err_num;
//          cout << " abs err_num(dF_" << i << ") = " << err_num << endl;
         err_num2 += SQ(err_num);
//          cout << " IntGrad-"<<idat<<"[" << i << "] = "<<IntGrad[i]<<endl;
      }

      err_Int[idat] = sqrt( covMatrix.Similarity(IntGrad) );
      err_num2 = sqrt( err_num2 / err_Int[idat] ); // abs numerical error
//       cout << " err_num(dBW-" << idat << ") = " << err_num2 << endl;
      if ( err_num2 > epsilon ) {
         cout << " WARNING: " << __func__ << " numerical error"
            " of integration of dF-" << idat << " is too big "
            << err_num2 << endl;
      }
   }

   vector<ValErr> ret { ValErr {Integ[0], err_Int[0]},
                        ValErr {Integ[1], err_Int[1]}  };
   return ret;
}

//-------------------------------------------------------------------------
void combine_Intfr_Ar(string fname09, string fname12, string pdf="") {
//-------------------------------------------------------------------------
   // Get un-binned data
   TH1D* hist[2];
   vector<double> mkk09 = get_mkk_hist( fname09, "mkk_09_cp", &hist[0] );
   vector<double> mkk12 = get_mkk_hist( fname12, "mkk_12_cp", &hist[1] );
   TH1D* h09 = hist[0];
   TH1D* h12 = hist[1];

   int n09 = mkk09.size();
   int n12 = mkk12.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Dat(dr, n09+n12);
   for ( int i = 0; i < n09; ++i ) {
      Dat.Add(-mkk09[i]);
   }
   for ( int i = 0; i < n12; ++i ) {
      Dat.Add(mkk12[i]);
   }

   //-----------------------------------------------------------------------
   // Fit data
   double sl09 = -1.78;
   double sl12 = -1.87;
   // function MUST be normalized to 1 on the fit range
   auto Lfit = [sl09,sl12](double* x,double* p) -> double {
      double m = x[0];
      if ( m < 0 ) { // 2009
         double pp[] = { p[0],p[1],p[5],p[2],p[3],p[4],sl09 };
         double intfr = IntfrBWARGN(-m,pp,0);
         double bg = RevArgusN(-m,0.,sl09);
         return (1-p[7])*intfr + p[7]*bg;
      } else { // 2012
         double pp[] = { p[0],p[1],p[6],p[2],p[3],p[4],sl12 };
         double intfr = IntfrBWARGN(m,pp,0);
         double bg = RevArgusN(m,0.,sl12);
         return (1-p[7])*intfr + p[7]*bg;
//          return (1-p[8])*intfr + p[8]*bg;
      }
   };
   vector<string> par_name { "M#phi", "G#phi", "A", "F", "#vartheta",
                             "#sigma{09}", "#sigma{12}", "Abg" };

//    vector<double> par_ini {Mphi,Gphi,0.,1.0,0.8,1.5e-3,1.2e-3,0.};  // -neg
//    vector<double> par_ini {Mphi,Gphi,0.,1.0,-0.8,1.5e-3,1.2e-3,0.}; // -pos
   vector<double> par_ini {Mphi,Gphi,0.,0.8,0.5,1.5e-3,1.2e-3,0.01};  // neg
//    vector<double> par_ini {Mphi,Gphi,0.,0.8,-0.5,1.5e-3,1.2e-3,0.01}; // pos

   bool neg_pos = par_ini[4] > 0; // true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 WFit( *Ffit );
   fitter.SetFunction( WFit, false); // false == no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   fitter.Config().ParSettings(0).SetValue(1.01953); // like MC-sig
   fitter.Config().ParSettings(0).Fix();                    // Mphi
   fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3, Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();                    // Gphi
   fitter.Config().ParSettings(2).Fix();                    // A
   fitter.Config().ParSettings(3).SetLimits(0.01, 10.);     // F
   fitter.Config().ParSettings(4).SetLimits(-M_PI, M_PI);   // vartheta
   fitter.Config().ParSettings(5).SetLimits(0.2e-3, 2.e-3); // sigma09
   fitter.Config().ParSettings(6).SetLimits(0.2e-3, 2.e-3); // sigma12
   fitter.Config().ParSettings(7).SetLimits(0., 1.);        // Abg

//    fitter.Config().ParSettings(4).SetValue(0.); // vartheta
//    fitter.Config().ParSettings(4).Fix();

   // == Fit
   fitter.LikelihoodFit( Dat, false ); // true == extended likelihood fit
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
//    res.PrintCovMatrix(cout); // print error matrix and correlations

   double Lmin = res.MinFcnValue();
   printf(" Lmin= %.6f\n", Lmin );
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   vector<string> names { "N(KK)", "Nphi", "Nnonphi", "Nifr"  };
   vector<ValErr> Nos;
   cout << " n09= " << n09 << "  n12= " << n12 << endl;
   vector<double> norm { double(n09), double(n12) };
   Nos.reserve(8);
   for ( int idx = 0; idx <= 3; ++idx ) {
      vector<ValErr> Vinteg = CombIntErAr( res,sl09,sl12,idx );
      for ( int j = 0; j < 2 ; ++j ) {
         const auto& Integ = Vinteg[j];
//          printf("Int[%i,%i] = %s\n",idx,j,Integ.prt(".4f")); // DEBUG
         double Num = norm[j] * Integ.val;
         double Err = fabs(Num) *
            sqrt( 1./norm[j] + SQ(Integ.err/Integ.val) );
         Nos.push_back( ValErr {Num,Err} );
      }
      printf("%s09 = %s %s12 = %s\n",
            names[idx].c_str(),Nos[2*idx].prt(".1f"),
            names[idx].c_str(),Nos[2*idx+1].prt(".1f") );
   }
   double Nphi09 = Nos[2].val, err_Nphi09 = Nos[2].err;
   double Nphi12 = Nos[3].val, err_Nphi12 = Nos[3].err;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof09 = [Fpar,sl09](double x) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[5],Fpar[2],Fpar[3],Fpar[4],sl09};
      double intfr = IntfrBWARGN(x,pp,0);
      double bg = RevArgusN(x,0.,sl09);
      return (1-Fpar[7])*intfr + Fpar[7]*bg;
   };
   rmath_fun< decltype(Lgof09) > ftest09(Lgof09);
   ROOT::Math::GoFTest* goftest09 = new ROOT::Math::GoFTest(
         mkk09.size(),mkk09.data(),ftest09,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS09 = goftest09 -> KolmogorovSmirnovTest();
   cout << " pvalueKS09= " << pvalueKS09 << endl;

   auto Lgof12 = [Fpar,sl12](double x) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[6],Fpar[2],Fpar[3],Fpar[4],sl12};
      double intfr = IntfrBWARGN(x,pp,0);
      double bg = RevArgusN(x,0.,sl12);
      return (1-Fpar[7])*intfr + Fpar[7]*bg;
//       return (1-Fpar[8])*intfr + Fpar[8]*bg;
   };
   rmath_fun< decltype(Lgof12) > ftest12(Lgof12);
   ROOT::Math::GoFTest* goftest12 = new ROOT::Math::GoFTest(
         mkk12.size(),mkk12.data(),ftest12,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS12 = goftest12 -> KolmogorovSmirnovTest();
   cout << " pvalueKS12= " << pvalueKS12 << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   auto Ldr09 = [n09,Fpar,sl09](double* x,double* p) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[5],Fpar[2],Fpar[3],Fpar[4],sl09};
      return bW*n09 * (1-Fpar[7])*IntfrBWARGN( x[0],pp,int(p[0]) );
   };
   TF1* f09 = new TF1("f09", Ldr09, dL, dU, 1);
   f09 -> SetNpx(500);

   auto Ldr09_bg = [n09,Fpar,sl09](double* x,double* p) -> double {
      return bW*n09 * Fpar[7]*RevArgusN(x[0],0.,sl09);
   };
   TF1* f09_bg = new TF1("f09_bg", Ldr09_bg, dL, dU, 0);

   auto Ldr12 = [n12,Fpar,sl12](double* x,double* p) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[6],Fpar[2],Fpar[3],Fpar[4],sl12};
      return bW*n12 * (1-Fpar[7])*IntfrBWARGN( x[0],pp,int(p[0]) );
   };
   TF1* f12 = new TF1("f12", Ldr12, dL, dU, 1);
   f12 -> SetNpx(500);

   auto Ldr12_bg = [n12,Fpar,sl12](double* x,double* p) -> double {
      return bW*n12 * Fpar[7]*RevArgusN(x[0],0.,sl12);
   };
   TF1* f12_bg = new TF1("f12_bg", Ldr12_bg, bL, dU, 0);

   // find max/min to draw
   f09 -> SetParameter(0, 1); // BW
   int mb09 = h09 -> GetMaximumBin();
   double hmax09 = h09 -> GetBinContent(mb09) + h09 -> GetBinError(mb09);
   double fmax09 = f09 -> GetMaximum( 1.01, 1.03, 1e-6);
   hmax09 = floor(1.1 * max( hmax09, fmax09 ));
   if ( hmax09 > 0 ) {
      h09 -> SetMaximum(hmax09);
   }

   f09 -> SetParameter(0, 3); // interference
   double hmin09 = f09 -> GetMinimum( 1.01, 1.03, 1e-6);
   hmin09 = floor(1.15 * min( 0., hmin09 ));
   if ( hmin09 < 0 ) {
      h09 -> SetMinimum(hmin09);
//       h09->SetMinimum(-15.);
   }

   f12 -> SetParameter(0, 1); // BW
   int mb12 = h12 -> GetMaximumBin();
   double hmax12 = h12 -> GetBinContent(mb12) + h12 -> GetBinError(mb12);
   double fmax12 = f12 -> GetMaximum( 1.01, 1.03, 1e-6);
   hmax12 = floor(1.1 * max( hmax12, fmax12 ));
   if ( hmax12 > 0 ) {
      h12 -> SetMaximum(hmax12);
   }

   f12 -> SetParameter(0, 3); // interference
   double hmin12 = f12 -> GetMinimum( 1.01, 1.03, 1e-6);
   hmin12 = floor(1.15 * min( 0., hmin12 ));
   if ( hmin12 < 0 ) {
      h12 -> SetMinimum(hmin12);
//       h12 -> SetMinimum(-35.);
   }

   //-----------------------------------------------------------------------
   // Draw results
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000);
   c1 -> Divide(1,2);

   c1 -> cd(1);
   gPad -> SetGrid();

   SetHstFace(h09);
   h09 -> GetYaxis() -> SetTitleOffset(1.);
   h09 -> SetLineWidth(2);
   h09 -> SetLineColor(kBlack);
   h09 -> SetMarkerStyle(20);

   h09->Draw("EP");

   TLegend* leg = new TLegend(0.51,0.45,0.89,0.89);
   leg -> SetHeader("#bf{2009(top)   2012(bottom)}","C");
   leg -> AddEntry( h09,"Data","LEP" );

   f09 -> SetParameter(0, 0); // Sum
   f09 -> SetLineWidth(2);
   f09 -> SetLineColor(kRed+1);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Combined fit", "L" );

   f09 -> SetParameter(0, 1); // BW
   f09 -> SetLineWidth(1);
   f09 -> SetLineStyle(kDashed);
   f09 -> SetLineColor(kGreen+2);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Breit-Wigner #phi#eta", "L");

   f09 -> SetParameter(0, 2); // Argus
   f09 -> SetLineColor(kBlue);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Non-#phi KK#eta", "L");

   f09 -> SetParameter(0, 3); // interference
   f09 -> SetLineColor(kMagenta+1);
   f09 -> DrawCopy("SAME");
   leg -> AddEntry( f09 -> Clone(), "Interference", "L");

   f09_bg -> SetLineWidth(2); // non-interfering bkg
   f09_bg -> SetLineColor(kBlue+2);
   f09_bg -> DrawCopy("SAME");
   leg -> AddEntry( f09_bg, "Non-interfering non-#phi KK#eta", "L" );

   leg -> Draw();
   gPad -> RedrawAxis();

   c1 -> cd(2);
   gPad -> SetGrid();
   SetHstFace(h12);
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
   f12 -> SetLineWidth(1);
   f12 -> SetLineStyle(kDashed);
   f12 -> SetLineColor(kGreen+2);
   f12 -> DrawCopy("SAME");

   f12 -> SetParameter(0, 2); // Argus
   f12 -> SetLineColor(kBlue);
   f12 -> DrawCopy("SAME");

   f12 -> SetParameter(0, 3); // interference
   f12 -> SetLineColor(kMagenta+1);
   f12 -> DrawCopy("SAME");

   f12_bg -> SetLineWidth(2); // non-interfering bkg
   f12_bg -> SetLineColor(kBlue+2);
   f12_bg -> DrawCopy("SAME");

   TPaveText* pt09 = new TPaveText(0.51,0.82,0.89,0.99,"NDC");
   pt09 -> SetTextAlign(12);
   pt09 -> SetTextFont(42);
   pt09 -> AddText( Form("#it{p-value(2009) = %.3f}", pvalueKS09) );
   pt09 -> AddText( Form("N_{#phi}(2009) = %.1f #pm %.1f",
                    Nphi09,err_Nphi09) );
   pt09 -> AddText( Form("#sigma(2009)= %s MeV", PE.Eform(5,".2f",1e3)) );
   pt09 -> Draw();

   TPaveText* pt12 = new TPaveText(0.51,0.65,0.89,0.82,"NDC");
   pt12 -> SetTextAlign(12);
   pt12 -> SetTextFont(42);
   pt12 -> AddText( Form("#it{p-value(2012) = %.3f}", pvalueKS12) );
   pt12 -> AddText( Form("N_{#phi}(2012) = %.1f #pm %.1f",
                    Nphi12,err_Nphi12) );
   pt12 -> AddText( Form("#sigma(2012)= %s MeV", PE.Eform(6,".2f",1e3)) );
   pt12 -> Draw();

   TPaveText* pt = new TPaveText(0.51,0.33,0.89,0.65,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
   pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt -> AddText( Form("a= %s",PE.Eform(2,".2f")) );
   pt -> AddText( Form("F= %s",PE.Eform(3,".2f")) );
   pt -> AddText( Form("#vartheta= %s",PE.Eform(4,".2f")) );
   pt -> AddText( Form("f_{s} = %s", PE.Eform(7,".4f")) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      pdf += (neg_pos ? "_n.pdf" : "_p.pdf");
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void combine_Intfr_Ar_scan(string fname09, string fname12, string pdf="") {
//-------------------------------------------------------------------------
   // Get un-binned data
   TH1D* hist[2];
   vector<double> mkk09 = get_mkk_hist( fname09, "mkk_09_cp", &hist[0] );
   vector<double> mkk12 = get_mkk_hist( fname12, "mkk_12_cp", &hist[1] );
   TH1D* h09 = hist[0];
   TH1D* h12 = hist[1];

   int n09 = mkk09.size();
   int n12 = mkk12.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Dat(dr, n09+n12);
   for ( int i = 0; i < n09; ++i ) {
      Dat.Add(-mkk09[i]);
   }
   for ( int i = 0; i < n12; ++i ) {
      Dat.Add(mkk12[i]);
   }

   //-----------------------------------------------------------------------
   // Fit data
   double sl09 = -1.78;
   double sl12 = -1.87;
   // function MUST be normalized to 1 on the fit range
   auto Lfit = [sl09,sl12](double* x,double* p) -> double {
      double m = x[0];
      if ( m < 0 ) { // 2009
         double pp[] = { p[0],p[1],p[5],p[2],p[3],p[4],sl09 };
         double intfr = IntfrBWARGN(-m,pp,0);
         double bg = RevArgusN(-m,0.,sl09);
         return (1-p[7])*intfr + p[7]*bg;
      } else { // 2012
         double pp[] = { p[0],p[1],p[6],p[2],p[3],p[4],sl12 };
         double intfr = IntfrBWARGN(m,pp,0);
         double bg = RevArgusN(m,0.,sl12);
         return (1-p[7])*intfr + p[7]*bg;
//          return (1-p[8])*intfr + p[8]*bg;
      }
   };
   vector<string> par_name { "M#phi", "G#phi", "A", "F", "#vartheta",
                             "#sigma{09}", "#sigma{12}", "Abg" };

//    vector<double> par_ini {Mphi,Gphi,0.,0.8,0.5,1.5e-3,1.2e-3,0.01};// 30deg
   vector<double> par_ini {Mphi,Gphi,0.,1.4,1.3,1.4e-3,1.1e-3,0.01};// -75deg

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
//    fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 WFit( *Ffit );
   fitter.SetFunction( WFit, false); // false == no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   fitter.Config().ParSettings(0).SetValue(1.01953); // like MC-sig
   fitter.Config().ParSettings(0).Fix();                    // Mphi
   fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3, Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();                    // Gphi
   fitter.Config().ParSettings(2).Fix();                    // A
   fitter.Config().ParSettings(3).SetLimits(0.01, 10.);     // F
//    fitter.Config().ParSettings(4).SetLimits(-M_PI, M_PI);   // vartheta
   fitter.Config().ParSettings(5).SetLimits(0.2e-3, 2.e-3); // sigma09
   fitter.Config().ParSettings(6).SetLimits(0.2e-3, 2.e-3); // sigma12
   fitter.Config().ParSettings(7).SetLimits(0., 1.);        // Abg

   //-----------------------------------------------------------------------
   // Functions to draw
   vector<double> Fpar;
   auto Ldr09 = [n09,&Fpar,sl09](double* x,double* p) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[5],Fpar[2],Fpar[3],Fpar[4],sl09};
      return bW*n09 * (1-Fpar[7])*IntfrBWARGN( x[0],pp,int(p[0]) );
   };
   TF1* f09 = new TF1("f09", Ldr09, dL, dU, 1);
   f09 -> SetNpx(500);

   auto Ldr09_bg = [n09,&Fpar,sl09](double* x,double* p) -> double {
      return bW*n09 * Fpar[7]*RevArgusN(x[0],0.,sl09);
   };
   TF1* f09_bg = new TF1("f09_bg", Ldr09_bg, dL, dU, 0);

   auto Ldr12 = [n12,&Fpar,sl12](double* x,double* p) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[6],Fpar[2],Fpar[3],Fpar[4],sl12};
      return bW*n12 * (1-Fpar[7])*IntfrBWARGN( x[0],pp,int(p[0]) );
   };
   TF1* f12 = new TF1("f12", Ldr12, dL, dU, 1);
   f12 -> SetNpx(500);

   auto Ldr12_bg = [n12,&Fpar,sl12](double* x,double* p) -> double {
      return bW*n12 * Fpar[7]*RevArgusN(x[0],0.,sl12);
   };
   TF1* f12_bg = new TF1("f12_bg", Ldr12_bg, bL, dU, 0);

   //-----------------------------------------------------------------------
   // loop for angles (degrees)
   int angmin = -75, angmax = 75, angstep = 5;
//    int angmin = 30, angmax = 95, angstep = 5;
//    int angmin = 30, angmax = -95, angstep = -5;

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000);
   c1 -> Divide(1,2);
   c1 -> cd();
   c1 -> Print((pdf+"[").c_str()); // open pdf-file

   SetHstFace(h09);
   h09 -> GetYaxis() -> SetTitleOffset(1.);
   h09 -> SetLineWidth(2);
   h09 -> SetLineColor(kBlack);
   h09 -> SetMarkerStyle(20);
   h09 -> SetMinimum(-15.);

   SetHstFace(h12);
   h12 -> GetYaxis() -> SetTitleOffset(1.);
   h12 -> SetLineWidth(2);
   h12 -> SetLineColor(kBlack);
   h12 -> SetMarkerStyle(20);
   h12 -> SetMinimum(-35.);

   vector<double> angl, Lmin;
   for( int ang = angmin; ang != angmax; ang+=angstep ) {

      fitter.Config().ParSettings(4).SetValue(M_PI*ang/180.); // vartheta
      fitter.Config().ParSettings(4).Fix();

      // == Fit
      fitter.LikelihoodFit( Dat, false ); // true == extended likelihood fit
//       fitter.CalculateMinosErrors();

      ROOT::Fit::FitResult res = fitter.Result();
//       res.Print(cout);
//       res.PrintCovMatrix(cout); // print error matrix and correlations
      ParAndErr PE(res);
      Fpar = PE.Fpar;

      double lmin = res.MinFcnValue();
      printf(" ang= %i degree =>  Lmin= %.6f\n", ang,lmin );
      angl.push_back(ang);
      Lmin.push_back(lmin);

      //--------------------------------------------------------------------
      // Draw results

      c1 -> cd(1);
      gPad -> SetGrid();
      h09->Draw("EP");

      TLegend* leg = new TLegend(0.51,0.45,0.89,0.89);
      leg -> SetHeader("#bf{2009(top)   2012(bottom)}","C");
      leg -> AddEntry( h09,"Data","LEP" );

      f09 -> SetParameter(0, 0); // Sum
      f09 -> SetLineWidth(2);
      f09 -> SetLineStyle(kSolid);
      f09 -> SetLineColor(kRed+1);
      f09 -> DrawCopy("SAME");
      leg -> AddEntry( f09 -> Clone(), "Combined fit", "L" );

      f09 -> SetParameter(0, 1); // BW
      f09 -> SetLineWidth(1);
      f09 -> SetLineStyle(kDashed);
      f09 -> SetLineColor(kGreen+2);
      f09 -> DrawCopy("SAME");
      leg -> AddEntry( f09 -> Clone(), "Breit-Wigner #phi#eta", "L");

      f09 -> SetParameter(0, 2); // Argus
      f09 -> SetLineColor(kBlue);
      f09 -> DrawCopy("SAME");
      leg -> AddEntry( f09 -> Clone(), "Non-#phi KK#eta", "L");

      f09 -> SetParameter(0, 3); // interference
      f09 -> SetLineColor(kMagenta+1);
      f09 -> DrawCopy("SAME");
      leg -> AddEntry( f09 -> Clone(), "Interference", "L");

      f09_bg -> SetLineWidth(2); // non-interfering bkg
      f09_bg -> SetLineStyle(kSolid);
      f09_bg -> SetLineColor(kBlue+2);
      f09_bg -> DrawCopy("SAME");
      leg -> AddEntry( f09_bg, "Non-interfering non-#phi KK#eta", "L" );

      leg -> Draw();
      gPad -> RedrawAxis();

      c1 -> cd(2);
      gPad -> SetGrid();
      h12 -> Draw("EP");

      f12 -> SetParameter(0, 0); // Sum
      f12 -> SetLineWidth(2);
      f12 -> SetLineStyle(kSolid);
      f12 -> SetLineColor(kRed+1);
      f12 -> DrawCopy("SAME");

      f12 -> SetParameter(0, 1); // BW
      f12 -> SetLineWidth(1);
      f12 -> SetLineStyle(kDashed);
      f12 -> SetLineColor(kGreen+2);
      f12 -> DrawCopy("SAME");

      f12 -> SetParameter(0, 2); // Argus
      f12 -> SetLineColor(kBlue);
      f12 -> DrawCopy("SAME");

      f12 -> SetParameter(0, 3); // interference
      f12 -> SetLineColor(kMagenta+1);
      f12 -> DrawCopy("SAME");

      f12_bg -> SetLineWidth(2); // non-interfering bkg
      f12_bg -> SetLineStyle(kSolid);
      f12_bg -> SetLineColor(kBlue+2);
      f12_bg -> DrawCopy("SAME");

      TPaveText* pt = new TPaveText(0.51,0.34,0.89,0.99,"NDC");
      pt -> SetTextAlign(12);
      pt -> SetTextFont(42);
      pt -> AddText( Form("#vartheta= %s deg",PE.Eform(4,".0f",180./M_PI)) );
      pt -> AddText( Form("#it{L_{min}} = %.1f",lmin) );
      pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
      pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
      pt -> AddText( Form("a= %s",PE.Eform(2,".2f")) );
      pt -> AddText( Form("F= %s",PE.Eform(3,".2f")) );
      pt -> AddText( Form("f_{s} = %s", PE.Eform(7,".4f")) );
      pt -> AddText( Form("#sigma(2009)= %s MeV", PE.Eform(5,".2f",1e3)) );
      pt -> AddText( Form("#sigma(2012)= %s MeV", PE.Eform(6,".2f",1e3)) );
      pt -> Draw();

      gPad -> RedrawAxis();
      c1 -> Update();
      c1 -> Print(pdf.c_str()); // add to pdf-file
   }
   //-----------------------------------------------------------------------
   // final draw

   c1 -> cd();
   c1 -> Clear();
   c1 -> SetCanvasSize(800,800); // resize
   c1 -> cd();
   gPad -> SetGrid();

   double minL = *(min_element(Lmin.begin(),Lmin.end()));
   for( auto &l : Lmin ) {
      l -= minL;
   }
   int nch = angl.size();
   auto gr = new TGraph( nch, angl.data(), Lmin.data() );
   gr -> SetTitle(";#vartheta, degrees;#it{-2log(L/L_{max})}");
   gr -> GetYaxis() -> SetMaxDigits(3);
   gr -> GetYaxis()->SetTitleOffset(1.2);
   gr -> SetMarkerColor(kBlue);
   gr -> SetMarkerStyle(21);
   gr -> SetLineWidth(2);

   gr -> Draw("APL");

   gPad -> RedrawAxis();
   c1 -> Update();
   c1 -> Print(pdf.c_str()); // add to pdf-file
   c1 -> Print((pdf+"]").c_str()); // close pdf-file

}

// {{{1 MAIN:
//-------------------------------------------------------------------------
void mass_kk_fit() {
//-------------------------------------------------------------------------
   gROOT -> Reset();
   gStyle -> SetOptStat(0);
   gStyle -> SetStatFont(62);
   gStyle -> SetLegendFont(42);

   if ( DEBUG ) {
      // BOOK HISTOGRAMMS
      htt[0] = new TH1D("eBWG","log10(BWG conv. error)",100,-6.,-1.);
      htt[1] = new TH1D("nBW","BreitWigner() norm",1000,0.,2.);
      htt[2] = new TH1D("eBW","log10(BW norm. error)",100,-6.,-1.);
      htt[3] = new TH1D("nAr","RevArgus() norm",1000,0.,10.);
      htt[4] = new TH1D("eAr","log10(RevArgus norm. error)",100,-6.,-1.);
      htt[5] = new TH1D("eInG","log10(Intfr. conv. error)",100,-6.,-1.);
   }

   vector<string> fnames = {
      "data_09psip_all.root",
      "data_12psip_all.root",
      "data_3650_all.root",
      "mcinc_09psip_all.root",
      "mcinc_12psip_all.root",
      "mcsig_kkmc_09.root",
      "mcsig_kkmc_12.root",
      "mckketa_kkmc_09.root",
      "mckketa_kkmc_12.root"
   };
   vector<string> titles = {
      "Data 2009",
      "Data 2012",
      "E=3.65 GeV",
      "MC incl. 2009",
      "MC incl. 2012",
      "MC signal #phi#eta 2009",
      "MC signal #phi#eta 2012",
      "MC non-#phi KK#eta 2009",
      "MC non-#phi KK#eta 2012"
   };
   // Numbers of items for date:
   int id = 0, ii = 3, is = 5, ib = 7; // 2009
//    int id = 1, ii = 5, is = 6, ib = 8; // 2012

   // = define GSL error handler which does nothing =
   gsl_set_error_handler_off();

   // = set handler =
//    gsl_set_error_handler( my_handler );

   // set integrator: ROOT::Math::GSLIntegrator adaptive method (QAG)
   ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Adaptive");

// ------------- signal ---------------
//    test_BreitWigner();
//    string pdf = string("mkk")+((id==0) ? "09" : "12")+"_sig_log.pdf";
//    sig_fit(fnames.at(is),titles.at(is),pdf); // fit phi-eta by BW x Gauss

//    sig_fit_mc("prod-10/mcsig_kkmc_09.root","mcmkk09_sig.pdf");
//    sig_fit_mc("prod-10/mcsig_kkmc_12.root","mcmkk12_sig.pdf");

// ------------- background ---------------
//    test_Agrus();
//    test_RevAgrus();

//    string pdf = string("mkk")+((id==0) ? "09" : "12")+"_fitbg.pdf";
//    bkg_fit(fnames.at(ib),titles.at(ib),pdf); // fit kketa by Argus

//    bkg_fit("mcrhoeta_kkmc_12.root","MC #rho#eta 2012",
//            "mkk_rhoeta12_argus.pdf");

   // sum 09+12
//    string pdf = string("mkk_sum_fitbg.pdf");
//    bkg_fit2(fnames[7],fnames[8],pdf); // fit kketa by Argus

// ------------- data: no interference ---------------
//    int Mphix = 1; // 1 - mphi is fixed at the PDG
//    string pdf = string("mkk") + ((id==0) ? "09" : "12") + "_fit"
//       + to_string(Mphix) + ".pdf";
//    data_fit(fnames.at(id),titles.at(id),2009+id*3,Mphix,pdf);

// ------------- data: interference ---------------
//    test_Intfr();
//    string pdf = string("mkk")+((id==0) ? "09" : "12")+"_ifr";
//    data_Intfr_fit(fnames.at(id),titles.at(id),pdf);

   // combined with sb
   string pdf = string("mkk")+((id==0) ? "09" : "12")+"_ifr_SB";
   data_Intfr_sb(fnames.at(id),titles.at(id),pdf);

// ------------- data: common fit ---------------
//    combine_Intfr(fnames[0],fnames[1],"mkk_cf_std");
//    combine_Intfr_bg(fnames[0],fnames[1],"mkk_cf_bkg");
//    combine_Intfr_sb(fnames[0],fnames[1],"mkk_cf_SB");

// ------------- data: interference + incoherent sum ---------------
//  DO NOT USE IT ---> see memo_intfr_plus.tex
//    string pdf = string("mkk")+((id==0) ? "09" : "12")+"_ifrbg";
//    data_Intfr_Ar(fnames.at(id),titles.at(id),pdf);

//    combine_Intfr_Ar(fnames[0],fnames[1],"mkk_cf_Abg");
//    combine_Intfr_Ar_scan(fnames[0],fnames[1],"mkk_cf_Abg_scan.pdf");

// -----------------------------------------------
}