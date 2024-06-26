// unbinned LH fit of M(K+K-) distributions
// after the cuts (see cuts.h): Mrec + chi^2(4C) + Mgg
// Formulas for:
// *) Breit-Wigner for psi -> K+K- (psi from J/Psi -> phi eta)
// *) convolution of BW with the Gauss distribution
// *) ARGUS function (and "reverse" Argus)
// *) interference
//    -> mkkYY_fit??.pdf
//    -> mkk_cf_???.pdf (combined fit)

#include <functional>        // for std::function
#include "gsl/gsl_errno.h"   // GSL error handler

#include "masses.h"

// {{{1 constants, TODO: struct ?
//--------------------------------------------------------------------
// the meson radial parameter in Blatt-Weisskopf form factor
static const double R_BW = 3.; // GeV^{-1}; +/-1 uncertainties study

static const double dL = 2*Mk; // the left cutoff = 0.987354

// binning for Breit-Wigner: bin width 1.0 MeV
static const int Nbins = 100;
static const double bL = 0.98; // first bin < dL
static const double dU = 1.08; // MUST BE < 1.0835 !!!
static const double bW = (dU-bL)/Nbins; // bin width

// {{{1 handler for GSL exceptions
//--------------------------------------------------------------------
static char gsl_integral_func[91] {"not used"};
// histograms for debug
static bool DEBUG = false;
static TH1D* htt[10];

//--------------------------------------------------------------------
void my_handler(const char* reason, const char* file, int line,
      int gsl_errno) {
//--------------------------------------------------------------------
   printf("-- GSL EXCEPTION: %s in func %s\n",
         reason,gsl_integral_func);
   return;
}

// {{{1 helper functions
//--------------------------------------------------------------------
// adapter class to use lambda functions with closures in ROOT::MATH
// Using:  rmath_fun< decltype(Lambda) > Functor(Lambda);
//--------------------------------------------------------------------
template< typename F >
class rmath_fun: public ROOT::Math::IBaseFunctionOneDim {
//--------------------------------------------------------------------
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

   // if the upper and lower errors differ no more than the
   // 'threshold', the maximum of them is used as a symmetric error
   double threshold = 0.1;

   // ctor
   //-----------------------------------------------------------------
   ParAndErr(const ROOT::Fit::FitResult& res,double thr = 0.1) {
   //-----------------------------------------------------------------
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
   const char* Eform (int i, string fm, double X = 1) {
   //-----------------------------------------------------------------
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
            + "^{#lower[0.15]{#kern[0.15]{#plus}#kern[0.05]{%"
            + fm + "}}" + "}"
            + "_{#lower[-0.15]{#kern[0.15]{#minus}#kern[0.05]{%"
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
   double Middle (int i) {
   //-----------------------------------------------------------------
      return Fpar[i] + 0.5*(Uerr[i]-Lerr[i]);
   }
};

//--------------------------------------------------------------------
struct ValErr {
   double val = 0.;
   double err = 0.;
   const char* prt(string fm) const {
      string fmt = "%" + fm + " +/- %" + fm;
      return Form(fmt.c_str(),val,err);
   }
};

//--------------------------------------------------------------------
// adapter function to spead up the Kolmogorov-Smirnov test
// linear extrapolation of Fun on a uniform grid of n points in the
// region [dL,dU] is used
//--------------------------------------------------------------------
double FastKSTest(const function<double (double)>& Fun,
      const vector<double>& Mkk, size_t NMkk = 0) {
//--------------------------------------------------------------------
   int n = 10000;
   vector<double> vec(n);
   double dx = (dU-dL)/(n-1);
   for ( int i = 0; i < n; ++i ) {
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
   double pvalueKS = goftest->KolmogorovSmirnovTest();
   return pvalueKS;
   // double pvalueAD = gofcr->AndersonDarlingTest();
   // cout << " pvalueAD= " << pvalueAD << endl;
}

//--------------------------------------------------------------------
constexpr double SQ(double x) {
//--------------------------------------------------------------------
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

// {{{1 data processing
//--------------------------------------------------------------------
vector<double> get_mkk_hist( string fname, string hname, TH1D* hst[],
      int type = 0 ) {
//--------------------------------------------------------------------
#include "cuts.h"
   bool isMC = (fname.find("mcsig") != string::npos);

   // name of folder with root files
   // static string dir("prod-12/");
   string dir("prod_v709/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   TCut c_here = c_Mrec+c_chi2;
   // c_here += TCut("chsq3g>ch2");// with the search for third gamma
   c_here += TCut(Form("%f<=Mkk&&Mkk<%f",dL,dU));
   if ( type == 0 ) {        // central part
      c_here += c_cpgg;
   } else if ( type == 1 ) { // side-band
      c_here += c_sbgg;
   }
   if ( isMC ) {
      cout << " THIS IS MONTE CARLO EVENTS" << endl;
      c_here += c_MCmkk;
   }

   int n = a4c->Draw("Mkk", c_here, "goff");
   // cout << " n= " << n << endl;
   double* buffer = a4c->GetVal(0);
   vector<double> mkk(buffer, buffer+n);

   string title(";M^{ inv}_{ K^{#plus}K^{#minus }}, GeV/c^{2}");
   int ibW = int(bW*1e5+0.1);
   title += string(";Entries / ") +
      ( (ibW%10 == 0) ? string(Form("%.0f MeV/c^{2}",bW*1e3))
        : string(Form("%.2f MeV/c^{2}",bW*1e3)) );
   hst[0] = new TH1D(hname.c_str(),title.c_str(),Nbins,bL,dU);

   hst[0]->Sumw2(true);
   for ( const auto& mk : mkk ) {
      hst[0]->Fill(mk);
   }
   return mkk;
}

//--------------------------------------------------------------------
TH1D* plot_Mkk( string fname, string hname, int type = 0 ) {
//--------------------------------------------------------------------
   TH1D* hist[1];
   get_mkk_hist( fname, hname, hist, type );
   return hist[0];
}

//--------------------------------------------------------------------
TH1D* get_mcMkk(string fname) {
//--------------------------------------------------------------------
   // name of folder with root files
   // static string dir("prod-12/");
   string dir("prod_v709/");

   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( !froot ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TH1D* mcmkk = (TH1D*) gDirectory->Get("mc_Mkk");
   if( !mcmkk ) {
      cerr << "can not find mc_Mkk" << endl;
      exit(EXIT_FAILURE);
   }

   string title(";M^{ inv}_{ K^{#plus}K^{#minus }}, GeV/c^{2}");
   title += string(";Entries / 1 MeV/c^{2}");
   mcmkk->SetTitle(title.c_str());
   // mcmkk->Sumw2(true);

   return mcmkk;
}

// {{{1 Breigt Wigner for phi -> KK
//--------------------------------------------------------------------
double BreitWigner(double m, double mphi, double gphi, double R=R_BW){
//--------------------------------------------------------------------
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
   double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k2; // == m*G(m)

   return (kap * SQ(fm)) / ( SQ(m2-mphi2) + SQ(GM) );
}

//--------------------------------------------------------------------
double BreitWignerGauss( double m,
                         double mphi, double gphi, // BW parameters
                         double sigma,             // Gauss
                         double slope              // eff(m)
                       ) {
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

//--------------------------------------------------------------------
double BreitWignerGaussN( double m,
                          double mphi, double gphi, double sigma,
                          double slope
                        ) {
//--------------------------------------------------------------------
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

//--------------------------------------------------------------------
void test_BreitWigner() {
//--------------------------------------------------------------------

   // 1) BreitWigner (internally normalized to 1)
   auto Lbw = [](const double* x,const double* p) -> double {
      return p[0]*BreitWigner(x[0],p[1],p[2],p[3]);
   };
   TF1* bw = new TF1("bw", Lbw, dL, dU, 4);
   bw->SetParNames("Norm","Mphi","Gphi","Rff");
   bw->SetParameters(1., Mphi, Gphi,R_BW);
   bw->SetLineWidth(1);
   bw->SetLineColor(kBlue);
   bw->SetNpx(500);

   double norm = 1. / bw->Integral(dL,dU,1e-8);
   printf("norm = %.7f\n",norm);
   bw->SetParameter(0, norm );

   // test of stability
   // for ( int i = 0; i < 100; i++ ) { // mphi = Mphi +/- 1MeV
      // double mphi = Mphi + (-1. + 0.02*i)*1e-3;
      // for ( int j = 0; j < 100; j++ ) { // sigma from 1.1 to 1.3MeV
         // double sigma = (1.1 + 0.02*j)*1e-3;
         // BreitWignerGaussN(mphi-4*sigma,mphi,Gphi,sigma);
      // }
   // }

   // 2) BreitWignerGaussN
   auto Lbwgn = [](const double* x,const double* p) -> double {
      return p[0]*BreitWignerGaussN(x[0],p[1],p[2],p[3],p[4]);
   };
   TF1* bwgn = new TF1("bwgn", Lbwgn, dL, dU, 5);
   bwgn->SetParNames("Norm","Mphi","Gphi","Sigma","Slope");
   bwgn->SetParameters(1., Mphi, Gphi, 1.2e-3, 0.);
   bwgn->SetLineWidth(2);
   bwgn->SetLineColor(kRed);
   bwgn->SetLineStyle(kDashed);
   bwgn->SetNpx(500);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(true);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader("Breit-Wigner (M_{#phi}, #Gamma_{#phi})","C");
   leg->AddEntry(bw->Clone(), "BW, Rff=3GeV^{-1}", "L");

   bw->DrawCopy("L");
   c1->Update();

   bw->Update();
   bw->SetParameters(1., Mphi, Gphi, 6.); // Blatt-Weisskopf x2
   bw->SetLineColor(kGreen);
   bw->DrawCopy("LSAME");
   leg->AddEntry(bw->Clone(), "BW, Rff=6GeV^{-1}", "L");

   bw->Update();
   bw->SetParameters(1., Mphi, Gphi, 1.); // Blatt-Weisskopf ff
   bw->SetLineColor(kRed);
   bw->DrawCopy("LSAME");
   leg->AddEntry(bw->Clone(), "BW, Rff=1GeV^{-1}", "L");

   // bwgn->DrawCopy("SAME");
   // leg->AddEntry( bwgn->Clone(), Form("#sigma=%.2f sl=%.2f",
            // bwgn->GetParameter(3)*1e3,
            // bwgn->GetParameter(4) ), "L" );

   // bwgn->SetParameter(4, -1.9 );  // slope
   // bwgn->SetLineColor(kGreen);
   // bwgn->DrawCopy("SAME");
   // leg->AddEntry(bwgn->Clone(), Form("#sigma=%.2f sl=%.2f",
            // bwgn->GetParameter(3)*1e3,
            // bwgn->GetParameter(4) ), "L");

   leg->Draw();
   gPad->RedrawAxis();
   c1->Update();
}

//--------------------------------------------------------------------
void sig_fit_mc(int date, string pdf="") {
//--------------------------------------------------------------------
   string fname( Form("mcsig_kkmc_%02i.root",date%100) );
   TH1D* hst = get_mcMkk(fname);
   double bin_width = hst->GetBinWidth(1);

   // Fit signal MC(phi,eta) by Breit-Wigner
   auto Lbwg = [bin_width](const double* x,const double* p) {
      return bin_width*p[0]*BreitWigner(x[0],p[1],p[2],p[3]);
   };
   TF1* bwg = new TF1("bwg", Lbwg, 0.98, 2., 4);
   bwg->SetParNames("Norm","Mphi","Gphi","r");
   bwg->SetParameters(1.e5, Mphi, Gphi, R_BW);

   // bwg->FixParameter(1, Mphi);
   // bwg->FixParameter(2, Gphi);
   // bwg->FixParameter(3, R_BW);
   bwg->SetParLimits(3, 0.5, 5.);

   bwg->SetLineWidth(2);
   bwg->SetLineColor(kRed);
   bwg->SetNpx(500);

   // gStyle->SetOptFit(111); // do not print fixed parameters
   gStyle->SetOptFit(112); // print all parameters
   gStyle->SetStatFont(42);
   gStyle->SetFitFormat(".7g"); // or .6f
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   SetHstFace(hst);
   hst->SetMinimum(1.);
   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetTitleOffset(1.25);
   hst->SetLineWidth(3);
   hst->SetLineColor(kBlack);

   hst->Draw("EP");
   double maxMkk=1.0835;
   hst->Fit("bwg","E","",dL,maxMkk);

   TLine* lB = new TLine;
   lB->SetLineColor(kBlue+1);
   lB->SetLineWidth(2);
   lB->SetLineStyle(kDashed);
   lB->DrawLine(maxMkk,1,maxMkk,hst->GetBinContent(103));

   TF1* bwg1 = (TF1*) bwg->Clone();
   bwg1->SetLineStyle(kDashed);
   bwg1->DrawF1(maxMkk,1.125,"SAME");

   double genEv = hst->GetEntries();
   double corr = bwg->GetParameter(0) / genEv;
   double dcorr = bwg->GetParError(0) / genEv;
   cout << " genEv= " << genEv << " corr= " << corr << endl;

   double intg108 = bwg->Integral(dL,1.08,1e-8) / bin_width;
   double intg12  = bwg->Integral(dL,1.2,1e-8) / bin_width;
   double norm108 = bwg->GetParameter(0) / intg108;
   double norm12  = bwg->GetParameter(0) / intg12;
   cout << " norm108= " << norm108 << endl;
   cout << " norm12= " << norm12 << endl;

   TLegend* leg = new TLegend(0.22,0.20,0.62,0.40);
   leg->SetTextSize(0.03);
   leg->SetHeader( Form("Cut-off correction: %.3f#pm%.3f",corr,dcorr),
         "C" );
   leg->AddEntry(hst, Form("MC truth %i",date),"LEP");
   leg->AddEntry(bwg,"Breit-Wigner","L");
   leg->AddEntry(lB,"Cut-off","L");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

//--------------------------------------------------------------------
void sig_fit(int date, string pdf="") {
//--------------------------------------------------------------------
   // Get un-binned data
   TH1D* hist[1];
   string fname( Form("mcsig_kkmc_%02i.root",date%100) );
   vector<double> mkk = get_mkk_hist( fname, "mkk_phi", hist );
   TH1D* hst = hist[0];
   double norm = hst->Integral();

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Sig(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Sig.Add(mkk[i]);
   }
   cout << " norm= " << norm << " n= " << n << endl;

   //-----------------------------------------------------------------
   // Fit signal MC(phi,eta) by Breit-Wigner convoluted with Gauss
   // function MUST be normalized to 1 on the fit range
   // double slope = -1.8; // prod-12
   double slope = -2.0; // v709
   auto Lbwg = [slope](const double* x,const double* p) -> double {
      return BreitWignerGaussN(x[0],p[0],p[1],p[2],slope);
   };
   vector<string> par_name { "Mphi", "Gphi", "sigma" };
   vector<double> par_ini  {  Mphi,   Gphi,  1.2e-3  };

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* bwg = new TF1("bwg", Lbwg, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 fbwg( *bwg );
   fitter.SetFunction(fbwg,false); // false=no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(Mphi-0.01,Mphi+0.01);
//    fitter.Config().ParSettings(0).Fix();
//    fitter.Config().ParSettings(1).SetLimits(Gphi-1e-4,Gphi+1e-4);
//    fitter.Config().ParSettings(1).SetValue(4.247e-3);
   fitter.Config().ParSettings(1).Fix();
   fitter.Config().ParSettings(2).SetLimits(0.2e-3, 2.e-3);

   fitter.LikelihoodFit(Sig,false); // true=extended likelihood fit
   fitter.CalculateHessErrors(); // without MINOS here
//    fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   // use limited statistics (larger than the data)
   size_t max_mkk = 0;
   if ( date == 2009 ) {  // ~800ev
      max_mkk = 0; // all
   } else if ( date == 2012 ) { // ~2600ev
      max_mkk = 9000;
   } else if ( date == 2021 ) { // ~16 000ev
      max_mkk = 20000;
   }

   auto Fun = [Lbwg,Fpar](double m) -> double {
      return Lbwg(&m,Fpar.data());
   };
   double pvalueKS = FastKSTest(Fun,mkk,max_mkk);
   printf(" pvalueKS= %.3g (%zu)\n", pvalueKS,max_mkk);

   //-----------------------------------------------------------------
   // Functions to draw
   double Wn = bW*n;
   auto Lfit = [Wn,Lbwg,Fpar](double* x,double* p) -> double {
      return  Wn * Lbwg(x,Fpar.data());
   };
   TF1* ffit = new TF1("ffit", Lfit, dL, dU, 0);
   ffit->SetLineWidth(2);
   ffit->SetLineColor(kRed);
   ffit->SetNpx(500);

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   SetHstFace(hst);
   hst->SetMinimum(1.);
   hst->GetYaxis()->SetTitleOffset(1.25);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
   ffit->Draw("SAME");

   TLegend* leg = new TLegend(0.55,0.81,0.89,0.89);
   leg->AddEntry(hst,Form("MC signal #phi#eta %i",date),"LEP");
   leg->AddEntry(bwg,"(#it{BW}#times#it{E}) #otimes #it{Gauss}","L");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.55,0.63,0.89,0.80,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   // pt->AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt->AddText( Form("#it{p-value(KS) = %.3f}",pvalueKS) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt->AddText( Form("#sigma = %s MeV", PE.Eform(2,".2f",1e3)) );
   pt->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 Argus functions
//--------------------------------------------------------------------
double Argus(double x, double a) {
//--------------------------------------------------------------------
//ARGUS distribution: https://en.wikipedia.org/wiki/ARGUS_distribution
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

//--------------------------------------------------------------------
double RevArgus(double m, double A) {
//--------------------------------------------------------------------
// Argus function with left cutoff (L)
//      U is the upper limit of fitting
//      L,U - must be const for fit

   static const double L = 2*Mk;
   static const double U = Mjpsi - Meta;

   double x = (U-m)/(U-L); // transformation: (L,U) -> (1,0)
   return Argus(x,A)/(U-L);
}

//--------------------------------------------------------------------
double RevArgusN(double m, double A, double slope) {
//--------------------------------------------------------------------
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
      ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-8, 1000);
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

//--------------------------------------------------------------------
void test_Argus() {
//--------------------------------------------------------------------
   auto Largus = [](const double* x,const double* p) -> double {
//       return p[0]*Argus(x[0],p[1]);
      return p[0]*RevArgus(x[0],p[1]);
   };

//    TF1* fbkg = new TF1("fbkg", Largus, 0., 1., 2);
   TF1* fbkg = new TF1("fbkg", Largus, 0.98, 2., 2);
   fbkg->SetParNames("N","A");
   fbkg->SetParameters(1., 2.);
   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad->SetGrid();
   fbkg->DrawCopy();

   // test normalization:
   printf(" norm= %.15f\n",fbkg->Integral(2*Mk,Mjpsi-Meta));

   fbkg->SetParameters(1., 1.);
   fbkg->SetLineColor(kRed);
   fbkg->DrawCopy("SAME");
}

//--------------------------------------------------------------------
void test_RevArgus() {
//--------------------------------------------------------------------
   auto Lrargus = [](const double* x,const double* p) -> double {
      return p[0]*RevArgusN(x[0],p[1],p[2]);
   };

   TF1* fbkg = new TF1("fbkg", Lrargus, 0.98, dU, 3);
   fbkg->SetParNames("N","A","Sl");
   fbkg->SetParameters(1., 0.98, 0.);
   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);

   TLegend* leg = new TLegend(0.35,0.20,0.75,0.40);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad->SetGrid();
   fbkg->DrawCopy();
   leg->AddEntry(fbkg->Clone(), Form("RevArgus A=%.2f Sl=%.2f",
            fbkg->GetParameter(1),fbkg->GetParameter(2)),"L");

   // test normalization:
   printf(" RevArgusN norm= %.15f\n",fbkg->Integral( dL,dU ) );

   fbkg->SetParameter(1, 0. );  // A
   fbkg->SetParameter(2, -2.0 ); // Slope
   fbkg->SetLineColor(kRed);
   fbkg->SetLineStyle(kDashed);
   fbkg->DrawCopy("SAME");
   leg->AddEntry(fbkg->Clone(), Form("RevArgus A=%.2f Sl=%.2f",
            fbkg->GetParameter(1),fbkg->GetParameter(2)),"L");
   printf(" RevArgusN norm2= %.15f\n",fbkg->Integral( dL,dU ) );

   fbkg->SetParameter(1, 0.98 );  // A
   fbkg->SetLineColor(kGreen);
   fbkg->DrawCopy("SAME");
   leg->AddEntry(fbkg->Clone(), Form("RevArgus A=%.2f Sl=%.2f",
            fbkg->GetParameter(1),fbkg->GetParameter(2)),"L");

   leg->Draw();
}

//--------------------------------------------------------------------
void bkg_fit(int date, int type=0) {
//--------------------------------------------------------------------
   string fname, title, pdf;
   if ( type == 0 ) { // KK-eta PHSP
      fname = string( Form("mckketa_kkmc_%02i.root",date%100) );
      title = string( Form("MC non-#phi KK#eta %i",date) );
      pdf   = string( Form("mkk%02i_fitbg.pdf",date%100) );
   } else if ( type == 1 ) { // side-band events in data
      fname = string( Form("data_%02ipsip_all.root",date%100) );
      title = string( Form("Data %i (side-band)",date) );
      pdf   = string( Form("mkk_sb%02i_fitar.pdf",date%100) );
   } else {
      cerr << "ERROR: unknown type = " << type << endl;
      exit(EXIT_FAILURE);
   }

   // Get un-binned data
   TH1D* hist[1];
   vector<double> mkk = get_mkk_hist( fname, "mkk", hist, type );
   TH1D* hst = hist[0];
   double norm = hst->Integral();

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Bkg(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Bkg.Add(mkk[i]);
   }
   cout << " norm= " << norm << " n= " << n << endl;

   //-----------------------------------------------------------------
   // Fit KKeta background by reversed Argus function
   // function MUST be normalized to 1 on the fit range
   // double slope = -1.8; // prod-12
   double slope = -2.0; // v709
   auto Lan = [slope](const double* x,const double* p) -> double {
      return RevArgusN(x[0],p[0],slope);
   };
   vector<string> par_name { "a" };
   vector<double> par_ini  {  0. };

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* bkg = new TF1("bkg", Lan, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 fbkg( *bkg );
   fitter.SetFunction(fbkg,false); // false=no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(-0.0001, 10.);

   fitter.LikelihoodFit(Bkg,false); // true=extended likelihood fit
   fitter.CalculateHessErrors(); // without MINOS here
   // fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   size_t max_mkk = min( 4000, int(mkk.size()) );
   auto Fun = [Lan,Fpar](double x) -> double {
      return Lan(&x,Fpar.data());
   };
   double pvalueKS = FastKSTest(Fun,mkk,max_mkk);
   printf(" pvalueKS= %.3g (%zu)\n", pvalueKS,max_mkk);

   //-----------------------------------------------------------------
   // Functions to draw
   double Wn = bW*n;
   auto Lfit = [Wn,Lan,Fpar](double* x,double* p) -> double {
      return Wn * Lan(x,Fpar.data());
   };
   TF1* ffit = new TF1("ffit", Lfit, dL-1e-2, dU, 0);
   ffit->SetLineWidth(2);
   ffit->SetLineColor(kBlue);
   ffit->SetNpx(200);

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst);
   double hst_max = hst->GetMaximum();
   if ( hst_max < 10 ) {
      hst->SetMaximum( hst_max+2);
   } else {
      hst->SetMaximum( 1.3*hst_max );
   }
   hst->GetYaxis()->SetMaxDigits(3);
   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetTitleOffset(1.25);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
   ffit->Draw("SAME");

   TLegend* leg = new TLegend(0.11,0.77,0.50,0.89);
   leg->AddEntry(hst,title.c_str(),"LEP");
   leg->AddEntry(ffit,"Argus function","L");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.55,0.78,0.89,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   // pt->AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt->AddText( Form("#it{p-value(KS) = %.3f}",pvalueKS) );
   if ( type == 1 ) { // low statistic
      pt->AddText( Form("a = %s",PE.Eform(0,".1f")) );
   } else {
      pt->AddText( Form("a = %s",PE.Eform(0,".2f")) );
   }
   pt->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 data_fit() no interference
//--------------------------------------------------------------------
void data_fit(int date, int Mphix, string pdf="") {
//--------------------------------------------------------------------
   string fname( Form("data_%02ipsip_all.root",date%100) );

   // Get un-binned data
   TH1D* hist[1];
   vector<double> mkk = get_mkk_hist( fname, "mkk_data_cp", hist );
   TH1D* hst = hist[0];
   double norm = hst->Integral();

   int n = mkk.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Dat(dr, n);
   for ( int i = 0; i < n; ++i ) {
      Dat.Add(mkk[i]);
   }
   cout << " norm= " << norm << " n= " << n << endl;

   //-----------------------------------------------------------------
   // Fit data
   // double slope = -1.8; // prod-12
   double slope = -2.0; // v709
   // function represents the total number of events (extended LH)
   auto Lsum = [slope](const double* x,const double* p) -> double {
      return p[4] * BreitWignerGaussN(x[0],p[0],p[1],p[2],slope) +
             p[5] * RevArgusN(x[0],fabs(p[3]),slope);
   };
   vector<string> par_name {"Mphi","Gphi","sigma","a","Nphi","Nbkg"};
   vector<double> par_ini  { Mphi, Gphi, 1.16e-3,
                             0., 0.94*norm, 0.06*norm };

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* fsum = new TF1("fsum", Lsum, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 WFsum( *fsum );
   fitter.SetFunction(WFsum,false); // false=no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   if ( Mphix == 1 ) {
      fitter.Config().ParSettings(0).SetValue(1.01952); // like MC
      fitter.Config().ParSettings(0).Fix();
   } else {
      fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   }
   fitter.Config().ParSettings(1).Fix(); // Gphi
   fitter.Config().ParSettings(2).SetLimits(0.5e-3, 2.e-3);// sigma
   fitter.Config().ParSettings(3).SetLimits(-5.,5.);       // Argus
   fitter.Config().ParSettings(3).Fix();                   // Ar
   // fitter.Config().ParSettings(4).SetLimits(0.,1.5*norm);  // Nphi
   // fitter.Config().ParSettings(5).SetLimits(0.,norm);      // Nbkg

   fitter.LikelihoodFit(Dat,true); // true=extended likelihood fit
   // fitter.CalculateHessErrors();
   fitter.CalculateMinosErrors();

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
   printf(" pvalueKS= %.3g\n", pvalueKS);

   //-----------------------------------------------------------------
   // Functions to draw
   auto Lfit = [Lsum,Fpar](double* x,double* p) -> double {
      return bW * Lsum(x,Fpar.data());
   };
   TF1* ffit = new TF1("ffit", Lfit, dL, dU, 0);
   ffit->SetLineWidth(2);
   ffit->SetLineColor(kRed);
   ffit->SetNpx(500);

   auto Lbkg = [Fpar,slope](double* x,double* p) -> double {
      return bW * Fpar[5] * RevArgusN(x[0],fabs(Fpar[3]),slope);
   };
   TF1* fbkg = new TF1("fbkg", Lbkg, dL, dU, 0);
   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);
   fbkg->SetLineStyle(kDashed);

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst);
   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetTitleOffset(1.25);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
   ffit->Draw("SAME");
   fbkg->Draw("SAME");

   TLegend* leg = new TLegend(0.55,0.76,0.89,0.89);
   leg->AddEntry(hst,Form("Data %i",date),"LEP");
   leg->AddEntry(ffit,"Fit result","L");
   leg->AddEntry(fbkg,"non-#phi KK#eta","L");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.55,0.45,0.89,0.75,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   // pt->AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt->AddText( Form("#it{p-value(KS) = %.3f}",pvalueKS) );
   pt->AddText( Form("N_{#phi}= %s",PE.Eform(4,".1f")) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt->AddText( Form("#sigma = %s MeV", PE.Eform(2,".2f",1e3)) );
   pt->AddText( Form("N_{non-#phi}= %s",PE.Eform(5,".1f")) );
   pt->AddText( Form("a = %s",PE.Eform(3,".2f")) );
   pt->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// dataSB_fit() no interference + side-band fit
//--------------------------------------------------------------------
struct myFCN_nointer_sb {
   // efficiency parameters:
   // const double sl = -1.8; // prod-12
   const double sl = -2.0; // v709
   // Function for side-band: 0 - constant, 1 - Argus
   const int funSB = 1;

   vector<double> mkk;     // data central part
   vector<double> sb;      // data side band

   // minimization function
   // the signature of this operator() MUST be exactly this:
   // THE EXTENDED MAXIMUM LIKELIHOOD ESTIMATOR
   // Kai-Feng Chen "FittingMinuitRooFit" p.35
   //-----------------------------------------------------------------
   double operator() (const double* p) {
   //-----------------------------------------------------------------
      double ar    = p[3];
      double Nphi  = p[4];
      double Nnphi = p[5];
      double Nsb   = p[6];
      double asb   = p[7];

      double res = 2*(Nphi+Nnphi+Nsb);
      for ( const auto& m : mkk ) {
         double L = Nphi  * BreitWignerGaussN(m,p[0],p[1],p[2],sl) +
                    Nnphi * RevArgusN(m,ar,sl);
         if ( funSB == 0 ) {
            L += Nsb/(dU-dL);
         } else {
            L += Nsb * RevArgusN(m,asb,sl);
         }
         if (L > 0.) {
            res -= 2*log(L);
         } else {
            res += FLT_MAX; // ~3e38
         }
      }

      // fit SB by constant
      res += 2*Nsb;
      if ( funSB == 0 ) {
         // by constant (mkk >= dL=2Mk)
         res += -double(sb.size()) * 2*log(Nsb/(dU-dL));
      } else {
         // by Argus
         for ( const auto& m : sb ) {
            double L = Nsb * RevArgusN(m,asb,sl);
            if (L > 0.) {
               res -= 2*log(L);
            } else {
               res += FLT_MAX; // ~3e38
            }
         }
      }

      return res;
   }
};

//--------------------------------------------------------------------
void dataSB_fit( int date, int Mphix, string pdf="") {
//--------------------------------------------------------------------
   bool NoNonPhi = false; // if true Nnonphi fix to 0.
   string fname( Form("data_%02ipsip_all.root",date%100) );

   myFCN_nointer_sb my_fcn;        // class for 'FitFCN'
   const double sl = my_fcn.sl;    // efficiency parameters
   const int FunSB = my_fcn.funSB; // 0 - constant, 1 - Argus

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
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   vector<string> par_name { "Mphi", "Gphi", "sig", "a",
                             "Nphi", "Nnon", "Nbkg", "asb" };
   vector<double> par_ini;
   if ( date == 2009 ) {
      par_ini = { Mphi,Gphi,1.3e-3,0., 800,70,7,7. };
      // par_ini = { Mphi,Gphi,0.2e-3,0., 685,2,3,0. }; // I/O
   } else if ( date == 2012 ) {
      par_ini = { Mphi,Gphi,1.1e-3,0., 2550,200,33,2.5};
      // par_ini = { Mphi,Gphi,1.1e-3,0.,2700,10,10,0.}; // I/O
   } else if ( date == 2021 ) {
      par_ini = { Mphi,Gphi,1.1e-3,0., 16200,1300,190.,5.};
   }

   const unsigned int Npar = par_name.size(); // number of parameters

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   if ( Mphix == 1 ) {
      fitter.Config().ParSettings(0).SetValue(1.01952); // like MC
      fitter.Config().ParSettings(0).Fix();
   } else {
      fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   }
   fitter.Config().ParSettings(1).Fix(); // Gphi
   fitter.Config().ParSettings(2).SetLimits(0.2e-3, 2.e-3); // sig
   fitter.Config().ParSettings(3).SetLimits(-5.,5.);        // a
   fitter.Config().ParSettings(3).Fix();
   // fitter.Config().ParSettings(4).SetLimits(1.,1.5*norm);   // Nphi
   // fitter.Config().ParSettings(5).SetLimits(1.,norm);       // Nnon
   // fitter.Config().ParSettings(6).SetLimits(1.,2.*nsb);     // Nbkg
   // fitter.Config().ParSettings(7).SetLimits(0.,15.);        // asb
   if ( NoNonPhi ) {  // no non-phi KK eta
      fitter.Config().ParSettings(5).SetValue(0.); // no Nnonphi
      fitter.Config().ParSettings(5).Fix();
      // fitter.Config().ParSettings(6).Fix();        // fix Nbkg
      // fitter.Config().ParSettings(7).Fix();        // fix a(sb)
   }

   // == Fit
   int Ndat = n + nsb;
   fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // false=likelihood
   fitter.CalculateHessErrors();
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   vector<double> Fpar = res.Parameters();
   if ( res.IsValid() ) {
      printf("\n=> FINAL RESULT <=\n");
   } else {
      printf("\n=> INVALID RESULT <=\n");
   }
   res.Print(cout);
   // res.PrintCovMatrix(cout);

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res,0.05); // ignore 5% upper/lower errors

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   auto Fun = [Fpar,sl,FunSB](double x) -> double {
      double ret =
         Fpar[4] * BreitWignerGaussN(x,Fpar[0],Fpar[1],Fpar[2],sl) +
         Fpar[5] * RevArgusN(x,Fpar[3],sl);
      if ( FunSB == 0 ) {
         ret += Fpar[6]/(dU-dL);
      } else {
         ret += Fpar[6] * RevArgusN(x,Fpar[7],sl);
      }
      // normalization on 1
      double Ntot = Fpar[4] + Fpar[5] + Fpar[6];
      return ret/Ntot;
   };
   double pvalueKS = FastKSTest(Fun,mkk);
   printf(" pvalueKS(cr)= %.3g\n", pvalueKS);

   // side-band
   auto Fsb = [Fpar,sl,FunSB](double x) -> double {
      double ret = 0;
      if ( FunSB == 0 ) {
         ret = 1./(dU-dL);
      } else {
         ret = RevArgusN(x,Fpar[7],sl);
      }
      return ret;
   };
   double pvalueKSsb = FastKSTest(Fsb,sb);
   printf(" pvalueKS(sb)= %.3g\n", pvalueKSsb);

   //-----------------------------------------------------------------
   // Functions to draw
   auto Ldr = [Fpar,sl,FunSB](double* x,double* p) -> double {
      int isw = int(p[0]+0.5);
      double m = x[0];
      double BW =
         Fpar[4] * BreitWignerGaussN(m,Fpar[0],Fpar[1],Fpar[2],sl);
      if ( isw == 1 ) { return bW * BW; }

      double Ar = Fpar[5] * RevArgusN(m,Fpar[3],sl);
      if ( isw == 2 ) { return bW * Ar; }

      double Bg = 0;
      if ( FunSB == 0 ) {
         Bg = Fpar[6]/(dU-dL);
      } else {
         Bg = Fpar[6] * RevArgusN(m,Fpar[7],sl);
      }
      if ( isw == 3 ) { return bW * Bg; }

      return bW * (BW + Ar + Bg);
   };
   TF1* fdr = new TF1("fdr", Ldr, bL, dU, 1);
   fdr->SetNpx(500);

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst);
   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetTitleOffset(1.25);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");

   TLegend* leg = new TLegend(0.55,0.745,0.89,0.89);
   leg->AddEntry(hst,Form("Data %i",date),"LEP");

   fdr->SetParameter(0, 0); // Sum
   fdr->SetLineWidth(2);
   fdr->SetLineColor(kRed+1);
   fdr->DrawCopy("SAME");
   leg->AddEntry(fdr->Clone(),"Result of fit","L");

   fdr->SetParameter(0, 1); // BW
   fdr->SetLineWidth(2);
   fdr->SetLineStyle(kDashed);
   fdr->SetLineColor(kGreen+2);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Breit-Wigner #phi#eta", "L");

   fdr->SetParameter(0, 2); // Ar
   fdr->SetLineWidth(2);
   fdr->SetLineStyle(kDashed);
   fdr->SetLineColor(kBlue);
   fdr->DrawCopy("SAME");
   leg->AddEntry(fdr->Clone(),"non-#phi KK#eta","L");

   // fdr->SetParameter(0, 3); // SB
   // fdr->SetLineWidth(1);
   // fdr->SetLineStyle(kDashed);
   // fdr->SetLineColor(kCyan+3);
   // fdr->DrawCopy("SAME");
   // leg->AddEntry(fdr->Clone(),"background","L");

   leg->Draw();

   TPaveText* pt = new TPaveText(0.55,0.46,0.89,0.74,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#it{p-value(KS) = %.3f}",pvalueKS) );
   pt->AddText( Form("N_{#phi}= %s",PE.Eform(4,".1f")) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt->AddText( Form("#sigma = %s MeV", PE.Eform(2,".2f",1e3)) );
   pt->AddText( Form("N_{non-#phi}= %s",PE.Eform(5,".1f")) );
   pt->AddText( Form("a = %s",PE.Eform(3,".2f")) );
   pt->Draw();

   TPaveText* ptsb = new TPaveText(0.55,0.34,0.89,0.46,"NDC");
   ptsb->SetTextAlign(12);
   ptsb->SetTextFont(42);
   ptsb->AddText( Form("#it{p-val(side-band) = %.3f}", pvalueKSsb) );
   ptsb->AddText( Form("#lower[-0.1]{Nbg = %s}", PE.Eform(6,".1f")) );
   ptsb->AddText( Form("#lower[-0.1]{a(sb) = %s}",PE.Eform(7,".1f")));
   ptsb->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      if ( NoNonPhi ) {
         auto pos = pdf.find(".pdf");
         if ( pos != string::npos) {
            pdf = pdf.substr(0,pos) + "_noNonphi.pdf";
         }
      }
      c1->Print(pdf.c_str());
   }
}

// {{{1 Interference BW with Argus bkg.
// {{{2 Common
//--------------------------------------------------------------------
vector<double> IntfrBWAR( double m,
      double mphi, double gphi,      // B-W
      double A, double F, double Ang // Argus & interference
      ) {
//--------------------------------------------------------------------
// see memo and BreitWigner() function
// The return vector contains:
//                   0    1     2        3
//        ret     { Sum, B-W, Argus, Interference }
//--------------------------------------------------------------------
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
      ) {
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
   ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);
   double result = gsl_int.Integral(Fint,-5.,+5.);

   return Ngauss * result;
}

//--------------------------------------------------------------------
double IntfrBWARGN( double m,
      const double p[], // {mphi,gphi,sigma,A,F,Ang,sl}
      int idx = 0       // what to return
      ) {
//--------------------------------------------------------------------
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
         return IntfrBWARG(x,p,0);  // MUST be idx=0
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

// {{{2 test_Intfr()
//--------------------------------------------------------------------
void test_Intfr() {
//--------------------------------------------------------------------
   double mphi = Mphi, gphi = Gphi;
   auto Lintfr = [mphi,gphi](const double* x, const double* p) {
      double pp[] { mphi, gphi, p[0],p[1],p[2],p[3],p[4] };
      return IntfrBWARGN( x[0], pp, int(p[5]) );
   };

   int npar = 6;
   TF1* fun = new TF1("fun", Lintfr, dL, dU, npar);
   fun->SetParNames("Sigma", "A", "F", "Ang",  "Sl", "T");
   fun->SetParameters(1.2e-3, 0., 0.8,   0.7,  -1.99,  0.);
   fun->SetLineWidth(2);
   fun->SetNpx(500);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   TH1D* tmp = new TH1D("tmp","",500,bL,dU);
   tmp->SetAxisRange(-20.,130.,"Y");
   tmp->Draw();

   fun->SetParameter(npar-1, 0); // draw all
   fun->SetLineColor(kRed);
   fun->DrawCopy("SAME");
   printf(" Norm= %.6f\n", fun->Integral(dL,dU,1e-4) );

   fun->SetParameter(npar-1, 2); // draw Argus
   fun->SetLineColor(kBlue);
   fun->SetLineStyle(kDashed);
   fun->DrawCopy("SAME");

   fun->SetParameter(npar-1, 1); // draw BW
   fun->SetLineColor(kGreen+2);
   fun->DrawCopy("SAME");

   fun->SetParameter(npar-1, 3); // draw interference
   fun->SetLineColor(kMagenta);
   fun->SetLineStyle(kSolid);
   fun->DrawCopy("SAME");

   c1->Update();
}

// {{{2 data_Intfr_fit()
//--------------------------------------------------------------------
ValErr CalcIntEr( const ROOT::Fit::FitResult& res,
      double sl, unsigned int idx ) {
//--------------------------------------------------------------------
// Numerical calculation of integrals and their errors
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

//--------------------------------------------------------------------
void data_Intfr_fit(int date, string pdf="") {
//--------------------------------------------------------------------
   string fname( Form("data_%02ipsip_all.root",date%100) );

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

   //-----------------------------------------------------------------
   // Fit data
   // efficiency parameters
   // double sl = -1.8; // prod-12
   double sl = -2.0; // v709
   // function MUST be normalized to 1 on the fit range
   auto Lfit = [sl](const double* x,const double* p) -> double {
      const double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl };
      return IntfrBWARGN( x[0], pp, 0 );
   };
   vector<string> par_name { "Mphi","Gphi","sig","a","F","angle" };

   vector<double> par_ini;
   if ( date == 2009 ) {
      par_ini = {Mphi,Gphi,1.4e-3,0.,1.1,1.1}; //n
      // par_ini = {Mphi,Gphi,1.4e-3,0.,1.2,-1.1}; //p
   } else if ( date == 2012 ) {
      par_ini = {Mphi,Gphi,1.1e-3,0.,1.0,0.7}; //n
      // par_ini = {Mphi,Gphi,1.1e-3,0.,1.0,-0.8}; //p
   } else if ( date == 2021 ) {
      par_ini = {Mphi,Gphi,1.1e-3,0.,1.0,0.8}; //n
      // par_ini = {Mphi,Gphi,1.1e-3,0.,1.1,-0.8}; //p
   }
   bool neg_pos = par_ini.back() > 0; //true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 WFit( *Ffit );
   fitter.SetFunction(WFit,false); // false=no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   fitter.Config().ParSettings(0).SetValue(1.01952);        // like MC
   fitter.Config().ParSettings(0).Fix();                    // Mphi
   fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3, Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();                    // Gphi
   fitter.Config().ParSettings(2).SetLimits(0.5e-3, 2.e-3); // sig
   fitter.Config().ParSettings(3).Fix();                    // a
   fitter.Config().ParSettings(4).SetLimits(0.01, 10.);     // F
   fitter.Config().ParSettings(5).SetLimits(-M_PI, M_PI);   // angle

   fitter.LikelihoodFit(Dat,false); // true=extended likelihood fit
   fitter.CalculateMinosErrors();
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
      // printf("Int[%i] = %s\n",idx,Integ.prt(".4f")); // DEBUG
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
   double fmax = ffit->GetMaximum( 1.01, 1.03, 1e-6);
   hmax = floor(1.1 * max( hmax, fmax ));
   if ( hmax > 0 ) {
      hst->SetMaximum(hmax);
   }

   ffit->SetParameter(0, 3); // interference
   double hmin = ffit->GetMinimum( 1.01, 1.03, 1e-6);
   hmin = floor(1.15 * min( 0., hmin ));
   if ( hmin < 0 ) {
      hst->SetMinimum(hmin);
   }

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst);

   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetTitleOffset(1.3);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");

   TLegend* leg = new TLegend(0.58,0.725,0.89,0.89);
   leg->AddEntry(hst,Form("Data %i",date),"LEP");

   ffit->SetParameter(0, 0); // SUM
   ffit->SetLineWidth(2);
   ffit->SetLineColor(kRed);
   ffit->DrawCopy("SAME");
   leg->AddEntry( ffit->Clone(), "Result of fit", "L");

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

   TPaveText* pt = new TPaveText(0.58,0.48,0.89,0.72,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   // pt->AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt->AddText( Form("#it{p-value(KS) = %.3f}",pvalueKS) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt->AddText( Form("#sigma = %s MeV", PE.Eform(2,".2f",1e3)) );
   pt->AddText( Form("a = %s",PE.Eform(3,".1f")) );
   pt->AddText( Form("F = %s",PE.Eform(4,".2f")) );
   pt->AddText( Form("#vartheta = %s",PE.Eform(5,".2f")) );
   pt->Draw();
   TPaveText* ptr = new TPaveText(0.58,0.45,0.89,0.48,"NDC");
   ptr->SetTextAlign(12);
   ptr->SetTextFont(42);
   ptr->AddText( Form("#lower[0.1]{N_{#phi} = %.1f #pm %.1f}",
            Nphi,err_Nphi) );
   ptr->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      pdf += (neg_pos ? "_n.pdf" : "_p.pdf");
      c1->Print(pdf.c_str());
   }
}

// {{{2 dataSB_Intfr(): combined with side-band
//--------------------------------------------------------------------
struct myFCN_inter {
   // efficiency parameter:
   // const double sl = -1.8; // prod-12
   const double sl = -2.0; // v709
   // Function for side-band:
   const int funSB = 1; // 0 - constant, 1 - Argus

   vector<double> mkk;  // data central part
   vector<double> sb;   // data side band

   // minimization function
   // the signature of this operator() MUST be exactly this:
   // THE EXTENDED MAXIMUM LIKELIHOOD ESTIMATOR
   // Kai-Feng Chen "FittingMinuitRooFit" p.35
   //-----------------------------------------------------------------
   double operator() (const double* p) {
   //-----------------------------------------------------------------
      double Nkk = p[6];
      double Nsb = p[7];
      double asb = p[8];

      double res = 2*(Nkk+Nsb);

      const double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl };
      for ( const auto& m : mkk ) {
         double L = Nkk * IntfrBWARGN( m,pp,0 );
         if ( funSB == 0 ) {
            L += Nsb/(dU-dL);
         } else {
            L += Nsb * RevArgusN(m,asb,sl);
         }

         if (L > 0.) {
            res -= 2*log(L);
         } else {
            res += FLT_MAX; // ~3e38
         }
      }

      // fit SB
      res += 2*Nsb;
      if ( funSB == 0 ) {
         // by constant (mkk >= dL=2Mk)
         res += -double(sb.size()) * 2*log(Nsb/(dU-dL));
      } else {
         // by Argus
         for ( const auto& m : sb ) {
            double L = Nsb * RevArgusN(m,asb,sl);
            if (L > 0.) {
               res -= 2*log(L);
            } else {
               res += FLT_MAX; // ~3e38
            }
         }
      }

      return res;
   }
};

//--------------------------------------------------------------------
ValErr CalcIntErSb( const ROOT::Fit::FitResult& res,
      double sl, unsigned int idx, bool noerrors=false) {
//--------------------------------------------------------------------
// Numerical calculation of integrals and their errors
   // desired errors:
   constexpr double eps_abs = 1.e-6;
   constexpr double eps_rel = 1.e-6;
   constexpr double epsilon = 1.e-3; // max numerical error

   // get parameters and covariation matrix of error
   int Npar = res.NPar();
   vector<double> Fpar = res.Parameters();

   if ( Npar != 9 ) {
      cerr << " FATAL ERROR in " << __func__
         << " Npar= " << Npar << endl;
      exit(EXIT_FAILURE);
   }
   Npar = 6; // remove NKK, Nbg and asb parameters!

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
      double q[] { p[0],p[1],p[2],p[3],p[4],p[5],sl }; // real slope
      double normI = IntfrBWARGN( m,q,100 );
      q[6] = 0.; // set slope to zero
      double intfr = IntfrBWARG( m,q,idx ) / normI;
      return intfr;
   };
   TF1 FC("FC", Lfc, dL, dU, Npar); // dL !
   FC.SetParameters( Fpar.data() );

   double err_num = 0;
   double Integ = FC.IntegralOneDim(dL,dU,eps_rel,eps_abs,err_num);
   // printf(" err_FC_Integ[%d] = %.2g\n",idx,err_num);
   if ( err_num > epsilon ) {
      printf(" WARNING: %s numerical error of integration of "
            "F_C[%d] is too big: %.2g\n", __func__,idx,err_num);
   }

   if ( noerrors ) { // skip error calculation
      return (ValErr {Integ,0});
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
      IntGrad[i] = dF.IntegralOneDim(dL,dU,eps_rel,eps_abs,err_num);
      err_num = covMatrix(i,i) * IntGrad[i] * err_num;
      err_num2 += SQ(err_num);
      // printf(" IntGrad(dF[%i]/d_%d)= %.4g -> err_num= %.4g\n",
            // idx,i, IntGrad[i], err_num);
   }

   double err_Int = sqrt( covMatrix.Similarity(IntGrad) );
   // printf(" Integral[%d] = %.6g +/- %.2g\n",idx,Integ,err_Int);

   err_num2 = sqrt( err_num2 / err_Int ); // abs numerical error
   if ( err_num2 > epsilon ) {
      printf(" WARNING: %s numerical error of integration of "
            "dF[%d] is too big: %.2g\n", __func__,idx,err_num2);
   }

   return (ValErr {Integ,err_Int});
}

//--------------------------------------------------------------------
void dataSB_Intfr(int date, string pdf, bool plotsb=true) {
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

   my_fcn.sb = get_mkk_hist( fname, "mkk_sb", &hist[1], 1 );
   TH1D* hsb = hist[1];
   const vector<double>& sb = my_fcn.sb;

   //-----------------------------------------------------------------
   // Fit data
   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   vector<string> par_name
      { "Mphi","Gphi","sig","a","F","angle","NKK","Nbg","asb" };

   vector<double> par_ini;
   if ( date == 2009 ) {
      par_ini = {Mphi,Gphi,1.4e-3,0., 1.1, 1.0, 890., 10., 7. };//n
      // par_ini = {Mphi,Gphi,1.4e-3,0., 1.1,-1.0, 890., 10., 7. };//p
   } else if ( date == 2012 ) {
      par_ini = {Mphi,Gphi,1.1e-3,0., 0.8, 0.3, 2750.,35., 4. };//n
      // par_ini = {Mphi,Gphi,1.1e-3,0., 0.8,-0.3, 2750.,35., 4. };//p
   } else if ( date == 2021 ) {
      par_ini = {Mphi,Gphi,1.1e-3,0., 0.9, 0.5, 17500.,195.,6.};//n
      // par_ini = {Mphi,Gphi,1.1e-3,0., 0.9,-0.5, 17500.,195.,6.};//p
   }
   bool neg_pos = par_ini[5] > 0; // true for negative interference

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
   fitter.Config().ParSettings(5).SetLimits(-M_PI, M_PI);   // angle
   // fitter.Config().ParSettings(6).SetLimits(0.,1.5*mkk.size());// NKK
   // fitter.Config().ParSettings(7).SetLimits(0., 2.*sb.size()); // Nbg
   if ( FunSB == 0 ) {                                      // asb
      fitter.Config().ParSettings(8).SetValue(0.);
      fitter.Config().ParSettings(8).Fix();
   } else {
      fitter.Config().ParSettings(8).SetLimits(0.,15.);
   }

   // == Fit
   int Ndat = mkk.size() + sb.size();
   fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // false=likelihood
   fitter.CalculateMinosErrors();
   fitter.CalculateHessErrors(); // for correct integral errors

   const ROOT::Fit::FitResult& res = fitter.Result();
   vector<double> Fpar = res.Parameters();
   if ( res.IsValid() ) {
      printf("\n=> FINAL RESULT <=\n");
   } else {
      printf("\n=> INVALID RESULT <=\n");
   }
   res.Print(cout);
   // res.PrintCovMatrix(cout);

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res,0.05); // ignore 5% upper/lower errors

   double Nfit = Fpar[6];
   double Nfit_err = max({PE.Perr[6],PE.Uerr[6],PE.Lerr[6]});
   cout << "NKK_fit(no slope corrections)= " << Nfit
        << " +/- " << Nfit_err << endl;
   vector<string> names { "NKK", "Nphi", "Nnonphi", "Nifr" };
   vector<ValErr> Nos( names.size() );
   for ( int idx = 0; idx <= 3; ++idx ) {
      ValErr Integ = CalcIntErSb( res, sl, idx );
      // printf("Int[%i] = %s\n",idx,Integ.prt(".6f")); // DEBUG
      double Num = Nfit * Integ.val;
      double Err = fabs(Num) *
         sqrt( SQ(Integ.err/Integ.val) + SQ(Nfit_err/Nfit) );
      Nos[idx] = {Num,Err};
      printf("%7s = %s\n",names[idx].c_str(),Nos[idx].prt("7.1f"));
   }
   double Nkk = Nos[0].val, err_Nkk = Nos[0].err;
   double Nphi = Nos[1].val, err_Nphi = Nos[1].err;

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   // central part
   auto Fcr = [Fpar,sl,FunSB](double x) -> double {
      double pp[]
      {Fpar[0],Fpar[1],Fpar[2],Fpar[3],Fpar[4],Fpar[5],sl};
      double ret = Fpar[6] * IntfrBWARGN( x,pp,0 );
      if ( FunSB == 0 ) {
         ret += Fpar[7]/(dU-dL);
      } else {
         ret += Fpar[7] * RevArgusN( x,Fpar[8],sl );
      }
      // normalization on 1
      double Ntot = Fpar[6] + Fpar[7];
      return ret/Ntot;
   };
   double pvalueKS = FastKSTest(Fcr,mkk);
   printf("pvalueKS(cr)= %.3g\n", pvalueKS);

   // side-band
   auto Fsb = [Fpar,sl,FunSB](double x) -> double {
      double ret = 0;
      if ( FunSB == 0 ) {
         ret = 1./(dU-dL);
      } else {
         ret = RevArgusN( x,Fpar[8],sl );
      }
      return ret;
   };
   double pvalueKSsb = FastKSTest(Fsb,sb);
   printf(" pvalueKS(sb)= %.3g\n", pvalueKSsb);

   //-----------------------------------------------------------------
   // Functions to draw
   auto Ldr = [Fpar,sl,FunSB](double* x,double* p) -> double {
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

   // find max/min to draw
   fdr->SetParameter(0, 1); // BW
   int maxbin = hst->GetMaximumBin();
   double hmax = hst->GetBinContent(maxbin)+hst->GetBinError(maxbin);
   double fmax = fdr->GetMaximum( 1.01, 1.03, 1e-6);
   hmax = floor(1.1 * max( hmax, fmax ));
   if ( hmax > 0 ) {
      hst->SetMaximum(hmax);
   }

   fdr->SetParameter(0, 3); // interference
   double hmin = fdr->GetMinimum( 1.01, 1.03, 1e-6);
   hmin = floor(1.15 * min( 0., hmin ));
   if ( hmin < 0 ) {
      hst->SetMinimum(hmin);
   }

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   pdf += (neg_pos ? "_n.pdf" : "_p.pdf");
   c1->Print((pdf+"[").c_str()); // open pdf-file

   SetHstFace(hst);
   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetTitleOffset(1.3);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");

   TLegend* leg = new TLegend(0.58,0.725,0.89,0.89);
   leg->AddEntry(hst,Form("Data %i",date),"LEP");

   fdr->SetParameter(0, 0); // SUM
   fdr->SetLineWidth(2);
   fdr->SetLineColor(kRed+1);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Result of fit", "L");

   fdr->SetParameter(0, 1); // BW
   fdr->SetLineWidth(2);
   fdr->SetLineStyle(kDashed);
   fdr->SetLineColor(kGreen+2);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Breit-Wigner #phi#eta", "L");

   fdr->SetParameter(0, 2); // Argus
   fdr->SetLineWidth(2);
   fdr->SetLineStyle(kDashed);
   fdr->SetLineColor(kBlue);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Non-#phi KK#eta", "L");

   fdr->SetParameter(0, 3); // interference
   fdr->SetLineWidth(2);
   fdr->SetLineStyle(kDashed);
   fdr->SetLineColor(kMagenta+1);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Interference", "L");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.58,0.48,0.89,0.72,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#it{p-value = %.3f}",pvalueKS) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt->AddText( Form("#sigma = %s MeV", PE.Eform(2,".2f",1e3)) );
   pt->AddText( Form("a = %s",PE.Eform(3,".1f")) );
   pt->AddText( Form("F = %s",PE.Eform(4,".2f")) );
   pt->AddText( Form("#vartheta = %s",PE.Eform(5,".2f")) );
   pt->AddText( Form("NKK(fit) = %s", PE.Eform(6,".1f")) );
   pt->Draw();

   TPaveText* ptsb = new TPaveText(0.58,0.39,0.89,0.48,"NDC");
   ptsb->SetTextAlign(12);
   ptsb->SetTextFont(42);
   ptsb->AddText( Form("#it{p-val(side-band) = %.3f}", pvalueKSsb) );
   ptsb->AddText( Form("#lower[-0.1]{Nbg = %s}",PE.Eform(7,".1f")) );
   if ( FunSB != 0 ) {
      ptsb->AddText( Form("#lower[-0.1]{a(sb) = %s}",
               PE.Eform(8,".1f")) );
   }
   ptsb->Draw();

   TPaveText* ptres = new TPaveText(0.58,0.33,0.89,0.39,"NDC");
   ptres->SetTextAlign(12);
   ptres->SetTextFont(42);
   ptres->AddText(Form("N_{#etaKK} = %.1f #pm %.1f", Nkk,err_Nkk));
   ptres->AddText(Form("N_{#phi } = %.1f #pm %.1f",Nphi,err_Nphi));
   ptres->Draw();

   gPad->RedrawAxis();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   if ( !plotsb ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
      return;
   }

   //-----------------------------------------------------------------
   // Side-band

   SetHstFace(hsb);
   maxbin = hsb->GetMaximumBin();
   hmax = hsb->GetBinContent(maxbin)+hsb->GetBinError(maxbin);
   if ( hmax < 10 ) {
      hsb->SetMaximum( hmax+2);
   } else {
      hsb->SetMaximum( 1.2*hmax );
   }
   hsb->GetXaxis()->SetTitleOffset(1.1);
   hsb->GetYaxis()->SetTitleOffset(1.3);
   hsb->SetLineWidth(2);
   hsb->SetLineColor(kBlack);
   hsb->SetMarkerStyle(20);

   hsb->Draw("EP");

   TLegend* leg2 = new TLegend(0.55,0.81,0.89,0.89);
   leg2->AddEntry(hsb, Form("Side-band: data %i",date),"LEP");

   fdr->SetParameter(0, 4); // SB
   fdr->SetLineWidth(2);
   fdr->SetLineStyle(kSolid);
   fdr->SetLineColor(kRed+1);
   fdr->DrawCopy("SAME");
   leg2->AddEntry( fdr->Clone(), "Fit by Argus function", "L");
   leg2->Draw();

   ptsb->SetX1NDC(0.55);
   ptsb->SetX2NDC(0.89);
   ptsb->SetY1NDC(0.70);
   ptsb->SetY2NDC(0.80);
   ptsb->Draw();

   gPad->RedrawAxis();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
}

// {{{2  dataSB_Intfr_scan()
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

// {{{2 dataSBBR_Intfr(): combined with side-band + fitBR
//--------------------------------------------------------------------
struct myFCN_sbbr {
   // efficiency parameter:
   // const double sl = -1.8; // prod-12
   const double sl = -2.0; // v709
   // Function for side-band:
   const int funSB = 1; // 0 - constant, 1 - Argus

   vector<double> mkk;  // data central part
   vector<double> sb;   // data side band

   const double mphi = 1.01952; // Mphi like in MC
   const double gphi = Gphi;
   const double ar = 0.; // Argus parameter
   double Br2Nkk = 0.;   // convertion Br(J/Psi -> KK eta) -> Nkk
   double Br2Nphi = 0.;  // convertion Br(J/Psi -> phi eta) -> Nphi

   // integrals
   double IeBW = 0, IeAr = 0, IeIc = 0, IeIs = 0;
   double IBW = 0, IAr = 0, IIc = 0, IIs = 0;

   // normalizations for Argus in side-band
   double normArsb = 1;

   //-----------------------------------------------------------------
   myFCN_sbbr(int date) { // set parameters for 2009/2012
   //-----------------------------------------------------------------
      // Br(eta->2gamma) = 39.36 +/- 0.18% PDG 2023
      const double breta = 0.3936;
      // Br(phi->K+K-) = 49.1 +/- 0.5% PDG 2023
      const double brphi = 0.491;

      // prod-12
      // const double NJpsi09 = 15740130; // 2009
      // const double NJpsi12 = 48674393; // 2012
      // const double e_phi09 = 0.3280; // eff chi2<60
      // const double e_phi12 = 0.3229; // eff chi2<60

      // prod v709 first
      const double NJpsi09 =  15990331; // 2009
      const double NJpsi12 =  49882436; // 2012
      const double NJpsi21 = 324842271; // 2021
      const double e_phi09 = 0.2841; // eff chi2<40
      const double e_phi12 = 0.3112; // eff chi2<40
      const double e_phi21 = 0.3163; // eff chi2<40

      if ( date == 2009 ) {
         Br2Nkk = NJpsi09 * e_phi09 * breta;
      } else if ( date == 2012 ) {
         Br2Nkk = NJpsi12 * e_phi12 * breta;
      } else if ( date == 2021 ) {
         Br2Nkk = NJpsi21 * e_phi21 * breta;
      }
      Br2Nphi = Br2Nkk * brphi;
   }

   // normalization for Argus SB
   //-----------------------------------------------------------------
   void calcNormAr(const double asb) {
   //-----------------------------------------------------------------
      // cache for calculated values
      static double asb_save = -1.;

      if ( funSB == 0 ) {
         return;
      }
      if ( asb_save == asb ) {
         return;
      }
      asb_save = asb;

      double par[] = {asb,sl};
      auto Lintar = [](double x, void* pp) -> double{
         double* p = static_cast<double*>(pp);
         return RevArgus(x,p[0]) * (1+p[1]*(x-1.02));
      };
      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);
      normArsb = gsl_int.Integral(Lintar,(void *)par,dL,dU);

      // debug print
      // printf(" asb= %.3f, normArsb= %.3g\n", asb, normArsb);
   }

   //-----------------------------------------------------------------
   void calcIntegrals(double sig) {
   //-----------------------------------------------------------------
      // cach old values
      static double sig_save = 0.;
      if ( sig_save == sig ) {
         return;
      }
      sig_save = sig;

      // integrand lambda function: mphi,gphi,sig,ar,F,ang,sl,idx
      const double F = 1.;
      const double ang = 0.; //0(PI/2) for cos(sin) members of Infr.
      double pp[] { mphi,gphi,sig,ar,F,ang,sl, 1 };

      auto Lint = [](double x, void* pp) -> double{
         const double* p = static_cast<const double*>(pp);
         int idx = int(p[7]);
         return IntfrBWARG(x,p,idx);
      };
      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);

      // function for check results
      auto checkInt = [](string t, double Int, double pp[]) -> void {
         if ( !isfinite(Int) ) {
            printf("%s: pp= %g,%g,%g,%g,%g,%g,%g,%g\n", t.data(),
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
         Ff = max(0.,-B/(2*A));
         penalty += 1e5*fabs(Dis);
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
   // THE EXTENDED MAXIMUM LIKELIHOOD ESTIMATOR
   // Kai-Feng Chen "FittingMinuitRooFit" p.35
   //-----------------------------------------------------------------
   double operator() (const double* p) {
   //-----------------------------------------------------------------
      if ( !isfinite(p[0]) || !isfinite(p[1]) || !isfinite(p[2]) ||
           !isfinite(p[3]) || !isfinite(p[4]) || !isfinite(p[5]) ) {
         printf("NAN: p[]= %g,%g,%g,%g,%g,%g\n",
               p[0],p[1],p[2],p[3],p[4],p[5]);
         return DBL_MAX;
      }

      double Brkk = p[0];
      double Brphi = p[1];
      double Nkk = Brkk * Br2Nkk;
      double Nphi = Brphi * Br2Nphi;

      double ang = p[2];
      double sigma = p[3];
      calcIntegrals(sigma);

      double Ff = 0., penalty = 0.;
      tie(Ff,penalty) = calcF(p);

      double Nsb = p[4];
      double asb = p[5];
      calcNormAr(asb);
      const double NsbNorm = Nsb / normArsb;

//+ this long calculation depends only on Br(eta->2gamma)...
//+       double normI = IBW+Ff*((IIc*cos(ang)+IIs*sin(ang))+Ff*IAr);
//+       double NFit = Nkk/normI;
//+ short computation is mathematically completely equivalent to long
      double NFit = Nphi/IBW; // === Nkk/normI;

      double normE = IeBW+Ff*((IeIc*cos(ang)+IeIs*sin(ang))+Ff*IeAr);
      double NkkFit = NFit*normE;

      double res = 2*(NkkFit+Nsb) + penalty;

      const double pp[] { mphi,gphi,sigma,ar,Ff,ang,sl };
      for ( const auto& m : mkk ) {
         double L = NFit * IntfrBWARG( m,pp,0 );
         if ( funSB == 0 ) {
            L += Nsb/(dU-dL);
         } else {
            L += NsbNorm * RevArgus(m,asb)*(1+sl*(m-1.02));
         }
         if (L > 0.) {
            res -= 2*log(L);
         } else {
            res += FLT_MAX; // ~3e38
         }
      }

      // fit SB
      res += 2*Nsb;
      if ( funSB == 0 ) {
         // by constant (mkk >= dL=2Mk)
         res += -double(sb.size()) * 2*log(Nsb/(dU-dL));
      } else {
         // by Argus
         for ( const auto& m : sb ) {
            double L = NsbNorm * RevArgus(m,asb)*(1+sl*(m-1.02));
            if (L > 0.) {
               res -= 2*log(L);
            } else {
               res += FLT_MAX; // ~3e38
            }
         }
      }

      // printf("RES: p[]= %g,%g,%g,%g,%g,%g -> %.3f\n",
              // p[0],p[1],p[2],p[3],p[4],p[5],res);
      return res;
   }
};

//--------------------------------------------------------------------
void dataSBBR_Intfr(int date, string pdf, bool plotsb=true) {
//--------------------------------------------------------------------
   string fname( Form("data_%02ipsip_all.root",date%100) );

   myFCN_sbbr my_fcn(date);  // class for 'FitFCN'
   const double sl = my_fcn.sl; // efficiency parameters
   const int FunSB = my_fcn.funSB;

   // Get un-binned data
   TH1D* hist[2];
   my_fcn.mkk = get_mkk_hist( fname, "mkk_cp", &hist[0] );
   TH1D* hst = hist[0];
   const vector<double>& mkk = my_fcn.mkk;

   my_fcn.sb = get_mkk_hist( fname, "mkk_sb", &hist[1], 1 );
   TH1D* hsb = hist[1];
   const vector<double>& sb = my_fcn.sb;

   //-----------------------------------------------------------------
   // Fit data
   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   vector<string> par_name {"Brkk","Brphi","angle","sig","Nbg","asb"};

   vector<double> par_ini;
   if ( date == 2009 ) {
      par_ini = {5.0e-4, 10.3e-4, 1.0, 1.4e-3, 7., 7.}; //n
      // par_ini = {5.0e-4, 8.5e-4,-1.0, 1.4e-3, 7., 7.}; //p
   } else if ( date == 2012 ) {
      par_ini = {4.5e-4, 8.8e-4, 0.3, 1.1e-3, 34., 4.}; //n
      // par_ini = {4.5e-4, 8.4e-4,-0.3, 1.1e-3, 34., 4.}; //p
   } else if ( date == 2021 ) {
      // par_ini = {4.4e-4, 8.5e-4, 0.5, 1.1e-3, 195., 6.}; //n
      par_ini = {4.4e-4, 7.9e-4,-0.5, 1.1e-3, 195., 6.}; //p
   }
   bool neg_pos = par_ini[2] > 0; // true for negative interference

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
   // fitter.Config().ParSettings(3).Fix();
   // fitter.Config().ParSettings(4).SetLimits(0.,2*sb.size()); // Nbg
   fitter.Config().ParSettings(5).SetLimits(0., 15.);     // asb

   // == Fit
   int Ndat = mkk.size() + sb.size();
   fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // false=likelihood
   // fitter.CalculateHessErrors();
   fitter.CalculateMinosErrors();

   const ROOT::Fit::FitResult& res = fitter.Result();
   vector<double> Fpar = res.Parameters();
   // variables for KS-test and draw-function
   double Ff,penalty;
   tie(Ff,penalty) = my_fcn.calcF(Fpar.data());
   my_fcn.calcNormAr(Fpar[5]); // asb

   // print result
   if ( !isfinite(penalty) || penalty > 0. ) {
      printf("\n=> INVALID RESULT: penalty= %f <=\n",penalty);
   } else {
      printf("\n=> FINAL RESULT <=\n");
   }
   res.Print(cout);
   // res.PrintCovMatrix(cout);

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res,0.05); // ignore 5% upper/lower errors

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   // central part
   double Nkk = Fpar[0] * my_fcn.Br2Nkk;
   double Nphi = Fpar[1] * my_fcn.Br2Nphi;
   double NFit = Nphi/my_fcn.IBW;

   double ang = Fpar[2];
   double sig = Fpar[3];

   // {mphi,gphi,sigma,A,F,Ang,sl}
   const double pp[] {
      my_fcn.mphi, my_fcn.gphi, sig, my_fcn.ar, Ff, ang, my_fcn.sl
   };

   // normalization on 1
   double Nsb = Fpar[4];
   double asb = Fpar[5];
   double Ntot = Nkk + Nsb;
   NFit /= Ntot;

   double NsbNorm = Nsb/Ntot;
   if ( FunSB != 0 ) {
      NsbNorm /= my_fcn.normArsb;
   }

   auto Fcr = [NFit,pp,FunSB,NsbNorm,asb,sl](double x) -> double {
      double ret = NFit*IntfrBWARG( x,pp,0 );
      if ( FunSB == 0 ) {
         ret += NsbNorm/(dU-dL);
      } else {
         ret += NsbNorm * RevArgus(x,asb)*(1+sl*(x-1.02));
      }
      return ret;
   };
   double pvalueKS = FastKSTest(Fcr,mkk);
   printf(" pvalueKS(cr)= %.3g\n", pvalueKS);

   // side-band
   NsbNorm = 1. / my_fcn.normArsb; // norm to 1 for Lsb
   auto Fsb = [FunSB,NsbNorm,asb,sl](double x) -> double {
      double ret = 0;
      if ( FunSB == 0 ) {
         ret = 1./(dU-dL);
      } else {
         ret = NsbNorm * RevArgus(x,asb)*(1+sl*(x-1.02));
      }
      return ret;
   };
   double pvalueKSsb = FastKSTest(Fsb,sb);
   printf(" pvalueKS(sb)= %.3g\n", pvalueKSsb);

   //-----------------------------------------------------------------
   // Functions to draw
   auto Ldr = [Fpar,my_fcn,Ff,asb,FunSB](double* x,double* p) {
      int isw = int(p[0]+0.5);

      double Nphi = Fpar[1] * my_fcn.Br2Nphi;
      double NFit = Nphi/my_fcn.IBW;
      double ang = Fpar[2];
      double sig = Fpar[3];
      const double pp[] {
         my_fcn.mphi, my_fcn.gphi, sig, my_fcn.ar, Ff, ang, my_fcn.sl
      };

      double BWARG = NFit * IntfrBWARG( x[0],pp,isw );
      if ( isw > 0 && isw < 4 ) {
         return bW * BWARG;
      }

      double Bg = 0;
      if ( FunSB == 0 ) {
         Bg = Fpar[4]/(dU-dL);
      } else {
         double NsbNorm = Fpar[4] / my_fcn.normArsb;
         Bg = NsbNorm *
            RevArgus(x[0],asb)*(1+my_fcn.sl*(x[0]-1.02));
      }
      if ( isw == 4 ) {
         return bW * Bg;
      }

      return bW * (BWARG + Bg);
   };
   TF1* fdr = new TF1("fdr", Ldr, bL, dU, 1);
   fdr->SetNpx(500);

   // find max/min to draw
   fdr->SetParameter(0, 1); // BW
   int maxbin = hst->GetMaximumBin();
   double hmax = hst->GetBinContent(maxbin)+hst->GetBinError(maxbin);
   double fmax = fdr->GetMaximum( 1.01, 1.03, 1e-6);
   hmax = floor(1.1 * max( hmax, fmax ));
   if ( hmax > 0 ) {
      hst->SetMaximum(hmax);
   }

   fdr->SetParameter(0, 3); // interference
   double hmin = fdr->GetMinimum( 1.01, 1.03, 1e-6);
   hmin = floor(1.15 * min( 0., hmin ));
   if ( hmin < 0 ) {
      hst->SetMinimum(hmin);
   }

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   pdf += (neg_pos ? "_n.pdf" : "_p.pdf");
   c1->Print((pdf+"[").c_str()); // open pdf-file

   SetHstFace(hst);
   hst->GetXaxis()->SetTitleOffset(1.1);
   hst->GetYaxis()->SetTitleOffset(1.3);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");

   TLegend* leg = new TLegend(0.58,0.70,0.89,0.89);
   leg->AddEntry(hst,Form("Data %i",date),"LEP");

   fdr->SetParameter(0, 0); // SUM
   fdr->SetLineWidth(2);
   fdr->SetLineColor(kRed+1);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Result of fit", "L");

   fdr->SetParameter(0, 1); // BW
   fdr->SetLineWidth(2);
   fdr->SetLineStyle(kDashed);
   fdr->SetLineColor(kGreen+2);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Breit-Wigner #phi#eta", "L");

   fdr->SetParameter(0, 2); // Argus
   fdr->SetLineWidth(2);
   fdr->SetLineStyle(kDashed);
   fdr->SetLineColor(kBlue);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Non-#phi KK#eta", "L");

   fdr->SetParameter(0, 3); // interference
   fdr->SetLineWidth(2);
   fdr->SetLineStyle(kDashed);
   fdr->SetLineColor(kMagenta+1);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Interference", "L");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.58,0.47,0.89,0.69,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#it{p-value = %.3f}",pvalueKS) );
   pt->AddText( Form("Br(KK#eta)= %s #times10^{-4 }",
            PE.Eform(0,".2f",1e4)) );
   pt->AddText( Form("Br(#phi#eta)= %s #times10^{-4}",
            PE.Eform(1,".2f",1e4)) );
   pt->AddText( Form("#vartheta = %s",PE.Eform(2,".2f")) );
   pt->AddText( Form("#sigma = %s MeV", PE.Eform(3,".2f",1e3)) );
   pt->Draw();

   TPaveText* ptsb = new TPaveText(0.58,0.35,0.89,0.47,"NDC");
   ptsb->SetTextAlign(12);
   ptsb->SetTextFont(42);
   ptsb->AddText( Form("#it{p-val(side-band) = %.3f}", pvalueKSsb) );
   ptsb->AddText( Form("#lower[-0.1]{Nbg = %s}",PE.Eform(4,".1f")) );
if ( FunSB != 0 ) {
      ptsb->AddText( Form("#lower[-0.1]{a(sb) = %s}",
               PE.Eform(5,".1f")) );
   }
   ptsb->Draw();

   gPad->RedrawAxis();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   if ( !plotsb ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
      return;
   }

   //-----------------------------------------------------------------
   // Side-band

   SetHstFace(hsb);
   maxbin = hsb->GetMaximumBin();
   hmax = hsb->GetBinContent(maxbin)+hsb->GetBinError(maxbin);
   if ( hmax < 10 ) {
      hsb->SetMaximum( hmax+2);
   } else {
      hsb->SetMaximum( 1.2*hmax );
   }
   hsb->GetXaxis()->SetTitleOffset(1.1);
   hsb->GetYaxis()->SetTitleOffset(1.3);
   hsb->SetLineWidth(2);
   hsb->SetLineColor(kBlack);
   hsb->SetMarkerStyle(20);

   hsb->Draw("EP");

   TLegend* leg2 = new TLegend(0.55,0.81,0.89,0.89);
   leg2->AddEntry(hsb, Form("Side-band: data %i",date),"LEP");

   fdr->SetParameter(0, 4); // SB
   fdr->SetLineWidth(2);
   fdr->SetLineStyle(kSolid);
   fdr->SetLineColor(kRed+1);
   fdr->DrawCopy("SAME");
   leg2->AddEntry( fdr->Clone(), "Fit by Argus function", "L");
   leg2->Draw();

   ptsb->SetX1NDC(0.55);
   ptsb->SetX2NDC(0.89);
   ptsb->SetY1NDC(0.70);
   ptsb->SetY2NDC(0.80);
   ptsb->Draw();

   gPad->RedrawAxis();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
}

// {{{1 Combined fit
// {{{2   combine_Intfr()
//--------------------------------------------------------------------
vector<ValErr> CombIntEr( const ROOT::Fit::FitResult& res,
      double sl, unsigned int idx ) {
//--------------------------------------------------------------------
// Numerical calculation of integrals and their errors
   // desired errors:
   constexpr double eps_abs = 1.e-6;
   constexpr double eps_rel = 1.e-6;
   constexpr double epsilon = 1.e-3; // max numerical error

   // get parameters and covariation matrix of error
   int Npar = res.NPar();
   vector<double> Fp = res.Parameters();

   if ( Npar != 8 ) {
      cerr << " FATAL ERROR in " << __func__
         << " Npar= " << Npar << endl;
      exit(EXIT_FAILURE);
   }

   int covMatrStatus = res.CovMatrixStatus();
   if ( covMatrStatus != 3 ) {
      cout << " WARNING: " << __func__ << " covariance matrix"
         " status code is " << covMatrStatus << endl;
   }

   // cycle by year
   vector<ValErr> ret(3); // return value
   for ( size_t idat = 0; idat < 3; ++idat ) {
      // parameters for one year
      vector<double> Fpar {Fp[0],Fp[1],Fp[5+idat],Fp[2],Fp[3],Fp[4]};
      int npar = Fpar.size();
      vector<double> cov_m; // covariation matrix for Fpar
      cov_m.resize(npar*npar,0);
      // map raw of CovMatrix -> row of cov_m
      vector<int> MAP {0,1,3,4,5,  -1,-1,-1};
      MAP[5+idat] = 2;

      for ( int i = 0; i < Npar; ++i ) {
         int ii = MAP[i];
         for ( int j = 0; j < Npar; ++j ) {
            int jj = MAP[j];
            if ( ii != -1 && jj != -1 ) {
               cov_m[ii*npar+jj] = res.CovMatrix(ii,jj);
            }
         }
      }

      // 1) calculate integrals of normalized component
      auto Lfc = [sl,idx](const double* x,const double* p) {
         double m = x[0];
         double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl };
         double normI = IntfrBWARGN( m,pp,100 );
         pp[6] = 0; // set slope to zero
         double intfr = IntfrBWARG( m,pp,idx ) / normI;
         return intfr;
      };
      TF1 FC = TF1("FC",Lfc, dL, dU, npar);
      FC.SetParameters( Fpar.data() );

      double err_num = 0;
      auto& Integ = ret[idat].val;
      Integ = FC.IntegralOneDim(dL,dU,eps_rel,eps_abs, err_num);
      // cout << " err_num(" << idat << ") = " << err_num << endl;
      if ( err_num > epsilon ) {
         cout << " WARNING: " << __func__ << " numerical error"
            " of integration of FC-" << idat << " is too big "
            << err_num << endl;
      }

      // 2) calculate error of this integral
      TMatrixDSym covMatrix(npar);
      covMatrix.Use(npar,cov_m.data());
      // printf("\ncovMatrix-%d : ",idat);
      // covMatrix.Print();

      TVectorD IntGrad(npar);
      double err_num2 = 0;
      for ( int i = 0; i < npar; ++i ) {
         // skip parameters with zero error
         if ( covMatrix(i,i) == 0 ) {
            continue;
         }

         auto L_dF = [FC,i](const double* x, const double* p) {
            TF1 tmp(FC); // may modify 'tmp' but not FC
            return tmp.GradientPar(i,x);
         };
         TF1 dF("dF",L_dF,dL,dU,0);
         double err_num = 0;
         IntGrad[i] =
            dF.IntegralOneDim(dL,dU,eps_rel,eps_abs,err_num);
         err_num = covMatrix(i,i) * IntGrad[i] * err_num;
         err_num2 += SQ(err_num);
         // cout << " abs err_num(dF_" << i << ") = " << err_num << endl;
         // cout << " IntGrad-" << idat << "[" << i << "] = "
              // << IntGrad[i] << endl;
      }

      auto& err_Int = ret[idat].err;
      err_Int = sqrt( covMatrix.Similarity(IntGrad) );
      err_num2 = sqrt( err_num2 / err_Int ); // abs error
      // cout << " err_num(dBW-" << idat << ") = " << err_num2 << endl;
      if ( err_num2 > epsilon ) {
         cout << " WARNING: " << __func__ << " numerical error"
            " of integration of dF-" << idat << " is too big "
            << err_num2 << endl;
      }
   }
   return ret;
}

//--------------------------------------------------------------------
void combine_Intfr(string pdf, int neg_pos=-1) {
//--------------------------------------------------------------------
   // Get un-binned data
   int date[3] {2009, 2012, 2021};
   TH1D* hist[3];
   vector<double> mkk[3];
   size_t ntot = 0;
   for ( size_t i = 0; i < 3; ++i ) {
      string fname( Form("data_%02ipsip_all.root",date[i]%100) );
      mkk[i] = get_mkk_hist(fname,Form("mkk_%i_cp",date[i]),&hist[i]);
      ntot += mkk[i].size();
   }

   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Dat(dr, ntot);
   const double shift[3] { 10., 20., 0. };
   for ( size_t i = 0; i < 3; ++i ) {
      size_t nj = mkk[i].size();
      double sh = shift[i];
      const auto& mkki = mkk[i];
      for ( size_t j = 0; j < nj; ++j ) {
         Dat.Add( sh + mkki[j]);
      }
   }

   //-----------------------------------------------------------------
   // Fit data
   // efficiency parameters
   // double sl = -1.8; // prod-12
   double sl = -2.0; // v709

   // function MUST be normalized to 1 on the fit range
   auto Lfit = [sl](double* x, double* p) -> double {
      double pp[] = { p[0],p[1], p[7], p[2],p[3],p[4],sl };
      double m = x[0];
      if ( m > 10  ) {
         if ( m > 20 ) { // 2012
            m -= 20.;
            pp[2] = p[6];
         } else  {      // 2009
            m -= 10.;
            pp[2] = p[5];
         }
      }
      return IntfrBWARGN( m,pp,0 );
   };

   vector<string> par_name { "Mphi", "Gphi", "a", "F", "angle",
      "sig09", "sig12", "sig21" };

   vector<double> par_ini {Mphi, Gphi, 0., 1.0, 0.8,
      1.5e-3, 1.1e-3, 1.1e-3 };

   // neg_pos is -/+ for destructive/constructive interference
   if ( neg_pos > 0 ) {
      par_ini[4] = -par_ini[4];
   }

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 WFit( *Ffit );
   fitter.SetFunction(WFit,false); // false=no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   // fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   fitter.Config().ParSettings(0).SetValue(1.01952); // like MC-sig
   fitter.Config().ParSettings(0).Fix();                    // Mphi
   fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3, Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();                    // Gphi
   // fitter.Config().ParSettings(2).SetLimits(-5., 5.);       // a
   fitter.Config().ParSettings(2).Fix();                    // a
   fitter.Config().ParSettings(3).SetLimits(0.01, 10.);     // F
   fitter.Config().ParSettings(4).SetLimits(-M_PI, M_PI);   // angle
   fitter.Config().ParSettings(5).SetLimits(0.2e-3, 2.e-3); // sig09
   fitter.Config().ParSettings(6).SetLimits(0.2e-3, 2.e-3); // sig12
   fitter.Config().ParSettings(7).SetLimits(0.2e-3, 2.e-3); // sig21

   // == Fit
   fitter.LikelihoodFit(Dat,false); // true=extended likelihood fit
   fitter.CalculateHessErrors(); // for correct integral errors
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   // res.PrintCovMatrix(cout);

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------
   // Nphi and others
   vector<string> names { "NKK", "Nphi", "Nnonphi", "Nifr" };
   vector<ValErr> Nos[4];
   const auto& Nphi = Nos[1];
   string sepline(70,'='); // separation line
   printf("%s\n",sepline.c_str());
   for ( size_t i = 0; i < 3; ++i ) {
      printf(" n%02i = %zu", date[i]%100, mkk[i].size());
   }
   printf("\n");
   for ( size_t idx = 0; idx < names.size(); ++idx ) {
      vector<ValErr> Vinteg = CombIntEr( res,sl,idx );

      Nos[idx].resize(3); // for each year
      for ( size_t i = 0; i < 3; ++i ) {
         double norm = mkk[i].size();
         const auto& Integ = Vinteg[i];
         double Num = norm * Integ.val;
         double Err = fabs(Num) *
            sqrt( 1./norm + SQ(Integ.err/Integ.val) );
         Nos[idx][i] = ValErr {Num,Err};
         printf(" %7s_%02i = %s", names[idx].c_str(),date[i]%100,
               Nos[idx][i].prt("7.1f"));
      }
      printf("\n");
   }
   printf("%s\n",sepline.c_str());

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   vector<double> pvalueKS(3,0.);
   for ( size_t i = 0; i < 3; ++i ) {
      auto Fun = [Fpar,sl,i](double x) -> double {
         double pp[]
         {Fpar[0],Fpar[1],Fpar[5+i],Fpar[2],Fpar[3],Fpar[4], sl};
         return IntfrBWARGN( x,pp,0 );
      };
      pvalueKS[i] = FastKSTest(Fun,mkk[i]);
      printf(" pvalueKS_%i= %.3g\n", date[i],pvalueKS[i]);
   }
   printf("%s\n",sepline.c_str());

   //-----------------------------------------------------------------
   // Functions to draw
   std::function<double (double*, double*)> Ldr[3];
   vector<TF1*> ff(3,nullptr);
   for ( size_t i = 0; i < 3; ++i ) {
      auto bWn = bW*mkk[i].size();
      Ldr[i] = [bWn,Fpar,sl,i](double* x,double* p) -> double {
         double pp[]
         {Fpar[0],Fpar[1],Fpar[5+i],Fpar[2],Fpar[3],Fpar[4],sl};
         return bWn * IntfrBWARGN( x[0],pp, int(p[0]) );
      };

      ff[i] = new TF1(Form("ff%i",date[i]), Ldr[i], dL,dU,1);
      ff[i]->SetNpx(500);
   }

   // find max/min to draw
   for ( size_t i = 0; i < 3; ++i ) {
      auto& hst = hist[i];
      int mbin = hst->GetMaximumBin();
      double hmax = hst->GetBinContent(mbin)+hst->GetBinError(mbin);

      auto& fun = ff[i];
      fun->SetParameter(0, 1); // BW
      double fmax = fun->GetMaximum( 1.01, 1.03, 1e-6);

      hmax = floor(1.1 * max(hmax,fmax));
      if ( hmax > 0 ) {
         hst->SetMaximum(hmax);
      }

      fun->SetParameter(0, 3); // interference
      double hmin = fun->GetMinimum( 1.01, 1.03, 1e-6);
      hmin = floor(1.15 * min(0.,hmin));
      if ( hmin < 0 ) {
         hst->SetMinimum(hmin);
      }
   }

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,700,1000);
   c1->Divide(1,3);

   TLegend* leg = new TLegend(0.11,0.50,0.35,0.89);
   vector<TPaveText*> pt(3,nullptr);

   for ( size_t i = 0; i < 3; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();

      auto& hst = hist[i];
      auto& fun = ff[i];

      SetHstFace(hst);
      hst->GetXaxis()->SetTitleOffset(1.1);
      hst->GetYaxis()->SetTitleOffset(1.0);
      hst->SetLineWidth(2);
      hst->SetLineColor(kBlack);
      hst->SetMarkerStyle(20);
      hst->SetMarkerSize(0.8);
      hst->Draw("EP");
      if ( i == 0 ) {
         leg->AddEntry( hist[0],"Data","LEP" );
      }

      fun->SetParameter(0, 0); // Sum
      fun->SetLineWidth(2);
      fun->SetLineColor(kRed+1);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         leg->AddEntry( fun->Clone(), "Combined fit", "L" );
      }

      fun->SetParameter(0, 1); // BW
      fun->SetLineWidth(2);
      fun->SetLineStyle(kDashed);
      fun->SetLineColor(kGreen+3);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         leg->AddEntry( fun->Clone(), "Breit-Wigner #phi#eta", "L");
      }

      fun->SetParameter(0, 2); // Argus
      fun->SetLineWidth(1);
      fun->SetLineColor(kBlue);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         leg->AddEntry( fun->Clone(), "Non-#phi KK#eta", "L");
      }

      fun->SetParameter(0, 3); // interference
      fun->SetLineWidth(1);
      fun->SetLineColor(kMagenta+1);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         leg->AddEntry( fun->Clone(), "Interference", "L");
      }

      pt[i] = new TPaveText(0.61,0.72,0.89,0.89,"NDC");
      pt[i]->SetTextAlign(12);
      pt[i]->SetTextFont(42);
      pt[i]->AddText( Form("#it{p-value(%i) = %.3f}",
               date[i],pvalueKS[i]) );
      pt[i]->AddText( Form("#sigma(%i) = %s MeV",
               date[i],PE.Eform(5+i,".2f",1e3)) );
      pt[i]->AddText( Form("N_{#phi}(%i) = %.1f #pm %.1f",
               date[i],Nphi[i].val,Nphi[i].err) );
   }

   TPaveText* ptc = new TPaveText(0.61,0.40,0.89,0.72,"NDC");
   ptc->SetTextAlign(12);
   ptc->SetTextFont(42);
   ptc->AddText( "#it{common parameters}");
   ptc->AddText( Form("#lower[0.1]{M_{#phi}= %s MeV}",
            PE.Eform(0,".2f",1e3)) );
   ptc->AddText( Form("#Gamma_{#phi}= %s MeV",
            PE.Eform(1,".3f",1e3)) );
   ptc->AddText( Form("a = %s",PE.Eform(2,".2f")) );
   ptc->AddText( Form("F = %s",PE.Eform(3,".2f")) );
   ptc->AddText( Form("#vartheta = %s",PE.Eform(4,".2f")) );

   c1->cd(1);
   leg->Draw();
   pt[0]->Draw();
   ptc->Draw();
   gPad->RedrawAxis();

   c1->cd(2);
   pt[1]->Draw();
   gPad->RedrawAxis();

   c1->cd(3);
   pt[2]->Draw();
   gPad->RedrawAxis();

   c1->Update();
   if ( !pdf.empty() ) {
      pdf += (neg_pos<0 ? "_n.pdf" : "_p.pdf");
      c1->Print(pdf.c_str());
   }
}

// {{{2 combineSB_Intfr(): combined plus side-band
//--------------------------------------------------------------------
struct myFCN_comb {
   // efficiency parameters
   // const double sl = -1.8; // prod-12
   const double sl = -2.0; // v709
   // Function for side-band:
   const int funSB = 0; // 0 - 'asb' for each year are independent
                        // 1 - one 'asb' for all years

   const vector<int> date { 2009, 2012, 2021 };
   const size_t ndate = 3;

   vector<double> mkk[3]; // data central part
   vector<double> sb[3];  // data side band

   // minimization function
   // the signature of this operator() MUST be exactly this:
   // THE EXTENDED MAXIMUM LIKELIHOOD ESTIMATOR
   // Kai-Feng Chen "FittingMinuitRooFit" p.35
   //-----------------------------------------------------------------
   double operator() (const double* p) {
   //-----------------------------------------------------------------
      double res = 0;
      for ( size_t idat = 0; idat < ndate; ++idat ) {
         size_t idx = 5 + 4*idat;
         double sig = p[idx];
         double Nkk = p[idx+1];
         double Nsb = p[idx+2];
         double asb = p[idx+3];
         if ( funSB == 1 ) {
            asb = p[5+3];
         }

         res += 2*(Nkk+Nsb);
         const double pp[] = { p[0],p[1],sig,p[2],p[3],p[4],sl };
         for ( const auto& m : mkk[idat] ) {
            double L = Nkk * IntfrBWARGN( m,pp,0 );
            L += Nsb * RevArgusN(m,asb,sl);
            if (L > 0.) {
               res -= 2*log(L);
            } else {
               res += FLT_MAX; // ~3e38
            }
         }

         // fit SB
         res += 2*Nsb;
         for ( const auto& m : sb[idat] ) {
            double L = Nsb * RevArgusN(m,asb,sl);
            if (L > 0.) {
               res -= 2*log(L);
            } else {
               res += FLT_MAX; // ~3e38
            }
         }
      }

      return res;
   }
};

//--------------------------------------------------------------------
vector<ValErr> CombIntErSB(const ROOT::Fit::FitResult& res,
      double sl, unsigned int idx, bool noerrors=false) {
//--------------------------------------------------------------------
// Numerical calculation of integrals and their errors
   // desired errors:
   constexpr double eps_abs = 1.e-6;
   constexpr double eps_rel = 1.e-6;
   constexpr double epsilon = 1.e-3; // max numerical error

   // get parameters and covariation matrix of error
   int Npar = res.NPar();
   vector<double> Fp = res.Parameters();

   if ( Npar != 17 ) {
      cerr << " FATAL ERROR in " << __func__
         << " Npar= " << Npar << endl;
      exit(EXIT_FAILURE);
   }

   int covMatrStatus = res.CovMatrixStatus();
   if ( covMatrStatus != 3 ) {
      cout << " WARNING: " << __func__ << " covariance matrix"
         " status code is " << covMatrStatus << endl;
   }

   // cycle by year
   const char* sdat[] { "09","12","21" }; // for printf
   vector<ValErr> ret(3); // return value
   for ( size_t idat = 0; idat < 3; ++idat ) {
      // parameters for one year
      double sig = Fp[5 + idat*4];
      vector<double> Fpar {Fp[0],Fp[1],sig,Fp[2],Fp[3],Fp[4]};
      int npar = Fpar.size();
      vector<double> cov_m; // covariation matrix for Fpar
      cov_m.resize(npar*npar,0);
      // map raw of CovMatrix -> row of cov_m
      vector<int> MAP {0,1,3,4,5,
         -1,-1,-1,-1,
         -1,-1,-1,-1,
         -1,-1,-1,-1
      };
      MAP[5+idat*4] = 2;

      for ( int i = 0; i < Npar; ++i ) {
         int ii = MAP[i];
         for ( int j = 0; j < Npar; ++j ) {
            int jj = MAP[j];
            if ( ii != -1 && jj != -1 ) {
               cov_m[ii*npar+jj] = res.CovMatrix(ii,jj);
            }
         }
      }

      // 1) calculate integrals of normalized component
      auto Lfc = [sl,idx](const double* x,const double* p) {
         double m = x[0];
         double pp[] { p[0],p[1],p[2],p[3],p[4],p[5],sl };
         double normI = IntfrBWARGN( m,pp,100 );
         pp[6] = 0; // set slope to zero
         double intfr = IntfrBWARG( m,pp,idx ) / normI;
         return intfr;
      };
      TF1 FC = TF1("FC",Lfc, dL, dU, npar);
      FC.SetParameters( Fpar.data() );

      double err_num = 0;
      auto& Integ = ret[idat].val;
      Integ = FC.IntegralOneDim(dL,dU,eps_rel,eps_abs, err_num);
      // printf(" Integ%s[%d] = %.6g; err_num= %.2g\n",
            // sdat[idat],idx,Integ,err_num);
      if ( err_num > epsilon ) {
         printf(" WARNING: %s numerical error of integration of "
               "FC%s[%d] is too big: %.2g\n",
               __func__,sdat[idat],idx,err_num);
      }

      if ( noerrors ) {
         continue; // skip error calculation
      }

      // 2) calculate error of this integral
      TMatrixDSym covMatrix(npar);
      covMatrix.Use(npar,cov_m.data());
      // printf("\ncovMatrix-%s : ",sdat[idat]);
      // covMatrix.Print();

      TVectorD IntGrad(npar);
      double err_num2 = 0;
      for ( int i = 0; i < npar; ++i ) {
         // skip parameters with zero error
         if ( covMatrix(i,i) == 0 ) {
            continue;
         }

         auto L_dF = [FC,i](const double* x, const double* p) {
            TF1 tmp(FC); // may modify 'tmp' but not FC
            return tmp.GradientPar(i,x);
         };
         TF1 dF("dF",L_dF,dL,dU,0);
         double err_num = 0;
         IntGrad[i] =
            dF.IntegralOneDim(dL,dU,eps_rel,eps_abs,err_num);
         err_num = covMatrix(i,i) * IntGrad[i] * err_num;
         err_num2 += SQ(err_num);
         // printf(" IntGrad(dF%s[%d]/d_%d)= %.4g; err_num= %.4g\n",
               // sdat[idat],idx,i, IntGrad[i], err_num);
      }

      auto& err_Int = ret[idat].err;
      err_Int = sqrt( covMatrix.Similarity(IntGrad) );
      err_num2 = sqrt( err_num2 / err_Int ); // abs error
      // printf(" --> Integral_%s[idx=%d] = %.6g +/- %.6g\n",
            // sdat[idat],idx,Integ,err_Int);
      if ( err_num2 > epsilon ) {
         printf(" WARNING: %s numerical error of integration of "
               "dF_%s[%d] is too big: %.2g\n",
               __func__,sdat[idat],idx,err_num2);
      }
   }
   return ret;
}

//--------------------------------------------------------------------
void combineSB_Intfr(string pdf, int neg_pos=-1, bool plotsb=true) {
//--------------------------------------------------------------------
   myFCN_comb my_fcn;              // class for 'FitFCN'
   const double sl = my_fcn.sl;    // efficiency parameters
   const int FunSB = my_fcn.funSB; // Function for side-band
   auto& date      = my_fcn.date;

   // Get un-binned data
   TH1D* hist[6];
   size_t ntot = 0;
   for ( size_t i = 0; i < 3; ++i ) {
      string fname( Form("data_%02ipsip_all.root",date[i]%100) );
      my_fcn.mkk[i] =
         get_mkk_hist(fname,Form("mkk_%i_cp",date[i]),&hist[i]);
      ntot += my_fcn.mkk[i].size();
      my_fcn.sb[i] =
         get_mkk_hist(fname,Form("mkk_%i_sb",date[i]),&hist[3+i], 1);
      ntot += my_fcn.sb[i].size();
   }

   //-----------------------------------------------------------------
   // Fit data
   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   vector<string> par_name {
      "Mphi", "Gphi", "a", "F", "angle",
      "sig09", "Nkk09", "Nsb09", "asb09",
      "sig12", "Nkk12", "Nsb12", "asb12",
      "sig21", "Nkk21", "Nsb21", "asb21"
   };

   vector<double> par_ini {
      Mphi, Gphi, 0., 0.85, 0.5,
      1.5e-3,   890.,   8., 7.,
      1.1e-3,  2750.,  33., 4.,
      1.1e-3, 17500., 195., 6.
   };

   // neg_pos is -/+ for destructive/constructive interference
   if ( neg_pos > 0 ) {
      par_ini[4] = -par_ini[4];
   }

   if ( FunSB == 1 ) {
      par_ini[8] = par_ini[12] = par_ini[16]; // common asb
   }

   const unsigned int Npar = par_name.size(); // number of parameters

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetValue(1.01952); // like MC-sig
   fitter.Config().ParSettings(0).Fix();                    // Mphi
   fitter.Config().ParSettings(1).Fix();                    // Gphi
   fitter.Config().ParSettings(2).Fix();                    // a
   fitter.Config().ParSettings(3).SetLimits(0.01, 10.);     // F
   fitter.Config().ParSettings(4).SetLimits(-M_PI, M_PI);   // angle

   for ( size_t i = 0; i < 3; ++i ) { // sig,Nkk,Nsb,asb per each year
      fitter.Config().ParSettings(5+4*i).SetLimits(0.2e-3, 5.e-3);
      // fitter.Config().ParSettings(6+4*i).SetLimits(0.,2.*Nmkk[i]);
      // fitter.Config().ParSettings(7+4*i).SetLimits(0.,1.5*Nsb[i]);
      fitter.Config().ParSettings(8+4*i).SetLimits(0.,15.);
      if ( FunSB == 1 && i > 0 ) {
         fitter.Config().ParSettings(8+4*i).Fix(); // asb
      }
   }

   // == Fit
   int Ndat = ntot;
   fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // false=likelihood
   fitter.CalculateMinosErrors();
   fitter.CalculateHessErrors(); // for calculation of errors of Nphi

   ROOT::Fit::FitResult res = fitter.Result();
   vector<double> Fpar = res.Parameters();
   if ( res.IsValid() ) {
      printf("\n=> FINAL RESULT <=\n");
   } else {
      printf("\n=> INVALID RESULT <=\n");
   }
   res.Print(cout);
   // res.PrintCovMatrix(cout);

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res,0.05); // ignore 5% upper/lower errors

   //-----------------------------------------------------------------
   // Nphi and others
   vector<string> names { "NKK", "Nphi", "Nnonphi", "Nifr" };
   vector<ValErr> Nos[4];
   const auto& Nekk = Nos[0];
   const auto& Nphi = Nos[1];
   string sepline(70,'='); // separation line
   printf("%s\n",sepline.c_str());

   // printf( " no slope corrections: ");
   printf("%21i%21i%21i\n",date[0],date[1],date[2]);
   printf("%8s","NKK_fit");
   for ( size_t i = 0; i < 3; ++i ) {
      size_t iy = 6+i*4;
      double Nkk = Fpar[iy];
      double Nkk_err = PE.Perr[iy];
      if ( PE.Uerr[iy]*PE.Lerr[iy] != 0 ) {
         Nkk_err = max(PE.Uerr[iy],PE.Lerr[iy]);
      }
      printf( "  %7.1f +/- %7.1f", Nkk,Nkk_err);
   }
   printf("\n");

   for ( size_t idx = 0; idx < names.size(); ++idx ) {
      vector<ValErr> Vinteg = CombIntErSB( res,sl,idx );

      printf("%8s",names[idx].c_str());
      Nos[idx].resize(3); // for each year
      for ( size_t i = 0; i < 3; ++i ) {
         size_t iy = 6+i*4;
         double Nkk = Fpar[iy];
         double Nkk_err = PE.Perr[iy];
         if ( PE.Uerr[iy]*PE.Lerr[iy] != 0 ) {
            Nkk_err = max(PE.Uerr[iy],PE.Lerr[iy]);
         }

         const auto& Integ = Vinteg[i];
         // printf("Int%02i[%zu] = %s\n",date[i],idx,Integ.prt(".6f"));
         double Num = Nkk * Integ.val;
         double Err = fabs(Num) *
            sqrt( SQ(Integ.err/Integ.val) + SQ(Nkk_err/Nkk) );
         Nos[idx][i] = ValErr {Num,Err};
         printf("  %s",Nos[idx][i].prt("7.1f"));
      }
      printf("\n");
   }
   printf("%s\n",sepline.c_str());

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   vector<double> pvalueKS(3,0.), pvalueKSsb(3,0.);
   for ( size_t i = 0; i < 3; ++i ) {
      // central part
      auto Fcr = [Fpar,sl,FunSB,i](double x) -> double {
         size_t iy = 5 + 4*i;
         double sig = Fpar[iy];
         double Nkk = Fpar[iy+1];
         double Nsb = Fpar[iy+2];
         double asb = Fpar[iy+3];
         if ( FunSB == 1 ) {
            asb = Fpar[5+3];
         }
         double pp[] {Fpar[0],Fpar[1],sig,Fpar[2],Fpar[3],Fpar[4],sl};
         double ret = Nkk * IntfrBWARGN( x,pp,0 );
         ret += Nsb * RevArgusN( x,asb,sl );
         return ret/(Nkk+Nsb); // normalization on 1
      };
      pvalueKS[i] = FastKSTest(Fcr,my_fcn.mkk[i]);

      // side-band
      auto Fsb = [Fpar,sl,FunSB,i](double x) -> double {
         size_t iy = 5 + 4*i;
         double asb = Fpar[iy+3];
         if ( FunSB == 1 ) {
            asb = Fpar[5+3];
         }
         return RevArgusN( x,asb,sl );
      };
      pvalueKSsb[i] = FastKSTest(Fsb,my_fcn.sb[i]);
      printf("%i: pvalueKS(cr)= %.3f pvalueKS(sb)= %.3f\n",
            date[i],pvalueKS[i],pvalueKSsb[i]);
   }
   printf("%s\n",sepline.c_str());

   //-----------------------------------------------------------------
   // Functions to draw
   std::function<double (double*, double*)> Ldr[3];
   vector<TF1*> ff(3,nullptr);
   for ( size_t i = 0; i < 3; ++i ) {
      Ldr[i] = [Fpar,sl,FunSB,i](double* x,double* p) -> double {
         int isw = int(p[0]+0.5);
         size_t iy = 5 + 4*i;
         double sig = Fpar[iy];
         double Nkk = Fpar[iy+1];
         double Nsb = Fpar[iy+2];
         double asb = Fpar[iy+3];
         if ( FunSB == 1 ) {
            asb = Fpar[5+3];
         }
         double pp[] {Fpar[0],Fpar[1],sig,Fpar[2],Fpar[3],Fpar[4],sl};
         double BWARG = Nkk * IntfrBWARGN(x[0],pp,isw);
         if ( isw > 0 && isw < 4 ) {
            return bW * BWARG;
         }
         double Bg = Nsb * RevArgusN( x[0],asb,sl );
         if ( isw == 4 ) {
            return bW * Bg;
         }
         return bW * (BWARG + Bg);
      };
      ff[i] = new TF1(Form("ff%i",date[i]), Ldr[i], dL,dU,1);
      ff[i]->SetNpx(500);
   }

   // find max/min to draw
   for ( size_t i = 0; i < 3; ++i ) {
      auto& hst = hist[i];
      int mbin = hst->GetMaximumBin();
      double hmax = hst->GetBinContent(mbin)+hst->GetBinError(mbin);

      auto& fun = ff[i];
      fun->SetParameter(0, 1); // BW
      double fmax = fun->GetMaximum( 1.01, 1.03, 1e-6);

      hmax = floor(1.1 * max(hmax,fmax));
      if ( hmax > 0 ) {
         hst->SetMaximum(hmax);
      }

      fun->SetParameter(0, 3); // interference
      double hmin = fun->GetMinimum( 1.01, 1.03, 1e-6);
      hmin = floor(1.15 * min(0.,hmin));
      if ( hmin < 0 ) {
         hst->SetMinimum(hmin);
      }

      if ( plotsb ) { // Side-band
         auto& hsb = hist[3+i];
         mbin = hsb->GetMaximumBin();
         hmax = hsb->GetBinContent(mbin)+hsb->GetBinError(mbin);
         if ( hmax < 10 ) {
            hsb->SetMaximum( hmax+2);
         } else {
            hsb->SetMaximum( 1.2*hmax );
         }
      }
   }

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000);
   c1->Divide(1,3);

   pdf += (FunSB == 1) ? "_ar1" : "_ar";
   pdf += (neg_pos<0 ? "_n.pdf" : "_p.pdf");
   c1->Print((pdf+"[").c_str()); // open pdf-file

   TLegend* leg = new TLegend(0.11,0.50,0.35,0.89);
   vector<TPaveText*> pt(3,nullptr), ptr(3,nullptr);

   for ( size_t i = 0; i < 3; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();

      auto& hst = hist[i];
      auto& fun = ff[i];

      SetHstFace(hst); // Data
      hst->GetXaxis()->SetTitleOffset(1.1);
      hst->GetYaxis()->SetTitleOffset(1.0);
      hst->SetLineWidth(2);
      hst->SetLineColor(kBlack);
      hst->SetMarkerStyle(20);
      hst->SetMarkerSize(0.8);
      hst->Draw("EP");
      if ( i == 0 ) {
         leg->AddEntry( hist[0],"Data","LEP" );
      }

      fun->SetParameter(0, 0); // Sum
      fun->SetLineWidth(2);
      fun->SetLineStyle(kSolid);
      fun->SetLineColor(kRed+1);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         leg->AddEntry( fun->Clone(), "Combined fit", "L" );
      }

      fun->SetParameter(0, 1); // BW
      fun->SetLineWidth(2);
      fun->SetLineStyle(kDashed);
      fun->SetLineColor(kGreen+2);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         leg->AddEntry( fun->Clone(), "Breit-Wigner #phi#eta", "L");
      }

      fun->SetParameter(0, 2); // Argus
      fun->SetLineWidth(1);
      fun->SetLineStyle(kDashed);
      fun->SetLineColor(kBlue);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         leg->AddEntry( fun->Clone(), "Non-#phi KK#eta", "L");
      }

      fun->SetParameter(0, 3); // interference
      fun->SetLineWidth(1);
      fun->SetLineStyle(kDashed);
      fun->SetLineColor(kMagenta+1);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         leg->AddEntry( fun->Clone(), "Interference", "L");
      }

      double yptm = 0.59 + 0.06*(FunSB==1);
      pt[i] = new TPaveText(0.61,yptm,0.89,0.89,"NDC");
      pt[i]->SetTextAlign(12);
      pt[i]->SetTextFont(42);
      pt[i]->AddText( Form("#it{p-value(%i) = %.3f}",
               date[i],pvalueKS[i]) );
      pt[i]->AddText( Form("#sigma(%i) = %s MeV",
               date[i],PE.Eform(5+i*4,".2f",1e3)) );
      pt[i]->AddText( Form("N_{KK}(%i) = %s",
               date[i],PE.Eform(6+i*4,".1f")) );
      pt[i]->AddText( Form("N_{sb}(%i) = %s",
               date[i],PE.Eform(7+i*4,".1f")) );
      if ( FunSB == 0 ) {
         pt[i]->AddText( Form("a_{sb}(%i) = %s",
                  date[i],PE.Eform(8+i*4,".1f")) );
      }

      ptr[i] = new TPaveText(0.61,yptm-0.12,0.89,yptm,"NDC");
      ptr[i]->SetTextAlign(12);
      ptr[i]->SetTextFont(42);
      ptr[i]->AddText( Form("N_{#etaKK}(%i) = %.1f #pm %.1f",
               date[i],Nekk[i].val,Nekk[i].err) );
      ptr[i]->AddText( Form("N_{#phi}(%i) = %.1f #pm %.1f",
               date[i],Nphi[i].val,Nphi[i].err) );
   }

   double yptm = 0.47 + 0.06*(FunSB==1);
   TPaveText* ptc = new TPaveText(0.61,0.29,0.89,yptm,"NDC");
   ptc->SetTextAlign(12);
   ptc->SetTextFont(42);
   ptc->AddText( "#it{common parameters}");
   // ptc->AddText( Form("#lower[0.1]{M_{#phi}= %s MeV}",
            // PE.Eform(0,".2f",1e3)) );
   // ptc->AddText( Form("#Gamma_{#phi}= %s MeV",
            // PE.Eform(1,".3f",1e3)) );
   // ptc->AddText( Form("a = %s",PE.Eform(2,".2f")) );
   ptc->AddText( Form("F = %s",PE.Eform(3,".2f")) );
   ptc->AddText( Form("#vartheta = %s",PE.Eform(4,".2f")) );
   if ( FunSB == 1 ) {
      ptc->AddText( Form("a_{sb} = %s", PE.Eform(8,".1f")) );
   }

   c1->cd(1);
   leg->Draw();
   pt[0]->Draw();
   ptr[0]->Draw();
   ptc->Draw();
   gPad->RedrawAxis();

   c1->cd(2);
   pt[1]->Draw();
   ptr[1]->Draw();
   gPad->RedrawAxis();

   c1->cd(3);
   pt[2]->Draw();
   ptr[2]->Draw();
   gPad->RedrawAxis();

   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   if ( !plotsb ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
      return;
   }

   //-----------------------------------------------------------------
   // Side-band
   TLegend* legsb = new TLegend(0.11,0.72,0.40,0.89);
   vector<TPaveText*> ptsb(3,nullptr);

   for ( size_t i = 0; i < 3; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();

      auto& hst = hist[3+i];
      auto& fun = ff[i];

      SetHstFace(hst);
      hst->GetXaxis()->SetTitleOffset(1.1);
      hst->GetYaxis()->SetTitleOffset(1.3);
      hst->SetLineWidth(2);
      hst->SetLineColor(kBlack);
      hst->SetMarkerStyle(20);
      hst->Draw("EP");
      if ( i == 0 ) {
         legsb->AddEntry(hst,"Side-band data","LEP");
      }

      fun->SetParameter(0, 4); // SB
      fun->SetLineWidth(2);
      fun->SetLineStyle(kSolid);
      fun->SetLineColor(kRed+1);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         legsb->AddEntry(fun->Clone(),"Fit by Argus function:","L");
      }

      ptsb[i] = new TPaveText(0.70,0.61,0.89,0.89,"NDC");
      ptsb[i]->SetTextAlign(12);
      ptsb[i]->SetTextFont(42);
      ptsb[i]->AddText( Form("#bf{%16i}", date[i]) );
      ptsb[i]->AddText( Form("p-value = %.3f",pvalueKSsb[i]) );
      ptsb[i]->AddText( Form("N_{sb} = %s", PE.Eform(7+i*4,".1f")) );
      size_t ia = 8 + (i*4)*(FunSB==0);
      ptsb[i]->AddText( Form("#lower[-0.01]{a_{sb} = %s}",
               PE.Eform(ia,".1f")) );
   }

   c1->cd(1);
   legsb->Draw();
   ptsb[0]->Draw();
   gPad->RedrawAxis();

   c1->cd(2);
   ptsb[1]->Draw();
   gPad->RedrawAxis();

   c1->cd(3);
   ptsb[2]->Draw();
   gPad->RedrawAxis();

   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
}

// {{{2   combineSB_Intfr_scan()
//--------------------------------------------------------------------
void combineSB_Intfr_scan(string sname) {
//--------------------------------------------------------------------

   myFCN_comb my_fcn;              // class for 'FitFCN'
   const double sl = my_fcn.sl;    // efficiency parameters
   const int FunSB = my_fcn.funSB; // Function for side-band
   auto& date      = my_fcn.date;

   // Get un-binned data
   TH1D* hist[6];
   size_t ntot = 0;
   for ( size_t i = 0; i < 3; ++i ) {
      string fname( Form("data_%02ipsip_all.root",date[i]%100) );
      my_fcn.mkk[i] =
         get_mkk_hist(fname,Form("mkk_%i_cp",date[i]),&hist[i]);
      ntot += my_fcn.mkk[i].size();
      my_fcn.sb[i] =
         get_mkk_hist(fname,Form("mkk_%i_sb",date[i]),&hist[3+i], 1);
      ntot += my_fcn.sb[i].size();
   }

   int Ndat = ntot;
   //-----------------------------------------------------------------
   // Fit data
   ROOT::Fit::Fitter fitter;
   // fitter.Config().MinimizerOptions().SetPrintLevel(3);

   vector<string> par_name {
      "Mphi", "Gphi", "a", "F", "angle",
      "sig09", "Nkk09", "Nsb09", "asb09",
      "sig12", "Nkk12", "Nsb12", "asb12",
      "sig21", "Nkk21", "Nsb21", "asb21"
   };

   vector<double> par_ini {
      Mphi, Gphi, 0., 1.5,  0.,
      1.5e-3,   890.,   8., 7.,
      1.1e-3,  2750.,  33., 4.,
      1.1e-3, 17500., 195., 6.
   };

   if ( FunSB == 1 ) {
      par_ini[8] = par_ini[12] = par_ini[16]; // common
   }

   // parameters of loop for angles (degrees)
   // int angmin = -90, angmax = 90, angstep = 181; // test
   int angmin = -90, angmax = 90, angstep = 1; // to file

   const unsigned int Npar = par_name.size(); // number of parameters

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetValue(1.01952); // like MC-sig
   fitter.Config().ParSettings(0).Fix();                    // Mphi
   fitter.Config().ParSettings(1).Fix();                    // Gphi
   fitter.Config().ParSettings(2).Fix();                    // Ar
   fitter.Config().ParSettings(3).SetLimits(0.01, 10.);     // F

   for ( size_t i = 0; i < 3; ++i ) { // sig,Nkk,Nbg,asb per each year
      fitter.Config().ParSettings(5+4*i).SetLimits(0.2e-3, 5.e-3);
      // fitter.Config().ParSettings(6+4*i).SetLimits(0.,2.*Nmkk[i]);
      // fitter.Config().ParSettings(7+4*i).SetLimits(0.,1.5*Nsb[i]);
      fitter.Config().ParSettings(8+4*i).SetLimits(0.,15.);
      if ( FunSB == 1 && i > 0 ) {
         fitter.Config().ParSettings(8+4*i).Fix(); // asb
      }
   }

   //-----------------------------------------------------------------
   // Functions to draw
   vector<double> Fpar;
   std::function<double (double*, double*)> Ldr[3];
   vector<TF1*> ff(3,nullptr);
   for ( size_t i = 0; i < 3; ++i ) {
      Ldr[i] = [&Fpar,sl,FunSB,i](double* x,double* p) -> double {
         int isw = int(p[0]+0.5);
         size_t iy = 5 + 4*i;
         double sig = Fpar[iy];
         double Nkk = Fpar[iy+1];
         double Nsb = Fpar[iy+2];
         double asb = Fpar[iy+3];
         if ( FunSB == 1 ) {
            asb = Fpar[5+3];
         }
         double pp[] {Fpar[0],Fpar[1],sig,Fpar[2],Fpar[3],Fpar[4],sl};
         double BWARG = Nkk * IntfrBWARGN(x[0],pp,isw);
         if ( isw > 0 && isw < 4 ) {
            return bW * BWARG;
         }
         double Bg = Nsb * RevArgusN( x[0],asb,sl );
         if ( isw == 4 ) {
            return bW * Bg;
         }
         return bW * (BWARG + Bg);
      };
      ff[i] = new TF1(Form("ff%i",date[i]), Ldr[i], dL,dU,1);
      ff[i]->SetNpx(500);
   }

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000);
   c1->Divide(1,3);

   sname += (FunSB == 1) ? "_ar1" : "_ar";
   string hfile = sname + ".root";
   TFile c_out(hfile.c_str(),"recreate"); // file for scan results
   c_out.cd();

   string pdf = sname + ".pdf";
   c1->Print((pdf+"[").c_str()); // open pdf-file

   vector<double> hmin { -15., -35., -200.};
   for ( size_t i = 0; i < 3; ++i ) {
      auto& hst = hist[i];

      SetHstFace(hst); // Data
      hst->GetXaxis()->SetTitleOffset(1.1);
      hst->GetYaxis()->SetTitleOffset(1.0);
      hst->SetLineWidth(2);
      hst->SetLineColor(kBlack);
      hst->SetMarkerStyle(20);
      hst->SetMarkerSize(0.8);
      hst->SetMinimum(hmin[i]);
   }

   //-----------------------------------------------------------------
   // start the cycle
   vector<double> angl, Lmin, Nphi[3];
   for( int ang = angmin; ang <= angmax; ang+=angstep ) {

      fitter.Config().ParSettings(4).SetValue(M_PI*ang/180.); // angle
      fitter.Config().ParSettings(4).Fix();

      // == Fit
      fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false);//false=likelihood
      // fitter.CalculateMinosErrors();
      // fitter.CalculateHessErrors();

      ROOT::Fit::FitResult res = fitter.Result();
      Fpar = res.Parameters();
      if ( angstep > 9 ) {
         res.Print(cout);
      }
      // res.PrintCovMatrix(cout);
      ParAndErr PE(res);

      double lmin = res.MinFcnValue();
      angl.push_back(ang);
      Lmin.push_back(lmin);

      // calculate Nphi only and without error
      string sepline(50,'='); // separation line
      printf("%s\n",sepline.c_str());
      printf("ang= %i degree =>  Lmin= %.2f\n",ang,lmin);
      vector<ValErr> Vinteg = CombIntErSB( res,sl,1,true );
      for ( size_t i = 0; i < 3; ++i ) {
         size_t iy = 6+i*4;
         double Nkk = Fpar[iy];
         const auto& Integ = Vinteg[i];
         double Num = Nkk * Integ.val;
         Nphi[i].push_back(Num);
         printf("Nphi(%02i)= %.1f, ",date[i]%100,Num);
      }
      printf("\n");
      printf("%s\n",sepline.c_str());

      //--------------------------------------------------------------
      // Draw intermidiate results

      TLegend* leg = new TLegend(0.11,0.50,0.35,0.89);
      vector<TPaveText*> pt(3,nullptr);

      for ( size_t i = 0; i < 3; ++i ) {
         c1->cd(i+1);
         gPad->SetGrid();

         auto& hst = hist[i];
         auto& fun = ff[i];

         hst->Draw("EP");
         if ( i == 0 ) {
            leg->AddEntry( hist[0],"Data","LEP" );
         }

         fun->SetParameter(0, 0); // Sum
         fun->SetLineWidth(2);
         fun->SetLineStyle(kSolid);
         fun->SetLineColor(kRed+1);
         fun->DrawCopy("SAME");
         if ( i == 0 ) {
            leg->AddEntry( fun->Clone(), "Combined fit", "L" );
         }

         fun->SetParameter(0, 1); // BW
         fun->SetLineWidth(2);
         fun->SetLineStyle(kDashed);
         fun->SetLineColor(kGreen+2);
         fun->DrawCopy("SAME");
         if ( i == 0 ) {
            leg->AddEntry(fun->Clone(), "Breit-Wigner #phi#eta", "L");
         }

         fun->SetParameter(0, 2); // Argus
         fun->SetLineWidth(1);
         fun->SetLineStyle(kDashed);
         fun->SetLineColor(kBlue);
         fun->DrawCopy("SAME");
         if ( i == 0 ) {
            leg->AddEntry( fun->Clone(), "Non-#phi KK#eta", "L");
         }

         fun->SetParameter(0, 3); // interference
         fun->SetLineWidth(1);
         fun->SetLineStyle(kDashed);
         fun->SetLineColor(kMagenta+1);
         fun->DrawCopy("SAME");
         if ( i == 0 ) {
            leg->AddEntry( fun->Clone(), "Interference", "L");
         }

         // pt[i] = new TPaveText(0.61,0.59,0.89,0.89,"NDC");
         pt[i] = new TPaveText(0.61,0.40,0.89,0.70,"NDC");
         pt[i]->SetTextAlign(12);
         pt[i]->SetTextFont(42);
         pt[i]->AddText( Form("#sigma(%i) = %s MeV",
                  date[i],PE.Eform(5+i*4,".2f",1e3)) );
         pt[i]->AddText( Form("N_{KK}(%i) = %s",
                  date[i],PE.Eform(6+i*4,".1f")) );
         pt[i]->AddText( Form("N_{sb}(%i) = %s",
                  date[i],PE.Eform(7+i*4,".1f")) );
         size_t ia = 8 + (i*4)*(FunSB==0);
         pt[i]->AddText( Form("a_{sb}(%i) = %s",
                  date[i],PE.Eform(ia,".1f")) );
         pt[i]->AddText( Form("N_{#phi}(%i) = %.1f",
                  date[i],Nphi[i].back()) );
      }

      // TPaveText* ptc = new TPaveText(0.11,0.31,0.35,0.49,"NDC");
      TPaveText* ptc = new TPaveText(0.61,0.71,0.89,0.89,"NDC");
      ptc->SetTextAlign(12);
      ptc->SetTextFont(42);
      ptc->AddText( Form("#vartheta= %s deg",
               PE.Eform(4,".0f",180./M_PI)) );
      ptc->AddText( Form("#it{L_{min}} = %.1f",lmin) );
      // ptc->AddText( Form("#lower[0.1]{M_{#phi}= %s MeV}",
      // PE.Eform(0,".2f",1e3)) );
      // ptc->AddText( Form("#Gamma_{#phi}= %s MeV",
      // PE.Eform(1,".3f",1e3)) );
      // ptc->AddText( Form("a = %s",PE.Eform(2,".2f")) );
      ptc->AddText( Form("F = %s",PE.Eform(3,".2f")) );

      c1->cd(1);
      leg->Draw();
      ptc->Draw();
      pt[0]->Draw();
      gPad->RedrawAxis();

      c1->cd(2);
      pt[1]->Draw();
      gPad->RedrawAxis();

      c1->cd(3);
      pt[2]->Draw();
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

   // * Nphi vs angl
   vector <TGraph*> grn(3,nullptr);
   for ( size_t i = 0; i < 3; ++i ) {
      c1->cd(i+1);
      grn[i] = new TGraph( nch, angl.data(), Nphi[i].data() );
      grn[i]->SetName( Form("Nphi_scan_%02i",date[i]%100) );
      grn[i]->SetTitle( Form("%i;#vartheta, degrees;N_{#phi}",
               date[i]) );
      grn[i]->SetMarkerColor(kBlue);
      grn[i]->SetMarkerStyle(21);
      grn[i]->SetLineWidth(2);
      grn[i]->GetYaxis()->SetTitleOffset(1.2);
      grn[i]->Draw("APL");
      grn[i]->Write();
   }
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   // * Lmin vs Nphi

   // subtract minimal value
   double minL = *(min_element(Lmin.begin(),Lmin.end()));
   for( auto &l : Lmin ) {
      l -= minL;
   }

   vector <TGraph*> grnl(3,nullptr);
   for ( size_t i = 0; i < 3; ++i ) {
      c1->cd(i+1);
      grnl[i] = new TGraph( nch, Nphi[i].data(), Lmin.data() );
      grnl[i]->SetName( Form("Nphi_Lmin_scan_%02i",date[i]%100) );
      grnl[i]->SetTitle( Form(";N_{#phi }(%i);#it{-2log(L/L_{max})}",
               date[i]) );
      grnl[i]->SetMarkerColor(kBlue);
      grnl[i]->SetMarkerStyle(20);
      grnl[i]->SetLineWidth(2);
      grnl[i]->GetYaxis()->SetTitleOffset(1.2);
      grnl[i]->Draw("APL");
      grnl[i]->Write();
   }
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   // * Lmin vs angle
   c1->cd();
   c1->Clear();
   c1->SetCanvasSize(800,800); // resize
   c1->cd();
   gPad->SetGrid();

   auto gr = new TGraph( nch, angl.data(), Lmin.data() );
   gr->SetName("Lmin_scan_cf");
   gr->SetTitle(";#vartheta, degrees;#it{-2log(L/L_{max})}");
   gr->GetYaxis()->SetMaxDigits(3);
   gr->GetYaxis()->SetTitleOffset(1.);
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(20);
   gr->SetLineWidth(2);
   gr->Draw("APL");
   gr->Write();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
   c_out.Close(); // close root-file with scan results
}

// {{{2 combineSBBR_Intfr(): combined plus side-band + fitBR
//--------------------------------------------------------------------
struct myFCN_combsbbr {
   // const double sl = -1.8; // prod-12
   const double sl = -2.0; // v709
   // Function for side-band:
   const int funSB = 1; // 0 - 'asb' for each year are independent
                        // 1 - one 'asb' for all years

   const vector<int> date { 2009, 2012, 2021 };
   const size_t ndate = 3;

   const double mphi = 1.01952; // Mphi like in MC
   const double gphi = Gphi;
   const double ar = 0.; // Argus parameter in interference

   vector<double> mkk[3];  // data central part
   vector<double> sb[3];   // data side band

   vector<double> Br2Nkk;  // conversion Br(J/Psi -> KK eta) -> Nkk
   vector<double> Br2Nphi; // conversion Br(J/Psi -> phi eta) -> Nphi

   // integrals
   vector<double> IeBW, IeAr, IeIc, IeIs;
   vector<double> IBW, IAr, IIc, IIs;

   vector<double> normArsb; // normalizations for Argus in side-band

   //-----------------------------------------------------------------
   myFCN_combsbbr() {
   //-----------------------------------------------------------------
      // Br(eta->2gamma) = 39.36 +/- 0.18% PDG 2023
      const double breta = 0.3936;
      // Br(phi->K+K-) = 49.1 +/- 0.5% PDG 2023
      const double brphi = 0.491;

      // prod-12
      // const double NJpsi09 = 15740130; // 2009
      // const double NJpsi12 = 48674393; // 2012
      // const double e_phi09 = 0.3131; // chi2<40
      // const double e_phi12 = 0.3075; // chi2<40

      // prod v709 first (2009,2012,2021)
      vector<double> NJpsi  { 15990331, 49882436, 324842271 };
      vector<double> e_phi  { 0.2841, 0.3112, 0.3163}; // eff chi2<40

      Br2Nkk.resize(ndate);
      Br2Nphi.resize(ndate);
      for ( size_t idat = 0; idat < ndate; ++idat ) {
         Br2Nkk[idat] = NJpsi[idat] * e_phi[idat] * breta;
         Br2Nphi[idat] = Br2Nkk[idat] * brphi;
      }

      IeBW.resize(ndate,0.);
      IeAr.resize(ndate,0.);
      IeIc.resize(ndate,0.);
      IeIs.resize(ndate,0.);
      IBW.resize(ndate,0.);
      IAr.resize(ndate,0.);
      IIc.resize(ndate,0.);
      IIs.resize(ndate,0.);

      normArsb.resize(ndate,1.);
   }

   // normalization for Argus SB
   //-----------------------------------------------------------------
   void calcNormAr(const vector<double>& asb) {
   //-----------------------------------------------------------------
      // cache for calculated values
      static vector<double> asb_save(3,-1); // ndate ?

      // integrand lambda function
      double par[] { 0., sl };
      auto Lintar = [](double x, void* pp) -> double{
         double* p = static_cast<double*>(pp);
         return RevArgus(x,p[0]) * (1+p[1]*(x-1.02));
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);

      // integrals
      for ( size_t idat = 0; idat < ndate; ++idat ) {
         if ( asb_save[idat] == asb[idat] ) {
            continue;
         }
         asb_save[idat] = asb[idat];

         if ( funSB == 1 && idat > 0 ) {
            normArsb[idat] = normArsb[0];
            continue;
         }

         par[0] = asb[idat];
         normArsb[idat] = gsl_int.Integral(Lintar,(void *)par,dL,dU);

         // debug print
         // printf(" %i: asb= %.3f, normArsb= %.3g\n",
               // date[idat], asb[idat], normArsb[idat]);
      }
   }

   //-----------------------------------------------------------------
   void calcIntegrals(const vector<double>& sig) {
   //-----------------------------------------------------------------
      // cache for calculated values
      static vector<double> sig_save(3,0.); // ndate ?

      // integrand lambda function: mphi,gphi,sig,ar,F,ang,sl,idx
      const double F = 1.;
      const double ang = 0.; //0(PI/2) for cos(sin) members of Infr.
      double pp[] { mphi,gphi,0.,ar,F,ang,sl, 1. };

      auto Lint = [](double x, void* pp) -> double{
         const double* p = static_cast<const double*>(pp);
         int idx = int(p[7]);
         return IntfrBWARG(x,p,idx);
      };

      // function for check results
      auto checkInt = [](string t, double Int, double pp[]) -> void {
         if ( !isfinite(Int) ) {
            printf("%s: pp= %g,%g,%g,%g,%g,%g,%g,%g\n", t.data(),
                  pp[0],pp[1],pp[2],pp[3],pp[4],pp[5],pp[6],pp[7]);
            exit(EXIT_FAILURE);
         }
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);

      // integrals
      for ( size_t idat = 0; idat < ndate; ++idat ) {
         if ( sig_save[idat] == sig[idat] ) {
            continue;
         }
         sig_save[idat] = sig[idat];
         pp[2] = sig[idat]; // sigma

         // calculate integrals
         pp[6] = sl; // for Ie..

         pp[7] = 1; // idx
         IeBW[idat] = gsl_int.Integral(Lint,(void *)pp,dL,dU);
         checkInt(Form("IeBW_%02i",date[idat]%100),IeBW[idat],pp);
         pp[7] = 2;
         IeAr[idat] = gsl_int.Integral(Lint,(void *)pp,dL,dU);
         checkInt(Form("IeAr_%02i",date[idat]%100),IeAr[idat],pp);
         pp[7] = 3;
         pp[5] = 0; // ang=0 for cos
         IeIc[idat] = gsl_int.Integral(Lint,(void *)pp,dL,dU);
         checkInt(Form("IeIc_%02i",date[idat]%100),IeIc[idat],pp);
         pp[5] = M_PI/2; // ang = pi/2 for sin
         IeIs[idat] = gsl_int.Integral(Lint,(void *)pp,dL,dU);
         checkInt(Form("IeIs_%02i",date[idat]%100),IeIs[idat],pp);

         pp[6] = 0; // sl=0 for I..

         pp[7] = 1; // idx
         IBW[idat] = gsl_int.Integral(Lint,(void *)pp,dL,dU);
         checkInt(Form("IBW_%02i",date[idat]%100),IBW[idat],pp);
         pp[7] = 2;
         IAr[idat] = gsl_int.Integral(Lint,(void *)pp,dL,dU);
         checkInt(Form("IAr_%02i",date[idat]%100),IAr[idat],pp);
         pp[7] = 3;
         pp[5] = 0; // ang=0 for cos
         IIc[idat] = gsl_int.Integral(Lint,(void *)pp,dL,dU);
         checkInt(Form("IIc_%02i",date[idat]%100),IIc[idat],pp);
         pp[5] = M_PI/2; // ang = pi/2 for sin
         IIs[idat] = gsl_int.Integral(Lint,(void *)pp,dL,dU);
         checkInt(Form("IIs_%02i",date[idat]%100),IIs[idat],pp);
      }

      // debug print
      // for ( size_t idat = 0; idat < ndate; ++idat ) {
         // printf(" %i: IeBW= %.6g IeAr= %.6g IeIc= %.6g IeIs= %.6g\n",
               // date[idat], IeBW[idat], IeAr[idar], IeIc[idat],
               // IeIs[idat]);
         // printf(" %i: IBW= %.6g IAr= %.6g IIc= %.6g IIs= %.6g\n",
               // date[idat], IBW[idat], IAr[idat], IIc[idat],
               // IIs[idat]);
      // }
   }

   // integrals MUST be calculated before calling this function
   //-----------------------------------------------------------------
   tuple<double,double> calcF(const double* p) const {
   //-----------------------------------------------------------------
   // NOTE: the ratio Nkk/Nphi = Brkk / (Brphi*Br(phi->KK)) does not
   // depend on the year, but the integrals do, though very slightly
   // so we use the integrals for the year with the largest statistics

      double Brkk = p[0];
      double Nkk = Brkk * Br2Nkk.back();

      double Brphi = p[1];
      double Nphi = Brphi * Br2Nphi.back();

      double ang = p[2];

      // calculate F:
      const size_t iy = 2; // 2021, 0 and 1 for systematic
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
         penalty = 1e5*fabs(Ff);
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
      double Brkk  = p[0];
      double Brphi = p[1];
      double ang   = p[2];

      // integrals:
      this->calcIntegrals({p[3],p[6],p[9]}); // sig09...sig21

      // normalization for Argus SB
      this->calcNormAr({p[5],p[8],p[11]}); // asb09...asb21

      // calculate F: must be after integrals
      double Ff = 0., penalty = 0.;
      tie(Ff,penalty) = this->calcF(p);

      double res = penalty;
      for ( size_t idat = 0; idat < ndate; ++idat ) {
         size_t idx = 3 + 3*idat;
         double sig = p[idx];
         double Nsb = p[idx+1];
         double asb = p[idx+2];
         if ( funSB == 1 ) {
            asb = p[3+2];
         }

         double Nkk  = Brkk  * Br2Nkk[idat];
         double Nphi = Brphi * Br2Nphi[idat];
         double NsbNorm = Nsb / normArsb[idat];

         double normI = IBW[idat]
            + Ff*( (IIc[idat]*cos(ang) + IIs[idat]*sin(ang))
                  + Ff*IAr[idat] );
         double NFit = Nkk/normI;

         double normE = IeBW[idat]
            + Ff*( (IeIc[idat]*cos(ang) + IeIs[idat]*sin(ang))
                  + Ff*IeAr[idat] );
         double NkkFit = NFit*normE;

         res += 2*(NkkFit+Nsb);

         const double pp[] { mphi,gphi,sig,ar,Ff,ang,sl };
         for ( const auto& m : mkk[idat] ) {
            double L = NFit * IntfrBWARG( m,pp,0 );
            L += NsbNorm * RevArgus(m,asb)*(1+sl*(m-1.02));
            if (L > 0.) {
               res -= 2*log(L);
            } else {
               res += FLT_MAX; // ~3e38
            }
         }

         // fit SB
         res += 2*Nsb;
         for ( const auto& m : sb[idat] ) {
            double L = NsbNorm * RevArgus(m,asb)*(1+sl*(m-1.02));
            if (L > 0.) {
               res -= 2*log(L);
            } else {
               res += FLT_MAX; // ~3e38
            }
         }
      }

      return res;
   }
};

//--------------------------------------------------------------------
void combineSBBR_Intfr(string pdf, int neg_pos=-1, bool plotsb=true) {
//--------------------------------------------------------------------
   myFCN_combsbbr my_fcn;          // class for 'FitFCN'
   const auto& date = my_fcn.date; // to shorten
   const size_t ndate = date.size();

   // Get un-binned data
   vector<TH1D*> hist(2*ndate,nullptr);
   size_t Ndat = 0;
   for ( size_t i = 0; i < ndate; ++i ) {
      string fname( Form("data_%02ipsip_all.root",date[i]%100) );
      my_fcn.mkk[i] =
         get_mkk_hist(fname,Form("cp_%i",date[i]),&hist[i]);
      Ndat += my_fcn.mkk[i].size();
      my_fcn.sb[i] =
         get_mkk_hist(fname,Form("sb_%i",date[i]),&hist[ndate+i], 1);
      Ndat += my_fcn.sb[i].size();
   }

   //-----------------------------------------------------------------
   // Fit data
   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   vector<string> par_name {
      "Brkk",  "Brphi", "angle",
      "sig09", "Nsb09", "asb09",
      "sig12", "Nsb12", "asb12",
      "sig21", "Nsb21", "asb21"
   };

   vector<double> par_ini {
      4.4e-4,  8.6e-4, 0.5,
      1.5e-3,      9.,  7.,
      1.1e-3,     34.,  4.,
      1.1e-3,    195.,  6.
   };

   // neg_pos is -/+ for destructive/constructive interference
   if ( neg_pos > 0 ) {
      par_ini[1] = 8.0e-4;
      par_ini[2] = -par_ini[2];
   }

   if ( my_fcn.funSB == 1 ) {
      par_ini[5] = par_ini[8] = par_ini[11]; // common asb
   }

   const size_t Npar = par_name.size(); // number of parameters

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(3e-4, 6e-4);   // Brkk
   fitter.Config().ParSettings(1).SetLimits(5e-4, 12e-4);  // Brphi
   fitter.Config().ParSettings(2).SetLimits(-M_PI, M_PI);  // angle

   for ( size_t i = 0; i < ndate; ++i ) { // sig,Nsb,asb for each year
      fitter.Config().ParSettings(3+3*i).SetLimits(0.5e-3,2.e-3);
      // fitter.Config().ParSettings(4+3*i).SetLimits(0.,2.*Nsb[i]);
      fitter.Config().ParSettings(5+3*i).SetLimits(0.,15.);
      if ( my_fcn.funSB == 1 && i > 0 ) {
         fitter.Config().ParSettings(5+3*i).Fix(); // asb
      }
   }

   // == Fit
   fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // likelihood
   // fitter.CalculateHessErrors();
   fitter.CalculateMinosErrors();

   const ROOT::Fit::FitResult res = fitter.Result();
   vector<double> Fpar = res.Parameters();
   // variables for KS-test and draw-function
   double Ff,penalty;
   tie(Ff,penalty) = my_fcn.calcF(Fpar.data());
   my_fcn.calcNormAr({Fpar[5],Fpar[8],Fpar[11]}); // asb09...asb21

   // print result
   if ( !isfinite(penalty) || penalty > 0. ) {
      printf("\n=> INVALID RESULT: penalty= %f <=\n",penalty);
   } else {
      printf("\n=> FINAL RESULT <=\n");
   }
   res.Print(cout);
   // res.PrintCovMatrix(cout); // print error matrix and correlations

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res,0.05); // ignore 5% upper/lower errors

   //-----------------------------------------------------------------
   // Print "middle points"
   // printf("++ Br(J/Psi->KKeta)= %7.3f *1e4\n",1e4*PE.Middle(0));
   // printf("++ Br(J/Psi->phi eta)= %6.2f *1e4\n",1e4*PE.Middle(1));

   //-----------------------------------------------------------------
   // "Goodness of fit" using KS-test
   string sepline(70,'='); // separation line
   printf("%s\n",sepline.c_str());
   vector<double> pvalueKS(ndate,0.), pvalueKSsb(ndate,0.);
   for ( size_t i = 0; i < ndate; ++i ) {
      // central part
      double Nkk = Fpar[0] * my_fcn.Br2Nkk[i];
      double Nphi = Fpar[1] * my_fcn.Br2Nphi[i];
      double NFit = Nphi/my_fcn.IBW[i]; // short formula?
      double ang = Fpar[2];
      double sig = Fpar[3+3*i];

      // {mphi,gphi,sigma,A,F,Ang,sl}
      const double pp[] {
         my_fcn.mphi, my_fcn.gphi, sig, my_fcn.ar, Ff, ang, my_fcn.sl
      };

      double Nsb = Fpar[4+3*i];
      double asb = Fpar[5+3*i];
      if ( my_fcn.funSB == 1 ) {
         asb = Fpar[5];
      }

      // normalization on 1
      double Ntot = Nkk + Nsb;
      NFit /= Ntot;

      double NsbNorm = Nsb / Ntot;
      NsbNorm /= my_fcn.normArsb[i];

      const auto& sl = my_fcn.sl;
      auto Fcr = [NFit,pp,NsbNorm,asb,sl](double x) -> double {
         double ret = NFit*IntfrBWARG( x,pp,0 );
         ret += NsbNorm * RevArgus(x,asb)*(1+sl*(x-1.02));
         return ret;
      };
      pvalueKS[i] = FastKSTest(Fcr,my_fcn.mkk[i]);

      // side-band
      NsbNorm = 1. / my_fcn.normArsb[i]; // norm to 1 for SB
      auto Fsb = [NsbNorm,asb,sl](double x) -> double {
         return NsbNorm * RevArgus(x,asb)*(1+sl*(x-1.02));
      };
      pvalueKSsb[i] = FastKSTest(Fsb,my_fcn.sb[i]);
      printf("%i: pvalueKS(cr)= %.3f pvalueKS(sb)= %.3f\n",
            date[i],pvalueKS[i],pvalueKSsb[i]);
   }
   printf("%s\n",sepline.c_str());

   //-----------------------------------------------------------------
   // Functions to draw
   std::function<double (double*, double*)> Ldr[3]; // ndate
   vector<TF1*> ff(ndate,nullptr);
   for ( size_t i = 0; i < ndate; ++i ) {
      Ldr[i] = [Fpar,my_fcn,Ff,i](double* x,double* p) {
         int isw = int(p[0]+0.5);
         double Nphi = Fpar[1] * my_fcn.Br2Nphi[i];
         double NFit = Nphi/my_fcn.IBW[i]; // short formula?
         double ang = Fpar[2];
         double sig = Fpar[3+3*i];

         // {mphi,gphi,sigma,A,F,Ang,sl}
         const double pp[] {
            my_fcn.mphi,my_fcn.gphi,sig,my_fcn.ar,Ff,ang,my_fcn.sl
         };
         double BWARG = NFit * IntfrBWARG( x[0],pp,isw );
         if ( isw > 0 && isw < 4 ) {
            return bW * BWARG;
         }

         double Nsb = Fpar[4+3*i];
         double asb = Fpar[5+3*i];
         if ( my_fcn.funSB == 1 ) {
            asb = Fpar[5];
         }
         double NsbNorm = Nsb / my_fcn.normArsb[i];
         double Bg = NsbNorm *
            RevArgus(x[0],asb)*(1+my_fcn.sl*(x[0]-1.02));
         if ( isw == 4 ) {
            return bW * Bg;
         }

         return bW * (BWARG + Bg);
      };
      ff[i] = new TF1(Form("ff%i",date[i]), Ldr[i], dL,dU,1);
      ff[i]->SetNpx(500);
   }

   // find max/min to draw
   for ( size_t i = 0; i < ndate; ++i ) {
      auto& hst = hist[i];
      int mbin = hst->GetMaximumBin();
      double hmax = hst->GetBinContent(mbin)+hst->GetBinError(mbin);

      auto& fun = ff[i];
      fun->SetParameter(0, 1); // BW
      double fmax = fun->GetMaximum( 1.01, 1.03, 1e-6);

      hmax = floor(1.1 * max(hmax,fmax));
      if ( hmax > 0 ) {
         hst->SetMaximum(hmax);
      }

      fun->SetParameter(0, 3); // interference
      double hmin = fun->GetMinimum( 1.01, 1.03, 1e-6);
      hmin = floor(1.15 * min(0.,hmin));
      if ( hmin < 0 ) {
         hst->SetMinimum(hmin);
      }

      if ( plotsb ) { // Side-band
         auto& hsb = hist[3+i];
         mbin = hsb->GetMaximumBin();
         hmax = hsb->GetBinContent(mbin)+hsb->GetBinError(mbin);
         if ( hmax < 10 ) {
            hsb->SetMaximum( hmax+2);
         } else {
            hsb->SetMaximum( 1.2*hmax );
         }
      }
   }

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000);
   c1->Divide(1,3);

   pdf += (my_fcn.funSB == 1) ? "_ar1" : "_ar";
   pdf += (neg_pos<0 ? "_n.pdf" : "_p.pdf");
   c1->Print((pdf+"[").c_str()); // open pdf-file

   TLegend* leg = new TLegend(0.11,0.50,0.35,0.89);
   vector<TPaveText*> pt(ndate,nullptr);

   for ( size_t i = 0; i < ndate; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();

      auto& hst = hist[i];
      auto& fun = ff[i];

      SetHstFace(hst); // Data
      hst->GetXaxis()->SetTitleOffset(1.1);
      hst->GetYaxis()->SetTitleOffset(1.0);
      hst->SetLineWidth(2);
      hst->SetLineColor(kBlack);
      hst->SetMarkerStyle(20);
      hst->SetMarkerSize(0.8);
      hst->Draw("EP");
      if ( i == 0 ) {
         leg->AddEntry( hist[0],"Data","LEP" );
      }

      fun->SetParameter(0, 0); // Sum
      fun->SetLineWidth(2);
      fun->SetLineStyle(kSolid);
      fun->SetLineColor(kRed+1);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         leg->AddEntry( fun->Clone(), "Combined fit", "L" );
      }

      fun->SetParameter(0, 1); // BW
      fun->SetLineWidth(2);
      fun->SetLineStyle(kDashed);
      fun->SetLineColor(kGreen+2);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         leg->AddEntry( fun->Clone(), "Breit-Wigner #phi#eta", "L");
      }

      fun->SetParameter(0, 2); // Argus
      fun->SetLineWidth(1);
      fun->SetLineStyle(kDashed);
      fun->SetLineColor(kBlue);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         leg->AddEntry( fun->Clone(), "Non-#phi KK#eta", "L");
      }

      fun->SetParameter(0, 3); // interference
      fun->SetLineWidth(1);
      fun->SetLineStyle(kDashed);
      fun->SetLineColor(kMagenta+1);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         leg->AddEntry( fun->Clone(), "Interference", "L");
      }

      double yptm = 0.59 + 0.06*(my_fcn.funSB==1);
      pt[i] = new TPaveText(0.61,yptm,0.89,0.89,"NDC");
      pt[i]->SetTextAlign(12);
      pt[i]->SetTextFont(42);
      pt[i]->AddText( Form("#it{p-value(%i) = %.3f}",
               date[i],pvalueKS[i]) );
      pt[i]->AddText( Form("#sigma(%i) = %s MeV",
               date[i],PE.Eform(3+i*3,".2f",1e3)) );
      pt[i]->AddText( Form("#it{p-sideband(%i) = %.3f}",
               date[i],pvalueKSsb[i]) );
      pt[i]->AddText( Form("N_{sb}(%i) = %s",
               date[i],PE.Eform(4+i*3,".1f")) );
      if ( my_fcn.funSB == 0 ) {
         pt[i]->AddText( Form("a_{sb}(%i) = %s",
                  date[i],PE.Eform(5+i*3,".1f")) );
      }
   }

   double yptm = 0.59 + 0.06*(my_fcn.funSB==1);
   TPaveText* ptc = new TPaveText(0.61,0.29,0.89,yptm,"NDC");
   ptc->SetTextAlign(12);
   ptc->SetTextFont(42);
   ptc->AddText( "#it{common parameters}");
   ptc->AddText( Form("Br(KK#eta)= %s #times10^{-4 }",
            PE.Eform(0,".2f",1e4)) );
   ptc->AddText( Form("Br(#phi#eta)= %s #times10^{-4 }",
            PE.Eform(1,".2f",1e4)) );
   ptc->AddText( Form("#vartheta = %s",PE.Eform(2,".2f")) );
   if ( my_fcn.funSB == 1 ) {
      ptc->AddText( Form("a_{sb} = %s", PE.Eform(5,".1f")) );
   }

   c1->cd(1);
   leg->Draw();
   pt[0]->Draw();
   ptc->Draw();
   gPad->RedrawAxis();

   c1->cd(2);
   pt[1]->Draw();
   gPad->RedrawAxis();

   c1->cd(3);
   pt[2]->Draw();
   gPad->RedrawAxis();

   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   if ( !plotsb ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
      return;
   }

   //-----------------------------------------------------------------
   // Side-band
   TLegend* legsb = new TLegend(0.11,0.72,0.40,0.89);
   vector<TPaveText*> ptsb(ndate,nullptr);

   for ( size_t i = 0; i < 3; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();

      auto& hst = hist[3+i];
      auto& fun = ff[i];

      SetHstFace(hst);
      hst->GetXaxis()->SetTitleOffset(1.1);
      hst->GetYaxis()->SetTitleOffset(1.3);
      hst->SetLineWidth(2);
      hst->SetLineColor(kBlack);
      hst->SetMarkerStyle(20);
      hst->Draw("EP");
      if ( i == 0 ) {
         legsb->AddEntry(hst,"Side-band data","LEP");
      }

      fun->SetParameter(0, 4); // SB
      fun->SetLineWidth(2);
      fun->SetLineStyle(kSolid);
      fun->SetLineColor(kRed+1);
      fun->DrawCopy("SAME");
      if ( i == 0 ) {
         legsb->AddEntry(fun->Clone(),"Fit by Argus function:","L");
      }

      ptsb[i] = new TPaveText(0.70,0.61,0.89,0.89,"NDC");
      ptsb[i]->SetTextAlign(12);
      ptsb[i]->SetTextFont(42);
      ptsb[i]->AddText( Form("#bf{%16i}", date[i]) );
      ptsb[i]->AddText( Form("p-value = %.3f",pvalueKSsb[i]) );
      ptsb[i]->AddText( Form("N_{sb} = %s", PE.Eform(4+i*3,".1f")) );
      size_t ia = 5 + (i*3)*(my_fcn.funSB==0);
      ptsb[i]->AddText( Form("#lower[-0.01]{a_{sb} = %s}",
               PE.Eform(ia,".1f")) );
   }

   c1->cd(1);
   legsb->Draw();
   ptsb[0]->Draw();
   gPad->RedrawAxis();

   c1->cd(2);
   ptsb[1]->Draw();
   gPad->RedrawAxis();

   c1->cd(3);
   ptsb[2]->Draw();
   gPad->RedrawAxis();

   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
}

// {{{2  combineSBBR_Intfr_scan()
//--------------------------------------------------------------------
void combineSBBR_Intfr_scan(string sname) {
//--------------------------------------------------------------------
   myFCN_combsbbr my_fcn;          // class for 'FitFCN'
   const auto& date = my_fcn.date; // to shorten
   const size_t ndate = date.size();

   // Get un-binned data
   vector<TH1D*> hist(2*ndate,nullptr);
   size_t Ndat = 0;
   for ( size_t i = 0; i < ndate; ++i ) {
      string fname( Form("data_%02ipsip_all.root",date[i]%100) );
      my_fcn.mkk[i] =
         get_mkk_hist(fname,Form("cp_%i",date[i]),&hist[i]);
      Ndat += my_fcn.mkk[i].size();
      my_fcn.sb[i] =
         get_mkk_hist(fname,Form("sb_%i",date[i]),&hist[ndate+i], 1);
      Ndat += my_fcn.sb[i].size();
   }

   //-----------------------------------------------------------------
   // Fit data
   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   vector<string> par_name {
      "Brkk",  "Brphi", "angle",
      "sig09", "Nsb09", "asb09",
      "sig12", "Nsb12", "asb12",
      "sig21", "Nsb21", "asb21"
   };

   vector<double> par_ini { // for funSB=0
      4.3e-4,  7.5e-4, -0.8,
      1.4e-3,     10.,   7.,
      1.1e-3,     35.,   4.,
      1.1e-3,    195.,   6.
   };

   if ( my_fcn.funSB == 1 ) {
      par_ini[5] = par_ini[8] = par_ini[11]; // common asb
   }

   // parameters of the loop for Brphi
   // double bp_min = 7.5e-4, bp_max = 9.21e-4, bp_step = 2.e-4;
   double bp_min = 7.5e-4, bp_max = 9.21e-4, bp_step = 0.02e-4;

   const size_t Npar = par_name.size(); // number of parameters

   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetLimits(3e-4, 6e-4);   // Brkk
   fitter.Config().ParSettings(2).SetLimits(-M_PI, M_PI);  // angle

   for ( size_t i = 0; i < ndate; ++i ) { // sig,Nsb,asb for each year
      fitter.Config().ParSettings(3+3*i).SetLimits(0.2e-3,5.e-3);
      // fitter.Config().ParSettings(4+3*i).SetLimits(0.,2.*Nsb[i]);
      fitter.Config().ParSettings(5+3*i).SetLimits(0.,15.);
      if ( my_fcn.funSB == 1 && i > 0 ) {
         fitter.Config().ParSettings(5+3*i).Fix(); // asb
      }
   }

   //-----------------------------------------------------------------
   // Functions to draw
   vector<double> Fpar;
   double Ff = 0;
   std::function<double (double*, double*)> Ldr[3]; // ndate?
   vector<TF1*> ff(ndate,nullptr);
   for ( size_t i = 0; i < ndate; ++i ) {
      Ldr[i] = [&Fpar,&my_fcn,&Ff,i](double* x,double* p) {
         int isw = int(p[0]+0.5);
         double Nphi = Fpar[1] * my_fcn.Br2Nphi[i];
         double NFit = Nphi/my_fcn.IBW[i]; // short formula?
         double ang = Fpar[2];
         double sig = Fpar[3+3*i];

         // {mphi,gphi,sigma,A,F,Ang,sl}
         const double pp[] {
            my_fcn.mphi,my_fcn.gphi,sig,my_fcn.ar,Ff,ang,my_fcn.sl
         };
         double BWARG = NFit * IntfrBWARG( x[0],pp,isw );
         if ( isw > 0 && isw < 4 ) {
            return bW * BWARG;
         }

         double Nsb = Fpar[4+3*i];
         double asb = Fpar[5+3*i];
         if ( my_fcn.funSB == 1 ) {
            asb = Fpar[5];
         }
         double NsbNorm = Nsb / my_fcn.normArsb[i];
         double Bg = NsbNorm *
            RevArgus(x[0],asb)*(1+my_fcn.sl*(x[0]-1.02));
         if ( isw == 4 ) {
            return bW * Bg;
         }

         return bW * (BWARG + Bg);
      };
      ff[i] = new TF1(Form("ff%i",date[i]), Ldr[i], dL,dU,1);
      ff[i]->SetNpx(500);
   }

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000);
   c1->Divide(1,3);
   c1->cd();

   sname += (my_fcn.funSB == 1) ? "_ar1" : "_ar";
   string hfile = sname + ".root";
   TFile c_out(hfile.c_str(),"recreate"); // file for scan results
   c_out.cd();

   string pdf = sname + ".pdf";
   c1->Print((pdf+"[").c_str()); // open pdf-file

   vector<double> hmin { -15., -35., -200.};
   for ( size_t i = 0; i < ndate; ++i ) {
      auto& hst = hist[i];

      SetHstFace(hst); // Data
      hst->GetXaxis()->SetTitleOffset(1.1);
      hst->GetYaxis()->SetTitleOffset(1.0);
      hst->SetLineWidth(2);
      hst->SetLineColor(kBlack);
      hst->SetMarkerStyle(20);
      hst->SetMarkerSize(0.8);
      hst->SetMinimum(hmin[i]);
   }

   //-----------------------------------------------------------------
   // start the cycle
   vector<double> br_phi, Lmin;
   for( double bp = bp_min; bp < bp_max; bp += bp_step ) {

      fitter.Config().ParSettings(1).SetValue(bp); // Brphi
      fitter.Config().ParSettings(1).Fix();

      // == Fit
      fitter.FitFCN(Npar,my_fcn,nullptr,Ndat,false); // likelihood
      // fitter.CalculateHessErrors();
      // fitter.CalculateMinosErrors();

      ROOT::Fit::FitResult res = fitter.Result();
      Fpar = res.Parameters();
      if ( bp_step > 0.5 ) {
         res.Print(cout);
      }
      ParAndErr PE(res);

      double lmin = res.MinFcnValue();
      br_phi.push_back(bp);
      Lmin.push_back(lmin);

      string sepline(50,'='); // separation line
      printf("%s\n",sepline.c_str());
      printf("bp= %e =>  Lmin= %.2f\n", bp,lmin);
      printf("%s\n",sepline.c_str());

      //--------------------------------------------------------------
      // Draw intermidiate results
      double penalty = 0;
      tie(Ff,penalty) = my_fcn.calcF(Fpar.data());

      TLegend* leg = new TLegend(0.11,0.50,0.35,0.89);
      vector<TPaveText*> pt(ndate,nullptr);

      for ( size_t i = 0; i < ndate; ++i ) {
         c1->cd(i+1);
         gPad->SetGrid();

         auto& hst = hist[i];
         auto& fun = ff[i];

         hst->Draw("EP");
         if ( i == 0 ) {
            leg->AddEntry( hist[0],"Data","LEP" );
         }

         fun->SetParameter(0, 0); // Sum
         fun->SetLineWidth(2);
         fun->SetLineStyle(kSolid);
         fun->SetLineColor(kRed+1);
         fun->DrawCopy("SAME");
         if ( i == 0 ) {
            leg->AddEntry( fun->Clone(), "Combined fit", "L" );
         }

         fun->SetParameter(0, 1); // BW
         fun->SetLineWidth(2);
         fun->SetLineStyle(kDashed);
         fun->SetLineColor(kGreen+2);
         fun->DrawCopy("SAME");
         if ( i == 0 ) {
            leg->AddEntry(fun->Clone(), "Breit-Wigner #phi#eta", "L");
         }

         fun->SetParameter(0, 2); // Argus
         fun->SetLineWidth(1);
         fun->SetLineStyle(kDashed);
         fun->SetLineColor(kBlue);
         fun->DrawCopy("SAME");
         if ( i == 0 ) {
            leg->AddEntry( fun->Clone(), "Non-#phi KK#eta", "L");
         }

         fun->SetParameter(0, 3); // interference
         fun->SetLineWidth(1);
         fun->SetLineStyle(kDashed);
         fun->SetLineColor(kMagenta+1);
         fun->DrawCopy("SAME");
         if ( i == 0 ) {
            leg->AddEntry( fun->Clone(), "Interference", "L");
         }

         pt[i] = new TPaveText(0.61,0.47,0.89,0.65,"NDC");
         pt[i]->SetTextAlign(12);
         pt[i]->SetTextFont(42);
         pt[i]->AddText( Form("#sigma(%i) = %s MeV",
                  date[i],PE.Eform(3+i*3,".2f",1e3)) );
         pt[i]->AddText( Form("N_{sb}(%i) = %s",
                  date[i],PE.Eform(4+i*3,".1f")) );
         size_t ia = 5 + (i*3)*(my_fcn.funSB==0);
         pt[i]->AddText( Form("a_{sb}(%i) = %s",
                  date[i],PE.Eform(ia,".1f")) );
      }

      TPaveText* ptc = new TPaveText(0.61,0.65,0.89,0.89,"NDC");
      ptc->SetTextAlign(12);
      ptc->SetTextFont(42);
      ptc->AddText( Form("Br(#phi#eta)= %s #times10^{-4 }",
               PE.Eform(1,".2f",1e4)) );
      ptc->AddText( Form("#it{L_{min}} = %.1f",lmin) );
      ptc->AddText( Form("Br(KK#eta)= %s #times10^{-4 }",
               PE.Eform(0,".2f",1e4)) );
      ptc->AddText( Form("#vartheta = %s",PE.Eform(2,".2f")) );

      c1->cd(1);
      leg->Draw();
      ptc->Draw();
      pt[0]->Draw();
      gPad->RedrawAxis();

      c1->cd(2);
      pt[1]->Draw();
      gPad->RedrawAxis();

      c1->cd(3);
      pt[2]->Draw();
      gPad->RedrawAxis();

      c1->Update();
      c1->Print(pdf.c_str()); // add to pdf-file
   }

   //-----------------------------------------------------------------
   // final draw
   int nch = br_phi.size();
   if ( nch < 3 ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
      return;
   }

   // subtract minimal value
   double minL = *(min_element(Lmin.begin(),Lmin.end()));
   for( auto &l : Lmin ) {
      l -= minL;
   }

   // * Br vs Lmin
   c1->cd();
   c1->Clear();
   c1->SetCanvasSize(800,800); // resize
   c1->cd();
   gPad->SetGrid();

   auto gr = new TGraph( nch, br_phi.data(), Lmin.data() );
   gr->SetTitle(";Br(#phi#eta);#it{-2log(L/L_{max})}");
   gr->GetYaxis()->SetMaxDigits(3);
   gr->GetYaxis()->SetTitleOffset(1.);
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(20);
   gr->SetLineWidth(2);
   gr->Draw("APL");
   gr->Write();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
   c_out.Close(); // close root-file with scan results
}

// {{{1 MAIN:
//--------------------------------------------------------------------
void mass_kk_fit(int date=2009) {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(42);
   gStyle->SetLegendFont(42);

   // check date
   if ( date != 2009 && date != 2012 && date != 2021 ) {
      cerr << "ERROR: unknown date= " << date << endl;
      exit(EXIT_FAILURE);
   }

   if ( DEBUG ) {
      // BOOK HISTOGRAMMS
      htt[0] = new TH1D("eBWG","log10(BWG conv. error)",100,-6.,-1.);
      htt[1] = new TH1D("nBW","BreitWigner norm",1000,0.,2.);
      htt[2] = new TH1D("eBW","log10(BW norm. error)",100,-6.,-1.);
      htt[3] = new TH1D("nAr","RevArg norm",1000,0.,10.);
      htt[4] = new TH1D("eAr","log10(RevArg norm.err.)",100,-6.,-1.);
      htt[5] = new TH1D("eInG","log10(Intfr.conv.err.)",100,-6.,-1.);
   }

//--------------------------------------------------------------------

   // = define GSL error handler which does nothing =
   gsl_set_error_handler_off();

   // = set handler =
   // gsl_set_error_handler( my_handler );

   // set integrator: ROOT::Math::GSLIntegrator adaptive method (QAG)
   ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Adaptive");
   ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-9);

// ------------- signal sec.6.1 --------------------------------------
   // test_BreitWigner();

   // string pdf_mcsig( Form("mcmkk%02i_sig.pdf",date%100) );
   // sig_fit_mc(date, pdf_mcsig);

   // fit MC phi-eta by BW x Gauss
   // string pdf_mc( Form("mkk%02i_sig_log.pdf",date%100) );
   // sig_fit(date,pdf_mc);

// ------------- background sec.6.2 ----------------------------------
   // test_Argus();
   // test_RevArgus(); // with normalization

   // bkg_fit(date, 0); // fit kk-eta by Argus
   // bkg_fit(date, 1); // fit data side-band by Argus

   // TODO: fit rho-eta by Argus
   // bkg_fit("mcrhoeta_kkmc_12.root","MC #rho#eta 2012",
           // "mkk_rhoeta12_argus.pdf");

// ------------- data: no interference sec.6.3.1 ---------------------
   int FixM = 1; // 0/1 - mphi is float/fixed
   // string pdf_dat( Form("mkk%02i_fit%i.pdf",date%100,FixM) );
   // data_fit(date,FixM,pdf_dat);

   // ++ combined with Side-Band (constant or Argus SB)
   // string pdf_dat( Form("mkk%02i_fitSB_%i.pdf",date%100,FixM) );
   // dataSB_fit(date,FixM,pdf_dat);

   // TODO: I/O check: use mcinc files as data
   // string pdf = string("mkk") + ((id==0) ? "09" : "12") +
      // "_chkIO_" + to_string(Mphix) + ".pdf";
   // dataSB_fit(fnames.at(id+3),titles.at(id+3),Mphix,pdf);

// ------------- data: interference sec 6.3.2 ------------------------
   // test_Intfr();

   // data_Intfr_fit(date, Form("mkk%02i_ifr",date%100) );

   // ++ combined with Side-Band
   // bool plotSB = false; // last argument in _SB funcs
   // dataSB_Intfr( date, Form("mkk%02i_ifrSB",date%100) );

   // string pdf_sc( Form("mkk%02i_ifrSB_scan",date%100) );
   // dataSB_Intfr_scan(date,pdf_sc);

   // ++ combined with Side-Band + fitBR
   // dataSBBR_Intfr( date, Form("mkk%02i_ifrSBBR",date%100) );

// ------------- data: common fit ------------------------------------
   // combine_Intfr("mkk_cf_ifr", -1); // +/- = neg/pos interference

   // combined with Side-Band
   // combineSB_Intfr("mkk_cfSB",-1); // +/- = neg/pos interfer.

   // combineSB_Intfr_scan("mkk_cfSB_scan");

   // combined with Side-Band + fitBR
   // combineSBBR_Intfr("mkk_cfSBBR",-1);

   combineSBBR_Intfr_scan("mkk_cfSBBR_scan");

//--------------------------------------------------------------------
}
