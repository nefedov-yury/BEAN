// fit distributions M(K+K-) for data
// NOTE: see also PsipJpsiPhiEta/mass_kk_fit.cc and mkk_fitch2.cc
//       memo_brf_v15a: phase_angle= 0.05^+0.52_-0.63 = [-0.58,+0.58]
// -> outputs:
//    'mkk_inter/' for interference model BW with Argus
//    'mkk_noint/' for model with sum BW + Argus


#include "gsl/gsl_errno.h" // GSL error handler

#include "masses.h"

// {{{1 constants
//--------------------------------------------------------------------
// the meson radial parameter in Blatt-Weisskopf form factor
static const double R_BW = 3.; // GeV^{-1}; +/-1 uncertainties study

static const double dL = 2*Mk; // the left cutoff = 0.987354

static const double MphiMC = 1.01953; // Mphi like in signal MC

// binning for Breit-Wigner: bin width 1.0 MeV
static const int Nbins = 50;
static const double bL = 0.98; // first bin < dL
static const double dU = 1.08; // MUST BE < 1.0835 !!!
static const double bW = (dU-bL)/Nbins; // bin width

// {{{1 helper function
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

   // pretty print for table
   //-----------------------------------------------------------------
   const char* Tform (int i,string fm, double X=1, bool asym=false) {
   //-----------------------------------------------------------------
      if ( Perr[i] > 0 ) {
         string fmt = "%" +fm+ " \\pm %" +fm;
         return Form(fmt.c_str(),Fpar[i]*X,Perr[i]*X);
      } else if ( Perr[i] < 0 ) {
         if ( asym ) { // print asymmetric error
            string fmt = "%" + fm
               + "^{%+" + fm + "}" + "_{%+" + fm + "}";
            return Form(fmt.c_str(),Fpar[i]*X,Uerr[i]*X,-Lerr[i]*X);
         } else { // print the symmetric maximum error
            double perr = max(Uerr[i],Lerr[i]);
            string fmt = "%" +fm+ " \\pm %" +fm;
            return Form(fmt.c_str(),Fpar[i]*X,perr*X);
         }
      }
      string fmt = "%" +fm+ " (fixed)";
      return Form(fmt.c_str(),Fpar[i]*X);
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

// {{{1 data processing (see mass_KK.cc)
//--------------------------------------------------------------------
vector<double> get_mkk_hist( string fname, string hname,
                              TH1D* hst[], int type = 0 ) {
//--------------------------------------------------------------------
#include "cuts.h"

   // name of folder with root files
   static string dir("Ntpls/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("SelectKKgg");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   TCut c_here = c_chi2+c_phi;
   if ( type == 0 ) {        // central part
      c_here += c_cpgg;
   } else if ( type == 1 ) { // side-band
      c_here += c_sbgg;
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
tuple<TH1D*,int,double,double> Subtract(string fname) {
//--------------------------------------------------------------------
// Side-band subtraction: return resulting histogram, number of bins
// with nonzero error and total number of events (integral)
   TH1D* hst1 = plot_Mkk(fname,string("mkk_data_cp"),0);
   TH1D* hst2 = plot_Mkk(fname,string("mkk_data_sb"),1);
   TH1D* sub = (TH1D*)hst1->Clone("sub");
   sub->Add(hst1,hst2,1.,-1.);

   size_t nbx = sub->GetNbinsX();

   // set negative bins to zero for fitting by likelihood
   for ( int ib = 1; ib <= nbx; ++ib ) {
      double bc = sub->GetBinContent(ib);
      if ( bc < -1e-3 ) {
         double bw = sub->GetBinWidth(ib);
         double mkk = sub->GetBinLowEdge(ib);
         sub->SetBinContent(ib,0.);
         if ( mkk+bw < dL ) { // ignore this bin completely
            sub->SetBinError(ib,0.);
            continue;
         }
         // subtract 'bc' from nearby (+/-10) bins
         for ( int j = 0; j < 20; j++ ) {
            int jbin = ib + (1+j/2)*(1-2*(j%2)); // 1,-1,2,-2...
               double bcp = sub->GetBinContent(jbin);
               while ( bcp > 1e-3 && bc < -1e-3 ) {
                  bcp-=1;
                  bc+=1;
               }
               sub->SetBinContent(jbin,bcp);
               if ( bc > -1e-3 ) {
                  break;
               }
         }
      }
   }

   // calculate numbers: nonzero bins, events after subtraction
   int Nbin = 0;
   double Nev = 0;
   for ( int ib = 1; ib <= nbx; ++ib ) {
      double bc = sub->GetBinContent(ib);
      Nev += bc;
      double be = sub->GetBinError(ib);
      if ( fabs(be) > 1e-3 ) {
         Nbin += 1;
      }
   }
   // number side-band events
   double Nsb = hst2->Integral();
   return make_tuple(sub,Nbin,Nev,Nsb);
}

//--------------------------------------------------------------------
TH1D* get_mcMkk(string fname) {
//--------------------------------------------------------------------
   // name of folder with root files
   static string dir("Ntpls/");

   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( !froot ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("SelectKKgg");
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
// Breigt Wigner for e+ e- -> phi eta -> K+K- eta
// Eee is the energy of (e+e-) != Mjpsi
//--------------------------------------------------------------------

//--------------------------------------------------------------------
double BreitWigner(double m, double Eee, double mphi, double gphi,
      double R = R_BW) {
//--------------------------------------------------------------------
   if ( m < dL ) { // phase space becomes zero (see p_k)
      return 0;
   }

   // static const double R   = R_BW; // Blatt-Weisskopf ff
   constexpr double Meta2  = SQ(Meta);
   constexpr double Mk2    = SQ(Mk);

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
   double BB_phi = sqrt( (1+SQ(R*p_phi0)) / (1+SQ(R*p_phi)) );
   double BB_k2  = (1+SQ(R*p_k0)) / (1+SQ(R*p_k));
   double BB_k = sqrt( BB_k2 );

   double fm = r_phi*r_k*BB_phi*BB_k;
   double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k2; // == m*G(m)

   return (kap * SQ(fm)) / ( SQ(m2-mphi2) + SQ(GM) );
}

//--------------------------------------------------------------------
double BreitWignerGauss( double m, double Eee,
                         double mphi, double gphi, // BW parameters
                         double sigma              // Gauss
                       ) {
//--------------------------------------------------------------------
// Function Fun(x) folded with a normal distribution:
//
//        1                        [                  (X'-X)^2    ]
// ---------------- *   Integral   [ Fun(X') * exp( - --------- ) ] dX'
// sqrt(2*pi)*sigma   (-oo<X'<+oo) [                  2*sigma^2   ]
//
//--------------------------------------------------------------------
// sigma ~ 1/3 gphi (bes3-exp)
// integration in +/-5 sigma interval
//--------------------------------------------------------------------
   static const double Ngauss = 1./sqrt(2*M_PI);

   // integrand lambda function
   double pp[] = { m, Eee, mphi, gphi, sigma };
   auto Lint = [](double t, void* pp) -> double{ // t == (X'-X)/sigma
      double* p = static_cast<double*>(pp);
      double x = p[0] + t*p[4];
      return exp(-0.5*t*t) * BreitWigner( x,p[1],p[2],p[3]);
   };

   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-7, 1.e-6, 1000);
   double result = gsl_int.Integral(Lint,pp,-5.,+5.);

   return Ngauss * result;
}

//--------------------------------------------------------------------
double BreitWignerGaussN( double m, double Eee,
                          double mphi, double gphi, double sigma
                        ) {
//--------------------------------------------------------------------
// Numeric normalisation on one for range [dL,dU]

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

//--------------------------------------------------------------------
void test_BreitWigner() {
//--------------------------------------------------------------------

   // 1) BreitWigner (internally normalized to 1)
   // double Eee = Mjpsi;
   double Eee = 3.08;
   auto Lbw = [Eee](const double* x,const double* p) -> double {
      return p[0]*BreitWigner(x[0],Eee,p[1],p[2]);
   };
   TF1* bw = new TF1("bw", Lbw, dL, dU, 3);
   bw->SetParNames("Norm","Mphi","Gphi");
   bw->SetParameters(1., Mphi, Gphi);
   bw->SetLineWidth(1);
   bw->SetLineColor(kBlue);
   bw->SetNpx(500);

   double norm = 1./bw->Integral(dL,dU,1e-8);
   printf("norm = %.7f\n",norm);
   bw->SetParameter(0, norm );

   // 2) BreitWignerGaussN
   auto Lbwgn = [Eee](const double* x,const double* p) -> double {
      return p[0]*BreitWignerGaussN(x[0],Eee,p[1],p[2],p[3]);
   };
   TF1* bwgn = new TF1("bwgn", Lbwgn, dL, dU, 4);
   bwgn->SetParNames("Norm","Mphi","Gphi","Sigma");
   bwgn->SetParameters(1., Mphi, Gphi, 1.2e-3);
   bwgn->SetLineWidth(2);
   bwgn->SetLineColor(kRed);
   bwgn->SetLineStyle(kDashed);
   bwgn->SetNpx(500);

   double normN = 1./bwgn->Integral(dL,dU,1e-8);
   printf("normN = %.6f\n",normN);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader("Breit-Wigner (M_{#phi}, #Gamma_{#phi})","C");
   leg->AddEntry(bw->Clone(), "BW, Rff=3GeV^{-1}", "L");
   leg->AddEntry( bwgn->Clone(), Form("BWG(#sigma=%.2f MeV)",
            bwgn->GetParameter(3)*1e3), "L" );

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(true);

   bw->Draw();
   bwgn->Draw("SAME");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
}

// {{{1 Fit signal MC by Breit-Wigner
//--------------------------------------------------------------------
void sig_fit_mc(string Ename, double Eee, string pdf) {
//--------------------------------------------------------------------
   string fname( Form("ntpl_mcgpj_%s.root",Ename.c_str()) );
   TH1D* hst = get_mcMkk(fname);
   double bin_width = hst->GetBinWidth(1);

   // Fit signal MC(phi,eta) by Breit-Wigner
   auto Lbwg = [Eee,bin_width](const double* x,const double* p) {
      return bin_width*p[0]*BreitWigner(x[0],Eee,p[1],p[2],p[3]);
   };
   TF1* bwg = new TF1("bwg", Lbwg, 0.98, 1.2, 4);
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

   TCanvas* c1 = new TCanvas("c1","...",0,0,1000,800);
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
   double maxMkk=1.083;
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
   double norm108 = bwg->GetParameter(0) / intg108;
   cout << " norm108= " << norm108 << endl;

   TLegend* leg = new TLegend(0.22,0.20,0.62,0.40);
   leg->SetTextSize(0.03);
   leg->SetHeader( Form("Cut-off correction: %.3f#pm%.3f",corr,dcorr),
         "C" );
   string title(Ename);
   if ( Ename == "3080_rs" ) {
      title = "3080MeV";
   } else if ( Ename == "J4" ) {
      title = "3097MeV";
   }
   leg->AddEntry(hst, Form("MC truth %s",title.c_str()),"LEP");
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
void sig_fit(string Ename, double Eee, string pdf) {
//--------------------------------------------------------------------
   // Use binned data
   string fname( Form("ntpl_mcgpj_%s.root",Ename.c_str()) );
   TH1D* hst = nullptr;
   int Nbin = 0;
   double Nev = 0, Nsb = 0;
   tie(hst,Nbin,Nev,Nsb) = Subtract(fname); // subtract side-band

   ROOT::Fit::DataOptions opt;
   // use integral of bin content instead of bin center (def: false)
   opt.fIntegral = true;
   // use empty bins (default is false) with a fixed error of 1
   opt.fUseEmpty = true; // must be even for LH

   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::BinData Sig(opt,dr);
   ROOT::Fit::FillData(Sig, hst);

   // Fit signal MC by Breit-Wigner convoluted with Gauss
   auto Lbwg = [Eee](const double* x,const double* p) -> double {
      return bW * p[3]*BreitWignerGaussN(x[0],Eee,p[0],p[1],p[2]);
   };

   vector<string> par_name { "Mphi", "Gphi", "sigma", "Nphi" };
   vector<double> par_ini  {  Mphi,   Gphi,  1.e-3, Nev-Nsb };
   const size_t Npar = par_name.size(); // number of parameters
   TF1* bwg = new TF1("bwg", Lbwg, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 fbwg( *bwg );
   fitter.SetFunction(fbwg,false); // false=no parameter derivatives

   // Set parameters, data must come before names
   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   //-----------------------------------------------------------------
   fitter.Config().ParSettings(0).SetLimits(Mphi-0.01,Mphi+0.01);
   // fitter.Config().ParSettings(0).Fix();
   fitter.Config().ParSettings(1).Fix();
   // fitter.Config().ParSettings(1).SetLimits(Gphi-1e-4,Gphi+1e-4);
   fitter.Config().ParSettings(2).SetLimits(0.1e-3, 3.e-3); // sigma
   // fitter.Config().ParSettings(4).SetLimits(0.8*Nev,1.2*Nev);// Nphi

   // Fit
   //-----------------------------------------------------------------
   fitter.LikelihoodFit( Sig );
   fitter.CalculateMinosErrors();
   // fitter.CalculateHessErrors(); // without MINOS here

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res,0.1);
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------
   // "Goodness of fit" evaluated by Chi2
   double chi2 = res.Chi2();
   unsigned int Ndf = res.Ndf();

   //-----------------------------------------------------------------
   // Functions to draw
   auto Ldr = [Lbwg,Fpar](double* x,double* p) -> double {
      return Lbwg(x,Fpar.data());
   };
   TF1* fdr = new TF1("fdr", Ldr, dL, dU, 0);
   fdr->SetLineWidth(2);
   fdr->SetLineColor(kRed);
   fdr->SetNpx(500);

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,1000,800);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   SetHstFace(hst);
   hst->SetMinimum(1.);
   hst->GetYaxis()->SetTitleOffset(1.2);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
   fdr->Draw("SAME");

   TLegend* leg = new TLegend(0.55,0.79,0.89,0.89);
   string title(Ename);
   if ( Ename == "3080_rs" ) {
      title = "3080MeV R-scan";
   } else if ( Ename == "J4" ) {
      title = "3097MeV (2018)";
   }
   leg->SetHeader(title.c_str(),"C");
   leg->AddEntry(hst,"MC signal #phi#eta","LEP");
   leg->AddEntry(bwg,"Breit-Wigner #otimes Gauss","L");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.55,0.63,0.89,0.78,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#chi^{2}/ndf = %.2f / %u",chi2,Ndf) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt->AddText( Form("#sigma = %s MeV", PE.Eform(2,".2f",1e3)) );
   pt->AddText( Form("N_{#phi}= %s",PE.Eform(3,".0f")) );
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
double RevArgus(double m, double A, double Eee) {
//--------------------------------------------------------------------
// Argus function with left cutoff (L)
//      U is the upper limit of fitting
//      L,U - must be constants for fitting

   static const double L = 2*Mk;
   const double U = Eee - Meta;

   double x = (U-m)/(U-L); // transformation: (L,U) -> (1,0)
   return Argus(x,A)/(U-L);
}

//----------------------------------------------------------------------
double RevArgusN(double m, double A, double Eee) {
//----------------------------------------------------------------------
// Numeric normalisation on one for range [dL,dU]

   double norm = 0;
   // cash parameters: Eee is the constant
   static double cacheN = 0;
   static double cacheA = 0;
   if ( cacheN > 0 && A == cacheA ) {
      norm = cacheN;
   } else {
      // integrand lambda function
      double p[] = {A,Eee};
      auto Lint = [](double x, void* pp) -> double{
         double* p = static_cast<double*>(pp);
         return RevArgus(x,p[0],p[1]);
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-6, 1.e-6, 1000);
      norm = gsl_int.Integral(Lint,p,dL,dU);
      cacheN = norm;
      cacheA = A;
   }

   return RevArgus(m,A,Eee) / norm;
}

//--------------------------------------------------------------------
void test_RevArgus() {
//--------------------------------------------------------------------

   double Eee = 3.08;
   // Create Unuran 1D distribution object
   auto Lar = [](const double* x,const double* p) -> double {
      return RevArgus(x[0],p[0],p[1]);
   };
   TF1* far = new TF1("far", Lar, dL,dU, 2);
   far->SetParameters(0.,Eee);
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

   // Generate
   TH1D* h1 = new TH1D("h1Ar","", 100,dL,dU);
   h1->SetLineWidth(2);
   h1->SetLineColor(kBlack);
   h1->SetMarkerStyle(20);

   int Ngen = 1000;
   for (int i = 0; i < Ngen; ++i) {
      double x = unr.Sample();
      h1->Fill(x);
   }

   // Draw
   double Wb = Ngen * (dU-dL) / h1->GetNbinsX();
   auto Larg = [Wb](const double* x,const double* p)->double {
      return Wb*RevArgusN(x[0],p[0],p[1]);
   };

   TF1* fbkg = new TF1("fbkg", Larg, dL, dU, 2);
   fbkg->SetParameter(0, far->GetParameter(0) );
   fbkg->SetParameter(1, Eee );
   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);

   TLegend* leg = new TLegend(0.13,0.74,0.43,0.89);
   leg->SetHeader( Form("RevArgus(Eee=%.2f) a=%.2f",
            Eee, fbkg->GetParameter(0)),"C" );
   leg->AddEntry( h1,"Unuran generated","PLE");
   leg->AddEntry( fbkg->Clone(), "RevArgusN function","L");

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   h1->Draw("EP");
   fbkg->DrawCopy("SAME");

   fbkg->SetParameter(0, 2.); // Argus parameter
   fbkg->SetLineColor(kRed+2);
   fbkg->DrawCopy("SAME");
   leg->AddEntry( fbkg->Clone(), "RevArgusN(A=2)","L");

   leg->Draw();
   gPad->RedrawAxis();
   c1->Update();

   // test normalization:
   printf(" RevArgusN norm= %.6f\n",fbkg->Integral(dL,dU)/Wb );
}

// {{{1 Interference BW with Argus bkg.
//--------------------------------------------------------------------
vector<double> IntfrBWAR( double m, double Eee,
      double mphi, double gphi,      // B-W
      double A, double F, double Ang // Argus & interference
      ) {
//--------------------------------------------------------------------
// see BreitWigner() function
// The return vector contains:
//                   0    1     2        3
//        ret     { Sum, B-W, Argus, Interference }
//--------------------------------------------------------------------
   if ( m < dL ) { // phase space becomes zero
      vector<double> ret(4,0.);
      return ret;
   }

   constexpr double Meta2  = SQ(Meta);
   constexpr double Mk2    = SQ(Mk);
   static const double R   = R_BW; // Blatt-Weisskopf ff

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
   double Ar_fun = RevArgus(m,A,Eee);
   double Ar = Ar_fun*SQ(F);

   // interference:
   double Intfr = 2 * tmp*F*sqrt(kap*Ar_fun) *
                  ( (m2-mphi2)*cos(Ang) - GM*sin(Ang) );

   double Sum = BW + Ar + Intfr;
   vector<double> ret { Sum, BW, Ar, Intfr };
   return ret;
}

//--------------------------------------------------------------------
double IntfrBWARG( double m, double Eee,
      const double p[], // {mphi,gphi,sigma,A,F,Ang}
      unsigned int idx = 0 // what to return
      ) {
//--------------------------------------------------------------------
// Function Fun(x) folded with a normal distribution:
//
//        1                        [                  (X'-X)^2    ]
// ---------------- *   Integral   [ Fun(X') * exp( - --------- ) ] dX'
// sqrt(2*pi)*sigma   (-oo<X'<+oo) [                  2*sigma^2   ]
//
//--------------------------------------------------------------------
// integration in +/-5 sigma interval
//--------------------------------------------------------------------
   static const double Ngauss = 1./sqrt(2*M_PI);

   // integrand lambda function
   auto Lint = [m,Eee,p,idx](double t) -> double{
      // t == (X'-X)/sigma
      double x = m+t*p[2];
      vector<double> Fun = IntfrBWAR(x,Eee,p[0],p[1],p[3],p[4],p[5]);
      return exp(-0.5*t*t) * Fun[idx];
   };

   rmath_fun< decltype(Lint) > Fint(Lint);

   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-6, 1.e-5, 1000);
   double result = gsl_int.Integral(Fint,-5.,+5.);

   return Ngauss * result;
}

//--------------------------------------------------------------------
double IntfrBWARGN( double m, double Eee,
      const double p[], //{mphi,gphi,sigma,A,F,Ang}
      int idx = 0       // what to return
      ) {
//--------------------------------------------------------------------
// Numerical normalisation to one for range [dL,dU]

   // cash parameters: Eee is the constant
   static double cacheN = 0;
   constexpr int npar = 6;
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
      auto Lint = [Eee,p](double x) -> double {
         return IntfrBWARG(x,Eee,p,0);  // MUST be zero!
      };
      rmath_fun< decltype(Lint) > Fint(Lint);

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-5, 1.e-4, 1000);
      norm = gsl_int.Integral(Fint,dL,dU);

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

   return IntfrBWARG(m,Eee,p,idx) / norm;
}

//--------------------------------------------------------------------
void test_Intfr() {
//--------------------------------------------------------------------
   double Eee = 3.097; // ~ M(J/Psi)

   // Create Unuran 1D distribution object
   auto Lgen = [Eee](const double* x,const double* p) -> double {
      // vector<double> res =
         // IntfrBWAR( x[0], Eee, p[0],p[1],p[2],p[3],p[4] );
      // return res[0];
      return IntfrBWARG( x[0],Eee, p, 0 );
   };
   TF1* fgen = new TF1("fgen", Lgen, dL,dU, 6); // 5 or 6
   double Sg=1.4e-3, Ag = 0, Fg = 1.0, Ang = 0.8;
   // fgen->SetParameters(Mphi,Gphi,Ag,Fg,Ang);
   fgen->SetParameters(Mphi,Gphi,Sg,Ag,Fg,Ang);

   double Xmax = fgen->GetMaximumX(Mphi-0.4e-3,Mphi+0.4e-3);
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

   // Generate
   TH1D* h1 = new TH1D("h1Intfr","", 100,dL,dU);
   h1->SetLineWidth(2);
   h1->SetLineColor(kBlack);
   h1->SetMarkerStyle(20);

   int Ngen = 1000;
   for (int i = 0; i < Ngen; ++i) {
      double x = unr.Sample();
      h1->Fill(x);
   }

   // Draw
   double Wb = Ngen * (dU-dL) / h1->GetNbinsX();
   auto Lintfr=[Eee,Wb,Sg,Ag,Fg,Ang](double* x, double* p) -> double {
      // vector<double> res = IntfrBWAR( x[0],Eee, Mphi,Gphi,Ag,Fg,Ang);
      // return Wb*res[int(p[0])];
      const double pp[] {Mphi,Gphi,Sg,Ag,Fg,Ang};
      return Wb * IntfrBWARG( x[0],Eee, pp, int(p[0]));
   };

   TF1* fun = new TF1("fun", Lintfr, dL, dU, 1);
   fun->SetLineWidth(2);
   fun->SetNpx(500);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   TLegend* leg = new TLegend(0.59,0.59,0.89,0.89);
   leg->SetHeader("Interference BW #otimes Argus","C");
   leg->AddEntry( h1,"Unuran generated","PLE");

   h1->SetMaximum(1.2*fun->Eval(Xmax));
   h1->SetMinimum(-25);
   h1->Draw("EP");

   fun->SetParameter(0, 0.); // draw all
   fun->SetLineColor(kRed);
   fun->SetLineStyle(kSolid);
   fun->DrawCopy("SAME");
   leg->AddEntry( fun->Clone(), "IntfrBWARG function", "L");

   fun->SetParameter(0, 1); // draw BW
   fun->SetLineColor(kGreen+2);
   fun->SetLineStyle(kDashed);
   fun->DrawCopy("SAME");
   leg->AddEntry( fun->Clone(), "Breit-Wigner", "L");

   fun->SetParameter(0, 2); // draw Argus
   fun->SetLineColor(kBlue);
   fun->SetLineStyle(kDashed);
   fun->DrawCopy("SAME");
   leg->AddEntry( fun->Clone(), "Argus", "L");

   fun->SetParameter(0, 3); // draw interference
   fun->SetLineColor(kMagenta);
   fun->SetLineStyle(kSolid);
   fun->DrawCopy("SAME");
   leg->AddEntry( fun->Clone(), "Interference", "L");

   leg->Draw();
   c1->Update();
}

// {{{1 Fitting without interference: BWG + Argus background
//--------------------------------------------------------------------
void do_fit(string fname, double Eee, bool bkgfit,
      string title, string pdf="") {
//--------------------------------------------------------------------
   // Use binned data
   TH1D* hst = nullptr;
   int Nbin = 0;
   double Nev = 0, Nsb = 0;
   tie(hst,Nbin,Nev,Nsb) = Subtract(fname); // subtract side-band

   ROOT::Fit::DataOptions opt;
   // use integral of bin content instead of bin center (def: false)
   opt.fIntegral = true;
   // use empty bins (default is false) with a fixed error of 1
   opt.fUseEmpty = true; // must be even for LH

   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::BinData Dat(opt,dr);
   ROOT::Fit::FillData(Dat, hst);

   // Fit by Breit-Wigner convoluted with Gauss + Argus bkg
   auto Lfit = [Eee](const double* x,const double* p) -> double {
      return bW * ( p[4]*BreitWignerGaussN(x[0],Eee,p[0],p[1],p[2])
            + p[5]*RevArgusN(x[0],p[3],Eee) );
   };

   vector<string> par_name {"Mphi","Gphi","Sigma","Ar","Nphi","Nbg"};
   vector<double> par_ini  {Mphi,Gphi,1.5e-3, 0., 0.97*Nev,0.03*Nev};
   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 Wfit( *Ffit );
   fitter.SetFunction(Wfit,false); // false=no parameter derivatives

   // Set parameters, data must come before names
   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   //-----------------------------------------------------------------
   fitter.Config().ParSettings(0).SetValue(MphiMC); // like MC
   fitter.Config().ParSettings(0).Fix();
   // fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);

   // fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3,Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();

   // fitter.Config().ParSettings(2).SetValue(1e-3);
   fitter.Config().ParSettings(2).SetLimits(0.1e-3, 3e-3); // Sigma

   fitter.Config().ParSettings(3).SetValue(0.); // Ar
   fitter.Config().ParSettings(3).Fix();

   fitter.Config().ParSettings(4).SetLimits(0.8*Nev,1.2*Nev);// Nphi

   if ( bkgfit ) {
      fitter.Config().ParSettings(5).SetLimits(0.,0.2*Nev); // Nbg
   } else {
      fitter.Config().ParSettings(5).SetValue(0.);
      fitter.Config().ParSettings(5).Fix();
   }

   // Fit
   //-----------------------------------------------------------------
   fitter.LikelihoodFit( Dat );
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res,0.1); // ignore 10% upper/lower errors
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------
   // "Goodness of fit" evaluated by Chi2
   double chi2 = res.Chi2();
   unsigned int Ndf = res.Ndf();
   // printf("Lmin=%g, chi2=%g, Ndof=%u\n", Lmin,chi2,Ndf);

   // do it yourself -> the same!
   // double ch2 = 0;
   // int ndof = -res.NFreeParameters();
   // for ( int ib = 1; ib <= hst->GetNbinsX(); ++ib ) {
      // double vd = hst->GetBinContent(ib);
      // // if ( vd < 1.e-6 ) continue;
      // double ed = hst->GetBinError(ib);
      // if ( fabs(ed) < 1.e-1 ) continue;
      // double bc[] { hst->GetBinCenter(ib) };
      // double fd = Lfit(bc,Fpar.data());
      // ch2 += SQ((vd-fd)/ed);
      // ndof += 1;
   // }
   // printf("my ch2 / Ndof =%g / %d\n",ch2,ndof);

   //-----------------------------------------------------------------
   // Functions to draw
   auto Ldr = [Lfit,Fpar](double* x,double* p) -> double {
      return Lfit(x,Fpar.data());
   };
   TF1* fdr = new TF1("fdr", Ldr, dL, dU, 0);
   fdr->SetLineWidth(2);
   fdr->SetLineColor(kRed);
   fdr->SetNpx(500);

   auto Lbkg = [Fpar,Eee](double* x,double* p) -> double {
      return bW * Fpar[5] * RevArgusN(x[0],Fpar[3],Eee);
   };
   TF1* fbkg = new TF1("fbkg", Lbkg, dL, dU, 0);
   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);
   fbkg->SetLineStyle(kDashed);

   // set min/max to draw
   hst->SetMinimum(-2.);
   int maxbin = hst->GetMaximumBin();
   double hmax = hst->GetBinContent(maxbin)+hst->GetBinError(maxbin);
   double fmax = fdr->GetMaximum( 1.01, 1.03, 1e-6);
   hmax = floor(1.1 * max( hmax, fmax ));
   if ( hmax > 0 ) {
      hst->SetMaximum(hmax);
   }

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,1000,800);
   // TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();
   // gPad->SetLogy(true);

   SetHstFace(hst);
   hst->GetYaxis()->SetMaxDigits(3);
   hst->GetYaxis()->SetTitleOffset(1.2);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
   fdr->Draw("SAME");

   double yleg = 0.80 - 0.04*int(bkgfit);
   TLegend* leg = new TLegend(0.55,yleg,0.89,0.89);
   leg->AddEntry(hst,title.c_str(),"LEP");
   leg->AddEntry(fdr,"Breit-Wigner #otimes Gauss","L");
   if ( bkgfit ) {
      fbkg->Draw("SAME");
      leg->AddEntry(fbkg,"flat background","L");
   }
   leg->Draw();

   double ypt = yleg - 0.04*(5+int(bkgfit));
   TPaveText* pt = new TPaveText(0.55,ypt,0.89,yleg-0.01,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#chi^{2}/ndf = %.2f / %u",chi2,Ndf) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",
            PE.Eform(1,".3f",1e3)) );
   pt->AddText( Form("#sigma= %s MeV", PE.Eform(2,".2f",1e3)) );
   pt->AddText( Form("N_{#phi}= %s",PE.Eform(4,".1f")) );
   if ( bkgfit ) {
      // pt->AddText( Form("a = %s",PE.Eform(3,".2f")) );
      pt->AddText( Form("N_{bkg}= %s",PE.Eform(5,".1f")) );
   }
   pt->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      if ( !bkgfit ) {
         pdf += "_nobg";
      }
      pdf += ".pdf";
      c1->Print(pdf.c_str());
   }

   // print for table
   printf(" No Interference: %s\n",title.c_str());
   printf(" Data: Ndat= %.1f Nsb= %.1f\n",Nev,Nsb);
   printf(" Fit: chi^2/Ndf= %.2f / %u\n", chi2, Ndf);
   printf("      sigma= %s MeV  F= %s\n",
         PE.Tform(2,".2f",1e3), PE.Tform(4,".2f") );
   printf("      Nphi= %s  Nbg= %s\n",
         PE.Tform(4,".1f"), PE.Tform(5,".1f") );
}

// {{{1  Fit Interference(BW + RevArgus)
// {{{2 class chi^2 for root fitter
//--------------------------------------------------------------------
class myChi2 {
   public:
      bool prtchi = false;
      //--------------------------------------------------------------
      myChi2(const TH1D* hst, TF1* f): func(f) {
      //--------------------------------------------------------------
         // save bins with nonzero error
         size_t Nb = hst->GetNbinsX();
         Hbc.reserve(Nb);
         Hbe.reserve(Nb);
         Hx1.reserve(Nb);
         Hx2.reserve(Nb);
         for ( int ib = 1; ib <= Nb; ++ib ) {
            double bx = hst->GetBinCenter(ib);
            if ( bx < dL || bx > dU ) {
               continue;
            }
            double bw = hst->GetBinWidth(ib);
            double bc = hst->GetBinContent(ib);
            double be = hst->GetBinError(ib);
            if ( fabs(be) > 1e-3 ) {
               Hbc.push_back(bc);
               Hbe.push_back(be);
               Hx1.push_back(bx-0.5*bw);
               Hx2.push_back(bx+0.5*bw);
               printf("myChi2: bc=%g, be=%g, bx=%g\n", bc,be,bx);
            } else {
               // debug print
               if ( fabs(bc) > 1e-3 ) {
                  printf("myChi2::ERROR: bc=%g, be=%g, bx=%g\n",
                        bc,be,bx);
               }
            }
         }
      }

      size_t Size() const {
         return Hbc.size();
      }

      // chi^2 - function
      // the signature of this operator() MUST be exactly this:
      //--------------------------------------------------------------
      double operator() (const double* p) {
      //--------------------------------------------------------------
         const double maxResValue = DBL_MAX / 1000;
         const double epsrel = 1e-6;

         func->SetParameters(p);

         double chi2 = 0;
         size_t Nb = Size();
         for ( int i = 0; i < Nb; ++i ) {
            double bc = Hbc[i];
            double be = Hbe[i];
            double x1 = Hx1[i];
            double x2 = Hx2[i];

            // integral over bin
            double err_int = 0;
            double fval = func->IntegralOneDim( x1, x2,
                  epsrel, be*1e-1, err_int ) / (x2-x1);

            double resval = SQ( (bc-fval) / be );

            if ( prtchi ) {
               printf("chi2(x=%.3f)=%.1f  [%.2f, %.2f]/{%.2f}\n",
                     0.5*(x1+x2),resval, bc,fval, be);
            }

            // avoid inifinity or nan in chi2
            if( resval < maxResValue ) {
               chi2 += resval;
            } else {
               chi2 += maxResValue;
            }
         }

         return chi2;
      }

   private:
      TF1* func;
      vector<double> Hbc;
      vector<double> Hbe;
      vector<double> Hx1;
      vector<double> Hx2;
};

// {{{2 do_fitI
//--------------------------------------------------------------------
void do_fitI(string fname, double Eee, string title, string pdf="") {
//--------------------------------------------------------------------
   bool is2012 = (title.find("2012") != string::npos);
   bool is2018 = (title.find("2018") != string::npos);
   bool isR = (title.find("R-scan") != string::npos);

   // Use binned data
   TH1D* hst = nullptr;
   int Nbin = 0;
   double Nev = 0, Nsb = 0;
   tie(hst,Nbin,Nev,Nsb) = Subtract(fname); // subtract side-band

   ROOT::Fit::DataOptions opt;
   // use integral of bin content instead of bin center (def: false)
   opt.fIntegral = true;
   // use empty bins (default is false) with a fixed error of 1
   opt.fUseEmpty = true; // for LH ?

   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::BinData Dat(opt,dr);
   ROOT::Fit::FillData(Dat, hst);

   //-----------------------------------------------------------------
   // Fit signal MC(phi,eta) by Breit-Wigner convoluted with Gauss
   // function is normalized to 1 on the fit range
   auto Lfit = [Eee](const double* x,const double* p) -> double {
      return bW * p[6]*IntfrBWARGN( x[0], Eee, p, 0 );
   };

   vector<string> par_name
      { "Mphi", "Gphi", "Sigma", "Ar", "F", "vartheta", "NKK" };
   vector<double> par_ini
      {  Mphi,   Gphi,   1e-3,    0.,   0.7,  0.,        Nev };

   // bool neg_pos = par_ini[5] > 0; // true for negative interference

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* Ffit = new TF1("Ffit", Lfit, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 Wfit( *Ffit );
   fitter.SetFunction( Wfit, false); // false=no parameter derivatives

   // Set parameters, data must come before names
   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   //-----------------------------------------------------------------
   // fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   fitter.Config().ParSettings(0).SetValue(MphiMC); // like MC
   fitter.Config().ParSettings(0).Fix();

   // fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3,Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();

   fitter.Config().ParSettings(2).SetLimits(0.2e-3, 5e-3); // Sigma

   fitter.Config().ParSettings(3).SetValue(0.);             // Ar
   fitter.Config().ParSettings(3).Fix();

   fitter.Config().ParSettings(4).SetLimits(0., 5.);       // F

   // fitter.Config().ParSettings(5).SetLimits(-M_PI,M_PI);//vartheta
   fitter.Config().ParSettings(5).SetValue(0.); // MEMO
   // fitter.Config().ParSettings(5).SetValue(0.58); // + 1sigma
   // fitter.Config().ParSettings(5).SetValue(-0.58); // - 1sigma
   fitter.Config().ParSettings(5).Fix();

   // NetaKK limits
   if ( Nev < 33 ) {
      fitter.Config().ParSettings(6).SetLimits(1.,50.);
   } else {
      fitter.Config().ParSettings(6).SetLimits(0.5*Nev,1.5*Nev);
   }

   // Fix sigma for small number of events
   if ( Nev < 15 ) { // fix sigma ~ 1e-3
      fitter.Config().ParSettings(2).SetValue(1e-3);
      fitter.Config().ParSettings(2).Fix();
      if ( Nbin < 3 ) { // fix F to 0
         fitter.Config().ParSettings(4).SetValue(0.);
         fitter.Config().ParSettings(4).Fix();
      }
   }

   // Fit
   //-----------------------------------------------------------------
   // minimization 'function' == LeastSquareFit( Dat )
   // myChi2 chi2_fit( hst, Ffit);
   // int Nd = chi2_fit.Size();
   // fitter.FitFCN(Npar, chi2_fit, nullptr, Nd, true);

   // fitter.LeastSquareFit( Dat );
   fitter.LikelihoodFit( Dat );
   fitter.CalculateHessErrors(); // for err_Nphi
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double Lmin = res.MinFcnValue();
   ParAndErr PE(res,0.1); // ignore 10% upper/lower errors
   vector<double>& Fpar = PE.Fpar;

   //-----------------------------------------------------------------
   // "Goodness of fit" evaluated by Chi2
   double chi2 = res.Chi2();
   unsigned int Ndf = res.Ndf();

   //-----------------------------------------------------------------
   // print chi^2 point by point
   // printf("\n");
   // chi2_fit.prtchi = true;
   // chi2_fit( PE.Fpar.data() );
   // chi2_fit.prtchi = false;
   // printf("\n");

   //-----------------------------------------------------------------
   // Calculate Nphi with error
   ValErr CalcIntEr( const ROOT::Fit::FitResult& res,
         double Eee, unsigned int idx ); // prototype

   ValErr Integ = CalcIntEr( res, Eee, 1 ); // BW
   double Nphi = Integ.val;
   double err_Nphi = fabs(Integ.err);
   if ( SQ(err_Nphi) < Nphi ) {
      printf(" WARNING: Too small error of Nphi: %.1f +/- %.1f\n",
            Nphi, err_Nphi);
      err_Nphi = sqrt( Nphi );
   }

   //-----------------------------------------------------------------
   // Functions to draw
   auto Ldr = [Eee,Fpar](double* x,double* p) -> double {
      int idx = int(p[0]);
      return bW * Fpar[6]*IntfrBWARGN( x[0],Eee,Fpar.data(),idx );
   };
   TF1* fdr = new TF1("fdr", Ldr, dL, dU, 1);
   fdr->SetLineWidth(2);
   fdr->SetLineColor(kRed);
   fdr->SetNpx(500);

   // set max/min to draw
   if ( Nev < 15 ) { // use defaults
      hst->SetMaximum();
   } else {
      fdr->SetParameter(0, 1); // BW
      int maxbin = hst->GetMaximumBin();
      double hmax =
         hst->GetBinContent(maxbin)+hst->GetBinError(maxbin);
      double fmax = fdr->GetMaximum( 1.01, 1.03, 1e-6);
      hmax = floor(1.1 * max( hmax, fmax ));
      hst->SetMaximum(hmax);
   }

   fdr->SetParameter(0, 3);   // interference
   double fmin = fdr->GetMinimum( 1.01, 1.03, 1e-6);
   double hmin = min( -1.5, 1.15*fmin );
   hst->SetMinimum(hmin);

   //-----------------------------------------------------------------
   // TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   TCanvas* c1 = new TCanvas("c1","...",0,0,1000,800);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst);
   hst->GetYaxis()->SetMaxDigits(3);
   hst->GetYaxis()->SetTitleOffset(1.2);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");

   double yleg = 0.89 - 0.03*5;
   // TLegend* leg = new TLegend(0.55,yleg,0.89,0.89);
   TLegend* leg = new TLegend(0.60,yleg,0.89,0.89);
   leg->AddEntry(hst,title.c_str(),"LEP");

   fdr->SetParameter(0, 0); // SUM
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Result of the fit", "L");

   // draw components
   fdr->SetLineWidth(1);
   fdr->SetLineStyle(kDashed);

   fdr->SetParameter(0, 1);   // BW
   fdr->SetLineColor(kGreen+2);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Breit-Wigner", "L");

   fdr->SetParameter(0, 2);   // Argus
   fdr->SetLineColor(kBlue);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Argus", "L");

   fdr->SetParameter(0, 3);   // interference
   fdr->SetLineColor(kMagenta+1);
   fdr->DrawCopy("SAME");
   leg->AddEntry( fdr->Clone(), "Interference", "L");

   leg->Draw();

   double ypt = yleg - 0.03*8;
   // TPaveText* pt = new TPaveText(0.55,ypt,0.89,yleg-0.01,"NDC");
   TPaveText* pt = new TPaveText(0.60,ypt,0.89,yleg-0.01,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#chi^{2}/ndf= %.2f / %u",chi2,Ndf) );
   pt->AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt->AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt->AddText( Form("#sigma= %s MeV", PE.Eform(2,".2f",1e3)) );
   // pt->AddText( Form("a= %s",PE.Eform(3,".1f")) );
   pt->AddText( Form("F= %s",PE.Eform(4,".2f")) );
   pt->AddText( Form("#vartheta= %s",PE.Eform(5,".2f")) );
   // pt->AddText( Form("N_{#etaKK}= %s",PE.Eform(6,".1f")) );
   pt->AddText( Form("N_{KK}= %s",PE.Eform(6,".1f")) ); // memo
   pt->AddText( Form("N_{#phi}= %.1f #pm %.1f",Nphi,err_Nphi) );
   pt->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      // pdf = string("mkk_pict/") + pdf +
         // (neg_pos ? "_n.pdf" : "_p.pdf");
      pdf += ".pdf";
      c1->Print(pdf.c_str());
   }

   // print for table
   printf(" Interference: %s\n",title.c_str());
   printf(" Data: Ndat= %.1f Nsb= %.1f\n",Nev,Nsb);
   printf(" Fit: chi^2/Ndf= %.2f / %u\n", chi2, Ndf);
   printf("      sigma= %s MeV  F= %s\n",
         PE.Tform(2,".2f",1e3), PE.Tform(4,".2f") );
   printf("      Nphi= %.1f \\pm %.1f  N_etaKK= %s\n",
         Nphi, err_Nphi, PE.Tform(6,".1f") );
}

//--------------------------------------------------------------------
ValErr CalcIntEr( const ROOT::Fit::FitResult& res,
                  double Eee, unsigned int idx ) {
//--------------------------------------------------------------------
// Numerical calculation of integrals and their errors
   // desired errors:
   constexpr double eps_abs = 1.e-5;
   constexpr double eps_rel = 1.e-5;
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
      return p[6]*IntfrBWARGN( x[0], Eee, p, idx );
   };
   TF1 FC("FC", Lfc, dL, dU, Npar);
   FC.SetParameters(Fpar.data());

   double err_num = 0;
   double Integ = FC.IntegralOneDim(dL,dU,eps_rel,eps_abs, err_num);
   // printf("++ Integ[%u]= %g, rel_err= %g\n",idx,Integ,err_num);
   if ( err_num > epsilon ) {
      cout << " WARNING: " << __func__ << " numerical error of"
         " integration of F_C(" << idx << ") is too big "
         << err_num << endl;
   }

   // DEBUG:
   // double errD = 0;
   // double IntD = FC.IntegralOneDim(1.01,1.03,eps_rel,eps_abs, errD);
   // printf("++ IntD[%u]= %g, relD= %g\n",idx,IntD,errD);

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
      double e_rel = ((i!=2) ? eps_rel : 1e-3); // sigma =~ 1e-3
      IntGrad[i] = dF.IntegralOneDim(dL,dU,e_rel,eps_abs, err_num);
      err_num = covMatrix(i,i) * IntGrad[i] * err_num;
      // printf("++++ IntGrad(%i)= %g +/- %g\n",i,IntGrad[i],err_num);
      err_num2 += SQ(err_num);
   }

   double err_Int = sqrt( covMatrix.Similarity(IntGrad) );
   err_num2 = sqrt( err_num2 / err_Int ); // abs numerical error
   // printf("++ err_Int[%u]= %g, err_num2= %g\n",idx,err_Int,err_num2);
   if ( err_num2 > epsilon ) {
      cout << " WARNING: " << __func__ << " numerical error "
         "of integration of dF is too big " << err_num2 << endl;
   }

   return (ValErr {Integ,err_Int});
}

// {{{1 MAIN:
//--------------------------------------------------------------------
void mass_KK_fit(int interference=1) {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(111); // do not print fixed parameters!
   // gStyle->SetOptFit(112); // print all parameters (fixed)
   gStyle->SetStatFont(62);
   gStyle->SetLegendFont(42);

   // gStyle->SetStatH(0.25);
   // gStyle->SetStatW(0.23);

   //-----------------------------------------------------------------
   // = define GSL error handler which does nothing =
   // gsl_set_error_handler_off();

   // set integrator: ROOT::Math::GSLIntegrator adaptive method (QAG)
   ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Adaptive");
   ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1e-9);

   //-----------------------------------------------------------------
   // -> test functions
   // test_BreitWigner();
   // test_RevArgus();
   // test_Intfr();

   //-----------------------------------------------------------------
   // -> fit signal MC

   string Ename("3080_rs");
   double Eee = 3.08;
   // string Ename("J4");
   // double Eee = 3.096986;
   // --------------
   // string Ename("J3");
   // double Eee = 3.096203;
   // string Ename("3097"); // 2012
   // double Eee = 3.097777;


   // string pdf_mcsig( Form("mcMkk_sig_%s.pdf",Ename.c_str()) );
   // sig_fit_mc(Ename,Eee,pdf_mcsig);

   // fit MC phi-eta by BW x Gauss
   // string pdf_sig( Form("mkk_sig_log_%s.pdf",Ename.c_str()) );
   // sig_fit(Ename,Eee,pdf_sig);

   // return;

   //-----------------------------------------------------------------
   // structure with data description
   //-----------------------------------------------------------------
   struct DataDescription {
      string Basename;
      double Energy;    // GeV
      string Title;
   };
   const vector<DataDescription> DD {
      // -- 2019 3080 no lum! 0
      // {"3080_2019",3.08,"2019: 3080MeV"},
      //-- R-scan 2015: 1-6
      // {"2900_rs",2.90, "R-scan: 2900MeV"},
      // {"2950_rs",2.95, "R-scan: 2950MeV"},
      // {"2981_rs",2.981,"R-scan: 2981MeV"},
      // {"3000_rs",3.0,  "R-scan: 3000MeV"},
      // {"3020_rs",3.02, "R-scan: 3020MeV"},
      {"3080_rs",3.08, "R-scan: 3080MeV"},
      //-- tau-scan 2018: 7-14
      {"J1", 3.087593, "2018: 3087.593MeV"},
      {"J2", 3.095726, "2018: 3095.726MeV"},
      {"J3", 3.096203, "2018: 3096.203MeV"},
      {"J4", 3.096986, "2018: 3096.986MeV"},
      {"J5", 3.097226, "2018: 3097.226MeV"},
      {"J6", 3.097654, "2018: 3097.654MeV"},
      {"J7", 3.098728, "2018: 3098.728MeV"},
      {"J8", 3.104,    "2018: 3104.000MeV"},
      //-- J/Psi-scan 2012: 15-30
      {"3050",3.049642,"2012: 3049.642MeV"},
      {"3060",3.058693,"2012: 3058.693MeV"},
      {"3080",3.079645,"2012: 3079.645MeV"},
      {"3083",3.082496,"2012: 3082.496MeV"},
      {"3090",3.088854,"2012: 3088.854MeV"},
      {"3093",3.091760,"2012: 3091.760MeV"},
      {"3094",3.094697,"2012: 3094.697MeV"},
      {"3095",3.095430,"2012: 3095.430MeV"},
      {"3096",3.095826,"2012: 3095.826MeV"},
      {"3097",3.097213,"2012: 3097.213MeV"},
      {"3098",3.098340,"2012: 3098.340MeV"},
      {"3099",3.099042,"2012: 3099.042MeV"},
      {"3102",3.101359,"2012: 3101.359MeV"},
      {"3106",3.105580,"2012: 3105.580MeV"},
      {"3112",3.112051,"2012: 3112.051MeV"},
      {"3120",3.119878,"2012: 3119.878MeV"},
   };

   //-----------------------------------------------------------------
   // Fit: Interference(BW + RevArgus)
   //-----------------------------------------------------------------
   if ( interference == 1 ) {
      for ( const auto& dd : DD ) {
         string Filename = "ntpl_" + dd.Basename + ".root";
         string outdir("mkk_inter/Std2/");
         // string outdir("mkk_inter/Ch2_60/");
         // string outdir("mkk_inter/Weta2/");
         // string outdir("mkk_inter/AngM/");
         string pdfout = outdir + "mkk_" + dd.Basename; //+".pdf"
         string txtout = outdir + "mkk_" + dd.Basename + ".txt";
         fflush(stdout);
         // save fd
         int fd_save1 = dup( fileno(stdout) );
         int fd_save2 = dup( fileno(stderr) );
         // redirection stdout,stderr -> file
         freopen(txtout.c_str(),"w",stdout);
         if ( -1 == dup2(fileno(stdout), fileno(stderr)) ) {
            perror("ERROR of redirection\n");
         }
         do_fitI( Filename, dd.Energy, dd.Title, pdfout );
         fflush(stdout); // do not close!
         // restore stdout, stderr
         dup2(fd_save1, fileno(stdout));
         close(fd_save1);
         dup2(fd_save2, fileno(stderr));
         close(fd_save2);
         printf(" The end of %s\n", dd.Basename.c_str());
      }
      return;

   //-----------------------------------------------------------------
   // Fitting without interference: BWG + Argus background
   //-----------------------------------------------------------------
   } else if ( interference == 0 ) {
      // bool bkgfit=true;
      bool bkgfit = false; // only BWG

      for ( const auto& dd : DD ) {
         string Filename = "ntpl_" + dd.Basename + ".root";
         string outdir("mkk_noint/");
         string pdfout = outdir + "mkk_" + dd.Basename; //+".pdf"
         string txtout = outdir + "mkk_" + dd.Basename + ".txt";
         // redirection stdout,stderr -> file
         FILE* fout = freopen(txtout.c_str(),"w",stdout);
         if ( -1 == dup2(fileno(fout), fileno(stderr)) ) {
            perror("ERROR of redirection\n");
         }
         do_fit( Filename, dd.Energy, bkgfit, dd.Title, pdfout );
         fclose(fout);
         printf(" The end of %s\n", dd.Basename.c_str());
      }
      return;
   }

   //-----------------------------------------------------------------
   // Debug: outputs in stdout, pdf in the curent dir
   //-----------------------------------------------------------------
   bool bkgfit = false;
   // do_fit( "ntpl_J4.root", 3.096986, bkgfit,
         // "2018: 3096.986MeV","mkk_J4_M" ); // 866
   // do_fit( "ntpl_3097.root", 3.097227, bkgfit,
         // "2012: 3097.227MeV","mkk_3097" ); // 559
   // do_fit( "ntpl_3080_2019.root", 3.08, bkgfit,
         // "2019: 3080MeV","mkk_3080_2019"); // 159
   // do_fit( "ntpl_3080_rs.root", 3.08, bkgfit,
         // "R-scan: 3080MeV","mkk_3080_rs"); // 224


   // do_fitI( "ntpl_3080_rs.root", 3.08,
         // "R-scan: 3080MeV","mkk_3080_rs");     // 224, sb=4
   // do_fitI( "ntpl_3080.root", 3.079645,
         // "2012: 3079.645MeV", "mkk_3080");     // 35, sb=0

   do_fitI( "ntpl_J1.root", 3.087659,
         "2018: 3087.659MeV","mkk_J1" ); // 8,sb=1
   // do_fitI( "ntpl_J3.root", 3.096203,
         // "2018: 3096.203MeV","mkk_J3" ); // 985,sb=31
   // do_fitI( "ntpl_J4.root", 3.096986,
         // "2018: 3096.986MeV","mkk_J4" ); // 866

   // do_fitI( "ntpl_3050.root", 3.049663,
         // "2012: 3049.663MeV","mkk_3050" ); // 22
   // do_fitI( "ntpl_3095.root", 3.095444,
         // "2012: 3095.444MeV","mkk_3095" ); // 100
   // do_fitI( "ntpl_3097.root", 3.097227,
         // "2012: 3097.227MeV","mkk_3097" ); // 559,sb=16
   // do_fitI( "ntpl_3099.root", 3.099056,
         // "2012: 3099.056MeV","mkk_3099" ); // 24
   // do_fitI( "ntpl_3106.root", 3.105594,
         // "2012: 3105.594MeV","mkk_3106" ); // 13
}
