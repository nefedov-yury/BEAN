// plot M(K+K-) for data
// using cuts (see cuts.h): Mrec + chi^2(4C) + Mgg
// and fit M(K+K-) by Breit-Wigner convoluted with Gauss and
// reversed ARGUS function[(arXiv:1805.05613v3) p.7] as background
// -> mkkYY_fit??.pdf

// {{{1 constants && helper function
//----------------------------------------------------------------------
// adapter class to use lambda functions with closures in gsl_int
// Using:  gsl_fun<decltype(lambda)> Functor(Lambda);
//----------------------------------------------------------------------
template< typename F >
class gsl_fun: public ROOT::Math::IBaseFunctionOneDim {
//----------------------------------------------------------------------
   public:
      gsl_fun(const F& lam) : _f(lam) {}
      gsl_fun* Clone() const {
         return nullptr; // do something more smart
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

static const double Mjpsi = 3.096916; // 3096.916 +/- 0.011 MeV
static const double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
static const double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV
static const double Meta  = 0.547862; // 547.862  +/- 0.017 MeV
static const double Mk    = 0.493677; // 493.677  +/- 0.016 MeV

static const double Mjpsi2 = SQ(Mjpsi);
static const double Meta2  = SQ(Meta);
static const double Mk2    = SQ(Mk);
static const double R      = 3.; // GeV^{-1} for Blatt-Weisskopf ff

static const double dL = 2*Mk; // the left cutoff = 0.987354

// binning
// --------------- old -------------------------
// static const double bL = 0.98; // first bin < dL
// static const double dU = 1.101;
// ---------------------------------------------
// ------------- background only ---------------
// static const double bL = 0.98; // first bin < dL
// static const double dU = 1.18;
// static const int Nbins = 100;
// ---------------------------------------------
//
static const double dU = 1.08; // MUST BE < 1.0835 !!!
// ----------------- bW = 1.1 MeV [standard] --------------
static const double bL = 0.981; // first bin < dL
static const int Nbins = 90;
// ---------------------------------------------
// ----------------- bW = 1.2 MeV --------------
// static const double bL = 0.9804; // first bin < dL
// static const int Nbins = 83;
// static const double bL = 0.98; // first bin < dL
// static const int Nbins = 80; // bW=1.25
// ---------------------------------------------
// ----------------- bW = 1.0 MeV --------------
// static const double bL = 0.98; // first bin < dL
// static const int Nbins = 100;
// ---------------------------------------------
// ----------------- bW = 0.9 MeV --------------
// static const double bL = 0.981; // first bin < dL
// static const int Nbins = 110;
// ---------------------------------------------
//
static const double bW = (dU-bL)/Nbins; // bin width

// histograms for debug
static bool DEBUG = false;
// static vector<TH1D*> htt(10,nullptr);
static TH1D* htt[10];

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

// {{{1 DATA
//-------------------------------------------------------------------------
// DATA
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
TH1D* plot_Mkk(string fname, string hname, int type = 0) {
//-------------------------------------------------------------------------
#include "cuts.h"

   // name of folder with root files
   static string dir("prod-9/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   string title(";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2}");
   title += string(";Entries/") + string(Form("%.2fMeV/c^{2}",bW*1e3));
   TH1D* mkk = new TH1D(hname.c_str(), title.c_str(), Nbins,bL,dU);

   mkk->Sumw2(true);
   TCut c_here = c_Mrec+c_chi2;
   if ( type == 0 ) {   // central part
      c_here += c_cpgg;
   } else {             // side-band
      c_here += c_sbgg;
   }
   string dr = string("Mkk>>") + hname;
   a4c->Draw(dr.c_str(),c_here,"goff");

   return mkk;
}

//-------------------------------------------------------------------------
TH1D* Subtruct(string fname) {
//-------------------------------------------------------------------------
// Side-band subtraction
   TH1D* hst1 = plot_Mkk(fname,string("mkk_data_cp"),0);
   TH1D* hst2 = plot_Mkk(fname,string("mkk_data_sb"),1);
   TH1D* sub = (TH1D*)hst1->Clone("sub");
   sub->Add(hst1,hst2,1.,-1.);

   // errors for zero and negative bins
   for ( int ib = 1; ib <= sub->GetNbinsX(); ++ib ) {
      if ( sub->GetBinCenter(ib) < dL-bW ) {
         continue;
      }
      if ( sub->GetBinContent(ib) > 1.e-3 ) {
         continue;
      }

      if( sub->GetBinError(ib) < 1e-3 ) {
         sub->SetBinError(ib,1.);
      }

      // print "not positive" bins for (hst1-hst2)
//       cout << " bin# " << ib
//            << " mkk= " << sub->GetBinCenter(ib) << " -> " << endl
//            << "\t h1: " << hst1->GetBinContent(ib)
//            << " +/- " << hst1->GetBinError(ib) << endl
//            << "\t h2: " << hst2->GetBinContent(ib)
//            << " +/- " << hst2->GetBinError(ib) << endl
//            << "\t sub:" << sub->GetBinContent(ib)
//            << " +/- " << sub->GetBinError(ib) << endl;
   }
   return sub;
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
   title += string(";Entries/1MeV/c^{2}");
   mcmkk->SetTitle(title.c_str());
//    mcmkk->Sumw2(true);

   return mcmkk;
}

// {{{1 Breigt Wigner for phi -> KK
//----------------------------------------------------------------------
// Breigt Wigner for phi -> KK
//----------------------------------------------------------------------

//----------------------------------------------------------------------
double BreitWigner(double m, double mphi, double gphi) {
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
   double BB_k = sqrt( (1+SQ(R*p_k0)) / (1+SQ(R*p_k)) );

   double fm = r_phi*r_k*BB_phi*BB_k;

   double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k; // == m*G(m)

   return (kap * SQ(fm)) / ( SQ(m2-mphi2) + SQ(GM) );
}

//----------------------------------------------------------------------
double BreitWignerGauss( double m,
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
   static const double Ngauss = 1./sqrt(2*M_PI);

   // integrand lambda function
   double pp[] = { m, mphi, gphi, sigma };
   auto Lint = [](double t, void* pp) -> double{ // t == (X'-X)/sigma
      double* p = static_cast<double*>(pp);
      return exp(-0.5*t*t) * BreitWigner( (p[0]+t* p[3]),p[1],p[2]);
   };

   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);
   double result = gsl_int.Integral(Lint,pp,-5.,+5.);

   // int status = gsl_int.Status();
   if ( DEBUG ) {
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
                          double mphi, double gphi, double sigma
                        ) {
//----------------------------------------------------------------------
// Numeric normalisation on one for range [dL,dU]

   double norm = 0;
   // cash parameters
   static double cacheN = 0;
   static double cacheM = 0;
   static double cacheG = 0;
   static double cacheS = 0;
   if ( cacheN > 0 &&
         mphi == cacheM && gphi == cacheG && sigma == cacheS ) {
      norm = cacheN;
   } else {
      // integrand lambda function
      double p[] = {mphi,gphi,sigma};
      auto Lint = [](double x, void* pp) -> double{
         double* p = static_cast<double*>(pp);
         return BreitWignerGauss(x,p[0],p[1],p[2]);
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
   }

   return BreitWignerGauss(m,mphi,gphi,sigma) / norm;
}

//-------------------------------------------------------------------------
void test_BreitWigner() {
//-------------------------------------------------------------------------

   // 1) BreitWigner (internally normalized to 1)
//    cout << " BW(1.02)= " << BreitWigner(1.02,Mphi,Gphi) << endl;
   auto Lbw = [](const double* x,const double* p) -> double {
      return p[0]*BreitWigner(x[0],p[1],p[2]);
   };
   TF1* bw = new TF1("bw", Lbw, 0.98, dU, 3);
   bw->SetParNames("Norm","Mphi","Gphi");
   bw->SetParameters(1., Mphi, Gphi);
   bw->SetLineWidth(1);
   bw->SetLineColor(kBlue);
   bw->SetNpx(500);

   double norm = 1./bw->Integral(dL,dU,1e-8);
   printf("norm = %.6f\n",norm);
   bw->SetParameter(0, norm );

   // test of stability
//    for ( int i = 0; i < 100; i++ ) { // mphi = Mphi +/- 1MeV
//       double mphi = Mphi + (-1. + 0.02*i)*1e-3;
//       for ( int j = 0; j < 100; j++ ) { // sigma from 1.1 to 1.3 MeV
//          double sigma = (1.1 + 0.02*j)*1e-3;
//          BreitWignerGaussN(mphi-4*sigma,mphi,Gphi,sigma);
//       }
//    }

   // 2) BreitWignerGaussN
   cout << " BWG(1.02)= " << BreitWignerGaussN(1.02,Mphi,Gphi,1e-3) << endl;
   auto Lbwg = [](const double* x,const double* p) -> double {
      return p[0]*BreitWignerGaussN(x[0],p[1],p[2],p[3]);
   };
   TF1* bwg = new TF1("bwg", Lbwg, 0.98, dU, 4);
   bwg->SetParNames("Norm","Mphi","Gphi","Sigma");
   bwg->SetParameters(1., Mphi, Gphi, 1.2e-3);
   bwg->SetLineWidth(2);
   bwg->SetLineColor(kRed);
   bwg->SetNpx(500);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad->SetGrid();
   bw->Draw();
   bwg->Draw("SAME");
   bw->Draw("SAME");
}

//-------------------------------------------------------------------------
void sig_fit(string fname, string title, string pdf="") {
//-------------------------------------------------------------------------
   TH1D* hst = plot_Mkk(fname,string("mkk_phi"));

   // Fit signal MC(phi,eta) by Breit-Wigner convoluted with Gauss
   auto Lbwg = [](const double* x,const double* p) -> double {
      return bW*p[0]*BreitWignerGaussN(x[0],p[1],p[2],p[3]);
   };
   TF1* bwg = new TF1("bwg", Lbwg, 0.98, dU, 4);
   bwg->SetParNames("Nphi","Mphi","Gphi","Sigma");
   bwg->SetParameters(1.e4, Mphi, Gphi, 1.2e-3);

   bwg->FixParameter(1, Mphi);
   bwg->FixParameter(2, Gphi);

   bwg->SetLineWidth(1);
   bwg->SetLineColor(kRed);
   bwg->SetNpx(500);

   gStyle->SetOptFit(111); // do not print fixed parameters!

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();
   gStyle->SetFitFormat(".3g");
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);

   SetHstFace(hst);
//    hst->GetYaxis()->SetMaxDigits(3);
   hst->GetYaxis()->SetTitleOffset(1.25);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
   hst->Fit("bwg","E","",dL,1.08);

   TLegend* leg = new TLegend(0.60,0.57,0.89,0.72);
   leg->AddEntry(hst,title.c_str(),"LEP");
   leg->AddEntry(bwg,"BW #otimes Gauss","L");
   leg->Draw();

//    double Nsig = bwg->GetParameter(0)/dx;
   double Nsig = bwg->Integral(dL,1.08,1e-6)/bW;
   double Nsig_err = Nsig * bwg->GetParError(0)/bwg->GetParameter(0);
   string Tsig( Form("N(#phi#eta) = %.0f #pm  %.0f",Nsig,Nsig_err) );
   cout << Tsig << endl;

//    TPaveText* pt = new TPaveText(0.65,0.48,0.95,0.58,"NDC");
//    TText* t1 = pt->AddText( Tsig.c_str() );
//    t1->SetTextColor(kRed);
//    pt->Draw();

   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void sig_fit_mc(string fname, string pdf="") {
//-------------------------------------------------------------------------
   TH1D* hst = get_mcMkk(fname);
   double bin_width = 1e-3;

   // Fit signal MC(phi,eta) by Breit-Wigner
   auto Lbwg = [bin_width](const double* x,const double* p) -> double {
      return bin_width*p[0]*BreitWigner(x[0],p[1],p[2]);
   };
   TF1* bwg = new TF1("bwg", Lbwg, 0.98, 2., 3);
   bwg->SetParNames("Norm","Mphi","Gphi");
   bwg->SetParameters(1.e5, Mphi, Gphi);

//    bwg->FixParameter(1, Mphi);
//    bwg->FixParameter(2, Gphi);

   bwg->SetLineWidth(2);
   bwg->SetLineColor(kRed);
   bwg->SetNpx(500);

   gStyle->SetOptFit(111); // do not print fixed parameters!

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();
   gStyle->SetFitFormat(".6g");
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

   double genEv = hst->GetEntries();
   double corr = bwg -> GetParameter(0) / genEv;
   double dcorr = bwg -> GetParError(0) / genEv;
   cout << " genEv= " << genEv << " corr= " << corr << endl;

//    gStyle->SetLegendFont(62);
   gStyle->SetLegendTextSize(0.03);
   TLegend* leg = new TLegend(0.22,0.20,0.62,0.40);
   if ( genEv > 8e5 ) {
      leg->AddEntry(hst,"MC truth 2012","L");
   } else {
      leg->AddEntry(hst,"MC truth 2009","L");
   }
   leg->AddEntry(bwg,"Breit-Wigner","L");
   leg->AddEntry(lB,"Cut-off","L");
   string txt(Form(" Cut-off corr.: %.3f #pm %.3f",corr,dcorr));
   leg -> SetHeader( txt.c_str(), "C" );
   leg->Draw();

   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
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

   // normalization constant
   static const double one_over_sqrt2pi = 0.5*M_2_SQRTPI*M_SQRT1_2;
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
double RevArgusN(double m, double A) {
//----------------------------------------------------------------------
// Numeric normalisation on one for range [dL,dU]

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
      ROOT::Math::GSLIntegrator gsl_int(1.e-12, 1.e-6, 1000);
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
   }

   return RevArgus(m,A) / norm;
}

//-------------------------------------------------------------------------
void test_Agrus() {
//-------------------------------------------------------------------------
   auto Largus = [](const double* x,const double* p) -> double {
      return p[0]*Argus(x[0],p[1]);
   };

   TF1* fbkg = new TF1("fbkg", Largus, 0., 1., 2);
   fbkg->SetParNames("N","A");
   fbkg->SetParameters(1., 1.4);
   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad->SetGrid();
   fbkg->DrawCopy();

   // test normalization:
   printf(" norm= %.15f\n",fbkg->Integral(0.,1.));

   fbkg->SetParameters(1., 1.);
   fbkg->SetLineColor(kRed);
   fbkg->DrawCopy("SAME");
}

//-------------------------------------------------------------------------
void test_RevAgrus() {
//-------------------------------------------------------------------------
   auto Lrargus = [](const double* x,const double* p) -> double {
//       return p[0]*RevArgus(x[0],p[1]);
      return p[0]*RevArgusN(x[0],p[1]);
   };

   TF1* fbkg = new TF1("fbkg", Lrargus, 0.98, dU, 2);
   fbkg->SetParNames("N","A");
   fbkg->SetParameters(1., 0.97);
   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);

   TLegend* leg = new TLegend(0.35,0.20,0.65,0.40);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad->SetGrid();
   fbkg->DrawCopy();
   leg -> AddEntry(fbkg->Clone(),
         Form("RevArgus A=%.2f",fbkg -> GetParameter(1)),"L");

   // test normalization:
//    printf(" RevArgus norm= %.15f\n",fbkg->Integral( dL,(Mjpsi-Meta) ) );
   printf(" RevArgusN norm= %.15f\n",fbkg->Integral( dL,dU ) );

   fbkg->SetParameters(1., 5.1);
   fbkg->SetLineColor(kRed);
   fbkg->DrawCopy("SAME");
   leg -> AddEntry(fbkg->Clone(),
         Form("RevArgus A=%.2f",fbkg -> GetParameter(1)),"L");

   leg -> Draw();
}

//-------------------------------------------------------------------------
void bkg_fit(string fname, string title, string pdf="") {
//-------------------------------------------------------------------------
   TH1D* hst = plot_Mkk(fname,string("mkk_kketa"));
   hst->SetMaximum( 1.3*hst->GetMaximum() );

   // Fit KKeta background by reversed Argus function
   auto Lan = [](const double* x,const double* p) -> double {
      return bW*p[0]*RevArgusN(x[0],p[1]);
   };

   TF1* fbkg = new TF1("fbkg", Lan, 0.98, dU, 2);
   fbkg->SetParNames("Nbkg","A");
   fbkg->SetParameters(1e4, 1.);
//    fbkg->FixParameter(1, 0.5); // A

   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();
   gStyle->SetFitFormat(".3g");
   gStyle->SetStatX(0.895);
   gStyle->SetStatY(0.895);
   gStyle->SetStatH(0.13);

   gStyle->SetOptFit(111); // do not print fixed parameters!

   SetHstFace(hst);
   hst->GetYaxis()->SetTitleOffset(1.3);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
   hst->Fit("fbkg","E","",dL-1e-3,dU); // L -- Loglikelihood method

   TLegend* leg = new TLegend(0.105,0.765,0.5,0.895);
   leg->AddEntry(hst,title.c_str(),"LEP");
   leg->AddEntry(fbkg,"Argus function","L");
   leg->Draw();

//    double Nbg = fbkg->GetParameter(0);
//    double Nbg_err = fbkg->GetParError(0);
//    string Tbg( Form("N(KK#eta) = %.0f #pm  %.0f",Nbg,Nbg_err) );
//    TPaveText* pt = new TPaveText(0.25,0.65,0.55,0.73,"NDC");
//    TText* t1 = pt->AddText( Tbg.c_str() );
//    t1->SetTextColor(kBlue);
//    pt->Draw();

   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void bkg_fit2(string fname09, string fname12, string pdf="") {
//-------------------------------------------------------------------------
   TH1D* hst = plot_Mkk(fname09,string("mkk_kketa09"));
   TH1D* hst2 = plot_Mkk(fname12,string("mkk_kketa12"));
   // normalization on number of Psi(2S) events:
   double norm09 = 107.0/1.2;
   double norm12 = 341.1/3.6;
   hst->Add(hst,hst2,1.,norm12/norm09);
   hst->SetMaximum( 1.3*hst->GetMaximum() );

   // Fit KKeta background by reversed Argus function
   auto Lan = [](const double* x,const double* p) -> double {
      return bW*p[0]*RevArgusN(x[0],p[1]);
   };

   TF1* fbkg = new TF1("fbkg", Lan, 0.98, dU, 2);
   fbkg->SetParNames("Nbkg","A");
   fbkg->SetParameters(2e4, 1.);
//    fbkg->FixParameter(1, 0.5); // A=0.5

   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();
   gStyle->SetFitFormat(".3g");
   gStyle->SetStatX(0.895);
   gStyle->SetStatY(0.895);
   gStyle->SetStatH(0.13);

   gStyle->SetOptFit(111); // do not print fixed parameters!

   SetHstFace(hst);
   hst->GetYaxis()->SetTitleOffset(1.3);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
//    TFitResultPtr pfit =
   hst->Fit("fbkg","ES","",dL-1e-3,dU); // L -- Loglikelihood method

//    pfit -> Print("V"); // correlation matrix

   TLegend* leg = new TLegend(0.105,0.70,0.5,0.895);
   leg -> SetHeader("#bf{2009 & 2012}", "C");
   leg -> AddEntry(hst,"MC KK#eta","LEP");
   leg -> AddEntry(fbkg,"Argus function","L");
   leg -> Draw();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 data_fit()
//-------------------------------------------------------------------------
void data_fit(string fname, string title, int date, int Mphix,
      string pdf="") {
//-------------------------------------------------------------------------
   TH1D* hst = Subtruct(fname); // subtract side-band

   if ( Mphix == 1 ) {
      pdf += "1.pdf";
   } else {
      pdf += "0.pdf";
   }

   // Fit data
   auto Lsum = [](double* x,double* p) -> double {
      return bW*p[0]*( BreitWignerGaussN(x[0],p[1],p[2],p[3]) +
                       p[4]*RevArgusN(x[0],p[5]) );
   };

   TF1* fall = new TF1("fall", Lsum, dL, dU, 6);
   fall->SetParNames("N#phi","M#phi","G#phi","#sigma",
                     "F_{bkg}","A");
   if ( date == 2009 ) {
      fall->SetParameters(830, Mphi,Gphi,1.4e-3, 0.05, 0.97);
   } else { // 2012
      fall->SetParameters(2.5e3, Mphi,Gphi,1.1e-3, 0.06, 1.);
   }

   if ( Mphix == 1 ) {
      fall->FixParameter(1, Mphi);
   }
   fall->FixParameter(2, Gphi);
   fall->FixParameter(5, 0.97);   // A

   fall->SetLineWidth(2);
   fall->SetLineColor(kRed);
   fall->SetNpx(500);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad->SetGrid();
   gStyle->SetFitFormat(".5f");
   gStyle->SetOptFit(111); // do not print fixed parameters!

   SetHstFace(hst);
   hst->GetYaxis()->SetTitleOffset(1.3);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
   hst->Fit("fall","E","",dL,dU);

   // plot background
   double* par = fall->GetParameters();
   auto Lbkg = [par](double* x,double* p) -> double {
      return bW*par[0]*par[4]*RevArgusN(x[0],par[5]);
   };
   TF1* fbkg = new TF1("fbkg", Lbkg, dL, dU, 0);
   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);
   fbkg->SetLineStyle(kDashed);
   fbkg->Draw("SAME");

   TLegend* leg = new TLegend(0.62,0.50+Mphix*0.03,0.95,0.68+Mphix*0.03);
   leg->AddEntry(hst,title.c_str(),"LEP");
   leg->AddEntry(fall,"Signal + Bkg.","L");
   leg->AddEntry(fbkg,"Background","L");
   leg->Draw();

//    double Nsig = fall->GetParameter(0)/dx;
//    double Nsig_err = fall->GetParError(0)/dx;
//    string Tsig( Form("N(#phi#eta) = %.0f #pm  %.0f",Nsig,Nsig_err) );
//    double Nbg = fall->GetParameter(4)/dx;
//    double Nbg_err = fall->GetParError(4)/dx;
//    string Tbg( Form("N(KK#eta) = %.0f #pm  %.0f",Nbg,Nbg_err) );

//    TPaveText* pt = new TPaveText(0.62,0.21,0.95,0.33,"NDC");
//    TText* t1 = pt->AddText( Tsig.c_str() );
//    t1->SetTextColor(kRed);
//    TText* t2 = pt->AddText( Tbg.c_str() );
//    t2->SetTextColor(kBlue);
//    pt->Draw();

   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 Interference BW with Argus bkg.
//----------------------------------------------------------------------
// Interference BW with Argus bkg.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
double IntfrBWAR( double m,
                  double mphi, double gphi,               // B-W
                  double A,                               // Argus
                  double Nr, double Fb, double Ang,       // interference
                  int ret_type = 0                        // what to return
                ) {
//----------------------------------------------------------------------
// see BAM-00117:  PhysRevD.91.112001.pdf (arXiv:1504.03194v2)
//
// The return value depends on the variable 'ret_type':
//      ret_type =    0      1     2       3
//        return ->  Sum    B-W   Argus   Interference
//----------------------------------------------------------------------
   double BWnorm = 0;
   static double cacheN = 0;
   static double cacheM = 0;
   static double cacheG = 0;
   if ( cacheN > 0 &&
         mphi == cacheM && gphi == cacheG ) {
      BWnorm = cacheN;
   } else {
      auto Lbw = [mphi,gphi](const double* x,const double* p) -> double {
         return BreitWigner(x[0],mphi,gphi);
      };
      TF1* bw = new TF1("bw", Lbw, 0.98, dU, 0);
      BWnorm = 1./bw->Integral(dL,dU,1e-8);
      printf("BWnorm = %.6f\n",BWnorm);

      cacheN = BWnorm;
      cacheM = mphi;
      cacheG = gphi;
   }

   if ( m < dL ) { // phase space becomes zero
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
   double BB_k = sqrt( (1+SQ(R*p_k0)) / (1+SQ(R*p_k)) );

   double fm = r_phi*r_k*BB_phi*BB_k;
   double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k; // == m*G(m)

   // common multiplier
   double tmp =  fm / ( SQ(m2-mphi2) + SQ(GM) );

   // signal
   double BW = Nr*BWnorm*kap * fm * tmp;

   // background:
   double Ar = Nr * Fb * RevArgusN(m,A);

   // interference:
   double Intfr = 2 * sqrt(Nr*BWnorm*kap * Ar) * tmp *
                  ( (m2-mphi2)*cos(Ang) + GM*sin(Ang) );

   double Sum = Ar + BW + Intfr;

   if ( ret_type == 0 ) {
      return Sum;
   } else if ( ret_type == 1 ) {
      return BW;
   } else if ( ret_type == 2 ) {
      return Ar;
   } else if ( ret_type == 3 ) {
      return Intfr;
   }

   return 0.;
}

//----------------------------------------------------------------------
double IntfrBWARG( double m,
                   double mphi, double gphi,         // B-W
                   double sigma,                     // Gauss
                   double A,                         // Argus
                   double Nr, double Fb, double Ang, // interference
                   int ret_type = 0                  // what to return
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
   static const double Ngauss = 1./sqrt(2*M_PI);

   // integrand lambda function
   auto Lint = [=](double t) -> double{ // t == (X'-X)/sigma
      return exp(-0.5*t*t) *
      IntfrBWAR( (m+t*sigma),mphi,gphi, A, Nr,Fb,Ang, ret_type );
   };
   gsl_fun< decltype(Lint) > Fint(Lint);

   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-9, 1.e-6, 1000);
   double result = gsl_int.Integral(Fint,-5.,+5.);

   // int status = gsl_int.Status();
   if ( DEBUG ) {
      // the estimate of the absolute Error
      if ( gsl_int.Error() > 1e-6 ) {
         htt[5]->Fill(log10(gsl_int.Error()));
      } else {
         htt[5]->Fill(-7.);
      }
   }
   return Ngauss * result;
}

//-------------------------------------------------------------------------
void test_Intfr() {
//-------------------------------------------------------------------------
   double mphi = Mphi, gphi = Gphi;
   auto Lintfr = [mphi,gphi](const double* x,const double* p) -> double {
      return IntfrBWARG( x[0],
                         mphi, gphi,
                         p[0],
                         p[1],
                         p[2],p[3],p[4],
                         int(p[5]) );
   };

   TF1* fun = new TF1("fun", Lintfr, bL, dU, 6);
   fun->SetParNames("Sigma",  "A", "Nphi","Fbkg","Ang","T");
   fun->SetParameters(1.2e-3, 0.97,  1.,   0.05,   0.,  0);

   fun->SetLineWidth(2);
   fun->SetLineColor(kRed);
   fun->SetNpx(500);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad->SetGrid();
   fun->SetParameter(5, 0); // draw all
   fun->DrawCopy();
//    ymax = fun->GetMaximum();

   fun->SetParameter(5, 2); // draw Argus
   fun->SetLineColor(kBlue);
   fun->SetLineStyle(kDashed);
   fun->DrawCopy("SAME");

   fun->SetParameter(5, 1); // draw BW
   fun->SetLineColor(kGreen+2);
   fun->DrawCopy("SAME");

   fun->SetParameter(5, 3); // draw interference
   fun->SetLineColor(kMagenta);
   fun->SetLineStyle(kSolid);
   fun->DrawCopy("SAME");

   c1->Update();
}

//-------------------------------------------------------------------------
void data_Intfr_fit(string fname, string title, int date, string pdf="") {
//-------------------------------------------------------------------------
   TH1D* hst = Subtruct(fname); // subtract side-band
   if( date == 2009 ) {
      hst->SetMinimum(-10.);
   } else {
      hst->SetMinimum(-25.);
   }

   double mphi = Mphi, gphi = Gphi;
   auto Lintfr = [mphi,gphi](const double* x,const double* p) -> double {
      return bW*IntfrBWARG( x[0],
                         mphi, gphi,
                         p[0],
                         p[1],
                         p[2],p[3],p[4],
                         int(p[5]) );
   };

   TF1* fun = new TF1("fun", Lintfr, bL, dU, 6);
   fun->SetParNames("#sigma", "A", "N#phi","F_{bkg}","#vartheta","T");
   if ( date == 2009 ) {
      fun->SetParameters(1.4e-3, 0.97,  840., 0.007, 0.,  0);
   } else { // 2012
      fun->SetParameters(1.1e-3, 0.97, 2540., 0.02, 0.,  0);
   }

   fun->FixParameter(5, 0); // MUST BE FIXED!
//    gStyle->SetOptFit(112);  // print fixed parameters!
   gStyle->SetOptFit(111); // do not print fixed parameters!
//    gStyle->SetFitFormat(".3g");

   fun->FixParameter(1, 0.97);  // A
//    fun->FixParameter(0, 1.13e-3); // Sigma
//    fun->FixParameter(4, 0.); // Ang

   fun->SetLineWidth(2);
   fun->SetLineColor(kRed);
   fun->SetNpx(500);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst);
   hst->GetYaxis()->SetTitleOffset(1.3);
   hst->SetLineWidth(2);
   hst->SetLineColor(kBlack);
   hst->SetMarkerStyle(20);

   hst->Draw("EP");
   hst->Fit("fun","E","",dL,dU);
//    TFitResultPtr res = hst->Fit("fun","SE","",dL,dU);
//    res->Print();

   // draw components of Interference
   double par[20];
   fun->GetParameters(par);
   auto Fdraw = [mphi,gphi,par](const double* x,const double* p) -> double {
      return bW*IntfrBWARG( x[0],
                         mphi, gphi,
                         par[0],
                         par[1],
                         par[2],par[3],par[4],
                         int(p[0]) );
   };
   TF1* fdr = new TF1("fdr", Fdraw, bL, dU, 1);
   fdr->SetLineWidth(1);
   fdr->SetLineStyle(kDashed);
   fdr->SetNpx(500);

   TLegend* leg = new TLegend(0.62,0.42,0.98,0.67);
   leg->AddEntry(hst,title.c_str(),"LEP");
   leg->AddEntry(fun,"Result of fit","L");

   fdr->SetParameter(0, 1); // BW
   fdr->SetLineColor(kGreen+2);
   fdr->DrawCopy("SAME");
   leg->AddEntry(fdr->Clone(),"Breit-Wigner","L");

   fdr->SetParameter(0, 2); // Argus
   fdr->SetLineColor(kBlue);
   fdr->DrawCopy("SAME");
   leg->AddEntry(fdr->Clone(),"Argus","L");
//    double Nkke = fdr->Integral(L,U,1e-6)/dx;
//    double tmp = fun->GetParError(3)/par[3]; // Nb
//    double Nkke_err = Nkke * tmp;
//    string Tkke( Form("N(KK#eta)= %.1f #pm %.1f",Nkke,Nkke_err) );

   fdr->SetParameter(0, 3); // interference
   fdr->SetLineColor(kMagenta+1);
   fdr->DrawCopy("SAME");
   leg->AddEntry(fdr->Clone(),"Interference","L");

   leg->Draw();

   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void data_Intfr_scan(string fname, string title, int date) {
//-------------------------------------------------------------------------
   TH1D* hst = Subtruct(fname); // subtract side-band
   hst->SetMinimum(-30.);

   string pdf = (date==2012) ? string("mkk12_scan_full.pdf")
                : string("mkk09_scan_full.pdf");
//-------------------------------------------------------------------------
   // Fit data
   double mphi = Mphi, gphi = Gphi;
   auto Lintfr = [mphi,gphi](const double* x,const double* p) -> double {
      return bW*IntfrBWARG( x[0],
                         mphi, gphi,
                         p[0],
                         p[1],
                         p[2],p[3],p[4],
                         int(p[5]) );
   };

   TF1* fun = new TF1("fun", Lintfr, bL, dU, 6);
   fun->SetParNames("#sigma", "A",   "N#phi","F_{bkg}","#vartheta","T");
   if ( date == 2009 ) {
      fun->SetParameters(1.4e-3, 0.97,  800.,  0.03,       0,  0);
   } else { // 2012
      fun->SetParameters(1.0e-3, 0.97, 2250.,  0.04,       0,  0);
   }

   fun->FixParameter(5, 0); // MUST BE FIXED!
   fun->FixParameter(1, 0.97);  // A

   fun->SetLineWidth(2);
   fun->SetLineColor(kRed);
   fun->SetNpx(500);

   // angles in degrees
   int angmin = 90, angmax = -95, angstep = -5;

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->Print((pdf+"[").c_str()); // open pdf-file

   gStyle->SetOptFit(111); // do not print fixed parameters!
   gStyle->SetFitFormat(".5g");

   TLegend* leg = new TLegend(0.50,0.50,0.88,0.70);
   leg->AddEntry(hst,title.c_str(),"LEP");
   TLegendEntry* la = leg->AddEntry(fun,"","L"); // template

   vector<double> angl, chi2;
   for( int ang = angmin; ang != angmax; ang+=angstep ) {
      fun->FixParameter(4, M_PI*ang/180.); // Ang -> radians

      c1->cd();
      gPad->SetGrid();

      SetHstFace(hst);
      hst->GetYaxis()->SetTitleOffset(1.3);
      hst->SetLineWidth(2);
      hst->SetLineColor(kBlack);
      hst->SetMarkerStyle(20);

      hst->Draw("EP");
      hst->Fit("fun","","",dL,dU);

      la->SetLabel(Form("#vartheta= %i deg",ang));
      leg->Draw();

      c1->Update();
      c1->Print(pdf.c_str()); // add to pdf-file

      double ch2 = fun->GetChisquare();
      ch2 /= fun->GetNDF();
      angl.push_back(ang);
      chi2.push_back(ch2);
   }

   c1->cd();
   gPad->SetGrid();

   int nch = chi2.size();
   auto gChi2 = new TGraph( nch,angl.data(),chi2.data() );
   gChi2->SetTitle(";#vartheta, degrees;#chi^{2} / ndf");
   gChi2->SetMarkerColor(kBlue);
   gChi2->SetMarkerStyle(21);
   gChi2->GetYaxis()->SetTitleOffset(1.4);
   gChi2->SetLineWidth(2);

   gChi2->Draw("APL");

   double width = 0.025*title.size();
   TLegend* leg2 = new TLegend(0.89-width,0.79,0.89,0.89);
   leg2 -> SetHeader(title.c_str(), "C");
   leg2 -> Draw();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
}

// {{{1 Combined fit
//-------------------------------------------------------------------------
// Combined fit
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
class myChi2 {
   public:
      myChi2(TH1D* hst1, TH1D* hst2) {
         int ntot = 200;
         mkk.reserve(ntot);
         data.reserve(ntot);
         er_data.reserve(ntot);

         for(int n = 1; n <= hst1->GetNbinsX(); ++n) {
            double m = hst1->GetBinCenter(n);
            if( m < dL ) {
               continue;
            }
            mkk.push_back(m);
            data.push_back(hst1->GetBinContent(n));
            er_data.push_back(SQ(hst1->GetBinError(n)));
         }
         nbins1 = mkk.size();

         for(int n = 1; n <= hst2->GetNbinsX(); ++n) {
            double m = hst2->GetBinCenter(n);
            if( m < dL ) {
               continue;
            }
            mkk.push_back(m);
            data.push_back(hst2->GetBinContent(n));
            er_data.push_back(SQ(hst2->GetBinError(n)));
         }
      }

      // Size function
      unsigned int Size() const {
         return mkk.size();
      }

      // chi^2 - function
      // the signature of this operator() MUST be exactly this:
      double operator() (const double* p) {
         static const double maxResValue = DBL_MAX / 1000;

         // A,Fb,Ang - common parameters
         double mphi = Mphi, gphi = Gphi;
         auto F1 = [mphi,gphi,p](double m) -> double {
            return bW*IntfrBWARG( m, mphi,gphi,
                   //sigma, A,   Nr,   Fb,   Ang, ret
                     p[0], p[2], p[3], p[5], p[6], 0 );
         };
         auto F2 = [mphi,gphi,p](double m) -> double {
            return bW*IntfrBWARG( m, mphi,gphi,
                   //sigma, A,   Nr,   Fb,   Ang, ret
                     p[1], p[2], p[4], p[5], p[6], 0 );
         };

         double chi2 = 0;
         unsigned int n = Size();
         for(unsigned int i = 0; i < n; i++) {
//          for(unsigned int i = 0; i < nbins1; i++) { // debug!
            double m = mkk[i];
            double f = ( i < nbins1 ) ? F1(m) : F2(m);

            double dif = SQ(data[i] - f);
            double er2 = er_data[i];

            double resval = dif / er2;
            // avoid inifinity or nan in chi2
            if( resval < maxResValue ) {
               chi2 += resval;
            } else {
               chi2 += maxResValue;
            }
         } // end of for()

         return chi2;
      }

   private:
      int nbins1; // number of bins in first set of data
      vector<double> mkk;
      vector<double> data, er_data;
};
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void combine_Intfr(string pdf="") {
//-------------------------------------------------------------------------
   TH1D* hst1 = Subtruct("data_12psip_all.root"); // subtract side-band
   TH1D* hst2 = Subtruct("data_09psip_all.root"); // subtract side-band
   hst1->SetMinimum(-25.);
   hst2->SetMinimum(-10.);

   myChi2 chi2_fit(hst1,hst2);

   cout << " chi2_fit: Size= " << chi2_fit.Size() << endl;
//    for ( auto m : chi2_fit.mkk ) {
//       cout << m << endl;
//    }

   // ========================= Fit with ROOT =========================
   // == fit configuration
   ROOT::Fit::Fitter fitter;
   // set parameters of fitter: (Minuit,Minuit2,Fumili,GSLMultiFit...)
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");
   fitter.Config().MinimizerOptions().SetPrintLevel(3);
   // print the default minimizer option values
   ROOT::Math::MinimizerOptions min_opt;
   min_opt.Print();

   // set names and start values for parameters
   vector<string> par_name
        {"sig12", "sig09",  "A", "N12", "N09", "Fbkg", "Ang" };
   vector<double> par_ini
        { 1.1e-3, 1.4e-3,  0.97, 2550.,  820.,   0.01,   0. };

   const unsigned int Npar = par_name.size(); // number of parameters
   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for(unsigned int i = 0; i < Npar; i++) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

//    fitter.Config().ParSettings(1).Fix(); // sig09 debug!
//    fitter.Config().ParSettings(2).Fix(); // sig12 debug!

   // fix  parameters
   fitter.Config().ParSettings(2).SetValue(0.97);              // A
   fitter.Config().ParSettings(2).Fix();
//    fitter.Config().ParSettings(2).SetStepSize(0.001);
//    fitter.Config().ParSettings(2).SetLimits(1.e-3,10.);
//    fitter.Config().ParSettings(6).SetValue(0.);                 // Ang
//    fitter.Config().ParSettings(6).Fix();


   // == Fit
   fitter.FitFCN(Npar, chi2_fit, nullptr, chi2_fit.Size(), true);

   // to obtain reliable errors:
//    fitter.CalculateHessErrors();
   fitter.CalculateMinosErrors();

   // == Fit result
   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double chi2   = res.Chi2();
   int ndf       = res.Ndf();
   string Tchi2( Form("#chi^{2}/ndf = %.1f / %i",chi2,ndf) );
   const vector<double> par   = res.Parameters();
   const vector<double> er_par = res.Errors();
   string Ts12( Form("#sigma(2012) = %.2f #pm  %.2f MeV",
                par[0]*1e3, er_par[0]*1e3) );
   string Ts09( Form("#sigma(2009) = %.2f #pm  %.2f MeV",
                par[1]*1e3, er_par[1]*1e3) );
   string Ta( Form("A = %.2f ",par[2]) );
   if ( er_par[2] == 0 ) {
      Ta += string( "(fixed)" );
   } else {
      Ta += string( Form("#pm  %.2f",er_par[2]) );
   }
   string Tn12( Form("N_{#phi}(2012) = %.1f #pm  %.1f",par[3],er_par[3]) );
   string Tn09( Form("N_{#phi}(2009) = %.1f #pm  %.1f",par[4],er_par[4]) );
   string Tbkg( Form("F_{bkg} = %.3f #pm  %.3f",par[5],er_par[5]) );
   string Tang( Form("#vartheta = %.3f #pm  %.3f",par[6],er_par[6]) );

   // ========================= Draw results =========================
//    TCanvas* c1 = new TCanvas("c1","...",0,0,900,900); // slide
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000); // memo
   c1->Divide(1,2);

   // draw components of Interference
   double mphi = Mphi, gphi = Gphi;
   auto F1 = [mphi,gphi,par](const double* x,const double* p) -> double {
      return bW*IntfrBWARG( x[0], mphi,gphi,
             //sigma,  A,     Nr,    Fb,   Ang,      ret
               par[0],par[2],par[3],par[5],par[6], int(p[0]) );
   };
   auto F2 = [mphi,gphi,par](const double* x,const double* p) -> double {
      return bW*IntfrBWARG( x[0], mphi,gphi,
             //sigma,  A,     Nr,    Fb,   Ang,      ret
               par[1],par[2],par[4],par[5],par[6], int(p[0]) );
   };

   TF1* fdr1 = new TF1("fdr1", F1, bL, dU, 1);
   fdr1->SetNpx(500);
   TF1* fdr2 = new TF1("fdr2", F2, bL, dU, 1);
   fdr2->SetNpx(500);

//    TLegend* leg = new TLegend(0.51,0.40,0.89,0.89); // slide
   TLegend* leg = new TLegend(0.51,0.45,0.89,0.89); // memo
   leg->SetHeader("#bf{2012(top)   2009(bottom)}","C");
   leg->AddEntry(hst1,"Data","LEP");

   c1->cd(1);
   gPad->SetGrid();
   SetHstFace(hst1);
   hst1->SetLineWidth(2);
   hst1->SetLineColor(kBlack);
   hst1->SetMarkerStyle(20);

   hst1->Draw("EP");
   fdr1->SetParameter(0, 0); // Sum
   fdr1->SetLineWidth(2);
   fdr1->SetLineColor(kRed+1);
   fdr1->DrawCopy("SAME");
   leg->AddEntry(fdr1->Clone(),"Combined fit","L");

   fdr1->SetParameter(0, 1); // BW
   fdr1->SetLineWidth(1);
   fdr1->SetLineStyle(kDashed);
   fdr1->SetLineColor(kGreen+2);
   fdr1->DrawCopy("SAME");
   leg->AddEntry(fdr1->Clone(),"Breit-Wigner","L");

   fdr1->SetParameter(0, 2); // Argus
   fdr1->SetLineColor(kBlue);
   fdr1->DrawCopy("SAME");
   leg->AddEntry(fdr1->Clone(),"Argus","L");

   fdr1->SetParameter(0, 3); // interference
   fdr1->SetLineColor(kMagenta+1);
   fdr1->DrawCopy("SAME");
   leg->AddEntry(fdr1->Clone(),"Interference","L");
   leg->Draw();

   c1->cd(2);
   gPad->SetGrid();
   SetHstFace(hst2);
   hst2->SetLineWidth(2);
   hst2->SetLineColor(kBlack);
   hst2->SetMarkerStyle(20);

   hst2->Draw("EP");
   fdr2->SetParameter(0, 0); // Sum
   fdr2->SetLineWidth(2);
   fdr2->SetLineColor(kRed+1);
   fdr2->DrawCopy("SAME");

   fdr2->SetParameter(0, 1); // BW
   fdr2->SetLineWidth(1);
   fdr2->SetLineStyle(kDashed);
   fdr2->SetLineColor(kGreen+2);
   fdr2->DrawCopy("SAME");

   fdr2->SetParameter(0, 2); // Argus
   fdr2->SetLineColor(kBlue);
   fdr2->DrawCopy("SAME");

   fdr2->SetParameter(0, 3); // interference
   fdr2->SetLineColor(kMagenta+1);
   fdr2->DrawCopy("SAME");

//    TPaveText* pt = new TPaveText(0.51,0.26,0.89,0.999,"NDC"); // slide
   TPaveText* pt = new TPaveText(0.51,0.40,0.89,0.99,"NDC"); // memo
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
   pt -> AddText( Tchi2.c_str() );
//    TLine* l0 = pt -> AddLine(0.,0.87,1.,0.87 );
//    l0 -> SetLineWidth(2);
   Ta = "#lower[0.4]{" + Ta + "}";
   pt -> AddText( Ta.c_str() );
   Tbkg = "#lower[0.3]{" + Tbkg + "}";
   pt -> AddText( Tbkg.c_str() );
   Tang = "#lower[0.3]{" + Tang + "}";
   pt -> AddText( Tang.c_str() );
   TLine* l1 = pt -> AddLine(0.,0.52,1.,0.52 );
   l1 -> SetLineWidth(2);
   Tn12 = "#lower[0.7]{" + Tn12 + "}";
   pt -> AddText( Tn12.c_str() );
   Ts12 = "#lower[0.7]{" + Ts12 + "}";
   pt -> AddText( Ts12.c_str() );
   TLine* l2 = pt -> AddLine(0.,0.26,1.,0.26 );
   l2 -> SetLineWidth(2);
   Tn09 = "#lower[1.1]{" + Tn09 + "}";
   pt -> AddText( Tn09.c_str() );
   Ts09 = "#lower[1.1]{" + Ts09 + "}";
   pt -> AddText( Ts09.c_str() );
   pt -> AddText( "" );
   pt -> Draw();

   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 MAIN:
//-------------------------------------------------------------------------
void mass_kk_fit() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);
   gStyle->SetLegendFont(42);
//    gStyle->SetLegendTextSize(0.04);

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
      "MC #phi#eta 2009",
      "MC #phi#eta 2012",
      "MC KK#eta 2009",
      "MC KK#eta 2012"
   };
   // Numbers of items for date:
   int id = 0, ii = 3, is = 5, ib = 7; // 2009
//    int id = 1, ii = 5, is = 6, ib = 8; // 2012

// ------------- signal ---------------
//    test_BreitWigner();
//    string pdf = string("mkk")+((id==0) ? "09" : "12")+"_sig_log.pdf";
//    sig_fit(fnames.at(is),titles.at(is),pdf); // fit phi-eta by BW x Gauss

//    sig_fit_mc("prod-9/mcsig_kkmc_12.root","mcmkk12_sig.pdf");
//    sig_fit_mc("prod-9/mcsig_kkmc_09.root","mcmkk09_sig.pdf");

// ------------- background ---------------
//    test_Agrus();
//    test_RevAgrus();

//    string pdf = string("mkk")+((id==0) ? "09" : "12")+"_fitbg.pdf";
//    bkg_fit(fnames.at(ib),titles.at(ib),pdf); // fit kketa by Argus

   // sum 09+12
//    string pdf = string("mkk_sum_fitbg.pdf");
//    bkg_fit2(fnames[7],fnames[8],pdf); // fit kketa by Argus

// ------------- data: no interference ---------------
//    string pdf = string("mkk")+((id==0) ? "09" : "12")+"_fit";
//    int Mphix = 1; // 1 - mphi is fixed at the PDG
//    data_fit(fnames.at(id),titles.at(id),2009+id*3,Mphix,pdf);

// ------------- data: interference ---------------
//    test_Intfr();
//    string pdf = string("mkk")+((id==0) ? "09" : "12")+"_intfr.pdf";
//    data_Intfr_fit(fnames.at(id),titles.at(id),2009+id*3,pdf);

//    data_Intfr_scan(fnames.at(id),titles.at(id),2009+id*3);

// ------------- data: common fit ---------------
   combine_Intfr("mkk_cf_test.pdf");

}
