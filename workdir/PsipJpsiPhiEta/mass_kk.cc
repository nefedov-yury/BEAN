// mass_kk.cc
// plot M(K+K-) for data and signal MC for Mgg in central and side
// bands.
//      -> mass_kk_{YEAR}.pdf
//      -> mass_kksb_{YEAR}.pdf

#include "masses.h"

// {{{1 helper functions
// GLOBAL: name of folder with root files
string Dir;

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

// {{{1 data processing
//--------------------------------------------------------------------
TH1D* get_Mkk(string fname, string hname, int type=0)
//--------------------------------------------------------------------
{
#include "cuts.h"

   fname = Dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   const double dU = 1.08;
   const double bL = 0.98; // first bin < dL
   const int Nbins = 100;
   const double bW = (dU-bL)/Nbins; // bin width
   string title(";M^{ inv}_{ K^{#plus}K^{#minus }}, GeV/c^{2}");
   title += string(";Entries/")+string(Form("%.0fMeV/c^{2}",bW*1e3));
   TH1D* mkk = new TH1D(hname.c_str(), title.c_str(), Nbins,bL,dU);

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

// {{{1 Argus functions
//--------------------------------------------------------------------
double Argus(double x, double a)
//--------------------------------------------------------------------
{
   // ARGUS distribution:
   // https://en.wikipedia.org/wiki/ARGUS_distribution
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

   static const double dL = 2*Mk; // the left cutoff = 0.987354
   static const double dU = 1.08; // MUST BE < 1.0835 !!!

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
      cacheN = norm;
      cacheA = A;
      cacheS = slope;
   }

   return RevArgus(m,A) * (1+slope*(m-1.02)) / norm;
}

//--------------------------------------------------------------------
void test_RevAgrus()
//--------------------------------------------------------------------
{
   static const double dL = 2*Mk; // the left cutoff = 0.987354
   static const double dU = 1.08; // MUST BE < 1.0835 !!!

   auto Lrargus = [](const double* x,const double* p) -> double {
      return p[0]*RevArgusN(x[0],p[1],p[2]);
   };

   TF1* fbkg = new TF1("RevArgus", Lrargus, 0.98, dU, 3);
   fbkg->SetParNames("N","A","Sl");
   fbkg->SetParameters(1., 9., -1.8);
   fbkg->SetLineWidth(2);
   fbkg->SetLineColor(kBlue);

   TLegend* leg = new TLegend(0.49,0.69,0.89,0.89);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad->SetGrid();
   fbkg->DrawCopy();
   leg->AddEntry(fbkg->Clone(), Form("RevArgus A=%.2f Sl=%.2f",
            fbkg->GetParameter(1),fbkg->GetParameter(2)),"L");

   // test normalization:
   printf(" RevArgusN norm= %.15f\n",fbkg->Integral( dL,dU ) );

   fbkg->SetParameter(1, 6.0 );  // A
   fbkg->SetLineColor(kRed);
   fbkg->SetLineStyle(kDashed);
   fbkg->DrawCopy("SAME");
   leg->AddEntry(fbkg->Clone(), Form("RevArgus A=%.2f Sl=%.2f",
            fbkg->GetParameter(1),fbkg->GetParameter(2)),"L");
   printf(" RevArgusN norm2= %.15f\n",fbkg->Integral( dL,dU ) );

   fbkg->SetParameter(1, 0.0 );  // A
   fbkg->SetLineColor(kGreen+4);
   fbkg->DrawCopy("SAME");
   leg->AddEntry(fbkg->Clone(), Form("RevArgus A=%.2f Sl=%.2f",
            fbkg->GetParameter(1),fbkg->GetParameter(2)),"L");

   leg->Draw();
}

// {{{1 Fit Side-Band
//-------------------------------------------------------------------------
TF1* Fit_Sb( TH1D* hst, int type=1, TH1D* mcsig=nullptr )
//-------------------------------------------------------------------------
{
   static const double dL = 2*Mk; // the left cutoff = 0.987354
   static const double dU = 1.08; // MUST BE < 1.0835 !!!
   TF1* fret = nullptr;
   double bW = hst->GetBinWidth(1); // bin width
   double norm = hst->Integral();
   if ( type == 1 ) { // fit side-band by constant
      if ( mcsig ) {
         // add MC signal
         double W = bW/(dU-2*Mk); // normalization on 1
         auto Fs = [W,mcsig](const double* x,const double* p) {
            int bin = mcsig->GetXaxis()->FindBin(x[0]);
            return W*p[0] + p[1]*mcsig->GetBinContent(bin);
         };
         fret = new TF1("fret", Fs, 0.98, dU, 2);
         fret->SetParNames("Nbg","norm");
         fret->SetParameters(norm, 1.);
         fret->FixParameter(1, 1.); // norm
      } else {
         fret = (TF1*)gROOT->GetFunction("pol0")->Clone();
      }
   }

   if ( type == 2 ) { // Reverse Argus
      if ( !mcsig ) {
         auto Lan = [bW](const double* x,const double* p) {
            return bW*p[0]*RevArgusN(x[0],p[1],p[2]);
         };
         fret = new TF1("fret", Lan, 0.98, dU, 3);
         fret->SetParNames("Nbg","a","Sl");
         fret->SetParameters(norm, 5., -1.8);
         fret->FixParameter(2, -1.8); // slope
      } else {
         auto Lan = [bW,mcsig](const double* x,const double* p) {
            int bin = mcsig->GetXaxis()->FindBin(x[0]);
            return bW*p[0]*RevArgusN(x[0],p[1],p[2]) +
               p[3]*mcsig->GetBinContent(bin);
         };
         fret = new TF1("fret", Lan, 0.98, dU, 4);
         fret->SetParNames("Nbg","a","Sl","norm");
         fret->SetParameters(norm, 5., -1.8,1.);
         fret->FixParameter(2, -1.8); // slope
         fret->FixParameter(3, 1.); // norm
      }
   }

   if ( fret ) {
      fret->SetLineColor(kRed);
      fret->SetLineWidth(2);
      // fret->SetLineStyle(kDashed);
      hst->Fit(fret,"LE","E",dL,dU);
   }
   return fret;
}

// {{{1 plot_mass_kk()
//--------------------------------------------------------------------
void plot_mass_kk(int date,  int Cx=500, int Cy=900)
//--------------------------------------------------------------------
{
#include "cuts.h"

   string fname[3], hn[3];
   // data
   fname[0] = string( Form("data_%02ipsip_all.root",date%100) );
   hn[0] = string( Form("mkk_data_%i",date) );
   // "MC inclusive"
   fname[1] = string( Form("mcinc_%02ipsip_all.root",date%100) );
   hn[1] = string( Form("mkk_mcinc_%i",date) );
   // "MC signal #phi#eta"
   fname[2] = string( Form("mcsig_kkmc_%02i.root",date%100) );
   hn[2] = string( Form("mkk_mcsig_%i",date) );

   TH1D* dat[2];
   TH1D* mcs[2];
   for ( int isb = 0; isb < 2; ++isb  )  {
      string ssb = (isb==1) ? "_sb" : "";
      dat[isb] = get_Mkk(fname[0],hn[0]+ssb,isb);
      mcs[isb] = get_Mkk(fname[2],hn[2]+ssb,isb);
   }

   for ( auto& d : dat ) { // data
      SetHstFace(d);
      d->SetLineWidth(2);
      d->SetLineColor(kBlack);
      d->SetMarkerStyle(20);
      d->SetMarkerSize(0.7);
   }

   // normalize MC on data
   double sig_scale =
      ( dat[0]->Integral() - dat[1]->Integral() ) /
      ( mcs[0]->Integral() - mcs[1]->Integral() );
   for ( auto& mc : mcs ) { // MC phi eta
      mc->Scale( sig_scale );
      mc->SetLineWidth(2);
      mc->SetLineColor(kGreen+2);
   }

   auto name = Form("c1_mkk_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   double ymax = max( dat[0]->GetMaximum(), mcs[0]->GetMaximum() );
   dat[0]->SetMaximum(1.2*ymax);

   dat[0]->GetXaxis()->SetTitleOffset(1.1);
   dat[0]->GetYaxis()->SetTitleOffset(1.1);
   dat[0]->GetYaxis()->SetMaxDigits(3);

   dat[0]->Draw("E");
   mcs[0]->Draw("HIST SAME");
   dat[0]->Draw("E,SAME");

   TLegend* leg = new TLegend(0.50,0.65,0.892,0.89);
   leg->SetTextSize(0.05);
   leg->SetHeader(
         "|M_{#lower[-0.2]{#gamma#gamma}}#minus M_{#eta}| < 0.024",
         "C");
   leg->AddEntry(dat[0],Form("Data %i",date),"PLE");
   leg->AddEntry(mcs[0],"MC signal #phi#eta","L");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf( Form("mass_kk_%02i.pdf",date%100) );
   c1->Print(pdf.c_str());

   auto name2 = Form("c1_mkk_sb_%i",date);
   TCanvas* c2 = new TCanvas(name2,name2,0,Cy/2,Cx,Cy);
   c2->cd();
   gPad->SetGrid();

   // if ( date == 2009 ) {
      // hst[10]->SetMaximum(4);
   // } else if ( date == 2012 ) {
      // hst[10]->SetMaximum(6);
   // } else if ( date == 2021 ) {
      // hst[10]->SetMaximum(14);
   // }

   dat[1]->GetXaxis()->SetTitleOffset(1.1);
   dat[1]->GetYaxis()->SetTitleOffset(1.1);

   TF1* fit = nullptr;
   int type = 1; // line
   // int type = 2; // rev-argus
   // fit = Fit_Sb(dat[1],type,mcs[1]); // take into account MC sig.
   if ( !fit ) {
      dat[1]->Draw("E");
   }

   mcs[1]->Draw("HIST SAME");
   auto Sig  = mcs[1]->Integral();
   auto Sdat = dat[1]->Integral();
   printf("%i: side-band Signal/Data= %.1f+/-%.1f%%\n",
         date,Sig/Sdat*100,sqrt(Sig)/Sdat*100);

   TPaveText* pt = new TPaveText(0.45,0.76,0.892,0.89,"NDC,TR");
   pt->SetTextAlign(22);
   pt->SetTextFont(42);
   pt->SetTextSize(0.05);
   pt->AddText("Side-Band");
   pt->AddText(
         Form( "%.3f<"
            "|M_{#lower[-0.2]{#gamma#gamma}}#minus M_{#eta}|"
            "<%.3f", shift_eta,shift_eta+weta)
         );
   pt->Draw();

   gPad->RedrawAxis();

   c2->Update();
   string pdfsb( Form("mass_kksb_%02i.pdf",date%100 ) );
   c2->Print(pdfsb.c_str());
}

// {{{1 MAIN:
//--------------------------------------------------------------------
void mass_kk() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(42);
   // gStyle->SetOptFit(0);
   gStyle->SetOptFit(111); // do not print fixed params
   gStyle->SetFitFormat(".1f");
   gStyle->SetStatX(0.892);
   gStyle->SetStatY(0.755);
   gStyle->SetStatW(0.245);


   //========================================================
   // set the name of the folder with the root files
   Dir = "prod_v709n4/";
   //========================================================

   size_t Cx = 800, Cy = 640; // canvas sizes, X/Y = 1.25

   // just test
   // test_RevAgrus();

   // fig.13
   // for ( int date : {2009} ) {
   for ( int date : {2009, 2012, 2021} ) {
      plot_mass_kk(date, Cx, Cy);
   }
}
