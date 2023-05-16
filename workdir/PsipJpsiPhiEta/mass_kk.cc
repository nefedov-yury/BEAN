// plot M(K+K-) for data, inclusive MC, signal MC and MC KKeta
// cuts (see cuts.h): Mrec + chi^2(4C) + Mgg
// -> mass_kk_YEAR.pdf

#include "masses.h"

// {{{1 helper functions
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
TH1D* get_Mkk(string fname, string hname, int type=0) {
//--------------------------------------------------------------------
#include "cuts.h"

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

   // const double dU = 1.101;
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
double Argus(double x, double a) {
//--------------------------------------------------------------------
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
void test_RevAgrus() {
//--------------------------------------------------------------------
   static const double dL = 2*Mk; // the left cutoff = 0.987354
   static const double dU = 1.08; // MUST BE < 1.0835 !!!

   auto Lrargus = [](const double* x,const double* p) -> double {
      return p[0]*RevArgusN(x[0],p[1],p[2]);
   };

   TF1* fbkg = new TF1("RevArgus", Lrargus, 0.98, dU, 3);
   fbkg -> SetParNames("N","A","Sl");
   fbkg -> SetParameters(1., 9., -1.8);
   fbkg -> SetLineWidth(2);
   fbkg -> SetLineColor(kBlue);

   TLegend* leg = new TLegend(0.49,0.69,0.89,0.89);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1 -> cd();
   gPad -> SetGrid();
   fbkg -> DrawCopy();
   leg -> AddEntry(fbkg -> Clone(), Form("RevArgus A=%.2f Sl=%.2f",
            fbkg -> GetParameter(1),fbkg -> GetParameter(2)),"L");

   // test normalization:
   printf(" RevArgusN norm= %.15f\n",fbkg -> Integral( dL,dU ) );

   fbkg -> SetParameter(1, 6.0 );  // A
   fbkg -> SetLineColor(kRed);
   fbkg -> SetLineStyle(kDashed);
   fbkg -> DrawCopy("SAME");
   leg -> AddEntry(fbkg -> Clone(), Form("RevArgus A=%.2f Sl=%.2f",
            fbkg -> GetParameter(1),fbkg -> GetParameter(2)),"L");
   printf(" RevArgusN norm2= %.15f\n",fbkg -> Integral( dL,dU ) );

   fbkg -> SetParameter(1, 0.0 );  // A
   fbkg -> SetLineColor(kGreen+4);
   fbkg -> DrawCopy("SAME");
   leg -> AddEntry(fbkg -> Clone(), Form("RevArgus A=%.2f Sl=%.2f",
            fbkg -> GetParameter(1),fbkg -> GetParameter(2)),"L");

   leg -> Draw();
}

// {{{1 Fit Side-Band
//-------------------------------------------------------------------------
TF1* Fit_Sb( TH1D* hst, int type = 1 ) {
//-------------------------------------------------------------------------
   static const double dL = 2*Mk; // the left cutoff = 0.987354
   static const double dU = 1.08; // MUST BE < 1.0835 !!!
   TF1* fret = nullptr;
   if ( type == 1 ) { // fit side-band by constant
      fret = (TF1*)gROOT -> GetFunction("pol0")->Clone();
   }

   if ( type == 2 ) { // Reverse Argus
      double bW = hst -> GetBinWidth(1); // bin width
      double norm = hst -> Integral();
      auto Lan = [bW](const double* x,const double* p) -> double {
         return bW*p[0]*RevArgusN(x[0],p[1],p[2]);
      };
      fret = new TF1("fret", Lan, 0.98, dU, 3);
      fret -> SetParNames("N","a","Sl");
      fret -> SetParameters(norm, 5., -1.8);
      fret -> FixParameter(2, -1.8); // slope
   }

   if ( fret ) {
      fret -> SetLineColor(kRed);
      fret -> SetLineWidth(2);
      fret -> SetLineStyle(kDashed);
      hst -> Fit(fret,"LE","E",dL,dU);
   }
   return fret;
}

// {{{1 plot_mass_kk()
//--------------------------------------------------------------------
void plot_mass_kk(int date, bool PR=false) { // defalt is memo
//--------------------------------------------------------------------
// PR = true -> for presentation
#include "cuts.h"

   vector<string> fnames = {
      "data_09psip_all.root",
      "data_12psip_all.root",
      "data_21psip_all.root",
      "data_3650_all.root",
      "mcinc_09psip_all.root",
      "mcinc_12psip_all.root",
      "mcinc_21psip_all.root",
      "mcsig_kkmc_09.root",
      "mcsig_kkmc_12.root",
      "mcsig_kkmc_21.root",
   };
      // "mckketa_kkmc_09.root",
      // "mckketa_kkmc_12.root"

   vector<string> titles = {
      "Data 2009",
      "Data 2012",
      "Data 2021",
      "E=3.65 GeV",
      "MC inclusive 2009",
      "MC inclusive 2012",
      "MC inclusive 2021",
      "MC signal #phi#eta 2009",
      "MC signal #phi#eta 2012",
      "MC signal #phi#eta 2021",
   };
      // "MC non-#phi KK#eta 2009",
      // "MC non-#phi KK#eta 2012"

   int id = 0, ii = 4, is = 7; // 2009
   if ( date == 2012 ) {
       id = 1, ii = 5, is = 8; // 2012
   }
   if ( date == 2021 ) {
       id = 2, ii = 6, is = 9; // 2021
   }

   TH1D* hst[30];
   // vector<int> tmp {id, is, ib};
   vector<int> tmp {id, is};
   for( auto i : tmp )  {
      string hname = string("mkk_") + to_string(i+1);
      hst[i] = get_Mkk(fnames.at(i),hname,0);
      hst[10+i] = get_Mkk(fnames.at(i),hname+string("_sb"),1);
   }

   TCanvas* c1 = nullptr;
   if ( PR ) {
      c1 = new TCanvas("c1","...",0,0,1200,500); // Pr
      c1->Divide(2,1);
   } else {
      c1 = new TCanvas("c1","...",0,0,500,900); // memo
      c1->Divide(1,2);
   }

   for ( int it : {0,10} ) {
      hst[id+it]->SetLineWidth(2);
      hst[id+it]->SetLineColor(kBlack); // data
      hst[id+it]->SetMarkerStyle(20);
      hst[id+it]->SetMarkerSize(0.7);

      hst[is+it]->SetLineWidth(2);
      hst[is+it]->SetLineColor(kGreen+2); // MC phi eta

      // hst[ib+it]->SetLineColor(kBlue+1);
      // hst[ib+it]->SetLineWidth(2);
   }

   // normalize MC on data
   double sig_scale =
      ( hst[id]->Integral() - hst[id+10]->Integral() ) /
      ( hst[is]->Integral() - hst[is+10]->Integral() );
   hst[is]->Scale( sig_scale );

   /*
   const double wphi = 5*Gphi;

   int n1 = hst[id]->FindBin(Mphi-wphi);
   int n2 = hst[id]->FindBin(Mphi+wphi);
   int n3 = hst[id]->GetNbinsX();
   double sig_scale = 1.0;
   double bg_scale = 1.0;
   for (int iter = 0; iter < 3; iter++) {
      double sig = 1.0;
      if ( iter==0 ) {
         sig = hst[id]->Integral(n1,n2) /
            hst[is]->Integral(n1,n2);
      } else {
         sig = (hst[id]->Integral(n1,n2) - hst[ib]->Integral(n1,n2))
            / hst[is]->Integral(n1,n2);
      }
      hst[is]->Scale( sig );
      sig_scale *= sig;

      double bg = (hst[id]->Integral(n2,n3) - hst[is]->Integral(n2,n3))
         / hst[ib]->Integral(n2,n3);
      hst[ib]->Scale( bg );
      bg_scale *= bg;
      cout << " iter# " << iter << " sig_scale= " << sig_scale
           << " bg_scale= " << bg_scale << endl;
   }
   */

   // normalize side-band on the same numbers
   hst[10+is]->Scale( sig_scale );
   // hst[10+ib]->Scale( bg_scale );

   double ymax = max( hst[id]->GetMaximum(), hst[is]->GetMaximum() );
   hst[id]->SetMaximum(1.2*ymax);

   // hst[21] = (TH1D*)hst[is]->Clone("mc_sig_clone");
   // hst[21]->Add(hst[ib]); // MC(sig) + MC(bg)
   // hst[21]->SetLineColor(kRed+1);

   // hst[22] = (TH1D*)hst[10+is]->Clone("mc_sig_sb_clone");
   // hst[22]->Add(hst[10+ib]); // MC(sig) + MC(bg)
   // hst[22]->SetLineColor(kRed+1);

   c1->cd(1);
   gPad->SetGrid();

   SetHstFace(hst[id]);
   hst[id]->GetXaxis()->SetTitleOffset(1.1);
   hst[id]->GetYaxis()->SetTitleOffset(1.3);

   hst[id]->Draw("E");
   hst[is]->Draw("HIST SAME");
   // hst[ib]->Draw("SAME HIST");
   // hst[21]->Draw("SAME HIST");
   hst[id]->Draw("E,SAME");

   TLegend* leg = nullptr;
   if ( PR ) {
      leg = new TLegend(0.48,0.55,0.89,0.89);
   } else {
      leg = new TLegend(0.49,0.53,0.89,0.89); // memo
      leg->SetTextSize(0.05);
   }
   leg->SetHeader(
         "|M_{#lower[-0.2]{#gamma#gamma}}#minus M_{#eta}| < 0.024",
         "C");
   leg->AddEntry(hst[id],titles[id].c_str(),"PLE");
   leg->AddEntry(hst[is],titles[is].c_str(),"L");
   // leg->AddEntry(hst[ib],titles[ib].c_str(),"L");
   // leg->AddEntry(hst[21],"Sum of MC","L");
   leg->Draw();

   gPad->RedrawAxis();

   c1->cd(2);
   gPad->SetGrid();

   if ( date == 2009 ) {
      hst[10+id]->SetMaximum(4);
   } else if ( date == 2012 ) {
      hst[10+id]->SetMaximum(6);
   } else if ( date == 2021 ) {
      hst[10+id]->SetMaximum(14);
   }

   SetHstFace(hst[10+id]);
   hst[10+id]->GetXaxis()->SetTitleOffset(1.1);

   int type = 2;
   TF1* fit = Fit_Sb(hst[10+id],type);
   if ( !fit ) {
      hst[10+id]->Draw("E");
   }

   hst[10+is]->Draw("HIST SAME");
   cout << " Integral signal= " << hst[10+is]->Integral() << endl;
   // hst[10+ib]->Draw("SAME HIST");
   // hst[10+id]->Draw("E,SAME");

   TPaveText* pt = new TPaveText(0.45,0.76,0.89,0.89,"NDC,TR");
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

   c1->Update();

   string pdf( Form("mass_kk_%02i%s.pdf",date,(PR) ? "_PR" : "" ) );
   c1->Print(pdf.c_str());
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
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.755);
   gStyle->SetStatW(0.245);


   // just test
   // test_RevAgrus();

   // presentation
   // plot_mass_kk(2009,true);
   // plot_mass_kk(2012,true);
   plot_mass_kk(2021,true);

   // memo
   // plot_mass_kk(2009);
   // plot_mass_kk(2012);
}
