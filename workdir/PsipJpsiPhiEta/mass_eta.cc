// plot M(gamma gamma) distributions
// after 4C kinematic fit and chi^2 cut
// -> mass_eta.pdf

#include "masses.h"

//--------------------------------------------------------------------
void SetHstFace(TH1* hst) {
//--------------------------------------------------------------------
   TAxis* X = hst -> GetXaxis();
   if ( X ) {
      X -> SetLabelFont(62);
      X -> SetLabelSize(0.04);
      X -> SetTitleFont(62);
      X -> SetTitleSize(0.04);
   }
   TAxis* Y = hst -> GetYaxis();
   if ( Y ) {
      Y -> SetLabelFont(62);
      Y -> SetLabelSize(0.04);
      Y -> SetTitleFont(62);
      Y -> SetTitleSize(0.04);
   }
   TAxis* Z = hst -> GetZaxis();
   if ( Z ) {
      Z -> SetLabelFont(62);
      Z -> SetLabelSize(0.04);
      Z -> SetTitleFont(62);
      Z -> SetTitleSize(0.04);
   }
}

//--------------------------------------------------------------------
TH1D* get_mass_hist(string fname, string hname) {
//--------------------------------------------------------------------
#include "cuts.h"

   // name of folder with root files
//    static string dir("prod-9/");
   static string dir("prod-11/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot -> cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory -> Get("a4c");

   double bound = 0.075; // shift_eta=6sigma
   if ( shift_eta+weta > bound ) {
      bound = 0.085;     // shift_eta=7sigma
   }
   int nbins = int(2000*bound+0.5);
   TH1D* mgg = new TH1D(hname.c_str(),
                        ";M^{ inv}_{ #gamma#gamma}, GeV/c^{2}"
                        ";Entries / 1 MeV/c^{2}",
                        nbins, Meta-bound,Meta+bound);


   string select = string("Mgg>>") + hname;
   TCut c_here = c_Mrec+c_chi2+c_phi;
   a4c -> Draw(select.c_str(),c_here,"goff");

   return mgg;
}

//--------------------------------------------------------------------
void GaussFit(TH1D* hist) {
//--------------------------------------------------------------------
   const double ns = 2.; // number of sigmas to fit
   TF1* gs = (TF1*)gROOT -> GetFunction("gaus");
   gs -> SetLineWidth(2);
   gs -> SetLineColor(kRed);
   hist -> Fit(gs,"Q0"); // 1st fit
   double gsmean = gs -> GetParameter(1);
   double gssig  = gs -> GetParameter(2);
   hist -> Fit(gs,"Q0","",gsmean-ns*gssig,gsmean+ns*gssig);
   gsmean = gs -> GetParameter(1);
   gssig  = gs -> GetParameter(2);
   hist -> Fit(gs,"","",gsmean-ns*gssig,gsmean+ns*gssig); // 2nd fit
   gsmean = gs -> GetParameter(1);
   gssig  = gs -> GetParameter(2);
   hist -> Fit(gs,"","",gsmean-ns*gssig,gsmean+ns*gssig); // 3d fit
   hist -> DrawCopy();
}

//--------------------------------------------------------------------
void plot_mass_eta() {
//--------------------------------------------------------------------
#include "cuts.h"

   vector<string> fnames = {
      "data_09psip_all.root",
      "data_12psip_all.root",
      "mcinc_09psip_all.root",
      "mcinc_12psip_all.root",
      "mcsig_kkmc_09.root",
      "mcsig_kkmc_12.root"
   };

   vector<string> date = {
      "2009",
      "2012"
   };
   vector<string> titles = {
      "Data",
      "MC inclusive",
      "MC signal #phi#eta",
   };

   TH1D* mgg[10];
   for (int i = 0; i < 6; i++) {
      string hn = string("mgg_")+to_string(i);
      mgg[i] = get_mass_hist(fnames[i],hn);

      SetHstFace(mgg[i]);
      mgg[i] -> GetXaxis() -> SetTitleOffset(1.1);
      mgg[i] -> GetYaxis() -> SetMaxDigits(3);
      mgg[i] -> SetLineWidth(2);
      mgg[i] -> SetLineColor(kBlack);

      if ( i < 2 ) {
         mgg[i] -> SetOption("E");
      }
   }

   TLine* lR = new TLine;
   lR -> SetLineColor(kRed+1);
   lR -> SetLineWidth(2);
   lR -> SetLineStyle(7);
   TLine* lB = new TLine;
   lB -> SetLineColor(kBlue+1);
   lB -> SetLineWidth(2);
   lB -> SetLineStyle(7);

   TCanvas* c1 = new TCanvas("c1","note",0,0,800,1050);
   int nx = 2;
   int ny = 3;
   c1 -> Divide(nx,ny);
   gStyle -> SetStatX(0.99);
   gStyle -> SetStatW(0.24);
   gStyle -> SetStatH(0.22);
   gStyle -> SetFitFormat(".3g");

   vector<TPaveText*> pt(nx*ny,nullptr);

   int i = 0;
   for ( int y = 0; y < ny; ++y ) {
      for ( int x = 0; x < nx; ++x ) {
         c1 -> cd(i+1);
         GaussFit(mgg[i]);

         double ymax=0.4*mgg[i] -> GetMaximum();
         lR -> DrawLine(Meta-weta,0,Meta-weta,ymax);
         lR -> DrawLine(Meta+weta,0,Meta+weta,ymax);

         double shft  = shift_eta;  // near edge of side-band
         double shftw = shift_eta+weta; // far edge of side-band
         lB -> DrawLine(Meta-shftw,0,Meta-shftw,ymax);
         lB -> DrawLine(Meta-shft, 0,Meta-shft, ymax);
         lB -> DrawLine(Meta+shft, 0,Meta+shft, ymax);
         lB -> DrawLine(Meta+shftw,0,Meta+shftw,ymax);

         string tit = titles[y] + " " + date[x];
         if ( tit.size() < 10 ) { // for data
            tit = "  " + tit + "  ";
         }

         pt[i] = new TPaveText(0.11,0.79,0.4,0.89,"NDC");
         pt[i] -> SetTextAlign(22);
         pt[i] -> SetTextFont(42);
         pt[i] -> AddText( tit.c_str() );
         pt[i] -> Draw();

         i = i+1;
      }
   }

   string pdf("mass_eta.pdf");
   c1 -> Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot_mass_eta_PR() {
//--------------------------------------------------------------------
#include "cuts.h"

   vector<string> fnames = {
      "data_12psip_all.root",
      "mcinc_12psip_all.root",
      "mcsig_kkmc_12.root"
   };

   vector<string> date = {
      "2012"
   };
   vector<string> titles = {
      "Data 2012",
      "MC inclusive",
      "MC signal #phi#eta",
   };

   TH1D* mgg[10];
   for (int i = 0; i < 3; i++) {
      string hn = string("mgg_")+to_string(i);
      mgg[i] = get_mass_hist(fnames[i],hn);
      SetHstFace(mgg[i]);
      mgg[i] -> GetYaxis() -> SetMaxDigits(3);
      mgg[i] -> GetYaxis() -> SetTitleOffset(1.2);
      mgg[i] -> SetLineWidth(2);
      mgg[i] -> SetLineColor(kBlack);
      if ( i == 0 ) {
         mgg[i] -> SetOption("E");
      }
   }

   TLatex* t1 = new TLatex;
   t1 -> SetTextSize(0.08);
   t1 -> SetTextColor(kMagenta+3);

   TLine* lR = new TLine;
   lR -> SetLineColor(kRed+1);
   lR -> SetLineWidth(2);
   lR -> SetLineStyle(7);
   TLine* lB = new TLine;
   lB -> SetLineColor(kBlue+1);
   lB -> SetLineWidth(2);
   lB -> SetLineStyle(7);

//    TCanvas* c1 = new TCanvas("c1","PR",0,0,900,650);
//    int nx = 2;
//    int ny = 2;
   TCanvas* c1 = new TCanvas("c1","PR",0,0,1300,400);
   int nx = 3;
   int ny = 1;
   c1 -> Divide(nx,ny);
//    gStyle -> SetStatX(0.99);
   gStyle -> SetStatW(0.24);
   gStyle -> SetStatH(0.22);
   gStyle -> SetFitFormat(".3g");

   vector<TPaveText*> pt(nx*ny,nullptr);

   int i = 0;
   for ( int y = 0; y < ny; ++y ) {
      for ( int x = 0; x < nx; ++x ) {
         c1 -> cd(i+1);
         if ( i < 3 ) {
            GaussFit(mgg[i]);
         } else {
            gPad -> Clear();
            continue;
         }

         double ymax=0.4*mgg[i] -> GetMaximum();
         lR -> DrawLine(Meta-weta,0,Meta-weta,ymax);
         lR -> DrawLine(Meta+weta,0,Meta+weta,ymax);

         double shft  = shift_eta;  // near edge of side-band
         double shftw = shift_eta+weta; // far edge of side-band
         lB -> DrawLine(Meta-shftw,0,Meta-shftw,ymax);
         lB -> DrawLine(Meta-shft, 0,Meta-shft, ymax);
         lB -> DrawLine(Meta+shft, 0,Meta+shft, ymax);
         lB -> DrawLine(Meta+shftw,0,Meta+shftw,ymax);

//          double yy=0.9*mgg[i] -> GetMaximum();
//          t1 -> DrawLatex(0.5+seta-shft_eta,yy,date[0].c_str());
//          t1 -> DrawLatex(0.5+seta-shft_eta,0.8*yy,titles[i].c_str());

         pt[i] = new TPaveText(0.11,0.79,0.4,0.89,"NDC");
         pt[i] -> SetTextAlign(22);
         pt[i] -> SetTextFont(42);
         pt[i] -> AddText( titles[x].c_str() );
         pt[i] -> Draw();

         i = i+1;
      }
   }

   string pdf("mass_eta_PR.pdf");
   c1 -> Print(pdf.c_str());
}

//--------------------------------------------------------------------
void mass_eta() {
//--------------------------------------------------------------------
   gROOT -> Reset();
   gStyle -> SetOptStat(0);
   gStyle -> SetOptFit(112); // print all parameters (fixed)
//    gStyle -> SetFitFormat(".3g");
   gStyle -> SetStatFont(42);
//    gStyle -> SetStatX(0.89);
//    gStyle -> SetStatY(0.89);

   plot_mass_eta(); // memo
//    plot_mass_eta_PR(); // presentation OLD
}
