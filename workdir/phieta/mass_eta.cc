// plot M(gamma gamma) distributions
// after chi^2 cut for Mphi region
// -> mass_eta_(PR).pdf

#include "masses.h"

// {{{1 helper functions and constants
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

//--------------------------------------------------------------------
void GaussFit(TH1D* hist) {
//--------------------------------------------------------------------
   const double ns = 2.0; // number of sigmas to fit
   TF1* gs = (TF1*)gROOT->GetFunction("gaus");
   gs->SetLineWidth(2);
   gs->SetLineColor(kRed);
   double mean = hist->GetMean();
   double sig  = hist->GetRMS();
   hist->Fit(gs,"Q0","",mean-ns*sig,mean+ns*sig);
   mean = gs->GetParameter(1);
   sig  = gs->GetParameter(2);
   hist->Fit(gs,"Q0","",mean-ns*sig,mean+ns*sig); // 1st fit
   mean = gs->GetParameter(1);
   sig  = gs->GetParameter(2);
   hist->Fit(gs,"Q0","",mean-ns*sig,mean+ns*sig); // 2nd fit
   mean = gs->GetParameter(1);
   sig  = gs->GetParameter(2);
   hist->Fit(gs,"","",mean-ns*sig,mean+ns*sig); // 3d fit
   hist->DrawCopy();
}

// {{{1 Fill histograms
//--------------------------------------------------------------------
TH1D* get_mass_hist(string fname, string hname, bool isMC=false) {
//--------------------------------------------------------------------
#include "cuts.h"

   // name of folder with root files
   static string dir("Ntpls/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("SelectKKgg");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   TCut c_here = c_chi2 + c_phi;
   // TCut c_here = c_chi2 + c_phiT;
   if (isMC) {
      c_here += c_xisr + c_MCmkk; // X_isr>0.9 && mc_Mkk<1.08
   }

   TH1D* mgg = new TH1D(hname.c_str(),
         ";M^{ inv}_{ #gamma#gamma}, GeV/c^{2}"
         ";Entries/0.002 GeV/c^{2}",
         75, Meta-0.075,Meta+0.075);

   a4c->Draw( Form("Mgg>>%s",hname.c_str()), c_here, "goff" );

   return mgg;
}

// {{{1 Plot mass_eta
//--------------------------------------------------------------------
void plot_mass_eta(bool presentation=true) {
//--------------------------------------------------------------------
#include "cuts.h"

   vector<string> fnames = {
        "ntpl_3080_rs.root",
        "ntpl_J4.root",
        "ntpl_3097.root",
        "ntpl_mcgpj_3080_rs.root",
        "ntpl_mcgpj_J4.root",
        "ntpl_mcgpj_3097.root",
   };
   vector<string> titles = {
      "3080MeV R-scan",
      "3097MeV (2018)",
      "3097MeV (2012)",
      "3080MeV MC #phi#eta",
      "3097MeV MC #phi#eta 2018",
      "3097MeV MC #phi#eta 2012",
   };

   size_t Nhst = fnames.size();
   vector<TH1D*> mgg( Nhst, nullptr );
   vector<TPaveText*> pt( Nhst, nullptr );
   for ( size_t i = 0; i < Nhst; ++i ) {
      string hn( Form("mgg_%zu",i) );
      bool isMC = (i >= 3);
      mgg[i] = get_mass_hist(fnames[i],hn,isMC);

      SetHstFace(mgg[i]);
      mgg[i]->GetYaxis()->SetMaxDigits(3);
      mgg[i]->GetXaxis()->SetTitleOffset(1.1);
      mgg[i]->SetLineWidth(1);
      mgg[i]->SetLineColor(kBlack);
      if ( !isMC ) { // data
         mgg[i]->SetOption("E");
         mgg[i]->SetMarkerStyle(20);
         mgg[i]->SetMarkerSize(0.5);
      }

      pt[i] = new TPaveText(0.11,0.79,0.4,0.89,"NDC NB");
      pt[i]->SetTextAlign(12);
      pt[i]->SetTextFont(42);
      pt[i]->AddText( titles[i].c_str() );
   }

   //------------------------------------------------------------------

   // lines for central part and side-band
   TLine* lR = new TLine;
   lR->SetLineColor(kRed+2);
   lR->SetLineWidth(2);
   lR->SetLineStyle(7);
   TLine* lB = new TLine;
   lB->SetLineColor(kBlue+2);
   lB->SetLineWidth(2);
   lB->SetLineStyle(7);

   // Draw
   TCanvas* c1 = nullptr;
   size_t nx = 3, ny = 2;
   vector<size_t> map_pads {0,1,2,3,4,5};
   string pdf("mass_eta");
   if ( presentation ) {
      c1 = new TCanvas("c1","PR",0,0,1100,500);
      pdf += "_PR";
   } else {
      c1 = new TCanvas("c1","MEMO",0,0,900,700);
      nx = 2;
      map_pads = {0,1,3,4};
   }
   c1->Divide(nx,ny);

   gStyle->SetStatX(0.99);
   gStyle->SetStatW(0.24);
   gStyle->SetStatH(0.22);
   gStyle->SetFitFormat(".3g");

   for ( size_t j = 0; j < nx*ny; ++j ) {
      c1->cd(j+1);

      auto i = map_pads[j];
      GaussFit( mgg[i] );
      pt[i]->Draw();

      double ymax=0.4*mgg[i]->GetMaximum();
      lR->DrawLine(Meta-weta,0,Meta-weta,ymax);
      lR->DrawLine(Meta+weta,0,Meta+weta,ymax);

      double shft  = shift_eta;  // near edge of side-band
      double shftw = shift_eta+weta; // far edge of side-band
      lB->DrawLine(Meta-shftw,0,Meta-shftw,ymax);
      lB->DrawLine(Meta-shft, 0,Meta-shft, ymax);
      lB->DrawLine(Meta+shft, 0,Meta+shft, ymax);
      lB->DrawLine(Meta+shftw,0,Meta+shftw,ymax);

      gPad->RedrawAxis();
   }

   c1->Update();
   // pdf += "T.pdf"; // for tight cuts
   pdf += ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Main
//--------------------------------------------------------------------
void mass_eta() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(112); // print all parameters (fixed)
   gStyle->SetStatFont(62);

   // plot_mass_eta();
   plot_mass_eta(false); // for memo
}
