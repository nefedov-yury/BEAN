// after 4C kinematic fit and chi^2 cut:
// 1) plot M(gamma,gamma) vs M(K+K-)
//    -> Mass2D_XXXX.pdf; 
//
// 2) Dalitz plots ( for >=prod-2 )
//    type = 1 M^2(K-eta) vs M^2(K+eta)
//    type = 2 M^2(K+K-) vs M^2(K+eta)
//    -> Dalitz_xxx.pdf

#include "masses.h"

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
void plot_2D(string fname, string title, string pdf="") {
//--------------------------------------------------------------------
#include "cuts.h"

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot -> cd("SelectKKgg");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   TH2D* m2d = new TH2D("m2d",
                        ";M^{ inv}_{ #gamma#gamma} , GeV/c^{2}"
                        ";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2}",
                        100, 0.0, 1.2,   // gg
                        100, 0.9, 2.0  );// KK

   a4c -> Draw("Mkk:Mgg>>m2d",c_chi2,"goff");

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,873,800);

   c1 -> cd();
   SetHstFace(m2d);
   m2d -> GetYaxis()->SetTitleOffset(1.2);
   m2d -> SetMarkerStyle(20);
   m2d -> SetMarkerSize(0.5);

   m2d -> Draw();

   TEllipse* l = new TEllipse;
   l -> SetLineColor(kRed+1);
   l -> SetLineWidth(3);
   double Reta=0.05;
   double Rphi=0.04;
   l -> SetFillStyle(0);
   l -> DrawEllipse(Meta,Mphi,Reta,Rphi,0.,360.,0.);

   TLatex* t = new TLatex;
   t -> SetTextFont(42);
   t -> SetTextColor(kRed+1);
   t -> DrawLatex(0.6,0.94,"#phi#eta signal");

   double width = 0.02*title.size();
   TLegend* leg = new TLegend(0.89-width,0.79,0.89,0.89);
   leg -> SetHeader(title.c_str(), "C");
   leg -> Draw();

   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// Dalitz plot for KKeta
// type = 1 M^2(K-eta) vs M^2(K+eta)
// type = 2 M^2(K+K-) vs M^2(K+eta)
//--------------------------------------------------------------------
void plot_Dalitz(string fname, int type, string title, string pdf="") {
//--------------------------------------------------------------------
#include "cuts.h"

   if ( type == 1 ) {
      pdf = string("Dalitz_") + pdf;
   } else if ( type == 2 ) {
      pdf = string("Dalitz2_") + pdf;
   } else {
      cerr << "wrong type= " << type << endl;
      exit(0);
   }

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot -> cd("SelectKKgg");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   TCut c_here = c_chi2; // common cuts
   c_here += c_cpgg;            // only eta
//    c_here += TCut("Mkk<1.08");  // only near phi

   int n = 0;
   if ( type == 1 ) {
      n = a4c -> Draw("M2kmeta:M2kpeta",c_here,"goff");
   } else if ( type == 2 ) {
      n = a4c -> Draw("(Mkk*Mkk):M2kpeta",c_here,"goff");
   }

   double* m2X = a4c -> GetVal(1);
   double* m2Y = a4c -> GetVal(0);
   cout << n << endl;

   //-----------------------------------------------------------------

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1 -> cd();
//   gPad -> SetGrid();

   TGraph* gr = new TGraph(n, m2X, m2Y);

   // set ranges: https://root-forum.cern.ch/t/how-to-set-ranges-on-axis/28254
   gr -> GetXaxis() -> SetLimits(0.,7.);
   if ( type == 1 ) {
      gr -> SetTitle(";M^{2}(K^{+}#eta), GeV^{2}/c^{4}"
                     ";M^{2}(K^{-}#eta), GeV^{2}/c^{4}");
      gr -> GetHistogram() -> SetMinimum(0.);
      gr -> GetHistogram() -> SetMaximum(7.);
   } else if ( type == 2 ) {
      gr -> SetTitle(";M^{2}(K^{+}#eta), GeV^{2}/c^{4}"
                     ";M^{2}(K^{+}K^{-}), GeV^{2}/c^{4}");
      gr -> GetHistogram() -> SetMinimum(0.);
      gr -> GetHistogram() -> SetMaximum(7.);
   }
//    gr -> GetYaxis() -> SetTitleOffset(1.2);

   gr -> Draw("AP");

//    double width = max(0.2,0.015*title.size());
   double width = 0.021*title.size();
   TLegend* leg = new TLegend(0.89-width,0.79,0.89,0.89);
   leg -> SetHeader(title.c_str(), "C");
   leg -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//--------------------------------------------------------------------
void mass_2D() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);
   gStyle->SetLegendFont(42);
   gStyle->SetLegendTextSize(0.04); // must be here

   vector<string> fnames  {
      "ntpl_3080_rs.root",
      "ntpl_3080_2019.root",
      "ntpl_3080_SUM.root",
      "ntpl_mcgpj_3080_rs.root",
      "ntpl_mcgpj_3080_2019.root",

      "ntpl_3097.root",
      "ntpl_J4.root",
      "ntpl_mcgpj_3097.root",
      "ntpl_mcgpj_J4.root",

      "ntpl_jpsi_incl.root",
      "ntpl_qq_kkmc_3080.root",
      "ntpl_mcgpj_PHSP_3080_2019.root",
      "ntpl_mcgpj_J4.root",
   };

   vector<string> titles = {
      "3080MeV R-scan",              // iset = 0
      "3080MeV (2019)",
      "3080MeV R-scan+2019",
      "MCGPJ #phi#eta at 3080MeV",      // R-scan: BOSS-704
      "MC #phi#eta at 3080MeV",         // tau-scan 2018: BOSS-704

      "3097MeV J/#Psi-scan", // iset = 5
      "3097MeV (2018)",
      "MCGPJ #phi#epa at 3097MeV",     // BOSS-664
      "MC #phi#epa at 3097MeV",        // BOSS-704

      "MC J/#Psi inclusive",           // iset = 8 BOSS-664
      "KKMC at 3080MeV",               // BOSS-664
      "MCGPJ non-#phi KK#eta at 3080", // Mkk < 1.082 in MC
      "MCGPJ non-#phi KK#eta at 3080", // Mkk < 1.082 in MC
   };

   vector<string> pdfs = {
      "3080_rs.pdf",
      "3080_2019.pdf",
      "3080_sum.pdf",
      "mcphieta_3080_rc.pdf",
      "mcphieta_3080_2019.pdf",

      "3097_2012.pdf",
      "3097_2018.pdf",
      "mcphieta_3097_2012.pdf",
      "mcphieta_3097_2018.pdf",

      "Jpsi_incl.pdf",
      "KKMC_3080.pdf",
      "kketa_3080_rc.pdf",
      "kketa_3080_2019.pdf",
   };

   bool Mass2D = true;          // <--- Mass2D/Dalitz
   vector<int> list {0,1,5,6};
   for (int iset : list) {
      string fn = fnames.at(iset);
      string fname = "Ntpls/" + fn;
      if (Mass2D) {
         string pdf = "Mass2D_" + pdfs[iset];
         plot_2D(fname, titles[iset], pdf);
      } else {
         string pdf = "Dalitz2_" + pdfs[iset];
         plot_Dalitz(fname, 2, titles[iset], pdfs[iset]);
      }
   }
}
