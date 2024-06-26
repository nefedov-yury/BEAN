// 1) plot M(gamma,gamma) vs M(K+K-)
//    -> Mass2D_XXXX.pdf;
//
// 2) Dalitz plots (for >=prod-9m)
//    after 4C kinematic fit and chi^2 cut
//      type = 1 M^2(K-eta) vs M^2(K+eta)
//      type = 2 M^2(K+K-) vs M^2(K+eta)
//    -> Dalitz_xxx.pdf

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
void plot_Mgg(string fname, string title, string pdf="") {
//--------------------------------------------------------------------
#include "cuts.h"

   gStyle->SetOptStat(0);

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");
   TH1D* mgg = new TH1D("mgg", ";M(#gamma#gamma)",
                        400, 0.0, 1.1);   // gg
   TCut c_here = c_Mrec+c_chi2;
   c_here += TCut("1.4<Mkk&&Mkk<1.85");

   a4c->Draw("Mgg>>mgg",c_here,"goff");

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1->cd();
   gPad->SetGrid();
   mgg->Draw();

   TLatex* t1 = new TLatex;
   t1->SetTextSize(0.05);
   t1->SetTextColor(kMagenta+3);
   double y = mgg->GetMaximum();
   t1->DrawLatex(0.6,0.8*y,title.c_str());

   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

//--------------------------------------------------------------------
void plot_Mkk(string fname, string title, string pdf="") {
//--------------------------------------------------------------------
#include "cuts.h"

   gStyle->SetOptStat(0);

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");
   TH1D* mkk = new TH1D("mkk",
                        ";M(K^{#plus}K^{#minus})",
                        100, 0.9, 2.1  );// KK
   TCut c_here = c_Mrec+c_chi2;
   c_here += c_cpgg;
   a4c->Draw("Mkk>>mkk",c_here,"goff");

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1->cd();
   gPad->SetGrid();
   mkk->Draw();

   TLatex* t1 = new TLatex;
   t1->SetTextSize(0.05);
   t1->SetTextColor(kMagenta+3);
   double y = mkk->GetMaximum();
   t1->DrawLatex(1.4,0.8*y,title.c_str());

   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 2D: M(KK) vs M(gg)
//--------------------------------------------------------------------
void plot_2D(string fname, string title, string pdf="") {
//--------------------------------------------------------------------
#include "cuts.h"

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   TH2D* m2d = new TH2D("m2d",
         ";M^{ inv}_{ #gamma#gamma}, GeV/c^{2}"
         ";M^{ inv}_{ K^{#plus}K^{#minus }}, GeV/c^{2}",
         100, 0.0, 1.2,   // gg
         100, 0.9, 2.0  );// KK

   a4c->Draw("Mkk:Mgg>>m2d",c_Mrec+c_chi2,"goff");

   TCanvas* c1 = new TCanvas("c1","...",0,0,873,800);

   c1->cd();
   SetHstFace(m2d);
   // m2d->GetXaxis()->SetTitleOffset(0.9);
   m2d->GetYaxis()->SetTitleOffset(1.2);
   m2d->SetMarkerStyle(20);
   m2d->SetMarkerSize(0.3);

   m2d->Draw();

   TEllipse* l = new TEllipse;
   l->SetLineColor(kRed+1);
   l->SetLineWidth(3);
   double Reta=0.05;
   double Rphi=0.04;
   l->SetFillStyle(0);
   l->DrawEllipse(Meta,Mphi,Reta,Rphi,0.,360.,0.);

   TLatex* t = new TLatex;
   t->SetTextFont(42);
   t->SetTextColor(kRed+1);
   // t->DrawLatex(0.6,0.94,"#phi#eta signal");
   t->DrawLatex(0.6,0.94,"#phi#eta");

   double width = 0.02*title.size();
   TLegend* leg = new TLegend(0.89-width,0.79,0.89,0.89);
   leg->SetTextSize(0.04);
   leg->SetHeader(title.c_str(), "C");
   leg->Draw();

   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 Dalitz plot for KKeta
// type = 1 M^2(K-eta) vs M^2(K+eta)
// type = 2 M^2(K+K-) vs M^2(K+eta)
//--------------------------------------------------------------------
void plot_Dalitz(string fname,int type,string title,string pdf="") {
//--------------------------------------------------------------------
#include "cuts.h"

   if ( type == 1 ) {
      pdf = string("Dalitz_") + pdf;
   } else if ( type == 2 ) {
      pdf = string("Dalitz2_") + pdf;
   } else {
      cerr << "wrong type= " << type << endl;
      exit(EXIT_FAILURE);
   }

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   TCut c_here = c_Mrec+c_chi2; // common cuts
   c_here += c_cpgg;            // only eta
   // c_here += TCut("Mkk<1.08");  // only near phi

   int n = 0;
   if ( type == 1 ) {
      n = a4c -> Draw("M2kmeta:M2kpeta",c_here,"goff");
   } else if ( type == 2 ) {
      n = a4c -> Draw("(Mkk*Mkk):M2kpeta",c_here,"goff");
   }

   double* m2X = a4c -> GetVal(1);
   double* m2Y = a4c -> GetVal(0);
   cout << n << endl;
   n = min(n,65000); // > 10*Ndata12

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1 -> cd();
   // gPad -> SetGrid();

   TGraph* gr = new TGraph(n, m2X, m2Y);

   // set ranges:
   // https://root-forum.cern.ch/t/how-to-set-ranges-on-axis/28254
   gr -> GetXaxis() -> SetLimits(0.,7.);
   if ( type == 1 ) {
      gr -> SetTitle(";M^{2}(K^{#plus}#eta), GeV^{2}/c^{4}"
                     ";M^{2}(K^{#minus}#eta), GeV^{2}/c^{4}");
      gr -> GetHistogram() -> SetMinimum(0.);
      gr -> GetHistogram() -> SetMaximum(7.);
   } else if ( type == 2 ) {
      gr -> SetTitle(";M^{2}(K^{#plus}#eta), GeV^{2}/c^{4}"
                     ";M^{2}(K^{#plus}K^{#minus}), GeV^{2}/c^{4}");
      gr -> GetHistogram() -> SetMinimum(0.);
      gr -> GetHistogram() -> SetMaximum(7.);
   }
   // gr -> GetYaxis() -> SetTitleOffset(1.2);

   gr -> Draw("AP");

   double width = max(0.2,0.015*title.size());
   TLegend* leg = new TLegend(0.89-width,0.79,0.89,0.89);
   leg -> SetTextSize(0.04);
   leg -> SetHeader(title.c_str(), "C");
   leg -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

// {{{1 Main
//--------------------------------------------------------------------
void mass_2D() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);
   gStyle->SetLegendFont(42);

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

   // string dir = string("prod-12/");
   string dir("prod_v709/");

   bool Mass2D = true;          // <--- Mass2D/Dalitz
   vector<int> list {0,1,2, 4,5,6};
   if ( !Mass2D ) {
      list = {0,1,2};
   }

   for (int iset : list) {
      string fn = fnames.at(iset);
      string::size_type idx = fn.find("_");
      string fname = dir + fn;

      if (Mass2D) {
         // Mass2D .. fig.9
         string pdf = string("Mass2D_") + fn.substr(0,idx)
                      + fn.substr(idx+1,2) + ".pdf";
         plot_2D(fname, titles[iset], pdf);
      } else {
         // Dalitz .. fig.14
         string pdf = fn.substr(0,idx);
         if ( iset < 5 ) {
            pdf += fn.substr(idx+1,2) + ".pdf";
         } else {
            pdf += fn.substr(idx+6,2) + ".pdf";
         }
         // type = 1 M^2(K-eta) vs M^2(K+eta)
         // type = 2 M^2(K+K-) vs M^2(K+eta)
         int type = 1;
         plot_Dalitz(fname, type, titles[iset], pdf);
      }

      gSystem->ProcessEvents(); //Process pending events
      if ( iset != list.back() ) {
         getchar(); // pause
      }
   }

   // plot_Mgg(fname, titles.at(iset)); //, "Mgg_09.pdf");
   // plot_Mkk(fname, titles.at(iset));
}
