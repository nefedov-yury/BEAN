// mass_2D.cc
// 1) plot M(gamma,gamma) vs M(K+K-)
//    -> Mass2D_{data|mcinc}{Year}.pdf;
//
// 2) Dalitz plots (for >=prod-9m)
//    after 4C kinematic fit and chi^2 cut
//      type = 1 M^2(K-eta) vs M^2(K+eta)
//      type = 2 M^2(K+K-) vs M^2(K+eta)
//    -> Dalitz{type}_{data|mcinc}{Year}.pdf

#include "masses.h"

// {{{1 helper functions and constants
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

// {{{1 plot 1D - M(gg) & M(KK)
//--------------------------------------------------------------------
void plot_Mgg(string fname, string title, string pdf="",
      int Cx=800, int Cy=800)
//--------------------------------------------------------------------
{
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

   auto name = Form("c1_gg_%s",pdf.c_str());
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);

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
void plot_Mkk(string fname, string title, string pdf="",
      int Cx=800, int Cy=800)
//--------------------------------------------------------------------
{
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

   auto name = Form("c1_kk_%s",pdf.c_str());
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);

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
void plot_2D(string fname, string title, string pdf, int Cx, int Cy)
//--------------------------------------------------------------------
{
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
         // 100, 0.0, 1.2,   // gg
         100, 0.0, 1.1,   // gg
         100, 0.9, 2.0  );// KK

   a4c->Draw("Mkk:Mgg>>m2d",c_Mrec+c_chi2,"goff");

   auto name = Form("c2_%s",pdf.c_str());
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);

   c1->cd();
   SetHstFace(m2d);
   // m2d->GetXaxis()->SetTitleOffset(0.9);
   m2d->GetXaxis()->SetTitleOffset(1.1);
   m2d->GetYaxis()->SetTitleOffset(1.1);
   m2d->SetMarkerStyle(20);
   m2d->SetMarkerSize(0.3);

   // m2d->Draw("SCAT"); // depricated in root6.30
   gPad->SetLogz();
   m2d->Draw("COLZ0"); // show 0 bins

   c1->Update();
   // scale width of palette
   TPaletteAxis * palette =
      (TPaletteAxis *) m2d->GetListOfFunctions()->FindObject("palette");
   if ( palette ) {
      palette->SetX2NDC(0.93);
      c1->Modified();
      c1->Update();
   }

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
void plot_Dalitz(string fname, string title, int type, string pdf="",
      int Cx=873, int Cy=800)
//--------------------------------------------------------------------
{
#include "cuts.h"

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
      n = a4c->Draw("M2kmeta:M2kpeta",c_here,"goff");
   } else if ( type == 2 ) {
      n = a4c->Draw("(Mkk*Mkk):M2kpeta",c_here,"goff");
   }

   double* m2X = a4c->GetVal(1);
   double* m2Y = a4c->GetVal(0);
   cout << n << endl;
   n = min(n,50000); // ~ 2021 data

   auto name = Form("c2d_%s",pdf.c_str());
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);

   c1->cd();
   // gPad->SetGrid();

   TGraph* gr = new TGraph(n, m2X, m2Y);

   // set ranges:
   // https://root-forum.cern.ch/t/how-to-set-ranges-on-axis/28254
   gr->GetXaxis()->SetLimits(0.,7.);
   if ( type == 1 ) {
      gr->SetTitle(";M^{2}(K^{#plus}#eta), GeV^{2}/c^{4}"
            ";M^{2}(K^{#minus}#eta), GeV^{2}/c^{4}");
      gr->GetHistogram()->SetMinimum(0.);
      gr->GetHistogram()->SetMaximum(7.);
   } else if ( type == 2 ) {
      gr->SetTitle(";M^{2}(K^{#plus}#eta), GeV^{2}/c^{4}"
            ";M^{2}(K^{#plus}K^{#minus}), GeV^{2}/c^{4}");
      gr->GetHistogram()->SetMinimum(0.);
      gr->GetHistogram()->SetMaximum(7.);
   }
   // gr->GetYaxis()->SetTitleOffset(1.2);

   gr->Draw("AP");

   double width = max(0.2,0.015*title.size());
   TLegend* leg = new TLegend(0.89-width,0.79,0.89,0.89);
   leg->SetTextSize(0.04);
   leg->SetHeader(title.c_str(), "C");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 Main
//--------------------------------------------------------------------
void mass_2D()
//--------------------------------------------------------------------
{
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);
   gStyle->SetLegendFont(42);

   //========================================================
   // set the name of the folder with the root files
   Dir = "prod_v709n4/";
   //========================================================

   size_t Cx = 800, Cy = 640; // canvas sizes, X/Y = 1.25

   // Mass2D, fig.8
   // for ( int date : {2009, 2012, 2021} ) {
   for ( int date : {2021}) {
      // DATA
      string fnd( Form("data_%02ipsip_all.root",date%100) );
      string titd( Form("Data %d",date) );
      string pdfd( Form("Mass2D_data%d.pdf",date%100) );
      // plot_2D(Dir+fnd, titd, pdfd, Cx, Cy);
      // MC incl.
      string fnm( Form("mcinc_%02ipsip_all.root",date%100) );
      string titm( Form("MC inclusive %d",date) );
      string pdfm( Form("Mass2D_mcinc%d.pdf",date%100) );
      // plot_2D(Dir+fnm, titm, pdfm, Cx, Cy);

      // OLD: 1D-plots
      // plot_Mgg(Dir+fnd, titd);
      // plot_Mkk(Dir+fnd, titd);
   }

   // Dalitz, fig.13 or 14
   // type = 1 M^2(K-eta) vs M^2(K+eta)
   // type = 2 M^2(K+K-) vs M^2(K+eta)
   int type = 1;
   // for ( int date : {2009, 2012, 2021} ) {
   for ( int date : {2021}) {
      // DATA
      string fnd( Form("data_%02ipsip_all.root",date%100) );
      string titd( Form("Data %d",date) );
      string pdfd( Form("Dalitz%i_data%d.pdf",type,date%100) );
      plot_Dalitz(Dir+fnd, titd, type, pdfd, Cx, Cy);
      // MC signal
      string fnm( Form("mcsig_kkmc_%02i.root",date%100) );
      string titm( Form("MC signal #phi#eta %d",date) );
      string pdfm( Form("Dalitz%i_mcsig%d.pdf",type,date%100) );
      plot_Dalitz(Dir+fnm, titm, type, pdfm, Cx, Cy);
   }
}
