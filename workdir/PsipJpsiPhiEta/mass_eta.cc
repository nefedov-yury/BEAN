// mass_eta.cc
// plot M(gamma gamma) distributions for data and MC
// fit central part, show side-bands
// -> mass_eta_{dat|mi|mc}_{YEAR}.pdf

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

//--------------------------------------------------------------------
void GaussFit(TH1D* hist)
//--------------------------------------------------------------------
{
   const double ns = 2.; // number of sigmas to fit
   TF1* gs = (TF1*)gROOT->GetFunction("gaus");
   gs->SetLineWidth(2);
   gs->SetLineColor(kRed);
   hist->Fit(gs,"Q0"); // 1st fit
   double gsmean = gs->GetParameter(1);
   double gssig  = gs->GetParameter(2);
   hist->Fit(gs,"Q0","",gsmean-ns*gssig,gsmean+ns*gssig);
   gsmean = gs->GetParameter(1);
   gssig  = gs->GetParameter(2);
   hist->Fit(gs,"","",gsmean-ns*gssig,gsmean+ns*gssig); // 2nd fit
   gsmean = gs->GetParameter(1);
   gssig  = gs->GetParameter(2);
   hist->Fit(gs,"","",gsmean-ns*gssig,gsmean+ns*gssig); // 3d fit
   hist->DrawCopy();
}

// {{{1 Fill histograms
//--------------------------------------------------------------------
TH1D* get_mass_hist(string fname, string hname)
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

   // calculate offset from Meta for histogram
   double offset = 0.005*floor(1000*(shift_eta+weta+seta)/5);
   int nbins = int(ceil(2*offset/0.001)); //bin width ~ 0.001
   TH1D* mgg = new TH1D(hname.c_str(),
         ";M^{ inv}_{ #gamma#gamma}, GeV/c^{2}"
         ";Entries / 1 MeV/c^{2}",
         nbins, Meta-offset,Meta+offset);

   string select = string("Mgg>>") + hname;
   TCut c_here = c_Mrec+c_chi2+c_phi;
   a4c->Draw(select.c_str(),c_here,"goff");

   return mgg;
}

// {{{1 M(gg)
//--------------------------------------------------------------------
void plot_mass_eta(int date, int mc, int Cx=600, int Cy=600)
//--------------------------------------------------------------------
{
#include "cuts.h"

   string fname, hn, title, pdf;
   switch (mc) {
      case 0: // data
         fname = string( Form("data_%02ipsip_all.root",date%100) );
         hn = string( Form("mgg_data_%i",date) );
         title = string( Form("Data %i",date) );
         pdf = string( Form("mass_eta_dat_%02i.pdf", date%100) );
         break;
      case 1: // MC inc
         fname = string( Form("mcinc_%02ipsip_all.root",date%100) );
         hn = string( Form("mgg_mcinc_%i",date) );
         // title = string( Form("MC inclusive %02i",date%100) );
         title = string( "MC inclusive" );
         pdf = string( Form("mass_eta_mi_%02i.pdf", date%100) );
         break;
      case 2: // MC sig
         fname = string( Form("mcsig_kkmc_%02i.root",date%100) );
         hn = string( Form("mgg_mcsig_%i",date) );
         // title = string( Form("MC signal #phi#eta %02i",date%100) );
         title = string( "MC signal #phi#eta" );
         pdf = string( Form("mass_eta_ms_%02i.pdf", date%100) );
         break;
      default:
         printf("ERROR: mc= %i\n",mc);
         exit(1);
   }

   TH1D* mgg = get_mass_hist(fname,hn);

   SetHstFace(mgg);
   mgg->GetXaxis()->SetTitleOffset(1.1);
   mgg->GetYaxis()->SetTitleOffset(1.2);
   mgg->GetYaxis()->SetMaxDigits(3);
   mgg->SetLineWidth(2);
   mgg->SetLineColor(kBlack);
   if ( mc == 0 ) {
      mgg->SetOption("E");
   }

   TLine* lR = new TLine;
   lR->SetLineColor(kRed+1);
   lR->SetLineWidth(2);
   lR->SetLineStyle(7);
   TLine* lB = new TLine;
   lB->SetLineColor(kBlue+1);
   lB->SetLineWidth(2);
   lB->SetLineStyle(7);

   TPaveText* pt = new TPaveText(0.11,0.79,0.4,0.89,"NDC");
   pt->SetTextAlign(22);
   pt->SetTextFont(42);
   pt->AddText( title.c_str() );

   auto name = Form("c1_%i_%i",mc,date);
   TCanvas* c1 = new TCanvas(name,name,20*mc,10*(date-2009),Cx,Cy);
   c1->cd();

   GaussFit(mgg);

   double ymax=0.4*mgg->GetMaximum();
   lR->DrawLine(Meta-weta,0,Meta-weta,ymax);
   lR->DrawLine(Meta+weta,0,Meta+weta,ymax);

   double shft  = shift_eta;  // near edge of side-band
   double shftw = shift_eta+weta; // far edge of side-band
   lB->DrawLine(Meta-shftw,0,Meta-shftw,ymax);
   lB->DrawLine(Meta-shft, 0,Meta-shft, ymax);
   lB->DrawLine(Meta+shft, 0,Meta+shft, ymax);
   lB->DrawLine(Meta+shftw,0,Meta+shftw,ymax);

   pt->Draw();
   gPad->RedrawAxis();

   c1->Update();
   c1->Print(pdf.c_str());
}

// {{{1 M(gg) presentation TODO: is it OLD?
//--------------------------------------------------------------------
void plot_mass_eta_PR()
//--------------------------------------------------------------------
{
#include "cuts.h"

   vector<string> fnames = {
      "data_21psip_all.root",
      "mcinc_21psip_all.root",
      "mcsig_kkmc_21.root"
   };

   vector<string> titles = {
      "Data 2021",
      "MC inclusive",
      "MC signal #phi#eta",
   };

   TH1D* mgg[10];
   size_t Nhst = fnames.size();
   for ( size_t i = 0; i < Nhst; ++i ) {
      string hn = string("mgg_")+to_string(i);
      mgg[i] = get_mass_hist(fnames[i],hn);

      SetHstFace(mgg[i]);
      mgg[i]->GetXaxis()->SetTitleOffset(1.1);
      mgg[i]->GetYaxis()->SetMaxDigits(3);
      mgg[i]->SetLineWidth(2);
      mgg[i]->SetLineColor(kBlack);

      if ( fnames[i].find("data") != string::npos ) {
         mgg[i]->SetOption("E");
      }
   }

   TLine* lR = new TLine;
   lR->SetLineColor(kRed+1);
   lR->SetLineWidth(2);
   lR->SetLineStyle(7);
   TLine* lB = new TLine;
   lB->SetLineColor(kBlue+1);
   lB->SetLineWidth(2);
   lB->SetLineStyle(7);

   TCanvas* c1 = new TCanvas("c1","PR",0,0,1300,400);
   int nx = 3;
   int ny = 1;
   c1->Divide(nx,ny);
   // gStyle->SetStatX(0.99);
   gStyle->SetStatW(0.24);
   gStyle->SetStatH(0.22);
   gStyle->SetFitFormat(".3g");

   vector<TPaveText*> pt(nx*ny,nullptr);

   int i = 0;
   for ( int y = 0; y < ny; ++y ) {
      for ( int x = 0; x < nx; ++x ) {
         c1->cd(i+1);
         GaussFit(mgg[i]);

         double ymax=0.4*mgg[i]->GetMaximum();
         lR->DrawLine(Meta-weta,0,Meta-weta,ymax);
         lR->DrawLine(Meta+weta,0,Meta+weta,ymax);

         double shft  = shift_eta;  // near edge of side-band
         double shftw = shift_eta+weta; // far edge of side-band
         lB->DrawLine(Meta-shftw,0,Meta-shftw,ymax);
         lB->DrawLine(Meta-shft, 0,Meta-shft, ymax);
         lB->DrawLine(Meta+shft, 0,Meta+shft, ymax);
         lB->DrawLine(Meta+shftw,0,Meta+shftw,ymax);

         pt[i] = new TPaveText(0.11,0.79,0.4,0.89,"NDC");
         pt[i]->SetTextAlign(22);
         pt[i]->SetTextFont(42);
         pt[i]->AddText( titles[x].c_str() );
         pt[i]->Draw();

         i = i+1;
      }
   }

   string pdf("mass_eta_PR.pdf");
   c1 -> Print(pdf.c_str());
}

// {{{1 Main
//--------------------------------------------------------------------
void mass_eta()
//--------------------------------------------------------------------
{
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);
   gStyle->SetOptFit(112); // print all parameters (fixed)

   gStyle->SetFitFormat(".3g");
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.18);
   gStyle->SetStatH(0.12);

   //========================================================
   // set the name of the folder with the root files
   Dir = "prod_v709n3/";
   //========================================================

   size_t Cx = 880, Cy = 760; // canvas sizes

   for ( int date : {2009, 2012, 2021} ) {
      for ( int mc : {0,2} ) { // data, MC-sig
         plot_mass_eta(date, mc, Cx, Cy); // memo
      }
   }
}
