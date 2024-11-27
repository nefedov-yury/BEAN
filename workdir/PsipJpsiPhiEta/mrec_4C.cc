// mrec_4C.cc
// We use special selection with 4C-kinematic constrints for data,
// continuum and inclusive MC events to estimate the contribution of
// Psi(2S) -> pi+ pi- (phi eta) (non J/Psi) events.
// The 'Mrec' is drawn and number of events in signal and side-band
// regions are calculated.
//      -> Mrec4C_{Year}.pdf

#include <algorithm>
#include "masses.h"

// {{{1 helper functions and constants
//--------------------------------------------------------------------
// GLOBAL: name of folder with root files
string Dir;

//--------------------------------------------------------------------
constexpr double SQ(double x)
//--------------------------------------------------------------------
{
   return x*x;
}

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

// {{{1 Mrec distribution and numbers
//-------------------------------------------------------------------------
tuple<TH1D*,vector<double> > get_Mrec(string fname, string hname)
//-------------------------------------------------------------------------
{
#include "cuts.h"

   fname = Dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* c4c = (TTree*)gDirectory->Get("c4C");

   TH1D* hst = new TH1D(hname.c_str(),
         ";M^{rec}_{#pi^{#plus}#pi^{#minus }}, GeV/c^{2}"
         ";Entries / 0.001GeV/c^{2}",
         200,3.0,3.2 // 1bin = 1MeV
         );
   hst->Sumw2(true);

   TCut c_here = c_chi2;
   c_here += TCut(Form("%f<=Mkk&&Mkk<%f",2*Mk,1.08));
   c_here += c_cpgg;
   // cout << " c_here= " << c_here << endl;

   string dr = string("Mrec>>") + hname;
   c4c->Draw( dr.c_str(), c_here, "goff" );

   // number of events in central part and in side-bands:
   TCut c_Mr0("abs(Mrec-3.097)<0.005"); // [3.092, 3.102]
   TCut c_MrL("Mrec>3.00&&Mrec<3.05");  // [3.00, 3.05] left SB
   TCut c_MrR("Mrec>3.144&&Mrec<3.194");// [3.144, 3.194] right SB
   double n0 = c4c->Draw( "Mrec", c_here+c_Mr0, "goff" );
   double nL = c4c->Draw( "Mrec", c_here+c_MrL, "goff" );
   double nR = c4c->Draw( "Mrec", c_here+c_MrR, "goff" );
   vector<double> vec {n0,nL,nR};

   return make_tuple(hst,vec);
}

// {{{1 draw
//-------------------------------------------------------------------------
void draw_Pict(int date, int Cx, int Cy, string pdf)
//-------------------------------------------------------------------------
{
#include "norm.h"

   printf(" --- %i ----\n", date);

   string fn_data( Form("data_%02ipsip_all.root",date%100) );
   auto [hstd,vec_d] = get_Mrec( fn_data, Form("d4c_%i",date) );
   printf("Data: Ncp=%f, Lsb=%f, Rsb=%f\n",
         vec_d[0], vec_d[1], vec_d[2]);

   string fn_dc("data_3650_2021.root");
   auto [hstc,vec_c] = get_Mrec( fn_dc, "d4c_3650" );
   printf("3650: Ncp=%f, Lsb=%f, Rsb=%f\n",
         vec_c[0], vec_c[1], vec_c[2]);

   string fn_mc( Form("mcinc_%02ipsip_all.root",date%100) );
   auto [hstm,vec_m] = get_Mrec( fn_mc, Form("m4c_%i",date) );
   printf("MCinc: Ncp=%f, Lsb=%f, Rsb=%f\n",
         vec_m[0], vec_m[1], vec_m[2]);

   // string fn_mcs( Form("mcsig_kkmc_%02i.root",date%100) );
   // auto [hstms,vec_ms] = get_Mrec( fn_mcs, Form("ms4c_%i",date) );
   // printf(" mc_sig: Ncp=%f, Lsb=%f, Rsb=%f\n",
         // vec_ms[0], vec_ms[1], vec_ms[2]);

   SetHstFace(hstd);
   hstd->SetLineColor(kBlack); // data
   hstd->SetLineWidth(2);
   hstd->SetMarkerStyle(20);
   hstd->SetMarkerSize(0.7);

   hstc->SetLineColor(kBlue+1); // Continuum data
   hstc->SetLineWidth(2);

   hstm->SetLineColor(kGreen+2); // inclusive MC
   hstm->SetLineWidth(2);

   // hstms->SetLineColor(kGreen+2); // MC phi eta
   // hstms->SetLineWidth(2);


   // data:
   double Ndat = vec_d[0];
   double Nsb = vec_d[1]+vec_d[2];
   double err_Nsb = sqrt(Nsb);
   double delta = Nsb/Ndat /10. *100;
   double err_delta = err_Nsb/Ndat /10. *100;
   printf("\n");
   printf("Data : $%.0f$  $%.0f\\pm%.1f$  $%.2f\\pm%.2f$\n",
         Ndat,Nsb,err_Nsb,delta,err_delta);

   // inclusive MC, normalization on data
   hstm->Scale( MC_Dat.at(date) );
   for ( int i = 0; i < 3; i++ ) {
      vec_m[i] *= MC_Dat.at(date);
   }
   double Nmc = vec_m[0];
   double Nmc_sb = vec_m[1]+vec_m[2];
   double err_Nmc_sb = sqrt(Nmc_sb*MC_Dat.at(date));
   double del_mc = Nmc_sb/Nmc /10. *100;
   double err_del_mc = err_Nmc_sb/Nmc /10. *100;
   printf("MCinc: $%.0f$  $%.0f\\pm%.1f$  $%.2f\\pm%.2f$\n",
         Nmc,Nmc_sb,err_Nmc_sb,del_mc,err_del_mc);

   // continuum data, normalization on data
   // (3.2-3.0)/(3.102-3.092) = 20
   double Ncont= hstc->GetEntries();
   double err_Ncont = sqrt(Ncont);
   Ncont *= C_Dat.at(date) /20;
   err_Ncont *= C_Dat.at(date) /20;
   hstc->Scale( C_Dat.at(date) );
   printf("CONT : $%.2f\\pm%.2f$\n", Ncont,err_Ncont);

   // Draw
   //-----------------------------------------------------------------------

   TLine* lR = new TLine;
   lR->SetLineColor(kRed+1);
   lR->SetLineWidth(2);
   lR->SetLineStyle(7);
   TLine* lB = new TLine;
   lB->SetLineColor(kMagenta);
   lB->SetLineWidth(2);
   lB->SetLineStyle(7);

   auto name = Form("c1_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   gPad->SetLogy(true);
   double hmin = 0.5;
   hstd->SetMinimum(hmin);
   double hmax = hstd->GetMaximum();

   hstd->GetYaxis()->SetTitleOffset(1.1);

   hstd->Draw("E");
   hstc->Draw("HIST SAME");
   hstm->Draw("HIST SAME");
   hstd->Draw("E SAME");

   lR->DrawLine(3.092,hmin,3.092,1.2*hmax);
   lR->DrawLine(3.102,hmin,3.102,1.2*hmax);
   lB->DrawLine(3.001,hmin,3.001,0.03*hmax);
   lB->DrawLine(3.050,hmin,3.050,0.03*hmax);
   lB->DrawLine(3.144,hmin,3.144,0.03*hmax);
   lB->DrawLine(3.194,hmin,3.194,0.03*hmax);

   TLegend* leg = new TLegend(0.60,0.67,0.892,0.89);
   leg->AddEntry( hstd,Form("Data %i",date), "PLE" );
   leg->AddEntry( hstc, "Off-resonance data", "L" );
   leg->AddEntry( hstm, "MC inclusive", "L" );
   leg->AddEntry( lR, "J/#Psi area", "L" );
   leg->AddEntry( lB, "side-band", "L" );
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   c1->Print(pdf.c_str());
}

// {{{1 Main
//-------------------------------------------------------------------------
void mrec_4C()
//-------------------------------------------------------------------------
{
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);

   //========================================================
   // set the name of the folder with the root files
   Dir = "prod_v709n4/";
   //========================================================

   size_t Cx = 800, Cy = 640; // canvas sizes, X/Y = 1.25

   for ( int date : { 2021,2012,2009 } ) {
      string pdf( Form("Mrec4C_%02i.pdf",date%100) );
      draw_Pict(date, Cx,Cy, pdf);
   }

}

