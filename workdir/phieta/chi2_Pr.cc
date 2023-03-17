// plot chi2 of the kinematic fit
// data (window for Mkk is used) for signal ans side-band events
// versus signal MC (plus another MC if any):
// -> chi2_sb_XXXXX.pdf

static double chi2_cut = 0.; // to keep value from cuts.h file

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
TH1D* get_chi2(string fname, string hname, int icuts=1) {
//--------------------------------------------------------------------
#include "cuts.h"
   chi2_cut = chi2M;

   // name of folder with root files
   static string dir("Ntpls/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if ( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("SelectKKgg");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   TCut c_here(""); // no cuts if icuts == 0
   if ( icuts != 0 ) {
      bool isMC = (icuts < 0);
      if (isMC) {
         c_here += c_xisr + c_MCmkk; // X_isr>0.9 && mc_Mkk<1.08
      }
      if ( abs(icuts) == 1 ) { // central part of Mgg
         c_here += c_cpgg;
      } else if ( abs(icuts) == 2 ) { // side-band
         c_here += c_sbgg;
      }
      // c_here += c_phi; // Mkk in [2*Mk, 1.08GeV]
      c_here += c_phiT; // Mkk in [1.01, 1.03GeV]
   }

   TH1D* chi2 = new TH1D(hname.c_str(),
         ";#chi^{2};Entries/2",
         100, 0.,200.);

   a4c->Draw( Form("ch2>>%s",hname.c_str()),c_here,"goff");

   return chi2;
}

//--------------------------------------------------------------------
void chi2_SB( string stype ) {
//--------------------------------------------------------------------
   double maxh = 0.;
   string Dname, MCname, MCinc, MCqq;
   string Tleg, pdf;
   if( stype == "3080_rs" ) {
      Dname = string("ntpl_3080_rs.root");
      MCname= string("ntpl_mcgpj_3080_rs.root");
      Tleg  = string("#chi^{2}(4C): R-scan at 3080MeV");
      pdf   = string("chi2_sb_3080_rs");
      MCqq  = string("ntpl_qq_kkmc_3080.root");
   } else if( stype == "3080_2019" ) {
      Dname = string("ntpl_3080_2019.root");
      MCname= string("ntpl_mcgpj_3080_2019.root");
      Tleg  = string("#chi^{2}(4C): 3080MeV in 2019");
      pdf   = string("chi2_sb_3080_2019");
      MCqq  = string("ntpl_qq_kkmc_3080.root");
   } else if ( stype == "3097" ) {
      Dname = string("ntpl_3097.root");
      MCname= string("ntpl_mcgpj_3097.root");
      Tleg  = string("#chi^{2}(4C): J/#Psi-scan at 3097MeV");
      pdf   = string("chi2_sb_3097");
      MCinc = string("ntpl_jpsi_incl.root");
   } else if ( stype == "3097J" ) {
      Dname = string("ntpl_J4.root");
      MCname= string("ntpl_mcgpj_J4.root");
      Tleg  = string("#chi^{2}(4C): 3097MeV in 2018");
      pdf   = string("chi2_sb_3097J");
      MCinc = string("ntpl_jpsi_incl.root");
   } else if ( stype == "2900_rs" ) {
      Dname = string("ntpl_2900_rs.root");
      MCname= string("ntpl_mcgpj_2900_rs.root");
      Tleg  = string("#chi^{2}(4C): R-scan at 2900MeV");
      pdf   = string("chi2_sb_2900_rs");
   }

   TH1D* chi2cp = get_chi2(Dname,"chi2cp",1);
   SetHstFace(chi2cp);
   chi2cp -> SetLineColor(kBlack);
   chi2cp -> SetLineWidth(2);
   chi2cp -> GetYaxis() -> SetTitleOffset(1.25);
   maxh = max(maxh, chi2cp->GetMaximum());

   TH1D* chi2sb = get_chi2(Dname,"chi2sb",2);
   chi2sb -> SetLineColor(kBlue+2);
   chi2sb -> SetFillStyle(3001);
   chi2sb -> SetFillColor(kBlue+1);

   TH1D* chi2_mcS = get_chi2(MCname,"chi2_mcS",-1);
   chi2_mcS -> SetLineColor(kGreen+2);
   chi2_mcS -> SetLineWidth(2);
   // normalize MC signal to data
   chi2_mcS -> Scale(
         (chi2cp->Integral() - chi2sb->Integral()) /
         chi2_mcS -> Integral()
         );
   maxh = max(maxh, chi2_mcS -> GetMaximum());

   TH1D* chi2mc_inc = nullptr;
   if ( !MCinc.empty() ) {
      chi2mc_inc = get_chi2(MCinc,"chi2_mcinc",-1);
      // normalize MC inclusive to data
      chi2mc_inc -> Scale(
         chi2cp->Integral() / chi2mc_inc -> Integral()
         );
      maxh = max(maxh, chi2mc_inc -> GetMaximum());
      chi2mc_inc -> SetLineColor(kRed+2);
      chi2mc_inc -> SetLineWidth(2);
      chi2mc_inc -> SetLineStyle(kDashed);
   }

   TH1D* chi2mc_qq = nullptr;
   if ( !MCqq.empty() ) {
      TH1D* chi2_0 = get_chi2(Dname,"chi2_0",0);
      TH1D* chi2qq_0 = get_chi2(MCqq,"chi2qq_0",0);
      // normalizition for MC-qq to data without Mqq cuts
      double norm = chi2_0->Integral() / chi2qq_0 -> Integral();

      chi2mc_qq = get_chi2(MCqq,"chi2_mcqq",-1);
      chi2mc_qq -> Scale(norm);
      maxh = max(maxh, chi2mc_qq -> GetMaximum());
      chi2mc_qq -> SetLineColor(kRed+2);
      chi2mc_qq -> SetLineWidth(2);
      chi2mc_qq -> SetLineStyle(kDashed);
   }

   double ymax = 1.15*maxh;
   chi2cp -> SetMaximum(ymax);

   // box to show rejected region
   TBox* box = new TBox;
   box -> SetFillStyle(3001);
   box -> SetFillColor(kRed-10);
   box -> SetLineColor(kRed-9);
   box -> SetLineWidth(1);

   // Draw
   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1 -> cd();
   gPad -> SetGrid();

   chi2cp -> SetMaximum(ymax);
   chi2cp -> DrawCopy("E1");

   double ymin = gPad->PixeltoY(0);
   box -> DrawBox(chi2_cut, ymin, 200., ymax);
   chi2cp -> DrawCopy("E1,SAME");

   chi2sb -> DrawCopy("SAME HIST");
   chi2_mcS -> DrawCopy("SAME HIST");
   if ( chi2mc_inc ) chi2mc_inc -> DrawCopy("SAME HIST");
   if ( chi2mc_qq ) chi2mc_qq -> DrawCopy("SAME HIST");

   TLegend* leg = new TLegend(0.46,0.64,0.89,0.89);
   leg->SetHeader(Tleg.c_str(),"C");
   leg->AddEntry(chi2cp,"Data: signal region","LE");
   leg->AddEntry(chi2sb,"Data: side-band","F");
   leg->AddEntry(chi2_mcS,"MC signal #phi#eta","L");
   if ( chi2mc_inc ) {
      leg->AddEntry(chi2mc_inc,"J/#Psi inclusive MC (2012)","L");
   }
   if ( chi2mc_qq ) {
      leg->AddEntry(chi2mc_qq,"MC QQ: KKMC at 3080MeV","L");
   }
   leg -> AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   // pdf += ".pdf";
   pdf += "_T.pdf";  // for tight Mkk cut
   c1 -> Print(pdf.c_str());
}

//--------------------------------------------------------------------
void chi2_Pr() {
//--------------------------------------------------------------------
   gROOT -> Reset();
   gStyle -> SetOptStat(0);
   gStyle -> SetLegendFont(42);

   // chi2_SB("3080_rs");
   // chi2_SB("3097J"); // J4
   // chi2_SB("3097");

//    chi2_SB("2900_rs");
//    chi2_SB("3080_2019");
//    chi2_SB("");

}

