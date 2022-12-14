// plot eff, purity and signal/background ration as function of chi^2
// cut estimated by Monte Carlo for inclusive decays of J/Psi
// -> chi2_eff_pur.pdf

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
void eff_pur(string pdf) {
//--------------------------------------------------------------------
   const char* fname = "Ntpls/ntpl_jpsi_incl.root";
   TFile* froot = TFile::Open(fname,"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("SelectKKgg");
   TTree *a4c = (TTree*)gDirectory->Get("a4c");

#include "cuts.h"
   TCut c_here = c_xisr + c_MCmkk; // X_isr>0.9 && mc_Mkk<1.08
   c_here += c_cpgg;               // central part of Mgg
   c_here += c_phi;                // Mkk in [2*Mk, 1.08GeV]

   TCut c_good("dec==68");
   TCut c_bad("dec!=68");

   TH1D* ch2_good = new TH1D("ch2_good",
         ";#chi^{2};Entries/2", 100, 0.,200.);
   TH1D* ch2_bad  = new TH1D("ch2_bad",
         ";#chi^{2};Entries/2", 100, 0.,200.);

   a4c->Draw("ch2>>ch2_good",c_here+c_good,"goff");
   a4c->Draw("ch2>>ch2_bad", c_here+c_bad, "goff");

   double *Ig = ch2_good->GetIntegral();
   double *Ib = ch2_bad->GetIntegral();
   int Nbins = 100;
   double Ng = Ig[Nbins+1]; // number of entries
   double Nb = Ib[Nbins+1];

   TH1D* eff = new TH1D("eff",
         ";cut on #chi^{2}",//;efficiency",
         Nbins,0.,200.);
   TH1D* pur = new TH1D("pur",
         ";cut on #chi^{2}",//;purity" Sig/(Sig+Bkg)
         Nbins,0.,200.);
   TH1D* rat = new TH1D("rat",
         ";cut on #chi^{2}",// Sig/sqrt(Sig+Bkg)
         Nbins,0.,200.);

   for(int i=1; i <= Nbins; i++) {
      double bc = ch2_good -> GetBinCenter(i);
      eff -> Fill( bc, 1-Ig[i] );
      pur -> Fill( bc, Ng*Ig[i]/(Ng*Ig[i]+Nb*Ib[i]) );
      rat -> Fill( bc, Ng*Ig[i]/sqrt(Ng*Ig[i]+Nb*Ib[i]) );
   }
   SetHstFace(eff);
   eff->SetLineWidth(2);
   SetHstFace(pur);
   pur->SetLineWidth(2);
   SetHstFace(rat);
   rat->SetLineWidth(2);

   double pur80 = pur->GetBinContent(81); // purity for ch2 < 80
   printf(" purity for ch2<80 is %.2f%%\n", pur80*100);

   // Draw
   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->Divide(1,2);

//    c1->cd(1);
//    gPad->SetGrid();
//    ch2_good->Draw("hist");
//    pur->SetAxisRange(0.95,1.0,"Y");
//    pur->Draw("hist");

   c1->cd(1);
   gPad->SetGrid();
   rat->SetAxisRange(0.,35.,"Y");
   rat->Draw("hist L");

   auto pt1 = new TPaveText(0.12,0.79,0.49,0.89,"NDC");
   pt1 -> SetTextAlign(12);
   pt1 -> SetTextFont(42);
   pt1 -> AddText( "Signal #lower[0.15]{#scale[1.6]{/}}"
         "#sqrt{Signal+Background}");
   pt1 -> Draw();


   c1->cd(2);
   gPad->SetGrid();
//    ch2_bad->Draw("hist");
//    eff->SetAxisRange(0.75,1.05,"Y");
   eff->Draw("hist L");

   auto pt2 = new TPaveText(0.32,0.79,0.89,0.89,"NDC");
   pt2 -> SetTextAlign(12);
   pt2 -> SetTextFont(42);
   pt2 -> AddText("Fraction of events lost due to #chi^{2} cut");
   pt2 -> Draw();

   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

//--------------------------------------------------------------------
void jpsi_incl() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   string pdf("chi2_eff_pur.pdf");
   eff_pur(pdf);
}
