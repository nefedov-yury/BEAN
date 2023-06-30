// selection efficiency of the J/Psi -> phi eta as a function
// of true M(K+K-). We fit this efficiency with constant and
// linear dependences on M(K+K-)
// -> eff_sig_XXX.pdf

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
constexpr double SQ(double x) {
//--------------------------------------------------------------------
   return x*x;
}

// {{{1 get_eff()
//--------------------------------------------------------------------
void get_eff(string fname, string title, string pdf="") {
//--------------------------------------------------------------------
#include "cuts.h"

   bool is_phsp = (fname.find("PHSP") != string::npos); // phase space
   bool is_sig = !is_phsp;

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cout << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }
   string sepline(70,'='); // separation line
   printf("%s\n",sepline.c_str());
   printf(" %s\n", title.c_str());
   printf(" file: %s\n",fname.c_str());
   froot->cd("SelectKKgg");

   //-----------------------------------------------------------------
   // get initial numbers:
   //-----------------------------------------------------------------
   TH1D* mc_Mkk_ini = (TH1D*)gROOT->FindObject("mc_Mkk_ini");
   if( !mc_Mkk_ini ) {
      cout << " can not find mc_Mkk_ini histo" << endl;
      exit(EXIT_FAILURE);
   }
   // Note: Xisr > 0.9 cut is applied in SelectKKgg.cxx
   double Nini = mc_Mkk_ini->Integral(1,100); // Mkk<1.08GeV
   printf(" number of generated 'e+e- -> phi eta'"
         " for Xisr>0.9 is %.0f\n",Nini);

   // 100bins [0.98,108]
   // int Nbins = 100;
   // TH1D* mkk_i = mc_Mkk_ini;

   // re-bin to 50 bins
   int Nbins = 50;
   TH1D* mkk_i = (TH1D*)mc_Mkk_ini->Rebin(2,"mkk_i");

   // re-bin to 25 bins
   // int Nbins = 25;
   // TH1D* mkk_i = (TH1D*)mc_Mkk_ini->Rebin(4,"mkk_i");

   if ( mkk_i->GetNbinsX() != Nbins ) {
      cout << " something wrong with mkk_i: Nbins="
           << mkk_i->GetNbinsX() << endl;
      exit(EXIT_FAILURE);
   }
   mkk_i->Sumw2(true);

   //-----------------------------------------------------------------
   // get final numbers:
   //-----------------------------------------------------------------
   TTree* a4c = (TTree*)gDirectory->Get("a4c");
   if ( !a4c ) {
      cout << " can not find a4c" << endl;
      exit(EXIT_FAILURE);
   }
#include "a4c_prod3.h"

   double Lphi = 0.98, Uphi = 1.08; // histograms
   TH1D* mkk_cp = new TH1D("mkk_cp","central part",Nbins,Lphi,Uphi);
   TH1D* mkk_sb = new TH1D("mkk_sb","sideband",Nbins,Lphi,Uphi);
   mkk_cp->Sumw2(true);
   mkk_sb->Sumw2(true);

   double Ncp = 0, eNcp = 0;
   double Nsb = 0, eNsb = 0;

   Long64_t nentries = a4c->GetEntries();
   for ( Long64_t i = 0; i < nentries; ++i ) {
      a4c->GetEntry(i);
      if ( !f_xisr(xisr) ) continue;
      if ( !f_MCmkk(mcmkk) ) continue;
      if ( !f_chi2(ch2) ) continue;
      if ( !f_phi(Mkk) ) continue;

      double w = 1; // no weights
      if ( f_cpgg(Mgg) ) {              // central part
         Ncp  += w;
         eNcp += SQ(w);
         mkk_cp->Fill(mcmkk,w);
      } else if ( f_sbgg(Mgg) ) {       // side-band
         Nsb  += w;
         eNsb += SQ(w);
         mkk_sb->Fill(mcmkk,w);
      }
   }

   double nCP = mkk_cp->Integral(1,Nbins);
   double nSB = mkk_sb->Integral(1,Nbins);
   if ( fabs(nCP-Ncp) > 1e-6 || fabs(nSB-Nsb) > 1e-6 ) {
      printf(" DEBUG: CP(Integ)= %f, Ncp= %f\n",nCP,Ncp);
      printf(" DEBUG: SB(Integ)= %f, Nsb= %f\n",nSB,Nsb);
   }

   // Efficiency with side-band subtraction
   double eff = (Ncp - Nsb) / Nini;
   double err = eff*sqrt( (eNcp+eNsb)/SQ(Ncp-Nsb) + 1./Nini );
   TH1D* heff = (TH1D*)mkk_i->Clone("heff");
   heff->Add(mkk_cp,mkk_sb,1.,-1.);    // side-band subtraction
   // heff->Divide(heff,mkk_i,1.,1.,"B"); // efficiency
   heff->Divide(heff,mkk_i,1.,1.); // efficiency

   // Efficiency by central part only
   // double eff = Ncp / Nini;
   // double err = eff*sqrt( eNcp/SQ(Ncp) + 1./Nini );
   // TH1D* heff = (TH1D*)mkk_i->Clone("heff");
   // heff->Divide(mkk_cp,mkk_i,1.,1.,"B"); // efficiency

   printf("\n integral efficiency: %.5f +/- %.5f\n",eff,err);

   //-----------------------------------------------------------------
   // Draw results + fit
   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,1000,800);
   c1->cd();
   gPad->SetGrid();

   heff->SetTitle( ";M^{ inv}_{ K^{#plus}K^{#minus}}, GeV/c^{2};"
         "efficiency"
         );
   heff->GetXaxis()->SetTitleOffset(1.1);
   heff->GetYaxis()->SetTitleOffset(1.3);
   heff->SetLineWidth(2);
   heff->SetLineColor(kBlack);
   heff->SetMarkerStyle(20);
   heff->SetAxisRange(0.05,0.25,"Y");
   heff->GetYaxis()->SetNdivisions(504);
   heff->Draw("EP");

   // left and right fit boundaries:
   double bin_width = heff->GetBinWidth(1);
   double Lfit = (int(2*Mk/bin_width)+1)*bin_width;
   // double Ufit = Uphi - bin_width;
   double Ufit = Uphi;

   auto Lf = [](double* x,double* p) -> double {
      return p[0]*(1+p[1]*(x[0]-1.02));
   };
   TF1* ff = new TF1("ff", Lf, Lfit, Ufit, 2);
   ff->SetParNames("e","k");
   ff->SetParameters(0.2, 0.);
   ff->FixParameter(1, 0.); // fit by constant
   ff->SetLineColor(kBlue);
   ff->SetLineWidth(1);

   heff->Fit("ff","EQ","",Lfit,Ufit);

   TF1* ff2 = new TF1("ff2", Lf, Lfit, Ufit, 2);
   ff2->SetParameters(0.2, -1.);
   ff2->SetLineColor(kRed);
   ff2->SetLineWidth(1);
   ff2->SetLineStyle(kDashed);

   heff->Fit("ff2","EQN","",Lfit,Ufit);
   ff2->DrawCopy("SAME");

   double chi2 = ff->GetChisquare();
   int ndf = ff->GetNDF();
   double a = ff->GetParameter(0);
   double ea = ff->GetParError(0);

   double Ch2 = ff2->GetChisquare();
   int Ndf = ff2->GetNDF();
   double a2 = ff2->GetParameter(0);
   double ea2 = ff2->GetParError(0);
   double b2 = ff2->GetParameter(1);
   double eb2 = ff2->GetParError(1);

   TPaveText* pt = new TPaveText(0.12,0.77,0.50,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   string header(title);
   if ( is_sig ) {
      header = "MC signal #phi#eta, " + header;
   } else {
      header = "MC non-#phi K^{+}K^{-}#eta, " + header;
   }
   pt->AddText( header.c_str() );
   pt->AddText( Form("#chi^{2}/ndf = %.1f / %d",chi2,ndf) );
   pt->AddText( Form("efficiency = %.4f #pm %.4f",a,ea) );
   pt->Draw();

   TPaveText* pt2 = new TPaveText(0.59,0.77,0.89,0.89,"NDC");
   pt2->SetTextAlign(12);
   pt2->SetTextFont(42);
   pt2->AddText("#color[2]{ eff = e(1 #plus k (M #minus 1.02))}");
   pt2->AddText(Form("#color[2]{#chi^{2}/ndf = %.1f / %d}",Ch2,Ndf));
   pt2->AddText(Form("#color[2]{e = %.4f #pm %.4f}",a2,ea2));
   pt2->AddText(Form("#color[2]{k = %.2f #pm %.2f}",b2,eb2));
   pt2->Draw();

   // printf("\n%s\n",sepline.c_str());
   printf("\n constant fit:\n");
   printf("     chi2/ndf = %.1f / %d\n", chi2,ndf);
   printf("     eff= %.5f +/- %.5f\n", a,ea);
   printf("\n fitting by e*(1+k*(mkk-1.02)):\n");
   printf("     chi2/ndf = %.1f / %d\n", Ch2,Ndf);
   printf("     e= %.4f +/- %.4f; k= %.2f +/- %.2f\n",a2,ea2,b2,eb2);
   printf("%s\n",sepline.c_str());

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 Main
//--------------------------------------------------------------------
void eff_MC() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetLegendFont(42);

   // get_eff("Ntpls/ntpl_mcgpj_J4.root",
           // "2018: 3097MeV", "eff_sig_3097J.pdf");
   // get_eff("Ntpls/ntpl_mcgpj_3080_rs.root",
           // "R-scan: 3080MeV ", "eff_sig_3080rc.pdf");

   // get_eff("Ntpls/ntpl_mcgpj_3097.root",
           // "3097MeV J/#Psi-scan", "eff_sig_3097.pdf");
   // get_eff("Ntpls/ntpl_mcgpj_3080_2019.root",
           // "3080 (2019)", "eff_sig_3080_2019.pdf");
}
