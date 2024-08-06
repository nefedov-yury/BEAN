// trk_eff_wts.cc
// plot weights for K+K- and pi+ pi- pairs:
// -> wts_{KK|PiPi}_{YEAR}.pdf

#include "RewTrkPiK.hpp"    // RewTrk functions with HC

// {{{1 helper functions
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

// {{{1 Pion weights, use nt1-tree
//--------------------------------------------------------------------
void plot_WPi(string fname, int date) {
//--------------------------------------------------------------------
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }
   cout << " file: " << fname << endl;

   froot->cd("PsipJpsiPhiEta");

   // Psi(2S) -> pi+ pi- J/Psi
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");
   if ( !nt1 ) {
      cerr << "ERROR in "<< __func__
         << " can not find nt1" << endl;
      exit(EXIT_FAILURE);
   }

   //Declaration of leaves types
   //        those not used here are commented out
   // vector<float> * Mrec = nullptr; // must be !
   Float_t         Mrs;
   Float_t         Ptsp;
   Float_t         Ptsm;
   // Float_t         Mrb;
   // Float_t         Ptp;
   // Float_t         Cpls;
   // Float_t         Ptm;
   // Float_t         Cmns;
   Int_t           dec;
   // Float_t         mcmkk;

   // Set branch addresses.
   // nt1->SetBranchAddress("Mrec",&Mrec);
   nt1->SetBranchAddress("Mrs",&Mrs);
   nt1->SetBranchAddress("Ptsp",&Ptsp);
   nt1->SetBranchAddress("Ptsm",&Ptsm);
   // nt1->SetBranchAddress("Mrb",&Mrb);
   // nt1->SetBranchAddress("Ptp",&Ptp);
   // nt1->SetBranchAddress("Cpls",&Cpls);
   // nt1->SetBranchAddress("Ptm",&Ptm);
   // nt1->SetBranchAddress("Cmns",&Cmns);
   nt1->SetBranchAddress("dec",&dec);
   // nt1->SetBranchAddress("mcmkk",&mcmkk);

   TH1D* hstWpiP = new TH1D("hstWpiP",
         ";correction factor;Entries/0.001",
         200,0.9,1.1);
   TH1D* hstWpiM = new TH1D("hstWpiM",
         ";correction factor;Entries/0.001",
         200,0.9,1.1);
   TH1D* hstWpi = new TH1D("hstWpi",
         ";correction factor;Entries/0.001",
         200,0.9,1.1);

   Long64_t nentries = nt1->GetEntries();

   for ( Long64_t i = 0; i < nentries; ++i ) {
      nt1->GetEntry(i);

      if ( dec == 64 ) {
         if ( Mrs > 1 ) {
            double wp = RewTrkPi( date, Ptsp, +1.);
            double wm = RewTrkPi( date, Ptsm, -1.);
            hstWpiP->Fill(wp);
            hstWpiM->Fill(wm);
            hstWpi->Fill(wp*wm);
         }
      }
   }

   double hmax = hstWpi->GetMaximum();
   hmax = max( hmax,hstWpiP->GetMaximum() );
   hmax = max( hmax,hstWpiM->GetMaximum() );
   // cout << "hmax= " << hmax << endl;

   auto cname = Form("c1_%i",date);
   TCanvas* c1 = new TCanvas(cname,cname,0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   gPad->SetLogy(true);
   hstWpi->SetMinimum(1.);
   hstWpi->SetMaximum(2*hmax);

   SetHstFace(hstWpi);
   hstWpi->SetLineWidth(3);
   hstWpi->GetYaxis()->SetMaxDigits(4);
   hstWpi->GetYaxis()->SetTitleOffset(1.2);
   hstWpi->Draw("HIST");

   hstWpiP->SetLineColor(kRed+1);
   hstWpiM->SetLineColor(kGreen+3);
   hstWpiP->SetLineStyle(kDashed);
   hstWpiM->SetLineStyle(kDashed);
   hstWpiP->SetLineWidth(2);
   hstWpiM->SetLineWidth(2);
   hstWpiP->Draw("SAME HIST");
   hstWpiM->Draw("SAME HIST");

   TLegend* leg = new TLegend(0.60,0.74,0.89,0.89);
   leg->SetHeader( Form("MC inclusive %i",date), "C" );
   leg->AddEntry( hstWpi, "pair of #pi^{#plus}#pi^{#minus}", "L");
   leg->AddEntry( hstWpiP, "#pi^{#plus} only", "L");
   leg->AddEntry( hstWpiM, "#pi^{#minus} only", "L");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "wts_PiPi_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Kaon weights
//--------------------------------------------------------------------
void plot_WK(string fname, int date) {
//--------------------------------------------------------------------
#include "masses.h"
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }
   cout << " file: " << fname << endl;

   froot->cd("PsipJpsiPhiEta");

   // J/Psi -> phi eta
   TTree* a4c = (TTree*)gDirectory->Get("a4c");
   if ( !a4c ) {
      cerr << "ERROR in "<< __func__
         << " can not find a4c" << endl;
      exit(EXIT_FAILURE);
   }

   // Declaration of leaves types
   // #include "a4c_v709.h"
   //        just the ones we use here
   Double_t        Mrec;
   a4c->SetBranchAddress("Mrec",&Mrec);
   Double_t        ch2;
   a4c->SetBranchAddress("ch2",&ch2);
   Double_t        Pkp;
   a4c->SetBranchAddress("Pkp",&Pkp);
   Double_t        Ckp;
   a4c->SetBranchAddress("Ckp",&Ckp);
   Double_t        Pkm;
   a4c->SetBranchAddress("Pkm",&Pkm);
   Double_t        Ckm;
   a4c->SetBranchAddress("Ckm",&Ckm);
   Double_t        Mkk;
   a4c->SetBranchAddress("Mkk",&Mkk);
   Double_t        Mgg;
   a4c->SetBranchAddress("Mgg",&Mgg);
   Double_t        mcmkk;
   a4c->SetBranchAddress("mcmkk",&mcmkk);

   // this is cut for a4c-tuple
   auto c_Mrec = [](double Mrec)->bool{
      return fabs(Mrec-3.097)<0.005;
   };

   // chi^2 cut
   auto c_chi2 = [](double ch2)->bool{
      return ch2<80;
   };

   // Mphi cut: [dL, dU]
   double dL = 2*Mk; // the left cutoff
   double dU = 1.08; // MUST BE < 1.0835 for signal
   auto c_phi = [dL,dU](double Mkk) -> bool{
      return ( Mkk > dL && Mkk < dU);
   };

   // Meta: central part
   const double seta = 0.008;
   const double weta = 3*seta;
   auto c_cpgg = [weta](double Mgg) -> bool{
      return fabs(Mgg-Meta) < weta;
   };

   TH1D* hstWKP = new TH1D("hst_WKP",
         ";correction factor;Entries/0.001",
         200,0.9,1.1);
   TH1D* hstWKM = new TH1D("hst_WKM",
         ";correction factor;Entries/0.001",
         200,0.9,1.1);
   TH1D* hstWK = new TH1D("hst_WK",
         ";correction factor;Entries/0.001",
         200,0.9,1.1);

   Long64_t nentries = a4c->GetEntries();
   for (Long64_t i=0; i<nentries;i++) {
      a4c->GetEntry(i);
      if ( !(mcmkk<1.08) ) continue;
      if ( !c_Mrec(Mrec) ) continue;
      if ( !c_chi2(ch2) ) continue;
      if ( !c_phi(Mkk) ) continue;
      if ( !c_cpgg(Mgg) ) continue;

      // correction for K+,K- eff:
      double Ptkp = Pkp*sqrt(1-Ckp*Ckp); // P -> Pt
      double Ptkm = Pkm*sqrt(1-Ckm*Ckm);
      double wp = RewTrk_K( date, Ptkp, +1.);
      double wm = RewTrk_K( date, Ptkm, -1.);
      hstWKP->Fill(wp);
      hstWKM->Fill(wm);
      hstWK->Fill(wp*wm);
   }

   double hmax = hstWK->GetMaximum();
   hmax = max( hmax,hstWKP->GetMaximum() );
   hmax = max( hmax,hstWKM->GetMaximum() );
   // cout << "hmax= " << hmax << endl;

   auto cname = Form("c1_%i",date);
   TCanvas* c1 = new TCanvas(cname,cname,0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   gPad->SetLogy(true);
   hstWK->SetMinimum(1.);
   hstWK->SetMaximum(2*hmax);
   // hstWK->SetMaximum(1.1*hmax); // liner

   SetHstFace(hstWK);
   hstWK->SetLineWidth(3);
   hstWK->GetYaxis()->SetMaxDigits(4);
   hstWK->GetYaxis()->SetTitleOffset(1.2);
   hstWK->Draw("HIST");

   hstWKP->SetLineColor(kRed+1);
   hstWKM->SetLineColor(kGreen+3);
   hstWKP->SetLineStyle(kDashed);
   hstWKM->SetLineStyle(kDashed);
   hstWKP->SetLineWidth(2);
   hstWKM->SetLineWidth(2);
   hstWKP->Draw("SAME HIST");
   hstWKM->Draw("SAME HIST");

   TLegend* leg = new TLegend(0.60,0.74,0.89,0.89);
   leg->SetHeader( Form("MC signal %i",date), "C" );
   leg->AddEntry( hstWK, "pair of K^{#plus}K^{#minus}", "L");
   leg->AddEntry( hstWKP, "K^{#plus} only", "L");
   leg->AddEntry( hstWKM, "K^{#minus} only", "L");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "wts_KK_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Main
//--------------------------------------------------------------------
void trk_eff_wts() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   // gStyle->SetOptStat(1100);
   // gStyle->SetStatX(0.89);
   // gStyle->SetStatY(0.79);
   // gStyle->SetStatW(0.3);
   // gStyle->SetStatFont(42);
   gStyle->SetLegendFont(42);

   const string Dir = "prod_v709n3/";

   for ( auto date : {2009, 2012, 2021} ) {
      // string mcincfile( Form("mcinc_%02ipsip_all.root",date%100) );
      // plot_WPi(Dir + mcincfile, date);

      // string mcsigfile( Form("mcsig_kkmc_%02i.root",date%100) );
      // plot_WK(Dir + mcsigfile, date);
   }
}
