// plot weights for K+K- and pi+ pi- pairs:
// -> wts_KK_${date}.pdf; wts_PiPi_${date}.pdf;

#include "ReWeightTrkPid_11.h"

//----------------------------------------------------------------------
void SetHstFace(TH1* hst) {
//----------------------------------------------------------------------
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

//-------------------------------------------------------------------------
void plot_WPi(string fname, int date) {
//-------------------------------------------------------------------------
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }
   cout << " file: " << fname << endl;

   froot->cd("PsipJpsiPhiEta");

   // Psi(2S) -> pi+ pi- J/Psi
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");

//Declaration of leaves types
   Float_t         Mrs;
   Float_t         mdp;
   Float_t         MrsW;
//    vector<float>   Mrec;
   Float_t         Mrb;
   Float_t         Ptp;
   Float_t         Cpls;
   Float_t         Ptm;
   Float_t         Cmns;
   UShort_t        nch;
   UShort_t        Nm;
   UShort_t        dec;
//    Float_t         mcmkk;
//    Float_t         mcmkpet;
//    Float_t         mccosT;

   // Set branch addresses.
   nt1->SetBranchAddress("Mrs",&Mrs);
   nt1->SetBranchAddress("mdp",&mdp);
   nt1->SetBranchAddress("MrsW",&MrsW);
//    nt1->SetBranchAddress("Mrec",&Mrec);
   nt1->SetBranchAddress("Mrb",&Mrb);
   nt1->SetBranchAddress("Ptp",&Ptp);
   nt1->SetBranchAddress("Cpls",&Cpls);
   nt1->SetBranchAddress("Ptm",&Ptm);
   nt1->SetBranchAddress("Cmns",&Cmns);
   nt1->SetBranchAddress("nch",&nch);
   nt1->SetBranchAddress("Nm",&Nm);
   nt1->SetBranchAddress("dec",&dec);
//    nt1->SetBranchAddress("mcmkk",&mcmkk);
//    nt1->SetBranchAddress("mcmkpet",&mcmkpet);
//    nt1->SetBranchAddress("mccosT",&mccosT);


   TH1D* hstMrsW = new TH1D("hstMrsW",
         ";correction factor;Entries/0.001",
         200,0.9,1.1);

   nt1->Draw("MrsW>>hstMrsW","dec==64","goff");

   TH1D* hstMrbW = new TH1D("hstMrbW",
         ";correction factor;Entries/0.001",
         200,0.9,1.1);

   TH1D* hstDelta = new TH1D("hstDelta",
         ";delta MrsW-MrbW;Entries/0.0001",
         200,-0.01,0.01);

   Long64_t nentries = nt1 -> GetEntries();
   for (Long64_t i=0; i<nentries;i++) {
      nt1 -> GetEntry(i);

      if ( dec != 64 ) continue;
      if ( Mrb < 3.0 || Mrb > 3.2 ) continue;

      // correction for pi+,pi- eff:
      double wp = ReWeightTrkPid(date,0,Ptp);
      double wm = ReWeightTrkPid(date,0,Ptm);
      double MrbW = wp*wm;

      hstMrbW -> Fill(MrbW);

      double delta = MrsW - MrbW;
      hstDelta -> Fill(delta);
   }

   TLegend* leg = new TLegend(0.59,0.79,0.89,0.89);
   leg->SetHeader(Form("#pi^{#plus}#pi^{#minus} pairs %i",date),"C");

   TCanvas* c1 = new TCanvas("c1","...",0,0,700,700);
   c1->cd();
   gPad->SetGrid();

   TH1D* hstWpi = hstMrbW;
   SetHstFace(hstWpi);
   hstWpi -> SetLineWidth(2);
//    hstWK -> SetLineColor(kBlack);
//    hstWpi -> GetYaxis() -> SetMaxDigits(3);
   hstWpi -> GetYaxis() -> SetTitleOffset(1.25);
   hstWpi -> Draw("HIST");

   // test
//    hstMrsW -> SetLineWidth(1);
//    hstMrsW -> SetLineColor(kRed);
//    hstMrsW -> Draw("HIST SAME");
//    hstDelta -> Draw("HIST");

   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("wts_PiPi_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_WK(string fname, int date) {
//-------------------------------------------------------------------------
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }
   cout << " file: " << fname << endl;

   froot->cd("PsipJpsiPhiEta");

   // J/Psi -> phi eta
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   // Declaration of leaves types
   Double_t        Mrec;
   Double_t        ch2;
   Double_t        Ptpip;
   Double_t        Cpip;
   Double_t        Ptpim;
   Double_t        Cpim;
   Double_t        Ptkp;
   Double_t        Ckp;
   Double_t        Ptkm;
   Double_t        Ckm;
   Double_t        Eg1;
   Double_t        Cg1;
   Double_t        Eg2;
   Double_t        Cg2;
   Double_t        Ptgg;
   Double_t        Cgg;
   Double_t        Mkk;
   Double_t        Mgg;
   Double_t        dec;
   Double_t        decj;

   // Set branch addresses.
   a4c->SetBranchAddress("Mrec",&Mrec);
   a4c->SetBranchAddress("ch2",&ch2);
   a4c->SetBranchAddress("Ptpip",&Ptpip);
   a4c->SetBranchAddress("Cpip",&Cpip);
   a4c->SetBranchAddress("Ptpim",&Ptpim);
   a4c->SetBranchAddress("Cpim",&Cpim);
   a4c->SetBranchAddress("Ptkp",&Ptkp);
   a4c->SetBranchAddress("Ckp",&Ckp);
   a4c->SetBranchAddress("Ptkm",&Ptkm);
   a4c->SetBranchAddress("Ckm",&Ckm);
   a4c->SetBranchAddress("Eg1",&Eg1);
   a4c->SetBranchAddress("Cg1",&Cg1);
   a4c->SetBranchAddress("Eg2",&Eg2);
   a4c->SetBranchAddress("Cg2",&Cg2);
   a4c->SetBranchAddress("Ptgg",&Ptgg);
   a4c->SetBranchAddress("Cgg",&Cgg);
   a4c->SetBranchAddress("Mkk",&Mkk);
   a4c->SetBranchAddress("Mgg",&Mgg);
   a4c->SetBranchAddress("dec",&dec);
   a4c->SetBranchAddress("decj",&decj);

   // this is cut for a4c-tuple
   auto c_Mrec = [](double Mrec)->bool{return fabs(Mrec-3.097)<0.005;};

   // Cuts for a4c: (see cuts.h)
   auto c_chi2 = [](double ch2)->bool{return ch2<80;}; // std
//    auto c_chi2 = [](double ch2)->bool{return ch2<100;}; // unc: 60, 100

   // Mphi cut: [dL, dU] (see mass_kk_fit.cc)
   const double mk   = 0.493677; // 493.677  +/- 0.016 MeV
   const double dL = 2*mk; // the left cutoff
   const double dU = 1.08; // MUST BE < 1.0835 !!!
   auto c_phi = [dL,dU](double Mkk)->bool{return (Mkk>dL && Mkk<dU);};

   // Meta: central part
   const double meta = 0.547862; //  547.862 +/- 0.017 MeV
   const double seta = 0.008;
   const double weta = 3*seta; // standard
//    const double weta = 4*seta; // uncertainties study: 2x, 4x
   auto c_cpgg = [meta,weta](double Mgg)->bool{return fabs(Mgg-meta)<weta;};

   TH1D* hstWK = new TH1D("hst_WK",
         ";correction factor;Entries/0.001",
         200,0.9,1.1);

   Long64_t nentries = a4c->GetEntries();
   for (Long64_t i=0; i<nentries;i++) {
      a4c->GetEntry(i);
      if ( !c_Mrec(Mrec) ) continue;
      if ( !c_chi2(ch2) ) continue;
      if ( !c_phi(Mkk) ) continue;
      if ( !c_cpgg(Mgg) ) continue;

      // correction for K+,K- eff:
      double wp = ReWeightTrkPid(date,1,Ptkp);
      double wm = ReWeightTrkPid(date,1,Ptkm);
      double w = wp*wm;
      hstWK -> Fill(w);
   }

   TLegend* leg = new TLegend(0.59,0.79,0.89,0.89);
   leg->SetHeader( Form("K^{#plus}K^{#minus} pairs %i",date),"C");

   TCanvas* c1 = new TCanvas("c1","...",0,0,700,700);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hstWK);
   hstWK -> SetLineWidth(2);
//    hstWK -> SetLineColor(kBlack);
   hstWK -> GetYaxis() -> SetMaxDigits(3);
   hstWK -> Draw("HIST");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("wts_KK_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void trk_eff_wts() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetStatFont(42);
   gStyle->SetOptStat(1100);
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.79);
   gStyle->SetStatW(0.3);
   gStyle->SetLegendFont(42);

//    plot_WK("prod-11/mcsig_kkmc_09.root", 2009);
//    plot_WK("prod-11/mcsig_kkmc_12.root", 2012);

//    plot_WPi("prod-11/mcsig_kkmc_09.root", 2009);
//    plot_WPi("prod-11/mcsig_kkmc_12.root", 2012);

}
