// plot weights for K+K- and pi+ pi- pairs:
// -> wts_KK_${date}.pdf; wts_PiPi_${date}.pdf;

//----------------------------------------------------------------------
void SetHstFace(TH1* hst) {
//----------------------------------------------------------------------
   TAxis* X = hst->GetXaxis();
   if ( X ) {
      X->SetLabelFont(62);
      X->SetLabelSize(0.04);
      X->SetTitleFont(62);
      X->SetTitleSize(0.045);
      X->SetTitleOffset(1.0);
   }
   TAxis* Y = hst->GetYaxis();
   if ( Y ) {
      Y->SetLabelFont(62);
      Y->SetLabelSize(0.04);
      Y->SetTitleFont(62);
      Y->SetTitleSize(0.04);
      Y->SetTitleOffset(1.3);
   }
   TAxis* Z = hst->GetZaxis();
   if ( Z ) {
      Z->SetLabelFont(62);
      Z->SetLabelSize(0.04);
      Z->SetTitleFont(62);
      Z->SetTitleSize(0.04);
   }
}

// copy from trk_eff_fit.cc
//----------------------------------------------------------------------
double ReWeightTrkPid(int DataPeriod, int Kp, double Pt) {
//----------------------------------------------------------------------
// This correction is based on "prod-12eff"
// and independent of cos(Theta)
// Kp = 1 for kaons and Kp = 0 for pions

   double W = 1.;

   auto CUBE = [](double x)-> double{return x*x*x;};
   if ( Kp == 1 ) {             // kaons
      Pt = max( 0.1, Pt );
      Pt = min( 1.4, Pt );
      if ( DataPeriod == 2009 ) {
         W = 1.00931 - 0.02363 * Pt;
         // OLD: prod-te6; -te10
//          W = 1.00703 - 0.01977 * Pt;
         if ( Pt < 0.2 ) {
            W = 0.9278;
//             W = 0.9311; // OLD: prod-te6; -te10
         }
      } else if ( DataPeriod == 2012 ) {
         static TF1* cK12 = nullptr;
         if ( !cK12 ) {
            int nch = 4;
            auto Lchb = [nch](const double* xx, const double* p) -> double {
               if (nch == 0) { return p[0]; }
               double x = xx[0];
               double sum = p[0] + x*p[1];
               if (nch == 1) { return sum; }
               double T0 = 1, T1 = x;
               for ( int i = 2; i <= nch; ++i ) {
                  double Tmp = 2*x*T1 - T0;
                  sum += p[i]*Tmp;
                  T0 = T1;
                  T1 = Tmp;
               }
               return sum;
            };
            cK12 = new TF1("cK12", Lchb, 0.1, 1.4,nch+1);
            cK12->SetParameters(1.82144,-1.41435,0.83606,-0.32437,0.05736);
            // OLD: prod-te6; -te10
//             cK12->SetParameters(2.04378,-1.78748,1.05229,-0.40293,0.07065);
         }
         W = cK12->Eval(Pt);
      }
   } else if ( Kp == 0 ) {      // pions
      Pt = max( 0.05, Pt );
      Pt = min( 0.4, Pt );
      if ( DataPeriod == 2009 ) {
         W = 0.9878 + CUBE(0.0219/Pt);
         // OLD: prod-te6; -te10
//          W = 0.9871 + CUBE(0.0224/Pt);
      } else if ( DataPeriod == 2012 ) {
         W = 0.9859 + CUBE(0.02974/Pt);
         // OLD: prod-te6; -te10
//          W = 0.9843 + CUBE(0.03015/Pt);
      }
   }
   return W;
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

   TH1D* hstWpi = new TH1D("hstWpi",";correction factor",200,0.9,1.1);

   nt1->Draw("MrsW>>hstWpi","dec==64","goff");

   TLegend* leg = new TLegend(0.49,0.79,0.89,0.89);
//    leg->SetHeader( (string("#bf{")+to_string(date)+"}").c_str(),"C");
   leg->SetHeader( Form("#pi^{+}#pi^{-} pair: %i",date),"C");

   TCanvas* c1 = new TCanvas("c1","...",0,0,700,700);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hstWpi);
   hstWpi -> SetLineWidth(2);
//    hstWK -> SetLineColor(kBlack);
   hstWpi -> GetYaxis() -> SetMaxDigits(3);
   hstWpi -> Draw("HIST");
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

   TH1D* hstWK = new TH1D("hst_WK",";correction factor",200,0.9,1.1);

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

   TLegend* leg = new TLegend(0.49,0.79,0.89,0.89);
//    leg->SetHeader( (string("#bf{")+to_string(date)+"}").c_str(),"C");
   leg->SetHeader( Form("K^{+}K^{-} pair: %i",date),"C");

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
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(62);
   gStyle->SetStatFont(62);


   plot_WK("prod-9/mcsig_kkmc_12.root", 2012);
//    plot_WK("prod-9/mcsig_kkmc_09.root", 2009);

//    plot_WPi("prod-9/mcinc_12psip_all.root", 2012);

//    plot_WPi("prod-9/mcsig_kkmc_12.root", 2012);
//    plot_WPi("prod-9/mcsig_kkmc_09.root", 2009);
}
