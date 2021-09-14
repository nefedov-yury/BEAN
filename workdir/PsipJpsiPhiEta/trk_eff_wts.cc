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

// copy from trk_eff_fit.cc
//----------------------------------------------------------------------
double ReWeightTrkPid_OLD(int DataPeriod, int Kp, double Pt) {
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

//----------------------------------------------------------------------
double ReWeightTrkPid(int DataPeriod, int Kp, double Pt) {
//----------------------------------------------------------------------
// The corrections are based on production-11(helix corrections for MC).
// They do not dependent of the sign of the particle and cos(Theta) of
// the track. Input parameters are following.
// Kp is the type of the particle: 1 for kaon and 0 for pion.
// Pt is the transverse momentum if the particle.
// The return value is the weight of MC event with such a particle. 

   // parameters
   static const vector<double> K09 {0.991,-0.018};
   static const double K09_first = 0.922;
   static const vector<double> pi09 {0.9872,0.0243};

   static const vector<double> K12 {0.9891,-0.0005,0.0074,0.0111,0.0102};
   static const vector<double> pi12 {0.9851,0.032};

   double W = 1.;

   if ( Kp == 1 ) {             // kaons
      const double Ptmin = 0.1, Ptmax = 1.4;
      Pt = max( Ptmin, Pt );
      Pt = min( Ptmax, Pt );

      int nch = (DataPeriod == 2009) ? 1 : 4;
      double xmin = Ptmin, xmax = Ptmax;
      auto Lchb = [nch,xmin,xmax](double xx, const double* p) {
         if (nch == 0) { return p[0]; }
         // [xmin,xmax] -> [-1,+1]
         double x = (2*xx-xmin-xmax)/(xmax-xmin);
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

      if ( DataPeriod == 2009 ) {
         W = Lchb(Pt,K09.data());
         if ( Pt < 0.2 ) {
            W = K09_first;
         }
      } else if ( DataPeriod == 2012 ) {
         W = Lchb( Pt, K12.data() );
      }
   } else if ( Kp == 0 ) {      // pions
      const double Ptmin = 0.05, Ptmax = 0.4;
      Pt = max( Ptmin, Pt );
      Pt = min( Ptmax, Pt );

      auto CUBE = [](double x)-> double{return x*x*x;};
      if ( DataPeriod == 2009 ) {
         W = pi09[0] + CUBE(pi09[1]/Pt);
      } else if ( DataPeriod == 2012 ) {
         W = pi12[0] + CUBE(pi12[1]/Pt);
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

// OLD
//    plot_WK("archive/prod-9/mcsig_kkmc_09.root", 2009);
//    plot_WK("archive/prod-9/mcsig_kkmc_12.root", 2012);
//
//    plot_WPi("archive/prod-9/mcsig_kkmc_09.root", 2009);
//    plot_WPi("archive/prod-9/mcsig_kkmc_12.root", 2012);
//    plot_WPi("archive/prod-9/mcinc_09psip_all.root", 2009);  // for testing

// NEW
//    plot_WK("prod-11/mcsig_kkmc_09.root", 2009);
//    plot_WK("prod-11/mcsig_kkmc_12.root", 2012);

//    plot_WPi("prod-11/mcsig_kkmc_09.root", 2009);
//    plot_WPi("prod-11/mcsig_kkmc_12.root", 2012);

}
