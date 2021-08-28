// efficiencies for selection of the J/Psi -> phi eta
// as a function of true MC M(K+K-) and
// fit it by a straight line
// -> eff_(sig|bkg)_YYYY.pdf

//----------------------------------------------------------------------
constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

//--------------------------------------------------------------------
double ReWeightEtaEff(int DataPeriod) {
//--------------------------------------------------------------------
// This correction is based on "prod-13eff"
   double W = 0.997; //0.9970 +/- 0.0044
   return W;
}

// copy from trk_eff_fit.cc
//----------------------------------------------------------------------
double ReWeightTrkPid(int DataPeriod, int Kp, double Pt) {
//----------------------------------------------------------------------
// This correction is based on "prod-12/13eff"
// and independent of cos(Theta)
// Kp = 1 for kaons and Kp = 0 for pions

   double W = 1.;

   auto CUBE = [](double x)-> double{return x*x*x;};
   if ( Kp == 1 ) {             // kaons
      Pt = max( 0.1, Pt );
      Pt = min( 1.4, Pt );
      if ( DataPeriod == 2009 ) {
         W = 1.00931 - 0.02363 * Pt;
         if ( Pt < 0.2 ) {
            W = 0.9278;
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
         }
         W = cK12->Eval(Pt);
      }
   } else if ( Kp == 0 ) {      // pions
      Pt = max( 0.05, Pt );
      Pt = min( 0.4, Pt );
      if ( DataPeriod == 2009 ) {
//          W = 0.9878 + CUBE(0.0219/Pt);
         W = 0.9863 + CUBE(0.02234/Pt); // Mar 2020
      } else if ( DataPeriod == 2012 ) {
//          W = 0.9859 + CUBE(0.02974/Pt);
         W = 0.9856 + CUBE(0.02967/Pt); // Mar 2020
      }
   }
   return W;
}

//-------------------------------------------------------------------------
void get_eff(string fname, string pdf="") {
//-------------------------------------------------------------------------
   bool is2009 = (fname.find("_09") != string::npos);
   bool is_sig = (fname.find("mcsig") != string::npos);

   int date = (is2009) ? 2009 : 2012;

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }
   cout << " file: " << fname << endl;

   froot->cd("PsipJpsiPhiEta");

   // initial number of Psi' -> pi+ pi- J/Psi -> phi eta
   TH1D* MCdec = (TH1D*)gROOT->FindObject("mc_dcj0");
   if ( !MCdec ) {
      cout << " can not find mc_dcj0" << endl;
      exit(0);
   }
   double Nini = MCdec->GetBinContent(69);
   cout << " number of generated J/Psi -> phi eta (dec# "
        << MCdec->GetBinCenter(69) << ")= " << Nini << endl;


   // number of pi+pi-J/Psi after Mrs cut [3.092, 3.102]
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");
   TCut c_Mrb = TCut("abs(Mrb-3.097)<0.005");
   cout << " Recoil mass cut: " << c_Mrb.GetTitle() << endl;
   // this is cut for a4c-tuple
   auto c_Mrec = [](double Mrec)->bool{return fabs(Mrec-3.097)<0.005;};

   double Lphi = 0.98, Uphi = 1.08; // histogramms
   int Nbins = 50;
   TCut c_MCmkk("mcmkk<1.08"); // cut for MC generated events
   TH1D* mkk_i = new TH1D("mkk_i","",Nbins,Lphi,Uphi);
   TH1D* mkk_f = new TH1D("mkk_f","",Nbins,Lphi,Uphi);
   mkk_i -> Sumw2(true);
   mkk_f -> Sumw2(true);

   double Nppj = nt1->Draw("mcmkk>>mkk_i",c_Mrb&&c_MCmkk,"goff");
   cout << " number of pi+pi-J/Psi after 'Mrb' cut = "
        << Nppj << endl;

   double nI = mkk_i -> Integral(1,Nbins);
   if ( fabs(nI-Nppj) > 1e-3 ) {
      cout << " DEBUG: mkk_i(Integ)= " << nI << " Nppj= " << Nppj << endl;
      double sc = Nppj/nI;
      cout << "  scale mkk_i = " << sc << endl;
   }

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
   Double_t        mcmkk;

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
   a4c->SetBranchAddress("mcmkk",&mcmkk);

   // Cuts for a4c: (see cuts.h)
   auto c_chi2 = [](double ch2)->bool{return ch2<80;}; // std
//    auto c_chi2 = [](double ch2)->bool{return ch2<60;}; // unc: 60, 100

   // Mphi cut: [dL, dU] (see mass_kk_fit.cc)
   const double mk   = 0.493677; // 493.677  +/- 0.016 MeV
   double dL = 2*mk; // the left cutoff
   double dU = 1.08; // MUST BE < 1.0835 for signal
   auto c_phi = [dL,dU](double Mkk)->bool{return (Mkk>dL && Mkk<dU);};

   // Meta: central part
   const double meta = 0.547862; //  547.862 +/- 0.017 MeV
   const double seta = 0.008;
   const double weta = 3*seta; // standard
//    const double weta = 2*seta; // uncertainties study: 2x, 4x
   auto c_cpgg = [meta,weta](double Mgg)->bool{return fabs(Mgg-meta)<weta;};

   // Meta: side-band
   double shft_eta = 6*seta - weta;
   auto c_sbgg = [meta,weta,shft_eta](double Mgg)->bool{
      return (fabs(Mgg-meta) > weta+shft_eta) &&
             (fabs(Mgg-meta) < 2*weta+shft_eta);
   };

   double Ncp = 0, eNcp = 0;
   double Nsb = 0, eNsb = 0;

   double Weta = ReWeightEtaEff(date);
   Long64_t nentries = a4c -> GetEntries();
   for (Long64_t i=0; i<nentries;i++) {
      a4c->GetEntry(i);
      if ( !(mcmkk<1.08) ) continue;
      if ( !c_Mrec(Mrec) ) continue;
      if ( !c_chi2(ch2) ) continue;
      if ( !c_phi(Mkk) ) continue;

      // correction for K+,K- eff:
      double wp = ReWeightTrkPid(date,1,Ptkp);
      double wm = ReWeightTrkPid(date,1,Ptkm);
      double w = wp*wm;

      // correction for eta eff:
      w *= Weta;

      if ( c_cpgg(Mgg) ) {              // central part
         Ncp  += w;
         eNcp += SQ(w);
         mkk_f -> Fill(mcmkk,w);
      } else if ( c_sbgg(Mgg) ) {       // side-band
         Nsb  += w;
         eNsb += SQ(w);
      }
   }

   // side-band subtraction (OLD)
//    double eff = (Ncp - Nsb) / Nppj;
//    double err = eff*sqrt( (eNcp+eNsb)/SQ(Ncp-Nsb) + 1./Nppj);
//    cout << " Ncp= " << Ncp << " Nsb= " << Nsb << endl;

   // only central part
   double nF = mkk_f -> Integral(1,Nbins);
   cout << " DEBUG: mkk_f(Integ)= " << nF << " Ncp= " << Ncp << endl;

   double eff = Ncp / Nppj;
   double err = eff*sqrt( eNcp/SQ(Ncp) + 1./Nppj);
   printf("\n integral eff(%d) = %.5f +/- %.5f\n",date,eff,err);

   //-----------------------------------------------------------------------
   // Draw results + fit
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1 -> cd();
   gPad -> SetGrid();

   TH1D* heff = (TH1D*)mkk_i -> Clone("heff");
   heff -> Divide(mkk_f,mkk_i,1.,1.,"B");

   string title(";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2};efficiency");
   heff -> SetTitle(title.c_str());
   heff -> GetYaxis() -> SetTitleOffset(1.4);
   heff -> SetLineWidth(2);
   heff -> SetLineColor(kBlack);
   heff -> SetMarkerStyle(20);
//    heff -> SetAxisRange(0.2,0.5,"Y");
   heff -> Draw("EP");

   auto Lf = [](double* x,double* p) -> double {
      return p[0]*(1+p[1]*(x[0]-1.02));
   };
   TF1* ffit = new TF1("ffit", Lf, Lphi, Uphi, 2);
   ffit -> SetParNames("e","k");
   ffit -> SetParameters(0.3, -1.9);
   ffit -> SetLineWidth(1);

   TF1* ff2 = new TF1("ff2", Lf, Lphi, Uphi, 2);
   ff2 -> SetParameters(0.3, -1.9);
   ff2 -> SetLineColor(kBlue);
   ff2 -> SetLineWidth(2);
   ff2 -> SetLineStyle(kDashed);
   if ( is2009 ) {
      ff2 -> FixParameter(1, -1.6);
   } else {
      ff2 -> FixParameter(1, -1.9);
   }

//    double Lfit = int(2*mk*500)/500., Ufit = Uphi - 0.002;
   double Lfit = int(2*mk*500)/500., Ufit = dU;
   heff -> Fit("ffit","EQ","",Lfit,Ufit);
   heff -> Fit("ff2","EQN","",Lfit,Ufit);
   ff2 -> DrawCopy("SAME");
   double chi2 = ffit -> GetChisquare();
   int ndf = ffit -> GetNDF();
   double a = ffit -> GetParameter(0);
   double ea = ffit -> GetParError(0);
   double b = ffit -> GetParameter(1);
   double eb = ffit -> GetParError(1);

   TPaveText* pt = new TPaveText(0.59,0.70,0.89,0.89,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
   string tit;
   if ( is_sig ) {
      tit = "MC signal #phi#eta, ";
   } else {
//       tit = "MC K^{+}K^{-}(not #phi) #eta, ";
      tit = "MC non-#phi K^{+}K^{-}#eta, ";
   }
   tit += to_string(date);
   pt -> AddText( tit.c_str() );
//    pt -> AddText( "#varepsilon = e (1 + k (M-1.02))" );
   pt -> AddText( "eff = e (1 + k (M-1.02))" );
   pt -> AddText( Form("#chi^{2}/ndf = %.1f / %d",chi2,ndf) );
   pt -> AddText( Form("e= %.4f #pm %.4f",a,ea) );
   pt -> AddText( Form("k= %.2f #pm %.2f",b,eb) );
   pt -> Draw();

   printf(" fit: eff(%d)=e*(1+k*(mkk-1.02)):",date);
   printf(" e= %.4f +/- %.4f; k= %.2f +/- %.2f\n\n",a,ea,b,eb);

   double a0 = ff2 -> GetParameter(0);
   double ea0 = ff2 -> GetParError(0);
   double b0 = ff2 -> GetParameter(1);
   double eb0 = ff2 -> GetParError(1);
   if ( fabs(eb0) < 1.e-4 ) {
      printf(" e= %.4f +/- %.4f; k= %.2f (fixed)\n\n",a0,ea0,b0);
   }

   gPad -> RedrawAxis();
   c1 -> Update();
//    if ( !pdf.empty() ) {
//       c1 -> Print(pdf.c_str());
//    }
}

//-------------------------------------------------------------------------
void eff_mc() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle -> SetOptStat(0);
   gStyle -> SetOptFit(0);
   gStyle -> SetStatFont(62);
   gStyle -> SetLegendFont(42);

   get_eff("prod-9/mcsig_kkmc_09.root","eff_sig_2009.pdf");
//    get_eff("prod-9/mcsig_kkmc_12.root","eff_sig_2012.pdf");

//    get_eff("prod-9/mckketa_kkmc_09.root","eff_bkg_2009.pdf");
//    get_eff("prod-9/mckketa_kkmc_12.root","eff_bkg_2012.pdf");
}
