// efficiencies for selection of the J/Psi -> phi eta
// -> eff_low.pdf

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
void plot_eff(string fname, string fsname="", string pdf="") {
//-------------------------------------------------------------------------
   const double mjpsi  = 3.096916; // 3096.916  +/- 0.011   MeV
   const double meta   = 0.547862; // 547.862   +/- 0.017   MeV
   const double mk     = 0.493677; // 493.677   +/- 0.016   MeV

   bool is2009 = (fname.find("_09") != string::npos);
   bool is_sig = (fname.find("mcsig") != string::npos);
   int date = (is2009) ? 2009 : 2012;

   if ( is2009 || is_sig ) {
      cerr << " this is " << date << " signal:" << is_sig
           << " STOP" << endl;
      exit(0);
   }

   string prod("prod-10/"); // <- MC MUST contain nt1
   fname = prod + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
   }
   cout << " file: " << fname << endl;
   froot->cd("PsipJpsiPhiEta");

   //----------------------------------------------------------------------
   // get initial numbers:
   //----------------------------------------------------------------------

   // initial number of Psi' -> pi+ pi- J/Psi
   TH1D* MCdec = (TH1D*)gROOT->FindObject("mc_dcj0");
   if ( !MCdec ) {
      cout << " can not find mc_dcj0" << endl;
      exit(0);
   }
   double Nini = 0;
   if ( is_sig ) {
      Nini = MCdec->GetBinContent(69);
   } else {
      Nini = MCdec->GetBinContent(1);
   }
   printf(" number of generated J/Psi -> phi eta (dec# %.0f)= %.0f\n",
          MCdec->GetBinCenter(69),Nini);


   // number of pi+pi-J/Psi after Mrs cut [3.092, 3.102]
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");
   if ( !nt1 ) {
      cout << " can not find nt1" << endl;
      exit(0);
   }

   TCut c_Mrb = TCut("abs(Mrb-3.097)<0.005");
   cout << " Recoil mass cut: " << c_Mrb.GetTitle() << endl;

   TCut c_MCmkk("mcmkk<1.18"); // cut for MC generated events ??

   int Nbins = 50;
   TH1D* h_ini = new TH1D("h_ini","",Nbins,-1.,1.);
   TH1D* h_fin = new TH1D("h_fin","",Nbins,-1.,1.);
   TH1D* h_fin2 = new TH1D("h_fin2","",Nbins,-1.,1.);
   h_ini -> Sumw2(true);
   h_fin -> Sumw2(true);
   h_fin2 -> Sumw2(true);

   double Nppj = nt1->Draw("mccosT>>h_ini",c_Mrb+c_MCmkk,"goff");
   cout << " number of pi+pi-J/Psi after 'Mrb' cut = "
        << Nppj << endl;

   // 2-D plot of efficiency
   int Nkk = 40, Nke = 40;
   TH2D* h2i = new TH2D("h2i","",Nkk,0.98,1.18,Nke,1.55,2.55);
   TH2D* h2f = new TH2D("h2f","",Nkk,0.98,1.18,Nke,1.55,2.55);
   h2i -> Sumw2(true);
   h2f -> Sumw2(true);

   nt1->Draw("mcmkpet:mcmkk>>h2i",c_Mrb,"goff");

   //----------------------------------------------------------------------
   // get final numbers:
   //----------------------------------------------------------------------

   // J/Psi -> phi eta
   TTree* a4c = (TTree*)gDirectory->Get("a4c");
   if ( !a4c ) {
      cout << " can not find a4c" << endl;
      exit(0);
   }
#include "p10a4cMC.h"

   // cut on recoil mass: [3.092, 3.102]
   auto c_Mrec = [](double Mrec)->bool{return fabs(Mrec-3.097)<0.005;};
   // chi^2 cut:
   auto c_chi2 = [](double ch2)->bool{return ch2<80;};
   // Meta: central part
   auto c_cpgg = [meta](double Mgg)->bool{return fabs(Mgg-meta)<0.024;};
   // Mkk cut
//    auto c_mkk = [mk](double Mkk)->bool{return (Mkk>2*mk&&Mkk<1.1);};
   auto c_mkk = [mk](double Mkk)->bool{return (Mkk>2*mk&&Mkk<1.18);};

   double Weta = ReWeightEtaEff(date);
   Long64_t nentries = a4c -> GetEntries();
   for (Long64_t i=0; i<nentries;i++) {
      a4c->GetEntry(i);

      // selection:
//       if ( !(mcmkk<1.1) ) continue;
      if ( !c_Mrec(Mrec) ) continue;
      if ( !c_chi2(ch2) ) continue;
      if ( !c_cpgg(Mgg) ) continue;
      if ( !c_mkk(Mkk) ) continue;

      // correction for K+,K- eff:
      double wp = ReWeightTrkPid(date,1,Ptkp);
      double wm = ReWeightTrkPid(date,1,Ptkm);
      double w = wp*wm;

      // correction for eta eff:
      w *= Weta;
      h_fin -> Fill(mccosT,w);

      TVector3 VPj; // J/Psi in Lab system
      VPj.SetMagThetaPhi(Pj,acos(Cj),phij);
      TLorentzVector LVPj;
      LVPj.SetVectM(VPj,mjpsi);
      TVector3 boostJ = LVPj.BoostVector();

      TVector3 VPkk; // K+K- system in Lab
      VPkk.SetMagThetaPhi(Pkk,acos(Ckk),phikk);
      TLorentzVector LVPkk;
      LVPkk.SetVectM(VPkk,Mkk);
      TVector3 boostKK = LVPkk.BoostVector();
      LVPkk.Boost(-boostJ);
      VPkk = LVPkk.Vect(); // KK in J/Psi rest system

      TVector3 VPkp; // K+ in Lab
      VPkp.SetMagThetaPhi(Pkp,acos(Ckp),phikp);
      TLorentzVector LVPkp;
      LVPkp.SetVectM(VPkp,mk);
      LVPkp.Boost(-boostKK);
      VPkp = LVPkp.Vect(); // K+ in KK rest system

      double cos_T = VPkp.Dot(VPkk)/(VPkp.Mag()*VPkk.Mag());
      h_fin2 -> Fill(cos_T,w);

      // 2-D hists
      h2f -> Fill(mcmkk,mcmkpet,w);
   }

   double nF = 0., nF_err = 0.;
   nF = h_fin  -> IntegralAndError(1,Nbins,nF_err);
   double eff = nF / Nppj;
   double err = eff*sqrt(SQ(nF_err/nF) + 1/Nppj);
   printf(" integral eff(%d) = %.5f +/- %.5f\n",date,eff,err);

   //----------------------------------------------------------------------
   // efficiencies
   //----------------------------------------------------------------------
   TH1D* heff = (TH1D*)h_ini -> Clone("heff");
   heff -> Divide(h_fin,h_ini,1.,1.,"B");

   TH1D* heff2 = (TH1D*)h_ini -> Clone("heff2");
   heff2 -> Divide(h_fin2,h_ini,1.,1.,"B");

   TH2D* h2eff = (TH2D*)h2i -> Clone("h2eff");
   h2eff -> Divide(h2f,h2i,1.,1.,"B");

   //----------------------------------------------------------------------
   // Draw results
   //----------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1 -> cd();
   gPad -> SetGrid();

   /*
   // 1-D
   string title(";cosT;efficiency");
   heff -> SetTitle(title.c_str());
   heff -> GetYaxis() -> SetTitleOffset(1.4);
   heff -> SetLineWidth(2);
   heff -> SetLineColor(kBlack);
   heff -> SetMarkerStyle(20);
//    heff -> SetAxisRange(0.,1.,"Y");
   heff -> Draw("EP");

   heff2 -> SetLineColor(kRed);
   heff2 -> Draw("HIST,SAME");
   */

   // 2-D
//    h2f -> Draw("Lego"); // Surf3
   h2eff -> Draw("COLZ");

   // save 2d-eff
   if ( !fsname.empty() ) {
      fsname = prod + fsname;
      TFile* fsr = TFile::Open(fsname.c_str(),"NEW");
      h2eff->Write();
      fsr->Close();
   }

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void eff_low() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle -> SetOptStat(0);
//    gStyle -> SetOptFit(0);
   gStyle -> SetStatFont(62);
   gStyle -> SetLegendFont(42);

//    plot_eff("mckketa_kkmc_12.root","eff_low.pdf");
   plot_eff("mckketa_kkmc_12.root","eff_low.root","eff_low.pdf");

}
