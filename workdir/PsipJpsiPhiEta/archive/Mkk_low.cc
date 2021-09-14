// study low M(K+K-) region: Mkk < 1.1 GeV/c^2

//----------------------------------------------------------------------
constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

// cosT: T is the angle between the K+ direction in K+K- rest frame
// and the prior direction of the K+K- system in J/Psi rest frame
//-------------------------------------------------------------------
struct Miaa {
   double mkk;     // inv mass M(K+K-)
   double cosT;
   double w;       // weight
};
//-------------------------------------------------------------------

// return values of Legendre polynomials in x for indexes 0..kmax
//-------------------------------------------------------------------
vector<double> Pleg(int kmax, double x) {
//-------------------------------------------------------------------
   vector<double> Pl(max(kmax+1,1),1.);
   if ( kmax > 0 ) {
      Pl[1] = x;
      for ( int k = 2; k <= kmax; ++k ) {
         Pl[k] = ((2*k-1)*x*Pl[k-1]-(k-1)*Pl[k-2])/double(k);
      }
   }
   return Pl;
}

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

//-------------------------------------------------------------------------
TH2D* get_loweff() {
//-------------------------------------------------------------------------
   string fname("prod-10/eff_low.root");
//    fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   TH2D* eff = (TH2D*)gROOT->FindObject("h2eff");
   if ( !eff ) {
      cout << " can not find h2eff" << endl;
      exit(0);
   }

   return eff;
}

//-------------------------------------------------------------------------
double Eff_weight(double mkk, double mkpeta) {
//-------------------------------------------------------------------------
   const double mk     = 0.493677; // 493.677   +/- 0.016   MeV
   static TH2D* eff = nullptr;

   if ( !eff ) {
      eff = get_loweff();
   }

   int ibin = eff->FindFixBin(mkk,mkpeta);
   if ( ibin < 0 ) {
      cout << "ERROR: ibin= " << ibin
           << " mkk= " << mkk << " mkpeta= " << mkpeta << endl;
      exit(0);
   }
   double Eff = eff->GetBinContent(ibin);

   // phase-spase
   double phsp = sqrt(1.-SQ(2*mk/mkk));

   // weight = 1 / (Eff*phsp)
//    double weight = Eff*phsp;
   double weight = phsp;
//    double weight = 1;
   if ( weight > 1e-6 ) {
      weight = 1. / weight;
   } else {
      cout << "ERROR: weight= " << weight
           << " mkk= " << mkk << " mkpeta= " << mkpeta << endl;
//       exit(0);
   }

   return weight;
}

//-------------------------------------------------------------------------
vector<Miaa> get_cosT(string fname) {
//-------------------------------------------------------------------------
   const double mjpsi  = 3.096916; // 3096.916  +/- 0.011   MeV
   const double meta   = 0.547862; // 547.862   +/- 0.017   MeV
   const double mk     = 0.493677; // 493.677   +/- 0.016   MeV
//    const double mphi   = 1.019461; //1019.461  +/- 0.019 MeV

   bool isMC = (fname.find("mcsig") != string::npos);
   if ( isMC ) {
      cout << " THIS IS MONTE CARLO EVENTS" << endl;
   }

   // name of folder with root files
   static string dir("prod-10/"); // MUST BE >= 10
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory -> Get("a4c");
#include "p10a4c.h"

   // cut on recoil mass: [3.092, 3.102]
   auto c_Mrec = [](double Mrec)->bool{return fabs(Mrec-3.097)<0.005;};
   // chi^2 cut:
   auto c_chi2 = [](double ch2)->bool{return ch2<80;};
   // Meta: central part
   auto c_cpgg = [meta](double Mgg)->bool{return fabs(Mgg-meta)<0.024;};
   // Mkk cut
   auto c_mkk = [mk](double Mkk)->bool{return (Mkk>2*mk&&Mkk<1.15);};

   vector<Miaa> cosT;
   cosT.reserve(4000); // 3027

   Long64_t nentries = a4c->GetEntries();
   for ( Long64_t i=0; i < nentries; ++i ) {
      a4c -> GetEntry(i);

      // selection:
      if ( !c_Mrec(Mrec) ) continue;
      if ( !c_chi2(ch2) ) continue;
      if ( !c_cpgg(Mgg) ) continue;
      if ( !c_mkk(Mkk) ) continue;

      TVector3 VPj; // J/Psi in Lab system
      VPj.SetMagThetaPhi(Pj,acos(Cj),phij);
      TLorentzVector LVPj;
      LVPj.SetVectM(VPj,mjpsi);
      TVector3 boostJ = LVPj.BoostVector();

      // check
      TLorentzVector LVtmp = LVPj;
      LVtmp.Boost(-boostJ);
      if( fabs(LVtmp.E() - mjpsi) > 1e-7 ) {
         cout << " ERROR: boostJ: delta= " << LVtmp.E() - mjpsi << endl;
         exit(0);
      }

      TVector3 VPkk; // K+K- system in Lab
      VPkk.SetMagThetaPhi(Pkk,acos(Ckk),phikk);
      TLorentzVector LVPkk;
      LVPkk.SetVectM(VPkk,Mkk);
      TVector3 boostKK = LVPkk.BoostVector();

      // check
      LVtmp = LVPkk;
      LVtmp.Boost(-boostKK);
      if( fabs(LVtmp.E() - Mkk) > 1e-7 ) {
         cout << " ERROR: boostKK: delta= " << LVtmp.E() - Mkk << endl;
         exit(0);
      }

      LVPkk.Boost(-boostJ);
      VPkk = LVPkk.Vect(); // KK in J/Psi rest system

      TVector3 VPkp; // K+ in Lab
      VPkp.SetMagThetaPhi(Pkp,acos(Ckp),phikp);
      TLorentzVector LVPkp;
      LVPkp.SetVectM(VPkp,mk);
      LVPkp.Boost(-boostKK);
      VPkp = LVPkp.Vect(); // K+ in KK rest system

      double cos_T = VPkp.Dot(VPkk)/(VPkp.Mag()*VPkk.Mag());
      double m_kpeta = sqrt(M2kpeta);
      double w = Eff_weight( Mkk, m_kpeta );
      cosT.push_back({Mkk,cos_T,w});

      // check
      double cos_t = VPkk.Dot(VPkp)/(VPkp.Mag()*VPkk.Mag());
      if( fabs(cos_t - cos_T) > 1e-7 ) {
         cout << " ERROR: cosT: delta= " << cos_t - cos_T << endl;
         exit(0);
      }
   }

   cerr << " select " << cosT.size() << " events " << endl;
   return cosT;
}

//-------------------------------------------------------------------------
TH1D* hst_cosT(const vector<Miaa>& mia, string hname) {
//-------------------------------------------------------------------------

   const int Nbins = 50;
   string title = string("angle of K^{+} in K^{+}K^{-} RF"
                         " wrt direction of K^{+}K^{-} in J/Psi RF");
   title += ";cos #Theta(K^{+} in KK)";
   title += ";Events/0.04";
   auto hst = new TH1D(hname.c_str(),title.c_str(),Nbins,-1.,1.);

   hst -> Sumw2(true);
   for ( const auto& ct : mia ) {
      hst -> Fill(ct.cosT,ct.w);
   }
   return hst;
}

//-------------------------------------------------------------------------
void plot_cosT(const vector<Miaa>& mia, string pdf) {
//-------------------------------------------------------------------------

   auto hst = hst_cosT(mia,"cosT");

   TCanvas* c1 = new TCanvas("c1","Mkk low",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst);
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
   hst -> SetMarkerStyle(20);
   hst -> GetYaxis() -> SetTitleOffset(1.25);
   hst -> Draw("EP");

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void get_Yi(const vector<Miaa>& mia, vector<double> Y[]) {
//-------------------------------------------------------------------------
   const double tmp = 1./sqrt(4*M_PI);
   const int Ny = 5;

//    const double p2y[] = {tmp, sqrt(3)*tmp, sqrt(5)*tmp };
   vector<double> p2y(Ny,tmp);
   for ( int i = 1; i < Ny; ++i ) {
      p2y[i] *= sqrt(double(2*i+1));
   }

   int n = mia.size();
   for ( int i = 0; i < Ny; ++i ) {
      Y[i].reserve(n);
   }

   for ( const auto& mc : mia ) {
//       vector<double> Pl = Pleg(2,mc.cosT); // 0,1,2
      vector<double> Pl = Pleg(Ny-1,mc.cosT);
      for ( int i = 0; i < Ny; ++i ) {
         Y[i].push_back(mc.w*p2y[i]*Pl[i]);
      }
   }
}

//-------------------------------------------------------------------------
void plot_Yi(const vector<Miaa>& mia,const vector<double> Y[],string pdf) {
//-------------------------------------------------------------------------
   const int Ny = 5;
   const int Nbins = 85;
   string title = ";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2}";
   title += ";Entries / 2 MeV/c^{2}";
   TH1D* hst[10];
   for ( int i = 0; i < Ny; ++i ) {
      string hname = "Y" + to_string(i);
      hst[i] = new TH1D(hname.c_str(),title.c_str(),Nbins,0.98,1.15);
      hst[i] -> Sumw2(true);
      const auto& y = Y[i];
      for ( unsigned int j = 0; j < mia.size(); ++j ) {
         hst[i] -> Fill(mia[j].mkk,y[j]);
      }
   }

   TPaveText* pt[10];
   for ( int i = 0; i < Ny; ++i ) {
      double dy = -0.1*(i==2);
      pt[i] = new TPaveText(0.75,0.75+dy,0.89,0.89+dy,"NDC");
      pt[i] -> SetTextAlign(12);
      pt[i] -> SetTextFont(42);
      pt[i] -> AddText( Form("#bf{#LTY_{%d}^{0}#GT}",i) );
   }

   TCanvas* c1 = new TCanvas("c1","Mkk low",0,0,1500,900);
   c1->cd();
//    c1->Divide(3,1);
   c1->Divide(3,2);

   for ( int i = 0; i < Ny; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();
      SetHstFace(hst[i]);
      hst[i] -> SetLineWidth(2);
      hst[i] -> SetLineColor(kBlack);
      hst[i] -> SetMarkerStyle(20);
      hst[i] -> GetYaxis() -> SetTitleOffset(1.25);
      hst[i] -> Draw("EP");
      pt[i] -> Draw();
   }

   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void plot_PSphi(const vector<Miaa>& mia,const vector<double> Y[],string pdf) {
//-------------------------------------------------------------------------
   const int Nbins = 60;
   string title = ";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2}";
   title += ";Entries / 2 MeV/c^{2}";
   TH1D* hst[3];
   for ( int i = 0; i < 3; ++i ) {
      string hname = "AA" + to_string(i);
      hst[i] = new TH1D(hname.c_str(),title.c_str(),Nbins,0.98,1.1);
      hst[i] -> Sumw2(true);
   }
   const auto& y0 = Y[0];
   const auto& y1 = Y[1];
   const auto& y2 = Y[2];
   for ( unsigned int j = 0; j < mia.size(); ++j ) {
      double s2 = sqrt(4*M_PI)*y0[j] - sqrt(5*M_PI)*y2[j];
      hst[0] -> Fill(mia[j].mkk,s2);
      double p2 = sqrt(5*M_PI)*y2[j];
      hst[1] -> Fill(mia[j].mkk,p2);
      double tmp = (2*y0[j]-sqrt(5)*y2[j])*sqrt(5)*y2[j];
      double cosphi = y1[j] / sqrt(fabs(tmp));
//       hst[2] -> Fill(mia[j].mkk,tmp); // DEBUG
      hst[2] -> Fill(mia[j].mkk,cosphi);
   }

   TPaveText* pt[3];
   for ( int i = 0; i < 3; ++i ) {
      pt[i] = new TPaveText(0.75,0.65,0.89,0.79,"NDC");
      pt[i] -> SetTextAlign(12);
      pt[i] -> SetTextFont(42);
   }
   pt[0] -> AddText("#bf{|S|^{2}}");
   pt[1] -> AddText("#bf{|P|^{2}}");
   pt[2] -> AddText("#bf{cos(#phi)}");

   TCanvas* c1 = new TCanvas("c1","Mkk low",0,0,1800,600);
   c1->cd();
   c1->Divide(3,1);

   for ( int i = 0; i < 3; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();
      SetHstFace(hst[i]);
      hst[i] -> SetLineWidth(2);
      hst[i] -> SetLineColor(kBlack);
      hst[i] -> SetMarkerStyle(20);
      hst[i] -> GetYaxis() -> SetTitleOffset(1.25);
      hst[i] -> Draw("EP");
      pt[i] -> Draw();
   }

   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void Mkk_low() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(62);
//    gStyle->SetLegendTextSize(0.03);
//    gStyle->SetStatFont(62);

//    vector<Miaa> mia = get_cosT("data_12psip_all.root");
//    plot_cosT(mia,"cosT_ini.pdf");
//    vector<Miaa> mia = get_cosT("mckketa_kkmc_12.root");
//    plot_cosT(mia,"cosT_kketa_ini.pdf");
   vector<Miaa> mia = get_cosT("mcsig_kkmc_12.root");
   plot_cosT(mia,"cosT_sig_ini.pdf");

//    vector<double> Y[5];
//    get_Yi(mia,Y);
//    plot_Yi(mia,Y,"avY.pdf");
//    plot_Yi(mia,Y,"avY_sig.pdf");
//    plot_Yi(mia,Y,"avY_kketa.pdf");

//    plot_PSphi(mia,Y,"PSphi.pdf");

}
