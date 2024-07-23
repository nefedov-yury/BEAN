// trk_eff_fit.cc - The reconstruction efficiency (data and MC)
// and their ratio (data/MC) for pions and kaons.
// Get corrections as function of P_t and Z.
//
// Kolmogorov–Smirnov  && Chi2 probability tests to
// compare the ratios for K+ and K- (pi+ and pi-)

// #include <filesystem>
// namespace fs = std::filesystem;

#include "RewTrkPiK.hpp"    // RewTrk functions with HC

// {{{1 Common parameters: Params
//--------------------------------------------------------------------
struct Params {
   Params(int dat, int kpi, int pm, int rew);

   TFile* OpenFile(int mc);
   TTree* GetEff(int mc);
   const char* Sdate() { return Form("%i",date); }

   // name of folder with root files
   // string Dir = "prod_v709/";
   const string Dir = "prod_v709n3/";
   const bool use_nohc = false;
   string datafile;
   string mcincfile;

   int date;    // 2009 ...
   int use_kpi; // 1: kaons; 2: pions
   int use_pm;  //  0: sum "+" and "-"
                //  1: "+" only
                // -1: "-" only
   int use_rew; // 0 - no weights;
                // 1 - calculate weights

   // cuts for ntuples
   TCut Cmcsig; // mc-signal
   TCut Ckaon;  // std cuts for kaons
   TCut CKPtC;  // cuts for Pt and cos(Theta) of kaons
   TCut Cpion;  // std cuts for pions
   TCut CPiPtC; // cuts for Pt and cos(Theta) of pions
};

// {{{2 > ctor
//--------------------------------------------------------------------
Params::Params(int dat, int kpi = 1, int pm = 0, int rew = 0) {
//--------------------------------------------------------------------
   date = dat;
   use_kpi = kpi;
   use_pm = pm;
   use_rew = rew;

   // set the names of data and mc files
   datafile  = string( Form("data_%02ipsip_all.root",date%100) );
   mcincfile = string( Form("mcinc_%02ipsip_all.root",date%100) );
   if ( use_nohc ) {
      mcincfile = "NoHC/" + mcincfile;
      printf("ATTENTION:"
            " You are using MC without helix corrections\n");
   }

   // CUTS ATTENTION: check also Fill_PIhst() and Fill_Khst()

   // mc-signal
   Cmcsig = TCut("good==1");

   // std cuts for kaons: mk**2 = 0.243717 (GeV**2)
   Ckaon = TCut("fabs(Mrec-3.097)<0.002&&"
         "fabs(Mrk2-0.243717)<0.025");
   // cuts for Pt and cos(Theta) of kaons
   CKPtC = TCut("0.1<Ptk&&Ptk<1.4&&fabs(Ck)<0.8");

   // std cuts for pions: mpi**2 = 0.0194797849
   Cpion = TCut("fabs(MppKK-3.097)<0.01&&"
         "fabs(Mrpi2-0.0194798)<0.025");
   // cuts for Pt and cos(Theta) of pions
   // CPiPtC = TCut("0.05<Ptpi&&Ptpi<0.4&&fabs(Cpi)<0.8");
   CPiPtC = TCut("0.025<Ptpi&&Ptpi<0.425&&fabs(Cpi)<0.8");
}

// {{{2 > OpenFile()
//--------------------------------------------------------------------
TFile* Params::OpenFile(int mc) {  // 1 for MC
//--------------------------------------------------------------------
   string dfname = Dir + ( (mc!=1) ? datafile : mcincfile );
   // cout " OpenFile: " << dfname << endl;
   TFile* froot = TFile::Open(dfname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << dfname << endl;
      exit(EXIT_FAILURE);
   }
   froot->cd("PipPimKpKm");
   return froot;
}

// {{{2 > GetEff() : root-tree for K or Pi
//--------------------------------------------------------------------
TTree* Params::GetEff(int mc = 0) {
//--------------------------------------------------------------------
   TFile* froot = this->OpenFile(mc);
   string name( ((use_kpi == 1) ? "eff_K" : "eff_Pi") );
   TTree* eff = (TTree*)gDirectory->Get(name.c_str());
   if ( !eff ) {
      cout << "can not find " << name << " in "
         << froot->GetName() << endl;
      exit(EXIT_FAILURE);
   }
   return eff;
}

// {{{1 helper functions
//--------------------------------------------------------------------
constexpr double SQ(double x) {
//--------------------------------------------------------------------
   return x*x;
}

//--------------------------------------------------------------------
void SetHstFace(TH1* hst) {
//--------------------------------------------------------------------
   TAxis* X = hst->GetXaxis();
   if ( X ) {
      X->SetLabelFont(62);
      X->SetLabelSize(0.04);
      X->SetTitleFont(42);
      X->SetTitleSize(0.05);
   }
   TAxis* Y = hst->GetYaxis();
   if ( Y ) {
      Y->SetLabelFont(62);
      Y->SetLabelSize(0.04);
      Y->SetTitleFont(42);
      Y->SetTitleSize(0.05);
   }
   TAxis* Z = hst->GetZaxis();
   if ( Z ) {
      Z->SetLabelFont(62);
      Z->SetLabelSize(0.04);
      Z->SetTitleFont(42);
      Z->SetTitleSize(0.04);
   }
}

// for small plots
//--------------------------------------------------------------------
template<class Tmp>
void SetHstFaceTbl(Tmp hst) {
//--------------------------------------------------------------------
   TAxis* X = hst->GetXaxis();
   if ( X ) {
      X->SetLabelFont(42);
      X->SetLabelSize(0.06);
      X->SetTitleFont(42);
      X->SetTitleSize(0.06);
      X->SetTickLength(0.05);
   }
   TAxis* Y = hst->GetYaxis();
   if ( Y ) {
      Y->SetLabelFont(42);
      Y->SetLabelSize(0.06);
      Y->SetTitleFont(42);
      Y->SetTitleSize(0.06);
   }
}

// {{{1 Fill Pi-histograms
//--------------------------------------------------------------------
void Fill_PIhst( Params* p, int mc, const vector<TH1D*>& hst,
      double Ptmin, double Ptmax ) {
//--------------------------------------------------------------------
   TTree* eff_Pi = p->GetEff(mc);

   //Declaration of leaves types
   //        those not used here are commented out
   Double_t        Zpi;
   Double_t        Ptpi;
   Double_t        Cpi;
   Double_t        fl;
   // Double_t        dP;
   // Double_t        dTh;
   Double_t        MppKK;
   Double_t        Mrpi2;
   // Double_t        Egsum;
   // Double_t        Egmax;
   // Double_t        good;

   // Set branch addresses.
   eff_Pi->SetBranchAddress("Zpi",&Zpi);
   eff_Pi->SetBranchAddress("Ptpi",&Ptpi);
   eff_Pi->SetBranchAddress("Cpi",&Cpi);
   eff_Pi->SetBranchAddress("fl",&fl);
   // eff_Pi->SetBranchAddress("dP",&dP);
   // eff_Pi->SetBranchAddress("dTh",&dTh);
   eff_Pi->SetBranchAddress("MppKK",&MppKK);
   eff_Pi->SetBranchAddress("Mrpi2",&Mrpi2);
   // eff_Pi->SetBranchAddress("Egsum",&Egsum);
   // eff_Pi->SetBranchAddress("Egmax",&Egmax);
   // eff_Pi->SetBranchAddress("good",&good);

   // List of cuts:
   auto c_pion = [](double MppKK, double Mrpi2)->bool{
      return (fabs(MppKK-3.097)<0.01) &&
             (fabs(Mrpi2-0.0194798)<0.025);
   };
   auto c_Ptbin = [Ptmin,Ptmax](double Pt)->bool {
      return (Ptmin<=Pt && Pt<Ptmax);
   };

   bool CosAndPt = (hst.size() == 4);

   int rew = 0;
   if ( mc == 1 ) {
      rew = p->use_rew;
   }

   Long64_t nentries = eff_Pi->GetEntries();
   for ( Long64_t i = 0; i < nentries; ++i ) {
      eff_Pi->GetEntry(i);
      if ( !c_pion(MppKK,Mrpi2) ) {
         continue;
      }
      if ( !c_Ptbin(Ptpi) ) {
         continue;
      }
      if ( CosAndPt ) {
         if ( !(fabs(Cpi) < 0.8) ) {
            continue;
         }
      }

      if ( p->use_pm == 1  && Zpi < 0 ) {
         continue;
      }
      if ( p->use_pm == -1  && Zpi > 0 ) {
         continue;
      }

      double W = 1; // weights: use only in numerator of efficiency
      if ( rew == 1 ) {
         W = RewTrkPi(p->date, Ptpi, Zpi);
      }

      hst[0]->Fill(Cpi);
      bool trk_pid_ok = (fl > 1.5); // found track and it is pion
      if ( trk_pid_ok ) {
         hst[1]->Fill(Cpi,W);
      }
      if ( CosAndPt ) {
         hst[2]->Fill(Ptpi);
         if ( trk_pid_ok ) {
            hst[3]->Fill(Ptpi,W);
         }
      }
   }
}

// {{{1 Fill K-histograms
//--------------------------------------------------------------------
void Fill_Khst( Params* p, int mc, const vector<TH1D*>& hst,
      double Ptmin, double Ptmax ) {
//--------------------------------------------------------------------
   TTree* eff_K = p->GetEff(mc);

   //Declaration of leaves types
   //        those not used here are commented out
   Double_t        Zk;
   Double_t        Ptk;
   Double_t        Ck;
   Double_t        fl;
   // Double_t        dP;
   // Double_t        dTh;
   Double_t        Mrec;
   Double_t        Mrk2;
   // Double_t        Egsum;
   // Double_t        Egmax;
   // Double_t        good;

   // Set branch addresses.
   eff_K->SetBranchAddress("Zk",&Zk);
   eff_K->SetBranchAddress("Ptk",&Ptk);
   eff_K->SetBranchAddress("Ck",&Ck);
   eff_K->SetBranchAddress("fl",&fl);
   // eff_K->SetBranchAddress("dP",&dP);
   // eff_K->SetBranchAddress("dTh",&dTh);
   eff_K->SetBranchAddress("Mrec",&Mrec);
   eff_K->SetBranchAddress("Mrk2",&Mrk2);
   // eff_K->SetBranchAddress("Egsum",&Egsum);
   // eff_K->SetBranchAddress("Egmax",&Egmax);
   // eff_K->SetBranchAddress("good",&good);

   // List of cuts:
   auto c_kaon = [](double Mrec, double Mrk2)->bool{
      return (fabs(Mrec-3.097)<0.002) && (fabs(Mrk2-0.243717)<0.025);
   };
   auto c_Ptbin = [Ptmin,Ptmax](double Pt)->bool {
      return (Ptmin<=Pt && Pt<Ptmax);
   };

   bool CosAndPt = (hst.size() == 4);

   int rew = 0;
   if ( mc == 1 ) {
      rew = p->use_rew;
   }

   Long64_t nentries = eff_K->GetEntries();
   for ( Long64_t i = 0; i < nentries; ++i ) {
      eff_K->GetEntry(i);
      if ( !c_kaon(Mrec,Mrk2) ) {
         continue;
      }
      if ( !c_Ptbin(Ptk) ) {
         continue;
      }
      if ( CosAndPt ) {
         if ( !(fabs(Ck) < 0.8) ) {
            continue;
         }
      }

      if ( p->use_pm == 1  && Zk < 0 ) {
         continue;
      }
      if ( p->use_pm == -1  && Zk > 0 ) {
         continue;
      }

      double W = 1; // weights: use only in numerator of efficiency
      if ( rew == 1 ) {
         W = RewTrk_K(p->date, Ptk, Zk);
      }

      hst[0]->Fill(Ck);
      bool trk_pid_ok = (fl > 1.5); // found track and it is kaon
      if ( trk_pid_ok ) {
         hst[1]->Fill(Ck,W);
      }
      if ( CosAndPt ) {
         hst[2]->Fill(Ptk);
         if ( trk_pid_ok ) {
            hst[3]->Fill(Ptk,W);
         }
      }
   }
}

// {{{1 Get histograms as a function of cos and Pt
//--------------------------------------------------------------------
 vector<TH1D*> get_hst(Params* p, int mc) {
//--------------------------------------------------------------------
   vector<TH1D*> hst(4,nullptr);

   auto hPtK = [](string nm) {
      vector<double> xbins {0.05, 0.15};
      for (double x = 0.2; x < 1.36; x += 0.05) {
         xbins.push_back(x);
      }
      xbins.push_back(1.45);
      return new TH1D(nm.c_str(),"",xbins.size()-1,xbins.data());
      // return new TH1D(nm.c_str(),"",28,0.05,1.45);
      // return new TH1D(nm.c_str(),"",26,0.1,1.4);
   };
   auto hPtPi= [](string nm) {
      vector<double> xbins {0.025, 0.075};
      for (double x = 0.1; x < 0.376; x += 0.025) {
         xbins.push_back(x);
      }
      xbins.push_back(0.425);
      return new TH1D(nm.c_str(),"",xbins.size()-1,xbins.data());
      // return new TH1D(nm.c_str(),"",16,0.025,0.425);
      // return new TH1D(nm.c_str(),"",14,0.05,0.4);
   };
   auto hC = [](string nm) {
      return (new TH1D(nm.c_str(),"",32,-0.8,0.8));
   };

   string hname( ((mc == 0) ? "data" : "mc") );
   hname += to_string(p->use_pm); // for KS-test
   hst[0] = hC( hname + "_C0" );
   hst[1] = hC( hname + "_C1" );
   if ( p->use_kpi == 1 ) {   // Kaons
      hst[2] = hPtK( hname+"_Pt0" );
      hst[3] = hPtK( hname+"_Pt1" );
   } else {                     // Pions
      hst[2] = hPtPi( hname+"_Pt0" );
      hst[3] = hPtPi( hname+"_Pt1" );
   }
   for ( auto& h : hst ) {
      h->Sumw2(true);
   }

   // Fill
   if ( p->use_kpi == 1 ) {   // Kaons
      // double Ptmin = 0.1, Ptmax = 1.4;
      double Ptmin = 0.05, Ptmax = 1.45;
      Fill_Khst( p, mc, hst, Ptmin, Ptmax);
   } else {                     // Pions
      // double Ptmin = 0.05, Ptmax = 0.4;
      double Ptmin = 0.025, Ptmax = 0.425;
      Fill_PIhst( p, mc, hst, Ptmin, Ptmax);
   }

   return hst;
}

// {{{1 Efficiency as a function of Pt and cos
//--------------------------------------------------------------------
vector<TH1D*> get_eff(Params* p, int mc) {
//--------------------------------------------------------------------
   vector<TH1D*> eff(2,nullptr);

   vector<TH1D*> hst = get_hst(p, mc);

   // attention change order: eff[0] is Pt and [1] is cos
   for( int i = 0; i < 2; i++ ) {
      int ih = 2*i;
      string name = string("eff_") + string(hst[ih]->GetName());
      eff[1-i]=(TH1D*)hst[ih]->Clone( name.c_str() );
      eff[1-i]->Divide(hst[ih+1],hst[ih],1.,1.,"B");
   }

   return eff;
}

// {{{1 Efficiency ratio DATA/MC
//--------------------------------------------------------------------
vector<TH1D*> get_ratio( const vector<TH1D*>& effdat,
      const vector<TH1D*>& effmc) {
//--------------------------------------------------------------------
   int Nh = effdat.size();
   vector<TH1D*> rat(Nh,nullptr);
   for ( int i = 0; i < Nh; ++i ) {
      string name = string("r_") + string( effdat[i]->GetName() );
      rat[i]=(TH1D*)effdat[i]->Clone( name.c_str() );
      rat[i]->Divide(effmc[i]);
   }
   return rat;
}

// {{{1 Test +/- identity
//--------------------------------------------------------------------
void test_PM(int date, int kpi, int test=1) {
//--------------------------------------------------------------------
// compare K+ and K- (pi+ and pi-) 1D distributions of the data/MC
// ratio to check their identity

   Params* parP = new Params(date,kpi,1,0); // date, K/pi, "+", rew
   vector<TH1D*> eff_d  = get_eff( parP, 0 );
   vector<TH1D*> eff_mc = get_eff( parP, 1 );
   vector<TH1D*> ratP   = get_ratio( eff_d, eff_mc );

   Params* parM = new Params(date,kpi,-1,0); // date, K/pi, "-", rew
   eff_d  = get_eff( parM, 0 );
   eff_mc = get_eff( parM, 1 );
   vector<TH1D*> ratM   = get_ratio( eff_d, eff_mc );

   if ( test == 1 ) {
      cout << "Kolmogorov–Smirnov";
   } else if ( test == 2 ) {
      cout << "chi^2";
   }
   cout << " probability test for data/MC efficiency ratio\n";

   cout << " data-" << date << ":\n";
   bool verbose = false;
   for (int i = 1; i >= 0; --i ) {
      string info = (kpi == 1) ?  " K+ <--> K- " : " pi+ <--> pi- ";
      info += (i==0) ? "; P_{t}:      " : "; cos(#Theta):";
      if ( test == 1 ) {
         // use the Kolmogorov–Smirnov test:
         // http://root.cern.ch/root/html/TH1.html#TH1:KolmogorovTest
         if ( verbose ) {
            cout << info << endl;
            ratP[i]->KolmogorovTest(ratM[i],"D");
            cout << endl;
         } else {
            printf("%s p-value(chi2)= %.3f Max Dist= %.4f\n",
                  info.c_str(), ratP[i]->KolmogorovTest(ratM[i]),
                  ratP[i]->KolmogorovTest(ratM[i],"M") );
         }
      } else if ( test == 2 ) {
         // use chi^2 test:
         // http://root.cern.ch/root/html/TH1.html#TH1:Chi2Test
         if ( verbose ) {
            cout << info << endl;
            ratP[i]->Chi2Test(ratM[i],"WWP");
            cout << endl;
         } else {
            int igood = 0, ndf = 0;
            double chi2 = 0;
            double pval =
               ratP[i]->Chi2TestX(ratM[i],chi2,ndf,igood,"WW");
            printf("%s p-value(chi2)= %.3f chi^2/NDF= %.1f/%d\n",
                  info.c_str(),pval,chi2,ndf);
         }
      }
   }
}

// {{{1 Plot Pi
//--------------------------------------------------------------------
void plot_pict_pi(int date, int pm, int rew, int Cx=600, int Cy=600) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,pm,rew); // pi, "+/-"
   vector<TH1D*> eff_d  = get_eff( par, 0 ); // data
   vector<TH1D*> eff_mc = get_eff( par, 1 ); // mc
   vector<TH1D*> rat0   = get_ratio( eff_d, eff_mc );

   int Nh = eff_d.size(); // must be 2

   // Attributes of draw
   auto setBlue = [](TH1* h) {
      h->SetLineColor(kBlue+1);
      h->SetMarkerColor(kBlue+2);
      h->SetMarkerStyle(22);
      h->SetMarkerSize(0.9);
   };
   auto setRed = [](TH1* h) {
      h->SetLineColor(kRed+1);
      h->SetMarkerColor(kRed+2);
      h->SetMarkerStyle(23);
      h->SetMarkerSize(0.9);
   };
   for (int i = 0; i < Nh; ++i ) {
      setBlue(eff_d[i]);
      setRed(eff_mc[i]);
   }

   // for fit
   TF1* pl0 = (TF1*)gROOT->GetFunction("pol0")->Clone();
   pl0->SetLineColor(kRed);

   // Draw:
   ////////////////////////////////////
   // double eff_min = 0.5, eff_max = 1.0;
   double rat_min = 0.85, rat_max = 1.15;
   double eff_min = 0.3, eff_max = 1.0;
   if ( rew == 1 ) {
      rat_min = 0.90; rat_max = 1.10;
   }
   ////////////////////////////////////
   string p_pdf("Pi"), p_leg("#pi");
   if ( pm == 1 ) { p_pdf += "p"; p_leg += "^{#kern[0.25]{#plus}}";}
   if ( pm == -1 ){ p_pdf += "m"; p_leg += "^{#kern[0.25]{#minus}}";}

   vector<TCanvas*> cc(4,nullptr);
   vector<TLegend*> leg(4,nullptr);
   for ( size_t i = 0; i < cc.size(); ++i ) {
      int x0 = 700*(i%2);
      int y0 = 500*(i/2);
      // printf("%zu -> x0=%i y0=%i\n",i,x0,y0);
      auto name = Form("c%zu_%s_%i",1+i,p_pdf.c_str(),date);
      cc[i] = new TCanvas( name,name, x0,y0,Cx,Cy);

      auto cdat = Form("%i  %s",date,p_leg.c_str());
      if ( i < 2 ) { // efficiency
         leg[i] = new TLegend(0.61,0.30,0.89,0.50);
         leg[i]->SetHeader( cdat,"C");
         leg[i]->AddEntry(eff_d[0],  "  data ", "LP");
         leg[i]->AddEntry(eff_mc[0], "  MC ",   "LP");
      } else { // retio
         leg[i] = new TLegend(0.14,0.81,0.42,0.89);
         leg[i]->SetHeader( cdat,"C");
      }
   }

   // plot efficiencies
   for (int i = 0; i < Nh; ++i ) {
      auto& ci = cc[0+i];
      ci->cd();
      gPad->SetGrid();
      string title =
         string( ((i==0) ? ";P_{t}, GeV/c" : ";cos(#Theta)") ) +
         string(";#epsilon(#pi)");
      eff_d[i]->SetTitle(title.c_str());
      SetHstFace(eff_d[i]);
      eff_d[i]->SetAxisRange(eff_min,eff_max,"Y");
      eff_d[i]->GetYaxis()->SetNdivisions(505);
      eff_d[i]->GetYaxis()->SetTitleOffset(0.9);
      eff_d[i]->GetXaxis()->SetTitleOffset(0.9);
      eff_d[i]->Draw("E");
      eff_mc[i]->Draw("E SAME");
      leg[0+i]->Draw();

      gPad->RedrawAxis();
      ci->Update();
      string pp = p_pdf + ( ( i == 0 ) ? "_pt" : "_cos" ) +
         ( ( rew == 1 ) ? "_w" : "" );
      ci->Print( Form("trkeff_%s_%i.pdf",pp.c_str(),date) );
   }

   // plot ratio
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.25);
   gStyle->SetFitFormat(".4f");

   for (int i = 0; i < Nh; ++i ) {
      auto& ci = cc[2+i];
      ci->cd();
      gPad->SetLeftMargin(gPad->GetLeftMargin()+0.02);
      gPad->SetRightMargin(gPad->GetRightMargin()-0.02);
      gPad->SetGrid();
      string title =
         string( (i==0) ? ";P_{t}, GeV/c" : ";cos(#Theta)" ) +
         string(";#epsilon(data) / #epsilon(MC)");
      rat0[i]->SetTitle(title.c_str());
      SetHstFace(rat0[i]);
      rat0[i]->SetAxisRange(rat_min,rat_max,"Y");
      rat0[i]->GetYaxis()->SetTitleOffset(1.15);
      rat0[i]->GetXaxis()->SetTitleOffset(0.9);
      rat0[i]->SetLineWidth(2);
      rat0[i]->SetMarkerStyle(20);
      rat0[i]->SetLineColor(kBlack);
      rat0[i]->Fit( pl0, "" );
      leg[2+i]->Draw();
      rat0[i]->Draw("SAME E0");

      gPad->RedrawAxis();
      ci->Update();
      string pp = p_pdf + ( ( i == 0 ) ? "_pt" : "_cos" ) +
         ( ( rew == 1 ) ? "_w" : "" );
      ci->Print( Form("trkeff_rat_%s_%i.pdf",pp.c_str(),date) );
   }
}

// {{{1 Plot K
//--------------------------------------------------------------------
void plot_pict_K(int date, int pm, int rew, int Cx=600, int Cy=600) {
//--------------------------------------------------------------------
   Params* par = new Params(date,1,pm,rew); // K, "+/-"
   vector<TH1D*> eff_d  = get_eff( par, 0 );
   vector<TH1D*> eff_mc = get_eff( par, 1 );
   vector<TH1D*> rat0   = get_ratio( eff_d, eff_mc );

   int Nh = eff_d.size(); // must be 2

   // Attributes of draw
   auto setBlue = [](TH1* h) {
      h->SetLineColor(kBlue+1);
      h->SetMarkerColor(kBlue+2);
      h->SetMarkerStyle(22);
      h->SetMarkerSize(0.9);
   };
   auto setRed = [](TH1* h) {
      h->SetLineColor(kRed+1);
      h->SetMarkerColor(kRed+2);
      h->SetMarkerStyle(23);
      h->SetMarkerSize(0.9);
   };
   for (int i = 0; i < Nh; ++i ) {
      setBlue(eff_d[i]);
      setRed(eff_mc[i]);
   }

   // for fit
   TF1* pl0 = (TF1*)gROOT->GetFunction("pol0")->Clone();
   pl0->SetLineColor(kRed);

   // Draw:
   ////////////////////////////////////
   double eff_min = 0.2, eff_max = 1.0;
   double rat_min = 0.85, rat_max = 1.15;
   if ( rew == 1 ) {
      rat_min = 0.90; rat_max = 1.10;
   }
   ////////////////////////////////////
   string p_pdf("K"), p_leg("K");
   if ( pm == 1 ) { p_pdf += "p"; p_leg += "^{#plus}";}
   if ( pm == -1 ){ p_pdf += "m"; p_leg += "^{#minus}";}

   vector<TCanvas*> cc(4,nullptr);
   vector<TLegend*> leg(4,nullptr);
   for ( size_t i = 0; i < cc.size(); ++i ) {
      int x0 = 700*(i%2);
      int y0 = 500*(i/2);
      // printf("%zu -> x0=%i y0=%i\n",i,x0,y0);
      auto name = Form("c%zu_%s_%i",1+i,p_pdf.c_str(),date);
      cc[i] = new TCanvas( name,name, x0,y0,Cx,Cy);

      auto cdat = Form("%i  %s",date,p_leg.c_str());
      if ( i < 2 ) { // efficiency
         leg[i] = new TLegend(0.61,0.30,0.89,0.50);
         leg[i]->SetHeader( cdat,"C");
         leg[i]->AddEntry(eff_d[0],  "  data ", "LP");
         leg[i]->AddEntry(eff_mc[0], "  MC ",   "LP");
      } else { // retio
         leg[i] = new TLegend(0.14,0.81,0.42,0.89);
         leg[i]->SetHeader( cdat,"C");
      }
   }

   // plot efficiencies
   for (int i = 0; i < Nh; ++i ) {
      auto& ci = cc[0+i];
      ci->cd();
      gPad->SetGrid();
      eff_d[i]->SetAxisRange(eff_min,eff_max,"Y");
      string title =
         string( ((i==0) ? ";P_{t}, GeV/c" : ";cos(#Theta)") ) +
         string(";#epsilon (K)");
      eff_d[i]->SetTitle(title.c_str());
      SetHstFace(eff_d[i]);
      eff_d[i]->GetYaxis()->SetTitleOffset(0.9);
      eff_d[i]->GetXaxis()->SetTitleOffset(0.9);
      eff_d[i]->Draw("E");
      eff_mc[i]->Draw("E SAME");
      leg[0+i]->Draw();

      gPad->RedrawAxis();
      ci->Update();
      string pp = p_pdf + ( ( i == 0 ) ? "_pt" : "_cos" ) +
         ( ( rew == 1 ) ? "_w" : "" );
      ci->Print( Form("trkeff_%s_%i.pdf",pp.c_str(),date) );
   }

   // plot ratio
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.25);
   gStyle->SetFitFormat(".4f");

   for (int i = 0; i < Nh; ++i ) {
      auto& ci = cc[2+i];
      ci->cd();
      gPad->SetLeftMargin(gPad->GetLeftMargin()+0.02);
      gPad->SetRightMargin(gPad->GetRightMargin()-0.02);
      gPad->SetGrid();
      rat0[i]->SetAxisRange(rat_min,rat_max,"Y");
         // string( Form("%s, data/MC %i",p_leg.c_str(),date) )
      string title =
         string( (i==0) ? ";P_{t}, GeV/c" : ";cos(#Theta)" ) +
         string(";#epsilon (data) / #epsilon (MC)");
      rat0[i]->SetTitle(title.c_str());
      SetHstFace(rat0[i]);
      rat0[i]->GetYaxis()->SetTitleOffset(1.15);
      rat0[i]->GetXaxis()->SetTitleOffset(0.9);
      rat0[i]->SetLineWidth(2);
      rat0[i]->SetMarkerStyle(20);
      rat0[i]->SetLineColor(kBlack);
      rat0[i]->Fit( pl0, "" );
      leg[2+i]->Draw();
      rat0[i]->Draw("SAME E0");

      gPad->RedrawAxis();
      ci->Update();
      string pp = p_pdf + ( ( i == 0 ) ? "_pt" : "_cos" ) +
         ( ( rew == 1 ) ? "_w" : "" );
      ci->Print( Form("trkeff_rat_%s_%i.pdf",pp.c_str(),date) );
   }
}

// {{{1 Get average Pt value in [Ptmin,Ptmax] and |cos(Th)|<0.8
//--------------------------------------------------------------------
double get_Pt_mean(Params* p, double Ptmin, double Ptmax) {
//--------------------------------------------------------------------
   // TTree* eff = p->GetEff(1); // use MC distribution, we correct MC
   TTree* eff = p->GetEff(0);

   int nbins = 20;
   TH1D* hst = new TH1D( "hpt","",nbins,Ptmin,Ptmax);
   hst->Sumw2(true);

   // Fill
   if ( p->use_kpi == 1 ) {   // Kaons
      TCut cutD = p->Ckaon + TCut("fabs(Ck)<0.8");
      if ( p->use_pm == 1 ) { // K+
         cutD += TCut("Zk>0");
      } else if ( p->use_pm == -1 ) { // K-
         cutD += TCut("Zk<0");
      }
      eff->Draw("Ptk>>hpt",cutD,"goff");
   } else {                   // Pions
      TCut cutD = p->Cpion + TCut("fabs(Cpi)<0.8");
      if ( p->use_pm == 1 ) { // pi+
         cutD += TCut("Zpi>0");
      } else if ( p->use_pm == -1 ) { // pi-
         cutD += TCut("Zpi<0");
      }
      eff->Draw("Ptpi>>hpt",cutD,"goff");
   }

   double ptmean = hst->GetMean();
   delete hst;

   return ptmean;
}

// {{{1 Get efficiency ratio as a function of cos(Theta) to fit
//--------------------------------------------------------------------
TH1D* get_eff_cos(Params* p, int mc, double Ptmin, double Ptmax) {
//--------------------------------------------------------------------
   vector<TH1D*> hst(2,nullptr);
   int nbins = ( p->date == 2009 ) ? 20 : 40;
   for ( int i = 0; i < 2; ++i ) {
      hst[i] = new TH1D( Form("C%i",i),"",nbins,-1.,1.);
      hst[i]->Sumw2(true);
   }

   // Fill
   if ( p->use_kpi == 1 ) {   // Kaons
      Fill_Khst( p, mc, hst, Ptmin, Ptmax);
   } else {                   // Pions
      Fill_PIhst( p, mc, hst, Ptmin, Ptmax);
   }

   // calculate efficiency
   auto name = Form("%s_%.0f", (mc==0) ? "data" : "mc", Ptmin*100);
   TH1D* heff = (TH1D*) hst[0]->Clone( name );
   heff->Divide(hst[1],hst[0],1.,1.,"B");
   // heff->Divide(hst[1],hst[0]); // poison errors

   delete hst[0];
   delete hst[1];

   return heff;
}

// {{{1 Fit Ratio
//--------------------------------------------------------------------
void FitRatio(int date, int kpi, int pm=0, int rew=0) {
//--------------------------------------------------------------------
   Params* par = new Params(date,kpi,pm,rew); // date,K/pi,"+/-",rew

   double Ptmin, Ptmax, Nbins; // binning for Pt
   if ( kpi == 1 ) {        // Kaons
      Ptmin = 0.05;
      Ptmax = 1.45;
      Nbins = 14;
   } else if ( kpi == 2 ) { // Pions
      Ptmin = 0.025;
      Ptmax = 0.425;
      Nbins = 8;
   }

   // common for name of pdf-files
   string p_pdf( ((kpi==1) ? "K" : "Pi") );
   if ( pm == +1 ) { p_pdf += "p"; }
   if ( pm == -1 ) { p_pdf += "m"; }
   if ( rew == 1 ) { p_pdf += "_w";}

   // First fit the ratio as the function of "cos(Theta)"
   // function to the first fit
   TF1* fit1 = (TF1*)gROOT->GetFunction("pol0")->Clone("fit1");
   fit1->SetLineColor(kRed);
   fit1->SetLineWidth(2);

   // Second fit the obtained parameters as the function of "Pt"
   TF1* fit2 = nullptr;
   // or use the spline if fit2 is zero
   TSpline3* sp = nullptr;
   if ( rew == 1 ) {
      fit2 = (TF1*)gROOT->GetFunction("pol0")->Clone("fit2");
   }
   // else {
      // if ( kpi == 1 ) {
         // // functions to the second fit of Kaons
         // int nch = 0; // nch is the order of Chebyshev polynomials
         // if ( date == 2009 ) {
            // nch = 1;  // straight line
         // } else if ( date == 2012 ) {
            // nch = 3;
         // } else if ( date == 2021 ) {
            // nch = 5;
         // }
         // double xmin = Ptmin, xmax = Ptmax;
         // auto Lchb =
            // [nch,xmin,xmax](const double* xx, const double* p) {
               // if ( nch == 0 ) {
                  // return p[0];
               // }
               // // [xmin,xmax]->[-1,+1]
               // double x = (2*xx[0]-xmin-xmax)/(xmax-xmin);
               // double sum = p[0] + x*p[1];
               // if ( nch == 1 ) {
                  // return sum;
               // }
               // double T0 = 1, T1 = x;
               // for ( int i = 2; i <= nch; ++i ) {
                  // double Tmp = 2*x*T1 - T0;
                  // sum += p[i]*Tmp;
                  // T0 = T1;
                  // T1 = Tmp;
               // }
               // return sum;
            // };

         // fit2 = new TF1("fit2", Lchb, Ptmin, Ptmax,nch+1);
         // for ( size_t ipar = 0; ipar <= nch; ++ipar ) {
            // double ini_val = ( ipar==0 ) ? 0.99 : 0.01 ;
            // fit2->SetParameter(ipar, ini_val);
         // }
      // }
   // }

   if ( fit2 ) {
      fit2->SetRange(Ptmin,Ptmax);
      fit2->SetLineColor(kRed);
      fit2->SetLineWidth(2);
   }

   // ++++ First fit ++++
   TCanvas* c1 = new TCanvas("c1","...",850,0,700,1000);
   int Nc1 = 14;
   if ( kpi == 1 ) {    // Kaons
      c1->Divide(2,7);
   } else {
      Nc1 = 8;
      c1->Divide(2,4);
   }
   int Ic1 = 1;

   string pdf1 = "rat_fit1_"+p_pdf+"_"+to_string(date)+".pdf";
   if ( rew == 0 ) {
      c1->Print( (pdf1+"[").c_str() ); // open pdf-file
   }

   vector<double> xx(Nbins), yy(Nbins);
   vector<double> ex(Nbins,1e-3), ey(Nbins);
   double dp = (Ptmax-Ptmin)/double(Nbins);
   for (int i = 0; i < Nbins; ++i ) {
      double ptmin = Ptmin + dp*i;
      double ptmax = ptmin + dp;
      // double ptav = 0.5*(ptmin + ptmax);
      double ptav = get_Pt_mean(par,ptmin,ptmax);
      xx[i] = ptav;

      TH1D* eff_dat = get_eff_cos(par,0,ptmin,ptmax);
      TH1D* eff_mc  = get_eff_cos(par,1,ptmin,ptmax);
      string name = string("r_") + string( eff_dat->GetName() );
      TH1D* rat = (TH1D*)eff_dat->Clone( name.c_str() );
      rat->Divide(eff_dat, eff_mc);

      c1->cd(Ic1);
      gPad->SetGrid();
      rat->SetAxisRange(0.9,1.1,"Y");
      if ( i == 0 ) { // 1st bin
         rat->SetAxisRange(0.5,1.5,"Y");
      }
      if ( kpi == 1 && i == Nbins-1 ) {    // K and last bin
         rat->SetAxisRange(0.5,1.5,"Y");
      }
      rat->GetYaxis()->SetNdivisions(1004);
      rat->SetTitle(Form(
               "P_{t} #in [%.3f, %.3f] GeV"
               "; cos(#Theta)"
               ";#epsilon(data) / #epsilon(MC)",
               ptmin,ptmax));
      SetHstFaceTbl(rat);
      if ( kpi == 1 ) {
         gStyle->SetTitleFontSize(0.1);
         rat->GetYaxis()->SetTitleOffset(0.5);
      } else {
         gStyle->SetTitleFontSize(0.05);
         rat->GetYaxis()->SetTitleOffset(0.8);
      }
      rat->GetXaxis()->SetTitleOffset(0.8);
      TFitResultPtr rs = rat->Fit(fit1,"SEQ","",-0.8,0.8);
      rat->Draw("SAME E0");
      // rs->Print();
      yy[i] = fit1->GetParameter(0);
      ey[i] = fit1->GetParError(0);
      auto ch2 = fit1->GetChisquare();
      auto ndf = fit1->GetNDF();
      printf("%2i <pt>=%.3f  r=%.4f+/-%.4f  ch2/NDF= %.1f/%i\n",
            i,ptav,yy[i],ey[i],ch2,ndf);

      Ic1 += 1;
      if ( i == Nbins-1 ) { // clear all other pads
         for (; Ic1 <= Nc1; ++Ic1 ) {
             c1->cd(Ic1);
             gPad->Clear();
         }
      }
      if ( Ic1 > Nc1 ) {
         Ic1 = 1;
         c1->Update();
         if ( rew == 0 ) {
            c1->Print( pdf1.c_str() ); // add to pdf-file
         }
      }
   } // end of for
   if ( rew == 0 ) {
      c1->Print( (pdf1+"]").c_str() ); // close pdf-file
   }

   // ++++ Second fit ++++
   TCanvas* c2 = new TCanvas("c2","...",0,0,800,750);
   c2->cd();
   gPad->SetGrid();
   //  gStyle->SetFitFormat(".4f");

   TGraphErrors* gX = new TGraphErrors(Nbins,xx.data(),yy.data(),
                                             ex.data(),ey.data() );
   double Ymin = 0.9;
   double Ymax = 1.1;
   if ( rew == 1 ) {
      Ymin = 0.95;
      Ymax = 1.05;
   }
   if ( yy[0] < Ymin ) {
      Ymin = 0.1 * std::floor( (yy[0]-0.01)*10 );
   }
   gX->GetHistogram()->SetMinimum(Ymin);
   if ( yy[0] > Ymax ) {
      Ymax = 0.1 * std::ceil( (yy[0]+0.01)*10 );
   }
   gX->GetHistogram()->SetMaximum(Ymax);

   gX->SetMarkerColor(kBlue);
   gX->SetMarkerStyle(21);
   gX->Draw("AP");
   gX->GetYaxis()->SetTitleOffset(1.4);
   gX->SetTitle(
         "; P_{t}, GeV/c"
         "; #epsilon(data)/#epsilon(MC)"
         );

   TFitResultPtr res(nullptr);
   TF1* Fsp = nullptr;
   if ( fit2 ) {
      if ( kpi == 1 && rew == 0 && date == 2009 ) {
         fit2->SetRange(0.2,Ptmax); // reduce fit range
      } else if ( kpi == 1 && rew == 0 && date == 2012 ) {
         fit2->SetRange(0.3,Ptmax); // reduce fit range
      } else if ( kpi == 1 && rew == 0 && date == 2021 ) {
         fit2->SetRange(0.3,Ptmax); // reduce fit range
      }
      res = gX->Fit("fit2","SER EX0");
      gX->Draw("SAME 0P");
      if ( rew == 0 ) {
         res->Print(); // print results
         // res->Print("V"); // + error and correlation matrices
      }
   } else {
      // Save points in file for python smooth splines
      string fname = p_pdf+"_"+to_string(date)+".dat";
      FILE* fp = fopen(fname.c_str(),"w");
      if ( !fp ) {
         cerr << "Error open file " << fname << endl;
         fp = stdout;
      }
      fprintf( fp, "fX = [ " );
      for ( size_t i = 0; i < Nbins; ++i ) {
         fprintf( fp, "%f%s",
               gX->GetPointX(i), i<Nbins-1 ? ", " : " ]\n" );
      }
      fprintf( fp, "fY = [ " );
      for ( size_t i = 0; i < Nbins; ++i ) {
         fprintf( fp, "%f%s",
               gX->GetPointY(i), i<Nbins-1 ? ", " : " ]\n" );
      }
      fprintf( fp, "dY = [ " );
      for ( size_t i = 0; i < Nbins; ++i ) {
         fprintf( fp, "%f%s",
               gX->GetErrorY(i), i<Nbins-1 ? ", " : " ]\n" );
      }
      if (fp != stdout ) {
         fclose(fp);
      }

      sp = new TSpline3("sp",gX,"e2b2",0.,0.);
      // string fsp = p_pdf+"_"+to_string(date)+"_spline.cc";
      // printf("Save spline to file: %s\n", fsp.c_str());
      // sp->SaveAs(fsp.c_str());

      // to draw within [Ptmin,Ptmax]
      auto Lsp = [sp](double* x, double* p) {
         return sp->Eval(x[0]);
      };
      Fsp = new TF1("Fsp",Lsp,Ptmin,Ptmax,0);
      Fsp->SetTitle("Cubic spline");
      Fsp->SetLineColor(kRed);
      Fsp->SetLineWidth(2);
      Fsp->Draw("L SAME");
      gX->Draw("SAME 0P");
   }

   string p_leg( ((kpi==1) ? "K" : "#pi") );
   if ( pm == 0 ) { p_leg += "^{#kern[0.25]{#pm}}"; }
   if ( pm == 1 ) { p_leg += "^{#kern[0.25]{#plus}}"; }
   if ( pm == -1 ){ p_leg += "^{#kern[0.25]{#minus}}"; }
   if ( par->use_nohc ) { p_leg += " noHC"; }
   TLegend* leg = nullptr;
   if ( fit2 ) {
      leg = new TLegend(0.20,0.82,0.50,0.89);
      leg->SetHeader( Form("%i  %s",date,p_leg.c_str()),"C" );
   } else {
      leg = new TLegend(0.55,0.77,0.89,0.89);
      leg->AddEntry(gX,Form("%i  %s",date,p_leg.c_str()),"PE");
      leg->AddEntry(Fsp,Fsp->GetTitle(),"L");
   }
   leg->Draw();

   gPad->RedrawAxis();
   c2->Update();
   string pdf2 = "rat_fit_"+p_pdf+"_"+to_string(date)+".pdf";
   c2->Print( pdf2.c_str() );
}

// {{{1 Main
//--------------------------------------------------------------------
void trk_eff_fit() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(42);
   gStyle->SetStatFont(62);
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   // gStyle->SetStatW(0.25);

   // 1) Reconstruction efficiency for data and Monte Carlo as
   //   a function of Pt and cos(θ) and eff.ratio eff(data)/eff(MC)
   // size_t Cx = 630, Cy = 600; // canvas sizes
   // for ( auto date : {2009, 2012, 2021} ) {
      // for ( auto sign : {+1, -1} ) {
         // const int rew = 1;     // 1 - use corrections
         // Kaons: fig.38-39 efficiency, 40-41 efficiency ratio
         // plot_pict_K(date,sign,rew,Cx,Cy);
         // Pions: fig.42-43 efficiency, 44-45 efficiency ratio
         // plot_pict_pi(date,sign,rew,Cx,Cy);
      // }
   // }

   // 2) checking the compatibility of efficiency ratios for
   //   K+ and K- (pi+/pi-)
   // int test = 1; // 1 - Kolmogorov–Smirnov; 2 - Chi2Test
   // for ( auto date : {2009, 2012, 2021} ) {
      // test_PM(date,1,test); // Kaons
      // test_PM(date,2,test); // Pions
   // }

   // 3) Get Efficiency Corrections
   const int Kaons = 1;
   const int Pions = 2;

   const int date = 2009; // 2009, 2012, 2021
   const int Z = +1;      // 1=+, -1=-, 0=+/-
   const int rew = 1;     // use corrections, weights

   // FitRatio(date,Pions,Z);
   // FitRatio(date,Pions,Z,rew);

   // FitRatio(date,Kaons,Z);
   // FitRatio(date,Kaons,Z,rew);

}
