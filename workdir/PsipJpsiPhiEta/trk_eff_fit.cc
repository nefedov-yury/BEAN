// The reconstruction efficiency (data and MC)
// and their ratio (data/MC) for pions and kaons.
// Get corrections as function of P_t
// by fitting the ratio.
// Kolmogorov–Smirnov probability test to compare
// the ratios for K+ and K- (pi+ and pi-)
// TODO: rewighting!!!

//----------------------------------------------------------------------
// Common parameters:
struct Params {
   // names of files to use
   vector<string> fnames;

   int date;    // 2009 or 2012
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

   //-----------------------------------------------------------------
   Params(int dat = 2012, int kpi = 1, int pm = 0, int rew = 0) {
   //-----------------------------------------------------------------
      date = dat;
      use_kpi = kpi;
      use_pm = pm;
      use_rew = rew;

      // name of folder with root-files
      string dir("prod-13eff/");
      fnames = {
         "data_09psip_all.root", "mcinc_09psip_all.root",
         "data_12psip_all.root", "mcinc_12psip_all.root"
      };
      for ( auto& fn : fnames ) {
         fn = dir + fn;
      }
      // mc-signal
      Cmcsig = TCut("good==1");

      // std cuts for kaons: mk**2 = 0.243717 (GeV**2)
      Ckaon = TCut("fabs(Mrec-3.097)<0.002&&"
                   "fabs(Mrk2-0.243717)<0.025");
      // cuts for Pt and cos(Theta) of kaons
      CKPtC = TCut("0.1<Ptk&&Ptk<1.4&&fabs(Ck)<0.8");

      // std cuts for pions: mpi**2 = 0.0194797849
      Cpion = TCut("fabs(MppKK-3.096)<0.009&&"
                   "fabs(Mrpi2-0.0194798)<0.025");
      // cuts for Pt and cos(Theta) of pions
      CPiPtC = TCut("0.05<Ptpi&&Ptpi<0.4&&fabs(Cpi)<0.8");
   }

   //-----------------------------------------------------------------
   TFile* OpenFile(int mc) {  // 1 for MC
   //-----------------------------------------------------------------
      // open file
      int idx = mc + ((date==2009) ? 0 : 2);
      string dfname = fnames[idx];
      TFile* froot = TFile::Open(dfname.c_str(),"READ");
      if( froot == 0 ) {
         cerr << "can not open " << dfname << endl;
         exit(0);
      }
//       cout << " file: " << dfname << endl;
      froot->cd("PipPimKpKm");
      return froot;
   }

   //-----------------------------------------------------------------
   TTree* GetEff(int mc = 0) { // get tree for K or Pi
   //-----------------------------------------------------------------
      TFile* froot = this->OpenFile(mc);
      string name( ((use_kpi == 1) ? "eff_K" : "eff_Pi") );
      TTree* eff = (TTree*)gDirectory->Get(name.c_str());
      if ( !eff ) {
         cout << "can not find " << name << " in "
            << froot->GetName() << endl;
         exit(0);
      }
      return eff;
   }

};
//----------------------------------------------------------------------

//----------------------------------------------------------------------
constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

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

//--------------------------------------------------------------------
void get_hst(Params* p, int mc, vector<TH1D*>& hst) {
//--------------------------------------------------------------------
// get histograms: mc=0,1 / data,MC

   TTree* eff = p -> GetEff(mc);

   TCut Crec;
   if ( p -> use_kpi == 1 ) {   // Kaons
      Crec += p -> Ckaon + p -> CKPtC;
      if ( p -> use_pm == 1 ) {         // "K+"
         Crec += TCut("Zk>0");
      } else if ( p -> use_pm == -1 ) { // "K-"
         Crec += TCut("Zk<0");
      }
   } else {                     // Pions
      Crec += p -> Cpion + p -> CPiPtC;
      if ( p -> use_pm == 1 ) {         // "pi+"
         Crec += TCut("Zpi>0");
      } else if ( p -> use_pm == -1 ) { // "pi-"
         Crec += TCut("Zpi<0");
      }
   }
//    cout << " DEBUG: Crec= " << Crec.GetTitle() << endl;
   TCut Crec0 = Crec + TCut("fl<1.5");
   TCut Crec1 = Crec + TCut("fl>1.5"); // trk && pid id!

   string hname( ((mc==0) ? "data" : "mc") );
   hname += to_string(p -> use_pm); // for KS-test
   vector<string> hns = {
      hname + "_P0", hname + "_P1",
      hname + "_C0", hname + "_C1"
   };
   auto hPtK = [](string nm) {return new TH1D(nm.c_str(),"",26,0.1,1.4);};
   auto hPtPi= [](string nm) {return new TH1D(nm.c_str(),"",14,0.05,0.4);};
   auto hC = [](string nm) {return (new TH1D(nm.c_str(),"",32,-0.8,0.8));};

   hst.clear();
   hst.resize(4,nullptr);
   if ( p -> use_kpi == 1 ) {   // Kaons
      hst[0] = hPtK(hns[0]);
      hst[1] = hPtK(hns[1]);
      hst[2] = hC(hns[2]);
      hst[3] = hC(hns[3]);
   } else {                     // Pions
      hst[0] = hPtPi(hns[0]);
      hst[1] = hPtPi(hns[1]);
      hst[2] = hC(hns[2]);
      hst[3] = hC(hns[3]);
   }
   for ( auto& h : hst ) {
      h -> Sumw2(true);
   }

   // Fill
   if ( p -> use_kpi == 1 ) {   // Kaons
      eff -> Draw( ("Ptk>>"+hns[0]).c_str(),Crec0,"goff");
      eff -> Draw( ("Ck>>"+hns[2]).c_str(),Crec0,"goff");
      eff -> Draw( ("Ptk>>"+hns[1]).c_str(),Crec1,"goff");
      eff -> Draw( ("Ck>>"+hns[3]).c_str(),Crec1,"goff");
   } else {                     // Pions
      eff -> Draw( ("Ptpi>>"+hns[0]).c_str(),Crec0,"goff");
      eff -> Draw( ("Cpi>>"+hns[2]).c_str(),Crec0,"goff");
      eff -> Draw( ("Ptpi>>"+hns[1]).c_str(),Crec1,"goff");
      eff -> Draw( ("Cpi>>"+hns[3]).c_str(),Crec1,"goff");
   }
}

//--------------------------------------------------------------------
void get_eff(Params* p, int mc, vector<TH1D*>& eff ) {
//--------------------------------------------------------------------
// Get histograms and calculate efficiencies

   vector<TH1D*> hst;
   get_hst(p, mc, hst);

   eff.clear();
   eff.resize(2,nullptr);
   for( int i = 0; i < 2; i++ ) {
      int ih = 2*i;
      string name = string("eff_") + string(hst[ih]->GetName());
      eff[i]=(TH1D*)hst[ih]->Clone( name.c_str() );
      eff[i]->Add(hst[ih+1]);
      eff[i]->Divide(hst[ih+1],eff[i],1.,1.,"B");
   }
}

//--------------------------------------------------------------------
void get_ratio( const vector<TH1D*>& effdat,
                const vector<TH1D*>& effmc,
                                vector<TH1D*>& rat ) {
//--------------------------------------------------------------------
// calculate ratio of efficiencies DATA/MC

   int Nh = effdat.size();
   rat.clear();
   rat.resize(Nh,nullptr);
   for ( int i = 0; i < Nh; ++i ) {
      string name = string("r_") + string( effdat[i]->GetName() );
      rat[i]=(TH1D*)effdat[i]->Clone( name.c_str() );
      rat[i]->Divide(effmc[i]);
   }
}

//--------------------------------------------------------------------
void plot_pict_K(int date, int pm = 0) {
//--------------------------------------------------------------------
   Params* par = new Params(date,1,pm,0); // date, K, "+/-",rew
   vector<TH1D*> eff_d, eff_mc, rat0;
   get_eff(par, 0, eff_d);
   get_eff(par, 1, eff_mc);
   get_ratio( eff_d, eff_mc, rat0 );

   int Nh = eff_d.size(); // must be 2

   // Attributes of draw
   auto setBlue = [](TH1* h) {
      h -> SetLineColor(kBlue+1);
      h -> SetMarkerColor(kBlue+2);
      h -> SetMarkerStyle(22);
      h -> SetMarkerSize(0.9);
   };
   auto setRed = [](TH1* h) {
      h -> SetLineColor(kRed+1);
      h -> SetMarkerColor(kRed+2);
      h -> SetMarkerStyle(23);
      h -> SetMarkerSize(0.9);
   };
   for (int i = 0; i < Nh; ++i ) {
      setBlue(eff_d[i]);
      setRed(eff_mc[i]);
   }

   // for fit
   TF1* pl0 = (TF1*)gROOT -> GetFunction("pol0")->Clone();
   pl0 -> SetLineColor(kRed);

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,400);
   c1 -> Divide(2,1);

///////////////////////////////////////////////////////////////////
   double eff_min = 0.2, eff_max = 1.0;
   double rat_min = 0.9, rat_max = 1.1;
///////////////////////////////////////////////////////////////////
   string p_pdf("K"), p_leg("K");
   if ( pm == 1 ) { p_pdf += "p"; p_leg += "^{+}";}
   if ( pm == -1 ){ p_pdf += "m"; p_leg += "^{-}";}

   TLegend* leg = new TLegend(0.65,0.40,0.89,0.60);
   leg->SetHeader( Form("%i  %s",date,p_leg.c_str()),"C");
   leg -> AddEntry(eff_d[0],  "  data ", "LP");
   leg -> AddEntry(eff_mc[0], "  MC ",   "LP");

   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.25);
   gStyle->SetFitFormat(".3f");

   // plot efficiencies
   for (int i = 0; i < Nh; ++i ) {
      c1 -> cd(i+1);
      gPad -> SetGrid();
      eff_d[i] -> SetAxisRange(eff_min,eff_max,"Y");
      string title =
         string( ((i==0) ? ";P_{t}, GeV/c" : ";cos(#Theta)") ) +
         string(";#epsilon(K)");
      eff_d[i] -> SetTitle(title.c_str());
      SetHstFace(eff_d[i]);
      eff_d[i] -> GetYaxis() -> SetTitleOffset(1.0);
      eff_d[i] -> Draw("E");
      eff_mc[i] -> Draw("E SAME");
      leg -> Draw();
   }

   gPad -> RedrawAxis();
   c1 -> Update();
   string pdf1 = string("eff_") + p_pdf + "_" + to_string(date) + ".pdf";
   c1 -> Print( pdf1.c_str() );
//    return;

   TCanvas* c2 = new TCanvas("c2","...",0,500,1200,400);
   c2 -> Divide(2,1);
   // plot ratio
   for (int i = 0; i < Nh; ++i ) {
      c2 -> cd(i+1);
      gPad -> SetGrid();
      rat0[i] -> SetAxisRange(rat_min,rat_max,"Y");
      string title = string(Form("%s, data/MC %i",p_leg.c_str(),date)) +
                     ( (i==0) ? ";P_{t}, GeV/c" : ";cos(#Theta)" ) +
                     string(";#epsilon(data) / #epsilon(MC)");
      rat0[i] -> SetTitle(title.c_str());
      SetHstFace(rat0[i]);
      rat0[i] -> GetYaxis() -> SetTitleOffset(1.15);
      rat0[i] -> SetLineWidth(2);
      rat0[i] -> SetMarkerStyle(20);
      rat0[i] -> SetLineColor(kBlack);
      rat0[i] -> Fit( pl0, "" );
   }

   gPad -> RedrawAxis();
   c2 -> Update();
   string pdf2 = string("rat_") + p_pdf + "_" + to_string(date) + ".pdf";
   c2 -> Print( pdf2.c_str() );
}

//--------------------------------------------------------------------
void plot_pict_pi(int date, int pm = 0) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,pm,0); // date, pi, "+/-",rew
   vector<TH1D*> eff_d, eff_mc, rat0;
   get_eff(par, 0, eff_d);
   get_eff(par, 1, eff_mc);
   get_ratio( eff_d, eff_mc, rat0 );

   int Nh = eff_d.size(); // must be 2

   // Attributes of draw
   auto setBlue = [](TH1* h) {
      h -> SetLineColor(kBlue+1);
      h -> SetMarkerColor(kBlue+2);
      h -> SetMarkerStyle(22);
      h -> SetMarkerSize(0.9);
   };
   auto setRed = [](TH1* h) {
      h -> SetLineColor(kRed+1);
      h -> SetMarkerColor(kRed+2);
      h -> SetMarkerStyle(23);
      h -> SetMarkerSize(0.9);
   };
   for (int i = 0; i < Nh; ++i ) {
      setBlue(eff_d[i]);
      setRed(eff_mc[i]);
   }

   // for fit
   TF1* pl0 = (TF1*)gROOT -> GetFunction("pol0")->Clone();
   pl0 -> SetLineColor(kRed);

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,400);
   c1 -> Divide(2,1);

///////////////////////////////////////////////////////////////////
   double eff_min = 0.2, eff_max = 1.0;
   double rat_min = 0.9, rat_max = 1.1;
///////////////////////////////////////////////////////////////////
   string p_pdf("Pi"), p_leg("#pi");
   if ( pm == 1 ) { p_pdf += "p"; p_leg += "^{+}";}
   if ( pm == -1 ){ p_pdf += "m"; p_leg += "^{-}";}

   TLegend* leg = new TLegend(0.65,0.40,0.89,0.60);
   leg->SetHeader( Form("%i  %s",date,p_leg.c_str()),"C");
   leg -> AddEntry(eff_d[0],  "  data ", "LP");
   leg -> AddEntry(eff_mc[0], "  MC ",   "LP");

   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.25);
   gStyle->SetFitFormat(".4f");

   // plot efficiencies
   for (int i = 0; i < Nh; ++i ) {
      c1 -> cd(i+1);
      gPad -> SetGrid();
      eff_d[i] -> SetAxisRange(eff_min,eff_max,"Y");
      string title =
         string( ((i==0) ? ";P_{t}, GeV/c" : ";cos(#Theta)") ) +
         string(";#epsilon(#pi)");
      eff_d[i] -> SetTitle(title.c_str());
      SetHstFace(eff_d[i]);
      eff_d[i] -> GetYaxis() -> SetTitleOffset(1.0);
      eff_d[i] -> Draw("E");
      eff_mc[i] -> Draw("E SAME");
      leg -> Draw();
   }

   gPad -> RedrawAxis();
   c1 -> Update();
   string pdf1 = string("eff_") + p_pdf + "_" + to_string(date) + ".pdf";
   c1 -> Print( pdf1.c_str() );
//    return;

   TCanvas* c2 = new TCanvas("c2","...",0,500,1200,400);
   c2 -> Divide(2,1);
   // plot ratio
   for (int i = 0; i < Nh; ++i ) {
      c2 -> cd(i+1);
      gPad -> SetGrid();
      rat0[i] -> SetAxisRange(rat_min,rat_max,"Y");
      string title = string(Form("%s, data/MC %i",p_leg.c_str(),date)) +
                     ( (i==0) ? ";P_{t}, GeV/c" : ";cos(#Theta)" ) +
                     string(";#epsilon(data) / #epsilon(MC)");
      rat0[i] -> SetTitle(title.c_str());
      SetHstFace(rat0[i]);
      rat0[i] -> GetYaxis() -> SetTitleOffset(1.15);
      rat0[i] -> SetLineWidth(2);
      rat0[i] -> SetMarkerStyle(20);
      rat0[i] -> SetLineColor(kBlack);
      rat0[i] -> Fit( pl0, "" );
   }

   gPad -> RedrawAxis();
   c2 -> Update();
   string pdf2 = string("rat_") + p_pdf + "_" + to_string(date) + ".pdf";
   c2 -> Print( pdf2.c_str() );
}

//-------------------------------------------------------------------------
void test_KS(int date, int kpi) {
//-------------------------------------------------------------------------
// use the Kolmogorov–Smirnov test to compare K+ and K- (pi+ and pi-)
// 1D distributions of the data/MC ratio
// see http://root.cern.ch/root/html/TH1.html#TH1:KolmogorovTest

   Params* parP = new Params(date,kpi,1,0); // date, K/pi, "+", rew
   vector<TH1D*> eff_d, eff_mc, ratP;
   get_eff(parP, 0, eff_d);
   get_eff(parP, 1, eff_mc);
   get_ratio( eff_d, eff_mc, ratP );

   Params* parM = new Params(date,kpi,-1,0); // date, K/pi, "-", rew
   vector<TH1D*> ratM;
   get_eff(parM, 0, eff_d);
   get_eff(parM, 1, eff_mc);
   get_ratio( eff_d, eff_mc, ratM );

   cout << " data-" << date << ": probability test for data/MC" << endl;
   for (int i = 0; i < 2; ++i ) {
      if ( kpi == 1 ) {
         cout << " K+ <--> K- ";
      } else {
         cout << " pi+ <--> pi- ";
      }
      cout << ( (i==0) ? ";P_{t} " : ";cos(#Theta) ");
//       cout << endl; ratP[i]->KolmogorovTest(ratM[i],"D"); // debug
//       cout << " Prob= " << ratP[i]->KolmogorovTest(ratM[i],"");
//       cout << " Prob(X)= " << ratP[i]->KolmogorovTest(ratM[i],"X");
      cout << " p-value(chi2)= " << ratP[i]->KolmogorovTest(ratM[i]);
      cout << endl;
   }
}

//--------------------------------------------------------------------
void Fill_Khst( Params* p, int mc, TH1D* hst[],
                double Ptmin, double Ptmax      ) {
//--------------------------------------------------------------------
   TTree* eff_K = p -> GetEff(mc);

   //Declaration of leaves types
   Double_t        Zk;
   Double_t        Ptk;
   Double_t        Ck;
   Double_t        fl;
   Double_t        dP;
   Double_t        dTh;
   Double_t        Mrec;
   Double_t        Mrk2;
   Double_t        Egsum;
   Double_t        Egmax;
   Double_t        good;

   // Set branch addresses.
   eff_K -> SetBranchAddress("Zk",&Zk);
   eff_K -> SetBranchAddress("Ptk",&Ptk);
   eff_K -> SetBranchAddress("Ck",&Ck);
   eff_K -> SetBranchAddress("fl",&fl);
   eff_K -> SetBranchAddress("dP",&dP);
   eff_K -> SetBranchAddress("dTh",&dTh);
   eff_K -> SetBranchAddress("Mrec",&Mrec);
   eff_K -> SetBranchAddress("Mrk2",&Mrk2);
   eff_K -> SetBranchAddress("Egsum",&Egsum);
   eff_K -> SetBranchAddress("Egmax",&Egmax);
   eff_K -> SetBranchAddress("good",&good);

   // List of cuts:
   auto c_kaon = [](double Mrec, double Mrk2) -> bool{
      return (fabs(Mrec-3.097)<0.002) && (fabs(Mrk2-0.243717)<0.025);
   };
//    auto c_KPtC = [](double Ptk, double Ck) -> bool {
//       return (0.1<Ptk && Ptk<1.4) && (fabs(Ck)<0.8);
//    };
   auto c_plus = [](double Z) -> bool{
      return Z>0;
   };
   auto c_Ptbin = [Ptmin,Ptmax](double Pt) -> bool {
      return (Ptmin<=Pt && Pt<Ptmax);
   };

   int rew = 0;
   if ( mc == 1 ) {
      rew = p -> use_rew;
   }

   Long64_t nentries = eff_K->GetEntries();
   for ( Long64_t i = 0; i < nentries; ++i ) {
      eff_K -> GetEntry(i);
      if ( !c_kaon(Mrec,Mrk2) ) continue;
      if ( p -> use_pm == 1  && !c_plus(Zk) ) continue;
      if ( p -> use_pm == -1  && c_plus(Zk) ) continue;
      if ( !c_Ptbin(Ptk) ) continue;

      // reweiting
      double W = 1;
      if ( rew == 1 ) {
         W = ReWeightTrkPid(p -> date, 1, Ptk);
      }

      hst[0] -> Fill(Ck);
      if ( fl > 1.5 ) {         // trk && pid id!
         hst[1] -> Fill(Ck,W);
      }
   }
}

//--------------------------------------------------------------------
void Fill_PIhst( Params* p, int mc, TH1D* hst[],
                double Ptmin, double Ptmax      ) {
//--------------------------------------------------------------------
   TTree* eff_Pi = p -> GetEff(mc);

   //Declaration of leaves types
   Double_t        Zpi;
   Double_t        Ptpi;
   Double_t        Cpi;
   Double_t        fl;
   Double_t        dP;
   Double_t        dTh;
   Double_t        MppKK;
   Double_t        Mrpi2;
   Double_t        Egsum;
   Double_t        Egmax;
   Double_t        good;

   // Set branch addresses.
   eff_Pi->SetBranchAddress("Zpi",&Zpi);
   eff_Pi->SetBranchAddress("Ptpi",&Ptpi);
   eff_Pi->SetBranchAddress("Cpi",&Cpi);
   eff_Pi->SetBranchAddress("fl",&fl);
   eff_Pi->SetBranchAddress("dP",&dP);
   eff_Pi->SetBranchAddress("dTh",&dTh);
   eff_Pi->SetBranchAddress("MppKK",&MppKK);
   eff_Pi->SetBranchAddress("Mrpi2",&Mrpi2);
   eff_Pi->SetBranchAddress("Egsum",&Egsum);
   eff_Pi->SetBranchAddress("Egmax",&Egmax);
   eff_Pi->SetBranchAddress("good",&good);

   // List of cuts:
   auto c_pion = [](double MppKK, double Mrpi2) -> bool{
      return (fabs(MppKK-3.096)<0.009) &&
             (fabs(Mrpi2-0.0194798)<0.025);
   };
//    auto c_PiPtC = [](double Ptpi, double Cpi) -> bool {
//       return (0.05<Ptpi && Ptpi<0.4) && (fabs(Cpi)<0.8);
//    };
   auto c_plus = [](double Z) -> bool{
      return Z>0;
   };
   auto c_Ptbin = [Ptmin,Ptmax](double Pt) -> bool {
      return (Ptmin<=Pt && Pt<Ptmax);
   };

   int rew = 0;
   if ( mc == 1 ) {
      rew = p -> use_rew;
   }

   Long64_t nentries = eff_Pi->GetEntries();
   for ( Long64_t i = 0; i < nentries; ++i ) {
      eff_Pi -> GetEntry(i);
      if ( !c_pion(MppKK,Mrpi2) ) continue;
      if ( p -> use_pm == 1  && !c_plus(Zpi) ) continue;
      if ( p -> use_pm == -1  && c_plus(Zpi) ) continue;
      if ( !c_Ptbin(Ptpi) ) continue;

      // reweiting
      double W = 1;
      if ( rew == 1 ) {
         W = ReWeightTrkPid(p -> date, 0, Ptpi);
      }

      hst[0] -> Fill(Cpi);
      if ( fl > 1.5 ) {         // trk && pid id!
         hst[1] -> Fill(Cpi,W);
      }
   }
}
//--------------------------------------------------------------------
TH1D* get_eff_cos(Params* p, int mc, double Ptmin, double Ptmax) {
//--------------------------------------------------------------------
   vector<string> hns = { "C0", "C1" };
   auto hC = [](string nm) {return (new TH1D(nm.c_str(),"",40,-1.,1.));};
//    auto hC = [](string nm) {return (new TH1D(nm.c_str(),"",20,-1.,1.));};

   vector<TH1D*> hst(2,nullptr);
   for ( int i = 0; i < 2; ++i ) {
      hst[i] = hC(hns[i]);
      hst[i] -> Sumw2(true);
   }

   // Fill
   if ( p -> use_kpi == 1 ) {   // Kaons
      Fill_Khst( p, mc, hst.data(), Ptmin, Ptmax);
   } else {                     // Pions
      Fill_PIhst( p, mc, hst.data(), Ptmin, Ptmax);
   }

   // calculate efficiency
   string name( ((mc==0) ? "data_" : "mc_") );
   name += to_string(int(Ptmin*100));
   TH1D* heff = (TH1D*) hst[0] -> Clone( name.c_str() );
   heff -> Divide(hst[1],hst[0],1.,1.,"B");

   delete hst[0];
   delete hst[1];

   return heff;
}

//-------------------------------------------------------------------------
void FitRatio(int date, int kpi, int pm=0, int rew=0) {
//-------------------------------------------------------------------------
   Params* par = new Params(date,kpi,pm,rew); // date, K/pi, "+/-", rew

   double Ptmin, Ptmax, Nbins; // binning for Pt
   if ( kpi == 1 ) {    // Kaons
      Ptmin = 0.1;
      Ptmax = 1.4;
      Nbins = 13;
//       Nbins = 26;
   } else {             // Pions
      Ptmin = 0.05;
      Ptmax = 0.4;
      Nbins = 7;
//       Nbins = 14;
   }

   // First fit the ratio as the function of "cos(Theta)"
   // function to the first fit
   TF1* fit1 = (TF1*)gROOT->GetFunction("pol0") -> Clone("fit1");
   fit1 -> SetLineColor(kRed);
   fit1 -> SetLineWidth(2);

   // Second fit the obtained parameters as the function of "Pt"
   TF1* fit2 = nullptr;
   if ( rew == 1 ) {
      fit2 = (TF1*)gROOT->GetFunction("pol0") -> Clone("fit2");
      fit2 -> SetRange(Ptmin,Ptmax);
   } else {
      if ( kpi == 1 ) {    // Kaons
         // functions to the second fit of kaons
         // nch - max order of Chebyshev polynomials
         int nch = 4; // 2012
         if ( date == 2009 ) {
            nch = 1;  // strait line
         }

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
         // fit2 = new TF1("fit2", "cheb4", Ptmin, Ptmax);
         fit2 = new TF1("fit2", Lchb, Ptmin, Ptmax,nch+1);
         if ( date == 2012 ) {
            fit2 -> SetParameters(1.77, -1.34, 0.80, -0.315, 0.057);
         }
      } else {             // Pions
         // functions to the second fit of pions
         auto fPi = [](const double* x, const double* p) -> double {
            double t = p[1]/x[0];
            return p[0] + t*t*t;
         };
         fit2 = new TF1("fit2", fPi, Ptmin, Ptmax, 2);
         fit2 -> SetParNames("a","b");
//          fit2 -> SetParameters(0.98, 0.02);
         fit2 -> SetParameters(1., 0.);
      }
   }
   fit2 -> SetLineColor(kRed);
   fit2 -> SetLineWidth(2);

   // First fit:
   TCanvas* c1 = new TCanvas("c1","...",850,0,700,1000);
   int Nc1 = 14;
   if ( kpi == 1 ) {    // Kaons
      c1 -> Divide(2,7);
   } else {
      Nc1 = 8;
      c1 -> Divide(2,4);
   }
   int Ic1 = 1;
   string p_pdf( ((kpi==1) ? "K" : "Pi") );
   if ( pm == 1 ) { p_pdf += "p"; }
   if ( pm == -1 ){ p_pdf += "m"; }

   vector<double> xx(Nbins), yy(Nbins);
   vector<double> ex(Nbins,1e-3), ey(Nbins);
   double dp = (Ptmax-Ptmin)/double(Nbins);
   for (int i = 0; i < Nbins; ++i ) {
      double ptmin = Ptmin + dp*i;
      double ptmax = ptmin + dp;
      double ptav = 0.5*(ptmin + ptmax);
      xx[i] = ptav;

      TH1D* eff_dat = get_eff_cos(par,0,ptmin,ptmax);
      TH1D* eff_mc  = get_eff_cos(par,1,ptmin,ptmax);
      string name = string("r_") + string( eff_dat -> GetName() );
      TH1D* rat = (TH1D*)eff_dat -> Clone( name.c_str() );
      rat -> Divide(eff_dat, eff_mc);

      c1 -> cd(Ic1);
      gPad -> SetGrid();
      rat -> SetAxisRange(0.8,1.2,"Y");
      rat -> GetYaxis() -> SetNdivisions(1004);
      rat -> GetYaxis() -> SetTitleOffset(0.9);
      rat -> SetTitle(Form(
               "<Pt> = %.3f GeV"
               "; cos(#Theta)"
               ";#epsilon(data) / #epsilon(MC)",
               ptav));
      SetHstFace(rat);
      rat -> Fit(fit1,"Q","",-0.799,0.799);
      yy[i] = fit1 -> GetParameter(0);
      ey[i] = fit1 -> GetParError(0);
      cout << " i= " << i << " pt= " << ptav
           << " r= " << yy[i] << " +/- " << ey[i] << endl;

      Ic1 += 1;
      if ( i == Nbins-1 ) { // clear all other pads
         for (; Ic1 <= Nc1; ++Ic1 ) {
             c1 -> cd(Ic1);
             gPad -> Clear();
         }
      }
      if ( Ic1 > Nc1 ) {
         Ic1 = 1;
         c1 -> Update();
         string pdf = string("rat_fit1_") + p_pdf + "_" +
            to_string(date) + ".pdf";
         if ( pm == 0 && rew == 0 ) {
            c1 -> Print( pdf.c_str() );
         }
//          gPad->WaitPrimitive(); // pause
      }
   } // end of for

   // second fit
   TCanvas* c2 = new TCanvas("c2","...",0,0,800,800);
   c2 -> cd();
   gPad -> SetGrid();
   gStyle->SetFitFormat(".4f");

   TGraphErrors* gX = new TGraphErrors(Nbins,xx.data(),yy.data(),
                                             ex.data(),ey.data() );
   gX->GetHistogram()->SetMaximum(1.1);
   gX->GetHistogram()->SetMinimum(0.9);
   gX->SetMarkerColor(kBlue);
   gX->SetMarkerStyle(21);
   gX->Draw("AP");
   gX->GetYaxis()->SetTitleOffset(1.4);
   gX->SetTitle(
         "; P_{t}, GeV/c"
         "; #epsilon(data)/#epsilon(MC)"
               );

   if ( kpi == 1 && date == 2009 && rew == 0 ) { // reduce fit range
      gX->Fit("fit2","EX0","", 0.2, Ptmax);
   } else {
      gX->Fit("fit2","REX0");
   }

   // --------------------------------------------
   TLegend* leg = new TLegend(0.20,0.82,0.50,0.89);
   string p_leg( ((kpi==1) ? "K" : "#pi") );
   if ( pm == 0 ) { p_leg += "^{#pm}"; }
   if ( pm == 1 ) { p_leg += "^{+}"; }
   if ( pm == -1 ){ p_leg += "^{-}"; }
   leg -> SetHeader( Form("%i  %s",date,p_leg.c_str()),"C" );
   leg -> Draw();

   gPad -> RedrawAxis();
   c2 -> Update();
   string pdf = string("rat_fit_") + p_pdf +
      ( (rew==1) ? "_w_" : "_" ) + to_string(date) + ".pdf";
   c2 -> Print( pdf.c_str() );
}

//-------------------------------------------------------------------------
void trk_eff_fit() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(62);
   gStyle->SetStatFont(62);
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
//    gStyle->SetStatW(0.25);

//    plot_pict_K(2012,1);
//    plot_pict_K(2012,-1);
//    plot_pict_K(2009,1);
//    plot_pict_K(2009,-1);

//    test_KS(2012,1); // Kaons
//    test_KS(2009,1); // Kaons

//    plot_pict_pi(2012,1);
//    plot_pict_pi(2012,-1);
//    plot_pict_pi(2009,1);
//    plot_pict_pi(2009,-1);

//    test_KS(2012,2); // Pions
//    test_KS(2009,2); // Pions

//    FitRatio(2012,1); // Kaons
//    FitRatio(2009,1); // Kaons
//    FitRatio(2012,2); // Pions
//    FitRatio(2009,2); // Pions


//    FitRatio(2012,1,1,1); // Kaons+ re-weighted
//    FitRatio(2012,1,-1,1); // Kaons- re-weighted
//    FitRatio(2012,1,0,1); // Kaons+/- re-weighted
//    FitRatio(2009,1,1,1); // Kaons+ re-weighted
//    FitRatio(2009,1,-1,1); // Kaons- re-weighted
//    FitRatio(2009,1,0,1); // Kaons+/- re-weighted

//    FitRatio(2012,2,1,1); // Pions+ re-weighted
//    FitRatio(2012,2,-1,1); // Pions- re-weighted
   FitRatio(2012,2,0,1); // Pions+/- re-weighted
//    FitRatio(2009,2,1,1); // Pions+ re-weighted
//    FitRatio(2009,2,-1,1); // Pions- re-weighted
//    FitRatio(2009,2,0,1); // Pions+/- re-weighted
}
