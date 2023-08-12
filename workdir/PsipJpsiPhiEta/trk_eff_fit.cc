// trk_eff_fit.cc - The reconstruction efficiency (data and MC)
// and their ratio (data/MC) for pions and kaons.
// Get corrections as function of P_t
// by fitting the ratio.
// Kolmogorov–Smirnov  && Chi2 probability tests to
// compare the ratios for K+ and K- (pi+ and pi-)

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
      Z->SetTitleFont(62);
      Z->SetTitleSize(0.04);
   }
}

// for very small hst
//--------------------------------------------------------------------
void SetHstFaceTbl(TH1* hst) {
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

// {{{1 RewTrk functions with HC
//--------------------------------------------------------------------
double RewTrkPi(int DataPeriod, double Pt, double Z) {
//--------------------------------------------------------------------
   // Corrections for the efficiency of reconstruction a pion having
   // transverse momentum Pt and sign Z.
   // The return value is the weight for the MC event.
   // v709, DelClonedTrk, helix corrections

   const double Ptmin = 0.05, Ptmax = 0.4;
   Pt = max( Ptmin, Pt );
   Pt = min( Ptmax, Pt );

   // 'spline-function' rewritten from ROOT-generated function
   auto Lsp = [](double x, const double* fY,
         const double* fB, const double* fC, const double* fD) {
      const int fNp = 7;
      const double fDelta = 0.050, fXmin = 0.075, fXmax = 0.375;
      // const double fX[7] = {
         // 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375 };

      // Equidistant knots
      int klow = int( (x-fXmin)/fDelta );
      klow = max(klow,0);
      klow = min(klow,fNp-1);

      // Evaluate now
      // double dx = x - fX[klow];
      double dx = x - (fXmin + klow * fDelta);
      return fY[klow] + dx*(fB[klow] + dx*(fC[klow] + dx*fD[klow]));
   };


   double W = 1.;
   if ( DataPeriod == 2021 ) {
      if ( Z > 0 ) {
         // W = Pip_2021_spline( Pt );
         static const double fY[7] = {
            0.981691, 0.987616, 0.985132, 0.978823, 0.981374,
            0.984894, 0.983846 };
         static const double fB[7] = {
            0.154886, 0.0457408, -0.131388, -0.0478018, 0.0971147,
            0.0236108, -0.0432394 };
         static const double fC[7] = {
            -5.55112e-16, -2.1829, -1.35968, 3.03141, -0.133076,
            -1.337, 0.05 };
         static const double fD[7] = {
            -14.5527, 5.48818, 29.2739, -21.0965, -8.02619,
            8.91336, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      } else {
         // W = Pim_2021_spline( Pt );
         static const double fY[7] = {
            0.986859, 0.987314, 0.985509, 0.989826, 0.987409,
            0.988992, 0.981344 };
         static const double fB[7] = {
            0.0332347, -0.0391776, 0.0424775, 0.0199988, -0.00851865,
            -0.0359787, -0.211467 };
         static const double fC[7] = {
            6.93889e-17, -1.44825, 3.08135, -3.53092, 2.96057,
            -3.50977, 0.05 };
         static const double fD[7] = {
            -9.65498, 30.1973, -44.0818, 43.2766, -43.1357,
            23.3985, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      }
   } else if ( DataPeriod == 2012 ) {
      if ( Z > 0 ) {
         // W = Pip_2012_spline( Pt );
         static const double fY[7] = {
            1.03169, 0.993442, 0.992989, 0.980393, 0.984639,
            0.982477, 0.987807 };
         static const double fB[7] = {
            -0.99227, -0.310335, -0.0883914, -0.119001, 0.0633865,
            -0.00948767, 0.164633 };
         static const double fC[7] = {
            4.44089e-15, 13.6387, -9.19985, 8.58765, -4.93989,
            3.48241, 0.05 };
         static const double fD[7] = {
            90.9248, -152.257, 118.583, -90.1836, 56.1487,
            -23.216, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      } else {
         // W = Pim_2012_spline( Pt );
         static const double fY[7] = {
            1.03121, 0.995748, 0.994734, 0.998636, 0.985478,
            0.993828, 0.988655 };
         static const double fB[7] = {
            -0.877667, -0.372377, 0.178602, -0.16878, -0.0588284,
            0.115605, -0.212984 };
         static const double fC[7] = {
            -2.22045e-15, 10.1058, 0.9138, -7.86144, 10.0605,
            -6.57179, 0.05 };
         static const double fD[7] = {
            67.3719, -61.2799, -58.5016, 119.479, -110.882,
            43.8119, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      }
   } else if ( DataPeriod == 2009 ) {
      if ( Z > 0 ) {
         // W = Pip_2009_spline( Pt );
         static const double fY[7] = {
            0.992996, 0.99419, 0.990618, 0.982758, 0.986737,
            0.981416, 0.995198 };
         static const double fB[7] = {
            0.0372614, -0.00290172, -0.168377, -0.00952206, -0.0263622,
            0.0344359, 0.396241 };
         static const double fC[7] = {
            -2.77556e-16, -0.803263, -2.50624, 5.68333, -6.02013,
            7.2361, 0.05 };
         static const double fD[7] = {
            -5.35509, -11.3532, 54.5971, -78.0231, 88.3749,
            -48.2407, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      } else {
         // W = Pim_2009_spline( Pt );
         static const double fY[7] = {
            1.00415, 0.975639, 0.983931, 0.995963, 0.99217,
            0.991772, 1.00196 };
         static const double fB[7] = {
            -0.755805, -0.198795, 0.338097, 0.0658337, -0.107085,
            0.11107, 0.249998 };
         static const double fC[7] = {
            2.22045e-15, 11.1402, -0.402369, -5.04289, 1.58452,
            2.77857, 0.05 };
         static const double fD[7] = {
            74.268, -76.9504, -30.9368, 44.1827, 7.96038,
            -18.5238, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      }
   }
   return W;
}

//--------------------------------------------------------------------
double RewTrk_K(int DataPeriod, double Pt, double Z) {
//--------------------------------------------------------------------
   // Corrections for the efficiency of reconstruction a kaon having
   // transverse momentum Pt and sign Z.
   // The return value is the weight for the MC event.
   // v709, DelClonedTrk, helix corrections

   const double Ptmin = 0.1, Ptmax = 1.4;
   Pt = max( Ptmin, Pt );
   Pt = min( Ptmax, Pt );

   // [Ptmin,Ptmax]->[-1,+1]
   double x = (2*Pt-Ptmin-Ptmax) / (Ptmax-Ptmin);
   auto Lchb = [](int nch, double x, const vector<double> p) {
      if ( nch == 0 ) {
         return p[0];
      }
      double sum = p[0] + x*p[1];
      if ( nch == 1 ) {
         return sum;
      }
      double T0 = 1, T1 = x;
      for ( int i = 2; i <= nch; ++i ) {
         double Tmp = 2*x*T1 - T0;
         sum += p[i]*Tmp;
         T0 = T1;
         T1 = Tmp;
      }
      return sum;
   };

   double W = 1.;
   if ( DataPeriod == 2009 ) {
      if ( Pt > 0.2 ) {
         if ( Z > 0 ) {
            W = Lchb(1, x, {0.9902,-0.013});
         } else {
            W = Lchb(1, x, {0.9930,-0.023});
         }
      } else {
         W = 0.92;
      }
   } else if ( DataPeriod == 2012 ) {
      if ( Z > 0 ) {
         if ( Pt > 0.3 ) {
            W = Lchb(3, x, {0.9677,0.022,-0.018,0.018});
         } else {
            W = ( Pt > 0.2 ) ? 0.97 : 1.01;
         }
      } else {
         if ( Pt > 0.3 ) {
            W = Lchb(3, x, {0.9723,0.028,-0.016,0.021});
         } else {
            W = ( Pt > 0.2 ) ? 1.01 : 0.99;
         }
      }
   } else if ( DataPeriod == 2021 ) {
      if ( Z > 0 ) {
         if ( Pt > 0.3 ) {
            W = Lchb(5,x,{0.9402,0.074,-0.056,0.043,-0.018,0.005});
         } else {
            W = ( Pt > 0.2 ) ? 0.966 : 0.953;
         }
      } else {
         if ( Pt > 0.3 ) {
            W = Lchb(5,x,{0.9524,0.043,-0.038,0.029,-0.007,0.001});
         } else {
            W = ( Pt > 0.2 ) ? 0.996 : 1.004;
         }
      }
   }
   return W;
}

// {{{1 RewTrk functions noHC
//--------------------------------------------------------------------
double RewTrkPi0(int DataPeriod, double Pt, double Z) {
//--------------------------------------------------------------------
   // Corrections for the efficiency of reconstruction a pion having
   // transverse momentum Pt and sign Z.
   // The return value is the weight for the MC event.
   // v709, DelClonedTrk, NO helix corrections

   const double Ptmin = 0.05, Ptmax = 0.4;
   Pt = max( Ptmin, Pt );
   Pt = min( Ptmax, Pt );

   double W = 1.;
   if ( DataPeriod == 2009 ) {
      if ( Z > 0 ) {
         W = 0.991 - 0.017 * Pt;
      } else {
         W = 0.975 + 0.064 * Pt;
      }
   } else if ( DataPeriod == 2012 ) {
      auto SQ = [](double x) -> double{return x*x;};
      if ( Z > 0 ) {
         W = 0.9806 + SQ(0.0150/Pt);
      } else {
         W = 0.9891 + SQ(0.0131/Pt);
      }
   } else if ( DataPeriod == 2021 ) {
      if ( Z > 0 ) {
         W = 0.9825;
      } else {
         W = 0.9874;
      }
   }
   return W;
}

// {{{1 Fill histograms
//--------------------------------------------------------------------
void get_hst(Params* p, int mc, vector<TH1D*>& hst) {
//--------------------------------------------------------------------
// get histograms: mc=0,1 / data,MC

   TTree* eff = p->GetEff(mc);

   TCut Crec;
   if ( p->use_kpi == 1 ) {   // Kaons
      Crec += p->Ckaon + p->CKPtC;
      if ( p->use_pm == 1 ) {         // "K+"
         Crec += TCut("Zk>0");
      } else if ( p->use_pm == -1 ) { // "K-"
         Crec += TCut("Zk<0");
      }
   } else {                     // Pions
      Crec += p->Cpion + p->CPiPtC;
      if ( p->use_pm == 1 ) {         // "pi+"
         Crec += TCut("Zpi>0");
      } else if ( p->use_pm == -1 ) { // "pi-"
         Crec += TCut("Zpi<0");
      }
   }
   // cout << " DEBUG: Crec= " << Crec.GetTitle() << endl;
   TCut Crec0 = Crec + TCut("fl<1.5");
   TCut Crec1 = Crec + TCut("fl>1.5"); // trk && pid id!

   string hname( ((mc==0) ? "data" : "mc") );
   hname += to_string(p->use_pm); // for KS-test
   vector<string> hns = {
      hname + "_P0", hname + "_P1",
      hname + "_C0", hname + "_C1"
   };
   auto hPtK = [](string nm) {
      return new TH1D(nm.c_str(),"",26,0.1,1.4);
   };
   auto hPtPi= [](string nm) {
      return new TH1D(nm.c_str(),"",14,0.05,0.4);
   };
   auto hC = [](string nm) {
      return (new TH1D(nm.c_str(),"",32,-0.8,0.8));
   };

   hst.clear();
   hst.resize(4,nullptr);
   if ( p->use_kpi == 1 ) {   // Kaons
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
      h->Sumw2(true);
   }

   // Fill
   if ( p->use_kpi == 1 ) {   // Kaons
      eff->Draw( ("Ptk>>"+hns[0]).c_str(),Crec0,"goff");
      eff->Draw( ("Ck>>"+hns[2]).c_str(),Crec0,"goff");
      eff->Draw( ("Ptk>>"+hns[1]).c_str(),Crec1,"goff");
      eff->Draw( ("Ck>>"+hns[3]).c_str(),Crec1,"goff");
   } else {                     // Pions
      eff->Draw( ("Ptpi>>"+hns[0]).c_str(),Crec0,"goff");
      eff->Draw( ("Cpi>>"+hns[2]).c_str(),Crec0,"goff");
      eff->Draw( ("Ptpi>>"+hns[1]).c_str(),Crec1,"goff");
      eff->Draw( ("Cpi>>"+hns[3]).c_str(),Crec1,"goff");
   }
}

// {{{1 Efficiency calculation
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

// {{{1 Ratio calculation
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

// {{{1 Plot K
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
   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,400);
   c1->Divide(2,1);

   ////////////////////////////////////
   double eff_min = 0.2, eff_max = 1.0;
   double rat_min = 0.9, rat_max = 1.1;
   ////////////////////////////////////
   string p_pdf("K"), p_leg("K");
   if ( pm == 1 ) { p_pdf += "p"; p_leg += "^{#plus}";}
   if ( pm == -1 ){ p_pdf += "m"; p_leg += "^{#minus}";}

   TLegend* leg = new TLegend(0.65,0.30,0.89,0.60);
   leg->SetHeader( Form("%i  %s",date,p_leg.c_str()),"C");
   leg->AddEntry(eff_d[0],  "  data ", "LP");
   leg->AddEntry(eff_mc[0], "  MC ",   "LP");

   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.25);
   gStyle->SetFitFormat(".3f");

   // plot efficiencies
   for (int i = 0; i < Nh; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();
      eff_d[i]->SetAxisRange(eff_min,eff_max,"Y");
      string title =
         string( ((i==0) ? ";P_{t}, GeV/c" : ";cos(#Theta)") ) +
         string(";#epsilon(K)");
      eff_d[i]->SetTitle(title.c_str());
      SetHstFace(eff_d[i]);
      eff_d[i]->GetYaxis()->SetTitleOffset(0.8);
      eff_d[i]->GetXaxis()->SetTitleOffset(0.9);
      eff_d[i]->Draw("E");
      eff_mc[i]->Draw("E SAME");
      leg->Draw();
   }

   gPad->RedrawAxis();
   c1->Update();
   string pdf1 = "eff_" + p_pdf + "_" + to_string(date) + ".pdf";
   c1->Print( pdf1.c_str() );

   TCanvas* c2 = new TCanvas("c2","...",0,500,1200,400);
   c2->Divide(2,1);
   // plot ratio
   for (int i = 0; i < Nh; ++i ) {
      c2->cd(i+1);
      gPad->SetGrid();
      rat0[i]->SetAxisRange(rat_min,rat_max,"Y");
      string title =
         string( Form("%s, data/MC %i",p_leg.c_str(),date) )
         + ( (i==0) ? ";P_{t}, GeV/c" : ";cos(#Theta)" )
         + string(";#epsilon(data) / #epsilon(MC)");
      rat0[i]->SetTitle(title.c_str());
      SetHstFace(rat0[i]);
      rat0[i]->GetYaxis()->SetTitleOffset(0.9);
      rat0[i]->GetXaxis()->SetTitleOffset(0.9);
      rat0[i]->SetLineWidth(2);
      rat0[i]->SetMarkerStyle(20);
      rat0[i]->SetLineColor(kBlack);
      rat0[i]->Fit( pl0, "" );
      rat0[i]->Draw("SAME E0");
   }

   gPad->RedrawAxis();
   c2->Update();
   string pdf2 = "rat_" + p_pdf + "_" + to_string(date) + ".pdf";
   c2->Print( pdf2.c_str() );
}

// {{{1 Plot Pi
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
   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,400);
   c1->Divide(2,1);

   ////////////////////////////////////
   double eff_min = 0.2, eff_max = 1.0;
   double rat_min = 0.9, rat_max = 1.1;
   ////////////////////////////////////
   string p_pdf("Pi"), p_leg("#pi");
   if ( pm == 1 ) { p_pdf += "p"; p_leg += "^{#kern[0.25]{#plus}}";}
   if ( pm == -1 ){ p_pdf += "m"; p_leg += "^{#kern[0.25]{#minus}}";}

   TLegend* leg = new TLegend(0.65,0.30,0.89,0.60);
   leg->SetHeader( Form("%i  %s",date,p_leg.c_str()),"C");
   leg->AddEntry(eff_d[0],  "  data ", "LP");
   leg->AddEntry(eff_mc[0], "  MC ",   "LP");

   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.25);
   gStyle->SetFitFormat(".4f");

   // plot efficiencies
   for (int i = 0; i < Nh; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();
      eff_d[i]->SetAxisRange(eff_min,eff_max,"Y");
      string title =
         string( ((i==0) ? ";P_{t}, GeV/c" : ";cos(#Theta)") ) +
         string(";#epsilon(#pi)");
      eff_d[i]->SetTitle(title.c_str());
      SetHstFace(eff_d[i]);
      eff_d[i]->GetYaxis()->SetTitleOffset(0.8);
      eff_d[i]->GetXaxis()->SetTitleOffset(0.9);
      eff_d[i]->Draw("E");
      eff_mc[i]->Draw("E SAME");
      leg->Draw();
   }

   gPad->RedrawAxis();
   c1->Update();
   string pdf1 = "eff_" + p_pdf + "_" + to_string(date) + ".pdf";
   c1->Print( pdf1.c_str() );
   // return;

   TCanvas* c2 = new TCanvas("c2","...",0,500,1200,400);
   c2->Divide(2,1);
   // plot ratio
   for (int i = 0; i < Nh; ++i ) {
      c2->cd(i+1);
      gPad->SetGrid();
      rat0[i]->SetAxisRange(rat_min,rat_max,"Y");
      // if ( date == 2012 && i == 0 ) {
         // rat0[i]->SetAxisRange(0.95,1.25,"Y");
      // }
      string title =
         string(Form("%s, data/MC %i",p_leg.c_str(),date)) +
         ( (i==0) ? ";P_{t}, GeV/c" : ";cos(#Theta)" ) +
         string(";#epsilon(data) / #epsilon(MC)");
      rat0[i]->SetTitle(title.c_str());
      SetHstFace(rat0[i]);
      rat0[i]->GetYaxis()->SetTitleOffset(0.9);
      rat0[i]->GetXaxis()->SetTitleOffset(0.9);
      rat0[i]->SetLineWidth(2);
      rat0[i]->SetMarkerStyle(20);
      rat0[i]->SetLineColor(kBlack);
      rat0[i]->Fit( pl0, "" );
      rat0[i]->Draw("SAME E0");
   }

   gPad->RedrawAxis();
   c2->Update();
   string pdf2 = "rat_" + p_pdf + "_" + to_string(date) + ".pdf";
   c2->Print( pdf2.c_str() );
}

// {{{1 Test +/- identity
//--------------------------------------------------------------------
void test_PM(int date, int kpi, int test=1) {
//--------------------------------------------------------------------
// compare K+ and K- (pi+ and pi-) 1D distributions of the data/MC
// ratio to check their identity

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

   cout << " data-" << date << ": ";
   if ( test == 1 ) {
      cout << "Kolmogorov–Smirnov";
   } else if ( test == 2 ) {
      cout << "chi^2";
   }
   cout << " probability test for data/MC efficiency ratio" << endl;

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

// {{{1 Fill K-histograms
//--------------------------------------------------------------------
void Fill_Khst( Params* p, int mc, TH1D* hst[],
                double Ptmin, double Ptmax      ) {
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
   // auto c_KPtC = [](double Ptk, double Ck)->bool {
      // return (0.1<Ptk && Ptk<1.4) && (fabs(Ck)<0.8);
   // };
   auto c_Ptbin = [Ptmin,Ptmax](double Pt)->bool {
      return (Ptmin<=Pt && Pt<Ptmax);
   };

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

      if ( p->use_pm == 1  && Zk < 0 ) {
         continue;
      }
      if ( p->use_pm == -1  && Zk > 0 ) {
         continue;
      }

      // reweiting
      double W = 1;
      if ( rew == 1 ) {
         W = RewTrk_K(p->date, Ptk, Zk);
      }

      hst[0]->Fill(Ck);
      if ( fl > 1.5 ) {         // trk && pid id!
         hst[1]->Fill(Ck,W);
      }
   }
}

// {{{1 Fill Pi-histograms
//--------------------------------------------------------------------
void Fill_PIhst( Params* p, int mc, TH1D* hst[],
                double Ptmin, double Ptmax      ) {
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
      return (fabs(MppKK-3.096)<0.009) &&
             (fabs(Mrpi2-0.0194798)<0.025);
   };
   // auto c_PiPtC = [](double Ptpi, double Cpi)->bool {
      // return (0.05<Ptpi && Ptpi<0.4) && (fabs(Cpi)<0.8);
   // };
   auto c_Ptbin = [Ptmin,Ptmax](double Pt)->bool {
      return (Ptmin<=Pt && Pt<Ptmax);
   };

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

      if ( p->use_pm == 1  && Zpi < 0 ) {
         continue;
      }
      if ( p->use_pm == -1  && Zpi > 0 ) {
         continue;
      }

      // reweiting
      double W = 1;
      if ( rew == 1 ) {
         W = RewTrkPi(p->date, Ptpi, Zpi);
      }

      hst[0]->Fill(Cpi);
      if ( fl > 1.5 ) {         // trk && pid id!
         hst[1]->Fill(Cpi,W);
      }
   }
}

// {{{1 Get efficiency as a finction of cos(Theta)
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
      Fill_Khst( p, mc, hst.data(), Ptmin, Ptmax);
   } else {                   // Pions
      Fill_PIhst( p, mc, hst.data(), Ptmin, Ptmax);
   }

   // calculate efficiency
   string name( ((mc==0) ? "data_" : "mc_") );
   name += to_string(int(Ptmin*100));
   TH1D* heff = (TH1D*) hst[0]->Clone( name.c_str() );
   heff->Divide(hst[1],hst[0],1.,1.,"B");

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
      Ptmin = 0.1;
      Ptmax = 1.4;
      Nbins = 13;
   } else if ( kpi == 2 ) { // Pions
      Ptmin = 0.05;
      Ptmax = 0.4;
      Nbins = 7;
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
   } else {
      if ( kpi == 1 ) {
         // functions to the second fit of Kaons
         int nch = 0; // nch is the order of Chebyshev polynomials
         if ( date == 2009 ) {
            nch = 1;  // straight line
         } else if ( date == 2012 ) {
            nch = 3;
         } else if ( date == 2021 ) {
            nch = 5;
         }
         double xmin = Ptmin, xmax = Ptmax;
         auto Lchb =
            [nch,xmin,xmax](const double* xx, const double* p) {
               if ( nch == 0 ) {
                  return p[0];
               }
               // [xmin,xmax]->[-1,+1]
               double x = (2*xx[0]-xmin-xmax)/(xmax-xmin);
               double sum = p[0] + x*p[1];
               if ( nch == 1 ) {
                  return sum;
               }
               double T0 = 1, T1 = x;
               for ( int i = 2; i <= nch; ++i ) {
                  double Tmp = 2*x*T1 - T0;
                  sum += p[i]*Tmp;
                  T0 = T1;
                  T1 = Tmp;
               }
               return sum;
            };

         fit2 = new TF1("fit2", Lchb, Ptmin, Ptmax,nch+1);
         for ( size_t ipar = 0; ipar <= nch; ++ipar ) {
            double ini_val = ( ipar==0 ) ? 0.99 : 0.01 ;
            fit2->SetParameter(ipar, ini_val);
         }
      }
   }

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
   c1->Print( (pdf1+"[").c_str() ); // open pdf-file

   vector<double> xx(Nbins), yy(Nbins);
   vector<double> ex(Nbins,1e-3), ey(Nbins);
   double dp = (Ptmax-Ptmin)/double(Nbins);
   for (int i = 0; i < Nbins; ++i ) {
      double ptmin = Ptmin + dp*i;
      double ptmax = ptmin + dp;
      double ptav = 0.5*(ptmin + ptmax);
      xx[i] = ptav;

      TH1D* eff_dat = get_eff_cos(par,0,ptmin,ptmax);
      double Ndat = eff_dat->GetEntries();
      TH1D* eff_mc  = get_eff_cos(par,1,ptmin,ptmax);
      double Nmc = eff_mc->GetEntries();
      string name = string("r_") + string( eff_dat->GetName() );
      TH1D* rat = (TH1D*)eff_dat->Clone( name.c_str() );
      rat->Divide(eff_dat, eff_mc);

      c1->cd(Ic1);
      gPad->SetGrid();
      rat->SetAxisRange(0.8,1.2,"Y");
      if ( i == 0 || i == Nbins-1 ) {
         rat->SetAxisRange(0.5,1.5,"Y");
      }
      rat->GetYaxis()->SetNdivisions(1004);
      rat->SetTitle(Form(
               "<Pt> = %.3f GeV"
               "; cos(#Theta)"
               ";#epsilon(data) / #epsilon(MC)",
               ptav));
      SetHstFaceTbl(rat);
      if ( kpi == 1 ) {
         gStyle->SetTitleFontSize(0.1);
         rat->GetYaxis()->SetTitleOffset(0.5);
      } else {
         gStyle->SetTitleFontSize(0.08);
         rat->GetYaxis()->SetTitleOffset(0.8);
      }
      rat->GetXaxis()->SetTitleOffset(0.8);
      // rat->Fit(fit1,"Q","",-0.799,0.799);
      TFitResultPtr rs = rat->Fit(fit1,"SEQ","",-0.8,0.8);
      rat->Draw("SAME E0");
      // rs->Print();
      yy[i] = fit1->GetParameter(0);
      ey[i] = fit1->GetParError(0);
      printf("%2i  <pt>=%.3f  r=%.4f+/-%.4f  N(dat,mc)=%.1f %.1f\n",
            i,ptav,yy[i],ey[i],Ndat,Nmc);

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
         c1->Print( pdf1.c_str() ); // add to pdf-file
      }
   } // end of for
   c1->Print( (pdf1+"]").c_str() ); // close pdf-file

   // ++++ Second fit ++++
   TCanvas* c2 = new TCanvas("c2","...",0,0,800,800);
   c2->cd();
   gPad->SetGrid();
   gStyle->SetFitFormat(".4f");

   TGraphErrors* gX = new TGraphErrors(Nbins,xx.data(),yy.data(),
                                             ex.data(),ey.data() );
   if ( rew == 0 ) {
      gX->GetHistogram()->SetMaximum(1.1);
      gX->GetHistogram()->SetMinimum(0.9);
   } else {
      if ( kpi == 1 ) {
         gX->GetHistogram()->SetMaximum(1.05);
         gX->GetHistogram()->SetMinimum(0.95);
      } else {
         gX->GetHistogram()->SetMaximum(1.025);
         gX->GetHistogram()->SetMinimum(0.975);
      }
      gX->GetYaxis()->SetNdivisions(1005);
   }
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
      // sp = new TSpline3("sp",gX,"e1b1",0.,0.);
      sp = new TSpline3("sp",gX,"e2b2",0.,0.);// natural cubic spline
      string fsp = p_pdf+"_"+to_string(date)+"_spline.cc";
      printf("Save spline to file: %s\n", fsp.c_str());
      sp->SaveAs(fsp.c_str());
      // to draw within [Ptmin,Ptmax]
      auto Lsp = [sp](double* x, double* p) {
         return sp->Eval(x[0]);
      };
      Fsp = new TF1("Fsp",Lsp,Ptmin,Ptmax,0);
      Fsp->SetLineColor(kRed);
      Fsp->SetLineWidth(2);
      Fsp->Draw("L SAME");
      gX->Draw("SAME 0P");
   }

   string p_leg( ((kpi==1) ? "K" : "#pi") );
   if ( pm == 0 ) { p_leg += "^{#kern[0.25]{#pm}}"; }
   if ( pm == 1 ) { p_leg += "^{#kern[0.25]{#plus}}"; }
   if ( pm == -1 ){ p_leg += "^{#kern[0.25]{#minus}}"; }
   TLegend* leg = nullptr;
   if ( fit2 ) {
      leg = new TLegend(0.20,0.82,0.50,0.89);
      leg->SetHeader( Form("%i  %s",date,p_leg.c_str()),"C" );
   } else {
      leg = new TLegend(0.55,0.77,0.89,0.89);
      leg->AddEntry(gX,Form("%i  %s",date,p_leg.c_str()),"PE");
      leg->AddEntry(Fsp,"Cubic Spline","L");
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

   // plot_pict_pi(2009,1);  // +
   // plot_pict_pi(2009,-1); // -
   // plot_pict_pi(2012,1);
   // plot_pict_pi(2012,-1);
   // plot_pict_pi(2021,1);
   // plot_pict_pi(2021,-1);

   // plot_pict_K(2009,1);  // +
   // plot_pict_K(2009,-1); // -
   // plot_pict_K(2012,1);
   // plot_pict_K(2012,-1);
   // plot_pict_K(2021,1);
   // plot_pict_K(2021,-1);

   // int test = 1; // 1 - Kolmogorov–Smirnov; 2 - Chi2Test
   // for ( auto date : {2009, 2012, 2021} ) {
      // test_PM(date,1,test); // Kaons
      // test_PM(date,2,test); // Pions
   // }

   // FitRatio(2009,2,1); // 2009, Pions, (1=+, -1=-, 0=+/-)
   // FitRatio(2009,1,1); // 2009, Kaons
   // FitRatio(2009,2,1,1); // re-weighted

   // FitRatio(2012,2,1); // 2012, Pions, (1=+, -1=-, 0=+/-)
   // FitRatio(2012,1,1); // 2012, Kaons
   // FitRatio(2012,2,1,1); // re-weighted

   // FitRatio(2021,2,1); // 2012, Pions, (1=+, -1=-, 0=+/-)
   // FitRatio(2021,1,1); // 2021, Kaons,
   // FitRatio(2021,2,1,1); // re-weighted

}
