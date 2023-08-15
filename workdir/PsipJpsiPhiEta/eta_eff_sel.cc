// eta_eff_sel.cc - Pictures for presentation/memo
// Study of the eta->2gamma reconstruction efficiency

// {{{1 Common parameters: Params
//--------------------------------------------------------------------
struct Params {
   Params(int dat, int slc, int rew);

   TFile* OpenFile(int mc);
   TTree* GetEffEta(int mc);
   double W_g_eta();    // wights for MC-gamma-eta2
   double W_phi_eta();  // wights for MC-phi-eta2
   const char* Sdate() { return Form("%i",date); }

   // name of folder with root files
   const string Dir = "prod_v709n3/";
   string datafile;
   string mcincfile;
   string mcsigf1;  // MC for gamma eta
   string mcsigf2;  // MC for phi eta

   int date;
   int slct;    // >0: g-eta; <0: phi-eta
                // 1 - photon eff; 2 - eta eff;
   int use_rew; // 0 - no weights;
                // 1 - calculate weights

   TCut Cmcsig; // for mc-signal
   TCut Cbg;    // cuts against the background
   TCut Cph;    // cuts for selection photons
   TCut Ceta;   // cuts for selection eta

   const double Br_eta_2gamma = 0.3936; // PDG 2023
   const double Br_phi_KK = 0.491; // PDG 2023
};

// {{{2 > ctor
//--------------------------------------------------------------------
Params::Params(int dat, int slc = 2, int rew = 0) {
//--------------------------------------------------------------------
   date = dat;
   slct = slc;
   use_rew = rew;

   // set the names:
   datafile  = string( Form("data_%02ipsip_all.root",date%100) );
   mcincfile = string( Form("mcinc_%02ipsip_all.root",date%100) );
   // mcsigf1 = string( Form("mcgammaeta2_kkmc_%02i.root",date%100) );
   // mcsigf2 = string( Form("mcphieta2_kkmc_%02i.root",date%100) );
   mcsigf1 = mcincfile; // debug
   mcsigf2 = string( Form("mcsig_kkmc_%02i.root",date%100) ); // debug

   // mc-signal
   if ( slct > 0 ) {    // gamma-eta
      Cmcsig += TCut("decj==22");
   } else {             // phi-eta
      Cmcsig += TCut("decj==68");
   }

   // cuts against the background
   if ( slct > 0 ) {         // gamma-eta
      Cbg += TCut("m2fr<0.002");
      Cbg += TCut("fabs(Cg0)<0.8");
      Cbg += TCut("fabs(Cg1)<0.8 && Eg1>0.2");
   } else {                  // phi-eta
      Cbg += TCut("m2fr<0.002");
   }

   // cuts for selection photons
   Cph += TCut("fabs(Cg2)<0.8||(fabs(Cg2)>0.85&&fabs(Cg2)<0.92)");
   Cph += TCut("fabs(Cg1)<0.8||(fabs(Cg1)>0.85&&fabs(Cg1)<0.92)");
   if ( slct > 0 ) {         // gamma-eta
      Cph += TCut("Eg2>0.1&&Eg2<1.4");
   } else {                  // phi-eta
      Cph += TCut("Eg2>0.05&&Eg2<1.45");
   }

   // cuts for selection eta
   Ceta += TCut("fabs(Ceta)<0.9");
   Ceta += TCut("fabs(Cg2)<0.8||(fabs(Cg2)>0.85&&fabs(Cg2)<0.92)");
   if ( slct > 0 ) {         // gamma-eta
      Ceta += TCut("Peta>1.3&&Peta<1.7");
   } else {                  // phi-eta
      Ceta += TCut("Peta>1.15&&Peta<1.5");
   }
}

// {{{2 > OpenFile(mc = 0: data, 1: MC inc, 2: MC-signal)
//--------------------------------------------------------------------
TFile* Params::OpenFile(int mc) {
//--------------------------------------------------------------------
   string dfname = Dir;
   if ( mc == 0 ) {
      dfname += datafile;
   } else if ( mc == 1 ) {
      dfname += mcincfile;
   } else if ( mc == 2 ) {
      if ( slct > 0 ) {
         dfname += mcsigf1; // gamma-eta-2
      } else {
         dfname += mcsigf2; // phi-eta-2
      }
   }
   // cout " OpenFile: " << dfname << endl;
   TFile* froot = TFile::Open(dfname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << dfname << endl;
      exit(EXIT_FAILURE);
   }
   if ( slct > 0 ) {
      froot->cd("PsipJpsiGammaEta");
   } else {
      froot->cd("PsipJpsiPhiEta");
   }
   return froot;
}

// {{{2 > GetEffEta : root-tree for eta
//--------------------------------------------------------------------
TTree* Params::GetEffEta(int mc = 0) {
//--------------------------------------------------------------------
   TFile* froot = this->OpenFile(mc);
   TTree* eff_eta = (TTree*)gDirectory->Get("eff_eta");
   if ( !eff_eta ) {
      cout << "can not find eff_eta in " << froot->GetName() << endl;
      exit(EXIT_FAILURE);
   }
   return eff_eta;
}

// {{{2 > W_g_eta() & W_phi_eta()
//--------------------------------------------------------------------
double Params::W_g_eta() { // wights for MC-gamma-eta
//--------------------------------------------------------------------
   // normalization on numbers in "official inclusive MC"
   // 664: ((date==2012) ? (137258./5e5) : (35578./1.5e5))
   double W = 1;
   // switch (date) {
      // case 2009: W =  41553./150e3; break;
      // case 2012: W = 130865./500e3; break;
      // case 2021: W = 883676./2.5e6; break;
   // }
   return W;
}

//--------------------------------------------------------------------
double Params::W_phi_eta() { // wights for MC-sig
//--------------------------------------------------------------------
   // normalization on numbers in "official inclusive MC"
   // 664: ((date==2012) ? (104950./5e5) : (27274./1.5e5));
   double W = 1;
   switch (date) {
      case 2009: W =  27547./200e3; break;
      case 2012: W =  89059./600e3; break;
      case 2021: W = 598360./3.0e6; break;
   }
   return W * Br_phi_KK;
}

// {{{1 helper functions and constants
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

// {{{1 Fill histograms
//--------------------------------------------------------------------
void get_hst( Params* p, int mc,
      const vector<string> hns, vector<TH1D*>& hst ) {
//--------------------------------------------------------------------
// read set of histograms to "hst" according of names "hns"

   int Nh = hns.size();
   hst.clear();
   hst.resize(Nh,nullptr);

   TFile* froot = p->OpenFile(mc);

   for ( int i = 0; i < Nh; ++i ) {
      hst[i] = (TH1D*)gROOT->FindObject(hns[i].c_str());
      if ( !hst[i] ) {
         cout << " can not find histo:" << hns[i] << endl;
         exit(EXIT_FAILURE);
      }
      hst[i]->Sumw2(true);
      SetHstFace(hst[i]);

      // rename:
      string name( Form("%s_%02i_%s", (mc==1) ? "MC_" : "Dat",
               p->date, hst[i]->GetName()) );
      hst[i]->SetName( name.c_str() );
   }
}

// {{{1 Plot pi0
//--------------------------------------------------------------------
void plot_pi0(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   vector<string> hndat {"Mg2_pi0"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   // TE: true eta->2gamma
   vector<string> hnmc {"Mg2_pi0", "mcMg2_pi0TE"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0]->Integral(25,150)/hmc[0]->Integral(25,150);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   // Data
   hdat[0]->SetAxisRange(0.005,0.03,"X");
   hdat[0]->SetTitle(";M^{2}_{#gamma#gamma} , GeV^{2}/c^{4}"
         ";Entries/0.0002 GeV^{2}/c^{4}"); // OK!
   if ( date == 2009 ) {
      hdat[0]->GetYaxis()->SetMaxDigits(2);
      hdat[0]->GetYaxis()->SetNdivisions(1005);
   } else {
      hdat[0]->GetYaxis()->SetMaxDigits(3);
   }
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.1 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(0.013,.0,0.022,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[1]->SetLineWidth(1);
   hmc[1]->SetLineColor(kBlue+1);
   hmc[1]->SetFillStyle(3001);
   hmc[1]->SetFillColor(kGreen+1);
   hmc[1]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC signal","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_pi0_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 ++ J/Psi -> gamma eta ++ Plot recoil mass pi+pi-gamma
//--------------------------------------------------------------------
void plot_Mpipig(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   vector<string> hndat {"Mrg2"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmcBG {"mcMrg2F"};
   vector<TH1D*> hmcBG;
   get_hst(par, 1, hnmcBG, hmcBG); // only bg !

   vector<string> hnmcS {"mcMrg2T", "mcMrg2TT"};
   vector<TH1D*> hmcS;
   get_hst(par, 2, hnmcS, hmcS); // only signal !

   // sum of signal and bg.
   vector<TH1D*> hmc(4);
   hmc[0] = (TH1D*)hmcS[0]->Clone("MC_total");
   hmc[1] = hmcS[0];
   hmc[2] = hmcBG[0];
   hmc[0]->Add( hmc[1], hmc[2], par->W_g_eta(), 1.);
   hmc[1]->Scale( par->W_g_eta() );
   hmc[3] = hmcS[1];
   hmc[3]->Scale( par->W_g_eta() );

   // rebining
   hdat[0]->Rebin();
   for ( auto& h : hmc ) {
      h->Rebin();
   }

   // normalization on DATA
   double scale = hdat[0]->Integral(30,70) / hmc[0]->Integral(30,70);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);
   const double meta   = 0.547862; // 547.862   +/- 0.017   MeV
   double M2l = SQ(meta)-0.2;
   double M2r = SQ(meta)+0.2;

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   // Data
   // hdat[0]->SetAxisRange(-0.1,0.7,"X");
   hdat[0]->SetTitle(
         ";M^{2}_{recoil}(#pi^{#plus}#pi^{#minus}#gamma),"
         " GeV^{2}/c^{4}"
         ";Entries/0.01 GeV^{2}/c^{4}");
   if ( date == 2009 ) {
      hdat[0]->GetYaxis()->SetNdivisions(1005);
   }
   hdat[0]->GetYaxis()->SetMaxDigits(2);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.07 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(-0.2,0.,M2l,ymax);
   box->DrawBox( M2r,0.,0.8,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[1]->SetLineWidth(1);
   hmc[1]->SetLineColor(kGreen+1);
   hmc[1]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   // hmc[2]->SetFillStyle(3001);
   // hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   hmc[3]->SetLineWidth(1);
   hmc[3]->SetLineStyle(kDashed);
   hmc[3]->SetLineColor(kGreen+1);
   // hmc[3]->Draw("HIST,SAME"); // true gamma from J/Psi decay!

   TLegend* leg = new TLegend(0.12,0.65,0.42,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[1], "MC signal","L");
   leg->AddEntry(hmc[2], "MC background","L");
   leg->AddEntry(hmc[0], "MC total","L");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_Mpipig_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot invariant mass of 2nd photon and «missing photon»
//--------------------------------------------------------------------
void plot_Minv2g(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   vector<string> hndat {"Mgg2_a"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmcBG {"mcMgg2_aF"};
   vector<TH1D*> hmcBG;
   get_hst(par, 1, hnmcBG, hmcBG); // only bg !

   vector<string> hnmcS {"mcMgg2_aT", "mcMgg2_aTT"}; // TT or TE
   vector<TH1D*> hmcS;
   get_hst(par, 2, hnmcS, hmcS); // only signal !

   // sum of signal and bg.
   vector<TH1D*> hmc(4);
   hmc[0] = (TH1D*)hmcS[0]->Clone("MC_total");
   hmc[1] = hmcS[0];
   hmc[2] = hmcBG[0];
   hmc[0]->Add( hmc[1], hmc[2], par->W_g_eta(), 1.);
   hmc[1]->Scale( par->W_g_eta() );
   hmc[3] = hmcS[1];
   hmc[3]->Scale( par->W_g_eta() );
   hmc[3]->Add( hmc[3], hmc[2] );

   double bg = hmc[2]->Integral(84,113); // [bin1,bin2]
   double sum = hmc[0]->Integral(84,113);
   printf(" %i %s background is %.2f%%\n",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral(84,113)/hmc[0]->Integral(84,113);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   // Data
   hdat[0]->SetAxisRange(0.198,0.4019,"X");
   hdat[0]->SetTitle(
         ";M^{2}_{inv}(#gamma_{2}#gamma_{missing} ), GeV^{2}/c^{4}"
         ";Entries/0.003 GeV^{2}/c^{4}");
   if ( date == 2009 ) {
      hdat[0]->GetYaxis()->SetNdivisions(1005);
   }
   hdat[0]->GetYaxis()->SetMaxDigits(2);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.07 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(0.198,0.,0.25, ymax);
   box->DrawBox(0.34, 0.,0.402,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   hmc[3]->SetLineWidth(1);
   hmc[3]->SetLineStyle(kDashed);
   hmc[3]->SetLineColor(kGreen+1);
   // hmc[3]->Draw("HIST,SAME"); // true gamma or true eta->2gamma

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_Minv2g_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot missing mass of pi+ pi- gamma gamma
//--------------------------------------------------------------------
void plot_M2mis(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   vector<string> hndat {"M2fr_min"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmcBG {"mcM2fr_minF"};
   vector<TH1D*> hmcBG;
   get_hst(par, 1, hnmcBG, hmcBG); // only bg !

   vector<string> hnmcS {"mcM2fr_minT", "mcM2fr_minTT"}; // TT or TE
   vector<TH1D*> hmcS;
   get_hst(par, 2, hnmcS, hmcS); // only signal !

   // sum of signal and bg.
   vector<TH1D*> hmc(4);
   hmc[0] = (TH1D*)hmcS[0]->Clone("MC_total");
   hmc[1] = hmcS[0];
   hmc[2] = hmcBG[0];
   hmc[0]->Add( hmc[1], hmc[2], par->W_g_eta(), 1.);
   hmc[1]->Scale( par->W_g_eta() );
   hmc[3] = hmcS[1];
   hmc[3]->Scale( par->W_g_eta() );

   double bg = hmc[2]->Integral(1,10); // [bin1,bin2]
   double sum = hmc[0]->Integral(1,10);
   printf(" %i %s background is %.2f%%\n",date,__func__,bg/sum*100);
   double bg_g = hmc[1]->Integral(1,10) - hmc[3]->Integral(1,10);
   printf(" %i %s bg(wrong gamma) is %.2f%%\n",
         date, __func__, bg_g/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral(1,10) / hmc[0]->Integral(1,10);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   hdat[0]->SetAxisRange(0.,0.0079,"X");
   hdat[0]->SetMinimum(0.5);
   hdat[0]->SetTitle(
         ";M^{2}_{missing}(#pi^{#plus}#pi^{#minus}#gamma#gamma),"
         " GeV^{2}/c^{4}"
         ";Entries/0.0002 GeV^{2}/c^{4}"); // OK!
   hdat[0]->GetXaxis()->SetNdivisions(1004);
   hdat[0]->GetYaxis()->SetMaxDigits(2);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(0.002,0.5,0.008, ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   hmc[3]->SetLineWidth(1);
   hmc[3]->SetLineStyle(kDashed);
   hmc[3]->SetLineColor(kGreen+1);
   // hmc[3]->Draw("HIST,SAME"); // true gamma or true eta->2gamma

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_M2miss_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot E of gamma from decay eta->2gamms
//--------------------------------------------------------------------
void plot_Eg(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   auto hst = [](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 34,0.,1.7);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut Cbg1; // see Cbg but without cut on Eg1
   Cbg1 += TCut("m2fr<0.002");
   Cbg1 += TCut("fabs(Cg0)<0.8");
   Cbg1 += TCut("fabs(Cg1)<0.8");

   TCut cutD = Cbg1 + par->Ceta;
   eff_etaD->Draw("Eg1>>hdat",cutD,"goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(4,nullptr);
   hmc[2] = hst("hmcB");
   eff_etaMC->Draw( "Eg1>>hmcB",cutD +!par->Cmcsig,"goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw( "Eg1>>hmcS",cutD + par->Cmcsig,"goff");
   hmc[3] = hst("hmcSTT");
   eff_etaMC2->Draw( "Eg1>>hmcSTT",cutD+par->Cmcsig+TCut("dgam!=1"),
         "goff");

   // sum of signal and bg.
   hmc[0] = hst("hmcT");
   hmc[0]->Add( hmc[1], hmc[2], par->W_g_eta(), 1.);
   hmc[1]->Scale( par->W_g_eta() );
   hmc[3]->Scale( par->W_g_eta() );

   double bg = hmc[2]->Integral(5,34); // [bin1,bin2]
   double sum = hmc[0]->Integral(5,34);
   printf(" %i %s background is %.2f%%\n",date,__func__,bg/sum*100);
   double bg_g = hmc[3]->Integral(5,34);
   printf(" %i %s bg(wrong gamma) is %.2f%%\n",
         date, __func__, bg_g/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral(5,34) / hmc[0]->Integral(5,34);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   double ymax = hmc[0]->GetMaximum();
   ymax = 1.08 * max(ymax,hdat[0]->GetMaximum());
   hdat[0]->SetMaximum(ymax);

   hdat[0]->SetTitle(
         ";E_{#gamma}, GeV"
         ";Entries/0.025 GeV"); // OK!
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   // hdat[0]->GetYaxis()->SetNdivisions(504);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   box->DrawBox(0.,0.,0.2, ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   hmc[3]->SetLineWidth(2);
   hmc[3]->SetLineStyle(kDashed);
   hmc[3]->SetLineColor(kBlue+3);
   // hmc[3]->Draw("HIST,SAME"); // wrong gamma

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_Eg_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot: match Energy and Theta
//--------------------------------------------------------------------
void plot_rE(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   vector<string> hndat {"rE"}; // bug in binning in prod-12!!!
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"rE", "mcrET", "mcrEF"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetLogy(true);

   // Data
   hdat[0]->SetMinimum(1.);
   hdat[0]->SetTitle(
         ";E_{#gamma}(pred) / E_{#gamma}(rec)"
         ";Entries / 0.01");
   // hdat[0]->GetXaxis()->SetNdivisions(1004);
   hdat[0]->SetAxisRange(0.0,2.0,"X");
   hdat[0]->GetYaxis()->SetMaxDigits(2);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(0., 1.,0.4, ymax);
   box->DrawBox(1.8,1.,2.01, ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_rE_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}


//--------------------------------------------------------------------
void plot_dTh(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   vector<string> hndat {"dTh"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"dTh", "mcdThT", "mcdThF"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetLogy(true);

   // Data
   hdat[0]->SetMinimum(1.);
   hdat[0]->SetTitle(
         ";#delta#Theta(#gamma), deg."
         ";Entries / 0.2 deg.");
   // hdat[0]->GetXaxis()->SetNdivisions(1004);
   hdat[0]->GetYaxis()->SetMaxDigits(2);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(10., 1.,20., ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_dTh_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot of Mgg
//--------------------------------------------------------------------
void plot_Mgg(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   auto hst = [](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 80,0.51,0.59);
      h->Sumw2(true);
      return h;
   };

   const double meta   = 0.547862; // 547.862   +/- 0.017   MeV
   const double weta   = 3*0.008; // see cuts.h

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cbg + par->Ceta + TCut("fl>0.5");
   eff_etaD->Draw("mggf>>hdat",cutD,"goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(3,nullptr);
   hmc[2] = hst("hmcB");
   eff_etaMC->Draw( "mggf>>hmcB",cutD +!par->Cmcsig,"goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw( "mggf>>hmcS",cutD + par->Cmcsig,"goff");

   // sum of signal and bg.
   hmc[0] = hst("hmcT");
   hmc[0]->Add( hmc[1], hmc[2], par->W_g_eta(), 1.);
   hmc[1]->Scale( par->W_g_eta() );

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // fit function
   TF1* gs = (TF1*)gROOT->GetFunction("gaus");
   gs->SetLineWidth(2);
   gs->SetLineColor(kGreen+2);

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.19);
   gStyle->SetStatH(0.12);
   gStyle->SetFitFormat(".3g");

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetTitle(
         ";M_{inv}(#gamma#gamma), GeV/c^{2}"
         ";Entries/0.001 GeV/c^{2}");
   hdat[0]->GetYaxis()->SetNdivisions(1005);
   hdat[0]->GetYaxis()->SetMaxDigits(2);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->SetMarkerStyle(20);
   double ymax = 1.12 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   double Mmin = meta-weta, Mmax = meta+weta;
   hdat[0]->Fit(gs,"QM","",Mmin,Mmax);
   hdat[0]->Draw("E");

   box->DrawBox( 0.51,      0., meta-weta, ymax);
   box->DrawBox( meta+weta, 0., 0.59,      ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(2);
   hmc[0]->SetLineColor(kRed+1);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.11,0.64,0.36,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(gs, "Fit of data","L");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_Mgg_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot predicted P(eta) and cos(Theta(eta))
//--------------------------------------------------------------------
void plot_Peta(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   auto hst = [](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 55,1.25,1.8);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut Ceta;   // cuts for selection eta ++ without Peta cut ++
   Ceta += TCut("fabs(Ceta)<0.9");
   Ceta += TCut("fabs(Cg2)<0.8||(fabs(Cg2)>0.85&&fabs(Cg2)<0.92)");

   TCut cutD = par->Cbg + Ceta;
   eff_etaD->Draw("Peta>>hdat",cutD,"goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(4,nullptr);
   hmc[2] = hst("hmcB");
   eff_etaMC->Draw( "Peta>>hmcB",cutD +!par->Cmcsig,"goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw( "Peta>>hmcS",cutD + par->Cmcsig,"goff");
   hmc[3] = hst("hmcSF");
   eff_etaMC2->Draw( "Peta>>hmcSF",cutD+par->Cmcsig+TCut("dgam!=1"),
         "goff");

   // sum of signal and bg.
   hmc[0] = hst("hmcT");
   hmc[0]->Add( hmc[1], hmc[2], par->W_g_eta(), 1.);
   hmc[1]->Scale( par->W_g_eta() );
   hmc[3]->Scale( par->W_g_eta() );

   double bg = hmc[2]->Integral(6,45); // [bin1,bin2]
   double sum = hmc[0]->Integral(6,45);
   printf(" %i %s background is %.2f%%\n",date,__func__,bg/sum*100);
   double bg_g = hmc[3]->Integral(6,45);
   printf(" %i %s bg(wrong gamma) is %.2f%%\n",
         date, __func__, bg_g/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetTitle(
         ";P_{#eta}, GeV/c"
         ";Entries/0.01 GeV/c");
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   hmc[3]->SetLineWidth(2);
   hmc[3]->SetLineStyle(kDashed);
   hmc[3]->SetLineColor(kBlue+3);
   // hmc[3]->Draw("HIST,SAME"); // wrong gamma

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_Peta_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot_Ceta(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   int Nbins = 40;
   auto hst = [Nbins](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", Nbins,-1.,1.);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut Ceta;   // cuts for selection eta ++ without Ceta cuts ++
   Ceta += TCut("Peta>1.3&&Peta<1.7");
   Ceta += TCut("fabs(Cg2)<0.8||(fabs(Cg2)>0.85&&fabs(Cg2)<0.92)");

   TCut cutD = par->Cbg + Ceta;
   eff_etaD->Draw("Ceta>>hdat",cutD,"goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(4,nullptr);
   hmc[2] = hst("hmcB");
   eff_etaMC->Draw( "Ceta>>hmcB",cutD +!par->Cmcsig,"goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw( "Ceta>>hmcS",cutD + par->Cmcsig,"goff");
   hmc[3] = hst("hmcSTT");
   eff_etaMC2->Draw( "Ceta>>hmcSTT",cutD+par->Cmcsig+TCut("dgam!=1"),
         "goff");

   // sum of signal and bg.
   hmc[0] = hst("hmcT");
   hmc[0]->Add( hmc[1], hmc[2], par->W_g_eta(), 1.);
   hmc[1]->Scale( par->W_g_eta() );
   hmc[3]->Scale( par->W_g_eta() );

   double bg = hmc[2]->Integral(3,38); // [bin1,bin2]
   double sum = hmc[0]->Integral(3,38);
   printf(" %i %s background is %.2f%%\n",date,__func__,bg/sum*100);
   double bg_g = hmc[3]->Integral(3,38);
   printf(" %i %s bg(wrong gamma) is %.2f%%\n",
         date, __func__, bg_g/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);

   hdat[0]->SetTitle( Form(
            ";cos(#Theta_{#eta}) "
            ";Entries/%g",2./Nbins));
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   // hdat[0]->GetXaxis()->SetNdivisions(504);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   hmc[3]->SetLineWidth(2);
   hmc[3]->SetLineStyle(kDashed);
   hmc[3]->SetLineColor(kBlue+3);
   // hmc[3]->Draw("HIST,SAME"); // wrong gamma

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_Ceta_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 ++ J/Psi -> phi eta ++ Plot invariant mass of K+K-
//--------------------------------------------------------------------
void plot2_MKK(int date) {
//--------------------------------------------------------------------
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   vector<string> hndat {"E_mkk2"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmcBG {"Emc_mkk2F"};
   vector<TH1D*> hmcBG;
   get_hst(par, 1, hnmcBG, hmcBG); // bg only

   vector<string> hnmcS {"Emc_mkk2T"};
   vector<TH1D*> hmcS;
   get_hst(par, 2, hnmcS, hmcS); // signal

   // sum of signal and bg.
   vector<TH1D*> hmc(3);
   hmc[0] = (TH1D*)hmcS[0]->Clone("MC_total");
   hmc[1] = hmcS[0];
   hmc[1]->Scale( par->W_phi_eta() );
   hmc[2] = hmcBG[0];
   hmc[0]->Add( hmc[1], hmc[2] );

   // normalization on DATA
   // double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   double scale = hdat[0]->Integral(49,69) / hmc[0]->Integral(49,69);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   // Data
   hdat[0]->SetAxisRange(0.98,1.0999,"X");
   hdat[0]->SetTitle(
         ";M^{2}_{inv}(K^{#plus}K^{#minus}), GeV^{2}/c^{4}"
         ";Entries/0.001 GeV^{2}/c^{4}");
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.1 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(0.98,0.,1.01,ymax);
   box->DrawBox(1.07,0.,1.1,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   // hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaphi_Mkk2_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot invariant mass of photon and «missing photon»
//--------------------------------------------------------------------
void plot2_Minv2g(int date) {
//--------------------------------------------------------------------
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   vector<string> hndat {"E_M2gg"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmcBG {"Emc_M2ggF"};
   vector<TH1D*> hmcBG;
   get_hst(par, 1, hnmcBG, hmcBG); // only bg !

   vector<string> hnmcS {"Emc_M2ggT", "Emc_M2ggTE"};
   vector<TH1D*> hmcS;
   get_hst(par, 2, hnmcS, hmcS); // only signal !

   // sum of signal and bg.
   vector<TH1D*> hmc(4);
   hmc[0] = (TH1D*)hmcS[0]->Clone("MC_total");
   hmc[1] = hmcS[0];
   hmc[1]->Scale( par->W_phi_eta() );
   hmc[2] = hmcBG[0];
   hmc[0]->Add( hmc[1], hmc[2] );
   hmc[3] = hmcS[1];
   hmc[3]->Scale( par->W_phi_eta() );

   // normalization on DATA
   double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   // double scale = hdat[0]->Integral(84,113)/hmc[0]->Integral(84,113);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   // Data
   hdat[0]->SetAxisRange(0.198,0.4019,"X");
   hdat[0]->SetTitle(
         ";M^{2}_{inv}(#gamma#gamma_{missing} ), GeV^{2}/c^{4}"
         ";Entries/0.003 GeV^{2}/c^{4}");
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.1 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(0.198,0.,0.25, ymax);
   box->DrawBox(0.34, 0.,0.402,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   hmc[3]->SetLineWidth(1);
   hmc[3]->SetLineStyle(kDashed);
   hmc[3]->SetLineColor(kGreen+1);
   // hmc[3]->Draw("HIST,SAME"); // true eta->2gamma decay

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaphi_Minv2g_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot missing mass of pi+ pi- gamma gamma
//--------------------------------------------------------------------
void plot2_M2mis(int date) {
//--------------------------------------------------------------------
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   vector<string> hndat {"E_M2fr"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmcBG {"Emc_M2frF"};
   vector<TH1D*> hmcBG;
   get_hst(par, 1, hnmcBG, hmcBG); // only bg !

   vector<string> hnmcS {"Emc_M2frT","Emc_M2frTE"};
   vector<TH1D*> hmcS;
   get_hst(par, 2, hnmcS, hmcS); // only signal !

   // sum of signal and bg.
   vector<TH1D*> hmc(4);
   hmc[0] = (TH1D*)hmcS[0]->Clone("MC_total");
   hmc[1] = hmcS[0];
   hmc[1]->Scale( par->W_phi_eta() );
   hmc[2] = hmcBG[0];
   hmc[0]->Add( hmc[1], hmc[2] );
   hmc[3] = hmcS[1];
   hmc[3]->Scale( par->W_phi_eta() );

   double bg = hmc[2]->Integral(1,10); // [bin1,bin2]
   double sum = hmc[0]->Integral(1,10);
   printf(" %i %s background is %.2f%%\n",date,__func__, bg/sum*100);
   double bg_eta = hmc[1]->Integral(1,10) - hmc[3]->Integral(1,10);
   printf(" %i %s bg(eta not 2g) is %.2f%%\n",
         date, __func__, bg_eta/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral(1,10) / hmc[0]->Integral(1,10);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   hdat[0]->SetAxisRange(0.,0.0079,"X");
   hdat[0]->SetMinimum(0.5);
   hdat[0]->SetTitle( ";M^{2}_{missing}"
         "(#pi^{#plus}#pi^{#minus}K^{#plus}K^{#minus}#gamma),"
         " GeV^{2}/c^{4}"
         ";Entries/0.0002 GeV^{2}/c^{4}");
   hdat[0]->GetXaxis()->SetNdivisions(1005);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(0.002,0.5,0.008, ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   hmc[3]->SetLineWidth(1);
   hmc[3]->SetLineStyle(kDashed);
   hmc[3]->SetLineColor(kGreen+1);
   // hmc[3]->Draw("HIST,SAME"); // true eta->2gamma decay

   TLegend* leg = new TLegend(0.59,0.67,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaphi_M2miss_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot: match Energy and Theta
//--------------------------------------------------------------------
void plot2_rE(int date) {
//--------------------------------------------------------------------
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   vector<string> hndat {"E_rE"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"E_rE", "Emc_rET", "Emc_rEF"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetLogy(true);

   // Data
   hdat[0]->SetMinimum(0.5);
   hdat[0]->SetTitle(
         ";E_{#gamma}(pred) / E_{#gamma}(rec)"
         ";Entries/0.01");
   hdat[0]->SetAxisRange(0.0,2.0,"X");
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(0., 0.5,0.4, ymax);
   box->DrawBox(1.8,0.5,2.01, ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.59,0.67,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaphi_rE_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot2_dTh(int date) {
//--------------------------------------------------------------------
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   vector<string> hndat {"E_dTh"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"E_dTh", "Emc_dThT", "Emc_dThF"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetLogy(true);

   // Data
   hdat[0]->SetMinimum(0.5);
   hdat[0]->SetTitle(
         ";#delta#Theta(#gamma), deg."
         ";Entries/0.2 deg.");
   // hdat[0]->GetXaxis()->SetNdivisions(1004);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(10.,0.5,20.,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.59,0.67,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaphi_dTh_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot predicted P(eta) and cos(Theta(eta))
//--------------------------------------------------------------------
void plot2_Peta(int date) {
//--------------------------------------------------------------------
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   auto hst = [](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 45,1.1,1.55);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   // cuts for selection eta ++ without Peta & Ceta cuts ++
   TCut Ceta;
   // Ceta += TCut("fabs(Ceta)<0.9");
   Ceta += TCut("fabs(Cg2)<0.8||(fabs(Cg2)>0.85&&fabs(Cg2)<0.92)");

   TCut cutD = par->Cbg + Ceta;
   eff_etaD->Draw("Peta>>hdat",cutD,"goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(4,nullptr);
   hmc[2] = hst("hmcB");
   eff_etaMC->Draw( "Peta>>hmcB",cutD +!par->Cmcsig,"goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw( "Peta>>hmcS",cutD + par->Cmcsig,"goff");
   hmc[3] = hst("hmcSF");
   eff_etaMC2->Draw( "Peta>>hmcSF",cutD+par->Cmcsig+TCut("deta!=1"),
         "goff");

   // sum of signal and bg.
   hmc[0] = hst("hmcT");
   hmc[1]->Scale( par->W_phi_eta() );
   hmc[0]->Add( hmc[1], hmc[2] );
   hmc[3]->Scale( par->W_phi_eta() );

   double bg = hmc[2]->Integral(); // [bin1,bin2]
   double sum = hmc[0]->Integral();
   printf(" %i %s background is %.2f%%\n",date,__func__, bg/sum*100);
   double bg_eta = hmc[3]->Integral();
   printf(" %i %s bg(eta not 2g) is %.2f%%\n",
         date, __func__, bg_eta/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetTitle(
         ";P_{#eta}, GeV/c"
         ";Entries/0.01 GeV/c");
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   hmc[0]->SetLineWidth(2);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   hmc[3]->SetLineWidth(2);
   hmc[3]->SetLineStyle(kDashed);
   hmc[3]->SetLineColor(kBlue+3);
   // hmc[3]->Draw("HIST,SAME"); // wrong eta!

   TLegend* leg = new TLegend(0.59,0.71,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaphi_Peta_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot2_Ceta(int date) {
//--------------------------------------------------------------------
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   int Nbins = 40;
   auto hst = [Nbins](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", Nbins,-1.,1.);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   // cuts for selection eta ++ without Peta & Ceta cuts ++
   TCut Ceta;
   // Ceta += TCut("fabs(Ceta)<0.9");
   Ceta += TCut("fabs(Cg2)<0.8||(fabs(Cg2)>0.85&&fabs(Cg2)<0.92)");

   TCut cutD = par->Cbg + Ceta;
   eff_etaD->Draw("Ceta>>hdat",cutD,"goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(4,nullptr);
   hmc[2] = hst("hmcB");
   eff_etaMC->Draw( "Ceta>>hmcB",cutD +!par->Cmcsig,"goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw( "Ceta>>hmcS",cutD + par->Cmcsig,"goff");
   hmc[3] = hst("hmcSF");
   eff_etaMC2->Draw( "Ceta>>hmcSF",cutD+par->Cmcsig+TCut("deta!=1"),
         "goff");

   // sum of signal and bg.
   hmc[0] = hst("hmcT");
   hmc[1]->Scale( par->W_phi_eta() );
   hmc[0]->Add( hmc[1], hmc[2] );
   hmc[3]->Scale( par->W_phi_eta() );

   double bg = hmc[2]->Integral(); // [bin1,bin2]
   double sum = hmc[0]->Integral();
   printf(" %i %s background is %.2f%%\n",date,__func__, bg/sum*100);
   double bg_eta = hmc[3]->Integral();
   printf(" %i %s bg(eta not 2g) is %.2f%%\n",
         date,__func__,bg_eta/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);

   hdat[0]->SetTitle( Form (
            ";cos(#Theta_{#eta}) "
            ";Entries/%1g",2./Nbins));
   hdat[0]->GetXaxis()->SetNdivisions(504);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   hmc[3]->SetLineWidth(2);
   hmc[3]->SetLineStyle(kDashed);
   hmc[3]->SetLineColor(kBlue+3);
   // hmc[3]->Draw("HIST,SAME"); // wrong eta!

   TLegend* leg = new TLegend(0.59,0.71,0.89,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaphi_Ceta_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Main
//--------------------------------------------------------------------
void eta_eff_sel() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(42);

   for ( auto date : {2009, 2012, 2021} ) {
   // for ( auto date : {2009} ) {
      // ++ pi0 rejection, fig.50
      // plot_pi0(date);

      // J/Psi -> gamma eta
      // +fig 51
      // plot_Mpipig(date);
      // +fig 52
      // plot_Minv2g(date);
      // +fig 53
      // plot_M2mis(date);
      // +fig 54
      // plot_Eg(date);
      // +fig 55
      // plot_rE(date);
      // plot_dTh(date);
      // +fig 56
      // plot_Mgg(date);
      // +fig 57
      // plot_Peta(date);
      // plot_Ceta(date);

      // J/Psi -> phi eta
      // +fig 60
      // plot2_MKK(date);
      // +fig 61
      // plot2_Minv2g(date);
      // +fig 62
      // plot2_M2mis(date);
      // +fig 63
      // plot2_rE(date);
      // plot2_dTh(date);
      // +fig 64
      // plot2_Peta(date);
      // plot2_Ceta(date);
   }

}
