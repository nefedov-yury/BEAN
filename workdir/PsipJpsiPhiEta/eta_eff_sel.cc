// eta_eff_sel.cc
// Study of the eta->2gamma reconstruction efficiency
// Pictures for event selection.

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
   // const string Dir = "prod_v709n3/";
   const string Dir = "prod_v709n4/";
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
   TCut Cbg0;   // Cbg without cut on Eg1
   TCut Cph;    // cuts for selection photons
   TCut Ceta;   // cuts for selection eta
   TCut Cfnd;   // predicted eta found

   const double meta   = 0.54786; // 547.862   +/- 0.017   MeV
   const double Br_eta_2gamma = 0.3936; // PDG
   const double Br_phi_KK = 0.491; // PDG
};

// {{{2 > ctor
//--------------------------------------------------------------------
Params::Params(int dat, int slc = 2, int rew = 0)
//--------------------------------------------------------------------
{
   date = dat;
   slct = slc;
   use_rew = rew;

   // set the names:
   datafile  = string( Form("data_%02ipsip_all.root",date%100) );
   mcincfile = string( Form("mcinc_%02ipsip_all.root",date%100) );
   mcsigf1 = string( Form("mcgammaeta2_kkmc_%02i.root",date%100) );
   mcsigf2 = string( Form("mcphieta2_kkmc_%02i.root",date%100) );

   // mc-signal
   if ( slct > 0 ) {    // gamma-eta
      // Cmcsig += TCut("decj==22");
      Cmcsig += TCut("abs(decj-22)<0.1"); // decj is float here
   } else {             // phi-eta
      Cmcsig += TCut("decj==68");
   }

   // cuts against the background
   if ( slct > 0 ) {         // gamma-eta
      Cbg += TCut("abs(Cg0)<0.8");
      Cbg += TCut("abs(Cg1)<0.8");
      Cbg += TCut("m2fr<0.002");
      Cbg0 = Cbg; // for study of the cut on Eg1
      Cbg += TCut("Eg1>0.25&&Eg1<1.7");
   } else {                  // phi-eta
      Cbg += TCut("m2fr<0.001");
   }

   // cuts for selection photons
   Cph += TCut("abs(Cg1)<0.8||(abs(Cg1)>0.85&&abs(Cg1)<0.92)");
   Cph += TCut("abs(Cg2)<0.8||(abs(Cg2)>0.85&&abs(Cg2)<0.92)");
   if ( slct > 0 ) {         // gamma-eta
      Cph += TCut("Eg2>0.1&&Eg2<1.4");
   } else {                  // phi-eta
      Cph += TCut("Eg2>0.05&&Eg2<1.45");
   }

   // cuts for selection eta
   if ( slct > 0 ) {         // gamma-eta
      Ceta += TCut("Peta>1.3&&Peta<1.7");
   } else {                  // phi-eta
      Ceta += TCut("Peta>1.15&&Peta<1.5");
   }
   // doubtful:
   // Ceta += TCut("abs(Ceta)<0.9");
   // Ceta += TCut("abs(Cg2)<0.8||(abs(Cg2)>0.85&&abs(Cg2)<0.92)");

   // predicted eta is found
   Cfnd += TCut("dTh<8.&&0.8<rE&&rE<1.4");
   Cfnd += TCut( Form("abs(mggf-%.6f)<0.024",meta) );
}

// {{{2 > OpenFile(mc = 0: data, 1: MC inc, 2: MC-signal)
//--------------------------------------------------------------------
TFile* Params::OpenFile(int mc)
//--------------------------------------------------------------------
{
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
TTree* Params::GetEffEta(int mc = 0)
//--------------------------------------------------------------------
{
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
double Params::W_g_eta() // wights for MC-gamma-eta
//--------------------------------------------------------------------
{
   // normalization on numbers in "official inclusive MC"
   // decay code for J/Psi -> gamma eta is 22;
   // BOSS-664: ((date==2012) ? (137258./5e5) : (35578./1.5e5))
   double W = 1;
   switch (date) {
      case 2009:
         W =  41553./150e3; // 0.28
         break;
      case 2012:
         W = 130865./500e3; // 0.26
         break;
      case 2021:
         W = 883676./2.5e6; // 0.35
         break;
   }
   return W;
}

//--------------------------------------------------------------------
double Params::W_phi_eta() // wights for MC-sig
//--------------------------------------------------------------------
{
   // normalization on numbers in "official inclusive MC"
   // decay code for J/Psi -> phi eta is 68;
   // BOSS-664: ((date==2012) ? (104950./5e5) : (27274./1.5e5));
   double W = 1;
   switch (date) {
      case 2009:
         W =  27547./150e3; // 0.18
         break;
      case 2012:
         W =  89059./500e3; // 0.18
         break;
      case 2021:
         W = 598360./2.5e6; // 0.24
         break;
   }
   return W * Br_phi_KK;
}

// {{{1 helper functions and constants
//--------------------------------------------------------------------
constexpr double SQ(double x)
//--------------------------------------------------------------------
{
   return x*x;
}

//--------------------------------------------------------------------
void SetHstFace(TH1* hst)
//--------------------------------------------------------------------
{
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
      const vector<string>& hns, vector<TH1D*>& hst )
//--------------------------------------------------------------------
{
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
void plot_pi0(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
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
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
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
   hdat[0]->GetYaxis()->SetTitleOffset(1.15);
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

   TLegend* leg = new TLegend(0.64,0.65,0.892,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   // leg->AddEntry(hmc[1], "MC signal","F");
   leg->AddEntry(hmc[1], "MC signal #gamma#eta","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_pi0_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 ++ J/Psi -> gamma eta ++ Plot recoil mass pi+pi-gamma
//--------------------------------------------------------------------
void plot_Mpipig(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   // Here we use histograms filled in the PsipJpsiGammaEta.
   // It is incorrect to use ntuple here.

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
   double M2l = SQ(par->meta)-0.2;
   double M2r = SQ(par->meta)+0.2;

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
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
   hdat[0]->GetXaxis()->SetTitleOffset(1.05);
   hdat[0]->GetYaxis()->SetTitleOffset(1.15);
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

   TLegend* leg = new TLegend(0.12,0.64,0.40,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   // leg->AddEntry(hmc[1], "MC signal","L");
   leg->AddEntry(hmc[1], "MC signal #gamma#eta","L");
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
void plot_Minv2g(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   // Here we use histograms filled in the PsipJpsiGammaEta.
   // It is incorrect to use ntuple here.

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

   // recalculation: Min,Max -> bin1,bin2
   double Minv2_min = 0.25, Minv2_max = 0.34;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Minv2_min+bW/2);
   int bin2 = hdat[0]->FindBin(Minv2_max-bW/2);
   // cout << "bW=" << bW
      // << " bin1=" << bin1 << " bin2=" << bin2 << endl;

   // normalization on DATA
   double scale = hdat[0]->Integral(bin1,bin2)/
      hmc[0]->Integral(bin1,bin2);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral(bin1,bin2);
   double sum = hmc[0]->Integral(bin1,bin2);
   printf("%i %s bkg: %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100);

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   // Data
   double Xmin = 0.198, Xmax = 0.4019;
   hdat[0]->SetAxisRange(Xmin,Xmax,"X");
   hdat[0]->SetTitle(
         ";M^{2}_{inv}(#gamma_{2}#gamma_{missing} ), GeV^{2}/c^{4}"
         ";Entries/0.003 GeV^{2}/c^{4}");
   if ( date == 2009 ) {
      hdat[0]->GetYaxis()->SetNdivisions(1005);
   }
   hdat[0]->GetYaxis()->SetMaxDigits(2);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->GetYaxis()->SetTitleOffset(1.15);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.07 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(Xmin,0.,      Minv2_min,ymax);
   box->DrawBox(Minv2_max,0., Xmax,ymax);
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

   TLegend* leg = new TLegend(0.61,0.65,0.892,0.89);
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

// {{{1 New Plot invariant mass of 2nd photon and «missing photon»
//--------------------------------------------------------------------
void plot_Minv2g_n(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   // Filling histograms from ntuple
   double Xmin = 0.2, Xmax = 0.4;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 100,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };
   TCut cutD = TCut("abs(Cg0)<0.8&&abs(Cg1)<0.8") +
      TCut("Peta>1.3&&Peta<1.7");
   eff_etaD->Draw("M2gg>>hdat", cutD, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(4,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("M2gg>>hmcB", cutD+!(par->Cmcsig), "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   // add to background g-eta but !(eta->2gamma)
   eff_etaMC2->Draw("M2gg>>hmcT", cutD+TCut("deta!=1"), "goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_g_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("M2gg>>hmcS", cutD+TCut("deta==1"), "goff");
   hmc[1]->Scale( par->W_g_eta() );

   hmc[3] = hst("hmcSTT"); // to estimate incorect gamma matching
   eff_etaMC2->Draw("M2gg>>hmcSTT", cutD+TCut("dgam!=1"), "goff");
   hmc[3]->Scale( par->W_g_eta() );

   // sum of signal and bg.
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Minv2_min = 0.25, Minv2_max = 0.34;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Minv2_min+bW/2);
   int bin2 = hdat[0]->FindBin(Minv2_max-bW/2);
   // cout << "bW=" << bW
      // << " bin1=" << bin1 << " bin2=" << bin2 << endl;

   // normalization on DATA
   double scale = hdat[0]->Integral(bin1,bin2)/
      hmc[0]->Integral(bin1,bin2);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral(bin1,bin2);
   double sum = hmc[0]->Integral(bin1,bin2);
   printf("%i %s bkg: %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100);

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   // Data
   hdat[0]->SetTitle(
         ";M^{2}_{inv}(#gamma_{2}#gamma_{missing} ), GeV^{2}/c^{4}"
         ";Entries/0.003 GeV^{2}/c^{4}");
   if ( date == 2009 ) {
      hdat[0]->GetYaxis()->SetNdivisions(1005);
   }
   hdat[0]->GetYaxis()->SetMaxDigits(2);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->GetYaxis()->SetTitleOffset(1.15);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.07 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(Xmin,0.,      Minv2_min,ymax);
   box->DrawBox(Minv2_max,0., Xmax,ymax);
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

   TLegend* leg = new TLegend(0.61,0.65,0.892,0.89);
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
void plot_M2mis(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   // Filling histograms from ntuple
   double Xmin = 0., Xmax = 0.008;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 80,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };
   TCut cutD = TCut("abs(Cg0)<0.8&&abs(Cg1)<0.8") + par->Ceta;
   eff_etaD->Draw("m2fr>>hdat", cutD, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(4,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("m2fr>>hmcB", cutD+!(par->Cmcsig), "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   // add to background g-eta but !(eta->2gamma)
   eff_etaMC2->Draw("m2fr>>hmcT", cutD+TCut("deta!=1"), "goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_g_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("m2fr>>hmcS", cutD+TCut("deta==1"), "goff");
   hmc[1]->Scale( par->W_g_eta() );

   hmc[3] = hst("hmcSTT"); // to estimate incorect gamma matching
   eff_etaMC2->Draw("m2fr>>hmcSTT", cutD+TCut("dgam!=1"), "goff");
   hmc[3]->Scale( par->W_g_eta() );

   // sum of signal and bg.
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Mmis2_min = 0., Mmis2_max = 0.002;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = 1;
   int bin2 = hdat[0]->FindBin(Mmis2_max-bW/2);
   // cout << "bW=" << bW
      // << " bin1=" << bin1 << " bin2=" << bin2 << endl;

   // normalization on DATA
   double scale = hdat[0]->Integral(bin1,bin2) /
      hmc[0]->Integral(bin1,bin2);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral(bin1,bin2);
   double sum = hmc[0]->Integral(bin1,bin2);
   double bg_g = hmc[3]->Integral(bin1,bin2);
   printf("%i %s, bkg: %.2f+/-%.2f%%, wrong gamma: %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100,
         bg_g/sum*100, sqrt(bg_g)/sum*100);

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   SetHstFace(hdat[0]);
   double Ymin = 0.5;
   hdat[0]->SetMinimum(Ymin);
   hdat[0]->SetTitle(
         ";M^{2}_{missing}(#pi^{#plus}#pi^{#minus}#gamma#gamma_{2}),"
         " GeV^{2}/c^{4}"
         ";Entries/0.0002 GeV^{2}/c^{4}"); // OK!
   hdat[0]->GetXaxis()->SetNdivisions(1004);
   hdat[0]->GetYaxis()->SetMaxDigits(2);
   hdat[0]->GetYaxis()->SetTitleOffset(1.15);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(Mmis2_max,Ymin,Xmax,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   // incorect gamma matching
   // hmc[3]->SetLineWidth(1);
   // hmc[3]->SetLineStyle(kDashed);
   // hmc[3]->SetLineColor(kGreen+1);
   // hmc[3]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.61,0.65,0.892,0.89);
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
void plot_Eg(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   // Filling histograms from ntuple
   auto hst = [](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 34,0.,1.7);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cbg0 + TCut("Peta>1.3&&Peta<1.7");
   eff_etaD->Draw("Eg1>>hdat", cutD, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(4,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("Eg1>>hmcB", cutD+!(par->Cmcsig), "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   // add to background g-eta but !(eta->2gamma)
   eff_etaMC2->Draw("Eg1>>hmcT", cutD+TCut("deta!=1"), "goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_g_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("Eg1>>hmcS", cutD+TCut("deta==1"), "goff");
   hmc[1]->Scale( par->W_g_eta() );

   hmc[3] = hst("hmcSTT"); // to estimate incorect gamma matching
   eff_etaMC2->Draw("Eg1>>hmcSTT", cutD+TCut("dgam!=1"), "goff");
   hmc[3]->Scale( par->W_g_eta() );

   // sum of signal and bg.
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Eg1_min = 0.25, Eg1_max = 1.7;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Eg1_min+bW/2);
   int bin2 = hdat[0]->FindBin(Eg1_max-bW/2);
   // cout << "bW=" << bW
      // << " bin1=" << bin1 << " bin2=" << bin2 << endl;

   // normalization on DATA
   double scale = hdat[0]->Integral(bin1,bin2) /
      hmc[0]->Integral(bin1,bin2);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral(bin1,bin2);
   double sum = hmc[0]->Integral(bin1,bin2);
   double bg_g = hmc[3]->Integral(bin1,bin2);
   printf("%i %s bkg: %.2f+/-%.2f%%, wrong gamma: %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100,
         bg_g/sum*100, sqrt(bg_g)/sum*100);

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   double ymax = hmc[0]->GetMaximum();
   ymax = 1.08 * max(ymax,hdat[0]->GetMaximum());
   hdat[0]->SetMaximum(ymax);

   hdat[0]->SetTitle(
         ";E(#gamma_{2}), GeV"
         ";Entries/0.025 GeV");
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   // hdat[0]->GetYaxis()->SetNdivisions(504);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   box->DrawBox(0.,0.,Eg1_min, ymax);
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

   TLegend* leg = new TLegend(0.61,0.65,0.892,0.89);
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
void plot_rE(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   // Filling histograms from ntuple
   double Xmin = 0.4, Xmax = 1.8;
   double Ymin = 0.5; // log-Y
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 70,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cbg + par->Ceta;
   TCut Covr( Form("rE<%f||%f<=rE",Xmin,Xmax) ); // under & overflows
   eff_etaD->Draw("rE>>hdat", cutD, "goff");
   eff_etaD->Draw("1.79>>+hdat", cutD+Covr, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("rE>>hmcB", cutD+!(par->Cmcsig), "goff");
   eff_etaMC->Draw("1.79>>+hmcB", cutD+!(par->Cmcsig)+Covr, "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   // add to background g-eta but !(eta->2gamma)
   eff_etaMC2->Draw("rE>>hmcT", cutD+TCut("deta!=1"), "goff");
   eff_etaMC2->Draw("1.79>>+hmcT",cutD+TCut("deta!=1")+Covr,"goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_g_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("rE>>hmcS", cutD+TCut("deta==1"), "goff");
   eff_etaMC2->Draw("1.79>>+hmcS",cutD+TCut("deta==1")+Covr,"goff");
   hmc[1]->Scale( par->W_g_eta() );

   // sum of signal and bg.
   hmc[0]->Add( hmc[1], hmc[2] );

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
   auto name = Form("c1_rE_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetLogy(true);

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetMinimum(Ymin);
   hdat[0]->SetTitle(
         ";E_{#gamma}#kern[0.1]{(pred)} / "
         "E_{#gamma}#kern[0.1]{(rec)}"
         ";Entries / 0.01");
   // hdat[0]->GetXaxis()->SetNdivisions(1004);
   hdat[0]->GetYaxis()->SetMaxDigits(2);
   hdat[0]->GetYaxis()->SetTitleOffset(1.15);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(Xmin,Ymin, 0.8, ymax);
   box->DrawBox(1.4, Ymin, Xmax,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   // TLegend* leg = new TLegend(0.61,0.65,0.892,0.89);
   TLegend* leg = new TLegend(0.60,0.65,0.882,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   // leg->AddEntry(box, "Rejection area","F");
   leg->AddEntry(box, "Mismatch area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_rE_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot_dTh(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   // Filling histograms from ntuple
   double Xmin = 0., Xmax = 15.;
   double Ymin = 1.; // log-Y
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 75,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cbg + par->Ceta;
   TCut Covr( Form("%f<=dTh",Xmax) ); // under & overflows
   eff_etaD->Draw("dTh>>hdat", cutD, "goff");
   eff_etaD->Draw("14.9>>+hdat", cutD+Covr, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("dTh>>hmcB", cutD+!(par->Cmcsig), "goff");
   eff_etaMC->Draw("14.9>>+hmcB", cutD+!(par->Cmcsig)+Covr, "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   // add to background g-eta but !(eta->2gamma)
   eff_etaMC2->Draw("dTh>>hmcT", cutD+TCut("deta!=1"), "goff");
   eff_etaMC2->Draw("14.9>>+hmcT", cutD+TCut("deta!=1")+Covr,"goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_g_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("dTh>>hmcS", cutD+TCut("deta==1"), "goff");
   eff_etaMC2->Draw("14.9>>+hmcS", cutD+TCut("deta==1")+Covr,"goff");
   hmc[1]->Scale( par->W_g_eta() );

   // sum of signal and bg.
   hmc[0]->Add( hmc[1], hmc[2] );

   // normalization on DATA
   double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral();
   double sum = hmc[0]->Integral();
   printf("%i %s, bkg: %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100);

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   auto name = Form("c1_dTh_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,Cy/2,Cx,Cy);
   c1->cd();
   gPad->SetLogy(true);

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetMinimum(Ymin);
   hdat[0]->SetTitle(
         ";#delta#Theta(#gamma), deg."
         ";Entries / 0.2 deg.");
   // hdat[0]->GetXaxis()->SetNdivisions(1004);
   hdat[0]->GetYaxis()->SetMaxDigits(2);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(8., Ymin, Xmax, ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   // TLegend* leg = new TLegend(0.61,0.65,0.892,0.89);
   TLegend* leg = new TLegend(0.60,0.65,0.882,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   // leg->AddEntry(box, "Rejection area","F");
   leg->AddEntry(box, "Mismatch area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_dTh_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot of Mgg
//--------------------------------------------------------------------
void plot_Mgg(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   // Filling histograms from ntuple
   double Xmin = 0.508, Xmax = 0.588;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 80,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   auto& meta = par->meta;
   const double weta = 3*0.008; // see cuts.h
   double Mgg_min = meta-weta, Mgg_max = meta+weta;

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cbg + par->Ceta + TCut("dTh<8.&&0.8<rE&&rE<1.4");
   eff_etaD->Draw("mggf>>hdat", cutD, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(4,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("mggf>>hmcB", cutD+!(par->Cmcsig), "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC signal
   // add to background g-eta but !(eta->2gamma)
   eff_etaMC2->Draw("mggf>>hmcT", cutD+TCut("deta!=1"), "goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_g_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("mggf>>hmcS", cutD+TCut("deta==1"), "goff");
   hmc[1]->Scale( par->W_g_eta() );

   hmc[3] = hst("hmcSTT"); // to estimate incorect gamma matching
   eff_etaMC2->Draw("mggf>>hmcSF", cutD+TCut("dgam!=1"), "goff");
   hmc[3]->Scale( par->W_g_eta() );

   // sum of signal and bg.
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Mgg_min+bW/2);
   int bin2 = hdat[0]->FindBin(Mgg_max-bW/2);
   // cout << "bW=" << bW
      // << " bin1=" << bin1 << " bin2=" << bin2 << endl;

   // normalization on DATA
   double scale = hdat[0]->Integral(bin1,bin2) /
      hmc[0]->Integral(bin1,bin2);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral(bin1,bin2);
   double sum = hmc[0]->Integral(bin1,bin2);
   double bg_g = hmc[3]->Integral(bin1,bin2);
   printf("%i %s, bkg(in eta): %.2f+/-%.2f%%, "
         "wrong gamma %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100,
         bg_g/sum*100, sqrt(bg_g)/sum*100);

   // fit function
   TF1* gs = (TF1*)gROOT->GetFunction("gaus");
   // TF1* gs = new TF1("gs","gaus+[3]",Mgg_min,Mgg_max);
   // gs->SetParameters(1.,meta,0.008,0.);
   // gs->SetParNames("Constant","Mean","Sigma","Bkg const");
   gs->SetLineWidth(2);
   gs->SetLineColor(kGreen+2);

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   gStyle->SetStatX(0.892);
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
   hdat[0]->GetYaxis()->SetTitleOffset(1.15);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->SetMarkerStyle(20);
   double ymax = 1.12 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Fit(gs,"QM","",meta-2*0.008,meta+2*0.008);
   hdat[0]->Draw("E");

   box->DrawBox( Xmin,    0., Mgg_min, ymax);
   box->DrawBox( Mgg_max, 0., Xmax,    ymax);
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
   // leg->AddEntry(box, "Rejection area","F");
   leg->AddEntry(box, "Mismatch area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaeff_Mgg_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot predicted P(eta) and cos(Theta(eta))
//--------------------------------------------------------------------
void plot_Peta(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   // Filling histograms from ntuple
   auto hst = [](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 55,1.25,1.8);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cbg;
   eff_etaD->Draw("Peta>>hdat", cutD, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(4,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("Peta>>hmcB", cutD+!(par->Cmcsig), "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC signal
   // add to background g-eta but !(eta->2gamma)
   eff_etaMC2->Draw("Peta>>hmcT", cutD+TCut("deta!=1"), "goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_g_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("Peta>>hmcS", cutD+TCut("deta==1"), "goff");
   hmc[1]->Scale( par->W_g_eta() );

   hmc[3] = hst("hmcSTT"); // to estimate incorect gamma matching
   eff_etaMC2->Draw("Peta>>hmcSTT", cutD+TCut("dgam!=1"), "goff");
   hmc[3]->Scale( par->W_g_eta() );

   // sum of signal and bg.
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Peta_min = 1.3, Peta_max = 1.7;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Peta_min+bW/2);
   int bin2 = hdat[0]->FindBin(Peta_max-bW/2);
   // cout << "bW=" << bW
      // << " bin1=" << bin1 << " bin2=" << bin2 << endl;

   // normalization on DATA
   double scale = hdat[0]->Integral(bin1,bin2) /
      hmc[0]->Integral(bin1,bin2);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral(bin1,bin2);
   double sum = hmc[0]->Integral(bin1,bin2);
   double bg_g = hmc[3]->Integral(bin1,bin2);
   printf("%i %s, bkg: %.2f+/-%.2f%%, "
         "wrong gamma %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100,
         bg_g/sum*100, sqrt(bg_g)/sum*100);

   // Draw:
   auto name = Form("c1_peta_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetTitle(
         ";P_{#eta}, GeV/c"
         ";Entries/0.01 GeV/c");
   hdat[0]->GetYaxis()->SetMaxDigits(3);
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

   // incorect gamma matching
   // hmc[3]->SetLineWidth(2);
   // hmc[3]->SetLineStyle(kDashed);
   // hmc[3]->SetLineColor(kBlue+3);
   // hmc[3]->Draw("HIST,SAME"); // wrong gamma

   TLegend* leg = new TLegend(0.61,0.69,0.892,0.89);
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
void plot_Ceta(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew

   // Filling histograms from ntuple
   int Nbins = 40;
   auto hst = [Nbins](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", Nbins,-1.,1.);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cbg + par->Ceta;
   eff_etaD->Draw("Ceta>>hdat", cutD, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(4,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("Ceta>>hmcB", cutD+!(par->Cmcsig), "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC signal
   // add to background g-eta but !(eta->2gamma)
   eff_etaMC2->Draw("Ceta>>hmcT", cutD+TCut("deta!=1"), "goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_g_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("Ceta>>hmcS", cutD+TCut("deta==1"), "goff");
   hmc[1]->Scale( par->W_g_eta() );

   hmc[3] = hst("hmcSTT"); // to estimate incorect gamma matching
   eff_etaMC2->Draw("Ceta>>hmcSTT", cutD+TCut("dgam!=1"), "goff");
   hmc[3]->Scale( par->W_g_eta() );

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Ceta_min = -1., Ceta_max = 1.;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Ceta_min+bW/2);
   int bin2 = hdat[0]->FindBin(Ceta_max-bW/2);
   // cout << "bW=" << bW
      // << " bin1=" << bin1 << " bin2=" << bin2 << endl;

   // normalization on DATA
   double scale = hdat[0]->Integral(bin1,bin2) /
      hmc[0]->Integral(bin1,bin2);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral(bin1,bin2);
   double sum = hmc[0]->Integral(bin1,bin2);
   double bg_g = hmc[3]->Integral(bin1,bin2);
   printf("%i %s, bkg: %.2f+/-%.2f%%, "
         "wrong gamma %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100,
         bg_g/sum*100, sqrt(bg_g)/sum*100);

   // Draw:
   auto name = Form("c1_ceta_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,Cy/2,Cx,Cy);
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

   // incorect gamma matching
   // hmc[3]->SetLineWidth(2);
   // hmc[3]->SetLineStyle(kDashed);
   // hmc[3]->SetLineColor(kBlue+3);
   // hmc[3]->Draw("HIST,SAME"); // wrong gamma

   TLegend* leg = new TLegend(0.61,0.69,0.892,0.89);
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
void plot2_MKK(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   // Here we use histograms filled in the PsipJpsiPhiEta: EtaEff()
   // It is incorrect to use ntuple here.

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
   double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   // double scale = hdat[0]->Integral(49,69) / hmc[0]->Integral(49,69);
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
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   // Data
   hdat[0]->SetAxisRange(0.98,1.0999,"X");
   hdat[0]->SetTitle(
         ";M^{2}_{inv}(K^{#plus}K^{#minus}), GeV^{2}/c^{4}"
         ";Entries/0.001 GeV^{2}/c^{4}");
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.1 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(0.98,0.,1.02,ymax);
   box->DrawBox(1.06,0.,1.10,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   // hmc[2]->Draw("HIST,SAME"); // bkg is phi+X, it is too big

   TLegend* leg = new TLegend(0.61,0.65,0.892,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   // leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaphi_Mkk2_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot invariant mass of photon and «missing photon»
//--------------------------------------------------------------------
void plot2_Minv2g(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   // Here we use histograms filled in the PsipJpsiPhiEta: EtaEff()
   // M2gg not in tuple yet

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
   // double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
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
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   // Data
   hdat[0]->SetAxisRange(0.198,0.4019,"X");
   hdat[0]->SetTitle(
         ";M^{2}_{inv}(#gamma#gamma_{missing} ), GeV^{2}/c^{4}"
         ";Entries/0.003 GeV^{2}/c^{4}");
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
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

   TLegend* leg = new TLegend(0.61,0.65,0.892,0.89);
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
void plot2_M2mis(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   // Filling histograms from ntuple
   double Xmin = 0., Xmax = 0.008;
   double Ymin = 0.5; // log
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 80,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Ceta;
   eff_etaD->Draw("m2fr>>hdat", cutD, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("m2fr>>hmcB", cutD+!(par->Cmcsig), "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-phi-eta: signal
   // add to background phi-eta but !(eta->2gamma)
   eff_etaMC2->Draw("m2fr>>hmcT", cutD+TCut("deta!=1"), "goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_phi_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("m2fr>>hmcS", cutD+TCut("deta==1"), "goff");
   hmc[1]->Scale( par->W_phi_eta() );

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Mmis2_min = 0., Mmis2_max = 0.001;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = 1;
   int bin2 = hdat[0]->FindBin(Mmis2_max-bW/2);
   // cout << "bW=" << bW
      // << " bin1=" << bin1 << " bin2=" << bin2 << endl;

   // normalization on DATA
   double scale = hdat[0]->Integral(bin1,bin2) /
      hmc[0]->Integral(bin1,bin2);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral(bin1,bin2);
   double sum = hmc[0]->Integral(bin1,bin2);
   printf("%i %s, bkg: %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100);

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetMinimum(Ymin);
   hdat[0]->SetTitle( ";M^{2}_{missing}"
         "(#pi^{#plus}#pi^{#minus}K^{#plus}K^{#minus}#gamma),"
         " GeV^{2}/c^{4}"
         ";Entries/0.0002 GeV^{2}/c^{4}");
   hdat[0]->GetXaxis()->SetNdivisions(1005);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(Mmis2_max,Ymin, Xmax,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.61,0.65,0.892,0.89);
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
void plot2_rE(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   // Filling histograms from ntuple
   double Xmin = 0.4, Xmax = 1.8;
   double Ymin = 0.5; // log-Y
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 70,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cbg + par->Ceta;
   TCut Covr( Form("rE<%f||%f<=rE",Xmin,Xmax) ); // under & overflows
   eff_etaD->Draw("rE>>hdat", cutD, "goff");
   eff_etaD->Draw("1.79>>+hdat", cutD+Covr, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("rE>>hmcB", cutD+!(par->Cmcsig), "goff");
   eff_etaMC->Draw("1.79>>+hmcB", cutD+!(par->Cmcsig)+Covr, "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   // add to background g-eta but !(eta->2gamma)
   eff_etaMC2->Draw("rE>>hmcT", cutD+TCut("deta!=1"), "goff");
   eff_etaMC2->Draw("1.79>>+hmcT",cutD+TCut("deta!=1")+Covr,"goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_phi_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("rE>>hmcS", cutD+TCut("deta==1"), "goff");
   eff_etaMC2->Draw("1.79>>+hmcS",cutD+TCut("deta==1")+Covr,"goff");
   hmc[1]->Scale( par->W_phi_eta() );

   // sum of signal and bg.
   hmc[0]->Add( hmc[1], hmc[2] );

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
   auto name = Form("c1_eE_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetLogy(true);

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetMinimum(Ymin);
   hdat[0]->SetTitle(
         ";E_{#gamma}(pred) / E_{#gamma}(rec)"
         ";Entries/0.01");
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.15);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(Xmin,Ymin, 0.8, ymax);
   box->DrawBox(1.4, Ymin, Xmax,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   // TLegend* leg = new TLegend(0.61,0.65,0.892,0.89);
   TLegend* leg = new TLegend(0.60,0.65,0.882,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   // leg->AddEntry(box, "Rejection area","F");
   leg->AddEntry(box, "Mismatch area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaphi_rE_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot2_dTh(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   // Filling histograms from ntuple
   double Xmin = 0., Xmax = 15.;
   double Ymin = 0.5; // log-Y
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 75,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cbg + par->Ceta;
   TCut Covr( Form("%f<=dTh",Xmax) ); // under & overflows
   eff_etaD->Draw("dTh>>hdat", cutD, "goff");
   eff_etaD->Draw("14.9>>+hdat", cutD+Covr, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("dTh>>hmcB", cutD+!(par->Cmcsig), "goff");
   eff_etaMC->Draw("14.9>>+hmcB", cutD+!(par->Cmcsig)+Covr, "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC-gamma-eta: signal
   // add to background g-eta but !(eta->2gamma)
   eff_etaMC2->Draw("dTh>>hmcT", cutD+TCut("deta!=1"), "goff");
   eff_etaMC2->Draw("14.9>>+hmcT", cutD+TCut("deta!=1")+Covr,"goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_phi_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("dTh>>hmcS", cutD+TCut("deta==1"), "goff");
   eff_etaMC2->Draw("14.9>>+hmcS", cutD+TCut("deta==1")+Covr,"goff");
   hmc[1]->Scale( par->W_phi_eta() );

   // sum of signal and bg.
   hmc[0]->Add( hmc[1], hmc[2] );

   // normalization on DATA
   double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral();
   double sum = hmc[0]->Integral();
   printf("%i %s, bkg: %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100);

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   auto name = Form("c1_dTh_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,Cy/2,Cx,Cy);
   c1->cd();
   gPad->SetLogy(true);

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetMinimum(Ymin);
   hdat[0]->SetTitle(
         ";#delta#Theta(#gamma), deg."
         ";Entries/0.2 deg.");
   // hdat[0]->GetXaxis()->SetNdivisions(1004);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(8., Ymin, Xmax, ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   // TLegend* leg = new TLegend(0.61,0.65,0.892,0.89);
   TLegend* leg = new TLegend(0.60,0.65,0.882,0.89);
   leg->SetHeader( Form("%i",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   // leg->AddEntry(box, "Rejection area","F");
   leg->AddEntry(box, "Mismatch area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaphi_dTh_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot of Mgg
//--------------------------------------------------------------------
void plot2_Mgg(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   Params* par = new Params(date,-2,0); // date, phieta, no_rew

   double Xmin = 0.508, Xmax = 0.588;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 80,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   auto& meta = par->meta;
   const double weta = 3*0.008; // see cuts.h
   double Mgg_min = meta-weta, Mgg_max = meta+weta;

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cbg + par->Ceta + TCut("dTh<8.&&0.8<rE&&rE<1.4");
   eff_etaD->Draw("mggf>>hdat", cutD, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("mggf>>hmcB", cutD+!par->Cmcsig, "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC signal
   // add to background !(eta->2gamma)
   eff_etaMC2->Draw("mggf>>hmcT", cutD+TCut("deta!=1"), "goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_phi_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("mggf>>hmcS", cutD+TCut("deta==1"), "goff");
   hmc[1]->Scale( par->W_phi_eta() );

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Mgg_min+bW/2);
   int bin2 = hdat[0]->FindBin(Mgg_max-bW/2);
   // cout << "bW=" << bW
      // << " bin1=" << bin1 << " bin2=" << bin2 << endl;

   // normalization on DATA
   double scale = hdat[0]->Integral(bin1,bin2) /
      hmc[0]->Integral(bin1,bin2);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral(bin1,bin2);
   double sum = hmc[0]->Integral(bin1,bin2);
   printf("%i %s, bkg(in eta): %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100);

   // fit function
   TF1* gs = (TF1*)gROOT->GetFunction("gaus");
   // TF1* gs = new TF1("gs","gaus+[3]",Mgg_min,Mgg_max);
   // gs->SetParameters(1.,meta,0.008,0.);
   // gs->SetParNames("Constant","Mean","Sigma","Bkg const");
   gs->SetLineWidth(2);
   gs->SetLineColor(kGreen+2);

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   gStyle->SetStatX(0.892);
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
   hdat[0]->GetYaxis()->SetTitleOffset(1.15);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->SetMarkerStyle(20);
   double ymax = 1.12 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Fit(gs,"QM","",meta-2*0.008,meta+2*0.008);
   hdat[0]->Draw("E");

   box->DrawBox( Xmin,    0., Mgg_min, ymax);
   box->DrawBox( Mgg_max, 0., Xmax,    ymax);
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
   // leg->AddEntry(box, "Rejection area","F");
   leg->AddEntry(box, "Mismatch area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "etaphi_Mgg_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot predicted P(eta) and cos(Theta(eta))
//--------------------------------------------------------------------
void plot2_Peta(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
   // date, eta_eff in phieta, no_rew
   Params* par = new Params(date,-2,0);

   auto hst = [](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 45,1.1,1.55);
      h->Sumw2(true);
      return h;
   };

   TTree* eff_etaD = par->GetEffEta(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cbg;
   eff_etaD->Draw("Peta>>hdat", cutD, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("Peta>>hmcB", cutD+!par->Cmcsig, "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC signal
   // add to background g-eta but !(eta->2gamma)
   eff_etaMC2->Draw("Peta>>hmcT", cutD+TCut("deta!=1"), "goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_phi_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("Peta>>hmcS", cutD+TCut("deta==1"), "goff");
   hmc[1]->Scale( par->W_phi_eta() );

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Peta_min = 1.15, Peta_max = 1.5;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Peta_min+bW/2);
   int bin2 = hdat[0]->FindBin(Peta_max-bW/2);
   // cout << "bW=" << bW
      // << " bin1=" << bin1 << " bin2=" << bin2 << endl;

   // normalization on DATA
   double scale = hdat[0]->Integral(bin1,bin2) /
      hmc[0]->Integral(bin1,bin2);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral(bin1,bin2);
   double sum = hmc[0]->Integral(bin1,bin2);
   printf("%i %s, bkg: %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100);

   // Draw:
   auto name = Form("c1_peta_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetTitle(
         ";P_{#eta}, GeV/c"
         ";Entries/0.01 GeV/c");
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
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


   TLegend* leg = new TLegend(0.61,0.71,0.892,0.89);
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
void plot2_Ceta(int date, int Cx = 800, int Cy = 800)
//--------------------------------------------------------------------
{
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

   TCut cutD = par->Cbg + par->Ceta;
   eff_etaD->Draw("Ceta>>hdat", cutD, "goff");

   TTree* eff_etaMC = par->GetEffEta(1); // MC: bg only!
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[2] = hst("hmcB");
   eff_etaMC->Draw("Ceta>>hmcB", cutD+!par->Cmcsig, "goff");

   TTree* eff_etaMC2 = par->GetEffEta(2); // MC signal
   // add to background !(eta->2gamma)
   eff_etaMC2->Draw("Ceta>>hmcT", cutD+TCut("deta!=1"), "goff");
   hmc[2]->Add( hmc[2], hmc[0], 1., par->W_phi_eta() );

   hmc[1] = hst("hmcS");
   eff_etaMC2->Draw("Ceta>>hmcS", cutD+TCut("deta==1"), "goff");
   hmc[1]->Scale( par->W_phi_eta() );

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Ceta_min = -1., Ceta_max = 1.;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Ceta_min+bW/2);
   int bin2 = hdat[0]->FindBin(Ceta_max-bW/2);
   // cout << "bW=" << bW
      // << " bin1=" << bin1 << " bin2=" << bin2 << endl;

   // normalization on DATA
   double scale = hdat[0]->Integral(bin1,bin2) /
      hmc[0]->Integral(bin1,bin2);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // background estimation
   double bg = hmc[2]->Integral(bin1,bin2);
   double sum = hmc[0]->Integral(bin1,bin2);
   printf("%i %s, bkg: %.2f+/-%.2f%%\n",
         date, __func__, bg/sum*100, sqrt(bg)/sum*100);

   // Draw:
   auto name = Form("c1_ceta_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,Cy/2,Cx,Cy);
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

   TLegend* leg = new TLegend(0.61,0.71,0.892,0.89);
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
void eta_eff_sel()
//--------------------------------------------------------------------
{
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(42);

   size_t Cx = 800, Cy = 640; // canvas sizes, X/Y = 1.25

   // for ( auto date : {2021} ) {
   for ( auto date : {2009, 2012, 2021} ) {
      // ++ pi0 rejection ++
      // plot_pi0(date,Cx,Cy); // fig.B2

      // ++ J/Psi -> gamma eta ++
      // plot_Mpipig(date,Cx,Cy); // fig B3
      // plot_Minv2g(date,Cx,Cy); // fig B4
      //--  plot_Minv2g_n(date,Cx,Cy); // fig B4 -> no M2gg in data ---
      // plot_M2mis(date,Cx,Cy);  // fig B5
      // plot_Eg(date,Cx,Cy);     // fig B6
      // plot_dTh(date,Cx,Cy);    // fig B7 (top)
      // plot_rE(date,Cx,Cy);     // fig B7 (bottom)
      // plot_Mgg(date,Cx,Cy);    // fig B8
      // plot_Peta(date,Cx,Cy);   // fig B9 (top)
      // plot_Ceta(date,Cx,Cy);   // fig B9 (bottom)

      // -- J/Psi -> phi eta ++
      // plot2_MKK(date,Cx,Cy);       // fig B12
      // plot2_Minv2g(date,Cx,Cy);    // fig B13
      // plot2_M2mis(date,Cx,Cy);     // fig B14
      // plot2_dTh(date,Cx,Cy);       // fig B15 (top)
      // plot2_rE(date,Cx,Cy);        // fig B15 (bottom)
      // plot2_Mgg(date,Cx,Cy);       // fig B16
      // plot2_Peta(date,Cx,Cy);      // fig B17 (top)
      // plot2_Ceta(date,Cx,Cy);      // fig B17 (bottom)
   }
}
