// trk_eff_sel.cc
// Study of the efficiency of reconstruction of pion and kaon tracks.
// Pictures for event selection.

// {{{1 Common parameters: Params
//--------------------------------------------------------------------
struct Params {
   Params(int dat, int kpi, int pm, int rew);

   TFile* OpenFile(int mc);
   TTree* GetEff(int mc);
   const char* Sdate() { return Form("%i",date); }

   // name of folder with root files
   const string Dir = "prod_v709n4/Eff/";
   string datafile;
   string mcincfile;

   int date;    // 2009 ...
   int use_kpi; // 1: kaons; 2: pions
   int use_pm;  //  0: sum "+" and "-"
                //  1: "+" only
                // -1: "-" only
   int use_rew; // 0 - no weights;
                // 1 - calculate weights

   const double mk   = 0.49368;  // 493.677   +/- 0.015 MeV
   const double mpi  = 0.13957;  // 139.57039 +/- 0.00018 MeV

   // cuts for ntuples
   TCut Cmcsig; // mc-signal

   TCut Ckaon;  // std cuts for kaons
   TCut CKPtC;  // cuts for Pt and cos(Theta) of kaons
   TCut CKf;    // predicted kaon found

   TCut Cpion;  // std cuts for pions
   TCut CPiPtC; // cuts for Pt and cos(Theta) of pions
   TCut CPif;   // predicted pion found
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

   // *** MC signal ***
   Cmcsig = TCut("good==1");
   // Cmcsig = TCut("abs(good-1)<0.1"); // < n4_eff

   // *** cuts for kaons *** (mk**2 = 0.243717 GeV**2)
   Ckaon = TCut("abs(Mrec-3.097)<0.002") +
      TCut( Form("abs(Mrk2-%.6f)<0.025",mk*mk) );
   // limits for Pt and cos(Theta) of kaons
   CKPtC = TCut("abs(Ck)<0.8") + TCut("0.05<Ptk&&Ptk<1.45");
   // predicted kaon found (we do not use it here)
   // CKf   = TCut("dTh<10&&(-0.12<dP&&dP<0.08)&&Nk==2");

   // *** cuts for pions *** (mpi**2 = 0.0194798 GeV**2)
   // Cpion = TCut("abs(MppKK-3.097)<0.01") +
   Cpion = TCut("abs(MpKK-3.097)<0.01") +
      TCut( Form("abs(Mrpi2-%.6f)<0.025",mpi*mpi) );
   // limits for Pt and cos(Theta) of pions
   CPiPtC = TCut("abs(Cpi)<0.8") + TCut("0.025<Ptpi&&Ptpi<0.425");
   // predicted pion found (we do not use it here)
   // CPif   = TCut("dTh<15&&(-0.12<dP&&dP<0.08)&&Npi==4");

}

// {{{2 > OpenFile()
//--------------------------------------------------------------------
TFile* Params::OpenFile(int mc) {  // 1 for MC
//--------------------------------------------------------------------
   string dfname = Dir + ( (mc!=1) ? datafile : mcincfile );
   TFile* froot = TFile::Open(dfname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << dfname << endl;
      exit(EXIT_FAILURE);
   }
   froot->cd("PipPimKpKm");
   return froot;
}

// {{{2 > GetEff()
//--------------------------------------------------------------------
TTree* Params::GetEff(int mc = 0) { // get tree for K or Pi
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
void plot_pi0(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date);

   vector<string> hndat {"S1b_Mg2"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S1b_Mg2", "MC_M2ggT", "MC_M2ggF"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0]->Integral(60,109)/hmc[0]->Integral(60,109);
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
   auto name = Form("c1_pi0_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0]->SetAxisRange(0.005,0.03,"X");
   double bW = hdat[0]->GetBinWidth(1);
   hdat[0]->SetTitle( Form(";M^{2}_{#gamma#gamma} , GeV^{2}/c^{4}"
         ";Entries/%.4f GeV^{2}/c^{4}",bW) ); // 0.0002
   if ( date > 2009 ) {
      hdat[0]->GetYaxis()->SetMaxDigits(3);
   }
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   double ymax = 1.05*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(0.012,.0,0.022,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[1]->SetLineWidth(1);
   hmc[1]->SetLineColor(kGreen+1);
   hmc[1]->SetFillStyle(3001);
   hmc[1]->SetFillColor(kGreen+1);
   hmc[1]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.59,0.69,0.892,0.89);
   leg->SetHeader( par->Sdate(),"C");
   leg->AddEntry(hdat[0], "Data","LEP`");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC signal","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_pi0_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 ++++ K-eff: Plot MmisK
//--------------------------------------------------------------------
void plot_MmisK(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date);

   // Filling histograms from ntuple
   double Xmin = SQ(par->mk)-0.06, Xmax = SQ(par->mk)+0.06;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 96,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = TCut("abs(Mrec-3.097)<0.002");
   effD->Draw("Mrk2>>hdat", cutD, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("Mrk2>>hmcS", cutD+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("Mrk2>>hmcB", cutD+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double M2l = SQ(par->mk)-0.025;
   double M2r = SQ(par->mk)+0.025;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(M2l+bW/2);
   int bin2 = hdat[0]->FindBin(M2r-bW/2);
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
   auto name = Form("c1_MmisK_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetTitle( Form(
         ";M^{2}_{recoil}(2(#pi^{+}#pi^{-})K^{#pm}), GeV^{2}/c^{4}"
         ";Entries/%.4f GeV^{2}/c^{4}",bW)
         ); // 0.0012
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   double ymax = 1.06*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(Xmin,0., M2l, ymax);
   box->DrawBox(M2r, 0., Xmax,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.64,0.67,0.892,0.89);
   leg->SetHeader( par->Sdate(),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_MmK_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot PtKp(Km) CKp(Km)
//--------------------------------------------------------------------
void plot_PtKp(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date,1); // Kaons

   // Filling histograms from ntuple
   double Xmin = 0., Xmax = 1.5;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 60,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Ckaon + TCut("Zk>0")
      + TCut("abs(Ck)<0.8");
   effD->Draw("Ptk>>hdat", cutD, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("Ptk>>hmcS", cutD+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("Ptk>>hmcB", cutD+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Ptmin = 0.05;
   double Ptmax = 1.45;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Ptmin+bW/2);
   int bin2 = hdat[0]->FindBin(Ptmax-bW/2);
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

   // Draw:
   auto name = Form("c1_ptkp_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetTitle( Form(";P_{t}(K), GeV/c"
         ";Entries/%.3f GeV/c",bW)
         ); // 0.025
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
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

   TLegend* leg = new TLegend(0.59,0.69,0.892,0.89);
   leg->SetHeader( Form("%i  K^{#plus}",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_PtKp_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot_PtKm(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date, 1); // Kaons

   // Filling histograms from ntuple
   double Xmin = 0., Xmax = 1.5;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 60,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Ckaon + TCut("Zk<0")
      + TCut("abs(Ck)<0.8");
   effD->Draw("Ptk>>hdat", cutD, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("Ptk>>hmcS", cutD+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("Ptk>>hmcB", cutD+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Ptmin = 0.05;
   double Ptmax = 1.45;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Ptmin+bW/2);
   int bin2 = hdat[0]->FindBin(Ptmax-bW/2);
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

   // Draw:
   auto name = Form("c1_ptkm_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,Cy/2,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetTitle( Form(";P_{t}(K), GeV/c"
         ";Entries/%.3f GeV/c",bW)
         ); // 0.025
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
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

   TLegend* leg = new TLegend(0.59,0.69,0.892,0.89);
   leg->SetHeader( Form("%i  K^{#minus}",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_PtKm_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot_CKp(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date,1); // Kaons

   // Filling histograms from ntuple
   double Xmin = -1., Xmax = 1.;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 40,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Ckaon + TCut("Zk>0")
      + TCut("0.05<Ptk&&Ptk<1.45");
   effD->Draw("Ck>>hdat", cutD, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("Ck>>hmcS", cutD+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("Ck>>hmcB", cutD+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Thmin = -0.8;
   double Thmax = +0.8;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Thmin+bW/2);
   int bin2 = hdat[0]->FindBin(Thmax-bW/2);
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

   // Draw:
   auto name = Form("c1_Ckp_%i",date);
   TCanvas* c1 = new TCanvas(name,name,Cx/2,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->SetTitle( Form(";cos(#Theta_{K}) "
         ";Entries/%.2f ", bW)
         ); // 0.05
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
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

   TLegend* leg = new TLegend(0.59,0.69,0.892,0.89);
   leg->SetHeader( Form("%i  K^{#plus}",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_CKp_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot_CKm(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date,1);

   // Filling histograms from ntuple
   double Xmin = -1., Xmax = 1.;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 40,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Ckaon + TCut("Zk<0")
      + TCut("0.05<Ptk&&Ptk<1.45");
   effD->Draw("Ck>>hdat", cutD, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("Ck>>hmcS", cutD+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("Ck>>hmcB", cutD+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Thmin = -0.8;
   double Thmax = +0.8;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Thmin+bW/2);
   int bin2 = hdat[0]->FindBin(Thmax-bW/2);
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

   // Draw:
   auto name = Form("c1_Ckm_%i",date);
   TCanvas* c1 = new TCanvas(name,name,Cx/2,Cy/2,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->SetTitle( Form(";cos(#Theta_{K}) "
         ";Entries/%.2f ", bW)
         ); // 0.05
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
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

   TLegend* leg = new TLegend(0.59,0.69,0.892,0.89);
   leg->SetHeader( Form("%i  K^{#minus}",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_CKm_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot dPk dThK
//--------------------------------------------------------------------
void plot_dPK(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date);

   // Filling histograms from ntuple
   double Xmin = -0.15, Xmax = +0.15;
   double Ymin = 0.5; // log-Y
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 60,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };
   double bW = hdat[0]->GetBinWidth(1);

   TCut cutD = par->Ckaon + par->CKPtC;
   TCut Covr( Form("dP<%f||%f<=dP",Xmin,Xmax) ); // under & overflows
   effD->Draw("dP>>hdat", cutD, "goff");
   effD->Draw(Form("%f>>+hdat",Xmax-bW/2), cutD+Covr, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("dP>>hmcS", cutD+par->Cmcsig, "goff");
   effMC->Draw(Form("%f>>+hmcS",Xmax-bW/2),
         cutD+Covr+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("dP>>hmcB", cutD+!par->Cmcsig, "goff");
   effMC->Draw(Form("%f>>+hmcB",Xmax-bW/2),
         cutD+Covr+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
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
   auto name = Form("c1_dPK_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetMinimum(Ymin);
   hdat[0]->SetTitle( Form(";#deltaP(K), GeV/c"
         ";Entries/%.3f GeV/c",bW)
         ); // 0.005
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   // double ymax=1.06*hdat[0]->GetMaximum(); // Lin
   double ymax = 1.9*hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(Xmin,Ymin,-0.12,ymax);
   box->DrawBox(0.08,Ymin, Xmax,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.59,0.69,0.882,0.89);
   leg->SetHeader( par->Sdate(), "C" );
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Mismatch area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_dPK_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot_dThK(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date);

   // Filling histograms from ntuple
   double Xmin = 0., Xmax = 15;
   double Ymin = 0.5; // log-Y
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 75,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };
   double bW = hdat[0]->GetBinWidth(1);

   TCut cutD = par->Ckaon + par->CKPtC;
   TCut Covr( Form("%f<=dTh",Xmax) ); // under & overflows
   effD->Draw("dTh>>hdat", cutD, "goff");
   effD->Draw(Form("%f>>+hdat",Xmax-bW/2), cutD+Covr, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("dTh>>hmcS", cutD+par->Cmcsig, "goff");
   effMC->Draw(Form("%f>>+hmcS",Xmax-bW/2),
         cutD+Covr+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("dTh>>hmcB", cutD+!par->Cmcsig, "goff");
   effMC->Draw(Form("%f>>+hmcB",Xmax-bW/2),
         cutD+Covr+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
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
   auto name = Form("c1_dThK_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,Cy/2,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetMinimum(Ymin);
   hdat[0]->SetTitle( Form(";#delta#Theta(K), deg"
         ";Entries/%.1f deg",bW)
         ); // 0.2
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   // double ymax=1.06*hdat[0]->GetMaximum(); // Lin
   double ymax=1.9*hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(10.,Ymin, Xmax,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.59,0.69,0.882,0.89);
   leg->SetHeader( par->Sdate(), "C" );
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Mismatch area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_dThK_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 ++++ Pi-eff: Plot MinvJ
//--------------------------------------------------------------------
void plot_MinvJ(int date, int Cx=800, int Cy=800, bool fitgs=false) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2); // Pions

   // Here we use histograms filled in the PsipPiPiKK.
   // It is incorrect to use ntuple here.

   vector<string> hndat {"S4_MppKKb"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S4_MppKKb", "MC_MppKK_T", "MC_MppKK_F"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // recalculation: Min,Max -> bin1,bin2
   // v709n4: [3.087,3.107]
   double Ml = 3.087;
   double Mr = 3.107;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Ml+bW/2);
   int bin2 = hdat[0]->FindBin(Mr-bW/2);
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
   auto name = Form("c1_MinvJ_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   double Xmin = 3.0617, Xmax = 3.1328;
   hdat[0]->SetAxisRange(Xmin,Xmax,"X");
   hdat[0]->SetTitle( Form(
         ";M^{ inv}_{#pi^{#plus}#pi^{#minus}K^{#plus}K^{#minus}},"
         " GeV/c^{2}"
         ";Entries/%.4f GeV/c^{2}",bW)
         ); // 0.0009
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   double ymax = 1.06 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);

   if ( fitgs ) {
      gStyle->SetStatX(0.892);
      gStyle->SetStatY(0.89);
      gStyle->SetStatW(0.16);
      gStyle->SetFitFormat(".4f");
      TF1* gs = (TF1*)gROOT->GetFunction("gaus");
      gs->SetLineWidth(2);
      gs->SetLineStyle(5);
      gs->SetLineColor(kMagenta);
      hdat[0]->Fit(gs,"Q","",3.092,3.102); // just at peak
   } else {
      hdat[0]->Draw("E");
   }

   box->DrawBox(Xmin,0., Ml,  ymax);
   box->DrawBox(Mr,  0., Xmax,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = nullptr;
   if ( fitgs ) {
      leg = new TLegend(0.12,0.69,0.38,0.89);
   } else {
      leg = new TLegend(0.61,0.69,0.892,0.89);
   }
   leg->SetHeader( par->Sdate(), "C" );
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_MiJ_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot MmisP
//--------------------------------------------------------------------
void plot_MmisP(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2); // Pions

   // Filling histograms from ntuple
   double Xmin = SQ(par->mpi)-0.02, Xmax = SQ(par->mpi)+0.02;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 80,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   // TCut cutD = TCut("abs(MppKK-3.097)<0.01");
   TCut cutD = TCut("abs(MpKK-3.097)<0.01");
   effD->Draw("Mrpi2>>hdat", cutD, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("Mrpi2>>hmcS", cutD+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("Mrpi2>>hmcB", cutD+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double M2l = SQ(par->mpi)-0.01;
   double M2r = SQ(par->mpi)+0.01;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(M2l+bW/2);
   int bin2 = hdat[0]->FindBin(M2r-bW/2);
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
   auto name = Form("c1_MmisP_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetTitle( Form(";M^{2}_{recoil}"
         "(#pi^{#pm}#pi^{#plus}#pi^{#minus}K^{#plus}K^{#minus}),"
         " GeV^{2}/c^{4}"
         ";Entries/%.4f GeV^{2}/c^{4}",bW)
         ); //0.0005
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   double ymax = 1.07*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(Xmin,0., M2l, ymax);
   box->DrawBox(M2r, 0., Xmax,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.59,0.69,0.892,0.89);
   leg->SetHeader( par->Sdate(), "C" );
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_MmP_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot PtPip(Pim) CPip(Pim)
//--------------------------------------------------------------------
void plot_PtPip(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2); // Pions

   // Filling histograms from ntuple
   double Xmin = 0., Xmax = 0.45;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 36,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cpion + TCut("Zpi>0")
      + TCut("abs(Cpi)<0.8");
   effD->Draw("Ptpi>>hdat", cutD, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("Ptpi>>hmcS", cutD+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("Ptpi>>hmcB", cutD+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Ptmin = 0.025;
   double Ptmax = 0.425;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Ptmin+bW/2);
   int bin2 = hdat[0]->FindBin(Ptmax-bW/2);
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

   // Draw:
   auto name = Form("c1_ptpip_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetTitle( Form(";P_{t}(#pi), GeV/c"
         ";Entries/%.4f GeV/c",bW)
         ); // 0.0125
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
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

   // TLegend* leg = new TLegend(0.11,0.69,0.38,0.89);
   TLegend* leg = new TLegend(0.62,0.69,0.892,0.89);
   leg->SetHeader( Form("%i  #pi^{#kern[0.25]{#plus}}",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_PtPip_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot_PtPim(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2); // pi

   // Filling histograms from ntuple
   double Xmin = 0., Xmax = 0.45;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 36,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cpion + TCut("Zpi<0")
      + TCut("abs(Cpi)<0.8");
   effD->Draw("Ptpi>>hdat", cutD, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("Ptpi>>hmcS", cutD+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("Ptpi>>hmcB", cutD+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Ptmin = 0.025;
   double Ptmax = 0.425;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Ptmin+bW/2);
   int bin2 = hdat[0]->FindBin(Ptmax-bW/2);
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

   // Draw:
   auto name = Form("c1_ptpim_%i",date);
   TCanvas* c1 = new TCanvas(name,name,Cx/2,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetTitle( Form(";P_{t}(#pi), GeV/c"
         ";Entries/%.4f GeV/c",bW)
         ); // 0.0125
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
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

   TLegend* leg = new TLegend(0.62,0.69,0.892,0.89);
   leg->SetHeader( Form("%i  #pi^{#kern[0.25]{#minus}}",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_PtPim_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot_CPip(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2);

   // Filling histograms from ntuple
   double Xmin = -1., Xmax = 1.;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 40,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cpion + TCut("Zpi>0")
      + TCut("0.025<Ptpi&&Ptpi<0.425");
   effD->Draw("Cpi>>hdat", cutD, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("Cpi>>hmcS", cutD+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("Cpi>>hmcB", cutD+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Thmin = -0.8;
   double Thmax = +0.8;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Thmin+bW/2);
   int bin2 = hdat[0]->FindBin(Thmax-bW/2);
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

   // Draw:
   auto name = Form("c1_Thpip_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,Cy/2,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->SetTitle( Form(";cos(#Theta_{#pi}) "
         ";Entries/%.2f ",bW)
         ); // 0.05
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
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

   TLegend* leg = new TLegend(0.62,0.69,0.892,0.89);
   leg->SetHeader( Form("%i  #pi^{#kern[0.25]{#plus}}",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_CPip_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot_CPim(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2);

   // Filling histograms from ntuple
   double Xmin = -1., Xmax = 1.;
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 40,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };

   TCut cutD = par->Cpion + TCut("Zpi<0")
      + TCut("0.025<Ptpi&&Ptpi<0.425");
   effD->Draw("Cpi>>hdat", cutD, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("Cpi>>hmcS", cutD+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("Cpi>>hmcB", cutD+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // recalculation: Min,Max -> bin1,bin2
   double Thmin = -0.8;
   double Thmax = +0.8;
   double bW = hdat[0]->GetBinWidth(1);
   int bin1 = hdat[0]->FindBin(Thmin+bW/2);
   int bin2 = hdat[0]->FindBin(Thmax-bW/2);
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

   // Draw:
   auto name = Form("c1_Thpim_%i",date);
   TCanvas* c1 = new TCanvas(name,name,Cx/2,Cy/2,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   SetHstFace(hdat[0]);
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->SetTitle( Form(";cos(#Theta_{#pi}) "
         ";Entries/%.2f ",bW)
         ); // 0.05
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
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

   TLegend* leg = new TLegend(0.62,0.69,0.892,0.89);
   leg->SetHeader( Form("%i  #pi^{#kern[0.25]{#minus}}",date),"C");
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_CPim_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot dPpi dThPi
//--------------------------------------------------------------------
void plot_dPpi(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2); // Pions

   // Filling histograms from ntuple
   double Xmin = -0.15, Xmax = +0.15;
   double Ymin = 0.5; // log-Y
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 60,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };
   double bW = hdat[0]->GetBinWidth(1);

   TCut cutD = par->Cpion + par->CPiPtC;
   TCut Covr( Form("dP<%f||%f<=dP",Xmin,Xmax) ); // under & overflows
   effD->Draw("dP>>hdat", cutD, "goff");
   effD->Draw(Form("%f>>+hdat",Xmax-bW/2), cutD+Covr, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("dP>>hmcS", cutD+par->Cmcsig, "goff");
   effMC->Draw(Form("%f>>+hmcS",Xmax-bW/2),
         cutD+Covr+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("dP>>hmcB", cutD+!par->Cmcsig, "goff");
   effMC->Draw(Form("%f>>+hmcB",Xmax-bW/2),
         cutD+Covr+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
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
   auto name = Form("c1_dPpi_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetMinimum(Ymin);
   hdat[0]->SetTitle( Form(";#deltaP(#pi), GeV/c"
         ";Entries/%.3f GeV/c",bW)
         ); // 0.005
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   // double ymax=1.06*hdat[0]->GetMaximum(); // Lin
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(Xmin,Ymin,-0.12,ymax);
   box->DrawBox(0.08,Ymin, Xmax,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.59,0.69,0.882,0.89);
   leg->SetHeader( par->Sdate(), "C" );
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Mismatch area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_dPpi_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot_dThPi(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2); // Pions

   // Filling histograms from ntuple
   double Xmin = 0., Xmax = 20;
   double Ymin = 0.5; // log-Y
   auto hst = [Xmin,Xmax](string nm) {
      TH1D* h = new TH1D(nm.c_str(),"", 100,Xmin,Xmax);
      h->Sumw2(true);
      return h;
   };

   TTree* effD = par->GetEff(0); // data
   vector<TH1D*> hdat { hst("hdat") };
   double bW = hdat[0]->GetBinWidth(1);

   TCut cutD = par->Cpion + par->CPiPtC;
   TCut Covr( Form("%f<=dTh",Xmax) ); // under & overflows
   effD->Draw("dTh>>hdat", cutD, "goff");
   effD->Draw(Form("%f>>+hdat",Xmax-bW/2), cutD+Covr, "goff");

   TTree* effMC = par->GetEff(1); // MC
   vector<TH1D*> hmc(3,nullptr);
   hmc[0] = hst("hmcT");

   hmc[1] = hst("hmcS");
   effMC->Draw("dTh>>hmcS", cutD+par->Cmcsig, "goff");
   effMC->Draw(Form("%f>>+hmcS",Xmax-bW/2),
         cutD+Covr+par->Cmcsig, "goff");

   hmc[2] = hst("hmcB");
   effMC->Draw("dTh>>hmcB", cutD+!par->Cmcsig, "goff");
   effMC->Draw(Form("%f>>+hmcB",Xmax-bW/2),
         cutD+Covr+!par->Cmcsig, "goff");

   // sum of signal and background
   hmc[0]->Add( hmc[1], hmc[2] );

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
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
   auto name = Form("c1_dThPi_%i",date);
   TCanvas* c1 = new TCanvas(name,name,0,Cy/2,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   SetHstFace(hdat[0]);
   hdat[0]->SetMinimum(Ymin);
   hdat[0]->SetTitle( Form(";#delta#Theta(#pi), deg"
         ";Entries/%.1f deg",bW)
         ); // 0.2
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.1);
   hdat[0]->GetXaxis()->SetTitleOffset(1.1);
   // double ymax=1.06*hdat[0]->GetMaximum(); // Lin
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(15.,Ymin, Xmax,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   TLegend* leg = new TLegend(0.59,0.69,0.882,0.89);
   leg->SetHeader( par->Sdate(), "C" );
   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Mismatch area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_dThPi_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot Weights: OLD use trk_eff_wts.cc
//--------------------------------------------------------------------
void plot_Weights(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date);

   vector<string> hnmc {"MC_3_WK", "MC_5_WPi"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   TLegend* leg = new TLegend(0.59,0.79,0.89,0.89);
   leg->SetHeader( (string("#bf{")+to_string(date)+"}").c_str(),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,600);
   c1->Divide(2,1);
   // gPad->SetGrid();

   // Kaons
   c1->cd(1);
   // hmc[0]->SetAxisRange(-0.15,0.15,"X");
   hmc[0]->SetTitle(";weights for K");
   hmc[0]->SetLineWidth(2);
   hmc[0]->SetLineColor(kBlue+1);
   hmc[0]->GetYaxis()->SetMaxDigits(3);
   hmc[0]->Draw("HIST");
   leg->Draw();

   // Pions
   c1->cd(2);
   hmc[1]->SetTitle(";weights for #pi");
   hmc[1]->SetLineWidth(2);
   hmc[1]->SetLineColor(kRed+1);
   hmc[1]->GetYaxis()->SetMaxDigits(3);
   hmc[1]->Draw("HIST");
   leg->Draw();

   // leg->AddEntry(hmc[0], "Weights for K","L");
   // leg->AddEntry(hmc[1], "Weights for #pi","L");

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_WW_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Main
//--------------------------------------------------------------------
void trk_eff_sel() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(42);

   size_t Cx = 800, Cy = 640; // canvas sizes, X/Y = 1.25

   for ( auto date : {2009, 2012, 2021} ) {
   // for ( auto date : {2021} ) {
      // -- pi0 rejection, fig.A1
      // plot_pi0(date,Cx,Cy);

      // ---- K+ K- ----
      // fig A2
      // plot_MmisK(date,Cx,Cy);

      // fig A3,A4
      // plot_PtKp(date,Cx,Cy);
      // plot_PtKm(date,Cx,Cy);
      // plot_CKp(date,Cx,Cy);
      // plot_CKm(date,Cx,Cy);

      // fig A5
      // plot_dThK(date,Cx,Cy);
      // plot_dPK(date,Cx,Cy);

      // ---- pi+ pi- ----
      // fig A6
      // plot_MinvJ(date,Cx,Cy);
      // plot_MinvJ(date,Cx,Cy,true); // fit J/Psi peak position

      // fig A7
      // plot_MmisP(date,Cx,Cy);

      // fig A8,A9
      // plot_PtPip(date,Cx,Cy);
      // plot_PtPim(date,Cx,Cy);
      // plot_CPip(date,Cx,Cy);
      // plot_CPim(date,Cx,Cy);

      // fig A10
      // plot_dThPi(date,Cx,Cy);
      // plot_dPpi(date,Cx,Cy);
   }

   // OLD: do not use it
   // plot_Weights(2012);
   // plot_Weights(2009);
}
