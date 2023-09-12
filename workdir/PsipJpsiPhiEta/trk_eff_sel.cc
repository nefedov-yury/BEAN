// trk_eff_sel.cc - Pictures for presentation/memo
// Study of the track reconstruction efficiency for pi and K

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

   // std cuts for pions: mpi**2 = 0.0194798 (GeV**2)
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

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( par->Sdate(),"C");

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0]->SetAxisRange(0.005,0.03,"X");
   hdat[0]->SetTitle(";M^{2}_{#gamma#gamma} , GeV^{2}/c^{4}"
         ";Entries/0.0002 GeV^{2}/c^{4}");
   if ( date > 2009 ) {
      hdat[0]->GetYaxis()->SetMaxDigits(3);
   }
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
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
   // hmc[1]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP`");
   leg->AddEntry(hmc[0], "MC","L");
   // leg->AddEntry(hmc[1], "MC signal","F");
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

   vector<string> hndat {"S3_M2K"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S3_M2K", "MC_M2KT", "MC_M2KF"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   int nm = hdat[0]->GetNbinsX() / 2;
   double scale = hdat[0]->Integral(nm-10,nm+10) /
      hmc[0]->Integral(nm-10,nm+10);
   // double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( par->Sdate(),"C");

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0]->SetAxisRange(0.15,0.35,"X");
   hdat[0]->SetTitle(
         ";M^{2}_{recoil}(2(#pi^{+}#pi^{-})K^{#pm}), GeV^{2}/c^{4}"
         ";Entries/0.0012 GeV^{2}/c^{4}"
         );
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   if ( date == 2009 ) {
      hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   }
   double ymax = 1.06*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   const double mk   = 0.493677; // 493.677  +/- 0.016 MeV
   double M2l = SQ(mk)-0.025;
   double M2r = SQ(mk)+0.025;
   box->DrawBox(0.15,0.,M2l,ymax);
   box->DrawBox(M2r,0.,0.35,ymax);
   hdat[0]->Draw("E,SAME");

   double dh  = hdat[0]->GetBinWidth(1);
   double bl  = hdat[0]->GetBinLowEdge(1);
   double sum =
      hmc[0]->Integral(int((M2l-bl)/dh)-1,int((M2r-bl)/dh)+1);
   double bg  =
      hmc[2]->Integral(int((M2l-bl)/dh)-1,int((M2r-bl)/dh)+1);
   printf("%i: background in the selection window is %.1f%%\n",
      date, bg/sum*100);

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

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
   auto hPt = [](string nm) {
      return (new TH1D(nm.c_str(),"",60,0.,1.5));
   };

   TTree* effKd = par->GetEff(0); // data
   vector<TH1D*> hdat { hPt("hdat") };
   TCut cutD = par->Ckaon  + TCut("Zk>0") + TCut("fabs(Ck)<0.8");
   effKd->Draw("Ptk>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
      // cout<<h->GetName()<<" Entries= "<<h->GetEntries()<<endl;
   }

   TTree* effKmc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hPt("hmc"), hPt("hmcBG") };
   effKmc->Draw("Ptk>>hmc",cutD,"goff");
   effKmc->Draw("Ptk>>hmcBG",cutD + !(par->Cmcsig),"goff");
   // cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1]->Integral(); // [bin1,bin2]
   double sum = hmc[0]->Integral();
   printf("%i %s background is %.1f%%\n",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i  K^{#plus}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0]->SetTitle(
         ";P_{t}(K), GeV/c"
         ";Entries/0.025 GeV/c"
         );
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[1]->SetLineWidth(1);
   hmc[1]->SetLineColor(kBlue+1);
   hmc[1]->SetFillStyle(3001);
   hmc[1]->SetFillColor(kBlue+1);
   hmc[1]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
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
   auto hPt = [](string nm) {
      return (new TH1D(nm.c_str(),"",60,0.,1.5));
   };

   TTree* effKd = par->GetEff(0); // data
   vector<TH1D*> hdat { hPt("hdat") };
   TCut cutD = par->Ckaon  + TCut("Zk<0") + TCut("fabs(Ck)<0.8");
   effKd->Draw("Ptk>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
   }

   TTree* effKmc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hPt("hmc"), hPt("hmcBG") };
   effKmc->Draw("Ptk>>hmc",cutD,"goff");
   effKmc->Draw("Ptk>>hmcBG",cutD + !(par->Cmcsig),"goff");
   // cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1]->Integral(); // [bin1,bin2]
   double sum = hmc[0]->Integral();
   printf("%i %s background is %.1f%%\n",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i  K^{#minus}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0]->SetTitle(
         ";P_{t}(K), GeV/c"
         ";Entries/0.025 GeV/c"
         );
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[1]->SetLineWidth(1);
   hmc[1]->SetLineColor(kBlue+1);
   hmc[1]->SetFillStyle(3001);
   hmc[1]->SetFillColor(kBlue+1);
   hmc[1]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
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
   auto hC = [](string nm) {
      return (new TH1D(nm.c_str(),"",40,-1.,1.));
   };

   TTree* effKd = par->GetEff(0); // data
   vector<TH1D*> hdat { hC("hdat") };
   TCut cutD = par->Ckaon + TCut("Zk>0") + TCut("0.1<Ptk&&Ptk<1.4");
   effKd->Draw("Ck>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
   }

   TTree* effKmc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hC("hmc"), hC("hmcBG") };
   effKmc->Draw("Ck>>hmc",cutD,"goff");
   effKmc->Draw("Ck>>hmcBG",cutD + !(par->Cmcsig),"goff");
   // cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1]->Integral(); // [bin1,bin2]
   double sum = hmc[0]->Integral();
   printf("%i %s background is %.1f%%\n",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral(5,36) / hmc[0]->Integral(5,36);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i  K^{#plus}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->SetTitle(
         ";cos(#Theta_{K}) "
         ";Entries/0.05 "
                       );
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[1]->SetLineWidth(1);
   hmc[1]->SetLineColor(kBlue+1);
   hmc[1]->SetFillStyle(3001);
   hmc[1]->SetFillColor(kBlue+1);
   hmc[1]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
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
   auto hC = [](string nm) {
      return (new TH1D(nm.c_str(),"",40,-1.,1.));
   };

   TTree* effKd = par->GetEff(0); // data
   vector<TH1D*> hdat { hC("hdat") };
   TCut cutD = par->Ckaon  + TCut("Zk<0") + TCut("0.1<Ptk&&Ptk<1.4");
   effKd->Draw("Ck>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
   }

   TTree* effKmc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hC("hmc"), hC("hmcBG") };
   effKmc->Draw("Ck>>hmc",cutD,"goff");
   effKmc->Draw("Ck>>hmcBG",cutD + !(par->Cmcsig),"goff");
   // cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1]->Integral(); // [bin1,bin2]
   double sum = hmc[0]->Integral();
   printf("%i %s background is %.1f%%\n",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral(5,36) / hmc[0]->Integral(5,36);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i  K^{#minus}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->SetTitle(
         ";cos(#Theta_{K}) "
         ";Entries/0.05 "
         );
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[1]->SetLineWidth(1);
   hmc[1]->SetLineColor(kBlue+1);
   hmc[1]->SetFillStyle(3001);
   hmc[1]->SetFillColor(kBlue+1);
   hmc[1]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
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

   vector<string> hndat {"S3_dPK"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S3_dPK", "MC_dPK_T", "MC_dPK_F"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   // double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( par->Sdate(), "C" );

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   hdat[0]->SetAxisRange(-0.15,0.15,"X");
   hdat[0]->SetMinimum(1.);
   hdat[0]->SetTitle(
         ";#deltaP(K), GeV/c"
         ";Entries/0.002 GeV/c"
         );
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   // double ymax=1.06*hdat[0]->GetMaximum(); // Lin
   double ymax = 1.9*hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(-0.15,0.,-0.12,ymax);
   box->DrawBox(0.08,0.,0.15,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "No match","F");
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

   vector<string> hndat {"S3_dThK"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S3_dThK", "MC_dThK_T", "MC_dThK_F"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   // double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( par->Sdate(), "C" );

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   // hdat[0]->SetAxisRange(-0.15,0.15,"X");
   hdat[0]->SetMinimum(1.);
   hdat[0]->SetTitle(
         ";#delta#Theta(K), deg"
         ";Entries/0.2 deg"
         );
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   // double ymax=1.06*hdat[0]->GetMaximum(); // Lin
   double ymax=1.9*hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(10.,0.,20.,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "No match","F");
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
   Params* par = new Params(date);

   vector<string> hndat {"S4_MppKKb"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S4_MppKKb", "MC_MppKK_T", "MC_MppKK_F"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   // double scale = hdat[0]->Integral() / hmc[0]->Integral();
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

   TLegend* leg = nullptr;
   if ( fitgs ) {
      leg = new TLegend(0.12,0.69,0.405,0.895);
   } else {
      leg = new TLegend(0.61,0.69,0.895,0.895);
   }
   leg->SetHeader( par->Sdate(), "C" );

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0]->SetAxisRange(3.06,3.13,"X");
   hdat[0]->SetTitle(
         ";M^{ inv}_{#pi^{#plus}#pi^{#minus}K^{#plus}K^{#minus}},"
         " GeV/c^{2}"
         ";Entries/0.0009 GeV/c^{2}"
         );
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   double ymax = 1.06 * hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);

   if ( fitgs ) {
      gStyle->SetStatY(0.895);
      TF1* gs = (TF1*)gROOT->GetFunction("gaus");
      gs->SetLineWidth(2);
      gs->SetLineStyle(5);
      gs->SetLineColor(kMagenta);
      hdat[0]->Fit(gs,"Q","",3.087,3.105);
   } else {
      hdat[0]->Draw("E");
   }

   // fabs(MppKK-3.096)<0.009
   box->DrawBox(3.06,0.,3.087,ymax);
   box->DrawBox(3.105,0.,3.13,ymax);
   // [3.08, 3.11]
   // box->DrawBox(3.06,0.,3.08,ymax);
   // box->DrawBox(3.11,0.,3.13,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

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
   Params* par = new Params(date);

   vector<string> hndat {"S5_M2Pi"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S5_M2Pi", "MC_M2piT", "MC_M2piF"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   // double scale = hdat[0]->Integral() / hmc[0]->Integral();
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

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( par->Sdate(), "C" );

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   // hdat[0]->SetAxisRange(3.06,3.13,"X");
   hdat[0]->SetTitle(
         ";M^{2}_{recoil}"
         "(#pi^{#pm}#pi^{#plus}#pi^{#minus}K^{#plus}K^{#minus}),"
         " GeV^{2}/c^{4}"
         ";Entries/0.0003 GeV^{2}/c^{4}"
         );
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   double ymax = 1.07*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   const double mpi  = 0.13957;  // 139.57018 +/- 0.00035 MeV
   double M2l = SQ(mpi)-0.01;
   double M2r = SQ(mpi)+0.01;
   box->DrawBox(-0.01,0.,M2l,ymax);
   box->DrawBox(M2r,0.,0.05,ymax);
   hdat[0]->Draw("E,SAME");

   double dh  = hdat[0]->GetBinWidth(1);
   double bl  = hdat[0]->GetBinLowEdge(1);
   double sum =
      hmc[0]->Integral(int((M2l-bl)/dh)-1,int((M2r-bl)/dh)+1);
   double bg  =
      hmc[2]->Integral(int((M2l-bl)/dh)-1,int((M2r-bl)/dh)+1);
   printf("%s %i: background in the selection window is %.1f%%\n",
      __func__, date, bg/sum*100);

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

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
   auto hPt = [](string nm) {
      return (new TH1D(nm.c_str(),"",40,0.,0.4));
   };

   TTree* effPid = par->GetEff(0); // data
   vector<TH1D*> hdat { hPt("hdat") };
   TCut cutD = par->Cpion  + TCut("Zpi>0") + TCut("fabs(Cpi)<0.8");
   effPid->Draw("Ptpi>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
   }

   TTree* effPimc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hPt("hmc"), hPt("hmcBG") };
   effPimc->Draw("Ptpi>>hmc",cutD,"goff");
   effPimc->Draw("Ptpi>>hmcBG",cutD + !(par->Cmcsig),"goff");
   // cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1]->Integral(); // [bin1,bin2]
   double sum = hmc[0]->Integral();
   printf("%i %s background is %.1f%%\n",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   TLegend* leg = new TLegend(0.11,0.69,0.38,0.89);
   leg->SetHeader( Form("%i  #pi^{#kern[0.25]{#plus}}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0]->SetTitle(
         ";P_{t}(#pi), GeV/c"
         ";Entries/0.01 GeV/c"
         );
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[1]->SetLineWidth(1);
   hmc[1]->SetLineColor(kBlue+1);
   hmc[1]->SetFillStyle(3001);
   hmc[1]->SetFillColor(kBlue+1);
   hmc[1]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
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
   auto hPt = [](string nm) {
      return (new TH1D(nm.c_str(),"",40,0.,0.4));
   };

   TTree* effPid = par->GetEff(0); // data
   vector<TH1D*> hdat { hPt("hdat") };
   TCut cutD = par->Cpion  + TCut("Zpi<0") + TCut("fabs(Cpi)<0.8");
   effPid->Draw("Ptpi>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
   }

   TTree* effPimc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hPt("hmc"), hPt("hmcBG") };
   effPimc->Draw("Ptpi>>hmc",cutD,"goff");
   effPimc->Draw("Ptpi>>hmcBG",cutD + !(par->Cmcsig),"goff");
   // cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1]->Integral(); // [bin1,bin2]
   double sum = hmc[0]->Integral();
   printf("%i %s background is %.1f%%\n",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   TLegend* leg = new TLegend(0.11,0.69,0.38,0.89);
   leg->SetHeader( Form("%i  #pi^{#kern[0.25]{#minus}}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0]->SetTitle(
         ";P_{t}(#pi), GeV/c"
         ";Entries/0.01 GeV/c"
                       );
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[1]->SetLineWidth(1);
   hmc[1]->SetLineColor(kBlue+1);
   hmc[1]->SetFillStyle(3001);
   hmc[1]->SetFillColor(kBlue+1);
   hmc[1]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
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
   auto hC = [](string nm) {
      return (new TH1D(nm.c_str(),"",40,-1.,1.));
   };

   TTree* effPid = par->GetEff(0); // data
   vector<TH1D*> hdat { hC("hdat") };
   TCut cutD = par->Cpion  + TCut("Zpi>0")
      + TCut("0.05<Ptpi&&Ptpi<0.4");
   effPid->Draw("Cpi>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
   }

   TTree* effPimc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hC("hmc"), hC("hmcBG") };
   effPimc->Draw("Cpi>>hmc",cutD,"goff");
   effPimc->Draw("Cpi>>hmcBG",cutD + !(par->Cmcsig),"goff");
   // cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1]->Integral(); // [bin1,bin2]
   double sum = hmc[0]->Integral();
   printf("%i %s background is %.1f%%\n",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral(5,36) / hmc[0]->Integral(5,36);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   TLegend* leg = new TLegend(0.11,0.69,0.38,0.89);
   leg->SetHeader( Form("%i  #pi^{#kern[0.25]{#plus}}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->SetTitle(
         ";cos(#Theta_{#pi}) "
         ";Entries/0.05 "
                       );
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[1]->SetLineWidth(1);
   hmc[1]->SetLineColor(kBlue+1);
   hmc[1]->SetFillStyle(3001);
   hmc[1]->SetFillColor(kBlue+1);
   hmc[1]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
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
   auto hC = [](string nm) {
      return (new TH1D(nm.c_str(),"",40,-1.,1.));
   };

   TTree* effPid = par->GetEff(0); // data, pi
   vector<TH1D*> hdat { hC("hdat") };
   TCut cutD = par->Cpion  + TCut("Zpi<0")
      + TCut("0.05<Ptpi&&Ptpi<0.4");
   effPid->Draw("Cpi>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
   }

   TTree* effPimc = par->GetEff(1); // MC, pi
   vector<TH1D*> hmc { hC("hmc"), hC("hmcBG") };
   effPimc->Draw("Cpi>>hmc",cutD,"goff");
   effPimc->Draw("Cpi>>hmcBG",cutD + !(par->Cmcsig),"goff");
   // cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1]->Integral(); // [bin1,bin2]
   double sum = hmc[0]->Integral();
   printf("%i %s background is %.1f%%\n",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0]->Integral(5,36) / hmc[0]->Integral(5,36);
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   TLegend* leg = new TLegend(0.11,0.69,0.38,0.89);
   leg->SetHeader( Form("%i  #pi^{#kern[0.25]{#minus}}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0]->SetMaximum(ymax);
   hdat[0]->SetTitle(
         ";cos(#Theta_{#pi}) "
         ";Entries/0.05 "
         );
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.25);
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->Draw("E");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[1]->SetLineWidth(1);
   hmc[1]->SetLineColor(kBlue+1);
   hmc[1]->SetFillStyle(3001);
   hmc[1]->SetFillColor(kBlue+1);
   hmc[1]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
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
   Params* par = new Params(date);

   vector<string> hndat {"S5_dPpi"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S5_dPpi", "MC_dPpi_T", "MC_dPpi_F"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   // double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( par->Sdate(), "C" );

   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   hdat[0]->SetAxisRange(-0.15,0.15,"X");
   hdat[0]->SetMinimum(1.);
   hdat[0]->SetTitle(
         ";#deltaP(#pi), GeV/c"
         ";Entries/0.002 GeV/c"
         );
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   // double ymax=1.06*hdat[0]->GetMaximum(); // Lin
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(-0.15,0.,-0.12,ymax);
   box->DrawBox(0.08,0.,0.15,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "No match","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_dPpi_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

//--------------------------------------------------------------------
void plot_dThPi(int date, int Cx = 800, int Cy = 800) {
//--------------------------------------------------------------------
   Params* par = new Params(date);

   vector<string> hndat {"S5_dThPi"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S5_dThPi", "MC_dThPi_T", "MC_dThPi_F"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0]->Integral() / hmc[0]->Integral();
   // double scale = hdat[0]->GetMaximum() / hmc[0]->GetMaximum();
   for ( auto& h : hmc ) {
      h->Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( par->Sdate(), "C" );


   // Draw:
   TCanvas* c1 = new TCanvas(par->Sdate(),par->Sdate(),0,0,Cx,Cy);
   c1->cd(0);
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   // hdat[0]->SetAxisRange(-0.15,0.15,"X");
   hdat[0]->SetMinimum(1.);
   hdat[0]->SetTitle(
         ";#delta#Theta(#pi), deg"
         ";Entries/0.2 deg"
         );
   hdat[0]->SetLineWidth(2);
   hdat[0]->SetMarkerStyle(20);
   hdat[0]->SetLineColor(kBlack);
   hdat[0]->GetYaxis()->SetMaxDigits(3);
   hdat[0]->GetYaxis()->SetTitleOffset(1.2);
   // double ymax=1.06*hdat[0]->GetMaximum(); // Lin
   double ymax = 1.9 * hdat[0]->GetMaximum(); // Log
   hdat[0]->SetMaximum(ymax);
   hdat[0]->Draw("E");

   box->DrawBox(15.,0.,20.,ymax);
   hdat[0]->Draw("E,SAME");

   hmc[0]->SetLineWidth(1);
   hmc[0]->SetLineColor(kRed+2);
   hmc[0]->Draw("HIST,SAME");

   hmc[2]->SetLineWidth(1);
   hmc[2]->SetLineColor(kBlue+1);
   hmc[2]->SetFillStyle(3001);
   hmc[2]->SetFillColor(kBlue+1);
   hmc[2]->Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "No match","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = "trkeff_dThPi_" + to_string(date) + ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Plot ...
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
   // gStyle->SetStatFont(62); // ?

   size_t Cx = 880, Cy = 800; // canvas sizes

   // for ( auto date : {2009, 2012} ) {
   for ( auto date : {2009, 2012, 2021} ) {
      // -- pi0 rejection, fig.28
      // plot_pi0(date,Cx,Cy);

      // ---- K+ K- ----
      // fig 29
      // plot_MmisK(date,Cx,Cy);

      // fig 30,31
      // plot_PtKp(date,Cx,Cy);
      // plot_PtKm(date,Cx,Cy);
      // plot_CKp(date,Cx,Cy);
      // plot_CKm(date,Cx,Cy);

      // fig 32
      // plot_dThK(date,Cx,Cy);
      // plot_dPK(date,Cx,Cy);

      // ---- pi+ pi- ----
      // fig 33
      // plot_MinvJ(date,Cx,Cy);
      // plot_MinvJ(date,Cx,Cy,true); // fit J/Psi peak position

      // fig 34
      // plot_MmisP(date,Cx,Cy);

      // fig 35,36
      // plot_PtPip(date,Cx,Cy);
      // plot_PtPim(date,Cx,Cy);
      // plot_CPip(date,Cx,Cy);
      // plot_CPim(date,Cx,Cy);

      // fig 37
      // plot_dThPi(date,Cx,Cy);
      // plot_dPpi(date,Cx,Cy);
   }

   // other: do not use
   // plot_Weights(2012);
   // plot_Weights(2009);
}
