// trk_eff_sel.cc - Pictures for presentation/memo
// Study of the track reconstruction efficiency for pi and K:

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

//-------------------------------------------------------------------------
void get_hst( Params* p, int mc,
      const vector<string> hns, vector<TH1D*>& hst ) {
//-------------------------------------------------------------------------
// get set of histograms from file "fname" according of "hns"

   int Nh = hns.size();
   hst.clear();
   hst.resize(Nh,nullptr);

   TFile* froot = p->OpenFile(mc);

   for ( int i = 0; i < Nh; ++i ) {
      hst[i] = (TH1D*)gROOT->FindObject(hns[i].c_str());
      if ( !hst[i] ) {
         cout << " can not find histo:" << hns[i] << endl;
         exit(1);
      }
      hst[i]->Sumw2(true);
      SetHstFace(hst[i]);

      // rename:
      string name = string((mc==1) ? "MC_" : "Dat") +
                    string((p->date==2009) ? "09_"  : "12_") +
                    string( hst[i]->GetName() );
      hst[i]->SetName( name.c_str() );
//       hst[i]->SetTitle("");
   }
}

//-------------------------------------------------------------------------
void plot_pi0(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date);

   vector<string> hndat {"S1b_Mg2"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S1b_Mg2", "MC_M2ggT", "MC_M2ggF"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
//    double scale = hdat[0] -> Integral(26,150) / hmc[0] -> Integral(26,150);
   double scale = hdat[0] -> Integral(60,109) / hmc[0] -> Integral(60,109);
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box -> SetFillStyle(3001);
   box -> SetFillColor(kRed-10);
   box -> SetLineColor(kRed-9);
   box -> SetLineWidth(1);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg -> SetHeader( Form("%i",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,600);
   c1 -> cd(0);
   gPad -> SetGrid();

   // Data
   hdat[0] -> SetAxisRange(0.005,0.03,"X");
   hdat[0] -> SetTitle(";M^{2}_{#gamma#gamma} , GeV^{2}/c^{4}"
                        ";Entries/0.0002 GeV^{2}/c^{4}");
   hdat[0] -> GetYaxis() -> SetTitleOffset(1.);
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> Draw("E");

   double ymax=1.05*hdat[0]->GetMaximum();
   box -> DrawBox(0.012,.0,0.022,ymax);
   hdat[0] -> Draw("E,SAME");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[1] -> SetLineWidth(1);
   hmc[1] -> SetLineColor(kGreen+1);
   hmc[1] -> SetFillStyle(3001);
   hmc[1] -> SetFillColor(kGreen+1);
//    hmc[1] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP`");
   leg->AddEntry(hmc[0], "MC","L");
//    leg->AddEntry(hmc[1], "MC signal","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_pi0_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_MmisK(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date);

   vector<string> hndat {"S3_M2K"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S3_M2K", "MC_M2KT", "MC_M2KF"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
//    double scale = hdat[0] -> Integral() / hmc[0] -> Integral();
   double scale = hdat[0] -> GetMaximum() / hmc[0] -> GetMaximum();
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);
   const double mk   = 0.493677; // 493.677  +/- 0.016 MeV
   double M2l = SQ(mk)-0.025;
   double M2r = SQ(mk)+0.025;
   double dh = 0.0012; // bin width
   double bg = hmc[2] -> Integral(int((M2l-0.125)/dh),int((M2r-0.125)/dh));
   double sum = hmc[0] -> Integral(int((M2l-0.125)/dh),int((M2r-0.125)/dh));
//    cout << " bg=" << bg << " sum= " << sum << " R= " << bg/sum << endl;
   printf(" %i background is %.1f%%",date,bg/sum*100);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( (string("#bf{")+to_string(date)+"}").c_str(),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0] -> SetAxisRange(0.15,0.35,"X");
   hdat[0] -> SetTitle(
         ";M^{2}_{recoil}(2(#pi^{+}#pi^{-})K^{#pm}), GeV^{2}/c^{4}"
         ";Entries/0.0012 GeV^{2}/c^{4}"
                       );
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   if ( date == 2012 ) {
      hdat[0] -> GetYaxis() -> SetTitleOffset(1.0);
   }
   hdat[0] -> Draw("E");

   double ymax=1.06*hdat[0]->GetMaximum();
   box->DrawBox(0.15,0.,M2l,ymax);
   box->DrawBox(M2r,0.,0.35,ymax);
   hdat[0] -> Draw("E,SAME");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[2] -> SetLineWidth(1);
   hmc[2] -> SetLineColor(kBlue+1);
   hmc[2] -> SetFillStyle(3001);
   hmc[2] -> SetFillColor(kBlue+1);
   hmc[2] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_MmK_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_dPK(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date);

   vector<string> hndat {"S3_dPK"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S3_dPK", "MC_dPK_T", "MC_dPK_F"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0] -> Integral() / hmc[0] -> Integral();
//    double scale = hdat[0] -> GetMaximum() / hmc[0] -> GetMaximum();
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( (string("#bf{")+to_string(date)+"}").c_str(),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd(0);
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   hdat[0] -> SetAxisRange(-0.15,0.15,"X");
   hdat[0] -> SetMinimum(1.);
   hdat[0] -> SetTitle(
         ";#deltaP(K), GeV/c"
         ";Entries/0.002 GeV/c"
                       );
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   hdat[0] -> Draw("E");

//    double ymax=1.06*hdat[0]->GetMaximum(); // Lin
   double ymax=1.8*hdat[0]->GetMaximum(); // Log
   box->DrawBox(-0.15,0.,-0.12,ymax);
   box->DrawBox(0.08,0.,0.15,ymax);
   hdat[0] -> Draw("E,SAME");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[2] -> SetLineWidth(1);
   hmc[2] -> SetLineColor(kBlue+1);
   hmc[2] -> SetFillStyle(3001);
   hmc[2] -> SetFillColor(kBlue+1);
   hmc[2] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "No match","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_dPK_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_dThK(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date);

   vector<string> hndat {"S3_dThK"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S3_dThK", "MC_dThK_T", "MC_dThK_F"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0] -> Integral() / hmc[0] -> Integral();
//    double scale = hdat[0] -> GetMaximum() / hmc[0] -> GetMaximum();
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( (string("#bf{")+to_string(date)+"}").c_str(),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd(0);
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
//    hdat[0] -> SetAxisRange(-0.15,0.15,"X");
   hdat[0] -> SetMinimum(1.);
   hdat[0] -> SetTitle(
         ";#delta#Theta(K), deg"
         ";Entries/0.2 deg"
                       );
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   hdat[0] -> Draw("E");

//    double ymax=1.06*hdat[0]->GetMaximum(); // Lin
   double ymax=1.8*hdat[0]->GetMaximum(); // Log
   box->DrawBox(10.,0.,20.,ymax);
   hdat[0] -> Draw("E,SAME");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[2] -> SetLineWidth(1);
   hmc[2] -> SetLineColor(kBlue+1);
   hmc[2] -> SetFillStyle(3001);
   hmc[2] -> SetFillColor(kBlue+1);
   hmc[2] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "No match","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_dThK_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_PtKp(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date,1); // Kaons
   auto hPt = [](string nm) {return (new TH1D(nm.c_str(),"",60,0.,1.5));};

   TTree* effKd = par -> GetEff(0); // data
   vector<TH1D*> hdat { hPt("hdat") };
   TCut cutD = par -> Ckaon  + TCut("Zk>0") + TCut("fabs(Ck)<0.8");
   effKd->Draw("Ptk>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
//       cout << h->GetName() << " Entries= " << h->GetEntries() << endl;
   }

   TTree* effKmc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hPt("hmc"), hPt("hmcBG") };
   effKmc->Draw("Ptk>>hmc",cutD,"goff");
   effKmc->Draw("Ptk>>hmcBG",cutD + !(par -> Cmcsig),"goff");
//    cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1] -> Integral(); // [bin1,bin2]
   double sum = hmc[0] -> Integral();
   printf(" %i %s background is %.1f%%",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0] -> Integral() / hmc[0] -> Integral();
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i  K^{+}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,600);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0] -> SetTitle(
         ";P_{t}(K), GeV/c"
         ";Entries/0.025 GeV/c"
                       );
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   hdat[0] -> GetYaxis() -> SetTitleOffset(1.);
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> Draw("E");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[1] -> SetLineWidth(1);
   hmc[1] -> SetLineColor(kBlue+1);
   hmc[1] -> SetFillStyle(3001);
   hmc[1] -> SetFillColor(kBlue+1);
   hmc[1] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_PtKp_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_PtKm(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date, 1); // Kaons
   auto hPt = [](string nm) {return (new TH1D(nm.c_str(),"",60,0.,1.5));};

   TTree* effKd = par -> GetEff(0); // data
   vector<TH1D*> hdat { hPt("hdat") };
   TCut cutD = par -> Ckaon  + TCut("Zk<0") + TCut("fabs(Ck)<0.8");
   effKd->Draw("Ptk>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
//       cout << h->GetName() << " Entries= " << h->GetEntries() << endl;
   }

   TTree* effKmc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hPt("hmc"), hPt("hmcBG") };
   effKmc->Draw("Ptk>>hmc",cutD,"goff");
   effKmc->Draw("Ptk>>hmcBG",cutD + !(par -> Cmcsig),"goff");
//    cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1] -> Integral(); // [bin1,bin2]
   double sum = hmc[0] -> Integral();
   printf(" %i %s background is %.1f%%",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0] -> Integral() / hmc[0] -> Integral();
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i  K^{-}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,600);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0] -> SetTitle(
         ";P_{t}(K), GeV/c"
         ";Entries/0.025 GeV/c"
                       );
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   hdat[0] -> GetYaxis() -> SetTitleOffset(1.);
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> Draw("E");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[1] -> SetLineWidth(1);
   hmc[1] -> SetLineColor(kBlue+1);
   hmc[1] -> SetFillStyle(3001);
   hmc[1] -> SetFillColor(kBlue+1);
   hmc[1] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_PtKm_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_CKp(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date,1); // Kaons
   auto hC = [](string nm) {return (new TH1D(nm.c_str(),"",40,-1.,1.));};

   TTree* effKd = par -> GetEff(0); // data
   vector<TH1D*> hdat { hC("hdat") };
   TCut cutD = par -> Ckaon  + TCut("Zk>0") + TCut("0.1<Ptk&&Ptk<1.4");
   effKd->Draw("Ck>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
//       cout << h->GetName() << " Entries= " << h->GetEntries() << endl;
   }

   TTree* effKmc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hC("hmc"), hC("hmcBG") };
   effKmc->Draw("Ck>>hmc",cutD,"goff");
   effKmc->Draw("Ck>>hmcBG",cutD + !(par -> Cmcsig),"goff");
//    cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1] -> Integral(); // [bin1,bin2]
   double sum = hmc[0] -> Integral();
   printf(" %i %s background is %.1f%%",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0] -> Integral(5,36) / hmc[0] -> Integral(5,36);
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i  K^{+}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,600);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0] -> SetMaximum(ymax);
   hdat[0] -> SetTitle(
         ";cos(#Theta_{K}) "
         ";Entries/0.05 "
                       );
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   hdat[0] -> GetYaxis() -> SetTitleOffset(1.);
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> Draw("E");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[1] -> SetLineWidth(1);
   hmc[1] -> SetLineColor(kBlue+1);
   hmc[1] -> SetFillStyle(3001);
   hmc[1] -> SetFillColor(kBlue+1);
   hmc[1] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_CKp_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_CKm(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date,1);
   auto hC = [](string nm) {return (new TH1D(nm.c_str(),"",40,-1.,1.));};

   TTree* effKd = par -> GetEff(0); // data
   vector<TH1D*> hdat { hC("hdat") };
   TCut cutD = par -> Ckaon  + TCut("Zk<0") + TCut("0.1<Ptk&&Ptk<1.4");
   effKd->Draw("Ck>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
//       cout << h->GetName() << " Entries= " << h->GetEntries() << endl;
   }

   TTree* effKmc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hC("hmc"), hC("hmcBG") };
   effKmc->Draw("Ck>>hmc",cutD,"goff");
   effKmc->Draw("Ck>>hmcBG",cutD + !(par -> Cmcsig),"goff");
//    cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1] -> Integral(); // [bin1,bin2]
   double sum = hmc[0] -> Integral();
   printf(" %i %s background is %.1f%%",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0] -> Integral(5,36) / hmc[0] -> Integral(5,36);
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( Form("%i  K^{-}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,600);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0] -> SetMaximum(ymax);
   hdat[0] -> SetTitle(
         ";cos(#Theta_{K}) "
         ";Entries/0.05 "
                       );
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   hdat[0] -> GetYaxis() -> SetTitleOffset(1.);
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> Draw("E");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[1] -> SetLineWidth(1);
   hmc[1] -> SetLineColor(kBlue+1);
   hmc[1] -> SetFillStyle(3001);
   hmc[1] -> SetFillColor(kBlue+1);
   hmc[1] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_CKm_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_MinvJ(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date);

   vector<string> hndat {"S4_MppKKb"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S4_MppKKb", "MC_MppKK_T", "MC_MppKK_F"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
//    double scale = hdat[0] -> Integral() / hmc[0] -> Integral();
   double scale = hdat[0] -> GetMaximum() / hmc[0] -> GetMaximum();
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   TLegend* leg = new TLegend(0.61,0.69,0.895,0.895);
   leg->SetHeader( (string("#bf{")+to_string(date)+"}").c_str(),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0] -> SetAxisRange(3.06,3.13,"X");
   hdat[0] -> SetTitle(
         ";M^{ inv}_{#pi^{+}#pi^{-}K^{+}K^{-}}, GeV/c^{2}"
         ";Entries/0.0009 GeV/c^{2}"
                       );
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   if ( date == 2012 ) {
      hdat[0] -> GetYaxis() -> SetTitleOffset(1.0);
   }
   hdat[0] -> Draw("E");

   double ymax=1.06*hdat[0]->GetMaximum();
   box->DrawBox(3.06,0.,3.087,ymax);
   box->DrawBox(3.105,0.,3.13,ymax);
   hdat[0] -> Draw("E,SAME");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[2] -> SetLineWidth(1);
   hmc[2] -> SetLineColor(kBlue+1);
   hmc[2] -> SetFillStyle(3001);
   hmc[2] -> SetFillColor(kBlue+1);
   hmc[2] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_MiJ_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_MmisP(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date);

   vector<string> hndat {"S5_M2Pi"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S5_M2Pi", "MC_M2piT", "MC_M2piF"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
//    double scale = hdat[0] -> Integral() / hmc[0] -> Integral();
   double scale = hdat[0] -> GetMaximum() / hmc[0] -> GetMaximum();
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);
   const double mpi  = 0.13957;  // 139.57018 +/- 0.00035 MeV
   double M2l = SQ(mpi)-0.01;
   double M2r = SQ(mpi)+0.01;
   double dh = 0.0003; // bin width
   double bg = hmc[2] -> Integral(int((M2l+0.01)/dh),int((M2r+0.01)/dh));
   double sum = hmc[0] -> Integral(int((M2l+0.01)/dh),int((M2r+0.01)/dh));
//    cout << " bg=" << bg << " sum= " << sum << " R= " << bg/sum << endl;
   printf(" %i background is %.1f%%",date,bg/sum*100);

//    TLegend* leg = new TLegend(0.61,0.69,0.895,0.895);
   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( (string("#bf{")+to_string(date)+"}").c_str(),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd(0);
   gPad->SetGrid();

   // Data
//    hdat[0] -> SetAxisRange(3.06,3.13,"X");
   hdat[0] -> SetTitle(
         ";M^{2}_{recoil}(#pi^{#pm}#pi^{+}#pi^{-}K^{+}K^{-}), GeV^{2}/c^{4}"
         ";Entries/0.0003 GeV^{2}/c^{4}"
                       );
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   if ( date == 2012 ) {
      hdat[0] -> GetYaxis() -> SetTitleOffset(1.0);
   }
   hdat[0] -> Draw("E");

   double ymax=1.06*hdat[0]->GetMaximum();
   box->DrawBox(-0.01,0.,M2l,ymax);
   box->DrawBox(M2r,0.,0.05,ymax);
   hdat[0] -> Draw("E,SAME");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[2] -> SetLineWidth(1);
   hmc[2] -> SetLineColor(kBlue+1);
   hmc[2] -> SetFillStyle(3001);
   hmc[2] -> SetFillColor(kBlue+1);
   hmc[2] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_MmP_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_dPpi(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date);

   vector<string> hndat {"S5_dPpi"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S5_dPpi", "MC_dPpi_T", "MC_dPpi_F"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0] -> Integral() / hmc[0] -> Integral();
//    double scale = hdat[0] -> GetMaximum() / hmc[0] -> GetMaximum();
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( (string("#bf{")+to_string(date)+"}").c_str(),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd(0);
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
   hdat[0] -> SetAxisRange(-0.15,0.15,"X");
   hdat[0] -> SetMinimum(1.);
   hdat[0] -> SetTitle(
         ";#deltaP(#pi), GeV/c"
         ";Entries/0.002 GeV/c"
                       );
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   hdat[0] -> Draw("E");

//    double ymax=1.06*hdat[0]->GetMaximum(); // Lin
   double ymax=1.8*hdat[0]->GetMaximum(); // Log
   box->DrawBox(-0.15,0.,-0.12,ymax);
   box->DrawBox(0.08,0.,0.15,ymax);
   hdat[0] -> Draw("E,SAME");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[2] -> SetLineWidth(1);
   hmc[2] -> SetLineColor(kBlue+1);
   hmc[2] -> SetFillStyle(3001);
   hmc[2] -> SetFillColor(kBlue+1);
   hmc[2] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "No match","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_dPpi_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_dThPi(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date);

   vector<string> hndat {"S5_dThPi"};
   vector<TH1D*> hdat;
   get_hst(par, 0, hndat, hdat);

   vector<string> hnmc {"S5_dThPi", "MC_dThPi_T", "MC_dThPi_F"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   // normalization on DATA
   double scale = hdat[0] -> Integral() / hmc[0] -> Integral();
//    double scale = hdat[0] -> GetMaximum() / hmc[0] -> GetMaximum();
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(1);

   TLegend* leg = new TLegend(0.59,0.69,0.89,0.89);
   leg->SetHeader( (string("#bf{")+to_string(date)+"}").c_str(),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd(0);
   gPad->SetGrid();
   gPad->SetLogy(true);

   // Data
//    hdat[0] -> SetAxisRange(-0.15,0.15,"X");
   hdat[0] -> SetMinimum(1.);
   hdat[0] -> SetTitle(
         ";#delta#Theta(#pi), deg"
         ";Entries/0.2 deg"
                       );
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   hdat[0] -> Draw("E");

//    double ymax=1.06*hdat[0]->GetMaximum(); // Lin
   double ymax=1.8*hdat[0]->GetMaximum(); // Log
   box->DrawBox(15.,0.,20.,ymax);
   hdat[0] -> Draw("E,SAME");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[2] -> SetLineWidth(1);
   hmc[2] -> SetLineColor(kBlue+1);
   hmc[2] -> SetFillStyle(3001);
   hmc[2] -> SetFillColor(kBlue+1);
   hmc[2] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[2], "MC background","F");
   leg->AddEntry(box, "No match","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_dThPi_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_PtPip(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date,2); // Pions
   auto hPt = [](string nm) {return (new TH1D(nm.c_str(),"",40,0.,0.4));};

   TTree* effPid = par -> GetEff(0); // data
   vector<TH1D*> hdat { hPt("hdat") };
   TCut cutD = par -> Cpion  + TCut("Zpi>0") + TCut("fabs(Cpi)<0.8");
   effPid->Draw("Ptpi>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
//       cout << h->GetName() << " Entries= " << h->GetEntries() << endl;
   }

   TTree* effPimc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hPt("hmc"), hPt("hmcBG") };
   effPimc->Draw("Ptpi>>hmc",cutD,"goff");
   effPimc->Draw("Ptpi>>hmcBG",cutD + !(par -> Cmcsig),"goff");
//    cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1] -> Integral(); // [bin1,bin2]
   double sum = hmc[0] -> Integral();
   printf(" %i %s background is %.1f%%",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0] -> Integral() / hmc[0] -> Integral();
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   TLegend* leg = new TLegend(0.11,0.69,0.38,0.89);
   leg->SetHeader( Form("%i  #pi^{+}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,600);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0] -> SetTitle(
         ";P_{t}(#pi), GeV/c"
         ";Entries/0.01 GeV/c"
                       );
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   hdat[0] -> GetYaxis() -> SetTitleOffset(1.);
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> Draw("E");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[1] -> SetLineWidth(1);
   hmc[1] -> SetLineColor(kBlue+1);
   hmc[1] -> SetFillStyle(3001);
   hmc[1] -> SetFillColor(kBlue+1);
   hmc[1] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_PtPip_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_PtPim(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date,2); // pi
   auto hPt = [](string nm) {return (new TH1D(nm.c_str(),"",40,0.,0.4));};

   TTree* effPid = par -> GetEff(0); // data
   vector<TH1D*> hdat { hPt("hdat") };
   TCut cutD = par -> Cpion  + TCut("Zpi<0") + TCut("fabs(Cpi)<0.8");
   effPid->Draw("Ptpi>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
//       cout << h->GetName() << " Entries= " << h->GetEntries() << endl;
   }

   TTree* effPimc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hPt("hmc"), hPt("hmcBG") };
   effPimc->Draw("Ptpi>>hmc",cutD,"goff");
   effPimc->Draw("Ptpi>>hmcBG",cutD + !(par -> Cmcsig),"goff");
//    cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1] -> Integral(); // [bin1,bin2]
   double sum = hmc[0] -> Integral();
   printf(" %i %s background is %.1f%%",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0] -> Integral() / hmc[0] -> Integral();
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   TLegend* leg = new TLegend(0.11,0.69,0.38,0.89);
   leg->SetHeader( Form("%i  #pi^{-}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,600);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   hdat[0] -> SetTitle(
         ";P_{t}(#pi), GeV/c"
         ";Entries/0.01 GeV/c"
                       );
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   hdat[0] -> GetYaxis() -> SetTitleOffset(1.);
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> Draw("E");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[1] -> SetLineWidth(1);
   hmc[1] -> SetLineColor(kBlue+1);
   hmc[1] -> SetFillStyle(3001);
   hmc[1] -> SetFillColor(kBlue+1);
   hmc[1] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_PtPim_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_CPip(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date,2);
   auto hC = [](string nm) {return (new TH1D(nm.c_str(),"",40,-1.,1.));};

   TTree* effPid = par -> GetEff(0); // data
   vector<TH1D*> hdat { hC("hdat") };
   TCut cutD = par -> Cpion  + TCut("Zpi>0") + TCut("0.05<Ptpi&&Ptpi<0.4");
   effPid->Draw("Cpi>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
//       cout << h->GetName() << " Entries= " << h->GetEntries() << endl;
   }

   TTree* effPimc = par->GetEff(1); // MC
   vector<TH1D*> hmc { hC("hmc"), hC("hmcBG") };
   effPimc->Draw("Cpi>>hmc",cutD,"goff");
   effPimc->Draw("Cpi>>hmcBG",cutD + !(par -> Cmcsig),"goff");
//    cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1] -> Integral(); // [bin1,bin2]
   double sum = hmc[0] -> Integral();
   printf(" %i %s background is %.1f%%",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0] -> Integral(5,36) / hmc[0] -> Integral(5,36);
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   TLegend* leg = new TLegend(0.11,0.69,0.38,0.89);
   leg->SetHeader( Form("%i  #pi^{+}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,600);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0] -> SetMaximum(ymax);
   hdat[0] -> SetTitle(
         ";cos(#Theta_{#pi}) "
         ";Entries/0.05 "
                       );
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   hdat[0] -> GetYaxis() -> SetTitleOffset(1.);
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> Draw("E");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[1] -> SetLineWidth(1);
   hmc[1] -> SetLineColor(kBlue+1);
   hmc[1] -> SetFillStyle(3001);
   hmc[1] -> SetFillColor(kBlue+1);
   hmc[1] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_CPip_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_CPim(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date,2);
   auto hC = [](string nm) {return (new TH1D(nm.c_str(),"",40,-1.,1.));};

   TTree* effPid = par -> GetEff(0); // data, pi
   vector<TH1D*> hdat { hC("hdat") };
   TCut cutD = par -> Cpion  + TCut("Zpi<0") + TCut("0.05<Ptpi&&Ptpi<0.4");
   effPid->Draw("Cpi>>hdat",cutD,"goff");
   for ( auto& h : hdat ) {
      SetHstFace(h);
//       cout << h->GetName() << " Entries= " << h->GetEntries() << endl;
   }

   TTree* effPimc = par->GetEff(1); // MC, pi
   vector<TH1D*> hmc { hC("hmc"), hC("hmcBG") };
   effPimc->Draw("Cpi>>hmc",cutD,"goff");
   effPimc->Draw("Cpi>>hmcBG",cutD + !(par -> Cmcsig),"goff");
//    cout << "MC Entries= " << hmc[0]->GetEntries() << endl;

   double bg = hmc[1] -> Integral(); // [bin1,bin2]
   double sum = hmc[0] -> Integral();
   printf(" %i %s background is %.1f%%",date,__func__,bg/sum*100);

   // normalization on DATA
   double scale = hdat[0] -> Integral(5,36) / hmc[0] -> Integral(5,36);
   for ( auto& h : hmc ) {
      h -> Scale(scale);
   }

   TLegend* leg = new TLegend(0.11,0.69,0.38,0.89);
   leg->SetHeader( Form("%i  #pi^{-}",date),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,600);
   c1->cd(0);
   gPad->SetGrid();

   // Data
   double ymax=1.5*hdat[0]->GetMaximum();
   hdat[0] -> SetMaximum(ymax);
   hdat[0] -> SetTitle(
         ";cos(#Theta_{#pi}) "
         ";Entries/0.05 "
                       );
   hdat[0] -> GetYaxis() -> SetMaxDigits(3);
   hdat[0] -> GetYaxis() -> SetTitleOffset(1.);
   hdat[0] -> SetLineWidth(2);
   hdat[0] -> SetMarkerStyle(20);
   hdat[0] -> SetLineColor(kBlack);
   hdat[0] -> Draw("E");

   hmc[0] -> SetLineWidth(1);
   hmc[0] -> SetLineColor(kRed+2);
   hmc[0] -> Draw("HIST,SAME");

   hmc[1] -> SetLineWidth(1);
   hmc[1] -> SetLineColor(kBlue+1);
   hmc[1] -> SetFillStyle(3001);
   hmc[1] -> SetFillColor(kBlue+1);
   hmc[1] -> Draw("HIST,SAME");

   leg->AddEntry(hdat[0], "Data","LEP");
   leg->AddEntry(hmc[0], "MC","L");
   leg->AddEntry(hmc[1], "MC background","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_CPim_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void plot_Weights(int date) {
//-------------------------------------------------------------------------
   Params* par = new Params(date);

   vector<string> hnmc {"MC_3_WK", "MC_5_WPi"};
   vector<TH1D*> hmc;
   get_hst(par, 1, hnmc, hmc);

   TLegend* leg = new TLegend(0.59,0.79,0.89,0.89);
   leg->SetHeader( (string("#bf{")+to_string(date)+"}").c_str(),"C");

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,600);
   c1 -> Divide(2,1);
//    gPad->SetGrid();

   // Kaons
   c1 -> cd(1);
//    hmc[0] -> SetAxisRange(-0.15,0.15,"X");
   hmc[0] -> SetTitle(
         ";weights for K"
                       );
   hmc[0] -> SetLineWidth(2);
   hmc[0] -> SetLineColor(kBlue+1);
   hmc[0] -> GetYaxis() -> SetMaxDigits(3);
   hmc[0] -> Draw("HIST");
   leg->Draw();

   // Pions
   c1 -> cd(2);
   hmc[1] -> SetTitle(
         ";weights for #pi"
                       );
   hmc[1] -> SetLineWidth(2);
   hmc[1] -> SetLineColor(kRed+1);
   hmc[1] -> GetYaxis() -> SetMaxDigits(3);
   hmc[1] -> Draw("HIST");
   leg->Draw();

//    leg->AddEntry(hmc[0], "Weights for K","L");
//    leg->AddEntry(hmc[1], "Weights for #pi","L");

   gPad->RedrawAxis();
   c1->Update();
   string pdf = string("trkeff_WW_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}
//-------------------------------------------------------------------------
void trk_eff_sel() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(62);
   gStyle->SetStatFont(62);

//    plot_pi0(2012);
//    plot_pi0(2009);

   // ---- K+ K- ----
//    plot_MmisK(2012);
//    plot_MmisK(2009);

//    plot_dPK(2012);
//    plot_dPK(2009);
//    plot_dThK(2012);
//    plot_dThK(2009);

//    plot_PtKp(2012);
//    plot_PtKp(2009);
//    plot_PtKm(2012);
//    plot_PtKm(2009);
//    plot_CKp(2012);
//    plot_CKp(2009);
//    plot_CKm(2012);
//    plot_CKm(2009);

   // ---- pi+ pi- ----
//    plot_MinvJ(2012);
//    plot_MinvJ(2009);

//    plot_MmisP(2012);
//    plot_MmisP(2009);

//    plot_dPpi(2012);
//    plot_dPpi(2009);
//    plot_dThPi(2012);
//    plot_dThPi(2009);

//    plot_PtPip(2012);
//    plot_PtPip(2009);
//    plot_PtPim(2012);
//    plot_PtPim(2009);
//    plot_CPip(2012);
//    plot_CPip(2009);
//    plot_CPim(2012);
//    plot_CPim(2009);

   // other
//    plot_Weights(2012); // do not use
//    plot_Weights(2009);
}