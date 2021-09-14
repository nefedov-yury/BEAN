// study the tail in M(K+K-) distribution
// 0) plot Mkk for data
// plot data vs MC for 1.1 < Mkk < 2.0 GeV/c^2

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
void plot_Mkk(int date, string pdf) {
//-------------------------------------------------------------------------
   vector<string> fnames = {
      "data_09psip_all.root",
      "data_12psip_all.root",
   };
   int idx = ((date == 2009) ? 0 : 1);

   // name of folder with root files
   static string dir("prod-10/");
   string fname = dir + fnames.at(idx);
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

//    int Nbins = 204;
//    double Vmin = 0.98, Vmax = 2.0;
   int Nbins = 184;
   double Vmin = 1.08, Vmax = 2.0;
   string title = ";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2}"
                  ";Entries / 5 MeV/c^{2}";
   TH1D* hst = new TH1D("mkk",title.c_str(),Nbins,Vmin,Vmax);

   // cut on recoil mass: [3.092, 3.102]
   TCut c_Mrec("abs(Mrec-3.097)<0.005");
   // chi^2 cut:
   TCut c_chi2("ch2<80");
   // Meta: central part
   // const double meta = 0.547862; //  547.862 +/- 0.017 MeV
   TCut c_cpgg("abs(Mgg-0.547862)<0.024");

   TCut c_here = c_Mrec+c_chi2+c_cpgg;

   string select = string("Mkk>>mkk");
   a4c->Draw(select.c_str(),c_here,"goff");

   TCanvas* c1 = new TCanvas("c1","tail",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   TLegend* leg = new TLegend(0.69,0.79,0.89,0.89);

   SetHstFace(hst);
   hst -> SetOption("E");
   hst -> SetLineWidth(2);
   hst -> SetLineColor(kBlack);
//    hst -> SetMarkerStyle(20);
   hst -> GetYaxis()->SetTitleOffset(1.25);

//    hst -> Draw("EP");
   hst -> Draw("HIST");

   leg->AddEntry(hst, (string("Data ")+to_string(date)).c_str(), "L");
   leg->Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
TH1D* get_hist(string fname, string hname, string var, bool mcsig) {
//-------------------------------------------------------------------------
   // name of folder with root files
   static string dir("prod-9/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   // cut on recoil mass: [3.092, 3.102]
   TCut c_Mrec("abs(Mrec-3.097)<0.005");
   // chi^2 cut:
   TCut c_chi2("ch2<80");
   // Meta: central part
   // const double meta = 0.547862; //  547.862 +/- 0.017 MeV
   TCut c_cpgg("abs(Mgg-0.547862)<0.024");

   // tail of Mkk: data & mc non-phi KKeta
   TCut c_tail1("Mkk>1.1&&Mkk<2.0");
   // tail of Mkk: mcsig
   TCut c_tail2("Mkk>1.04&&Mkk<1.08");
//    TCut c_tail2("Mkk>1.01&&Mkk<1.03"); // phi-peak

   TCut c_tail = ( (!mcsig) ? c_tail1 : c_tail2 );

   TCut c_here = c_Mrec+c_chi2+c_cpgg+c_tail;
   string title("1.1<Mkk<2.0");
   int Nbins = 100;
   double Vmin = 0, Vmax = 0;
   if ( var == "Cgg" ) {
      Nbins = 50; Vmin = -1; Vmax = 1;
      title += ";cos #Theta(#eta)";
      title += ";Events/0.04";
   } else if ( var == "Mgg" ) {
      const double meta = 0.547862; //  547.862 +/- 0.017 MeV
      Nbins = 50; Vmin = meta-0.075; Vmax = meta+0.075;
      title = string("M^{ inv}_{ #gamma#gamma} :") + title;
      title += ";M^{ inv}_{ #gamma#gamma}, GeV/c^{2}";
      title += ";Entries/3MeV/c^{2}";
      c_here = c_Mrec+c_chi2+c_tail; // no cut in Mgg
   } else if ( var == "Ckp" ) {
      Nbins = 50; Vmin = -1; Vmax = 1;
      title = string("K^{+} :") + title;
      title += ";cos #Theta(K^{+})";
      title += ";Events/0.04";
   } else if ( var == "Ckm" ) {
      Nbins = 50; Vmin = -1; Vmax = 1;
      title = string("K^{-} :") + title;
      title += ";cos #Theta(K^{-})";
      title += ";Events/0.04";
   } else if ( var == "Ptkp" ) {
      Nbins = 50; Vmin = 0.; Vmax = 1.5;
      title = string("K^{+} :") + title;
      title += ";P_{t}(K^{+}), GeV/c";
      title += ";Events/30MeV/c";
   } else if ( var == "Ptkm" ) {
      Nbins = 50; Vmin = 0.; Vmax = 1.5;
      title = string("K^{-} :") + title;
      title += ";P_{t}(K^{-}), GeV/c";
      title += ";Events/30MeV/c";
   } else {
      cerr << " unknown var=" << var << endl;
      exit(0);
   }

   TH1D* hst = new TH1D(hname.c_str(),title.c_str(),Nbins,Vmin,Vmax);
   string select = var + ">>" + hname;
   a4c->Draw(select.c_str(),c_here,"goff");

   return hst;
}

//-------------------------------------------------------------------------
TH1D* get_hist10(string fname, string hname, string var, bool mcsig) {
//-------------------------------------------------------------------------
   const double mjpsi  = 3.096916; // 3096.916  +/- 0.011   MeV
   const double meta   = 0.547862; // 547.862   +/- 0.017   MeV
   const double mk     = 0.493677; // 493.677   +/- 0.016   MeV

   // name of folder with root files
   static string dir("prod-10/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");
#include "p10a4c.h"

   // cut on recoil mass: [3.092, 3.102]
   auto c_Mrec = [](double Mrec)->bool{return fabs(Mrec-3.097)<0.005;};
   // chi^2 cut:
   auto c_chi2 = [](double ch2)->bool{return ch2<80;};
   // Meta: central part
   auto c_cpgg = [meta](double Mgg)->bool{return fabs(Mgg-meta)<0.024;};
   // tail of Mkk: data & mc non-phi KKeta
//    auto c_tail1 = [](double Mkk)->bool{return (Mkk>1.1&&Mkk<2.0);};
   auto c_tail1 = [](double Mkk)->bool{return (Mkk>1.08&&Mkk<1.3);};
//    auto c_tail1 = [](double Mkk)->bool{return (Mkk>0.98&&Mkk<1.1);}; // ALL!
   // tail of Mkk: mcsig
//    auto c_tail2 = [](double Mkk)->bool{return (Mkk>1.04&&Mkk<1.08);};
   // phi-peak
   auto c_tail2 = [](double Mkk)->bool{return (Mkk>1.01&&Mkk<1.03);};

   bool isCj = (var == "Cj");
   bool isCgg = (var == "Cgg");
   bool isCkk = (var == "Ckk");
   bool isMgg = (var == "Mgg");
   bool isCkp = (var == "Ckp");
   bool isCkm = (var == "Ckm");
   bool isPkp = (var == "Pkp");
   bool isPkm = (var == "Pkm");
   bool isCkpkm = (var == "Ckpkm");
   bool isCkpKK = (var == "CkpKK");

//    string title("1.1<Mkk<2.0");
   string title("1.08<Mkk<1.3");
//    string title("0.98<Mkk<1.1");
   int Nbins = 100;
   double Vmin = 0, Vmax = 0;
   if ( isCj ) {
      Nbins = 50; Vmin = -1; Vmax = 1;
      title = string("J/#psi: ") + title;
      title += ";cos #Theta(J/#psi)";
      title += ";Events/0.04";
   } else if ( isCgg ) {
      Nbins = 50; Vmin = -1; Vmax = 1;
      title = string("#eta(in J/#psi RF): ") + title;
      title += ";cos #Theta(#eta)";
      title += ";Events/0.04";
   } else if ( isCkk ) {
      Nbins = 50; Vmin = -1; Vmax = 1;
      title = string("angle of K^{+}K^{-} system in J/#psi RF: ") + title;
      title += ";cos #Theta(K^{+}K^{-})";
      title += ";Events/0.04";
   } else if ( isMgg ) {
      const double meta = 0.547862; //  547.862 +/- 0.017 MeV
      Nbins = 50; Vmin = meta-0.075; Vmax = meta+0.075;
      title = string("M^{ inv}_{ #gamma#gamma} :") + title;
      title += ";M^{ inv}_{ #gamma#gamma}, GeV/c^{2}";
      title += ";Entries/3MeV/c^{2}";
   } else if ( isCkp ) {
      Nbins = 50; Vmin = -1; Vmax = 1;
      title = string("K^{+}(in J/#psi RF): ") + title;
      title += ";cos #Theta(K^{+})";
      title += ";Events/0.04";
   } else if ( isCkm ) {
      Nbins = 50; Vmin = -1; Vmax = 1;
      title = string("K^{-}(in J/#psi RF): ") + title;
      title += ";cos #Theta(K^{-})";
      title += ";Events/0.04";
   } else if ( isPkp ) {
      Nbins = 30; Vmin = 0.; Vmax = 1.2;
      title = string("K^{+}(in J/#psi RF): ") + title;
      title += ";P(K^{+}), GeV/c";
      title += ";Events/40MeV/c";
   } else if ( isPkm ) {
      Nbins = 30; Vmin = 0.; Vmax = 1.2;
      title = string("K^{-}(in J/#psi RF): ") + title;
      title += ";P(K^{-}), GeV/c";
      title += ";Events/40MeV/c";
   } else if ( isCkpkm ) {
//       Nbins = 50; Vmin = -1; Vmax = 1;
      Nbins = 50; Vmin = 0; Vmax = 1;
      title = string("angle #alpha between K^{+} and K^{-} in J/#psi RF: ") 
            + title;
      title += ";cos #alpha(K^{+}K^{-})";
//       title += ";Events/0.04";
      title += ";Events/0.02";
   } else if ( isCkpKK ) {
      Nbins = 50; Vmin = -1; Vmax = 1;
      title = string("angle of K^{+} in K^{+}K^{-} RF"
                     " wrt direction of K^{+}K^{-} in J/Psi RF: ") + title;
      title += ";cos #Theta(K^{+} in KK)";
      title += ";Events/0.04";
   } else {
      cerr << " unknown var=" << var << endl;
      exit(0);
   }

   TH1D* hst = new TH1D(hname.c_str(),title.c_str(),Nbins,Vmin,Vmax);

   bool prt = false;
   Long64_t nentries = a4c->GetEntries();
   for ( Long64_t i=0; i < nentries; ++i ) {
      a4c -> GetEntry(i);

      // selection:
      if ( !c_Mrec(Mrec) ) continue;
      if ( !c_chi2(ch2) ) continue;
      if ( (!isMgg) && !c_cpgg(Mgg) ) continue;
      if ( !mcsig ) {
         if ( !c_tail1(Mkk) ) continue;
      } else {
         if ( !c_tail2(Mkk) ) continue;
      }

      TVector3 VPj; // J/Psi
      VPj.SetMagThetaPhi(Pj,acos(Cj),phij);
      TLorentzVector LVPj;
      LVPj.SetVectM(VPj,mjpsi);
      TVector3 boostJ = LVPj.BoostVector();
      if ( prt ) {
         cout << " Pj,Cj,phij= " << Pj << ", " << Cj << ", " << phij << endl;
         cout << " VPj= "; VPj.Print(); cout << endl;
         cout << " boostJ= "; boostJ.Print(); cout << endl;
      }

      double xvar = 0.;
      if ( isCj ) {
         xvar = Cj;
      } else if ( isMgg ) {
         xvar = Mgg;
      } else if ( isCgg ) {
         TVector3 VPgg;
         VPgg.SetMagThetaPhi(Pgg,acos(Cgg),phigg);
         TLorentzVector LVPgg;
         LVPgg.SetVectM(VPgg,meta);
         LVPgg.Boost(-boostJ);
         xvar = LVPgg.CosTheta();
      } else if ( isCkk ) {
         TVector3 VPkk; // K+K- system
         VPkk.SetMagThetaPhi(Pkk,acos(Ckk),phikk);
         TLorentzVector LVPkk;
         LVPkk.SetVectM(VPkk,Mkk);
         LVPkk.Boost(-boostJ);
         xvar = LVPkk.CosTheta();

      } else if ( isCkp || isPkp ) {
         TVector3 VPkp;
         VPkp.SetMagThetaPhi(Pkp,acos(Ckp),phikp);
         TLorentzVector LVPkp;
         LVPkp.SetVectM(VPkp,mk);
         LVPkp.Boost(-boostJ);
         if ( isCkp ) {
            xvar = LVPkp.CosTheta();
         } else {
            xvar = LVPkp.P();
         }
         if ( prt ) {
            cout << " Pkp,Ckp,phikp= " << Pkp << ", " << Ckp
                 << ", " << phikp << endl;
            cout << " VPkp= "; VPkp.Print(); cout << endl;
            cout << " cos= " << xvar << endl;
         }
      } else if ( isCkm || isPkm ) {
         TVector3 VPkm;
         VPkm.SetMagThetaPhi(Pkm,acos(Ckm),phikm);
         TLorentzVector LVPkm;
         LVPkm.SetVectM(VPkm,mk);
         LVPkm.Boost(-boostJ);
         if ( isCkm ) {
            xvar = LVPkm.CosTheta();
         } else {
            xvar = LVPkm.P();
         }
      } else if ( isCkpkm ) {
         TVector3 VPkp;
         VPkp.SetMagThetaPhi(Pkp,acos(Ckp),phikp);
         TLorentzVector LVPkp;
         LVPkp.SetVectM(VPkp,mk);
         LVPkp.Boost(-boostJ);
         VPkp = LVPkp.Vect();

         TVector3 VPkm;
         VPkm.SetMagThetaPhi(Pkm,acos(Ckm),phikm);
         TLorentzVector LVPkm;
         LVPkm.SetVectM(VPkm,mk);
         LVPkm.Boost(-boostJ);
         VPkm = LVPkm.Vect();

         xvar = VPkp.Dot(VPkm)/(VPkp.Mag()*VPkm.Mag());
      } else if ( isCkpKK ) {
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

         xvar = VPkp.Dot(VPkk)/(VPkp.Mag()*VPkk.Mag());
      }

      hst -> Fill(xvar);
      prt=false;
   }

   return hst;
}

//-------------------------------------------------------------------------
void plot_tail(int date, string var, string pdf) {
//-------------------------------------------------------------------------
   vector<string> fnames = {
      "data_09psip_all.root",
      "data_12psip_all.root",
      "mckketa_kkmc_09.root",
      "mckketa_kkmc_12.root",
      "mcsig_kkmc_09.root",
      "mcsig_kkmc_12.root",
//       "mcinc_09psip_all.root",
//       "mcinc_12psip_all.root",
   };

   TH1D* hst[10];
   int idx = ((date == 2009) ? 0 : 1);
   for( int i = 0; i < 3; ++i ) {
      string fn = fnames.at(idx+2*i);
      string::size_type idx = fn.find(".");
      string hname=fn.substr(0,idx);
//       hst[i] = get_hist(fn,hname,var,i==2);
      hst[i] = get_hist10(fn,hname,var,i==2);
   }

   // normaliza MC on data:
   double scale1 = 1., scale2 = 1.;
   if ( var == "Cj" || var == "CkpKK" ) {
      double Ndat = hst[0] -> Integral();
      double Nmc  = hst[1] -> Integral();
      double NmcS = hst[2] -> Integral();
      scale1 = Ndat / Nmc;
      scale2 = Ndat / NmcS;
   } else if ( var == "Cgg" || var == "Ckk" ) {
      double Ndat = hst[0] -> Integral();
      double Nmc  = hst[1] -> Integral();
      scale1 = Ndat / Nmc;
      double Ndat2 = hst[0] -> Integral(6,45);  // [bin1,bin2]
      double NmcS  = hst[2] -> Integral(6,45);
      scale2 = Ndat2 / NmcS;
   } else if ( var == "Mgg" ||
               var == "Pkp" || var == "Pkm" ||
               var == "Ckpkm"  ) {
      double Ndat = hst[0] -> GetMaximum();
      double Nmc  = hst[1] -> GetMaximum();
      double NmcS = hst[2] -> GetMaximum();
      scale1 = Ndat / Nmc;
      scale2 = Ndat / NmcS;
   } else if ( var == "Ckp" || var == "Ckm" ) {
      double Ndat = hst[0] -> Integral(6,45);  // [bin1,bin2]
      double Nmc  = hst[1] -> Integral(6,45);
      double NmcS = hst[2] -> Integral(6,45);
      scale1 = Ndat / Nmc;
      scale2 = Ndat / NmcS;
   }
   hst[1]->Scale( scale1 );
   hst[2]->Scale( scale2 );
   printf(" scale(KKeta)= %g,",scale1);
   printf(" scale(phieta)=%g\n",scale2);

   // legenda
   TLegend* leg = new TLegend(0.54,0.77,0.89,0.89);
   if ( var == "Ckp" || var == "Ckm" ) {
      delete leg;
      leg = new TLegend(0.33,0.77,0.68,0.89);
   } else if ( var == "Pkp" || var == "Pkm" || var == "Ckpkm" ) {
      delete leg;
      leg = new TLegend(0.11,0.77,0.46,0.89);
   }

   TCanvas* c1 = new TCanvas("c1","tail",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   // 1) data
   SetHstFace(hst[0]);
   hst[0] -> SetOption("E");
   hst[0] -> SetLineWidth(2);
   hst[0] -> SetLineColor(kBlack);
   hst[0] -> SetMarkerStyle(20);
   hst[0] ->GetYaxis()->SetTitleOffset(1.25);

   hst[0] -> Draw("EP");
   leg->AddEntry(hst[0], (string("Data ")+to_string(date)).c_str(), "EP");

   // 2) MC KKeta
   hst[1] -> SetLineColor(kBlue+1);
   hst[1] -> SetLineWidth(2);
   hst[1] -> Draw("SAME HIST");
   leg->AddEntry(hst[1], "MC non-#phi KK#eta (PHSP)", "L");

   // 2) MC phi eta
   hst[2] -> SetLineColor(kRed+1);
   hst[2] -> SetLineWidth(2);
   hst[2] -> Draw("SAME HIST");
//    leg->AddEntry(hst[2], "MC #phi#eta Mkk in [1.04,1.08]", "L");
   leg->AddEntry(hst[2], "MC #phi#eta Mkk in [1.01,1.03]", "L");

   leg->Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void Mkk_tail() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(62);
//    gStyle->SetLegendTextSize(0.03);
//    gStyle->SetStatFont(62);

   //----------------------------------------------------------------
//    plot_Mkk(2012, "Mkk_12.pdf");
//    plot_Mkk(2009, "Mkk_09.pdf");
   //----------------------------------------------------------------
//    plot_tail(2012, "Mgg", "Mgg_tail_12.pdf");
//    plot_tail(2009, "Mgg", "Mgg_tail_09.pdf");
   //----------------------------------------------------------------
//    plot_tail(2012, "Cj", "Cjpsi_tail_12.pdf");
   //----------------------------------------------------------------
//    plot_tail(2012, "Cgg", "CggJs_tail_12.pdf");
//    plot_tail(2009, "Cgg", "Cgg_tail_09.pdf");
//    plot_tail(2012, "Ckk", "CkkJs_tail_12.pdf");
   //----------------------------------------------------------------
//    plot_tail(2012, "Ckp", "CkpJs_tail_12.pdf");
//    plot_tail(2012, "Ckm", "CkmJs_tail_12.pdf");
//    plot_tail(2012, "Pkp", "Pkp_tail_12.pdf");
//    plot_tail(2012, "Pkm", "Pkm_tail_12.pdf");
   //----------------------------------------------------------------
//    plot_tail(2012, "Ckpkm", "CkpkmJs_tail_12.pdf");
   plot_tail(2012, "CkpKK", "CkpKK_tail_12.pdf");
   //----------------------------------------------------------------
//    plot_tail(2012, "CkpKK", "CkpKK_all_12.pdf");
}
