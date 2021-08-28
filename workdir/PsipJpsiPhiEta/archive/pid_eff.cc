// Study of PID reconstruction efficiency for pi and K
// -> pideff_[date].pdf

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
   TAxis* Y = hst->GetXaxis();
   if ( Y ) {
      Y->SetLabelFont(62);
      Y->SetLabelSize(0.04);
      Y->SetTitleFont(62);
      Y->SetTitleSize(0.04);
   }
   TAxis* Z = hst->GetXaxis();
   if ( Z ) {
      Z->SetLabelFont(62);
      Z->SetLabelSize(0.04);
      Z->SetTitleFont(62);
      Z->SetTitleSize(0.04);
   }
}

//----------------------------------------------------------------------
TH1D* HisDiff(TH1D* A, TH1D* B) {
//----------------------------------------------------------------------
// calculate difference = (A-B) with sigma^2 = sigma^2(A)-sigma^2(B)

   TH1D* diff = static_cast<TH1D*>(A->Clone("diff"));
   int nx = A->GetNbinsX();
   for ( int i = 0; i <= nx+1; ++i ) {
      diff->SetBinContent( i,
            A->GetBinContent(i) - B->GetBinContent(i)
                         );
      diff->SetBinError( i,
            sqrt( fabs(SQ(A->GetBinError(i)) - SQ(B->GetBinError(i))) )
                       );
   }
   return diff;

}

//-------------------------------------------------------------------------
void get_pid_eff( string fname, const string& ver,
                  vector<TH2D*>& hst, vector<TH1D*>& eff ) {
//-------------------------------------------------------------------------
   // name of histograms
   vector<string> hisn {
      "S6_Kp0",  "S6_Kp1",  "S6_Km0",  "S6_Km1",
      "S6_Pip0", "S6_Pip1", "S6_Pim0", "S6_Pim1",
   };
   if ( ver.size() != 0 ) { // version of selections
      for ( auto& hn : hisn ) { hn += ver; }
   }
   int Nh = hisn.size();

   hst.clear();
   hst.resize(Nh,nullptr);

   // DATA/MC:
   int isdm = 0; // 0,1 for data,MC
   if ( fname.find("mcinc") != string::npos ) {
      isdm = 1;
   }

   // name of folder with root files
   static const string dir("prod-te6/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }
   cout << " File: " << fname  << "\t Version: " << ver << endl;
   froot->cd("PipPimKpKm");

   for ( int i = 0; i < Nh; ++i ) {
      hst[i] = (TH2D*)gROOT->FindObject(hisn[i].c_str());
      if ( !hst[i] ) {
         cout << " can not find histo:" << hisn[i] << endl;
         return;
      }
      hst[i]->Sumw2(true);
   }

   // Get projections and calculate efficiencies
   eff.clear();
   eff.resize(Nh,nullptr);
   vector<TH1D*> tmp(2*Nh,nullptr); // proj. here
   for ( int i = 0; i < Nh; ++i ) {
      int isKp = (i/4)%2; // 0,1 for K,pi
      string pxn = hisn[i].substr(3) + string("_ct");
      string pyn = hisn[i].substr(3) + string("_pt");
//       cout << " pxn[" << i << "]= " << pxn << endl;
      int nx = hst[i]->GetNbinsX();
      int ny = hst[i]->GetNbinsY();

      tmp[2*i] = (TH1D*)hst[i]->ProjectionX(pxn.c_str(),1,ny);
      tmp[2*i+1] = (TH1D*)hst[i]->ProjectionY(pyn.c_str(),2,nx-1);

      if ( i%2 == 1 ) {
         eff[i-1] = (TH1D*)tmp[2*i]->Clone((string("eff_")+pxn).c_str());
         eff[i-1]->Add(tmp[2*i-2]);
         eff[i-1]->Divide(tmp[2*i],eff[i-1],1.,1.,"B");

         eff[i] = (TH1D*)tmp[2*i+1]->Clone((string("eff_")+pyn).c_str());
         eff[i]->Add(tmp[2*i-1]);
         eff[i]->Divide(tmp[2*i+1],eff[i],1.,1.,"B");
      }
   }

   // Attributes of histograms
   string DataMc[2] = {"Data", "MC"};
   string SXY[2] = {"PID-0","PID-1"};
   string Zs[2] = {"^{+}","^{-}"};
   string Kpi[2] = {"K","#pi"};
   for ( int i = 0; i < Nh; i++) {
      int isXY = i%2;     // 0,1 cos,Pt for eff or PID: 0/1
      int isZ  = (i/2)%2; // 0,1 for plus,minus
      int isKp = (i/4)%2; // 0,1 for K,pi
      string title = Kpi[isKp] + Zs[isZ] + " " +DataMc[isdm]
                   + SXY[isXY]
                   + "; cos(#Theta); P_{t}, GeV/c^{2}";
//       cout << " i= " << i << " title: " << title << endl;
      hst[i]->SetTitle(title.c_str());
      SetHstFace(hst[i]);
      hst[i]->GetXaxis()->SetTitleOffset(1.);
      hst[i]->GetYaxis()->SetTitleOffset(1.);

      string teff = Kpi[isKp] + Zs[isZ] + " " + DataMc[isdm]
                  + ((isXY==0) ? "; cos(#Theta)" : "; P_{t}, GeV/c^{2}")
                  + "; efficiency";
//       cout << " i= " << i << " teff: " << teff << endl;
      eff[i]->SetTitle(teff.c_str());
      SetHstFace(eff[i]);
      eff[i]->GetXaxis()->SetTitleOffset(1.2);
      eff[i]->GetYaxis()->SetTitleOffset(1.2);
   }
}

//-------------------------------------------------------------------------
void get_eff_rat( const vector<TH1D*>& effdat, const vector<TH1D*>& effmc,
                  vector<TH1D*>& rat) {
//-------------------------------------------------------------------------
   // calculate ratio DATA/MC
   int Nh = effdat.size();
   rat.clear();
   rat.resize(Nh,nullptr);
   for ( int i = 0; i < Nh; i++) {
      string effn(effdat[i]->GetName());
      string ratn = string("rat") + effn.substr(3);
//       cout << " ratn[" << i << "]= " << ratn << endl;
      rat[i]=(TH1D*)effdat[i]->Clone(ratn.c_str());
      rat[i]->Divide(effmc[i]);
   }

   // Attributes of ratio-histograms
   string Zs[2] = {"^{+}","^{-}"};
   string Kpi[2] = {"K","#pi"};
   for ( int i = 0; i < Nh; i++) {
      int isXY = i%2;     // 0,1 for cos,Pt for eff
      int isZ  = (i/2)%2; // 0,1 for plus,minus
      int isKp = (i/4)%2; // 0,1 for K,pi
      string trat = Kpi[isKp] + Zs[isZ] + " data/MC"
                  + ((isXY==0) ? "; cos(#Theta)" : "; P_{t}, GeV/c^{2}")
                  + "; data/MC";
//       cout << " i= " << i << " trat: " << trat << endl;

      rat[i]->SetTitle(trat.c_str());
      SetHstFace(rat[i]);
      rat[i]->GetXaxis()->SetTitleOffset(1.2);
      rat[i]->GetYaxis()->SetTitleOffset(1.2);
   }
}

//-------------------------------------------------------------------------
void plot_pid_ratio(int date, string ver, string pdf) {
//-------------------------------------------------------------------------
   static const vector<string> fnames {
        "data_09psip.root", "mcinc_09psip.root",
        "data_12psip.root", "mcinc_12psip.root"
   };

   vector<TH2D*> piddat, pidmc;
   vector<TH1D*> effdat, effmc;

   int idx = (date==2009) ? 0 : 2;
   get_pid_eff(fnames[idx], ver, piddat, effdat);
   get_pid_eff(fnames[idx+1], ver, pidmc, effmc);

   vector<TH1D*> rat; // data/MC
   get_eff_rat( effdat, effmc, rat );
   int Nh = rat.size();

   // Attributes of draw
   for (int i = 0; i < Nh; ++i ) {
      effdat[i]->SetLineColor(kBlue+1);
      effdat[i]->SetMarkerColor(kBlue+2);
      effdat[i]->SetMarkerStyle(22);
      effdat[i]->SetMarkerSize(0.9);

      effmc[i]->SetLineColor(kRed+1);
      effmc[i]->SetMarkerColor(kRed+2);
      effmc[i]->SetMarkerStyle(23);
      effmc[i]->SetMarkerSize(0.9);

      rat[i]->SetMarkerStyle(20);
      rat[i]->SetLineColor(kBlack);
   }

   TLegend* leg1 = new TLegend(0.40,0.20,0.72,0.40);
   leg1->AddEntry(effdat[0], Form("#color[%i]{data %i}",
                  effdat[0]->GetLineColor(),date),"LP");
   leg1->AddEntry(effmc[0], Form("#color[%i]{MC %i}",
                  effmc[0]->GetLineColor(),date),"LP");

   TLine* lC = new TLine(-0.80,1.,0.80,1.);
   TLine* lPtK = new TLine(0.1,1.,1.4,1.);
   TLine* lPtPi = new TLine(0.05,1.,0.4,1.);
   lC->SetLineColor(kRed+1);
   lC->SetLineWidth(2);
   lC->SetLineStyle(kSolid); // kDashed
   lPtK->SetLineColor(kRed+1);
   lPtK->SetLineWidth(2);
   lPtK->SetLineStyle(kSolid); // kDashed
   lPtPi->SetLineColor(kRed+1);
   lPtPi->SetLineWidth(2);
   lPtPi->SetLineStyle(kSolid); // kDashed

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->Divide(2,2);

   bool ispdf = pdf.size() > 0;
   if ( ispdf ) {
      c1->Print((pdf+"[").c_str()); // open pdf-file
   }

   int Nh2 = Nh/2;
   // efficiencies K
   for ( int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();
      if ( i%2 == 0 ) {
         effdat[i]->SetAxisRange(-0.799,0.799,"X");
         effdat[i]->SetAxisRange(0.8,1.0,"Y");
      } else {
         effdat[i]->SetAxisRange(0.1,1.399,"X");
         effdat[i]->SetAxisRange(0.8,1.0,"Y");
      }
      effdat[i]->Draw("E");
      effmc[i]->Draw("E SAME");
      leg1->Draw();
   }
   c1->Update();
   if ( ispdf ) {
      c1->Print(pdf.c_str()); // add to pdf-file
   } else {
      gPad->WaitPrimitive(); // pause
   }

   // ratio K
   TLine* line[2] {lC, lPtK};
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();
      if ( i%2 == 0 ) {
         rat[i]->SetAxisRange(-0.799,0.799,"X");
         rat[i]->SetAxisRange(0.90,1.1,"Y");
      } else {
         rat[i]->SetAxisRange(0.1,1.399,"X");
         rat[i]->SetAxisRange(0.90,1.1,"Y");
      }
      rat[i]->Draw("E");
      line[i%2]->Draw();
   }
   c1->Update();
   if ( ispdf ) {
      c1->Print(pdf.c_str()); // add to pdf-file
   } else {
      gPad->WaitPrimitive(); // pause
   }

   // efficiencies Pi
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();
      if ( i%2 == 0 ) {
         effdat[Nh2+i]->SetAxisRange(-0.799,0.799,"X");
         effdat[Nh2+i]->SetAxisRange(0.8,1.0,"Y");
      } else {
         effdat[Nh2+i]->SetAxisRange(0.05,0.399,"X");
         effdat[Nh2+i]->SetAxisRange(0.8,1.0,"Y");
      }
      effdat[Nh2+i]->Draw("E");
      effmc[Nh2+i]->Draw("E SAME");
      leg1->Draw();
   }
   c1->Update();
   if ( ispdf ) {
      c1->Print(pdf.c_str()); // add to pdf-file
   } else {
      gPad->WaitPrimitive(); // pause
   }

   // ratio Pi
   line[1] = lPtPi;
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();
      if ( i%2 == 0 ) {
         rat[Nh2+i]->SetAxisRange(-0.799,0.799,"X");
         rat[Nh2+i]->SetAxisRange(0.95,1.05,"Y");
      } else {
         rat[Nh2+i]->SetAxisRange(0.05,0.399,"X");
         rat[Nh2+i]->SetAxisRange(0.95,1.05,"Y");
      }
      rat[Nh2+i]->Draw("E");
      line[i%2]->Draw();
   }
   c1->Update();
   if ( ispdf ) {
      c1->Print(pdf.c_str()); // add to pdf-file
   }

   if ( ispdf ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
   }
}

//-------------------------------------------------------------------------
void plot_diff_ver(string fname) {
//-------------------------------------------------------------------------

   vector<TH2D*> hst; // unused here
   vector<TH1D*> effA, effB;

   get_pid_eff(fname, "_A", hst, effA);
   get_pid_eff(fname, "", hst, effB);
   int Nh = effA.size();

   vector<TH1D*> diff(Nh,nullptr);
   for ( int i = 0; i < Nh; ++i ) {
      diff[i] = HisDiff(effA[i],effB[i]);
      string name( effA[i]->GetName() );
      string sdif = string("dif") + name.substr(3);
      cout << " sdif[" << i << "]= " << sdif << endl;
      diff[i]->SetName( sdif.c_str() ); // rename
   }

   // Attributes of dif-histograms
   string Zs[2] = {"^{+}","^{-}"};
   string Kpi[2] = {"K","#pi"};
   for ( int i = 0; i < Nh; i++) {
      int isXY = i%2;     // 0,1 for cos,Pt for eff
      int isZ  = (i/2)%2; // 0,1 for plus,minus
      int isKp = (i/4)%2; // 0,1 for K,pi
      string tdiff = Kpi[isKp] + Zs[isZ] + " #Delta(efficiencies)"
                  + ((isXY==0) ? "; cos(#Theta)" : "; P_{t}, GeV/c^{2}")
                  + "; effA - effB";

      diff[i]->SetTitle(tdiff.c_str());
      SetHstFace(diff[i]);
//       rat[i]->GetXaxis()->SetTitleOffset(1.2);
//       rat[i]->GetYaxis()->SetTitleOffset(1.2);
   }

   TLine* lC = new TLine(-0.80,0.,0.80,0.);
   TLine* lPtK = new TLine(0.1,0.,1.4,0.);
   TLine* lPtPi = new TLine(0.05,0.,0.4,0.);
   lC->SetLineColor(kRed+1);
   lC->SetLineWidth(1);
   lC->SetLineStyle(kSolid); // kDashed
   lPtK->SetLineColor(kRed+1);
   lPtK->SetLineWidth(1);
   lPtK->SetLineStyle(kSolid); // kDashed
   lPtPi->SetLineColor(kRed+1);
   lPtPi->SetLineWidth(1);
   lPtPi->SetLineStyle(kSolid); // kDashed

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->Divide(2,2);

   int Nh2 = Nh/2;
   // efficiencies K
   TLine* line[2] {lC, lPtK};
   for ( int i = 0; i < Nh2; ++i ) {
      c1->cd(i%4+1);
      gPad->SetGrid();
      if ( i%2 == 0 ) {
         diff[i]->SetAxisRange(-0.799,0.799,"X");
         diff[i]->SetAxisRange(-0.02,0.02,"Y");
      } else {
         diff[i]->SetAxisRange(0.1,1.399,"X");
         diff[i]->SetAxisRange(-0.02,0.02,"Y");
      }
      diff[i]->Draw("E");
      line[i%2]->Draw();
   }
   c1->Update();
   gPad->WaitPrimitive(); // pause

   // efficiencies Pi
   line[1] = lPtPi;
   for ( int i = 0; i < Nh2; ++i ) {
      c1->cd(i%4+1);
      gPad->SetGrid();
      if ( i%2 == 0 ) {
         diff[Nh2+i]->SetAxisRange(-0.799,0.799,"X");
         diff[Nh2+i]->SetAxisRange(-0.01,0.01,"Y");
      } else {
         diff[Nh2+i]->SetAxisRange(0.05,0.399,"X");
         diff[Nh2+i]->SetAxisRange(-0.01,0.01,"Y");
      }
      diff[Nh2+i]->Draw("E");
      line[i%2]->Draw();
   }
   c1->Update();
}

//-------------------------------------------------------------------------
void pid_eff() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(62);
//    gStyle->SetLegendTextSize(0.03);

/*
   int date = 2012;
//    int date = 2009;
//    string ver = "";
   string ver = "_A";

//    string pdf = ""; // debug
   string pdf = string("pid_rat_") + to_string(date) + ver
              + string(".pdf");

   plot_pid_ratio(date,ver,pdf);
*/

  plot_diff_ver("data_12psip.root");

}
