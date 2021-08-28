// Plot the track reconstruction efficiency for pi and K
//

//-------------------------------------------------------------------------
void get_hst_eff( string fname, int subSB,
                  vector<TH2D*>& hst, vector<TH1D*>& eff ) {
//-------------------------------------------------------------------------
   // name of histograms
   vector<string> hisn {
      "S3_Kp5",  "S3_Kp6",  "S3_Km5",  "S3_Km6",
      "S5_Pip5", "S5_Pip6", "S5_Pim5", "S5_Pim6",
      "S3_Kp5sb",  "S3_Kp6sb",  "S3_Km5sb",  "S3_Km6sb",
      "S5_Pip5sb", "S5_Pip6sb", "S5_Pim5sb", "S5_Pim6sb",
   };
   for ( auto& hn : hisn ) { // second set of selections
      hn += string("_A");
//       hn += string("W_A"); // re-weighted
   }
   int Nh = hisn.size();
   if ( subSB == 0 ) {
      Nh = Nh/2; // do not read Side-Band histograms
   }

   hst.clear();
   hst.resize(Nh,nullptr);

   // DATA/MC:
   int isdm = 0; // 0,1 for data,MC
   if ( fname.find("mcinc") != string::npos ) {
      isdm = 1;
   }

   // name of folder with root files
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }
   cout << " file: " << fname << endl;
   froot->cd("PipPimKpKm");

   for ( int i = 0; i < Nh; ++i ) {
      hst[i] = (TH2D*)gROOT->FindObject(hisn[i].c_str());
      if ( !hst[i] ) {
         cout << " can not find histo:" << hisn[i] << endl;
         return;
      }
      hst[i]->Sumw2(true);
   }

   if ( subSB != 0 ) { // Subtract side-band
      int Nh2 = Nh/2;
      for ( int i = 0; i < Nh2; ++i ) {
         hst[i]->Add(hst[i+Nh2],-1.);
         // check negative content ???
         int nx = hst[i]->GetNbinsX();
         int ny = hst[i]->GetNbinsY();
         for ( int jx = 1; jx <= nx; jx++ ) {
            for ( int jy = 1; jy <= ny; jy++ ) {
               double t = hst[i]->GetBinContent(jx,jy);
               if ( t < 0 ) {
                  cout << " NEGATIVE BIN in: "<<hst[i]->GetTitle() << endl
                       << " x= " << hst[i]->GetXaxis()->GetBinCenter(jx)
                       << " y= " << hst[i]->GetYaxis()->GetBinCenter(jy)
                       << " content= " << t << endl;
               }
            }
         }
      }
      Nh = Nh2;
      hst.resize(Nh);
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

//       if ( isKp == 1 ) { // pions
         tmp[2*i] = (TH1D*)hst[i]->ProjectionX(pxn.c_str(),1,ny);
//       } else { // Kaons
//          tmp[2*i] = (TH1D*)hst[i]->ProjectionX(pxn.c_str(),2,ny);
//       }
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
   string S56[2] = {"5 tracks","6 tracks"};
   string Zs[2] = {"^{+}","^{-}"};
   string Kpi[2] = {"K","#pi"};
   string Sb = "(side-band subtracted) ";
   for ( int i = 0; i < Nh; i++) {
      int is56 = i%2;     // 0,1 for 5,6 tracks or cos,Pt for eff
      int isZ  = (i/2)%2; // 0,1 for plus,minus
      int isKp = (i/4)%2; // 0,1 for K,pi
      string title = Kpi[isKp] + Zs[isZ] + " " +DataMc[isdm]
                   + ( (subSB != 0) ? Sb : " " )
                   + S56[is56]
                   + "; cos(#Theta); P_{t}, GeV/c^{2}";
//       cout << " i= " << i << " title: " << title << endl;
      hst[i]->SetTitle(title.c_str());

      hst[i]->GetXaxis()->SetLabelFont(62);
      hst[i]->GetYaxis()->SetLabelFont(62);
      hst[i]->GetZaxis()->SetLabelFont(62);
      hst[i]->GetXaxis()->SetLabelSize(0.04);
      hst[i]->GetYaxis()->SetLabelSize(0.04);
      hst[i]->GetZaxis()->SetLabelSize(0.04);
      hst[i]->GetXaxis()->SetTitleFont(62);
      hst[i]->GetYaxis()->SetTitleFont(62);
      hst[i]->GetXaxis()->SetTitleSize(0.04);
      hst[i]->GetYaxis()->SetTitleSize(0.04);
      hst[i]->GetXaxis()->SetTitleOffset(1.);
      hst[i]->GetYaxis()->SetTitleOffset(1.);

      string teff = Kpi[isKp] + Zs[isZ] + " " + DataMc[isdm]
                  + ( (subSB != 0) ? Sb : " " )
                  + ((is56==0) ? "; cos(#Theta)" : "; P_{t}, GeV/c^{2}")
                  + "; efficiency";
//       cout << " i= " << i << " teff: " << teff << endl;
      eff[i]->SetTitle(teff.c_str());

      eff[i]->GetXaxis()->SetLabelFont(62);
      eff[i]->GetYaxis()->SetLabelFont(62);
      eff[i]->GetXaxis()->SetLabelSize(0.04);
      eff[i]->GetYaxis()->SetLabelSize(0.04);
      eff[i]->GetXaxis()->SetTitleFont(62);
      eff[i]->GetYaxis()->SetTitleFont(62);
      eff[i]->GetXaxis()->SetTitleSize(0.04);
      eff[i]->GetYaxis()->SetTitleSize(0.04);
      eff[i]->GetXaxis()->SetTitleOffset(1.2);
      eff[i]->GetYaxis()->SetTitleOffset(1.2);
   }
}

//-------------------------------------------------------------------------
void plot_trk_eff() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(62);
//    gStyle->SetLegendTextSize(0.03);

   vector<TH2D*> hst;
   vector<TH1D*> eff;
   int subSB = 0; // 1 - subtract Side-Band

//    string fname("mcincl_09.root");
   string fname("mcincl_12.root");

   get_hst_eff(fname, subSB, hst, eff);
   int Nh = eff.size();

   // a4 210x297
   TCanvas* c1 = new TCanvas("c1","...",0,0,700,1000);
   c1->Divide(2,4);

   // efficiencies K
   int Nh2 = Nh/2;
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      gPad->SetGrid();
      if ( i%2 == 0 ) {
         eff[i]->SetAxisRange(-0.8,0.8,"X");
         eff[i]->SetAxisRange(0.9,1.0,"Y");
      } else {
         eff[i]->SetAxisRange(0.9,1.0,"Y");
      }
      eff[i]->SetLineColor(kBlue+1);
      eff[i]->SetMarkerColor(kBlue+2);
      eff[i]->SetMarkerStyle(22);
      eff[i]->SetMarkerSize(0.9);
      eff[i]->Draw("E");
   }
//    c1->Update();

   // efficiencies Pi
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+5);
      gPad->SetGrid();
      if ( i%2 == 0 ) {
         eff[Nh2+i]->SetAxisRange(-0.8,0.8,"X");
         eff[Nh2+i]->SetAxisRange(0.9,1.0,"Y");
      } else {
         eff[Nh2+i]->SetAxisRange(0.9,1.0,"Y");
      }
      eff[Nh2+i]->SetLineColor(kBlue+1);
      eff[Nh2+i]->SetMarkerColor(kBlue+2);
      eff[Nh2+i]->SetMarkerStyle(22);
      eff[Nh2+i]->SetMarkerSize(0.9);
      eff[Nh2+i]->Draw("E");
   }
   c1->Update();

}
