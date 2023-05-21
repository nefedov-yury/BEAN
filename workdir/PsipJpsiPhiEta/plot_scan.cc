// script to plot scan results:
// 'χ2 = -2log(L/L_{max})' as a function of the vartheta-angle

// {{{1 helper function
//--------------------------------------------------------------------
void SetFace(TGraph* gr) {
//--------------------------------------------------------------------
   TAxis* X = gr->GetXaxis();
   if ( X ) {
      X->SetLabelFont(62);
      X->SetLabelSize(0.04);
      X->SetTitleFont(62);
      X->SetTitleSize(0.04);
      X->CenterTitle();
   }
   TAxis* Y = gr->GetYaxis();
   if ( Y ) {
      Y->SetLabelFont(62);
      Y->SetLabelSize(0.04);
      Y->SetTitleFont(62);
      Y->SetTitleSize(0.04);
      Y->CenterTitle();
   }
}

// {{{1 plot scan in given range
//--------------------------------------------------------------------
void plot_3_scans(string pdf) {
//--------------------------------------------------------------------
   // name of folder with root files
   string dir("prod_v709/");

   vector<int> Years {2009,2012,2021};
   vector<TGraph*> Vgr( Years.size(), nullptr );
   for ( size_t i = 0; i < Years.size(); ++i  ) {
      int date = Years[i];
      string fname( Form("mkk%02i_ifrSB_scan_ar.root",date%100) );
      fname = dir + fname;
      cout << fname << endl;
      TFile* froot = TFile::Open(fname.c_str(),"READ");
      if( froot == 0 ) {
         printf("can not open %s\n",fname.c_str());
         exit(EXIT_FAILURE);
      }
      froot->cd();
      TGraph* tmp = (TGraph*) gDirectory->Get("Graph"); // *TODO
      if( tmp == 0 ) {
         printf("can not read graph in file %s\n",fname.c_str());
         gDirectory->ls();
         exit(EXIT_FAILURE);
      }
      Vgr[i] = (TGraph*) tmp->Clone( Form("Scan_%02i",date%100) );
      froot->Close();
   }

   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","phase scan",0,0,900,900);
   c1->cd();
   gPad->SetGrid();

   auto& gr = Vgr[0];
   gr->SetTitle(";#vartheta, degrees;#it{#minus2#kern[0.2]"
         "{l}og(L#kern[0.5]{/}#kern[0.1]{L_{max}}#kern[0.2]{)}}");
   SetFace(gr);
   gr->SetMarkerSize(0.7);
   gr->SetMarkerColor(kGreen+2);
   gr->SetLineColor(kGreen+2);
   // gr->GetYaxis()->SetTitleOffset(1.2);

   gr->GetXaxis()->SetLimits(-90,90);
   gr->SetMinimum(-0.1);
   gr->SetMaximum(9.9);

   gr->Draw("APL");

   auto& gr1 = Vgr[1];
   gr1->SetMarkerSize(0.7);
   gr1->SetMarkerColor(kBlue);
   gr1->SetLineColor(kBlue);
   gr1->Draw("SAME,PL");

   auto& gr2 = Vgr[2];
   gr2->SetMarkerSize(0.7);
   gr2->SetMarkerColor(kRed+1);
   gr2->SetLineColor(kRed+1);
   gr2->Draw("SAME,PL");

   TLegend* leg = new TLegend(0.4,0.75,0.6,0.89);
   for ( size_t i = 0; i < Years.size(); ++i  ) {
      leg->AddEntry( Vgr[i],Form("%i",Years[i]),"LP" );
   }
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 main
//--------------------------------------------------------------------
void plot_scan() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   //   gStyle->SetOptFit(112);
   //   gStyle->SetStatX(0.48);
   //   gStyle->SetStatY(0.88);

   plot_3_scans("Scan_ifrSB.pdf");
}

