// plot Data vs MC distributions after all cuts applied
// -> data_vs_mc_XXXX.pdf

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

//--------------------------------------------------------------------
vector<TH1D*> get_hists(string name, string var, bool isMC) {
//--------------------------------------------------------------------
   string fname("Ntpls/");
   if( !isMC ) {
      fname += string("ntpl_") + name + string(".root");
   } else {
      fname += string("ntpl_mcgpj_") + name + string(".root");
   }
//    cout << " file: " << fname << endl;

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("SelectKKgg");
   TTree *a4c = (TTree*)gDirectory->Get("a4c");

   string title;
   int Nbins = 100;
   double Vmin = 0, Vmax = 0;
   if ( var == "Pkp" ) {
      Nbins = 40; Vmin = 0.; Vmax = 1.;
      title = ";P(K^{#plus}), GeV/c;Events/25 MeV/c";
   } else if ( var == "Pkm" ) {
      Nbins = 40; Vmin = 0.; Vmax = 1.;
      title = ";P(K^{#minus}), GeV/c;Events/25 MeV/c";
   } else if ( var == "Ckp" ) {
      Nbins = 40; Vmin = -1.; Vmax = 1.;
      title = ";cos(#Theta(K^{#plus}));Events/0.05";
   } else if ( var == "Ckm" ) {
      Nbins = 40; Vmin = -1.; Vmax = 1.;
      title = ";cos(#Theta(K^{#minus}));Events/0.05";
   } else if ( var == "Pgg" ) {
      Nbins = 50; Vmin = 1.; Vmax = 1.5;
      title = ";P(#eta), GeV/c;Events/10 MeV/c";
   } else if ( var == "Cgg" ) {
      Nbins = 40; Vmin = -1.; Vmax = 1.;
      title = ";cos(#Theta(#eta));Events/0.05";
   } else {
      cerr << " unknown var=" << var << endl;
      exit(0);
   }

#include "cuts.h"
   TCut c_here = c_chi2;          // ch2 < 80
   if (isMC) {
      c_here += c_xisr + c_MCmkk; // X_isr>0.9 && mc_Mkk<1.08
   }
   c_here += c_phi;               // Mkk in [2*Mk, 1.08GeV]
//    c_here += c_cpgg;              // central part of Mgg
//    c_here += c_sbgg;              // side-band

   vector<TH1D*> hst(2,nullptr);
   string hncp = var + ((isMC) ? "_mc" : "") + "_cp";
   string hnsb = var + ((isMC) ? "_mc" : "") + "_sb";
   hst[0] = new TH1D(hncp.c_str(),title.c_str(), Nbins,Vmin,Vmax);
   hst[1] = new TH1D(hnsb.c_str(),"Side-band", Nbins,Vmin,Vmax);
   a4c->Draw( (var+">>"+hncp).c_str(), c_here+c_cpgg, "goff");
   a4c->Draw( (var+">>"+hnsb).c_str(), c_here+c_sbgg, "goff");

   return hst;
}

//--------------------------------------------------------------------
void fill_hists(string name,vector<TH1D*>& Dhst,vector<TH1D*>& Mhst) {
//--------------------------------------------------------------------
   vector<string> vars {"Pkp", "Pkm", "Ckp", "Ckm", "Pgg", "Cgg" };
   for ( auto var : vars ) {
      auto dh = get_hists(name,var,false); // data
      auto dm = get_hists(name,var,true);  // MC

      // normaliza MC on data:
      double Ndat = dh[0]->Integral() - dh[1]->Integral();
      double Nmc  = dm[0]->Integral() - dm[1]->Integral();
      double scale = Ndat/Nmc;
      dm[0]->Scale( scale );
      dm[1]->Scale( scale );

      Dhst.insert( end(Dhst), begin(dh), end(dh) );
      Mhst.insert( end(Mhst), begin(dm), end(dm) );
   }
}

//--------------------------------------------------------------------
void PlotDataMc(string name, string title) {
//--------------------------------------------------------------------
   string pdf = "data_vs_mc_" + name +".pdf";

   vector<TH1D*> Dhst, Mhst;
   fill_hists(name,Dhst,Mhst);
   int Ncp = Dhst[0] -> GetEntries();
   int Nsb = Dhst[1] -> GetEntries();

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   c1->Print((pdf+"[").c_str()); // just open pdf-file

   TLegend* leg = new TLegend(0.12,0.75,0.52,0.92);
   string head = string("Data vs MC ") + title;
   leg->SetHeader(head.c_str(),"C");
   leg->AddEntry(Dhst[0],Form("Data: %i events",Ncp),"EP");
   leg->AddEntry(Dhst[1],Form("Side-band: %i events",Nsb),"F");
   leg->AddEntry(Mhst[0],"MC signal #phi#eta","L");

   int Nhst = Dhst.size();
   for ( int j = 0; j < Nhst; j += 2 ) {
      Dhst[j]->SetMarkerStyle(21);
      Dhst[j]->SetMarkerColor(kBlue+2);
      Dhst[j]->Draw("E1,P");

      Dhst[j+1]->SetLineColor(kRed);
      Dhst[j+1]->SetFillStyle(3001);
      Dhst[j+1]->SetFillColor(kRed);
      Dhst[j+1]->Draw("SAME");

      Mhst[j]->SetLineWidth(2);
      Mhst[j]->SetLineColor(kGreen+2);
      Mhst[j]->DrawCopy("HIST SAME");

      leg->Draw();

      c1->Update();
      c1->Print(pdf.c_str()); // add to pdf-file
   }

   c1->Print((pdf+"]").c_str()); // just close pdf-file
}

//--------------------------------------------------------------------
void plot_data_mc() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
//    gStyle->SetOptStat(111110);
//    gStyle->SetOptFit(0);

//    PlotDataMc("3097","3097 MeV");
//    PlotDataMc("J4","3097 MeV (2018)");
//    PlotDataMc("3080_rs","3080 MeV (R-scan)");
//    PlotDataMc("3080_2019","3080 MeV (2019)");
}