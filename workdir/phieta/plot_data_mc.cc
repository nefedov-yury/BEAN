// plot Data vs MC distributions after all cuts applied
// -> data_vs_mc_XXXX.pdf

// {{{1 helper functions and constants
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
vector<TH1D*> get_hists(string name, string var, bool isMC) {
//--------------------------------------------------------------------
   string fname("Ntpls/");
   if( !isMC ) {
      fname += string("ntpl_") + name + string(".root");
   } else {
      fname += string("ntpl_mcgpj_") + name + string(".root");
   }
   // cout << " file: " << fname << endl;

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("SelectKKgg");
   TTree *a4c = (TTree*)gDirectory->Get("a4c");

   string title;
   int Nbins = 100;
   double Vmin = 0, Vmax = 0;
   string var2;
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
   } else if ( var == "Pk" ) {
      Nbins = 80; Vmin = -1.0; Vmax = 1.0;
      title = ";P(K^{#minus}) and P(K^{#plus}), GeV/c;"
         "Events/25 MeV/c";
      var = "Pkp";
      var2 = "-Pkm";
   } else if ( var == "Ck" ) {
      Nbins = 40; Vmin = -1.; Vmax = 1.;
      title = ";cos(#Theta(K^{#pm}));Events/0.05";
      var = "Ckp";
      var2 = "Ckm";
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
   // c_here += c_phiT;              // Mkk in [1.01, 1.03GeV]
   // c_here += c_cpgg;              // central part of Mgg
   // c_here += c_sbgg;              // side-band

   vector<TH1D*> hst(2,nullptr);
   string hncp = var + ((isMC) ? "_mc" : "") + "_cp";
   string hnsb = var + ((isMC) ? "_mc" : "") + "_sb";
   hst[0] = new TH1D(hncp.c_str(),title.c_str(), Nbins,Vmin,Vmax);
   hst[1] = new TH1D(hnsb.c_str(),"Side-band", Nbins,Vmin,Vmax);
   a4c->Draw( (var+">>"+hncp).c_str(), c_here+c_cpgg, "goff");
   a4c->Draw( (var+">>"+hnsb).c_str(), c_here+c_sbgg, "goff");
   if ( !var2.empty() ) {
      a4c->Draw( (var2+">>+"+hncp).c_str(), c_here+c_cpgg, "goff");
      a4c->Draw( (var2+">>+"+hnsb).c_str(), c_here+c_sbgg, "goff");
   }

   return hst;
}

//--------------------------------------------------------------------
tuple<vector<TH1D*>,vector<TH1D*> > HstDM(string name, string var) {
//--------------------------------------------------------------------
   auto hd = get_hists(name,var,false); // data
   auto hm = get_hists(name,var,true);  // MC

   // normaliza MC on data:
   double Ndat = hd[0]->Integral() - hd[1]->Integral();
   double Nmc  = hm[0]->Integral() - hm[1]->Integral();
   double scale = Ndat/Nmc;
   hm[0]->Scale( scale );
   hm[1]->Scale( scale );

   return make_tuple(hd,hm);
}

// {{{1 Plot for slides
//--------------------------------------------------------------------
void PlotDataMc(string name, string title) {
//--------------------------------------------------------------------
   string pdf = "data_vs_mc_" + name + ".pdf";
   // string pdf = "data_vs_mc_" + name + "_T.pdf"; // tight Mkk cut

   vector<TH1D*> Dhst, Mhst;
   vector<string> vars {"Pkp", "Pkm", "Ckp", "Ckm", "Pgg", "Cgg" };
   for ( auto var : vars ) {
      vector<TH1D*> hd, hm;
      tie(hd,hm) = HstDM(name,var);
      Dhst.insert( end(Dhst), begin(hd), end(hd) );
      Mhst.insert( end(Mhst), begin(hm), end(hm) );
   }

   int Ncp = Dhst[0]->GetEntries();
   int Nsb = Dhst[1]->GetEntries();

   TCanvas* c1 = new TCanvas("c1","...",0,0,1100,500);
   c1->Divide(3,2);

   TLegend* leg = new TLegend(0.12,0.65,0.52,0.92);
   string head(Form("#bf{Data vs MC %s}",title.c_str()));
   leg->SetHeader(head.c_str(),"C");
   leg->AddEntry(Dhst[0],Form("Data: %i events",Ncp),"EP");
   leg->AddEntry(Dhst[1],Form("Side-band: %i events",Nsb),"F");
   leg->AddEntry(Mhst[0],"MC signal #phi#eta","L");

   vector<int> npad { 1, 2, 4, 5, 3, 6 };
   vector<TPaveText*> pt(6,nullptr);
   for ( int i = 0; i < 6; ++i ) {
      pt[i] = new TPaveText(0.80,0.77,0.89,0.89,"NDC");
      pt[i]->SetTextAlign(12);
      pt[i]->SetTextFont(42);
      if ( i%3 == 0 ) {
         pt[i]->AddText(Form("#color[%i]{#bf{K^{#plus}}}",kRed+1));
      } else if ( i%3 == 1 ) {
         pt[i]->AddText(Form("#color[%i]{#bf{K^{#minus}}}",kRed+1));
      } else if ( i%3 == 2 ) {
         pt[i]->AddText(Form("#color[%i]{#bf{ #eta}}",kRed+1));
      }
   }
   vector<double> Hmax;
   if ( name == "3080_rs" ) {
         Hmax = { 37, 37, 16, 16, 83, 16 };
   } else if ( name == "3097" ) {
         Hmax = { 70, 70, 30, 30, 175, 30 };
   } else if ( name == "J4" ) {
         Hmax = { 105, 105, 40, 40, 300, 40 };
   }

   int Nhst = Dhst.size();
   for ( int j = 0; j < Nhst; j += 2 ) {
      int np = npad[j/2];
      c1->cd(np);
      gPad->SetGrid();
      Dhst[j]->SetMarkerStyle(20);
      Dhst[j]->SetMarkerSize(0.7);
      Dhst[j]->SetMarkerColor(kBlue+2);
      if ( !Hmax.empty() ) {
         Dhst[j]->SetMaximum( Hmax[j/2] );
      }
      SetHstFace(Dhst[j]);
      Dhst[j]->Draw("E1,P");

      Dhst[j+1]->SetLineColor(kRed);
      Dhst[j+1]->SetFillStyle(3001);
      Dhst[j+1]->SetFillColor(kRed);
      Dhst[j+1]->Draw("SAME");

      Mhst[j]->SetLineWidth(2);
      Mhst[j]->SetLineColor(kGreen+2);
      Mhst[j]->DrawCopy("HIST SAME");

      if ( np == 1 ) {
         leg->Draw();
      }
      pt[np-1]->Draw();

      c1->Update();
   }

   c1->cd();
   c1->Update();
   c1->Print(pdf.c_str());
}

// {{{1 Plot for memo
//--------------------------------------------------------------------
void PlotDataMcMEMO(string name, string var, string pdf) {
//--------------------------------------------------------------------
   vector<TH1D*> hd, hm; // data & mc
   tie(hd,hm) = HstDM(name,var);

   int maxbin = hd[0]->GetMaximumBin();
   double hmax = 0;
   hmax = hd[0]->GetBinContent(maxbin)+hd[0]->GetBinError(maxbin);
   hmax *= 1.15;

   int Ncp = hd[0]->GetEntries();
   int Nsb = hd[1]->GetEntries();

   double X1=0.30, Y1=0.75, X2=0.70, Y2=0.89; // legend
   if ( var == "Pgg" || var == "Cgg" ) {
      X1=0.12; X2=0.52;
   }
   TLegend* leg = new TLegend(X1,Y1,X2,Y2);
   string title;
   if ( name == "3080_rs" ) {
      title = "3080MeV R-scan";
   } else if ( name == "3097" ) {
      title = "3097MeV (2012)";
   } else if ( name == "J4" ) {
      title = "3097MeV (2018)";
   }
   leg->SetHeader(title.c_str(),"C");
   leg->AddEntry(hd[0],Form("Data: %i events",Ncp),"E");
   leg->AddEntry(hd[1],Form("Side-band: %i events",Nsb),"F");
   leg->AddEntry(hm[0],"MC signal #phi#eta","L");

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(hd[0]);
   hd[0]->SetLineColor(kBlack);
   hd[0]->SetLineWidth(2);
   hd[0]->GetYaxis()->SetTitleOffset(1.25);
   if ( hmax > 1 ) {
      hd[0]->SetMaximum(hmax);
   }
   hd[0]->Draw("E1");

   hd[1]->SetLineColor(kBlue+2);
   hd[1]->SetFillStyle(3001);
   hd[1]->SetFillColor(kBlue+1);
   hd[1]->Draw("SAME");

   hm[0]->SetLineWidth(2);
   hm[0]->SetLineColor(kGreen+2);
   hm[0]->DrawCopy("HIST SAME");

   leg->Draw();

   c1->Update();
   c1->Print(pdf.c_str());
}

// {{{1 Main
//--------------------------------------------------------------------
void plot_data_mc() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(42);
   // gStyle->SetOptStat(111110);
   // gStyle->SetOptFit(0);

   // PlotDataMc("3080_rs","3080MeV R-scan");
   // PlotDataMc("J4","3097MeV (2018)");
   // PlotDataMc("3097","3097MeV (2012)");

   // PlotDataMcMEMO("3080_rs","Pgg","Pgg_3080rs.pdf");
   // PlotDataMcMEMO("3080_rs","Cgg","Cgg_3080rs.pdf");
   // PlotDataMcMEMO("3080_rs","Pk","Pk_3080rs.pdf");
   // PlotDataMcMEMO("3080_rs","Ck","Ck_3080rs.pdf");

   // PlotDataMcMEMO("J4","Pgg","Pgg_3097J.pdf");
   // PlotDataMcMEMO("J4","Cgg","Cgg_3097J.pdf");
   // PlotDataMcMEMO("J4","Pk","Pk_3097J.pdf");
   // PlotDataMcMEMO("J4","Ck","Ck_3097J.pdf");
}
