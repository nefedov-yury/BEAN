// data vs MC for  Jpsi -> phi eta  selection (see cuts.h)
// -> variable_YEAR.pdf

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
TH1D* get_hist(string fname, string var, string hname, int type=0) {
//-------------------------------------------------------------------------
#include "cuts.h"

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

   string title;
   int Nbins = 100;
   double Vmin = 0, Vmax = 0;
   if ( var == "Cgg" ) {
      Nbins = 50; Vmin = -1; Vmax = 1;
      title = ";cos #Theta(#eta);Events/0.04";
   } else if ( var == "Pgg" ) {
      Nbins = 55; Vmin = 1.05; Vmax = 1.6;
      title = ";P(#eta), GeV/c;Events/0.01GeV/c";
      var="Ptgg/sqrt(1-Cgg*Cgg)";
   } else if ( var == "Ptk" ) {
      Nbins = 110; Vmin = -1.1; Vmax = 1.1;
      title = ";P_{t}(K), GeV/c;Events/0.02GeV/c";
      var="Ptkp";
   } else if ( var == "Ckk" ) {
      Nbins = 40; Vmin = -1; Vmax = 1;
      title = ";cos #Theta(K^{#pm});Events/0.05";
      var="Ckp";
   } else {
      cerr << " unknown var=" << var << endl;
      exit(0);
   }

   TH1D* hst = new TH1D(hname.c_str(),title.c_str(),Nbins,Vmin,Vmax);

   TCut c_here = c_Mrec+c_chi2+c_phi;
   if ( type == 0 ) {           // central part
      c_here += c_cpgg;
   } else if ( type == 1) {     // side-band
      c_here += c_sbgg;
   } else if ( type == -1) {     // MC: background
      c_here += c_cpgg+TCut("decj!=68");
   }

   string select = var + string(">>") + hname;
   a4c->Draw(select.c_str(),c_here,"goff");

   if ( var == "Ptkp" ) {
      select = "-Ptkm" + string(">>+") + hname;
      a4c->Draw(select.c_str(),c_here,"goff");
   } else if ( var == "Ckp" ) {
      select = "Ckm" + string(">>+") + hname;
      a4c->Draw(select.c_str(),c_here,"goff");
   }

   return hst;
}

//-------------------------------------------------------------------------
void fill_hists(string var, int date, TH1D* hst[]) {
//-------------------------------------------------------------------------
   vector<string> fnames = {
      "data_09psip_all.root",
      "data_12psip_all.root",
      "mcinc_09psip_all.root",
      "mcinc_12psip_all.root",
      "mcsig_kkmc_09.root",
      "mcsig_kkmc_12.root"
   };

   int idx = ((date == 2009) ? 0 : 1);
   for( int i = 0; i < 3; i++ ) {
      string fn = fnames.at(idx+2*i);
      string::size_type idx = fn.find(".");
      string hname=fn.substr(0,idx);
      hst[i] = get_hist(fn,var,hname);

      if ( i == 0 ) { // side-band:
         hname += "_SB";
         hst[3] = get_hist(fn,var,hname,1);
      } else if ( i == 1 ) { // bkg from inclusive MC
         hname += "_BKG";
         hst[4] = get_hist(fn,var,hname,-1);
      }
   }

   // normaliza MC on data:
   double Ndat = hst[0] -> Integral();
   double NmcI = hst[1] -> Integral();
   hst[1]->Scale( Ndat / NmcI );
   hst[4]->Scale( Ndat / NmcI ); // bkg from inclusive MC
   double NmcS = hst[2] -> Integral();
   hst[2]->Scale( Ndat / NmcS );

   // debug print out:
   double Nsb = hst[3] -> Integral();
   printf( " N(side-band)= %.1f (%.2f%%)\n", Nsb, 100*Nsb/Ndat );
   double Nbkg = hst[4] -> Integral();
   printf( " N(bkg)= %.1f (%.2f%%)\n", Nbkg, 100*Nbkg/Ndat );

}

//-------------------------------------------------------------------------
void plot_var(string var, int date, string pdf) {
//-------------------------------------------------------------------------
   TH1D* hst[10];
   fill_hists( var, date, hst);

//    double max_dat = hst[0] -> GetMaximum();
//    double max_mc = hst[1] -> GetMaximum();
//    double max_sig = hst[2] -> GetMaximum();

   TLegend* leg = new TLegend(0.59,0.74,0.89,0.89);
   if ( var == "Cgg" ) {
      hst[0] -> SetMaximum( 1.35 * hst[0] -> GetMaximum() );
   } else if ( var == "Ptk" ) {
      leg = new TLegend(0.35,0.74,0.65,0.89);
   } else if ( var == "Ckk" ) {
      hst[0] -> SetMaximum( 1.35 * hst[0] -> GetMaximum() );
   }

   TCanvas* c1 = new TCanvas("c1","note",0,0,900,900);

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

   // 2) signal MC
   hst[2] -> SetLineColor(kGreen+2);
   hst[2] -> SetLineWidth(2);
   hst[2] -> Draw("SAME HIST");
   leg->AddEntry(hst[2], "MC signal #phi#eta", "L");

   // 3) side-band
   hst[3] -> SetLineColor(kBlue+2);
   hst[3] -> SetFillColor(kBlue+1);
   hst[3] -> SetFillStyle(3001);
   hst[3] -> Draw("SAME HIST");
   leg->AddEntry(hst[3],"Data: side-band","F");

   // ?) all inclusive MC ???
//    hst[1] -> SetLineColor(kRed+1);
//    hst[1] -> SetLineWidth(2);
//    hst[1] -> Draw("SAME HIST");
//    leg->AddEntry(hst[1], "MC incl.", "L");

   // 4) background from inclusive MC
   hst[4] -> SetLineColor(kRed+1);
   hst[4] -> SetLineWidth(2);
   hst[4] -> SetLineStyle(kDashed);
   hst[4] -> Draw("SAME HIST");
   leg->AddEntry(hst[4], "Bkg. from incl. MC", "L");

   leg->Draw();

   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void kkgg_dataMC() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
//    gStyle->SetOptFit(112); // print all parameters (fixed)
//    gStyle->SetFitFormat(".3g");
//    gStyle->SetLegendFont(62);
   gStyle->SetLegendTextSize(0.03);
   gStyle->SetStatFont(62);

   //----------------------------------------------------------------
//    plot_var("Pgg",2009, "Pgg_datMC09.pdf");
//    plot_var("Cgg", 2009, "Cgg_datMC09.pdf");
//    plot_var("Pgg",2012, "Pgg_datMC12.pdf");
//    plot_var("Cgg", 2012, "Cgg_datMC12.pdf");
//
//    plot_var("Ptk",2009, "Ptk_datMC09.pdf");
//    plot_var("Ckk",2009, "Ckk_datMC09.pdf");
//    plot_var("Ptk",2012, "Ptk_datMC12.pdf");
//    plot_var("Ckk",2012, "Ckk_datMC12.pdf");
   //----------------------------------------------------------------
}
