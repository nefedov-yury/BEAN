// pipi_dataMC.cc - plot data distributions vs MC
// for pi+ pi- in selected Psi(2S) -> Jpsi pi+pi-
//      -> VarName_YEAR.pdf

// pointer on function wich fill one histo
typedef TH1D* (*FILL_H)(string, string, string);

// {{{1 helper functions and constants
//--------------------------------------------------------------------
// global name of folder with root files
string Dir;

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

// {{{1 plot_hist() - functions
//--------------------------------------------------------------------
TH1D* get_hst(string fname, string hname) {
//--------------------------------------------------------------------
   fname = Dir + fname;
   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "ERROR in "<< __func__
         << ": can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TH1D* hst = (TH1D*)gROOT->FindObject(hname.c_str());
   if ( !hst ) {
      cout << " can not find " << hname << endl;
      exit(EXIT_FAILURE);
   }

   return hst;
}

//--------------------------------------------------------------------
void get_hist(int date, string hname, TH1D* hst[]) {
//--------------------------------------------------------------------
#include "norm.h"

   string datafile( Form("data_%02ipsip_all.root",date%100) );
   string mcincfile( Form("mcinc_%02ipsip_all.root",date%100) );

   hst[0]=get_hst(datafile, hname);
   hst[0]->SetMarkerStyle(20);
   hst[0]->SetMarkerSize(0.7);
   hst[0]->GetYaxis()->SetMaxDigits(3);
   // cout << " Integral(data)= " << hst[0]->Integral() << endl;

   hst[1]=get_hst("data_3650_2021.root", hname);
   hst[1]->Scale(C_Dat.at(date));
   hst[1]->SetLineColor(kBlue+1);
   hst[1]->SetLineWidth(2);
   // cout << " Integral(3650)= " << hst[1]->Integral() << endl;

   hst[2]=get_hst(mcincfile, hname);
   hst[2]->Scale(MC_Dat.at(date));
   hst[2]->SetLineColor(kRed+1);
   hst[2]->SetLineWidth(1);
   // cout << " Integral(mc)= " << hst[2]->Integral() << endl;
}

//--------------------------------------------------------------------
void plot_hist(string hname, int date) {
//--------------------------------------------------------------------
   TCanvas* c1 = new TCanvas(Form("c1_%s_%i",hname.c_str(),date),
         Form("%i %s",date,hname.c_str()),0,0,800,800);

   TH1D* hst[10];
   get_hist(date, hname, hst);
   hst[3]=(TH1D*)hst[2]->Clone("Sum");
   hst[3]->Add(hst[1]); // sum MC inclusive and continuum
   double hst_max = max(hst[0]->GetMaximum(),hst[3]->GetMaximum());
   double win_max = 1.15 * hst_max;

   // relative normalisation
   hst[3]->Scale( hst[0]->Integral()/hst[3]->Integral() );

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(2);

   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst[0]);
   hst[0]->SetAxisRange(0.,win_max,"Y");

   int pos = 0;
   if ( hname == "cosPM" ) {
      hst[0]->SetTitle(
            ";cos #Theta_{#pi^{#plus}#pi^{#minus}}"
            ";Events/0.02" );
      if ( date == 2009 ) {
         hst[0]->GetYaxis()->SetTitleOffset(1.2);
      }
   } else if ( hname == "invPM" ) {
      hst[0]->SetTitle(
            ";M^{ inv}_{#pi^{#plus}#pi^{#minus}}, GeV/c^{2}"
            ";Events/0.005 GeV/c^{2}" );
      if ( date == 2009 ) {
         hst[0]->GetYaxis()->SetTitleOffset(1.2);
      }
   } else if ( hname == "Pip_P" ) {
      pos = 1;
      box = nullptr;
      hst[0]->SetTitle(
            ";Momentum of #pi^{#plus}, GeV/c"
            ";Events/0.005 GeV/c" );
      hst[0]->GetYaxis()->SetTitleOffset(1.2);
   } else if ( hname == "Pim_P" ) {
      pos = 1;
      box = nullptr;
      hst[0]->SetTitle(
            ";Momentum of #pi^{#minus}, GeV/c"
            ";Events/0.005 GeV/c" );
      hst[0]->GetYaxis()->SetTitleOffset(1.2);
   } else if ( hname == "Pip_C") {
      pos = 2;
      hst[0]->SetTitle(
            ";cos #Theta_{#pi^{#plus}}"
            ";Events/0.02" );
      hst[0]->GetYaxis()->SetTitleOffset(1.2);
   } else if ( hname == "Pim_C") {
      pos = 2;
      hst[0]->SetTitle(
            ";cos #Theta_{#pi^{#minus}}"
            ";Events/0.02" );
      hst[0]->GetYaxis()->SetTitleOffset(1.2);
   } else if ( hname == "Bnch") {
      box = nullptr;
      hst[0]->SetTitle(
            ";Number of charged tracks"
            ";Events");
      hst[0]->SetAxisRange(1.,17.,"X");
      hst[0]->SetMinimum(1.);
      hst[0]->SetMaximum(); // default
      gPad->SetLogy(true);
      hst[0]->GetYaxis()->SetTitleOffset(1.2);
   }

   hst[0]->Draw("E"); // data
   if ( hname == "cosPM" ) {
      box->DrawBox(0.8,.0,1.0,win_max);
   } else if ( hname == "invPM" ) {
      double mk0    = 0.497611;
      double d = 0.008;
      box->DrawBox(mk0-d,.0,mk0+d,win_max);
   } else if ( hname == "Pip_C" || hname == "Pim_C" ) {
      gPad->SetLogy(false);
      if ( date == 2009 ) {
         hst[0]->GetYaxis()->SetNdivisions(1004);
      }
      box->DrawBox(-1.,0.,-0.8,win_max);
      box->DrawBox(0.8,0.,1.,win_max);
      // gPad->SetLogy(true);
      // hst[0]->SetAxisRange(1.,win_max,"Y");
      // box->DrawBox(-1.,1.,-0.8,win_max);
      // box->DrawBox(0.8,1.,1.,win_max);
   }
   hst[0]->Draw("E, SAME"); // data
   hst[3]->Draw("SAME,HIST"); // sum
   hst[1]->Draw("SAME,HIST"); // cont.
   gPad->RedrawAxis();

   TLegend* leg;
   if ( pos == 0 ) { // default
      leg = new TLegend(0.53,0.65+0.07*(box==nullptr),0.89,0.89);
   } else if ( pos == 1 ) {
      leg = new TLegend(0.60,0.75,0.92,0.92);
   } else if ( pos == 2 ) {
      leg = new TLegend(0.32,0.26,0.68,0.50);
   }
   leg->AddEntry(hst[0],
         (string("Data ")+to_string(date)).c_str(), "EP");
   leg->AddEntry(hst[3],
         Form("#color[%i]{MC + Continuum}",hst[2]->GetLineColor()),
         "L");
   leg->AddEntry(hst[1],
         Form("#color[%i]{Continuum}",hst[1]->GetLineColor()),
         "L");
   if ( box ) {
      leg->AddEntry(box, "Rejection area","F");
   }
   leg->Draw();

   c1->Update();
   string pdf = hname + "_" + to_string(date)+".pdf";
   c1->Print( pdf.c_str() );
}

// {{{1 plot_One and functions
//--------------------------------------------------------------------
TH1D* fill_nch(string fname, string hname, string hcut) {
//--------------------------------------------------------------------
   fname = Dir + fname;
   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "ERROR in "<< __func__
         << ": can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");

   TH1D* hst = new TH1D(hname.c_str(),
         ";Number of charged tracks"
         ";Events",
         20,0.5,20.5
         );
   hst->Sumw2();

   string dr = string("nch>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//--------------------------------------------------------------------
TH1D* fill_Nm(string fname, string hname, string hcut) {
//--------------------------------------------------------------------
   fname = Dir + fname;
   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "ERROR in "<< __func__
         << ": can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");

   TH1D* hst = new TH1D(hname.c_str(),
         ";Number of M^{rec}_{#pi^{#plus}#pi^{#minus}} combinations"
         ";Events",
         20,0.5,20.5
         );
   hst->Sumw2();

   string dr = string("@Mrec.size()+(Mrs>0)>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//--------------------------------------------------------------------
TH1D* fill_ppip(string fname, string hname, string hcut) {
//--------------------------------------------------------------------
   fname = Dir + fname;
   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "ERROR in "<< __func__
         << ": can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");

   TH1D* hst = new TH1D(hname.c_str(),
         ";Momentum of #pi^{#plus}, GeV/c"
         ";Events/0.005GeV/c",
         100,0.,0.5
         );
   hst->Sumw2();

   string dr = string("Ppls>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//--------------------------------------------------------------------
TH1D* fill_ppim(string fname, string hname, string hcut) {
//--------------------------------------------------------------------
   fname = Dir + fname;
   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "ERROR in "<< __func__
         << ": can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");

   TH1D* hst = new TH1D(hname.c_str(),
         ";Momentum of #pi^{#minus}, GeV/c"
         ";Events/0.005GeV/c",
         100,0.,0.5
         );
   hst->Sumw2();

   string dr = string("Pmns>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//--------------------------------------------------------------------
TH1D* fill_cpip(string fname, string hname, string hcut) {
//--------------------------------------------------------------------
   fname = Dir + fname;
   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "ERROR in "<< __func__
         << ": can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");

   TH1D* hst = new TH1D(hname.c_str(),
         ";cos #Theta(#pi^{#plus})"
         ";Events/0.02",
         100,-1.,1.
         );
   hst->Sumw2();

   string dr = string("Cpls>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//--------------------------------------------------------------------
TH1D* fill_cpim(string fname, string hname, string hcut) {
//--------------------------------------------------------------------
   fname = Dir + fname;
   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "ERROR in "<< __func__
         << ": can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");

   TH1D* hst = new TH1D(hname.c_str(),
         ";cos #Theta(#pi^{#minus})"
         ";Events/0.02",
         100,-1.,1.
         );
   hst->Sumw2();

   string dr = string("Cmns>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//--------------------------------------------------------------------
TH1D* fill_cosPM(string fname, string hname, string hcut) {
//--------------------------------------------------------------------
   fname = Dir + fname;
   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "ERROR in "<< __func__
         << ": can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");

   TH1D* hst = new TH1D(hname.c_str(),
         ";cos #Theta(#pi^{#plus}#pi^{#minus})"
         ";Events/0.01",
         200,-1.,1.
         );
   hst->Sumw2();

   string dr = string("cosPM>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//--------------------------------------------------------------------
TH1D* fill_invPM(string fname, string hname, string hcut) {
//--------------------------------------------------------------------
   fname = Dir + fname;
   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "ERROR in "<< __func__
         << ": can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");

   TH1D* hst = new TH1D(hname.c_str(),
         ";M^{inv}(#pi^{#plus}#pi^{#minus}), GeV/c^{2}"
         ";Events/0.005GeV/c^{2}",
         100,0.2,0.7
         );
   hst->Sumw2();

   string dr = string("invPM>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//--------------------------------------------------------------------
void fill_hist(int date, string var, FILL_H fill_h, TH1D* hst[]) {
//--------------------------------------------------------------------
#include "norm.h"

   string datafile( Form("data_%02ipsip_all.root",date%100) );
   string mcincfile( Form("mcinc_%02ipsip_all.root",date%100) );

   string hn0 = var + string(Form("_%02i",date%100));
   hst[0]=fill_h(datafile, hn0, "");
   hst[0]->SetMarkerStyle(20);
   hst[0]->SetMarkerSize(0.7);
   hst[0]->GetYaxis()->SetMaxDigits(3);
   cout << " Integral(data)= " << hst[0]->Integral() << endl;

   string hn1 = var + string("_3650");
   hst[1]=fill_h("data_3650_2021.root", hn1, "");
   hst[1]->Scale(C_Dat.at(date));
   hst[1]->SetLineColor(kBlue+1);
   hst[1]->SetLineWidth(2);
   cout << " Integral(3650)= " << hst[1]->Integral() << endl;

   string hn2 = var + string(Form("_%02iMC",date%100));
   hst[2]=fill_h(mcincfile, hn2, "");
   hst[2]->Scale(MC_Dat.at(date));
   hst[2]->SetLineColor(kRed+1);
   hst[2]->SetLineWidth(1);
   cout << " Integral(mc)= " << hst[2]->Integral() << endl;

   string hn3 = hn2 + "bg";
   hst[3]=fill_h(mcincfile, hn3, "dec!=64");
   hst[3]->Scale(MC_Dat.at(date));
   hst[3]->SetLineColor(kMagenta+1);
   hst[3]->SetLineWidth(2);
   cout << " Integral(mcbg)= " << hst[3]->Integral() << endl;
}

//--------------------------------------------------------------------
void plot_One(int date, string var, FILL_H fill_h) {
//--------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   string pdf = var+string("_")+to_string(date)+string(".pdf");
   c1->Print((pdf+"[").c_str()); // open pdf-file

   TH1D* hst[10];
   fill_hist(date,var,fill_h,hst);
   hst[2]->Add(hst[1]); // sum MC inclusive and continuum
   // double hst_max = max(hst[0]->GetMaximum(),hst[2]->GetMaximum());
   // double win_max = 1.15 * hst_max;

   c1->cd();
   gPad->SetGrid();
   int log = 1;
   if ( log == 1 ) {
      gPad->SetLogy();
      hst[0]->SetMinimum(1.);
   }

   SetHstFace(hst[0]);
   // hst[0]->SetAxisRange(0.,win_max,"Y");
   hst[0]->GetYaxis()->SetTitleOffset(1.2);
   hst[0]->Draw("E"); // data
   hst[2]->Draw("SAME,HIST"); // inc. MC
   hst[1]->Draw("SAME,HIST"); // cont.
   hst[3]->Draw("SAME,HIST"); // inc. MC not 64

   int pos = 0;
   TLegend* leg;
   if ( pos == 0 ) { // default
      leg = new TLegend(0.53,0.65,0.89,0.89);
   } else {
      leg = new TLegend(0.11,0.65,0.47,0.89);
   }
   leg->AddEntry(hst[0],
         (string("Data ")+to_string(date)).c_str(), "EP");
   leg->AddEntry(hst[2],
         Form("#color[%i]{MC + Continuum}",hst[2]->GetLineColor()),
         "L");
   leg->AddEntry(hst[1],
         Form("#color[%i]{Continuum}",hst[1]->GetLineColor()),
         "L");
   leg->AddEntry(hst[3],
         Form("#color[%i]{non #pi^{#plus}#pi^{#minus}J/#Psi decay}",
            hst[3]->GetLineColor() ), "L");
   leg->Draw();

   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
}

// {{{1 main()
//--------------------------------------------------------------------
void pipi_dataMC() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetStatFont(62);
   gStyle->SetLegendFont(42);
   // gStyle->SetLegendTextSize(0.03);

   //========================================================
   // set the name of the folder with the root files
   Dir = "prod_v709n3/";
   //========================================================

   for ( auto date : {2009, 2012, 2021} ) {
      plot_hist("cosPM",date);
      plot_hist("invPM",date);

      plot_hist("Pip_P",date);
      plot_hist("Pim_P",date);

      // plot_hist("Bnch",date); // Nch

      // slow
      // plot_One(date,"Nm",fill_Nm);
   }

   // ugly
   // plot_hist("Pip_C",2021); // or log see fun()
}
