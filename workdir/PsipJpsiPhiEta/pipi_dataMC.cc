// plot data distributions vs MC
// for pi+ pi- in selected Psi(2S) -> Jpsi pi+pi-
//      -> VarName_YEAR.pdf

// pointer on function wich fill one histo
typedef TH1D* (*FILL_H)(string, string, string);

// global name of folder with root files
static string dir("prod-11/");

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
TH1D* fill_nch(string fname, string hname, string hcut) {
//-------------------------------------------------------------------------
   fname = dir + fname;
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
                        ";Number of charged tracks"
                        ";Events",
                        20,0.5,20.5
                       );
   hst->Sumw2();

   string dr = string("nch>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//-------------------------------------------------------------------------
TH1D* fill_Nm(string fname, string hname, string hcut) {
//-------------------------------------------------------------------------
   fname = dir + fname;
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
         ";Number of M^{rec}_{#pi^{#plus}#pi^{#minus}} combinations"
         ";Events",
         20,0.5,20.5
//       30,0.5,30.5
         );
   hst->Sumw2();

   string dr = string("Nm+(Mrs>0)>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//-------------------------------------------------------------------------
TH1D* fill_ppip(string fname, string hname, string hcut) {
//-------------------------------------------------------------------------
   fname = dir + fname;
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
                        ";Momentum of #pi^{+}, GeV/c"
                        ";Events/0.005GeV/c",
                        100,0.,0.5
                       );
   hst->Sumw2();

   string dr = string("Ppls>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//-------------------------------------------------------------------------
TH1D* fill_ppim(string fname, string hname, string hcut) {
//-------------------------------------------------------------------------
   fname = dir + fname;
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
                        ";Momentum of #pi^{-}, GeV/c"
                        ";Events/0.005GeV/c",
                        100,0.,0.5
                       );
   hst->Sumw2();

   string dr = string("Pmns>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//-------------------------------------------------------------------------
TH1D* fill_cpip(string fname, string hname, string hcut) {
//-------------------------------------------------------------------------
   fname = dir + fname;
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
                        ";cos #Theta(#pi^{+})"
                        ";Events/0.02",
                        100,-1.,1.
                       );
   hst->Sumw2();

   string dr = string("Cpls>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//-------------------------------------------------------------------------
TH1D* fill_cpim(string fname, string hname, string hcut) {
//-------------------------------------------------------------------------
   fname = dir + fname;
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
                        ";cos #Theta(#pi^{-})"
                        ";Events/0.02",
                        100,-1.,1.
                       );
   hst->Sumw2();

   string dr = string("Cmns>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//-------------------------------------------------------------------------
TH1D* fill_cosPM(string fname, string hname, string hcut) {
//-------------------------------------------------------------------------
   fname = dir + fname;
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
                        ";cos #Theta(#pi^{+}#pi^{-})"
                        ";Events/0.01",
                        200,-1.,1.
                       );
   hst->Sumw2();

   string dr = string("cosPM>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//-------------------------------------------------------------------------
TH1D* fill_invPM(string fname, string hname, string hcut) {
//-------------------------------------------------------------------------
   fname = dir + fname;
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
                        ";M^{inv}(#pi^{+}#pi^{-}), GeV/c^{2}"
                        ";Events/0.005GeV/c^{2}",
                        100,0.2,0.7
                       );
   hst->Sumw2();

   string dr = string("invPM>>") + hname;
   nt1->Draw(dr.c_str(),hcut.c_str(),"goff");

   return hst;
}

//-------------------------------------------------------------------------
void fill_hist2009(FILL_H fill_h, string var, TH1D* hst[]) {
//-------------------------------------------------------------------------
#include "norm.h"

   string hn0 = var + string("09");
   hst[0]=fill_h("data_09psip_all.root", hn0, "");
   hst[0]->SetMarkerStyle(20);
   hst[0]->SetMarkerSize(0.7);
   hst[0]->GetYaxis()->SetMaxDigits(3);

   string hn1 = var + string("3650A");
   hst[1]=fill_h("data_3650_all.root", hn1, "");
   hst[1]->Scale(Cto09);
   hst[1]->SetLineColor(kBlue+1);
   hst[1]->SetLineWidth(2);

   string hn2 = var + string("09MC");
   hst[2]=fill_h("mcinc_09psip_all.root", hn2, "");
   hst[2]->Scale(Ito09);
   hst[2]->SetLineColor(kRed+1);
   hst[2]->SetLineWidth(1);

   string hn3 = var + string("09MCbg");
   hst[3]=fill_h("mcinc_09psip_all.root", hn3, "dec!=64");
   hst[3]->Scale(Ito09);
   hst[3]->SetLineColor(kMagenta+1);
   hst[3]->SetLineWidth(2);
}

//-------------------------------------------------------------------------
void fill_hist2012(FILL_H fill_h, string var, TH1D* hst[]) {
//-------------------------------------------------------------------------
#include "norm.h"

   string hn0 = var + string("12");
   hst[0]=fill_h("data_12psip_all.root", hn0, "");
   hst[0]->SetMarkerStyle(20);
   hst[0]->SetMarkerSize(0.7);
   hst[0]->GetYaxis()->SetMaxDigits(3);

   string hn1 = var + string("3650B");
   hst[1]=fill_h("data_3650_all.root", hn1, "");
   hst[1]->Scale(Cto12);
   hst[1]->SetLineColor(kBlue+1);
   hst[1]->SetLineWidth(2);

   string hn2 = var + string("12MC");
   hst[2]=fill_h("mcinc_12psip_all.root", hn2, "");
   hst[2]->Scale(Ito12);
   hst[2]->SetLineColor(kRed+1);
   hst[2]->SetLineWidth(1);

   string hn3 = var + string("12MCbg");
   hst[3]=fill_h("mcinc_12psip_all.root", hn3, "dec!=64");
   hst[3]->Scale(Ito12);
   hst[3]->SetLineColor(kMagenta+1);
   hst[3]->SetLineWidth(2);
}

//-------------------------------------------------------------------------
void plot_One(int date, string var, FILL_H fill_h, int pos=0, int log=0) {
//-------------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   string pdf = var+string("_")+to_string(date)+string(".pdf");
   c1->Print((pdf+"[").c_str()); // open pdf-file

   TH1D* hst[10];
   if( date==2009 ) {
      fill_hist2009(fill_h,var,hst);
   } else if(date==2012) {
      fill_hist2012(fill_h,var,hst);
   }
   hst[2]->Add(hst[1]);

   c1->cd();
   gPad->SetGrid();
   if ( log == 1 ) {
      gPad->SetLogy();
//       hst[0]->SetMinimum(1.);
   }
   SetHstFace(hst[0]);
   hst[0]->GetYaxis()->SetTitleOffset(1.2);
   hst[0]->Draw("E"); // data
   hst[2]->Draw("SAME,HIST"); // inc. MC
   hst[1]->Draw("SAME,HIST"); // cont.
   hst[3]->Draw("SAME,HIST"); // inc. MC not 64

   TLegend* leg;
   if ( pos == 0 ) { // default
     leg = new TLegend(0.53,0.65,0.89,0.89);
   } else {
     leg = new TLegend(0.11,0.65,0.47,0.89);
   }
   leg->AddEntry(hst[0], (string("Data ")+to_string(date)).c_str(), "EP");
   leg->AddEntry(hst[2], Form("#color[%i]{MC + Continuum}",
                              hst[2]->GetLineColor() ),"L");
   leg->AddEntry(hst[1], Form("#color[%i]{Continuum}",
                              hst[1]->GetLineColor() ),"L");
   leg->AddEntry(hst[3],
         Form("#color[%i]{non #pi^{#plus}#pi^{#minus}J/#Psi decay}",
         hst[3]->GetLineColor() ),"L");
   leg->Draw();

   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
}

//-------------------------------------------------------------------------
void plot_Ppi(int date) {
//-------------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   string pdf = string("Ppi_")+to_string(date)+string(".pdf");
   c1->Print((pdf+"[").c_str()); // open pdf-file

   TH1D* hst[10];
   if( date==2009 ) {
      fill_hist2009(fill_ppip,"ppls",&hst[0]);
      hst[0]->SetMaximum(1e6);
   } else if(date==2012) {
      fill_hist2012(fill_ppip,"ppls",&hst[0]);
      hst[0]->SetMaximum(3.2e6);
   }
   hst[2]->Add(hst[1]);

   hst[0]->Draw("E"); // data
   hst[2]->Draw("SAME,HIST"); // inc. MC
   hst[1]->Draw("SAME,HIST"); // cont.
   hst[3]->Draw("SAME,HIST"); // inc. MC not 64

   TLegend* leg = new TLegend(0.62,0.75,0.94,0.94);
   leg->AddEntry(hst[0], (string("Data ")+to_string(date)).c_str(), "EP");
   leg->AddEntry(hst[2], Form("#color[%i]{MC + Continuum}",
                              hst[2]->GetLineColor() ),"L");
   leg->AddEntry(hst[1], Form("#color[%i]{Continuum}",
                              hst[1]->GetLineColor() ),"L");
   leg->AddEntry(hst[3], Form("#color[%i]{non #pi^{+}#pi^{-}J/#Psi decay}",
                              hst[3]->GetLineColor() ),"L");
   leg->Draw();

   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   if( date==2009 ) {
      fill_hist2009(fill_ppim,"pmns",&hst[4]); // Pi-
      hst[4]->SetMaximum(1e6);
   } else if(date==2012) {
      fill_hist2012(fill_ppim,"pmns",&hst[4]); // Pi-
      hst[4]->SetMaximum(3.2e6);
   }
   hst[6]->Add(hst[5]);

   hst[4]->Draw("E"); // data
   hst[6]->Draw("SAME,HIST"); // inc. MC
   hst[5]->Draw("SAME,HIST"); // cont.
   hst[7]->Draw("SAME,HIST"); // inc. MC not 64
   leg->Draw();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
}

//-------------------------------------------------------------------------
void plot_Cpi(int date) {
//-------------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   string pdf = string("Cpi_")+to_string(date)+string(".pdf");
   c1->Print((pdf+"[").c_str()); // open pdf-file

   TH1D* hst[10];
   if( date==2009 ) {
      fill_hist2009(fill_cpip,"cpls",&hst[0]);
   } else if(date==2012) {
      fill_hist2012(fill_cpip,"cpls",&hst[0]);
   }
   hst[2]->Add(hst[1]);

   hst[0]->Draw("E"); // data
   hst[2]->Draw("SAME,HIST"); // inc. MC
   hst[1]->Draw("SAME,HIST"); // cont.
   hst[3]->Draw("SAME,HIST"); // inc. MC not 64

   TLegend* leg = new TLegend(0.34,0.50,0.66,0.70);
   leg->AddEntry(hst[0], (string("Data ")+to_string(date)).c_str(), "EP");
   leg->AddEntry(hst[2], Form("#color[%i]{MC + Continuum}",
                              hst[2]->GetLineColor() ),"L");
   leg->AddEntry(hst[1], Form("#color[%i]{Continuum}",
                              hst[1]->GetLineColor() ),"L");
   leg->AddEntry(hst[3], Form("#color[%i]{non #pi^{+}#pi^{-}J/#Psi decay}",
                              hst[3]->GetLineColor() ),"L");
   leg->Draw();

   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   if( date==2009 ) {
      fill_hist2009(fill_cpim,"cmns",&hst[4]); // Pi-
   } else if(date==2012) {
      fill_hist2012(fill_cpim,"cmns",&hst[4]); // Pi-
   }
   hst[6]->Add(hst[5]);

   hst[4]->Draw("E"); // data
   hst[6]->Draw("SAME,HIST"); // inc. MC
   hst[5]->Draw("SAME,HIST"); // cont.
   hst[7]->Draw("SAME,HIST"); // inc. MC not 64
   leg->Draw();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
}

//-------------------------------------------------------------------------
TH1D* get_hst(string fname, string hname) {
//-------------------------------------------------------------------------
   fname = dir + fname;
   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "ERROR in "<< __func__
           << ": can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TH1D* hst = (TH1D*)gROOT->FindObject(hname.c_str());
   if ( !hst ) {
      cout << " can not find " << hname << endl;
      exit(0);
   }

   return hst;
}

//-------------------------------------------------------------------------
void get_hist2009(string hname, TH1D* hst[]) {
//-------------------------------------------------------------------------
#include "norm.h"

   hst[0]=get_hst("data_09psip_all.root", hname);
   hst[0]->SetMarkerStyle(20);
   hst[0]->SetMarkerSize(0.7);
   hst[0]->GetYaxis()->SetMaxDigits(3);

   hst[1]=get_hst("data_3650_all.root", hname);
   hst[1]->Scale(Cto09);
   hst[1]->SetLineColor(kBlue+1);
   hst[1]->SetLineWidth(2);

   hst[2]=get_hst("mcinc_09psip_all.root", hname);
   hst[2]->Scale(Ito09);
   hst[2]->SetLineColor(kRed+1);
   hst[2]->SetLineWidth(1);
}

//-------------------------------------------------------------------------
void get_hist2012(string hname, TH1D* hst[]) {
//-------------------------------------------------------------------------
#include "norm.h"

   hst[0]=get_hst("data_12psip_all.root", hname);
   hst[0]->SetMarkerStyle(20);
   hst[0]->SetMarkerSize(0.7);
   hst[0]->GetYaxis()->SetMaxDigits(3);

   hst[1]=get_hst("data_3650_all.root", hname);
   hst[1]->Scale(Cto12);
   hst[1]->SetLineColor(kBlue+1);
   hst[1]->SetLineWidth(2);

   hst[2]=get_hst("mcinc_12psip_all.root", hname);
   hst[2]->Scale(Ito12);
   hst[2]->SetLineColor(kRed+1);
   hst[2]->SetLineWidth(1);
}

//-------------------------------------------------------------------------
void plot_hist(string hname, int date) {
//-------------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   TH1D* hst[10];
   if ( date==2009 ) {
      get_hist2009(hname,hst);
   } else if ( date==2012 ) {
      get_hist2012(hname,hst);
   }
   hst[3]=(TH1D*)hst[2]->Clone("Sum");
   hst[3]->Add(hst[1]); // sum MC inclusive and continuum

   // relative normalisation
//    hst[3]->Scale( hst[0]->Integral()/hst[3]->Integral() );

   // box to show cut
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-9);
   box->SetLineWidth(2);

   c1->cd();
   gPad->SetGrid();

   SetHstFace(hst[0]);

   double hst_max = hst[3]->GetMaximum();
   int pos = 0;
   if ( hname == "cosPM" ) {
      hst[0]->SetTitle( ";cos #Theta_{#pi^{#plus}#pi^{#minus}}"
                        ";Events/0.02" );
      if ( date == 2009 ) {
         hst[0]->GetYaxis()->SetTitleOffset(1.2);
      }
   } else if ( hname == "invPM" ) {
      hst[0]->SetTitle( ";M^{ inv}_{#pi^{#plus}#pi^{#minus}}, GeV/c^{2}"
                        ";Events/0.005 GeV/c^{2}" );
      hst[0]->SetAxisRange(0.,1.15*hst_max,"Y");
      if ( date == 2009 ) {
         hst[0]->GetYaxis()->SetTitleOffset(1.2);
      }
   } else if ( hname == "Pip_P" ) {
      pos = 1;
      box = nullptr;
      hst[0]->SetTitle( ";Momentum of #pi^{#plus}, GeV/c"
                        ";Events/0.005 GeV/c" );
      hst[0]->SetAxisRange(0.,1.15*hst_max,"Y");
      hst[0]->GetYaxis()->SetTitleOffset(1.2);
   } else if ( hname == "Pim_P" ) {
      pos = 1;
      box = nullptr;
      hst[0]->SetTitle( ";Momentum of #pi^{#minus}, GeV/c"
                        ";Events/0.005 GeV/c" );
      hst[0]->SetAxisRange(0.,1.15*hst_max,"Y");
      hst[0]->GetYaxis()->SetTitleOffset(1.2);
   } else if ( hname == "Pip_C") {
      pos = 2;
      hst[0]->SetTitle( ";cos #Theta_{#pi^{#plus}}"
                        ";Events/0.02" );
   } else if ( hname == "Pim_C") {
      pos = 2;
      hst[0]->SetTitle( ";cos #Theta_{#pi^{#minus}}"
                        ";Events/0.02" );
   }

   hst[0]->Draw("E"); // data
   if ( hname == "cosPM" ) {
      box->DrawBox(0.8,.0,1.0,hst_max);
   } else if ( hname == "invPM" ) {
      double mk0    = 0.497611;
      double d = 0.008;
      box->DrawBox(mk0-d,.0,mk0+d,hst_max);
   } else if ( hname == "Pip_C" || hname == "Pim_C" ) {
      box->DrawBox(0.8,.0,1.0,hst_max);
      box->DrawBox(-1.0,.0,-0.8,hst_max);
      gPad->SetLogy(true);
   }
   hst[0]->Draw("E, SAME"); // data
   hst[3]->Draw("SAME,HIST"); // sum
   hst[1]->Draw("SAME,HIST"); // cont.
   gPad->RedrawAxis();

   TLegend* leg;
   if ( pos == 0 ) { // default
      leg = new TLegend(0.53,0.65,0.89,0.89);
   } else if ( pos == 1 ) {
      leg = new TLegend(0.60,0.75,0.92,0.92);
   } else if ( pos == 2 ) {
      leg = new TLegend(0.32,0.26,0.68,0.50);
   }
   leg->AddEntry(hst[0], (string("Data ")+to_string(date)).c_str(), "EP");
   leg->AddEntry(hst[3], Form("#color[%i]{MC + Continuum}",
                              hst[2]->GetLineColor() ),"L");
   leg->AddEntry(hst[1], Form("#color[%i]{Continuum}",
                              hst[1]->GetLineColor() ),"L");
   if ( box ) {
      leg->AddEntry(box, "Rejection area","F");
   }
   leg->Draw();

   c1->Update();

   string pdf = hname + "_" + to_string(date)+".pdf";
   c1->Print( pdf.c_str() );
}

//-------------------------------------------------------------------------
void pipi_dataMC() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetStatFont(62);
   gStyle->SetLegendFont(42);
//    gStyle->SetLegendTextSize(0.03);

//========================================================

//    plot_hist("cosPM",2009);
//    plot_hist("cosPM",2012);
//    plot_hist("invPM",2009);
//    plot_hist("invPM",2012);

//    plot_hist("Pip_P",2009);
//    plot_hist("Pim_P",2009);
//    plot_hist("Pip_P",2012);
//    plot_hist("Pim_P",2012);

//    plot_One(2009,"Nm",fill_Nm,0,1);
//    plot_One(2009,"Nch",fill_nch,0,1);
//    plot_One(2012,"Nm",fill_Nm,0,1);
//    plot_One(2012,"Nch",fill_nch,0,1);

// OLD:
//    plot_hist("Pip_C",2012); // log; lin - ???

//    plot_Ppi(2009); // ?
//    plot_Ppi(2012); // ?
//    plot_Cpi(2009); // ?
//    plot_Cpi(2012); // ?

}
