// plot chi2 of 4C kinematic fit
// 1) data vs inclusive MC after 4C kinematic fit
//    -> chi2_YEAR.pdf
// 2) data vs signal MC in the window for M(K+K-) and M(gg)
//    -> chi2_sb_YEAR.pdf

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
TH1D* get_chi2(string fname, string hname, int icut=0) {
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

   TCut c_here = c_Mrec + c_phi;
   if ( icut == 1 ) { // CP
      c_here += c_cpgg;
   } else if ( icut == 2 ) { // SB
      c_here += c_sbgg;
   }

   TH1D* chi2 = new TH1D(hname.c_str(),
                         ";#chi^{2};Entries/2",
                         100, 0., 200.);
   string select = string("ch2>>") + hname;
   a4c->Draw(select.c_str(), c_here, "goff");

   return chi2;
}

//-------------------------------------------------------------------------
void chi2_Plot(int date) {
//-------------------------------------------------------------------------
   string dname = (date==2009) ? "09" : "12";
   double maxh = 0.;
//    cout << " DEBUG: " << __func__ << " dname= " << dname << endl;

   string fn1 = string("data_") + dname + string("psip_all.root");
   TH1D* chi2_data = get_chi2(fn1, "chi2_data");
//    double ndata = chi2_data->Integral();
   maxh = max(maxh, chi2_data->GetMaximum());

   // MC inclusive: arbitrary normalization
   string fn2 = string("mcinc_") + dname + string("psip_all.root");
   TH1D* chi2_mc = get_chi2(fn2, "chi2_mc");
   chi2_mc->SetLineColor(kRed);
//    double nmc = chi2_mc->Integral();
//    chi2_mc->Scale(ndata/nmc);
   int nbin = 10;
   chi2_mc->Scale( chi2_data->GetBinContent(nbin)/
                   chi2_mc->GetBinContent(nbin) );
   maxh = max(maxh, chi2_mc->GetMaximum());

   // MC phi eta: arbitrary normalization
   string fn3 = string("mcsig_kkmc_") + dname + string(".root");
   TH1D* chi2_mcS = get_chi2(fn3.c_str(), "chi2_mcS");
//    double nmcS = chi2_mcS->Integral();
   chi2_mcS->Scale( chi2_data->GetBinContent(nbin)/
                    chi2_mcS->GetBinContent(nbin) );
   chi2_mcS->SetLineColor(kGreen+2);
   chi2_mcS->SetLineStyle(kDashed);

   // MC KKeta: arbitrary  normalization
   string fn4 = string("mckketa_kkmc_") + dname + string(".root");
   TH1D* chi2_mcB = get_chi2(fn4.c_str(), "chi2_mcB");
   chi2_mcB->Scale( chi2_data->GetBinContent(nbin)/
                    chi2_mcB->GetBinContent(nbin) );
   chi2_mcB->SetLineColor(kBlue+1);
   chi2_mcB->SetLineStyle(kDashed);

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1->cd();
   gPad->SetGrid();
   chi2_data->SetMaximum(1.1*maxh);
   chi2_data->DrawCopy("E1");
   chi2_mc->DrawCopy("SAME HIST");
   chi2_mcS->DrawCopy("SAME HIST");
//    chi2_mcB->DrawCopy("SAME HIST");

   double xx = 80.; // ----> draw chi2_cut
   double yy = 0.25*maxh;
   TLine* l1 = new TLine(xx,0.,xx,yy);
   TLatex* t1 = new TLatex(xx+5,0.9*yy,
                           Form("cut: #chi^{2}<%.0f",xx));
   l1->SetLineColor(kMagenta+2);
   l1->SetLineWidth(3);
   l1->Draw();
   t1->SetTextColor(kMagenta+2);
   t1->Draw();

   TLegend* leg = new TLegend(0.50,0.65,0.88,0.88);
   leg->AddEntry(chi2_data,(string("Data ")+to_string(date)).c_str(), "LE");
   leg->AddEntry(chi2_mc,"MC inclusive","L");
   leg->AddEntry(chi2_mcS,"MC #phi#eta","L");
//    leg->AddEntry(chi2_mcB,"MC KK#eta","L");
   leg->Draw();
   c1->Update();

   string pdf = string("chi2_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void chi2_SB(int date) {
//-------------------------------------------------------------------------
   string dname = (date==2009) ? "09" : "12";
   double maxh = 0.;
//    cout << " DEBUG: " << __func__ << " dname= " << dname << endl;

   string fn1 = string("data_") + dname + string("psip_all.root");
   TH1D* chi2cp = get_chi2(fn1.c_str(), "chi2_cp", 1);
   chi2cp->SetLineColor(kBlack);
   TH1D* chi2sb = get_chi2(fn1.c_str(), "chi2_sb", 2);
   chi2sb->SetLineColor(kBlue+2);
   chi2sb->SetFillStyle(3001);
   chi2sb->SetFillColor(kBlue+1);
   maxh = max(maxh, chi2cp->GetMaximum());

   // MC phi eta: arbitrary normalization
   string fn3 = string("mcsig_kkmc_") + dname + string(".root");
   TH1D* chi2_mcS = get_chi2(fn3.c_str(), "chi2_mcS", 1);
   chi2_mcS->SetLineColor(kRed+1);

   chi2_mcS->Scale(
      (chi2cp->Integral()-chi2sb->Integral())
      / chi2_mcS->Integral()
   );

   maxh = max(maxh, chi2_mcS->GetMaximum());

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1->cd();
   gPad->SetGrid();

   SetHstFace(chi2cp);
//    chi2cp->GetYaxis()->SetMaxDigits(3);
   chi2cp->SetLineWidth(2);
   chi2cp->GetYaxis()->SetTitleOffset(1.25);

   chi2cp->SetMaximum(1.15*maxh);
   chi2cp->DrawCopy("E1");
   chi2sb->DrawCopy("SAME");
   chi2_mcS->DrawCopy("SAME HIST");

   double xx = 80.; // ----> draw chi2_cut
   double yy = 0.25*maxh;
   TLine* l1 = new TLine(xx,0.,xx,yy);
   TLatex* t1 = new TLatex(xx+5,0.9*yy,
                           Form("reject #chi^{2} > %.0f",xx));
   l1->SetLineColor(kMagenta+2);
   l1->SetLineWidth(3);
   l1->Draw();
   t1->SetTextColor(kMagenta+2);
   t1->Draw();

   TLegend* leg = new TLegend(0.46,0.64,0.89,0.89);
   leg->SetHeader(
      (string("#bf{#chi^{2}(5C): ")+to_string(date)+"}").c_str(),"C");
   leg->AddEntry(chi2cp,"Data: signal region","LE");
   leg->AddEntry(chi2sb,"Data: side-band","F");
   leg->AddEntry(chi2_mcS,"MC signal #phi#eta","L");
   leg->Draw();
   c1->Update();

   string pdf = string("chi2_sb_") + to_string(date) + string(".pdf");
   c1->Print(pdf.c_str());
}

//-------------------------------------------------------------------------
void chi2_Pr()
//-------------------------------------------------------------------------
{
   gROOT->Reset();

   gStyle->SetOptStat(0);
//    gStyle->SetLegendFont(62); // not use it!
   gStyle->SetLegendTextSize(0.04);

//    chi2_Plot(2009);
//    chi2_Plot(2012);

   chi2_SB(2009);
//    chi2_SB(2012);
}

