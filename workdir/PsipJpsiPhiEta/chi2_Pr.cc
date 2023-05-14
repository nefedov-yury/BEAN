// plot chi2 of the kinematic fit
// 1) data vs signal MC in the window for M(K+K-) and M(gg)
//    -> chi2_sb_YEAR.pdf
//
// 2) data vs inclusive MC
//    -> chi2_YEAR.pdf

// {{{1 helper functions and constants
static double chi2_cut = 0.; // to keep value from cuts.h file

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
TH1D* get_chi2(string fname, string hname, int icp=1) {
//--------------------------------------------------------------------
#include "cuts.h"
   chi2_cut = chi2M;

   // name of folder with root files
   // static string dir("prod-12/");
   string dir("prod_v709/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   // 1) cut on recoil mass: [3.092, 3.102]
   TCut c_here = c_Mrec;

   // 2) cut on Mkk:
   if ( icp < 10 ) {
      c_here += c_phi; // Mkk in [2*Mk, 1.08GeV]
   } else {  // for chi2opt
      c_here += TCut( "Mkk>1.025&&Mkk<1.08" ); // [1.025,1.08 GeV]
   }

   // 3) central part and side-band of Mgg:
   if ( icp%10 == 1 ) { // CP
      c_here += c_cpgg;
   } else if ( icp%10 == 2 ) { // SB
      c_here += c_sbgg;
   }

   TH1D* chi2 = new TH1D(hname.c_str(),
         ";#chi^{2};Entries/2",
         100, 0., 200.);
   string select = string("ch2>>") + hname;
   a4c->Draw(select.c_str(), c_here, "goff");

   return chi2;
}

// {{{1 data vs inclusive MC
//--------------------------------------------------------------------
void chi2_Plot(int date) {
//--------------------------------------------------------------------
   double maxh = 0.;

   string datafile( Form("data_%02ipsip_all.root",date%100) );
   TH1D* chi2_data = get_chi2(datafile, "chi2_data");
   // double ndata = chi2_data->Integral();
   int nbin = chi2_data->GetMaximumBin();
   maxh = max(maxh, 1.2*chi2_data->GetMaximum());
   SetHstFace(chi2_data);
   chi2_data->GetYaxis()->SetMaxDigits(3);
   chi2_data->GetYaxis()->SetTitleOffset(1.2);

   // MC inclusive: arbitrary normalization
   string mcincfile( Form("mcinc_%02ipsip_all.root",date%100) );
   TH1D* chi2_mc = get_chi2(mcincfile, "chi2_mc");
   chi2_mc->SetLineColor(kRed);
   double nmc = chi2_mc->Integral();
   // chi2_mc->Scale(ndata/nmc);
   chi2_mc->Scale( chi2_data->GetBinContent(nbin)/
                   chi2_mc->GetBinContent(nbin) );
   maxh = max(maxh, chi2_mc->GetMaximum());
   chi2_mc->SetLineWidth(2);
   chi2_mc->SetLineColor(kRed);
   chi2_mc->SetLineStyle(kDashed);

   // MC phi eta: arbitrary normalization
   string mcsigfile( Form("mcsig_kkmc_%02i.root",date%100) );
   TH1D* chi2_mcS = get_chi2(mcsigfile, "chi2_mcS");
   // double nmcS = chi2_mcS->Integral();
   chi2_mcS->Scale( chi2_data->GetBinContent(nbin)/
                    chi2_mcS->GetBinContent(nbin) );
   maxh = max(maxh, chi2_mcS->GetMaximum());
   chi2_mcS->SetLineWidth(2);
   chi2_mcS->SetLineColor(kGreen+2);

   // MC KKeta: arbitrary  normalization
   // string mckketafile( Form("mckketa_kkmc_%02i.root",date%100) );
   // TH1D* chi2_mcB = get_chi2(mckketafile, "chi2_mcB");
   // chi2_mcB->Scale( chi2_data->GetBinContent(nbin)/
                    // chi2_mcB->GetBinContent(nbin) );
   // chi2_mcB->SetLineColor(kBlue+1);
   // chi2_mcB->SetLineStyle(kDashed);

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1->cd();
   gPad->SetGrid();
   chi2_data->DrawCopy("E1");
   chi2_mc->DrawCopy("SAME HIST");
   chi2_mcS->DrawCopy("SAME HIST");
   // chi2_mcB->DrawCopy("SAME HIST");

   // double xx = chi2_cut;
   // double yy = 0.25*maxh;
   // TLine* l1 = new TLine(xx,0.,xx,yy);
   // TLatex* t1 = new TLatex(xx+5,0.9*yy,
         // Form("cut: #chi^{2}<%.0f",xx));
   // l1->SetLineColor(kMagenta+2);
   // l1->SetLineWidth(3);
   // l1->Draw();
   // t1->SetTextColor(kMagenta+2);
   // t1->Draw();

   TLegend* leg = new TLegend(0.50,0.65,0.88,0.88);
   leg->AddEntry(chi2_data, Form("Data %i",date), "LE");
   leg->AddEntry(chi2_mc,"MC inclusive","L");
   leg->AddEntry(chi2_mcS,"MC signal #phi#eta","L");
   // leg->AddEntry(chi2_mcB,"MC KK#eta","L");
   leg->Draw();
   c1->Update();

   string pdf( Form("chi2_%i.pdf",date) );
   c1->Print(pdf.c_str());
}

// {{{1 optimization of chi2-cut
// plot Sig/sqrt(Sig+Bkg) ratio (type=0)
// or 1/(relative error) (type=1)
//--------------------------------------------------------------------
TGraph* optSigBkg(TH1D* data, TH1D* sb, int type = 0) {
//--------------------------------------------------------------------
   int n = data->GetXaxis()->GetNbins();
   if ( sb->GetXaxis()->GetNbins() != n ) {
      return nullptr;
   }

   vector<double> x(n);
   // data->GetXaxis()->GetCenter( x.data() );
   data->GetXaxis()->GetLowEdge( x.data() );
   double wx = data->GetBinWidth(1);
   for_each( begin(x),end(x),[wx](double& v){v+=wx;} );

   // data = sig + bkg
   const double* y = data->GetArray();
   vector<double> d(y+1,y+1+n); // exclude underflow/overflows bins
   partial_sum(begin(d), end(d), begin(d));

   // sideband = bkg
   const double* t = sb->GetArray();
   vector<double> b(t+1,t+1+n); // exclude underflow/overflows bins
   partial_sum(begin(b), end(b), begin(b));

   vector<double> SBrat(n);
   double rat = 0.;
   for ( int i = 0; i < n; i++ ) {
      if ( type == 0 ) { // default
         rat = (d[i]>0) ? (d[i]-b[i])/sqrt(d[i]) : 0;
      } else if (type == 1 ) {
         rat = (d[i]-b[i])/sqrt(d[i]+b[i]);
      }
      SBrat[i] = rat;
   }

   TGraph* gr = new TGraph( n, x.data(), SBrat.data() );
   gr->SetTitle(";#chi^{2};#it{Sig/#sqrt{Sig+Bkg}}");
   gr->GetYaxis()->SetMaxDigits(3);
   gr->GetYaxis()->SetTitleOffset(1.2);
   gr->SetMarkerColor(kViolet);
   gr->SetMarkerStyle(20);

   return gr;
}

//--------------------------------------------------------------------
void chi2opt(int date) {
//--------------------------------------------------------------------
   double maxh = 0.;

   string datafile( Form("data_%02ipsip_all.root",date%100) );
   TH1D* chi2cp = get_chi2(datafile, "chi2_cp", 11);
   SetHstFace(chi2cp);
   maxh = max(maxh, 1.2*chi2cp->GetMaximum());
   chi2cp->SetLineColor(kBlack);
   chi2cp->SetLineWidth(1);
   chi2cp->GetYaxis()->SetMaxDigits(3);
   chi2cp->GetYaxis()->SetTitleOffset(1.25);

   TH1D* chi2sb = get_chi2(datafile, "chi2_sb", 12);
   chi2sb->SetLineColor(kBlue+2);
   chi2sb->SetFillStyle(3001);
   chi2sb->SetFillColor(kBlue+1);

   // MC signal phi eta: arbitrary normalization
   string mcsigfile( Form("mcsig_kkmc_%02i.root",date%100) );
   TH1D* chi2_mcS = get_chi2(mcsigfile, "chi2_mcS", 11);
   chi2_mcS->SetLineColor(kGreen+2);
   chi2_mcS->SetLineWidth(1);
   chi2_mcS->Scale(
         (chi2cp->Integral() - chi2sb->Integral())
         / chi2_mcS->Integral()
         );
   maxh = max(maxh, chi2_mcS->GetMaximum());

   TGraph* gr = optSigBkg(chi2cp,chi2sb,0);
   gr->GetXaxis()->SetLimits(20.,200.);
   gr->GetYaxis()->SetNdivisions(505);
   if ( date == 2009 ) {
      gr->SetMinimum(10.);
   } else if ( date == 2012 ) {
      gr->SetMinimum(21.);
   } else if ( date == 2021 ) {
      gr->SetMinimum(50.);
   }

   TCanvas* c1 = new TCanvas("c1","...",0,0,1500,800);
   c1->Divide(2,1);

   c1->cd(1);
   gPad->SetGrid();

   chi2cp->SetMaximum(maxh);
   chi2cp->DrawCopy("E1");
   chi2sb->DrawCopy("SAME");
   chi2_mcS->DrawCopy("SAME HIST");

   TLegend* leg = new TLegend(0.46,0.64,0.89,0.89);
   leg->SetHeader( Form("#chi^{2}(5C): %i",date),"C");
   leg->AddEntry(chi2cp,"Data: signal region","LE");
   leg->AddEntry(chi2sb,"Data: side-band","F");
   leg->AddEntry(chi2_mcS,"MC signal #phi#eta","L");
   leg->Draw();

   gPad -> RedrawAxis();

   c1->cd(2);
   gPad->SetGrid();

   gr->Draw("APL");

   c1->Update();
   string pdf( Form("chi2opt_%i.pdf",date) );
   c1 -> Print(pdf.c_str());
}

// {{{1 data vs signal MC
//--------------------------------------------------------------------
void chi2_SB(int date) {
//--------------------------------------------------------------------
   double maxh = 0.;

   string datafile( Form("data_%02ipsip_all.root",date%100) );
   TH1D* chi2cp = get_chi2(datafile, "chi2_cp", 1);
   maxh = max(maxh, 1.2*chi2cp->GetMaximum());
   SetHstFace(chi2cp);
   chi2cp->SetLineColor(kBlack);
   chi2cp->SetLineWidth(2);
   chi2cp->GetYaxis()->SetMaxDigits(3);
   chi2cp->GetYaxis()->SetTitleOffset(1.2);

   TH1D* chi2sb = get_chi2(datafile, "chi2_sb", 2);
   chi2sb->SetLineColor(kBlue+2);
   chi2sb->SetFillStyle(3001);
   chi2sb->SetFillColor(kBlue+1);

   // MC signal phi eta: arbitrary normalization
   string mcsigfile( Form("mcsig_kkmc_%02i.root",date%100) );
   TH1D* chi2_mcS = get_chi2(mcsigfile, "chi2_mcS", 1);
   chi2_mcS->SetLineWidth(2);
   chi2_mcS->SetLineColor(kGreen+2);
   chi2_mcS -> Scale(
         (chi2cp->Integral() - chi2sb->Integral())
         / chi2_mcS->Integral()
         );
   maxh = max(maxh, chi2_mcS->GetMaximum());

   // box to show rejected region
   double ch2_cut = chi2_cut;
   double ch2_max = 200.;
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-10);
   box->SetLineWidth(2);

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   chi2cp->SetMaximum(maxh);
   chi2cp->DrawCopy("E1");
   chi2sb->DrawCopy("SAME");
   chi2_mcS->DrawCopy("SAME HIST");

   // double ymin = gPad->PixeltoY(0);
   double ymin = 0;
   box->DrawBox(ch2_cut,ymin,ch2_max,maxh);

   TLegend* leg = new TLegend(0.46,0.64,0.89,0.89);
   leg->SetHeader( Form("#chi^{2}(5C): %i",date),"C");
   leg->AddEntry(chi2cp,"Data: signal region","LE");
   leg->AddEntry(chi2_mcS,"MC signal #phi#eta","L");
   leg->AddEntry(chi2sb,"Data: side-band","F");
   leg->AddEntry(box, "Rejection area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();

   string pdf( Form("chi2_sb_%i.pdf",date) );
   c1->Print(pdf.c_str());
}

// {{{1 Main
//--------------------------------------------------------------------
void chi2_Pr() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(42);
   // gStyle->SetLegendTextSize(0.04);

   // 0) optimization of chi2-cut
   // chi2opt(2009);
   // chi2opt(2012);
   // chi2opt(2021);

   // 1) data vs signal MC in the window for M(K+K-) and M(gg) .. fig.11
   chi2_SB(2009);
   // chi2_SB(2012);
   // chi2_SB(2021);

   // 2) data vs inclusive MC
   // chi2_Plot(2009);
   // chi2_Plot(2012);
   // chi2_Plot(2021);
}

