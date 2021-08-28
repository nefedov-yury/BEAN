// plot M(K+K-) distributions
// after 4C kinematic fit and chi^2 cut (see cuts.h)
// -> mass_phi.pdf

// #include<TH1.h>
// #include<TF1.h>
// #include<TCanvas.h>

#include <Math/GSLIntegrator.h>

// #include<TRandom.h>

// #include <iostream>
// #include <cmath>

//----------------------------------------------------------------------
double BWGauss_integrand(double u, void* params) {
//----------------------------------------------------------------------
// integrand for BWGauss -> u = (X'-X)/(sqrt2*sigma)
//----------------------------------------------------------------------
  double* p = (double *) params;
  double a = u + p[0];
  double b = p[1];
  double f = exp(-u*u)/(a*a+b*b);
  return f;
}

//----------------------------------------------------------------------
double BWGauss(double x, void* params) {
//----------------------------------------------------------------------
// Breit-Wigner distribution folded with a normal distribution
//
//                               [  1           G/2
//   Integral ( -oo < X' < +oo ) [  -- * -------------------- *
//                               [  pi   (G/2)^2 + (X' - M)^2
//
//                       1                   (X'-X)^2    ]
//              * ----------------  * exp( - --------- ) ] dX'
//                sqrt(2*pi)*sigma           2*sigma^2   ]
//
// Use GSL-QAGI adaptive integration on infinite intervals
//----------------------------------------------------------------------
  static const double Norm = 1./(M_PI*sqrt(M_PI));

  // get parameters:
  double* p    = (double *) params;
  double M     = p[0];
  double G     = p[1]/2;
  double sigma = p[2];

  double tmp  = 1./(M_SQRT2*sigma);
  double ab[] = {(x-M)*tmp, G*tmp};

  // desired errors:                abs    rel
  ROOT::Math::GSLIntegrator gsl_int(1.e-9, 1.e-6, 1000);
  double result = gsl_int.Integral(BWGauss_integrand,ab); // used qagi
  // double error = gsl_int.Error();
  // int status = gsl_int.Status();

  return Norm*G*tmp*tmp * result;
}

//----------------------------------------------------------------------
double BWG_fit_func(double* x, double* p) {
//----------------------------------------------------------------------
// function for fitting bu Breit-Wigner distribution
// folded with a normal distribution

  return p[0]*BWGauss(x[0],&p[1]);
}

//----------------------------------------------------------------------
void BWGaussFit(TH1D* hist) {
//----------------------------------------------------------------------
  const double Mphi = 1.019460; //1019.460 +/- 0.016 MeV
  const double Gphi = 4.247e-3; //   4.247 +/- 0.016 MeV

  TF1* bwg = new TF1("bwg", BWG_fit_func, Mphi-0.04, Mphi+0.04, 4);
  bwg->SetLineWidth(2);
  bwg->SetLineColor(kRed);
  bwg->SetParNames("Norm","Mphi","Gphi","Sigma");
  bwg->SetParameters(1., Mphi, Gphi, 1.7e-3);

  bwg->SetRange(Mphi-5*Gphi, Mphi+5*Gphi);
//   bwg->FixParameter(1, Mphi);    // Mass
  bwg->FixParameter(2, Gphi);  // Gamma

  hist->Fit("bwg","R");
//   hist->DrawCopy("E1");
  hist->DrawCopy();
}

//-------------------------------------------------------------------------
void get_mass_hist(string fname, string hname, TH1D* mkk[])
//-------------------------------------------------------------------------
{
#include "cuts.h"

  TFile* froot = TFile::Open(fname.c_str(),"READ");
  if( froot == 0 ) {
    cerr << "can not open " << fname << endl;
    exit(0);
  }

  froot->cd("PsipJpsiPhiEta");
  TTree *a4c = (TTree*)gDirectory->Get("a4c");

  // central part of Mgg
  mkk[0] = new TH1D(hname.c_str(),";M(K^{+}K^{-})",
                    100, mphi-0.04, mphi+0.04);
  string select = string("Mkk>>") + hname;
  a4c->Draw(select.c_str(), c_Mrec+c_chi2+c_cpgg, "goff");

  // side band of Mgg
  string hnsb = hname+string("_sb");
  mkk[1] = new TH1D(hnsb.c_str(),";M(K^{+}K^{-})",
                     100, mphi-0.04,mphi+0.04);
  string select2 = string("Mkk>>") + hnsb;
  a4c->Draw(select2.c_str(), c_Mrec+c_chi2+c_sbgg, "goff");

  mkk[0]->SetLineWidth(2);
  mkk[0]->GetXaxis()->SetTitleSize(0.05);
  mkk[0]->GetXaxis()->SetLabelFont(62);
  mkk[0]->GetXaxis()->SetLabelSize(0.04);
  mkk[0]->GetYaxis()->SetTitleSize(0.05);
  mkk[0]->GetYaxis()->SetLabelFont(62);
  mkk[0]->GetYaxis()->SetLabelSize(0.04);

  mkk[1]->SetLineColor(kGreen+3);
  mkk[1]->SetLineWidth(2);
}

//-------------------------------------------------------------------------
void plot_mass_phi()
//-------------------------------------------------------------------------
{
  gStyle->SetOptStat(0);
//   gStyle->SetOptStat(10);
  gStyle->SetOptFit(112); // print all parameters (fixed)
  gStyle->SetStatFont(62);
  gStyle->SetLegendTextSize(0.08);
  gStyle->SetStatH(0.25);
  gStyle->SetStatW(0.23);

#include "cuts.h"

  vector<string> fnames = {
        "data_09psip_all.root",
        "data_12psip_all.root",
        "mcinc_09psip_all.root",
        "mcinc_12psip_all.root",
        "ntpl_kkmc_09.root",
        "ntpl_kkmc_12.root"
  };

  vector<string> date = {
        "2009",
        "2012"
  };
  vector<string> titles = {
        "date",
        "MC inclusive",
        "MC signal",
  };

  TH1D* mkk[20];
  for (int i = 0; i < 6; i++) {
      string hn = string("mkk_")+to_string(i);
      get_mass_hist(fnames[i].c_str(), hn, &mkk[2*i]);
      if ( i < 2 ) mkk[2*i]->SetOption("E");
  }
  //-----------------------------------------------------------------------

  TLatex* t1 = new TLatex;
  t1->SetTextSize(0.08);
  t1->SetTextColor(kMagenta+3);

  TLine* lR = new TLine;
  lR->SetLineColor(kRed+2);
  lR->SetLineWidth(2);
  lR->SetLineStyle(7);

  TLine* lin = new TLine;
  lin->SetLineColor(kRed);
  TLegend* leg = new TLegend(0.56,0.6,0.975,0.93);
  leg->AddEntry(mkk[0],"data","LE");
  leg->AddEntry(mkk[1]," side-band","L");
  leg->AddEntry(lin,"BW+Gauss fit","L");

  TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
  int nx = 2;
  int ny = 3;
  c1->Divide(nx,ny);

  int i = 0;
  for ( int y = 0; y < ny; ++y ) {
     for ( int x = 0; x < nx; ++x ) {
        c1->cd(i+1);
//         gPad->SetGrid();
        BWGaussFit(mkk[2*i]);
        mkk[2*i+1]->DrawCopy("SAME");

        double ymax=0.4*mkk[2*i]->GetMaximum();
        lR->DrawLine(mphi-wphi,0,mphi-wphi,ymax);
        lR->DrawLine(mphi+wphi,0,mphi+wphi,ymax);

        double yy=0.9*mkk[2*i]->GetMaximum();
        t1->DrawLatex(0.985,yy,date[x].c_str());
        t1->DrawLatex(0.985,0.8*yy,titles[y].c_str());

        i = i+1;
     }
  }

  leg->Draw();

  c1->Update();
  c1->Print("mass_phi.pdf");
}

//-------------------------------------------------------------------------
void mass_phi()
//-------------------------------------------------------------------------
{
  gROOT->Reset();
  plot_mass_phi();
}
