// MrecFitSB.cc
// * MrecDraw() : draws a picture of recoil mass of pi+pi-
//    -> Mrec{YEAR}[_zoom].pdf
// * DoFitSB() : calculate number of Psi(2S) -> J/Psi pi+ pi- decays
//               in data using side-band method
//    -> Mrec{YEAR}_fsb[_sys{Sys}]_T{Npol}.pdf
//
// We fit the distribution of recoil mass (Mrec) of two pions in the
// area far from the peak of J/Psi by rescaling inclusive MC.
// We rescale only background from non pi+ pi- J/Psi by a polynomial,
// and correct signal only by «normalization» constant.
// The contribution of the continuum is taken into account.
// Number of signal events is calculated by subtracting the estimated
// background from the data.

#include <filesystem>
// namespace fs = std::filesystem;

#include "RewTrkPiK.hpp"    // RewTrkPi(), RewTrk_K() functions

// {{{1 helper functions and constants
//--------------------------------------------------------------------
// GLOBAL: name of folder with root files
string Dir;

//--------------------------------------------------------------------
constexpr double SQ(double x)
//--------------------------------------------------------------------
{
   return x*x;
}

//--------------------------------------------------------------------
void SetHstFace(TH1* hst)
//--------------------------------------------------------------------
{
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
void set_draw_opt(vector<TH1D*>& hst) {
//--------------------------------------------------------------------
   // data
   hst[0]->SetMarkerStyle(20);
   hst[0]->SetMarkerSize(0.7);
   // data 3650
   hst[1]->SetLineColor(kBlue+1);
   hst[1]->SetLineWidth(2);
   // MC signal
   hst[2]->SetLineColor(kGreen+3);
   hst[2]->SetLineWidth(2);
   // MC bg from pi+pi-J/Psi
   hst[3]->SetLineColor(kCyan+2);
   hst[3]->SetLineWidth(2);
   // MC bg not pi+pi-J/Psi
   hst[4]->SetLineColor(kMagenta+1);
   hst[4]->SetLineWidth(2);
}

// {{{1 Fill histograms
//--------------------------------------------------------------------
vector<TH1D*> fill_mrec(string fname, string hname,
      int date, bool isMC=false, bool sys_NoHC=false)
//--------------------------------------------------------------------
{
   if ( sys_NoHC && isMC ) { // sys
      fname = "NoHC/" + fname;
   }

   fname = Dir + fname;
   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if ( froot == 0 ) {
      cerr << "ERROR in "<< __func__
         << ": can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");
   if ( !nt1 ) {
      cerr << "ERROR in "<< __func__
         << " can not find nt1" << endl;
      exit(EXIT_FAILURE);
   }

   //Declaration of leaves types
   //        those not used here are commented out
   vector<float> * Mrec = nullptr; // must be !
   Float_t         Mrs;
   Float_t         Ptsp;
   Float_t         Ptsm;
   // Float_t         Mrb;
   // Float_t         Ptp;
   // Float_t         Cpls;
   // Float_t         Ptm;
   // Float_t         Cmns;
   Int_t           dec;
   // Float_t         mcmkk;

   // Set branch addresses.
   nt1->SetBranchAddress("Mrec",&Mrec);
   nt1->SetBranchAddress("Mrs",&Mrs);
   nt1->SetBranchAddress("Ptsp",&Ptsp);
   nt1->SetBranchAddress("Ptsm",&Ptsm);
   // nt1->SetBranchAddress("Mrb",&Mrb);
   // nt1->SetBranchAddress("Ptp",&Ptp);
   // nt1->SetBranchAddress("Cpls",&Cpls);
   // nt1->SetBranchAddress("Ptm",&Ptm);
   // nt1->SetBranchAddress("Cmns",&Cmns);
   nt1->SetBranchAddress("dec",&dec);
   // nt1->SetBranchAddress("mcmkk",&mcmkk);

   TH1D* hst = new TH1D(hname.c_str(),
         ";M^{rec}_{#pi^{#plus}#pi^{#minus }}, GeV/c^{2}"
         ";Entries/1 MeV/c^{2}",
         200,3.0,3.2 // 1bin = 1MeV
         // 400,3.0,3.2 // 1bin = 0.5MeV
         );
   hst->Sumw2(true);

   Long64_t nentries = nt1->GetEntries();

   if ( !isMC ) { // Data: return all recoil masses
      for ( Long64_t i = 0; i < nentries; ++i ) {
         nt1->GetEntry(i);
         for (const auto& mr : *Mrec ) {
            hst->Fill(mr);
         }
      }
      return vector<TH1D*> { hst };
   }

   // MC return:
   // 0 recoil mass of true pi+pi- pair (MC)
   // 1 background for dec==64 (pi+pi-J/Psi)
   // 2 background for dec!=64 (other Psi' decays)

   vector<TH1D*> hist(3,nullptr);
   vector<string> hn = { hname+"sig", hname+"bg1", hname+"bg2" };
   for ( size_t i = 0; i < hist.size(); ++i ) {
      hist[i] = (TH1D*)hst->Clone(hn[i].c_str());
   }
   delete hst;
   hst = nullptr;

   for ( Long64_t i = 0; i < nentries; ++i ) {
      nt1->GetEntry(i);
      // cerr<<" size Mrec="<<Mrec->size()<<" dec="<<dec<<endl;
      // cerr<<" Mrs= "<<Mrs<<" Pt="<<Ptsp<<","<<Ptsm<<endl;
      if ( dec == 64 ) {
         if ( Mrs > 1 ) {
            double wp = RewTrkPi( date, Ptsp, +1.);
            double wm = RewTrkPi( date, Ptsm, -1.);
            hist[0]->Fill(Mrs,wp*wm);
         }
         for (const auto& mr : *Mrec ) {
            hist[1]->Fill(mr);
         }
      } else {
         for (const auto& mr : *Mrec ) {
            hist[2]->Fill(mr);
         }
      }
   }
   return hist;
}

//--------------------------------------------------------------------
vector<TH1D*> fill_hist(int date, bool sys_NoHC=false)
//--------------------------------------------------------------------
{
#include "norm.h"

   string sd(Form("%02i",date%100));
   vector<string> hname = {
      "data"+sd, "data3650_"+sd,
      "mc"+sd+"sig", "mc"+sd+"bg1", "mc"+sd+"bg2"
   };

   vector<TH1D*> hst(hname.size(),nullptr);

   // read from cache file
   string cachef = Dir +
      ((sys_NoHC) ? "NoHC/" : "") +
      string(Form("MrecFitSB_%i.root",date));
   if ( std::filesystem::exists( cachef )  ) {
      TFile* froot = TFile::Open(cachef.c_str(),"READ");
      if ( !froot ) {
         cerr << "ERROR in "<< __func__
            << ": unreadable file: " << cachef << endl;
         exit(EXIT_FAILURE);
      }
      for ( size_t i = 0; i < hname.size(); ++i ) {
         hst[i]=(TH1D*)froot->Get(hname[i].c_str());
      }
      return hst;
   }

   string datafile( Form("data_%02ipsip_all.root",date%100) );
   string mcincfile( Form("mcinc_%02ipsip_all.root",date%100) );

   auto hs = fill_mrec(datafile, hname[0], date);
   hst[0] = hs[0];

   hs = fill_mrec("data_3650_2021.root", hname[1], 2021);
   hst[1] = hs[0];
   hst[1]->Scale(C_Dat.at(date));

   hs = fill_mrec(mcincfile, "mc"+sd, date, true, sys_NoHC);
   hst[2]=hs[0];
   hst[2]->Scale(MC_Dat.at(date));

   hst[3]=hs[1];
   hst[3]->Scale(MC_Dat.at(date));

   hst[4]=hs[2];
   hst[4]->Scale(MC_Dat.at(date));

   // save histos in cache file
   TFile* froot = TFile::Open(cachef.c_str(),"NEW");
   for ( const auto& h : hst ) {
      h->Write();
   }
   froot->Close();
   return hst;
}

/*
//--------------------------------------------------------------------
vector<TH1D*> fill_hist_IOcheck(int date) {
//--------------------------------------------------------------------
   string sd(Form("%02i",date%100));
   vector<string> hname = {
      "MCdata"+sd, "empty"+sd,
      "mc"+sd+"sig", "mc"+sd+"bg1", "mc"+sd+"bg2"
   };

   vector<TH1D*> hst(hname.size(),nullptr);

   // check cache file
   string cachef = Dir + string(Form("MrecFitSB_%i_IO.root",date));
   if ( std::filesystem::exists( cachef )  ) {
      TFile* froot = TFile::Open(cachef.c_str(),"READ");
      if ( !froot ) {
         cerr << "ERROR in "<< __func__
            << ": unreadeble file: " << cachef << endl;
         exit(EXIT_FAILURE);
      }
      for ( size_t i = 0; i < hname.size(); ++i ) {
         hst[i]=(TH1D*)froot->Get(hname[i].c_str());
      }
      return hst;
   }

   string mcincfile( Form("mcinc_%02ipsip_all.root",date%100) );

   hst[0]=fill_mrec(mcincfile, hname[0], 10);

   hst[1]=(TH1D*)hst[0]->Clone(hname[1].c_str());
   hst[1]->Reset(); // must be empty

   hst[2]=fill_mrec(mcincfile,hname[2],11);

   hst[3]=fill_mrec(mcincfile,hname[3],2);
   hst[4]=fill_mrec(mcincfile,hname[4],3);

   // save histos in cache file
   TFile* froot = TFile::Open(cachef.c_str(),"NEW");
   for ( const auto& h : hst ) {
      h->Write();
   }
   froot->Close();
   return hst;
}
*/

// {{{1 Function for fitting
//--------------------------------------------------------------------
double ChebN(int nch, double Xmin, double Xmax,
      double x, const double* p)
//--------------------------------------------------------------------
{
   // calculate value of the Chebyshev polynomial of 'nch' order
   // [Xmin Xmax] is the range of polynomial orthogonality
   if (nch == 0) {
      return p[0];
   }
   x = (2*x-Xmin-Xmax)/(Xmax-Xmin); // map [Xmin,Xmax] -> [-1,+1]
   double sum = p[0] + x*p[1];
   if (nch == 1) {
      return sum;
   }
   double T0 = 1, T1 = x;
   for ( int i = 2; i <= nch; ++i ) {
      double Tmp = 2*x*T1 - T0;
      sum += p[i]*Tmp;
      T0 = T1;
      T1 = Tmp;
   }
   return sum;
}

//--------------------------------------------------------------------
double PolN(double x, const double* p, int n)
//--------------------------------------------------------------------
{
   // calculate value of the polynomial by Horner's method:
   // here 'n' is the polynomial order, size of coefficients is p[n+1]
   double res = 0;
   for ( int i = n; i >= 0; --i ) {
      res = res*x + p[i];
   }
   return res;
}

// {{{1 Class for fitting
//--------------------------------------------------------------------
class myChi2 {
   public:
      myChi2(const vector<TH1D*>& hst, double ExMin, double ExMax);

      // operator() => chi^2 function
      double operator() (const double* p) { return Chi2(p); }
      double Chi2(const double* p);

      unsigned int Size() const {
         return data.size();
      }

      void SetNpol(int n) {
         npol = n;
      }

   private:
      vector<double> en; // energy
      vector<double> data, er_data;
      vector<double> cont, er_cont;
      vector<double> sig, er_sig;
      vector<double> bg, er_bg;
      vector<double> bgn, er_bgn;

      int npol = 3;
      const double Emin_chb = 3.0, Emax_chb = 3.2;
};

// {{{2 ctor:
//--------------------------------------------------------------------
myChi2::myChi2(const vector<TH1D*>& hst, double ExMin, double ExMax)
//--------------------------------------------------------------------
{
   // ExMin/Max -> boundary of excluded region
   int nbins = hst[0]->GetNbinsX();
   en.reserve(nbins);
   data.reserve(nbins);
   er_data.reserve(nbins);
   cont.reserve(nbins);
   er_cont.reserve(nbins);
   sig.reserve(nbins);
   er_sig.reserve(nbins);
   bg.reserve(nbins);
   er_bg.reserve(nbins);
   bgn.reserve(nbins);
   er_bgn.reserve(nbins);
   for ( int n = 1; n <= nbins; ++n ) {
      double x = hst[0]->GetBinCenter(n);
      if ( x > ExMin && x < ExMax ) {
         continue;
      }
      // en.push_back(x-3.1); // subtract midpoint
      en.push_back(x); // for Chebyshev polynomial

      data.push_back(hst[0]->GetBinContent(n));
      er_data.push_back(SQ(hst[0]->GetBinError(n)));

      cont.push_back(hst[1]->GetBinContent(n));
      er_cont.push_back(SQ(hst[1]->GetBinError(n)));

      sig.push_back(hst[2]->GetBinContent(n));
      er_sig.push_back(SQ(hst[2]->GetBinError(n)));

      bg.push_back(hst[3]->GetBinContent(n));
      er_bg.push_back(SQ(hst[3]->GetBinError(n)));

      bgn.push_back(hst[4]->GetBinContent(n));
      er_bgn.push_back(SQ(hst[4]->GetBinError(n)));
   }
}

// {{{2 Chi2() function
//--------------------------------------------------------------------
double myChi2::Chi2(const double* p)
//--------------------------------------------------------------------
{
   static const double maxResValue = DBL_MAX / 1000;

   double chi2 = 0;
   size_t n = Size();
   double A = p[0];

   for ( size_t i = 0; i < n; i++ ) {
      double x = en[i];
      // double bgscl = PolN(x,p,npol);
      double bgscl = ChebN(npol,Emin_chb,Emax_chb,x,p);
      if ( bgscl < 0 ) {
         chi2 += maxResValue;
         continue;
      }
      double bgscl2 = bgscl*bgscl;

      double dif = data[i] -
         ( cont[i] +
           A*sig[i] +
           bgscl * (bg[i] + bgn[i])
         );
      double er2 = er_data[i] +
         er_cont[i] +
         SQ(A)*er_sig[i] +
         SQ(bgscl) * (er_bg[i] + er_bgn[i]);

      double resval = SQ(dif) / er2;
      // avoid inifinity or nan in chi2
      if( resval < maxResValue ) {
         chi2 += resval;
      } else {
         chi2 += maxResValue;
      }
   }
   return chi2;
}

//--------------------------------------------------------------------

// {{{1 Corrections for histograms and print final numbers
//--------------------------------------------------------------------
TH1D* rescale_bg(const TH1D* hst, const vector<double>& par)
//--------------------------------------------------------------------
{
   // rescale ignoring over and underflow bins

   // polynomial parameters
   const double Emin_chb = 3.0, Emax_chb = 3.2;
   int npol = par.size() - 1;
   const double* p = par.data();

   // create copy of hst
   string new_name = string(hst->GetName()) + string("_Rescaled");
   TH1D* hst_n = (TH1D*)hst->Clone(new_name.c_str());

   int nbins = hst->GetNbinsX();
   for(int n = 1; n <= nbins; ++n) {
      double x = hst->GetBinCenter(n);
      // double w = PolN(x-3.1,p,npol);
      double w = ChebN(npol,Emin_chb,Emax_chb,x,p);
      double d = hst->GetBinContent(n);
      double e = hst->GetBinError(n);
      hst_n->SetBinContent(n,d*w);
      hst_n->SetBinError(n,e*w);
   }
   return hst_n;
}

//--------------------------------------------------------------------
void print_Numbers(const vector<TH1D*>& hst, const TH1D* SumBG,
      double Emin, double Emax, double N0=0., bool verbose=false)
//--------------------------------------------------------------------
{
   // ATTENTION: here I assume:
   // hst[0] - data
   // hst[2] - correct signal MC (for I/O check only)
   // sumBG - CONT(norm to Lum_data) + BG(after rescale_bg)

   double ndata   = 0;
   double er_data = 0;
   double nbg     = 0;
   double er_bg   = 0;
   double ntrue   = 0; // I/O check

   double Emin_hst = 0, Emax_hst = 0;
   int nbins = hst[0] -> GetNbinsX();
   double wbin = hst[0] -> GetBinWidth(1);
   for(int n = 1; n <= nbins; ++n) {
      double x = hst[0]->GetBinCenter(n);
      if( x < Emin ) {
         continue;
      }
      if( Emin_hst < 1. ) {
         Emin_hst = x - wbin/2;
      }
      if( x > Emax ) {
         Emax_hst = hst[0] -> GetBinCenter(n-1) + wbin/2;
         break;
      }
      ndata   += hst[0]->GetBinContent(n);
      er_data += SQ(hst[0]->GetBinError(n));
      nbg     += SumBG->GetBinContent(n);
      er_bg   += SQ(SumBG->GetBinError(n));
      ntrue   += hst[2]->GetBinContent(n); // I/O check
   }

   double n1   = ndata-nbg;
   double er1  = sqrt(er_data+er_bg);
   er_data     = sqrt(er_data);
   er_bg       = sqrt(er_bg);

   if ( N0 < 1 ) {
      printf("[%.3f,%.3f]: ", Emin_hst,Emax_hst);
      printf("NJpsi(SB)= $(%.3f \\pm %.3f)\\times10^6$",
            n1/1e6,er1/1e6);
      printf(" = %.3fe6\n",n1/1e6);
   } else {
      double del = 100*fabs(N0 - n1) / N0;
      printf("[%.3f, %.3f]: ", Emin_hst,Emax_hst);
      printf("NJpsi(data-bg)= $%.3f\\times10^6$ & %.2f %%\n",
            n1/1e6,del);
   }
   if ( !verbose ) {
      return;
   }

   printf("\n");
   printf(" Number of events in [%.3f, %.3f]:\n", Emin_hst,Emax_hst);
   printf(" Data:            %9.0f +/- %4.0f\n", ndata,er_data);
   printf(" MC(bg):          %9.0f +/- %4.0f\n", nbg,er_bg);
   printf(" NJpsi(data-bg):  %9.0f +/- %4.0f\n", n1,er1);
   printf("\n");
   printf(" The correct number of J/Psi is %9.0f\n",ntrue);
}

// {{{1 Fit
//--------------------------------------------------------------------
void DoFitSB(int date, int Cx=800, int Cy=800, int Sys=0)
//--------------------------------------------------------------------
{
   // IO check: trivially get the exact result!
   // bool IOcheck = true;
   bool IOcheck = false;

   bool sys_NoHC = Sys==5;
   vector<TH1D*> hst = fill_hist(date,sys_NoHC);

   // limits for drawing
   double maxWin = 0, minWin = 0;
   if ( date==2009 ) {
      maxWin = 7.e6;
      minWin = 0.99e5;
   } else if(date==2012) {
      maxWin = 2.e7;
      minWin = 3e5;
   } else if(date==2021) {
      maxWin = 1.3e8;
      minWin = 2e6;
   }
   hst[0]->SetMaximum(maxWin);
   hst[0]->SetMinimum(minWin);

   // exclude (ExMin; ExMax) from fit
   double ExMin = 3.06;
   double ExMax = 3.14;
   if ( Sys == 3 || Sys == 4 ) {
      // systematic study (shift +/- 0.01) : sysSB-Range
      double ddE = 0.01 * (1-2*(Sys%2)); //  -/+ 0.01
      ExMin += ddE;
      ExMax -= ddE;
   }

   myChi2 chi2_fit(hst,ExMin,ExMax);

   // ========================= Fit with ROOT ========================
   // == fit configuration
   ROOT::Fit::Fitter fitter;
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   // Possible printing levels are from 0 (minimal) to 3 (maximum)
   fitter.Config().MinimizerOptions().SetPrintLevel(int(Sys==0));
   // print the default minimizer option values
   ROOT::Math::MinimizerOptions min_opt;
   if ( !Sys ) {
      min_opt.Print();
   }

   // Npol is the order of polynomial, systematic study: +/-1
   size_t Npol = (date==2009) ? 4 : 5;
   if ( Sys == 1 || Sys == 2  ) {
      Npol += 1 - 2 * (Sys%2); // -/+ 1
   }
   chi2_fit.SetNpol(Npol);

   // set start values for parameters,
   const size_t Npar = Npol+1; // number of parameters
   vector<string> par_name(Npar);
   for( size_t i = 0; i <= Npol; ++i ) {
      par_name[i] = string( Form("p%zu",i) );
   }
   vector<double> par_ini;
   if ( date == 2009 ) {
      // par_ini = { 1.088, 0.009, 0.005, 0.002, 0.002 }; // n3
      par_ini = { 1.086, 0.009, 0.004, 0.002, 0.002 };
   } else if ( date == 2012 ) {
      // par_ini = { 1.106, 0.016, 0.006, 0.002, 0.002, 0.002 }; // n3
      par_ini = { 1.104, 0.016, 0.005, 0.002, 0.002, 0.002 };
   } else if ( date == 2021 ) {
      // par_ini = { 1.071, 0.001, 0.006, 0.002, 0.002, 0.001 }; // n3
      par_ini = { 1.069, 0.0003, 0.006, 0.002, 0.002, 0.001 };
   }
   if ( par_ini.size() != Npar ) {
      par_ini.resize(Npar,0.);
      cout << "WARNING: set size of par_ini to " << Npar << endl;
   }
   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for( size_t i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix  parameters
   // fitter.Config().ParSettings(0).SetValue(1.);
   // fitter.Config().ParSettings(0).Fix();
   // fitter.Config().ParSettings(2).SetValue(0.);
   // fitter.Config().ParSettings(2).Fix();

   // == Fit
   fitter.FitFCN(Npar, chi2_fit, nullptr, chi2_fit.Size(), true);

   // to obtain reliable errors:
   fitter.CalculateHessErrors();
   // fitter.CalculateMinosErrors(); // check

   // == Fit result
   ROOT::Fit::FitResult res = fitter.Result();
   if ( !Sys ) {
      res.Print(cout);
   }
   // double chi2   = res.Chi2(); // = -1 in root-6.30
   double chi2   = res.MinFcnValue(); // chi2 or likelihood
   int ndf       = res.Ndf();
   const vector<double> par   = res.Parameters();
   const vector<double> er_par = res.Errors();

   string sepline(70,'='); // separation line
   printf("%s\n",sepline.c_str());
   // printf("Error matrix and correlations\n");
   // res.PrintCovMatrix(cout); // print error matrix and correlations
   // printf("%s\n",sepline.c_str());

   // ========================= draw results ========================
   auto name = (!Sys) ? Form("c_%i",date) : Form("c_%i_%i",date,Sys);
   TCanvas* c1 = new TCanvas(name,name,0,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   set_draw_opt(hst);
   SetHstFace(hst[0]);
   hst[0]->GetXaxis()->SetTitleOffset(1.1);
   hst[0]->GetYaxis()->SetTitleOffset(1.1);

   hst[0]->Draw("E"); // data

   // box to show not fitted region
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-10);
   box->SetLineWidth(2);
   box->DrawBox(ExMin,minWin,ExMax,maxWin);

   TBox* box2 = new TBox;
   box2->SetFillStyle(3001);
   box2->SetFillColor(kGreen-6);
   box2->SetLineColor(kGreen-6);
   // box2->DrawBox(3.092,minWin,3.102,maxWin);

   hst[0]->Draw("E, SAME"); // redraw data

   TH1D* SumBG = (TH1D*)hst[1]->Clone("SumBG"); // clone Continuum
   // rescale BG:
   TH1D* hstBG1 = rescale_bg( hst[3], par );
   TH1D* hstBG2 = rescale_bg( hst[4], par );
   SumBG->Add(hstBG1);
   SumBG->Add(hstBG2);
   SumBG->Draw("SAME,HIST");

   TH1D* MCsig =(TH1D*)hst[2]->Clone("MCsig"); // clone MC Signal
   // rescale Signal
   MCsig->Scale(par[0]);
   MCsig->SetLineStyle(kDashed);
   MCsig->SetLineWidth(2);
   MCsig->Draw("SAME,HIST");

   TH1D* SumMC =(TH1D*)MCsig->Clone("SumMC");
   SumMC->Add(SumBG);
   SumMC->SetLineColor(kRed+1);
   SumMC->SetLineStyle(kSolid);
   SumMC->Draw("SAME,HIST");

   // hstBG2->Draw("SAME,HIST"); // Rescaled Bg2

   TLegend* leg = new TLegend(0.58,0.69,0.892,0.89);
   if ( Sys ) {
      leg->SetHeader( Form("Systematics #%i",Sys), "C" );
   }
   if ( IOcheck ) {
      leg->SetHeader( Form("Pseudo-Data %i",date), "C" );
   }

   leg->AddEntry(hst[0],Form("Data %i",date), "EP");

   leg->AddEntry(SumBG,
         Form("#color[%i]{Sum of all backgrounds}",
         SumBG->GetLineColor() ),"L");

   const char* Lpp = "#pi#lower[-0.8]{#scale[0.8]{#plus}}"
      "#pi#lower[-0.8]{#scale[0.8]{#minus}}";
   leg->AddEntry(MCsig,
         Form("#color[%i]{MC signal}",MCsig->GetLineColor()), "L" );

   leg->AddEntry(SumMC,
         Form("#color[%i]{Full fit}",SumMC->GetLineColor()), "L");

   // leg->AddEntry(hstBG2,
         // Form("#color[%i]{MC bg non %sJ/#Psi}",
            // hstBG2->GetLineColor(),Lpp ),"L");

   leg->AddEntry(box, "area excluded from fitting","F");
   leg->Draw();

   double Ypt = 0.68 - 0.04*(Npol-3);
   TPaveText* pt = new TPaveText(0.11,Ypt,0.40,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#bf{Fit in [%.2f,%.2f] & [%.2f,%.2f] }",
            3.0,ExMin,ExMax,3.2));
   pt->AddText( Form("#chi^{2} / ndf   %.0f / %i",chi2,ndf) );
   for ( size_t i = 0; i < par.size(); i++ ) {
      pt->AddText( Form("%s =       %.4f #pm %.4f",
               par_name[i].c_str(),par[i],er_par[i]) );
   }
   pt->Draw();
   gPad->RedrawAxis();
   c1->Update();

   // print for E in [Emin,Emax]
   double Ns = 0., Nw = 0.;
   if ( !Sys ) {
      printf("* T_%zu: Chi2 / NDf = %.0f / %i = %.1f\n",
            Npol,chi2,ndf,chi2/ndf );
   } else {
      printf("* T_%zu, [%.2f,%.2f]&[%.2f,%.2f]:"
            " Chi2 / NDf = %.0f / %i = %.1f\n",
            Npol,3.0,ExMin,ExMax,3.2,chi2,ndf,chi2/ndf );
      if ( date == 2009 ) {
         // Ns = ; Nw = ;
         Ns = 15.998e6; Nw = 18.066e6;
      } else if ( date == 2012 ) {
         Ns = 49.903e6; Nw = 56.728e6;
      } else if ( date == 2021 ) {
         Ns = 324.947e6; Nw = 370.197e6;
      }
   }
   printf("%s\n",sepline.c_str());
   print_Numbers(hst,SumBG,3.092,3.102,Ns);
   print_Numbers(hst,SumBG,3.06,3.14,Nw);
   printf("%s\n",sepline.c_str());

   string suff;
   if ( Sys ) {
      suff = string( Form("_sys%i",Sys) );
   }
   string pdf( Form("Mrec%d_fsb%s_T%zu",date,suff.c_str(),Npol) );
   if ( IOcheck ) {
      pdf += "_IO";
   }
   pdf += ".pdf";
   c1->Print(pdf.c_str());
}

// {{{1 Draw Mrec
//--------------------------------------------------------------------
void MrecDraw(int date, bool zoom=false, int Cx=800, int Cy=800)
//--------------------------------------------------------------------
{
   vector<TH1D*> hst = fill_hist(date);
   if ( date==2009 ) {
      hst[0]->SetMaximum(7.e6);
      hst[0]->SetMinimum(3.e4);
      if( zoom ) {
         hst[0]->SetMaximum(0.65e6);
         hst[0]->SetMinimum(0.);
         hst[0]->GetYaxis()->SetMaxDigits(3);
      }
   } else if(date==2012) {
      hst[0]->SetMaximum(2.e7);
      hst[0]->SetMinimum(1.e5);
      if( zoom ) {
         hst[0]->SetMaximum(1.9e6);
         hst[0]->SetMinimum(0.);
         hst[0]->GetYaxis()->SetMaxDigits(3);
      }
   } else if(date==2021) {
      hst[0]->SetMaximum(1.3e8);
      hst[0]->SetMinimum(5.e5);
      if( zoom ) {
         hst[0]->SetMaximum(1.3e7);
         hst[0]->SetMinimum(0.);
         hst[0]->GetYaxis()->SetMaxDigits(3);
      }
   }

   auto name = Form("c1_%i_%i",int(zoom),date);
   TCanvas* c1 = new TCanvas(name,name, 0, int(zoom)*Cy/2, Cx, Cy);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(!zoom);

   set_draw_opt(hst);
   SetHstFace(hst[0]);
   hst[0]->GetXaxis()->SetTitleOffset(1.1);
   hst[0]->GetYaxis()->SetTitleOffset(1.1);

   hst[0]->Draw("E"); // data

   hst[2]->SetLineStyle(kDashed);
   hst[2]->Draw("SAME,HIST"); // Signal

   TH1D* hstSUM=(TH1D*)hst[2]->Clone("SumAll");
   hstSUM->Add(hst[1]);
   hstSUM->Add(hst[3]);
   hstSUM->Add(hst[4]);
   hstSUM->SetLineColor(kRed+1);
   hstSUM->Draw("SAME,HIST"); // Sum

   hst[3]->Draw("SAME,HIST"); // Bg1
   hst[4]->Draw("SAME,HIST"); // Bg2

   TH1D* hstC10 = (TH1D*)hst[1]->Clone("Cont10");
   hstC10->Scale(10.);
   hstC10->Draw("SAME,HIST"); // Cont

   TLegend* leg = new TLegend(0.55,0.60,0.892,0.89);
   leg->AddEntry(hst[0], Form("Data %i",date), "EP");

   // Form("#color[%i]{MCsig #plus MCbg #plus Cont.}",
   leg->AddEntry(hstSUM,
         Form("#color[%i]{MC #plus Non-resonant bg}",
            hstSUM->GetLineColor() ),"L");

   const char* Lpp = "#pi#lower[-0.8]{#scale[0.8]{#plus}}"
      "#pi#lower[-0.8]{#scale[0.8]{#minus}}";
   // Form("#color[%i]{MC %sJ/#Psi correct %s}",
   leg->AddEntry(hst[2],
         Form("#color[%i]{MC signal}",hst[2]->GetLineColor()),"L");

   // Form("#color[%i]{MC %sJ/#Psi wrong %s}",
   leg->AddEntry(hst[3],
         Form("#color[%i]{MC bg, incorrect %s}",
            hst[3]->GetLineColor(),Lpp ),"L");

   // Form("#color[%i]{MC bg non %sJ/#Psi}",
   leg->AddEntry(hst[4],
         Form("#color[%i]{MC bg, non %sJ/#Psi}",
            hst[4]->GetLineColor(),Lpp ),"L");

   // Form("#color[%i]{Continuum bg #times 10}",
   leg->AddEntry(hstC10,
         Form("#color[%i]{Non-resonant bg #times 10}",
            hstC10->GetLineColor() ),"L");
   leg->Draw();

   gPad->RedrawAxis();

   c1->Update();
   string pdf(Form("Mrec%i%s.pdf",date,(zoom ? "_zoom" : "")));
   c1->Print(pdf.c_str());
}

// {{{1 Main
//--------------------------------------------------------------------
void MrecFitSB()
//--------------------------------------------------------------------
{
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(112);
   gStyle->SetStatFont(62);
   gStyle->SetLegendFont(42);

   //========================================================
   // set the name of the folder with the root files
   // Dir = "prod_v709n3/";
   Dir = "prod_v709n4/";
   //========================================================

   size_t Cx = 800, Cy = 640; // canvas sizes, X/Y = 1.25

   // Fig.4
   for ( int date : { 2009,2012,2021 } ) {
      bool zoom = true;
      // MrecDraw(date, !zoom, Cx, Cy);
      // MrecDraw(date,  zoom, Cx, Cy);
   }


   // main fit, Fig.6
   // for ( int date : { 2009,2012,2021 } ) {
      // cout << " >>> " << date << " <<<" << endl;
      // DoFitSB(date, Cx, Cy);
   // }

   // SB systematic: NT && Range
   // for ( int date : { 2009,2012,2021 } ) {
      // cout << " >>> " << date << " <<<" << endl;
      // for ( int sys : {1,2,3,4} ) {
         // DoFitSB(date, Cx, Cy, sys);
      // }
   // }

   // SB systematic: No HC (slow)
   // for ( int date : { 2009,2012,2021 } ) {
      // cout << " >>> " << date << " <<<" << endl;
      // DoFitSB(date, Cx, Cy, 5);
   // }

}
