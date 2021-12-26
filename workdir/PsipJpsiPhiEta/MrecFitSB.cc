// Calculate number of Psi(2S) -> J/Psi pi+ pi- decays in data
// using side-band method.
// We fit the distribution of recoil mass (Mrec) of two pions
// in the area far from the peak of J/Psi by rescaling inclusive MC.
// Here We rescale only background from non pi+ pi- J/Psi
// by a polynomial.
// Number of signal events is calculated by subtracting the estimated
// background from the data.  (The contribution of the continuum
// is taken into account.)
//                MrecDraw() -> Mrec_YEAR.pdf
//                DoFitSB()  -> Mrec_YEAR_fsb.pdf

// for prod-12 it should be commented
// #include "ReWeightTrkPid_11.h"

// {{{1 helper functions and constants
//--------------------------------------------------------------------
// GLOBAL: name of folder with root files
static const string dir("prod-12/");  // must be the same as prod-11

//--------------------------------------------------------------------
constexpr double SQ(double x) {
//--------------------------------------------------------------------
   return x*x;
}

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
void set_draw_opt(TH1D* hst[]) {
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
   hst[3]->SetLineColor(kBlue+3);
   hst[3]->SetLineWidth(2);
// MC bg not pi+pi-J/Psi
   hst[4]->SetLineColor(kMagenta+1);
   hst[4]->SetLineWidth(2);
}

// {{{1 Fill histograms
//--------------------------------------------------------------------
TH1D* fill_mrec(string fname, string hname, int type=0) {
//--------------------------------------------------------------------
   // type = 0 all recoil masses (default)
   // type = 1 recoil mass of true pi+pi- pair (MC)
   // type = 2 background for dec==64 (pi+pi-J/Psi)
   // type = 3 background for dec!=64 (other Psi' decays)

   fname = dir + fname;
   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if ( froot == 0 ) {
      cerr << "ERROR in "<< __func__
           << ": can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");

   //Declaration of leaves types
   Float_t         Mrs;
   Float_t         MrsW;
   Float_t         Mrb;
   Float_t         Ptp;
   Float_t         Ptm;
   UShort_t        dec;
   // Set branch addresses.
   nt1->SetBranchAddress("Mrs",&Mrs);
   nt1->SetBranchAddress("MrsW",&MrsW);
   nt1->SetBranchAddress("Mrb",&Mrb);
   nt1->SetBranchAddress("Ptp",&Ptp);
   nt1->SetBranchAddress("Ptm",&Ptm);
   nt1->SetBranchAddress("dec",&dec);

   TH1D* hst = new TH1D(hname.c_str(),
         ";M^{rec}_{#pi^{#plus}#pi^{#minus }}, GeV/c^{2}"
         ";Entries/0.001 GeV/c^{2}",
         200,3.0,3.2 // 1bin = 1MeV
//          400,3.0,3.2 // 1bin = 0.5MeV
                       );
   hst->Sumw2(true);

   string dr = string("Mrec>>") + hname;
   TCut cut;
   if ( type == 1 ) {
      dr = string("Mrs>>") + hname;
      cut = TCut("MrsW*(dec==64)"); // MrsW pi+pi- corrections

      /* prod-11
      bool is2009 = (fname.find("_09") != string::npos);
      int date = (is2009) ? 2009 : 2012;

      Long64_t nentries = nt1 -> GetEntries();
      for (Long64_t i=0; i<nentries;i++) {
         nt1 -> GetEntry(i);

         if ( dec != 64 ) continue;
         if ( Mrs < 3.0 || Mrs > 3.2 ) continue;

         // correction for pi+,pi- eff:
         double wp = ReWeightTrkPid(date,0,Ptp);
         double wm = ReWeightTrkPid(date,0,Ptm);
         double W = wp*wm;
         hst -> Fill(Mrs,W);
      }
      return hst;
      */
   } else if ( type == 2 ) {
      cut = TCut("dec==64");
   } else if ( type == 3 ) {
      cut = TCut("dec!=64");
   }

   nt1->Draw(dr.c_str(),cut,"goff");

   return hst;
}

//--------------------------------------------------------------------
void fill_hist2009(TH1D* hst[]) {
//--------------------------------------------------------------------
#include "norm.h"

   const string hname[5] = {
      "data09", "data3650_09",
      "mc09sig", "mc09bg1", "mc09bg2"
   };

   // check cache file
   const string cachef = dir + string("MrecFitSB_2009.root");
   TFile* froot = TFile::Open(cachef.c_str(),"READ");
   if ( froot != 0 ) {
      for ( int i = 0; i < 5; ++i ) {
         hst[i]=(TH1D*)froot->Get(hname[i].c_str());
      }
      return;
   }

   hst[0]=fill_mrec("data_09psip_all.root", hname[0]);

   hst[1]=fill_mrec("data_3650_all.root", hname[1]);
   hst[1]->Scale(Cto09);

   hst[2]=fill_mrec("mcinc_09psip_all.root",hname[2],1);
   hst[2]->Scale(Ito09);

   hst[3]=fill_mrec("mcinc_09psip_all.root",hname[3],2);
   hst[3]->Scale(Ito09);

   hst[4]=fill_mrec("mcinc_09psip_all.root",hname[4],3);
   hst[4]->Scale(Ito09);

   // save histos in cache file
   froot = TFile::Open(cachef.c_str(),"NEW");
   for ( int i = 0; i < 5; ++i ) {
      hst[i]->Write();
   }
   froot->Close();
}

//--------------------------------------------------------------------
void fill_hist2012(TH1D* hst[]) {
//--------------------------------------------------------------------
#include "norm.h"

   const string hname[5] = {
      "data12", "data3650_12",
      "mc12sig", "mc12bg1", "mc12bg2"
   };

   // check cache file
   const string cachef = dir + string("MrecFitSB_2012.root");
   TFile* froot = TFile::Open(cachef.c_str(),"READ");
   if ( froot != 0 ) {
      for ( int i = 0; i < 5; ++i ) {
         hst[i]=(TH1D*)froot->Get(hname[i].c_str());
      }
      return;
   }

   hst[0]=fill_mrec("data_12psip_all.root", hname[0]);

   hst[1]=fill_mrec("data_3650_all.root", hname[1]);
   hst[1]->Scale(Cto12);

   hst[2]=fill_mrec("mcinc_12psip_all.root",hname[2],1);
   hst[2]->Scale(Ito12);

   hst[3]=fill_mrec("mcinc_12psip_all.root",hname[3],2);
   hst[3]->Scale(Ito12);

   hst[4]=fill_mrec("mcinc_12psip_all.root",hname[4],3);
   hst[4]->Scale(Ito12);

   // save histos in cache file
   froot = TFile::Open(cachef.c_str(),"NEW");
   for ( int i = 0; i < 5; ++i ) {
      hst[i]->Write();
   }
   froot->Close();
}

// {{{1 Function for fitting
//-------------------------------------------------------------------------
double ChebN(int nch, double Xmin, double Xmax, double x, const double* p) {
//-------------------------------------------------------------------------
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

//-------------------------------------------------------------------------
double PolN(double x, const double* p, int n) {
//-------------------------------------------------------------------------
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
      // ctor:
      //--------------------------------------------------------------
      myChi2(TH1D* hst[], double ExMin, double ExMax) {
      //--------------------------------------------------------------
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
//             en.push_back(x-3.1); // subtract midpoint
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

      // Size function
      //--------------------------------------------------------------
      unsigned int Size() const {
      //--------------------------------------------------------------
         return data.size();
      }

      //--------------------------------------------------------------
      void SetNpol(int n) {
      //--------------------------------------------------------------
         npol = n;
      }

      // chi^2 - function
      // the signature of this operator() MUST be exactly this:
      //--------------------------------------------------------------
      double operator() (const double* p) {
      //--------------------------------------------------------------
         static const double maxResValue = DBL_MAX / 1000;

         double chi2 = 0;
         unsigned int n = Size();
         for(unsigned int i = 0; i < n; i++) {
            double x = en[i];
//             double bgscl = PolN(x,p,npol);
            double bgscl = ChebN(npol,Emin_chb,Emax_chb,x,p);
            if ( bgscl < 0 ) {
               return maxResValue;
            }
            double bgscl2 = bgscl*bgscl;

            double dif = data[i] -
                         ( sig[i] +
                           cont[i] +
                           bgscl * (bg[i] + bgn[i]) );
            double er2 = er_data[i] +
                         er_sig[i] +
                         er_cont[i] +
                         bgscl2 * (er_bg[i] + er_bgn[i]);

            double resval = dif*dif / er2;
            // avoid inifinity or nan in chi2
            if( resval < maxResValue ) {
               chi2 += resval;
            } else {
               return maxResValue;
            }
         } // end of for()

         return chi2;
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
//--------------------------------------------------------------------

// {{{1 Corrections for histograms and print final numbers
//--------------------------------------------------------------------
TH1D* rescale_bg(TH1D* hst, const double* p, int npol) {
//--------------------------------------------------------------------
   // rescale ignoring over and underflow bins
   const double Emin_chb = 3.0, Emax_chb = 3.2;

   // create copy of hst
   string new_name = string(hst->GetName()) + string("_Rescaled");
   TH1D* hst_n = (TH1D*)hst->Clone(new_name.c_str());

   int nbins = hst->GetNbinsX();
   for(int n = 1; n <= nbins; ++n) {
      double x = hst->GetBinCenter(n);
//       double w = PolN(x-3.1,p,npol);
      double w = ChebN(npol,Emin_chb,Emax_chb,x,p);
      double d = hst->GetBinContent(n);
      double e = hst->GetBinError(n);
      hst_n->SetBinContent(n,d*w);
      hst_n->SetBinError(n,e*w);
   }
   return hst_n;
}

//--------------------------------------------------------------------
void print_Numbers(TH1D* hst[], double Emin, double Emax) {
//--------------------------------------------------------------------
   // ATTENTION: here I assume:
   // [0] - data
   // [7] - CONT(norm to Lum_data) + BG(after rescale_bg)

   double ndata   = 0;
   double er_data = 0;
   double nbg     = 0;
   double er_bg   = 0;

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
      nbg   += hst[7]->GetBinContent(n);
      er_bg += SQ(hst[7]->GetBinError(n));
   }

   double n1   = ndata-nbg;
   double er1  = sqrt(er_data+er_bg);
   er_data     = sqrt(er_data);
   er_bg       = sqrt(er_bg);

   printf("\n");
   printf(" Number of events in [%.3f, %.3f]:\n", Emin_hst,Emax_hst);
   printf(" Data:            %9.0f +/- %4.0f\n", ndata,er_data);
   printf(" MC(bg):          %9.0f +/- %4.0f\n", nbg,er_bg);
   printf(" NPsip(data-bg):  %9.0f +/- %4.0f\n", n1,er1);
   printf("\n");
}

// {{{1 Fit
//--------------------------------------------------------------------
void DoFitSB(int date) {
//--------------------------------------------------------------------
   TH1D* hst[20];
   if ( date==2009 ) {
      fill_hist2009(hst);
      hst[0]->SetMaximum(7.e6);
      hst[0]->SetMinimum(0.99e5);
   } else if(date==2012) {
      fill_hist2012(hst);
      hst[0]->SetMaximum(2.e7);
      hst[0]->SetMinimum(3.e5);
   }

   // exclude (ExMin; ExMax) from fit
   double ExMin = 3.06;
   double ExMax = 3.14;
   // systematic study (shift +/- 0.01)
//    double ddE = -0.01;
//    ExMin -= ddE;
//    ExMax += ddE;

   myChi2 chi2_fit(hst,ExMin,ExMax);

   // ========================= Fit with ROOT ========================
   // == fit configuration
   ROOT::Fit::Fitter fitter;
   // set parameters of fitter: (Minuit,Minuit2,Fumili,GSLMultiFit...)
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");
   fitter.Config().MinimizerOptions().SetPrintLevel(3);
   // print the default minimizer option values
   ROOT::Math::MinimizerOptions min_opt;
   min_opt.Print();

   const unsigned int Npol = 3; // order of polynomial (base=3)
                                // systematic study: 2 && 4
   chi2_fit.SetNpol(Npol);

   // set names and start values for parameters
   vector<string> par_name;
   for ( int i = 0; i <= Npol; ++i ) {
      par_name.push_back( string("p") + to_string(i) );
   }
   vector<double> par_ini(Npol+1,1.); // for Chebyshev
   if ( Npol == 3 ) {
      if ( date == 2009 ) {
         par_ini = { 1.0343, 0.0069, 0.0083, 0.0011 };
      } else if ( date == 2012 ) {
         par_ini = { 1.0478, 0.0044, 0.0076, 0.0008 };
      }
   }
   if ( par_ini.size() != Npol+1 ) {
      cout << " WARNING: something wrong for Npol= " << Npol << endl;
      getchar();
      par_ini.resize(Npol+1);
   }

   const unsigned int Npar = par_name.size(); // number of parameters
   // must be first
   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for(unsigned int i = 0; i < Npar; i++) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix  parameters
//   fitter.Config().ParSettings(0).SetValue(1.);
//   fitter.Config().ParSettings(0).Fix();
//   fitter.Config().ParSettings(2).SetValue(0.);
//   fitter.Config().ParSettings(2).Fix();

   // == Fit
   fitter.FitFCN(Npar, chi2_fit, nullptr, chi2_fit.Size(), true);

   // to obtain reliable errors:
//    fitter.CalculateHessErrors();
   fitter.CalculateMinosErrors();

   // == Fit result
   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double chi2   = res.Chi2();
   int ndf       = res.Ndf();
   const vector<double> par   = res.Parameters();
   const vector<double> er_par = res.Errors();

   // rescale BG:
   hst[13] = rescale_bg( hst[3], par.data(), Npol );
   hst[14] = rescale_bg( hst[4], par.data(), Npol );

   // ========================= draw results ========================
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   set_draw_opt(hst);
   SetHstFace(hst[0]);
   hst[0]->GetXaxis()->SetTitleOffset(1.1);
   hst[0]->GetYaxis()->SetTitleOffset(1.25);

   hst[0]->Draw("E"); // data

   // box to show not fitted region
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-10);
   box->SetLineWidth(2);
   double ymin = gPad->PixeltoY(0);
   double ymax = hst[0]->GetMaximum();
   box->DrawBox(ExMin,ymin,ExMax,ymax);
   hst[0]->Draw("E, SAME"); // data

   hst[7]=(TH1D*)hst[1]->Clone("SumBG");
   hst[7]->Add(hst[13]);
   hst[7]->Add(hst[14]);
   hst[7]->Draw("SAME,HIST"); // SumBG

   hst[8]=(TH1D*)hst[2]->Clone("SumAll");
   hst[8]->Add(hst[7]);
   hst[8]->SetLineColor(kRed+1);
   hst[8]->Draw("SAME,HIST"); // Sum

   hst[2]->SetLineStyle(kDashed);
   hst[2]->Draw("SAME,HIST"); // Signal

//    hst[14]->Draw("SAME,HIST"); // Rescaled Bg

   TLegend* leg = new TLegend(0.55,0.65,0.89,0.89);
//    leg->SetHeader("#bf{Recoil Mass of #pi^{#plus}#pi^{#minus}}", "C");
   leg->AddEntry(hst[0],
         (string("Data ")+to_string(date)).c_str(), "EP");

   leg->AddEntry(hst[7],
         Form("#color[%i]{sum of backgrounds}",
         hst[7]->GetLineColor() ),"L");

   leg->AddEntry(hst[2],
         Form("#color[%i]{MC #pi^{#plus}#pi^{#minus}J/#Psi "
         "correct #pi^{#plus}#pi^{#minus}}",
         hst[2]->GetLineColor() ),"L");

   leg->AddEntry(hst[8],
         Form("#color[%i]{MCsignal + bgs}",
         hst[8]->GetLineColor() ),"L");

//    leg->AddEntry(hst[14],
//          Form("#color[%i]{MC bg non #pi^{#plus}#pi^{#minus}J/#Psi}",
//          hst[14]->GetLineColor() ),"L");

   leg->AddEntry(box, "excluded from fit","F");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.11,0.65,0.44,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText( Form("#bf{Fit in [%.2f,%.2f] & [%.2f,%.2f] }",
                     3.0,ExMin,ExMax,3.2));
   pt->AddText( Form("#chi^{2}/ndf = %.1f / %i",chi2,ndf) );
   for(unsigned int i = 0; i < par.size(); i++) {
      pt->AddText( Form("%s =       %.4f #pm %.4f",
                   par_name[i].c_str(),par[i],er_par[i]) );
   }
   pt->Draw();

   gPad->RedrawAxis();

   c1->Update();
   string pdf = string("Mrec") + to_string(date) + "_fsb.pdf";
   c1->Print(pdf.c_str());

   // print for E in [Emin,Emax]
   print_Numbers(hst,3.092,3.102); // see cuts.h !
   print_Numbers(hst,3.055,3.145); // BAM-42
}

// {{{1 Draw Mrec
//--------------------------------------------------------------------
void MrecDraw(int date, bool zoom=false) {
//--------------------------------------------------------------------
   TH1D* hst[10];
   string pdf = string("Mrec") + to_string(date) +
      ((zoom) ? "_zoom.pdf" : ".pdf");
   if ( date==2009 ) {
      fill_hist2009(hst);
      hst[0]->SetMaximum(7.e6);
      hst[0]->SetMinimum(3.e4);
      if( zoom ) {
         hst[0]->SetMaximum(0.65e6);
         hst[0]->SetMinimum(0.);
         hst[0]->GetYaxis()->SetMaxDigits(3);
      }
   } else if(date==2012) {
      fill_hist2012(hst);
      hst[0]->SetMaximum(2.e7);
      hst[0]->SetMinimum(1.e5);
      if( zoom ) {
         hst[0]->SetMaximum(1.9e6);
         hst[0]->SetMinimum(0.);
         hst[0]->GetYaxis()->SetMaxDigits(3);
      }
   }

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();
   if ( !zoom ) {
      gPad->SetLogy();
   }

   set_draw_opt(hst);
   SetHstFace(hst[0]);
   hst[0]->GetXaxis()->SetTitleOffset(1.1);
   hst[0]->GetYaxis()->SetTitleOffset(1.25);

   hst[0]->Draw("E"); // data

   hst[5]=(TH1D*)hst[1]->Clone("Cont10");
   hst[5]->Scale(10.);
   hst[5]->Draw("SAME,HIST"); // Cont

   hst[3]->Draw("SAME,HIST"); // Bg1
   hst[4]->Draw("SAME,HIST"); // Bg2

   hst[6]=(TH1D*)hst[2]->Clone("SumAll");
   hst[6]->Add(hst[1]);
   hst[6]->Add(hst[3]);
   hst[6]->Add(hst[4]);
   hst[6]->SetLineColor(kRed+1);
   hst[6]->Draw("SAME,HIST"); // Sum

   hst[2]->SetLineStyle(kDashed);
   hst[2]->Draw("SAME,HIST"); // Signal

   TLegend* leg = new TLegend(0.55,0.60,0.89,0.89);
//    leg->SetHeader("Recoil Mass of #pi^{#plus}#pi^{#minus}", "C");
   leg->AddEntry(hst[0],
         (string("Data ")+to_string(date)).c_str(), "EP");

   leg->AddEntry(hst[6],
         Form("#color[%i]{MCsig + MCbg + Cont.}",
         hst[6]->GetLineColor() ),"L");

   leg->AddEntry(hst[2],
         Form("#color[%i]{MC #pi^{#plus}#pi^{#minus}J/#Psi "
         "correct #pi^{#plus}#pi^{#minus}}",
         hst[2]->GetLineColor() ),"L");

   leg->AddEntry(hst[3],
         Form("#color[%i]{MC #pi^{#plus}#pi^{#minus}J/#Psi "
         "wrong #pi^{#plus}#pi^{#minus}}",
         hst[3]->GetLineColor() ),"L");

   leg->AddEntry(hst[4],
         Form("#color[%i]{MC bg non #pi^{#plus}#pi^{#minus}J/#Psi}",
         hst[4]->GetLineColor() ),"L");

   leg->AddEntry(hst[5],
         Form("#color[%i]{Continuum bg #times 10}",
         hst[5]->GetLineColor() ),"L");
   leg->Draw();
   c1->Update();
   c1->Print(pdf.c_str());
}

// {{{1 Main
//--------------------------------------------------------------------
void MrecFitSB() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(112);
   gStyle->SetStatFont(62);
   gStyle->SetLegendFont(42);

   int date=2009;
//    int date=2012;

   bool zoom = true;
//    MrecDraw(date, !zoom);
//    MrecDraw(date,  zoom);

   DoFitSB(date);
}
