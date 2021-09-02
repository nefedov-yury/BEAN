// FOR MC WITHOUT HELIX CORRECTIONS
// - fit data by convolution of MC inclusive signal
//   with Gauss and scaling of MC inclusive background;
// - calculate number of J/Psi in Psi'->J/Psi pi+ pi- decay
//      -> Mrec_YEAR_fit.pdf
//

// for cuts.h:
#define CONSTANTS_ONLY

void print_Numbers(TH1D* hst[], double Emin, double Emax);
void print_eff(int date, TH1D* hst[],
      const vector<double>& par, const vector<double>& er_par,
      double Emin, double Emax);

//----------------------------------------------------------------------
// GLOBAL:
// name of folder with root files
static const string dir("prod-11/");

//----------------------------------------------------------------------
constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

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
// Fill histograms
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
TH1D* fill_mrec(string fname, string hname, int type=0, double shift=0.) {
//-------------------------------------------------------------------------
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

   TH1D* hst = new TH1D(hname.c_str(),
                        ";M_{rec}(#pi^{+}#pi^{-}), GeV/c^{2}"
//                         ";Entries/0.001GeV/c^{2}",
                        ";Entries/0.0001, GeV/c^{2}",
//                         200,3.0,3.2 // 1bin = 1MeV
                        2000,3.0,3.2 // 1bin = 0.1MeV
                       );
   hst->Sumw2(true);

   string dr = string("Mrec>>") + hname;
   TCut cut;
   if ( type == 1 ) {
      dr = string("Mrs>>") + hname;
      cut = TCut("MrsW*(dec==64)");
   } else if ( type == 2 ) {
      cut = TCut("dec==64");
   } else if ( type == 3 ) {
      cut = TCut("dec!=64");
//    } else if ( type == 10 ) {
//       dr = string("Mrs>>") + hname;
//       cut = TCut("dec==64");
   }

   // only for MC: small energy shift
   if ( type > 0 && fabs(shift) > 1e-5 ) {
      dr = to_string(shift) + string("+") + dr;
      // to_string() converts only 6 digits after the point - check it!
      cout << " fill_mrec::INFO dr= " << dr << endl;
   }

   nt1->Draw(dr.c_str(),cut,"goff");

   return hst;
}

//-------------------------------------------------------------------------
void fill_hist2009(TH1D* hst[]) {
//-------------------------------------------------------------------------
#include "cuts.h"

   // optimization of Mrec shift(data-MC)
//    const double shift = 0.000221; // 0.221MeV chi2 = 3391.4 ***
                                  // 0.220MeV chi2 = 3458.7
                                  // 0.222MeV chi2 = 3414.9

   const double shift = 0.000030; // after helix corrections
                                  //   MeV      chi^2
                                  // 0.0        106815
                                  // +0.040     105639
                                  // +0.030     105504 ***
                                  // +0.020     105643

   vector<string> hname {
      "data09", "data3650_09",
      "mc09sig", "mc09bg1", "mc09bg2",
      "mc09sig10"
   };

   // check cache file
   const string cachef = dir + string("MrecFit_2009.root");
   const bool use_cache = true;
   if ( use_cache ) {
      TFile* froot = TFile::Open(cachef.c_str(),"READ");
      if ( froot != 0 ) {
         for ( unsigned int i = 0; i < hname.size(); ++i ) {
            hst[i]=(TH1D*)froot->Get(hname[i].c_str());
         }
         return;
      }
   }

   hst[0]=fill_mrec("data_09psip_all.root", hname[0]);

   hst[1]=fill_mrec("data_3650_all.root", hname[1]);
   hst[1]->Scale(Cto09);

   hst[2]=fill_mrec("mcinc_09psip_all.root",hname[2],1,shift);
   hst[5] = (TH1D*)hst[2]->Clone( hname[5].c_str() );
   hst[2]->Scale(Ito09);

   hst[3]=fill_mrec("mcinc_09psip_all.root",hname[3],2,shift);
   hst[3]->Scale(Ito09);

   hst[4]=fill_mrec("mcinc_09psip_all.root",hname[4],3,shift);
   hst[4]->Scale(Ito09);

   // save histos in cache file
   if ( use_cache ) {
      TFile* froot = TFile::Open(cachef.c_str(),"NEW");
      for ( unsigned int i = 0; i < hname.size(); ++i ) {
         hst[i]->Write();
      }
      froot->Close();
   }
}

//-------------------------------------------------------------------------
void fill_hist2012(TH1D* hst[]) {
//-------------------------------------------------------------------------
#include "cuts.h"

   // optimization of Mrec shift(data-MC)
//    const double shift = 0.000381; // 0.381MeV chi2=6980.65 ***
                                  // 0.380MeV chi2=7054.27
                                  // 0.382MeV chi2=7179.46

   const double shift =-0.000160; // after helix corrections
                                  //   MeV      chi^2
                                  // 0.0        278971
                                  // -0.160     109815 <- !
                                  // -0.170     106767
                                  // -0.180     104651
                                  // -0.190     103430
   vector<string> hname {
      "data12", "data3650_12",
      "mc12sig", "mc12bg1", "mc12bg2",
      "mc12sig10"
   };

   // check cache file
   const string cachef =  dir + string("MrecFit_2012.root");
   const bool use_cache = true;
   if ( use_cache ) {
      TFile* froot = TFile::Open(cachef.c_str(),"READ");
      if ( froot != 0 ) {
         for ( unsigned int i = 0; i < hname.size(); ++i ) {
            hst[i]=(TH1D*)froot->Get(hname[i].c_str());
         }
         return;
      }
   }

   hst[0]=fill_mrec("data_12psip_all.root", hname[0]);

   hst[1]=fill_mrec("data_3650_all.root", hname[1]);
   hst[1]->Scale(Cto12);

   hst[2]=fill_mrec("mcinc_12psip_all.root",hname[2],1,shift);
   hst[5] = (TH1D*)hst[2]->Clone( hname[5].c_str() );
   hst[2]->Scale(Ito12);

   hst[3]=fill_mrec("mcinc_12psip_all.root",hname[3],2,shift);
   hst[3]->Scale(Ito12);

   hst[4]=fill_mrec("mcinc_12psip_all.root",hname[4],3,shift);
   hst[4]->Scale(Ito12);

   // save histos in cache file
   if ( use_cache ) {
      TFile* froot = TFile::Open(cachef.c_str(),"NEW");
      for ( unsigned int i = 0; i < hname.size(); ++i ) {
         hst[i]->Write();
      }
      froot->Close();
   }
}

//-------------------------------------------------------------------------
void set_draw_opt(TH1D* hst[]) {
//-------------------------------------------------------------------------
// data
   hst[0]->SetMarkerStyle(20);
   hst[0]->SetMarkerSize(0.7);
// data 3650
   hst[1]->SetLineColor(kBlue+1);
   hst[1]->SetLineWidth(2);
// MC signal
   hst[2]->SetLineColor(kGreen+3);
   hst[2]->SetLineWidth(1);
// MC bg from pi+pi-J/Psi
   hst[3]->SetLineColor(kBlue+3);
   hst[3]->SetLineWidth(2);
// MC bg not pi+pi-J/Psi
   hst[4]->SetLineColor(kMagenta+1);
   hst[4]->SetLineWidth(2);
}


//-------------------------------------------------------------------------
// FIT
//-------------------------------------------------------------------------

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

//-------------------------------------------------------------------------
vector<double> Vsmear(double wbin, const double p[]) {
//-------------------------------------------------------------------------
   // prepare weight vector for Gauss convolution with signal
   // wbin - width of bin

   vector<double> gw;

   // decode parameters
   double sig1 = p[0];

   if ( sig1 < 0 ) {
      return gw;
   }

   gw.reserve(100);
   double wstep1 = wbin/(sqrt(2.)*sig1);
   double x01 = 0.5*wstep1;
   double erf01 = erf(-x01);
   double wsum = 0;
   for ( double x1 = x01; ; x1 += wstep1 ) {
      double erf1 = erf(x1);
      double w = 0.5*(erf1-erf01);
      wsum += 2*w ;
      if ( w < 1e-9 ) {
         break;
      }
      gw.push_back(w);
      erf01=erf1;
   }

   return gw;
}

//-------------------------------------------------------------------------
vector<double> Vsmear2(double wbin, const double p[]) {
//-------------------------------------------------------------------------
   // prepare weight vector for Gauss convolution with signal
   // wbin - width of bin

   vector<double> gw;

   // decode parameters
   double sig1 = p[0];
   double sig2 = p[1];
   double part = p[2]; // must be in [0,1]

   if ( sig1 < 0 || sig2 < 0 || part < 0 || part > 1 ) {
      return gw;
   }

   gw.reserve(100);
   double wstep1 = wbin/(sqrt(2.)*sig1);
   double wstep2 = wbin/(sqrt(2.)*sig2);
   double x01 = 0.5*wstep1;
   double x02 = 0.5*wstep2;
   double erf01 = erf(-x01);
   double erf02 = erf(-x02);
   double wsum = 0;
   for ( double x1 = x01, x2 = x02; ; x1+=wstep1,x2+=wstep2 ) {
      double erf1 = erf(x1);
      double erf2 = erf(x2);
      double w1 = 0.5*(erf1-erf01);
      double w2 = 0.5*(erf2-erf02);
      double w = w1*part + w2*(1-part);
      wsum += 2*w ;
      if ( w < 1e-9 ) {
         break;
      }
      gw.push_back(w);
      erf01=erf1;
      erf02=erf2;
   }

//    printf(" DEBUG: gw.size= %zd, sig1= %.3g, sig2= %.3g, part= %.3g",
//           gw.size(),sig1,sig2,part);
//    if ( gw.size() > 0 )  {
//       printf(" wsum-gw[0] = %.9f\n", wsum-gw[0]);
//    } else {
//       printf("\n");
//    }

   return gw;
}

//-------------------------------------------------------------------------
class myChi2 {
   public:
      // ctor:
      //-------------------------------------------------------------------
      myChi2(TH1D* hst[],double EMin,double EMax) : Emin(EMin),Emax(EMax) {
      //-------------------------------------------------------------------
         // ExMin/Max -> region to fit signal

         Ndata = 0;
         wbin = hst[0]->GetBinWidth(1);

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
            en.push_back(x);
            if(x >= Emin && x <= Emax) {
               Ndata++;
            }

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
      //-------------------------------------------------------------------
      unsigned int Size() const {
      //-------------------------------------------------------------------
         return Ndata;
      }

      //-------------------------------------------------------------------
      void SetNN(int ng,int np) {
      //-------------------------------------------------------------------
         nsm = ng;
         npol = np;
      }

      // chi^2 - function
      // the signature of this operator() MUST be exactly this:
      //-------------------------------------------------------------------
      double operator() (const double* p) {
      //-------------------------------------------------------------------
         static const double maxResValue = DBL_MAX / 1000;

         int n = Size();

         double Sc = p[0];  // scale constant for signal
         double Sc2 = SQ(Sc);
         vector<double> gw;
         int nsh = 0;
         // smoothing by Gaussian
         if ( nsm == 1 ) {
            gw = Vsmear(wbin, &p[1]); // one Gaussian
            nsh = 2;
         } else if ( nsm == 2 ) {
            gw = Vsmear2(wbin, &p[1]); // two Gaussian
            nsh = 4;
         }
         int ngw = gw.size();
         if ( ngw == 0 ) {
            return maxResValue;
         }

         // calculate "new" MC signal by convolution with Gauss
         // ATTENTION: errors must be multiplied on sqrt(weight)
         vector<double> sig_n(n), er_sig_n(n);
         for ( int i = 0; i < n; i++ ) {
            double d = sig[i];
            double e = er_sig[i];
            for ( int j = 0; j < ngw; j++ ) {
               double dd = d*gw[j];
               double ed = e*gw[j];

               if ( i+j < n ) { // to forward
                  sig_n[i+j] += dd;
                  er_sig_n[i+j] += ed;
               } else {
                  sig_n[i] += dd;
                  er_sig_n[i] += ed;
               }

               if( j==0 ) {
                  continue;
               }

               if ( i-j >=0 ) { // to backward
                  sig_n[i-j] += dd;
                  er_sig_n[i-j] += ed;
               } else {
                  sig_n[i] += dd;
                  er_sig_n[i] += ed;
               }
            }
         }

         double chi2 = 0;
         for ( int i = 0; i < n; i++ ) {
            double x = en[i];
            if ( x < Emin || x > Emax ) {
               continue;
            }

            // re-weight MC background by polynomial
            double xm = x-3.1; // subtract midpoint
            double bgscl = PolN(xm, &p[nsh], npol); // p[nsh,nsh+npol]
            double bgscl2 = bgscl*bgscl;

            double dif = data[i] -
                         ( Sc*sig_n[i] + cont[i] +
                           bgscl * (bg[i] + bgn[i]) );
//                            bgscl * bg[i] + bgn[i] );
//                            bgscl * bgn[i] + bg[i] );
            double er2 = er_data[i] +
                         Sc2*er_sig_n[i] +
                         er_cont[i] +
                         bgscl2 * (er_bg[i] + er_bgn[i]);
//                          bgscl2*er_bg[i] + er_bgn[i];
//                          bgscl2*er_bgn[i] + er_bg[i];

            double resval = dif*dif / er2;
            // avoid inifinity or nan in chi2
            if ( resval < maxResValue ) {
               chi2 += resval;
            } else {
               return maxResValue;
            }
         } // end of for()

         return chi2;
      }

   private:
      double Emin, Emax;
      int Ndata;

      double wbin;        // bin width
      vector<double> en;  // energy
      vector<double> data, er_data;
      vector<double> cont, er_cont;
      vector<double> sig, er_sig;
      vector<double> bg, er_bg;
      vector<double> bgn, er_bgn;

      int npol = 3;
      int nsm = 2;
};
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
TH1D* smear_sig(const TH1D* hst, const double* p) {
//-------------------------------------------------------------------------
   // ignore over and underflow bins

   // create copy of hst
   string new_name = string(hst->GetName()) + string("_smear");
   TH1D* hst_n = (TH1D*)hst->Clone(new_name.c_str());
//   return hst_n;

   int n = hst->GetNbinsX();
   double wbin = hst->GetBinWidth(1);

   double Sc = p[0];  // scale constant for signal

   vector<double> gw = Vsmear(wbin, &p[1]); // p[1-3]
//    vector<double> gw = Vsmear2(wbin, &p[1]); // p[1-3]
   int ngw = gw.size();
   if ( ngw == 0 ) {
      cout << " ERROR in smear_sig: ngw= " << ngw << endl;
      return hst_n;
   }

   // calculate "new" MC signal
   // ATTENTION: errors must be multiplied on sqrt(weight)
   vector<double> sig_n(n), er_sig_n(n);
   for ( int i = 0; i < n; i++ ) {
      double d = hst->GetBinContent(1+i);
      double e = SQ(hst->GetBinError(1+i));
      for ( int j = 0; j < ngw; j++ ) {
         double dd = d*gw[j];
         double ed = e*gw[j];

         if ( i+j < n ) { // to forward
            sig_n[i+j] += dd;
            er_sig_n[i+j] += ed;
         } else {
            sig_n[i] += dd;
            er_sig_n[i] += ed;
         }

         if ( j==0 ) {
            continue;
         }

         if ( i-j >=0 ) { // to backward
            sig_n[i-j] += dd;
            er_sig_n[i-j] += ed;
         } else {
            sig_n[i] += dd;
            er_sig_n[i] += ed;
         }
      }
   }

   for ( int i = 0; i < n; i++ ) {
      hst_n->SetBinContent(i+1,Sc*sig_n[i]);
      hst_n->SetBinError(i+1,Sc*sqrt(er_sig_n[i]));
   }

   return hst_n;
}

//-------------------------------------------------------------------------
void rescale_bg(TH1D* hst, const double* p, int npol) {
//-------------------------------------------------------------------------
   // rescale ignoring over and underflow bins

   int nbins = hst->GetNbinsX();
   for ( int n = 1; n <= nbins; ++n ) {
      double x = hst->GetBinCenter(n) - 3.1;
      double w = PolN(x,p,npol);
      double d = hst->GetBinContent(n);
      double e = hst->GetBinError(n);
      hst->SetBinContent(n,d*w);
      hst->SetBinError(n,e*w);
   }
}

//-------------------------------------------------------------------------
void DoFit(int date) {
//-------------------------------------------------------------------------
   TH1D* hst[20];
   if ( date==2009 ) {
      fill_hist2009(hst);
      hst[0]->SetMaximum(7.e5);
      hst[0]->SetMinimum(0.99e4);
   } else if(date==2012) {
      fill_hist2012(hst);
      hst[0]->SetMaximum(2.e6);
      hst[0]->SetMinimum(3.e4);
   }

   // fit signal in region (ExMin; ExMax)
   double ExMin = 3.01;
   double ExMax = 3.19;
   myChi2 chi2_fit(hst,ExMin,ExMax);

   // ========================= Fit with ROOT =========================
   // == fit configuration
   ROOT::Fit::Fitter fitter;
   // set parameters of fitter: (Minuit,Minuit2,Fumili,GSLMultiFit...)
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");
   fitter.Config().MinimizerOptions().SetPrintLevel(3);
   // print the default minimizer option values
   ROOT::Math::MinimizerOptions min_opt;
   min_opt.Print();

   const unsigned int Nsm = 1;  // number of gausses for smoothing
   const unsigned int Npol = 3; // order of polynomial (base=3)
   chi2_fit.SetNN(Nsm,Npol);

   // set names and start values for parameters
   vector<string> par_name {"S ","#sigma"};
//    vector<string> par_name {"S ","#sigma_{1}","#sigma_{2}","f  "};
   for ( int i = 0; i <= Npol; ++i ) {
      par_name.push_back( string("p") + to_string(i) );
   }

   vector<double> par_ini;
   if ( date == 2009 ) {
//       par_ini = { 0.9825,6.41e-4,3.22e-3,0.866,  1.0394,0.11,-1.4,-16. };
      // helix corr
//       par_ini = { 1.0,1.e-9,  1.0261,0.036,1.66,4.4 };
      par_ini = { 1.0,1.e-4,  1.0252,0.054,1.90,1.7 };
   } else if ( date == 2012 ) {
//       par_ini = { 0.9958,7.27e-4,3.42e-3,0.854,  1.0543,0.181,-1.41,-31. };
      // helix corr
//       par_ini = { 1.,1e-9,1.e-9,0.5,  1.0402,0.021,1.51,3.1 };
   }

   const unsigned int Npar = par_name.size(); // number of parameters
   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for(unsigned int i = 0; i < Npar; i++) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // Signal scale:
   fitter.Config().ParSettings(0).SetStepSize(1e-3);  // S
   fitter.Config().ParSettings(1).SetLimits(0.,1.); // sig
   fitter.Config().ParSettings(1).SetStepSize(1e-3);
//    fitter.Config().ParSettings(1).Fix();
//    fitter.Config().ParSettings(2).SetLimits(0.,1e-2); // sig2
//    fitter.Config().ParSettings(2).SetStepSize(1e-4);
//    fitter.Config().ParSettings(3).SetLimits(0.,1.);    // fraction: [0,1]

//    for ( int i = 1; i < 2*Nsm; ++i ) {      // fix parameters of signal
//    for ( int i = 2*Nsm; i <= 2*Nsm+Npol; ++i ) { // fix pars of background
//       fitter.Config().ParSettings(i).Fix();
//    }

   // == Fit
   fitter.FitFCN(Npar, chi2_fit, nullptr, chi2_fit.Size(), true);

   // to obtain reliable errors:
//   fitter.CalculateHessErrors();
//   fitter.CalculateMinosErrors();

   // == Fit result
   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double chi2   = res.Chi2();
   int ndf       = res.Ndf();
   const vector<double> par    = res.Parameters();
   const vector<double> er_par = res.Errors();

   // smear MC signal
   hst[12] = smear_sig(hst[2],par.data()); // 2 -> 12 Mc(sig) new

   // rescale BG:
   hst[13]=(TH1D*)hst[3]->Clone("RescaledBG");
   rescale_bg( hst[13], &par[2*Nsm], Npol );
   hst[14]=(TH1D*)hst[4]->Clone("RescaledBGn");
   rescale_bg( hst[14], &par[2*Nsm], Npol );

   // ========================= draw results ========================
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   string pdf = string("Mrec")+to_string(date)+string("_fit.pdf");

   set_draw_opt(hst);
   SetHstFace(hst[0]);
//    hst[0]->GetXaxis()->SetTitleOffset(0.9);
   hst[0]->GetYaxis()->SetTitleOffset(1.2);

   hst[0]->Draw("E"); // data

   hst[7]=(TH1D*)hst[1]->Clone("SumBG");
   hst[7]->Add(hst[13]);
   hst[7]->Add(hst[14]);
   hst[7]->Draw("SAME,HIST"); // SumBG

   hst[8]=(TH1D*)hst[12]->Clone("SumAll");
   hst[8]->Add(hst[7]);
   hst[8]->SetLineColor(kRed+1);
   hst[8]->Draw("SAME,HIST"); // Sum

   hst[12]->SetLineStyle(kDashed);
   hst[12]->Draw("SAME,HIST"); // Smeared signal

//    hst[14]->Draw("SAME,HIST"); // Rescaled Bg not pi+pi-J/Psi

//    TLegend* leg = new TLegend(0.60,0.67,0.89,0.89);
//    leg->SetHeader("#bf{Recoil Mass of #pi^{+}#pi^{-}}", "C");
   TLegend* leg = new TLegend(0.60,0.70,0.89,0.89);
   leg->AddEntry(hst[0],
         (string("Data ")+to_string(date)).c_str(), "EP");
   leg->AddEntry(hst[8], Form("#color[%i]{Sig + Bkg + Cont}",
                              hst[8]->GetLineColor() ),"L");
   leg->AddEntry(hst[12],Form("#color[%i]{Sig}",
                              hst[12]->GetLineColor() ),"L");
   leg->AddEntry(hst[7], Form("#color[%i]{Bkg + Cont}",
                              hst[7]->GetLineColor() ),"L");
//    leg->AddEntry(hst[14],Form("#color[%i]{MC bg non #pi^{+}#pi^{-}J/#Psi}",
//                               hst[14]->GetLineColor() ),"L");
   leg->Draw();


//    TPaveText* pt = new TPaveText(0.11,0.60,0.42,0.89,"NDC");
   TPaveText* pt = new TPaveText(0.11,0.57,0.37,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->AddText(
      Form("#bf{Fit in [%.2f,%.2f]}",ExMin,ExMax)
   );
   pt->AddText(
      Form("#chi^{2}/ndf  %.1f / %i",chi2,ndf)
   );
   for ( unsigned int i = 0; i < par.size(); i++ ) {
      if ( i < 3 ) {
         pt->AddText( Form("%s  %.4f #pm %.4f",
                      par_name[i].c_str(),par[i],er_par[i]) );
      } else {
         pt->AddText( Form("%s  %.3f #pm %.3f",
                      par_name[i].c_str(),par[i],er_par[i]) );
      }
   }
   pt->Draw();
   gPad->RedrawAxis();

   c1->Update();
   c1->Print(pdf.c_str());

   // numbers for E in [Emin,Emax]
   print_Numbers(hst,3.092,3.102);// see cuts.h !
   print_Numbers(hst,3.055,3.145);// BAM-42

   print_eff(date,hst,par,er_par,3.092,3.102);
   print_eff(date,hst,par,er_par,3.055,3.145);
}


//-------------------------------------------------------------------------
// Draw, print ...
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void print_Numbers(TH1D* hst[], double Emin, double Emax) {
//-------------------------------------------------------------------------
   // ATTENTION: here I assume that re-weighting was done and
   // [0] - data
   // [7] - sum of backgrounds
   // [12] - MC(sig)

   bool fbin=true;

   double ndata   = 0;
   double er_data = 0;
   double nsig    = 0;
   double er_sig  = 0;
   double nbg     = 0;
   double er_bg   = 0;

   int nbins = hst[0]->GetNbinsX();
   for(int n = 1; n <= nbins; ++n) {
      double x = hst[0]->GetBinCenter(n);
      if( x < Emin ) {
         continue;
      }
      if( fbin ) {
         fbin=false;
//          cout << " first bin# " << n << " E= " << x
//               << " (Emin= " << Emin << ")" << endl;
      }
      if( x > Emax ) {
//          cout << " last bin# " << n-1 << " E= " << hst[0]->GetBinCenter(n-1)
//               << " (Emax= " << Emax << ")" << endl;
         break;
      }
      ndata   += hst[0]->GetBinContent(n);
      er_data += SQ(hst[0]->GetBinError(n));
      nsig    += hst[12]->GetBinContent(n);
      er_sig  += SQ(hst[12]->GetBinError(n));
      nbg     += hst[7]->GetBinContent(n);
      er_bg   += SQ(hst[7]->GetBinError(n));
   }

   double n1   = ndata-nbg;
   double er1  = sqrt(er_data+er_bg);
   er_data     = sqrt(er_data);
   er_sig      = sqrt(er_sig);
   er_bg       = sqrt(er_bg);

   printf("\n");
   printf(" Number of events in [%.3f, %.3f]:\n", Emin,Emax);
   printf(" Data:            %9.0f +/- %4.0f\n", ndata,er_data);
   printf(" Sum(bg):         %9.0f +/- %4.0f\n", nbg,er_bg);
   printf(" NPsip(data-bg):  %9.0f +/- %4.0f\n", n1,er1);
   printf(" NPsip(MCsig):    %9.0f +/- %4.0f\n", nsig,er_sig);
   double diff =  100*fabs(n1-nsig)/nsig;
   cout << " " << string(50,'-') << endl;
   printf(" dif= %.2f%%\n",diff);
}

//-------------------------------------------------------------------------
void print_eff(int date, TH1D* hst[],
      const vector<double>& par, const vector<double>& er_par,
      double Emin, double Emax) {
//-------------------------------------------------------------------------
// calculate efficiency of the Psi' -> pi+ pi- J/Psi selection
   // ATTENTION: here I assume:
   // [5] - MC(sig) without smearing and without normalization

   static string fnames[] = {
      "mcinc_09psip_all.root",
      "mcinc_12psip_all.root"
   };

   string fname = dir + fnames[(date > 2010)];
//    cout << endl << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");

   // initial number of pi+ pi- J/Psi
   TH1D* MCdec = (TH1D*)gROOT->FindObject("mc_dec0");
   if ( !MCdec ) {
      cout << " can not find mc_dec0" << endl;
      exit(0);
   }
   double Nppj = MCdec->GetBinContent(65);
   cout << " Initial number of pi+ pi- J/Psi (dec# "
        << MCdec->GetBinCenter(65) << ")= " << Nppj << endl;

   // to study systematic vary parameters of fit function
   // on one sigma:
   vector<double> eff(7,0.), err(7,0);
   double sys_err = 0;
   for (int it = 0; it < 7; ++it ) {
      vector<double> par1(par);
      par1[0] = 1.; // absolute scaling should be 1
      switch(it) {
         case 0: break;
         case 1: par1[1] = par[1] - er_par[1]; break;
         case 2: par1[1] = par[1] + er_par[1]; break;
         case 3: par1[2] = par[2] - er_par[2]; break;
         case 4: par1[2] = par[2] + er_par[2]; break;
         case 5: par1[3] = par[3] - er_par[3]; break;
         case 6: par1[3] = par[3] + er_par[3]; break;
      }

//       TH1D* hMC = hst[5]; // debug
      TH1D* hMC = smear_sig(hst[5],par1.data()); // 5(no MrecW) -> hMC
      double Nmc = 0.;
      double er_Nmc = 0.;

      int nbins = hMC->GetNbinsX();
      for ( int n = 1; n <= nbins; ++n ) {
         double x = hMC->GetBinCenter(n);
         if( x < Emin ) { continue; }
         if( x > Emax ) { break; }
         Nmc    += hMC->GetBinContent(n);
         er_Nmc += SQ(hMC->GetBinError(n));
      }
      er_Nmc = sqrt(er_Nmc);

      eff[it] = Nmc / Nppj;
      err[it] = eff[it] * sqrt(SQ(er_Nmc/Nmc) + 1./Nppj);
//       printf(" it=%i eff: %.5f +/- %.5f\n",it,eff[it],err[it]);

      if ( it > 0 ) {
         double dif = fabs(eff[it]-eff[0]);
         if( dif > sys_err ) sys_err = dif;
      }
   }

   printf("\n");
   printf(" Efficiency in [%.3f, %.3f]:\n", Emin,Emax);
   printf(" eff(pi+pi-J/Psi): %.3f +/- %.3f +/- %.3f %%\n",
         100*eff[0],100*err[0],100*sys_err);
}

//-------------------------------------------------------------------------
void MrecFit() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(112);
//    gStyle->SetLegendTextSize(0.03);
   gStyle->SetLegendFont(42);

   int date=2009;
//    int date=2012;

   DoFit(date);
}

