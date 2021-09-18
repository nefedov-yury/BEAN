// FOR MC AFTER HELIX CORRECTIONS
// - data fitting by correcting of MC signal and scaling of MC
//   background
//   -> Mrec_YEAR_fitMODEL.pdf
//   -> diffYEAR_MODEL.pdf for delta(data-MC)
// - calculate number of J/Psi in Psi' -> J/Psi pi+ pi- decay
//   and efficiency of selection;
// - estimate systematic associated with a fit model

#include "ReWeightTrkPid_11.h"

// {{{1 helper functions and constants
//--------------------------------------------------------------------
// GLOBAL: name of folder with root files
static const string dir("prod-11/");

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
   hst[2]->SetLineWidth(1);
// MC bg from pi+pi-J/Psi
   hst[3]->SetLineColor(kBlue+3);
   hst[3]->SetLineWidth(2);
// MC bg not pi+pi-J/Psi
   hst[4]->SetLineColor(kMagenta+1);
   hst[4]->SetLineWidth(2);
}

// {{{1 Fill histograms
//--------------------------------------------------------------------
TH1D* fill_mrec(string fname,string hname,int type=0,double shift=0.){
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
         ";Entries/0.0001 GeV/c^{2}",
         2000,3.0,3.2 // 1bin = 0.1MeV
         );
   hst->Sumw2(true);

   string dr = string("Mrec>>") + hname;
   TCut cut;
   if ( type == 1 ) {
//       dr = string("Mrs>>") + hname;
//       cut = TCut("MrsW*(dec==64)"); // MrsW pi+pi- corrections

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
         hst -> Fill(shift+Mrs,W);
      }
      return hst;
   } else if ( type == 2 ) {
      cut = TCut("dec==64");
   } else if ( type == 3 ) {
      cut = TCut("dec!=64");
   }

   // only for MC: small energy shift
   if ( type > 0 && fabs(shift) > 1e-5 ) {
      dr = to_string(shift) + string("+") + dr;
   }
   // to_string() converts only 6 digits after the point,
   // check output!
   cout << " fill_mrec::INFO dr= " << dr << endl;

   nt1->Draw(dr.c_str(),cut,"goff");

   return hst;
}

//--------------------------------------------------------------------
void fill_hist2009(TH1D* hst[]) {
//--------------------------------------------------------------------
#include "norm.h"

   // optimization of Mrec shift(data-MC) !after helix corrections!
   const double shift = 0.000026; //   100 iterations
                                  // MeV      chi^2(3M,fix-bg)
                                  //  0.        9214.39
                                  // +0.024     8139.61   6204.53
                                  // +0.025     8057.32   6122.76
                                  // +0.026**   8048.78   6113.98
                                  // +0.027     8126.97   6194.08

   vector<string> hname {
      "data09", "data3650_09",
      "mc09sig", "mc09bg1", "mc09bg2",
      "mc09sig10"
   };

   // check cache file
   string cachef = dir + string("MrecFit_2009");
   const bool use_cache = true;
   if ( use_cache ) {
      cachef += (shift<0) ? "_m" : "_";
      cachef += to_string(int(fabs(shift*1e6))) + ".root";
      cout << " INFO: cachef= " << cachef << endl;
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

//--------------------------------------------------------------------
void fill_hist2012(TH1D* hst[]) {
//--------------------------------------------------------------------
#include "norm.h"

   // optimization of Mrec shift(data-MC) !after helix corrections!
   const double shift =-0.000201; //  100 iterations
                                  //   MeV    chi^2(M1,fix-bg)
                                  //  0.      221199
                                  // -0.160    49651.3 <data>=<MC>
                                  // -0.199    42361     33334.1
                                  // -0.200    42191.5   33157.6
                                  // -0.201**  42092.1   33055.3
                                  // -0.202    42250.4   33208.9
   vector<string> hname {
      "data12", "data3650_12",
      "mc12sig", "mc12bg1", "mc12bg2",
      "mc12sig10"
   };

   // check cache file
   string cachef =  dir + string("MrecFit_2012");
   const bool use_cache = true;
   if ( use_cache ) {
      cachef += (shift<0) ? "_m" : "_";
      cachef += to_string(int(fabs(shift*1e6))) + ".root";
      cout << " INFO: cachef= " << cachef << endl;
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

// {{{1 Function for fitting
//--------------------------------------------------------------------
double ChebN(int nch, double Xmin, double Xmax, double x,
      const double* p) {
//--------------------------------------------------------------------
   // calculate value of the Chebyshev polynomial of 'nch' order
   // [Xmin Xmax] is the range of polynomial orthogonality
   // https://en.wikipedia.org/wiki/Chebyshev_polynomials
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
double PolN(double x, const double* p, int n) {
//--------------------------------------------------------------------
   // calculate value of the polynomial by Horner's method:
   // here 'n' is the polynomial order, size of coefficients is p[n+1]
   double res = 0;
   for ( int i = n; i >= 0; --i ) {
      res = res*x + p[i];
   }
   return res;
}

// integrated normal distribution in bin-width of histogram
// expressed in terms of CDF (erf function)
//--------------------------------------------------------------------
vector<double> dCDFnorm(double wbin, double sig) {
//--------------------------------------------------------------------
   // return weight vector for convolution with Gauss(sigma = sig)
   // wbin - width of bin

   vector<double> gw; // helper vector
   if ( sig < 0 ) {
      return gw;
   }

   gw.reserve(100);
   double wstep = wbin/(sqrt(2.)*sig);
   double x0 = 0.5*wstep;
   double erf0 = erf(-x0);
   double wsum = 0;
   for ( double x1 = x0; ; x1 += wstep ) {
      double erf1 = erf(x1);
      double w = 0.5*(erf1-erf0);
      gw.push_back(w);
      wsum += 2*w ;
      if ( w < 1e-9 ) {
         break;
      }
      erf0 = erf1;
   }

//    printf(" DEBUG: gw.size= %zu, wsum-gw[0]= %.9f (must be 1)\n",
//             gw.size(), wsum-gw[0]);

   int ngw = gw.size();
   vector<double> res(2*ngw-1); // result
   int m = ngw-1; // midle point of res
   res[m] = gw[0];
   for ( int i = 1; i < ngw ; ++i ) {
      res[m-i] = res[m+i] = gw[i];
   }

   // test: shift on one bin to right -> the same result!
//    res.resize(2*ngw);
//    for ( int i = 2*ngw-1; i > 0 ; --i ) {
//       res[i] = res[i-1];
//    }
//    res[0] = 0;

   return res;
}

// integrated generalized normal distribution II in bin-width of
// histogram
// https://en.wikipedia.org/wiki/Generalized_normal_distribution
//--------------------------------------------------------------------
vector<double> dCDFGnorm(double wbin, double alpha, double kappa) {
//--------------------------------------------------------------------
// ksi == 0; alpha > 0;  kappa = 0 => normal distribution
   const double eps = 1e-9;

   vector<double> res;
   if ( alpha <= 0 ) {
      return res;
   }

   res.reserve(128);

   bool norm = fabs(kappa) < eps;
   auto cdf = [norm,alpha,kappa](double x) -> double {
      double y = (norm) ? x/alpha : -log(1-kappa*x/alpha)/kappa;
      return 0.5*( 1. + erf(y*M_SQRT1_2) );
   };

   double x0 = (norm) ? -6.5*alpha : alpha*(1-exp(6.5*kappa))/kappa;
   double cdf0 = cdf(x0);
   double wsum = 0;
   for ( double x1 = x0+wbin; ; x1 += wbin ) {
      double cdf1 = cdf(x1);
      double w = cdf1 - cdf0;
      res.push_back(w);
      wsum += w ;
      if ( x1 > 0 && w < 1e-9 ) {
         break;
      }
      cdf0 = cdf1;
   }

//    printf(" DEBUG: res.size=%zu, wsum=%.9f(must be 1), res[0]=%g\n",
//             res.size(), wsum, res[0]);

   return res;
}

// CDF for sum of two normal distributions
//--------------------------------------------------------------------
vector<double> dCDFnorm2(double wbin,double s1,double s2,double f) {
//--------------------------------------------------------------------
   // return weight vector for convolution with mixture of two Gauss
   // sigma = s1 & s2 and fraction of the second Gause is f
   // wbin - width of bin

   vector<double> gw; // helper vector
   if ( s1 < 0 || s2 < 0 || f < 0 || f > 1 ) {
      return gw;
   }

   gw.reserve(200);
   double wstep1 = wbin/(sqrt(2.)*s1);
   double wstep2 = wbin/(sqrt(2.)*s2);
   double x01 = 0.5*wstep1;
   double x02 = 0.5*wstep2;
   double erf01 = erf(-x01);
   double erf02 = erf(-x02);
   double wsum = 0;
   for ( double x1 = x01, x2 = x02; ; x1 += wstep1, x2 += wstep2 ) {
      double erf1 = erf(x1);
      double erf2 = erf(x2);
      double w1 = 0.5*(erf1-erf01);
      double w2 = 0.5*(erf2-erf02);
      double w = w1*(1-f) + w2*f;
      gw.push_back(w);
      wsum += 2*w ;
      if ( w < 1e-9 ) {
         break;
      }
      erf01 = erf1;
      erf02 = erf2;
   }

//    printf(" DEBUG: gw.size= %zu, wsum-gw[0]= %.9f (must be 1)\n",
//             gw.size(), wsum-gw[0]);

   int ngw = gw.size();
   vector<double> res(2*ngw-1); // result
   int m = ngw-1; // midle point of res
   res[m] = gw[0];
   for ( int i = 1; i < ngw ; ++i ) {
      res[m-i] = res[m+i] = gw[i];
   }

   return res;
}

// convolution
//--------------------------------------------------------------------
bool conv_vec( vector<double>& vec, vector<double>& er2,
               double wbin, double sig) {
//--------------------------------------------------------------------
   vector<double> res = dCDFnorm(wbin, sig);
   if ( res.empty() ) {
      return false;
   }
   int nres = res.size();
   int mres = (nres-1)/2; // midle point of res

   int n = vec.size();
   vector<double> vec_n(n), er2_n(n);
   for ( int i = 0; i < n; ++i ) {
      double d = vec[i];
      double e = er2[i];
//       int jmin = max(0,mres-i);
//       int jmax = min(nres,n+mres-i);
      for ( int j = 0; j < nres; ++j ) {
         int idx = i - mres + j;
         if ( idx < 0 ) continue; // j >= mres-i
         if ( idx > n-1 ) break;  // j < n+mres-i
         vec_n[idx] += d*res[j];
         er2_n[idx] += e*res[j];
      }
   }

   // replace vec & er2
   vec = vec_n;
   er2 = er2_n;
   return true;
}

// deconvolution
// ATTENTION: er2 is sqares of errors
//--------------------------------------------------------------------
bool decon_vec( vector<double>& vec, vector<double>& er2,
                double wbin, double sig) {
//--------------------------------------------------------------------
   vector<double> res = dCDFnorm(wbin, sig);
   if ( res.empty() ) {
      return false;
   }
   int n = vec.size();
   res.resize(n,0.);

   int niter = 100;   // number of iterations
   int nrep = 1;      // number of repetitions, for boosted dec.
   double boost = 1.; // boosting coefficient
   static TSpectrum* s = new TSpectrum();
   const char* dec_err =
   s -> Deconvolution(vec.data(), res.data(), n, niter,nrep,boost);
//    s -> DeconvolutionRL(vec.data(), res.data(), n, niter,nrep,boost);

   if ( dec_err ) {
      cerr << " ERROR in decon_vec:  " << endl;
      cerr << " vec.size= " << vec.size() << " wbin= " << wbin
           << " sig= " << sig << endl;
      cerr <<  dec_err << endl;
      return false;
   }

   // TODO: I do not know how to propagate er2 to newer errors
   // remove all negative values from vec
   for_each( vec.begin(),vec.end(), [](double& v){ v=max(0.,v); } );
   // and replace er2 with vec (er2 = errors^2)
   er2 = vec;

   return true;
}

//--------------------------------------------------------------------
bool deconG_vec( vector<double>& vec, vector<double>& er2,
                 double wbin, double sigma, double kappa ) {
//--------------------------------------------------------------------
   vector<double> res = dCDFGnorm(wbin, sigma, kappa);
   if ( res.empty() ) {
      return false;
   }
   int n = vec.size();
   res.resize(n,0.);

   int niter = 100;   // number of iterations
   int nrep = 1;      // number of repetitions, for boosted dec.
   double boost = 1.; // boosting coefficient
   static TSpectrum* s = new TSpectrum();
   const char* dec_err =
   s -> Deconvolution(vec.data(), res.data(), n, niter,nrep,boost);
//    s -> DeconvolutionRL(vec.data(), res.data(), n, niter,nrep,boost);

   if ( dec_err ) {
      cerr << " ERROR in deconG_vec:  " << endl;
      cerr << " vec.size= " << vec.size() << " wbin= " << wbin
           << " sig= " << sigma << " kappa= " << kappa << endl;
      cerr <<  dec_err << endl;
      return false;
   }

   // TODO: I do not know how to propagate er2 to newer errors
   // remove all negative values from vec
   for_each( vec.begin(),vec.end(), [](double& v){ v=max(0.,v); } );
   // and replace er2 with vec (er2 = errors^2)
   er2 = vec;

   return true;
}

//
//--------------------------------------------------------------------
bool decon_vec2( vector<double>& vec, vector<double>& er2,
                 double wbin, double s1, double s2, double f ) {
//--------------------------------------------------------------------
   vector<double> res = dCDFnorm2(wbin, s1,s2,f);
   if ( res.empty() ) {
      return false;
   }
   int n = vec.size();
   res.resize(n,0.);

   int niter = 100;   // number of iterations
   int nrep = 1;      // number of repetitions, for boosted dec.
   double boost = 1.; // boosting coefficient
   static TSpectrum* s = new TSpectrum();
   const char* dec_err =
   s -> Deconvolution(vec.data(), res.data(), n, niter,nrep,boost);
//    s -> DeconvolutionRL(vec.data(), res.data(), n, niter,nrep,boost);

   if ( dec_err ) {
      cerr << " ERROR in decon_vec2:  " << endl;
      cerr << " vec.size= " << vec.size() << " wbin= " << wbin
           << " s1= " << s1 << " s2= " << s2 << " f= " << f << endl;
      cerr <<  dec_err << endl;
      return false;
   }

   // TODO: I do not know how to propagate er2 to newer errors
   // remove all negative values from vec
   for_each( vec.begin(),vec.end(), [](double& v){ v=max(0.,v); } );
   // and replace er2 with vec (er2 = errors^2)
   er2 = vec;

   return true;
}

//--------------------------------------------------------------------
bool cor_sig( vector<double>& vec, vector<double>& er2, double wbin,
              int model, double s1, double s2, double f ) {
//--------------------------------------------------------------------
   if ( model == 1 ) {
      return decon_vec(vec,er2,wbin,s1);
   } else if ( model == 2 ) {
      return deconG_vec(vec,er2,wbin,s1,s2); // s2 = kappa
   } else if ( model == 3 ) {
      return decon_vec2(vec,er2,wbin,s1,s2,f);
   }
   return true;
}

// {{{1 Class for fitting
//--------------------------------------------------------------------
class myChi2 {
   public:
      // ctor:
      //--------------------------------------------------------------
      myChi2(TH1D* hst[],double EMin,double EMax) :
         Emin(EMin),Emax(EMax) {
      //--------------------------------------------------------------
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

            // MC signal: [2] - normalized, [5] - not normalized.
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
         return Ndata;
      }

      //--------------------------------------------------------------
      void SetNpol(int n) {
      //--------------------------------------------------------------
         npol = n;
      }

      //--------------------------------------------------------------
      void SetModel(int m) {
      //--------------------------------------------------------------
         model = m;
      }

      // chi^2 - function
      // the signature of this operator() MUST be exactly this:
      //--------------------------------------------------------------
      double operator() (const double* p) {
      //--------------------------------------------------------------
         static const double maxResValue = DBL_MAX / 1000;

         double Sc = p[0];   // scale constant for signal
         double Sc2 = SQ(Sc);
         // correction of the MC-signal shape
         double Gsig1 = p[1]; // sigma1 for deconvolution
         double Gsig2 = p[2]; // kappa/sigma2
         double frac =  p[3]; // fraction of the sigma2
         vector<double> sig_n(sig), er_sig_n(er_sig);
         if (!cor_sig(sig_n,er_sig_n,wbin, model, Gsig1,Gsig2,frac)) {
            return maxResValue;
         }

         double chi2 = 0;
         int nbins = en.size();
         for ( int i = 0; i < nbins; i++ ) {
            double x = en[i];
            if ( x < Emin || x > Emax ) {
               continue;
            }

            // re-weight MC background by polynomial
            int ish = model+1;
            double bgscl = ChebN(npol,Emin_chb,Emax_chb,x,&p[ish]);
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

      int model = 1;

      int npol = 3;
      const double Emin_chb = 3.0, Emax_chb = 3.2;
};
//--------------------------------------------------------------------

// {{{1 Corrections for histograms and print final numbers
//--------------------------------------------------------------------
TH1D* cor_sig_hst(const TH1D* hst, int model, const double* p) {
//--------------------------------------------------------------------
   // ignore over and underflow bins

   // create copy of hst
   string new_name = string(hst->GetName()) + string("_smear");
   TH1D* hst_n = (TH1D*)hst->Clone(new_name.c_str());
//   return hst_n;

   int n = hst->GetNbinsX();
   double wbin = hst->GetBinWidth(1);
   vector<double> sig(n),er2(n);
   for ( int i = 0; i < n; ++i ) {
      sig[i] = hst->GetBinContent(1+i);
      er2[i] = SQ(hst->GetBinError(1+i));
   }

   double Sc = p[0];   // scale constant for signal
   double Gsig1 = p[1]; // sigma1 for deconvolution
   double Gsig2 = p[2]; // sigma2
   double frac =  p[3]; // fraction of the sigma2
   if ( !cor_sig(sig,er2,wbin, model, Gsig1,Gsig2,frac) ) {
      return hst_n;
   }

   // fill hst_n
   for ( int i = 0; i < n; ++i ) {
      hst_n->SetBinContent(i+1,Sc*sig[i]);
      hst_n->SetBinError(i+1,Sc*sqrt(er2[i]));
   }
   return hst_n;
}

//--------------------------------------------------------------------
TH1D* rescale_bg( TH1D* hst, const double* p, int npol ) {
//--------------------------------------------------------------------
   // rescale ignoring over and underflow bins
   const double Emin_chb = 3.0, Emax_chb = 3.2;

   // create copy of hst
   string new_name = string(hst->GetName()) + string("_rescl");
   TH1D* hst_n = (TH1D*)hst->Clone(new_name.c_str());

   int nbins = hst->GetNbinsX();
   for ( int n = 1; n <= nbins; ++n ) {
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
   // [12] - MC_signal(after cor_sig_hst)

   double ndata   = 0;
   double er_data = 0;
   double nsig    = 0;
   double er_sig  = 0;
   double nbg     = 0;
   double er_bg   = 0;

   double Emin_hst = 0, Emax_hst = 0;
   int nbins = hst[0] -> GetNbinsX();
   double wbin = hst[0] -> GetBinWidth(1);
   for(int n = 1; n <= nbins; ++n) {
      double x = hst[0] -> GetBinCenter(n);
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
      ndata   += hst[0] -> GetBinContent(n);
      er_data += SQ( hst[0] -> GetBinError(n) );
      nsig    += hst[12] -> GetBinContent(n);
      er_sig  += SQ( hst[12] -> GetBinError(n) );
      nbg     += hst[7] -> GetBinContent(n);
      er_bg   += SQ( hst[7] -> GetBinError(n) );
   }

   double n1   = ndata-nbg;
   double er1  = sqrt(er_data+er_bg);
   er_data     = sqrt(er_data);
   er_sig      = sqrt(er_sig);
   er_bg       = sqrt(er_bg);
   double diff = 100*fabs(n1-nsig)/nsig;

   printf("\n");
   printf(" Number of events in [%.3f, %.3f]:\n",Emin_hst,Emax_hst);
   printf(" Data:            %9.0f +/- %4.0f\n", ndata,er_data);
   printf(" Sum(bg):         %9.0f +/- %4.0f\n", nbg,er_bg);
   printf(" NPsip(data-bg):  %9.0f +/- %4.0f --> diff=  %.2f%%\n",
         n1,er1,diff);
   printf(" NPsip(MCsig):    %9.0f +/- %4.0f\n", nsig,er_sig);
   printf("\n");
}

//--------------------------------------------------------------------
void print_eff(int date, int Model, TH1D* hst[],
      const ROOT::Fit::FitResult& res, double Emin, double Emax) {
//--------------------------------------------------------------------
// calculate efficiency of the Psi' -> pi+ pi- J/Psi selection
   // ATTENTION: here I assume:
   // [5] - MC(sig) without corrections and without normalization

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

   froot -> cd("PsipJpsiPhiEta");

   // initial number of pi+ pi- J/Psi
   TH1D* MCdec = (TH1D*)gROOT -> FindObject("mc_dec0");
   if ( !MCdec ) {
      cerr << " can not find mc_dec0" << endl;
      exit(0);
   }
   double Nppj = MCdec -> GetBinContent(65);
   printf("\n");
   printf(" Initial number of pi+ pi- J/Psi (dec# %i) = %.1f\n",
         int(MCdec -> GetBinCenter(65)), Nppj);

   // to study systematic vary parameters of MC signal (Model)
   // on one sigma:
   const vector<double> par = res.Parameters();
   const vector<double> er_par = res.Errors();
   const int niter = 1 + 2*Model; // first is nominal parameters
   vector<double> eff(niter,0.), err(niter,0);
   double sys_err = 0;
   for (int it = 0; it < niter; ++it ) {
      vector<double> par1(par);
      par1[0] = 1.; // absolute scaling should be 1
      if ( it != 0 ) {
         int ipar = (it+1)/2; // 1,2 -> 1; 3,4 -> 2 ...
         if ( res.HasMinosError(ipar) ) {
            if ( it%2 ) {
               par1[ipar] -= fabs(res.LowerError(ipar));
            } else {
               par1[ipar] += fabs(res.UpperError(ipar));
            }
         } else {
            par1[ipar] += (1-2*(it%2))*er_par[ipar]; // -/+ err
         }
      }

//       TH1D* hMC = hst[5]; // debug
      TH1D* hMC = cor_sig_hst(hst[5],Model,par1.data());
      double Nmc = 0.;
      double er_Nmc = 0.;

      int nbins = hMC -> GetNbinsX();
      for ( int n = 1; n <= nbins; ++n ) {
         double x = hMC -> GetBinCenter(n);
         if( x < Emin ) { continue; }
         if( x > Emax ) { break; }
         Nmc    += hMC -> GetBinContent(n);
         er_Nmc += SQ( hMC -> GetBinError(n) );
      }
      er_Nmc = sqrt(er_Nmc);

      eff[it] = Nmc / Nppj;
      err[it] = eff[it] * sqrt(SQ(er_Nmc/Nmc) + 1./Nppj);
//       printf(" it=%i eff: %.5f +/- %.5f\n",it,eff[it],err[it]);

      if ( it > 0 ) {
         double diff = fabs(eff[it]-eff[0]);
         if( diff > sys_err ) sys_err = diff;
      }
   }

   // error on normalization constant
   double sys_norm = 0;
   if ( res.HasMinosError(0) ) {
      sys_norm = max( fabs(res.LowerError(0)),
                      fabs(res.UpperError(0)) ) / par[0];

   } else {
      sys_norm = er_par[0]/par[0];
   }
   sys_err = sqrt(SQ(sys_err)+SQ(eff[0]*sys_norm));

   printf("\n");
   printf(" Efficiency in [%.3f, %.3f]:\n", Emin,Emax);
   printf(" eff(pi+pi-J/Psi): %.3f +/- %.3f +/- %.3f %%\n",
         100*eff[0],100*err[0],100*sys_err);
}

// {{{1 Fit
//--------------------------------------------------------------------
void DoFit(int date) {
//--------------------------------------------------------------------
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

   // ========================= Fit with ROOT ========================
   // == fit configuration
   ROOT::Fit::Fitter fitter;
   // set parameters of fitter: (Minuit,Minuit2,Fumili,GSLMultiFit...)
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");
   fitter.Config().MinimizerOptions().SetPrintLevel(3);
   // print the default minimizer option values
   ROOT::Math::MinimizerOptions min_opt;
   min_opt.Print();

   // type of correction and the number of correction parameters
   // in the function cor_sig()
   const unsigned int Model = (date == 2009) ? 3 : 1;
//    const unsigned int Model = 3;
   chi2_fit.SetModel(Model);

   // order of polynomial for background correction
   const unsigned int Npol = 3;
   chi2_fit.SetNpol(Npol);

   // set names and start values for parameters
   vector<string> par_name;
   if ( Model == 1 ) {
      par_name = {"#it{N}","#sigma"};
   } else if ( Model == 2 ) {
      par_name = {"#it{N}","#sigma", "#kappa"};
   } else if ( Model == 3 ) {
      par_name = {"#it{N}","#sigma_{1}","#sigma_{2}","frac"};
   } else {
      cerr << " FATAL: Model= " << Model << endl;
      return;
   }
   for ( int i = 0; i <= Npol; ++i ) {
      par_name.push_back( string("p") + to_string(i) );
   }

   vector<double> par_ini;
   if ( date == 2009 ) {
      if ( Model == 1 ) {
//          par_ini = { 0.9788,0.868e-3,
//                      1.0343,0.0069,0.0083,0.0011 };
         par_ini = { 0.9823,0.862e-3,
                     1.0343,0.0069,0.0083,0.0011 }; // fix-bg
      } else if ( Model == 2 ) {
         par_ini = { 0.9786,0.873e-3,-0.0015,
                     1.0343,0.0069,0.0083,0.0011 };
      } else if ( Model == 3 ) {
         par_ini = { 0.9756, 1.24e-3, 0.44e-3, 0.48,
                     1.0340,0.0047,-0.0092,-0.0005 }; // final
//          par_ini = { 0.9800, 1.18e-3, 0.42e-3, 0.44,
//                      1.0343,0.0069,0.0083,0.0011 }; // fix bg
      }
   } else if ( date == 2012 ) {
      if ( Model == 1 ) {
         par_ini = { 0.9824,0.627e-3,
                     1.0482,0.0015,-0.0135,-0.0030 }; // final
//          par_ini = { 0.9875,0.611e-3,
//                      1.0478,0.0044,0.0076,0.0008 }; // fix-bg
      } else if ( Model == 2 ) {
         par_ini = { 0.9874,0.610e-3,-0.0646,       // 40272.8 ???
                     1.0478,0.0044,0.0076,0.0008 }; // fix-bg
      } else if ( Model == 3 ) {
         par_ini = { 0.9874, 1.e-3, 0.2e-3, 0.5,
                     1.0478,0.0044,0.0076,0.0008 }; // == M1
      }
   }

   const unsigned int Npar = par_name.size(); // number of parameters
   // must be first
   fitter.Config().SetParamsSettings(Npar,par_ini.data());
   for(unsigned int i = 0; i < Npar; i++) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
   fitter.Config().ParSettings(0).SetStepSize(1e-3);    // Norm
   fitter.Config().ParSettings(0).SetLimits(0.1,10.);
//    fitter.Config().ParSettings(0).Fix();
   fitter.Config().ParSettings(1).SetStepSize(1e-4);    // sig1
   fitter.Config().ParSettings(1).SetLimits(0.,1.e-1);
//    fitter.Config().ParSettings(1).Fix();

   if ( Model == 2 ) {
      fitter.Config().ParSettings(2).SetStepSize(1e-4); // kappa
//       fitter.Config().ParSettings(2).Fix();
   }

   if ( Model == 3 ) {
      fitter.Config().ParSettings(2).SetStepSize(1e-4); // sig2
      fitter.Config().ParSettings(2).SetLimits(0.,1.e-1);
//       fitter.Config().ParSettings(2).Fix();
      fitter.Config().ParSettings(2).SetStepSize(1e-2); // frac [0,1]
      fitter.Config().ParSettings(3).SetLimits(0.,1.);
//       fitter.Config().ParSettings(3).Fix();
   }

   bool fixbg = false;
   if ( fixbg ) {
      for ( int i = 0; i <= Npol; ++i ) { // fix pars of background
         int ii = 1 + Model + i;
         fitter.Config().ParSettings(ii).Fix();
      }
   }

   // == Fit
   fitter.FitFCN(Npar, chi2_fit, nullptr, chi2_fit.Size(), true);

   // to obtain reliable errors:
//   fitter.CalculateHessErrors();
  fitter.CalculateMinosErrors();

   // == Fit result
   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double chi2   = res.Chi2();
   int ndf       = res.Ndf();
   const vector<double> par    = res.Parameters();
   const vector<double> er_par = res.Errors();

   // corrections for MC signal
   hst[12] = cor_sig_hst(hst[2], Model, par.data());

   // corrections for background
   int ish = 1 + Model;
   hst[13] = rescale_bg( hst[3], &par[ish], Npol);
   hst[14] = rescale_bg( hst[4], &par[ish], Npol);

   // ========================= draw results ========================
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   set_draw_opt(hst);
   SetHstFace(hst[0]);
   hst[0]->SetAxisRange(ExMin,ExMax,"X");
//    hst[0]->SetAxisRange(3.085,3.11,"X");
   hst[0]->GetXaxis()->SetTitleOffset(1.1);
   hst[0]->GetYaxis()->SetTitleOffset(1.25);

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
   hst[12]->Draw("SAME,HIST"); // signal

//    hst[14]->Draw("SAME,HIST"); // Rescaled Bg not pi+pi-J/Psi

   TLegend* leg = new TLegend(0.60,0.70,0.89,0.89);
//    leg->SetHeader("Recoil Mass of #pi^{#plus}#pi^{#minus}", "C");
   leg->AddEntry(hst[0],
         (string("Data ")+to_string(date)).c_str(), "EP");
   leg->AddEntry(hst[8],
         Form("#color[%i]{Sig + Bkg + Cont}",
         hst[8]->GetLineColor() ),"L");
   leg->AddEntry(hst[12],
         Form("#color[%i]{Sig}",
         hst[12]->GetLineColor() ),"L");
   leg->AddEntry(hst[7],
         Form("#color[%i]{Bkg + Cont}",
         hst[7]->GetLineColor() ),"L");

//    leg->AddEntry(hst[14],
//          Form("#color[%i]{MC bg non #pi^{#plus}#pi^{#minus}J/#Psi}",
//          hst[14]->GetLineColor() ),"L");

   leg->Draw();

   TPaveText* pt = new TPaveText(0.11,0.57,0.40,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
//    pt->AddText( Form("#bf{Fit in [%.2f,%.2f]}",ExMin,ExMax) );
   pt->AddText( Form("#chi^{2}/ndf =  %.0f / %i",chi2,ndf) );
   pt->AddText( Form("%s = %.4f #pm %.4f",
                par_name[0].c_str(),par[0],er_par[0]) );
   pt->AddText( Form("%s = %.3f #pm %.3f MeV ",
                par_name[1].c_str(),1e3*par[1],1e3*er_par[1]) );
   if ( Model == 2 ) {
      pt->AddText( Form("%s = %.4f #pm %.4f",
                   par_name[2].c_str(),par[2],er_par[2]) );
   }
   if ( Model == 3 ) {
      pt->AddText( Form("%s = %.3f #pm %.3f MeV ",
                   par_name[2].c_str(),1e3*par[2],1e3*er_par[2]) );
      pt->AddText( Form("%s = %.3f #pm %.3f",
                   par_name[3].c_str(),par[3],er_par[3]) );
   }
   for ( int i = 1+Model; i <= 1+Model+Npol; ++i ) {
      pt->AddText( Form("%s =  %.4f #pm %.4f",
                   par_name[i].c_str(),par[i],er_par[i]) );
   }
   pt->Draw();
   gPad->RedrawAxis();

   c1->Update();

   string pdf = string("Mrec")+to_string(date)+string("_fitM")
      +to_string(Model)+((fixbg) ? "_fixbg" : "")+string(".pdf");
   c1->Print(pdf.c_str());

   // diff data-MC_sum
   hst[17] = (TH1D*)hst[0] -> Clone("Diff");
   hst[17] -> Add(hst[8],-1.);
   hst[17] -> SetAxisRange(3.08,3.12,"X");
   hst[17] -> GetYaxis() -> SetTitle("#delta(data-MC)");
   hst[17] -> GetYaxis() -> SetMaxDigits(3);
   TCanvas* c2 = new TCanvas("c2","...",900,0,800,800);
   c2 -> cd();
   hst[17] -> Draw("E");
   gPad -> RedrawAxis();
   c2 -> Update();
   string pdf2 = string("diff")+to_string(date)+string("_m")
      +to_string(Model)+((fixbg) ? "_fixbg" : "")+string(".pdf");
   c2 -> Print(pdf2.c_str());


   // numbers for E in [Emin,Emax]
   print_Numbers(hst,3.092,3.102);// see cuts.h !
   print_Numbers(hst,3.055,3.145);// BAM-42

   print_eff(date,Model,hst,res,3.092,3.102);
   print_eff(date,Model,hst,res,3.055,3.145);
}

// {{{1 Main
//--------------------------------------------------------------------
void MrecFit() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(112);
//    gStyle->SetLegendTextSize(0.03);
   gStyle->SetLegendFont(42);
//    gStyle->SetFitFormat(".8g"); // DEBUG

   int date=2009;
//    int date=2012;

   DoFit(date);
}

