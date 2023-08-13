// MrecFit.cc - Data fitting by correcting of MC signal and scaling of
// MC background
//   -> Mrec_{YEAR}_M{MODEL}_T{Npol}_[hohc].pdf
// - calculate number of J/Psi in Psi' -> J/Psi pi+ pi- decay
//   and efficiency of selection;
// - estimate systematic associated with a fit model

#include "RewTrkPiK.hpp"    // RewTrk functions with HC

// {{{1 helper functions and constants
//--------------------------------------------------------------------
// GLOBAL: name of folder with root files
string Dir;
bool USE_NOHC_SIGNAL_MC = true; // use MC without Helix Corrections

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
   hst[2]->SetLineWidth(1);
   // MC bg from pi+pi-J/Psi
   // hst[3]->SetLineColor(kBlue+3);
   hst[3]->SetLineColor(kCyan+2);
   hst[3]->SetLineWidth(2);
   // MC bg not pi+pi-J/Psi
   hst[4]->SetLineColor(kMagenta+1);
   hst[4]->SetLineWidth(2);
}

// {{{1 Fill histograms
//--------------------------------------------------------------------
vector<TH1D*> fill_mrec(string fname, string hname,
      int date, bool isMC = false, double shift = 0.) {
//--------------------------------------------------------------------
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
         ";Entries/0.1 MeV/c^{2}",
         2000,3.0,3.2 // 1bin = 0.1MeV
         );
   hst->Sumw2(true);

   Long64_t nentries = nt1->GetEntries();

   if ( !isMC ) { // data: return all recoil masses
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
            hist[0]->Fill(Mrs+shift,wp*wm);
         }
         for (const auto& mr : *Mrec ) {
            hist[1]->Fill(mr+shift);
         }
      } else {
         for (const auto& mr : *Mrec ) {
            hist[2]->Fill(mr+shift);
         }
      }
   }
   return hist;
}

//--------------------------------------------------------------------
vector<TH1D*> fill_hist(int date) {
//--------------------------------------------------------------------
#include "norm.h"
   double shift = 0.0;
   if ( date == 2009 ) {
      // shift = 0.000225; // chi2(M2) = 3106?
      shift = 0.000226; // chi2(M2) = 3144 ? 3057 ***
      // shift = 0.000227; // chi2(M2) = 3077?
   } else if ( date == 2012 ) {
      // shift = 0.000385; // chi2(M2) = 6183?
      shift = 0.000386; // chi2(M2) = 6185 ? 6013 ***
      // shift = 0.000387; // chi2(M2) = 6168?
   } else if ( date == 2021 ) {
      // shift = 0.000606; // chi2(M2) = 23 481
      shift = 0.000607; // chi2(M2) = 22 947 ?22 637***
      // shift = 0.000608; // chi2(M2) = 23 122
   }

   string sd(Form("%02i",date%100));
   vector<string> hname = {
      "data"+sd, "data3650_"+sd,
      "mc"+sd+"sig", "mc"+sd+"bg1", "mc"+sd+"bg2",
      "mc"+sd+"sig10"
   };

   vector<TH1D*> hst(hname.size(),nullptr);

   // read from cache file
   string cachef = Dir + (USE_NOHC_SIGNAL_MC ? "NoHC/" : "")
      + string(Form("MrecFit_%i",date));
   if ( fabs(shift) > 1e-7 ) {
      cachef += ((shift<0) ? "_m" : "_")
         + to_string(int(fabs(shift*1e6)));
   }
   cachef += ".root";
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

   string datafile( Form("data_%02ipsip_all.root",date%100) );
   string mcincfile( Form("mcinc_%02ipsip_all.root",date%100) );
   if ( USE_NOHC_SIGNAL_MC ) {
      mcincfile = "NoHC/" +  mcincfile;
   }

   auto hs = fill_mrec(datafile, hname[0], date);
   hst[0] = hs[0];

   hs = fill_mrec("data_3650_2021.root", hname[1], 2021);
   hst[1] = hs[0];
   hst[1]->Scale(C_Dat.at(date));

   hs = fill_mrec(mcincfile, "mc"+sd, date, true, shift);
   hst[2]=hs[0];
   hst[5] = (TH1D*)hst[2]->Clone( hname[5].c_str() );
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

// {{{1 Functions for fitting
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

   // printf(" DEBUG: gw.size= %zu, wsum-gw[0]= %.9f (must be 1)\n",
            // gw.size(), wsum-gw[0]);

   int ngw = gw.size();
   vector<double> res(2*ngw-1); // result
   int m = ngw-1; // midle point of res
   res[m] = gw[0];
   for ( int i = 1; i < ngw ; ++i ) {
      res[m-i] = res[m+i] = gw[i];
   }

   // test: shift on one bin to right -> the same result!
   // res.resize(2*ngw);
   // for ( int i = 2*ngw-1; i > 0 ; --i ) {
      // res[i] = res[i-1];
   // }
   // res[0] = 0;

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

   // printf(" DEBUG: res.size=%zu, wsum=%.9f(must be 1), res[0]=%g\n",
            // res.size(), wsum, res[0]);

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

   // printf(" DEBUG: gw.size= %zu, wsum-gw[0]= %.9f (must be 1)\n",
            // gw.size(), wsum-gw[0]);

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
      double wbin, double s1, double s2, double f, int NG ) {
//--------------------------------------------------------------------
   vector<double> res = ( NG == 2 ) ?
      dCDFnorm2(wbin, s1,s2,f) //two gausses
      :
      dCDFnorm(wbin, s1); // one gauss

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
      // int jmin = max(0,mres-i);
      // int jmax = min(nres,n+mres-i);
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
      s->Deconvolution(vec.data(), res.data(), n, niter,nrep,boost);
      // s->DeconvolutionRL(vec.data(), res.data(), n, niter,nrep,boost);

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
      s->Deconvolution(vec.data(), res.data(), n, niter,nrep,boost);
      // s->DeconvolutionRL(vec.data(), res.data(), n, niter,nrep,boost);

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
      s->Deconvolution(vec.data(), res.data(), n, niter,nrep,boost);
      // s->DeconvolutionRL(vec.data(), res.data(), n, niter,nrep,boost);

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

// {{{1 Class for fitting
//--------------------------------------------------------------------
class myChi2 {
   public:
      myChi2(const vector<TH1D*>& hst, double EMin, double EMax);

      // interface to convolution/deconvolution functions
      bool cor_sig( vector<double>& vec, vector<double>& er2,
            double s1, double s2, double f ) const;

      // operator() => chi^2 function
      double operator() (const double* p) { return Chi2(p); }
      double Chi2(const double* p);

      size_t Size() const {
         return Ndata;
      }

      void SetModel(int m) {
         model = m;
         switch ( model ) {
            case 1: m_par = 2; break;
            case 2: m_par = 4; break;
            default:
               cerr << "FATAL: SetModel: " << model << endl;
               exit(EXIT_FAILURE);
         }
      }

      int GetMPar() const {
         return m_par;
      }

      void SetNpol(int n) {
         npol = n;
      }

   private:
      double Emin=3.01, Emax=3.19; // limits to fit signal
      size_t Ndata=0;

      double wbin=0.;    // bin width
      vector<double> en; // energy
      vector<double> data, er_data;
      vector<double> cont, er_cont;
      vector<double> sig, er_sig;
      vector<double> bg, er_bg;
      vector<double> bgn, er_bgn;

      // models
      int model = 1; // see cor_sig() function
      int m_par = 2; // number parameters of the model

      int npol = 3;
      const double Emin_chb = 3.0, Emax_chb = 3.2;
};

// {{{2 ctor:
//--------------------------------------------------------------------
myChi2::myChi2(const vector<TH1D*>& hst, double EMin, double EMax) {
//--------------------------------------------------------------------
// Emin/Emax -> region to fit signal
   Emin = EMin;
   Emax = EMax;

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

// {{{2 cor_sig(): wbin & model from class
//--------------------------------------------------------------------
bool myChi2::cor_sig( vector<double>& vec, vector<double>& er2,
      double s1, double s2, double f ) const {
//--------------------------------------------------------------------
   switch ( model ) {
      case 1:
         return conv_vec(vec,er2,wbin,s1,s2,f,1); // gauss
      case 2:
         return conv_vec(vec,er2,wbin,s1,s2,f,2); // gauss(+)gauss
      // case 1:
         // return decon_vec(vec,er2,wbin,s1);
      // case 2:
         // return deconG_vec(vec,er2,wbin,s1,s2); // s2 = kappa
      // case 3:
         // return decon_vec2(vec,er2,wbin,s1,s2,f);
   }
   return true;
}

// {{{2 Chi2() function
//--------------------------------------------------------------------
double myChi2::Chi2(const double* p) {
//--------------------------------------------------------------------
   static const double maxResValue = DBL_MAX / 1000;

   double Sc = p[0];   // scale constant for signal
   double Sc2 = SQ(Sc);
   // correction of the MC-signal shape
   double Gsig1 = p[1]; // sigma1 for deconvolution/convolution
   double Gsig2 = p[2]; // sigma2/(kappa)
   double frac =  p[3]; // fraction of the sigma2
   vector<double> sig_n(sig), er_sig_n(er_sig);
   if ( !cor_sig(sig_n,er_sig_n, Gsig1,Gsig2,frac) ) {
      cerr << "ERROR:: cor_sig failed: model=" << model
         << " Gsig1=" << Gsig1 << " Gsig2=" << Gsig2
         << " frac=" << frac << endl;
      return maxResValue*Ndata;
   }

   double chi2 = 0;
   int nbins = en.size();
   for ( int i = 0; i < nbins; i++ ) {
      double x = en[i];
      if ( x < Emin || x > Emax ) {
         continue;
      }

      // re-weight MC background by polynomial
      double bgscl = ChebN(npol,Emin_chb,Emax_chb,x,&p[m_par]);
      double bgscl2 = bgscl*bgscl;

      double dif = data[i] -
         ( Sc*sig_n[i] +
           cont[i] +
           bgscl * (bg[i] + bgn[i]) );
      double er2 = er_data[i] +
         Sc2*er_sig_n[i] +
         er_cont[i] +
         bgscl2 * (er_bg[i] + er_bgn[i]);

      double resval = dif*dif / er2;
      // avoid inifinity or nan in chi2
      if ( resval < maxResValue ) {
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
TH1D* cor_sig_hst(const TH1D* hst, const myChi2& ch2,
      const vector<double>& p) {
//--------------------------------------------------------------------
   // create copy of hst
   string new_name = string(hst->GetName()) + string("_smear");
   TH1D* hst_n = (TH1D*)hst->Clone(new_name.c_str());
   // return hst_n;

   int n = hst->GetNbinsX();
   vector<double> sig(n),er2(n);
   for ( int i = 0; i < n; ++i ) {
      sig[i] = hst->GetBinContent(1+i);
      er2[i] = SQ(hst->GetBinError(1+i));
   }

   double Sc = p[0];    // scale constant for signal
   double Gsig1 = p[1]; // sigma1
   double Gsig2 = p[2]; // sigma2
   double frac =  p[3]; // fraction of the sigma2
   if ( !ch2.cor_sig(sig,er2, Gsig1,Gsig2,frac) ) {
      return hst_n;
   }

   // fill hst_n, ignore over and underflow bins
   for ( int i = 0; i < n; ++i ) {
      hst_n->SetBinContent(i+1,Sc*sig[i]);
      hst_n->SetBinError(i+1,Sc*sqrt(er2[i]));
   }
   return hst_n;
}

//--------------------------------------------------------------------
TH1D* rescale_bg(const TH1D* hst, const vector<double>& par) {
//--------------------------------------------------------------------
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
void print_Numbers(const vector<TH1D*>& hst, const TH1D* MCsig_cor,
      const TH1D* SumBG, double Emin, double Emax) {
//--------------------------------------------------------------------
   // ATTENTION: here I assume:
   // [0] - data
   // sumBG - CONT(norm to Lum_data) + BG(after rescale_bg)
   // MCsig_cor - MC_signal(after cor_sig_hst)

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
      nsig    += MCsig_cor -> GetBinContent(n);
      er_sig  += SQ( MCsig_cor -> GetBinError(n) );
      nbg     += SumBG -> GetBinContent(n);
      er_bg   += SQ( SumBG -> GetBinError(n) );
   }

   double n1   = ndata-nbg;
   double er1  = sqrt(er_data+er_bg);
   er_data     = sqrt(er_data);
   er_sig      = sqrt(er_sig);
   er_bg       = sqrt(er_bg);
   double diff = 100*fabs(n1-nsig)/nsig;

   bool short_print = true;
   if ( short_print ) {
      printf(" [%.3f, %.3f]: ", Emin_hst,Emax_hst);
      printf(" NJpsi(MCsig):    %9.0f +/- %4.0f\n", nsig,er_sig);
      return;
   }

   printf("\n");
   printf(" Number of events in [%.3f, %.3f]:\n",Emin_hst,Emax_hst);
   printf(" Data:            %9.0f +/- %4.0f\n", ndata,er_data);
   printf(" Sum(bg):         %9.0f +/- %4.0f\n", nbg,er_bg);
   printf(" NJpsi(data-bg):  %9.0f +/- %4.0f --> diff=  %.2f%%\n",
         n1,er1,diff);
   printf(" NJpsi(MCsig):    %9.0f +/- %4.0f\n", nsig,er_sig);
   printf("\n");
}

// {{{1 print efficiency and error on it
//--------------------------------------------------------------------
void print_eff(int date, const TH1D* MCsig, const myChi2& ch2,
      const vector<double>& par, const vector<double> er_par,
      double Emin, double Emax) {
//--------------------------------------------------------------------
// calculate efficiency of the Psi' -> pi+ pi- J/Psi selection
   // ATTENTION: here I assume:
   // MCsig - MC(sig) without corrections and without normalization

   string mcincfile( Form("mcinc_%02ipsip_all.root",date%100) );

   string fname = Dir + (USE_NOHC_SIGNAL_MC ? "NoHC/" : "")
      + mcincfile;
   // cout << endl << " Eff file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "ERROR in "<< __func__
         << ": can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }
   froot->cd("PsipJpsiPhiEta");

   // initial number of pi+ pi- J/Psi
   TH1D* MCdec = (TH1D*)gROOT -> FindObject("mc_dec0");
   if ( !MCdec ) {
      cerr << "ERROR in "<< __func__
         << " can not find mc_dec0" << endl;
      exit(EXIT_FAILURE);
   }
   double Nppj = MCdec->GetBinContent(65);

   // to study systematic vary 'model'-parameters of MC signal
   // on one sigma:
   const int m_par = ch2.GetMPar();
   const int niter = 2*m_par-1;
   vector<double> eff(niter,0.), err(niter,0);
   double sys_err = 0;
   for ( int it = 0; it < niter; ++it ) {
      vector<double> par1(par);
      par1[0] = 1.; // absolute scaling should be 1
      int ipar = it;
      if ( it != 0 ) {
         ipar = (it+1)/2; // 1,2 -> 1; 3,4 -> 2 ...
         par1[ipar] += (1-2*(it%2))*er_par[ipar]; // -/+ err
      }
      TH1D* MCsig_cor = cor_sig_hst(MCsig, ch2, par1);

      double Nmc = 0.;
      double er_Nmc = 0.;
      int nbins = MCsig_cor->GetNbinsX();
      for ( int n = 1; n <= nbins; ++n ) {
         double x = MCsig_cor->GetBinCenter(n);
         if( x < Emin ) { continue; }
         if( x > Emax ) { break; }
         Nmc    += MCsig_cor->GetBinContent(n);
         er_Nmc += SQ( MCsig_cor->GetBinError(n) );
      }
      er_Nmc = sqrt(er_Nmc);

      eff[it] = Nmc / Nppj;
      err[it] = eff[it] * sqrt(SQ(er_Nmc/Nmc) + 1./Nppj);
      // printf(" ipar=%i par=%g par1=%g -> eff: %.5f +/- %.5f\n",
            // ipar, par[ipar], par1[ipar], eff[it], err[it]);

      if ( it > 0 ) {
         double diff = fabs(eff[it]-eff[0]);
         if( diff > sys_err ) sys_err = diff;
      }
   }

   double sum_err = sqrt(SQ(err[0])+SQ(sys_err));
   double sum_rel_err = sum_err/eff[0];

   printf("\n");
   printf(" Initial number of pi+ pi- J/Psi (dec# %i) = %.1f\n",
         int(MCdec->GetBinCenter(65)), Nppj);
   printf(" Efficiency in [%.3f, %.3f]:\n", Emin,Emax);
   printf(" eff(pi+pi-J/Psi): %.3f +/- %.3f +/- %.3f %%",
         100*eff[0],100*err[0],100*sys_err);
   printf("  (sum_err= %.3f or %.2f%%)\n",
         100*sum_err,100*sum_rel_err);
   printf("\n");
}

// {{{1 diff data-MC_sum
//--------------------------------------------------------------------
void plot_diff(const TH1D* Data, const TH1D* SumMC,
      TPaveText* pt, string pdf) {
//--------------------------------------------------------------------
   TH1D* Diff = (TH1D*)Data->Clone("Diff");
   Diff->Add(SumMC,-1.);
   Diff->SetAxisRange(3.08,3.12,"X");
   Diff->GetYaxis()->SetTitle("#delta(data-MC)");
   Diff->GetYaxis()->SetMaxDigits(3);
   TCanvas* c2 = new TCanvas("diff","...",900,0,800,800);
   c2 -> cd();
   Diff->Draw("E");
   pt->Draw();
   gPad->RedrawAxis();

   c2 -> Update();
   c2 -> Print(pdf.c_str());
}

// {{{1 Fit
//--------------------------------------------------------------------
void DoFit(int date) {
//--------------------------------------------------------------------
   vector<TH1D*> hst = fill_hist(date);

   // limits for drawing, bin size ten times smaller than FitSB
   double maxWin = 0, minWin = 0;
   if ( date==2009 ) {
      maxWin = 7.e5;
      minWin = 0.99e4;
   } else if(date==2012) {
      maxWin = 2.e6;
      minWin = 3e4;
   } else if(date==2021) {
      maxWin = 1.3e7;
      minWin = 2e5;
   }
   hst[0]->SetMaximum(maxWin);
   hst[0]->SetMinimum(minWin);

   // fit signal in region (EMin; EMax)
   double EMin = 3.01;
   double EMax = 3.19;
   myChi2 chi2_fit(hst,EMin,EMax);

   // order of polynomial for background, the same as in MrecFitSB
   const size_t Npol = (date==2009) ? 4 : 5;
   chi2_fit.SetNpol(Npol);

   // set type of corrections, see myChi2::cor_sig()
   const int Model = 2;
   chi2_fit.SetModel(Model);
   const int m_par = chi2_fit.GetMPar();

   // ========================= Fit with ROOT ========================
   // == fit configuration
   ROOT::Fit::Fitter fitter;
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");
   fitter.Config().MinimizerOptions().SetPrintLevel(3);
   // print the default minimizer option values
   ROOT::Math::MinimizerOptions min_opt;
   min_opt.Print();

   // set names and start values for parameters
   vector<string> par_name;
   switch ( Model ) {
      case 1:
         par_name = {"#it{N}","#sigma"};
         break;
      case 2:
         par_name = {"#it{N}","#sigma_{1}","#sigma_{2}","frac"};
         break;
      // case 3:
         // par_name = {"#it{N}","#sigma", "#kappa"};
   }
   for( size_t i = 0; i <= Npol; ++i ) {
      par_name.push_back( string( Form("p%zu",i) ) );
   }
   const size_t Npar = par_name.size(); // number of parameters

   vector<double> par_ini;
   if ( date == 2009 ) {
      if ( Model == 1 ) {
         par_ini = { 0.969, 0.86e-3,
            1.094,0.005,-0.027,-0.001 }; // old +0.23MeV old
      } else if ( Model == 2 ) {
         par_ini = { 0.984, 3.20e-3, 0.62e-3, 0.87,
            1.098,0.008,0.007,0.000,0.013 }; // n3 +0.226MeV
      }
   } else if ( date == 2012 ) {
      if ( Model == 1 ) {
         par_ini = { 0.982, 0.98e-3,
            1.118,0.008,-0.030,-0.006 }; // old +0.38MeV old
      } else if ( Model == 2 ) {
         par_ini = { 0.992, 3.40e-3, 0.70e-3, 0.86,
            1.116,0.026,0.008,0.008,0.014,0.008 }; // n3 +0.386MeV
      }
   } else if ( date == 2021 ) {
      if ( Model == 1 ) {
         par_ini = { 0.945, 0.89-3,
            1.066,-0.006,-0.020,-0.002 }; // old +0.60MeV old
      } else if ( Model == 2 ) {
         par_ini = { 0.994, 3.05e-3, 0.62e-3, 0.86,
            1.078,0.006,0.009,0.006,0.010,0.005 }; // n3 +0.607MeV
      }
   }
   if ( par_ini.size() != Npar ) {
      par_ini.resize(Npar,0.);
      cout << " WARNING: set size of par_ini to " << Npar << endl;
   }

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
      fitter.Config().ParSettings(2).SetStepSize(1e-4); // sig2
      fitter.Config().ParSettings(2).SetLimits(0.,1.e-1);
//       fitter.Config().ParSettings(2).Fix();
      fitter.Config().ParSettings(2).SetStepSize(1e-2); // frac [0,1]
      fitter.Config().ParSettings(3).SetLimits(0.,1.);
//       fitter.Config().ParSettings(3).Fix();
   }

   bool fixbg = false;
   if ( fixbg ) {
      for ( size_t i = 0; i <= Npol; ++i ) { // fix pars of background
         int ii = m_par + i;
         fitter.Config().ParSettings(ii).Fix();
      }
   }

   // == Fit
   fitter.FitFCN(Npar, chi2_fit, nullptr, chi2_fit.Size(), true);

   // to obtain reliable errors:
   fitter.CalculateHessErrors();
   // fitter.CalculateMinosErrors();

   // == Fit result
   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
   double chi2   = res.Chi2();
   int ndf       = res.Ndf();
   const vector<double> par    = res.Parameters();
   const vector<double> er_par = res.Errors();

   // copy bg-params to a new vector for convenience
   vector<double> par_bg(begin(par)+m_par, end(par));
   vector<double> err_bg(begin(er_par)+m_par, end(er_par));

   // ========================= draw results ========================
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   set_draw_opt(hst);
   SetHstFace(hst[0]);
   hst[0]->SetAxisRange(EMin,EMax,"X");
   hst[0]->GetXaxis()->SetTitleOffset(1.1);
   hst[0]->GetYaxis()->SetTitleOffset(1.25);

   hst[0]->Draw("E"); // data

   TH1D* SumBG = (TH1D*)hst[1]->Clone("SumBG"); // clone Continuum
   // corrections for background
   TH1D* hstBG1 = rescale_bg( hst[3], par_bg);
   TH1D* hstBG2 = rescale_bg( hst[4], par_bg);
   SumBG->Add(hstBG1);
   SumBG->Add(hstBG2);
   SumBG->Draw("SAME,HIST");

   // corrections for MC signal
   TH1D* MCsig_cor = cor_sig_hst(hst[2], chi2_fit, par);
   MCsig_cor->SetLineStyle(kDashed);
   MCsig_cor->Draw("SAME,HIST");

   TH1D* SumMC =(TH1D*)MCsig_cor->Clone("SumMC");
   SumMC->Add(SumBG);
   SumMC->SetLineStyle(kSolid);
   SumMC->SetLineColor(kRed+1);
   SumMC->Draw("SAME,HIST");

   // hstBG2->Draw("SAME,HIST"); // Rescaled Bg2

   TLegend* leg = new TLegend(0.60,0.70,0.89,0.89);
   // leg->SetHeader("Recoil Mass of #pi^{#plus}#pi^{#minus}", "C");
   leg->AddEntry(hst[0], Form("Data %i",date), "EP");
   leg->AddEntry(SumMC,
         Form("#color[%i]{Sig + Bkg + Cont}",
            SumMC->GetLineColor() ),"L");
   leg->AddEntry(MCsig_cor,
         Form("#color[%i]{Signal}",
            MCsig_cor->GetLineColor() ),"L");
   leg->AddEntry(SumBG,
         Form("#color[%i]{Bkg + Cont}",
            SumBG->GetLineColor() ),"L");

   // leg->AddEntry(hstBG2,
         // Form("#color[%i]{MC bg non #pi^{#plus}#pi^{#minus}J/#Psi}",
         // hstBG2->GetLineColor() ),"L");

   leg->Draw();

   double Ypt = 0.63 - 0.03*(m_par-4) - 0.03*(Npol-3);
   TPaveText* pt = new TPaveText(0.11,0.57,0.40,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   // pt->AddText( Form("#bf{Fit in [%.2f,%.2f]}",ExMin,ExMax) );
   pt->AddText( Form("#chi^{2} / ndf = %.0f / %i",chi2,ndf) );
   pt->AddText( Form("%s = %.4f #pm %.4f",
            par_name[0].c_str(),par[0],er_par[0]) );
   pt->AddText( Form("%s = %.3f #pm %.3f MeV ",
            par_name[1].c_str(),1e3*par[1],1e3*er_par[1]) );
   if ( Model == 2 ) {
      pt->AddText( Form("%s = %.3f #pm %.3f MeV ",
               par_name[2].c_str(),1e3*par[2],1e3*er_par[2]) );
      pt->AddText( Form("%s = %.3f #pm %.3f",
               par_name[3].c_str(),par[3],er_par[3]) );
   }
   for ( size_t i = 0; i <= Npol; ++i ) {
      int ii = m_par + i;
      pt->AddText( Form("%s =   %.4f #pm %.4f",
               par_name[ii].c_str(),par_bg[i],err_bg[i]) );
   }
   pt->Draw();
   gPad->RedrawAxis();

   c1->Update();
   string pdf( Form("Mrec%d_fit_M%d_T%zu",date,Model,Npol) );
   if ( fixbg ) {
      pdf += "_fixbg";
   }
   if ( USE_NOHC_SIGNAL_MC ) {
      pdf += "_nohc";
   }
   pdf += ".pdf";
   c1->Print(pdf.c_str());

   // plot_diff(hst[0],SumMC,pt,string("Diff_")+pdf);

   string sepline(70,'='); // separation line
   printf("%s\n",sepline.c_str());
   // numbers for E in [Emin,Emax]
   print_Numbers(hst,MCsig_cor,SumBG,3.092,3.102);// see cuts.h !
   print_Numbers(hst,MCsig_cor,SumBG,3.055,3.145);// BAM-42
   printf("%s\n",sepline.c_str());

   print_eff(date, hst[5], chi2_fit, par, er_par, 3.092,3.102);
   print_eff(date, hst[5], chi2_fit, par, er_par, 3.055,3.145);
   printf("%s\n",sepline.c_str());
}

// {{{1 Main
//--------------------------------------------------------------------
void MrecFit() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(112);
   gStyle->SetLegendFont(42);
   // gStyle->SetFitFormat(".8g"); // DEBUG

   //========================================================
   // set the name of the folder with the root files
   // Dir = "prod-12/";  // must be the same as prod-11
   Dir = "prod_v709n3/";
   USE_NOHC_SIGNAL_MC = true; // use MC without Helix Corrections
   //========================================================

   // int date=2009;
   // int date=2012;
   int date=2021;

   DoFit(date);
}
