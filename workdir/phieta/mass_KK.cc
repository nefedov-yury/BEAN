// plot M(K+K-) distributions
// cuts: chi^2(4C) + Mgg
// -> mass_kk.pdf

#include "masses.h"

// {{{1 helper function
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

// {{{1 Breigt Wigner for phi -> KK
//----------------------------------------------------------------------
// Breigt Wigner for phi -> KK
//----------------------------------------------------------------------

//----------------------------------------------------------------------
double BreitWigner(double m, double mphi, double gphi) {
//----------------------------------------------------------------------
   constexpr double Mjpsi = 3.096916; // 3096.916 +/- 0.011 MeV
   constexpr double Meta  = 0.547862; // 547.862  +/- 0.017 MeV
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV

   constexpr double Mjpsi2 = SQ(Mjpsi);
   constexpr double Meta2  = SQ(Meta);
   constexpr double Mk2    = SQ(Mk);

   constexpr double dL = 2*Mk; // the left cutoff = 0.987354
   if ( m < dL ) { // phase space becomes zero (see p_k)
      return 0;
   }

   double m2    = SQ(m);
   double mphi2 = SQ(mphi);
   double gphi2 = SQ(gphi);

   double gam = mphi*sqrt(mphi2+gphi2);
   double kap = 2*M_SQRT2*mphi*gphi*gam / (M_PI*sqrt(mphi2+gam));

//?    double p_phi0 = sqrt(SQ(Mjpsi2-mphi2-Meta2)-4*mphi2*Meta2) / (2*Mjpsi);
//?    double p_phi  = sqrt(SQ(Mjpsi2-m2-Meta2)-4*m2*Meta2) / (2*Mjpsi);
//?    double r_phi  = p_phi / p_phi0;

   double p_k0 = 0.5*sqrt(mphi2-4*Mk2);
   double p_k  = 0.5*sqrt(m2-4*Mk2);
   double r_k  = p_k / p_k0;

   // B(m)/B(m0)
   constexpr double R  = 3.; // GeV^{-1} for Blatt-Weisskopf ff
//?    double BB_phi = sqrt( (1+SQ(R*p_phi0)) / (1+SQ(R*p_phi)) );
   double BB_k2  = (1+SQ(R*p_k0)) / (1+SQ(R*p_k));
   double BB_k = sqrt( BB_k2 );

//?    double fm = r_phi*r_k*BB_phi*BB_k;
   double fm = r_k*BB_k;

//    double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k; // == m*G(m)
   double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k2; // == m*G(m)

   return (kap * SQ(fm)) / ( SQ(m2-mphi2) + SQ(GM) );
}

//----------------------------------------------------------------------
double BreitWignerGauss( double m,
                         double mphi, double gphi, // BW parameters
                         double sigma              // Gauss
                       ) {
//----------------------------------------------------------------------
// Function Fun(x) folded with a normal distribution:
//
//        1                        [                  (X'-X)^2    ]
// ---------------- *   Integral   [ Fun(X') * exp( - --------- ) ] dX'
// sqrt(2*pi)*sigma   (-oo<X'<+oo) [                  2*sigma^2   ]
//
//----------------------------------------------------------------------
// sigma ~ 1/3 gphi => +/-5 sigma integration
//----------------------------------------------------------------------
   constexpr double one_over_sqrt2pi = 0.5*M_2_SQRTPI*M_SQRT1_2;

   // integrand lambda function
   double pp[] = { m, mphi, gphi, sigma };
   auto Lint = [](double t, void* pp) -> double{ // t == (X'-X)/sigma
      double* p = static_cast<double*>(pp);
      return exp(-0.5*t*t) * BreitWigner( (p[0]+t* p[3]),p[1],p[2]);
   };

   // desired errors:                abs    rel
   ROOT::Math::GSLIntegrator gsl_int(1.e-8, 1.e-6, 1000);
   double result = gsl_int.Integral(Lint,pp,-5.,+5.);

   // int status = gsl_int.Status();
//    if ( DEBUG ) {
      // the estimate of the absolute Error
//       if ( gsl_int.Error() > 1e-6 ) {
//          htt[0]->Fill(log10(gsl_int.Error()));
//       } else {
//          htt[0]->Fill(-7.);
//       }
//    }

   return one_over_sqrt2pi * result;
}

//----------------------------------------------------------------------
double BreitWignerGaussN( double m,
                          double mphi, double gphi, double sigma
                        ) {
//----------------------------------------------------------------------
// Numeric normalisation on one for range [dL,dU]
   constexpr double Mk    = 0.493677; // 493.677  +/- 0.016 MeV
   constexpr double dL = 2*Mk; // the left cutoff = 0.987354
   constexpr double dU = 1.08; // upper limit

   double norm = 0;
   // cash parameters
   static double cacheN = 0;
   static double cacheM = 0;
   static double cacheG = 0;
   static double cacheS = 0;
   if ( cacheN > 0 &&
         mphi == cacheM && gphi == cacheG && sigma == cacheS ) {
      norm = cacheN;
   } else {
      // integrand lambda function
      double p[] = {mphi,gphi,sigma};
      auto Lint = [](double x, void* pp) -> double{
         double* p = static_cast<double*>(pp);
         return BreitWignerGauss(x,p[0],p[1],p[2]);
      };

      // desired errors:                abs    rel
      ROOT::Math::GSLIntegrator gsl_int(1.e-7, 1.e-6, 1000);
      norm = gsl_int.Integral(Lint,p,dL,dU);
//       if ( DEBUG ) {
//          htt[1]->Fill(norm);
//          if ( gsl_int.Error() > 1e-6 ) {
//             htt[2]->Fill(log10(gsl_int.Error()));
//          } else {
//             htt[2]->Fill(-7.);
//          }
//       }
      cacheN = norm;
      cacheM = mphi;
      cacheG = gphi;
      cacheS = sigma;
   }

   return BreitWignerGauss(m,mphi,gphi,sigma) / norm;
}

// {{{1 test Breigt Wigner
//-------------------------------------------------------------------------
void test_BreitWigner() {
//-------------------------------------------------------------------------
   constexpr double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
   constexpr double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV

   constexpr double dL = 0.98;
   constexpr double dU = 1.08;

   // 1) BreitWigner (internally normalized to 1)
//    cout << " BW(1.02)= " << BreitWigner(1.02,Mphi,Gphi) << endl;
   auto Lbw = [](const double* x,const double* p) -> double {
      return p[0]*BreitWigner(x[0],p[1],p[2]);
   };
   TF1* bw = new TF1("bw", Lbw, dL, dU, 3);
   bw -> SetParNames("Norm","Mphi","Gphi");
   bw -> SetParameters(1., Mphi, Gphi);
   bw -> SetLineWidth(1);
   bw -> SetLineColor(kBlue);
   bw -> SetNpx(500);

   double norm = 1./bw -> Integral(dL,dU,1e-8);
   printf("norm = %.7f\n",norm);
   bw -> SetParameter(0, norm );

   // 2) BreitWignerGaussN
   auto Lbwgn = [](const double* x,const double* p) -> double {
      return p[0]*BreitWignerGaussN(x[0],p[1],p[2],p[3]);
   };
   TF1* bwgn = new TF1("bwgn", Lbwgn, dL, dU, 4);
   bwgn -> SetParNames("Norm","Mphi","Gphi","Sigma");
   bwgn -> SetParameters(1., Mphi, Gphi, 1.2e-3);
   bwgn -> SetLineWidth(2);
   bwgn -> SetLineColor(kRed);
   bwgn -> SetLineStyle(kDashed);
   bwgn -> SetNpx(500);

   double normN = 1./bwgn->Integral(dL,dU,1e-8);
   printf("normN = %.6f\n",normN);

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   bw -> Draw();
   bwgn -> Draw("SAME");
}

// {{{1  Data & Fit
//--------------------------------------------------------------------
 vector<double> get_mkk_hist( string fname, string hname,
       TH1D* hst[], int type=0 ) {
//--------------------------------------------------------------------
#include "cuts.h"

   constexpr double dL = 0.98;
   constexpr double dU = 1.08;
   constexpr int Nbins = 50;
   constexpr double bW = (dU-dL)/Nbins; // bin width

   // name of folder with root files
   static string dir("Ntpls/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot -> cd("SelectKKgg");
   TTree* a4c = (TTree*)gDirectory -> Get("a4c");

   TCut c_here = c_chi2+c_phi;
   if ( type == 0 ) {        // central part
      c_here += c_cpgg;
   } else if ( type == 1 ) { // side-band
      c_here += c_sbgg;
   }
   int n = a4c -> Draw("Mkk", c_here, "goff");
//    cout << " n= " << n << endl;
   double* buffer = a4c -> GetVal(0);
   vector<double> mkk(buffer, buffer+n);

   string title(";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2}");
   int ibW = int(bW*1e5);
   title += string(";Entries / ") +
      ( (ibW%10 == 0) ? string(Form("%.0f MeV/c^{2}",bW*1e3))
                      : string(Form("%.2f MeV/c^{2}",bW*1e3)) );
   hst[0] = new TH1D(hname.c_str(),title.c_str(),Nbins,dL,dU);

   hst[0] -> Sumw2(true);
   for ( const auto& mk : mkk ) {
      hst[0] -> Fill(mk);
   }
   return mkk;
}

//--------------------------------------------------------------------
TH1D* get_mkk_hist( string fname, string hname, int type=0 ) {
//--------------------------------------------------------------------
   TH1D* hist[1];
   get_mkk_hist( fname, hname, hist, type );
   return hist[0];
}

//--------------------------------------------------------------------
void plot_mass_phi() {
//--------------------------------------------------------------------
#include "cuts.h"

   // get histos
   vector<string> fnames = {
        "ntpl_3080_rs.root",
        "ntpl_3097.root",
        "ntpl_J4.root",
        "ntpl_mcgpj_3080_rs.root",
        "ntpl_mcgpj_3097.root",
        "ntpl_mcgpj_J4.root",
   };
        // "ntpl_3080_2019.root",
        // "ntpl_mcgpj_3080_2019.root",
        // "ntpl_jpsi_incl.root",

   vector<string> titles = {
      "3080MeV R-scan",
      "3097MeV (2012)",
      "3097MeV (2018)",
      "3080MeV MC #phi#eta",
      "3097MeV MC #phi#eta",
      "3097MeV MC #phi#eta 2018",
   };
      // "3080MeV (2019)",
      // "3080MeV MC(2019) #phi#eta",
      // "J/#Psi inclusive MC",
      // "3080MeV MC #phi#eta",
      // "3097MeV MC #phi#eta",
      // "3097MeV MC #phi#eta 2018",
      // "3097MeV J/#Psi-scan",

   int Nhst = 2*fnames.size(); // cp & sb
   vector<TH1D*> mkk( Nhst, nullptr );
   for ( int i = 0; i < Nhst; ++i ) {
      string hn = string("mkk") + to_string(i);
      mkk[i] = get_mkk_hist( fnames[i/2], hn, i%2 );

      if ( i%2 == 0 ) { // signal part
         SetHstFace(mkk[i]);
         mkk[i] -> GetYaxis() -> SetMaxDigits(3);
         mkk[i] -> GetXaxis() -> SetTitleOffset(1.1);
         if ( i/2 < 3 ) { // data
            mkk[i] -> SetLineWidth(1);
            mkk[i] -> SetLineColor(kBlack);
            mkk[i] -> SetOption("E");
            mkk[i] -> SetMarkerStyle(20);
            mkk[i] -> SetMarkerSize(0.5);
         } else { // MC
            mkk[i] -> SetLineWidth(1);
            mkk[i] -> SetLineColor(kGreen+2);
         }
      } else { // sideband
         // auto kCol = kBlue+2;
         auto kCol = kRed; // the same as data vs MC
         mkk[i] -> SetLineWidth(1);
         mkk[i] -> SetLineColor(kCol);
         mkk[i] -> SetFillStyle(3001);
         mkk[i] -> SetFillColor(kCol);
      }
   }

   //-----------------------------------------------------------------
   // Fit signal MC(phi,eta) by Breit-Wigner convoluted with Gauss
   // + constant bg
   constexpr double dL = 0.98;
   constexpr double dU = 1.08;
   double bW = (dU-dL) / mkk[0] -> GetNbinsX();
   auto Lbwg = [bW](const double* x,const double* p) -> double {
      return bW*p[0]*BreitWignerGaussN(x[0],p[1],p[2],p[3]) + p[4];
   };

   // line for BG
   TF1* pl0 = (TF1*)gROOT -> GetFunction("pol0")->Clone();
   pl0 -> SetRange(2*Mk,1.08);
   pl0 -> SetLineWidth(1);
   pl0 -> SetLineColor(kBlue);
   pl0 -> SetLineStyle(kDashed);

   TF1* bwg = new TF1("bwg", Lbwg, dL, dU, 5);
   bwg -> SetParNames("Nphi","Mphi","Gphi","Sigma","Bg");
   bwg -> SetParameters(1.e3, Mphi, Gphi, 1.8e-3, 1.);
   bwg -> SetLineWidth(1);
   bwg -> SetLineColor(kRed);
   bwg -> SetNpx(500);

//    bwg -> FixParameter(1, Mphi);
   bwg -> FixParameter(1, Mphi);
//    bwg -> SetParLimits(1, Mphi-0.02, Mphi+0.02);
   bwg -> FixParameter(2, Gphi);
//    bwg -> SetParLimits(2, Gphi-0.1e-3, Gphi+0.1e-3);
   bwg -> SetParLimits(3, 0.5e-3, 2.5e-3 );
   bwg -> SetParLimits(4, 0., 10. );

   //-----------------------------------------------------------------
   // Draw
   // TCanvas* c1 = new TCanvas("c1","PR",0,0,1100,500); // presentation
   TCanvas* c1 = new TCanvas("c1","PR",0,0,1100,350); // presentation
   int ny = 1;
   // int nx = Nhst/2/ny;
   int nx = 3;

   c1 -> Divide(nx,ny);
   gStyle -> SetOptFit(0);
//    gStyle -> SetOptFit(111); // do not print fixed parameters!
//    gStyle -> SetOptFit(112); // print all parameters (fixed)
//    gStyle -> SetStatX(0.89);
//    gStyle -> SetStatY(0.73);
//    gStyle -> SetStatW(0.22);
//    gStyle -> SetStatH(0.22);
//    gStyle -> SetFitFormat(".3g");

   // vector<TPaveText*> pt( Nhst/2, nullptr );
   vector<TLegend*> leg( nx, nullptr );

   int ih = 0;
   for ( int y = 0; y < ny; ++y ) {
      for ( int x = 0; x < nx; ++x ) {

         int i = y*nx+x;
         c1->cd(i+1);
         mkk[ih] -> Draw("EP");
         // mkk[ih] -> Fit("bwg","QLE","",2*Mk,1.08);

         mkk[ih+1] -> Draw("SAME HIST"); // side-band

         double Ncp = mkk[ih]->GetEntries();
         double Nsb = mkk[ih+1]->GetEntries();

         // normalize MC on data and draw
         int imc = ih + 6;
         double Ndat = Ncp - Nsb;
         double Nmc = mkk[imc]->GetEntries();
         // -mkk[imc+1]->GetEntries();
         double scale = Ndat/Nmc;
         mkk[imc]->Scale(scale);
         mkk[imc]->Draw("SAME HIST"); // MC


         /* for fit only
         // background
//          pl0 -> SetParameter( 0, bwg -> GetParameter(4) );
//          pl0 -> DrawCopy("SAME");

         pt[i] = new TPaveText(0.52,0.59,0.89,0.89,"NDC,NB");
         pt[i] -> SetTextAlign(12);
         pt[i] -> SetTextFont(42);
         pt[i] -> AddText( Form("#scale[1.2]{%s}",
                  titles[i].c_str()) );
         double Nphi = Ncp - Nsb;
         double eNphi = sqrt(Ncp + Nsb);
         if ( x==1 ) { // MC
            Nphi = Ncp;
            eNphi = sqrt(Ncp);
            pt[i] -> AddText( Form(
                     "#color[4]{Nphi(sb) = %.0f #pm %.0f}",
                     Nphi,eNphi) );
            pt[i] -> AddText( Form(
                     "#color[2]{Nphi(fit) = %.0f #pm %.0f}",
                     bwg->GetParameter(0),bwg->GetParError(0)) );
         } else {
            pt[i] -> AddText( Form(
                     "#color[4]{Nphi(sb) = %.1f #pm %.1f}",
                     Nphi,eNphi) );
            pt[i] -> AddText( Form(
                     "#color[2]{Nphi(fit) = %.1f #pm %.1f}",
                     bwg->GetParameter(0),bwg->GetParError(0)) );
         }
//          pt[i] -> AddText( Form("Fit: #chi^{2}/ndf = %.1f #pm %d",
//                   bwg->GetChisquare(),bwg->GetNDF()) );
//          pt[i] -> AddText( Form("Mphi= %.3f #pm %.3f GeV",
//                   bwg->GetParameter(1),bwg->GetParError(1)) );
//          pt[i] -> AddText( Form("Sigma= %.2f #pm %.2f MeV",
//                   1e3*bwg->GetParameter(3),1e3*bwg->GetParError(3)) );
//          pt[i] -> AddText( Form("Bkg= %.2g #pm %.2g",
//                   bwg->GetParameter(4),bwg->GetParError(4)) );
         pt[i] -> Draw();
         */

         leg[i] = new TLegend(0.52,0.69,0.89,0.89);
         string head(Form("#scale[1.2]{#bf{%s}}",titles[i].c_str()));
         leg[i]->SetHeader(head.c_str(),"C");
         leg[i]->AddEntry(mkk[ih],
               Form("Data: %.0f events",Ncp),"EP");
         leg[i]->AddEntry(mkk[ih+1],
               Form("Side-band: %.0f events",Nsb),"F");
         leg[i]->AddEntry(mkk[ih+6],"MC signal #phi#eta","L");
         leg[i]->Draw();

         gPad -> RedrawAxis();

         ih+=2;
      }
   }

   c1 -> Update();
   string pdf("mass_kk_PR.pdf");
   c1 -> Print(pdf.c_str());
}

// {{{1 MAIN:
//-------------------------------------------------------------------------
void mass_KK() {
//-------------------------------------------------------------------------
   gROOT -> Reset();
   gStyle -> SetOptStat(0);
   gStyle -> SetStatFont(62);
   gStyle -> SetLegendFont(42);

   //--------------------------------------------
//    test_BreitWigner();

   //--------------------------------------------
   plot_mass_phi();
}
