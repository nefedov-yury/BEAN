// selection efficiency of the J/Psi -> phi eta
// as a function of true MC M(K+K-)
// fit it by a straight line
// -> eff_sig.pdf

#include "masses.h"

//--------------------------------------------------------------------
constexpr double SQ(double x) {
//--------------------------------------------------------------------
   return x*x;
}

//--------------------------------------------------------------------
void get_eff(string fname, string sdat, string pdf="") {
//--------------------------------------------------------------------
#include "cuts.h"

   bool is_phsp = (fname.find("PHSP") != string::npos); // phase space
   bool is_sig = !is_phsp;

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cout << "can not open " << fname << endl;
      exit(0);
   }
   cout << " file: " << fname << endl;
   froot -> cd("SelectKKgg");

   //-----------------------------------------------------------------
   // get initial numbers:
   //-----------------------------------------------------------------
   TH1D* mc_Mkk_ini = (TH1D*)gROOT -> FindObject("mc_Mkk_ini");
   if( !mc_Mkk_ini ) {
      cout << " can not find mc_Mkk_ini histo" << endl;
      exit(0);
   }
   // 100bins [0.98,108] (see SelectKKgg.cxx)
   double Nini = mc_Mkk_ini -> Integral(1,100);
   cout << " initial number events is " << Nini << endl;

   // re-bin to 50 bins
//    int Nbins = 50;
//    TH1D* mkk_i = dynamic_cast<TH1D*>(mc_Mkk_ini -> Rebin(2,"mkk_i"));

   // re-bin to 25 bins
   int Nbins = 25;
   TH1D* mkk_i = dynamic_cast<TH1D*>(mc_Mkk_ini -> Rebin(4,"mkk_i"));

   double Lphi = 0.98, Uphi = 1.08; // histograms
   if ( mkk_i -> GetNbinsX() != Nbins ) {
      cout << " something wrong with mkk_i: Nbins="
           << mkk_i -> GetNbinsX() << endl;
      exit(0);
   }
   TH1D* mkk_f = new TH1D("mkk_f","",Nbins,Lphi,Uphi);
   mkk_i -> Sumw2(true);
   mkk_f -> Sumw2(true);

   //-----------------------------------------------------------------
   // get final numbers:
   //-----------------------------------------------------------------
   TTree* a4c = (TTree*)gDirectory -> Get("a4c");
   if ( !a4c ) {
      cout << " can not find a4c" << endl;
      exit(0);
   }
#include "a4c_prod3.h"

   double Ncp = 0, eNcp = 0;
   double Nsb = 0, eNsb = 0;

   Long64_t nentries = a4c->GetEntries();
   for ( Long64_t i = 0; i < nentries; ++i ) {
      a4c->GetEntry(i);
      if ( !f_xisr(xisr) ) continue;
      if ( !f_MCmkk(mcmkk) ) continue;
      if ( !f_chi2(ch2) ) continue;
      if ( !f_phi(Mkk) ) continue;

      double w = 1; // no weights
      if ( f_cpgg(Mgg) ) {              // central part
         Ncp  += w;
         eNcp += SQ(w);
         mkk_f -> Fill(mcmkk,w);
      } else if ( f_sbgg(Mgg) ) {       // side-band
         Nsb  += w;
         eNsb += SQ(w);
      }
   }

   // side-band subtraction (OLD)
//    double eff = (Ncp - Nsb) / Nini;
//    double err = eff*sqrt( (eNcp+eNsb)/SQ(Ncp-Nsb) + 1./Nini );
//    cout << " Ncp= " << Ncp << " Nsb= " << Nsb << endl;

   // only central part

//    double nF = mkk_f -> Integral(1,Nbins);
//    cout << " DEBUG: mkk_f(Integ)= " << nF << " Ncp= " << Ncp
//       << endl;

   double eff = Ncp / Nini;
   double err = eff*sqrt( eNcp/SQ(Ncp) + 1./Nini );
   printf("\n integral eff(%s) = %.5f +/- %.5f\n",
         sdat.c_str(),eff,err);

   //-----------------------------------------------------------------
   // Draw results + fit
   //-----------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1 -> cd();
   gPad -> SetGrid();

   TH1D* heff = (TH1D*)mkk_i -> Clone("heff");
   heff -> Divide(mkk_f,mkk_i,1.,1.,"B");

   string title(";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2};efficiency");
   heff -> SetTitle(title.c_str());
   heff -> GetXaxis() -> SetTitleOffset(1.1);
   heff -> GetYaxis() -> SetTitleOffset(1.4);
   heff -> SetLineWidth(2);
   heff -> SetLineColor(kBlack);
   heff -> SetMarkerStyle(20);
   heff -> SetAxisRange(0.05,0.25,"Y");
   heff -> GetYaxis() -> SetNdivisions(504);
   heff -> Draw("EP");

   auto Lf = [](double* x,double* p) -> double {
      return p[0]*(1+p[1]*(x[0]-1.02));
   };
   TF1* ffit = new TF1("ffit", Lf, Lphi, Uphi, 2);
   ffit -> SetParNames("e","k");
   ffit -> SetParameters(0.2, -1.1);
   ffit -> SetLineWidth(1);

   double bin_width = heff -> GetBinWidth(1);
   // left and right fit boundaries:
   double Lfit = (int(2*Mk/bin_width)+1)*bin_width;
   double Ufit = Uphi - bin_width;
   heff -> Fit("ffit","EQ","",Lfit,Ufit);

   TF1* ff2 = new TF1("ff2", Lf, Lphi, Uphi, 2);
   ff2 -> SetParameters(0.2, 0.);
   ff2 -> SetLineColor(kBlue);
   ff2 -> SetLineWidth(2);
   ff2 -> SetLineStyle(kDashed);
   ff2 -> FixParameter(1, 0.);

   heff -> Fit("ff2","EQN","",Lfit,Ufit);
   ff2 -> DrawCopy("SAME");

   double chi2 = ffit -> GetChisquare();
   int ndf = ffit -> GetNDF();
   double a = ffit -> GetParameter(0);
   double ea = ffit -> GetParError(0);
   double b = ffit -> GetParameter(1);
   double eb = ffit -> GetParError(1);

   TPaveText* pt = new TPaveText(0.55,0.70,0.89,0.89,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
   string tit;
   if ( is_sig ) {
      tit = "MC signal #phi#eta, ";
   } else {
      tit = "MC non-#phi K^{+}K^{-}#eta, ";
   }
   tit += sdat;
   pt -> AddText( tit.c_str() );
//    pt -> AddText( "eff = e (1 + k (M-1.02))" );
   pt -> AddText( Form("#chi^{2}/ndf = %.1f / %d",chi2,ndf) );
   pt -> AddText( Form("eff(#phi) = %.4f #pm %.4f",a,ea) );
   pt -> AddText( Form("k= %.2f #pm %.2f",b,eb) );
   pt -> Draw();

   printf(" ++ fit: eff(%s)=e*(1+k*(mkk-1.02)):\n",sdat.c_str());
   printf(" ++ e= %.4f +/- %.4f; k= %.2f +/- %.2f\n",a,ea,b,eb);

   double chi2_0 = ff2 -> GetChisquare();
   int ndf0 = ff2 -> GetNDF();
   double a0 = ff2 -> GetParameter(0);
   double ea0 = ff2 -> GetParError(0);
   double b0 = ff2 -> GetParameter(1);
   double eb0 = ff2 -> GetParError(1);
   if ( fabs(eb0) < 1.e-4 ) {
      printf(" -> e= %.4f +/- %.4f; k= %.2f (fixed)",a0,ea0,b0);
      printf(" chi^2/ndf = %.1f / %d\n\n",chi2_0,ndf0);
   }

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//--------------------------------------------------------------------
void eff_MC() {
//--------------------------------------------------------------------
   gROOT -> Reset();
   gStyle -> SetOptStat(0);
   gStyle -> SetOptFit(0);
   gStyle -> SetLegendFont(42);

//    get_eff("Ntpls/ntpl_mcgpj_3080_rs.root",
//            "3080 R-scan", "eff_sig_3080rc.pdf");
//    get_eff("Ntpls/ntpl_mcgpj_3080_2019.root",
//            "3080 (2019)", "eff_sig_3080_2019.pdf");

//    get_eff("Ntpls/ntpl_mcgpj_3097.root",
//            "3097MeV J/#Psi-scan", "eff_sig_3097.pdf");
//    get_eff("Ntpls/ntpl_mcgpj_J4.root",
//            "3097MeV (2018)", "eff_sig_3097J.pdf");

}
