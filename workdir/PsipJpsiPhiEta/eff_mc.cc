// efficiencies for selection of the J/Psi -> phi eta
// as a function of true MC M(K+K-) and
// fit it by a straight line
// -> eff_(sig|bkg)_YYYY.pdf

#include "ReWeightTrkPid_11.h"
#include "ReWeightEtaEff.h"

// {{{1 helper functions and constants
//--------------------------------------------------------------------
constexpr double SQ(double x) {
//--------------------------------------------------------------------
   return x*x;
}

// {{{1 Efficiency: MC-root files MUST contain nt1-tree
//--------------------------------------------------------------------
void get_eff(string fname, string pdf="", bool chkIO = false) {
//--------------------------------------------------------------------
#include "masses.h"

   // date detection from name of file
   int date = 0;
   if ( fname.find("_09") != string::npos ) {
      date = 2009;
   } else if ( fname.find("_12") != string::npos ) {
      date = 2012;
   } else if ( fname.find("_21") != string::npos ) {
      date = 2021;
   } else {
      cerr << "can not detect date in fname: " << fname << endl;
      exit(EXIT_FAILURE);
   }
   bool is_sig = (fname.find("mcsig") != string::npos);

   // name of folder with root files
   // string prod("prod-12/");
   string dir("prod_v709/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }
   printf(" file: %s\n",fname.c_str());
   froot->cd("PsipJpsiPhiEta");

   //-----------------------------------------------------------------
   // get initial numbers:
   //-----------------------------------------------------------------

   // initial number of Psi' -> pi+ pi- J/Psi -> phi eta
   TH1D* MCdec = (TH1D*)gROOT->FindObject("mc_dcj0");
   if ( !MCdec ) {
      cout << " can not find mc_dcj0" << endl;
      exit(EXIT_FAILURE);
   }
   double Nini = MCdec->GetBinContent(69);
   printf(" number of generated 'J/Psi -> phi eta'"
         " (dec# %.0f) is %.0f\n",
          MCdec->GetBinCenter(69),Nini);

   // number of pi+pi-J/Psi after Mrs cut [3.092, 3.102]
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");
   if ( !nt1 ) {
      cout << " can not find nt1" << endl;
      exit(EXIT_FAILURE);
   }

   TCut c_Mrb = TCut("abs(Mrb-3.097)<0.005");

   double Lphi = 0.98, Uphi = 1.08; // histogramms
   int Nbins = 50;
   TCut c_MCmkk("mcmkk<1.08"); // cut for MC generated events
   TH1D* mkk_i = new TH1D("mkk_i","",Nbins,Lphi,Uphi);
   TH1D* mkk_f = new TH1D("mkk_f","",Nbins,Lphi,Uphi);
   TH1D* mkk_sb= new TH1D("mkk_sb","",Nbins,Lphi,Uphi);
   mkk_i -> Sumw2(true);
   mkk_f -> Sumw2(true);
   mkk_sb -> Sumw2(true);

   double Nppj = nt1->Draw("mcmkk>>mkk_i",c_Mrb&&c_MCmkk,"goff");
   printf(" number of pi+pi-J/Psi after cut '%s' is %.0f\n",
         c_Mrb.GetTitle(),Nppj);

   double nI = mkk_i->Integral(1,Nbins);
   if ( fabs(nI-Nppj) > 1e-3 ) {
      printf(" DEBUG: mkk_i(Integ)= %f, Nppj= %f\n",nI,Nppj);
      printf("        scale mkk_i= %f\n",Nppj/nI);
   }

   //-----------------------------------------------------------------
   // get final numbers:
   //-----------------------------------------------------------------

   // J/Psi -> phi eta
   TTree* a4c = (TTree*)gDirectory->Get("a4c");
   if ( !a4c ) {
      cout << " can not find a4c" << endl;
      exit(EXIT_FAILURE);
   }

// #include "a4c_v709.h"
   Double_t        Mrec;
   a4c->SetBranchAddress("Mrec",&Mrec);
   Double_t        ch2;
   a4c->SetBranchAddress("ch2",&ch2);
   Double_t        Pkp;
   a4c->SetBranchAddress("Pkp",&Pkp);
   Double_t        Ckp;
   a4c->SetBranchAddress("Ckp",&Ckp);
   Double_t        Pkm;
   a4c->SetBranchAddress("Pkm",&Pkm);
   Double_t        Ckm;
   a4c->SetBranchAddress("Ckm",&Ckm);
   Double_t        Mkk;
   a4c->SetBranchAddress("Mkk",&Mkk);
   Double_t        Mgg;
   a4c->SetBranchAddress("Mgg",&Mgg);
   Double_t        mcmkk;
   a4c->SetBranchAddress("mcmkk",&mcmkk);

   // Cuts for a4c: (MUST BE THE SAME AS IN "cuts.h"!)
   auto c_Mrec = [](double Mrec) -> bool{
      return fabs(Mrec-3.097) < 0.005;
   };

   // chi^2 cut
   auto c_chi2 = [](double ch2) -> bool{
      return ch2 < 40;
   };

   // Mphi cut: [dL, dU]
   double dL = 2*Mk; // the left cutoff
   double dU = 1.08; // MUST BE < 1.0835 for signal
   auto c_phi = [dL,dU](double Mkk) -> bool{
      return ( Mkk > dL && Mkk < dU);
   };

   // Meta: central part
   const double seta = 0.008;
   const double weta = 3*seta; // standard: 3x,
                               // uncertainties study: 2x, 4x
   auto c_cpgg = [weta](double Mgg) -> bool{
      return fabs(Mgg-Meta) < weta;
   };

   // 'shift_eta' is the start of the side-band
   double shift_eta = 7*seta;
   auto c_sbgg = [weta,shift_eta](double Mgg) -> bool{
      return (fabs(Mgg-Meta) > shift_eta) &&
             (fabs(Mgg-Meta) < shift_eta+weta);
   };

   double Ncp = 0, eNcp = 0;
   double Nsb = 0, eNsb = 0;

   double Weta = 1;
   if ( date <= 2012
         && !chkIO   // disable corrections for I/O check
      ) {
      Weta = ReWeightEtaEff(date);
   }
   Long64_t nentries = a4c -> GetEntries();
   for ( Long64_t i = 0; i < nentries; ++i ) {
      a4c->GetEntry(i);
      if ( !(mcmkk<1.08) ) continue;
      if ( !c_Mrec(Mrec) ) continue;
      if ( !c_chi2(ch2) ) continue;
      if ( !c_phi(Mkk) ) continue;

      // correction for K+,K- eff:
      double wp = 1;
      double wm = 1;
      if ( date <= 2012
            && !chkIO   // disable corrections for I/O check
         ) {
         double Ptkp = Pkp*sin(Ckp);
         wp = ReWeightTrkPid(date,1,Ptkp);
         double Ptkm = Pkm*sin(Ckm);
         wm = ReWeightTrkPid(date,1,Ptkm);
      }
      double w = wp*wm;

      // correction for eta eff:
      w *= Weta;

      if ( c_cpgg(Mgg) ) {              // central part
         Ncp  += w;
         eNcp += SQ(w);
         mkk_f->Fill(mcmkk,w);
      } else if ( c_sbgg(Mgg) ) {       // side-band
         Nsb  += w;
         eNsb += SQ(w);
         mkk_sb->Fill(mcmkk,w);
      }
   }

   //----------------------------------------------------------------------
   // Draw results + fit
   //----------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,1000,800);
   c1->cd();
   gPad->SetGrid();

   TH1D* heff = (TH1D*)mkk_i->Clone("heff");
//    mkk_f -> Add(mkk_f,mkk_sb,1.,-1.); // side-band subtraction ?
   heff->Divide(mkk_f,mkk_i,1.,1.,"B");

   string title(";M^{ inv}_{ K^{#plus}K^{#minus }}, GeV/c^{2}"
         ";efficiency");
   heff->SetTitle(title.c_str());
   heff->GetXaxis()->SetTitleOffset(1.1);
   heff->GetYaxis()->SetTitleOffset(1.4);
   heff->SetLineWidth(2);
   heff->SetLineColor(kBlack);
   heff->SetMarkerStyle(20);
   heff->SetAxisRange(0.1,0.5,"Y");
   heff->Draw("EP");

   auto Lf = [](double* x,double* p) -> double {
      return p[0]*(1+p[1]*(x[0]-1.02));
   };
   TF1* ffit = new TF1("ffit", Lf, Lphi, Uphi, 2);
   ffit->SetParNames("e","k");
   ffit->SetParameters(0.3, -1.9);
   ffit->SetLineWidth(1);

   TF1* ff2 = nullptr;
   if ( pdf.empty() ) {
      // Systematics ONLY:
      ff2= new TF1("ff2", Lf, Lphi, Uphi, 2);
      ff2->SetParameters(0.3, -1.9);
      ff2->SetLineColor(kBlue);
      ff2->SetLineWidth(2);
      ff2->SetLineStyle(kDashed);

      // ff2->FixParameter(1, -1.8); // sys: -1.8 +/- 0.2
      ff2->FixParameter(1, -2.0); // +/- 0.13
   }

   double bin_width = heff->GetBinWidth(1);
   // left and right fit boundaries:
   double Lfit = (int(2*Mk/bin_width)+1)*bin_width;
   double Ufit = Uphi - bin_width;
   heff->Fit("ffit","EQ","",Lfit,Ufit);
   if ( ff2 ) {
      heff->Fit("ff2","EQN","",Lfit,Ufit);
      ff2->DrawCopy("SAME");
   }
   double chi2 = ffit->GetChisquare();
   int ndf   = ffit->GetNDF();
   double a  = ffit->GetParameter(0);
   double ea = ffit->GetParError(0);
   double b  = ffit->GetParameter(1);
   double eb = ffit->GetParError(1);

   TPaveText* pt = new TPaveText(0.59,0.70,0.89,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   string tit;
   if ( is_sig ) {
      tit = "MC signal #phi#eta, ";
   } else {
      tit = "MC non-#phi K^{+}K^{-}#eta, ";
   }
   tit += to_string(date);
   pt->AddText( tit.c_str() );
   // pt->AddText( "eff = eff(#phi) (1 + k (M-1.02))" );
   pt->AddText( Form("#chi^{2}/ndf = %.1f / %d",chi2,ndf) );
   pt->AddText( Form("eff(#phi) = %.4f #pm %.4f",a,ea) );
   pt->AddText( Form("k = %.2f #pm %.2f",b,eb) );
   pt->Draw();

   string sepline(70,'='); // separation line
   printf("%s\n",sepline.c_str());
   printf(" fit: efficiency(%d) = e*(1+k*(mkk-1.02))\n", date);
   printf(" e = %.4f +/- %.4f; k = %.2f +/- %.2f\n", a,ea, b,eb);
   printf("%s\n",sepline.c_str());

   if ( ff2 ) {
      double a0  = ff2->GetParameter(0);
      double ea0 = ff2->GetParError(0);
      double b0  = ff2->GetParameter(1);
      double eb0 = ff2->GetParError(1);
      printf(" %i) e= %.4f +/- %.4f; k= %.2f (fixed)\n",
            date, a0, ea0, b0);
      printf("%s\n",sepline.c_str());
   }

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 Main
//-------------------------------------------------------------------------
void eff_mc() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetLegendFont(42);

   // fig.16:
   // get_eff("mcsig_kkmc_09.root","eff_sig_2009.pdf");
   // get_eff("mcsig_kkmc_12.root","eff_sig_2012.pdf");
   // get_eff("mcsig_kkmc_21.root","eff_sig_2021.pdf");

   // no pdf-files, fit with k-fixed
   // get_eff("mcsig_kkmc_09.root");
   // get_eff("mcsig_kkmc_12.root");
   // get_eff("mcsig_kkmc_21.root");

   // Systematic study => NoHC: no helix corrections
   // get_eff("NoHC/mcsig_kkmc_09.root");
   // get_eff("NoHC/mcsig_kkmc_12.root");

   // Systematic study => I/O check: disable corrections
   // get_eff("mcsig_kkmc_09.root","",true);
   // get_eff("mcsig_kkmc_12.root","",true);

   // KKeta phase space MC ???
   // get_eff("mckketa_kkmc_09.root");
   // get_eff("mckketa_kkmc_12.root");
}
