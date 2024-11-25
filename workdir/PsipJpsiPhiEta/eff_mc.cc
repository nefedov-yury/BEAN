// eff_mc.cc
// selection efficiency of the  J/Psi -> phi eta  as a function
// of true M(K+K-)
//      -> eff_sig_{date}.pdf

#include "RewEtaEff.hpp"    // RewEtaEff() func
#include "RewTrkPiK.hpp"    // RewTrkPi(), RewTrk_K() functions
#include "masses.h"

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

// {{{1 Get initial and final hst: MC-root files MUST contain nt1-tree
//--------------------------------------------------------------------
vector<TH1D*> get_hst(int date)
//--------------------------------------------------------------------
{
#include "cuts.h"
   // Systematic study => I/O check: disable corrections
   bool NoW = false;

   // signal MC files
   string mcsigfile( Form("mcsig_kkmc_%02i.root",date%100) );
   string fname = Dir + mcsigfile;

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }
   printf("file: %s\n",fname.c_str());
   froot->cd("PsipJpsiPhiEta");

   //-----------------------------------------------------------------
   // get initial histograms
   //-----------------------------------------------------------------

   // initial number of Psi' -> pi+ pi- J/Psi -> phi eta
   TH1D* MCdec = (TH1D*)gROOT->FindObject("mc_dcj0");
   if ( !MCdec ) {
      cout << "can not find mc_dcj0" << endl;
      exit(EXIT_FAILURE);
   }
   double Nini = MCdec->GetBinContent(69);
   printf("number of generated 'J/Psi -> phi eta'"
         " (dec# %.0f) is %.0f\n",
          MCdec->GetBinCenter(69),Nini);

   // number of pi+pi-J/Psi after Mrs cut [3.092, 3.102]
   TTree* nt1 = (TTree*)gDirectory->Get("nt1");
   if ( !nt1 ) {
      cout << "can not find nt1" << endl;
      exit(EXIT_FAILURE);
   }

   double Lphi = 0.98, Uphi = 1.08; // histogramms
   int Nbins = 50;
   TH1D* mkk_i = new TH1D(Form("mkk_i_%i",date),"",Nbins,Lphi,Uphi);
   TH1D* mkk_f = new TH1D(Form("mkk_f_%i",date),"",Nbins,Lphi,Uphi);
   TH1D* mkk_sb= new TH1D(Form("mkk_sb_%i",date),"",Nbins,Lphi,Uphi);
   mkk_i->Sumw2(true);
   mkk_f->Sumw2(true);
   mkk_sb->Sumw2(true);

   // use Mrb here
   TCut c_Mrb = TCut("abs(Mrb-3.097)<0.005");
   nt1->Draw( Form("mcmkk>>%s",mkk_i->GetName()), c_Mrb, "goff" );
   double Nppj = mkk_i->Integral();
   printf("number of pi+pi-J/Psi after cut '%s' is %.0f\n",
         c_Mrb.GetTitle(),Nppj);

   //-----------------------------------------------------------------
   // get final histograms
   //-----------------------------------------------------------------

   // J/Psi -> phi eta
   TTree* a4c = (TTree*)gDirectory->Get("a4c");
   if ( !a4c ) {
      cout << " can not find a4c" << endl;
      exit(EXIT_FAILURE);
   }

   // see "a4c_v709.h"
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
   auto cl_Mrec = [](double Mrec) -> bool{
      return fabs(Mrec-3.097) < 0.005;
   };

   // chi^2 cut
   auto cl_chi2 = [](double ch2) -> bool{
      return ch2 < 100;
   };

   // Mphi cut: [dL, dU]
   double dL = 2*Mk; // the left cutoff
   double dU = 1.08; // MUST BE < 1.0835 for signal
   auto cl_phi = [dL,dU](double Mkk) -> bool{
      return ( Mkk > dL && Mkk < dU);
   };

   // Meta: central part
   auto cl_cpgg = [weta](double Mgg) -> bool{
      return fabs(Mgg-Meta) < weta;
   };

   auto cl_sbgg = [weta,shift_eta](double Mgg) -> bool{
      return (fabs(Mgg-Meta) > shift_eta) &&
             (fabs(Mgg-Meta) < shift_eta+weta);
   };

   double Ncp = 0, eNcp = 0;
   double Nsb = 0, eNsb = 0;

   double Weta = 1;
   if ( !NoW ) { // disable corrections for I/O check
      Weta = RewEtaEff(date);
   }

   Long64_t nentries = a4c -> GetEntries();
   for ( Long64_t i = 0; i < nentries; ++i ) {
      a4c->GetEntry(i);
      // if ( !(mcmkk<1.08) ) continue;
      if ( !cl_Mrec(Mrec) ) continue;
      if ( !cl_chi2(ch2) ) continue;
      if ( !cl_phi(Mkk) ) continue;

      // correction for K+,K- eff:
      double w = 1.;
      if ( !NoW ) { // disable corrections for I/O check
         double the_p = acos(Ckp);
         double Ptkp = Pkp*sin(the_p);
         double wp = RewTrk_K(date,Ptkp,+1);
         double the_m = acos(Ckm);
         double Ptkm = Pkm*sin(the_m);
         double wm = RewTrk_K(date,Ptkm,-1);
         double w = wp*wm;
      }
      // correction for eta eff:
      w *= Weta;

      if ( cl_cpgg(Mgg) ) {              // central part
         Ncp  += w;
         eNcp += SQ(w);
         mkk_f->Fill(mcmkk,w);
      } else if ( cl_sbgg(Mgg) ) {       // side-band
         Nsb  += w;
         eNsb += SQ(w);
         mkk_sb->Fill(mcmkk,w);
      }
   }

   vector<TH1D*> result {mkk_i,mkk_f,mkk_sb};
   return result;
}

// {{{1 Efficiency
//--------------------------------------------------------------------
vector<double> get_eff(int date, string pdf, int Cx, int Cy,
      bool fix=false, int shift = 0)
//--------------------------------------------------------------------
{
   vector<TH1D*> hst;
   if ( date != 0 ) {
      hst = get_hst(date);
   } else {
      for ( int d : {2021,2012,2009} ) { // order must be this
         auto res = get_hst(d);
         hst.insert(end(hst),begin(res),end(res));
      }
      for ( size_t ih = 0; ih < 3; ++ih ) {
         hst[ih]->Add(hst[ih+3], 0.764); //2012: 345.4/2259.3*(3./0.6)
         hst[ih]->Add(hst[ih+6], 0.715); //2009: 107.7/2259.3*(3./0.2)
      }
   }

   auto& mkk_i = hst[0];
   auto& mkk_f = hst[1];
   auto& mkk_sb = hst[2];

   TH1D* heff = (TH1D*)mkk_i->Clone("heff");
   // mkk_f->Add(mkk_f,mkk_sb,1.,-1.); // side-band subtraction ?
   heff->Divide(mkk_f,mkk_i,1.,1.,"B");

   // Draw
   auto name = Form("c1_%i",date);
   TCanvas* c1 = new TCanvas(name,name,shift,0,Cx,Cy);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(heff);
   string title(
         ";M^{ inv}_{ K^{#plus}K^{#minus }}, GeV/c^{2}"
         ";efficiency"
         );
   heff->SetTitle(title.c_str());
   heff->GetXaxis()->SetTitleOffset(1.1);
   heff->GetYaxis()->SetNdivisions(504);
   heff->GetYaxis()->SetTitleOffset(1.1);
   heff->SetLineWidth(2);
   heff->SetLineColor(kBlack);
   heff->SetMarkerStyle(20);
   heff->SetAxisRange(0.15,0.55,"Y");
   heff->Draw("EP");

   auto Lf = [](double* x,double* p) -> double {
      return p[0]*(1+p[1]*(x[0]-1.02));
   };

   // left and right fit boundaries
   double bW = heff->GetBinWidth(1);
   int bin1 = heff->FindBin(2*Mk+bW/2);
   int bin2 = heff->FindBin(1.08-bW/2);
   double Lfit = heff->GetBinCenter(bin1) - bW/2;
   double Ufit = heff->GetBinCenter(bin2) - bW/2;

   TF1* ffit = new TF1("ffit", Lf, Lfit, Ufit, 2);
   ffit->SetParNames("e","k");
   ffit->SetParameters(0.3, -1.9);
   ffit->SetLineWidth(1);

   heff->Fit("ffit","EQ","",Lfit,Ufit);
   double chi2 = ffit->GetChisquare();
   int ndf   = ffit->GetNDF();
   double a  = ffit->GetParameter(0);
   double ea = ffit->GetParError(0);
   double b  = ffit->GetParameter(1);
   double eb = ffit->GetParError(1);

   TPaveText* pt = new TPaveText(0.53,0.67,0.892,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   if ( date != 0 ) {
      pt->AddText( Form("MC signal #phi#eta, %i",date) );
   } else {
      pt->AddText( "MC signal #phi#eta" );
   }
   // pt->AddText( "eff = eff(#phi) (1 + k (M-1.02))" );
   pt->AddText( Form("#chi^{2}/ndf = %.1f / %d",chi2,ndf) );
   // pt->AddText( Form("eff(#phi) = %.4f #pm %.4f",a,ea) );
   pt->AddText( Form("#scale[1.2]{#varepsilon}_{0} = %.4f #pm %.4f",
            a,ea) );
   pt->AddText( Form("k = %.2f #pm %.2f",b,eb) );
   pt->Draw();

   string sepline(70,'='); // separation line
   printf("%s\n",sepline.c_str());
   printf(" fit: efficiency(%d) = E_0*(1+k*(mkk-1.02))\n", date);
   printf(" E_0 = %.4f +/- %.4f; k = %.2f +/- %.2f\n", a,ea, b,eb);
   printf("%s\n",sepline.c_str());

   if ( fix ) {
      // show a line with one or both parameters fixed
      TF1* ff2= new TF1("ff2", Lf, Lfit, Ufit, 2);
      ff2->SetParameters(0.3, -1.9);
      ff2->SetLineColor(kBlue);
      ff2->SetLineWidth(2);
      ff2->SetLineStyle(kDashed);

      ff2->FixParameter(1, -1.97); // +/- 0.13
      ff2->FixParameter(0, 0.3380); // +/- 0.0004
      // ff2->FixParameter(0, 0.3385); // NoHC +/- 0.0004

      if ( ff2->GetNumberFreeParameters() > 0 ) {
         heff->Fit("ff2","EQN","",Lfit,Ufit);

         double a0  = ff2->GetParameter(0);
         double ea0 = ff2->GetParError(0);
         double b0  = ff2->GetParameter(1);
         double eb0 = ff2->GetParError(1);
         double chi20 = ff2->GetChisquare();
         int ndf0     = ff2->GetNDF();
         printf("%i) e=%.4f+/-%.4f; k=%.2f(fixed); ch2/ndf=%.1f/%d\n",
               date, a0, ea0, b0,chi20,ndf0);
      } else {
         // calculate chi2
         double ch2 = 0;
         int ndf0 = 0;
         for ( int i=bin1; i<bin2; ++i ) {
            double y = heff->GetBinContent(i);
            double e = heff->GetBinError(i);
            double f = ff2->Eval( heff->GetBinCenter(i) );
            ch2 += SQ((y-f)/e);
            ndf0 += 1;
         }
         printf("%i) e=%.4f & k=%.2f(fixed); ch2/ndf=%.1f/%d\n",
               date,ff2->GetParameter(0),ff2->GetParameter(1),
               ch2,ndf0);
      }

      ff2->DrawCopy("SAME");

      printf("%s\n",sepline.c_str());
   }

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }

   // return vector
   vector<double> ret { double(date),chi2,double(ndf),a,ea,b,eb };
   return ret;
}

// {{{1 Print Table
//-------------------------------------------------------------------------
void PrtTable(int Cx, int Cy)
//-------------------------------------------------------------------------
{
   const vector<int> date {2021, 2012, 2009};
   vector<double> eff_d;
   // bool fix_par = false;
   bool fix_par = true; // show a line with fixed parameters
   for ( size_t i = 0; i < date.size(); ++i ) {
      string pdf( Form("eff_sig_%i.pdf",date[i]) );
      // string pdf( Form("eff_sig_%i_nohc.pdf",date[i]) );
      auto res = get_eff( date[i], pdf, Cx, Cy, fix_par, i*Cx/4 );
      eff_d.insert(end(eff_d),begin(res),end(res));
   }

   // print table: calculate average of k
   string dml(64,'='); // demarcation line
   printf("%s\n",dml.c_str());
   size_t Nd = size(eff_d) / 7;
   if ( Nd == 0 ) {
      return;
   }
   double avk = 0, sigk = 0;
   double ave = 0, sige = 0;
   printf("\\    & $k$            & $\\mathcal{E}_0$ \\\\\n");
   for ( size_t i=0; i < Nd; ++i ) {
      size_t j = 7*i;
      // date & k & eff
      printf("%.0f & $%.2f\\pm%.2f$ & $%.4f\\pm%.4f$ \\\\\n",
            eff_d[j],eff_d[j+5],eff_d[j+6],eff_d[j+3],eff_d[j+4]);
      avk += eff_d[j+5]/SQ(eff_d[j+6]);
      sigk += 1./SQ(eff_d[j+6]);
      ave += eff_d[j+3]/SQ(eff_d[j+4]);
      sige += 1./SQ(eff_d[j+4]);
   }
   avk = avk/sigk;
   sigk = 1./sqrt(sigk);
   ave = ave/sige;
   sige = 1./sqrt(sige);
   printf("%s\n",dml.c_str());
   printf("<k> = %.2f +/- %.2f\n",avk,sigk);
   printf("<eff> = %.4f +/- %.4f\n",ave,sige);
   printf("%s\n",dml.c_str());
}

// {{{1 Main
//-------------------------------------------------------------------------
void eff_mc()
//-------------------------------------------------------------------------
{
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetLegendFont(42);

   //========================================================
   // set the name of the folder with the root files
   Dir = "prod_v709n4/";
   //========================================================

   size_t Cx = 800, Cy = 640; // canvas sizes, X/Y = 1.25

   // Table 6 and Fig.16
   // 1) combined result:
   // get_eff(0,"eff_sig_sum.pdf",Cx,Cy);
   // 2) fitting for each year separately
   // PrtTable(Cx,Cy);

   // 3) Systematic study
   //    3a) => NoHC: no helix corrections
   // Dir = "prod_v709n4/NoHC/";
   // get_eff(0,"eff_sig_sum_nohc.pdf",Cx,Cy);
   //    3b) => weta = N*seta (see cuts.h)
   get_eff(0,"eff_sig_sum_w35.pdf",Cx,Cy);
}
