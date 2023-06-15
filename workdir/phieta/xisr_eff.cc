// plot efficiency as a function of X_isr
// -> Xisr_eff_XXXX.pdf

// {{{1 helper functions and constants
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

// {{{1 Fill histograms
//--------------------------------------------------------------------
tuple<TH1D*,TH1D*> get_xisr_histos(string name) {
//--------------------------------------------------------------------
   string fname = "Ntpls/ntpl_mcgpj_" + name + ".root";
   cout << " file: " << fname << endl;

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("SelectKKgg");

   // ini
   TH1D* Xisr_ini = (TH1D*)gROOT->FindObject("mcisr_xisr");
   if( !Xisr_ini ) {
      cerr << " can not find Xisr_ini" << endl;
      exit(EXIT_FAILURE);
   }
   Xisr_ini->SetTitle(";X_{ISR};Events/0.01");

   // fin
   TTree *a4c = (TTree*)gDirectory->Get("a4c");
#include "cuts.h"
   TCut c_here = c_chi2;          // ch2 < 80
   c_here += c_xisr + c_MCmkk;    // X_isr>0.9 && mc_Mkk<1.08
   c_here += c_cpgg;              // central part of Mgg
   c_here += c_phi;               // Mkk in [2*Mk, 1.08GeV]

   TH1D* xf = new TH1D("xf", "after selection", 100,0.,1.);
   a4c->Draw("xisr>>xf",c_here,"goff");

   return make_tuple(Xisr_ini,xf);
}

// {{{1  PlotEffXisr()
//--------------------------------------------------------------------
void PlotEffXisr(string name, string title) {
//--------------------------------------------------------------------
   string pdf = "Xisr_eff_" + name + ".pdf";

   TH1D* xi = nullptr;
   TH1D* xf = nullptr;
   tie(xi,xf) = get_xisr_histos(name);

   // get commulative integrals normalized to 1
   double* xini = xi->GetIntegral();
   double* xfin = xf->GetIntegral();
   double Ni = xini[101]; // number of entries
   double Nf = xfin[101];
   cout << " Ni=" << Ni << " Nf=" << Nf << endl;

   // efficiency (%) as a function of cut: [X_{ISR},1];
   TH1D* heff = new TH1D("heff", ";X_{ISR};Efficiency",
         10,0.9,1.0);

   for ( int i=91; i < 100; ++i ) { // do not use bin 100
      double bc = xi->GetBinCenter(i);
      double si = (1-xini[i])*Ni;
      double sf = (1-xfin[i])*Nf;
      double eff = sf/si;
      heff->Fill( bc, eff );

      // double err_eff = eff * sqrt(1/sf - 1/si); // binomial error
      // heff->SetBinError(i-90,err_eff);

      cout << " bc= " << bc << " ini= " << si
         << " fin= " << sf << " eff= " << eff
         // << " +/- " << err_eff
         << endl;
   }

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   SetHstFace(heff);
   heff->GetYaxis()->SetTitleOffset(1.25);
   heff->GetYaxis()->SetNdivisions(505);
   heff->SetAxisRange(0.9,0.985,"X");
   heff->SetAxisRange(0.1,0.2,"Y");
   heff->DrawCopy("HISTO L");

   TLegend* leg = new TLegend(0.35,0.81,0.89,0.89);
   leg->SetHeader(title.c_str(),"C");
   leg->Draw();

   c1->Update();
   c1->Print(pdf.c_str());
}

// {{{1 Main
//--------------------------------------------------------------------
void xisr_eff() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);

   PlotEffXisr("3080_rs","R-scan: 3080 MeV");
}
