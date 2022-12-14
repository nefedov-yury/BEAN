// plot efficiency as a function of X_isr
// -> Xisr_eff_XXXX.pdf

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
vector<TH1D*> get_xisr_histos(string name) {
//--------------------------------------------------------------------
   string fname = "Ntpls/ntpl_mcgpj_" + name + ".root";
   cout << " file: " << fname << endl;

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("SelectKKgg");
   vector<TH1D*> hst;

   // 1)
   TH1D* Xisr_ini = (TH1D*)gROOT->FindObject("mcisr_xisr");
   if( !Xisr_ini ) {
      cerr << " can not find Xisr_ini" << endl;
      exit(0);
   }
   Xisr_ini -> SetTitle(";X_{ISR};Events/0.01");
   SetHstFace(Xisr_ini);
   hst.push_back(Xisr_ini);

   // 2)
   TTree *a4c = (TTree*)gDirectory->Get("a4c");
#include "cuts.h"
   TCut c_here = c_chi2;          // ch2 < 80
   c_here += c_xisr + c_MCmkk;    // X_isr>0.9 && mc_Mkk<1.08
   c_here += c_cpgg;              // central part of Mgg
   c_here += c_phi;               // Mkk in [2*Mk, 1.08GeV]

   TH1D* xf = new TH1D("xf", "after selection", 100,0.,1.);
   a4c->Draw("xisr>>xf",c_here,"goff");
   SetHstFace(xf);
   hst.push_back(xf);

   return hst;
}

//--------------------------------------------------------------------
void PlotEffXisr(string name, string title) {
//--------------------------------------------------------------------
   string pdf = "Xisr_eff_" + name + ".pdf";

   auto hst = get_xisr_histos(name);
   double* xini = hst[0] -> GetIntegral();
   double* xfin = hst[1] -> GetIntegral();
   double Ni = xini[101];
   double Nf = xfin[101];
   cout << " Ni=" << Ni << " Nf=" << Nf << endl;

   // efficiency (%) as a function of cut: [X_{ISR},1];
   TH1D* heff = new TH1D("heff", ";X_{ISR};Efficiency",
         10,0.9,1.0);

   for ( int i=91; i < 100; ++i ) { // do not use bin 100
      double bc = hst[0] -> GetBinCenter(i);
      double si = (1-xini[i])*Ni;
      double sf = (1-xfin[i])*Nf;
      double eff = sf/si;
      heff -> Fill( bc, eff );

//       double err_eff = eff * sqrt(1/sf - 1/si); // binomial error
//       cout << " bc= " << bc << " ini= " << si << " fin= " << sf
//          << " eff= " << eff << " +/- " << err_eff << endl;
//       heff -> SetBinContent(i-90,eff);
//       heff -> SetBinError(i-90,err_eff);
  }

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1 -> cd();
   gPad -> SetGrid();

   SetHstFace(heff);
   heff -> GetYaxis() -> SetTitleOffset(1.25);
   heff -> GetYaxis() -> SetNdivisions(505);
   heff -> SetAxisRange(0.9,0.99,"X");
   heff -> SetAxisRange(0.1,0.2,"Y");
   heff -> DrawCopy("HISTO L");

   TLegend* leg = new TLegend(0.35,0.81,0.89,0.89);
   leg -> SetHeader(title.c_str(),"C");
   leg -> Draw();
   
   c1 -> Update();
   c1 -> Print(pdf.c_str());
}

//----------------------------------------------------------
void xisr_eff()
//----------------------------------------------------------
{
   gROOT->Reset();
   gStyle->SetOptStat(0);

   PlotEffXisr("3080_rs","3080 MeV (R-scan)");
}
