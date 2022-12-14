// To illustratt why X_{ISR} > 0.9 is chosen as the signal definition
// area, we plot X_{ISR} distributions (from MCGPJ) before and after
// selection
// -> Xisr_XXXX.pdf

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
//    Xisr_ini -> SetTitle(";X_{ISR} = S'/S;Events/0.01");
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
void PlotXisr(string name, string title) {
//--------------------------------------------------------------------
   string pdf = "Xisr_" + name + ".pdf";

   auto hst = get_xisr_histos(name);

   // initial
   hst[0] -> GetYaxis() -> SetTitleOffset(1.2);
   hst[0] -> SetLineWidth(2);
   double ymax = 1.15 * hst[0] -> GetMaximum();
   hst[0] -> SetAxisRange(1.,ymax,"Y");
   // final
   hst[1] -> SetLineWidth(2);
   hst[1] -> SetLineColor(kRed+1);

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   // box to show selection region:
   TBox* box = new TBox;
   box -> SetFillStyle(3001);
   box -> SetFillColor(kGreen-8);
   box -> SetLineColor(kGreen-6);
   box -> SetLineWidth(1);

   hst[0] -> DrawCopy("");
   box -> DrawBox(0.9,1.0,1.0,ymax);
   hst[0] -> DrawCopy("SAME");
   hst[1] -> DrawCopy("SAME");

//    TLegend* leg = new TLegend(0.12,0.6,0.80,0.88);
//    leg->SetHeader("The distribution of X_{ISR} = S'/S","C");
//    leg->AddEntry((TObject*)0,"S  - c.m. energy squared","");
//    leg->AddEntry((TObject*)0,"S' - after ISR process","");
   TLegend* leg = new TLegend(0.12,0.7,0.80,0.88);
   string head = string("X_{ISR} = S'/S for ") + title;
   leg->SetHeader(head.c_str(),"C");
   leg->AddEntry(hst[0],"MCGPJ generated events","L");
   leg->AddEntry(hst[1],"Reconstructed MC events","L");
   leg->AddEntry(box, "Signal definition region","F");
   leg->Draw();

   c1 -> Update();
   c1 -> Print(pdf.c_str());
}

//--------------------------------------------------------------------
void xisr_Pr() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);

   PlotXisr("3080_rs","3080 MeV (R-scan)");
//    PlotXisr("3080","3080 MeV (J/#Psi-scan)");
//    PlotXisr("3080_2019","3080 MeV (2019)");
}
