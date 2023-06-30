// plot true distributions for MCGPJ generator:
// 1) beam-spread       -> mcgpj_bs_Exxxx.pdf
// 2) cos Theta (eta)   -> mcgpj_ang_Exxxx.pdf
// 3) X_{ISR} distributions before and after selection
//    -> Xisr_XXXX.pdf

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

//--------------------------------------------------------------------
TF1* GaussFit(TH1D* hist) {
//--------------------------------------------------------------------
   const double ns = 2.0; // number of sigmas to fit
   TF1* gs = (TF1*)gROOT->GetFunction("gaus");
   gs->SetLineWidth(2);
   gs->SetLineColor(kRed);
   double mean = hist->GetMean();
   double sig  = hist->GetRMS();
   hist->Fit(gs,"Q0","",mean-ns*sig,mean+ns*sig);
   mean = gs->GetParameter(1);
   sig  = gs->GetParameter(2);
   hist->Fit(gs,"Q0","",mean-ns*sig,mean+ns*sig); // 1st fit
   mean = gs->GetParameter(1);
   sig  = gs->GetParameter(2);
   hist->Fit(gs,"Q0","",mean-ns*sig,mean+ns*sig); // 2nd fit
   mean = gs->GetParameter(1);
   sig  = gs->GetParameter(2);
   hist->Fit(gs,"","",mean-ns*sig,mean+ns*sig); // 3d fit
   hist->DrawCopy();
   return gs;
}

// {{{1 get hst
//--------------------------------------------------------------------
TH1D* get_hst(string fname, string hname) {
//--------------------------------------------------------------------
   // name of folder with root files
   static string dir("Ntpls/");

   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( !froot ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("SelectKKgg");
   TH1D* hst = (TH1D*) gDirectory->Get(hname.c_str());
   if( !hst ) {
      cerr << "can not find hst with name " << hname << endl;
      exit(EXIT_FAILURE);
   }

   return hst;
}

//--------------------------------------------------------------------
TH1D* get_xisr_rec(string fname) {
//--------------------------------------------------------------------
#include "cuts.h"
   // name of folder with root files
   static string dir("Ntpls/");

   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( !froot ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("SelectKKgg");
   TTree *a4c = (TTree*)gDirectory->Get("a4c");

   TCut c_here = c_chi2;          // ch2 < 80
   // c_here += c_xisr + c_MCmkk;    // X_isr>0.9 && mc_Mkk<1.08
   c_here += c_MCmkk;             // mc_Mkk<1.08
   c_here += c_cpgg;              // central part of Mgg
   c_here += c_phi;               // Mkk in [2*Mk, 1.08GeV]

   TH1D* xf = new TH1D("xf", "after selection", 100,0.,1.);
   a4c->Draw("xisr>>xf",c_here,"goff");

   return xf;
}

// {{{1 beam spread
//--------------------------------------------------------------------
void BeamSpread(string Ename, double Eee, string pdf) {
//--------------------------------------------------------------------
   string fname( Form("ntpl_mcgpj_%s.root", Ename.c_str()) );
   TH1D* hbs = get_hst(fname,"mcisr_bspr");
   SetHstFace(hbs);
   hbs->SetTitle(
         Form(";(E_{cm} #minus %.3f) MeV; Events/0.01",Eee*1e3) );
   hbs->GetYaxis()->SetMaxDigits(3);
   hbs->GetYaxis()->SetTitleOffset(1.2);
   hbs->SetAxisRange(-3.,3.,"X");

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   gStyle->SetStatX(0.98);
   gStyle->SetStatY(0.87);
   c1->cd();
   auto gs = GaussFit(hbs);

   TLegend* leg = new TLegend(0.62,0.88,0.98,0.98);
   leg->AddEntry(hbs, Form("MC %.3f MeV",Eee*1e3),"LEP");
   leg->AddEntry(gs,"Gaussian fit","L");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 cos(Theta) for eta
//--------------------------------------------------------------------
void AngEta(string Ename, double Eee, string pdf) {
//--------------------------------------------------------------------
   string fname( Form("ntpl_mcgpj_%s.root", Ename.c_str()) );
   TH1D* hst = get_hst(fname,"mc_EtaC");
   SetHstFace(hst);
   hst->SetTitle(";cos(#Theta(#eta));Events/0.02");
   hst->GetYaxis()->SetMaxDigits(3);
   hst->GetYaxis()->SetTitleOffset(1.2);

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetGrid();

   c1->cd();
   auto Lfit = [](const double* x, const double* p) {
      return p[0]*(1+p[1]*x[0]*x[0]);
   };
   TF1* ff = new TF1("ff",Lfit,-1.,1.,2);
   ff->SetParNames("p0","p1");
   ff->SetParameters(1.e3,1.);

   gStyle->SetStatX(0.68);
   gStyle->SetStatY(0.81);
   // hst->Draw("E");
   hst->Draw("E");
   hst->Fit("ff","E");

   TLegend* leg = new TLegend(0.32,0.82,0.68,0.89);
   leg->AddEntry(hst, Form("MC %.3f MeV",Eee*1e3),"LEP");
   leg->AddEntry(ff,"p0(1 #plus p1 cos^{2}#Theta)","L");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 Xisr
//--------------------------------------------------------------------
void Xisr(string Ename, double Eee, string pdf) {
//--------------------------------------------------------------------
   string fname( Form("ntpl_mcgpj_%s.root", Ename.c_str()) );

   // initial
   TH1D* hi = get_hst(fname,"mcisr_xisr");
   SetHstFace(hi);
   hi->SetTitle(";X_{ISR};Events/0.01");
   hi->GetYaxis()->SetMaxDigits(3);
   hi->GetYaxis()->SetTitleOffset(1.2);
   hi->SetLineWidth(2);
   // double ymax = 1.2 * hi->GetMaximum();
   double ymax = 2.5e5; // 250000 ev in MC
   hi->SetAxisRange(1.,ymax,"Y");

   // final
   TH1D* hf = get_xisr_rec(fname);
   SetHstFace(hf);
   hf->SetLineWidth(2);
   hf->SetLineColor(kRed+1);

   TCanvas* c1 = new TCanvas("c1","...",0,0,1000,800);
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   // box to show selection region:
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kGreen-8);
   box->SetLineColor(kGreen-6);
   box->SetLineWidth(1);

   hi->DrawCopy();
   box->DrawBox(0.9,1.0,1.0,ymax);
   hi->DrawCopy("SAME");
   hf->DrawCopy("SAME");

   TLegend* leg = new TLegend(0.12,0.7,0.80,0.88);
   leg->SetHeader(Form("X_{ISR} = S'/S, MC %.3f MeV",Eee*1e3),"C");
   leg->AddEntry(hi, "MCGPJ generated events","L");
   leg->AddEntry(hf, "Reconstructed MC events","L");
   leg->AddEntry(box,"Signal definition region","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 MAIN:
//--------------------------------------------------------------------
void mcgpj_plots() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(112); // print all parameters (fixed)
   // gStyle->SetStatFont(62);


   // string Ename("3080_rs");
   // double Eee = 3.08;
   string Ename("J4");
   double Eee = 3.096986;

   // string pdf_bs( Form("mcgpj_bs_%s.pdf",Ename.c_str()) );
   // BeamSpread(Ename,Eee,pdf_bs);

   // string pdf_ang( Form("mcgpj_ang_%s.pdf",Ename.c_str()) );
   // AngEta(Ename,Eee,pdf_ang);

   string pdf_xisr( Form("Xisr_%s.pdf",Ename.c_str()) );
   Xisr(Ename,Eee,pdf_xisr);
}
