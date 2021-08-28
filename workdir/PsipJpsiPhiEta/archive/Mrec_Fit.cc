// calculate number of J/Psi for Psi'->J/Psi pi+ pi-
// by fitting of Gaussian plus Polinomial functions
//      -> NFJpsi_???.pdf
//

bool BG_reject;

//----------------------------------------------------------------------
Double_t Gauss(Double_t* x, Double_t* par) {
//----------------------------------------------------------------------
  double Norm  = par[0];
  double Mean  = par[1];
  double Sigma = par[2];
  if( Sigma < 2e-10 ) return 0.;
  double t = (x[0]-Mean)/Sigma;
  return Norm * exp(-0.5*t*t);
}
//----------------------------------------------------------------------
Double_t Background(Double_t* x, Double_t* par) {
//----------------------------------------------------------------------
  double xx = x[0];
  if (BG_reject && fabs(xx-3.097)<0.01) {
    TF1::RejectPoint();
    return 0;
  }

  // 1-st order polinom
//   return par[0] + xx*par[1];

  // 2-nd order polinom
//   return par[0] + xx*(par[1] + xx*par[2]);

  // 3-d order polinom
  return par[0] + xx*(par[1] + xx*(par[2] + xx*par[3]));
}
//----------------------------------------------------------------------
Double_t FitFunction(Double_t *x, Double_t *par) {
//----------------------------------------------------------------------
  // Sum of background and signal
  return Gauss(x,par) + Background(x,&par[3]);
}


//-------------------------------------------------------------------------
TH1D* Mrec_fit(string fname, int midx, int mctrue,  TCanvas* cv)
//-------------------------------------------------------------------------
{
  int MIDX = midx; // MUST BE DEFINED BEFORE "cuts.h"
#include "cuts.h"

  cout << " file: " << fname << endl;
  TFile* froot = TFile::Open(fname.c_str(),"READ");
  if( froot == 0 ) {
    cerr << "can not open " << fname << endl;
    exit(0);
  }

  double Nall = 0;

  froot->cd("PsipJpsiPhiEta");
  if( mctrue == 1 ) {
    TH1D* mc_dec0 = (TH1D*)gDirectory->Get("mc_deccode0");
    Nall = mc_dec0->GetBinContent(66);
    printf(" MC True: Nall= %.0f\n",Nall);
  }

  TTree* nt1 = (TTree*)gDirectory->Get("nt1");

  // 1wjpsi = 10bins => 1 bin =0.2MeV
  TH1D* hmr = new TH1D("hmr",
                        ";M_{rec}(#pi^{+}#pi^{-}), GeV",
                        200,mjpsi-10*wjpsi,mjpsi+10*wjpsi
                     );
  hmr->SetLineWidth(2);
  hmr->SetLineColor(kBlack); // ???

  TCut c_all;
  if( mctrue ) c_all = c_MCtrue;
  nt1->Draw("Mrec>>hmr",c_all,"goff");

  // fit central part by Gauss
  TF1* fgn = new TF1("fgn", Gauss, mjpsi-7*wjpsi,mjpsi+7*wjpsi,3);
  fgn->SetParNames("Norm","Mean","Sigma");
  fgn->SetLineColor(kMagenta+1);
  fgn->SetLineWidth(2);
  fgn->SetParameters(hmr->GetMaximum(), 3.097,1.5e-3);
  hmr->Fit("fgn","0","", mjpsi-wjpsi,mjpsi+wjpsi);

  // background function
  TF1* fbg = new TF1("fbg", Background, mjpsi-10*wjpsi,mjpsi+10*wjpsi,4);
  fbg->SetLineColor(kBlue+1);
  fbg->SetLineWidth(2);
  fbg->SetLineStyle(2);
  BG_reject = true;
  hmr->Fit("fbg","0");

  // fit Gauss + polinom
//   TF1* ft = new TF1("ft", FitFunction, mjpsi-wjpsi,mjpsi+wjpsi,5);
//   ft->SetParNames("Norm","Mean","Sigma","p0","p1");
//   TF1* ft = new TF1("ft", FitFunction, mjpsi-wjpsi*5,mjpsi+wjpsi*5,6);
//   ft->SetParNames("Norm","Mean","Sigma","p0","p1","p2");
  TF1* ft = new TF1("ft", FitFunction, mjpsi-wjpsi*5,mjpsi+wjpsi*5,7);
  ft->SetParNames("Norm","Mean","Sigma","p0","p1","p2","p3");
  double par[] = {1,1,1,1,1,1,1,1,1,1};
  fgn->GetParameters(&par[0]);
  fbg->GetParameters(&par[3]);
  ft->SetParameters(par);
  ft->FixParameter(3,par[3]);
  ft->FixParameter(4,par[4]);
  ft->FixParameter(5,par[5]);
  ft->FixParameter(6,par[6]);

  BG_reject = false;
  hmr->Fit("ft", "SRE","", mjpsi-wjpsi,mjpsi+wjpsi);

  // plot Gauss and background functions
  ft->GetParameters(&par[0]);
  fgn->SetParameters(&par[0]);
  fbg->SetParameters(&par[3]);
  fgn->DrawCopy("SAME");
  fbg->DrawCopy("SAME");

 /*
  double nd = fgn->Integral(mjpsi-wjpsi,mjpsi+wjpsi);
  double err = sqrt(nd); // FIXME ???
//   printf(" N(J/Psi)= %.0f +/- %.0f\n",nd,err);
  printf(" N(J/Psi)= %g +/- %g\n",nd,err);

  if( mctrue == 1 ) {
    double eff = nd / Nall;
    double err_eff = eff*sqrt( (err*err)/(nd*nd) + 1./Nall);
    printf(" eff(J/Psi)= %f +/- %f\n",eff,err_eff);
    hmr->SetTitle(Form(" MC efficiency(J/#Psi)= %.4f +/- %.4f",eff,err_eff) );
  } else {
    hmr->SetTitle(Form(" Number of events N(J/#Psi)= %.0f +/- %.0f",nd,err) );
  }

  // Draw
  cv->cd();
  hmr->DrawCopy();


  TLegend* leg = new TLegend(0.55,0.7,0.89,0.89);
  if( mctrue == 1 ) {
    leg->SetHeader(Form("MC true N(J/#Psi)= %.0f",nd), "C");
  }
  leg->AddEntry(lR, Form("#color[%i]{Central N= %.0f}",kRed+1,Ncp), "L");
  leg->AddEntry(lB, Form("#color[%i]{SideBand N= %.0f}",kBlue+1,Nsb), "L");
  leg->Draw();
*/
  cv->Update();
  return hmr;
}

//-------------------------------------------------------------------------
void Mrec_Fit()
//-------------------------------------------------------------------------
{
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(112);
  gStyle->SetLegendTextSize(0.03);

  string rfiles[] = {
                        "mcinc_09psip_all.root",
                        "ntpl_kkmc_09.root",
                        "data_09psip_all.root",
                        "mcinc_12psip_all.root",
                        "ntpl_kkmc_12.root",
                        "data_12psip_all.root"
                    };
  string pdfs[]     {
                        "NFJpsi_mcinc_09.pdf",
                        "NFJpsi_mcsig_09.pdf",
                        "NFJpsi_data_09.pdf",
                        "NFJpsi_mcinc_12.pdf",
                        "NFJpsi_mcsig_12.pdf",
                        "NFJpsi_data_12.pdf"
                    };

  int isw = 0;      // <----------------------------------

  int midx = 2;            // mc
  if( isw == 2 ) midx = 0; // data09
  if( isw == 5 ) midx = 1; // data12

  bool mc = ( isw % 3 != 2 );
  bool mcsig = ( isw % 3 == 1 );

  TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
  c1->cd();
  gPad->SetGrid();

  string pdf = pdfs[isw];
  c1->Print((pdf+"[").c_str()); // open pdf-file

//   if( mc ) {
//     Mrec_sb_draw(rfiles[isw], midx, 1, c1);
//     c1->Print(pdf.c_str()); // add to pdf-file
//   }

//   if( !mcsig ) {
//     TH1D* hmr = Mrec_sb_draw(rfiles[isw], midx, 0, c1);
//     c1->Print(pdf.c_str()); // add to pdf-file
//   }

//   Mrec_fit(rfiles[isw], midx, 0, c1);
  Mrec_fit(rfiles[isw], midx, 1, c1);
  c1->Print(pdf.c_str()); // add to pdf-file

  c1->Print((pdf+"]").c_str()); // close pdf-file
}
