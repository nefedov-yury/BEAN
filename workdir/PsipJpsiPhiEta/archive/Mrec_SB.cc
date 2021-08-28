// calculate number of J/Psi for Psi'->J/Psi pi+ pi-
// for side-band method
//      -> NJpsi_???.pdf
//
//-------------------------------------------------------------------------
TH1D* Mrec_sb_draw(string fname, int midx, int mctrue,  TCanvas* cv)
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
                        140,mjpsi-7*wjpsi,mjpsi+7*wjpsi
                     );
  hmr->SetLineWidth(2);
  hmr->SetLineColor(kBlack); // ???

  TCut c_all;
  if( mctrue ) c_all = c_MCtrue;
  nt1->Draw("Mrec>>hmr",c_all,"goff");

  double Ncp = nt1->Draw("Mrec",c_all&&c_cp,"goff");
//   cout << " cp cut: " << str_cp << endl;
  printf(" Ncp= %.0f\n",Ncp);

  double Nsb = nt1->Draw("Mrec",c_all&&(c_left||c_right),"goff");
//   cout << " left cut: " << str_left << endl;
//   cout << " right cut: " << str_right << endl;
  printf(" Nsb= %.0f\n",Nsb);

  double nd = Ncp - Nsb;
  double err = sqrt(Ncp + Nsb);
  printf(" N(J/Psi)= %.0f +/- %.0f\n",nd,err);
  if( mctrue == 1 ) {
    double eff = nd / Nall;
    double err_eff = eff*sqrt( (Ncp+Nsb)/(nd*nd) + 1./Nall);
    printf(" eff(J/Psi)= %f +/- %f\n",eff,err_eff);
    hmr->SetTitle(Form(" MC efficiency(J/#Psi)= %.4f +/- %.4f",eff,err_eff) );
  } else {
    hmr->SetTitle(Form(" Number of events N(J/#Psi)= %.0f +/- %.0f",nd,err) );
  }

  // Draw
  cv->cd();
  hmr->DrawCopy();

  TLine* lR = new TLine;
  lR->SetLineColor(kRed+1);
  lR->SetLineWidth(2);
  TLine* lB = new TLine;
  lB->SetLineColor(kBlue+1);
  lB->SetLineWidth(2);

  double ymax=0.7*hmr->GetMaximum();
  lR->DrawLine(mjpsi-wjpsi,0,mjpsi-wjpsi,ymax);
  lR->DrawLine(mjpsi+wjpsi,0,mjpsi+wjpsi,ymax);

  const int n=5;
  lB->DrawLine(mjpsi-n*wjpsi,0,mjpsi-n*wjpsi,ymax);
  lB->DrawLine(mjpsi-(n+1)*wjpsi,0,mjpsi-(n+1)*wjpsi,ymax);
  lB->DrawLine(mjpsi+n*wjpsi,0,mjpsi+n*wjpsi,ymax);
  lB->DrawLine(mjpsi+(n+1)*wjpsi,0,mjpsi+(n+1)*wjpsi,ymax);

  TLegend* leg = new TLegend(0.55,0.7,0.89,0.89);
  if( mctrue == 1 ) {
    leg->SetHeader(Form("MC true N(J/#Psi)= %.0f",nd), "C");
  }
  leg->AddEntry(lR, Form("#color[%i]{Central N= %.0f}",kRed+1,Ncp), "L");
  leg->AddEntry(lB, Form("#color[%i]{SideBand N= %.0f}",kBlue+1,Nsb), "L");
  leg->Draw();

  cv->Update();
  return hmr;
}

//-------------------------------------------------------------------------
void Mrec_SB()
//-------------------------------------------------------------------------
{
  gROOT->Reset();
  gStyle->SetOptStat(0);
//   gStyle->SetOptFit(0);
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
                        "NsbJpsi_mcinc_09.pdf",
                        "NsbJpsi_mcsig_09.pdf",
                        "NsbJpsi_data_09.pdf",
                        "NsbJpsi_mcinc_12.pdf",
                        "NsbJpsi_mcsig_12.pdf",
                        "NsbJpsi_data_12.pdf"
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

  if( mc ) {
    Mrec_sb_draw(rfiles[isw], midx, 1, c1);
    c1->Print(pdf.c_str()); // add to pdf-file
  }

  if( !mcsig ) {
    TH1D* hmr = Mrec_sb_draw(rfiles[isw], midx, 0, c1);
    c1->Print(pdf.c_str()); // add to pdf-file
  }

  c1->Print((pdf+"]").c_str()); // close pdf-file
}
