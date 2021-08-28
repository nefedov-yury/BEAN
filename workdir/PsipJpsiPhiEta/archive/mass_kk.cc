// plot M(K+K-) for data, inclusive MC, signal MC and MC KKeta
// cuts (see cuts.h): Mrec + chi^2(4C) + Mgg
// -> mass_kk.pdf

// global name of folder with root files
static string dir("prod-6/");

//-------------------------------------------------------------------------
TH1D* plot_Mkk(string fname, string hname, int type=0)
//-------------------------------------------------------------------------
{
#include "cuts.h"

   gStyle->SetOptStat(0);

   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");
   TH1D* mkk = new TH1D(hname.c_str(),
                        ";M^{inv}_{K^{+}K^{-}}, GeV",
                        100, 0.98, 1.1  ); //2012

//                         ";M(K^{+}K^{-}) (GeV)",
//                         50, 0.98, 1.1  );// KK 2009

   TCut c_here = c_Mrec+c_chi2;
   if ( type == 0 ) {   // central part
      c_here += c_cpgg;
   } else {             // side-band
      c_here += c_sbgg;
   }
   string dr = string("Mkk>>") + hname;
   a4c->Draw(dr.c_str(),c_here,"goff");

   mkk->GetXaxis()->SetTitleSize(0.04);
   mkk->GetXaxis()->SetTitleOffset(0.9);
   mkk->GetXaxis()->SetLabelFont(62);
   mkk->GetXaxis()->SetLabelSize(0.04);
   mkk->GetYaxis()->SetTitleSize(0.04);
   mkk->GetYaxis()->SetTitleOffset(0.9);
   mkk->GetYaxis()->SetLabelFont(62);
   mkk->GetYaxis()->SetLabelSize(0.04);

   return mkk;
}

//-------------------------------------------------------------------------
void mass_kk()
//-------------------------------------------------------------------------
{
   gROOT->Reset();
   gStyle->SetStatFont(62);
   gStyle->SetLegendTextSize(0.06);

   vector<string> fnames = {
      "data_09psip_all.root",
      "data_12psip_all.root",
      "data_3650_all.root",
      "mcinc_09psip_all.root",
      "mcinc_12psip_all.root",
      "mcsig_kkmc_09.root",
      "mcsig_kkmc_12.root",
      "mckketa_kkmc_09.root",
      "mckketa_kkmc_12.root"
   };
   vector<string> titles = {
      "Data 2009",
      "Data 2012",
      "E=3.65 GeV",
      "MC incl. 2009",
      "MC incl. 2012",
      "MC #phi#eta 2009",
      "MC #phi#eta 2012",
      "MC KK#eta 2009",
      "MC KK#eta 2012"
   };

//    int id = 0, ii = 3, is = 5, ib = 7; // 2009
   int id = 1, ii = 4, is = 6, ib = 8; // 2012

   TH1D* hst[20];
//    vector<int> tmp {id, ii, is, ib};
   vector<int> tmp {id, is, ib};
   for( auto i : tmp )  {
      string hname = string("mkk_") + to_string(i+1);
      hst[i] = plot_Mkk(fnames.at(i),hname,0);
      hst[10+i] = plot_Mkk(fnames.at(i),hname+string("_sb"),1);
   }

   string pdf = (id==0) ? string("mass_kk09.pdf")
                : string("mass_kk12.pdf");

   //-----------------------------------------------------------------------

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->Divide(1,2);
   c1->Print((pdf+"[").c_str()); // open pdf-file

   TLine* lR = new TLine;
   lR->SetLineColor(kBlue+3);
   lR->SetLineWidth(2);
   lR->SetLineStyle(7);

   const double mphi = 1.019461; // 1019.461 +/- 0.019   MeV
   const double Gphi = 4.247e-3; //   4.247 +/- 0.016 MeV
   const double wphi = 5*Gphi; // standard

   for ( int it = 0; it < 11; it+=10 ) {
      hst[id+it]->SetLineWidth(2);
      hst[id+it]->SetLineColor(kBlack);
      hst[id+it]->SetMarkerStyle(20);

      hst[is+it]->SetLineWidth(2);
      hst[is+it]->SetLineColor(kBlue+1);

      hst[ib+it]->SetLineColor(kGreen+2);
      hst[ib+it]->SetLineWidth(2);
   }

   // normalize MC on data
   int n1 = hst[id]->FindBin(mphi-wphi);
   int n2 = hst[id]->FindBin(mphi+wphi);
   int n3 = hst[id]->GetNbinsX();
   double sig_scale = 1.0;
   double bg_scale = 1.0;
   for (int iter = 0; iter < 3; iter++) {
      double sig = 1.0;
      if ( iter==0 ) {
         sig = hst[id]->Integral(n1,n2) / hst[is]->Integral(n1,n2);
      } else {
         sig = (hst[id]->Integral(n1,n2) - hst[ib]->Integral(n1,n2))
                        / hst[is]->Integral(n1,n2);
      }
      hst[is]->Scale( sig );
      sig_scale *= sig;

      double bg = (hst[id]->Integral(n2,n3) - hst[is]->Integral(n2,n3))
                        / hst[ib]->Integral(n2,n3);
      hst[ib]->Scale( bg );
      bg_scale *= bg;
      cout << " iter# " << iter << " sig_scale= " << sig_scale
           << " bg_scale= " << bg_scale << endl;
   }

   // normalize side-band on the same numbers
   hst[10+is]->Scale( sig_scale );
   hst[10+ib]->Scale( bg_scale );

   double ymax = max( hst[id]->GetMaximum(), hst[is]->GetMaximum() );
   hst[id]->SetMaximum(1.2*ymax);

   hst[9]=(TH1D*)hst[is]->Clone("mc_sig_clone");
   hst[9]->Add(hst[ib]); // MC(sig) + MC(bg)
   hst[9]->SetLineColor(kRed+1);

   hst[19]=(TH1D*)hst[10+is]->Clone("mc_sig_sb_clone");
   hst[19]->Add(hst[10+ib]); // MC(sig) + MC(bg)
   hst[19]->SetLineColor(kRed+1);

//   c1->cd(1);
//   gPad->SetGrid();
//   gPad->SetLogy();
//   hst[id]->Draw("E");
//   hst[ii]->Scale(hst[id]->Integral(1,40)/hst[ii]->Integral(1,40));
//   hst[ii]->SetLineColor(kRed+1);
//   hst[ii]->SetLineWidth(2);
//   hst[ii]->Draw("HIST SAME");

//   TLegend* leg1 = new TLegend(0.56,0.6,0.975,0.93);
//   leg1->AddEntry(hst[id],titles[id].c_str(),"LE");
//   leg1->AddEntry(hst[ii],titles[ii].c_str(),"L");
//   leg1->Draw();

//   lR->DrawLine(mphi-wphi,0,mphi-wphi,ymax);
//   lR->DrawLine(mphi+wphi,0,mphi+wphi,ymax);

   c1->cd(1);
   gPad->SetGrid();
   hst[id]->Draw("E");
   hst[is]->Draw("HIST SAME");
   hst[ib]->Draw("SAME HIST");
   hst[9]->Draw("SAME HIST");

//    TLegend* leg = new TLegend(0.56,0.6,0.975,0.93);
   TLegend* leg = new TLegend(0.56,0.40,0.89,0.88);
   leg->SetHeader("| M(#gamma#gamma) - M(#eta) | < 0.024GeV", "C");
   leg->AddEntry(hst[id],titles[id].c_str(),"PLE");
   leg->AddEntry(hst[is],titles[is].c_str(),"L");
   leg->AddEntry(hst[ib],titles[ib].c_str(),"L");
   leg->AddEntry(hst[9],"Sum of MC","L");
   leg->Draw();

//    lR->DrawLine(mphi-wphi,0,mphi-wphi,ymax);
//    lR->DrawLine(mphi+wphi,0,mphi+wphi,ymax);

   c1->cd(2);
   gPad->SetGrid();
   hst[10+id]->Draw("E");
   hst[10+is]->Draw("HIST SAME");
   hst[10+ib]->Draw("SAME HIST");
   hst[19]->Draw("SAME HIST");

   TPaveText* pt = new TPaveText(0.50,0.70,0.88,0.88,"NDC");
   pt->AddText("Side-Band");
   pt->AddText("0.024 < | M(#gamma#gamma) - M(#eta) | < 0.048 GeV");
//    pt->AddText("0.024 < | M_{#gamma#gamma} - M_{#eta} | < 0.048 GeV");
   pt->Draw();

   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   /*
   c1->cd(1);
   gPad->SetGrid();
   gPad->SetLogy();
   hst[id]->Draw("E");
   hst[is]->Draw("HIST SAME");
   hst[ib]->Draw("SAME HIST");
   hst[9]->Draw("SAME HIST");

   lR->DrawLine(mphi-wphi,0,mphi-wphi,ymax);
   lR->DrawLine(mphi+wphi,0,mphi+wphi,ymax);
   leg->Draw();

   c1->cd(2);
   gPad->SetGrid();
   gPad->SetLogy();
   hst[10+id]->Draw("E");
   hst[10+is]->Draw("HIST SAME");
   hst[10+ib]->Draw("SAME HIST");
   hst[19]->Draw("SAME HIST");
   legSB->Draw();

   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file
   */

   c1->Print((pdf+"]").c_str()); // close pdf-file
}
