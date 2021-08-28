// plot M(K+K-) for data, inclusive MC, signal MC and MC KKeta
// cuts (see cuts.h): Mrec + chi^2(4C) + Mgg
// -> mass_kk.pdf

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

//-------------------------------------------------------------------------
TH1D* plot_Mkk(string fname, string hname, int type=0) {
//-------------------------------------------------------------------------
#include "cuts.h"

   // name of folder with root files
   static string dir("prod-9/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

//    const double dU = 1.101;
   const double dU = 1.08;
   const double bL = 0.98; // first bin < dL
   const int Nbins = 100;
   const double bW = (dU-bL)/Nbins; // bin width
   string title(";M^{ inv}_{ K^{+}K^{-}} , GeV/c^{2}");
   title += string(";Entries/") + string(Form("%.0fMeV/c^{2}",bW*1e3));
   TH1D* mkk = new TH1D(hname.c_str(), title.c_str(), Nbins,bL,dU);

   TCut c_here = c_Mrec+c_chi2;
   if ( type == 0 ) {   // central part
      c_here += c_cpgg;
   } else {             // side-band
      c_here += c_sbgg;
   }
   string dr = string("Mkk>>") + hname;
   a4c->Draw(dr.c_str(),c_here,"goff");

   return mkk;
}

//-------------------------------------------------------------------------
void mass_kk() {
//-------------------------------------------------------------------------
// #include "cuts.h"
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);
   gStyle->SetLegendTextSize(0.05);

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
      "MC signal #phi#eta",
      "MC signal #phi#eta",
      "MC non-#phi KK#eta",
      "MC non-#phi KK#eta"
   };

//    int id = 0, ii = 3, is = 5, ib = 7; // 2009
   int id = 1, ii = 4, is = 6, ib = 8; // 2012

   bool PR = true;      // presentation
//    bool PR = false;     // memo

   TH1D* hst[20];
//    vector<int> tmp {id, ii, is, ib};
   vector<int> tmp {id, is, ib};
   for( auto i : tmp )  {
      string hname = string("mkk_") + to_string(i+1);
      hst[i] = plot_Mkk(fnames.at(i),hname,0);
      hst[10+i] = plot_Mkk(fnames.at(i),hname+string("_sb"),1);
   }

   string pdf = (id==0) ? string("mass_kk09") : string("mass_kk12");
   if ( PR ) {
      pdf += "_PR";
   }
   pdf += string(".pdf");

   //-----------------------------------------------------------------------

   TCanvas* c1 = nullptr;
   if ( PR ) {
      c1 = new TCanvas("c1","...",0,0,1000,500); // Pr
      c1->Divide(2,1);
//       c1 = new TCanvas("c1","...",0,0,800,800); // Pr-zoom
   } else {
      c1 = new TCanvas("c1","...",0,0,500,900); // note
      c1->Divide(1,2);
   }

   c1->Print((pdf+"[").c_str()); // open pdf-file

   const double mphi = 1.019461; // 1019.461 +/- 0.019 MeV
   const double Gphi = 4.247e-3; //    4.247 +/- 0.016 MeV
   const double wphi = 5*Gphi;

   for ( int it = 0; it < 11; it+=10 ) {
      hst[id+it]->SetLineWidth(2);
      hst[id+it]->SetLineColor(kBlack); // data
      hst[id+it]->SetMarkerStyle(20);
      hst[id+it]->SetMarkerSize(0.7);

      hst[is+it]->SetLineWidth(2);
      hst[is+it]->SetLineColor(kGreen+2); // MC phi eta

      hst[ib+it]->SetLineColor(kBlue+1);
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

   c1->cd(1);
   gPad->SetGrid();

   SetHstFace(hst[id]);
   hst[id]->GetYaxis()->SetTitleOffset(1.25);

   hst[id]->Draw("E");
   hst[is]->Draw("HIST SAME");
   hst[ib]->Draw("SAME HIST");
   hst[9]->Draw("SAME HIST");
   hst[id]->Draw("E,SAME");

   TLegend* leg = nullptr;
   if ( PR ) {
      leg = new TLegend(0.45,0.45,0.89,0.89);
   } else {
      leg = new TLegend(0.49,0.45,0.89,0.89); // memo
   }
   leg->SetHeader("|M(#gamma#gamma)-M(#eta)| < 0.024", "C");
   leg->AddEntry(hst[id],titles[id].c_str(),"PLE");
   leg->AddEntry(hst[is],titles[is].c_str(),"L");
   leg->AddEntry(hst[ib],titles[ib].c_str(),"L");
   leg->AddEntry(hst[9],"Sum of MC","L");
   leg->Draw();
   gPad->RedrawAxis();

   c1->cd(2);
   gPad->SetGrid();

   SetHstFace(hst[10+id]);

//    TF1* pl0 = (TF1*)gROOT -> GetFunction("pol0")->Clone();
//    pl0 -> SetLineColor(kMagenta);
//    pl0 -> SetLineWidth(2);
//    pl0 -> SetLineStyle(kDashed);

   gStyle -> SetOptFit(0);
//    hst[10+id]->Fit(pl0,"L","E");
   hst[10+id]->Draw("E");

   hst[10+is]->Draw("HIST SAME");
   hst[10+ib]->Draw("SAME HIST");
   hst[19]->Draw("SAME HIST");
   hst[10+id]->Draw("E,SAME");

//    TPaveText* pt = new TPaveText(0.45,0.65,0.89,0.89,"NDC");
   TPaveText* pt = new TPaveText(0.45,0.73,0.89,0.89,"NDC");
   pt -> SetTextAlign(22);
   pt -> SetTextFont(42);
   pt -> AddText("Side-Band");
   pt -> AddText("0.048<|M(#gamma#gamma)-M(#eta)|<0.072");
   pt -> Draw();

//    TLegend* legSB = new TLegend(0.42,0.60,0.89,0.89);
//    legSB -> SetTextFont(62);
//    legSB -> SetTextSize(0.04);
//    legSB -> SetHeader("0.048<|M(#gamma#gamma)-M(#eta)|<0.072", "C");
//    legSB -> AddEntry( hst[10+id], "Side-Band", "PLE" );
//    legSB -> AddEntry(pl0,"LH fit of side-band","L");
//    legSB -> AddEntry( (TObject*)0, Form("p0= %.2f #pm %.2f",
//             pl0 -> GetParameter(0),pl0 -> GetParError(0)),"");
//    legSB -> Draw();

   gPad->RedrawAxis();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // close pdf-file
}

