// compare cross-sections
//   -> cmpr_cs_XXXX.pdf

// {{{1~ container for measured energies and cross sections
//--------------------------------------------------------------------
struct Data {
   Data() {
      Ebeam.reserve(100);
      Sig.reserve(100);
      ErrEbeam.reserve(100);
      ErrSig.reserve(100);
   }

   void CheckConsistensy() const {
      // check correctness of vectors:
      unsigned int Np = Ebeam.size();
      if( Sig.size() != Np
            // || ErrEbeam.size() != Np || ErrSig.size() != Np
        ) {
         cerr << "FATAL ERROR::Data() "
            "sizes of vectors are different" << endl
            << " Ebeam.size    = " << Ebeam.size() << endl
            << " Sig.size      = " << Sig.size() << endl
            << " ErrEbeam.size = " << ErrEbeam.size() << endl
            << " ErrSig.size   = " << ErrSig.size() << endl
            ;
         exit(EXIT_FAILURE);
      }
   }

   int Size() const {
      return Ebeam.size();
   }

   std::string         version;
   double              Xisr_min;
   std::vector<double> Ebeam;  // GeV
   std::vector<double> Sig;    // pb
   std::vector<double> ErrEbeam;
   std::vector<double> ErrSig;
};

// 29Jan23: KK_fit, latest MC : -> old
//--------------------------------------------------------------------
void SetData29J(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "Picts/cs_29Jan23_results.h"

   dt.version = "29Jan23";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

// 17Mar23: mass_KK_fit (LH), more accurate parameters
//--------------------------------------------------------------------
void SetData17M(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "Picts/cs_17Mar23_results.h"

   dt.version = "17Mar23";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

// 05Feb23: just side-band subtraction, latest MC : -> old
//--------------------------------------------------------------------
void SetData05F(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "Picts/cs_05Feb23_results.h"

   dt.version = "05Feb23";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

// 10Feb23: side-band subtraction, tight cut: 1.01 < Mkk < 1.03 GeV
//--------------------------------------------------------------------
void SetData10F(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "Picts/cs_10Feb23_results.h"

   dt.version = "10Feb23";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

// {{{1~ Main
//--------------------------------------------------------------------
void cmpr_cs() {
//--------------------------------------------------------------------
   gROOT->Reset();

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   // gStyle->SetLegendTextSize(0.05);

   Data dt[2];

   // compare old-new mass_KK_fit (LH)
   // SetData29J(dt[0]); // old
   // SetData17M(dt[1]);
   // string pdf( "cmpr_cs_29J_17M.pdf" ); // <- PDF

   // compare tight cut for Mkk with std range [2*Mk, 1.08GeV]
   // SetData05F(dt[0]); // [2*Mk, 1.08GeV]
   // SetData10F(dt[1]); // 1.01 < Mkk < 1.03 GeV
   // string pdf( "cmpr_cs_tight_cut.pdf" ); // <- PDF

   SetData17M(dt[0]); // mass_KK_fit (LH)
   dt[0].version = "LH-fit of M_{KK} with interference";
   SetData10F(dt[1]); // 1.01 < Mkk < 1.03 GeV
   dt[1].version = "Side-band subtraction";
   string pdf( "cmpr_cs_fit_vs_SB.pdf" ); // <- PDF
   // SetData05F(dt[1]); // Mkk < 1.08 GeV
   // dt[1].version = "Side-band subtraction (std)";
   // string pdf( "cmpr_cs_fit_vs_SB0.pdf" ); // <- PDF


   int Nd = dt[0].Size();
   if ( dt[1].Size() != Nd ) {
      cerr << "FATAL ERROR:: sizes of Data are different" << endl
         << " Size1= " << dt[0].Size()
         << " Size2= " << dt[1].Size() << endl;
      exit(EXIT_FAILURE);
   }

   vector<double> R(Nd);
   double Rmax=10., Rmin=-50.;
   vector<double> Ro(Nd);
   for(unsigned int i = 0; i < Nd; i++) {
      R[i] = 100*(1-dt[0].Sig[i]/dt[1].Sig[i]);
      if ( R[i] > Rmax ) {
         Ro[i] = Rmax-0.1;
      } else if ( R[i] < Rmin ) {
         Ro[i] = Rmin+0.1;
      } else {
         Ro[i] = 1e3; // something invisible
      }
   }

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1->Print((pdf+"[").c_str()); // just open pdf-file

   // 1) crosssections on one graph
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   TGraphErrors* cs[2];
   for ( int i = 0; i < 2; ++i ) {
      cs[i] = new TGraphErrors( Nd,
            dt[i].Ebeam.data(), dt[i].Sig.data(),
            dt[i].ErrEbeam.data(), dt[i].ErrSig.data() );
   }

   TGraph* gr = cs[0];
   string title = string(";center of mass energy (GeV)"
         ";cross-section (pb)");
   gr->SetTitle(title.c_str());
   gr->GetXaxis()->CenterTitle();
   gr->GetYaxis()->CenterTitle();
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   // gr->GetXaxis()->SetTitleOffset(1.1);
   gr->GetXaxis()->SetLabelSize(0.03);
   gr->GetYaxis()->SetLabelSize(0.03);

   gr->SetMaximum(1e4);
   gr->SetMinimum(5.);

   gr->SetMarkerStyle(22);
   gr->SetMarkerColor(kRed+1);
   // gr->SetMarkerSize(0.7);
   gr->Draw("AP");

   gr = cs[1];
   gr->SetMarkerStyle(23);
   gr->SetMarkerColor(kBlue+1);
   // gr->SetMarkerSize(0.7);
   gr->Draw("P"); // on top of pict


   TLegend* leg = new TLegend(0.12,0.72,0.50,0.88);
   leg->AddEntry( cs[0],
         Form("%s",dt[0].version.c_str()), "EP");
         // Form("cross-section %s",dt[0].version.c_str()), "EP");
   leg->AddEntry( cs[1],
         Form("%s",dt[1].version.c_str()), "EP");
         // Form("cross-section %s",dt[1].version.c_str()), "EP");
   leg->Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(false);

   TGraph* rat = new TGraph( Nd, dt[0].Ebeam.data(), R.data() );

   gr = rat;
   string titleR = string(";center of mass energy (GeV)"
         ";relative difference (%)");
   gr->SetTitle(titleR.c_str());
   gr->GetXaxis()->CenterTitle();
   gr->GetYaxis()->CenterTitle();
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   // gr->GetYaxis()->SetTitleOffset(1);
   gr->GetXaxis()->SetLabelSize(0.03);
   gr->GetYaxis()->SetLabelSize(0.03);

   rat->SetMaximum(10.);
   rat->SetMinimum(-50.);

   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kBlue);
   // gr->SetMarkerSize(0.9);
   gr->Draw("AP");

   TGraph* rat1 = new TGraph( Nd, dt[0].Ebeam.data(), Ro.data() );
   rat1->SetMarkerStyle(20);
   rat1->SetMarkerColor(kRed);
   rat1->Draw("P"); // on top of pict

   gPad -> RedrawAxis();
   c1 -> Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // just close pdf-file
}

