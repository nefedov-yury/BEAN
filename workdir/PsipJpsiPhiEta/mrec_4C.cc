// Here we use special selection with 4C-kinematic constrints
// for data (2009,2012 and continuum) and inclusive MC events.
// The 'Mrec' is drawn and the contribution of
// Psi(2S) -> pi+ pi- phi eta (non J/Psi) is calculated.
// -> mrec_4Cyear.pdf

#include "masses.h"

//--------------------------------------------------------------------
void SetHstFace(TH1* hst) {
//--------------------------------------------------------------------
   TAxis* X = hst -> GetXaxis();
   if ( X ) {
      X -> SetLabelFont(62);
      X -> SetLabelSize(0.04);
      X -> SetTitleFont(62);
      X -> SetTitleSize(0.04);
   }
   TAxis* Y = hst -> GetYaxis();
   if ( Y ) {
      Y -> SetLabelFont(62);
      Y -> SetLabelSize(0.04);
      Y -> SetTitleFont(62);
      Y -> SetTitleSize(0.04);
   }
   TAxis* Z = hst -> GetZaxis();
   if ( Z ) {
      Z -> SetLabelFont(62);
      Z -> SetLabelSize(0.04);
      Z -> SetTitleFont(62);
      Z -> SetTitleSize(0.04);
   }
}

//-------------------------------------------------------------------------
TH1D* get_Mrec(string fname, string hname, vector<double>& num) {
//-------------------------------------------------------------------------
#include "cuts.h"
  // name of folder with root files
   static string dir("prod-12/FourC/");
   fname = dir + fname;
//    cout << " open: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot -> cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory -> Get("a4c");

   TH1D* hst = new TH1D(hname.c_str(),
         ";M^{rec}_{#pi^{#plus}#pi^{#minus }}, GeV/c^{2}"
         ";Entries / 0.001GeV/c^{2}",
         200,3.0,3.2 // 1bin = 1MeV
         );
   hst -> Sumw2(true);

   TCut c_here = c_chi2;
//    c_here += TCut("chsq3g>ch2"); // the same cut as in prod-9
   c_here += c_cpgg;
   c_here += TCut(Form("%f<=Mkk&&Mkk<%f",2*Mk,1.08));
//    cout << " c_here= " << c_here << endl;

   string dr = string("Mrec>>") + hname;
   a4c -> Draw( dr.c_str(), c_here, "goff" );

   // number of events in central part and in side-bands:
   TCut c_Mr0("abs(Mrec-3.097)<0.005"); // [3.092, 3.102]
   TCut c_MrL("Mrec>3.00&&Mrec<3.05");  // [3.00, 3.05] left SB
   TCut c_MrR("Mrec>3.144&&Mrec<3.194");// [3.144, 3.194] right SB
   double n0 = a4c -> Draw( "Mrec", c_here+c_Mr0, "goff" );
   double nL = a4c -> Draw( "Mrec", c_here+c_MrL, "goff" );
   double nR = a4c -> Draw( "Mrec", c_here+c_MrR, "goff" );
   num.push_back(n0);
   num.push_back(nL);
   num.push_back(nR);

   return hst;
}

//-------------------------------------------------------------------------
void mrec_4C() {
//-------------------------------------------------------------------------
#include "norm.h"

   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);

   vector<string> fnames = {
      "data_09psip_all.root",
      "data_12psip_all.root",
      "data_3650_all.root",
      "mcinc_09psip_all.root",
      "mcinc_12psip_all.root",
   };
   vector<string> titles = {
      "Data 2009",
      "Data 2012",
      "E=3.65 GeV",
      "MC incl. 2009",
      "MC incl. 2012",
   };

   int id = 0, ii = 3, is = 5; // 2009
//    int id = 1, ii = 4, is = 6; // 2012
   int ic = 2;  // Continuum data

   string pdf = (id==0) ? string("mrec_4C09") : string("mrec_4C12");
   pdf += string(".pdf");

   TH1D* hst[20];
   vector<double> n_id;
   hst[id] = get_Mrec(fnames.at(id),Form("mrec4c_%i",id),n_id);
   // debug print:
//    cout << " DATA CP: " << n_id[0] 
//       << " SB: L-R " << n_id[1] << " " << n_id[2] << endl;  

   vector<double> n_ii;
   hst[ii] = get_Mrec(fnames.at(ii),Form("mrec4c_%i",ii),n_ii);
//    cout << " MC_I CP: " << n_ii[0]
//       << " SB: L-R " << n_ii[1] << " " << n_ii[2] << endl;  

//    vector<double> n_is;
//    hst[is] = get_Mrec(fnames.at(is),Form("mrec4c_%i",is),n_is);
//    cout << " MC_S SB: L-R " << n_is[1] << " " << n_is[2] << endl;  

   vector<double> n_ic;
   hst[ic] = get_Mrec(fnames.at(ic),Form("mrec4c_%i",ic),n_ic);
//    cout << " CONT  CP: " << n_ic[0]
//       << " SB: L-R " << n_ic[1] << " " << n_ic[2] << endl;  

   hst[id] -> SetLineColor(kBlack); // data
   hst[id] -> SetLineWidth(2);
   hst[id] -> SetMarkerStyle(20);
   hst[id] -> SetMarkerSize(0.7);

   hst[ii] -> SetLineColor(kGreen+2); // inclusive MC
   hst[ii] -> SetLineWidth(2);

//    hst[is] -> SetLineColor(kGreen+2); // MC phi eta
//    hst[is] -> SetLineWidth(2);

   hst[ic] -> SetLineColor(kBlue+1); // Continuum data
   hst[ic] -> SetLineWidth(2);

   // normalize inclusive MC on data
   double ninc = n_ii[1]+n_ii[2];
   double einc = sqrt(ninc);
   double IIscl = ( id==0 ) ? Ito09 : Ito12;
   hst[ii] -> Scale( IIscl );
   for ( int i = 0; i < 3; i++ ) {
      n_ii[i] *= IIscl;
   }
   ninc *= IIscl;
   einc *= IIscl;
//    cout << " norm II: " << ninc << " +/- " << einc << endl;

   // normalize signal MC on data ??
   /*
   double nsig = nums[7]+nums[8];
   double esig = sqrt(nsig);
   double dmax = hst[id]->GetMaximum();
   double mcs_max = hst[is]->GetMaximum();
   hst[is] -> Scale( dmax/mcs_max );
   for ( int i = 6; i < 9; i++) {
         nums[i] *= dmax/mcs_max;
   }
   nsig *= dmax/mcs_max;
   esig *= dmax/mcs_max;
   */

   // normaliza continuum data to data
   double n_cont= hst[ic] -> GetEntries(); 
   double e_cont = sqrt(n_cont);
   bool is_cont = n_cont > 0;
   if ( is_cont ) {
      double ICscl = ( id==0 ) ? Cto09 : Cto12;
      hst[ic] -> Scale( ICscl );
      for ( int i = 0; i < 3; i++ ) {
         n_ic[i] *= ICscl;
      }
      n_cont *= ICscl;
      e_cont *= ICscl;
   }
//    cout << " norm IC: " << n_cont << " +/- " << e_cont << endl;

   // data:
   double ndat = n_id[1]+n_id[2];
   double edat = sqrt(ndat);

   cout << " DATA CP: " << n_id[0] << endl;  
   cout << " DATA SB: " << ndat/10 << " +/- " << edat/10 << endl;  
   cout << " MC_I SB: " << ninc/10 << " +/- " << einc/10 << endl;  
//    cout << " MC_S SB: " << nsig/10 << " +/- " << esig/10 << endl;  
   cout << " CONT ALL: " << n_cont/20 << " +/- " << e_cont/20 << endl;  
   cout << " CONT SB: " << n_ic[1]+n_ic[2] << endl;  

   // draw
   //-----------------------------------------------------------------------

   TLine* lR = new TLine;
   lR -> SetLineColor(kRed+1);
   lR -> SetLineWidth(2);
   lR -> SetLineStyle(7);
   TLine* lB = new TLine;
   lB -> SetLineColor(kMagenta);
   lB -> SetLineWidth(2);
   lB -> SetLineStyle(7);

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1 -> cd();
   gPad -> SetGrid();

   gPad -> SetLogy(true);
   double hmin = 0.5;
   hst[id] -> SetMinimum(hmin);
   double hmax = hst[id] -> GetMaximum();

   SetHstFace(hst[id]);
   hst[id] -> GetYaxis() ->SetTitleOffset(1.25);

   hst[id] -> Draw("E");
   hst[ii] -> Draw("HIST SAME");
//    hst[is] -> Draw("HIST SAME");
   if ( is_cont ) {
      hst[ic] -> Draw("HIST SAME");
   }
   hst[id] -> Draw("E SAME");

   lR -> DrawLine(3.092,hmin,3.092,hmax);
   lR -> DrawLine(3.102,hmin,3.102,hmax);
   lB -> DrawLine(3.001,hmin,3.001,0.1*hmax);
   lB -> DrawLine(3.050,hmin,3.050,0.1*hmax);
   lB -> DrawLine(3.144,hmin,3.144,0.1*hmax);
   lB -> DrawLine(3.194,hmin,3.194,0.1*hmax);

   TLegend* leg = new TLegend(0.58,0.70,0.89,0.89);
   leg -> AddEntry(hst[id],titles[id].c_str(),"PLE");
   if ( is_cont ) {
      leg -> AddEntry(hst[ic],"Continuum data","L");
   }
   leg -> AddEntry(hst[ii],titles[ii].c_str(),"L");
//    leg -> AddEntry(hst[is],titles[is].c_str(),"L");
   leg -> AddEntry(lR,"signal region","L");
   leg -> AddEntry(lB,"side-band","L");
   leg -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   c1 -> Print(pdf.c_str());
}

