// estimating the cross-section values from the data:
// I ) simple subtraction of the side-band
// II) after mass_KK_fit.cc
//  *) plot efficiency of phi-eta selection on the base of MCGPJ
//     -> eff_12/18/R/all.pdf
//  *) plot/print cross-section for each energy point
//    -> cp_12/18/R/all.pdf
//    -> cs_results.h   : cpp-code with cross-section and energies
//    -> cs_results.tex or .txt : final tables in TeX format

#include <cstdio>
#include <time.h>   // see man strftime
#include "masses.h"

// {{{1 helper functions
//--------------------------------------------------------------------
constexpr double SQ(double x) {
//--------------------------------------------------------------------
   return x*x;
}

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

// {{{1 print tables functions
//--------------------------------------------------------------------
void print_eff( const vector<double>& eff, string title ) {
//--------------------------------------------------------------------
   int n = eff.size();
   printf("   // Efficiency for %s\n",title.c_str());
   printf("   vector<double> eff {\n");
   for(int i = 0; i < n; i += 5) {
      printf("         ");
      int j = i;
      for(; j < (i+5<n ? i+5 : n); j++) {
         printf(" %8.5f,",eff[j]);
      }
      printf("\n");
   }
   printf("   };\n");
}

//--------------------------------------------------------------------
void prtXcpp( FILE* fp, string Tinfo ) {
//--------------------------------------------------------------------
#include "cuts.h"
   const char* shift = "   "; // shift = 3 spaces
   fprintf(fp,"%s// %s",shift,Tinfo.c_str());
   fprintf(fp,"\n");
   fprintf(fp,"%sXisr_min = %.3f;\n", shift,Xisr_min);
}

//--------------------------------------------------------------------
void prtCScpp( FILE* fp, string title, string suf,
      const vector<double>& E, const vector<double>& erE,
      const vector<double> cs[] ) {
//--------------------------------------------------------------------
   const char* shift = "   "; // shift = 3 spaces
   const char* comma = ", ";
   int n = E.size();
   fprintf(fp,"\n");
   fprintf(fp,"%sstring ver%s(\"%s\");\n",
         shift,suf.c_str(),title.c_str());
   fprintf(fp,"%svector<double> Ebeam%s = { // GeV\n",
         shift,suf.c_str());
   for(int i = 0; i < n; i += 5) {
      fprintf(fp,"%s%s",shift,shift);
      for(int j = i; j < (i+5<n ? i+5 : n); j++) {
         fprintf(fp,"%.6f%s",E[j]*1e-3, (j!=n-1) ? comma : "" );
      }
      fprintf(fp,"\n");
   }
   fprintf(fp,"%s};\n",shift);

   fprintf(fp,"%svector<double> ErrEbeam%s = { // GeV\n",
         shift,suf.c_str());
   for(int i = 0; i < n; i += 5) {
      fprintf(fp,"%s%s",shift,shift);
      for(int j = i; j < (i+5<n ? i+5 : n); j++) {
         fprintf(fp,"%.6f%s",erE[j]*1e-3, (j!=n-1) ? comma : "" );
      }
      fprintf(fp,"\n");
   }
   fprintf(fp,"%s};\n",shift);

   fprintf(fp,"%svector<double> Sig%s = { // cross-section, pb\n",
         shift,suf.c_str());
   for(int i = 0; i < n; i += 5) {
      fprintf(fp,"%s%s",shift,shift);
      for(int j = i; j < (i+5<n ? i+5 : n); j++) {
         fprintf(fp,"%.5e%s",cs[0][j], (j!=n-1) ? comma : "" );
      }
      fprintf(fp,"\n");
   }
   fprintf(fp,"%s};\n",shift);

   fprintf(fp,"%svector<double> ErrSig%s = { //error, pb\n",
         shift,suf.c_str());
   for(int i = 0; i < n; i += 5) {
      fprintf(fp,"%s%s",shift,shift);
      for(int j = i; j < (i+5<n ? i+5 : n); j++) {
         fprintf(fp,"%.5e%s",cs[1][j], (j!=n-1) ? comma : "" );
      }
      fprintf(fp,"\n");
   }
   fprintf(fp,"%s};\n",shift);
   fprintf(fp,"\n");
}

//--------------------------------------------------------------------
void prtTable( FILE* fp, bool TeX, string title,
      const vector<double>& e, const vector<double>& lumi,
      const vector<double>& eff,
      const vector<double> nd[], const vector<double> cs[] ) {
//--------------------------------------------------------------------
   int n = e.size();
   fprintf(fp,"\n");
   fprintf(fp,"  %s\n", title.c_str());
   if ( TeX ) {
      const char* pm = "$\\pm$";
      const char* ca = "&";
      const char* ce = "\\\\";
      fprintf(fp,
         "  E(MeV)& Lum.(pb$^{-1}$)& Signal&    Eff.(\\%%)&"
         "  Cross Section(nb)\\\\\n");
      for ( int i = 0; i < n; ++i ) {
         fprintf(fp,
               " %-9.3f%s %8.6g%s  %5.1f%s%4.1f%s  %6.2f%s"
               " %6.3f%s%5.3f%s\n",
               e[i],ca,lumi[i],ca,nd[0][i],pm,nd[1][i],ca,
               eff[i]*100,ca,cs[0][i]*1e-3,pm,cs[1][i]*1e-3,ce);
      }
   } else {
      const char* pm = " +/- ";
      const char* ca = " ";
      fprintf(fp,
         "  E(GeV)    Lum.(pb-1)     Signal      Eff.(%%)"
         "   Cross Section(nb)\n");
      for ( int i = 0; i < n; ++i ) {
         fprintf(fp,
               " %-9.6f%s %8.6g%s %5.1f%s%4.1f%s  %6.2f%s"
               "   %6.4f%s%6.4f\n",
               e[i]*1e-3,ca,lumi[i],ca,nd[0][i],pm,nd[1][i],ca,
               eff[i]*100,ca,cs[0][i]*1e-3,pm,cs[1][i]*1e-3);
      }
      fprintf(fp,"\n");
   }
}

// {{{1 efficiency && cross-sections
//--------------------------------------------------------------------
double getMCini(string fname) {
//--------------------------------------------------------------------
// return number of MC events generated for Xisr > 0.9 and Mkk < 1.08
// NOTE: there is no cut Mkk < 1.08 in "mcisr_xisr" histo
//       so we use "mc_Mkk_ini"

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cout << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("SelectKKgg");
   TH1D* mc_Mkk_ini = (TH1D*)gROOT->FindObject("mc_Mkk_ini");
   if( !mc_Mkk_ini ) {
      cout << " can not find mc_Mkk_ini histo" << endl;
      exit(EXIT_FAILURE);
   }
   // Note: Xisr > 0.9 cut is applied in SelectKKgg.cxx
   double Nini = mc_Mkk_ini->Integral(1,100); // Mkk<1.08GeV
   // printf(" number of generated 'e+e- -> phi eta'"
         // " for Xisr>0.9 is %.0f\n",Nini);

   return Nini;
}

//--------------------------------------------------------------------
tuple<double,double> getNumHst( string fname, string hname,
      bool isMC, bool UseFitMkk, TH1D* hst[] ) {
//--------------------------------------------------------------------
// fill Mkk histograms for the central and side-band regions and
// return number of events in them
#include "cuts.h"

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cout << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("SelectKKgg");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   TCut c_here = c_chi2;
   if ( isMC ) {
      c_here += c_xisr + c_MCmkk; // X_isr>0.9 && mc_Mkk<1.08
   }
   if ( UseFitMkk ) {
      c_here += c_phi;     // [2*Mk, 1.08GeV]
   } else {
      c_here += c_phiT;     // [1.01, 1.03GeV]
   }

   // binning
   double dU = 1.08;
   double bL = 0.98; // first bin < dL
   int Nbins = 50;
   if ( isMC ) {
      Nbins = 100;
   }
   const double bW = (dU-bL)/Nbins; // bin width
   string title(";M^{ inv}_{ K^{#plus}K^{#minus}}, GeV/c^{2}");
   title += ";Entries/" + string(Form("%.0fMeV/c^{2}",bW*1e3));
   string hncp = hname + "_cp";
   string hnsb = hname + "_sb";
   hst[0] = new TH1D(hncp.c_str(), title.c_str(), Nbins,bL,dU);
   hst[1] = new TH1D(hnsb.c_str(), title.c_str(), Nbins,bL,dU);

   // fill MKK histograms and count number of events in them
   double Ncp = a4c->Draw( Form("Mkk>>%s",hncp.c_str()),
         c_here+c_cpgg, "goff" );
   double Nsb = a4c->Draw( Form("Mkk>>%s",hnsb.c_str()),
         c_here+c_sbgg, "goff" );

   return make_tuple(Ncp,Nsb);
}

//--------------------------------------------------------------------
void getEfficiency( bool UseFitMkk, const vector<string>& names,
      vector<double> Eff[], vector<TH1D*>& HstMc ) {
//--------------------------------------------------------------------
// calculate efficiency and fill Nmc & mkk
   int N = names.size();
   Eff[0].resize(N,0.);
   Eff[1].resize(N,0.);
   HstMc.resize(2*N,nullptr);

   for ( int i = 0; i < N; ++i ) {
      string fname = "Ntpls/ntpl_mcgpj_" + names[i] + ".root";
      double Nini = getMCini(fname);

      string hname = "mkk_mc_" + names[i];
      double Ncp = 0, Nsb = 0;
      tie(Ncp,Nsb) = getNumHst( fname, hname,
            true, UseFitMkk, &HstMc[2*i]);

      double eff = (Ncp - Nsb) / Nini;
      double err = eff*sqrt( (Ncp+Nsb)/SQ(Ncp-Nsb) + 1./Nini);
      // printf("%s efficiency: %.5f +/- %.5f\n",
            // names[i].c_str(),eff,err);

      Eff[0][i] = eff;
      Eff[1][i] = err;
   }
}

//--------------------------------------------------------------------
void getCrossSection( const vector<string>& names,
      const vector<double>& lumi,const vector<double>& eff,
      vector<double> CS[],vector<double> Nd[],vector<TH1D*>& HstD ) {
//--------------------------------------------------------------------
// calculate cross sections by couning events and subtracting
// the side-band
   const double Br_phiKK = 0.492; // PDG2018: phi->K+K-   (49.2±0.5)%
   // PDG2018: eta->2gamma (39.41±0.20)%
   // MC generated decay of eta

   int N = names.size();
   CS[0].resize(N,0.);
   CS[1].resize(N,0.);
   Nd[0].resize(N,0.);
   Nd[1].resize(N,0.);
   HstD.resize(2*N,nullptr);

   for ( int i = 0; i < N; ++i ) {
      string fname = "Ntpls/ntpl_" + names[i] + ".root";
      string hname = "mkk_dat_" + names[i];
      double Ncp = 0, Nsb = 0;
      // here UseFitMkk must be false
      tie(Ncp,Nsb) = getNumHst(fname,hname,false,false,&HstD[2*i]);

      double nd      = Ncp - Nsb;
      double err_nd  = sqrt( Ncp + Nsb );
      double cs      = nd / (lumi[i]*eff[i]*Br_phiKK);
      double err_cs  = cs * (err_nd/nd);
      // printf("%s cross section: %.3g +/- %.3g\n",
            // names[i].c_str(),cs,err_cs);

      Nd[0][i] = nd;
      Nd[1][i] = err_nd;
      CS[0][i] = cs;
      CS[1][i] = err_cs;
   }
}

//--------------------------------------------------------------------
void read_mkk_file(const vector<string>& names,
      vector<double> Nd[], string DIR) {
//--------------------------------------------------------------------
// read Mkk fitting results: Nphi and errNphi
   // directory to read files
   // const string DIR("mkk_inter/");
   int N = names.size();
   Nd[0].resize(N,0.);
   Nd[1].resize(N,0.);
   for ( int i = 0; i < N; ++i ) {
      string fname = DIR + "mkk_" + names[i] + ".txt";
      FILE* fmkk = fopen(fname.c_str(),"r");
      if ( !fmkk ) {
         cout << " can not open file " << fname << endl;
         exit(EXIT_FAILURE);
      }
      fseek(fmkk, 0, SEEK_END);
      auto pos = ftell(fmkk);
      int count = 0;
      while ( pos ) { // search for \n
         fseek(fmkk, --pos, SEEK_SET);
         if (fgetc(fmkk) == '\n') {
            if (count++ == 1) break;
         }
      }
      // gets last line
      char line[100];
      while ( fgets(line, sizeof(line), fmkk ) != NULL ) {
        string sl(line);
        if ( sl.find("Nphi") != string::npos ) {
           double Nphi = 0, err_Nphi = 0;
           int ret =
              sscanf(line,"      Nphi= %lf \\pm %lf",&Nphi,&err_Nphi);
           if ( ret != 2 ) {
              cout << "-- ERR: " << sl << endl;
              exit(1);
           }
           // cout << fname << " ++OK Nphi= " << Nphi << endl;
           Nd[0][i] = Nphi;
           Nd[1][i] = err_Nphi;
           break;
        } else {
           cout << "- ERR: " << sl << endl;
           exit(1);
        }
      }
      fclose(fmkk);
   }
}

//--------------------------------------------------------------------
void CrossSection( const vector<double> Nd[],
      const vector<double>& lumi, const vector<double>& eff,
      vector<double> CS[] ) {
//--------------------------------------------------------------------
// calculate cross sections
   const double Br_phiKK = 0.492; // PDG2018: phi->K+K-   (49.2±0.5)%
   // PDG2018: eta->2gamma (39.41±0.20)%
   // MC generated decay of eta

   int N = Nd[0].size();
   CS[0].resize(N,0.);
   CS[1].resize(N,0.);

   for ( int i = 0; i < N; ++i ) {
      double nd      = Nd[0][i];
      double err_nd  = Nd[1][i];
      double cs      = nd / (lumi[i]*eff[i]*Br_phiKK);
      double err_cs  = cs * (err_nd/nd);
      // cout << " cross section: " << cs << " +/- " << err_cs << endl;

      CS[0][i] = cs;
      CS[1][i] = err_cs;
   }
}

// {{{1  draw the histos in pdf
//--------------------------------------------------------------------
TGraphErrors* get_graph( const vector<double>& E,
      bool isMC, vector<double> res[] ){
//--------------------------------------------------------------------
   int N = E.size();
   vector<double> err_E(N,0.); // for graph
   string title("e^{#plus}e^{#minus} #rightarrow #phi#eta;E, MeV");
   if ( isMC ) {
      title += string(";efficiency");
   } else {
      title += string(";cross section, pb");
   }

   TGraphErrors* geff = new TGraphErrors(N,E.data(),res[0].data(),
         err_E.data(),res[1].data());
   geff->SetTitle( title.c_str() );
   return geff;
}

//--------------------------------------------------------------------
void drawHstRes( string pdf, const vector<double>& E,
      const vector<TH1D*>& Hst, bool isMC, vector<double> res[] ){
//--------------------------------------------------------------------
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->Print((pdf+"[").c_str()); // just open pdf-file
   c1->Divide(2,2);
   int j = 1, jmax = 4;

   int N = E.size();
   for( int i = 0; i < N; ++i ) {
      c1->cd(j);
      TH1D* hst = Hst[2*i];
      SetHstFace(hst);
      hst->GetYaxis()->SetTitleOffset(1.3);
      hst->SetLineWidth(1);
      hst->SetLineColor(kBlack);
      // hst->SetMarkerStyle(20);
      hst->Draw("E");

      hst = Hst[2*i+1];
      SetHstFace(hst);
      hst->SetLineWidth(1);
      // hst->SetLineColor(kBlue+2);
      hst->SetLineColor(kRed);
      hst->SetFillStyle(3001);
      hst->SetFillColor(kRed);
      hst->Draw("SAME");

      TLegend* leg = new TLegend(0.50,0.65,0.89,0.89);
      leg->SetTextFont(42);
      leg->SetHeader(Form("E = %.3f MeV",E[i]),"C");

      if ( isMC ) {
         leg->AddEntry(Hst[2*i],"MC signal #phi#eta","PLE");
         leg->AddEntry(Hst[2*i+1],"Side band", "F");
         leg->AddEntry( (TObject*)0,
               Form("#varepsilon = %.4f #pm %.4f",
                  res[0][i],res[1][i]),"");
      } else {
         leg->AddEntry(Hst[2*i],"Data","PLE");
         leg->AddEntry(Hst[2*i+1],"Side band", "F");
         leg->AddEntry( (TObject*)0,
               Form("#sigma = %.3f #pm %.3f pb",
                  res[0][i],res[1][i]), "");
      }
      leg->Draw();
      gPad->RedrawAxis();

      j += 1;
      if ( j > jmax ) {
         c1->Update();
         c1->Print(pdf.c_str()); // add to pdf-file
         j = 1;
      }
   }

   // clear empty pads
   for ( ;1 < j && j <= jmax; ++j ) {
      c1->cd(j);
      gPad->Clear();
      c1->Update();
      c1->Print(pdf.c_str()); // add to pdf-file
   }

   c1->cd();
   gPad->Clear();
   gPad->SetGrid();
   if ( !isMC ) {
      gPad->SetLogy();
   }

   TGraphErrors* g = get_graph( E, isMC, res );
   g->GetYaxis()->SetTitleOffset(1.);
   g->SetMarkerStyle(20);
   g->Draw("APL");

   gPad->RedrawAxis();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // just close pdf-file
}

//--------------------------------------------------------------------
void drawGraph( string pdf, vector<TGraphErrors*>& g, bool isMC ) {
//--------------------------------------------------------------------
   string cn = pdf.substr(0,pdf.find('.'));
   TCanvas* c1 = new TCanvas(cn.c_str(),"...",0,0,800,800);
   c1->Print((pdf+"[").c_str()); // just open pdf-file

   c1->cd();
   gPad->SetGrid();
   if ( !isMC ) {
      gPad->SetLogy();
   }

   // g[0]->GetXaxis()->SetLimits(2880.,3140.);
   g[0]->GetXaxis()->SetLimits(3030.,3140.);
   if ( isMC ) {
      g[0]->SetMinimum(0.14);
      g[0]->SetMaximum(0.18);
      g[0]->GetYaxis()->SetNdivisions(1005);
      g[0]->GetYaxis()->SetTitleOffset(1.3);
   } else {
      g[0]->SetMinimum(5.);
      g[0]->SetMaximum(5000.);
      g[0]->GetYaxis()->SetTitleOffset(1.2);
   }

   g[0]->SetMarkerStyle(20);
   g[0]->SetMarkerColor(kBlue+2);
   g[0]->SetLineColor(kBlue+2);
   g[0]->SetLineWidth(1);
   g[0]->SetLineStyle(kSolid); // kDashed;
   g[0]->Draw("APL");

   g[1]->SetMarkerStyle(33);
   g[1]->SetMarkerColor(kGreen+1);
   g[1]->SetLineColor(kGreen+1);
   g[1]->SetLineWidth(1);
   g[1]->SetLineStyle(kSolid); // kDashed;
   g[1]->Draw("PL SAME");

   g[2]->SetMarkerStyle(22);
   g[2]->SetMarkerColor(kRed+2);
   g[2]->SetLineColor(kRed+2);
   g[2]->SetLineWidth(1);
   g[2]->SetLineStyle(kSolid); // kDashed;
   g[2]->Draw("PL SAME");

   g[0]->Draw("P SAME");

   TLegend* leg = new TLegend(0.15,0.65,0.50,0.89);
   if ( isMC ) {
      leg->SetHeader("Efficiency","C");
   } else {
      leg->SetHeader("Cross Section","C");
   }
   leg->AddEntry(g[0], "J/#Psi scan 2012", "PLE");
   leg->AddEntry(g[1], "2018", "PLE");
   leg->AddEntry(g[2], "R-scan", "PLE");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // just close pdf-file
}

// {{{1  J/Psi scan 2012
//--------------------------------------------------------------------
void setJpsiScan12( vector<string>& names, vector<double>& ebeam,
      vector<double>& er_eb, vector<double>& lumi ) {
//--------------------------------------------------------------------
   names = {
      "3050", "3060", "3080", "3083", "3090",
      "3093", "3094", "3095", "3096", "3097",
      "3098", "3099", "3102", "3106", "3112",
      "3120"
   };
   int N = names.size();

   // OLD:
   // ebeam = {
      // 3050.213, 3059.257, 3080.195, 3083.060, 3089.418,
      // 3092.324, 3095.261, 3095.994, 3096.390, 3097.777,
      // 3098.904, 3099.606, 3101.923, 3106.144, 3112.615,
      // 3120.442
   // };
   // for ( auto& e : ebeam ) {
      // e -=0.55;
   // }

   // arXiv:2206.13674v1
   ebeam = {
      3049.642, 3058.693, 3079.645, 3082.496, 3088.854,
      3091.760, 3094.697, 3095.430, 3095.826, 3097.213,
      3098.340, 3099.042, 3101.359, 3105.580, 3112.051,
      3119.878
   };
   if( ebeam.size() != N ) {
      cout << " ERROR in " << __func__ << " size of ebeam" << endl;
      exit(EXIT_FAILURE);
   }

   er_eb = {
      0.026, 0.028, 0.023, 0.023, 0.022,
      0.025, 0.084, 0.081, 0.075, 0.076,
      0.075, 0.093, 0.106, 0.090, 0.093,
      0.115
   };
   if( er_eb.size() != N ) {
      cout << " ERROR in " << __func__ << " size of er_eb" << endl;
      exit(EXIT_FAILURE);
   }
   // additional error due to calibration: common for all...
   // double er_shift = 0.033;
   // for ( auto& er : er_eb ) {
      // er = sqrt( SQ(er) + SQ(er_shift) );
   // }

   lumi = {
      14.919, 15.060, 17.393, 4.769, 15.558,
      14.910,  2.143,  1.816, 2.135,  2.069,
       2.203,  0.756,  1.612, 2.106,  1.720,
       1.264
   };
   if( lumi.size() != N ) {
      cout << " ERROR in " << __func__ << " size of lumi" << endl;
      exit(EXIT_FAILURE);
   }
}

// {{{1  tau scan: 2018
//--------------------------------------------------------------------
void setScan18( vector<string>& names, vector<double>& ebeam,
      vector<double>& er_eb, vector<double>& lumi ) {
//--------------------------------------------------------------------
   names = {
      "J1", "J2", "J3", "J4",
      "J5", "J6", "J7", "J8"
   };
   int N = names.size();

   // bes3memo-v1.1.pdf
   ebeam = {
      3087.593, 3095.726, 3096.203, 3096.986,
      3097.226, 3097.654, 3098.728, 3104.000
   };
   if( ebeam.size() != N ) {
      cout << " ERROR in " << __func__ << " size of ebeam" << endl;
      exit(EXIT_FAILURE);
   }

   er_eb = {
      0.125, 0.077, 0.069, 0.083, 0.099,
      0.082, 0.078, 0.082
   };
   if( er_eb.size() != N ) {
      cout << " ERROR in " << __func__ << " size of er_eb" << endl;
      exit(EXIT_FAILURE);
   }

   lumi = {
      2.467, 2.918, 4.978, 3.103,
      1.681, 4.656, 5.638, 5.716
   };
   if( lumi.size() != N ) {
      cout << " ERROR in " << __func__ << " size of lumi" << endl;
      exit(EXIT_FAILURE);
   }
}

// {{{1  R-scan 2015
//--------------------------------------------------------------------
void setRscan15( vector<string>& names, vector<double>& ebeam,
      vector<double>& er_eb, vector<double>& lumi ) {
//--------------------------------------------------------------------
   names = { "3080_rs" };
   int N = names.size();

   ebeam =  { 3080.  };
   if( ebeam.size() != N ) {
      cout << " ERROR in " << __func__ << " size of ebeam" << endl;
      exit(EXIT_FAILURE);
   }

   er_eb = vector<double>(N,1.0);

   lumi = { 126.21 };
   if( lumi.size() != N ) {
      cout << " ERROR in " << __func__ << " size of lumi" << endl;
      exit(EXIT_FAILURE);
   }
}

//--------------------------------------------------------------------
void setRscan15all( vector<string>& names, vector<double>& ebeam,
      vector<double>& er_eb, vector<double>& lumi ) {
//--------------------------------------------------------------------
   names = {
      "2900_rs", "2950_rs", "2981_rs", "3000_rs", "3020_rs", "3080_rs"
           };
   int N = names.size();

   ebeam =  {
       2900., 2950., 2981., 3000., 3020., 3080.
   };
   if( ebeam.size() != N ) {
      cout << " ERROR in " << __func__ << " size of ebeam" << endl;
      exit(EXIT_FAILURE);
   }

   er_eb = vector<double>(N,1.0);

   lumi = {
      105.53, 15.960, 16.046, 15.849, 17.315, 126.21
   };
   if( lumi.size() != N ) {
      cout << " ERROR in " << __func__ << " size of lumi" << endl;
      exit(EXIT_FAILURE);
   }
}

// {{{1 get_cross_section()
//--------------------------------------------------------------------
void get_cross_section( bool UseFitMkk,
      string pdfeff="eff", string pdfcs="cs", string prtCS="cs",
      bool TeX=false, string dir="" ) {
//--------------------------------------------------------------------
   // set names, energies and lumi
   vector<string> name12,  name18,  nameR;
   vector<double> ebeam12, ebeam18, ebeamR;
   vector<double> er_eb12, er_eb18, er_ebR;  // for tables
   vector<double> lumi12,  lumi18,  lumiR;

   setJpsiScan12(name12,ebeam12,er_eb12,lumi12);
   setScan18(name18,ebeam18,er_eb18,lumi18);
   setRscan15(nameR,ebeamR,er_ebR,lumiR);

   //-----------------------------------------------------------------
   // obtain efficiency
   vector<double> Eff12[2], Eff18[2], EffR[2];
   vector<TH1D*> HstMc12, HstMc18, HstMcR;

   getEfficiency( UseFitMkk, name12, Eff12, HstMc12 );
   getEfficiency( UseFitMkk, name18, Eff18, HstMc18 );
   getEfficiency( UseFitMkk, nameR,  EffR,  HstMcR );

   if ( !pdfeff.empty() ) {
      // draw efficiency
      // string eff12 = pdfeff + "_12.pdf";
      // string eff18 = pdfeff + "_18.pdf";
      // string effR  = pdfeff + "_R.pdf";
      // drawHstRes( eff12, ebeam12, HstMc12, true, Eff12 );
      // drawHstRes( eff18, ebeam18, HstMc18, true, Eff18 );
      // drawHstRes( effR,  ebeamR,  HstMcR,  true, EffR );

      string effall= pdfeff + "_all.pdf";
      vector<TGraphErrors*> geff {
         get_graph( ebeam12, true, Eff12 ),
         get_graph( ebeam18, true, Eff18 ),
         get_graph( ebeamR,  true, EffR )
      };
      drawGraph( effall, geff, true );
   }

   //-----------------------------------------------------------------
   // calculate cross section
   vector<double> Sig12[2], Nd12[2];
   vector<double> Sig18[2], Nd18[2];
   vector<double> SigR[2],  NdR[2];
   vector<TH1D*> HisD12, HisD18, HisDR;


   if ( !UseFitMkk ) { // OLD just side-band
      getCrossSection( name12,lumi12,Eff12[0], Sig12,Nd12,HisD12);
      getCrossSection( name18,lumi18,Eff18[0], Sig18,Nd18,HisD18);
      getCrossSection( nameR, lumiR, EffR[0],  SigR, NdR, HisDR);
      // if ( !pdfcs.empty() ) {
         // // draw cross-section separatly
         // string cs12 = pdfcs + "_12.pdf";
         // string cs18 = pdfcs + "_18.pdf";
         // string csR  = pdfcs + "_R.pdf";
         // drawHstRes( cs12, ebeam12, HisD12, false, Sig12 );
         // drawHstRes( cs18, ebeam18, HisD18, false, Sig18 );
         // drawHstRes( csR,  ebeamR,  HisDR,  false, SigR );
      // }
   } else { // NEW: read after mass_KK_fit()
      string mkkdir("mkk_inter/");
      if ( !dir.empty() ) {
         mkkdir += dir + "/";
      }
      read_mkk_file(name12,Nd12,mkkdir);
      read_mkk_file(name18,Nd18,mkkdir);
      read_mkk_file(nameR, NdR,mkkdir);
      CrossSection( Nd12,lumi12,Eff12[0], Sig12 );
      CrossSection( Nd18,lumi18,Eff18[0], Sig18 );
      CrossSection( NdR, lumiR, EffR[0],  SigR );
   }

   if ( !pdfcs.empty() ) { // draw cross-section
      string csall= pdfcs + "_all.pdf";
      vector<TGraphErrors*> gcs {
         get_graph( ebeam12, false, Sig12 ),
         get_graph( ebeam18, false, Sig18 ),
         get_graph( ebeamR,  false, SigR )
      };
      drawGraph( csall, gcs, false );
   }

   // print tables with cross-section
   string Tinfo("Side-Band subtruction only");
   if ( UseFitMkk ) {
      Tinfo="Fitting the Mkk distributions";
   }
   string t12("J/Psi scan 2012");
   string t18("J/Psi scan 2018");
   string tR("R-scan 2015");
   FILE* fp = stdout;  // by default
   if ( !prtCS.empty() ) {
      string ftxt = prtCS + ((TeX) ? ".tex" : ".txt"); // to file
      if ( !dir.empty() ) {
         ftxt = "mkk_inter/"+dir+"/" + ftxt;
      }
      fp = fopen(ftxt.c_str(),"w");
      if ( !fp ) {
         printf("Error open file %s\n",ftxt.c_str());
         fp = stdout;
      }
   } else {
      time_t temp = time(NULL);
      struct tm * timeptr = localtime(&temp);
      char Cdate[100];
      strftime(Cdate,sizeof(Cdate),"%B %d %Y",timeptr);
      printf("\nResults obtained on %s\n",Cdate);
   }
   fprintf(fp,"\n  %s\n", Tinfo.c_str());
   prtTable( fp, TeX, t12, ebeam12, lumi12, Eff12[0], Nd12, Sig12 );
   prtTable( fp, TeX, t18, ebeam18, lumi18, Eff18[0], Nd18, Sig18 );
   prtTable( fp, TeX, tR,  ebeamR,  lumiR,  EffR[0],  NdR,  SigR );
   if ( fp != stdout ) {
      fclose(fp);
   }

   if ( !prtCS.empty() ) {
      // print cpp-code with cross-section in file
      string fcpp = prtCS + "_results.h";
      if ( !dir.empty() ) {
         fcpp = "mkk_inter/"+dir+"/" + fcpp;
      }
      FILE* fp = fopen(fcpp.c_str(),"w");
      if ( !fp ) {
         printf("Error open file %s\n", fcpp.c_str());
         fp = stdout;
      }
      prtXcpp( fp, Tinfo );
      prtCScpp( fp, t12, "12", ebeam12, er_eb12, Sig12 );
      prtCScpp( fp, t18, "18", ebeam18, er_eb18, Sig18 );
      prtCScpp( fp, tR,  "R",  ebeamR,  er_ebR,   SigR  );
      if ( fp != stdout ) {
         fclose(fp);
      }
   }
}

// {{{1 MAIN:
//--------------------------------------------------------------------
void cross_section() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   gStyle->SetStatFont(62);
   gStyle->SetLegendFont(42);

   time_t temp = time(NULL);
   struct tm* timeptr = localtime(&temp);
   char buf[32];
   strftime(buf,sizeof(buf),"%d%b%y",timeptr);
   string dat(buf);

   bool UseFitMkk = true;
   // bool UseFitMkk = false;  // SB only

   if ( !UseFitMkk ) {
      dat="SB_"+dat;
   }

   get_cross_section(UseFitMkk,"eff_"+dat,"cs_"+dat,"cs_"+dat);

   // for TeX tables
   // bool TeX = true;
   // get_cross_section(UseFitMkk,"eff_"+dat,"cs_"+dat,"cs_"+dat,TeX);

   // === systematic:
   //     tables only, txt format, variation param in name,
   //     last argument is the directory in mkk_inter/

   // dat="Ch2_100";
   // dat="Ch2_60";

   // dat="Weta4";
   // dat="Weta2";

   // dat="AngP";
   // dat="AngM";

   // get_cross_section(UseFitMkk,"","","cs_"+dat,false,dat);
}
