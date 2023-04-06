// study of systematic uncertainties in efficiency due to the
// variation of parameters A and phi within their errors
// in MCGPJ generator
//    -> sys_mcgpj_DATE.pdf

#include <time.h>   // see man strftime

// {{{1 helper functions
//--------------------------------------------------------------------
constexpr double SQ(double x) {
//--------------------------------------------------------------------
   return x*x;
}

// {{{1 Get histograms
// 1) "nominal" histograms
//--------------------------------------------------------------------
tuple<TH1D*,TH1D*> getIniFin(string name) {
//--------------------------------------------------------------------
#include "cuts.h"

   string fname = "Ntpls/ntpl_mcgpj_" + name + ".root";
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cout << "can not open " << fname << endl;
      exit(0);
   }

   froot -> cd("SelectKKgg");
   TH1D* Xisr_ini = (TH1D*)gROOT -> FindObject("mcisr_1x");
   if( !Xisr_ini ) {
      cout << " file: " << fname << endl;
      cout << " can not find mcisr_1x histo" << endl;
      exit(0);
   }

   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   // MC cuts
   TCut c_here = c_xisr + c_MCmkk; // X_isr>0.9 && mc_Mkk<1.08
   // selection cuts
   c_here += c_chi2;
   c_here += c_cpgg;
   c_here += c_phi;     // [2*Mk, 1.08GeV]

   TH1D* Xisr_fin = new TH1D("Xisr_fin", ";1-s'/s", 100,0.,1.);
   Xisr_fin -> Sumw2();
   a4c -> Draw("(1.-xisr)>>Xisr_fin",c_here,"goff");

   return make_tuple(Xisr_ini, Xisr_fin);
}

// 2) Tested parameters histogram
//--------------------------------------------------------------------
TH1D* getHst(string dir, string name) {
//--------------------------------------------------------------------

   if ( name[4] == '_' ) {
      name = name.substr(0,4) + "R"; // replace '_rs' -> 'R'
   }
   string fname = "Ntpls/" + dir + "/phieta_" + name + "mc.root";
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cout << "can not open " << fname << endl;
      exit(0);
   }
   froot -> cd();
   TH1D* X_ini = (TH1D*)gROOT -> FindObject("h_x");

   return X_ini;
}

// {{{1 Calculate efficiency ratio
//--------------------------------------------------------------------
tuple<double,double> getReff( string dir, string name,
      vector<TH1D*>& Xhst ) {
//--------------------------------------------------------------------
   const double Xisr_min = 0.9;
   const double X1_max = 1 - Xisr_min;

   TH1D* Xini = nullptr;
   TH1D* Xfin = nullptr;
   tie(Xini,Xfin) = getIniFin(name);

   TH1D* Xh_ini = getHst(dir,name);
   TH1D* Xh_div = (TH1D*)Xh_ini->Clone("Xh_div");
   Xh_div -> Divide(Xh_ini,Xini);

   Xhst = {Xini, Xfin, Xh_ini, Xh_div}; // return vector

   int nb = Xini -> GetNbinsX();
   // double Nini = Xini -> Integral();
   // int ib_max = nb+1;

   // Nini in [0,0.1] interval
   int ib_max = int(rint(X1_max*nb))+1;
   double Xlow = Xini -> GetBinLowEdge(ib_max);
   if ( fabs(Xlow-X1_max) > 1e-6 ) {
      cout << " name: " << name << endl;
      cout << " wrong bin for X1_max: ib_max = " << ib_max
         << " LowEdge= " << Xlow << " X1_max= " << X1_max
         << endl;
      exit(0);
   }
   double Nini = Xini -> Integral(0,ib_max-1);
   // double Nini_h = Xh_ini -> Integral(0,ib_max-1);
   // printf("ib_max= %d, Xlow= %.2f, Nini= %f, Nini_h= %f\n",
         // ib_max,Xlow,Nini,Nini_h);

   double eff0 = 0, err0 = 0;
   double eff1 = 0, err1 = 0;
   for( int i = 0; i < ib_max; ++i ) {
      double a  = Xfin -> GetBinContent(i);
      double b  = Xh_div -> GetBinContent(i);
      if( a == 0 || b == 0 ) {
         continue;
      }

      eff0 += a;
      eff1 += a*b;

      double ea = Xfin -> GetBinError(i);
      double eb = Xh_div -> GetBinError(i);
      err0 += SQ(ea);
      err1 += SQ(a*eb) + SQ(b*ea);
   }

   eff0 /= Nini;
   eff1 /= Nini;
   err0  = sqrt(err0)/Nini;
   err1  = sqrt(err1)/Nini;

   double ratio = eff1/eff0;
   double err_r = ratio * sqrt( SQ(err0/eff0) + SQ(err1/eff1) );

   // printf(" eff0= %.2f +/- %.2f %%", eff0*100,err0*100);
   // printf(" eff1= %.2f +/- %.2f %%\n", eff1*100,err1*100);
   // printf(" 1-eff1/eff0= %.2f +/- %.2f %%\n",
         // (1-ratio)*100, err_r*100);

   return make_tuple(ratio, err_r);
}

// {{{1 Names and Energies
//--------------------------------------------------------------------
vector<size_t> getNamesEng(vector<string>& Names,
      vector<double>& Energy ) {
//--------------------------------------------------------------------

   vector<string> Names12 {
      "3050", "3060", "3080", "3083", "3090",
      "3093", "3094", "3095", "3096", "3097",
      "3098", "3099", "3102", "3106", "3112",
      "3120"
   };
   // Note: 0.55MeV is already subtracted
   vector<double> Energy12 {
      3.049663, 3.058707, 3.079645, 3.082510, 3.088868,
      3.091774, 3.094711, 3.095444, 3.095840, 3.097227,
      3.098354, 3.099056, 3.101373, 3.105594, 3.112065,
      3.119892
   };

   vector<string> Names18 {
      "J1", "J2", "J3", "J4",
      "J5", "J6", "J7", "J8"
   };
   vector<double> Energy18 {
      3.087659, 3.095726, 3.096203, 3.096986,
      3.097226, 3.097654, 3.098728, 3.104000
   };

   vector<string> NamesRs {
      "3080_rs","3020_rs","3000_rs",
      "2981_rs","2950_rs","2900_rs"
   };
   vector<double> EnergyRs {
      3.080, 3.020,  3.000,
      2.981,  2.950, 2.900
   };

   Names.clear();
   Names.insert(end(Names), begin(Names12), end(Names12));
   Names.insert(end(Names), begin(Names18), end(Names18));
   Names.push_back(NamesRs[0]);

   Energy.clear();
   Energy.insert(end(Energy), begin(Energy12), end(Energy12));
   Energy.insert(end(Energy), begin(Energy18), end(Energy18));
   Energy.push_back(EnergyRs[0]);

   vector<size_t> ret { Names12.size(), Names18.size(), 1 };
   return ret;
}

// {{{1 Graph for one dir
//--------------------------------------------------------------------
TGraphErrors* getGraph(string dir, string pdf) {
//--------------------------------------------------------------------

   vector<string> Names;
   vector<double> Energy;
   getNamesEng(Names,Energy);
   int Nd = Names.size();
   vector<double> ErrE(Nd,0.01);

   TCanvas* c1 = nullptr;
   if ( !pdf.empty() ) {
      c1 = new TCanvas("c1","...",0,0,800,800);
      c1 -> Divide(1,2);
      c1 -> cd(1);
      gPad -> SetGrid();
      gPad -> SetLogy();
      c1 -> cd(2);
      gPad -> SetGrid();
      c1 -> Print((pdf+"[").c_str()); // just open pdf-file
   }

   vector<TH1D*> Xhst;
   vector<double> vReff, vErrR;
   for ( int i = 0; i < Nd; ++i ) {
      double Reff = 0., errR = 0.;
      string name = Names[i];
      tie(Reff,errR) = getReff(dir, name, Xhst);
      vReff.push_back(100*(1-Reff));
      vErrR.push_back(100*errR);

      if ( c1 ) {
         c1 -> cd(1);
         auto& Xini = Xhst[0];
         auto& Xh_ini = Xhst[2];
         Xini -> SetTitle(";1 - s'/s");
         // Xini -> SetAxisRange(1.,2.e5,"Y");
         Xini -> SetLineColor(kRed-7);
         Xini -> SetFillColor(kRed-7); // for "E" options
         Xini -> DrawCopy("E3");
         Xh_ini -> SetLineColor(kBlue-7);
         Xh_ini -> SetFillColor(kBlue-7); // for "E" options
         Xh_ini -> DrawCopy("E3,SAME");
         TLegend* leg = new TLegend(0.60,0.70,0.89,0.88);
         leg -> SetHeader( Form("%s",name.c_str()),"C");
         leg -> AddEntry( Xini,"Nominal","L");
         leg -> AddEntry( Xh_ini,"MCGPJ","L");
         leg -> Draw();

         c1 -> cd(2);
         auto& Xh_div = Xhst[3];
         Xh_div -> SetTitle(";1 - s'/s");
         Xh_div -> SetAxisRange(0.,0.19,"X");
         // Xh_div -> SetAxisRange(0.5,1.5,"Y");
         Xh_div -> SetAxisRange(0.,2,"Y");
         Xh_div -> SetLineColor(kBlue+1);
         Xh_div -> DrawCopy("E1");
         TPaveText* pt = new TPaveText(0.60,0.79,0.89,0.88,"NDC");
         pt -> SetTextAlign(12);
         pt -> SetTextFont(42);
         pt -> AddText(Form("Reff = %.2f #pm %.2f %%",
                  Reff*100, errR*100));
         pt -> Draw();

         c1 -> Update();
         c1 -> Print(pdf.c_str()); // add to pdf-file

         delete leg;
         delete pt;
      }
   }

   int Ng = vReff.size();
   TGraphErrors* gr = new TGraphErrors(Ng,
         Energy.data(), vReff.data(), ErrE.data(), vErrR.data());
   string title = string(";center of mass energy (GeV)"
         ";relative difference (%)");
   gr->SetTitle(title.c_str());
   gr->GetXaxis()->CenterTitle();
   gr->GetYaxis()->CenterTitle();
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   gr->GetXaxis()->SetLabelSize(0.03);
   gr->GetYaxis()->SetLabelSize(0.03);
   // gr->GetXaxis()->SetTitleOffset(1.1);
   gr->GetYaxis()->SetTitleOffset(1.1);

   gr->SetMinimum(-5.);
   gr->SetMaximum(+5.);

   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kMagenta+1);

   if ( c1 ) {
      c1 -> cd();
      c1 -> Clear();
      gPad -> SetGrid();
      gr -> Draw("AP");
      c1 -> Update();
      c1 -> Print(pdf.c_str()); // add to pdf-file

      c1->Print((pdf+"]").c_str()); // just close pdf-file
      delete c1;
   }

   return gr;
}

// {{{1 Plot all graphs and print relative errors
//--------------------------------------------------------------------
void prtSysRel(FILE* fp, string title, const vector<double>& Rel) {
//--------------------------------------------------------------------
   const char* shift = "   "; // shift = 3 spaces
   const char* comma = ", ";
   vector<string> Names;
   vector<double> Energy;
   vector<size_t> Ndate = getNamesEng(Names,Energy);
   if ( size(Ndate) != 3 ) {
      cerr << " FATAL ERROR: prtSysRel: size(Ndate)= "
         << size(Ndate) << endl;
      exit(1);
   }

   fprintf(fp,"\n");
   fprintf(fp,"%s// %s\n",shift,title.c_str());
   string sufs[] {"12","18","R"};
   for ( int k = 0, l = 0; k < 3; k++) {
      int n =  Ndate[k];
      const string& suf = sufs[k];
      fprintf(fp,"%svector<double> RelSysMCGPJ%s = {\n",
            shift,suf.c_str());
      for(int i = 0; i < n; i += 5) {
         fprintf(fp,"%s%s",shift,shift);
         for(int j = i; j < (i+5<n ? i+5 : n); j++) {
            fprintf(fp,"%.2e%s",Rel[l+j], (j!=n-1) ? comma : "" );
         }
         fprintf(fp,"\n");
      }
      l += n;
      fprintf(fp,"%s};\n",shift);
   }
}

//--------------------------------------------------------------------
void PlotAllGraph(string pdf, string header) {
//--------------------------------------------------------------------

   // phi+sigma(phi) A+sigma(A), phi-sigma(phi), A-sigma(A)
   string DIR="MCGPJ_01/";
   vector<string> dirs { "un_pp", "un_ap", "un_pm", "un_am" };
   int mrks[]          {  22,      26,       23,      32    };
   int clrs[]          { kRed+1,  kGreen+2, kRed+1, kGreen+2};
   // int clrs[]          { kGreen+2, kRed+1, kYellow+3,kOrange-3};

   TGraphErrors* gr[4];
   for ( int i = 0; i < 4; ++i ) {
      gr[i] = getGraph(DIR+dirs[i],"");
      gr[i] -> SetMarkerStyle(mrks[i]);
      gr[i] -> SetMarkerColor(clrs[i]);
      gr[i] -> SetLineColor(clrs[i]);
   }

   TCanvas* c = new TCanvas("c","...",0,0,900,900);
   c -> cd();
   gPad -> SetGrid();

   TLegend* leg = new TLegend(0.14,0.67,0.89,0.89);
   leg -> SetNColumns(2);
   leg -> SetHeader("varying MCGPJ parameters","C");
   leg -> AddEntry( gr[0], "#varphi #plus #sigma(#varphi)", "PE" );
   leg -> AddEntry( gr[1], "A #plus #sigma(A)", "PE" );
   leg -> AddEntry( gr[2], "#varphi #minus #sigma(#varphi)", "PE" );
   leg -> AddEntry( gr[3], "A #minus #sigma(A)", "PE" );

   gr[0] -> SetMaximum(8);
   gr[0] -> Draw("AP");
   gr[1] -> Draw("P");
   gr[2] -> Draw("P");
   gr[3] -> Draw("P");
   leg -> Draw();

   gPad -> RedrawAxis();
   c -> Update();
   c -> Print(pdf.c_str());

   // print to file or to stdout
   FILE* fp = stdout;
   if ( !header.empty() ) {
      fp = fopen(header.c_str(),"w");
      if ( !fp ) {
         cerr << "Error open file " << header << endl;
         fp = stdout;
      }
   }

   int Nd = gr[0] -> GetN();
   vector<double> maxR(Nd,0.);
   for ( int i = 0; i < 4; ++i ) {
      const double* R = gr[i] -> GetY();
      for ( int j = 0; j < Nd; ++j ) {
         auto r = fabs(R[j] * 1e-2 ); // % in gr[]
         if ( r > maxR[j] ) {
            maxR[j] = r;
         }
      }
   }
   prtSysRel(fp,"varying MCGPJ parameters",maxR);
}

// {{{1~ Main
//--------------------------------------------------------------------
void sys_mcgpj() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   // gStyle->SetLegendTextSize(0.05);

   // plot each variation separatly
   string DIR="MCGPJ_01/";
   // getGraph(DIR+"un_0", "sys_mcgpj_un_0.pdf"); // no variation
   // getGraph(DIR+"un_am", "sys_mcgpj_un_am.pdf");
   // getGraph(DIR+"un_ap", "sys_mcgpj_un_ap.pdf");
   // getGraph(DIR+"un_pm", "sys_mcgpj_un_pm.pdf");
   // getGraph(DIR+"un_pp", "sys_mcgpj_un_pp.pdf");

   time_t temp = time(NULL);
   struct tm * timeptr = localtime(&temp);
   char buf[32];
   strftime(buf,sizeof(buf),"%d%b%y",timeptr);
   string dat(buf);

   string pdf = "sys_mcgpj_" + dat + ".pdf";
   string header = "sys_mcgpj_" + dat + ".h";
   PlotAllGraph(pdf,header);
}
