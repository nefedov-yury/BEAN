// plot pictures for ToyMC

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
TH1D* get_hist(string fname, string var) {
//-------------------------------------------------------------------------
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }
   TTree* tmc = (TTree*)gDirectory -> Get("tmc");

   string hname = string("toymc_") + var;
   string title;
   int Nbins = 100;
   double Vmin = 0, Vmax = 0;
   string select = var + string(">>") + hname;
   if ( var == "lmin" ) {
      Vmin = -10300; Vmax = -9800;
//       Vmin = -10350; Vmax = -9850;
      title = ";#it{-2log(L_{max})};Events/5";
   } else if (var == "pv" ) {
      Vmin = 0.; Vmax = 1.;
      title = ";#it{p-value(K-S)};Events/0.01";
   } else if (var == "sig" ) {
      Nbins = 80; Vmin = 0.8; Vmax = 1.6;
      title = ";#sigma(M^{ inv}_{ K^{+}K^{-}}), MeV/c^{2}"
              ";Events/0.01MeV/c^{2}";
      select = var + string("*1e3>>") + hname;
   } else if (var == "fb" ) {
      Nbins = 80; Vmin = 0.6; Vmax = 1.4;
      title = ";F ;Events/0.01";
   } else if (var == "ang" ) {
      Vmin = -1.3; Vmax = 1.2;
      title = ";#vartheta, rad;Events/0.025";
   } else if (var == "nphi" ) {
      Vmin = 2200; Vmax = 2900;
      title = ";N_{#phi};Events/7";
   } else {
      cerr << " unknown var=" << var << endl;
      exit(0);
   }

   TH1D* hst = new TH1D(hname.c_str(),title.c_str(),Nbins,Vmin,Vmax);

   tmc -> Draw(select.c_str(),"","goff");

   return hst;
}

//-------------------------------------------------------------------------
TF1* GaussFit(TH1D* hist) {
//-------------------------------------------------------------------------
   const double ns = 2.; // number of sigmas to fit
   TF1* gs = (TF1*)gROOT -> GetFunction("gaus");
   gs -> SetLineWidth(2);
   gs -> SetLineColor(kRed);
   hist -> Fit(gs,"Q0");
   double gsmean = gs -> GetParameter(1);
   double gssig  = gs -> GetParameter(2);
   hist -> Fit(gs,"Q0","",gsmean-ns*gssig,gsmean+ns*gssig);
   gsmean = gs -> GetParameter(1);
   gssig  = gs -> GetParameter(2);
   hist -> Fit(gs,"","",gsmean-ns*gssig,gsmean+ns*gssig);
//    hist -> DrawCopy();
   return gs;
}

//-------------------------------------------------------------------------
void plot_var(string fname, string var, string pdf) {
//-------------------------------------------------------------------------
   bool isPosAng = (fname.find("_p_") != string::npos);
   TH1D* hst = get_hist(fname, var);

   TPaveText* pt = nullptr;
   if ( var == "ang" || var == "nphi" ) {
      pt = new TPaveText(0.40,0.80,0.60,0.89,"NDC");
   } else {
      pt = new TPaveText(0.15,0.80,0.35,0.89,"NDC");
   }
   pt -> SetTextAlign(22);
   pt -> SetTextFont(42);
   pt -> AddText(" Toy MC ");

   if ( var == "lmin" ) {
      hst -> GetXaxis() -> SetNdivisions(1005);
   } else if ( var == "sig" ) {
      GaussFit(hst);
      hst -> GetXaxis() -> SetTitleOffset(1.1);
      pt -> AddText("#sigma = 1.2 MeV/c^{2}");
   } else if ( var == "fb" ) {
      GaussFit(hst);
      pt -> AddText("F = 1.0");
   } else if ( var == "ang" ) {
      if ( isPosAng ) {
         pt -> AddText("#vartheta = 0.8");
      } else {
         pt -> AddText("#vartheta =-0.8");
      }
   } else if ( var == "nphi" ) {
      if ( isPosAng ) {
         pt -> AddText("N_{#phi} = 2781");
      } else {
         pt -> AddText("N_{#phi} = 2412");
      }
   }

   TCanvas* c1 = new TCanvas("c1","note",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   SetHstFace(hst);
   hst -> SetLineWidth(2);
   hst -> GetYaxis() -> SetTitleOffset(1.3);
//    hst -> SetLineColor(kBlack);
//    hst -> SetMarkerStyle(20);

//    hst -> Draw("EP");
   hst -> Draw();

   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      pdf += (isPosAng) ? "_p.pdf" : "_n.pdf";
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void plot_pvdif(string fname, string pdf) {
//-------------------------------------------------------------------------
   bool isPosAng = (fname.find("_p_") != string::npos);

   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }
   TTree* tmc = (TTree*)gDirectory -> Get("tmc");

   int n1 = tmc -> Draw("pv","np==-1","goff");
   double* pv1 = tmc -> GetVal(0);
   vector<double> Vpm(pv1, pv1+n1);

   int n2 = tmc -> Draw("pv","np==+1","goff");
   double* pv2 = tmc -> GetVal(0);
   vector<double> Vpp(pv2, pv2+n2);

   if ( n1 != n2 ) {
      cout << " ERROR: n1= " << n1 << " n2= " << n2 << endl;
      exit(0);
   }

   TCanvas* c1 = new TCanvas("c1","note",0,0,900,900);
   c1 -> cd();

/*
   TGraph* gr = new TGraph(n1, Vpm.data(), Vpp.data());
   gr -> SetTitle(";#it{p-value for negative interference}"
                  ";#it{p-value for positive interference}");
   gr -> Draw("AP");

   TLine* l = new TLine(0.,0.,1.1,1.1);
   l -> SetLineColor(kRed);
   l -> SetLineStyle(kDashed); // kSolid
   l -> Draw();
*/

   TH1D* hst = new TH1D("hst","",100,-0.05,0.05);
   hst -> SetTitle(";#it{p-value(-) - p-value(+)}");

   for ( int i = 0; i < n1; ++i ) {
      hst -> Fill (Vpm[i]-Vpp[i]);
   }
   gPad -> SetGrid();
   hst -> Draw();


   TPaveText* pt = new TPaveText(0.15,0.80,0.35,0.89,"NDC");
   pt -> SetTextAlign(22);
   pt -> SetTextFont(42);
   pt -> AddText(" Toy MC ");
   if ( isPosAng ) {
      pt -> AddText("#vartheta = 0.8");
   } else {
      pt -> AddText("#vartheta =-0.8");
   }
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      pdf += (isPosAng) ? "_p.pdf" : "_n.pdf";
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void picts() {
//-------------------------------------------------------------------------
   gROOT -> Reset();
   gStyle -> SetOptStat(0);
//    gStyle -> SetOptFit(0);
//    gStyle -> SetLegendFont(42);
//    gStyle -> SetLegendTextSize(0.03);

   //----------------------------------------------------------------
   // file-name:
//    string fname("mctoy_p_5000.root"); // positive theta = +0.8; 5000ev
//    string fname("mctoy_n_5000.root"); // negative theta = -0.8; 5000ev
//    string fname("mctoy_p_10K.root"); // positive theta = +0.8; 10 000ev
   string fname("mctoy_n_10K.root"); // negative theta = -0.8; 10 000ev

//    plot_var(fname,"lmin","Lmin_toyMC");
//    plot_var(fname,"pv","Pv_toyMC");
//
//    plot_var(fname,"sig","Sig_toyMC");
//    plot_var(fname,"fb","F_toyMC");
//
//    plot_var(fname,"ang","Ang_toyMC");
   plot_var(fname,"nphi","Nphi_toyMC");
//

//    plot_pvdif(fname,"Pvdiff_toyMC.pdf");
   //----------------------------------------------------------------
}
