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
   if ( var == "Lmin" ) {
//       Vmin = -82.5e3; Vmax = -72.5e3; // ch60
      Vmin = -76.e3; Vmax = -66.5e3;
      title = ";#it{-2log(L_{max})};Events/100";
   } else if (var == "pv09" || var == "pv12") {
      Vmin = 0.; Vmax = 1.;
      title = string(";#it{p-value(") +
         ((var == "pv09") ? "2009" : "2012") +
         ")};Events/0.01";
   } else if (var == "bkk") {
      select = var + string("*1e4>>") + hname;
      Nbins = 100; Vmin = 4.25; Vmax = 4.75;
      title = ";Br(J/#psi #rightarrow KK#eta) #times 10^{-4}";
      title += ";Events/5e-7";
   } else if (var == "bphi") {
      select = var + string("*1e4>>") + hname;
      Nbins = 100; Vmin = 7.5; Vmax = 9.5;
//       Nbins = 100; Vmin = 4.5; Vmax = 12.5;
      title = ";Br(J/#psi #rightarrow #phi#eta) #times 10^{-4}";
      title += ";Events/2e-6";
   } else if (var == "ang" ) {
      Vmin = -1.; Vmax = 1.;
      title = ";#vartheta, rad;Events/0.01";
   } else if (var == "sig09" || var == "sig12") {
      select = var + string("*1e3>>") + hname;
      title = ";#sigma(M^{ inv}_{ K^{+}K^{-}}), MeV/c^{2}";
      if (var == "sig09") {
         Nbins = 140; Vmin = 0.7; Vmax = 2.1;
      } else {
         Nbins = 80; Vmin = 0.7; Vmax = 1.5;
      }
      title += ";Events/0.01MeV/c^{2}";
   } else if (var == "nbg09" ) {
//       Vmin = 5.; Vmax = 30.;Nbins = 25; // ch60
      Vmin = 0.; Vmax = 20.;Nbins = 20;
      title = ";N_{bg}(2009);Events";
   } else if (var == "nbg12" ) {
//       Vmin = 30.; Vmax = 80.;Nbins = 50; //ch60
      Vmin = 15.; Vmax = 45.;Nbins = 30;
      title = ";N_{bg}(2012);Events";
   } else if (var == "ar09" ) {
//       Vmin = 4.5; Vmax = 11.5;Nbins = 70; // ch60
      Vmin = 2.5; Vmax = 9.5;Nbins = 70;
//       title = ";a_{side-band}(2009);Entries/0.1";
      title = ";a(side-band);Entries/0.1"; // common for 09 & 12
   } else if (var == "ar12" ) {
      Vmin = 1.5; Vmax = 8.5;Nbins = 70;
      title = ";a_{side-band}(2012);Entries/0.1";
   } else {
      cerr << " unknown var=" << var << endl;
      exit(0);
   }

   TH1D* hst = new TH1D(hname.c_str(),title.c_str(),Nbins,Vmin,Vmax);

   tmc -> Draw(select.c_str(),"st>0","goff");

   return hst;
}

//-------------------------------------------------------------------------
tuple<TH2D*,double> get_hist2D(string fname, string var) {
//-------------------------------------------------------------------------
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }
   TTree* tmc = (TTree*)gDirectory -> Get("tmc");

   string hname = string("toymc_2D") + var;
   string select;
   string cut("st>0");
   string cutTrue;
   string title;
   int Nbins = 100;
   double Vmin = 0, Vmax = 0;
   if ( var == "bkk" ) {
      title ="Toy MC: Br(J/#psi #rightarrow KK#eta) #times 10^{-4}";
      title+=";lower limit";
      title+=";upper limit";
      select = string("(bkk+ubkk)*1e4:(bkk-ubkk)*1e4>>") + hname;
      cut += string("&&ubkk*lbkk>0");
      cutTrue = cut + "&&(bkk+ubkk)*1e4>4.5&&(bkk-lbkk)*1e4<4.5";
      Vmin = 4.1; Vmax = 4.9;
   } else if (var == "bphi") {
      title ="Toy MC: Br(J/#psi #rightarrow #phi#eta) #times 10^{-4}";
      title+=";lower limit";
      title+=";upper limit";
      select = string("(bphi+ubphi)*1e4:(bphi-lbphi)*1e4>>") + hname;
      cut += string("&&ubphi*lbphi>0");
      cutTrue = cut + "&&(bphi+ubphi)*1e4>8.5&&(bphi-lbphi)*1e4<8.5";
      Vmin = 7.; Vmax = 10.;
   } else if (var == "ang" ) {
      title ="Toy MC: mixing angle #vartheta";
      title+=";lower limit";
      title+=";upper limit";
      select = string("(ang+uang):(ang-lang)>>") + hname;
      cut += string("&&uang*lang>0");
      cutTrue = cut + "&&(ang+uang)>0&&(ang-lang)<0";
      Vmin = -1.2; Vmax = 1.2;
   } else {
      cerr << " unknown var=" << var << endl;
      exit(0);
   }

   TH2D* h2d = new TH2D(hname.c_str(),title.c_str(),
         Nbins,Vmin,Vmax, Nbins,Vmin,Vmax);

   int Nall = tmc -> Draw(select.c_str(),cut.c_str(),"goff");
   int Ntrue = tmc -> Draw("st",cutTrue.c_str(),"goff");

   double ptrue = double(Ntrue) / double(Nall);
//    cout << " ptrue= " << ptrue << endl;

   return make_tuple(h2d,ptrue);
}

//-------------------------------------------------------------------------
TF1* GaussFit(TH1D* hist) {
//-------------------------------------------------------------------------
   const double ns = 3.; // number of sigmas to fit
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
   TH1D* hst = get_hist(fname, var);

   TPaveText* pt = nullptr;
   pt = new TPaveText(0.15,0.80,0.40,0.89,"NDC");
   pt -> SetTextAlign(22);
   pt -> SetTextFont(42);
   pt -> AddText(" Toy MC ");

   hst -> GetXaxis() -> SetTitleOffset(1.1);
   hst -> GetYaxis() -> SetTitleOffset(1.3);

   // must be here for GaussFit
   TCanvas* c1 = new TCanvas("c1","note",0,0,900,900);
   c1 -> cd();
   gPad -> SetGrid();

   if ( var == "Lmin" ) {
      hst -> GetXaxis() -> SetNdivisions(1005);
      pt -> SetX2(0.35);
   } else if ( var == "pv09" || var == "pv12" ) {
      pt -> SetX2(0.35);
   } else if ( var == "bkk" ) {
      GaussFit(hst);
   } else if ( var == "bphi" ) {
//       gStyle -> SetOptStat(1110);
//       gStyle->SetStatH(0.1);
      GaussFit(hst);
   } else if (var == "ang" ) {
      pt -> AddText("#vartheta = 0");
      hst -> GetXaxis() -> SetNdivisions(1005);
      hst -> GetYaxis() -> SetMaxDigits(3);
   } else if ( var == "sig09" || var == "sig12" ) {
      GaussFit(hst);
      if (var == "sig09") {
         pt -> AddText("#sigma(2009) = 1.4MeV/c^{2}");
      } else {
         pt -> AddText("#sigma(2012) = 1.1MeV/c^{2}");
      }
   } else if ( var == "nbg09" || var == "nbg12" ) {
      GaussFit(hst);
      pt -> SetX1(0.11);
      pt -> SetX2(0.36);
      if (var == "nbg09") {
//          pt -> AddText("N_{bg}(2009) = 17"); // ch 60
         pt -> AddText("N_{bg}(2009) = 10");
      } else {
//          pt -> AddText("N_{bg}(2012) = 54"); // ch 60
         pt -> AddText("N_{bg}(2012) = 30");
      }
   } else if ( var == "ar09" || var == "ar12" ) {
      GaussFit(hst);
      pt -> SetX1(0.11);
      pt -> SetX2(0.36);
      if (var == "ar09") {
//          pt -> AddText("a(sb2009) = 8."); // ch60
         pt -> AddText("a(sb) = 6."); // common for 09 & 12
      } else {
         pt -> AddText("a(sb2012) = 5.");
      }
   }

   SetHstFace(hst);
   hst -> SetLineWidth(2);

   hst -> Draw();

   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void plot_2D(string fname, string var, string pdf) {
//-------------------------------------------------------------------------
   TH2D* hst = nullptr;
   double ptrue = 0;
   tie(hst,ptrue) = get_hist2D(fname, var);

   TBox* box = new TBox;
   box -> SetFillStyle(3001);
   box -> SetFillColor(kRed-10);
   box -> SetLineColor(kRed-10);
   box -> SetLineWidth(2);

   TLatex Tl;
   Tl.SetTextAlign(12);
   Tl.SetTextSize(0.035);
   Tl.SetTextColor(kRed+3);

   hst -> GetXaxis() -> SetTitleOffset(1.2);
   hst -> GetYaxis() -> SetTitleOffset(1.3);

   TCanvas* c1 = new TCanvas("c1","2D",0,0,900,900);
   c1 -> cd();

//    SetHstFace(hst);
   gStyle->SetTitleFontSize(0.04);
   hst -> Draw();

   if ( var == "bkk" ) {
      box -> DrawBox(4.1,4.5,4.5,4.9);
      Tl.DrawLatex(4.14,4.85,"True value is included");
      Tl.DrawLatex(4.25,4.80,Form("%.0f%%",ptrue*100));
   } else if ( var == "bphi" ) {
      box -> DrawBox(7.0,8.5,8.5,10.);
      Tl.DrawLatex(7.12,9.8,"True value is included");
      Tl.DrawLatex(7.6,9.6,Form("%.0f%%",ptrue*100));
   } else if (var == "ang" ) {
      box -> DrawBox(-1.2,0.,0.,1.2);
      Tl.DrawLatex(-1.1,1.1,"True value is included");
      Tl.DrawLatex(-0.5,0.9,Form("%.0f%%",ptrue*100));
   }

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

//-------------------------------------------------------------------------
void picts() {
//-------------------------------------------------------------------------
   gROOT -> Reset();
   gStyle -> SetOptStat(0);
//    gStyle -> SetOptFit(0);

   //----------------------------------------------------------------
   // file-name:
   //----------------------------------------------------------------
   // ch60:
   // Br(KKeta)=4.5e-4; Br(phi eta)=8.5e-4; ang=0;
   // sig09=1.4e-3; Nbg09=17; arsb09=8.;
   // sig12=1.1e-3; Nbg12=54; arsb12=5.;
//    string fname("ToyMC_cf_4K.root");
//    string fname("ToyMC_cf_D4K.root");
   //----------------------------------------------------------------
   // ch40:
   // Br(KKeta)=4.5e-4; Br(phi eta)=8.5e-4; ang=0;
   // sig09=1.4e-3; Nbg09=10; arsb09=6.;
   // sig12=1.1e-3; Nbg12=30; arsb12=6.;
   string fname("ToyMC_cf_4K.root");
//    string fname("ToyMC_cfD_2K.root");
   //----------------------------------------------------------------

//    plot_var(fname,"Lmin","ToyMC_Lmin.pdf");
//    plot_var(fname,"pv09","ToyMC_pv09.pdf");
//    plot_var(fname,"pv12","ToyMC_pv12.pdf");

//    plot_var(fname,"bkk","ToyMC_bkk.pdf");
   plot_var(fname,"bphi","ToyMC_bphi.pdf");
//    plot_2D(fname,"bkk","ToyMC_2d_bkk.pdf");
//    plot_2D(fname,"bphi","ToyMC_2d_bphi.pdf");

//    plot_var(fname,"ang","ToyMC_ang.pdf");
//    plot_2D(fname,"ang","ToyMC_2d_ang.pdf");

//    plot_var(fname,"sig09","ToyMC_sig09.pdf");
//    plot_var(fname,"sig12","ToyMC_sig12.pdf");

//    plot_var(fname,"nbg09","ToyMC_nbg09.pdf");
//    plot_var(fname,"nbg12","ToyMC_nbg12.pdf");

//    plot_var(fname,"ar09","ToyMC_arsb09.pdf");
//    plot_var(fname,"ar12","ToyMC_arsb12.pdf");

   //----------------------------------------------------------------
}
