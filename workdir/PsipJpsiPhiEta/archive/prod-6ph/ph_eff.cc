// Study of the photon reconstruction efficiency:
// plot 1D distributions  ->

//----------------------------------------------------------------------
// GLOBAL:
// name of folder with root files
// static const string dir("local/");
const string dir("prod-6ph/");

// names of files to use
static vector<string> fnames;
void InitFnames() { // must be in main()
   fnames = {
      "data_09psip_ph.root", "mumu_09.root",
      "data_12psip_ph.root", "mumu_12.root"
//       "data_09psip_ph.root", "mcinc_09psip_ph.root",
//       "data_12psip_ph.root", "mcinc_12psip_ph.root"
   };
}

// use weighted histograms:
// 0 - no weights;
// 1 - calculate weights here; } plot only TRK*PID eff
// 2 - use "W" histograms.     }
// const int use_rew = 0;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

//----------------------------------------------------------------------
double ReWeightTrkPid(int DataPeriod, int Kp, double Pt) {
//----------------------------------------------------------------------
// This correction is based on "prod-te6"
// and independent of cos(Theta)
// Kp = 1 for kaons and Kp = 0 for pions

   double W = 1.;

   auto CUBE = [](double x)-> double{return x*x*x;};
   if ( Kp == 1 ) {             // kaons
      Pt = max( 0.1, Pt );
      Pt = min( 1.4, Pt );
      if ( DataPeriod == 2009 ) {
         W = 1.00703 - 0.01977 * Pt;
         if ( Pt < 0.2 ) {
            W = 0.9311;
         }
      } else if ( DataPeriod == 2012 ) {
         static TF1* cK12 = nullptr;
         if ( !cK12 ) {
            int nch = 4;
            auto Lchb = [nch](const double* xx, const double* p) -> double {
               if (nch == 0) { return p[0]; }
               double x = xx[0];
               double sum = p[0] + x*p[1];
               if (nch == 1) { return sum; }
               double T0 = 1, T1 = x;
               for ( int i = 2; i <= nch; ++i ) {
                  double Tmp = 2*x*T1 - T0;
                  sum += p[i]*Tmp;
                  T0 = T1;
                  T1 = Tmp;
               }
               return sum;
            };
            cK12 = new TF1("cK12", Lchb, 0.1, 1.4,nch+1);
            cK12->SetParameters(2.04378,-1.78748,1.05229,-0.40293,0.07065);
         }
         W = cK12->Eval(Pt);
      }
   } else if ( Kp == 0 ) {      // pions
      Pt = max( 0.05, Pt );
      Pt = min( 0.4, Pt );
      if ( DataPeriod == 2009 ) {
         W = 0.9871 + CUBE(0.0224/Pt);
//          W = 0.9870 + CUBE(0.0232/Pt); //second iteration
      } else if ( DataPeriod == 2012 ) {
         W = 0.9843 + CUBE(0.03015/Pt);
//          W = 0.9837 + CUBE(0.0315/Pt); //second iteration
      }
   }
   return W;
}

//----------------------------------------------------------------------
void ReWeightTrkPid(int data, int Kp, TH2D* hst) {
//----------------------------------------------------------------------

   int nx = hst->GetNbinsX();
   int ny = hst->GetNbinsY();
   for ( int iy = 0; iy <= ny+1; ++iy ) {
      double pt = hst->GetYaxis()->GetBinCenter(iy);
      double w = ReWeightTrkPid(data,Kp,pt);
      for ( int ix = 0; ix <= nx+1; ++ix ) {
         double d  = hst->GetBinContent(ix,iy);
         double ed = hst->GetBinError(ix,iy);

         hst->SetBinContent(ix,iy,d*w);
         hst->SetBinError(ix,iy,ed*w);
      }
   }
}

//----------------------------------------------------------------------
void PrintDeb(TH2D* hst) {  // DEBUG
//----------------------------------------------------------------------
   string ll(50,'-');
   cout << ll << endl
        << " Content of " << hst->GetName() << endl
        << ll << endl;
   int nx = hst->GetNbinsX();
   int ny = hst->GetNbinsY();
   for ( int ix = 0; ix <= nx+1; ++ix ) {
      double ct = hst->GetXaxis()->GetBinCenter(ix);
      for ( int iy = 0; iy <= ny+1; ++iy ) {
         double pt = hst->GetYaxis()->GetBinCenter(iy);
         double d  = hst->GetBinContent(ix,iy);
         double ed = hst->GetBinError(ix,iy);
         cout << "<ct=" << ct << "> <pt=" << pt << "> ==>"
              << d << " +/- " << ed << endl;
      }
   }
   cout << ll << endl;
}

//----------------------------------------------------------------------
void PrintDeb(TH1D* hst) {  // DEBUG
//----------------------------------------------------------------------
   string ll(50,'-');
   cout << ll << endl
        << " Content of " << hst->GetName() << endl
        << ll << endl;
   int nx = hst->GetNbinsX();
   for ( int ix = 0; ix <= nx+1; ++ix ) {
      double x = hst->GetXaxis()->GetBinCenter(ix);
      double d  = hst->GetBinContent(ix);
      double ed = hst->GetBinError(ix);
      cout << "<X=" << x << "> ==>" << d << " +/- " << ed << endl;
   }
   cout << ll << endl;
}

//-------------------------------------------------------------------------
void RmFromTitle(TH1* h, string rm) {
//-------------------------------------------------------------------------
   string title(h->GetTitle());
   string::size_type idx = title.find(rm);
   if ( idx == string::npos ) return;
   title = title.substr(0,idx) + title.substr(idx+rm.size());
   h->SetTitle(title.c_str());
}

//----------------------------------------------------------------------
TH1D* HstDiff(TH1D* A, TH1D* B) {
//----------------------------------------------------------------------
// calculate difference = (A-B) with sigma^2 = sigma^2(A)-sigma^2(B)

   TH1D* diff = static_cast<TH1D*>(A->Clone("diff"));
   int nx = A->GetNbinsX();
   for ( int i = 0; i <= nx+1; ++i ) {
      diff->SetBinContent( i,
            A->GetBinContent(i) - B->GetBinContent(i)
                         );
      diff->SetBinError( i,
            sqrt( fabs(SQ(A->GetBinError(i)) - SQ(B->GetBinError(i))) )
                       );
   }
   return diff;

}

//----------------------------------------------------------------------
void SetHstFace(TH1* hst) {
//----------------------------------------------------------------------
   TAxis* X = hst->GetXaxis();
   if ( X ) {
      X->SetLabelFont(62);
      X->SetLabelSize(0.04);
      X->SetTitleFont(62);
      X->SetTitleSize(0.04);
   }
   TAxis* Y = hst->GetXaxis();
   if ( Y ) {
      Y->SetLabelFont(62);
      Y->SetLabelSize(0.04);
      Y->SetTitleFont(62);
      Y->SetTitleSize(0.04);
   }
   TAxis* Z = hst->GetXaxis();
   if ( Z ) {
      Z->SetLabelFont(62);
      Z->SetLabelSize(0.04);
      Z->SetTitleFont(62);
      Z->SetTitleSize(0.04);
   }
}

//-------------------------------------------------------------------------
void get_hst( string fname, int type, vector<TH2D*>& hst ) {
//-------------------------------------------------------------------------
// get set of histograms from file "fname" according of "type"

   // names of histograms in MuMuGamma
   const vector<string> n_ph {  // photon efficiency
      "D_g0",  "D_g1"
   };

   bool isMC = (fname.find("data") == string::npos);

   vector<string> hisn;
   hisn = n_ph;
   if ( type == 1 ) {            // narrow region of impact
      for ( auto& hn : hisn ) {
         hn += string("S");
      }
   } else if ( type == 2 ) {     // Side-Band for impact parameters
      for ( auto& hn : hisn ) {
         hn += string("B");
      }
   } else if ( type != 0 ) {
      cout << " ERROR:: get_hst: type= " << type << endl;
      exit(1);
   }
   int Nh = hisn.size();

   hst.clear();
   hst.resize(Nh,nullptr);

   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }
   cout << " file: " << fname << endl;
   froot->cd("MuMuGamma");

   for ( int i = 0; i < Nh; ++i ) {
      hst[i] = (TH2D*)gROOT->FindObject(hisn[i].c_str());
      if ( !hst[i] ) {
         cout << " can not find histo:" << hisn[i] << endl;
         exit(1);
      }
      hst[i]->Sumw2(true);
   }

   // Names and Titles of histograms
   string str_dd("rec. #gamma"); // data description
   if ( !isMC ) {
      str_dd += " data";
   } else {
      str_dd += " MC";
   }
   if ( type == 1 ) {
      str_dd += "(narrow)";
   }
   if ( type == 2 ) {
      str_dd += "(SB)";
   }
   string SNY[2] = {" No "," Yes"};

   for ( int i = 0; i < Nh; i++) {
      string name = ((isMC) ? "M" : "D")
                  + string( hst[i]->GetName() ).substr(1);
//       cout << " name= " << name << endl; // DEB
      hst[i]->SetName( name.c_str() );

      string title = str_dd + SNY[i%2]
                   + "; cos(#Theta);"
                   + " E, GeV/c";
      hst[i]->SetTitle(title.c_str());
//       cout << " title= " << title << endl; // DEB

      SetHstFace(hst[i]);
   }

   // re-weighting PID hst
//    if ( isMC  && use_rew == 1 && type >= 10 ) {
//       int data = 2012;
//       if ( fname.find("_09") != string::npos ) {
//          data = 2009;
//       }
//       ReWeightTrkPid(data,1,hst[1]); // pid K+
//       ReWeightTrkPid(data,1,hst[3]); // pid K-
//       ReWeightTrkPid(data,0,hst[5]); // pid pi+
//       ReWeightTrkPid(data,0,hst[7]); // pid pi+
//       cout << " ===> re-weighting " << fname << endl;
//    }
}

//-------------------------------------------------------------------------
void get_eff_1D( string fname, int type, vector<TH1D*>& eff ) {
//-------------------------------------------------------------------------
// Get projections and calculate efficiencies

   vector<TH2D*> hst;
   get_hst(fname, type, hst);

   int Nh = hst.size();
   eff.clear();
   eff.resize(Nh,nullptr);

   for ( int i = 0; i < Nh/2; ++i ) {
      int j = 2*i; // index for hst

      int nx = hst[j]->GetNbinsX();
      int ny = hst[j]->GetNbinsY();

      string name = string( hst[j]->GetName() );
      string pxn = name + "_ct";
      string pyn = name + "_e";

      eff[j] = (TH1D*)hst[j]->ProjectionX(pxn.c_str(),1,ny);
      TH1D* tmp = (TH1D*)hst[j+1]->ProjectionX("tmp",1,ny);
      eff[j]->Add(tmp);
      eff[j]->Divide(tmp,eff[j],1.,1.,"B");

      eff[j+1] = (TH1D*)hst[j]->ProjectionY(pyn.c_str(),2,nx-1);
      tmp = (TH1D*)hst[j+1]->ProjectionY("tmp",2,nx-1);
      eff[j+1]->Add(tmp);
      eff[j+1]->Divide(tmp,eff[j+1],1.,1.,"B");
   }

   // Titles of histograms
   for ( int i = 0; i < Nh; i++) {
      string title( hst[i]->GetTitle() );
      title = "Eff. of " + title.substr(0,title.size()-3);
//       cout << " eff title= " << title << endl; // DEB
      eff[i]->SetTitle(title.c_str());
      eff[i]->GetYaxis()->SetTitle("efficiency");
      SetHstFace(eff[i]);
      eff[i]->GetXaxis()->SetTitleOffset(1.2);
      eff[i]->GetYaxis()->SetTitleOffset(1.2);
   }

}

//-------------------------------------------------------------------------
void get_ratio( const vector<TH1D*>& effdat, const vector<TH1D*>& effmc,
                  vector<TH1D*>& rat) {
//-------------------------------------------------------------------------
// calculate ratio of efficiencies DATA/MC

   int Nh = effdat.size();
   rat.clear();
   rat.resize(Nh,nullptr);
   for ( int i = 0; i < Nh; ++i ) {
      string name(effdat[i]->GetName());
      name = "r_" + name;
      rat[i]=(TH1D*)effdat[i]->Clone( name.c_str() );
      rat[i]->Divide(effmc[i]);
   }

   // Titles of histograms
   for ( int i = 0; i < Nh; ++i ) {
      string title( effmc[i]->GetTitle() );
      string::size_type idx = title.find("MC");
      title = title.substr(0,idx) + "data/" + title.substr(idx);
      rat[i]->SetTitle(title.c_str());
      rat[i]->GetYaxis()->SetTitle("data/MC");
      SetHstFace(rat[i]);
      rat[i]->GetXaxis()->SetTitleOffset(1.2);
      rat[i]->GetYaxis()->SetTitleOffset(1.2);
      rat[i]->SetMarkerStyle(20);
      rat[i]->SetLineColor(kBlack);
   }
}

//-------------------------------------------------------------------------
void get_diff( const vector<TH1D*>& hA, const vector<TH1D*>& hB,
               vector<TH1D*>& diff) {
//-------------------------------------------------------------------------
// return diff histograms calculated according HstDiff

   int Nh = hA.size();
   diff.clear();
   diff.resize(Nh,nullptr);

   for ( int i = 0; i < Nh; ++i ) {
      diff[i] = HstDiff(hA[i],hB[i]);
      string name( hA[i]->GetName() );
      name = "d_" + name;
      diff[i]->SetName( name.c_str() ); // rename
   }

   // Titles of histograms (#Delta in xpdf does not work)
   for ( int i = 0; i < Nh; i++) {
      string title( hA[i]->GetTitle() );
      title = title.substr(0,title.find("eff")) + "Delta";
      diff[i]->SetTitle(title.c_str());
      diff[i]->GetYaxis()->SetTitle("Delta");
      SetHstFace(diff[i]);
      diff[i]->GetXaxis()->SetTitleOffset(1.2);
      diff[i]->GetYaxis()->SetTitleOffset(1.2);
   }

//    DebPrint<TH1D*>(diff, "get_diff: ");
}

//-------------------------------------------------------------------------
void plot_pict(int date, string pdf) {
//-------------------------------------------------------------------------
   int idx = (date==2009) ? 0 : 2;

   vector<TH1D*> eff_d, eff_mc, rat0;
   get_eff_1D(fnames[idx], 0, eff_d);
   get_eff_1D(fnames[idx+1], 0, eff_mc);
   get_ratio( eff_d, eff_mc, rat0 );

//    vector<TH1D*> eff_dA, eff_mcA, ratA;
//    vector<TH1D*> diff_d, diff_mc, diff_rat;
//    get_eff_1D(fnames[idx], 1, eff_dA);
//    get_eff_1D(fnames[idx+1], 1, eff_mcA);
//    get_ratio( eff_dA, eff_mcA, ratA );
//    get_diff( eff_d, eff_dA, diff_d );
//    get_diff( eff_mc, eff_mcA, diff_mc );
//    get_diff( rat0, ratA, diff_rat );

   int Nh = eff_d.size();

   // Attributes of draw
   TLine* line0[3]  { new TLine(-0.80,0.,0.80,0.),
                      new TLine(0.1,0.,1.4,0.),
                      new TLine(0.05,0.,0.4,0.) };
   TLine* line1[3]  { new TLine(-0.80,1.,0.80,1.),
                      new TLine(0.1,1.,1.4,1.),
                      new TLine(0.05,1.,0.4,1.) };
   for ( int i = 0; i < 3; ++i ) {
      line0[i]->SetLineColor(kGreen+2);
      line0[i]->SetLineWidth(2);
      line0[i]->SetLineStyle(kSolid); // kDashed
      line1[i]->SetLineColor(kRed+1);
      line1[i]->SetLineWidth(2);
      line1[i]->SetLineStyle(kSolid); // kDashed
   }

   auto setC_X = [](TH1D* h) { h->SetAxisRange(-0.799,0.799,"X"); };
   auto setE_X = [](TH1D* h) { h->SetAxisRange(0.,1.899,"X"); };

   auto setBlue = [](TH1* h) {
      h->SetLineColor(kBlue+1);
      h->SetMarkerColor(kBlue+2);
      h->SetMarkerStyle(22);
      h->SetMarkerSize(0.9);
   };
   auto setRed = [](TH1* h) {
      h->SetLineColor(kRed+1);
      h->SetMarkerColor(kRed+2);
      h->SetMarkerStyle(23);
      h->SetMarkerSize(0.9);
   };

   for (int i = 0; i < Nh; ++i ) {
      RmFromTitle(eff_d[i],"data ");

      setBlue(eff_d[i]);
      setRed(eff_mc[i]);

      if ( i%2 == 0 ) {
         setC_X(eff_d[i]);
         setC_X(rat0[i]);
      } else {
         setE_X(eff_d[i]);
         setE_X(rat0[i]);
      }
   }

   TLegend* leg1 = new TLegend(0.40,0.20,0.72,0.40);
   leg1->AddEntry(eff_d[0], Form("#color[%i]{data %i}",
                  eff_d[0]->GetLineColor(),date),"LP");
   leg1->AddEntry(eff_mc[0], Form("#color[%i]{MC %i}",
                  eff_mc[0]->GetLineColor(),date),"LP");


   // for diff
   TF1* pl0 = (TF1*)gROOT->GetFunction("pol0")->Clone();
   pl0->SetLineColor(kGreen+2);

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->Divide(2,2);
   for (int i = 1; i <= 4; ++i ) {
      c1->cd(i);
      gPad->SetGrid();
   }

   bool ispdf = pdf.size() > 0;
   if ( ispdf ) {
      c1->Print((pdf+"[").c_str()); // open pdf-file
   }
   auto addPDF = [c1,ispdf,pdf](TVirtualPad* pad) {
      c1->Update();
      if ( ispdf ) {
         c1->Print(pdf.c_str()); // add to pdf-file
      } else {
         pad->WaitPrimitive(); // pause
      }
   };

///////////////////////////////////////////////////////////////////
   double eff_min = 0, eff_max = 1.;
//    double rat_min = 0.9, rat_max = 1.1;
   double rat_min = 0., rat_max = 2.;
///////////////////////////////////////////////////////////////////

   // Gamma efficiencies
   for (int i = 0; i < Nh; ++i ) {
      c1->cd(i+1);
      eff_d[i]->SetAxisRange(eff_min,eff_max,"Y");
      eff_d[i]->Draw("E");
      eff_mc[i]->Draw("E SAME");
      leg1->Draw();
   }
//    addPDF(gPad);

   // ratio
   for (int i = 0; i < Nh; ++i ) {
      c1->cd(i+3);
      rat0[i]->SetAxisRange(rat_min,rat_max,"Y");
      rat0[i]->Draw("E");
//       line1[i%2]->Draw();
   }
   addPDF(gPad);

///////////////////////////////////////////////////////////////////


   if ( ispdf ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
   }
}

//-------------------------------------------------------------------------
void ph_eff() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(62);
   InitFnames();

   // TESTS ---------------------------------------------------------
//    vector<TH2D*> hst;
//    get_hst(fnames[0], 1, hst);

//    vector<TH1D*> eff;
//    get_eff_1D(fnames[1], 0, eff);

//    int idx = 0;
//    int type = 0;
//    vector<TH1D*> effdat, effmc;
//    get_eff_1D(fnames[idx], type, effdat);
//    get_eff_1D(fnames[idx+1], type, effmc);
//    vector<TH1D*> rat; // data/MC
//    get_ratio( effdat, effmc, rat );
   // TESTS ---------------------------------------------------------

   int date = 2009;
//    int date = 2012;

//    string pdf = ""; // debug
   string pdf = string("Pheff_")
              + to_string(date) + string(".pdf");
//               + ( (use_rew==2) ? "W_" : ((use_rew==1) ? "_TW_" : "_") )
   plot_pict(date,pdf);
}
