// Study of the track reconstruction efficiency for pi and K:
// fit the ratio data/MC as function of cos(Theta) and P_t
// -> trk_fit_[date].pdf
// pdftk trk_fitPM_2012.pdf cat 5 output trfit_12_K.pdf

//----------------------------------------------------------------------
// GLOBAL:
// name of folder with root files
static string dir;
// names of files to use
static vector<string> fnames;

const int date = 2009;

// use weighted histograms:
// 0 - no weights;
// 1 - calculate weights here; } plot only TRK*PID eff
// 2 - use "W" histograms.     }
const int use_rew = 1;

// sum "+" and "-":
bool sumPM = true;
// bool sumPM = false;

void InitFnames() { // must be in main()
   dir = string("prod-12eff/");
   fnames = {
      "data_09psip_all.root", "mcinc_09psip_all.root",
      "data_12psip_all.root", "mcinc_12psip_all.root"
   };
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

//----------------------------------------------------------------------
double ReWeightTrkPid(int DataPeriod, int Kp, double Pt) {
//----------------------------------------------------------------------
// This correction is based on "prod-12eff"
// and independent of cos(Theta)
// Kp = 1 for kaons and Kp = 0 for pions

   double W = 1.;

   auto CUBE = [](double x)-> double{return x*x*x;};
   if ( Kp == 1 ) {             // kaons
      Pt = max( 0.1, Pt );
      Pt = min( 1.4, Pt );
      if ( DataPeriod == 2009 ) {
         W = 1.00931 - 0.02363 * Pt;
         // OLD: prod-te6; -te10
//          W = 1.00703 - 0.01977 * Pt;
         if ( Pt < 0.2 ) {
            W = 0.9278;
//             W = 0.9311; // OLD: prod-te6; -te10
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
            cK12->SetParameters(1.82144,-1.41435,0.83606,-0.32437,0.05736);
            // OLD: prod-te6; -te10
//             cK12->SetParameters(2.04378,-1.78748,1.05229,-0.40293,0.07065);
         }
         W = cK12->Eval(Pt);
      }
   } else if ( Kp == 0 ) {      // pions
      Pt = max( 0.05, Pt );
      Pt = min( 0.4, Pt );
      if ( DataPeriod == 2009 ) {
         W = 0.9878 + CUBE(0.0219/Pt);
         // OLD: prod-te6; -te10
//          W = 0.9871 + CUBE(0.0224/Pt);
      } else if ( DataPeriod == 2012 ) {
         W = 0.9859 + CUBE(0.02974/Pt);
         // OLD: prod-te6; -te10
//          W = 0.9843 + CUBE(0.03015/Pt);
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

//----------------------------------------------------------------------
string RmDig(const string& str) {
//----------------------------------------------------------------------
   static regex re_dig("[0-9]");
   return regex_replace(str,re_dig,"");
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

//----------------------------------------------------------------------
void DebPrint(const vector<TH2D*> hst, string msg="") {
//----------------------------------------------------------------------
   cout << msg << endl;
   for ( auto h : hst ) {
      cout << "       name: " << h->GetName()
           << " title: " << h->GetTitle()
           << " X: " << h->GetXaxis()->GetTitle()
           << " Y: " << h->GetYaxis()->GetTitle() << endl;
   }
}

//-------------------------------------------------------------------------
void get_hst( string fname, int type, vector<TH2D*>& hst ) {
//-------------------------------------------------------------------------
// get set of histograms from file "fname" according of "type"

   // names of histograms in PsipPiPiKK
   const vector<string> n_trk {  // track efficiency
      "S3_Kp5",  "S3_Kp6",  "S3_Km5",  "S3_Km6",
      "S5_Pip5", "S5_Pip6", "S5_Pim5", "S5_Pim6"
   };
   const vector<string> n_pid {  // PID efficiency
      "S6_Kp0",  "S6_Kp1",  "S6_Km0",  "S6_Km1",
      "S6_Pip0", "S6_Pip1", "S6_Pim0", "S6_Pim1"
   };

   bool isMC = (fname.find("data") == string::npos);

   vector<string> hisn;
   if ( type == 0 ) {            // "trk"
      hisn = n_trk;
   } else if ( type == 1 ) {     // "trk_A"
      hisn = n_trk;
      for ( auto& hn : hisn ) {
         hn += string("_A");
      }
   } else if ( type == 2 ) {     // "trkSB"
      hisn = n_trk;
      for ( auto& hn : hisn ) {
         hn += string("sb");
      }
   } else if ( type == 3 ) {     // "trkSB_A"
      hisn = n_trk;
      for ( auto& hn : hisn ) {
         hn += string("sb_A");
      }
   } else if ( type == 10 ) {    // "pid"
      hisn = n_pid;
      if ( isMC && use_rew == 2 ) {
         for ( int i : {1,3,5,7}  ) {
            hisn[i] += string("W");
         }
      }
   } else if ( type == 11 ) {    // "pid_A"
      hisn = n_pid;
      if ( isMC && use_rew == 2 ) {
         for ( int i : {1,3,5,7}  ) {
            hisn[i] += string("W");
         }
      }
      for ( auto& hn : hisn ) {
         hn += string("_A");
      }
   } else {
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
   froot->cd("PipPimKpKm");

   for ( int i = 0; i < Nh; ++i ) {
      hst[i] = (TH2D*)gROOT->FindObject(hisn[i].c_str());
      if ( !hst[i] ) {
         cout << " can not find histo:" << hisn[i] << endl;
         exit(1);
      }
      hst[i]->Sumw2(true);
   }

   // Names and Titles of histograms
   string Kpi[2] = {"K","#pi"};
   string Zs[2] = {"^{+}","^{-}"};
   string SNY[2] = {" No "," Yes"};
   string str_dd; // data description
   if ( type < 10 ) {
      str_dd += "trk";
   } else {
      str_dd += "pid";
   }
   if ( !isMC ) {
      str_dd += " data";
   } else {
      str_dd += " MC";
   }
   if ( type == 1 || type == 3 || type == 11 ) {
      str_dd += "(A)";
   }
   if ( type == 2 || type == 3 ) {
      str_dd += " sideband";
   }
   for ( int i = 0; i < Nh; i++) {
      string name = ((isMC) ? "M" : "D")
                  + string( hst[i]->GetName() ).substr(1);
      hst[i]->SetName( name.c_str() );

      int isKp = (i/4)%2; // 0,1 for K,pi
      int isZ  = (i/2)%2; // 0,1 for plus,minus
      int isNY = i%2;     // 0,1 for No/Yes selection (5/6 tracks)
      string title = Kpi[isKp] + Zs[isZ] + " " + str_dd + SNY[isNY]
                   + "; cos(#Theta);"
                   + " P_{t}, GeV/c";

      hst[i]->SetTitle(title.c_str());
      SetHstFace(hst[i]);
   }

   // re-weighting PID hst
   if ( isMC  && use_rew == 1 && type >= 10 ) {
      int dataP = 2012;
      if ( fname.find("_09") != string::npos ) {
         dataP = 2009;
      }
      ReWeightTrkPid(dataP,1,hst[1]); // pid K+
      ReWeightTrkPid(dataP,1,hst[3]); // pid K-
      ReWeightTrkPid(dataP,0,hst[5]); // pid pi+
      ReWeightTrkPid(dataP,0,hst[7]); // pid pi+
      cout << " ===> re-weighting " << fname
           << " dataP= " << dataP << endl;
   }

//    DebPrint(hst, "get_hst: ");
}

//-------------------------------------------------------------------------
void get_hst_SBsub( string fname, int type, vector<TH2D*>& hst ) {
//-------------------------------------------------------------------------
// The same as get_hst() but subtract side-band from central region

   vector<TH2D*> hstCR, hstSB;
   if ( type == 0 || type == 1 ) {
      get_hst(fname, type, hstCR);
      get_hst(fname, type+2, hstSB);
   } else {
      cout << " ERROR:: get_hst_SBsub: type= " << type << endl;
      exit(1);
   }

   int Nh = hstCR.size();
   hst.clear();
   hst.resize(Nh,nullptr);

   for ( int i = 0; i < Nh; ++i ) {
      string name = string(hstCR[i]->GetName()) + "_SBS";
      string title = string(hstCR[i]->GetTitle())
                   + string(" sideband subtructed");

      hst[i] = (TH2D*)hstCR[i]->Clone( name.c_str() );
      hst[i]->SetTitle( title.c_str() );

      hst[i]->Add(hstSB[i],-1.);

      // check negative content: What we should do whis this ?
      int nx = hst[i]->GetNbinsX();
      int ny = hst[i]->GetNbinsY();
      for ( int jx = 1; jx <= nx; jx++ ) {
         for ( int jy = 1; jy <= ny; jy++ ) {
            double t = hst[i]->GetBinContent(jx,jy);
            if ( t < 0. ) {
               cout << " NEGATIVE BIN in: " << hst[i]->GetTitle() << endl
                    << " x= " << hst[i]->GetXaxis()->GetBinCenter(jx)
                    << " y= " << hst[i]->GetYaxis()->GetBinCenter(jy)
                    << " content= " << t << endl;
            }
         }
      }
   } // end of for (i)

//    DebPrint(hst, "get_hst_SBsub: ");
}

//-------------------------------------------------------------------------
void SumPM( vector<TH2D*>& hst ) {
//-------------------------------------------------------------------------
// Summarize histograms for positive and negative particles
   // kaons 0,1 and pions 4,5
   for ( int i : {0,1,4,5} ) {
      hst[i]->Add(hst[i+2]);
      string title(hst[i]->GetTitle());
      string::size_type idx = title.find("+");
      title = title.substr(0,idx) + "#pm" + title.substr(idx+1);
      hst[i]->SetTitle(title.c_str());
   }
   // shift and resize
   hst[2] = hst[4];
   hst[3] = hst[5];
   hst.resize(4);
//    DebPrint(hst, "SumPM: ");
}

//-------------------------------------------------------------------------
void get_eff( string fname, int type_eff, bool isSum, vector<TH2D*>& eff ) {
//-------------------------------------------------------------------------
// calculate efficiencies in 2D
// 1) if type_eff < 20: type_eff == type in get_hst()
// 2) if type_eff = 20,21 => eff(trk)*eff(PID) (A)
// 3) isSum => call SumPM

   type_eff = abs(type_eff);
   bool isTP = (type_eff >= 20);
   int type = (isTP) ? type_eff - 20 : type_eff;
   string prefix = (isTP) ? "ept_" : ( (type < 10) ? "etrk_" : "epid_" );

   vector<TH2D*> hst, hpid;
   if ( type == 2 || type == 3 ) {
      get_hst_SBsub(fname, type, hst);
   } else {
      get_hst(fname, type, hst);
      if ( isTP ) {
         get_hst(fname, type+10, hpid);
      }
   }

   if ( isSum ) {
      SumPM( hst );
      if ( isTP ) {
         SumPM( hpid );
      }
   }

   int Nh = hst.size();
   eff.clear();
   eff.resize(Nh/2,nullptr);

   for ( int i = 0; i < Nh/2; ++i ) {
      int j = 2*i; // index for hst

      string name = prefix + RmDig(string( hst[j]->GetName() ));

      eff[i] = (TH2D*)hst[j]->Clone( name.c_str() );
      eff[i]->Add(hst[j+1]);
      if ( !isTP ) {
        eff[i]->Divide(hst[j+1],eff[i],1.,1.,"B");
      } else {
        eff[i]->Divide(hpid[j+1],eff[i],1.,1.,"B");
      }

      // Titles of histograms
      string title( hst[j]->GetTitle() );
      title = title.substr(0,title.size()-3)
            + (( !isTP ) ? "eff" : "eff(trk*pid)");
      eff[i]->SetTitle(title.c_str());
      eff[i]->GetZaxis()->SetTitle("efficiency");
      SetHstFace(eff[i]);
   }
//    DebPrint(eff, "get_eff_1D: ");
}

//-------------------------------------------------------------------------
void get_ratio( const vector<TH2D*>& effdat, const vector<TH2D*>& effmc,
                   vector<TH2D*>& rat) {
//-------------------------------------------------------------------------
// calculate ratio DATA/MC

   int Nh = effdat.size();
   rat.clear();
   rat.resize(Nh,nullptr);
   for ( int i = 0; i < Nh; ++i ) {
      string name(effdat[i]->GetName());
      name = "r_" + name;
      rat[i]=(TH2D*)effdat[i]->Clone( name.c_str() );
      rat[i]->Divide(effmc[i]);
   }

   // Titles of histograms
   for ( int i = 0; i < Nh; ++i ) {
      string title( effmc[i]->GetTitle() );
      string::size_type idx = title.find("MC");
      title = title.substr(0,idx) + "data/" + title.substr(idx);
      rat[i]->SetTitle(title.c_str());
      rat[i]->GetZaxis()->SetTitle("data/MC");
      SetHstFace(rat[i]);
   }

//    DebPrint(rat, "get_ratio: ");
}

//-------------------------------------------------------------------------
int ResetCanvas(TCanvas* c, int nx, int ny) {
//-------------------------------------------------------------------------
   c->cd();
   c->Clear();
   int n = nx*ny;
   if ( n == 1 ) {
      c -> SetCanvasSize(700,700);
      gPad -> SetGrid();
      gStyle -> SetStatY(0.89);
   } else {
      c -> SetCanvasSize(700,900);
      c -> Divide(nx,ny);
      gStyle -> SetStatY(0.87);
      for (int i = 1; i <= n; ++i ) {
         c->cd(i);
         gPad->SetGrid();
      }
   }
   return n;
}

//-------------------------------------------------------------------------
void Fit_cos_Pt(const vector<TH2D*>& rat, string pdf) {
//-------------------------------------------------------------------------
// First fit the ratio as the function of "cos(Theta)"
// and then fit the obtained parameters as the function of "Pt"

//    TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   TCanvas* c1 = new TCanvas("c1","...",0,0,700,900);

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

   // function to the first fit
   TF1* pol0 = (TF1*)gROOT->GetFunction("pol0")->Clone();
   pol0->SetLineColor(kRed);

   // functions to the second fit kaons
   // nch - max order of Chebyshev polynomials
   int nch = 4; // 2012
   if ( date == 2009 ) {
      nch = 1;
   }
   if ( use_rew > 0 ) nch = 0;

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
   double xmin = 0.1, xmax = 1.4;
//    TF1* tfK = new TF1("tfK", "cheb4", xmin, xmax);
   TF1* tfK = new TF1("tfK", Lchb, xmin, xmax,nch+1);
   tfK->SetLineColor(kRed);

   TLegend* legK = new TLegend(0.20,0.82,0.50,0.89);
   if ( sumPM ) {
      legK->SetHeader( Form("Kaons: %i",date),"C" );
   }

   auto fPi = [](const double* x, const double* p) -> double {
      double t = p[1]/x[0];
      return p[0] + t*t*t;
   };
   TF1* tfPi = new TF1("tfPi", fPi, 0.06, 0.4, 2);
   tfPi->SetParNames("a","b");
   tfPi->SetParameters(0.98, 0.02);

   TLegend* legPi = new TLegend(0.20,0.82,0.50,0.89);
   if ( sumPM ) {
      legPi->SetHeader( Form("Pions: %i",date),"C" );
   }

   int Nh = rat.size();
   for (int ih = 0; ih < Nh; ++ih ) {
//       cout << rat[ih]->GetTitle() << endl;
      int nx = rat[ih]->GetNbinsX();
      int ny = rat[ih]->GetNbinsY();

      bool isK = (ih < Nh/2);
      int ncd = ( (isK) ? ResetCanvas(c1,2,7) : ResetCanvas(c1,2,4) );

      if ( !sumPM ) {
         if ( isK ) {
            if( ih < Nh/4 ) {
               legK->SetHeader( Form("K^{+} : %i",date),"C" );
            } else {
               legK->SetHeader( Form("K^{-} : %i",date),"C" );
            }
         } else {
            if( ih < 3*Nh/4 ) {
               legPi->SetHeader( Form("#pi^{+} : %i",date),"C" );
            } else {
               legPi->SetHeader( Form("#pi^{-} : %i",date),"C" );
            }
         }
      }

      vector<double> xx(ny), yy(ny);
      vector<double> ex(ny,1e-3), ey(ny);
      int j = 1;
      for ( int iy = 1; iy <= ny; ++iy ) {
         c1->cd(j);
         string name = "xx_" + to_string(iy);
         TH1D* tmp = rat[ih]->ProjectionX(name.c_str(),iy,iy);
         if ( isK ) {
            tmp->SetAxisRange(0.8,1.2,"Y"); // ** RANGE-1 **
         } else {
            tmp->SetAxisRange(0.8,1.2,"Y");
         }
         double pt = rat[ih]->GetYaxis()->GetBinCenter(iy);
         tmp->SetTitle(Form("#bf{<Pt> = %.2f GeV}",pt));
         SetHstFace(tmp);

         tmp->Fit(pol0,"Q","",-0.799,0.799);
         xx[iy-1] = pt;
         yy[iy-1] = pol0->GetParameter(0);
         ey[iy-1] = pol0->GetParError(0);
//          cout << " iy= " << iy << " pt= " << pt
//               << " r= " << yy[iy-1] << " +/- " << ey[iy-1] << endl;

         j = j + 1;
         if ( iy == ny ) {
            for (; j<=ncd; ++j ) {
                c1->cd(j);
                gPad->Clear();
            }
         }
         if ( j > ncd ) {
            j = 1;
            addPDF(gPad);
         }
      }

      // second fit
      ResetCanvas(c1,1,1);
      TGraphErrors* gX = new TGraphErrors(ny,xx.data(),yy.data(),
                                             ex.data(),ey.data() );
      gX->GetHistogram()->SetMaximum(1.1); // ** RANGE-2 **
      gX->GetHistogram()->SetMinimum(0.9);
      gX->SetMarkerColor(kBlue);
      gX->SetMarkerStyle(21);
//       SetHstFace(gX -> GetHistogram());
      gX->Draw("AP");
      gX->GetYaxis()->SetTitleOffset(1.4);

      if ( isK ) {
//          gX->SetTitle("Kaons: chebyshev polynomials"
         gX->SetTitle(
                      "; P_{t}, GeV/c"
                      "; #epsilon(DATA)/#epsilon(MC)"
         );
         if ( date == 2012  || use_rew > 0 ) {
            gX->Fit("tfK","REX0");
         } else if ( date == 2009 ) {
            gX->Fit("tfK","EX0","", 0.2, 1.4);
         }
//          gPad->WaitPrimitive(); // pause
//          tfK->Draw("SAME");

         legK -> Draw();
      } else {
//          gX->SetTitle("Pions: a + (b/P_{t})^{3}"
         gX->SetTitle(
                      "; P_{t}, GeV/c"
                      "; #epsilon(DATA)/#epsilon(MC)"
         );
         if ( use_rew == 0 ) {
            gX->Fit("tfPi","REX0");
         } else {
            gX->Fit("pol0","REX0");
         }
         legPi -> Draw();
      }

      addPDF(gPad);
   }

   if ( ispdf ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
   }
}

//-------------------------------------------------------------------------
void Fit_Pt_cos(const vector<TH2D*>& rat, string pdf) {
//-------------------------------------------------------------------------
// First fit the ratio as the function of "Pt"
// and then fit the obtained parameters as the function of "cos(Theta)"

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   int ncx = 2;
   int ncy = 2;

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

   // function to the first fit
   TF1* pol0 = (TF1*)gROOT->GetFunction("pol0")->Clone();
   pol0->SetLineColor(kRed);

   int Nh = rat.size();
   for (int ih = 0; ih < Nh; ++ih ) {
//       cout << rat[ih]->GetTitle() << endl;
      int nx = rat[ih]->GetNbinsX();
      int ny = rat[ih]->GetNbinsY();

      int ncd = ResetCanvas(c1,ncx,ncy);

      bool isK = (ih < Nh/2);

      vector<double> xx(nx), yy(nx);
      vector<double> ex(nx,1e-3), ey(nx);
      int j = 1;
      for ( int ix = 1; ix <= nx; ++ix ) {
         c1->cd(j);
         string name = "yy_" + to_string(ix);
         TH1D* tmp = rat[ih]->ProjectionY(name.c_str(),ix,ix);
         tmp->SetAxisRange(0.7,1.3,"Y");
         double ct = rat[ih]->GetXaxis()->GetBinCenter(ix);
         string title = " <cos>= " + to_string(ct);
         tmp->SetTitle(title.c_str());

         tmp->Fit(pol0,"Q");
         xx[ix-1] = ct;
         yy[ix-1] = pol0->GetParameter(0);
         ey[ix-1] = pol0->GetParError(0);

//          cout << " ix= " << ix << " ct= " << ct
//               << " r= " << yy[iy-1] << " +/- " << ey[iy-1] << endl;

         j = j + 1;
         if ( ix == nx ) {
            for (; j<=ncd; ++j ) {
                c1->cd(j);
                gPad->Clear();
            }
         }
         if ( j > ncd ) {
            j = 1;
            addPDF(gPad);
         }
      }

      // second fit
      ResetCanvas(c1,1,1);
      TGraphErrors* gY = new TGraphErrors(nx,xx.data(),yy.data(),
                                             ex.data(),ey.data() );
//       gY->GetHistogram()->SetMaximum(1.1);
//       gY->GetHistogram()->SetMinimum(0.9);
      gY->SetMarkerColor(kBlue);
      gY->SetMarkerStyle(21);
      gY->Draw("ACP");

      gY->Fit("pol0","Q","",-0.8,0.8);

      addPDF(gPad);
   }

   if ( ispdf ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
   }
}

//-------------------------------------------------------------------------
void FitRatio(bool isSum = true) {
//-------------------------------------------------------------------------
   int idx = (date==2009) ? 0 : 2; // date is global

   vector<TH2D*> eff_d, eff_mc, rat;
   int type = 20;
   get_eff(fnames[idx], type, isSum, eff_d);
   get_eff(fnames[idx+1], type, isSum, eff_mc);
   get_ratio( eff_d, eff_mc, rat );

   // debug print K-
//    PrintDeb(eff_d[1]);
//    PrintDeb(eff_mc[1]);
//    PrintDeb(rat[1]);

//    string pdf = ""; // debug
   string pdf = string("trk_")
              + ( (use_rew==2) ? "OW_" : ((use_rew==1) ? "TW_" : "fit") )
              + ( (isSum) ? "PM_" : "" )
              + to_string(date) + string(".pdf");

   Fit_cos_Pt(rat,pdf);
//    Fit_Pt_cos(rat,pdf);
}

//-------------------------------------------------------------------------
void Fit2_cos_Pt(const vector<TH2D*>& hst, string pdf) {
//-------------------------------------------------------------------------
// First fit the hst as the function of "cos(Theta)"
// and then fit the obtained parameters as the function of "Pt"

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   int ncx = 2;
   int ncy = 2;

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

   // function to the first fit
   TF1* pol0 = (TF1*)gROOT->GetFunction("pol0")->Clone();
   pol0->SetLineColor(kRed);
//    auto mf = [](const double* x, const double* p) -> double {
//       return p[0] + p[1]*x[0]*x[0];
//    };
//    TF1* tmf = new TF1("tmf", mf, 0.06, 0.4, 2);
//    tmf->SetParNames("p0","p2");
//    tmf->SetParameters(0.98, 0.);
//    tmf->SetLineColor(kRed);

   int Nh = hst.size();
   for (int ih = 0; ih < Nh; ++ih ) {
//       cout << hst[ih]->GetTitle() << endl;
      int nx = hst[ih]->GetNbinsX();
      int ny = hst[ih]->GetNbinsY();

      int ncd = ResetCanvas(c1,ncx,ncy);

      bool isK = (ih < Nh/2);

      vector<double> xx(ny), yy(ny);
      vector<double> ex(ny,1e-3), ey(ny);
//       vector<double> y2(ny), ey2(ny);
      int j = 1;
      for ( int iy = 1; iy <= ny; ++iy ) {
         c1->cd(j);
         string name = "xx_" + to_string(iy);
         TH1D* tmp = hst[ih]->ProjectionX(name.c_str(),iy,iy);
         tmp->SetAxisRange(0.,2.,"Y");
         double pt = hst[ih]->GetYaxis()->GetBinCenter(iy);
         string title = "<Pt> = " + to_string(pt);
         tmp->SetTitle(title.c_str());

         tmp->Fit(pol0,"Q","",-0.799,0.799);
//          tmp->Fit(tmf,"Q","",-0.799,0.799);
         xx[iy-1] = pt;
         yy[iy-1] = pol0->GetParameter(0);
         ey[iy-1] = pol0->GetParError(0);
//          yy[iy-1] = tmf->GetParameter(0);
//          ey[iy-1] = tmf->GetParError(0);
//          y2[iy-1] = tmf->GetParameter(1);
//          ey2[iy-1] = tmf->GetParError(1);

//          cout << " iy= " << iy << " pt= " << pt
//               << " r= " << yy[iy-1] << " +/- " << ey[iy-1] << endl;

         j = j + 1;
         if ( iy == ny ) {
            for (; j<=ncd; ++j ) {
                c1->cd(j);
                gPad->Clear();
            }
         }
         if ( j > ncd ) {
            j = 1;
            addPDF(gPad);
         }
      }

      // second fit
      ResetCanvas(c1,1,1);
      TGraphErrors* gX = new TGraphErrors(ny,xx.data(),yy.data(),
                                             ex.data(),ey.data() );
//       gX->GetHistogram()->SetMaximum(1.3);
//       gX->GetHistogram()->SetMinimum(0.3);
      gX->SetMarkerColor(kBlue);
      gX->SetMarkerStyle(21);
      gX->Draw("ACP");

      gX->SetTitle("p0 + p1 #times P_{t}; P_{t}, GeV/c^{2}");
      gX->Fit("pol1","Q");
//       gX->Fit("pol2","Q");

      addPDF(gPad);

//       TGraphErrors* gX2 = new TGraphErrors(ny,xx.data(),y2.data(),
//                                              ex.data(),ey2.data() );
//       gX2->SetMarkerColor(kBlue);
//       gX2->SetMarkerStyle(21);
//       gX2->Draw("ACP");
//       gX2->Fit("pol1","Q");
//       addPDF(gPad);

   }

   if ( ispdf ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
   }
}

/* OLD - wrong!
//-------------------------------------------------------------------------
void FitCr(int date) {
//-------------------------------------------------------------------------
// calculate and fit ratio:
//              [eff(data)/(1-eff(data))] / [eff(mc)/(1-eff(mc))]

   int idx = (date==2009) ? 0 : 2;

   vector<TH2D*> eff_d, eff_mc;
   bool isSum = true;
   int type = 20;
   get_eff(fnames[idx], type, isSum, eff_d);
   get_eff(fnames[idx+1], type, isSum, eff_mc);

   int Nh = eff_d.size();
   vector<TH2D*> Cr(Nh,nullptr);
   for (int ih = 0; ih < Nh; ++ih ) {
      int nx = eff_d[ih]->GetNbinsX();
      int ny = eff_d[ih]->GetNbinsY();

      string name = "Cor_" + string( eff_d[ih]->GetName() );
      Cr[ih] = (TH2D*)eff_d[ih]->Clone( name.c_str() );
      for ( int ix = 0; ix <= nx+1; ++ix ) {
         for ( int iy = 0; iy <= ny+1; ++iy ) {
            double ed  = eff_d[ih]->GetBinContent(ix,iy);
            double eed = eff_d[ih]->GetBinError(ix,iy);
            double em  = eff_mc[ih]->GetBinContent(ix,iy);
            double eem = eff_mc[ih]->GetBinError(ix,iy);

            double ed1 = 1 - ed;
            double em1 = 1 - em;

            // prevention of infinity
            if ( fabs(em)  < 1e-4 ||
                 fabs(em1) < 1e-4 ||
                 fabs(ed)  < 1e-4 ||
                 fabs(ed1) < 1e-4
               ) {
              Cr[ih]->SetBinContent(ix,iy,0);
              Cr[ih]->SetBinError(ix,iy,0);
              continue;
            }

            double cr = (ed*em1)/(ed1*em);
            double ecr = cr * sqrt( SQ(eed/(ed*ed1)) + SQ(eem/(em*em1)) );
            Cr[ih]->SetBinContent(ix,iy,cr);
            Cr[ih]->SetBinError(ix,iy,ecr);
         }
      }
   }

//    string pdf = ""; // debug
   string pdf = string("trk_corrPM_") + to_string(date)
              + string(".pdf");

   Fit2_cos_Pt(Cr, pdf);
//    Fit_Pt_cos(Cr, pdf);
}
*/

//-------------------------------------------------------------------------
void trk_eff_fit() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(62);
   gStyle->SetStatFont(62);
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
//    gStyle->SetStatW(0.25);
   InitFnames();

   FitRatio(sumPM);
}
