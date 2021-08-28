// Study of the track reconstruction efficiency for pi and K:
// I)  plot 1D distributions  ->  Trkeff_[date].pdf
// II) Kolmogorov–Smirnov test to compare K+ and K- (pi+ and pi-)
// pdftk Trkeff__2012.pdf cat 9 output treff_12_Keff.pdf

// #include <regex>
//----------------------------------------------------------------------
// GLOBAL:
// name of folder with root files
static string dir;
// names of files to use
static vector<string> fnames;

const int date = 2009;

// use weighted histograms:
// <0 - do KolmogorovTest and exit
// 0 - no weights;
// 1 - calculate weights here; } plot only TRK*PID eff
// 2 - use "W" histograms.     }
const int use_rew = 0;

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
      X->SetTitleSize(0.045);
      X->SetTitleOffset(0.9);
   }
   TAxis* Y = hst->GetYaxis();
   if ( Y ) {
      Y->SetLabelFont(62);
      Y->SetLabelSize(0.04);
      Y->SetTitleFont(62);
      Y->SetTitleSize(0.04);
      Y->SetTitleOffset(1.35);
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
/*
template<typename T>
void DebPrint(const vector<T> hst, string msg="") {
//----------------------------------------------------------------------
   cout << msg << endl;
   for ( auto h : hst ) {
      cout << "       name: " << h->GetName()
           << " title: " << h->GetTitle()
           << " X: " << h->GetXaxis()->GetTitle()
           << " Y: " << h->GetYaxis()->GetTitle() << endl;
   }
}*/

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

//    DebPrint<TH2D*>(hst, "get_hst: ");
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

//    DebPrint<TH2D*>(hst, "get_hst_SBsub: ");
}

//-------------------------------------------------------------------------
void get_eff_1D( string fname, int type_eff, vector<TH1D*>& eff ) {
//-------------------------------------------------------------------------
// Get projections and calculate efficiencies
// if type_eff < 20: type_eff == type in get_hst()
// 20 => eff(trk)*eff(PID); 21=>(A)

   bool isTP = (type_eff >= 20);

   vector<TH2D*> hst, hpid;
   int type = (isTP) ? type_eff - 20 : type_eff;
   if ( type == 2 || type == 3 ) {
      get_hst_SBsub(fname, type, hst);
   } else {
      get_hst(fname, type, hst);
      if ( isTP ) {
         get_hst(fname, type+10, hpid);
      }

      // debug print content (K-)
//       if ( type_eff == 20 ) {
//          PrintDeb(hpid[3]);
//          PrintDeb(hst[2]);
//          PrintDeb(hst[3]);
//       }
   }

   string prefix = (isTP) ? "ept_" : ( (type < 10) ? "etrk_" : "epid_" );

   int Nh = hst.size();
   eff.clear();
   eff.resize(Nh,nullptr);

   for ( int i = 0; i < Nh/2; ++i ) {
      int j = 2*i; // index for hst

      int nx = hst[j]->GetNbinsX();
      int ny = hst[j]->GetNbinsY();

      string name = prefix + RmDig(string( hst[j]->GetName() ));
      string pxn = name + "_ct";
      string pyn = name + "_pt";

      eff[j] = (TH1D*)hst[j]->ProjectionX(pxn.c_str(),1,ny);
      TH1D* tmp = (TH1D*)hst[j+1]->ProjectionX("tmp",1,ny);
      eff[j]->Add(tmp);
      if ( isTP ) {
         tmp = (TH1D*)hpid[j+1]->ProjectionX("tmp",1,ny);
      }
      eff[j]->Divide(tmp,eff[j],1.,1.,"B");

      eff[j+1] = (TH1D*)hst[j]->ProjectionY(pyn.c_str(),2,nx-1);
      tmp = (TH1D*)hst[j+1]->ProjectionY("tmp",2,nx-1);
      eff[j+1]->Add(tmp);
      if ( isTP ) {
         tmp = (TH1D*)hpid[j+1]->ProjectionY("tmp",2,nx-1);
      }
      eff[j+1]->Divide(tmp,eff[j+1],1.,1.,"B");
   }

   // Titles of histograms
   for ( int i = 0; i < Nh; i++) {
      string title( hst[i]->GetTitle() );
      if ( isTP ) {
         string rm = "trk ";
         string::size_type idx = title.find(rm);
         if ( idx != string::npos ) {
            title = title.substr(0,idx) + title.substr(idx+rm.size());
         }
      }
      title = title.substr(0,title.size()-3)
            + (( !isTP ) ? "eff" : "eff(trk*pid)");
      eff[i]->SetTitle(title.c_str());
      eff[i]->GetYaxis()->SetTitle("efficiency");
      SetHstFace(eff[i]);
   }

//    DebPrint<TH1D*>(eff, "get_eff_1D: ");
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
      rat[i]->GetYaxis()->
         SetTitle("#epsilon (data) / #epsilon (#scale[0.7]{MC})");
      SetHstFace(rat[i]);
      rat[i]->SetMarkerStyle(20);
      rat[i]->SetLineColor(kBlack);
   }

//    DebPrint<TH1D*>(rat, "get_ratio: ");
}

//-------------------------------------------------------------------------
void get_Cr( const vector<TH1D*>& effdat, const vector<TH1D*>& effmc,
             vector<TH1D*>& Cr) {
//-------------------------------------------------------------------------
// calculate ratio [eff(data)/(1-eff(data))] / [eff(mc)/(1-eff(mc))]

   int Nh = effdat.size();
   Cr.clear();
   Cr.resize(Nh,nullptr);
   for ( int i = 0; i < Nh; ++i ) {
      int nx = effdat[i]->GetNbinsX();
      int ny = effdat[i]->GetNbinsY();

      string name(effdat[i]->GetName());
      name = "Cr_" + name;
      Cr[i]=(TH1D*)effdat[i]->Clone( name.c_str() );
      for ( int ix = 0; ix <= nx+1; ++ix ) {
         for ( int iy = 0; iy <= ny+1; ++iy ) {
            double ed  = effdat[i]->GetBinContent(ix,iy);
            double eed = effdat[i]->GetBinError(ix,iy);
            double em  = effmc[i]->GetBinContent(ix,iy);
            double eem = effmc[i]->GetBinError(ix,iy);

            double ed1 = 1 - ed;
            double em1 = 1 - em;

            // prevention of infinity
            if ( fabs(em)  < 1e-4 ||
                 fabs(em1) < 1e-4 ||
                 fabs(ed)  < 1e-4 ||
                 fabs(ed1) < 1e-4
               ) {
              Cr[i]->SetBinContent(ix,iy,0);
              Cr[i]->SetBinError(ix,iy,0);
              continue;
            }

            double cr  = (ed*em1)/(ed1*em);
            double ecr = cr*sqrt( SQ(eed/(ed*ed1)) + SQ(eem/(em*em1)) );
            Cr[i]->SetBinContent(ix,iy,cr);
            Cr[i]->SetBinError(ix,iy,ecr);
         }
      } // end for(ix)
   }

   // Titles of histograms
   for ( int i = 0; i < Nh; ++i ) {
      string title( effmc[i]->GetTitle() );
      string::size_type idx = title.find("MC");
      title = title.substr(0,idx) + "Cr" + title.substr(idx+2);
      Cr[i]->SetTitle(title.c_str());
      Cr[i]->GetYaxis()->SetTitle("Cr");
      SetHstFace(Cr[i]);
      Cr[i]->SetMarkerStyle(20);
      Cr[i]->SetLineColor(kBlack);
   }

//    DebPrint<TH1D*>(Cr, "get_Cr: ");
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
   }

//    DebPrint<TH1D*>(diff, "get_diff: ");
}

//-------------------------------------------------------------------------
void plot_pict(string pdf) {
//-------------------------------------------------------------------------
   int idx = (date==2009) ? 0 : 2; // date is global

   vector<TH1D*> eff_d, eff_mc, rat0;
   get_eff_1D(fnames[idx], 0, eff_d);
   get_eff_1D(fnames[idx+1], 0, eff_mc);
   get_ratio( eff_d, eff_mc, rat0 );
//    get_Cr( eff_d, eff_mc, rat0 );

   vector<TH1D*> eff_dA, eff_mcA, ratA;
   vector<TH1D*> diff_d, diff_mc, diff_rat;
   get_eff_1D(fnames[idx], 1, eff_dA);
   get_eff_1D(fnames[idx+1], 1, eff_mcA);
   get_ratio( eff_dA, eff_mcA, ratA );
//    get_Cr( eff_dA, eff_mcA, ratA );
   get_diff( eff_d, eff_dA, diff_d );
   get_diff( eff_mc, eff_mcA, diff_mc );
   get_diff( rat0, ratA, diff_rat );

   vector<TH1D*> eff_dp, eff_mcp, rat_p;
   get_eff_1D(fnames[idx], 10, eff_dp);
   get_eff_1D(fnames[idx+1], 10, eff_mcp);
   get_ratio( eff_dp, eff_mcp, rat_p );
//    get_Cr( eff_dp, eff_mcp, rat_p );

   vector<TH1D*> eff_dpA, eff_mcpA, rat_pA;
   vector<TH1D*> diff_dp, diff_mcp, diff_ratp;
   get_eff_1D(fnames[idx], 11, eff_dpA);
   get_eff_1D(fnames[idx+1], 11, eff_mcpA);
   get_ratio( eff_dpA, eff_mcpA, rat_pA );
//    get_Cr( eff_dpA, eff_mcpA, rat_pA );
   get_diff( eff_dp, eff_dpA, diff_dp );
   get_diff( eff_mcp, eff_mcpA, diff_mcp );
   get_diff( rat_p, rat_pA, diff_ratp );

   vector<TH1D*> etp_d, etp_mc, rat_tp;
   get_eff_1D(fnames[idx], 20, etp_d);
   get_eff_1D(fnames[idx+1], 20, etp_mc);
   get_ratio( etp_d, etp_mc, rat_tp );
//    get_Cr( etp_d, etp_mc, rat_tp );

   // debug print K-
//    PrintDeb(etp_d[2]);
//    PrintDeb(etp_d[3]);
//    PrintDeb(etp_mc[2]);
//    PrintDeb(etp_mc[3]);
//    PrintDeb(rat_tp[2]);
//    PrintDeb(rat_tp[3]);

   vector<TH1D*> etp_dA, etp_mcA, rat_tpA;
   vector<TH1D*> diff_tp;
   get_eff_1D(fnames[idx], 21, etp_dA);
   get_eff_1D(fnames[idx+1], 21, etp_mcA);
   get_ratio( etp_dA, etp_mcA, rat_tpA );
//    get_Cr( etp_dA, etp_mcA, rat_tpA );
   get_diff( rat_tp, rat_tpA, diff_tp );

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
   auto setPtK_X = [](TH1D* h) { h->SetAxisRange(0.1,1.399,"X"); };
   auto setPtPi_X = [](TH1D* h) { h->SetAxisRange(0.05,0.399,"X"); };
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
      RmFromTitle(eff_dp[i],"data ");
      RmFromTitle(etp_d[i],"data ");

      RmFromTitle(etp_d[i],"eff(trk*pid)");
      RmFromTitle(rat_tp[i],"eff(trk*pid)");

      setBlue(eff_d[i]);
      setBlue(diff_d[i]);
      setBlue(eff_dp[i]);
      setBlue(diff_dp[i]);
      setBlue(etp_d[i]);

      setRed(eff_mc[i]);
      setRed(diff_mc[i]);
      setRed(eff_mcp[i]);
      setRed(diff_mcp[i]);
      setRed(etp_mc[i]);

      if ( i%2 == 0 ) {
         setC_X(eff_d[i]);
         setC_X(diff_d[i]);
         setC_X(rat0[i]);
         setC_X(diff_rat[i]);
         setC_X(eff_dp[i]);
         setC_X(diff_dp[i]);
         setC_X(rat_p[i]);
         setC_X(diff_ratp[i]);
         setC_X(etp_d[i]);
         setC_X(rat_tp[i]);
         setC_X(diff_tp[i]);
      } else if ( i < Nh/2 ) {
         setPtK_X(eff_d[i]);
         setPtK_X(diff_d[i]);
         setPtK_X(rat0[i]);
         setPtK_X(diff_rat[i]);
         setPtK_X(eff_dp[i]);
         setPtK_X(diff_dp[i]);
         setPtK_X(rat_p[i]);
         setPtK_X(diff_ratp[i]);
         setPtK_X(etp_d[i]);
         setPtK_X(rat_tp[i]);
         setPtK_X(diff_tp[i]);
      } else {
         setPtPi_X(eff_d[i]);
         setPtPi_X(diff_d[i]);
         setPtPi_X(rat0[i]);
         setPtPi_X(diff_rat[i]);
         setPtPi_X(eff_dp[i]);
         setPtPi_X(diff_dp[i]);
         setPtPi_X(rat_p[i]);
         setPtPi_X(diff_ratp[i]);
         setPtPi_X(etp_d[i]);
         setPtPi_X(rat_tp[i]);
         setPtPi_X(diff_tp[i]);
      }
   }

//    TLegend* leg1 = new TLegend(0.40,0.20,0.72,0.40);
   TLegend* leg1 = new TLegend(0.35,0.25,0.65,0.45);
   leg1->AddEntry(eff_d[0], Form("#color[%i]{data %i}",
                  eff_d[0]->GetLineColor(),date),"LP");
   leg1->AddEntry(eff_mc[0], Form("#color[%i]{MC %i}",
                  eff_mc[0]->GetLineColor(),date),"LP");


   // for diff
   TF1* pl0 = (TF1*)gROOT->GetFunction("pol0")->Clone();
   pl0->SetLineColor(kGreen+2);

   // Draw:
//    TCanvas* c1 = new TCanvas("c1","...",0,0,1100,900);
   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,800);
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

   int Nh2 = Nh/2;

///////////////////////////////////////////////////////////////////
//    double eff_min = 0, eff_max = 1.;
   double eff_min = 0.4, eff_max = 1.;
   double rat_min = 0.9, rat_max = 1.1;
//    double rat_min = 0.8, rat_max = 1.2;
//    bool plot_delta_eff = false, plot_delta_rat = true;
   bool plot_delta_eff = false, plot_delta_rat = false;
   double diff_min = -0.05, diff_max = 0.05;
///////////////////////////////////////////////////////////////////

   if ( use_rew == 0 ) {

   // TRK efficiencies K
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      eff_d[i]->SetAxisRange(eff_min,eff_max,"Y");
      eff_d[i]->Draw("E");
      eff_mc[i]->Draw("E SAME");
      leg1->Draw();
   }
   addPDF(gPad);

   // TRK delta efficiencies K
   if ( plot_delta_eff ) {
      for (int i = 0; i < Nh2; ++i ) {
         c1->cd(i+1);
         diff_d[i]->SetAxisRange(diff_min,diff_max,"Y");
         diff_d[i]->Draw("E");
         diff_mc[i]->Draw("E SAME");
         line0[i%2]->Draw();
      }
      addPDF(gPad);
   }

   // TRK ratio K
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      rat0[i]->SetAxisRange(rat_min,rat_max,"Y");
      rat0[i]->Draw("E");
      line1[i%2]->Draw();
   }
   addPDF(gPad);

   // TRK delta ratio K
   if ( plot_delta_rat ) {
      for (int i = 0; i < Nh2; ++i ) {
         c1->cd(i+1);
         diff_rat[i]->SetAxisRange(diff_min,diff_max,"Y");
         diff_rat[i]->Draw("E");
         line0[i%2]->Draw();
      }
      addPDF(gPad);
   }

   // TRK efficiencies Pi
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      int j = Nh2+i;
      eff_d[j]->SetAxisRange(eff_min,eff_max,"Y");
      eff_d[j]->Draw("E");
      eff_mc[j]->Draw("E SAME");
      leg1->Draw();
   }
   addPDF(gPad);

   // TRK delta efficiencies Pi
   if ( plot_delta_eff ) {
      for (int i = 0; i < Nh2; ++i ) {
         c1->cd(i+1);
         int j = Nh2+i;
         diff_d[j]->SetAxisRange(diff_min,diff_max,"Y");
         diff_d[j]->Draw("E");
         diff_mc[j]->Draw("E SAME");
         line0[(i%2)*2]->Draw();
         leg1->Draw();
      }
      addPDF(gPad);
   }

   // TRK ratio Pi
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      int j = Nh2+i;
      rat0[j]->SetAxisRange(rat_min,rat_max,"Y");
      rat0[j]->Draw("E");
      line1[(i%2)*2]->Draw();
   }
   addPDF(gPad);

   // TRK delta ratio Pi
   if ( plot_delta_rat ) {
      for (int i = 0; i < Nh2; ++i ) {
         c1->cd(i+1);
         int j = Nh2+i;
         diff_rat[j]->SetAxisRange(diff_min,diff_max,"Y");
         diff_rat[j]->Draw("E");
         line0[(i%2)*2]->Draw();
      }
      addPDF(gPad);
   }

///////////////////////////////////////////////////////////////////

   // PID efficiencies K
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      eff_dp[i]->SetAxisRange(eff_min,eff_max,"Y");
      eff_dp[i]->Draw("E");
      eff_mcp[i]->Draw("E SAME");
      leg1->Draw();
   }
   addPDF(gPad);

   // delta PID efficiencies K
   if ( plot_delta_eff ) {
      for (int i = 0; i < Nh2; ++i ) {
         c1->cd(i+1);
         diff_dp[i]->SetAxisRange(diff_min,diff_max,"Y");
         diff_dp[i]->Draw("E");
         diff_mcp[i]->Draw("E SAME");
         line0[i%2]->Draw();
      }
      addPDF(gPad);
   }

   // PID ratio K
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      rat_p[i]->SetAxisRange(rat_min,rat_max,"Y");
      rat_p[i]->Draw("E");
      line1[i%2]->Draw();
   }
   addPDF(gPad);

   // PID delta ratio K
   if ( plot_delta_rat ) {
      for (int i = 0; i < Nh2; ++i ) {
         c1->cd(i+1);
         diff_ratp[i]->SetAxisRange(diff_min,diff_max,"Y");
         diff_ratp[i]->Draw("E");
         line0[i%2]->Draw();
      }
      addPDF(gPad);
   }

   // PID efficiencies Pi
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      int j = Nh2+i;
      eff_dp[j]->SetAxisRange(eff_min,eff_max,"Y");
      eff_dp[j]->Draw("E");
      eff_mcp[j]->Draw("E SAME");
      leg1->Draw();
   }
   addPDF(gPad);

   // PID delta efficiencies Pi
   if ( plot_delta_eff ) {
      for (int i = 0; i < Nh2; ++i ) {
         c1->cd(i+1);
         int j = Nh2+i;
         diff_dp[j]->SetAxisRange(diff_min,diff_max,"Y");
         diff_dp[j]->Draw("E");
         diff_mcp[j]->Draw("E SAME");
         line0[(i%2)*2]->Draw();
         leg1->Draw();
      }
      addPDF(gPad);
   }

   // PID ratio Pi
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      int j = Nh2+i;
      rat_p[j]->SetAxisRange(rat_min,rat_max,"Y");
      rat_p[j]->Draw("E");
      line1[(i%2)*2]->Draw();
   }
   addPDF(gPad);

   // PID delta ratio Pi
   if ( plot_delta_rat ) {
      for (int i = 0; i < Nh2; ++i ) {
         c1->cd(i+1);
         int j = Nh2+i;
         diff_ratp[j]->SetAxisRange(diff_min,diff_max,"Y");
         diff_ratp[j]->Draw("E");
         line0[(i%2)*2]->Draw();
      }
      addPDF(gPad);
   }

   } // end of if (use_rew == 0 )

///////////////////////////////////////////////////////////////////

   // TRK*PID efficiencies K
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      etp_d[i]->SetAxisRange(eff_min,eff_max,"Y");
      etp_d[i]->Draw("E");
      etp_mc[i]->Draw("E SAME");
      leg1->Draw();
   }
   addPDF(gPad);

   // TRK*PID ratio K
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      rat_tp[i]->SetAxisRange(rat_min,rat_max,"Y");
//       rat_tp[i]->Draw("E");
//       line1[i%2]->Draw();
      rat_tp[i]->Fit("pol0","Q");
   }
   addPDF(gPad);

   // TRK*PID delta ratio K
   if ( plot_delta_rat ) {
      for (int i = 0; i < Nh2; ++i ) {
         c1->cd(i+1);
         diff_tp[i]->SetAxisRange(diff_min,diff_max,"Y");
//          diff_tp[i]->Draw("E");
//          line0[i%2]->Draw();
         diff_tp[i]->Fit(pl0,"Q");
      }
      addPDF(gPad);
   }

   // TRK*PID efficiencies Pi
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      int j = Nh2+i;
      etp_d[j]->SetAxisRange(eff_min,eff_max,"Y");
      etp_d[j]->Draw("E");
      etp_mc[j]->Draw("E SAME");
      leg1->Draw();
   }
   addPDF(gPad);

   // TRK*PID ratio Pi
   for (int i = 0; i < Nh2; ++i ) {
      c1->cd(i+1);
      int j = Nh2+i;
      rat_tp[j]->SetAxisRange(rat_min,rat_max,"Y");
//       rat_tp[j]->Draw("E");
//       line1[(i%2)*2]->Draw();
      rat_tp[j]->Fit("pol0","Q");
   }
   addPDF(gPad);

   // TRK*PID delta ratio Pi
   if ( plot_delta_rat ) {
      for (int i = 0; i < Nh2; ++i ) {
         c1->cd(i+1);
         int j = Nh2+i;
         diff_tp[j]->SetAxisRange(diff_min,diff_max,"Y");
//          diff_tp[j]->Draw("E");
//          line0[(i%2)*2]->Draw();
         diff_tp[j]->Fit(pl0,"Q");
      }
      addPDF(gPad);
   }

   if ( ispdf ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
   }
}

//-------------------------------------------------------------------------
void test_KS() {
//-------------------------------------------------------------------------
// use the Kolmogorov–Smirnov test to compare K+ and K- (pi+ and pi-)
// 1D distributions of the data/MC ratio
// see http://root.cern.ch/root/html/TH1.html#TH1:KolmogorovTest

   int idx = (date==2009) ? 0 : 2;

   vector<TH1D*> eff_d, eff_mc, rat;
   get_eff_1D(fnames[idx], 20, eff_d);
   get_eff_1D(fnames[idx+1], 20, eff_mc);
   get_ratio( eff_d, eff_mc, rat );

   int Nh = rat.size();
   for (int ih : {0,1,4,5} ) {
      string title1( rat[ih]->GetTitle() );
      string title2( rat[ih+2]->GetTitle() );
      string::size_type idx = title1.find("data");
      string common = title1.substr(idx);
      string particle1 = title1.substr(0,idx);
      string particle2 = title2.substr(0,title2.find("data"));
      cout << " data-" << date << ": probability test for "
           << common << endl;
      cout << "       " << particle1 << " <--> " << particle2
           << ((ih%2==0) ? " cos(Theta)" : " P_t")
           << endl;
      rat[ih]->KolmogorovTest(rat[ih+2],"D");
      cout << endl;
   }
}

//-------------------------------------------------------------------------
void trk_eff() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(62);
   gStyle->SetLegendTextSize(0.05);
   gStyle->SetStatFont(62);
//    gStyle->SetStatFontSize(0.07);
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.3);
   gStyle->SetFitFormat(".4f");
   InitFnames();

   if ( use_rew < 0 ) {
      test_KS();
      return;
   }

//    string pdf = ""; // debug
   string pdf = string("Trkeff_")
              + ( (use_rew==2) ? "OW_" : ((use_rew==1) ? "TW_" : "") )
              + to_string(date)
              + string(".pdf");

   plot_pict(pdf);
}
