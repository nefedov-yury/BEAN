// Study of the eta -> eta reconstruction efficiency &
//                     single photon rec.efficiency
// -> Eff_[eta,ph]_[date]_[gamma,phi]eta.pdf

//--------------------------------------------------------------------
// Common parameters:
struct Params {
   // names of files to use
   vector<string> fnames;

   int date;
   int slct;    // >0: g-eta; <0: phi-eta
                // 1 - photon eff; 2 - eta eff;
   int use_rew; // 0 - no weights;
                // 1 - calculate weights

   TCut Cmcsig; // for mc-signal
   TCut Cbg;    // cuts against the background
   TCut Cph;    // cuts for selection photons
   TCut Ceta;   // cuts for selection eta

   //-----------------------------------------------------------------
   Params(int dat = 2012, int slc = 2, int rew = 0) {
   //-----------------------------------------------------------------
      date = dat;
      slct = slc;
      use_rew = rew;

      // name of folder with root-files
      string dir("prod-12eff/");
      fnames = {
         "data_09psip_all.root", "mcinc_09psip_all.root",
         "data_12psip_all.root", "mcinc_12psip_all.root"
      };
      for ( auto& fn : fnames ) {
         fn = dir + fn;
      }

      // mc-signal
      if ( slct > 0 ) {    // gamma-eta
         Cmcsig += TCut("decj==22");
      } else {             // phi-eta
         Cmcsig += TCut("decj==68");
      }

      // cuts against the background
      if ( slct > 0 ) {    // gamma-eta
         Cbg += TCut("m2fr<0.002");
         Cbg += TCut("fabs(Cg0)<0.8");
         Cbg += TCut("fabs(Cg1)<0.8 && Eg1>0.2");
      } else {             // phi-eta
         Cbg += TCut("m2fr<0.002");
      }

      // cuts for selection photons
      Cph += TCut("fabs(Cg2)<0.8||(fabs(Cg2)>0.85&&fabs(Cg2)<0.92)");
      Cph += TCut("fabs(Cg1)<0.8||(fabs(Cg1)>0.85&&fabs(Cg1)<0.92)");
      if ( slct > 0 ) {  // gamma-eta
         Cph += TCut("Eg2>0.1&&Eg2<1.4");
      } else {             // phi-eta
         Cph += TCut("Eg2>0.05&&Eg1<1.45");
      }

      // cuts for selection eta
      Ceta += TCut("fabs(Ceta)<0.9");
      Ceta += TCut("fabs(Cg2)<0.8||(fabs(Cg2)>0.85&&fabs(Cg2)<0.92)");
      if ( slct > 0 ) {  // gamma-eta
         Ceta += TCut("Peta>1.3&&Peta<1.7");
      } else {             // phi-eta
         Ceta += TCut("Peta>1.15&&Peta<1.5");
      }
   }

   //-----------------------------------------------------------------
   TFile* OpenFile(int mc = 0) { // 1 for MC
   //-----------------------------------------------------------------
      // open file
      int idx = mc + ((date==2009) ? 0 : 2);
      string dfname = fnames[idx];
      TFile* froot = TFile::Open(dfname.c_str(),"READ");
      if( froot == 0 ) {
         cerr << "can not open " << dfname << endl;
         exit(0);
      }
      if ( slct > 0 ) {
         froot->cd("PsipJpsiGammaEta");
      } else {
         froot->cd("PsipJpsiPhiEta");
      }
      return froot;
   }
};
//----------------------------------------------------------------------

//----------------------------------------------------------------------
constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

//----------------------------------------------------------------------
double ReWeightEtaEff(int DataPeriod) {
//----------------------------------------------------------------------
// This correction is based on "prod-12eff"

   double W = 1.;
   if ( DataPeriod == 2012 ) {
      W = 0.993;
   } else if ( DataPeriod == 2009 ) {
      W = 0.982;
   }
//    W = 0.990;
   return W;
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

//-------------------------------------------------------------------------
void get_hst(Params* p, int mc, vector<TH1D*>& hst) {
//-------------------------------------------------------------------------
// get histograms: mc=0,1 / data,MC

   TFile* froot = p->OpenFile(mc);

   TTree* eff_eta = (TTree*)gDirectory->Get("eff_eta");

   TCut Crec0 = p->Cbg;
   TCut Crec1 = p->Cbg;
   if ( abs(p->slct) == 1 ) {      // photon
      Crec0 += TCut("fl<0.5");
      Crec1 += TCut("fl>0.5");
      Crec0 += p->Cph;
      Crec1 += p->Cph;
   } else {                     // eta
      Crec0 += TCut("fl<1.5");
      Crec1 += TCut("fl>1.5");
      Crec0 += p->Ceta;
      Crec1 += p->Ceta;
   }
//    cout << " Crec0: " << Crec0.GetTitle() << endl;
//    cout << " Crec1: " << Crec1.GetTitle() << endl;

   hst.clear();
   hst.resize(4,nullptr);
   string hname = ((mc==0) ? "data_" : "mc_") + to_string(p->date);
   vector<string> hns = {
      hname + "_P", hname + "_P1",
      hname + "_C", hname + "_C1"
   };

   // variable binning for cos(Theta):
   double dc=0.2;
   vector<double> Cbins {-1.0,-0.92,-0.85};
   for(double c = -0.8; c < 0.85; c += dc ) {
     Cbins.push_back(c);
   }
   Cbins.push_back(0.85);
   Cbins.push_back(0.92);
   Cbins.push_back(1.0);
//    for( auto c : Cbins ) cout << c << " : ";
//    cout << endl;

   if ( p->slct == 1 ) {  // gamma-eta => one photon efficiency
      hst[0] = new TH1D( hns[0].c_str(), "#gamma ;E (GeV/c)", 14,0.05,1.45);
      hst[1] = new TH1D( hns[1].c_str(), "#gamma ;E (GeV/c)", 14,0.05,1.45);
      hst[2] = new TH1D( hns[2].c_str(), "#gamma ;cos(#Theta)",
                         Cbins.size()-1, Cbins.data() );
      hst[3] = new TH1D( hns[3].c_str(), "#gamma ;cos(#Theta)",
                         Cbins.size()-1, Cbins.data() );
   } else if ( p->slct == -1 ) {  // phi-eta => one photon efficiency
      if ( p->date == 2012 ) {
         hst[0] = new TH1D( hns[0].c_str(),"#gamma ;E (GeV/c)",14,0.05,1.45);
         hst[1] = new TH1D( hns[1].c_str(),"#gamma ;E (GeV/c)",14,0.05,1.45);
      } else {
         hst[0] = new TH1D( hns[0].c_str(),"#gamma ;E (GeV/c)",7,0.05,1.45);
         hst[1] = new TH1D( hns[1].c_str(),"#gamma ;E (GeV/c)",7,0.05,1.45);
      }
      hst[2] = new TH1D( hns[2].c_str(), "#gamma ;cos(#Theta)",
                         Cbins.size()-1, Cbins.data() );
      hst[3] = new TH1D( hns[3].c_str(), "#gamma ;cos(#Theta)",
                         Cbins.size()-1, Cbins.data() );
   } else if ( p->slct == 2 ) { // g-eta   => eta eff
      if ( p->date == 2012 ) {
         hst[0] = new TH1D( hns[0].c_str(),
               "#epsilon(#eta) ;P (GeV/c)", 8,1.3,1.7);
         hst[1] = new TH1D( hns[1].c_str(), "#eta ;P (GeV/c)", 8,1.3,1.7);
         hst[2] = new TH1D( hns[2].c_str(),
               "#epsilon(#eta) ;cos(#Theta)", 9,-0.9,0.9);
         hst[3] = new TH1D( hns[3].c_str(), "#eta ;cos(#Theta)", 9,-0.9,0.9);
      } else {
         hst[0] = new TH1D( hns[0].c_str(),
               "#epsilon(#eta);P (GeV/c)", 6,1.3,1.7);
         hst[1] = new TH1D( hns[1].c_str(), "#eta ;P (GeV/c)", 6,1.3,1.7);
         hst[2] = new TH1D( hns[2].c_str(),
               "#epsilon(#eta);cos(#Theta)", 9,-0.9,0.9);
         hst[3] = new TH1D( hns[3].c_str(), "#eta ;cos(#Theta)", 9,-0.9,0.9);
      }
   } else if ( p->slct == -2 ) { // phi-eta => eta eff
      hst[0] = new TH1D( hns[0].c_str(),
            "#epsilon(#eta);P (GeV/c)", 7,1.15,1.5);
      hst[1] = new TH1D( hns[1].c_str(), "#eta ;P (GeV/c)", 7,1.15,1.5);
      hst[2] = new TH1D( hns[2].c_str(),
            "#epsilon(#eta);cos(#Theta)", 6,-0.9,0.9);
      hst[3] = new TH1D( hns[3].c_str(), "#eta ;cos(#Theta)", 6,-0.9,0.9);
   }

   for ( auto& h : hst ) {
      h->Sumw2(true);
      SetHstFace(h);
   }

   if ( abs(p->slct) == 1 ) {  // photon eff
      eff_eta->Draw( ("Eg2>>"+hns[0]).c_str(),Crec0,"goff");
      eff_eta->Draw( ("Cg2>>"+hns[2]).c_str(),Crec0,"goff");

      eff_eta->Draw( ("Eg2>>"+hns[1]).c_str(),Crec1,"goff");
      eff_eta->Draw( ("Eg1>>+"+hns[1]).c_str(),Crec1,"goff");
      eff_eta->Draw( ("Cg2>>"+hns[3]).c_str(),Crec1,"goff");
      eff_eta->Draw( ("Cg1>>+"+hns[3]).c_str(),Crec1,"goff");

   }  else if ( abs(p->slct) == 2 ) {  // eta eff
      eff_eta->Draw( ("Peta>>"+hns[0]).c_str(),Crec0,"goff");
      eff_eta->Draw( ("Ceta>>"+hns[2]).c_str(),Crec0,"goff");

      eff_eta->Draw( ("Peta>>"+hns[1]).c_str(),Crec1,"goff");
      eff_eta->Draw( ("Ceta>>"+hns[3]).c_str(),Crec1,"goff");
   }
}

//-------------------------------------------------------------------------
void get_eff(Params* p, int mc, vector<TH1D*>& eff ) {
//-------------------------------------------------------------------------
// Get histograms and calculate efficiencies

   vector<TH1D*> hst;
   get_hst(p, mc, hst);

   eff.clear();
   eff.resize(2,nullptr);
   for( int i = 0; i < 2; i++ ) {
      int ih = 2*i;
      string name = string("eff_") + string(hst[ih]->GetName());
      eff[i]=(TH1D*)hst[ih]->Clone( name.c_str() );
      eff[i]->Add(hst[ih+1]);
      eff[i]->Divide(hst[ih+1],eff[i],1.,1.,"B");
      eff[i]->GetYaxis()->SetTitle("efficiency");
      eff[i]->GetYaxis()->SetTitleOffset(1.2);
   }

   // Re-weighting efficiency
   if ( mc == 1 && p->use_rew == 1 ) {
      for ( int i = 0; i < 2; i++ ) {
         double w = ReWeightEtaEff(p->date);
         auto& hst = eff[i];
         int nx = hst->GetNbinsX();
         for ( int ix = 0; ix <= nx+1; ++ix ) {
//             double par = hst->GetXaxis()->GetBinCenter(ix);
            double d  = hst->GetBinContent(ix);
            double ed = hst->GetBinError(ix);

            hst->SetBinContent(ix,d*w);
            hst->SetBinError(ix,ed*w);
      }
   }
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
      string name = string("r_") + string( effdat[i]->GetName() );
      rat[i]=(TH1D*)effdat[i]->Clone( name.c_str() );
      rat[i]->Divide(effmc[i]);
      rat[i]->GetYaxis()->SetTitle("data/MC");
   }
}

//-------------------------------------------------------------------------
void plot_pict(Params* p, string pdf) {
//-------------------------------------------------------------------------
   vector<TH1D*> eff_d, eff_mc, rat0;
   get_eff(p, 0, eff_d); // data
   get_eff(p, 1, eff_mc);// MC
   get_ratio( eff_d, eff_mc, rat0 );

   int Nh = eff_d.size();

   // Attributes of draw
   vector<TLine*> line1(2,nullptr);
   for ( int i = 0; i < 2; i++ ) {
      line1[i] = new TLine(
            rat0[i]->GetXaxis()->GetXmin(), 1.,
            rat0[i]->GetXaxis()->GetXmax(), 1.  );
   }
   for (auto& l : line1 ) {
      l->SetLineColor(kGreen+2);
      l->SetLineWidth(2);
      l->SetLineStyle(kSolid); // or kDashed
   }

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
      setBlue(eff_d[i]);
      setRed(eff_mc[i]);
   }

   TLegend* leg1 = new TLegend(0.35,0.20,0.65,0.40);
   leg1->AddEntry(eff_d[0], Form("#color[%i]{data %i}",
                  eff_d[0]->GetLineColor(),p->date),"LP");
   leg1->AddEntry(eff_mc[0], Form("#color[%i]{MC %i}",
                  eff_mc[0]->GetLineColor(),p->date),"LP");


   // for fit
   TF1* pl0 = (TF1*)gROOT->GetFunction("pol0")->Clone();
   pl0->SetLineColor(kRed);

   double xmid = 1.3;
   auto Fpl1 = [xmid](const double* xx, const double* p) -> double {
      double x = xx[0] - xmid;
      double pl1 = p[0] + x*p[1];
      return pl1;
   };
   TF1* pl1 = new TF1("pl1", Fpl1, 1.2, 1.7,2);
   pl1->SetLineColor(kRed);


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
   double eff_min = 0., eff_max = 1.;
//    double rat_min = 0.9, rat_max = 1.1;
   double rat_min = 0.8, rat_max = 1.2;
//    double rat_min = 0.7, rat_max = 1.3;
///////////////////////////////////////////////////////////////////

   // plot efficiencies
   for (int i = 0; i < Nh; ++i ) {
      c1->cd(i+1);
      eff_d[i]->SetAxisRange(eff_min,eff_max,"Y");
      eff_d[i]->Draw("E");
      eff_mc[i]->Draw("E SAME");
      leg1->Draw();
   }

   // plot ratio
   for (int i = 0; i < Nh; ++i ) {
      c1->cd(i+3);
      rat0[i]->SetAxisRange(rat_min,rat_max,"Y");
         rat0[i]->SetTitle("#epsilon(DATA)/#epsilon(MC)");
         rat0[i]->SetLineWidth(2);
         rat0[i]->Fit("pol0","Q");
//       line1[i%2]->Draw();
//       rat0[i]->Draw("SAME");
   }
   addPDF(gPad);

///////////////////////////////////////////////////////////////////

   if ( ispdf ) {
      c1->Print((pdf+"]").c_str()); // close pdf-file
   }
}

//-------------------------------------------------------------------------
void eta_eff() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);
   gStyle->SetLegendFont(62);
   gStyle->SetLegendTextSize(0.05);
//    gStyle->SetStatFontSize(0.07);
   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.35);
   gStyle->SetFitFormat(".3f");

   // date, slct (>0: g-eta | <0: phi-eta; 2: eta-eff | 1:g-eff), weights
   Params* par = new Params(2012,2,0);

//    string pdf = ""; // debug
   string pdf = string("Eff_")
      + ( (abs(par->slct) == 2 ) ? "eta_" : "ph_" )
      + to_string(par->date)
      + ( (par->slct > 0) ? "_gammaeta" : "_phieta" )
      + ((par->use_rew > 0) ? "RW" : "")
      + string(".pdf");

   plot_pict(par,pdf);

}
