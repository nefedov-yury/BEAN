// Study of the eta -> 2gamma reconstruction efficiency &
//                     single photon rec.efficiency
// -> Eff_[eta,ph]_[date]_[gamma,phi]eta.pdf

#include "ReWeightEtaEff.h"

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

   const double Br_eta_2gamma = 0.3941; // PDG
   const double Br_phi_KK = 0.492; // PDG

   //-----------------------------------------------------------------
   Params(int dat = 2012, int slc = 2, int rew = 0) {
   //-----------------------------------------------------------------
      date = dat;
      slct = slc;
      use_rew = rew;

      // names of root-files
//       string Dir("prod-13eff/");
      string Dir("prod-11/");
      fnames = {
         Dir+"data_09psip_all.root", Dir+"mcinc_09psip_all.root",
         Dir+"data_12psip_all.root", Dir+"mcinc_12psip_all.root",
         // signal
         Dir+"mcgammaeta_kkmc_09.root", Dir+"mcgammaeta_kkmc_12.root",
         Dir+"mcphieta2_kkmc_09.root",  Dir+"mcphieta2_kkmc_12.root"
         // signal prog-13eff
         // Dir+"mcgammaeta2_kkmc_09.root", Dir+"mcgammaeta2_kkmc_12.root",
         // Dir+"mcphieta2_kkmc_09.root",   Dir+"mcphieta2_kkmc_12.root"
         // eta -> 2gamma decay only prog-13eff
         // Dir+"mcgammaeta_kkmc_09.root", Dir+"mcgammaeta_kkmc_12.root"
         // Dir+"mcsig_kkmc_09.root",      Dir+"mcsig_kkmc_12.root"
      };

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
         Cph += TCut("Eg2>0.05&&Eg2<1.45");
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
   TFile* OpenFile(int mc = 0) { // 1: MC inc, 2: MC-signal
   //-----------------------------------------------------------------
      int idx = mc + ((date==2009) ? 0 : 2);
      if ( mc == 2 ) {
         if ( slct > 0 ) {
            idx = ((date==2009) ? 4 : 5); // gamma-eta signal
         } else {
            idx = ((date==2009) ? 6 : 7); // phi-eta signal
         }
      }
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

   //-----------------------------------------------------------------
   TTree* GetEffEta(int mc = 0) { // get tree
   //-----------------------------------------------------------------
      TFile* froot = this->OpenFile(mc);
      TTree* eff_eta = (TTree*)gDirectory->Get("eff_eta");
      if ( !eff_eta ) {
         cout << "can not find eff_eta in " << froot->GetName() << endl;
         exit(0);
      }
      return eff_eta;
   }

   //-----------------------------------------------------------------
   double W_g_eta() { // wights for MC-gamma-eta
   //-----------------------------------------------------------------
      // normalization on numbers in "official inclusive MC"
      double W = ((date==2012) ? (137258./5e5) : (35578./1.5e5));
      return W;
   }

   //-----------------------------------------------------------------
   double W_phi_eta() { // wights for MC-sig
   //-----------------------------------------------------------------
      // normalization on numbers in "official inclusive MC"
      double W = ((date==2012) ? (104950./5e5) : (27274./1.5e5));
      return W * Br_phi_KK;
   }

};
//--------------------------------------------------------------------

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
      X->SetTitleFont(42);
      X->SetTitleSize(0.05);
   }
   TAxis* Y = hst->GetYaxis();
   if ( Y ) {
      Y->SetLabelFont(62);
      Y->SetLabelSize(0.04);
      Y->SetTitleFont(42);
      Y->SetTitleSize(0.05);
   }
   TAxis* Z = hst->GetZaxis();
   if ( Z ) {
      Z->SetLabelFont(62);
      Z->SetLabelSize(0.04);
      Z->SetTitleFont(62);
      Z->SetTitleSize(0.04);
   }
}

//--------------------------------------------------------------------
void get_hst(Params* p, int mc, vector<TH1D*>& hst, int sigbg = 0) {
//--------------------------------------------------------------------
// get histograms: mc=0,1 / data,MC

   TTree* eff_eta = p -> GetEffEta(mc);

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
   if ( mc > 0 && sigbg == 1 )          { // signal only
      Crec0 += p -> Cmcsig;
      Crec1 += p -> Cmcsig;
   } else if ( mc > 0 && sigbg == 2 )   { // bg only
      Crec0 += !(p -> Cmcsig);
      Crec1 += !(p -> Cmcsig);
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
      hst[0] = new TH1D( hns[0].c_str(), "", 14,0.05,1.45);
      hst[1] = new TH1D( hns[1].c_str(), "", 14,0.05,1.45);
      hst[2] = new TH1D( hns[2].c_str(), "", Cbins.size()-1, Cbins.data() );
      hst[3] = new TH1D( hns[3].c_str(), "", Cbins.size()-1, Cbins.data() );
   } else if ( p->slct == -1 ) {  // phi-eta => one photon efficiency
      if ( p->date == 2012 ) {
         hst[0] = new TH1D( hns[0].c_str(),"",14,0.05,1.45);
         hst[1] = new TH1D( hns[1].c_str(),"",14,0.05,1.45);
      } else {
         hst[0] = new TH1D( hns[0].c_str(),"",7,0.05,1.45);
         hst[1] = new TH1D( hns[1].c_str(),"",7,0.05,1.45);
      }
      hst[2] = new TH1D( hns[2].c_str(), "", Cbins.size()-1, Cbins.data() );
      hst[3] = new TH1D( hns[3].c_str(), "", Cbins.size()-1, Cbins.data() );
   } else if ( p->slct == 2 ) { // g-eta   => eta eff
      hst[0] = new TH1D( hns[0].c_str(), "", 8,1.3,1.7);
      hst[1] = new TH1D( hns[1].c_str(), "", 8,1.3,1.7);
      hst[2] = new TH1D( hns[2].c_str(), "", 9,-0.9,0.9);
      hst[3] = new TH1D( hns[3].c_str(), "", 9,-0.9,0.9);
   } else if ( p->slct == -2 ) { // phi-eta => eta eff
      hst[0] = new TH1D( hns[0].c_str(), "", 7,1.15,1.5);
      hst[1] = new TH1D( hns[1].c_str(), "", 7,1.15,1.5);
      hst[2] = new TH1D( hns[2].c_str(), "", 9,-0.9,0.9);
      hst[3] = new TH1D( hns[3].c_str(), "", 9,-0.9,0.9);
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

//--------------------------------------------------------------------
void get_eff_data(Params* p, vector<TH1D*>& eff ) {
//--------------------------------------------------------------------
// Get histograms and calculate efficiencies

   vector<TH1D*> hst;
   get_hst(p, 0, hst);

   TGraphAsymmErrors* tmp = nullptr;
//    tmp = new TGraphAsymmErrors();

   eff.clear();
   eff.resize(2,nullptr);
   for( int i = 0; i < 2; i++ ) {
      int ih = 2*i;
      string name = string("eff_") + string(hst[ih]->GetName());
      eff[i]=(TH1D*)hst[ih]->Clone( name.c_str() );
      hst[ih]->Add(hst[ih+1]);
      if ( !tmp ) {
         eff[i]->Divide(hst[ih+1],hst[ih],1.,1.,"B");
         int N = eff[i] -> GetNbinsX();
         for ( int j = 1; j <= N; ++j ) {
            int m = int(hst[ih] -> GetBinContent(j));
            if ( int(hst[ih+1] -> GetBinContent(j)) == m ) {
               double le = pow((1-0.683)/2,1./double(m));
               eff[i] -> SetBinError(j,1-le);
               cout << " ++ m= " << m << " le= " << le
                  << " err= " << 1-le << endl;
            }
         }
      } else {
         tmp -> Divide(hst[ih+1],hst[ih],"cp v");
         int N = tmp -> GetN();
         for(int j = 0; j < N; ++j) {
            double x,y;
            tmp -> GetPoint(j,x,y);
            double dy = tmp -> GetErrorY(j); // sym error?
            eff[i] -> SetBinContent(j+1,y);
            eff[i] -> SetBinError(j+1,dy);
         }
      }
   }

   delete tmp;
}

//--------------------------------------------------------------------
void get_eff_mc(Params* p, vector<TH1D*>& eff ) {
//--------------------------------------------------------------------
// Get signal from MC-"signal" and bg from inclusive MC
// and calculate efficiencies

   vector<TH1D*> hst, hstBG;
   get_hst(p, 2, hst, 1);   // MC-"signal"
   get_hst(p, 1, hstBG, 2); // MC bg only

   // sum of signal and bg.
   double W = 1;
   if ( p -> slct > 0 ) {       // gamma-eta
      W = p -> W_g_eta();
   } else {                     // phi-eta
      W = p -> W_phi_eta();
   }
   for ( int i = 0; i < hst.size(); i++ ) {
      hst[i] -> Add( hst[i], hstBG[i], W, 1.);
   }

   eff.clear();
   eff.resize(2,nullptr);
   for( int i = 0; i < 2; i++ ) {
      int ih = 2*i;
      string name = string("eff_") + string(hst[ih]->GetName());
      eff[i]=(TH1D*)hst[ih]->Clone( name.c_str() );
      eff[i]->Add(hst[ih+1]);
      eff[i]->Divide(hst[ih+1],eff[i],1.,1.,"B");
   }

   // Re-weighting efficiency
   if ( p->use_rew == 1 ) {
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

//--------------------------------------------------------------------
void get_ratio( const vector<TH1D*>& effdat,
                const vector<TH1D*>& effmc,
                                vector<TH1D*>& rat ) {
//--------------------------------------------------------------------
// calculate ratio of efficiencies DATA/MC

   int Nh = effdat.size();
   rat.clear();
   rat.resize(Nh,nullptr);
   for ( int i = 0; i < Nh; ++i ) {
      string name = string("r_") + string( effdat[i]->GetName() );
      rat[i]=(TH1D*)effdat[i]->Clone( name.c_str() );
      rat[i]->Divide(effmc[i]);
   }
}

//--------------------------------------------------------------------
void plot_pict_gamma_eta(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date,2,0); // date, eta_eff, no_rew
//    Params* par = new Params(date,1,0); // date, gamma_eff, no_rew

   vector<TH1D*> eff_d, eff_mc, rat0;
   get_eff_data(par, eff_d);
   get_eff_mc(par, eff_mc);
   get_ratio( eff_d, eff_mc, rat0 );

   int Nh = eff_d.size();

   // Attributes of draw
   auto setBlue = [](TH1* h) {
      h -> SetLineColor(kBlue+1);
      h -> SetMarkerColor(kBlue+2);
      h -> SetMarkerStyle(22);
      h -> SetMarkerSize(0.9);
   };
   auto setRed = [](TH1* h) {
      h -> SetLineColor(kRed+1);
      h -> SetMarkerColor(kRed+2);
      h -> SetMarkerStyle(23);
      h -> SetMarkerSize(0.9);
   };
   for (int i = 0; i < Nh; ++i ) {
      setBlue(eff_d[i]);
      setRed(eff_mc[i]);
   }

   // for fit
   TF1* pl0 = (TF1*)gROOT -> GetFunction("pol0")->Clone();
   pl0 -> SetLineColor(kRed);

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,400);
   c1 -> Divide(2,1);

///////////////////////////////////////////////////////////////////
   double eff_min = 0.6, eff_max = 1.0;
   double rat_min = 0.9, rat_max = 1.1;
///////////////////////////////////////////////////////////////////

//    TLegend* leg = new TLegend(0.11,0.70,0.30,0.89);
   TLegend* leg = new TLegend(0.11,0.71,0.50,0.89);
   leg -> SetHeader( Form("%i",date),"C");
   leg -> SetNColumns(2);
   leg -> AddEntry(eff_d[0],  "  data ", "LP");
   leg -> AddEntry(eff_mc[0], "  MC ",   "LP");

   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.25);
   gStyle->SetFitFormat(".3f");

   // plot efficiencies
   for (int i = 0; i < Nh; ++i ) {
      c1 -> cd(i+1);
      gPad -> SetGrid();
      eff_d[i] -> SetAxisRange(eff_min,eff_max,"Y");
      string title;
      if (par->slct == 2 ) {
         title = string( ((i==0) ? ";P, GeV/c" : ";cos(#Theta)") ) +
                     string(";#epsilon(#eta)");
      } else {
         title = string( ((i==0) ? ";E_{#gamma}, GeV" : ";cos(#Theta)") ) +
                     string(";#epsilon(#gamma)");
      }
      eff_d[i] -> SetTitle(title.c_str());
      SetHstFace(eff_d[i]);
      eff_d[i] -> GetXaxis() -> SetTitleOffset(0.9);
      eff_d[i] -> GetYaxis() -> SetTitleOffset(1.);
      eff_d[i] -> GetYaxis() -> SetNdivisions(1005);
      eff_d[i] -> Draw("E");
      eff_mc[i] -> Draw("E SAME");
      leg -> Draw();
   }

   gPad -> RedrawAxis();
   c1 -> Update();
   string pdf1;
   if (par->slct == 2 ) {
      pdf1 = string("eta_getaeff_") + to_string(date) + string(".pdf");
   } else {
      pdf1 = string("gamma_getaeff_") + to_string(date) + string(".pdf");
   }
   c1 -> Print( pdf1.c_str() );

   TCanvas* c2 = new TCanvas("c2","...",0,500,1200,400);
   c2 -> Divide(2,1);
   // plot ratio
   for (int i = 0; i < Nh; ++i ) {
      c2 -> cd(i+1);
      gPad -> SetGrid();
      rat0[i] -> SetAxisRange(rat_min,rat_max,"Y");
      string title;
      if (par->slct == 2 ) {
         title = string(Form("#eta, data/MC %i",date)) +
                     ( (i==0) ? ";P, GeV/c" : ";cos(#Theta)" ) +
                     string(";#epsilon(data) / #epsilon(MC)");
      } else {
         title = string(Form("#gamma, data/MC %i",date)) +
                     ( (i==0) ? ";E_{#gamma}, GeV" : ";cos(#Theta)" ) +
                     string(";#epsilon(data) / #epsilon(MC)");
      }
      rat0[i] -> SetTitle(title.c_str());
      SetHstFace(rat0[i]);
      rat0[i] -> GetXaxis() -> SetTitleOffset(0.9);
      rat0[i] -> GetYaxis() -> SetTitleOffset(1.);
      rat0[i] -> SetLineWidth(2);
      rat0[i] -> SetMarkerStyle(20);
      rat0[i] -> SetLineColor(kBlack);
      rat0[i] -> Fit( pl0, "" );
      rat0[i] -> Draw("SAME E0");
   }

   gPad -> RedrawAxis();
   c2 -> Update();
   string pdf2;
   if (par->slct == 2 ) {
      pdf2 = string("eta_getarat_") + to_string(date) + string(".pdf");
   } else {
      pdf2 = string("gamma_getarat_") + to_string(date) + string(".pdf");
   }
   c2 -> Print( pdf2.c_str() );
}

//--------------------------------------------------------------------
void plot_pict_phi_eta(int date) {
//--------------------------------------------------------------------
   Params* par = new Params(date,-2,0); // date, eta_eff in phieta, no_rew

   vector<TH1D*> eff_d, eff_mc, rat0;
   get_eff_data(par, eff_d);
   get_eff_mc(par, eff_mc);
   get_ratio( eff_d, eff_mc, rat0 );

   int Nh = eff_d.size();

   // Attributes of draw
   auto setBlue = [](TH1* h) {
      h -> SetLineColor(kBlue+1);
      h -> SetMarkerColor(kBlue+2);
      h -> SetMarkerStyle(22);
      h -> SetMarkerSize(0.9);
   };
   auto setRed = [](TH1* h) {
      h -> SetLineColor(kRed+1);
      h -> SetMarkerColor(kRed+2);
      h -> SetMarkerStyle(23);
      h -> SetMarkerSize(0.9);
   };
   for (int i = 0; i < Nh; ++i ) {
      setBlue(eff_d[i]);
      setRed(eff_mc[i]);
   }

   // for fit
   TF1* pl0 = (TF1*)gROOT -> GetFunction("pol0")->Clone();
   pl0 -> SetLineColor(kRed);

   // Draw:
   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,400);
   c1 -> Divide(2,1);

///////////////////////////////////////////////////////////////////
   double eff_min = 0.6, eff_max = 1.1;
//    double rat_min = 0.9, rat_max = 1.1;
   double rat_min = 0.8, rat_max = 1.2;
///////////////////////////////////////////////////////////////////

//    TLegend* leg = new TLegend(0.12,0.70,0.35,0.89);
   TLegend* leg = new TLegend(0.11,0.71,0.50,0.89);
   leg -> SetHeader( Form("%i",date),"C");
   leg -> SetNColumns(2);
   leg -> AddEntry(eff_d[0], " data ", "LP");
   leg -> AddEntry(eff_mc[0], " MC ", "LP");

   gStyle->SetStatX(0.89);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.25);
   gStyle->SetFitFormat(".3f");

   // plot efficiencies
   for (int i = 0; i < Nh; ++i ) {
      c1 -> cd(i+1);
      gPad -> SetGrid();
      eff_d[i] -> SetAxisRange(eff_min,eff_max,"Y");
      string title = string( ((i==0) ? ";P, GeV/c" : ";cos(#Theta)") ) +
                     string(";#epsilon(#eta)");
      eff_d[i] -> SetTitle(title.c_str());
      SetHstFace(eff_d[i]);
      eff_d[i] -> GetXaxis() -> SetTitleOffset(0.9);
      eff_d[i] -> GetYaxis() -> SetTitleOffset(1.);
//       eff_d[i] -> GetYaxis() -> SetNdivisions(1005);
      eff_d[i] -> Draw("E");
      eff_mc[i] -> Draw("E SAME");
      leg -> Draw();
   }

   gPad -> RedrawAxis();
   c1 -> Update();
   string pdf1 = string("eta_phietaeff_") + to_string(date) + string(".pdf");
   c1 -> Print( pdf1.c_str() );

   TCanvas* c2 = new TCanvas("c2","...",0,500,1200,400);
   c2 -> Divide(2,1);
   // plot ratio
   for (int i = 0; i < Nh; ++i ) {
      c2 -> cd(i+1);
      gPad -> SetGrid();
      rat0[i] -> SetAxisRange(rat_min,rat_max,"Y");
      string title = string(Form("#eta, data/MC %i",date)) +
                     ( (i==0) ? ";P, GeV/c" : ";cos(#Theta)" ) +
                     string(";#epsilon(data) / #epsilon(MC)");
      rat0[i] -> SetTitle(title.c_str());
      SetHstFace(rat0[i]);
      rat0[i] -> GetXaxis() -> SetTitleOffset(0.9);
      rat0[i] -> GetYaxis() -> SetTitleOffset(1.);
      rat0[i] -> SetLineWidth(2);
      rat0[i] -> SetMarkerStyle(20);
      rat0[i] -> SetLineColor(kBlack);
      rat0[i] -> Fit( pl0, "" );
   }

   gPad -> RedrawAxis();
   c2 -> Update();
   string pdf2 = string("eta_phietarat_") + to_string(date) + string(".pdf");
   c2 -> Print( pdf2.c_str() );
}

//--------------------------------------------------------------------
void eta_eff() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(42);
//    gStyle->SetStatFont(62);

   // -- J/Psi -> gamma eta, fig 57,58
//    plot_pict_gamma_eta(2009);
//    plot_pict_gamma_eta(2012);

   // -- J/Psi -> phi eta, fig 64,65
//    plot_pict_phi_eta(2009);
//    plot_pict_phi_eta(2012);
}
