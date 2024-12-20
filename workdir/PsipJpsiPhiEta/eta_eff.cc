// eta_eff.cc
// Study of the eta->2gamma reconstruction efficiency
// Maybe: single photon reconstruction efficiency: need tests
// -> etaeff_{geta|phieta}_{efficiency or ratio variable}_{YEAR}.pdf

#include "RewEtaEff.hpp"    // RewEtaEff() func
#include "RewTrkPiK.hpp"    // RewTrkPi(), RewTrk_K() functions

// {{{1 Common parameters: Params
//--------------------------------------------------------------------
struct Params {
   Params(int dat, int slc, int rew);

   TFile* OpenFile(int mc);
   TTree* GetEffEta(int mc);
   double W_g_eta();    // wights for MC-gamma-eta2
   double W_phi_eta();  // wights for MC-phi-eta2
   const char* Sdate() { return Form("%i",date); }

   // name of folder with root files
   const string Dir = "prod_v709n4/Eff/";
   string datafile;
   string mcincfile;
   string mcsigf1;  // MC for gamma eta
   string mcsigf2;  // MC for phi eta

   int date;
   int slct;    // >0: g-eta; <0: phi-eta
                // 1 - photon eff; 2 - eta eff;
   int use_rew; // 0 - no weights;
                // 1 - calculate weights

   TCut Cmcsig; // select mc-signal
   TCut Cbg;    // background suppression
   TCut Cph;    // selection photons
   TCut Ceta;   // limits for fitting the ratio
   TCut Cfnd;   // predicted eta found
   TCut Cfnd_g; // predicted gamma found

   const double meta   = 0.54786; // 547.862   +/- 0.017   MeV
   const double Br_eta_2gamma = 0.3936; // PDG 2023
   const double Br_phi_KK = 0.491; // PDG 2023
};

// {{{2 > ctor
//--------------------------------------------------------------------
Params::Params(int dat, int slc = 2, int rew = 0)
//--------------------------------------------------------------------
{
   date = dat;
   slct = slc;
   use_rew = rew;

   // set the names:
   datafile  = string( Form("data_%02ipsip_all.root",date%100) );
   mcincfile = string( Form("mcinc_%02ipsip_all.root",date%100) );
   mcsigf1 = string( Form("mcgammaeta2_kkmc_%02i.root",date%100) );
   mcsigf2 = string( Form("mcphieta2_kkmc_%02i.root",date%100) );

   // select mc-signal
   if ( slct > 0 ) {    // gamma-eta
      Cmcsig += TCut("decj==22");
   } else {             // phi-eta
      Cmcsig += TCut("decj==68");
   }

   // cuts against the background
   if ( slct > 0 ) {         // gamma-eta
      Cbg += TCut("abs(Cg0)<0.8");
      Cbg += TCut("abs(Cg1)<0.8");
      Cbg += TCut("m2fr<0.002");
      Cbg += TCut("Eg1>0.25&&Eg1<1.7");
   } else {                  // phi-eta
      Cbg += TCut("1.02<Mkk2&&Mkk2<1.06");
      Cbg += TCut("m2fr<0.001");
   }

   // background suppression
   Cph += TCut("abs(Cg1)<0.8||(abs(Cg1)>0.85&&abs(Cg1)<0.92)");
   Cph += TCut("abs(Cg2)<0.8||(abs(Cg2)>0.85&&abs(Cg2)<0.92)");
   if ( slct > 0 ) {         // gamma-eta
      Cph += TCut("Eg2>0.1&&Eg2<1.4");
   } else {                  // phi-eta
      Cph += TCut("Eg2>0.05&&Eg2<1.45");
   }

   // selection photons
   if ( slct > 0 ) {         // gamma-eta
      Ceta += TCut("Peta>1.3&&Peta<1.7");
   } else {                  // phi-eta
      Ceta += TCut("Peta>1.15&&Peta<1.5");
   }
   Ceta += TCut("abs(Ceta)<0.9");
   // doubtful:
   // Ceta += TCut("abs(Cg2)<0.8||(abs(Cg2)>0.85&&abs(Cg2)<0.92)");

   // predicted eta is found
   Cfnd += TCut("dTh<8.&&0.8<rE&&rE<1.4");
   Cfnd += TCut( Form("abs(mggf-%.6f)<0.024",meta) );

   // predicted gamma is found
   Cfnd_g += TCut("dTh<8.&&0.8<rE&&rE<1.4");
}

// {{{2 > OpenFile(mc = 0: data, 1: MC inc, 2: MC-signal)
//--------------------------------------------------------------------
TFile* Params::OpenFile(int mc)
//--------------------------------------------------------------------
{
   string dfname = Dir;
   if ( mc == 0 ) {
      dfname += datafile;
   } else if ( mc == 1 ) {
      dfname += mcincfile;
   } else if ( mc == 2 ) {
      if ( slct > 0 ) {
         dfname += mcsigf1; // gamma-eta-2
      } else {
         dfname += mcsigf2; // phi-eta-2
      }
   }
   // cout " OpenFile: " << dfname << endl;
   TFile* froot = TFile::Open(dfname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << dfname << endl;
      exit(EXIT_FAILURE);
   }
   if ( slct > 0 ) {
      froot->cd("PsipJpsiGammaEta");
   } else {
      froot->cd("PsipJpsiPhiEta");
   }
   return froot;
}

// {{{2 > GetEffEta : root-tree for eta
//--------------------------------------------------------------------
TTree* Params::GetEffEta(int mc = 0)
//--------------------------------------------------------------------
{
   TFile* froot = this->OpenFile(mc);
   TTree* eff_eta = (TTree*)gDirectory->Get("eff_eta");
   if ( !eff_eta ) {
      cout << "can not find eff_eta in " << froot->GetName() << endl;
      exit(EXIT_FAILURE);
   }
   return eff_eta;
}

// {{{2 > W_g_eta() & W_phi_eta()
//--------------------------------------------------------------------
double Params::W_g_eta() // wights for MC gammaeta2
//--------------------------------------------------------------------
{
   // normalization on numbers in "official inclusive MC"
   // decay code for J/Psi -> gamma eta is 22;
   // BOSS-664: ((date==2012) ? (137258./5e5) : (35578./1.5e5))
   double W = 1;
   switch (date) {
      case 2009:
         W =  41553./150e3; // 0.28
         break;
      case 2012:
         W = 130865./500e3; // 0.26
         break;
      case 2021:
         W = 883676./2.5e6; // 0.35
         break;
   }
   return W;
}

//--------------------------------------------------------------------
double Params::W_phi_eta() // wights for MC phieta2
//--------------------------------------------------------------------
{
   // normalization on numbers in "official inclusive MC"
   // decay code for J/Psi -> phi eta is 68;
   // BOSS-664: ((date==2012) ? (104950./5e5) : (27274./1.5e5));
   // we have to use events with J/Psi -> phi eta, phi->K+K-
   // stored in histogram mc_MKK.
   double W = 1;
   switch (date) {
      case 2009:
         // W =  27547./150e3; // decj0 = 68
         W =  13519./150e3; // mc_MKK
         break;
      case 2012:
         // W =  89059./500e3; // decj0 = 68
         W =  44287./500e3; // mc_MKK
         break;
      case 2021:
         // W = 598360./2.5e6; // decj0 = 68
         W = 297126./2.5e6; // mc_MKK
         break;
   }
   // W *= Br_phi_KK; // for decj0 = 68
   return W;
}

// {{{1 helper functions and constants
//--------------------------------------------------------------------
constexpr double SQ(double x)
//--------------------------------------------------------------------
{
   return x*x;
}

//--------------------------------------------------------------------
void SetHstFace(TH1* hst)
//--------------------------------------------------------------------
{
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

// {{{1 OLD: fill histograms: mc=0,1 / data,MC
//--------------------------------------------------------------------
void get_hst(Params* p, int mc, vector<TH1D*>& hst, int sigbg = 0)
//--------------------------------------------------------------------
{
   TTree* eff_eta = p->GetEffEta(mc);

   TCut Crec0 = p->Cbg;
   TCut Crec1 = p->Cbg;
   if ( abs(p->slct) == 1 ) {      // photon
      Crec0 += !(p->Cfnd_g);
      Crec1 += p->Cfnd_g;
      Crec0 += p->Cph;
      Crec1 += p->Cph;
   } else {                     // eta
      Crec0 += !(p->Cfnd);
      Crec1 += p->Cfnd;
      Crec0 += p->Ceta;
      Crec1 += p->Ceta;
   }
   if ( mc > 0 && sigbg == 1 )          { // signal only
      Crec0 += p->Cmcsig;
      Crec1 += p->Cmcsig;
   } else if ( mc > 0 && sigbg == 2 )   { // bg only
      Crec0 += !(p->Cmcsig);
      Crec1 += !(p->Cmcsig);
   }
   // cout << " Crec0: " << Crec0.GetTitle() << endl;
   // cout << " Crec1: " << Crec1.GetTitle() << endl;

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
   // for( auto c : Cbins ) {
      // cout << c << " : ";
   // }
   // cout << endl;

   if ( p->slct == 1 ) {  // gamma-eta => one photon efficiency
      hst[0] = new TH1D( hns[0].c_str(), "", 14,0.05,1.45);
      hst[1] = new TH1D( hns[1].c_str(), "", 14,0.05,1.45);
      hst[2] = new TH1D( hns[2].c_str(), "",
            Cbins.size()-1, Cbins.data() );
      hst[3] = new TH1D( hns[3].c_str(), "",
            Cbins.size()-1, Cbins.data() );
   } else if ( p->slct == -1 ) {  // phi-eta => one photon efficiency
      if ( p->date == 2012 ) {
         hst[0] = new TH1D( hns[0].c_str(),"",14,0.05,1.45);
         hst[1] = new TH1D( hns[1].c_str(),"",14,0.05,1.45);
      } else {
         hst[0] = new TH1D( hns[0].c_str(),"",7,0.05,1.45);
         hst[1] = new TH1D( hns[1].c_str(),"",7,0.05,1.45);
      }
      hst[2] = new TH1D( hns[2].c_str(), "",
            Cbins.size()-1, Cbins.data() );
      hst[3] = new TH1D( hns[3].c_str(), "",
            Cbins.size()-1, Cbins.data() );
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

// {{{1 Fill gamma-eta histograms: mc=0,1 / data,MC
//--------------------------------------------------------------------
void FillGEtaHst(Params* p, int mc, vector<TH1D*>& hst, int sigbg=0)
//--------------------------------------------------------------------
{
   // cout << __func__ << " mc= " << mc << " sigbg= " << sigbg<<endl;
   if ( p->slct != 2 ) { // gamma-eta => eta eff
      cout << __func__ << " function, incorrect using: p->slct="
         << p->slct << " Stop!" << endl;
      exit(1);
   }

   TTree* eff_eta = p->GetEffEta(mc);

   double Ptp,Ptm,Mrec; // Pt(pi+), Pt(pi-), Mrec(pi+pi-)
   double Eg0,Cg0;      // E,cos(Theta) of gamma in J/Psi -> g0 eta
   double Peta,Ceta;    // P,cos(Theta) of eta
   double Eg1,Cg1;      // E,cos(Theta) of rec gamma in eta -> g1 g2
   double Eg2,Cg2;      // E,cos(Theta) of predicted gamma (g2)
   // double Egr,Cgr;      // E and cos(Theta) of found gamma (gr)
   double rE,dTh;       // relations btw predicted(g2) and found(gr)
                        // * variables used in selection:
   // double M2gr;         // Mrec^2(pi+pi-g0)
   // double M2gg;         // Minv^2(g1 g2)
   double m2fr;         // Mrec^2(pi+pi-g0g1) or g2 missing mass^2
   double mggf;         // Minv(g1,gr) after 4C kinematic constraints
   // double ch2;          // chi^2 of 4C kinematic constraints
                        // * MC-variables
   int    decj;         // MC: decay codes of J/Psi
   // int    dgam;         // MC: 1 if g0 found correctly
   // int    deta;         // MC: 1 if eta decays into two gammas

   // Set branch addresses.
   eff_eta->SetBranchAddress("Ptp", &Ptp );
   eff_eta->SetBranchAddress("Ptm", &Ptm );
   eff_eta->SetBranchAddress("Mrec",&Mrec);
   eff_eta->SetBranchAddress("Eg0",  &Eg0);
   eff_eta->SetBranchAddress("Cg0",  &Cg0);
   eff_eta->SetBranchAddress("Peta",&Peta);
   eff_eta->SetBranchAddress("Ceta",&Ceta);
   eff_eta->SetBranchAddress("Eg1",  &Eg1);
   eff_eta->SetBranchAddress("Cg1",  &Cg1);
   eff_eta->SetBranchAddress("Eg2",  &Eg2);
   eff_eta->SetBranchAddress("Cg2",  &Cg2);

   eff_eta->SetBranchAddress("rE",   &rE );
   eff_eta->SetBranchAddress("dTh",  &dTh);

   eff_eta->SetBranchAddress("m2fr",&m2fr);
   eff_eta->SetBranchAddress("mggf",&mggf);

   eff_eta->SetBranchAddress("decj",&decj);

   // see Params::Params() for the selection parameters
   auto c_geta = [](double Cg0, double Cg1, double m2fr, double Eg1) {
      return ( (fabs(Cg0)<0.8 && fabs(Cg1)<0.8) && m2fr<0.002
         && (Eg1>0.25&&Eg1<1.7) );
   };

   auto c_eta = [](double Peta, double Ceta)->bool{
      // return (Peta>1.3&&Peta<1.7);
      return (Peta>1.3&&Peta<1.7) && (fabs(Ceta)<0.9);
   };
   // doubtful:
   auto c_g2 = [](double Cg2)->bool{
      return (fabs(Cg2)<0.8) || (fabs(Cg2)>0.85&&fabs(Cg2)<0.92);
   };

   auto& meta = p->meta;
   auto c_fnd = [meta](double dTh, double rE, double mggf)->bool{
      return (dTh<8.&&(0.8<rE&&rE<1.4) && fabs(mggf-meta)<0.024);
   };

   hst.clear();
   hst.resize(4,nullptr);
   auto hn = Form("%s_%i",((mc==0) ? "data" : "mc"),p->date);
   hst[0] = new TH1D( Form("%s_P",hn),  "", 8,1.3,1.7);
   hst[1] = new TH1D( Form("%s_P1",hn), "", 8,1.3,1.7);
   hst[2] = new TH1D( Form("%s_C",hn),  "", 9,-0.9,0.9);
   hst[3] = new TH1D( Form("%s_C1",hn), "", 9,-0.9,0.9);
   for ( auto& h : hst ) {
      h->Sumw2(true);
   }

   Long64_t nentries = eff_eta->GetEntries();
   for ( Long64_t i = 0; i < nentries; ++i ) {
      eff_eta->GetEntry(i);

      if ( !c_geta(Cg0,Cg1,m2fr,Eg1) ) { continue; }
      if ( !c_eta(Peta,Ceta) ) { continue; }

      // doubtful:
      // if ( !c_g2(Cg2) ) { continue; }

      double W = 1; // weights: use only in numerator of efficiency
      if ( mc > 0 ) {
         if ( sigbg == 1 && decj != 22 )         { // signal only
            continue;
         } else if ( sigbg == 2 &&  decj == 22 ) { // bg only
            continue;
         }
         // W *= RewTrkPi(p->date, Ptp, 1);
         // W *= RewTrkPi(p->date, Ptm, -1);
      }

      hst[0]->Fill(Peta);
      hst[2]->Fill(Ceta);
      if ( c_fnd(dTh,rE,mggf) ) { // found eta
         hst[1]->Fill(Peta,W);
         hst[3]->Fill(Ceta,W);
      }
   }
}

// {{{1 Fill phi-eta histograms: mc=0,1 / data,MC
//--------------------------------------------------------------------
void FillPhiEtaHst(Params* p, int mc, vector<TH1D*>& hst, int sigbg=0)
//--------------------------------------------------------------------
{
   // cout << __func__ << " mc= " << mc << " sigbg= " << sigbg<<endl;
   if ( p->slct != -2 ) { // phi-eta => eta eff
      cout << __func__ << " function, incorrect using: p->slct="
         << p->slct << " Stop!" << endl;
      exit(1);
   }

   TTree* eff_eta = p->GetEffEta(mc);

   // Declaring only the leaf types we used here
   // double Ptp,Ptm,Mrec; // Pt(pi+), Pt(pi-), Mrec(pi+pi-)
   double Ptkp,Ptkm;    // Pt(K+), Pt(K-)
   double Peta,Ceta;    // P,cos(Theta) of eta
   // double Eg1,Cg1;      // E,cos(Theta) of rec gamma in eta -> g1 g2
   double Eg2,Cg2;      // E,cos(Theta) of predicted gamma (g2)
   // double Egr,Cgr;      // E and cos(Theta) of found gamma (gr)
   double rE,dTh;       // relations btw predicted(g2) and found(gr)
                        // * variables used in selection:
   double Mkk2;         // Minv2(K+K-)
   // double M2gg;         // ?Minv^2(g1 g2)
   double m2fr;         // Mrec^2(pi+pi-g0g1) or g2 missing mass^2
   double mggf;         // Minv(g1,gr) after 4C kinematic constraints
   // double ch2;          // chi^2 of 4C kinematic constraints
                        // * MC-variables
   int    decj;         // MC: decay codes of J/Psi
   // int    deta;         // MC: 1 if eta decays into two gammas
   // int    dphi;         // MC: 1 if phi decays into K+K-

   // Set branch addresses.
   eff_eta->SetBranchAddress("Ptkp",&Ptkp);
   eff_eta->SetBranchAddress("Ptkm",&Ptkm);
   eff_eta->SetBranchAddress("Peta",&Peta);
   eff_eta->SetBranchAddress("Ceta",&Ceta);
   // eff_eta->SetBranchAddress("Eg1",  &Eg1);
   // eff_eta->SetBranchAddress("Cg1",  &Cg1);
   eff_eta->SetBranchAddress("Eg2",  &Eg2);
   eff_eta->SetBranchAddress("Cg2",  &Cg2);

   eff_eta->SetBranchAddress("rE",  &rE);
   eff_eta->SetBranchAddress("dTh", &dTh);

   eff_eta->SetBranchAddress("Mkk2",&Mkk2);
   eff_eta->SetBranchAddress("m2fr",&m2fr);
   eff_eta->SetBranchAddress("mggf",&mggf);

   eff_eta->SetBranchAddress("decj",&decj);

   // see Params::Params() for the selection parameters
   auto c_phieta = [](double Mkk2, double m2fr)->bool{
      return (1.02<Mkk2&&Mkk2<1.06) && (m2fr<0.001);
   };

   auto c_eta = [](double Peta, double Ceta)->bool{
      // return (Peta>1.15&&Peta<1.5);
      return (Peta>1.15&&Peta<1.5) && (fabs(Ceta)<0.9);
   };
   // doubtful:
   auto c_g2 = [](double Cg2)->bool{
      return (fabs(Cg2)<0.8) || (fabs(Cg2)>0.85&&fabs(Cg2)<0.92);
   };

   auto& meta = p->meta;
   auto c_fnd = [meta](double dTh, double rE, double mggf)->bool{
      return (dTh<8.&&(0.8<rE&&rE<1.4) && fabs(mggf-meta)<0.024);
   };

   hst.clear();
   hst.resize(4,nullptr);
   auto hn = Form("%s_%i",((mc==0) ? "data" : "mc"),p->date);
   hst[0] = new TH1D( Form("%s_P",hn),  "", 7,1.15,1.5);
   hst[1] = new TH1D( Form("%s_P1",hn), "", 7,1.15,1.5);
   hst[2] = new TH1D( Form("%s_C",hn),  "", 9,-0.9,0.9);
   hst[3] = new TH1D( Form("%s_C1",hn), "", 9,-0.9,0.9);
   for ( auto& h : hst ) {
      h->Sumw2(true);
   }

   Long64_t nentries = eff_eta->GetEntries();
   for ( Long64_t i = 0; i < nentries; ++i ) {
      eff_eta->GetEntry(i);

      if ( !c_phieta(Mkk2,m2fr) ) { continue; }

      // if ( !(Eg1>0.25&&Eg1<1.7) ) { continue; } // == geta
      if ( !c_eta(Peta,Ceta) ) { continue; }

      // doubtful:
      // if ( !c_g2(Cg2) ) { continue; }

      double W = 1; // weights: use only in numerator of efficiency
      if ( mc > 0 ) {
         if ( sigbg == 1 && decj != 68 )         { // signal only
            continue;
         } else if ( sigbg == 2 &&  decj == 68 ) { // bg only
            continue;
         }

         // W *= RewTrk_K(p->date, Ptkp, 1);
         // W *= RewTrk_K(p->date, Ptkm, -1);
      }

      hst[0]->Fill(Peta);
      hst[2]->Fill(Ceta);
      if ( c_fnd(dTh,rE,mggf) ) { // found eta
         hst[1]->Fill(Peta,W);
         hst[3]->Fill(Ceta,W);
      }
   }
}

// {{{1 Get histograms and calculate efficiencies: data and MC
//--------------------------------------------------------------------
void get_eff_data(Params* p, vector<TH1D*>& eff )
//--------------------------------------------------------------------
{
   vector<TH1D*> hst;
   if ( p->slct == +2 ) { // gamma-eta => eta eff
      FillGEtaHst(p, 0, hst);
   } else if ( p->slct == -2 ) { // phi-eta => eta eff
      FillPhiEtaHst(p, 0, hst);
   } else {
      get_hst(p, 0, hst);
      hst[0]->Add(hst[1]);
      hst[2]->Add(hst[3]);
   }

   TGraphAsymmErrors* tmp = nullptr;
   // tmp = new TGraphAsymmErrors();

   eff.clear();
   eff.resize(2,nullptr);
   for( int i = 0; i < 2; i++ ) {
      int ih = 2*i;
      string name = string("eff_") + string(hst[ih]->GetName());
      eff[i]=(TH1D*)hst[ih]->Clone( name.c_str() );
      if ( !tmp ) {
         eff[i]->Divide(hst[ih+1],hst[ih],1.,1.,"B");
         int N = eff[i]->GetNbinsX();
         for ( int j = 1; j <= N; ++j ) {
            int m = int(hst[ih]->GetBinContent(j));
            if ( int(hst[ih+1]->GetBinContent(j)) == m ) {
               double le = pow((1-0.683)/2,1./double(m));
               eff[i]->SetBinError(j,1-le);
               cout << " ++ m= " << m << " le= " << le
                  << " err= " << 1-le << endl;
            }
         }
      } else {
         tmp->Divide(hst[ih+1],hst[ih],"cp v");
         int N = tmp->GetN();
         for(int j = 0; j < N; ++j) {
            double x,y;
            tmp->GetPoint(j,x,y);
            double dy = tmp->GetErrorY(j); // sym error?
            eff[i]->SetBinContent(j+1,y);
            eff[i]->SetBinError(j+1,dy);
         }
      }
   }

   delete tmp;
}

//--------------------------------------------------------------------
void get_eff_mc(Params* p, vector<TH1D*>& eff )
//--------------------------------------------------------------------
{
   // take signal from MC-"signal" and bg from inclusive MC
   vector<TH1D*> hst, hstBG;
   if ( p->slct == +2 ) {            // gamma-eta => eta eff
      FillGEtaHst(p, 2, hst, 1);        // MC signal
      FillGEtaHst(p, 1, hstBG, 2);      // MC bg only
   } else if ( p->slct == -2 ) {     // phi-eta => eta eff
      FillPhiEtaHst(p, 2, hst, 1);      // MC signal
      FillPhiEtaHst(p, 1, hstBG, 2);    // MC bg only
   } else {
      get_hst(p, 2, hst, 1);   // MC signal
      get_hst(p, 1, hstBG, 2); // MC bg only
      hst[0]->Add(hst[1]);
      hst[2]->Add(hst[3]);
   }

   // sum of signal and background
   double W = 1;
   if ( p->slct > 0 ) {       // gamma-eta
      W = p->W_g_eta();
   } else {                   // phi-eta
      W = p->W_phi_eta();
   }
   for ( int i = 0; i < hst.size(); i++ ) {
      hst[i]->Add( hst[i], hstBG[i], W, 1.);
   }

   eff.clear();
   eff.resize(2,nullptr);
   for( int i = 0; i < 2; i++ ) {
      int ih = 2*i;
      string name = string("eff_") + string(hst[ih]->GetName());
      eff[i]=(TH1D*)hst[ih]->Clone( name.c_str() );
      eff[i]->Divide(hst[ih+1],hst[ih],1.,1.,"B");
   }

   // Re-weighting efficiency
   if ( p->use_rew == 1 ) {
      for ( int i = 0; i < 2; i++ ) {
         double w = 1.; // DEBUG
         // double w = RewEtaEff( p->date ); // Note: date is not used
         auto& hst = eff[i];
         int nx = hst->GetNbinsX();
         for ( int ix = 0; ix <= nx+1; ++ix ) {
            // double par = hst->GetXaxis()->GetBinCenter(ix);
            double d  = hst->GetBinContent(ix);
            double ed = hst->GetBinError(ix);

            hst->SetBinContent(ix,d*w);
            hst->SetBinError(ix,ed*w);
         }
      }
   }

}

// {{{1 Ratio of efficiencies DATA/MC
//--------------------------------------------------------------------
void get_ratio( const vector<TH1D*>& effdat,
      const vector<TH1D*>& effmc,
      vector<TH1D*>& rat )
//--------------------------------------------------------------------
{
   int Nh = effdat.size();
   rat.clear();
   rat.resize(Nh,nullptr);
   for ( int i = 0; i < Nh; ++i ) {
      string name = string("r_") + string( effdat[i]->GetName() );
      rat[i]=(TH1D*)effdat[i]->Clone( name.c_str() );
      rat[i]->Divide(effmc[i]);
   }
}

// {{{1 Plot gamma-eta
//--------------------------------------------------------------------
vector<double> plot_pict_gamma_eta(int date, int rew,
      int Cx=600, int Cy=600)
//--------------------------------------------------------------------
{
   Params* par = new Params(date,2,rew); // date, eta_eff, rew
   // Params* par = new Params(date,1,rew); // date, gamma_eff, no_rew

   vector<TH1D*> eff_d, eff_mc, rat0;
   get_eff_data(par, eff_d);
   get_eff_mc(par, eff_mc);
   get_ratio( eff_d, eff_mc, rat0 );

   int Nh = eff_d.size();

   // Attributes of draw
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

   // for fit
   TF1* pl0 = (TF1*)gROOT->GetFunction("pol0")->Clone();
   // TF1* pl0 = (TF1*)gROOT->GetFunction("pol1")->Clone(); // test
   pl0->SetLineColor(kRed);

   // Draw:
   ///////////////////////////////////////////////////////////////////
   double eff_min = 0.45, eff_max = 0.9;
   double rat_min = 0.9, rat_max = 1.1;
   ///////////////////////////////////////////////////////////////////
   string p_pdf = (par->slct == 2) ?  "etaeff_geta" : "Geff_geta";

   vector<TCanvas*> cc(4,nullptr);
   vector<TLegend*> leg(4,nullptr);
   for ( size_t i = 0; i < cc.size(); ++i ) {
      int x0 = (Cx/2)*(i%2);
      int y0 = (Cy/2)*(i/2);
      // printf("%zu -> x0=%i y0=%i\n",i,x0,y0);
      auto name = Form("c%zu_%i",1+i,date);
      cc[i] = new TCanvas( name,name, x0,y0,Cx,Cy);

      auto cdat = Form("%i",date);
      if ( i < 2 ) { // efficiency
         leg[i] = new TLegend(0.12,0.73,0.40,0.88);
         leg[i]->SetHeader(cdat,"C");
         leg[i]->SetNColumns(2);
         leg[i]->AddEntry(eff_d[0], " Data", "LP");
         leg[i]->AddEntry(eff_mc[0], " MC", "LP");
      } else { // retio
         leg[i] = new TLegend(0.14,0.81,0.42,0.89);
         if ( par->slct == 2 ) {
            leg[i]->SetHeader(cdat,"C");
         } else {
            leg[i]->SetHeader( Form("#gamma, %i",date), "C" );
         }
      }
   }

   // plot efficiencies
   for (int i = 0; i < Nh; ++i ) {
      auto& ci = cc[0+i];
      ci->cd();
      gPad->SetGrid();
      if ( par->slct == 2 ) {
         eff_d[i]->SetTitle( Form(";%s;#epsilon(#eta)",
                  ( i == 0 ? "P, GeV/c" : "cos(#Theta)")) );
      } else {
         eff_d[i]->SetTitle( Form(";%s;#epsilon(#gamma)",
                  ( i == 0 ? "E_{#gamma}, GeV" : "cos(#Theta)")) );
      }
      SetHstFace(eff_d[i]);
      eff_d[i]->SetAxisRange(eff_min,eff_max,"Y");
      // eff_d[i]->GetXaxis()->SetTitleOffset(1.);
      eff_d[i]->GetYaxis()->SetTitleOffset(1.);
      eff_d[i]->GetYaxis()->SetNdivisions(1005);
      eff_d[i]->Draw("E");
      leg[0+i]->Draw();
      eff_mc[i]->Draw("E SAME");

      gPad->RedrawAxis();
      ci->Update();
      string pp = p_pdf + ( ( i == 0 ) ? "_p" : "_cos" ) +
         ( ( rew == 1 ) ? "_w" : "" );
      ci->Print( Form("%s_%i.pdf",pp.c_str(),date) );
   }

   // plot ratio
   gStyle->SetStatX(0.91);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.25);
   gStyle->SetFitFormat(".4f");
   if ( date == 2009 ) {
      gStyle->SetFitFormat(".3f");
   }

   vector<double> fit_res;
   for (int i = 0; i < Nh; ++i ) {
      auto& ci = cc[2+i];
      ci->cd();
      gPad->SetLeftMargin(gPad->GetLeftMargin()+0.02);
      gPad->SetRightMargin(gPad->GetRightMargin()-0.02);
      gPad->SetGrid();
      if ( par->slct == 2 ) {
         rat0[i]->SetTitle( Form(";%s;#epsilon(data) / #epsilon(MC)",
                  ( i == 0 ? "P, GeV/c" : "cos(#Theta)")) );
      } else {
         rat0[i]->SetTitle( Form(";%s;#epsilon(data) / #epsilon(MC)",
                  ( i == 0 ? "E_{#gamma}, GeV" : "cos(#Theta)" )) );
      }
      SetHstFace(rat0[i]);
      rat0[i]->SetAxisRange(rat_min,rat_max,"Y");
      // rat0[i]->GetYaxis()->SetNdivisions(1004);
      // rat0[i]->GetXaxis()->SetTitleOffset(1.);
      rat0[i]->GetYaxis()->SetTitleOffset(1.2);
      rat0[i]->SetLineWidth(2);
      rat0[i]->SetMarkerStyle(20);
      rat0[i]->SetLineColor(kBlack);
      TFitResultPtr rs = rat0[i]->Fit( pl0, "S" );
      leg[2+i]->Draw();
      rat0[i]->Draw("SAME E0");

      double ch2 = pl0->GetChisquare();
      double ndf = pl0->GetNDF();
      double rr0 = pl0->GetParameter(0);
      double er0 = pl0->GetParError(0);
      fit_res.push_back(double(date));
      fit_res.push_back(ch2);
      fit_res.push_back(ndf);
      fit_res.push_back(rr0);
      fit_res.push_back(er0);

      gPad->RedrawAxis();
      ci->Update();
      string pp = p_pdf + "_rat" + ( ( i == 0 ) ? "_p" : "_cos" ) +
         ( ( rew == 1 ) ? "_w" : "" );
      ci->Print( Form("%s_%i.pdf",pp.c_str(),date) );
   }
   return fit_res;
}

// {{{1 Plot phi-eta
//--------------------------------------------------------------------
vector<double> plot_pict_phi_eta(int date, int rew,
      int Cx=600, int Cy=600)
//--------------------------------------------------------------------
{
   Params* par = new Params(date,-2,rew); // date, phi-eta, rew

   vector<TH1D*> eff_d, eff_mc, rat0;
   get_eff_data(par, eff_d);
   get_eff_mc(par, eff_mc);
   get_ratio( eff_d, eff_mc, rat0 );

   int Nh = eff_d.size();

   // Attributes of draw
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

   // for fit
   TF1* pl0 = (TF1*)gROOT->GetFunction("pol0")->Clone();
   pl0->SetLineColor(kRed);

   // Draw:
   ///////////////////////////////////////////////////////////////////
   double eff_min = 0.5, eff_max = 1.0;
   double rat_min = 0.9, rat_max = 1.1;
   if( date == 2009 ) {
      rat_min = 0.85, rat_max = 1.15;
   }
   ///////////////////////////////////////////////////////////////////
   string p_pdf = "etaeff_phieta";

   vector<TCanvas*> cc(4,nullptr);
   vector<TLegend*> leg(4,nullptr);
   for ( size_t i = 0; i < cc.size(); ++i ) {
      int x0 = (Cx/2)*(i%2);
      int y0 = (Cy/2)*(i/2);
      // printf("%zu -> x0=%i y0=%i\n",i,x0,y0);
      auto name = Form("c%zu_%i",1+i,date);
      cc[i] = new TCanvas( name,name, x0,y0,Cx,Cy);

      auto cdat = Form("%i",date);
      if ( i < 2 ) { // efficiency
         leg[i] = new TLegend(0.12,0.73,0.40,0.88);
         leg[i]->SetHeader( cdat,"C");
         leg[i]->SetNColumns(2);
         leg[i]->AddEntry(eff_d[0], "data", "LP");
         leg[i]->AddEntry(eff_mc[0], "MC", "LP");
      } else { // retio
         leg[i] = new TLegend(0.14,0.81,0.42,0.89);
         leg[i]->SetHeader( cdat,"C");
      }
   }
   // plot efficiencies
   for (int i = 0; i < Nh; ++i ) {
      auto& ci = cc[0+i];
      ci->cd();
      gPad->SetGrid();
      eff_d[i]->SetTitle( Form(";%s;#epsilon(#eta)",
               ( i == 0 ? "P, GeV/c" : "cos(#Theta)")) );
      SetHstFace(eff_d[i]);
      eff_d[i]->SetAxisRange(eff_min,eff_max,"Y");
      eff_d[i]->GetYaxis()->SetTitleOffset(1.1);
      eff_d[i]->GetYaxis()->SetNdivisions(1005);
      eff_d[i]->Draw("E");
      leg[0+i]->Draw();
      eff_mc[i]->Draw("E SAME");

      gPad->RedrawAxis();
      ci->Update();
      string pp = p_pdf + ( ( i == 0 ) ? "_p" : "_cos" ) +
         ( ( rew == 1 ) ? "_w" : "" );
      ci->Print( Form("%s_%i.pdf",pp.c_str(),date) );
   }

   // plot ratio
   gStyle->SetStatX(0.91);
   gStyle->SetStatY(0.89);
   gStyle->SetStatW(0.25);
   gStyle->SetFitFormat(".3f");

   vector<double> fit_res;
   for (int i = 0; i < Nh; ++i ) {
      auto& ci = cc[2+i];
      ci->cd();
      gPad->SetLeftMargin(gPad->GetLeftMargin()+0.02);
      gPad->SetRightMargin(gPad->GetRightMargin()-0.02);
      gPad->SetGrid();
      rat0[i]->SetTitle( Form(";%s;#epsilon(data) / #epsilon(MC)",
               ( i == 0 ? "P, GeV/c" : "cos(#Theta)")) );
      SetHstFace(rat0[i]);
      rat0[i]->SetAxisRange(rat_min,rat_max,"Y");
      // rat0[i]->GetXaxis()->SetTitleOffset(1.0);
      rat0[i]->GetYaxis()->SetTitleOffset(1.2);
      rat0[i]->SetLineWidth(2);
      rat0[i]->SetMarkerStyle(20);
      rat0[i]->SetLineColor(kBlack);
      TFitResultPtr rs = rat0[i]->Fit( pl0, "S" );
      leg[2+i]->Draw();
      rat0[i]->Draw("SAME E0");

      double ch2 = pl0->GetChisquare();
      double ndf = pl0->GetNDF();
      double rr0 = pl0->GetParameter(0);
      double er0 = pl0->GetParError(0);
      fit_res.push_back(double(date));
      fit_res.push_back(ch2);
      fit_res.push_back(ndf);
      fit_res.push_back(rr0);
      fit_res.push_back(er0);

      gPad->RedrawAxis();
      ci->Update();
      string pp = p_pdf + "_rat" + ( ( i == 0 ) ? "_p" : "_cos" ) +
         ( ( rew == 1 ) ? "_w" : "" );
      ci->Print( Form("%s_%i.pdf",pp.c_str(),date) );
   }
   return fit_res;
}

// {{{1 Main
//--------------------------------------------------------------------
void eta_eff()
//--------------------------------------------------------------------
{
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(42);
   gStyle->SetStatFont(42);

   size_t Cx = 800, Cy = 640; // canvas sizes, X/Y = 1.25

   vector<double> g_eta;
   vector<double> p_eta;
   // for ( auto date : {2021} ) {
   for ( auto date : {2009, 2012, 2021} ) {
      const int rew = 0;     // 1 - use corrections

      // I. J/Psi->gamma eta, fig B10,B11
      // auto res = plot_pict_gamma_eta(date,rew,Cx,Cy);
      // g_eta.insert(end(g_eta),begin(res),end(res));

      // II. J/Psi->phi eta, fig B18,B19
      // auto res2 = plot_pict_phi_eta(date,rew,Cx,Cy);
      // p_eta.insert(end(p_eta),begin(res2),end(res2));
   }

   string dml(64,'='); // demarcation line

   // I. J/Psi->gamma eta, fig B10,B11
   int Ni = g_eta.size() / 10;

   double av = 0, sig = 0;
   for ( int i = 0; i < Ni; ++i ) {
      int j = 10*i; // P(eta)
      printf("%.0f     P, chi^2/NDF = %4.1f /%2.0f,  "
            "p0 = %.6f +/- %.6f\n",
            g_eta[j],g_eta[j+1],g_eta[j+2],g_eta[j+3],g_eta[j+4]);
      av += g_eta[j+3]/SQ(g_eta[j+4]);
      sig += 1./SQ(g_eta[j+4]);
   }
   if ( Ni > 0 ) {
      av = av/sig;
      sig = 1./sqrt(sig);
      printf("%s\n",dml.c_str());
      printf("average for   P: %.4f +/- %.4f\n\n",av,sig);
   }

   av = 0, sig = 0;
   for ( int i = 0; i < Ni; ++i ) {
      int j = 5 + 10*i; // cos(eta)
      printf("%.0f   cos, chi^2/NDF = %4.1f /%2.0f,  "
            "p0 = %.6f +/- %.6f\n",
            g_eta[j],g_eta[j+1],g_eta[j+2],g_eta[j+3],g_eta[j+4]);
      av += g_eta[j+3]/SQ(g_eta[j+4]);
      sig += 1./SQ(g_eta[j+4]);
   }
   if ( Ni > 0 ) {
      av = av/sig;
      sig = 1./sqrt(sig);
      printf("%s\n",dml.c_str());
      printf("average for cos: %.4f +/- %.4f\n",av,sig);
   }

   // II. J/Psi->phi eta, fig B18,B19
   int Mi = p_eta.size() / 10;

   av = 0, sig = 0;
   for ( int i = 0; i < Mi; ++i ) {
      int j = 10*i; // P(eta)
      printf("%.0f     P, chi^2/NDF = %4.1f /%2.0f,  "
            "p0 = %.6f +/- %.6f\n",
            p_eta[j],p_eta[j+1],p_eta[j+2],p_eta[j+3],p_eta[j+4]);
      av += p_eta[j+3]/SQ(p_eta[j+4]);
      sig += 1./SQ(p_eta[j+4]);
   }
   if ( Mi > 0 ) {
      av = av/sig;
      sig = 1./sqrt(sig);
      printf("%s\n",dml.c_str());
      printf("average for   P: %.4f +/- %.4f\n\n",av,sig);
   }

   av = 0, sig = 0;
   for ( int i = 0; i < Mi; ++i ) {
      int j = 5 + 10*i; // cos(eta)
      printf("%.0f   cos, chi^2/NDF = %4.1f /%2.0f,  "
            "p0 = %.6f +/- %.6f\n",
            p_eta[j],p_eta[j+1],p_eta[j+2],p_eta[j+3],p_eta[j+4]);
      av += p_eta[j+3]/SQ(p_eta[j+4]);
      sig += 1./SQ(p_eta[j+4]);
   }
   if ( Mi > 0 ) {
      av = av/sig;
      sig = 1./sqrt(sig);
      printf("%s\n",dml.c_str());
      printf("average for cos: %.4f +/- %.4f\n",av,sig);
   }
}
