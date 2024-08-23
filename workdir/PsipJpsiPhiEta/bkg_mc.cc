// bkg.mc
// tables of decay Psi(2S) -> J/Psi -> final state,
// to study background; use inclusive MC or sigbal MC

// {{{1 helper functions and constants
//--------------------------------------------------------------------
// GLOBAL: name of folder with root files
string Dir;

//--------------------------------------------------------------------
constexpr double SQ(double x)
//--------------------------------------------------------------------
{
   return x*x;
}

//--------------------------------------------------------------------
void PrtMap(string TblName, const map<string,size_t>& decays)
//--------------------------------------------------------------------
{
   // print sorted by values (in decreasing order of decay frequency)
   vector<pair<string,size_t>> vprt( begin(decays),end(decays) );
   sort( begin(vprt), end(vprt),
         []( pair<string,size_t> elem1, pair<string,size_t> elem2 ) {
         return elem1.second > elem2.second;
         }
       );

   size_t Ntot = 0;
   for ( const auto& p : vprt ) {
      Ntot += p.second;
   }
   double R = 100./double(Ntot);

   printf("%-20.20s: %zu\n",TblName.c_str(),Ntot);
   for ( const auto& p : vprt ) {
      const string& name = p.first;
      size_t ev = p.second;
      double pc = ev*R;
      printf("%-35.35s %9zu (%10.4f %%)\n",name.c_str(),ev,pc);
   }
}

// {{{1 Table
//--------------------------------------------------------------------
void PrtTbl(int date)
//--------------------------------------------------------------------
{
#include "masses.h"

   string mcincfile( Form("mcinc_%02ipsip_all.root",date%100) );
   string fname = Dir + mcincfile;

   // test signal MC
   // string mcsigfile( Form("mcsig_kkmc_%02i.root",date%100) );
   // string fname = Dir + mcsigfile;

   cout << " file: " << fname << endl;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if ( froot == 0 ) {
      cerr << "ERROR in "<< __func__
         << ": can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("PsipJpsiPhiEta");
   // J/Psi -> phi eta
   TTree* a4c = (TTree*)gDirectory->Get("a4c");
   if ( !a4c ) {
      cout << " can not find a4c" << endl;
      exit(EXIT_FAILURE);
   }

   // Declaration of leaves types from "a4c_v709.h"
   // just the ones we use here
   Double_t        Mrec;
   a4c->SetBranchAddress("Mrec",&Mrec);
   Double_t        ch2;
   a4c->SetBranchAddress("ch2",&ch2);
   Double_t        Pkp;
   a4c->SetBranchAddress("Pkp",&Pkp);
   Double_t        Ckp;
   a4c->SetBranchAddress("Ckp",&Ckp);
   Double_t        Pkm;
   a4c->SetBranchAddress("Pkm",&Pkm);
   Double_t        Ckm;
   a4c->SetBranchAddress("Ckm",&Ckm);
   Double_t        Mkk;
   a4c->SetBranchAddress("Mkk",&Mkk);
   Double_t        Mgg;
   a4c->SetBranchAddress("Mgg",&Mgg);
   Double_t        mcmkk;
   a4c->SetBranchAddress("mcmkk",&mcmkk);

   // for table
   Int_t           dec;
   a4c->SetBranchAddress("dec",&dec);
   Int_t           decj;
   a4c->SetBranchAddress("decj",&decj);
   string*         sdecp = nullptr;
   a4c->SetBranchAddress("sdecp",&sdecp);
   string*         sdecj = nullptr;
   a4c->SetBranchAddress("sdecj",&sdecj);

   // Cuts MUST BE THE SAME AS IN "cuts.h"
   auto c_Mrec = [](double Mrec) -> bool{
      return fabs(Mrec-3.097) < 0.005;
   };

   // chi^2 cut
   auto c_chi2 = [](double ch2) -> bool{
      return ch2 < 100;
   };

   // Mphi cut: [dL, dU]
   double dL = 2*Mk; // the left cutoff
   double dU = 1.08; // MUST BE < 1.0835 for signal
   auto c_phi = [dL,dU](double Mkk) -> bool{
      return ( Mkk > dL && Mkk < dU);
   };

   // Meta: central part
   const double seta = 0.008;
   const double weta = 3*seta; // standard: 3x,
                               // uncertainties study: 2x, 4x
   auto c_cpgg = [weta](double Mgg) -> bool{
      return fabs(Mgg-Meta) < weta;
   };

   // 'shift_eta' is the start of the side-band
   double shift_eta = 7*seta;
   auto c_sbgg = [weta,shift_eta](double Mgg) -> bool{
      return (fabs(Mgg-Meta) > shift_eta) &&
             (fabs(Mgg-Meta) < shift_eta+weta);
   };

   // map<string,size_t> Psi2STbl;
   map<string,size_t> JpsiTbl;
   map<string,size_t> JpsiTbl_sb;

   Long64_t nentries = a4c -> GetEntries();
   for ( Long64_t i = 0; i < nentries; ++i ) {
      a4c->GetEntry(i);
      if ( !(mcmkk<1.08) ) continue;
      if ( !c_Mrec(Mrec) ) continue;
      if ( !c_chi2(ch2) ) continue;
      if ( !c_phi(Mkk) ) continue;

      string Psi2S_decay(*sdecp);
      string Jpsi_decay(*sdecj);
      if ( Jpsi_decay == "unknown" ) {
         Jpsi_decay = "Psi(2s)->" + Psi2S_decay;
      } else if ( Jpsi_decay == "9020221 gamma" ) {
         Jpsi_decay = "eta(1405) gamma";
      } else if ( Jpsi_decay == "9000223 gamma" ) {
         Jpsi_decay = "f1(1510) gamma";
      } else if ( Jpsi_decay == "9010221 phi" ) {
         Jpsi_decay = "f0(980) phi";
      }

      if ( c_cpgg(Mgg) ) {              // central part
         // Psi2STbl[Psi2S_decay] += 1;
         JpsiTbl[Jpsi_decay] += 1;
      } else if ( c_sbgg(Mgg) ) {       // side-band
         JpsiTbl_sb[Jpsi_decay] += 1;
      }
   }

   // PrtMap("Psi(2S)", Psi2STbl);
   // printf("\n");
   PrtMap("J/Psi sig",JpsiTbl);
   printf("\n");
   PrtMap("J/Psi sb",JpsiTbl_sb);
}

// {{{1 Main
//--------------------------------------------------------------------
void bkg_mc(int date=2021)
//--------------------------------------------------------------------
{
   gROOT->Reset();
   // gStyle->SetOptStat(0);
   // gStyle->SetOptFit(0);
   // gStyle->SetLegendFont(42);

   //========================================================
   // set the name of the folder with the root files
   Dir = "prod_v709n4/";
   //========================================================

   PrtTbl(date);
}
