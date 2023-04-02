// compare cross-sections
//   -> cmpr_cs_XXXX.pdf

#include <cstdio>
#include <time.h>   // see man strftime

// {{{1~ container for measured energies and cross sections
//--------------------------------------------------------------------
struct Data {
   Data() {
      Ebeam.reserve(100);
      Sig.reserve(100);
      ErrEbeam.reserve(100);
      ErrSig.reserve(100);
   }

   void CheckConsistensy() const {
      // check correctness of vectors:
      unsigned int Np = Ebeam.size();
      if( Sig.size() != Np
            // || ErrEbeam.size() != Np || ErrSig.size() != Np
        ) {
         cerr << "FATAL ERROR::Data() "
            "sizes of vectors are different" << endl
            << " Ebeam.size    = " << Ebeam.size() << endl
            << " Sig.size      = " << Sig.size() << endl
            << " ErrEbeam.size = " << ErrEbeam.size() << endl
            << " ErrSig.size   = " << ErrSig.size() << endl
            ;
         exit(EXIT_FAILURE);
      }
   }

   int Size() const {
      return Ebeam.size();
   }

   std::string         version;
   double              Xisr_min;
   std::vector<double> Ebeam;  // GeV
   std::vector<double> Sig;    // pb
   std::vector<double> ErrEbeam;
   std::vector<double> ErrSig;

   // numbers of energy points for 2012, 2018 and R-scan data
   std::vector<size_t> Ndate;
};

// {{{1~ cs_25Mar23: mass_KK_fit (LH without zero bins)
//--------------------------------------------------------------------
void SetData25M(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "Picts/CS_OLD/cs_25Mar23_results.h"

   dt.version = "25Mar23";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ndate.push_back(Ebeam12.size());
   dt.Ndate.push_back(Ebeam18.size());
   dt.Ndate.push_back(EbeamR.size()-imR);

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

// {{{1~ 01Apr23: mass_KK_fit (LH fix bug with zero bins)
//--------------------------------------------------------------------
void SetData01A(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "Picts/cs_01Apr23_results.h"

   dt.version = "01Apr23";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ndate.push_back(Ebeam12.size());
   dt.Ndate.push_back(Ebeam18.size());
   dt.Ndate.push_back(EbeamR.size()-imR);

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

// {{{1~ SB_25Mar23: side-band subtraction: 1.01 < Mkk < 1.03 GeV
//--------------------------------------------------------------------
void SetDataSB_25M(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "Picts/CS_OLD/cs_SB_25Mar23_results.h"

   dt.version = "SB_25Mar23";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ndate.push_back(Ebeam12.size());
   dt.Ndate.push_back(Ebeam18.size());
   dt.Ndate.push_back(EbeamR.size()-imR);

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

// {{{1~ SB_01Apr23: side-band subtraction: 1.01 < Mkk < 1.03 GeV
//--------------------------------------------------------------------
void SetDataSB_01A(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "Picts/cs_SB_01Apr23_results.h"

   dt.version = "SB_01Apr23";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ndate.push_back(Ebeam12.size());
   dt.Ndate.push_back(Ebeam18.size());
   dt.Ndate.push_back(EbeamR.size()-imR);

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

// {{{1~ Systematic: ANGL variation MINUS - PLUS
//--------------------------------------------------------------------
void SetDataAnglM(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "mkk_inter/LH_angM/cs_angM_results.h"

   dt.version = "#vartheta=-0.58";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ndate.push_back(Ebeam12.size());
   dt.Ndate.push_back(Ebeam18.size());
   dt.Ndate.push_back(EbeamR.size()-imR);

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

//--------------------------------------------------------------------
void SetDataAnglP(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "mkk_inter/LH_angP/cs_angP_results.h"

   dt.version = "#vartheta=+0.58";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ndate.push_back(Ebeam12.size());
   dt.Ndate.push_back(Ebeam18.size());
   dt.Ndate.push_back(EbeamR.size()-imR);

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

// {{{1~ Systematic: chi2 variation MIN - MAX
//--------------------------------------------------------------------
void SetDataCh2MIN(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "mkk_inter/LH_ch2_60/cs_ch2_60_results.h"

   dt.version = "chi2 MIN";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ndate.push_back(Ebeam12.size());
   dt.Ndate.push_back(Ebeam18.size());
   dt.Ndate.push_back(EbeamR.size()-imR);

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

//--------------------------------------------------------------------
void SetDataCh2MAX(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "mkk_inter/LH_ch2_100/cs_ch2_100_results.h"

   dt.version = "chi2 MAX";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ndate.push_back(Ebeam12.size());
   dt.Ndate.push_back(Ebeam18.size());
   dt.Ndate.push_back(EbeamR.size()-imR);

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

// {{{1~ Systematic: Eta window variation MIN - MAX
//--------------------------------------------------------------------
void SetDataWETAMIN(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "mkk_inter/LH_weta2/cs_weta_2_results.h"

   dt.version = "weta MIN";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ndate.push_back(Ebeam12.size());
   dt.Ndate.push_back(Ebeam18.size());
   dt.Ndate.push_back(EbeamR.size()-imR);

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

//--------------------------------------------------------------------
void SetDataWETAMAX(Data& dt) {
//--------------------------------------------------------------------
   double Xisr_min;
#include "mkk_inter/LH_weta4/cs_weta_4_results.h"

   dt.version = "chi2 MAX";
   dt.Xisr_min = Xisr_min;

   int imR = 5; // imR=0: all, 1: skip 2900, etc...

   dt.Ndate.push_back(Ebeam12.size());
   dt.Ndate.push_back(Ebeam18.size());
   dt.Ndate.push_back(EbeamR.size()-imR);

   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam12), end(Ebeam12));
   dt.Ebeam.insert(end(dt.Ebeam), begin(Ebeam18), end(Ebeam18));
   dt.Ebeam.insert(end(dt.Ebeam), begin(EbeamR)+imR, end(EbeamR));

   dt.Sig.insert(end(dt.Sig), begin(Sig12), end(Sig12));
   dt.Sig.insert(end(dt.Sig), begin(Sig18), end(Sig18));
   dt.Sig.insert(end(dt.Sig), begin(SigR)+imR, end(SigR));

   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam12), end(ErrEbeam12));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeam18), end(ErrEbeam18));
   dt.ErrEbeam.insert(end(dt.ErrEbeam),
         begin(ErrEbeamR)+imR, end(ErrEbeamR));

   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig12), end(ErrSig12));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSig18), end(ErrSig18));
   dt.ErrSig.insert(end(dt.ErrSig), begin(ErrSigR)+imR, end(ErrSigR));

   dt.CheckConsistensy();
}

// {{{1~ print relative error for cross-section
//--------------------------------------------------------------------
void prtSysCS(string title, string fname, string tname,
      const vector<size_t>& Ndate, const vector<double>& Rat ) {
//--------------------------------------------------------------------
   const char* shift = "   "; // shift = 3 spaces
   const char* comma = ", ";
   if ( size(Ndate) != 3 ) {
      cerr << " FATAL ERROR: prtSysCS: size(Ndate)= "
         << size(Ndate) << endl;
      exit(1);
   }

   time_t temp = time(NULL);
   struct tm * timeptr = localtime(&temp);
   char buf[32];
   strftime(buf,sizeof(buf),"%d%b%y",timeptr);
   string dat(buf);
   string filename = fname + '_' + dat + ".h";
   FILE* fp = fopen(filename.c_str(),"w");
   if ( !fp ) {
      cerr << "Error open file " << filename << endl;
      fp = stdout;
   }

   fprintf(fp,"\n");
   fprintf(fp,"%s// %s\n",shift,title.c_str());
   string sufs[] {"12","18","R"};
   for ( int k = 0, l = 0; k < 3; k++) {
      int n =  Ndate[k];
      const string& suf = sufs[k];
      fprintf(fp,"%svector<double> RelSys%s%s = {\n",
            shift,tname.c_str(),suf.c_str());
      for(int i = 0; i < n; i += 5) {
         fprintf(fp,"%s%s",shift,shift);
         for(int j = i; j < (i+5<n ? i+5 : n); j++) {
            fprintf(fp,"%.2e%s",Rat[l+j], (j!=n-1) ? comma : "" );
         }
         fprintf(fp,"\n");
      }
      l += n;
      fprintf(fp,"%s};\n",shift);
   }
}

// {{{1~ systematic FIT vs SB
//--------------------------------------------------------------------
void cmpr_cs_FIT() {
//--------------------------------------------------------------------
   Data dt[2];
   bool save_header = false;

   //* compare old-new mass_KK_fit (LH)
   // SetData01A(dt[0]); // new
   // SetData25M(dt[1]); // old
   // double Rmax=10., Rmin=-40.;
   // string pdf( "cmpr_cs_01A_25M.pdf" ); // <- PDF

   SetData01A(dt[0]); // mass_KK_fit (LH)
   dt[0].version = "LH-fit of M_{KK} with interference";
   SetDataSB_01A(dt[1]); // SB: 1.01 < Mkk < 1.03 GeV
   dt[1].version = "Side-band subtraction";
   double Rmax=15., Rmin=-15.;
   string fname = "sys_FIT";
   string pdf = fname + "_01Apr23.pdf";
   save_header = true;

   int Nd = dt[0].Size();
   if ( dt[1].Size() != Nd ) {
      cerr << "FATAL ERROR:: sizes of Data are different" << endl
         << " Size1= " << dt[0].Size()
         << " Size2= " << dt[1].Size() << endl;
      exit(EXIT_FAILURE);
   }

   // Rat - relative error
   // R - percent
   // Ro - under and over [Rmin,Rmax]
   vector<double> Rat(Nd), R(Nd), Ro(Nd);
   for(unsigned int i = 0; i < Nd; i++) {
      Rat[i] = 1. - dt[1].Sig[i] / dt[0].Sig[i];
      R[i] = 100*Rat[i];
      if ( R[i] > Rmax ) {
         Ro[i] = Rmax-0.1;
      } else if ( R[i] < Rmin ) {
         Ro[i] = Rmin+0.1;
      } else {
         Ro[i] = 1e3; // something invisible
      }
   }

   if ( save_header ) {
      string titleH = dt[0].version + " vs " + dt[1].version;
      prtSysCS(titleH, fname, "CS", dt[0].Ndate, Rat);
   }

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1->Print((pdf+"[").c_str()); // just open pdf-file

   // 1) cross-sections on one graph
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   TGraphErrors* cs[2];
   for ( int i = 0; i < 2; ++i ) {
      cs[i] = new TGraphErrors( Nd,
            dt[i].Ebeam.data(), dt[i].Sig.data(),
            dt[i].ErrEbeam.data(), dt[i].ErrSig.data() );
   }

   TGraph* gr = cs[0];
   string title = string(";center of mass energy (GeV)"
         ";cross-section (pb)");
   gr->SetTitle(title.c_str());
   gr->GetXaxis()->CenterTitle();
   gr->GetYaxis()->CenterTitle();
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   gr->GetXaxis()->SetLabelSize(0.03);
   gr->GetYaxis()->SetLabelSize(0.03);
   // gr->GetXaxis()->SetTitleOffset(1.1);
   gr->GetYaxis()->SetTitleOffset(1.1);

   gr->SetMaximum(1e4);
   gr->SetMinimum(5.);

   gr->SetMarkerStyle(22);
   gr->SetMarkerColor(kRed+1);
   // gr->SetMarkerSize(0.7);
   gr->Draw("AP");

   gr = cs[1];
   gr->SetMarkerStyle(23);
   gr->SetMarkerColor(kBlue+1);
   // gr->SetMarkerSize(0.7);
   gr->Draw("P"); // on top of pict


   TLegend* leg = new TLegend(0.12,0.72,0.50,0.88);
   leg->AddEntry( cs[0],
         Form("%s",dt[0].version.c_str()), "EP");
         // Form("cross-section %s",dt[0].version.c_str()), "EP");
   leg->AddEntry( cs[1],
         Form("%s",dt[1].version.c_str()), "EP");
         // Form("cross-section %s",dt[1].version.c_str()), "EP");
   leg->Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(false);

   TGraph* rat = new TGraph( Nd, dt[0].Ebeam.data(), R.data() );

   gr = rat;
   string titleR = string(";center of mass energy (GeV)"
         ";relative difference (%)");
   gr->SetTitle(titleR.c_str());
   gr->GetXaxis()->CenterTitle();
   gr->GetYaxis()->CenterTitle();
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   gr->GetXaxis()->SetLabelSize(0.03);
   gr->GetYaxis()->SetLabelSize(0.03);
   // gr->GetXaxis()->SetTitleOffset(1.1);
   gr->GetYaxis()->SetTitleOffset(1.1);

   rat->SetMaximum(Rmax);
   rat->SetMinimum(Rmin);

   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kBlue);
   // gr->SetMarkerSize(0.9);
   gr->Draw("AP");

   TGraph* rat1 = new TGraph( Nd, dt[0].Ebeam.data(), Ro.data() );
   rat1->SetMarkerStyle(20);
   rat1->SetMarkerColor(kRed);
   rat1->Draw("P"); // on top of pict

   gPad -> RedrawAxis();
   c1 -> Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // just close pdf-file
}

// {{{1~ Systematic ANGL
//--------------------------------------------------------------------
void cmpr_cs_SYSANGL() {
//--------------------------------------------------------------------
   Data dt[3];

   SetData01A(dt[0]);
   dt[0].version = "std: #vartheta=0";
   SetDataAnglM(dt[1]);
   SetDataAnglP(dt[2]);
   string fname = "sys_ANGL";
   string pdf = fname + "_01Apr23.pdf";

   int Nd = dt[0].Size();
   if ( dt[1].Size() != Nd || dt[2].Size() != Nd  ) {
      cerr << "FATAL ERROR:: sizes of Data are different" << endl
         << " Size0= " << dt[0].Size()
         << " Size1= " << dt[1].Size()
         << " Size2= " << dt[1].Size() << endl;
      exit(EXIT_FAILURE);
   }

   // Rat - max relative error
   // R - percent
   vector<double> Rat(Nd), R1(Nd), R2(Nd);
   double Rmax=15., Rmin=-15.;
   for(unsigned int i = 0; i < Nd; i++) {
      double rat1 = 1. - dt[1].Sig[i] / dt[0].Sig[i];
      double rat2 = 1. - dt[2].Sig[i] / dt[0].Sig[i];
      R1[i] = 100*rat1;
      R2[i] = 100*rat2;
      if ( fabs(rat1) > fabs(rat2) ) {
         Rat[i] = rat1;
      } else {
         Rat[i] = rat2;
      }
   }
   string titleH = dt[0].version + " vs " + dt[1].version
      + " and " + dt[2].version;
   prtSysCS(titleH, fname, "ANGL", dt[0].Ndate, Rat);

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1->Print((pdf+"[").c_str()); // just open pdf-file

   // 1) cross-sections on one graph
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   TGraphErrors* cs[3];
   for ( int i = 0; i < 3; ++i ) {
      cs[i] = new TGraphErrors( Nd,
            dt[i].Ebeam.data(), dt[i].Sig.data(),
            dt[i].ErrEbeam.data(), dt[i].ErrSig.data() );
   }

   TGraph* gr = cs[0];
   string title = string(";center of mass energy (GeV)"
         ";cross-section (pb)");
   gr->SetTitle(title.c_str());
   gr->GetXaxis()->CenterTitle();
   gr->GetYaxis()->CenterTitle();
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   gr->GetXaxis()->SetLabelSize(0.03);
   gr->GetYaxis()->SetLabelSize(0.03);
   // gr->GetXaxis()->SetTitleOffset(1.1);
   gr->GetYaxis()->SetTitleOffset(1.1);

   gr->SetMaximum(1e4);
   gr->SetMinimum(5.);

   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kRed+1);
   // gr->SetMarkerSize(0.7);
   gr->Draw("AP");

   cs[1]->SetMarkerStyle(23);
   cs[1]->SetMarkerColor(kBlue+1);
   cs[1]->Draw("P"); // on top of pict
   cs[2]->SetMarkerStyle(22);
   cs[2]->SetMarkerColor(kGreen-1);
   cs[2]->Draw("P"); // on top of pict


   TLegend* leg = new TLegend(0.12,0.72,0.50,0.88);
   for ( int i = 0; i < 3; ++i ) {
      leg->AddEntry( cs[i],
            Form("%s",dt[i].version.c_str()), "EP");
   }
   leg->Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(false);

   TGraph* rat1 = new TGraph( Nd, dt[0].Ebeam.data(), R1.data() );

   gr = rat1;
   string titleR = string(";center of mass energy (GeV)"
         ";relative difference (%)");
   gr->SetTitle(titleR.c_str());
   gr->GetXaxis()->CenterTitle();
   gr->GetYaxis()->CenterTitle();
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   gr->GetXaxis()->SetLabelSize(0.03);
   gr->GetYaxis()->SetLabelSize(0.03);
   // gr->GetXaxis()->SetTitleOffset(1.1);
   gr->GetYaxis()->SetTitleOffset(1.1);

   rat1->SetMaximum(Rmax);
   rat1->SetMinimum(Rmin);

   gr->SetMarkerStyle(23);
   gr->SetMarkerColor(kBlue+1);
   // gr->SetMarkerSize(0.9);
   gr->Draw("AP");

   TGraph* rat2 = new TGraph( Nd, dt[0].Ebeam.data(), R2.data() );
   rat2->SetMarkerStyle(22);
   rat2->SetMarkerColor(kGreen-1);
   rat2->Draw("P"); // on top of pict

   TLegend* legD = new TLegend(0.12,0.80,0.52,0.88);
   legD->AddEntry( rat1, "1 #minus "
         "#sigma(#vartheta=-0.58) / #sigma(#vartheta=0)", "EP" );
   legD->AddEntry( rat2, "1 #minus "
         "#sigma(#vartheta=+0.58) / #sigma(#vartheta=0)", "EP" );
   legD->Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // just close pdf-file
}

// {{{1~ Systematic CH2
//--------------------------------------------------------------------
void cmpr_cs_SYSCH2() {
//--------------------------------------------------------------------
   Data dt[3];

   SetData01A(dt[0]); // mass_KK_fit (LH)
   dt[0].version = "std: #chi^{2} < 80";
   SetDataCh2MIN(dt[1]);
   // dt[1].version = "#chi^{2}<70";
   dt[1].version = "#chi^{2} < 60";
   SetDataCh2MAX(dt[2]);
   // dt[2].version = "#chi^{2} < 90";
   dt[2].version = "#chi^{2} < 100";
   string pdf( "sys_CH_01apr23.pdf" ); // <- PDF

   int Nd = dt[0].Size();
   if ( dt[1].Size() != Nd || dt[2].Size() != Nd  ) {
      cerr << "FATAL ERROR:: sizes of Data are different" << endl
         << " Size0= " << dt[0].Size()
         << " Size1= " << dt[1].Size()
         << " Size2= " << dt[1].Size() << endl;
      exit(EXIT_FAILURE);
   }

   // points with biggest statistic (R, 2012, 2018)
   vector<double> Ems {  3.08,
      3.095840, 3.097227, 3.098354,
      3.095726, 3.096203, 3.096986, 3.097226, 3.097654, 3.098728
   };
   // Rat - max relative error
   // R - percent
   vector<double> Rat(Nd), R1(Nd), R2(Nd);
   double Rmax=10., Rmin=-10.;
   vector<double> Rms(Ems.size());
   for( size_t i = 0; i < Nd; i++) {
      double rat1 = 1. - dt[1].Sig[i] / dt[0].Sig[i];
      double rat2 = 1. - dt[2].Sig[i] / dt[0].Sig[i];
      R1[i] = 100*rat1;
      R2[i] = 100*rat2;
      if ( fabs(rat1) > fabs(rat2) ) {
         Rat[i] = rat1;
      } else {
         Rat[i] = rat2;
      }

      for( size_t ie = 0; ie < Ems.size(); ++ie ) {
         if ( fabs(Ems[ie] - dt[0].Ebeam[i]) < 1e-6 ) {
            Rms[ie] = 100*Rat[i];
            printf(" E= %.6f, rat= %.2f %%\n", Ems[ie], Rms[ie]);
            break;
         }
      }
   }

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1->Print((pdf+"[").c_str()); // just open pdf-file

   // 1) cross-sections on one graph
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   TGraphErrors* cs[3];
   for ( int i = 0; i < 3; ++i ) {
      cs[i] = new TGraphErrors( Nd,
            dt[i].Ebeam.data(), dt[i].Sig.data(),
            dt[i].ErrEbeam.data(), dt[i].ErrSig.data() );
   }

   TGraph* gr = cs[0];
   string title = string(";center of mass energy (GeV)"
         ";cross-section (pb)");
   gr->SetTitle(title.c_str());
   gr->GetXaxis()->CenterTitle();
   gr->GetYaxis()->CenterTitle();
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   gr->GetXaxis()->SetLabelSize(0.03);
   gr->GetYaxis()->SetLabelSize(0.03);
   // gr->GetXaxis()->SetTitleOffset(1.1);
   gr->GetYaxis()->SetTitleOffset(1.1);

   gr->SetMaximum(1e4);
   gr->SetMinimum(5.);

   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kRed+1);
   // gr->SetMarkerSize(0.7);
   gr->Draw("AP");

   cs[1]->SetMarkerStyle(23);
   cs[1]->SetMarkerColor(kBlue+1);
   cs[1]->Draw("P"); // on top of pict
   cs[2]->SetMarkerStyle(22);
   cs[2]->SetMarkerColor(kGreen-1);
   cs[2]->Draw("P"); // on top of pict


   TLegend* leg = new TLegend(0.12,0.72,0.50,0.88);
   for ( int i = 0; i < 3; ++i ) {
      leg->AddEntry( cs[i],
            Form("%s",dt[i].version.c_str()), "EP");
   }
   leg->Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(false);

   TGraph* rat1 = new TGraph( Nd, dt[0].Ebeam.data(), R1.data() );

   gr = rat1;
   string titleR = string(";center of mass energy (GeV)"
         ";relative difference (%)");
   gr->SetTitle(titleR.c_str());
   gr->GetXaxis()->CenterTitle();
   gr->GetYaxis()->CenterTitle();
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   gr->GetXaxis()->SetLabelSize(0.03);
   gr->GetYaxis()->SetLabelSize(0.03);
   // gr->GetXaxis()->SetTitleOffset(1.1);
   gr->GetYaxis()->SetTitleOffset(1.1);

   rat1->SetMaximum(Rmax);
   rat1->SetMinimum(Rmin);

   gr->SetMarkerStyle(23);
   gr->SetMarkerColor(kBlue+1);
   // gr->SetMarkerSize(0.9);
   gr->Draw("AP");

   TGraph* rat2 = new TGraph( Nd, dt[0].Ebeam.data(), R2.data() );
   rat2->SetMarkerStyle(22);
   rat2->SetMarkerColor(kGreen-1);
   rat2->Draw("P"); // on top of pict

   TGraph* rat3 = new TGraph( Ems.size(), Ems.data(), Rms.data() );
   rat3->SetMarkerStyle(41);
   // rat3->SetMarkerSize(0.4);
   rat3->SetMarkerColor(kRed);
   rat3->Draw("P"); // on top of pict

   TLegend* legD = new TLegend(0.35,0.74,0.88,0.88);
   legD->SetHeader("variation of #chi^{2} cut:","C");
   legD->AddEntry( rat1, "1 #minus "
         "#sigma(#chi^{2} < 60) / #sigma(#chi^{2} < 80)", "EP" );
   legD->AddEntry( rat2, "1 #minus "
         "#sigma(#chi^{2} < 100) / #sigma(#chi^{2} < 80)", "EP" );
   legD->Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // just close pdf-file
}

// {{{1~ Systematic Weta
//--------------------------------------------------------------------
void cmpr_cs_SYSWETA() {
//--------------------------------------------------------------------
   Data dt[3];

   SetData01A(dt[0]); // mass_KK_fit (LH)
   dt[0].version = "std: |M(#gamma#gamma) #minus M(#eta)| < 24 MeV";
   SetDataWETAMIN(dt[1]);
   dt[1].version = "|M(#gamma#gamma) #minus M(#eta)| < 16 MeV";
   SetDataWETAMAX(dt[2]);
   dt[2].version = "|M(#gamma#gamma) #minus M(#eta)| < 32 MeV";
   string pdf( "sys_Weta_01apr23.pdf" ); // <- PDF

   int Nd = dt[0].Size();
   if ( dt[1].Size() != Nd || dt[2].Size() != Nd  ) {
      cerr << "FATAL ERROR:: sizes of Data are different" << endl
         << " Size0= " << dt[0].Size()
         << " Size1= " << dt[1].Size()
         << " Size2= " << dt[1].Size() << endl;
      exit(EXIT_FAILURE);
   }

   // points with biggest statistic (R, 2012, 2018)
   vector<double> Ems {  3.08,
      3.095840, 3.097227, 3.098354,
      3.095726, 3.096203, 3.096986, 3.097226, 3.097654, 3.098728
   };
   // Rat - max relative error
   // R - percent
   vector<double> Rat(Nd), R1(Nd), R2(Nd);
   double Rmax=10., Rmin=-10.;
   vector<double> Rms(Ems.size());
   for( size_t i = 0; i < Nd; i++) {
      double rat1 = 1. - dt[1].Sig[i] / dt[0].Sig[i];
      double rat2 = 1. - dt[2].Sig[i] / dt[0].Sig[i];
      R1[i] = 100*rat1;
      R2[i] = 100*rat2;
      if ( fabs(rat1) > fabs(rat2) ) {
         Rat[i] = rat1;
      } else {
         Rat[i] = rat2;
      }

      for( size_t ie = 0; ie < Ems.size(); ++ie ) {
         if ( fabs(Ems[ie] - dt[0].Ebeam[i]) < 1e-6 ) {
            Rms[ie] = 100*Rat[i];
            printf(" E= %.6f, rat= %.2f %%\n", Ems[ie], Rms[ie]);
            break;
         }
      }
   }

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1->Print((pdf+"[").c_str()); // just open pdf-file

   // 1) cross-sections on one graph
   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy();

   TGraphErrors* cs[3];
   for ( int i = 0; i < 3; ++i ) {
      cs[i] = new TGraphErrors( Nd,
            dt[i].Ebeam.data(), dt[i].Sig.data(),
            dt[i].ErrEbeam.data(), dt[i].ErrSig.data() );
   }

   TGraph* gr = cs[0];
   string title = string(";center of mass energy (GeV)"
         ";cross-section (pb)");
   gr->SetTitle(title.c_str());
   gr->GetXaxis()->CenterTitle();
   gr->GetYaxis()->CenterTitle();
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   gr->GetXaxis()->SetLabelSize(0.03);
   gr->GetYaxis()->SetLabelSize(0.03);
   // gr->GetXaxis()->SetTitleOffset(1.1);
   gr->GetYaxis()->SetTitleOffset(1.1);

   gr->SetMaximum(1e4);
   gr->SetMinimum(5.);

   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kRed+1);
   // gr->SetMarkerSize(0.7);
   gr->Draw("AP");

   cs[1]->SetMarkerStyle(23);
   cs[1]->SetMarkerColor(kBlue+1);
   cs[1]->Draw("P"); // on top of pict
   cs[2]->SetMarkerStyle(22);
   cs[2]->SetMarkerColor(kGreen-1);
   cs[2]->Draw("P"); // on top of pict


   TLegend* leg = new TLegend(0.12,0.72,0.50,0.88);
   for ( int i = 0; i < 3; ++i ) {
      leg->AddEntry( cs[i],
            Form("%s",dt[i].version.c_str()), "EP");
   }
   leg->Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(false);

   TGraph* rat1 = new TGraph( Nd, dt[0].Ebeam.data(), R1.data() );

   gr = rat1;
   string titleR = string(";center of mass energy (GeV)"
         ";relative difference (%)");
   gr->SetTitle(titleR.c_str());
   gr->GetXaxis()->CenterTitle();
   gr->GetYaxis()->CenterTitle();
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   gr->GetXaxis()->SetLabelSize(0.03);
   gr->GetYaxis()->SetLabelSize(0.03);
   // gr->GetXaxis()->SetTitleOffset(1.1);
   gr->GetYaxis()->SetTitleOffset(1.1);

   rat1->SetMaximum(Rmax);
   rat1->SetMinimum(Rmin);

   gr->SetMarkerStyle(23);
   gr->SetMarkerColor(kBlue+1);
   // gr->SetMarkerSize(0.9);
   gr->Draw("AP");

   TGraph* rat2 = new TGraph( Nd, dt[0].Ebeam.data(), R2.data() );
   rat2->SetMarkerStyle(22);
   rat2->SetMarkerColor(kGreen-1);
   rat2->Draw("P"); // on top of pict

   TGraph* rat3 = new TGraph( Ems.size(), Ems.data(), Rms.data() );
   rat3->SetMarkerStyle(41);
   // rat3->SetMarkerSize(0.4);
   rat3->SetMarkerColor(kRed);
   rat3->Draw("P"); // on top of pict

   TLegend* legD = new TLegend(0.35,0.74,0.88,0.88);
   legD->SetHeader("variation of M(#gamma#gamma) window","C");
   legD->AddEntry( rat1, "1 #minus "
         "#sigma(tight cut) / #sigma(std)", "EP" );
   legD->AddEntry( rat2, "1 #minus "
         "#sigma(loose cut) / #sigma(std)", "EP" );
   legD->Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   c1->Print(pdf.c_str()); // add to pdf-file

   c1->Print((pdf+"]").c_str()); // just close pdf-file
}

// {{{1~ Main
//--------------------------------------------------------------------
void cmpr_cs() {
//--------------------------------------------------------------------
   gROOT->Reset();

   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   // gStyle->SetLegendTextSize(0.05);

   // cmpr_cs_FIT();
   // cmpr_cs_SYSANGL();

   // cmpr_cs_SYSCH2();
   // cmpr_cs_SYSWETA();
}
