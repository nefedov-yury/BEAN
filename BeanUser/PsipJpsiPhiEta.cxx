//======================================================================//
//                                                                      //
// PsipJpsiPhiEta:                                                      //
// e+ e- -> Psi(2S) -> pi+ pi- J/Psi                                    //
//                             |-> phi       eta                        //
//                                 |-> K+K-  |-> 2gamma                 //
//                                                                      //
//======================================================================//

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TNtupleD.h>
#include <TDatabasePDG.h>

// #include "CLHEP/Units/PhysicalConstants.h"
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Geometry/Point3D.h>
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "AbsCor/AbsCor.h"
#include "VertexFit/VertexDbSvc.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include "MagneticField/MagneticFieldSvc.h"
#include "EventTag/EventTagSvc.h"
#include "EventTag/DecayTable.h"

#include "DstEvtRecTracks.h"
#include "ReadDst.h"


using namespace std;

//-------------------------------------------------------------------------
// Structure to save variables for a single event
//-------------------------------------------------------------------------
struct Select_PsipJpsiPhiEta {
   // Run-info
   int runNo;              // run-number
   int event;              // event-number
   HepLorentzVector LVcms; // Momentum in center of mass system
   Hep3Vector xorig;       // interaction point from DB

   // MC information:
   int decPsip, decJpsi;      // decay codes for MC events
   int dec_eta;               // 1 if eta->2gamma
   Hep3Vector mcPip, mcPim;   // true pi+/pi- momentums from Psi(2S) decay
   double mc_mkk;             // true inv. mass of K+K- from phi decay

   long int pip_pdg, pim_pdg; // pdg of the best selected candidates
   int      pip_ok,  pim_ok;

   double Mrec_sig;           // recoil mass of the true decay
   double Mrec_sig_mindp;
   RecMdcKalTrack* Pip_sig;   // pi+ Mrec_sig
   RecMdcKalTrack* Pim_sig;   // pi- Mrec_sig

   // soft pions
   set<RecMdcKalTrack*> good_pip;   // positive tracks
   set<RecMdcKalTrack*> good_pim;   // negative tracks

   // recoil masses of pi+ pi-
   vector<double> Mrec;     // all good recoil masses
   double Mrec_best;        // the best recoil mass - closest to M(J/Psi)
   double cosPM, invPM;     // cos(Theta pi+ pi-) & inv.mass of the best
   RecMdcKalTrack* trk_Pip; // pi+ best
   RecMdcKalTrack* trk_Pim; // pi- best

   // Kaons candidate
   vector<RecMdcKalTrack*> trk_Kp;  // positive tracks
   vector<RecMdcKalTrack*> trk_Km;  // negative tracks

   // gamma tracks
   vector<RecEmcShower*> gtrk;
   vector<double> angt;             // angle with closest charged track
   vector<HepLorentzVector> Pg;     // 4-momentum of gammas
   vector<RecEmcShower*> g4f;       // two best photons after 4C-fit

   Select_PsipJpsiPhiEta() {
      decPsip = decJpsi = -1;
      dec_eta = 0;
      mc_mkk = 0;
      pip_pdg = 0; pim_pdg = 0;
      pip_ok = -1; pim_ok = -1;
      Mrec_sig = 0;
      Mrec_sig_mindp = 0;
      Pip_sig = Pim_sig = nullptr;

      trk_Pip = trk_Pim = nullptr;

      g4f.resize(2,nullptr);
   }
};

typedef Select_PsipJpsiPhiEta Select;

//-------------------------------------------------------------------------
// Structure for tree
//-------------------------------------------------------------------------
struct Xnt1 {
   float Mrs;            // recoil mass of true signal pi+pi- for decPsip==64
   float Mrs_mindp;      // min delta momentum to find Mrs
   float MrsW;           // weight of Mrs
   vector<float> Mrec;   // recoil masses of all pairs except true signal MC
   float Mrbest;         // recoil mass closest to M(J/Psi)
   float Ptp;            // Pt(pi+) for best
   float Cpls;           // cos(Theta(pi+)) for best
   float Ptm;            // Pt(pi-) for best
   float Cmns;           // cos(Theta(pi-)) for best
   UShort_t nch;         // N charged tracks
   UShort_t Nm;          // number of els in Mrec: use @Mrec.size()
   UShort_t decPsip;     // decay code of Psi(2S)
   float mcmkk;          // true inv. mass of K+K- from phi decay
} xnt1;

//    float cosPM;          // cos(Theta(pi+ pi-)) for best
//    float invPM;          // invariant mass of (pi+, pi-) for best
//    float Ppls;           // momentum of pi+ for best
//    float Pmns;           // momentum of pi- for best

//-------------------------------------------------------------------------
// Global variables
//-------------------------------------------------------------------------
static const double beam_angle = 0.011; // 11 mrad

// masses of particles (GeV)           from PDG:
static const double mpsip  = 3.686097; // 3686.097  +/- 0.025   MeV
static const double mjpsi  = 3.096916; // 3096.916  +/- 0.011   MeV
static const double mpi    = 0.13957;  // 139.57018 +/- 0.00035 MeV
static const double mpi0   = 0.13498;  // 134.9766  +/- 0.0006  MeV
static const double meta   = 0.547862; // 547.862   +/- 0.017   MeV
static const double momega = 0.78265;  // 782.65    +/- 0.12    MeV
static const double mk     = 0.493677; // 493.677   +/- 0.016   MeV
static const double mk0    = 0.497611; // 497.611   +/- 0.013   MeV
static const double mphi   = 1.019461; //1019.461   +/- 0.019   MeV

static AbsCor* m_abscor = 0;
static EventTagSvc* m_EventTagSvc = 0;

// histograms
static vector<TH1*> hst;
static vector<TNtupleD*> m_tuple;
static TTree* m_nt1 = nullptr;

// decays classification
static DecayTable PsipTbl;  // narrow window [3.092, 3.102]
static DecayTable PsipTbl2; // wide window [3.055, 3.145]
static DecayTable PsipTbl4C;// side-band for 4C KF
static DecayTable JpsiTbl;  // final phi-eta selection
static DecayTable JpsiTbl2; // phi-eta: Mrec in [3.055, 3.145]
static DecayTable JpsiTblE; // for effitiency selection

// container for warnings
static map<string,int> warning_msg;

// expect one data taking period for all events
static int DataPeriod = 0;

static bool isMC = false;

//-------------------------------------------------------------------------
// Functions: use C-linkage names
//-------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif

//-------------------------------------------------------------------------
inline void Warning(const string& msg) {
   warning_msg[msg] += 1;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
inline double RtoD(double ang) {
   return ang*180/M_PI;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
inline double SQ(double x) {
   return x*x;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
static bool SelectPM(double cosPM, double invPM) {
//-------------------------------------------------------------------------
   // criteria for pair pions (+/-) for selection good Mrec
   bool ret = true;
   if ( (cosPM > 0.80) ||         // flying in one direction
         (abs(invPM-mk0) < 0.008 ) // pions from K^0_s decay
      ) {
      ret = false;
   }
   return ret;
}

//-------------------------------------------------------------------------
static double ReWeightTrkPid(int Kp, double Pt) {
//-------------------------------------------------------------------------
// This correction is based on "prod-12/13eff"
// and independent of cos(Theta)
// Kp = 1 for kaons and Kp = 0 for pions

   double W = 1.;

   auto CUBE = [](double x)-> double{return x*x*x;};
   if ( Kp == 1 ) {             // kaons
      Pt = max( 0.1, Pt );
      Pt = min( 1.4, Pt );
      if ( DataPeriod == 2009 ) {
         W = 1.00931 - 0.02363 * Pt;
         if ( Pt < 0.2 ) {
            W = 0.9278;
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
         }
         W = cK12->Eval(Pt);
      }
   } else if ( Kp == 0 ) {      // pions
      Pt = max( 0.05, Pt );
      Pt = min( 0.4, Pt );
      if ( DataPeriod == 2009 ) {
//          W = 0.9878 + CUBE(0.0219/Pt);
         W = 0.9863 + CUBE(0.02234/Pt); // Mar 2020
      } else if ( DataPeriod == 2012 ) {
//          W = 0.9859 + CUBE(0.02974/Pt);
         W = 0.9856 + CUBE(0.02967/Pt); // Mar 2020
      }
   }
   return W;
}

//-------------------------------------------------------------------------
void PsipJpsiPhiEtaStartJob(ReadDst* selector) {
//-------------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
   }

   hst.resize(500,nullptr);
   m_tuple.resize(10,nullptr);

   // init Absorption Correction
   m_abscor = new AbsCor(selector->AbsPath("Analysis/AbsCor"));

   // initialize EventTag
   m_EventTagSvc = EventTagSvc::instance();
   if ( !m_EventTagSvc->IsInitialized() ) {
      // set paths to pdg & decayCodes files:
      m_EventTagSvc->setPdtFile(
         selector->AbsPath( "Analysis/EventTag/share/pdt_bean.table" )
                               );
      m_EventTagSvc->setDecayTabsFile(
         selector->AbsPath(
            "Analysis/EventTag/share/DecayCodes/dcode_charmonium.txt"
         )                           );
      m_EventTagSvc->setIgnorePhotons(false); // for "dcode_charmonium.txt"
      if ( selector->Verbose() ) {
         m_EventTagSvc->setVerbose(1);
      }
      m_EventTagSvc->initialize();
   } else {
      cout << " WARNING in " << __func__ << ": "
           << "EventTagSvc has already been initialized" << endl;
      Warning("EventTagSvc has already been initialized");
   }

   // We have to initialize DatabaseSvc -----------------------------------
   DatabaseSvc* dbs = DatabaseSvc::instance();
   if ( (dbs->GetDBFilePath()).empty() ) {
      // set path to directory with databases:
      dbs->SetDBFilePath(selector->AbsPath("Analysis/DatabaseSvc/dat"));
   }

   // We have to initialize Magnetic field --------------------------------
   MagneticFieldSvc* mf = MagneticFieldSvc::instance();
   if ( (mf->GetPath()).empty() ) {
      // set path to directory with magnetic fields tables
      mf->SetPath(selector->AbsPath("Analysis/MagneticField"));
      mf->UseDBFlag(false); // like in the boss program
      mf->RunMode(3); // like in the boss program
   } else {
      cout << " WARNING:"
           << " MagneticFieldSvc has already been initialized" << endl
           << "                         path = " << mf->GetPath() << endl;
      Warning("MagneticFieldSvc has already been initialized");
   }

   // set path for ParticleID algorithm
   ParticleID* pid = ParticleID::instance();
#if (BOSS_VER < 700)
   pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
#else
   pid->set_path(selector->AbsPath("Analysis/ParticleID"));
#endif

   //--------- Book histograms --------------------------------------------

   hst[1] = new TH1D("All_cuts","selections cuts", 20,-0.5,19.5);

   //  ChargedTracksPiPi:
   hst[11] = new TH1D("Rxy","R_{xy}", 200,-2.,2.);
   hst[12] = new TH1D("Rz","R_{z}", 200,-20.,20.);
   hst[13] = new TH1D("cos_theta","cos(#theta)", 200,-1.,1.);
   hst[14] = new TH1D("theta","#theta", 180,0.,180.);
   hst[15] = new TH1D("Pid_clpi","lg(CL_{#pi})", 100,-4.,0.);
   hst[16] = new TH1D("Pid_ispi","1 - #pi, 0 - another particle", 2,-0.5,1.5);
   hst[17] = new TH1D("noKalTrk", "0,1 - no/yes mdcKalTrk", 2,-0.5,1.5);
   hst[18] = new TH1D("pi_QP","Charged momentum (Q*P)", 200,-1.,1.);
   hst[19] = new TH1D("pi_ct","soft pions cos(#theta)", 200,-1.,1.);
   hst[20] = new TH2D("pi_pm","N(+) vs N(-)", 10,-0.5,9.5, 10,-0.5,9.5);

   hst[21] = new TH1D("Mrec","recol mass of pi+ pi-", 100,3.,3.2);
   hst[22] = new TH2D("MrecNch","Mrec(pi+ pi-) vs Nch",
                      100,3.,3.2,10,1.5,11.5);
   hst[23] = new TH1D("cosPM", "cos(Theta_pi+_pi-)", 100,-1.0,1.0);
   hst[24] = new TH2D("cosPMNch", "cos(Theta_pi+_pi-) vs Nch",
                      100,-1.0,1.0,10,1.5,11.5);
   hst[25] = new TH1D("invPM", "inv.mass(pi+ pi-)", 100,0.25,0.75);
   hst[26] = new TH2D("invPMNch", "inv.mass(pi+ pi-) vs Nch",
                      100,0.25,0.75,10,1.5,11.5);
   hst[27] = new TH1D("NMrec","N good pi+ pi- pairs", 10,-0.5,9.5);

   hst[31] = new TH1D("BMrec","the best Mrec(pi+ pi-)", 100,3.05,3.15);
   hst[32] = new TH1D("BcosPM", "the best cos(Theta_pi+_pi-)", 100,-1.0,1.0);
   hst[33] = new TH1D("BinvPM", "the best inv.mass(pi+ pi-)", 100,0.25,0.75);
   hst[34] = new TH1D("Bnch","number of good trk", 20,-0.5,19.5);

   //  ChargedTracksKK:
   hst[41] = new TH1D("Pid_clK","lg(CL_{K})", 100,-4.,0.);
   hst[42] = new TH1D("Pid_isK","1 - K, 0 - another particle", 2,-0.5,1.5);
   hst[43] = new TH1D("K_QP","Charged momentum (Q*P)for K", 300,-1.5,1.5);
   hst[44] = new TH2D("K_pm","N(K+) vs N(K-)", 10,-0.5,9.5, 10,-0.5,9.5);
   hst[45] = new TH1D("N_K","number of good kaons", 10,-0.5,9.5);
   hst[46] = new TH1D("N_o","number of other good part.", 10,-0.5,9.5);

   // NeutralTracks:
   hst[51] = new TH1D("G_ang","angle with closest chg.trk", 180,0.,180.);
   hst[52] = new TH1D("G_n","N_{#gamma} in event", 11,-0.5,10.5);
   hst[53] = new TH1D("K_PM","momentum K+ and K-", 240,-1.2,1.2);
   hst[54] = new TH1D("G_P","Momentum of gammas", 200,0.,2.);

   // VertKinFit:
   hst[55] = new TH1D("vtx_fit", "vertex fit: 0/1 - bad/good", 2,-0.5,1.5);
   hst[56] = new TH1D("vtx_chi2", "vertex fit #chi^2", 100,0.,100.);

   hst[61] = new TH1D("fit_chi2","fit: min #chi^{2}", 250,0.,250.);
   hst[62] = new TH1D("fit_ch2_3delta",
                      "fit: #chi^{2}(2#gamma) - #chi^{2}(3#gamma)",
                      200,-100.,100.);
   hst[63] = new TH1D("fit_ch2_3","fit: #chi^{2}(3#gamma)", 250,0.,250.);
   hst[64] = new TH1D("fit_ch2_2","fit: #chi^{2}(2#gamma)", 250,0.,250.);

   hst[71] = new TH1D("fit_ang","fit: min angle gamma trk", 180,0.,180.);
   hst[72] = new TH1D("fit_eg","fit: E(#gamma) rejected", 500,0.,0.5);
   hst[73] = new TH1D("fit_egm","fit: E_{max}(#gamma) rejected", 100,0.,1.);

   hst[81] = new TH1D("mgg0","Minv(gg) befor fit", 120,meta-0.06,meta+0.06);
   hst[82] = new TH1D("mkk0","Minv(K+K-) befor fit", 100,0.98,1.1);
   hst[83] = new TH1D("fit_Mrec","fit: Mrec(pi+pi-)", 100,3.,3.2);
   hst[84] = new TH1D("fit_cosPM", "fit: cos Theta(pi+pi-)", 100,-1.0,1.0);
   hst[85] = new TH1D("fit_invPM", "fit: Minv(pi+pi-)", 100,0.25,0.75);
   hst[86] = new TH1D("Cphi_rest","cos #Theta(#phi) in J/Psi", 200,-1.0,1.0);

   //  ChargedTracksPiPi:
   hst[91] = new TH1D("Pip_P","P(#pi+)", 100,0.,0.5);
   hst[92] = new TH1D("Pip_T","Pt(#pi+)", 100,0.,0.5);
   hst[93] = new TH1D("Pip_C","cos #Theta(#pi+)", 200,-1.0,1.0);
   hst[94] = new TH1D("Pim_P","P(#pi-)", 100,0.,0.5);
   hst[95] = new TH1D("Pim_T","Pt(#pi-)", 100,0.,0.5);
   hst[96] = new TH1D("Pim_C","cos #Theta(#pi-)", 200,-1.0,1.0);
   hst[99] = new TH1D("Pi_W","weights for Mr(sig)", 200,0.9,1.1);

   // Monte Carlo histograms:
   hst[100] = new TH1D("mc_dec0", "decPsi(2S) nocut",256,-0.5,255.5);
   hst[101] = new TH1D("mc_dec2", "decPsi(2S) cut#2",256,-0.5,255.5);
   hst[102] = new TH1D("mc_dec4", "decPsi(2S) cut#4",256,-0.5,255.5);
   hst[103] = new TH1D("mc_dec6", "decPsi(2S) cut#6",256,-0.5,255.5);
   hst[104] = new TH1D("mc_dec7", "decPsi(2S) cut#7",256,-0.5,255.5);

   hst[105] = new TH1D("mc_dcj0", "dec J/psi nocut",300,-0.5,299.5);
   hst[106] = new TH1D("mc_dcj2", "dec J/psi cut#2",300,-0.5,299.5);
   hst[107] = new TH1D("mc_dcj4", "dec J/psi cut#4",300,-0.5,299.5);
   hst[108] = new TH1D("mc_dcj6", "dec J/psi cut#6",300,-0.5,299.5);
   hst[109] = new TH1D("mc_dcj7", "dec J/psi cut#7",300,-0.5,299.5);

   // FillHistoMC:
   hst[111] = new TH1D("mc_pdg", "PDG codes of all particles",
                       2001,-1000.5,1000.5);
   hst[112] = new TH1D("mc_pdg0", "PDG of particles from primary vertex",
                       2001,-1000.5,1000.5);

   hst[121] = new TH1D("mc_PsipPt", "Pt of #Psi(2S)", 100,0.,1.);
   hst[122] = new TH1D("mc_PsipC", "cos(#Theta) of Psi(2S)", 100,-1.,1.);

   hst[123] = new TH1D("mc_JpsiP", "Momentum of J/#Psi", 100,0.,1.);
   hst[124] = new TH1D("mc_JpsiPt","Pt of J/#Psi", 100,0.,1.);
   hst[125] = new TH1D("mc_JpsiC", "cos(#Theta) of J/#Psi", 100,-1.,1.);

   hst[126] = new TH1D("mc_PipP", "Momentum of #pi^{+}", 200,0.,1.);
   hst[127] = new TH1D("mc_PipC", "cos(#Theta) of #pi^{+}", 100,-1.,1.);
   hst[128] = new TH1D("mc_PimP", "Momentum of #pi^{-}", 200,0.,1.);
   hst[129] = new TH1D("mc_PimC", "cos(#Theta) of #pi^{-}", 100,-1.,1.);

   hst[131] = new TH1D("mc_EtaP", "Momentum of #eta", 1000,0.,2.);
   hst[132] = new TH1D("mc_EtaPt","Pt of #eta", 1000,0.,2.);
   hst[133] = new TH1D("mc_EtaC", "cos(#Theta) of #eta", 100,-1.,1.);
   hst[135] = new TH1D("mc_PhiP", "Momentum of #phi", 1000,0.,2.);
   hst[136] = new TH1D("mc_PhiPt","Pt of #phi", 1000,0.,2.);
   hst[137] = new TH1D("mc_PhiC", "cos(#Theta) of #phi", 100,-1.,1.);
   hst[138] = new TH1D("mc_KP", "Momentum of K^{#pm}", 1000,-2.,2.);
   hst[139] = new TH1D("mc_KC", "cos(#Theta) of K^{#pm}", 100,-1.,1.);

   hst[141] = new TH1D("mc_EtaRC",
              "cos(#Theta) of #eta in the rest frame of J/#Psi", 100,-1.,1.);
   hst[142] = new TH1D("mc_PhiRC",
              "cos(#Theta) of #phi in the rest frame of J/#Psi", 100,-1.,1.);
   hst[143] = new TH1D("mc_test",
              "angle(#eta,#phi) in the rest frame of J/#Psi", 100,179.,181.);
   hst[145] = new TH1D("mc_Mkk", "Minv(K^{+}K^{-})", 140,0.98,1.12);
   hst[146] = new TH1D("mc_KpCRC",
              "cos(#Theta) of K^{+} in the rest frame of J/#Psi",100,-1.,1.);
   hst[147] = new TH1D("mc_KmCRC",
              "cos(#Theta) of K^{-} in the rest frame of J/#Psi",100,-1.,1.);

   // MatchRecMcTrks:
   hst[151] = new TH1D("mc_mdp_p","REC-MC ok min dP #pi^{+}",100,-0.05,0.05);
   hst[152] = new TH1D("mc_mdp_np","REC-MC no min dP #pi^{+}",100,-0.05,0.05);
   hst[153] = new TH1D("mc_mnPp","no REC-MC #pi^{+} P", 100,0.,0.5);
   hst[154] = new TH1D("mc_mnTp","no REC-MC #pi^{+} cos(#Theta)",100,-1.,1.);
   hst[155] = new TH1D("mc_mdp_m","REC-MC ok min dP #pi^{-}",100,-0.05,0.05);
   hst[156] = new TH1D("mc_mdp_nn","REC-MC no min dP #pi^{-}",100,-0.05,0.05);
   hst[157] = new TH1D("mc_mnPn","no REC-MC #pi^{-} P", 100,0.,0.5);
   hst[158] = new TH1D("mc_mnTn","no REC-MC #pi^{-} cos(#Theta)",100,-1.,1.);

   // MatchMcRecMr:
   hst[161] = new TH1D("mc_dPx1","MC-REC dPx #pi^{+}", 120,-0.06,0.06);
   hst[162] = new TH1D("mc_dPy1","MC-REC dPy #pi^{+}", 120,-0.06,0.06);
   hst[163] = new TH1D("mc_dPz1","MC-REC dPz #pi^{+}", 120,-0.06,0.06);
   hst[164] = new TH1D("mc_dPa1","MC-REC delta angle #pi^{+}", 100,0.,0.1);
   hst[165] = new TH1D("mc_dP1m","MC-REC min dP #pi^{+}", 200,-0.1,0.1);
   hst[166] = new TH1D("mc_dPx2","MC-REC dPx #pi^{-}", 120,-0.06,0.06);
   hst[167] = new TH1D("mc_dPy2","MC-REC dPy #pi^{-}", 120,-0.06,0.06);
   hst[168] = new TH1D("mc_dPz2","MC-REC dPz #pi^{-}", 120,-0.06,0.06);
   hst[169] = new TH1D("mc_dPa2","MC-REC delta angle #pi^{-}", 100,0.,0.1);
   hst[170] = new TH1D("mc_dP2m","MC-REC min dP #pi^{-}", 200,-0.1,0.1);
   hst[171] = new TH1D("mc_mcP1","no MC-REC #pi^{+} P", 100,0.,0.5);
   hst[172] = new TH1D("mc_mcC1","no MC-REC #pi^{+} cos(#Theta)",100,-1.,1.);
   hst[173] = new TH1D("mc_mcP2","no MC-REC #pi^{-} P", 100,0.,0.5);
   hst[174] = new TH1D("mc_mcC2","no MC-REC #pi^{-} cos(#Theta)",100,-1.,1.);
   hst[175] = new TH1D("mc_mch","0(no) 1(yes) 2(select), ", 3,-0.5,2.5);
   hst[177] = new TH1D("mc_Mrs","signal Mrec(pi+pi-)", 100,3.,3.2);
   hst[178] = new TH1D("mc_cosPMs", "signal cos Theta(pi+pi-)", 100,-1.0,1.0);
   hst[179] = new TH1D("mc_invPMs", "signal Minv(pi+pi-)", 100,0.25,0.75);

   // data-EtaEff:
   hst[201] = new TH1D("E_mrec","recol mass of pi+ pi-", 200,3.055,3.145);
   hst[203] = new TH1D("E_M2pi0","M^{2}(2#gamma)", 200,0.,0.04);
   hst[204] = new TH1D("E_M2all","M^{2}(2#gamma)", 500,0.,2.);
   hst[205] = new TH1D("E_Npi0","N(#pi^{0})", 5,-0.5,4.5);
   hst[207] = new TH1D("E_mkk2","Minv^{2}(K+K-)", 170,0.98,1.15);
   hst[211] = new TH1D("E_M2gg","M^{2}(2#gamma)", 200,0.,0.6);
   hst[215] = new TH1D("E_M2fr","M^2(fr) gamma min", 100,0.,0.02);

   hst[221] = new TH1D("E_rE","Eg(pred)/Eg(rec)", 250,0.,2.5);
   hst[222] = new TH1D("E_dTh","#delta#Theta (pred-rec)", 100,0.,20.);
   hst[223] = new TH2D("E_2D","ratE vs dTheta;#Theta",
                      50,0.,20.,50,0.5,2.0);
   hst[225] = new TH1D("E_M2rl_min","M^2(real) min", 200,-0.02,0.02);
   hst[226] = new TH1D("E_Mgg2_rl","M^{2}(2#gamma) real", 200,0.,0.6);

   hst[231] = new TH1D("E_vtx_ch2", "vertex fit #chi^2", 100,0.,100.);
   hst[232] = new TH1D("E_fit_ch2","kin-fit: #chi^{2}", 250,0.,250.);
   hst[233] = new TH1D("E_Mgg_fit","M(2#gamma) fit", 100,0.5,0.6);

   // MC-EtaEff:
   hst[301] = new TH1D("Emc_M2pi0T","M^{2}(2#gamma) T", 200,0.,0.04);
   hst[302] = new TH1D("Emc_M2pi0F","M^{2}(2#gamma) F", 200,0.,0.04);
   hst[303] = new TH1D("Emc_M2pi0TE","M^{2}(2#gamma) TE", 200,0.,0.04);
   hst[304] = new TH1D("Emc_Npi0T","N(#pi^{0}) T", 5,-0.5,4.5);
   hst[305] = new TH1D("Emc_Npi0F","N(#pi^{0}) F", 5,-0.5,4.5);
   hst[306] = new TH1D("Emc_Npi0TE","N(#pi^{0}) TE", 5,-0.5,4.5);
   hst[307] = new TH1D("Emc_dcj", "dec J/psi if Npi0>0",300,-0.5,299.5);
   hst[308] = new TH1D("Emc_mkk2T","Minv^{2}(K+K-) T", 170,0.98,1.15);
   hst[309] = new TH1D("Emc_mkk2F","Minv^{2}(K+K-) F", 170,0.98,1.15);
   hst[311] = new TH1D("Emc_M2ggT","M^{2}(2#gamma) T", 200,0.,0.6);
   hst[312] = new TH1D("Emc_M2ggF","M^{2}(2#gamma) F", 200,0.,0.6);
   hst[313] = new TH1D("Emc_M2ggTE","M^{2}(2#gamma) TE", 200,0.,0.6);
   hst[315] = new TH1D("Emc_M2frT","M^2(fr) min T", 100,0.,0.02);
   hst[316] = new TH1D("Emc_M2frF","M^2(fr) min F", 100,0.,0.02);
   hst[317] = new TH1D("Emc_M2frTE","M^2(fr) min TE", 100,0.,0.02);

   hst[321] = new TH1D("Emc_rET","Eg(pred)/Eg(rec) T", 250,0.,2.5);
   hst[322] = new TH1D("Emc_dThT","#delta#Theta (pred-rec) T", 100,0.,20.);
   hst[323] = new TH1D("Emc_rEF","Eg(pred)/Eg(rec) F", 250,0.,2.5);
   hst[324] = new TH1D("Emc_dThF","#delta#Theta (pred-rec) F", 100,0.,20.);
   hst[325] = new TH2D("Emc_2DT","ratE vs dTheta T;#Theta",
                      50,0.,20.,50,0.5,2.0);
   hst[326] = new TH2D("Emc_2DF","ratE vs dTheta F;#Theta",
                      50,0.,20.,50,0.5,2.0);
   hst[331] = new TH1D("Emc_M2rl_minT","M^2(real) T", 200,-0.02,0.02);
   hst[332] = new TH1D("Emc_M2rl_minF","M^2(real) F", 200,-0.02,0.02);
   hst[335] = new TH1D("Emc_Mgg2_rlT","M^{2}(2#gamma) real T", 200,0.,0.6);
   hst[336] = new TH1D("Emc_Mgg2_rlF","M^{2}(2#gamma) real F", 200,0.,0.6);
   hst[338] = new TH1D("Emc_Mgg_fitT","M(2#gamma) fit T", 100,0.5,0.6);
   hst[339] = new TH1D("Emc_Mgg_fitF","M(2#gamma) fit F", 100,0.5,0.6);

   // TTree for pi+ pi- J/Psi selection:
   m_nt1 = new TTree("nt1","pi+ pi- J/Psi selection");
   m_nt1->Branch("Mrs",&xnt1.Mrs);
   m_nt1->Branch("mdp",&xnt1.Mrs_mindp);
   m_nt1->Branch("MrsW",&xnt1.MrsW);
   m_nt1->Branch("Mrec",&xnt1.Mrec);
   m_nt1->Branch("Mrb",&xnt1.Mrbest);
   m_nt1->Branch("Ptp",&xnt1.Ptp);
   m_nt1->Branch("Cpls",&xnt1.Cpls);
   m_nt1->Branch("Ptm",&xnt1.Ptm);
   m_nt1->Branch("Cmns",&xnt1.Cmns);
   m_nt1->Branch("nch",&xnt1.nch);
   m_nt1->Branch("Nm",&xnt1.Nm);  // just for convenience,
   m_nt1->Branch("dec",&xnt1.decPsip);
   m_nt1->Branch("mcmkk",&xnt1.mcmkk);

   // final ntuple for
   //   Psi(2S) -> pi+ pi- J/Psi
   //                      J/Psi -> phi eta
   //                               phi -> K+ K-
   //                                   eta -> 2gamma
   m_tuple[1] = new TNtupleD("a4c","after 4C kinematic fit",
                "Mrec:"           // recoil mass as in 'nt1'
                "ch2:"            // chi^2 of 5C fit
                "Ptpip:Cpip:"     // Pt and cos(Theta) of pi+
                "Ptpim:Cpim:"     // Pt and cos(Theta) of pi-
                "Ptkp:Ckp:"       // Pt and cos(Theta) of K+
                "Ptkm:Ckm:"       // Pt and cos(Theta) of K-
                "Eg1:Cg1:"        // E and cos(Theta) of gamma-1
                "Eg2:Cg2:"        // E and cos(Theta) of gamma-2
                "Ptgg:Cgg:"       // Pt and cos(Theta) of eta (2gamma)
                "Mkk:Mgg:"        // invariant masses of K+K- and 2gammas
                "M2kpg1:M2kpg2:M2kmg1:M2kmg2:" // M_inv^2( K(+/-)g(1/2) )
                "M2kpeta:M2kmeta"              // M_inv^2( K(+/-) eta )
                ":dec:decj"       // MC: decay codes of Psi(2S) and J/Psi
                ":mcmkk"          // MC: invariant masses of MCtrue K+K-
                            );

   // ntuple for eta efficiency study
   m_tuple[2] = new TNtupleD("eff_eta","eta reconstruction efficiency",
                "Pkk:Ckk:"      // P and cos(Theta) of K+K-
                "Peta:Ceta:"    // P and cos(Theta) of eta
                "Eg1:Cg1:"      // E and cos(Theta) of rec gamma
                "Eg2:Cg2:"      // E and cos(Theta) of rec gamma
                "Egr:Cgr:"      // E and cos(Theta) of found gamma (fl>0)
                "fl:"           // 0/1/2 - only predicted/gamma/eta-found
                "rE:dTh:"       // predicted/found relation
                "m2fr:m2rl:"    // final recoil M^2 predicted/found
                "mggpr:mggf"    // mass(eta) predicted/fitted
                ":decj"         // MC: decay codes of J/Psi
                ":deta"         // MC: 1 if eta decay to 2 gamma
                            );

   // register in selector to save in given directory
   const char* SaveDir = "PsipJpsiPhiEta";
   VecObj hsto(hst.begin(),hst.end());
   selector->RegInDir(hsto,SaveDir);
   VecObj ntuples(m_tuple.begin(),m_tuple.end());
   selector->RegInDir(ntuples,SaveDir);

   VecObj tt;
   tt.push_back(m_nt1);
   selector->RegInDir(tt,SaveDir);
}

//-------------------------------------------------------------------------
static Hep3Vector getVertexOrigin(int runNo, bool verbose = false) {
//-------------------------------------------------------------------------
   static int save_runNo = 0;
   static Hep3Vector xorigin;

   if ( runNo == save_runNo ) {
      return xorigin;
   }

   // update vertex for new run
   xorigin.set(0.,0.,0.);
   VertexDbSvc* vtxsvc = VertexDbSvc::instance();

   int run = abs(runNo);
   if (
      (run >= 8093 && run <= 9025)  // 2009 Psi(2S)
      ||  (run >= 9613 && run <= 9779)  // 2009 3650 data
      ||  (run >= 33725 && run <= 33772)  // 2013 3650 data
   ) {
      vtxsvc->SetBossVer("6.6.4");
   } else if (
      (run >= 25338 && run <= 27090)  // 2012 Psi(2S)
   ) {
      vtxsvc->SetBossVer("6.6.4.p03");
   }
   vtxsvc->handle(runNo);

   if ( vtxsvc->isVertexValid() ) {
      double* dbv = vtxsvc->PrimaryVertex();
      xorigin.set(dbv[0],dbv[1],dbv[2]);
      if( verbose ) {
         cout << " sqlite-db vertex: (x,y,z)= " << xorigin << endl;
      }
   } else {
      cout << " FATAL ERROR:"
           " Cannot obtain vertex information for run#" << runNo << endl;
      exit(1);
   }

   save_runNo = runNo;
   return xorigin;
}

//-------------------------------------------------------------------------
static void FillHistoMC(const ReadDst* selector, Select& Slct) {
//-------------------------------------------------------------------------
   if ( !isMC ) {
      return;
   }

   const TMcEvent*  m_TMcEvent  = selector->GetMcEvent();
   const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol();
   if ( !mcParticles ) {
      return;
   }

   // EventTag
   m_EventTagSvc->setMcEvent(const_cast<TMcEvent*>(m_TMcEvent));
   unsigned int evTag = m_EventTagSvc->getEventTag();
   // check general event type
   if ( (evTag & 0xF) != 5 ) { // 5 is Psi(2S) event
      printf(" WARNING: MC is not Psi(2s) evTag= 0x%08x\n",evTag);
      Warning("MC is not Psi(2s)");
      evTag = 0;
   }

   // decay code of Psi(2S):
   int decPsip  = (evTag >> 8) & 0xFF;
   Slct.decPsip = decPsip;
   hst[100]->Fill(Slct.decPsip);

   // decay code of J/Psi from Psi(2S) -> pi+ pi- J/Psi
   int decJpsi = 0;
   if( decPsip == 64 ) {
      decJpsi = ( evTag >> 16 ) & 0xFF;
      Slct.decJpsi = decJpsi;
      hst[105]->Fill(Slct.decJpsi);
   }

   int idx_psip=-1;
   int idx_jpsi=-1;
   int idx_phi=-1;
   int idx_eta=-1;
   PsipTbl.Reset();
   JpsiTbl.Reset();

   // momentum eta & phi in rest frame of J/Psi
   Hep3Vector beta_jpsi;
   HepLorentzVector LVeta, LVphi;
   // momentum K+ & K- from phi decay
   vector<HepLorentzVector> LVK;

   TIter mcIter(mcParticles);
   while( auto part = static_cast<TMcParticle*>(mcIter.Next()) ) {
      long int part_pdg = part->getParticleID ();

      Hep3Vector Vp( part->getInitialMomentumX(),
                     part->getInitialMomentumY(),
                     part->getInitialMomentumZ() );

      hst[111]->Fill(part_pdg);
      if ( part->getMother() == -99 ) { // primary vertex
         hst[112]->Fill(part_pdg);
      }

      if ( part_pdg == 100443 ) { // Psi(2S)
         hst[121]->Fill(Vp.rho());
         hst[122]->Fill(Vp.cosTheta());
         idx_psip = part->getTrackIndex();
      }

      if ( part->getMother() == idx_psip ) {      // decays of Psi(2S)
         if ( decPsip == 64 ) {                    // -> pi+ pi- J/Psi
            if ( part_pdg == 443 ) {                // J/Psi
               hst[123]->Fill(Vp.mag());
               hst[124]->Fill(Vp.rho());
               hst[125]->Fill(Vp.cosTheta());
               idx_jpsi = part->getTrackIndex();
               HepLorentzVector LVjpsi( Vp, sqrt(Vp.mag2() + SQ(mjpsi)) );
               beta_jpsi = LVjpsi.boostVector();
            } else if ( part_pdg == 211 ) {         // pi+
               hst[126]->Fill(Vp.mag());
               hst[127]->Fill(Vp.cosTheta());
               Slct.mcPip = Vp;
            } else if ( part_pdg == -211 ) {        // pi-
               hst[128]->Fill(Vp.mag());
               hst[129]->Fill(Vp.cosTheta());
               Slct.mcPim = Vp;
            }
         }

         // collect decays of Psi(2S)
         PsipTbl.vecdec.push_back(part_pdg);

      } // end Psi(2S) decays

      if( part->getMother() == idx_jpsi ) {  // decays of J/Psi
         // collect decays of J/psi
         JpsiTbl.vecdec.push_back(part_pdg);
      }

      if ( part->getMother() == idx_jpsi           // Psip -> pi+pi-J/Psi
            && decJpsi == 68                ) {    // J/Psi -> phi eta
         if ( part_pdg == 221 ) {                  // eta
            hst[131]->Fill(Vp.mag());
            hst[132]->Fill(Vp.rho());
            hst[133]->Fill(Vp.cosTheta());
            idx_eta = part->getTrackIndex();

            HepLorentzVector LV( Vp, sqrt(Vp.mag2() + SQ(meta)) );
            LV.boost(-beta_jpsi);
            LVeta=LV;
            hst[141]->Fill( LV.cosTheta() );

         } else if ( part_pdg == 333 ) {           // phi
            hst[135]->Fill(Vp.mag());
            hst[136]->Fill(Vp.rho());
            hst[137]->Fill(Vp.cosTheta());
            idx_phi = part->getTrackIndex();

            HepLorentzVector LV( Vp, sqrt(Vp.mag2() + SQ(mphi)) );
            LV.boost(-beta_jpsi);
            LVphi=LV;
            hst[142]->Fill( LV.cosTheta() );
         }
      }

      if ( part->getMother() == idx_phi ) { // J/Psi -> *phi* eta
                                            //            |-> K+ K-
         if ( abs(part_pdg) == 321 ) {      // K+ or K-
            int signK = ( part_pdg > 0 ) ? 1 : -1;
            hst[138]->Fill( signK*Vp.mag() );
            hst[139]->Fill(Vp.cosTheta());
            LVK.push_back( HepLorentzVector( Vp,sqrt(Vp.mag2()+SQ(mk)) ) );
         }
      }

      if ( part->getMother() == idx_eta ) { // J/Psi -> phi eta
                                            //               |-> 2gamma
         if ( part_pdg != 22 ) {            // it's not 2gamma decay
            Slct.dec_eta = -1;
         }
      }
   } // end of while

   if ( decJpsi == 68 ) {     // J/Psi -> phi eta
      // check angle(phi,eta) in the rest frame of J/Psi
      double alpha= LVeta.angle(LVphi.vect());
      hst[143]->Fill( RtoD(alpha) );

      // study cut-off in invariant mass of K+ K- from phi decays
      if ( LVK.size() == 2 ) {
         Slct.mc_mkk = (LVK[0]+LVK[1]).m();
         hst[145]->Fill( Slct.mc_mkk );
      }

      // set Slct.dec_eta
      if ( Slct.dec_eta == 0 ) {
         Slct.dec_eta = 1;  // eta -> 2 gamma decay
      } else {
         Slct.dec_eta = 0;
      }
   }

   if ( Slct.decPsip == 64 && Slct.decJpsi == 0 ) {
      //  add small corrections to Slct.decJpsi
      string strdec = JpsiTbl.StrDec();
//       cout << " DEBUG: 64 -> strdec= " << strdec << endl;
      if ( !strdec.empty() ) {
         int decJpsi = 0;

         // DEBUG
//          cout << " CHECK-0: JpsiTbl= :" << strdec << ":" << endl;

         // sortied by absolute value of PDG code and positive first
         if ( strdec == string("K+ K- eta") ||
              strdec == string("K+ K- eta -22") ) {
            decJpsi = 261;
         } else if ( strdec == string("K+ K- eta pi+ pi-") ) {
            decJpsi = 262;
         } else if ( strdec == string("K+ K- pi0 pi0") ) {
            decJpsi = 263;
         } else if ( strdec == string("K*+ K- gamma")
                  || strdec == string("K*- K+ gamma") ) {
            decJpsi = 264;
         }

         if ( decJpsi > 0 ) {
            Slct.decJpsi = decJpsi;
         }
      }
   }

//    if ( decPsip == 64 ) { // DEBUG
//       if (  Slct.decJpsi == 118 || Slct.decJpsi == 83 ) {
//          cout << " CHECK: Slct.decJpsi= " << Slct.decJpsi
//               << " JpsiTbl= " << JpsiTbl.StrDec() << endl;
//          selector->PrintMcDecayTree(-99,0);
//       }
//    }

   // fill mc_mkk for J/Psi -> K+ K- eta
   if ( Slct.decPsip == 64 && Slct.decJpsi == 261 ) {
      TIter mcIter(mcParticles);
      while( auto part = static_cast<TMcParticle*>(mcIter.Next()) ) {
         long int part_pdg = part->getParticleID ();
         Hep3Vector Vp( part->getInitialMomentumX(),
                        part->getInitialMomentumY(),
                        part->getInitialMomentumZ() );

         if ( part->getMother() == idx_jpsi ) {
            if ( abs(part_pdg) == 321 ) {      // K+ or K-
               int signK = ( part_pdg > 0 ) ? 1 : -1;
               hst[138]->Fill( signK*Vp.mag() );
               hst[139]->Fill(Vp.cosTheta());
               HepLorentzVector LV( Vp, sqrt(Vp.mag2()+SQ(mk)) );
               LVK.push_back(LV);

               LV.boost(-beta_jpsi);
               if ( signK > 0 ) {
                  hst[146]->Fill( LV.cosTheta() );
               } else {
                  hst[147]->Fill( LV.cosTheta() );
               }

            } else if ( part_pdg == 221 ) {    // eta
               hst[131]->Fill(Vp.mag());
               hst[132]->Fill(Vp.rho());
               hst[133]->Fill(Vp.cosTheta());

               HepLorentzVector LV( Vp, sqrt(Vp.mag2() + SQ(meta)) );
               LV.boost(-beta_jpsi);
               hst[141]->Fill( LV.cosTheta() );
            }
         } // end J/Psi daughters
      } // end of while()

      // invariant mass of K+ K-
      if ( LVK.size() == 2 ) {
         Slct.mc_mkk = (LVK[0]+LVK[1]).m();
      }
   }
}

//-------------------------------------------------------------------------
#if __cplusplus >= 201103L
// C++11 and above
[[gnu::unused]]
#endif
static void MatchRecMcTrks(ReadDst* selector, Select& Slct) {
//-------------------------------------------------------------------------
   // Match charged reconstructed tracks (with the
   // closest to M(J/Psi) recoil mass) with MC particles

   const static double maxdp = 0.100;      // 100MeV

   if ( !isMC ) {
      return;
   }

   const TMcEvent* m_TMcEvent = selector->GetMcEvent();
   const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol();
   if ( !mcParticles ) {
      return;
   }

   Hep3Vector Pip(Slct.trk_Pip->px(), Slct.trk_Pip->py(), Slct.trk_Pip->pz());
   Hep3Vector Pim(Slct.trk_Pim->px(), Slct.trk_Pim->py(), Slct.trk_Pim->pz());
   double Pip_mag = Pip.mag(), Pim_mag = Pim.mag();

   // pdg of the best selected candidates
   long int pip_pdg = 0, pim_pdg = 0;
   int      pip_ok = -1, pim_ok = -1;

   double mindp_pip = maxdp, mindp_pim = maxdp;

   int idx_psip=-1;
   TIter mcIter(mcParticles);
   while( auto part = static_cast<TMcParticle*>(mcIter.Next()) ) {

      long int pdg = part->getParticleID();
      if ( pdg == 100443 ) { // Psi(2S)
         idx_psip = part->getTrackIndex();
         continue;
      }

      TParticlePDG* pdgParticle = TDatabasePDG::Instance()->GetParticle(pdg);
      if( pdgParticle == 0  || pdgParticle->Charge() == 0 ) {
         continue;
      }

      Hep3Vector Vp( part->getInitialMomentumX(),
                     part->getInitialMomentumY(),
                     part->getInitialMomentumZ() );
      double Vp_mag = Vp.mag();

      // positive pion
      if (    fabs(Pip[0]-Vp[0]) < maxdp
              && fabs(Pip[1]-Vp[1]) < maxdp
              && fabs(Pip[2]-Vp[2]) < maxdp ) {
         if ( fabs(Pip_mag - Vp_mag) < fabs(mindp_pip) ) {
            mindp_pip = Pip_mag - Vp_mag;
            pip_pdg = pdg;
            pip_ok = int(part->getMother() == idx_psip);
         }
      }

      // negative pion
      if (    fabs(Pim[0]-Vp[0]) < maxdp
              && fabs(Pim[1]-Vp[1]) < maxdp
              && fabs(Pim[2]-Vp[2]) < maxdp ) {
         if ( fabs(Pim_mag - Vp_mag) < fabs(mindp_pim) ) {
            mindp_pim = Pim_mag - Vp_mag;
            pim_pdg = pdg;
            pim_ok = int(part->getMother() == idx_psip);
         }
      }

   } // end of while

   // check quality of matching
   if ( pip_ok == 1 ) {
      hst[151]->Fill( mindp_pip );
   } else if ( pip_ok == 0 ) {
      hst[152]->Fill( mindp_pip );
   } else {
      hst[153]->Fill(Pip_mag);
      hst[154]->Fill(Pip.cosTheta());
   }
   if ( pim_ok == 1 ) {
      hst[155]->Fill( mindp_pim );
   } else if ( pim_ok == 0 ) {
      hst[156]->Fill( mindp_pim );
   } else {
      hst[157]->Fill(Pim_mag);
      hst[158]->Fill(Pim.cosTheta());
   }

   // save in Select
   Slct.pip_pdg = pip_pdg;
   Slct.pip_ok  = pip_ok;
   Slct.pim_pdg = pim_pdg;
   Slct.pim_ok  = pim_ok;

   return;
}

//-------------------------------------------------------------------------
static void MatchMcRecMr(ReadDst* selector, Select& Slct) {
//-------------------------------------------------------------------------
   // Match MC pi+ pi- from decay of Psi(2S) -> J/Psi pi+ pi-
   // with reconstructed tracks and calculate their recoil mass
   //
   // IMPORTANT: it is necessary to check that this pair of pions
   //            passes all the selection criteria

   const static double maxdp = 0.1;   // 100MeV

   if ( !isMC ) {
      return;
   }
   if ( Slct.decPsip != 64 ) {
      return;
   }

   double mcPip_mag = Slct.mcPip.mag();
   RecMdcKalTrack* trp = nullptr;
   double mindp_pip = maxdp;
   for ( const auto& tr : Slct.good_pip ) {
      Hep3Vector Pip(tr->px(), tr->py(), tr->pz());
      double Pip_mag = Pip.mag();
      hst[161]->Fill( Pip[0]-Slct.mcPip[0] );
      hst[162]->Fill( Pip[1]-Slct.mcPip[1] );
      hst[163]->Fill( Pip[2]-Slct.mcPip[2] );
      if (    fabs(Pip[0]-Slct.mcPip[0]) < maxdp
              && fabs(Pip[1]-Slct.mcPip[1]) < maxdp
              && fabs(Pip[2]-Slct.mcPip[2]) < maxdp ) {
         hst[164]->Fill( Slct.mcPip.angle(Pip) );
         if ( fabs(Pip_mag - mcPip_mag) < fabs(mindp_pip) ) {
            mindp_pip = Pip_mag - mcPip_mag;
            trp   = tr;
         }
      }
   }
   if ( trp ) { // matching for pi+
      hst[165]->Fill( mindp_pip );
   } else {
      hst[171]->Fill( mcPip_mag );
//       hst[172]->Fill( RtoD(Slct.mcPip.theta()) );
      hst[172]->Fill( Slct.mcPip.cosTheta() );
   }

   double mcPim_mag = Slct.mcPim.mag();
   RecMdcKalTrack* trm = nullptr;
   double mindp_pim = maxdp;
   for ( const auto& tr : Slct.good_pim ) {
      Hep3Vector Pim(tr->px(), tr->py(), tr->pz());
      double Pim_mag = Pim.mag();
      hst[166]->Fill( Pim[0]-Slct.mcPim[0] );
      hst[167]->Fill( Pim[1]-Slct.mcPim[1] );
      hst[168]->Fill( Pim[2]-Slct.mcPim[2] );
      if (    fabs(Pim[0]-Slct.mcPim[0]) < maxdp
              && fabs(Pim[1]-Slct.mcPim[1]) < maxdp
              && fabs(Pim[2]-Slct.mcPim[2]) < maxdp ) {
         hst[169]->Fill( Slct.mcPim.angle(Pim) );
         if ( fabs(Pim_mag - mcPim_mag) < fabs(mindp_pim) ) {
            mindp_pim = Pim_mag - mcPim_mag;
            trm   = tr;
         }
      }
   }
   if ( trm ) { // matching for pi-
      hst[170]->Fill( mindp_pim );
   } else {
      hst[173]->Fill( mcPim_mag );
//       hst[174]->Fill( RtoD(Slct.mcPim.theta()) );
      hst[174]->Fill( Slct.mcPim.cosTheta() );
   }

   if ( !trp || !trm ) { // no matching for pair
      hst[175]->Fill( 0. );
      return;
   }
   hst[175]->Fill(1.);

   // pass these pi+ pi- through the selection
   Hep3Vector Vp(trp->px(), trp->py(), trp->pz());;
   HepLorentzVector LVp( Vp, sqrt( Vp.mag2() + SQ(mpi) ) );
   Hep3Vector Vm(trm->px(), trm->py(), trm->pz());;
   HepLorentzVector LVm( Vm, sqrt( Vm.mag2() + SQ(mpi) ) );

   double Mrec2 = (Slct.LVcms - LVp - LVm).mag2();
   if ( Mrec2 < 0 ) {
      return;
   }
   double Mrec = sqrt(Mrec2);
   if( Mrec <= 3.0 || Mrec >= 3.2 ) {
      return;
   }

   double cosPM = Vm.cosTheta(Vp); // Theta(pi+ pi-)
   double invPM = (LVp+LVm).m();   // Minv( pi+ pi-)

   hst[177]->Fill( Mrec );
   hst[178]->Fill( cosPM );
   hst[179]->Fill( invPM );

   // check Mrec candidate
   if ( !SelectPM( cosPM, invPM ) ) {
      return;
   }

   Slct.Mrec_sig = Mrec;
   Slct.Mrec_sig_mindp = max(fabs(mindp_pip),fabs(mindp_pim));
   Slct.Pip_sig = trp;
   Slct.Pim_sig = trm;

   hst[175]->Fill(2.);
}

//-------------------------------------------------------------------------
// Charged tracks
//-------------------------------------------------------------------------

// search for pi+ pi- J/Psi
//-------------------------------------------------------------------------
static bool ChargedTracksPiPi(ReadDst* selector, Select& Slct) {
//-------------------------------------------------------------------------
   static const double Rvxy0_max = 1.0;
   static const double Rvz0_max = 10.0;
   static const double cosTheta_max = 0.80;  // barrel only
//    static const double cosTheta_max = 0.93;

   const TEvtRecObject* m_TEvtRecObject = selector->GetEvtRecObject();
   const TEvtRecEvent* evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();
   ParticleID* pid = ParticleID::instance();

   // 1) select soft pion candidates (check PID information)

   int Ncharged = 0; // counter for all good charged tracks

   // temporary vectors for soft pions
   vector<RecMdcKalTrack*> trk_p; // plus
   vector<RecMdcKalTrack*> trk_m; // minus

   for(int i = 0; i < evtRecEvent->totalCharged(); i++) {

      DstEvtRecTracks* itTrk =
         static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(i));
      if( !itTrk->isMdcTrackValid() ) {
         continue;
      }

      RecMdcTrack* mdcTrk = itTrk->mdcTrack();

      double theta = mdcTrk->theta();
      double cosTheta = cos(theta);

      HepVector a = mdcTrk->helix();
      HepSymMatrix Ea = mdcTrk->err();
      HepPoint3D point0(0.,0.,0.);   // initial point for MDC recosntruction
      HepPoint3D IP(Slct.xorig[0],Slct.xorig[1],Slct.xorig[2]);
      VFHelix helixip(point0,a,Ea);
      helixip.pivot(IP);
      HepVector vecipa = helixip.a();
      double Rvxy0 = vecipa[0]; // the nearest distance to IP in xy plane
      double Rvz0  = vecipa[3]; // ... in z direction

      hst[11]->Fill(Rvxy0);
      hst[12]->Fill(Rvz0);
      hst[13]->Fill(cosTheta);
      hst[14]->Fill(RtoD(theta));

      if( fabs(Rvxy0) >= Rvxy0_max ) {
         continue;
      }
      if( fabs(Rvz0) >= Rvz0_max ) {
         continue;
      }
      if ( fabs(cosTheta) >= cosTheta_max ) {
         continue;
      }

      Ncharged++;

      // PID information
      pid->init();
      pid->setMethod(pid->methodProbability());
      pid->setChiMinCut(4);
      pid->setRecTrack(itTrk);

      // list of systems used for PID
      pid->usePidSys( pid->useDedx() |
                      pid->useTof1() | pid->useTof2() |
                      pid->useTofE() );

      pid->identify(
         pid->onlyPion()    |
         pid->onlyKaon()    |
         pid->onlyProton()  |
         pid->onlyMuon()    |
         pid->onlyElectron()
      );

      pid->calculate(Slct.runNo);
      int check_pion = 0;
      if ( pid->IsPidInfoValid() ) {
         if( pid->probPion() > 0 ) {
            hst[15]->Fill( log10(pid->probPion()) );
         } else {
            hst[15]->Fill(-5.);
         }

         // check that track is pion:
         if ( (pid->probPion() > 0.001)             &&
               (pid->probPion() > pid->probKaon())   &&
               (pid->probPion() > pid->probProton())
            ) {
            check_pion = 1;
         }
      } // end IsPidInfoValid
      hst[16]->Fill( double(check_pion) );
      if( check_pion == 0 ) {
         continue;
      }

      // require Kalman fit
      RecMdcKalTrack* mdcKalTrk = itTrk->mdcKalTrack();
      hst[17]->Fill( double(mdcKalTrk != 0) );
      if ( !mdcKalTrk ) {
         cout << " WARNING: No Kalman track! EventNo= "
              << Slct.event << endl;
         Warning("No Kalman track");
         continue;
      }
      mdcKalTrk->setPidType(RecMdcKalTrack::pion);

      hst[18]->Fill( mdcKalTrk->charge()*mdcKalTrk->p() );
      hst[19]->Fill( cos(mdcKalTrk->theta()) );
      if( mdcKalTrk->p() > 0.45 ) {
         continue;
      }

      if( mdcKalTrk->charge() > 0 ) {
         trk_p.push_back(mdcKalTrk);
      } else {
         trk_m.push_back(mdcKalTrk);
      }
   } //-------------------------------------------------End for 1)

   int Ntrkp = trk_p.size();
   int Ntrkm = trk_m.size();
   hst[20]->Fill(Ntrkp,Ntrkm);
   if ( Ntrkp < 1 || Ntrkm < 1 ) {
      return false;
   }
   hst[1]->Fill(1); // "cuts"


   // 2) calculate recoil mass of all pairs of pi+ pi-
   double delta_min = 1e3;
   for ( const auto& tp : trk_p ) {
      Hep3Vector Vp(tp->px(), tp->py(), tp->pz());
      HepLorentzVector LVp( Vp, sqrt( Vp.mag2() + SQ(mpi) ) );

      for ( const auto& tm : trk_m ) {
         Hep3Vector Vm(tm->px(), tm->py(), tm->pz());
         HepLorentzVector LVm( Vm, sqrt( Vm.mag2() + SQ(mpi) ) );

         // Mrec(pi+ pi-)
         double Mrec2 = (Slct.LVcms - LVp - LVm).m2();
         if( Mrec2 < 0 ) {
            continue;
         }
         double Mrec = sqrt(Mrec2);
         if( Mrec <= 3.0 || Mrec >= 3.2 ) { // preselection
            continue;
         }
         hst[21]->Fill( Mrec );
         hst[22]->Fill( Mrec, double(Ncharged) );

         // Theta(pi+ pi-)
         double cosPM = Vm.cosTheta(Vp);
         hst[23]->Fill( cosPM );
         hst[24]->Fill( cosPM, double(Ncharged) );

         // Minv( pi+ pi-)
         double invPM = (LVp+LVm).m();
         hst[25]->Fill( invPM );
         hst[26]->Fill( invPM, double(Ncharged) );

         // check Mrec candidate
         if ( !SelectPM( cosPM, invPM ) ) {
            continue;
         }

         // good candidate
         Slct.Mrec.push_back(Mrec);
         Slct.good_pip.insert(tp);
         Slct.good_pim.insert(tm);

         // the best candidate (closest to M(J/Psi))
         double delta = fabs(Mrec-mjpsi);
         if( delta < delta_min ) {
            delta_min = delta;
            Slct.Mrec_best = Mrec;
            Slct.trk_Pip = tp;
            Slct.trk_Pim = tm;
         }
      }
   }  //-------------------------------------------------End for 2)

   hst[27]->Fill( Slct.Mrec.size() );
   if( Slct.Mrec.size() == 0 ) {
      return false;
   }
   hst[1]->Fill(2); // "cuts"

   if ( isMC ) {
      hst[101]->Fill(Slct.decPsip);
      hst[106]->Fill(Slct.decJpsi);

      // fill decay table of Psi(2S):
      PsipTbl2.vecdec = PsipTbl.vecdec; // propagate for second table
      if ( fabs(Slct.Mrec_best - 3.097) < 0.005 ) { // [3.092, 3.102]
         PsipTbl.Add();
      }
      if ( fabs(Slct.Mrec_best - 3.1) < 0.045 ) { // [3.055, 3.145]
         PsipTbl2.Add();
      }

      // fill MC related fields in Slct
//       MatchRecMcTrks(selector, Slct);
      MatchMcRecMr(selector, Slct);
   }

   //----------------------------------------------------------------------
   // histo and ntuple for pi+ pi- J/Psi
   hst[31]->Fill( Slct.Mrec_best );

   Hep3Vector Vp(Slct.trk_Pip->px(), Slct.trk_Pip->py(), Slct.trk_Pip->pz());
   Hep3Vector Vm(Slct.trk_Pim->px(), Slct.trk_Pim->py(), Slct.trk_Pim->pz());
   Slct.cosPM = Vm.cosTheta(Vp);
   hst[32]->Fill( Slct.cosPM );

   HepLorentzVector LVp( Vp, sqrt( Vp.mag2() + SQ(mpi) ) );
   HepLorentzVector LVm( Vm, sqrt( Vm.mag2() + SQ(mpi) ) );
   Slct.invPM = (LVp+LVm).m();
   hst[33]->Fill( Slct.invPM );

   hst[34]->Fill( double(Ncharged) );

   for ( const auto& tp : Slct.good_pip ) {
      Hep3Vector Vp(tp->px(), tp->py(), tp->pz());
      hst[91]->Fill(Vp.mag());
      hst[92]->Fill(Vp.perp());
      hst[93]->Fill(Vp.cosTheta());
   }
   for ( const auto& tm : Slct.good_pim ) {
      Hep3Vector Vm(tm->px(), tm->py(), tm->pz());
      hst[94]->Fill(Vm.mag());
      hst[95]->Fill(Vm.perp());
      hst[96]->Fill(Vm.cosTheta());
   }

   // signal (only MC)
   xnt1.Mrs = Slct.Mrec_sig;
   xnt1.Mrs_mindp = Slct.Mrec_sig_mindp;
   if ( Slct.Pip_sig && Slct.Pim_sig ) {
      double ptp = Slct.Pip_sig->pxy();
      double wp  = ReWeightTrkPid(0,ptp);
      hst[99]->Fill(wp);
      double ptm = Slct.Pim_sig->pxy();
      double wm  = ReWeightTrkPid(0,ptm);
      hst[99]->Fill(wm);
      xnt1.MrsW = wp*wm;
   } else {
      xnt1.MrsW = 0;
   }

   // all Mrec, but in case of MC do not save Mrs in xnt1.Mrec
   bool remove_sig = isMC;
   xnt1.Mrec.clear();
   for ( const auto& mr : Slct.Mrec ) {
      if ( remove_sig ) {
         if ( fabs(mr - Slct.Mrec_sig) < 1e-6 ) {
            remove_sig = false;
            continue;
         }
      }
      xnt1.Mrec.push_back( mr );
   }
   xnt1.Nm = xnt1.Mrec.size();

   // the best: Mrec and momentums
   xnt1.Mrbest = Slct.Mrec_best;
   xnt1.Ptp = Vp.perp();
   xnt1.Cpls = Vp.cosTheta();
   xnt1.Ptm = Vm.perp();
   xnt1.Cmns = Vm.cosTheta();

   xnt1.nch = Ncharged;

   // MC decay
   xnt1.decPsip =
      (Slct.decPsip > 0) ? static_cast<UShort_t>(Slct.decPsip) : 0;
   xnt1.mcmkk = Slct.mc_mkk;

//    m_nt1->Fill(); // ATTENTION!

   return true;
}

// select K+K- candidates, no other tracks.
//-------------------------------------------------------------------------
static bool ChargedTracksKK(ReadDst* selector, Select& Slct) {
//-------------------------------------------------------------------------
   static const double Rvxy0_max = 1.0;
   static const double Rvz0_max = 10.0;
   static const double cosTheta_max = 0.80;  // barrel only
//    static const double cosTheta_max = 0.93;

   if( Slct.Mrec.size() != 1 ) {
      return false;
   }
   int Nother = 0; // other good tracks (not K; not pi+ pi- sel. pair)

   const TEvtRecObject* m_TEvtRecObject = selector->GetEvtRecObject();
   const TEvtRecEvent* evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();
   ParticleID* pid = ParticleID::instance();

   for(int i = 0; i < evtRecEvent->totalCharged(); i++) {

      DstEvtRecTracks* itTrk =
         static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(i));
      if( !itTrk->isMdcTrackValid() ) {
         continue;
      }

      RecMdcTrack* mdcTrk = itTrk->mdcTrack();

      double theta = mdcTrk->theta();
      double cosTheta = cos(theta);

      HepVector a = mdcTrk->helix();
      HepSymMatrix Ea = mdcTrk->err();
      HepPoint3D point0(0.,0.,0.);   // initial point for MDC reconstruction
      HepPoint3D IP(Slct.xorig[0],Slct.xorig[1],Slct.xorig[2]);
      VFHelix helixip(point0,a,Ea);
      helixip.pivot(IP);
      HepVector vecipa = helixip.a();
      double Rvxy0 = vecipa[0]; // the nearest distance to IP in xy plane
      double Rvz0  = vecipa[3]; // ... in z direction

      if( fabs(Rvxy0) >= Rvxy0_max ) {
         continue;
      }
      if( fabs(Rvz0) >= Rvz0_max ) {
         continue;
      }
      if ( fabs(cosTheta) >= cosTheta_max ) {
         continue;
      }

      // require Kalman fit
      RecMdcKalTrack* mdcKalTrk = itTrk->mdcKalTrack();
      if ( !mdcKalTrk ) {
         continue;
      }
      // skip pi+ pi- candiadate
      if( mdcKalTrk == Slct.trk_Pip || mdcKalTrk == Slct.trk_Pim ) {
         continue;
      }

      // PID information
      pid->init();
      pid->setMethod(pid->methodProbability());
      pid->setChiMinCut(4);
      pid->setRecTrack(itTrk);

      // list of systems used for PID
      pid->usePidSys( pid->useDedx() |
                      pid->useTof1() | pid->useTof2() |
                      pid->useTofE() );

      pid->identify(
         pid->onlyPion()    |
         pid->onlyKaon()    |
         pid->onlyProton()  |
         pid->onlyMuon()    |
         pid->onlyElectron()
      );

      pid->calculate(Slct.runNo);
      int check_kaon = 0;
      if ( pid->IsPidInfoValid() ) {
         if( pid->probKaon() > 0 ) {
            hst[41]->Fill( log10(pid->probKaon()) );
         } else {
            hst[41]->Fill(-5.);
         }

         // check that track is kaon:
         if ( (pid->probKaon() > 0.001)             &&
               (pid->probKaon() > pid->probPion())   &&
               (pid->probKaon() > pid->probProton())
            ) {
            check_kaon = 1;
         }
      } // end IsPidInfoValid
      hst[42]->Fill( double(check_kaon) );
      if( check_kaon == 0 ) {
         Nother++;
         continue;
      }

      mdcKalTrk->setPidType(RecMdcKalTrack::kaon);
      hst[43]->Fill( mdcKalTrk->charge()*mdcKalTrk->p() );
      if( mdcKalTrk->charge() > 0 ) {
         Slct.trk_Kp.push_back(mdcKalTrk);
      } else {
         Slct.trk_Km.push_back(mdcKalTrk);
      }
   } //-------------------------------------------------End for(i)

   int np = Slct.trk_Kp.size();
   int nm = Slct.trk_Km.size();
   hst[44]->Fill(np,nm);
   hst[45]->Fill(np+nm);
   hst[46]->Fill(Nother);

   // require exactly one "+" and one "-"
   if ( np != 1 || nm != 1 ) {
      return false;
   }
   // no other good tracks from primary vertex
   if ( Nother > 0 ) {
      return false;
   }
   hst[1]->Fill(3); // "cuts"

   return true;
}

// select gammas candidates
//-------------------------------------------------------------------------
static int NeutralTracks(ReadDst* selector, Select& Slct) {
//-------------------------------------------------------------------------
   // parameters of reconstruction
   static const double min_angle = 10 * M_PI/180; // 10 grad

   const TEvtRecObject* m_TEvtRecObject = selector->GetEvtRecObject();
   const TEvtRecEvent* evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

   for ( int i = evtRecEvent->totalCharged();
         i < evtRecEvent->totalTracks(); i++ ) {
      DstEvtRecTracks* itTrk =
         static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(i));
      if ( !itTrk->isEmcShowerValid() ) {
         continue;
      }
      RecEmcShower* emcTrk = itTrk->emcShower();

      // *) good EMC time:
      if ( emcTrk->time() < 0 || emcTrk->time() > 14 ) {
         continue;
      }

      // *) good EMC energy deposited in the barrel (endcap) part of EMC
      double eraw = emcTrk->energy();
      double absCosTheta = fabs(  cos(emcTrk->theta()) );

      bool GoodCluster=false;
      if ( absCosTheta < 0.8 ) {  //barrel
         GoodCluster = eraw > 25E-3;
      } else if ( absCosTheta > 0.85 && absCosTheta < 0.92 ) { //endcap
         GoodCluster = eraw > 50E-3;
      }
      if ( !GoodCluster ) {
         continue;
      }

      // *) the nearest charged track is far from cluster
      Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

      double tang = 200.; // min angle between cluster and track
      for ( int j = 0; j < evtRecEvent->totalCharged(); j++ ) {
         DstEvtRecTracks* jtTrk =
            static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(j));
         if ( !jtTrk->isExtTrackValid() ) {
            continue;
         }
         RecExtTrack* extTrk = jtTrk->extTrack();
         if ( extTrk->emcVolumeNumber() == -1 ) {
            continue;   //track does not hit EMC
         }
         Hep3Vector extpos = extTrk->emcPosition();
         double angd = fabs(extpos.angle(emcpos)); // [0,pi]
         if ( angd < tang ) {
            tang = angd;
         }
      } //--------------------------------------------End for(j)

      hst[51]->Fill(RtoD(tang));
      if ( tang < min_angle ) {
         continue;
      }

      Slct.gtrk.push_back(emcTrk);
      Slct.angt.push_back(tang);

      // (E,Px,Py,Pz) for good gammas
      Hep3Vector p3 = emcpos - Slct.xorig;
      p3 *= eraw / p3.mag();
      Slct.Pg.push_back( HepLorentzVector(p3,eraw) );
   } //-----------------------------------------------End for (i)

   int ng = Slct.gtrk.size();
   hst[52]->Fill(ng);

   if ( ng < 2 ) {
      return ng;
   }

   hst[1]->Fill(4); // "cuts"
   if ( isMC ) {
      hst[102]->Fill(Slct.decPsip);
      hst[107]->Fill(Slct.decJpsi);
   }

   // momentum of selected pair of Kaons and Gammas
   hst[53]->Fill( Slct.trk_Kp[0]->p() );
   hst[53]->Fill( -Slct.trk_Km[0]->p() );
   for(int i = 0; i < ng; i++) {
      hst[54]->Fill(Slct.Pg[i].e());
   }

   return ng;
}

// Vertex & Kinematic Fit
//-------------------------------------------------------------------------
static bool VertKinFit(ReadDst* selector, Select& Slct) {
//-------------------------------------------------------------------------
   // parameters of reconstruction
   bool do5C = true; // 5C (true) or 4C (false) fit
//    bool do5C = false; // 4C fit to check Psip -> pi+pi-phi eta (no J/Psi)

   if ( !do5C ) {
      Warning("VertKinFit: 4C fitN");
   }

   // search for a good vertex:
   WTrackParameter wp[4] = {
      WTrackParameter(mpi,Slct.trk_Pip->getZHelix(),
                      Slct.trk_Pip->getZError() ),
      WTrackParameter(mpi,Slct.trk_Pim->getZHelix(),
                      Slct.trk_Pim->getZError() ),
      WTrackParameter(mk, Slct.trk_Kp[0]->getZHelix(),
                      Slct.trk_Kp[0]->getZError() ),
      WTrackParameter(mk, Slct.trk_Km[0]->getZHelix(),
                      Slct.trk_Km[0]->getZError() )
   };

   HepPoint3D vx(0., 0., 0.);
   HepSymMatrix Evx(3, 0);
   double bx = 1E+6;
   double by = 1E+6;
   double bz = 1E+6;
   Evx[0][0] = bx*bx;
   Evx[1][1] = by*by;
   Evx[2][2] = bz*bz;

   VertexParameter vxpar;
   vxpar.setVx(vx);
   vxpar.setEvx(Evx);

   VertexFit* vtxfit = VertexFit::instance();
   vtxfit->init();
   for(int k = 0; k < 4; k++) {
      vtxfit->AddTrack(k,wp[k]);
   }
   vtxfit->AddVertex(0, vxpar,0, 1);
   bool ok = vtxfit->Fit(0);

   hst[55]->Fill( double(ok) );
   if ( !ok ) {
      return false;
   }
   hst[1]->Fill(5); // "cuts"
   hst[56]->Fill(vtxfit->chisq(0));

   // 6) 4C or 5C kinematic fit: select two gammas with the best chisq
   vtxfit->BuildVirtualParticle(0);
   vtxfit->Swim(0); // translate parameters of tracks to the vertex# 0
   for(int k = 0; k < 4; k++) {
      wp[k] = vtxfit->wtrk(k);
   }

   KinematicFit* kmfit = KinematicFit::instance();

   double chisq = 9999.;
   int ng = Slct.gtrk.size();
   for(int i = 0; i < ng-1; i++) {
      RecEmcShower* g1Trk = Slct.gtrk[i];
      for(int j = i+1; j < ng; j++) {
         RecEmcShower* g2Trk = Slct.gtrk[j];

         kmfit->init();
         for(int k = 0; k < 4; k++) {
            kmfit->AddTrack(k, wp[k]);
         }
         kmfit->AddTrack(4, 0.0, g1Trk);
         kmfit->AddTrack(5, 0.0, g2Trk);
         if ( do5C ) {
            kmfit->AddResonance(0, mjpsi, 2, 3, 4, 5);
            kmfit->AddFourMomentum(1, Slct.LVcms);
            if ( !kmfit->Fit(0) ) {
               continue;
            }
            if ( !kmfit->Fit(1) ) {
               continue;
            }
         } else {
            kmfit->AddFourMomentum(0, Slct.LVcms);
         }

         bool oksq = kmfit->Fit();
         if ( oksq ) {
            double chi2 = kmfit->chisq();
            if ( chi2 < chisq && chi2 > 0. ) {
               chisq = chi2;
               Slct.g4f[0] = g1Trk;
               Slct.g4f[1] = g2Trk;
            }
         }
      } //-----------------------------------------------End for(j)
   } //-----------------------------------------------End for(i)

   hst[61]->Fill(chisq);
   if ( !Slct.g4f[0] ) {
      return false;   // can not find good pair of gammas
   }

   hst[1]->Fill(6); // "cuts"
   if ( isMC ) {
      hst[103]->Fill(Slct.decPsip);
      hst[108]->Fill(Slct.decJpsi);
   }

   // ++ check that one more gamma does not give better fit ++
   bool f3good = false;
   for(int i = 0; i < ng; i++ ) {
      auto gt = Slct.gtrk[i];
      if( gt == Slct.g4f[0] || gt == Slct.g4f[1] ) {
         continue;
      }
      kmfit->init();
      for(int k = 0; k < 4; k++) {
         kmfit->AddTrack(k, wp[k]);
      }
      kmfit->AddTrack(4, 0.0, Slct.g4f[0]);
      kmfit->AddTrack(5, 0.0, Slct.g4f[1]);
      kmfit->AddTrack(6, 0.0, gt); // additional gamma
      if ( do5C ) {
         kmfit->AddResonance(0, mjpsi, 2, 3, 4, 5, 6);
         kmfit->AddFourMomentum(1, Slct.LVcms);
         if ( !kmfit->Fit(0) ) {
            continue;
         }
         if ( !kmfit->Fit(1) ) {
            continue;
         }
      } else {
         kmfit->AddFourMomentum(0, Slct.LVcms);
      }
      bool oksq = kmfit->Fit();
      if ( oksq ) {
         hst[62]->Fill(chisq - kmfit->chisq());
      }
      if ( !oksq || kmfit->chisq() > chisq ) {
         continue;
      }
      f3good = true;
      hst[63]->Fill(kmfit->chisq());
      hst[64]->Fill( chisq );
   } // end for( additional gamma )
   if ( f3good ) {
      return false;
   }

   hst[1]->Fill(7); // "cuts"
   if ( isMC ) {
      hst[104]->Fill(Slct.decPsip);
      hst[109]->Fill(Slct.decJpsi);
   }

   // ++ repeat fit for two best photons ++
   kmfit->init();
   for ( int k = 0; k < 4; k++ ) {
      kmfit->AddTrack(k, wp[k]);
   }
   kmfit->AddTrack(4, 0.0, Slct.g4f[0]);
   kmfit->AddTrack(5, 0.0, Slct.g4f[1]);
   if ( do5C ) {
      kmfit->AddResonance(0, mjpsi, 2, 3, 4, 5);
      kmfit->AddFourMomentum(1, Slct.LVcms);
   } else {
      kmfit->AddFourMomentum(0, Slct.LVcms);
   }
   bool oksq = kmfit->Fit();
   // check correctness of second fit
   if ( !oksq || fabs(kmfit->chisq()-chisq) > 0.1 ) {
      cout << " WARNING: Bad second fit";
      if ( oksq ) {
         cout << "chisq= "<<chisq<<" new fit chisq= "<<kmfit->chisq();
      }
      cout << endl;
      Warning("Bad second fit");
      return false;
   }
   chisq=kmfit->chisq();

   //----------------------------------------------------------------------
   // Final histograms and ntuples
   //----------------------------------------------------------------------
   double Eg_max = 0.; // max momentum of gammas not selected by fit
   double mang = 200.; // min angle wrt of charge trk for selected gammas
   HepLorentzVector Pgg0; // sum of momentums of good photons
   for(int i = 0; i < ng; i++ ) {
      auto gt = Slct.gtrk[i];
      if( gt == Slct.g4f[0] || gt == Slct.g4f[1] ) {
         // two selected photons
         Pgg0 += Slct.Pg[i];
         double deg = RtoD(Slct.angt[i]);
         hst[71]->Fill(deg);
         if( deg < mang ) {
            mang = deg;
         }
         continue;
      }
      // rejected photons
      double Eg = Slct.Pg[i].e();
      hst[72]->Fill(Eg);
      if( Eg > Eg_max ) {
         Eg_max = Eg;
      }
   }
   if ( Eg_max > 0.01 ) {
      hst[73]->Fill(Eg_max);
   }

   // invariant mass Mgg0 and Mkk0 using momentums before kin.fit.
   double Mgg0 = Pgg0.m();
   Hep3Vector p3p(Slct.trk_Kp[0]->px(),
                  Slct.trk_Kp[0]->py(),
                  Slct.trk_Kp[0]->pz());
   HepLorentzVector Pkp0( p3p, sqrt(p3p.mag2() + SQ(mk)) );
   Hep3Vector p3m(Slct.trk_Km[0]->px(),
                  Slct.trk_Km[0]->py(),
                  Slct.trk_Km[0]->pz());
   HepLorentzVector Pkm0( p3m, sqrt(p3m.mag2() + SQ(mk)) );
   double Mkk0 = (Pkp0+Pkm0).m();
   if ( chisq < 80 ) {
      if ( Mkk0 > 2*mk && Mkk0 < 1.1 ) {
         hst[81]->Fill(Mgg0);
      }
      if ( fabs(Mgg0-meta) < 0.06 ) {
         hst[82]->Fill(Mkk0);
      }
   }

   // Momentums after kinematic corrections:
   // pi+ and pi-
   HepLorentzVector Ppip = kmfit->pfit(0);
   HepLorentzVector Ppim = kmfit->pfit(1);
   double cosPMfit = Hep3Vector(Ppim).cosTheta(Hep3Vector(Ppip));
   double invPMfit = (Ppip+Ppim).m();
   double Mrec_fit = (Slct.LVcms - Ppip - Ppim).m();
   hst[83]->Fill(Mrec_fit);
   hst[84]->Fill(cosPMfit);
   hst[85]->Fill(invPMfit);

   // K+ and K-
   HepLorentzVector Pkp = kmfit->pfit(2);
   HepLorentzVector Pkm = kmfit->pfit(3);
   HepLorentzVector Pkk = Pkp + Pkm;
   double Mkk  = Pkk.m();

   // two gammas
   HepLorentzVector Pg1 = kmfit->pfit(4);
   HepLorentzVector Pg2 = kmfit->pfit(5);
   HepLorentzVector Pgg = Pg1 + Pg2;
   double Mgg  = Pgg.m();

   // J/Psi
   HepLorentzVector Pjpsi = Pkk + Pgg;
   Hep3Vector beta_jpsi = Pjpsi.boostVector();
   HepLorentzVector Pphi = Pkk;
   Pphi.boost(-beta_jpsi); // momentum phi in J/Psi rest system
   if ( chisq < 80 &&
         ( Mkk > 2*mk && Mkk < 1.08 ) &&
         fabs(Mgg-meta) < 0.024
      ) {
      hst[86]->Fill( Pphi.cosTheta() );
   }

   // M^2(Kg)
   double M2kpg1 = (Pkp + Pg1).m2();
   double M2kpg2 = (Pkp + Pg2).m2();
   double M2kmg1 = (Pkm + Pg1).m2();
   double M2kmg2 = (Pkm + Pg2).m2();
   // M^2(K eta)
   double M2kpeta = (Pkp + Pgg).m2();
   double M2kmeta = (Pkm + Pgg).m2();

   // HepLorentzVector.perp == HepLorentzVector.vect.rho
   // HepLorentzVector.rho == HepLorentzVector.vect.mag
   Double_t xfill[] = {
      Slct.Mrec_best,
      chisq,
      Ppip.perp(), Ppip.cosTheta(),
      Ppim.perp(), Ppim.cosTheta(),
      Pkp.perp(), Pkp.cosTheta(),
      Pkm.perp(), Pkm.cosTheta(),
      Pg1.e(), Pg1.cosTheta(),
      Pg2.e(), Pg2.cosTheta(),
      Pgg.perp(), Pgg.cosTheta(),
      Mkk, Mgg,
      M2kpg1,M2kpg2, M2kmg1,M2kmg2, M2kpeta,M2kmeta,
      double(Slct.decPsip), double(Slct.decJpsi),
      Slct.mc_mkk
   };
//       Mrec_fit,cosPMfit,invPMfit,
//       Mkk0, Mgg0,
//       Pphi.cosTheta(),
//       mang,
   m_tuple[1]->Fill( xfill );

   if ( isMC ) {
      static const double seta = 0.008;

      // fill decay table of Jpsi
      JpsiTbl2.vecdec = JpsiTbl.vecdec; // propagate for second table
      if ( chisq < 80
            && Mkk > 2*mk && Mkk < 1.08
            && fabs(Mgg-meta) < 3*seta
         ) {
         if ( fabs(Slct.Mrec_best - 3.097) < 0.005 ) { // [3.092, 3.102]
            JpsiTbl.Add();
         }
         if ( fabs(Slct.Mrec_best - 3.1) < 0.045 ) { // [3.055, 3.145]
            JpsiTbl2.Add();
         }
      }

      if ( !do5C ) { // side-band for 4C
         if (     (Slct.Mrec_best>3.0 && Slct.Mrec_best < 3.05)
               || (Slct.Mrec_best>3.144 && Slct.Mrec_best < 3.194) ) {
            if ( chisq < 80
               && fabs(Mgg-meta) < 3*seta
               && Mkk > 2*mk && Mkk < 1.08      ) {

               // fill decay table of Psi(2S):
               PsipTbl4C.vecdec = PsipTbl.vecdec; // propagate for 4C
               PsipTbl4C.Add();
               cout << "<<< EVENT IN SIDE-BAND>>>"
                    << " Mrec= " << Slct.Mrec_best
                    << " Mgg= " << Mgg
                    << " Mkk= " << Mkk
                    << endl;
               selector->PrintMcDecayTree(-99,0);
            }
         }
      }
   } // end if (isMC)

   return true;
}

//-------------------------------------------------------------------------
static void EtaEff(ReadDst* selector, Select& Slct) {
//-------------------------------------------------------------------------
   // window for selection of eta: see cuts.h (in PsipJpsiPhiEta)
   static const double seta = 0.008;
   static const double weta = 3*seta; // standard

   // cut for Mrec
   double Mrec = Slct.Mrec_best;
   hst[201]->Fill( Mrec );
   if( Mrec <= 3.092 || Mrec >= 3.102 ) { // narrow CUT on M(J/Psi)
//    if( Mrec <= 3.055 || Mrec >= 3.145 ) { // wide CUT on M(J/Psi)
      return;
   }
   Hep3Vector pi3p(Slct.trk_Pip->px(),
                   Slct.trk_Pip->py(),
                   Slct.trk_Pip->pz());
   HepLorentzVector LVpip( pi3p, sqrt(pi3p.mag2() + SQ(mpi)) );
   Hep3Vector pi3m(Slct.trk_Pim->px(),
                   Slct.trk_Pim->py(),
                   Slct.trk_Pim->pz());
   HepLorentzVector LVpim( pi3m, sqrt(pi3m.mag2() + SQ(mpi)) );
   HepLorentzVector LVjpsi = Slct.LVcms - LVpip - LVpim;

   int Ng = Slct.gtrk.size();

   // analyse gg pairs: reject pi0->2gamma
   int Npi0 = 0;
   for ( int i = 0; i < Ng-1; ++i ) {
      const auto& LVgi = Slct.Pg[i];
      for ( int j = i+1; j < Ng; ++j ) {
         const auto& LVgj = Slct.Pg[j];

         double Mgg2 = (LVgi+LVgj).m2();
         hst[203]->Fill(Mgg2);
         hst[204]->Fill(Mgg2);

         if ( isMC ) {
            if ( Slct.decJpsi == 68 ) {
               hst[301]->Fill(Mgg2);
               if ( Slct.dec_eta == 1 ) {
                  hst[303]->Fill(Mgg2);
               }
            } else {
               hst[302]->Fill(Mgg2);
            }
         }

         if ( 0.013 < Mgg2 && Mgg2 < 0.022 ) { // CUT pi0
            Npi0 += 1;
         }
      } // end of for(j)
   } // end of for(i)
   hst[205]->Fill(Npi0);
   if ( isMC ) {
      if ( Slct.decJpsi == 68 ) {
         hst[304]->Fill(Npi0);
         if ( Slct.dec_eta == 1 ) {
            hst[306]->Fill(Npi0);
         }
      } else {
         hst[305]->Fill(Npi0);
      }
      if ( Npi0 > 0 ) {
         hst[307]->Fill(Slct.decJpsi);
      }
   }

   // reject events with Npi0 > 0:
   if ( Npi0 > 0 ) {
      return;
   }

   // cut for K+K- invariant mass:
   Hep3Vector p3p(Slct.trk_Kp[0]->px(),
                  Slct.trk_Kp[0]->py(),
                  Slct.trk_Kp[0]->pz());
   HepLorentzVector LVkp( p3p, sqrt(p3p.mag2() + SQ(mk)) );
   Hep3Vector p3m(Slct.trk_Km[0]->px(),
                  Slct.trk_Km[0]->py(),
                  Slct.trk_Km[0]->pz());
   HepLorentzVector LVkm( p3m, sqrt(p3m.mag2() + SQ(mk)) );
   HepLorentzVector LVkk = LVkp+LVkm;
   double Mkk2 = LVkk.m2();
   hst[207]->Fill(Mkk2);
   if ( isMC ) {
      if ( Slct.decJpsi == 68 ) {
         hst[308]->Fill(Mkk2);
      } else {
         hst[309]->Fill(Mkk2);
      }
   }
   if ( Mkk2 <= 1.01 || Mkk2 >= 1.07 ) {  // CUT Mphi^2 = 1.039 +/- 0.006
      return;
   }

   // search for << J/Psi -> phi+eta  >> such that the
   // recoil mass of full system is the smallest
   double M2fr_min = 10;
   int i_min = -1;
   vector<HepLorentzVector> LVg_min(3); // save gammas here

   for ( int ig = 0; ig < Ng; ++ig ) {
      const auto& LVgi = Slct.Pg[ig];

      // calculate recoil moment
      HepLorentzVector LVrec = LVjpsi - LVkk - LVgi;
      // it must be gamma (from eta decay)
      HepLorentzVector LVgj;
      LVgj.setVectM(LVrec.vect(), 0.); // set zero mass

      // check that these gammas satisfy the decay of eta
      HepLorentzVector LVgg = LVgi + LVgj;
      double Mgg2 = LVgg.m2();
      hst[211]->Fill(Mgg2);
      if ( isMC ) {
         if ( Slct.decJpsi == 68 ) {
            hst[311]->Fill(Mgg2);
            if ( Slct.dec_eta == 1 ) {
               hst[313]->Fill(Mgg2);
            }
         } else {
            hst[312]->Fill(Mgg2);
         }
      }

      // Mgg^2 ~ 0.295 +/- 0.021
      if ( fabs(Mgg2-0.295) < 0.045 ) { // Select eta CUT
         // total recoil mass
         HepLorentzVector LVfr = LVrec - LVgj;
         double M2fr = LVfr.m2();
         if ( fabs(M2fr) < fabs(M2fr_min) ) {
            M2fr_min = M2fr;
            i_min = ig;
            LVg_min[0] = LVgi;
            LVg_min[1] = LVgj;
         }
      }
   } // end of gammas loop

   hst[215]->Fill(M2fr_min);
   if ( isMC ) {
      if ( Slct.decJpsi == 68 ) {
         hst[315]->Fill(M2fr_min);
         if ( Slct.dec_eta == 1 ) {
            hst[317]->Fill(M2fr_min);
         }
      } else {
         hst[316]->Fill(M2fr_min);
      }
   }

   if ( fabs(M2fr_min) > 0.01 ) { // CUT minimal total recoil mass
      return;
   }

   // search for real gamma close to predicted
   // with minimal recoil mass of full system
   const auto& LVgp = LVg_min[1]; // predicted
   HepLorentzVector LVrec = LVjpsi - LVkk - LVg_min[0];

   double M2real_min = 10;
   int j_min = -1;
   double rE_min = 0, dTh_min = 100.;
   for ( int jg = 0; jg < Ng; ++jg ) {
      if ( jg == i_min ) {
         continue;
      }
      const auto& LVgj = Slct.Pg[jg];

      double rE = LVgp.e()/LVgj.e();
      double cosTh = LVgp.vect().cosTheta( LVgj.vect() );
      double dTh = RtoD( acos(cosTh) );
      hst[221]->Fill(rE);
      hst[222]->Fill(dTh);
      if ( isMC ) {
         if ( Slct.decJpsi == 68 ) {
            hst[321]->Fill(rE);
            hst[322]->Fill(dTh);
         } else {
            hst[323]->Fill(rE);
            hst[324]->Fill(dTh);
         }
      }

      if ( 0.4 < rE && rE < 1.8 && dTh < 10 ) { // CUT
         HepLorentzVector LVfr = LVrec - LVgj;
         double M2fr = LVfr.m2();
         if ( fabs(M2fr) < fabs(M2real_min) ) {
            M2real_min = M2fr;
            j_min = jg;
            rE_min = rE;
            dTh_min = dTh;
         }
      } //
   } // end of for(jg)
   int found = ( j_min != -1 ) ? 1 : 0;

   double Mgg_found = 0.;
   if ( found ) {
      // save found real gamma in LVg_min
      const auto& LVgj = Slct.Pg[j_min];
      LVg_min[2] = LVgj;

      hst[223]->Fill(dTh_min,rE_min);
      hst[225]->Fill(M2real_min);
      if ( isMC ) {
         if ( Slct.decJpsi == 68 ) {
            hst[325]->Fill(dTh_min,rE_min);
         } else {
            hst[326]->Fill(dTh_min,rE_min);
         }
         if ( Slct.decJpsi == 68 ) {
            hst[331]->Fill(M2real_min);
         } else {
            hst[332]->Fill(M2real_min);
         }
      }

      // check eta-mass
      HepLorentzVector LVgg = LVg_min[0] + LVgj;
      double Mgg2 = LVgg.m2();
      hst[226]->Fill(Mgg2);
      if ( isMC ) {
         if ( Slct.decJpsi == 68 ) {
            hst[335]->Fill(Mgg2);
         } else {
            hst[336]->Fill(Mgg2);
         }
      }

      // here we must apply the same selection criteria  for "ete"
      // as in the main analysis

      // Vertex fit
      WTrackParameter wp[4] = {
         WTrackParameter(mpi,Slct.trk_Pip->getZHelix(),
                         Slct.trk_Pip->getZError() ),
         WTrackParameter(mpi,Slct.trk_Pim->getZHelix(),
                         Slct.trk_Pim->getZError() ),
         WTrackParameter(mk, Slct.trk_Kp[0]->getZHelix(),
                         Slct.trk_Kp[0]->getZError() ),
         WTrackParameter(mk, Slct.trk_Km[0]->getZHelix(),
                         Slct.trk_Km[0]->getZError() )
      };

      HepPoint3D vx(0., 0., 0.);
      HepSymMatrix Evx(3, 0);
      double bx = 1E+6;
      double by = 1E+6;
      double bz = 1E+6;
      Evx[0][0] = bx*bx;
      Evx[1][1] = by*by;
      Evx[2][2] = bz*bz;

      VertexParameter vxpar;
      vxpar.setVx(vx);
      vxpar.setEvx(Evx);

      VertexFit* vtxfit = VertexFit::instance();
      vtxfit->init();
      for(int k = 0; k < 4; k++) {
         vtxfit->AddTrack(k,wp[k]);
      }
      vtxfit->AddVertex(0, vxpar,0, 1);
      bool ok = vtxfit->Fit(0);

      if ( ok ) {
         hst[231]->Fill(vtxfit->chisq(0));

         // Kinematik fit
         vtxfit->BuildVirtualParticle(0);
         vtxfit->Swim(0); // translate parameters of tracks to the vertex# 0
         for(int k = 0; k < 4; k++) {
            wp[k] = vtxfit->wtrk(k);
         }
         KinematicFit* kmfit = KinematicFit::instance();
         kmfit->init();
         for ( int k = 0; k < 4; k++ ) {
            kmfit->AddTrack(k, wp[k]);
         }
         kmfit->AddTrack(4, 0.0, Slct.gtrk[i_min]);
         kmfit->AddTrack(5, 0.0, Slct.gtrk[j_min]);
         kmfit->AddFourMomentum(0, Slct.LVcms);
         bool oksq = kmfit->Fit();
         if ( oksq ) {
            hst[232]->Fill(kmfit->chisq());
            // Momentums after kinematic corrections:
            HepLorentzVector Pg1 = kmfit->pfit(4);
            HepLorentzVector Pg2 = kmfit->pfit(5);
            HepLorentzVector Pgg = Pg1 + Pg2;
            Mgg_found  = Pgg.m();

            hst[233]->Fill(Mgg_found);
            if ( isMC ) {
               if ( Slct.decJpsi == 68 ) {
                  hst[338]->Fill(Mgg_found);
               } else {
                  hst[339]->Fill(Mgg_found);
               }
            }

            // set found flag to 2 if we found good eta
            if ( fabs(Mgg_found-meta) < weta ) { // see cuts.h
               found = 2;
            }
         } // end kinematic fit
      } // end vertex fit
   } // end if( found )

   // save momentum and angle of eta && photons in ntuple
   HepLorentzVector LVeta = LVg_min[0] + LVg_min[1];
   double Mgg_pred = LVeta.m();
   Double_t xfill[] = {
         LVkk.vect().mag(), LVkk.cosTheta(),
         LVeta.vect().mag(), LVeta.cosTheta(),
         LVg_min[0].e(), LVg_min[0].cosTheta(),
         LVg_min[1].e(), LVg_min[1].cosTheta(),
         LVg_min[2].e(), LVg_min[2].cosTheta(),
         double(found),
         rE_min,dTh_min,
         M2fr_min, M2real_min,
         Mgg_pred,Mgg_found,
         double(Slct.decJpsi),
         double(Slct.dec_eta)
   };
   m_tuple[2]->Fill( xfill );

   if ( isMC ) {
      JpsiTblE.vecdec = JpsiTbl.vecdec; // propagate for Eff table
      JpsiTblE.Add(); // save in table
//       if ( Slct.decJpsi != 68 ) { // DEBUG
//          cout << " CHECK-Eff: decJpsi= " << Slct.decJpsi
//               << " JpsiTbl= :" << JpsiTbl.StrDec() << ":"
//               << endl;
//          selector->PrintMcDecayTree(-99,0);
//       }
   }

}

//-------------------------------------------------------------------------
bool PsipJpsiPhiEtaEvent( ReadDst*       selector,
                          TEvtHeader*    m_TEvtHeader,
                          TDstEvent*     m_TDstEvent,
                          TEvtRecObject* m_TEvtRecObject,
                          TMcEvent*      m_TMcEvent,
                          TTrigEvent*    m_TTrigEvent,
                          TDigiEvent*    m_TDigiEvent,
                          THltEvent*     m_THltEvent      ) {
//-------------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " start " << __func__ << "()" << endl;
   }

   m_abscor->AbsorptionCorrection(selector);

   Select Slct; // information for the current event

   //----------------------------------------------------------------------
   //-- Get event information --
   //----------------------------------------------------------------------
   int runNo   = m_TEvtHeader->getRunId();
   int eventNo = m_TEvtHeader->getEventId();
   Slct.runNo  = runNo;
   Slct.event  = eventNo;
   hst[1]->Fill(0); // "cuts"

   // define data taking period:
   if ( DataPeriod == 0 ) {
      // expect one data taking period for all runs in job
      int run = abs(runNo);
      if ( (run >= 8093 && run <= 9025)  ) {         // 2009 Psi(2S)
         DataPeriod = 2009;
      } else if ( (run >= 25338 && run <= 27090) ) { // 2012 Psi(2S)
         DataPeriod = 2012;
      } else {
         DataPeriod = -1;
         cout << " WARNING: Data taking period undefined for runNo= "
              << runNo << endl;
         Warning("Undefined data taking period");
      }
   }

   double Ecms = 3.686; // (GeV) energy in center of mass system
   Slct.LVcms = HepLorentzVector(Ecms*sin(beam_angle), 0, 0, Ecms);

   isMC = (runNo < 0);
   if ( (isMC && m_TMcEvent->getMcParticleCol()->IsEmpty() ) ||
         ( !isMC &&
           m_TMcEvent != 0 &&
           !m_TMcEvent->getMcParticleCol()->IsEmpty() )
      ) {
      cout << " WARNING: something wrong: isMC= " << isMC
           << " McParticles= " << m_TMcEvent->getMcParticleCol()->GetEntries()
           << endl;
      Warning("Incorrect number of MC particles");
      return false;
   }

   FillHistoMC(selector,Slct); // MC histo

   Hep3Vector xorigin = getVertexOrigin(runNo);
   Slct.xorig = xorigin;

   //----------------------------------------------------------------------
   if ( !ChargedTracksPiPi(selector, Slct) ) {
      return false;
   }
   if ( !ChargedTracksKK(selector, Slct) ) {
      return false;
   }

   //----------------------------------------------------------------------
   int ng = NeutralTracks(selector, Slct); // number of selected gammas

   //----------------------------------------------------------------------
   // normal selection
   bool fret = false;
   if ( ng >= 2 ) {
      fret = VertKinFit(selector, Slct);
   }

   //----------------------------------------------------------------------
   // study eta efficiency
   if ( ng >= 1 ) {
      EtaEff(selector, Slct);
   }

   //----------------------------------------------------------------------
   return fret;
}

//-------------------------------------------------------------------------
void PsipJpsiPhiEtaEndJob(ReadDst* selector) {
//-------------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << __func__ << "()" << endl;
   }

   if ( isMC ) {
      // print tables of decays
      cout << string(65,'#') << endl;
      cout << "Decays of Psi(2S) 1" << endl
           << "       Mrec in [3.092, 3.102] "
           << PsipTbl.ntot << " decays" << endl
           << "       size of table is " << PsipTbl.Size() << endl;
      PsipTbl.Print(0.1); // do not print decays with P<0.1% of all
      cout << "Enddecay" << endl << endl;

      cout << string(65,'#') << endl;
      cout << "Decays of Psi(2S) 2" << endl
           << "       Mrec in [3.055, 3.145] "
           << PsipTbl2.ntot << " decays" << endl
           << "       size of table is " << PsipTbl2.Size() << endl;
      PsipTbl2.Print(0.1);
      cout << "Enddecay" << endl << endl;

      // FIXME: how to save cuts automatically?
      string my_cuts(
         " List of cuts for J/Psi -> phi eta:\n "
         "      chisq<80 && 2*mk<Mkk<1.08 && abs(Mgg-meta)<3*seta"
      );
      cout << my_cuts << endl << endl;

      cout << string(65,'#') << endl;
      cout << "Decays of J/psi 1" << endl
           << "       Mrec in [3.092, 3.102] "
           << JpsiTbl.ntot << " decays" << endl
           << "       size of table is " << JpsiTbl.Size() << endl;
      JpsiTbl.Print(0.1);
      cout << "Enddecay" << endl << endl;

      cout << string(65,'#') << endl;
      cout << "Decays of J/psi 2" << endl
           << "       Mrec in [3.055, 3.145] "
           << JpsiTbl2.ntot << " decays" << endl
           << "       size of table is " << JpsiTbl2.Size() << endl;
      JpsiTbl2.Print(0.1);
      cout << "Enddecay" << endl << endl;

      cout << string(65,'#') << endl;
      cout << "Decays of J/psi Eff Phi Eta" << endl
           << "   reconstructed events " << JpsiTblE.ntot << endl
           << "       size of table is " << JpsiTblE.Size() << endl;
      JpsiTblE.Print(0.01); // do not print decays with P<0.01% of all
      cout << "Enddecay" << endl << endl;

      cout << string(65,'#') << endl;
      cout << "Decays of Psi(2S) 4C" << endl
           << "       Mrec in [3.00, 3.05] & [3.144, 3.194] "
           << PsipTbl4C.ntot << " decays" << endl
           << "       size of table is " << PsipTbl4C.Size() << endl;
      PsipTbl4C.Print(0.01); // do not print decays with P<0.01% of all
      cout << "Enddecay" << endl << endl;
   }

   string module = string(__func__);
   module = module.substr(0,module.size()-6);
   if ( warning_msg.empty() ) {
      cout << " There are no warnings in " << module << endl;
   } else {
      cout << " Check output for WARNINGS in " << module << endl;
      for(auto it = warning_msg.begin(); it != warning_msg.end(); ++it) {
         cout << it->first << " : " << it->second << endl;
      }
   }
}

//-------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
