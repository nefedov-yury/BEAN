//======================================================================//
//                                                                      //
// PsipMuMuGamma:                                                       //
// Study of the photon detection efficiency in process:                 //
//              e+ e- -> gamma(ISR) + mu+ mu-                           //
// using the data on energy Psi(2S)                                     //
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
// #include <TNtupleD.h>
// #include <TF1.h>
#include <TDatabasePDG.h>

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

#include "DstEvtRecTracks.h"
#include "ReadDst.h"

using namespace std;

//-----------------------------------------------------------------------------
// Struct to save variables per one event
//-----------------------------------------------------------------------------
struct MuMuGamma {
   // Run-info
   int runNo;              // run-number
   HepLorentzVector LVcms; // Momentum in center of mass system
   Hep3Vector xorig;       // interaction point from DB

   // MC information:
   int decPsip;            // decay code for MC events

   // tracks
   vector<int> Trk; // indexes of good tracks
   vector<double> ImpR; // impact parameters in XY plane
   vector<double> ImpZ; // Z-plane
   vector<RecMdcKalTrack*> trk_Mu;  // muon tracks
   vector<HepLorentzVector> LVmu;   // muon momentums

   HepLorentzVector LVisr; // predicted ISR photon (1C-fit);
   double chi2mu;
   double chisq;

   vector<RecEmcShower*>   trk_g;   // gamma tracks
   vector<HepLorentzVector> LVg;    // photons momentum

   MuMuGamma() : decPsip{0}, chi2mu(0),chisq(0) {}
};

//-----------------------------------------------------------------------------
// Global variables
//-----------------------------------------------------------------------------
static const double beam_angle = 0.011; // 11 mrad

// masses of particles (GeV)           from PDG:
static const double mpsip  = 3.686097; // 3686.097  +/- 0.025   MeV
static const double mjpsi  = 3.096916; // 3096.916  +/- 0.011   MeV
static const double mchic0 = 3.41471;  // 3414.71   +/- 0.30    MeV
static const double mchic1 = 3.51067;  // 3510.67   +/- 0.05    MeV
static const double mchic2 = 3.55617;  // 3556.17   +/- 0.07    MeV
static const double mmu    = 0.105658; // 105.6583745 +/- 0.0000024 MeV
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
// static vector<TNtupleD*> m_tuple;

// container for warnings
static map<string,int> warning_msg;

// expect one data taking period for all events
static int DataPeriod = 0;

static bool isMC = false;

//-----------------------------------------------------------------------------
// Functions: use C-linkage names
//-----------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
inline void Warning(const string& msg) {
   warning_msg[msg] += 1;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
inline double RtoD(double ang) {
   return ang*180/M_PI;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
inline double SQ(double x) {
   return x*x;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void PsipMuMuGammaStartJob(ReadDst* selector) {
//-----------------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
   }

   hst.resize(300,nullptr);
//    m_tuple.resize(10,static_cast<TNtupleD*>(0));

   // init Absorption Correction
   m_abscor = new AbsCor(selector->AbsPath("Analysis/AbsCor"));

   //initialize EventTag
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
      Warning("EventTagSvc has already been initialized");
   }

   // We have to initialize DatabaseSvc ---------------------------------------
   DatabaseSvc* dbs = DatabaseSvc::instance();
   if ( (dbs->GetDBFilePath()).empty() ) {
      // set path to directory with databases:
      dbs->SetDBFilePath(selector->AbsPath("Analysis/DatabaseSvc/dat"));
   }

   // We have to initialize Magnetic field ------------------------------------
   MagneticFieldSvc* mf = MagneticFieldSvc::instance();
   if ( (mf->GetPath()).empty() ) {
      // set path to directory with magnetic fields tables
      mf->SetPath(selector->AbsPath("Analysis/MagneticField"));
      mf->UseDBFlag(false); // like in the boss program
      mf->RunMode(3); // like in the boss program
   } else {
      cout << " WARNING:"
           << " MagneticFieldSvc has already been initialized" << endl
           << "          path = " << mf->GetPath() << endl;
      Warning("MagneticFieldSvc has already been initialized");
   }

   //--------- Book histograms ------------------------------------------------

   // DATA
   hst[1] = new TH1D("C1_Rxy", "1) R_{xy}", 400,-2.,2.);
   hst[2] = new TH1D("C1_Rz",  "1) R_{z}", 200,-20.,20.);
   hst[3] = new TH1D("C1_CT",  "1) cos(#theta)", 200,-1.,1.);
   hst[4] = new TH1D("C1_Ntrk","1) N_{trk}", 10,-0.5,9.5);
   hst[5] = new TH1D("C1_Ztrk","1) Z_{trk}", 9,-4.5,4.5);
   hst[6] = new TH1D("C1_Eemc","1) E(emc)", 100,0,1.);
   hst[7] = new TH1D("C1_Ep",  "1) E(emc)/P(mdc)", 110,0,1.1);
   hst[8] = new TH2D("C1_2D",  "1) 2D;P(mdc);E(emc)", 100,0,2, 100,0,2);
   hst[9] = new TH1D("C1_Nmu", "1) N_{mu} after Emc", 3,-0.5,2.5);

   hst[11] = new TH1D("C1_depth","1) MUC depth (cm)", 80,-15.,65.);
   hst[12] = new TH1D("C1_Nmuc", "1) N_{mu} after MuC", 3,-0.5,2.5);

   hst[15] = new TH1D("C1_tofN","1) Ntof per track", 10,-0.5,9.5);
   hst[16] = new TH1D("C1_tof", "1) TOF", 100, 1.,5.);
   hst[17] = new TH1D("C1_Ntof","1) N(tof)", 3,-0.5,2.5);
   hst[18] = new TH1D("C1_dTof","1) tof1-tof2", 200,-10,10);

   hst[21] = new TH1D("C1_Pmum","1) P(#mu^{-})", 300,0.,3.0);
   hst[22] = new TH1D("C1_Pmup","1) P(#mu^{+})", 300,0.,3.0);
   hst[23] = new TH1D("C1_Cmum","1) cos #theta(#mu^{-})", 200,-1.,1.);
   hst[24] = new TH1D("C1_Cmup","1) cos #theta(#mu^{+})", 200,-1.,1.);

   hst[31] = new TH1D("C2a_ok", "2a) vertex fit: bad(0)/good(1)", 2,-0.5,1.5);
   hst[32] = new TH1D("C2a_ch2","2a) vertex fit #chi^2", 100,0.,100.);
   hst[33] = new TH1D("C2b_2mu","2b) 4C fit 2#mu #chi^2", 100,0,100);
   hst[34] = new TH1D("C2c_ok", "2c) 1C fit: bad(0)/good(1)", 2,-0.5,1.5);
   hst[35] = new TH1D("C2c_ch2","2c) 1C fit #chi^2", 100,0,100);
   hst[38] = new TH1D("C2d_CdT","2d) cos #Theta(#gamma^{p} #mu)", 200,-1.,1.);
   hst[39] = new TH1D("C2d_dT", "2d) |#Theta_{#gamma^{p} #mu}|", 180,0.,90.);

   hst[41] = new TH1D("C2_inv","2) Minv2(mu+ mu-)", 500,5.,15.);
   hst[42] = new TH1D("C2_rec","2) before fit Mrec2(mu+ mu-)", 400,-0.2,0.2);
   hst[43] = new TH1D("C2_pEg","2) E(#gamma) predicted", 300,0,3.0);
   hst[44] = new TH1D("C2_pCg","2) cos #theta(#gamma) predicted",200,-1,1);

   // search for "hot" chanels in EMC
   hst[51] = new TH2D("C3_mapB","3) barrel Z vs phi",
                      270,-135.,+135., 360,-180.,180.);
   hst[52] = new TH2D("C3_mapC1","3) endcap Z<0 rho vs phi",
                      40, 50.,90., 360,-180.,180.);
   hst[53] = new TH2D("C3_mapC2","3) endcap Z>0 rho vs phi",
                      40, 50.,90., 360,-180.,180.);

   hst[55] = new TH1D("C3_Ng",   "3) N(#gamma)", 10,-0.5,9.5);
   hst[56] = new TH1D("C3a_Mg2", "3a) M^2(#gamma#gamma)", 200,0.,0.04);
   hst[57] = new TH1D("C3a_Mg",  "3a) M(#gamma#gamma)", 200,0.,1.);
   hst[58] = new TH1D("C3a_Npi0","3a) N(#pi^{0})", 10,-0.5,9.5);

   hst[61] = new TH1D("C3b_M2mu","3b) Minv2(#mu^{+} #mu^{-})", 150,8.8,10.3);
   hst[62] = new TH1D("C3b_inv", "3b) Minv2(J/Psi gamma)", 500,9.,14.);
   hst[63] = new TH1D("C3b_Nchi","3b) N(#chi_{CJ})", 10,-0.5,9.5);
   hst[64] = new TH2D("C3b_NN",  "3b) N(#chi_{CJ});N(#gamma)",
                      10,-0.5,9.5, 10,-0.5,9.5);

   hst[71] = new TH1D("C4_pEg",  "4) E(#gamma) predicted", 200,0,2.0);
   hst[72] = new TH1D("C4_pCg",  "4) cos #theta(#gamma) predicted",200,-1,1);
   hst[73] = new TH1D("C4_pPhig","4) #phi(#gamma) predicted",200,-M_PI,M_PI);
   hst[74] = new TH2D("C4_pCPhi","4) predicted #gamma;cos #theta;#phi",
                      100,-1,1, 100,-M_PI,M_PI);

   hst[81] = new TH1D("C4_rE","4) E(#gamma)/E(#gamma^{p})", 200,0.,2.);
   hst[82] = new TH1D("C4_dT","4) #Theta_{#gamma, #gamma^{p}}", 200,0.,100.);
   hst[83] = new TH2D("C4_2D","4) E/Ep vs dTheta;#Theta", 60,0,30,80,0.6,1.4);

   // study matching predicted-real
   hst[91] = new TH1D("C4_dEall","4) E(#gamma)-E(#gamma^{p})", 100,-0.2,0.2);
   hst[92] = new TH1D("C4_dE1","4) E(#gamma)-E(#gamma^{p}) 1b", 100,-0.1,0.1);
   hst[93] = new TH1D("C4_rE1","4) E(#gamma)/E(#gamma^{p}) 1b", 100,0.,2.);
   hst[94] = new TH1D("C4_dE2","4) E(#gamma)-E(#gamma^{p}) 2-5b", 200,-0.1,0.1);
   hst[95] = new TH1D("C4_rE2","4) E(#gamma)/E(#gamma^{p}) 2-5b", 200,0.,2.);
   hst[96] = new TH1D("C4_dE3","4) E(#gamma)-E(#gamma^{p}) >5b", 200,-0.1,0.1);
   hst[97] = new TH1D("C4_rE3","4) E(#gamma)/E(#gamma^{p}) >5b", 200,0.,2.);
   hst[98] = new TH1D("C4_dThe","4) #Theta_{#gamma, #gamma^{p}}", 60,0.,30.);

   // results:
   hst[201] = new TH2D("D_g1","E vs cos(Theta) ok",
                      18,-0.9,0.9, 18,0.05,1.85);
   hst[202] = new TH2D("D_g0","E vs cos(Theta) no",
                      18,-0.9,0.9, 18,0.05,1.85);

   hst[205] = new TH2D("D_imp","Impact parameters;Rxy;Rz",100,-1,1,100,-10,10);
   hst[206] = new TH2D("D_g1S","E vs cos(Theta) ok",
                      18,-0.9,0.9, 18,0.05,1.85);
   hst[207] = new TH2D("D_g0S","E vs cos(Theta) no",
                      18,-0.9,0.9, 18,0.05,1.85);
   hst[208] = new TH2D("D_g1B","E vs cos(Theta) ok",
                      18,-0.9,0.9, 18,0.05,1.85);
   hst[209] = new TH2D("D_g0B","E vs cos(Theta) no",
                      18,-0.9,0.9, 18,0.05,1.85);

   // Monte Carlo histograms:
   hst[101] = new TH1D("mc0_dec", "dec Psi(2S) all",256,-0.5,255.5);
   hst[102] = new TH1D("mc2_dec", "dec Psi(2S) all",256,-0.5,255.5);
   hst[103] = new TH1D("mc3a_dec","dec Psi(2S) all",256,-0.5,255.5);
   hst[104] = new TH1D("mc3b_dec","dec Psi(2S) all",256,-0.5,255.5);

   hst[108] = new TH1D("mc5_dec1","dec Psi(2S) all",256,-0.5,255.5);
   hst[109] = new TH1D("mc5_dec0","dec Psi(2S) all",256,-0.5,255.5);

   hst[111] = new TH1D("mc0_pdg", "PDG codes of all particles",
                       2001,-1000.5,1000.5);
   hst[112] = new TH1D("mc0_pdg0","PDG of particles from primary vertex",
                       2001,-1000.5,1000.5);
   // for e+ e- -> (ISR gamma) mu+ mu-
   hst[114] = new TH1D("mc0_Pmum", "P(#mu^{-})", 100,0.,2.);
   hst[115] = new TH1D("mc0_Cmum", "cos #Theta(#mu^{-})", 100,-1.,1.);
   hst[116] = new TH1D("mc0_Pmup", "P(#mu^{+})", 100,0.,2.);
   hst[117] = new TH1D("mc0_Cmup", "cos #Theta(#mu^{+})", 100,-1.,1.);
   hst[118] = new TH1D("mc0_Eisr", "E(#gamma ISR)", 100,0.,2.);
   hst[119] = new TH1D("mc0_Cisr", "cos #Theta(#gamma ISR)", 100,-1.,1.);

   // chi_CJ
   hst[121] = new TH1D("mc3b_Ng", "N(#gamma) for chi_CJ", 10,-0.5,9.5);
   hst[122] = new TH1D("mc3b_M2mu","Minv2(#mu^{+} #mu^{-})", 150,8.8,10.3);
   hst[123] = new TH1D("mc3b_Minv2", "Minv2(J/Psi gamma)", 500,9.,14.);

   // dec vs Eg(predicted)
   hst[131] = new TH1D("mc4_decP", "pik 0.48 < E < 0.56",256,-0.5,255.5);
   hst[132] = new TH1D("mc4_dec",  "not in pik",256,-0.5,255.5);


   // register in selector to save in given directory
   const char* SaveDir = "MuMuGamma";
   VecObj hsto(hst.begin(),hst.end());
   selector->RegInDir(hsto,SaveDir);
//    VecObj ntuples(m_tuple.begin(),m_tuple.end());
//    selector->RegInDir(ntuples,SaveDir);

}

//-----------------------------------------------------------------------------
static Hep3Vector getVertexOrigin(int runNo, bool verbose = false) {
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
static void FillHistoMC(const ReadDst* selector, MuMuGamma& mmg) {
//-----------------------------------------------------------------------------
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
      static int nprt = 0;
      if ( nprt < 1 ) {
         nprt += 1;
         printf(" WARNING: MC is not Psi(2S) evTag= 0x%08x\n",evTag);
      }
      Warning("MC is not Psi(2s)");
      evTag = 0;
   }

   // decay code of Psi(2S):
   int decPsip  = (evTag >> 8) & 0xFF;
   mmg.decPsip = decPsip;
   hst[101]->Fill(mmg.decPsip);

   double E_isr = 0.; // sum of all ISR photons
   TIter mcIter(mcParticles);
   while( auto part = static_cast<TMcParticle*>(mcIter.Next()) ) {
      int p_pdg = part->getParticleID();
//       int p_idx = part->getTrackIndex();
      int p_mother = part->getMother();

      Hep3Vector Vp( part->getInitialMomentumX(),
                     part->getInitialMomentumY(),
                     part->getInitialMomentumZ() );

      hst[111]->Fill(p_pdg);
      if ( p_mother == -99 ) { // primary vertex
         hst[112]->Fill(p_pdg);

         if ( abs(p_pdg) == 22 ) {      // ISR gamma
            E_isr += part->getInitialMomentumE();
            hst[119]->Fill(Vp.cosTheta());
         } else if ( p_pdg == 13 ) {    // mu-
            hst[114]->Fill(Vp.mag());
            hst[115]->Fill(Vp.cosTheta());
         } else if ( p_pdg == -13 ) {   // mu+
            hst[116]->Fill(Vp.mag());
            hst[117]->Fill(Vp.cosTheta());
         }
      }
      hst[118]->Fill(E_isr);

   } // end of while
}

//-----------------------------------------------------------------------------
static bool ChargedTracks(ReadDst* selector, MuMuGamma& mmg) {
//-----------------------------------------------------------------------------
   static const double Rvxy0_max = 1.0;
   static const double Rvz0_max = 10.0;
//    static const double cosTheta_max = 0.80;  // barrel only
   static const double cosTheta_max = 0.93;

   const TEvtRecObject* m_TEvtRecObject = selector->GetEvtRecObject();
   const TEvtRecEvent* evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

   //--------------------------------------------------------------------------
   // 1) select two opposite charged muons
   //--------------------------------------------------------------------------
   mmg.Trk.reserve(min(1,evtRecEvent->totalCharged()));
   int Nz = 0; // sum of charges
   for ( int i = 0; i < evtRecEvent->totalCharged(); i++ ) {
      auto itTrk = static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(i));
      if( !itTrk->isMdcTrackValid() ) {
         continue;
      }

      RecMdcTrack* mdcTrk = itTrk->mdcTrack();

      double theta = mdcTrk->theta();
      double cosTheta = cos(theta);

      HepVector a = mdcTrk->helix();
      HepSymMatrix Ea = mdcTrk->err();
      HepPoint3D point0(0.,0.,0.);   // initial point for MDC recosntruction
      HepPoint3D IP(mmg.xorig[0],mmg.xorig[1],mmg.xorig[2]);
      VFHelix helixip(point0,a,Ea);
      helixip.pivot(IP);
      HepVector vecipa = helixip.a();
      double Rvxy0 = vecipa[0]; // the nearest distance to IP in xy plane
      double Rvz0  = vecipa[3]; // ... in z direction

      hst[1]->Fill(Rvxy0);
      hst[2]->Fill(Rvz0);
      hst[3]->Fill(cosTheta);

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
      mdcKalTrk->setPidType(RecMdcKalTrack::muon);
      mmg.trk_Mu.push_back(mdcKalTrk);
      mmg.ImpR.push_back(Rvxy0);
      mmg.ImpZ.push_back(Rvz0);
      Nz += mdcKalTrk->charge();

      mmg.Trk.push_back(i);

   } // end for charged tracks loop

   int Ntrk = mmg.Trk.size();
   hst[4]->Fill(Ntrk);
   hst[5]->Fill(Nz);
   if ( Ntrk != 2 ) {
      return false;
   }
   if ( Nz != 0 ) {
      return false;
   }

   // eliminate electrons by Emc
   int Nmu = 0;
   for ( int k = 0; k < 2; k++ ) {
      auto itTrk = static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(mmg.Trk[k]));
      if ( !itTrk->isEmcShowerValid() ) {
         continue;
      }
      RecEmcShower* emcTrk = itTrk->emcShower();
      double Eemc = emcTrk->energy();
      double Pmdc = mmg.trk_Mu[k]->p();
      double E_p =  Eemc / Pmdc;
      hst[6]->Fill(Eemc);
      hst[7]->Fill(E_p);
      hst[8]->Fill(Pmdc,Eemc);

//       if ( E_p < 0.8 ) { // 1-st variant
//          Nmu += 1;
//       }
      if ( Eemc < 0.3 ) { // 2-nd var
         Nmu += 1;
      }
   }
   hst[9]->Fill(Nmu);
   if ( Nmu != 2 ) {
      return false;
   }

   // select muons by Muon Chambers signal
   int Nmuc = 0;
   for ( int k = 0; k < 2; k++ ) {
      auto itTrk = static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(mmg.Trk[k]));
      if ( !itTrk->isMucTrackValid() ) {
         continue;
      }
      RecMucTrack* mucTrk = itTrk->mucTrack();
      double mucdepth = mucTrk->depth();
      hst[11]->Fill(mucdepth);
      if ( mucdepth > 35. ) {
         Nmuc += 1;
      }
   }
   hst[12]->Fill(Nmuc);
   if ( Nmuc != 2 ) {
      return false;
   }

   // remove cosmic muons by signal in TOF
   vector<double> Ttof;
   for ( int k = 0; k < 2; k++ ) {
      auto itTrk = static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(mmg.Trk[k]));
      if ( !itTrk->isTofTrackValid() ) {
         continue;
      }
      const vector<RecTofTrack* >& tof_trk = itTrk->tofTrack();
      hst[15]->Fill( tof_trk.size() );
      Ttof.push_back( tof_trk[0]->tof() );
      for ( const auto& tt : tof_trk ) {
         hst[16]->Fill( tt->tof() );
      }
   }
   int Ntof = Ttof.size();
   hst[17]->Fill( Ntof );
   double dT = 100.;
   if ( Ntof == 2 ) {
      dT = Ttof[0] - Ttof[1];
      hst[18]->Fill( dT );
   }
   if ( fabs(dT) > 2. ) {
      return false;
   }

   // fill LVmu and histograms
   mmg.LVmu.reserve(2);
   for ( auto trk : mmg.trk_Mu ) {
      Hep3Vector Vmu( trk->px(), trk->py(), trk->pz());
      mmg.LVmu.push_back(HepLorentzVector( Vmu, sqrt(Vmu.mag2()+SQ(mmu)) ));
      int hZ = ( trk->charge() > 0 );
      hst[21+hZ]->Fill(trk->p());
      hst[23+hZ]->Fill(Vmu.cosTheta());
   }

   //--------------------------------------------------------------------------
   // 2) Vertex & Kinematic Fit
   //--------------------------------------------------------------------------
   // 2a) search for a good vertex

   HepPoint3D vx(0., 0., 0.);
   HepSymMatrix Evx(3, 0);
   Evx[0][0] = SQ(1e+6);
   Evx[1][1] = SQ(1e+6);
   Evx[2][2] = SQ(1e+6);

   VertexParameter vxpar;
   vxpar.setVx(vx);
   vxpar.setEvx(Evx);

   VertexFit* vtxfit = VertexFit::instance();
   vtxfit->init();
   for ( int k = 0; k < 2; k++ ) {
      const auto& trk = mmg.trk_Mu[k];
      vtxfit->AddTrack( k,
                        WTrackParameter(mmu,trk->getZHelix(),trk->getZError())
                      );
   }
   vtxfit->AddVertex(0, vxpar,0, 1);
   bool ok = vtxfit->Fit(0);

   hst[31]->Fill( double(ok) );
   if ( !ok ) {
      return false;
   }
   hst[32]->Fill(vtxfit->chisq(0));

   vtxfit->BuildVirtualParticle(0);
   vtxfit->Swim(0); // translate parameters of tracks to the vertex# 0

   // 2b) 4C kinematic fit of two muons only
   KinematicFit* kmfit = KinematicFit::instance();
   kmfit->init();
   for ( int k = 0; k < 2; k++ ) {
      kmfit->AddTrack( k, vtxfit->wtrk(k) );
   }
   kmfit->AddFourMomentum(0, mmg.LVcms);
   bool ok2mu = kmfit->Fit();
   double chi2mu = ( (ok2mu) ? kmfit->chisq() : 200. );
   hst[33]->Fill(chi2mu);

   if ( chi2mu < 30. ) { // we do not need ISR photon TODO:
      return false;
   }

   // 2c) 4C kinematic fit with fake ISR photon
   // recoil momentum of muons
   HepLorentzVector LVrec = mmg.LVcms - mmg.LVmu[0] - mmg.LVmu[1];

   // make fake ISR photon:
   //   the direction along recoil momentum and zero mass
   //   and set big errors for its parameters
   HepLorentzVector LVisr( LVrec.vect(), LVrec.vect().mag());
   double dphi = 1.;  // 1rad
   double dthe = 1.;  // 1rad
   double dE   = 3.;  // 3GeV
   WTrackParameter wISR( mmg.xorig, LVisr, dphi,dthe,dE );

   kmfit->init();
   for ( int k = 0; k < 2; k++ ) {
      kmfit->AddTrack( k, vtxfit->wtrk(k) );
   }
   kmfit->AddTrack(2, wISR);
   kmfit->AddFourMomentum(0, mmg.LVcms);
   bool oksq = kmfit->Fit();

   hst[34]->Fill( double(oksq) );
   if ( !oksq ) {
      return false;
   }
   double chisq = kmfit->chisq();
   hst[35]->Fill(chisq);

   if ( chisq > 25. ) {
      return false;
   }

   // results of fit
   for ( int k = 0; k < 2; k++ ) {
      mmg.LVmu.push_back( kmfit->pfit(k) );
   }
   mmg.LVisr = kmfit->pfit(2);
   mmg.chi2mu = chi2mu;
   mmg.chisq = chisq;

   // 2d) ignore events if the predicted gamma collinear to muons
   int flag_ang = 0;
   for ( const auto& LVmu : mmg.LVmu ) {
      double cosGmu = Hep3Vector(mmg.LVisr).cosTheta(Hep3Vector(LVmu));
      hst[38]->Fill(cosGmu);
      double ang = RtoD(acos(fabs(cosGmu)));
      hst[39]->Fill(ang);
      if ( ang < 5 ) { // degree!
         flag_ang += 1;
      }
   }

   if ( flag_ang > 0 ) {
      return false;
   }

   // invariant mass and recoil mass
   double inv2mu = (mmg.LVmu[0]+mmg.LVmu[1]).m2();
   hst[41]->Fill(inv2mu);
   hst[42]->Fill(LVrec.m2());

   double pEg = mmg.LVisr.e();
   double pCg = mmg.LVisr.cosTheta();
   hst[43]->Fill(pEg);
   hst[44]->Fill(pCg);

   // ignore events with a predicted gamma in an insensitive area
   if ( pEg < 0.02 ||
         fabs(pCg) > 0.95
      ) {
      return false;
   }

   if ( isMC ) {
      hst[102]->Fill(mmg.decPsip);
   }

   return true;
}

//-----------------------------------------------------------------------------
static bool NeutralTracks(ReadDst* selector, MuMuGamma& mmg) {
//-----------------------------------------------------------------------------
   // parameters of reconstruction
   static const double min_angle = 10 * M_PI/180; // 10 grad

   const TEvtRecObject* m_TEvtRecObject = selector->GetEvtRecObject();
   const TEvtRecEvent* evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

   //--------------------------------------------------------------------------
   // 3) neutral tracks
   //--------------------------------------------------------------------------
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
      // map of gammas
      if ( absCosTheta < 0.8 ) {  //barrel
         hst[51]->Fill(emcpos.z(), RtoD(emcpos.phi()) );
      } else { // endcap
         int hz = (emcpos.z() > 0);
         hst[52+hz]->Fill(emcpos.rho(), RtoD(emcpos.phi()) );
      }

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
      if ( tang < min_angle ) {
         continue;
      }

      mmg.trk_g.push_back( emcTrk );
      Hep3Vector p3 = emcpos - mmg.xorig;
      p3 *= eraw / p3.mag();
      mmg.LVg.push_back( HepLorentzVector(p3,eraw) );

   } // end of neutrals loop

   int Ng = mmg.trk_g.size(); // number of good photons
   hst[55]->Fill(Ng);

   // 3a) search for pi0
   int Npi0 = 0;
   for ( int i = 0; i < Ng-1; ++i ) {
      auto& LVgi = mmg.LVg[i];
      for ( int j = i+1; j < Ng; ++j ) {
         auto& LVgj = mmg.LVg[j];

         double Mgg2 = (LVgi+LVgj).m2(); // invariant mass of (gg)
         hst[56]->Fill(Mgg2);
         double Mgg = (Mgg2>0) ? sqrt(Mgg2) : 0.;
         hst[57]->Fill(Mgg);

         // mpi0^2=0.0182; sigma=0.002
         if ( fabs(Mgg2-SQ(mpi0)) < 0.006 ) { // strong cut!
            Npi0 += 1;
         }
      }
   }
   hst[58]->Fill(Npi0);

   if ( Npi0 > 0 ) {
      return false;
   }

   if ( isMC ) {
      hst[103]->Fill(mmg.decPsip);
   }

   // 3b) search for chi_Cj -> gamma J/Psi
   int NchiCj = 0;
   bool isMC_chiCJ = ( mmg.decPsip >=8 && mmg.decPsip <= 10 );
   HepLorentzVector LV2mu = mmg.LVmu[0] + mmg.LVmu[1];
   double inv2mu = LV2mu.m2();
   hst[61]->Fill(inv2mu);
   // sigma = 0.1; SQ(mjpsi) = 9.59
   if ( fabs(inv2mu-SQ(mjpsi)) < 0.5 ) { // J/Psi -> mu+mu-
      for ( const auto& LVg : mmg.LVg ) {
         double Minv2 = (LV2mu+LVg).m2();
         hst[62]->Fill(Minv2);

         if ( fabs(Minv2-12.45) < 0.55 ||
               fabs(Minv2-10.6) < 0.5 ) {  // TODO ??
            NchiCj += 1;
         }

         if ( isMC_chiCJ ) {
            hst[123]->Fill(Minv2);
         }
      }
   }
   hst[63]->Fill(NchiCj);
   hst[64]->Fill(Ng,NchiCj);

   if ( NchiCj > 0 ) {
      return false;
   }

   if ( isMC ) {
      hst[104]->Fill(mmg.decPsip);

      if ( mmg.decPsip >=8 && mmg.decPsip <= 10 ) {
         hst[121]->Fill(Ng);
         hst[122]->Fill(inv2mu);
      }
   }

   return true;
}

//-----------------------------------------------------------------------------
static void GammaEff(ReadDst* selector,MuMuGamma& mmg) {
//-----------------------------------------------------------------------------
//    const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

   //--------------------------------------------------------------------------
   // 4) the photon detection efficiency
   //--------------------------------------------------------------------------

   // predicted photon
   double pEg = mmg.LVisr.e();
   double pCg = mmg.LVisr.cosTheta();
   double pPhig = mmg.LVisr.phi();
   Hep3Vector pVg = mmg.LVisr.vect();

   hst[71]->Fill(pEg);
   hst[72]->Fill(pCg);
   hst[73]->Fill(pPhig);
   hst[74]->Fill(pCg,pPhig);
   if ( isMC ) {
      if ( pEg > 0.48 && pEg < 0.56 ) {
         hst[131]->Fill(mmg.decPsip);

//          if ( (mmg.decPsip >= 8 && mmg.decPsip <= 10)
//               ||  mmg.decPsip == 67
//             ) { // DEBUG
//             cout << " CHECK: mmg.decPsip= " << mmg.decPsip
//                  << " pEg= " << pEg << endl;
//             selector->PrintMcDecayTree(-99,0);
//          }
      } else {
         hst[132]->Fill(mmg.decPsip);
      }
   }

   // search for real photon
   bool find = false;
   for ( const auto& lvg : mmg.LVg ) {
      // ratio of energies
      double ratE = lvg.e()/pEg;
      hst[81]->Fill(ratE);

      // Theta(gamma gamma_p)
      double cosTT = lvg.vect().cosTheta(pVg);
      double dTheta = RtoD( acos(cosTT) );
      hst[82]->Fill(dTheta);
      hst[83]->Fill( dTheta, ratE );

      // study matching predicted-real
      if ( dTheta < 15 ) {
         double dE = lvg.e() - pEg;
         hst[91]->Fill(dE);
         if ( 0.05 < pEg && pEg < 0.15 ) {
            hst[92]->Fill(dE);
            hst[93]->Fill(ratE);
         } else if ( pEg < 0.55 ) {
            hst[94]->Fill(dE);
            hst[95]->Fill(ratE);
         } else {
            hst[96]->Fill(dE);
            hst[97]->Fill(ratE);
         }
      }
      if ( fabs(ratE-1) < 0.3 ) {
         hst[98]->Fill(dTheta);
      }

      // matching
      find = (fabs(ratE-1) < 0.3) && (dTheta < 15.);
   }

   if ( find ) {
      hst[201]->Fill(pCg,pEg);
      if ( isMC ) {
         hst[108]->Fill(mmg.decPsip);
      }
   } else {
      hst[202]->Fill(pCg,pEg);
      if ( isMC ) {
         hst[109]->Fill(mmg.decPsip);
      }
   }

   // study cosmic background by impact parameter
   bool fsmall_IP = true;
   for ( int k = 0; k < 2; ++k ) {
      hst[205]->Fill( mmg.ImpR[k], mmg.ImpZ[k] );
      if ( fabs(mmg.ImpR[k]) > 0.1 || fabs(mmg.ImpZ[k]) > 2 ) {
         fsmall_IP = false;
      }
   }

   if ( fsmall_IP ) {
      if ( find ) {
         hst[206]->Fill(pCg,pEg);
      } else {
         hst[207]->Fill(pCg,pEg);
      }
   }

   bool fbig_IP =
      ( fabs(mmg.ImpR[0]) > 0.5 || fabs(mmg.ImpZ[0]) > 5 ) &&
      ( fabs(mmg.ImpR[1]) > 0.5 || fabs(mmg.ImpZ[1]) > 5 );
   if ( fbig_IP ) {
      if ( find ) {
         hst[208]->Fill(pCg,pEg);
      } else {
         hst[209]->Fill(pCg,pEg);
      }
   }

}

// Loop for each event
//-----------------------------------------------------------------------------
bool PsipMuMuGammaEvent( ReadDst*       selector,
                         TEvtHeader*    m_TEvtHeader,
                         TDstEvent*     m_TDstEvent,
                         TEvtRecObject* m_TEvtRecObject,
                         TMcEvent*      m_TMcEvent,
                         TTrigEvent*    m_TTrigEvent,
                         TDigiEvent*    m_TDigiEvent,
                         THltEvent*     m_THltEvent       ) {
//-----------------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " start " << __func__ << "()" << endl;
   }

   m_abscor->AbsorptionCorrection(selector);

   MuMuGamma mmg;

   //--------------------------------------------------------------------------
   //-- Get event information --
   //--------------------------------------------------------------------------
   int runNo   = m_TEvtHeader->getRunId();
   mmg.runNo  = runNo;
//    int eventNo = m_TEvtHeader->getEventId();

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
   mmg.LVcms = HepLorentzVector(Ecms*sin(beam_angle), 0, 0, Ecms);

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

   FillHistoMC(selector,mmg); // MC histo

   Hep3Vector xorigin = getVertexOrigin(runNo);
   mmg.xorig = xorigin;

   // Study of the photon reconstruction efficiency
   if ( !ChargedTracks(selector, mmg) ) {
      return false;
   }
   if ( !NeutralTracks(selector, mmg) ) {
      return false;
   }
   GammaEff(selector,mmg);

   return true;
}

//-----------------------------------------------------------------------------
void PsipMuMuGammaEndJob(ReadDst* selector) {
//-----------------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << __func__ << "()" << endl;
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

//-----------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
