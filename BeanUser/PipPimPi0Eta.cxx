//======================================================================//
//                                                                      //
// SelectOmegaEta - search for e+ e- -> pi+ pi- pi0 Eta                 //
//                                                   |-> 2 gammas       //
//                                                                      //
//======================================================================//

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <vector>
// #include <map>

#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TVector3.h>

#include "CLHEP/Units/PhysicalConstants.h"
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/TwoVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Geometry/Point3D.h>

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
   typedef HepGeom::Point3D<double> HepPoint3D;
#endif

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
using CLHEP::pi;
using CLHEP::twopi;


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
#include "TofHitStatus.h"
#include "ReadDst.h"

using namespace std;

struct Select {
  double Ecms;                  // energy in center of mass system

  // charge tracks:
  vector<int> qtrks;            // trk# in EvtRecTrkCol
  int ipip, ipim;               // trk# for charge +/-
  HepLorentzVector ppip, ppim;  // 4-momentum of trks
  int ipip_mc, ipim_mc;         // correspondent mcParticle
  int pip_pdg, pim_pdg;         // PDG code of these mc-particles (or 0)

  HepLorentzVector P_missing;   // missing 4-momentum
  double SoverS0;               // for check MCGPJ

  // neutral (gammas) tracks
  vector<int> gammas;           // trk# in EvtRecTrkCol
  vector<HepLorentzVector> Pg;  // 4-momentum of gammas
  int ig1, ig2, ig3, ig4;       // the best gammas after 4C kin.fit
  HepLorentzVector ppi01, ppi02;// 4-momentum of pairs (g+g) ???

  HepLorentzVector Ptot;        // 4-momentum of the overall system

  // momentums after kinematic fit
  HepLorentzVector Ppip, Ppim, Pgg1,Pgg2, Pomega, Pgg;

  Select() {
    Ecms = -1;
    ipip = ipim = -1;
    ipip_mc = ipim_mc = -1;
    pip_pdg = pim_pdg = 0;
    ig1 = ig2 = ig3 = ig4 = -1;
    SoverS0 = -1;
  }

  int nGood() { return qtrks.size(); }
  int mc_trk(int i) { return (qtrks[i] == ipip) ? ipip_mc : ipim_mc; }
  int isPiPlus(int i) { return qtrks[i] == ipip; }
  int isPiMinus(int i) { return qtrks[i] == ipim; }
  int nGam() { return gammas.size(); }
  int gamma1() { return (ig1 >= 0) ? gammas[ig1] : -1; }
  int gamma2() { return (ig2 >= 0) ? gammas[ig2] : -1; }
  int gamma3() { return (ig3 >= 0) ? gammas[ig3] : -1; }
  int gamma4() { return (ig4 >= 0) ? gammas[ig4] : -1; }
  bool GoodGammas()
        {return (ig1 != -1 && ig2 != -1 && ig3 != -1 && ig4 != -1); }
};

//-----------------------------------------------------------------------------
// Global file variables
//-----------------------------------------------------------------------------
const static bool init_selection = true;

const static double beam_angle = 0.011; // 11 mrad

// masses of particles (GeV)           from PDG:
const static double mpi    = 0.13957;  // 139.57018 +/- 0.00035 MeV
const static double mpi0   = 0.13498;  // 134.9766  +/- 0.0006  MeV
const static double meta   = 0.547862; // 547.862   +/- 0.017   MeV
const static double momega = 0.78265;  // 782.65    +/- 0.12    MeV

static AbsCor* m_abscor = 0;
static EventTagSvc* m_EventTagSvc = 0;

static std::vector<TH1D*> his1;
static std::vector<TH2D*> his2;
static std::vector<TNtuple* > m_tuple;

static int ERROR_WARNING = 0;
static bool isMC = false;

static bool omega_to_pis = false;
static bool omega_best_mass = false;
static bool one_ISR = false;
static bool one_ISR_less50 = false;
static bool one_ISR_more50 = false;
static bool two_ISR = false;
static bool two_ISR_less50 = false;
static bool two_ISR_more50 = false;
static bool SoverS0 = false;
static bool EtaTo2Gammas = false;

// Functions: use C-linkage names
#ifdef __cplusplus
extern "C" {
#endif


inline Double_t RtoD(Double_t ang) {return ang*180/M_PI;}


void PipPimPi0EtaStartJob (ReadDst* selector)
{
  if ( selector->Verbose() ) cout << " Start: " << __func__ << "()" << endl;

  if ( init_selection ) {
    cout << " ===== INIT SELECTION =====" << endl;
  }

  his1.resize(500,(TH1D*)0);
  his2.resize(500,(TH2D*)0);
  m_tuple.resize(10,0);

  // init Absorption Correction
  m_abscor = new AbsCor(selector->AbsPath("Analysis/AbsCor"));

  //initialize EventTag
  m_EventTagSvc = EventTagSvc::instance();

  // set paths to pdg & decayCodes files:
  m_EventTagSvc->setPdtFile(
        selector->AbsPath("Analysis/EventTag/share/pdt.table") );
  m_EventTagSvc->setDecayTabsFile(
        selector->AbsPath("Analysis/EventTag/share/decay.codes") );

  if ( selector->Verbose() ) m_EventTagSvc->setVerbose(1);

  m_EventTagSvc->initialize();

  // We have to initialize DatabaseSvc ----------------------------------------
  DatabaseSvc* dbs = DatabaseSvc::instance();
  if ( (dbs->GetDBFilePath()).empty() ) {
    // set path to directory with databases:
    dbs->SetDBFilePath(selector->AbsPath("Analysis/DatabaseSvc/dat"));
  }

  // We have to initialize Magnetic field--------------------------------------
  MagneticFieldSvc* mf = MagneticFieldSvc::instance();
  if ( !(mf->GetPath()).empty() ) {
    cerr << " RhopiStartJob::WARNING:"
         << " MagneticFieldSvc has already been initialized" << endl
         << "                         path = " << mf->GetPath() << endl;
  }
  // set path to directory with magnetic fields tables
  mf->SetPath(selector->AbsPath("Analysis/MagneticField"));
  mf->UseDBFlag(false); // like in the boss program
  mf->RunMode(3); // like in the boss program


  his1[0] = new TH1D("InvMBef","Invariant mass before kinematic fit", 500, 0, 5);
  his1[1] = new TH1D("InvMAft","Invariant mass after kinematic fit", 500, 0, 5);
  his1[2] = new TH1D("InvMAft2","Invariant mass after kinematic fit2", 500, 0, 5);
  his1[3] = new TH1D("Pi0","Invariant mass for #pi_{0}", 500, 0, 5);
  his1[4] = new TH1D("Eta","Invariant mass for #eta", 500, 0, 5);
  his1[5] = new TH1D("Phi","Invariant mass for #phi", 500, 0, 5);

  m_tuple[0] = new TNtuple("a4C","after 4C fit",
                                "ch2:Eisr:PZisr:PTisr:"
                                "Etot:Ptot:Esum:PZsum:PTsum:"
                                "mcNg:mcEisr:mcEmisr:mcThisr:"
                                "dc:PZmis:PTmis:Mmis");

  VecObj his1o(his1.begin(),his1.end());
  selector->RegInDir(his1o, "one");

}

static Hep3Vector getVertexOrigin(int runNo)
{
  static int save_runNo = 0;
  static Hep3Vector xorigin;

  if ( runNo == save_runNo ) return xorigin;

  // update vertex for new run
  xorigin.set(0.,0.,0.);
  VertexDbSvc* vtxsvc = VertexDbSvc::instance();

  int run = abs(runNo);
  if (
          (run >= 28241 && run <= 28266)  // 3080-old
       || (run >=  9947 && run <= 10878)  // J/Psi 2009
     ) {
    vtxsvc->SetBossVer("6.6.3");
  } else if ( run < 30000 ) {   // here: J/Psi scan
    vtxsvc->SetBossVer("6.6.4");
  } else {
    vtxsvc->SetBossVer("6.6.5");  // here: 3080 from R-scan (2015)
  }

  vtxsvc->handle(runNo);

  if ( vtxsvc->isVertexValid() ) {
    double* dbv = vtxsvc->PrimaryVertex();
    xorigin.set(dbv[0],dbv[1],dbv[2]);
    // cout << " sqlite-db vertex: (x,y,z)= " << xorigin << endl;
  } else {
    cout << "Cannot obtain vertex information for run#" << runNo << endl;
    ERROR_WARNING++;
    exit(1);
  }

  save_runNo = runNo;
  return xorigin;
}

//-----------------------------------------------------------------------------
static double GetEcms(int run, double& Lumi)
//-----------------------------------------------------------------------------
{
  struct runInfo {
    int runS, runE; // first and last runs in period
    double Ebeam;   // energy of beam (MeV)
    double Spread;  // beam spread (KeV)
    double Lumi;    // luminosity (pb^-1)
  };

  static runInfo ListRuns[] =  {
        { 28312, 28346,     3050.254,    915,     14.918 },
        { 28347, 28381,     3059.273,    842,     15.059 },
        { 28241, 28266,     3080.224,    839,     17.513 },
        { 28382, 28387,     3082.491,    883,      2.711 },
        { 28466, 28469,     3084.351,   1019,      2.057 },
        { 28388, 28416,     3089.213,    905,     13.934 },
        { 28472, 28475,     3091.284,    777,      1.624 },
        { 28417, 28449,     3092.104,    864,     12.354 },
        { 28476, 28478,     3093.716,    837,      1.843 },
        { 28479, 28482,     3095.299,    819,      2.143 },
        { 28487, 28489,     3096.017,   1077,      1.816 },
        { 28490, 28492,     3096.432,    894,      2.134 },
        { 28493, 28495,     3097.783,    710,      2.069 },
        { 28496, 28498,     3098.943,    833,      2.203 },
        { 28499, 28501,     3099.572,   1170,      0.756 },
        { 28504, 28505,     3101.914,    960,      1.612 },
        { 28506, 28509,     3106.146,    981,      2.106 },
        { 28510, 28511,     3112.609,    667,      1.720 },
        { 28512, 28513,     3120.294,   1116,      1.264 },
        { 33490, 33556,     3810.951 + 0.55,   1088,     50.54 },
        { 33572, 33657,     3899.608 + 0.55,   1025,     52.61 },
        { 33659, 33719,     4090.238 + 0.55,    722,     52.63 },
        { 30372, 30437,     4190.557 + 0.55,    644,     43.09 },
        { 32046, 32140,     4218.961 + 0.55,   1497,     54.13 },
        { 30438, 30491,     4229.893 + 0.55,   1025,     44.40 },
        { 32141, 32226,     4245.102 + 0.55,   1583,     55.59 },
        { 30492, 30557,     4310.388 + 0.55,   1071,     44.90 },
        { 31281, 31325,     4389.833 + 0.55,   1577,     55.18 },
                        };

  int Np = sizeof(ListRuns)/sizeof(ListRuns[0]);

  int absrun = abs(run);

  // we need some values any way
  double Ecms   = 3.097;
  double Spread = 1000.;
  Lumi          = 0;

  int i = 0;
  for(; i < Np; i++) {
     if ( absrun >= ListRuns[i].runS && absrun <= ListRuns[i].runE ) {
       Ecms   = ListRuns[i].Ebeam  * 1.e-3;  // MeV -> GeV
       Spread = ListRuns[i].Spread * 1.e-6;  // KeV -> GeV
       Lumi   = ListRuns[i].Lumi;
       break;
     }
  }

  if ( i != Np ) {
    // correct energy according Yadi Wang fit of J/Psi peak"
    Ecms -= 0.55e-3; // -0.55MeV
  } else {
    // run not in the list
    if ( absrun >= 39355 && absrun <= 39618 ) { // 3080 from R-scan (2015)
      Ecms =  3.080; // no BEMS information
      Lumi = 120.11483;
    } else  if ( absrun >= 9947 && absrun <= 10878 ) { // J/Psi 2009
      Ecms =  3.097; // for inclusive MC
      Lumi = 0.;
    } else  if ( absrun >= 31327 && absrun <= 31390 ) {
      Ecms =  4.420; // no BEMS information
      Lumi =  44.67;
    } else if ( absrun >= 31983 && absrun <= 32045 ) {
      Ecms =  4.210; // no BEMS information
      Lumi =  54.55;
    } else if ( absrun >= 27147 && absrun <= 27233 ) {
      Ecms =  3.0827;
      Lumi = 13.5158;
      cout << " GetEcms::WARNING run# " << run
           << " this is old not recommended run" << endl;
    } else {
      cout << " GetEcms::WARNING unknown run# " << run
           << " use Ecms= " << Ecms << endl;
      ERROR_WARNING++;
    }
  }

  return Ecms;
}

bool PipPimPi0EtaEvent (ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
{
  if ( selector->Verbose() ) cout << " start " << __func__ << "()" << endl;

  m_abscor->AbsorptionCorrection(selector);

  Select Slct;

  int runNo   = m_TEvtHeader->getRunId();
  int eventNo = m_TEvtHeader->getEventId();

  Hep3Vector xorigin = getVertexOrigin(runNo);

  const TEvtRecEvent*  evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
  const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

  int Ncharge = 0;
  for(int i = 0; i < evtRecEvent->totalCharged(); i++) {

    DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
    if( !itTrk->isMdcTrackValid() ) continue;

    RecMdcTrack *mdcTrk = itTrk->mdcTrack();

    double theta = mdcTrk->theta();
    double cosTheta = cos(theta);

    HepVector a = mdcTrk->helix();
    HepSymMatrix Ea = mdcTrk->err();
    HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
    HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
    VFHelix helixip(point0,a,Ea);
    helixip.pivot(IP);
    HepVector vecipa = helixip.a();
    double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
    double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction

//    his1[11]->Fill(Rvxy0);
//    his1[12]->Fill(Rvz0);
//    his1[13]->Fill(fabs(cosTheta));

    if(fabs(Rvz0) >= 10.0) continue;
    if(fabs(Rvxy0) >= 1.0) continue;
    if(fabs(cosTheta) >= 0.93) continue;

    Slct.qtrks.push_back(i);
    Ncharge += mdcTrk->charge();
//    his1[14]->Fill( mdcTrk->charge()*mdcTrk->p() );
  }

  int nGood = Slct.nGood();

  if ( nGood != 2 || Ncharge != 0 ) return false;

  ParticleID *pid = ParticleID::instance();
#if (BOSS_VER < 700)
  pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
#else   
  pid->set_path(selector->AbsPath("Analysis/ParticleID"));
#endif

  int npp = 0, npm = 0; // number of +/- in event
  for(int i = 0; i < nGood; i++) {
    DstEvtRecTracks* itTrk=(DstEvtRecTracks*) evtRecTrkCol->At(Slct.qtrks[i]);

    pid->init();
    pid->setMethod(pid->methodProbability());
    pid->setChiMinCut(4);
    pid->setRecTrack(itTrk);

    // list of systems used for PID
    pid->usePidSys( pid->useDedx() |
                    pid->useTof1() | pid->useTof2() |
                    pid->useTofE() );
    //                 pid->useEmc()  |
    //                 pid->useMuc()

    pid->identify(
                pid->onlyPion()    |
                pid->onlyKaon()    |
                pid->onlyProton()  |
                pid->onlyMuon()    |
                pid->onlyElectron()
        );

    pid->calculate(runNo);

    if ( !pid->IsPidInfoValid() ) return false;

    if ( !init_selection ) {
      // select Pions only
      if ( pid->probPion() <  0.001 ) return false;
    }

    RecMdcKalTrack* mdcKalTrk = itTrk->mdcKalTrack();
//    his1[25]->Fill( (double)(mdcKalTrk != 0) );
    if ( !mdcKalTrk ) {
      cout << " No Kalman track ! EventNo= "<< eventNo << endl;
      return false;
    }

    // define charge by Kalman fit
    mdcKalTrk->setPidType(RecMdcKalTrack::pion);

    //-- E,Px,Py,Pz for good charged tracks
    Hep3Vector p3pi(mdcKalTrk->px(), mdcKalTrk->py(), mdcKalTrk->pz());
    double epi = sqrt( p3pi.mag2() + mpi*mpi );
    if ( mdcKalTrk->charge() > 0 ) {
      npp++;
      Slct.ipip = Slct.qtrks[i];
      Slct.ppip = HepLorentzVector( p3pi, epi );
    } else {
      npm++;
      Slct.ipim = Slct.qtrks[i];
      Slct.ppim = HepLorentzVector( p3pi, epi );
    }
  }

  if ( !init_selection ) {
    // require exactly one "+" and one "-"
    if ( npp != 1 || npm != 1 ) return false;
  }

  for(int i = evtRecEvent->totalCharged();
          i < evtRecEvent->totalTracks(); i++) {
    DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
    if(!itTrk->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = itTrk->emcShower();

    // 1) good EMC time:
    if ( emcTrk->time() < 0 || emcTrk->time() > 14 ) continue;

    // 2) good EMC energy deposited in the barrel (endcap) part of EMC
    double eraw = emcTrk->energy();
    double absCosTheta = fabs(  cos(emcTrk->theta()) );

    bool GoodCluster=false;
    if ( absCosTheta < 0.8 ) {  //barrel
      GoodCluster = eraw > 25E-3;
    } else if ( absCosTheta > 0.85 && absCosTheta < 0.92 ) { //endcap
      GoodCluster = eraw > 50E-3;
    }
    if ( !GoodCluster ) continue;

    // 3) the nearest charged track is far from cluster
    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

    double tang = 200.; // min angle between cluster and track
    for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
      DstEvtRecTracks* jtTrk = (DstEvtRecTracks*) evtRecTrkCol->At(j);
      if ( !jtTrk->isExtTrackValid() ) continue;
      RecExtTrack *extTrk = jtTrk->extTrack();
      if ( extTrk->emcVolumeNumber() == -1 )
        continue; //track does not hit EMC
      Hep3Vector extpos = extTrk->emcPosition();
      double angd = fabs(extpos.angle(emcpos)); // [0,pi]
      if ( angd < tang ) tang = angd;
    }

    static const double min_angle = 10 * M_PI/180; // 10 grad
    if ( tang < min_angle ) continue;

    Slct.gammas.push_back(i);

    // (E,Px,Py,Pz) for good gammas
    Hep3Vector p3 = emcpos - xorigin; // position emc cluster wrt of vertex
    p3 *= eraw / p3.mag();            // normalization on energy
    Slct.Pg.push_back( HepLorentzVector(p3,eraw) ); 
  }

  int nGam = Slct.nGam();
  
  if ( nGam != 4 ) return false;

  his1[0]->Fill( ( Slct.ppip + Slct.ppim + Slct.Pg[0] + Slct.Pg[1] + Slct.Pg[2] + Slct.Pg[3] ).m() );
  
  // Vertex & Kinematic Fit
  RecMdcKalTrack *pipTrk = ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.ipip))->mdcKalTrack();
  RecMdcKalTrack *pimTrk = ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.ipim))->mdcKalTrack();

  WTrackParameter wvpipTrk, wvpimTrk;
  wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());
  wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());

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
  vtxfit->AddTrack(0,  wvpipTrk);
  vtxfit->AddTrack(1,  wvpimTrk);
  vtxfit->AddVertex(0, vxpar,0, 1);
  bool ok = vtxfit->Fit(0);
  //his1[14]->Fill( double(ok) );
  if ( !ok )  return false; //skip event

  vtxfit->BuildVirtualParticle(0);
  //his1[15]->Fill(vtxfit->chisq(0));
  vtxfit->Swim(0); // translate parameters of tracks to the vertex# 0

  WTrackParameter wpip = vtxfit->wtrk(0);
  WTrackParameter wpim = vtxfit->wtrk(1);

  KinematicFit * kmfit = KinematicFit::instance();
  //--------------------------------------------------------------------------------------------------------------------------------Ecms =  GetEcms(runNo, lum);
   double lum = 0;
   double Ecms =  GetEcms(abs(runNo), lum);

   // double Ecms =  3.097;  // for inclusive MC


  HepLorentzVector ecms(Ecms*sin(beam_angle), 0, 0, Ecms);
  //------------------------------------------------------4C Fit
  // make fake ISR photon: Pt component must be small non zero
  // because of the error (bug) in the WTrackParameter
  double tmp_xy=1.e-6;
  double tmp_z=0.01;
  double tmp_e = sqrt(2*tmp_xy*tmp_xy+tmp_z*tmp_z);
  HepLorentzVector pisr(tmp_xy,tmp_xy,tmp_z,tmp_e);
  double dphi = 1.e-3;  // 1mrad
  double dthe = 1.e-3;  // 1mrad
  double dE = 3.;        // 3 GeV
  WTrackParameter wISR(xorigin,pisr,dphi,dthe,dE);


  int ig[4] = {-1,-1,-1,-1};
  double chisq = 9999.;
  double chisq4C = 9999., chisq5C = 9999.;
  for(int i = 0; i < nGam-3; i++) {
    RecEmcShower *g1Trk =
         ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.gammas[i]))->emcShower();
    for(int j = i+1; j < nGam-2; j++) {
      RecEmcShower *g2Trk =
         ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.gammas[j]))->emcShower();
      for(int k = j+1; k < nGam-1; k++){
        RecEmcShower *g3Trk =
          ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.gammas[k]))->emcShower();
        for(int l = k+1; l < nGam; l++){
          RecEmcShower *g4Trk =
            ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.gammas[l]))->emcShower();

          kmfit->init();
          kmfit->AddTrack(0, wpip);
          kmfit->AddTrack(1, wpim);
          kmfit->AddTrack(2, 0.0, g1Trk);
          kmfit->AddTrack(3, 0.0, g2Trk);
          kmfit->AddTrack(4, 0.0, g3Trk);
          kmfit->AddTrack(5, 0.0, g4Trk);
#ifdef  add_ISR
          kmfit->AddTrack(6, wISR);
#endif
          kmfit->AddFourMomentum(0, ecms);
          bool oksq = kmfit->Fit();
          //           cout << " i,j,k,l="<< i << j << k << l << " oksq= " << oksq
          //                << " chi2= " << kmfit->chisq() << endl
          //                << "    P_isr= " << kmfit->pfit(6) << endl;
          if ( oksq ) {
            double chi2 = kmfit->chisq();
            // his1[12]->Fill(chi2);
            if ( chi2 < chisq ) {
              chisq = chi2;
              ig[0] = i;
              ig[1] = j;
              ig[2] = k;
              ig[3] = l;
            }
          }
        }//------------------------------------------End for(l)
      }//--------------------------------------------End for(k)
    }//-----------------------------------------------End for(j)
  } //-----------------------------------------------End for(i)

  //his1[31]->Fill(chisq);
  //  if ( !Slct.GoodGammas() ) return false; // can not find 4  good gammas
 if ( ig[0]==-1 || ig[1]==-1 ||  ig[2]==-1 || ig[3]==-1) return false; // can not find 4  good gammas
 /*//his1[107]->Fill(3);
 	if(isMC && omega_best_mass) his1[280]->Fill(3);
 if(isMC){
    if(one_ISR) his1[300]->Fill(3);
    if(two_ISR) his1[301]->Fill(3);
    if(one_ISR && omega_best_mass) his1[302]->Fill(3);
    if(two_ISR && omega_best_mass) his1[303]->Fill(3);
  }
 chisq4C = chisq;
 if (chisq > 100) return false;
 //his1[107]->Fill(4);
 if(isMC && omega_best_mass) his1[280]->Fill(4);
 if(isMC){
    if(one_ISR) his1[300]->Fill(4);
    if(two_ISR) his1[301]->Fill(4);
    if(one_ISR && omega_best_mass) his1[302]->Fill(4);
    if(two_ISR && omega_best_mass) his1[303]->Fill(4);
  }


 //P_missing for E_mis & M_mis
 Slct.P_missing = ecms - Slct.ppip - Slct.ppim - Slct.Pg[ig[0]] - Slct.Pg[ig[1]] - Slct.Pg[ig[2]] - Slct.Pg[ig[3]];

    HepLorentzVector P_isr = kmfit->pfit(6);

  float xfill[] = {
                kmfit->chisq(),P_isr.e(),P_isr.z(),P_isr.perp(),
                P_tot.e(),P_tot.rho(),P_sum.e(),P_sum.z(),P_sum.perp(),
                Slct.Ng,Slct.Eisr,Slct.Emax_isr,Slct.Th_isr,
                decCode,Slct.Pmiss.z(),Slct.Pmiss.perp(),Slct.Pmiss.m()
                  };
  m_tuple[0]->Fill( xfill );*/

  his1[1]->Fill( ( Slct.ppip + Slct.ppim + Slct.Pg[0] + Slct.Pg[1] + Slct.Pg[2] + Slct.Pg[3] ).m() );  
  his1[2]->Fill( ( Slct.ppip + Slct.ppim + Slct.Pg[ig[0]] + Slct.Pg[ig[1]] + Slct.Pg[ig[2]] + Slct.Pg[ig[3]] ).m() );
  for ( int j = 0; j < 4; j++ )
    for ( int i = j + 1; i < 4; i++ )
      if ( (Slct.Pg[j] + Slct.Pg[j+i]).m() > 0.135 - 0.010 && (Slct.Pg[j] + Slct.Pg[j+i]).m() < 0.135 + 0.010 )
        {
          his1[3]->Fill( ( Slct.Pg[j] + Slct.Pg[j+i] ).m() );
          his1[5]->Fill( ( Slct.ppip + Slct.ppim + Slct.Pg[j] + Slct.Pg[j+i] ).m() );
        }
  for ( int j = 0; j < 4; j++ )
    for ( int i = j + 1; i < 4; i++ )
      if ( (Slct.Pg[j] + Slct.Pg[j+i]).m() > 0.547 - 0.010 && (Slct.Pg[j] + Slct.Pg[j+i]).m() < 0.547 + 0.010 )
        his1[4]->Fill( ( Slct.Pg[j] + Slct.Pg[j+i] ).m() );


  return false;
  }

}

void PipPimPi0EtaEndJob (ReadDst* selector)
{
  if ( selector->Verbose() ) cout << __func__ << "()" << endl;

  if ( ERROR_WARNING != 0 ) {
    cout << " ERROR_WARNING= " << ERROR_WARNING << " Check output!" << endl;
  }

  if ( init_selection ) {
    cout << " ===== INIT SELECTION =====" << endl;
  }
  
  his1[1]->Fit("gausn");
  
}
