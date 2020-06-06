//-----------------------------------------------------------------------------
// AbsCor: <=> boss/Analysis/PhotonCor/AbsCor
//
//         Nefedov: The original program from the BOSS can not be used
//         directly in the BEAN:
//         1) Corrections 'edgecor' use functions from the
//            detector description package (see EmcRecGeoSvc)
//            This is definitely not for BEAN.
//         2) The number of modifications to be applied exceeds
//            reasonable limits. The program becomes unreadable.
//
//-----------------------------------------------------------------------------

#include "AbsCor/AbsCor.h"   // must be first!

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"

#include "DstEvtRecTracks.h"

#ifdef  AbsCor_00_00_36
#include "Identifier/EmcID.h"
#endif

using namespace std;

//-----------------------------------------------------------------------------
AbsCor::AbsCor( const string& path,
                bool _usetof,
                bool _dodatacor,
                bool _edgecor
#ifdef  AbsCor_00_00_36
                ,
                bool _hotcellmask,
                bool _dopi0Cor,
                bool _MCuseTof
#endif
              )
//-----------------------------------------------------------------------------
{
  usetof = _usetof;
  dodatacor = _dodatacor;

  edgecor = _edgecor;
  if( edgecor ) {
      cout << " This type of correction can not be done" << endl;
      exit(1);
  }

#ifdef  AbsCor_00_00_36
  hotcellmask = _hotcellmask;
  dopi0Cor = _dopi0Cor;
  MCuseTof = _MCuseTof;
#endif

  string DataPathc = path;
  if( usetof ) {
#ifdef  AbsCor_00_00_28
    DataPathc += "/dat/00-00-28/c3ptof.txt";
#elif defined( AbsCor_00_00_36 )
    DataPathc += "/dat/00-00-36/c3ptof.txt";
#endif
  } else {
#ifdef  AbsCor_00_00_28
    DataPathc += "/dat/00-00-28/c3p.txt";
#elif defined( AbsCor_00_00_36 )
    DataPathc += "/dat/00-00-36/c3p.txt";
#endif
  }
  ifstream inc3p(DataPathc.c_str(),ios::in);
  if ( !inc3p.is_open() ) {
     cout << " can not open: " << DataPathc << endl;
     exit(1);
  }

  for(int i = 0; i < 4; i++){
    double am,ame;
    inc3p >> am;
    inc3p >> ame;
    ai[i] = am;
    // cerr << " ai[" << i << "]= " << am << endl;
  }
  if( inc3p.fail() ) {
    cout << " error while reading file " << DataPathc << endl;
    exit(1);
  }

#ifdef  AbsCor_00_00_36
  // Read energy correction parameters from PhotonCor/McCor
  string paraPath = path;
  if ( MCuseTof ) {
    paraPath += "/dat/00-00-36/evsetTof.txt";
  } else {
    paraPath += "/dat/00-00-36/evset.txt";
  }
  ifstream in2(paraPath.c_str(),ios::in);
  if ( !in2.is_open() ) {
     cout << " can not open: " << paraPath << endl;
     exit(1);
  }
  double energy,thetaid,peak1,peakerr1,res,reserr;
  dt = new TGraph2DErrors();
  dtErr = new TGraph2DErrors();
  //for(int i=0;i<560;i++){
  for(int i=0;i<1484;i++){  //53*28
    in2>>energy;
    in2>>thetaid;
    in2>>peak1;
    in2>>peakerr1;
    in2>>res;
    in2>>reserr;
    dt->SetPoint(i,energy,thetaid,peak1);
    dt->SetPointError(i,0,0,peakerr1);
    dtErr->SetPoint(i,energy,thetaid,res);
    dtErr->SetPointError(i,0,0,reserr);
    if(i<28) e25min[int(thetaid)]=energy;
    if(i>=1484-28) e25max[int(thetaid)]=energy;
    // if(i>=560-28) e25max[int(thetaid)]=energy;
  }
  if( in2.fail() ) {
    cout << " error while reading file " << paraPath << endl;
    exit(1);
  }
  in2.close();
#endif

  //-----------------------------------------------------------------
  // Suppression of hot crystals
  // Reading the map from hotcry.txt (Hajime, Jul 2013)
  for(int ih=0;ih<10;ih++){
    hrunstart[ih]=-1;
    hrunend[ih]=-1;
    hcell[ih]=-1;
  }
  int numhots=4; // numbers of hot crystals
  int dumring,dumphi,dummod,dumid;
  string HotList = path + "/dat/hotcry.txt";

  ifstream hotcrys;
  hotcrys.exceptions( ifstream::failbit | ifstream::badbit );
  hotcrys.open(HotList.c_str(),ios::in);
  for(int il=0; il<numhots; il++){
    hotcrys>>hrunstart[il];
    hotcrys>>hrunend[il];
    hotcrys>>hcell[il];
    hotcrys>>dumring;
    hotcrys>>dumphi;
    hotcrys>>dummod;
    hotcrys>>dumid;
  }
  hotcrys.close();
  //-----------------------------------------------------------------

}


// Nefedov: I keep this function for back compatibility
//          see BeanUser/ DaubleDTag & RadBhabha
//-----------------------------------------------------------------------------
void AbsCor::SuppressHotCrystals(ReadDst* selector)
//-----------------------------------------------------------------------------
{
  if( selector->Verbose() ) cout << " SuppressHotCrystals() " << endl;

  int runNo = const_cast<TEvtHeader* >(selector->GetEvtHeader())->getRunId();

  const TEvtRecEvent* evtRecEvent =
                                selector->GetEvtRecObject()->getEvtRecEvent();
  const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

  if( evtRecEvent->totalTracks() > evtRecTrkCol->GetSize() ) return;
  if( evtRecEvent->totalTracks() > 50 ) return;

  for(int i = 0; i < evtRecEvent->totalTracks(); i++) {
    DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
    if( !itTrk->isEmcShowerValid() ) continue;

    RecEmcShower* emcTrk = itTrk->emcShower();

    // If it is "hot", return "9999" (Hajime, Jul 2013)

    for ( int ih = 0; ih < 10; ih++ ) {
      if ( (hrunstart[ih] == -1) ||
           (hrunend[ih] == -1) ||
           (hcell[ih] == -1) ) continue;
      if ( (abs(runNo) < hrunstart[ih]) ||
           (abs(runNo) > hrunend[ih]) ) continue;

      if ( emcTrk->cellId() == hcell[ih] ) {
        emcTrk->setStatus(9999);

        if (selector->Verbose()) {
          cout << " AbsCor::SuppressHotCrystals(): suppressed cell="
               << emcTrk->cellId()
               << " energy=" << emcTrk->energy() << endl;
        }
      }
    }
  }

}


//-----------------------------------------------------------------------------
void AbsCor::AbsorptionCorrection(ReadDst* selector)
//-----------------------------------------------------------------------------
{
  if( selector->Verbose() ) cout << " AbsorptionCorrection() " << endl;

  int runNo = const_cast<TEvtHeader* >(selector->GetEvtHeader())->getRunId();

  const TEvtRecEvent* evtRecEvent =
                                selector->GetEvtRecObject()->getEvtRecEvent();
  const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

  if( evtRecEvent->totalTracks() > evtRecTrkCol->GetSize() ) return;
  if( evtRecEvent->totalTracks() > 50 ) return;

  for(int i = 0; i < evtRecEvent->totalTracks(); i++) {
    DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
    if( !itTrk->isEmcShowerValid() ) continue;

    RecEmcShower* emcTrk = itTrk->emcShower();

#ifdef  AbsCor_00_00_28
    // this piece of code prevents that AbsCor corrections
    // from being applied twice
    int st = emcTrk->status();
    if( selector->Verbose() )
      cout << " AbsCor:: EMC status= " << emcTrk->status();
    if( st > 10000 ) {
      if( selector->Verbose() ) cout << " skip ... " << endl;
      continue;
    }
    emcTrk->setStatus(st+20000);
    if( selector->Verbose() ) cout << " set to " << emcTrk->status() << endl;
#endif

#ifdef  AbsCor_00_00_28
    if( emcTrk->e5x5() < 0.015 ) continue;
#endif

    double etof=0;
    if( usetof && itTrk->isTofTrackValid() ) {
      const std::vector<RecTofTrack* >& recTofTrackVec = itTrk->tofTrack();
      if( !recTofTrackVec.empty() ) etof = recTofTrackVec[0]->energy();
      if( etof>100. ) etof = 0;
    }

    if( selector->Verbose() ) {
      cout << " AbsCor:: etof = " << etof << endl;
    }

    double energyC;
#ifdef  AbsCor_00_00_28
    energyC = emcTrk->energy() + etof;
#endif

#ifdef  AbsCor_00_00_36
    double e5x5=emcTrk->e5x5();

    Identifier id(emcTrk->cellId());

    unsigned int module, ntheta, nphi;
    module = EmcID::barrel_ec(id);
    ntheta = EmcID::theta_module(id);
    nphi = EmcID::phi_module(id);

    // id=EmcID::crystal_id(module,ntheta,nphi);

    unsigned int thetaModule = EmcID::theta_module(id);
    unsigned int phiModule = EmcID::phi_module(id);

    int thetaId;
    if (module==0||module==2) thetaId = thetaModule;
    if (module==1 && thetaModule<=21) thetaId = thetaModule + 6;
    if (module==1 && thetaModule>21)  thetaId = 43 - thetaModule + 6;

    if ( MCuseTof ) {
      energyC=ECorrMC(e5x5+etof,thetaId);
    } else {
      energyC=ECorrMC(e5x5,thetaId);
    }
#endif

    if( selector->Verbose() ) {
      cout << " AbsCor:: energyC= " << energyC << endl;
    }

    double lnEcor=1.0;
#ifdef  AbsCor_00_00_36
    if ( dopi0Cor ) {
#endif
      if ( runNo > 0 && dodatacor ) {
        double lnE = std::log(energyC);
#ifdef  AbsCor_00_00_36
        if ( energyC>1.0 ) lnE=std::log(1.0);
        if ( energyC<0.07 ) lnE=std::log(0.07);
#endif
        lnEcor = ai[0] + lnE*(ai[1] + lnE*(ai[2]+lnE*ai[3]));
      }
#ifdef  AbsCor_00_00_36
    }
#endif
    if( lnEcor < 0.5 ) continue;

#ifdef  AbsCor_00_00_36
    // Nefedov: function SuppressHotCrystals() is available in BEAN
    //          regardless of version
    //
    // If it is "hot", return "9999" (Hajime, Jul 2013)
    if ( hotcellmask ) {
      for ( int ih = 0; ih < 10; ih++ ) {
        if ( (hrunstart[ih] == -1) ||
             (hrunend[ih] == -1) ||
             (hcell[ih] == -1) ) continue;
        if ( (abs(runNo) < hrunstart[ih]) ||
             (abs(runNo) > hrunend[ih]) ) continue;

        if ( emcTrk->cellId() == hcell[ih] ) {
          emcTrk->setStatus(9999);
        }
      }
    }
#endif

    double enecor=1.;
    // remove part for "edgecor"

    double energyCC = energyC/(lnEcor*enecor);
    emcTrk->setEnergy(energyCC);

    if( selector->Verbose() ) {
      cout << " AbsCor:: energyCC = " << energyCC << endl;
    }
  }
}


#ifdef  AbsCor_00_00_36
//-----------------------------------------------------------------------------
//The following function is copied from PhotonCor/McCor
double AbsCor::ECorrMC(double eg, double theid) const
//-----------------------------------------------------------------------------
{
  double Energy5x5=eg;
  if(eg<E25min(int(theid))) eg=E25min(int(theid));
  if(eg>E25max(int(theid))) eg=E25max(int(theid))-0.001;

  if(theid<=0)theid=0.001;
  if(theid>=27)theid=26.999;
  Float_t einter = eg + 0.00001;
  Float_t tinter = theid+0.0001;
  //cout<<"inter="<< einter<<"   "<<tinter<<endl;
  double ecor=dt->Interpolate(einter,tinter);
  // cout<<"ecor="<<ecor<<endl;
  if(!(ecor))return Energy5x5;
  if(ecor<0.5)return Energy5x5;
  double EnergyCor=Energy5x5/ecor;
  return EnergyCor;
}


//-----------------------------------------------------------------------------
// Get energy error
double AbsCor::ErrMC(double eg, double theid) const
//-----------------------------------------------------------------------------
{
  if(eg<E25min(int(theid))) eg=E25min(int(theid));
  if(eg>E25max(int(theid))) eg=E25max(int(theid))-0.001;
  if(theid<=0)theid=0.001;
  if(theid>=27)theid=26.999;
  Float_t einter = eg + 0.00001;
  Float_t tinter = theid+0.0001;
  double err=dtErr->Interpolate(einter,tinter);
  return err;
}
#endif
