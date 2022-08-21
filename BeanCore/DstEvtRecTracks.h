#ifndef DstEvtRecTracks_H
#define DstEvtRecTracks_H

//
// This is transitional class to use Dst-tracks instead of
// "EvtRecEvent/EvtRecTrack.h"
//

#include <vector>

#include <TObject.h>

#include "RootEventData/TDstEvent.h"

#include "BMdcTrack.h"
#include "BMdcKalTrack.h"
#include "BTofTrack.h"
#include "BExtTrack.h"
#include "BEmcTrack.h"


class TEvtRecTrack;

typedef         BMdcTrack       RecMdcTrack;
typedef         TMdcDedx        RecMdcDedx;
typedef         BMdcKalTrack    RecMdcKalTrack;
typedef         BTofTrack       RecTofTrack;
typedef         BExtTrack       RecExtTrack;
typedef         BEmcTrack       RecEmcShower;
typedef         TMucTrack       RecMucTrack;

//
// SmartRefVector<RecTofTrack> -> std::vector<RecTofTrack*>& ;
//

class DstEvtRecTracks : public TObject
{
 public:
                DstEvtRecTracks(TEvtRecTrack* rec_trk,
                                TDstEvent* dst_event);
               ~DstEvtRecTracks();

  int           trackId() const;

  // TMdcTrack
  bool          isMdcTrackValid() const {return (mdc_trk != 0);}
  RecMdcTrack*  mdcTrack() const        {return mdc_trk;}

  // TMdcDedx
  bool          isMdcDedxValid() const  {return (mdc_dedx != 0);}
  RecMdcDedx*   mdcDedx() const         {return mdc_dedx;}

  // TMdcKalTrack
  bool          isMdcKalTrackValid() const
                                        {return (mdc_kal_trk != 0);}
  RecMdcKalTrack* mdcKalTrack() const   {return mdc_kal_trk;}

  // TofTrack
  bool          isTofTrackValid() const {return (!tof_trk.empty());}
  const std::vector<RecTofTrack* >&
                tofTrack() const        {return tof_trk;}

  // TExtTrack
  bool          isExtTrackValid() const {return (ext_trk != 0);}
  RecExtTrack*  extTrack() const        {return ext_trk;}

  // TEmcTrack
  bool          isEmcShowerValid() const {return (emc_trk != 0);}
  RecEmcShower*  emcShower() const       {return emc_trk;}

  // TMucTrack
  bool          isMucTrackValid() const {return (muc_trk != 0);}
  RecMucTrack*  mucTrack() const        {return muc_trk;}

 private:
  TEvtRecTrack*                 m_rec_trk;

  RecMdcTrack*                  mdc_trk;
  RecMdcDedx*                   mdc_dedx;
  RecMdcKalTrack*               mdc_kal_trk;
  RecExtTrack*                  ext_trk;
  std::vector<RecTofTrack* >    tof_trk;
  RecEmcShower*                 emc_trk;
  RecMucTrack*                  muc_trk;

};

typedef DstEvtRecTracks EvtRecTrack;

#endif
