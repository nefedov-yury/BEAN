#ifndef _ReadDst_h
#define _ReadDst_h

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// ReadDst                                                              //
//                                                                      //
// Primary class to read DST                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "DstFormat.h"
#include "Bean.h"

// class Bean;
class ReadDst;
class TProofOutputFile;
class TObject;
class TNamed;
class TEntryList;

class TMergeableMap;

// pointers on user functions:
// - start/end job
typedef void (*ssfn) (ReadDst*);
// - event function
typedef bool (*pufn) (ReadDst*,
                      TEvtHeader*,TDstEvent*,TEvtRecObject*,
                      TMcEvent*,TTrigEvent*,TDigiEvent*,THltEvent*);

typedef std::vector<TNamed* > VecObj;

class ReadDst : public DstFormat
{
public :
                        ReadDst();
   virtual             ~ReadDst();

   bool                 LoadConfig(Bean* _bean = 0);
   bool                 Verbose() const override;
   void                 SetVerbose() {bean->SetVerbose();}
   void                 SetSilent()  {bean->SetSilent();}
   std::string          GetBaseDir() const;
   std::string          AbsPath(std::string rel_path) const;

   const TObjArray*     GetEvtRecTrkCol() const {return m_evtRecTrkCol;}

   void                 SetEntryList(TEntryList* el);
   Long64_t             GetEntryNumber();
   void                 SaveEntryInList(TEntryList* el);

   // TSelector functions:
   void                 Begin(TTree* ) override;
   void                 SlaveBegin(TTree* ) override;
   void                 Init(TTree *tree) override;
   Bool_t               Notify() override;
   Bool_t               Process(Long64_t entry) override;
   void                 SlaveTerminate() override;
   void                 Terminate() override;

   // -- call user functions
   void                 UserStartJob();
   bool                 UserEvent(TTree* T);
   void                 UserEndJob();
   void                 CreateEvtRecTrkCol();

   // -- histogramming
   void                 RegInDir(const VecObj* hst, const char* dir = 0);
   void                 RegInDir(const VecObj& hst, const char* dir = 0){
      RegInDir(&hst,dir);
   }
   void                 RegInDir(const TList* hst, const char* dir = 0);
   void                 RegInDir(const TList& hst, const char* dir = 0){
      RegInDir(&hst,dir);
   }
   void                 Save_histo() const;

private :
   Bean*                bean; //-> all configuration parameters are here

   Long64_t             current_entry;
   TEntryList*          selected_entries;

   // -- track linking
   TObjArray*           m_evtRecTrkCol;

   // -- for selected events
   TTree*               T_select;
   TFile*               f_select;
   TProofOutputFile*    fp_select;
   int                  n_select_events;

   // -- time measuring
   time_t               start_time;
   Int_t                n_events;

   // histogram to directory map
   TMergeableMap*       dirMap;

   // -- internal functions:
   void                 CheckDupName(TObject *obj);
   void                 WriteJobInfo();

   // ClassVersionID=0 because we don't need object I/O
//    ClassDef(ReadDst,0); // Primary class to read DST
   // if class definition use `override` keyword
   ClassDefOverride(ReadDst,0); // Primary class to read DST
};
#endif
