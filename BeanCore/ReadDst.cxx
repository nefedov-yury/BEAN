//////////////////////////////////////////////////////////////////////////
//                                                                      //
// ReadDst                                                              //
//                                                                      //
// Primary class to read DST                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cstdlib>
#include <csignal>
#include <ctime>

volatile sig_atomic_t bean_termination;  // signal handler

#include <RVersion.h> // Root version
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <TProofOutputFile.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TH1.h>
#include <TH2.h>
#include <TProof.h>
#include <TObjectTable.h>
#include <TMacro.h>
#include <TFileCollection.h>
#include <TEntryList.h>

// if not, gproof->print redirection will be used.
//~ #define USE_MASTER_INFOHACK

#ifdef USE_MASTER_INFOHACK
    #include "MasterInfoHack.h"
#endif

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"

#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"
#include "RootEventData/TJobInfo.h"

#include "DstEvtRecTracks.h"
#include "TMergeableMap.h"
#include "ReadDst.h"

using namespace std;

//--------------------------------------------------------------------
ReadDst::ReadDst(): DstFormat()
//--------------------------------------------------------------------
{
   bean = 0;

   current_entry = -1;
   selected_entries = 0;

   m_evtRecTrkCol = new TObjArray();

   T_select = 0; // to save selected events
   f_select = 0;
   fp_select = 0;
   n_select_events = 0;

   n_events = 0;
   start_time = time(0);
   dirMap = new TMergeableMap();
   return;
}

//--------------------------------------------------------------------
ReadDst::~ReadDst()
//--------------------------------------------------------------------
{
   delete m_evtRecTrkCol;
}

//--------------------------------------------------------------------
bool ReadDst::LoadConfig(Bean* _bean)
//--------------------------------------------------------------------
{
   if( bean ) {
     cout << "ReadDst::LoadConfig: Warning: config is already loaded."
          << endl;
   }

   if( _bean ) {
     bean = _bean;
   } else {
     bean = (Bean *) fInput->FindObject("Bean");
     if( !bean ) {
       cout << "ReadDst::LoadConfig: ERROR: No Bean in fInput!"
            << endl;
     }
   }

   return (bean != 0);
}

//--------------------------------------------------------------------
bool ReadDst::Verbose() const
//--------------------------------------------------------------------
{
   if( !bean ) {
     cout << "ReadDst::Verbose: Warning: config is not yet loaded."
          << endl;
     return true;
   }

   return bean->Verbose();
}

//--------------------------------------------------------------------
string ReadDst::GetBaseDir() const
//--------------------------------------------------------------------
{
   if( !bean ) {
     cout << "ReadDst::GetBaseDir: ERROR: config is not yet loaded."
          << endl;
     exit(1);
   }

   return bean->GetBaseDir();
}

//--------------------------------------------------------------------
string ReadDst::AbsPath(string rel_path) const
//--------------------------------------------------------------------
{
   return (this->GetBaseDir() + "/" + rel_path);
}

//--------------------------------------------------------------------
void ReadDst::SetEntryList(TEntryList* el)
//--------------------------------------------------------------------
{
  if( selected_entries ) {
    cout << " ReadDst::SetEntryList WARNING: you will rewrite list"
            " of selected entries" << endl;
  }
  selected_entries = el;

  if( fChain ) {
    if( fChain->GetEntryList() ) {
      cout << " ReadDst::SetEntryList WARNING: "
         "list of selected entries of the fChain will be overridden"
         << endl;
    }
    fChain->SetEntryList(selected_entries);
  }
}

//--------------------------------------------------------------------
Long64_t ReadDst::GetEntryNumber()
//--------------------------------------------------------------------
{
  return fChain->GetChainEntryNumber(current_entry);
}

//--------------------------------------------------------------------
void ReadDst::SaveEntryInList(TEntryList* el)
//--------------------------------------------------------------------
{
  // current_entry ??
  el->Enter(this->GetEntryNumber(),fChain);
}

//--------------------------------------------------------------------
void ReadDst::Begin(TTree* )
//--------------------------------------------------------------------
{
   // This method is called before looping on the events in the Tree.
   // The user can create his histograms in this function.
   // When using PROOF Begin() is called on the client only. Histogram
   // creation should preferable be done in SlaveBegin() in that case.

   if( !bean ) LoadConfig();

   if( Verbose() ) cout << " ReadDst::Begin() " << endl;

   return;
}

//--------------------------------------------------------------------
void ReadDst::SlaveBegin(TTree* )
//--------------------------------------------------------------------
{
   // The SlaveBegin() function is called after the Begin() function.
   // This method is called on each PROOF worker node.
   // The user can create his histograms in this method.
   // In local mode this method is called on the client too.

   // The Init() function is called after the SlaveBegin() function
   // therefore the fChain is not accessible here on PROOF-workers

   if( !bean && !LoadConfig() ) {
      exit(1);
   }

   if( Verbose() ) {
      cout << " ReadDst::SlaveBegin() " << endl;
//      if(fChain ) fChain->Print();
      cout << " ... " << endl;
   }

   bean->LoadUserFcns();

   string new_dst_file = bean->DstFile();
   if( new_dst_file.size() > 0 ) {
      if( !bean->IsProof() ) {
         f_select = TFile::Open(new_dst_file.c_str(),"RECREATE");
         if( !f_select ) {
            cout  << " ReadDst::SlaveBegin:: ERROR can not open file: "
                  << new_dst_file << endl;
            exit(1);
         }
      } else {

         if ( bean->DstFileIsDataset() ) {
            fp_select = new TProofOutputFile(bean->DstFileName().c_str(),
                  TProofOutputFile::kDataset,
                  TProofOutputFile::kRegister |
                  TProofOutputFile::kVerify |
                  TProofOutputFile::kOverwrite );

            fp_select->GetFileCollection()->SetDefaultTreeName("/Event");

         } else {
            // here we are using old constructor
            fp_select =
               new TProofOutputFile(bean->DstFileName().c_str(), "M" );

            if (! bean->DstFileIsLocal() ) {
               fp_select->SetOutputFileName(bean->DstFile().c_str());
            }
         }


         f_select = fp_select->OpenFile("RECREATE");

         if( !f_select ) {
            cout << " ReadDst::SlaveBegin:: ERROR can not open "
                 << fp_select->GetDir() << "/" << fp_select->GetFileName()
                 << endl;
            exit(1);
         }
      }
   }

   UserStartJob();

   return;
}

//--------------------------------------------------------------------
void ReadDst::Init(TTree *tree)
//--------------------------------------------------------------------
{
   // The Init() function is called when the selector needs to
   // initialize a new tree or chain. Typically here the branch
   // addresses and branch pointers of the tree will be set.

   DstFormat::Init(tree);

   if( selected_entries ) {
      // file names (in fChain and in "selected_entries")
      // are expanded  to full path names
      fChain->SetEntryList(selected_entries);
   }
}

//--------------------------------------------------------------------
Bool_t ReadDst::Notify()
//--------------------------------------------------------------------
{
   // The Notify() function is called when a new file is opened.
   // This can be either for a new TTree in a TChain or when
   // a new TTree is started when using PROOF.
   // The return value is currently not used.

   // if( fChain ) {
     // fChain->SetParallelUnzip(kFALSE);
     // fChain->SetCacheSize(-1);
   // }

   DstFormat::Notify();

   if( Verbose() ) {
     cout << " ReadDst::Notify() " << endl;
   }

   if( f_select && T_select == 0 ) {
     // Note that only active branches are copied !
     // MUST be called DstFormat::Init() before this line
     T_select = fChain->CloneTree(0);
     T_select->SetDirectory(0);
     T_select->SetDirectory(f_select);
     // T_select->SetAutoSave(); // use default !
     // T_select->SetBasketSize("*",1024);
     T_select->SetAutoFlush(100);

     T_select->AutoSave();
   }

   // The cloned tree stays connected with this tree until this tree
   // is deleted. The ProofPlayer::Proccess(dataset,selector)
   // creates new tree for each element of the dataset. Therefore
   // we MUST reassign all branch addresses of T_select every new file:

   if( bean->IsProof() && T_select ) SetBranches(T_select);

   return kTRUE;
}

//--------------------------------------------------------------------
Bool_t ReadDst::Process(Long64_t entry)
//--------------------------------------------------------------------
{
   // This method is called to process an event. It is the user's
   // responsibility to read the corresponding entry in memory
   // (may be just a partial read).
   // Once the entry is in memory one can apply a selection and if the
   // event is selected histograms can be filled.
   //
   // The processing can be stopped by calling Abort()
   // in this case (root6.26) termination functions are not called

   current_entry = entry;

   ClearClasses(); // must be before "GetEntry"

   Int_t nbytes=0;

   // WARNING: entry is always the local entry number
   //          in the current tree!
   nbytes += fChain->GetTree()->GetEntry(entry);

   // debug print
   if( Verbose() ) {
     cout << "ReadDst::Process() Entry# " << entry
          << " nbytes= " << nbytes << endl;
     if( entry == 0 ) Show(entry);
   }

   if( bean->IsDump() ) {
     PrintHeader();
     PrintDst();
     PrintEvtRec();
     PrintMcEvent();
     PrintTrigEvent();
     PrintDigiEvent();
     PrintHltEvent();
#if BOSS_VER >= 661
     PrintEvtNavigator();
#endif
   }

   CreateEvtRecTrkCol();

   if( UserEvent(T_select) ) {
     T_select->Fill();
     n_select_events++;
     if( Verbose() ) {
       cout << " + n_select_events= " << n_select_events << endl;
     }
   }

   n_events += 1;
   if( bean_termination != 0 ) {
      printf("User signal '%s' had been received\n",
            strsignal(bean_termination));
      this->Abort("Stop the event loop...");
   }

   return kTRUE;
}

//--------------------------------------------------------------------
void ReadDst::SlaveTerminate()
//--------------------------------------------------------------------
{
   // This method is called at the end of the loop on all PROOF worker
   // nodes. In local mode this method is called on the client too.

   if( Verbose() )  {
     cout << " ReadDst::SlaveTerminate() " << endl;
     gObjectTable->Print();
   }

   UserEndJob();

   if( fp_select ) {            // Proof only
     Bool_t cleanup = kFALSE;
     TDirectory *savedir = gDirectory;
     if( Verbose() ) {
       cout << " ReadDst::SlaveTerminate n_select_events= "
            << n_select_events << endl;
     }
     if( n_select_events > 0 ) {
       f_select->cd();
       WriteJobInfo();
       T_select->Write();
       if( Verbose() ) {
         cout << " events saved in file: " << f_select->GetName() << endl;
       }
       fOutput->Add(fp_select);
     } else {
       cleanup = kTRUE;
     }
// #warning  T_select->SetDirectory(0); leads to crash.

     gDirectory = savedir;
     f_select->Close();

     // Cleanup, if needed
     if( cleanup ) {
       TUrl uf(*(f_select->GetEndpointUrl()));
       SafeDelete(f_select);
       gSystem->Unlink(uf.GetFile());
       SafeDelete(fp_select);
       if( Verbose() ) {
         cout << " cleanup file: " << uf.GetFile() << endl;
       }
     }
   }

    #ifdef USE_MASTER_INFOHACK
       if( bean->DstFileName().size() > 0 &&   bean->DstFileIsLocal() ){
                    fOutput->Add(new MasterInfoHack());
       }
    #endif

   dirMap->SetName("DirectoryMap");
   fOutput->Add( dirMap);
   return;
}

//--------------------------------------------------------------------
void ReadDst::Terminate()
//--------------------------------------------------------------------
{
   // This method is called at the end of the loop on all events.
   // When using PROOF Terminate() is call on the client only.
   // Typically one performs the fits on the produced histograms
   // or write the histograms to file in this method.

   if( Verbose() ) {
     cout << " ReadDst::Terminate() " << endl;
     if( fChain ) fChain->PrintCacheStats("a");
   }

   if( !bean->IsProof() && T_select ) {
     f_select->cd();
     T_select->Write();
     cout << " ReadDst::Terminate " << n_select_events
          << " events saved in file: " << f_select->GetName() << endl;

     WriteJobInfo();

     f_select->Close();
     SafeDelete(f_select);
   }

   if (!bean->IsProof() ) {
       time_t time_taken = time(0)  - start_time ;
       cout << " Performance summary: " << n_events
       << " events processed in " << time_taken << "s: ~" ;

       if (time_taken > 0) {
           cout << n_events / ( time_taken );
           }  else {
                   cout << "inf";
           }

       cout  << " events/sec" << endl;
   }

   Save_histo();


   // retrieving output DST file from master
   if( bean->IsProof() && !bean->Proof()->IsLite()
       && bean->DstFileName().size() > 0 && bean->DstFileIsLocal() ) {

       string serverDstLocation;
       bool use_xrootd = 0;

       #ifdef USE_MASTER_INFOHACK
           MasterInfoHack* masterInfoHack =
              (MasterInfoHack*) fOutput->FindObject("MasterInfoHack" );

       if (!( masterInfoHack->MasterWorkdir().size() == 0)) {
            string workdir = masterInfoHack->MasterWorkdir();
       #else
            string workdir = bean->GetProofMasterWorkdir();
       #endif


       // getfile return code, -1 for error (getfile is not supported
       // by server)
       int getfile_rc = 0;

       #if ROOT_VERSION_CODE >= ROOT_VERSION(5,25,2)
       if (!bean->ProofXrdOutput()) {
                           // download file using GetFile (default)
           string relative_path =  workdir + "/"+ bean->DstFileName();
           getfile_rc = bean->Proof()->GetManager()->
                GetFile( relative_path.c_str(), bean->DstFile().c_str() );
       }
       #else
           getfile_rc = -1; // there is no getfile functionality on client
       #endif

       if ( (bean->ProofXrdOutput())  || (getfile_rc == -1)  ) {
          // construct URL using master host, master workdir and Dst
          // file name
          serverDstLocation = "root://" + bean->ProofClr() + "/"
                  + workdir + "/"+ bean->DstFileName();
          use_xrootd = 1;
       }

       #ifdef USE_MASTER_INFOHACK
       } else {
         // in some rare cases (one worker) merging doesn't happen.
         // so we use worker dir to get file

           // FIXME: this will not work when master node has no worker
           // nodes

           TProofOutputFile* fProofFile = (TProofOutputFile*)
                        fOutput->FindObject(bean->DstFileName().c_str() );

           string fProofFileDir = fProofFile->GetDir();
           serverDstLocation =  fProofFileDir + "/" + bean->DstFileName();
           use_xrootd = 1;
           fProofFile->Print();
           cout << "WARNING: Merging of output file didn't occur."
                   " Will try to use worker dir to fetch file."
                   " This expected not to work if master node has"
                   " no workers nodes" << endl;
       }
       #endif


       if (use_xrootd) {
          cout << serverDstLocation << endl;
          //   sleep(100);
          TFile::Cp(serverDstLocation.c_str(), bean->DstFile().c_str() );
       }
   }
}

//--------------------------------------------------------------------
void ReadDst::UserStartJob()
//--------------------------------------------------------------------
{
   // User Start-Job functions:
   const VecUF& Ufn_start = bean->GetStartJobFns();
   VecUF::const_iterator iuf = Ufn_start.begin();

   // run all user functions
   for( ;iuf != Ufn_start.end(); iuf++) {
      ssfn user_func = (ssfn) (*iuf);
      user_func(this);
   }
}

//--------------------------------------------------------------------
bool ReadDst::UserEvent(TTree* T)
//--------------------------------------------------------------------
{
   const VecUF& Ufn_event = bean->GetUserEventFns();
   VecUF::const_iterator iuf = Ufn_event.begin();

   bool save_this = false;

   // run all user functions
   for( ;iuf != Ufn_event.end(); iuf++) {
      pufn user_func = (pufn) (*iuf);
      bool ret = user_func(this,
            m_TEvtHeader,m_TDstEvent,m_TEvtRecObject,
            m_TMcEvent,m_TTrigEvent,m_TDigiEvent,m_THltEvent);
      if( T && ret ) {
        save_this = true;
      }
   }

   if (T && Ufn_event.empty())
      save_this = true;
   return save_this;
}

//--------------------------------------------------------------------
void ReadDst::UserEndJob()
//--------------------------------------------------------------------
{
   const VecUF& Ufn_end = bean->GetEndJobFns();
   VecUF::const_iterator iuf = Ufn_end.begin();

   // run all user functions
   for( ;iuf != Ufn_end.end(); iuf++) {
      ssfn user_func = (ssfn) (*iuf);
      user_func(this);
   }
}

//--------------------------------------------------------------------
void ReadDst::CreateEvtRecTrkCol()
//--------------------------------------------------------------------
{
   if( Verbose() ) {
     cout << " CreateEvtRecTrkCol " << m_evtRecTrkCol->GetEntriesFast()
          << endl;
   }
   m_evtRecTrkCol->Delete();
   if (m_TEvtRecObject) {
      const TObjArray* m_evtRecTrackCol =
         m_TEvtRecObject->getEvtRecTrackCol();
      TIter evtRecTrackIter(m_evtRecTrackCol);
      TEvtRecTrack* evtRecTrack = 0;
      while(
            (evtRecTrack = (TEvtRecTrack*)evtRecTrackIter.Next()) != 0
           ) {
        m_evtRecTrkCol->AddLast(
              new DstEvtRecTracks(evtRecTrack,m_TDstEvent) );
      }
   }
}

// To add objects SHOULD be implemented the method 'Merge()'
//--------------------------------------------------------------------
void ReadDst::RegInDir(const VecObj* hst, const char* dir)
//--------------------------------------------------------------------
{
   if( Verbose() ) cout << " ReadDst::RegInDir(VecObj*) " << endl;

   // add the prefix (=="$dir_") to the name
   // so that the names is unique in the current directory
   string prefix;
   if ( dir ) {
      prefix = string(dir) + "_";
   }
   VecObj::const_iterator it1 = hst->begin();
   for( ; it1 != hst->end(); it1++) {
      if(*it1) {
        if( Verbose() ) cout << " + " << (*it1)->GetName() << endl;
        if ( dir ) {
           string nn = prefix+string((*it1)->GetName());
           (*it1)->SetName(nn.c_str());
        }
        CheckDupName(*it1);
        fOutput->Add(*it1);
        if( dir ) {
          dirMap->Add(
                new TObjString((*it1)->GetName()), new TObjString(dir));
        }
      }
   }
}

//--------------------------------------------------------------------
void  ReadDst::RegInDir(const TList* hst, const char* dir )
//--------------------------------------------------------------------
{
   if( Verbose() ) cout << " ReadDst::RegInDir(TList*) " << endl;
   TIter next(hst);
   VecObj hst_vec;
   TNamed * obj;
   while ((obj = (TNamed*) next())) {
       hst_vec.push_back(obj);
   }

   RegInDir(hst_vec, dir);
}

//--------------------------------------------------------------------
void ReadDst::Save_histo() const
//--------------------------------------------------------------------
{
   if( fOutput->GetSize() == 0 ) return;

   string hst_file = bean->HstFile();
   TFile* c_out = TFile::Open(hst_file.c_str(),"RECREATE");
   if( !c_out ) {
     cout << " ReadDst::Save_histo::ERROR can not open " << hst_file
          << endl;
     return;
   }
   c_out->cd();

   int nhisto = 0;
   int ntrees = 0;

   TMergeableMap* directoryMap =
                (TMergeableMap*) fOutput->FindObject("DirectoryMap");
   //~ cout << "directoryMap:" << directoryMap << endl;

   TIter next(fOutput);
   TObject* obj = 0;
   while( (obj = next()) ) {
      if( obj->IsA() == TProofOutputFile::Class() ) {
         ; // this object is for internal use
      #ifdef USE_MASTER_INFOHACK
      } else if( obj->IsA() == MasterInfoHack::Class()) {
         ; // this object is for internal use
      #endif
      } else if( obj == directoryMap ) {
      ; // this object is for internal use
      } else if( obj->InheritsFrom(TObject::Class()) &&
            obj->IsA()->GetMethodWithPrototype("GetName","") )
      {
         TNamed* hst  = (TNamed*) obj;

         if( obj->InheritsFrom("TH1") ) { // this is histograms
            if( ((TH1*) obj) -> GetEntries() < 0.01 ) {
               if( Verbose() ) {
                  cout << " skip empty histogram:  " << hst->GetName()
                       << endl;
               }
               continue;
            }
            nhisto++;
         } else if( obj->InheritsFrom("TTree") ) { // this is tree
            if( ((TTree*) obj) -> GetEntries() < 0.01 ) {
               if( Verbose() ) {
                  cout << " skip empty tree:  " << hst->GetName()
                       << endl;
               }
               continue;
            }
            ntrees++;
         }

         // save "objects" in the directory if it was specified
         // in this case remove the prefix (=="$dir_")
         // from the names (see RegInDir(VecObj*))
         if (directoryMap) {
            TObjString* dir_name =
               (TObjString*) directoryMap->GetValue(hst->GetName());
            size_t pref_size = string(dir_name->GetName()).size();
            if ( pref_size > 0 ) {
               string nn = string(hst->GetName()).substr(pref_size+1);
               hst->SetName(nn.c_str());
            }

            if (dir_name) {
               if ( !c_out->GetDirectory( dir_name->GetString() ) ) {
                  c_out->mkdir( dir_name->GetString() );
               }
               c_out->cd( dir_name->GetString() );
            } else {
               c_out->cd();
            }
         } else {
            c_out->cd();
         }

         // to store TMap/TCollection/etc as single object
         obj->Write(NULL, TObject::kSingleKey);

     } else {

       cout << " ReadDst::Save_histo::Warning: "
            << " an object in fOutput doesn't have a name "
            << typeid(*obj).name() << endl;
       TDirectory* dir = gDirectory;
       c_out->cd();
       obj->Write();
       dir->cd();

     }

   }

   c_out->Close();
   if( Verbose() ) {
     cout << " Save " << nhisto << " histograms " << endl
          << "  and " << ntrees << " trees in file " << hst_file
          << endl;
   }
}

//--------------------------------------------------------------------
void ReadDst::CheckDupName(TObject *obj)
//--------------------------------------------------------------------
{
   if( obj ) {
     TObject *org = fOutput->FindObject(obj->GetName());
     if( org && (org != obj) ) {
       cout << " ReadDst::CheckDupName:: FATAL ERROR " << endl
            << " object with name: " << obj->GetName()
            << " already in the fOutput list " << endl;
       exit(1);
     }
   }
}

//--------------------------------------------------------------------
void ReadDst::WriteJobInfo()
//--------------------------------------------------------------------
{
   // copy&pasted from BOSS
   TJobInfo* jobInfo = new TJobInfo;
   TTree* m_jobInfoTree = new TTree("JobInfoTree","Job info");
   m_jobInfoTree->Branch("JobInfo",&jobInfo);

   string m_bossVer = "BEAN";
   jobInfo->setBossVer(m_bossVer);

   // fill all members of TJobInfo with something
//    jobInfo->setDecayOptions(string(10,'+'));
//    jobInfo->addJobOptions(string(10,'+'));
//    vector<int> m_totEvtNo(2,1);
//    jobInfo->setTotEvtNo(m_totEvtNo);

   m_jobInfoTree->Fill();
   f_select->cd();
   m_jobInfoTree->Write();
}

