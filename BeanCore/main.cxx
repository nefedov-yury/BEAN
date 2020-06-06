#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <csignal>
#include <cstring>
#include <cerrno>
#include <climits>

#if defined (__unix__) || defined (__APPLE__)
// __unix__ is defined by compiler for Unix and __APPLE__ for MacOS
#  include <getopt.h>
#elif defined _WIN32    // _Win32 is defined for 32 or 64 bit Windows
#  include <windows.h>
#  include "win_getopt/getopt.h"
#endif

#include <TEnv.h>
#include <TSystem.h>
#include <TProof.h>
#include <TROOT.h>
#include <TChain.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TEnv.h>
#include <TProofLog.h>
#include <TFileCollection.h>

// #include <TObjectTable.h> // debug memeory leak

#include "Bean.h"
#include "ReadDst.h"
#if defined (__unix__)
   #include "DatabaseSvc/DatabaseSvc.h"
   #include "unix_cp.h"
#endif

using namespace std;

#if defined _WIN32
void win_exit(){::ExitProcess(0);}
#endif

// local:
static Bean* bean = 0;
#define TREE_CACHE_SIZE 10000000  //10 Mbytes

static void termination_handler(int isig);
static void set_user_termination();
static void segfault_handler(int isig);
static void set_proof_termination();
static void check_upload_enable(TProof* proof, const char* pname);

static bool save_proof_logs;

//------------------------------------------------------------------------------
void graceful_proof_exit()
//------------------------------------------------------------------------------
{
   if (gProof) {
      if (save_proof_logs) {
         gProof->GetManager()->GetSessionLogs()->Save("*", "bean_proof.log");
      }

      gProof->Close();
      delete gProof->GetManager();
   }
}

//------------------------------------------------------------------------------
void Usage(int argc, char **argv)
//------------------------------------------------------------------------------
{
   bool verbose = (bean) ? bean->Verbose() : false;

   cout << endl;
   cout << " Usage: " << argv[0] << " [ -option(s)] dst_file(s)" << endl;

   if( !verbose ) cout
        << "  Note: This is a short note, run \""
        << argv[0] << " -v\" for more information" << endl;
   else cout
        << "  Note: The program has been compiled for BOSS-"
        << BOSS_VER << " version" << endl;

#if defined _WIN32
   if( verbose ) cout
        << "    Important: Note that PROOF & PROOF-Lite is not yet available"
           " for Windows OS" << endl;
#endif

   cout << endl;
   cout << " Arguments:" << endl;

   cout << "  dst_file(s)   input ROOT files (mandatory";
   if( verbose ) cout
        << ", see also -L flag):" << endl;
   else cout
        << "):" << endl;
   cout << "     local_file.root - path to the local file" << endl
        << "     root://user@host/path/to/file.root"
           " - use remote file from xrootd" << endl;
   if( verbose ) cout
        << "     ds://DatasetName - use dataset registered on PROOF cluster"
        << endl;

   cout << endl;
   cout << " Options: " << endl;

   cout << "  -h hst_file   change default name for file of histograms"
        << " (bean_histo.root)" << endl;

   cout << "  -o out_file   name of output ROOT tree file"
        << " (default: no output is written)" << endl;
   if( verbose ) cout
        << "                In PROOF mode the file name could be " << endl
        << "                  1) out_file.root - the file will be"
                                 " merged on master node" << endl
        << "                     and then fetched back to client." << endl
        << "                  2) root://user@host/path/to/file.root"
                                 " - to save output" << endl
        << "                     at remote xrootd server" << endl
        << "                  3) ds://DatasetName"
                                 " - to register output as dataset " << endl
        << "                     on PROOF cluster (experts prefer) " << endl;

   cout << "  -u Uname      add user function user/Uname.cxx" << endl;
   if( verbose ) cout
        << "                (could be specified more than once) " << endl;

   cout << "  -v            set verbosity on (try \"" << argv[0] << " -v\")"
        << endl;

   cout << "  -D            event dump (detailed printout of DST content) "
        << endl;

   cout << "  -N num        process not more than \"num\" events in total"
        << endl;

   cout << "  -l            use PROOF-Lite (local PROOF optimized for"
        << " multi-core systems)" << endl;

   if( verbose ) cout
        << "  -p proof_clr  use the specified PROOF cluster" << endl;

#if defined (__unix__)
   cout << "  -S[path]      inicialize Sqlite database. The optional"
           " argument [path]" << endl
        << "                is used to specify the location where you want"
           " to create" << endl
        << "                a copy of the database to avoid locking"
           " problems on NFS." << endl
        << "                For example: \"-S/tmp\" or \"-S/dev/shm\""
        << endl;
#endif

   if( verbose ) {
   cout << endl;
   cout << " Advanced options:" << endl;

   cout << " * for input files:" << endl;
   cout << "  -L filelist   read dst_file(s) names from filelist, one per line"
        << endl;
   cout << "  -r prefix     add prefix to each file name" << endl;

   cout << " For output files:" << endl;
   cout << "  -x            use xrootd to fetch output file from master"
           " (with -o)" << endl
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,25,2)
        << "                (default: PROOF sandbox access will be used)"
#else
        << "                (This is the only avilable method."
           " In order to use PROOF" << endl
        << "                sandbox access you must install ROOT >= 5.22/2)"
#endif
        << endl;

   cout << " * for PROOF management: " << endl;
   cout << "  -a params     set PROOF parameter "
           "(\"valgrind\", \"workers=42\", etc)" << endl;
   cout << "  -A key=value  set PROOF INPUT parameter (TProof::SetParameter()) "
        << endl;

   cout << "  -g            save proof error logs to bean_proof.log file"
        << endl;
   cout << "  -d workers    disable some workers, comma-separated" << endl;
   cout << "  -E var=value  add variable to PROOF environment,"
           " for example:" << endl
        << "                PROOF_WRAPPERCMD=valgrind_opts:--leak-check=full"
        << endl;
   cout << "  -V level      set PROOF log level" << endl;

   cout << " * others: " << endl;
   cout << "  -e var=value  add variable to gEnv,"
           " for example: XProof.Debug=2" << endl;
   cout << endl;
   } // end of if( verbose )

   graceful_proof_exit();
   exit(0);
}

bool str2int(const char * str, long & val) {
   char* endptr;
   errno = 0;    // To distinguish success/failure after call
   val = strtol(str, &endptr, 10);

   // Check for various possible errors
   if( ( errno == ERANGE && (val == LONG_MAX||val == LONG_MIN))
       || (errno != 0 && val == 0))
   {
       return false;
   }

   if (endptr == str || *endptr != '\0') {
      return false;
   }

   return true;
}

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
//-----------------------------------------------------------------------------
{
   // Switch on synchronization with the standard C streams.
   ios_base::sync_with_stdio(true);

#if defined _WIN32
   // Enable two-digit exponent format
   _set_output_format(_TWO_DIGIT_EXPONENT);

   // suppress the abort message
   _set_abort_behavior( 0, _WRITE_ABORT_MSG );

   // Disable the message box for errors, unrecoverable problems, and so on
   _CrtSetReportMode( _CRT_ERROR, 0 );

   // Disable assertion failures
   _CrtSetReportMode( _CRT_ASSERT, 0 );

   // function win_exit() will run at exit()
   atexit( win_exit );
#endif

   bean = new Bean;

#ifdef BEANBASE
   bean->SetBaseDir(BEANBASE);
#else
   cout << " ERROR: BEANBASE not defined " << endl;
   exit(1);
#endif

   string list_fname = "";
   string prefix ="";

// ====================== PARSE COMMAND OPTIONS ======================
   bool is_error = false;
   save_proof_logs = false;
   bool copy_sqlite = false;
   (void)copy_sqlite; // prevent compiler to complain it is unused

   TString workers_disable_str;
   TString * str_arg;
   TString * str_env_var;
   int proof_log_level = 0;

   std::multimap<std::string, std::string>  proof_params;

   int oc; // option
   while( (oc = getopt(argc,argv,":d:r:L:xh:o:la:A:p:u:vDN:e:gS::E:V:")) != -1 ) {
     switch( oc ) {

     case 'r':  // input files  prefix
                prefix  = optarg;
                break;

     case 'L':  // list of input files
                list_fname  = optarg;
                break;
     case 'h':  // change default file name for histograms
                bean->SetHstFile(optarg);
                break;

     case 'o':  // define output ROOT tree file name
                bean->SetDstFile(optarg);
                break;

     case 'l':  // use PROOF-Lite cluster
                bean->SetProofClr("lite://");
                break;
     case 'x':  // use xrootd
                bean->SetProofXrdOutput(1);
                break;

     case 'p':  // use PROOF cluster
                bean->SetProofClr(optarg);
                break;
     case 'a':  // use PROOF param
                bean->SetProofParam(optarg);
                break;

     case 'u':  // add user function
                bean->AddUserFcn(optarg);
                break;

     case 'v':  // set verbosity on
                bean->SetVerbose();
                break;
     case 'D':  // event dump
                bean->SetEventDump();
                break;

     case 'N':  // process N events
                {
                   long val;
                   if (!str2int(optarg, val) || val < 0) {
                       cout << " option -N requires positive number " << endl;
                       cout << " incorrect argument: " << optarg << endl;
                       is_error = true;
                       break;
                   }

                   bean->SetMaxNumberEvents(val);
                }
                break;

     case 'd':  // disable some workers
                workers_disable_str = optarg;
                break;
     case 'e':
                str_arg = new TString(optarg);
                str_env_var = new TString(((*str_arg)(0, str_arg->First('='))));
                gEnv->SetValue(str_env_var->Data(),
                  (*str_arg)(str_arg->First('=')+1, str_arg->Length() ).Data());
                delete str_arg;
                delete str_env_var;

                break;
     case 'E':
                str_arg = new TString(optarg);
                str_env_var = new TString(((*str_arg)(0, str_arg->First('='))));
                TProof::AddEnvVar(str_env_var->Data(),
                (*str_arg)(str_arg->First('=') + 1, str_arg->Length() ).Data());

                delete str_arg;
                delete str_env_var;

                break;

     case 'A':
                str_arg = new TString(optarg);
                str_env_var = new TString(((*str_arg)(0, str_arg->First('='))));

                proof_params.insert( std::pair<std::string, std::string> (
                   str_env_var->Data(),
                   (*str_arg)(str_arg->First('=') + 1, str_arg->Length() ).Data()
                ));

                delete str_arg;
                delete str_env_var;

                break;


     case 'g':
                save_proof_logs = true;
                break;

     case 'V': // set PROOF log level
               {
                  long val;
                  if (!str2int(optarg, val) || val < 0) {
                     cout << " option -V requires positive number " << endl;
                     cout << " incorrect argument: " << optarg << endl;
                     is_error = true;
                     break;
                  }

                  proof_log_level = val;
               }
               break;

#if defined (__unix__)
     case 'S':  // initialize Sqlite database
                {
                  // cerr << " S-option with optarg= "
                  //     << ((optarg) ? string(optarg) : "not set") << endl;
                  DatabaseSvc* dbs = DatabaseSvc::instance();
                  char dir_name[256];
                  snprintf(dir_name, sizeof(dir_name),
                                "%s/Analysis/DatabaseSvc/dat",
                                (bean->GetBaseDir()).c_str()   );
                  if( optarg ) { // with argument
                    const char* dir_new = copy_dir_temp(dir_name,optarg);
                    dbs->SetDBFilePath(dir_new);
                    copy_sqlite = true;
                  } else {       //  without argument
                    dbs->SetDBFilePath(dir_name);
                  }

                  break;
                }
#endif
     case ':':  // no argument (first character of optstring MUST be ':')
                is_error = true;
                printf(" option `-%c' requires an argument\n",optopt);
                break;

     case '?':  // errors
     default:
                is_error = true;
                printf(" invalid option `-%c'\n",optopt);
                break;
     }
   }

   if( is_error ) {
     Usage(argc,argv);
   }

   vector<char*> file_names;
   for(int iarg = optind; iarg < argc; iarg++) {
     file_names.push_back(argv[iarg]);
   }

   // Read names of files from list_fname
   if( list_fname != "" ) {
      char * cin_fname;
      char cin_fname_buffer[1024];
      int cin_fname_size;

      istream* list_s;
      if( list_fname == "-" ) {
              list_s = &cin;
      } else {
              list_s = new ifstream(list_fname.c_str());
      }

      while (! list_s->eof() && !list_s->fail()) {
            list_s->getline(cin_fname_buffer,1024);
            cin_fname_size = strlen(cin_fname_buffer);
            if( cin_fname_size > 0 ) {
              cin_fname = (char*) malloc( sizeof(char) * (cin_fname_size+1));
              strcpy(cin_fname, cin_fname_buffer);
              file_names.push_back(cin_fname);
            }
      }

      if (list_s != &cin){
              delete list_s;
      }
   }

   if( file_names.empty() ) {
     if( argc > 2 || (argc==2 && !bean->Verbose()) ) {
       // skip this print if no options or only one option "-v"
       cerr << "ERROR: list of input ROOT files is required " << endl;
     }
     Usage(argc,argv);
   }

   //~ if( bean->NUserFns() == 0 && !bean->IsDump() ) {
     //~ cerr << " WARNING: you do not specify user function (option -u): "
          //~ << " use user/UserTest.cxx " << endl;
     //~ bean->AddUserFcn("UserTest");
   //~ }

   if( (!bean->IsProof() ) && ( bean->DstFileIsDataset() )) {
     cout << " ERROR: dataset output is not allowed in no-PROOF mode " << endl;
     Usage(argc,argv);
   }

   if( (!bean->IsProof() ) && (! bean->DstFileIsLocal() )) {
     cerr << " ERROR: remote file output is not allowed in no-PROOF mode " << endl;
     Usage(argc,argv);
   }

   bean->PrintOptions();

// ================== HANDLE THE TERMINATION SIGNALS =================
   if( bean->IsProof() ) set_proof_termination();
   else                  set_user_termination();

// ============================ EXECUTE ==============================

   TChain chain("Event");
   Long64_t nentries = bean->MaxNumberEvents();

   vector<string> dataset_names;

   for(unsigned int nf = 0; nf < file_names.size(); nf++) {
      string filename = prefix + file_names[nf];

      if( bean->Verbose() ) {
        cout << " add file: " << filename << endl;
      }

      // If we are working in proof mode and filename is definetely dataset
      string dataset_name = bean->ParseDatasetName(filename.c_str());
      if ( dataset_name.size() ) {
         if (!bean->IsProof()) {
            cerr << "Dataset input is not supported without PROOF - ignored."
                 << endl;
         } else {
            dataset_names.push_back(dataset_name);
         }
      } else {
         // found ordinary file
         chain.Add(filename.c_str());
      }
   }

   TProof* proof = 0;
   if( bean->IsProof() ) {
      gROOT->SetBatch(true); // Prevents ROOT from enabling graphics

      const string& clr = bean->ProofClr();
      const string& param = bean->ProofParam();

      // It's time for ugly hacks!

      // We use manager to check whether this will be Proof-Lite
      TProofMgr * manager = TProofMgr::Create( clr.c_str() );

      if (manager->IsLite()) {
         #if ROOT_VERSION_CODE < ROOT_VERSION(5,29,1)
            TProof::AddEnvVar("ROOTPROOFLITE", "1");
         #endif

         // If BEAN is built with xlinked ROOT there is no ROOT libraries in
         //                                           the ld.so search PATH.
         // But proofserv.exe need this libraries to work.
         // So in case of ProofLite we should set ld_library_path to
         // appropriate one
         #ifdef ROOTLIBDIR
            string new_library_path;
            new_library_path += ROOTLIBDIR;
            new_library_path +=":$LD_LIBRARY_PATH";
            TProof::AddEnvVar("LD_LIBRARY_PATH", new_library_path.c_str() );
         #endif
      }

      string clr_open = clr + "/?N"; // to force creation of a new session
      //~ TProof::Reset(clr.c_str());

      proof = TProof::Open(clr_open.c_str(), param.c_str());

      gEnv->SetValue("Proof.StatsHist",1);
      gEnv->SetValue("Proof.StatsTrace",1);
      gEnv->SetValue("Root.Stacktrace","no"); // do not print stack trace


      proof->SetLogLevel(proof_log_level);

      // set default proof parameters
      proof->SetParameter("PROOF_SavePartialResults", "1");

      // set user-specified proof parameters

      for (std::multimap<std::string, std::string>::const_iterator it = proof_params.begin();
           it != proof_params.end(); ++it )
      {
         proof->SetParameter(it->first.c_str(), it->second.c_str());
      }


      if ( bean->DstFileIsDataset() ) {
         // TODO: make some option to force overwrite
         if ( proof->GetDataSets(bean->DstFile().c_str())->GetSize() != 0) {
            cerr << " WARNING: dataset " << bean->DstFile() <<
               " already exists on cluster and will be overwriten" << endl;

         }
      }

      // handle datasets
      string dataset_names_string = "";
      if (chain.GetNtrees() == 0 ) {
         //all the data is datasets
         for (unsigned int i = 0; i < dataset_names.size(); ++i) {
            if (i!= 0) {
               dataset_names_string += "|";
            }

            dataset_names_string += dataset_names[i];
            dataset_names_string += "#Event";
         }
      } else {
         // if not, try to retrieve information from proof cluster
         for (unsigned int i = 0; i < dataset_names.size(); ++i) {
            TFileCollection* collection =
                                proof->GetDataSet(dataset_names[i].c_str());
            if ( collection ) {
               chain.AddFileInfoList((TCollection*) collection->GetList());
            } else {
               cerr << "Dataset " << dataset_names[i] << "does not exist"
                    << endl;
            }

         }
      }

      // Disable some workers
      TIter next(workers_disable_str.Tokenize(","));
      TObjString* worker;
      while ((worker = (TObjString*) next())) {
         proof->DeactivateWorker(worker->GetName());
      }


      check_upload_enable(proof, "par/Analysis.par");
      check_upload_enable(proof, "par/BeanUser.par");
      check_upload_enable(proof, 
            ("par/BeanLib_" + to_string(BOSS_VER) + ".par").c_str()
                         );


      bean->SetProof(proof);
      proof->AddInput(bean);





      if ( dataset_names_string.empty() ) {
         // use chain
         chain.SetProof();
         if( !nentries ) {
            chain.Process("ReadDst");
         } else {
            chain.Process("ReadDst","",nentries);
         }
      } else {
         //use datasets string
         if( !nentries ) {
            proof->Process(dataset_names_string.c_str(), "ReadDst");
         } else {
            proof->Process(dataset_names_string.c_str(), "ReadDst","",nentries);
         }
      }

   } else {

      ReadDst* selector = new ReadDst;

      TList* inputList  = new TList;
      inputList->Add(bean);
      selector->SetInputList(inputList);

      if( !nentries ) {
         chain.Process(selector);
      } else {
         chain.Process(selector,"",nentries);
      }
      delete selector;

   }

   if (bean->IsProof()) {
      graceful_proof_exit();
   }

// ================== CLEAN UP AFTER TERMINATION ========================
#if defined (__unix__)
   if( copy_sqlite ) { // remove copy of database
     DatabaseSvc* dbs = DatabaseSvc::instance();
     string new_dir = dbs->GetDBFilePath();
     rm_whole_dir(new_dir.c_str());
   }
#endif

// ========================= DEBUG MEMORY LEAKS =========================
//    // do not forget put following in .rootrc file:
//    // Root.MemStat: 1
//    // Root.MemStat.size: -1
//    // Root.MemStat.cnt: -1
//    // Root.ObjectStat: 1
//    // display the contents of the memory table:
//    gObjectTable->Print();

   return 0;
}

//------------------------------------------------------------------------------
static void termination_handler(int isig)
//------------------------------------------------------------------------------
{
   cout << endl <<" User signal \""
#if defined (__unix__) || defined (__APPLE__)
        << strsignal(isig)
#elif defined _WIN32
        << isig
#endif
        << "\" had been received." << endl << flush;

   static int UserSignal = 0;
   switch( isig ) {
     case SIGINT: // "program interrupt" (the user types CTRL-C )
                  if( UserSignal == SIGINT ) {  // 2-d CTRL-C
                    cout << endl
                         <<" second CRTL-C had been detected. Abort!"
                         << endl << flush;
                    abort();
                  }

     case SIGTERM:// politely ask a program to terminate.
#if defined (__unix__) || defined (__APPLE__)
     case SIGHUP: // "hang up" (the user's terminal is disconnected)
#endif
                  cout << " Normal job termination." << endl;
                  UserSignal=isig;
                  break;

     // SIGSEGV is the signal sent to a process when it makes an invalid
     //         memory reference, or segmentation fault.
     case SIGSEGV:
                  cout << endl
                       <<" SIGSEGV signal had been received."
                       << endl << flush
                       << "Stack trace: " << endl << flush;

                  gSystem->StackTrace();

                  abort();

     default:
                  break;
   }

   if( UserSignal && bean ) bean->SetUserSignal(UserSignal);
}

#if defined _WIN32
//------------------------------------------------------------------------------
BOOL CtrlHandler( DWORD fdwCtrlType )
//------------------------------------------------------------------------------
{
   cout << endl
        <<" Ctrl signal \"" << fdwCtrlType << "\" had been received."
        << endl << flush;

   static int CtrlSignal = -1; // CTRL_C_EVENT is 0
   switch( fdwCtrlType ) {
     case CTRL_C_EVENT: // Handle the CTRL-C signal.
                  cout << "Ctrl-C event" << endl;
                  if( CtrlSignal == CTRL_C_EVENT ) {
                    cout << endl
                         <<" second CRTL-C had been detected. Abort!"
                         << endl << flush;
                    abort();
                  }
                  break;

     case CTRL_CLOSE_EVENT: // CTRL-CLOSE: confirm that the user wants to exit.
                  cout << "Ctrl-Close event" << endl;
                  break;

     default:
                  return FALSE;
   }

  CtrlSignal = fdwCtrlType;

  // this is user interupt
  if( (CtrlSignal != -1) && bean ) bean->SetUserSignal(SIGINT);

  return( TRUE );
}
#endif

//------------------------------------------------------------------------------
static void set_user_termination()
//------------------------------------------------------------------------------
{
#if defined (__unix__) || defined (__APPLE__)
   struct sigaction new_action, old_action;

   // Set up the structure to specify the new action
   new_action.sa_handler = &termination_handler;
   sigemptyset(&new_action.sa_mask);
   new_action.sa_flags = 0;

   // avoid handling signals previously set to be ignored
   sigaction(SIGINT, NULL, &old_action);
   if( old_action.sa_handler != SIG_IGN ) {
     sigaction(SIGINT, &new_action, NULL);
   }
   sigaction(SIGTERM, NULL, &old_action);
   if( old_action.sa_handler != SIG_IGN ) {
     sigaction(SIGTERM, &new_action, NULL);
   }
   sigaction(SIGHUP, NULL, &old_action);
   if( old_action.sa_handler != SIG_IGN ) {
     sigaction(SIGHUP, &new_action, NULL);
   }

   // ignoring SIGSEGV results in undefined behavior
   sigaction(SIGSEGV, &new_action, NULL);
#elif defined _WIN32
  signal(SIGINT,  termination_handler);
  signal(SIGTERM, termination_handler);
  signal(SIGSEGV, termination_handler);

  if( !SetConsoleCtrlHandler( (PHANDLER_ROUTINE) CtrlHandler, TRUE ) ) {
    cout << " set_user_termination ERROR: Could not set control handler"
         << endl;
    exit(1);
  }
  if( bean->Verbose() ) {
    cout << " set_user_termination: The Control Handler is installed" << endl;
  }
#endif
}

//------------------------------------------------------------------------------
static void segfault_handler(int isig)
//------------------------------------------------------------------------------
{
   // SIGSEGV is the signal sent to a process when it makes an invalid
   //         memory reference, or segmentation fault.

   cout << endl <<" SIGSEGV signal \""
#if defined (__unix__) || defined (__APPLE__)
        << strsignal(isig)
#elif defined _WIN32
        << isig
#endif
        << "\" had been received." << endl << flush;

   bean->Proof()->GetManager()->GetSessionLogs()->Save("*", "bean_proof.log");

   abort();
}

//------------------------------------------------------------------------------
static void set_proof_termination()
//------------------------------------------------------------------------------
{
#if defined (__unix__) || defined (__APPLE__)
   struct sigaction new_action;

   // Set up the structure to specify the new action
   new_action.sa_handler = &segfault_handler;
   sigemptyset(&new_action.sa_mask);
   new_action.sa_flags = 0;

   // ignoring SIGSEGV results in undefined behavior
   sigaction(SIGSEGV, &new_action, NULL);
#elif defined _WIN32
  signal(SIGSEGV, segfault_handler);
#endif
}

//------------------------------------------------------------------------------
static void check_upload_enable(TProof* proof, const char* pname)
//------------------------------------------------------------------------------
{
   if( proof->UploadPackage(pname) != 0 ) {
     cerr << "Cannot upload package "<< pname <<" on PROOF. Exiting."
          << endl;
     graceful_proof_exit();
     exit(1);
   }
   if( bean->Verbose() )
     cout << " Package " << pname << " successfully uploaded" << endl;

   Bool_t notOnClient = kTRUE; // to enable packages also on the client
   if( proof->EnablePackage(pname, notOnClient) != 0 ) {
      cerr << "Cannot enable package "<< pname <<" on PROOF. Exiting."
         << endl;
      graceful_proof_exit();
      exit(1);
   }
   if( bean->Verbose() )
     cout << " Package " << pname << " successfully enabled" << endl;
}

