//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Bean                                                                 //
//                                                                      //
//  * Primary class to read DST                                         //
//  * Dynamic load and call of user functions                           //
//  * Save user histogram                                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>

// #include <typeinfo> // for debug

#include <TProof.h>
#include <TMacro.h>
#include <TObjString.h>
#include <TUrl.h>

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "Bean.h"

// These includes (espesially windows.h) MUST be after
// "RootEventData/TMcEvent.h" include
#if defined (__unix__) || defined (__APPLE__)
// __unix__ is defined by compiler for Unix and __APPLE__ for MacOS
#  include <dlfcn.h>    // interface to dynamic linking loader
#elif defined _WIN32    // _Win32 is defined for 32 or 64 bit Windows
#  include <windows.h>
#endif

// The type of "usrlib_handle" file variable depends from OS
#if defined (__unix__) || defined (__APPLE__)
static void* usrlib_handle = 0; //! this is not-persistent member
#elif defined _WIN32
static HINSTANCE usrlib_handle = 0; //! this is not-persistent member
#endif

ClassImp(Bean);

using namespace std;

//--------------------------------------------------------------------
Bean::Bean() : TObject()
//--------------------------------------------------------------------
{
   verbose = false;
   event_dump = false;
   hst_file = "bean_histo.root";
   max_number_events = 0;
   proof_xrd_output = 0;
   dst_file_is_local = 1;
   proof_lite = 0;
   dst_file_is_dataset = 0;
}

//--------------------------------------------------------------------
void Bean::SetDstFile(const char* opt)
//--------------------------------------------------------------------
{
   // parse output file path
   TString path(opt);
   dst_file_is_dataset = 0;
   dst_file_is_local = 0;

   if( path.BeginsWith("root://") ) {
         //~ dst_file_is_local = 0;
   } else {
      std::string dataset_name = ParseDatasetName(opt);
      if( dataset_name.size() ) {
         dst_file_is_dataset = 1;
         dst_file = dataset_name;
      } else {
         // In case of PROOF Lite the output file is *local*
         // but it is not *sandboxed*
         // So, DstFileIsLocal() will return false in such a case

         dst_file_is_local = 1;
      }
   }

   //path.Tokenize("/")->Last()->GetName()

   dst_file_name = path(path.Last('/')+1,path.Length()).Data();

   dst_file = opt;
}

//--------------------------------------------------------------------
void Bean::SetProof(TProof * val)
//--------------------------------------------------------------------
{
   proof = val;
   if (proof) {
       proof_lite = proof->IsLite();
   }
}

//--------------------------------------------------------------------
bool Bean::DstFileIsLocal() const
//--------------------------------------------------------------------
{
// Proof lite should be able to directly write file to proper
// destination on local PC So, the output file won't be 'local' i.e.
// it will be stored outside the sandbox

//~ local !lite  isLocal()
//~ 0      1     0
//~ 0      0     0
//~ 1      1     1
//~ 1      0     0
    // TODO: set PROOF first

   return dst_file_is_local && (!proof_lite);
}

//--------------------------------------------------------------------
void Bean::PrintOptions() const
//--------------------------------------------------------------------
{
   if( verbose ) {
     cout << " Bean::PrintOptions() " << endl;
     cout << " event_dump= " << ((event_dump) ? "true" : "false") << endl;
     if( max_number_events ) {
       cout << " process not more than " << max_number_events
            << " events in total " << endl;
     }
     cout << " hst_file= " << hst_file << endl;
     cout << " dst_file= " << dst_file << endl;
     if (DstFileIsLocal() ) cout << "   dst_file is local "<< endl;
     cout << " dst_file_name= " << dst_file_name << endl;

     if( IsProof() ) {
       cout << " Use PROOF cluster: " << proof_clr << endl;
     }

     cout << " base dir is: " << base_dir << endl;

     cout << " User functions: " << endl;
     VecName::const_iterator it = Ufn_names.begin();
     for( ;it != Ufn_names.end(); it++) {
        cout << "     " << *it << endl;
     }
   }
}

//--------------------------------------------------------------------
void Bean::AddUserFcn(const char* name)
//--------------------------------------------------------------------
{
   Ufn_names.push_back(string(name));
}

//--------------------------------------------------------------------
void Bean::LoadUserLib()
//--------------------------------------------------------------------
{
   if( verbose ) cout << " Bean::LoadUserLib()" << endl;

#if defined (__unix__) || defined (__APPLE__)
   // string usrlib_name = "libUser.so";
   // if( IsProof() ) usrlib_name = "BeanUser/" + usrlib_name;
   // usrlib_handle = dlopen(usrlib_name.c_str(), RTLD_NOW | RTLD_LOCAL);
   // RTLD_NOW - all undefined symbols in the library are resolved
   //            before dlopen()  returns;
   // RTLD_LOCAL - symbols defined in this library are not made available
   //            to resolve references in subsequently loaded libraries.

   usrlib_handle = dlopen(0, RTLD_NOW | RTLD_GLOBAL);
#elif defined _WIN32
   if( ! SetDllDirectory("lib") ) {
     cout << " Bean::LoadUserLib() ==> cannot SetDllDirectory(): "
          << GetLastError() << endl;
     exit(1);
   }
   usrlib_handle = LoadLibrary("BeanUser");
#endif

   if( !usrlib_handle ) {
     cout << " Bean::LoadUserLib() ==> cannot ";
#if defined (__unix__) || defined (__APPLE__)
     // cout << "dlopen(" << usrlib_name << ")"
     cout << "dlopen(0)" << endl;
     cout << dlerror() << endl;
     exit(1);
#elif defined _WIN32
     cout << "LoadLibrary(\"lib\\BeanUser.dll\")" << endl;
     cout << GetLastError() << endl;
#ifdef _DEBUG   // defined automatically in debug conf. of Visual Studio
     cout << "Try to LoadLibrary(\"lib/Debug/BeanUser.dll\")";
     if( SetDllDirectory("lib/Debug") )
       usrlib_handle = LoadLibrary("BeanUser");
     if( usrlib_handle ) {
       cout << " ==> success" << endl;
     } else {
       cout << " ==> failure" << GetLastError() << endl;
       exit(1);
     }
#else
     exit(1);
#endif  // _DEBUG
#endif  // _WIN32
   }

   if( verbose ) {
      cout << " User library is loaded" << endl;
   }
}

//--------------------------------------------------------------------
void Bean::LoadUserFcn(const char* name)
//--------------------------------------------------------------------
{
   if( !usrlib_handle ) LoadUserLib();

   string sname(name);
   string nameJobStart = sname + "StartJob";
   string nameEvent = sname + "Event";
   string nameJobEnd = sname + "EndJob";

   if( verbose ) {
     cout << " Load user function(s): ";
   }

#if defined (__unix__) || defined (__APPLE__)
   void* sym = dlsym(usrlib_handle, nameEvent.c_str());
#elif defined _WIN32
   FARPROC sym = GetProcAddress(usrlib_handle, nameEvent.c_str());
#endif
   if( sym ) {
     Ufn_event.push_back(sym);
     if( verbose ) cout << nameEvent << " ";
   } else {
     if( verbose ) cout << endl;
     cout << " User function " << nameEvent << " does not exist!" << endl;
     cout << " This is mandatory function. Stop!" << endl;
     exit(1);
   }

#if defined (__unix__) || defined (__APPLE__)
   sym = dlsym(usrlib_handle, nameJobStart.c_str());
#elif defined _WIN32
   sym = GetProcAddress(usrlib_handle, nameJobStart.c_str());
#endif
   if( sym ) {
     Ufn_start.push_back(sym);
     if( verbose ) cout << nameJobStart << " ";
   }

#if defined (__unix__) || defined (__APPLE__)
   sym = dlsym(usrlib_handle, nameJobEnd.c_str());
#elif defined _WIN32
   sym = GetProcAddress(usrlib_handle, nameJobEnd.c_str());
#endif
   if( sym ) {
     Ufn_end.push_back(sym);
     if( verbose ) cout << nameJobEnd << " ";
   }

   if( verbose ) cout << endl;
}

//--------------------------------------------------------------------
void Bean::LoadUserFcns()
//--------------------------------------------------------------------
{
   // load all user functions
   if( verbose ) cout << " Bean::LoadUserFcns() " << endl;

   Ufn_event.clear();
   Ufn_start.clear();
   Ufn_end.clear();

   VecName::const_iterator ufn = Ufn_names.begin();
   for( ;ufn != Ufn_names.end(); ufn++) {
      LoadUserFcn( ufn->c_str() );
   }
}

//--------------------------------------------------------------------
std::string  Bean::ParseDatasetName(const char * filename )
//--------------------------------------------------------------------
{
   TString name(filename);
   TString prefix("ds://");

   if (name.BeginsWith(prefix)) {
      int len = prefix.Length();
      return  name(len, name.Length() - len + 1).Data();
   } else {
      return "";
   }

}

//--------------------------------------------------------------------
std::string Bean::GetProofMasterWorkdir()
//--------------------------------------------------------------------
{
   if (proof) {
       RedirectHandle_t rh;
       TString tempFileName;
       gSystem->TempFileName(tempFileName);
       gSystem->RedirectOutput(tempFileName, "a", &rh);
       proof->Print();
       gSystem->RedirectOutput(0, 0, &rh);

       TMacro macro;
       macro.ReadFile(tempFileName);
       gSystem->Unlink(tempFileName);


       TObject* line = macro.GetLineWith("Working directory:");
       if (line) {
          TString s = line->GetName();
          TString ss(s(s.First(":")+1, s.Length()));
          TString sss(ss.Strip(TString::kBoth));
          return sss.Data();
       } else {
          return "";
       }

   } else {
       return "";
   }
}
