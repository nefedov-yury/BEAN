#ifndef _Bean_h
#define _Bean_h

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Bean                                                                 //
//                                                                      //
//  * Primary class to read DST                                         //
//  * Dynamic load and call of user functions                           //
//  * Save user histogram                                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <TObject.h>

class TProof;

typedef std::vector<void* > VecUF;
typedef std::vector<std::string> VecName;

class Bean : public TObject
{
public :
                Bean();

   virtual     ~Bean(){};

   // -- set/get options functions:
   void         SetVerbose()                    {verbose = true;}
   void         SetSilent()                     {verbose = false;}
   void         SetUserSignal(int sig)          {user_signal = sig;}
   void         SetEventDump()                  {event_dump = true;}
   void         SetHstFile(const char* opt)     {hst_file = opt;}
   void         SetDstFile(const char* opt);
   void         SetProofClr(const char* opt)    {proof_clr= (opt)? opt : " ";}
   void         SetProofParam(const char* opt)  {proof_param= opt;}
   void         SetMaxNumberEvents(int val)     {max_number_events = val;}

   void         SetProof(TProof * val);

   void         SetProofXrdOutput(bool val)     {proof_xrd_output = val;}

   bool         Verbose() const                 {return verbose;}
   int          UserSignal() const              {return user_signal;}
   bool         IsDump() const                  {return event_dump;}
   std::string  HstFile() const                 {return hst_file;}
   std::string  DstFile() const                 {return dst_file;}
   std::string  DstFileName() const             {return dst_file_name;}

// 'is local' means 'is in local sandbox'
// i.e. output should be stored in current directory.
// may be rename it to DstFileIsSandbox ?
   bool         DstFileIsLocal() const;
   bool         DstFileIsDataset() const        {return dst_file_is_dataset;}


   std::string  ProofClr() const                {return proof_clr;}
   std::string  ProofParam() const              {return proof_param;}
   bool         IsProof() const                 {return !proof_clr.empty();}
   TProof*      Proof() const                   {return proof;}
   bool         IsProofLite() const             {return proof_lite;}

   long         MaxNumberEvents() const         {return max_number_events;}
   void         PrintOptions() const;
   bool         ProofXrdOutput() const          {return proof_xrd_output;}

   // -- user functions
   void         AddUserFcn(const char* name);
   unsigned int NUserFns() const                {return Ufn_names.size();}
   const VecUF& GetStartJobFns() const          {return Ufn_start;}
   const VecUF& GetUserEventFns() const         {return Ufn_event;}
   const VecUF& GetEndJobFns() const            {return Ufn_end;}

   void         LoadUserLib();
   void         LoadUserFcn(const char* name);
   void         LoadUserFcns();

   std::string  GetProofMasterWorkdir();
   std::string  ParseDatasetName(const char * name);

   void         SetBaseDir(std::string dir)     {base_dir = dir;}
   std::string  GetBaseDir() const              {return base_dir;}

private :
   bool         verbose;         // by default false => silent mode
   int          user_signal;
   bool         event_dump;
   std::string  hst_file;
   std::string  dst_file; // full path
   std::string  proof_clr;
   std::string  proof_param;
   bool         proof_lite;
   long         max_number_events;

   TProof*      proof; //! this is not-persistent member

   VecUF        Ufn_start; //! vector of user functions StartJob
   VecUF        Ufn_event; //! vector of user Event functions
   VecUF        Ufn_end;   //! vector of user functions EndJob
   VecName      Ufn_names; // names of loaded functions

   int          proof_xrd_output;

   std::string  dst_file_name;
   bool         dst_file_is_local;
   bool         dst_file_is_dataset;

   std::string  base_dir; // directory with Analysis, BeanUser, etc.

   // ClassVersionID=1 because we have need object I/O
   ClassDef(Bean,1); // Primary class to read DST
};

#endif
