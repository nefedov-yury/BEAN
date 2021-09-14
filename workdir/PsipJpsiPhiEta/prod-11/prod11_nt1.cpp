{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Sat Sep 11 15:30:49 2021 by ROOT version6.24/02)
//   from TTree nt1/pi+ pi- J/Psi selection
//   found on file: mcsig_kkmc_09.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("mcsig_kkmc_09.root");
   if (!f) {
      f = new TFile("mcsig_kkmc_09.root");
   }
    TDirectory * dir = (TDirectory*)f->Get("mcsig_kkmc_09.root:/PsipJpsiPhiEta");
    dir->GetObject("nt1",tree);

//Declaration of leaves types
   Float_t         Mrs;
   Float_t         mdp;
   Float_t         MrsW;
   vector<float>   Mrec;
   Float_t         Mrb;
   Float_t         Ptp;
   Float_t         Cpls;
   Float_t         Ptm;
   Float_t         Cmns;
   UShort_t        nch;
   UShort_t        Nm;
   UShort_t        dec;
   Float_t         mcmkk;
   Float_t         mcmkpet;
   Float_t         mccosT;

   // Set branch addresses.
   nt1->SetBranchAddress("Mrs",&Mrs);
   nt1->SetBranchAddress("mdp",&mdp);
   nt1->SetBranchAddress("MrsW",&MrsW);
   nt1->SetBranchAddress("Mrec",&Mrec);
   nt1->SetBranchAddress("Mrb",&Mrb);
   nt1->SetBranchAddress("Ptp",&Ptp);
   nt1->SetBranchAddress("Cpls",&Cpls);
   nt1->SetBranchAddress("Ptm",&Ptm);
   nt1->SetBranchAddress("Cmns",&Cmns);
   nt1->SetBranchAddress("nch",&nch);
   nt1->SetBranchAddress("Nm",&Nm);
   nt1->SetBranchAddress("dec",&dec);
   nt1->SetBranchAddress("mcmkk",&mcmkk);
   nt1->SetBranchAddress("mcmkpet",&mcmkpet);
   nt1->SetBranchAddress("mccosT",&mccosT);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// nt1->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Long64_t nentries = nt1->GetEntries();

   Long64_t nbytes = 0;
//   for (Long64_t i=0; i<nentries;i++) {
//      nbytes += nt1->GetEntry(i);
//   }
}
