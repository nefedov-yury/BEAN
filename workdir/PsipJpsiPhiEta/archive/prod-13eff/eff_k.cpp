{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Thu Apr  2 10:49:08 2020 by ROOT version6.18/04)
//   from TTree eff_K/K reconstruction efficiency
//   found on file: mcinc_12psip_all.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("mcinc_12psip_all.root");
   if (!f) {
      f = new TFile("mcinc_12psip_all.root");
   }
    TDirectory * dir = (TDirectory*)f->Get("mcinc_12psip_all.root:/PipPimKpKm");
    dir->GetObject("eff_K",tree);

//Declaration of leaves types
   Double_t        Zk;
   Double_t        Ptk;
   Double_t        Ck;
   Double_t        fl;
   Double_t        dP;
   Double_t        dTh;
   Double_t        Mrec;
   Double_t        Mrk2;
   Double_t        Egsum;
   Double_t        Egmax;
   Double_t        good;

   // Set branch addresses.
   eff_K->SetBranchAddress("Zk",&Zk);
   eff_K->SetBranchAddress("Ptk",&Ptk);
   eff_K->SetBranchAddress("Ck",&Ck);
   eff_K->SetBranchAddress("fl",&fl);
   eff_K->SetBranchAddress("dP",&dP);
   eff_K->SetBranchAddress("dTh",&dTh);
   eff_K->SetBranchAddress("Mrec",&Mrec);
   eff_K->SetBranchAddress("Mrk2",&Mrk2);
   eff_K->SetBranchAddress("Egsum",&Egsum);
   eff_K->SetBranchAddress("Egmax",&Egmax);
   eff_K->SetBranchAddress("good",&good);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// eff_K->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Long64_t nentries = eff_K->GetEntries();

   Long64_t nbytes = 0;
//   for (Long64_t i=0; i<nentries;i++) {
//      nbytes += eff_K->GetEntry(i);
//   }
}
