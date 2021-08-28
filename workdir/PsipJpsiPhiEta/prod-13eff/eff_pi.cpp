{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Thu Apr  2 10:50:08 2020 by ROOT version6.18/04)
//   from TTree eff_Pi/Pi reconstruction efficiency
//   found on file: mcinc_12psip_all.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("mcinc_12psip_all.root");
   if (!f) {
      f = new TFile("mcinc_12psip_all.root");
   }
    TDirectory * dir = (TDirectory*)f->Get("mcinc_12psip_all.root:/PipPimKpKm");
    dir->GetObject("eff_Pi",tree);

//Declaration of leaves types
   Double_t        Zpi;
   Double_t        Ptpi;
   Double_t        Cpi;
   Double_t        fl;
   Double_t        dP;
   Double_t        dTh;
   Double_t        MppKK;
   Double_t        Mrpi2;
   Double_t        Egsum;
   Double_t        Egmax;
   Double_t        good;

   // Set branch addresses.
   eff_Pi->SetBranchAddress("Zpi",&Zpi);
   eff_Pi->SetBranchAddress("Ptpi",&Ptpi);
   eff_Pi->SetBranchAddress("Cpi",&Cpi);
   eff_Pi->SetBranchAddress("fl",&fl);
   eff_Pi->SetBranchAddress("dP",&dP);
   eff_Pi->SetBranchAddress("dTh",&dTh);
   eff_Pi->SetBranchAddress("MppKK",&MppKK);
   eff_Pi->SetBranchAddress("Mrpi2",&Mrpi2);
   eff_Pi->SetBranchAddress("Egsum",&Egsum);
   eff_Pi->SetBranchAddress("Egmax",&Egmax);
   eff_Pi->SetBranchAddress("good",&good);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// eff_Pi->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Long64_t nentries = eff_Pi->GetEntries();

   Long64_t nbytes = 0;
//   for (Long64_t i=0; i<nentries;i++) {
//      nbytes += eff_Pi->GetEntry(i);
//   }
}
