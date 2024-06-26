{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Mon Oct 31 16:43:00 2022 by ROOT version6.26/06)
//   from TTree a4c/after 4C kinematic fit
//   found on file: ntpl_3080_rs.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ntpl_3080_rs.root");
   if (!f) {
      f = new TFile("ntpl_3080_rs.root");
   }
    TDirectory * dir = (TDirectory*)f->Get("ntpl_3080_rs.root:/SelectKKgg");
    dir->GetObject("a4c",tree);

//Declaration of leaves types
   Double_t        ch2;
   Double_t        chsq3g;
   Double_t        Pkp;
   Double_t        Ckp;
   Double_t        phikp;
   Double_t        Pkm;
   Double_t        Ckm;
   Double_t        phikm;
   Double_t        Eg1;
   Double_t        Cg1;
   Double_t        phig1;
   Double_t        Eg2;
   Double_t        Cg2;
   Double_t        phig2;
   Double_t        Pgg;
   Double_t        Cgg;
   Double_t        phigg;
   Double_t        Mkk;
   Double_t        Mgg;
   Double_t        M2kpeta;
   Double_t        M2kmeta;
   Double_t        dec;
   Double_t        xisr;
   Double_t        mcmkk;
   Double_t        mcmkpet;

   // Set branch addresses.
   a4c->SetBranchAddress("ch2",&ch2);
   a4c->SetBranchAddress("chsq3g",&chsq3g);
   a4c->SetBranchAddress("Pkp",&Pkp);
   a4c->SetBranchAddress("Ckp",&Ckp);
   a4c->SetBranchAddress("phikp",&phikp);
   a4c->SetBranchAddress("Pkm",&Pkm);
   a4c->SetBranchAddress("Ckm",&Ckm);
   a4c->SetBranchAddress("phikm",&phikm);
   a4c->SetBranchAddress("Eg1",&Eg1);
   a4c->SetBranchAddress("Cg1",&Cg1);
   a4c->SetBranchAddress("phig1",&phig1);
   a4c->SetBranchAddress("Eg2",&Eg2);
   a4c->SetBranchAddress("Cg2",&Cg2);
   a4c->SetBranchAddress("phig2",&phig2);
   a4c->SetBranchAddress("Pgg",&Pgg);
   a4c->SetBranchAddress("Cgg",&Cgg);
   a4c->SetBranchAddress("phigg",&phigg);
   a4c->SetBranchAddress("Mkk",&Mkk);
   a4c->SetBranchAddress("Mgg",&Mgg);
   a4c->SetBranchAddress("M2kpeta",&M2kpeta);
   a4c->SetBranchAddress("M2kmeta",&M2kmeta);
   a4c->SetBranchAddress("dec",&dec);
   a4c->SetBranchAddress("xisr",&xisr);
   a4c->SetBranchAddress("mcmkk",&mcmkk);
   a4c->SetBranchAddress("mcmkpet",&mcmkpet);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// a4c->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Long64_t nentries = a4c->GetEntries();

   Long64_t nbytes = 0;
//   for (Long64_t i=0; i<nentries;i++) {
//      nbytes += a4c->GetEntry(i);
//   }
}
