{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Thu Dec 10 22:03:56 2020 by ROOT version6.22/02)
//   from TTree a4c/after 4C kinematic fit
//   found on file: data_12psip_all.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data_12psip_all.root");
   if (!f) {
      f = new TFile("data_12psip_all.root");
   }
    TDirectory * dir = (TDirectory*)f->Get("data_12psip_all.root:/PsipJpsiPhiEta");
    dir->GetObject("a4c",tree);

//Declaration of leaves types
   Double_t        Mrec;
   Double_t        ch2;
   Double_t        Ppip;
   Double_t        Cpip;
   Double_t        phipip;
   Double_t        Ppim;
   Double_t        Cpim;
   Double_t        phipim;
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
   Double_t        Pkk;
   Double_t        Ckk;
   Double_t        phikk;
   Double_t        Pgg;
   Double_t        Cgg;
   Double_t        phigg;
   Double_t        Pj;
   Double_t        Cj;
   Double_t        phij;
   Double_t        Mkk;
   Double_t        Mgg;
   Double_t        M2kpeta;
   Double_t        M2kmeta;
   Double_t        Ptpip;
   Double_t        Ptpim;
   Double_t        Ptkp;
   Double_t        Ptkm;
   Double_t        Ptgg;
   Double_t        dec;
   Double_t        decj;
   Double_t        mcmkk;

   // Set branch addresses.
   a4c->SetBranchAddress("Mrec",&Mrec);
   a4c->SetBranchAddress("ch2",&ch2);
   a4c->SetBranchAddress("Ppip",&Ppip);
   a4c->SetBranchAddress("Cpip",&Cpip);
   a4c->SetBranchAddress("phipip",&phipip);
   a4c->SetBranchAddress("Ppim",&Ppim);
   a4c->SetBranchAddress("Cpim",&Cpim);
   a4c->SetBranchAddress("phipim",&phipim);
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
   a4c->SetBranchAddress("Pkk",&Pkk);
   a4c->SetBranchAddress("Ckk",&Ckk);
   a4c->SetBranchAddress("phikk",&phikk);
   a4c->SetBranchAddress("Pgg",&Pgg);
   a4c->SetBranchAddress("Cgg",&Cgg);
   a4c->SetBranchAddress("phigg",&phigg);
   a4c->SetBranchAddress("Pj",&Pj);
   a4c->SetBranchAddress("Cj",&Cj);
   a4c->SetBranchAddress("phij",&phij);
   a4c->SetBranchAddress("Mkk",&Mkk);
   a4c->SetBranchAddress("Mgg",&Mgg);
   a4c->SetBranchAddress("M2kpeta",&M2kpeta);
   a4c->SetBranchAddress("M2kmeta",&M2kmeta);
   a4c->SetBranchAddress("Ptpip",&Ptpip);
   a4c->SetBranchAddress("Ptpim",&Ptpim);
   a4c->SetBranchAddress("Ptkp",&Ptkp);
   a4c->SetBranchAddress("Ptkm",&Ptkm);
   a4c->SetBranchAddress("Ptgg",&Ptgg);
   a4c->SetBranchAddress("dec",&dec);
   a4c->SetBranchAddress("decj",&decj);
   a4c->SetBranchAddress("mcmkk",&mcmkk);

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
