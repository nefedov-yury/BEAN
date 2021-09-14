// print simple statistic about production:
// 1) initial number of events for data and MC
// 2) number of good MC events and 
//    number of Psi(2S)->pi+pi-J/Psi (N64)

//-------------------------------------------------------------------------
void num_print(string fname)
//-------------------------------------------------------------------------
{
  TFile* froot = TFile::Open(fname.c_str(),"READ");
  if( froot == 0 ) {
    cout << "can not open " << fname << endl;
    return;
  }

  froot->cd("PsipJpsiPhiEta");

  TH1D* cuts = (TH1D*)gDirectory->Get("All_cuts");
  if ( cuts ) {
    printf("%s: (0_cuts(1)=%.1f)\n",fname.c_str(),cuts->GetBinContent(1));
  } else {
    cout << " no 0_cuts ?" << endl;
    return;
  }

  TH1D* mc_dec0 = (TH1D*)gDirectory->Get("mc_dec0");
  if ( mc_dec0 ) {
    double Nall = mc_dec0->Integral();
    double N64 = mc_dec0->GetBinContent(65);
    printf(" MC: Nall=%.1f N64=%.1f Br=%.4f%%\n",
             Nall,N64,100*N64/Nall);
  }

}

//-------------------------------------------------------------------------
void num_prints()
//-------------------------------------------------------------------------
{
  gROOT->Reset();
//   string dir("prod-4/");
  string dir("prod-5/");

  vector<string> fnames = {
        "data_3650_all.root",
//         "data_09psip_all.root",
//         "data_12psip_all.root",
//         "mcinc_09psip_all.root",
//         "mcinc_12psip_all.root",
//         "mcsig_kkmc_09.root",
//         "mcsig_kkmc_12.root"
  };

  for ( auto fn : fnames ) num_print(dir+fn);

}
