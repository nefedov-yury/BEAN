// print simple statistic about production:
// root -q -b 'num_prod.cc("psip_bean_5")' > prod5.test

//-------------------------------------------------------------------------
void num_print(string dir, string fname)
//-------------------------------------------------------------------------
{
  string f = dir + string("/") + fname;
  TFile* froot = TFile::Open(f.c_str(),"READ");
  if( froot == 0 ) {
    cout << "can not open " << f << endl;
    return;
  }

  froot->cd("PsipJpsiPhiEta");

  TH1D* cuts = (TH1D*)gDirectory->Get("All_cuts");
  if ( !cuts ) {
    // try old name
    cuts = (TH1D*)gDirectory->Get("0_cuts");
  }
  if ( cuts ) {
    printf("%s: (0_cuts(1)=%.1f)\n",fname.c_str(),cuts->GetBinContent(1));
  } else {
    cout << fname << ": no 0_cuts ?" << endl;
    return;
  }

}

//-------------------------------------------------------------------------
void num_prod(string dir)
//-------------------------------------------------------------------------
{
  gROOT->Reset();
//   string dir("psip_bean_5/");

  vector<string> fnames = {
        "data_3650_",
        "data_09psip_",
        "data_12psip_",
        "mcinc_09psip_",
        "mcinc_12psip_",
  };

  vector<int> Nf = {
        12,
        60,
        92,
        16,
        57
  };

  for ( int iset = 0; iset < 5; iset++ ) {
     for ( int i = 0; i < Nf[iset]; ++i ) {
        string i_str = to_string(i);
        if ( i < 10 ) i_str = string("0")+i_str;
        num_print(dir,fnames[iset]+i_str+string(".root"));
     }
  }

}
