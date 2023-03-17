// compare cross-sections for varying of chi^2 cut:  80 +/- 20
//   -> cmp_cs_chi2.pdf

/*
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -> chi2=60:
  // Efficiency from MCGPJ:
  double Xisr_min=0.9;
  double eff[17] = {
           0.12240,  0.12305,  0.12411,  0.12504,  0.12753,
           0.12895,  0.13452,  0.13713,  0.13578,  0.13663,
           0.13521,  0.13543,  0.13568,  0.13224,  0.12772,
           0.11060,
           0.12411 // this is 3080_new
                   };

  // Sig, ErrSig - cross-sections in picobarn
  double Xisr_min=0.9;
  Xmax = 1. - Xisr_min;
  vector<double> Sig =      { // pb
        2.35174e+01, 2.75883e+01, 3.31572e+01, 1.71469e+01, 3.09205e+01,
        4.04181e+01, 2.27003e+02, 7.30858e+02, 2.25034e+03, 3.55194e+03,
        1.73695e+03, 4.59391e+02, 1.12200e+02, 6.60865e+01, 9.30902e+01,
        4.38844e+01, 2.54581e+01,
                            };
  vector<double> ErrSig =   { // pb
        5.13192e+00, 5.51766e+00, 5.60459e+00, 9.07327e+00, 6.00986e+00,
        6.55668e+00, 4.01289e+01, 7.91925e+01, 1.29116e+02, 1.63849e+02,
        1.09632e+02, 9.57896e+01, 3.73999e+01, 2.43538e+01, 2.94377e+01,
        2.53367e+01, 1.84169e+00,
                            };

  Efficiency from MCGPJ: s(after ISR)/s(ini) in [0.9,1.0]
 =========================================================

  E(GeV)    Lum.(pb-1)     Signal      Eff.(%)   Cross Section(nb)
 3.049663     14.919     21 +/-  4.6    12.24    0.0235 +/- 0.0051
 3.058707      15.06     25 +/-  5.0    12.31    0.0276 +/- 0.0055
 3.079645     17.393     35 +/-  5.9    12.41    0.0332 +/- 0.0056
 3.082510      4.769      5 +/-  2.6    12.50    0.0171 +/- 0.0091
 3.088868     15.558     30 +/-  5.8    12.75    0.0309 +/- 0.0060
 3.091774      14.91     38 +/-  6.2    12.90    0.0404 +/- 0.0066
 3.094711      2.143     32 +/-  5.7    13.45    0.2270 +/- 0.0401
 3.095444      1.816     89 +/-  9.6    13.71    0.7309 +/- 0.0792
 3.095840      2.135    319 +/- 18.3    13.58    2.2503 +/- 0.1291
 3.097227      2.069    491 +/- 22.6    13.66    3.5519 +/- 0.1638
 3.098354      2.203    253 +/- 16.0    13.52    1.7370 +/- 0.1096
 3.099056      0.756     23 +/-  4.8    13.54    0.4594 +/- 0.0958
 3.101373      1.612     12 +/-  4.0    13.57    0.1122 +/- 0.0374
 3.105594      2.106      9 +/-  3.3    13.22    0.0661 +/- 0.0244
 3.112065       1.72     10 +/-  3.2    12.77    0.0931 +/- 0.0294
 3.119892      1.264      3 +/-  1.7    11.06    0.0439 +/- 0.0253
 3.080000     126.21    195 +/- 14.1    12.41    0.0255 +/- 0.0018

 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
/*
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -> chi2=100:
  // Efficiency from MCGPJ:
  double Xisr_min=0.9;
  double eff[17] = {
           0.12609,  0.12730,  0.12790,  0.12883,  0.13080,
           0.13238,  0.13752,  0.13971,  0.13831,  0.13943,
           0.13808,  0.13855,  0.13866,  0.13632,  0.13411,
           0.12357,
           0.12790 // this is 3080_new
                   };

  // Sig, ErrSig - cross-sections in picobarn
  double Xisr_min=0.9;
  Xmax = 1. - Xisr_min;
  vector<double> Sig =      { // pb
        2.50034e+01, 2.88006e+01, 3.21747e+01, 1.99709e+01, 3.01475e+01,
        4.14430e+01, 2.15112e+02, 7.57662e+02, 2.24380e+03, 3.60112e+03,
        1.71430e+03, 4.49046e+02, 1.09788e+02, 6.41085e+01, 8.86547e+01,
        3.92783e+01, 2.50837e+01,
                            };
  vector<double> ErrSig =   { // pb
        5.21356e+00, 5.54268e+00, 5.59174e+00, 9.41438e+00, 5.85962e+00,
        6.71453e+00, 3.98621e+01, 7.97923e+01, 1.28446e+02, 1.63505e+02,
        1.09024e+02, 9.36325e+01, 3.65961e+01, 2.36249e+01, 2.80351e+01,
        2.26773e+01, 1.81828e+00,
                            };

  Efficiency from MCGPJ: s(after ISR)/s(ini) in [0.9,1.0]
 =========================================================

  E(GeV)    Lum.(pb-1)     Signal      Eff.(%)   Cross Section(nb)
 3.049663     14.919     23 +/-  4.8    12.61    0.0250 +/- 0.0052
 3.058707      15.06     27 +/-  5.2    12.73    0.0288 +/- 0.0055
 3.079645     17.393     35 +/-  6.1    12.79    0.0322 +/- 0.0056
 3.082510      4.769      6 +/-  2.8    12.88    0.0200 +/- 0.0094
 3.088868     15.558     30 +/-  5.8    13.08    0.0301 +/- 0.0059
 3.091774      14.91     40 +/-  6.5    13.24    0.0414 +/- 0.0067
 3.094711      2.143     31 +/-  5.7    13.75    0.2151 +/- 0.0399
 3.095444      1.816     94 +/-  9.9    13.97    0.7577 +/- 0.0798
 3.095840      2.135    324 +/- 18.5    13.83    2.2438 +/- 0.1284
 3.097227      2.069    508 +/- 23.1    13.94    3.6011 +/- 0.1635
 3.098354      2.203    255 +/- 16.2    13.81    1.7143 +/- 0.1090
 3.099056      0.756     23 +/-  4.8    13.86    0.4490 +/- 0.0936
 3.101373      1.612     12 +/-  4.0    13.87    0.1098 +/- 0.0366
 3.105594      2.106      9 +/-  3.3    13.63    0.0641 +/- 0.0236
 3.112065       1.72     10 +/-  3.2    13.41    0.0887 +/- 0.0280
 3.119892      1.264      3 +/-  1.7    12.36    0.0393 +/- 0.0227
 3.080000     126.21    198 +/- 14.4    12.79    0.0251 +/- 0.0018

 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
//-------------------------------------------------------------------------
void cmpr_cs_chi2()
//-------------------------------------------------------------------------
{
  gROOT->Reset();

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLegendTextSize(0.05);

  int Nbeam = 17;
  double ebeam[] = {
                 3050.213, 3059.257, 3080.195, 3083.060, 3089.418,
                 3092.324, 3095.261, 3095.994, 3096.390, 3097.777,
                 3098.904, 3099.606, 3101.923, 3106.144, 3112.615,
                 3120.442,
                 3080.
                   };
  for(int i=0; i < Nbeam-1; i++) {// all except for 3080_new
     ebeam[i] -=0.55;
  }

  // draw errors for energy
  double ErrE[17];
  for(int i = 0; i < Nbeam; i++) ErrE[i] = 0.;

// chi2_cut = 80 +++ standard
  vector<double> Sig =      { // pb
        2.52517e+01, 2.80780e+01, 3.25385e+01, 1.68147e+01, 3.04337e+01,
        4.08263e+01, 2.23957e+02, 7.54028e+02, 2.23143e+03, 3.58859e+03,
        1.72718e+03, 4.52804e+02, 1.10674e+02, 6.47354e+01, 8.99491e+01,
        4.05821e+01, 2.49830e+01,
                            };
  vector<double> ErrSig =   { // pb
        5.26534e+00, 5.50655e+00, 5.50001e+00, 8.89751e+00, 5.91524e+00,
        6.53745e+00, 3.95903e+01, 7.98528e+01, 1.28580e+02, 1.63469e+02,
        1.09425e+02, 9.44162e+01, 3.68914e+01, 2.38559e+01, 2.84444e+01,
        2.34301e+01, 1.81638e+00,
                            };
// chi2_cut = 60 +
  vector<double> SigM =      { // pb
        2.35174e+01, 2.75883e+01, 3.31572e+01, 1.71469e+01, 3.09205e+01,
        4.04181e+01, 2.27003e+02, 7.30858e+02, 2.25034e+03, 3.55194e+03,
        1.73695e+03, 4.59391e+02, 1.12200e+02, 6.60865e+01, 9.30902e+01,
        4.38844e+01, 2.54581e+01,
                            };
  vector<double> ErrSigM =   { // pb
        5.13192e+00, 5.51766e+00, 5.60459e+00, 9.07327e+00, 6.00986e+00,
        6.55668e+00, 4.01289e+01, 7.91925e+01, 1.29116e+02, 1.63849e+02,
        1.09632e+02, 9.57896e+01, 3.73999e+01, 2.43538e+01, 2.94377e+01,
        2.53367e+01, 1.84169e+00,
                            };
// chi2_cut = 100 +
  vector<double> SigP =      { // pb
        2.50034e+01, 2.88006e+01, 3.21747e+01, 1.99709e+01, 3.01475e+01,
        4.14430e+01, 2.15112e+02, 7.57662e+02, 2.24380e+03, 3.60112e+03,
        1.71430e+03, 4.49046e+02, 1.09788e+02, 6.41085e+01, 8.86547e+01,
        3.92783e+01, 2.50837e+01,
                            };
  vector<double> ErrSigP =   { // pb
        5.21356e+00, 5.54268e+00, 5.59174e+00, 9.41438e+00, 5.85962e+00,
        6.71453e+00, 3.98621e+01, 7.97923e+01, 1.28446e+02, 1.63505e+02,
        1.09024e+02, 9.36325e+01, 3.65961e+01, 2.36249e+01, 2.80351e+01,
        2.26773e+01, 1.81828e+00,
                            };

  vector<double> RM(Nbeam), RP(Nbeam);
  string serr("  double err_chi2[] = {");
  for(unsigned int i = 0; i < Nbeam; i++) {
//     RM[i] = 100*fabs(1-SigM[i]/Sig[i]);
//     RP[i] = 100*fabs(1-SigP[i]/Sig[i]);
    RM[i] = 100*(1-SigM[i]/Sig[i]);
    RP[i] = 100*(1-SigP[i]/Sig[i]);
    cout << ebeam[i]
         << " ratio(-:+)= " << RM[i] << " : " << RP[i] << endl;
    serr += string(Form(" %.2f,",max(fabs(RM[i]),fabs(RP[i]))));
    if( (i+1)%5 == 0 ) serr += string(1,'\n')+string(23,' ');
  }
  serr[serr.size()-1]=' ';
  serr += "};";
  cout << serr << endl;

  TCanvas* c1 = new TCanvas("c1","...",0,0,800,600);

  string pdf( "cmp_cs_chi2.pdf" ); // <- PDF
  c1->Print((pdf+"[").c_str()); // just open pdf-file
  int DEBUG = 1;

  if( DEBUG==1 ) {
    c1->cd();
    gPad->SetGrid();
    gPad->SetLogy();

    TGraphErrors* grW1 = new TGraphErrors(
        Nbeam, ebeam, Sig.data(), ErrE, ErrSig.data());

    string titW1 = string(";center of mass energy (MeV)")+
                   string(";cross-section (pb)");
    grW1->SetTitle(titW1.c_str());
    grW1->SetMaximum(10000.);
    grW1->SetMinimum(5.);
    grW1->SetMarkerColor(kRed+2);
    grW1->SetMarkerStyle(20);
    grW1->SetMarkerSize(0.9);
    grW1->GetXaxis()->CenterTitle();
    grW1->GetYaxis()->CenterTitle();
    grW1->GetXaxis()->SetTitleSize(0.04);
    grW1->GetYaxis()->SetTitleSize(0.04);
    grW1->GetXaxis()->SetTitleOffset(1.1);
//     grW1->GetYaxis()->SetTitleOffset(1.0);
    grW1->GetXaxis()->SetLabelFont(62);
    grW1->GetYaxis()->SetLabelFont(62);
    grW1->GetXaxis()->SetLabelSize(0.03);
    grW1->GetYaxis()->SetLabelSize(0.03);

    grW1->Draw("AP");

    TGraphErrors* grW2 = new TGraphErrors(
        Nbeam, ebeam, SigM.data(), ErrE, ErrSigM.data());
    grW2->SetMarkerColor(kBlue+2);
    grW2->SetMarkerStyle(20);
    grW2->SetMarkerSize(1.1);
    grW2->Draw("P");

    TGraphErrors* grW3 = new TGraphErrors(
        Nbeam, ebeam, SigP.data(), ErrE, ErrSigP.data());
    grW3->SetMarkerColor(kGreen+2);
    grW3->SetMarkerStyle(20);
    grW3->SetMarkerSize(1.1);
    grW3->Draw("P");

    grW1->Draw("P"); // on top of pict

    TLegend* legD = new TLegend(0.12,0.52,0.50,0.88);
    legD->SetHeader("variation of #chi^{2} cut:","C");
    legD->AddEntry( grW1, "#sigma for #chi^{2}<80", "EP");
    legD->AddEntry( grW2, "#sigma for #chi^{2}<60", "EP");
    legD->AddEntry( grW3, "#sigma for #chi^{2}<100", "EP");
    legD->Draw();

    c1->Print(pdf.c_str()); // add to pdf-file
  } // end of if(DEBUG)

  c1->cd();
  gPad->SetGrid();
  gPad->SetLogy(false);

  TGraph* rat1 = new TGraph(Nbeam, ebeam,RM.data());
  TGraph* rat2 = new TGraph(Nbeam, ebeam,RP.data());
  string tit1 = string(";center of mass energy (MeV)")+
                string(";relative difference (%)");
  rat1->SetTitle(tit1.c_str());
  rat1->SetMaximum(20.);
  rat1->SetMinimum(-10.);
  rat1->SetMarkerColor(kBlue+2);
  rat1->SetMarkerStyle(20);
  rat1->SetMarkerSize(1.1);
  rat1->GetXaxis()->CenterTitle();
  rat1->GetYaxis()->CenterTitle();
  rat1->GetXaxis()->SetTitleSize(0.04);
  rat1->GetYaxis()->SetTitleSize(0.04);
  rat1->GetXaxis()->SetTitleOffset(1.1);
  rat1->GetYaxis()->SetTitleOffset(1.0);
  rat1->GetXaxis()->SetLabelFont(62);
  rat1->GetYaxis()->SetLabelFont(62);
  rat1->GetXaxis()->SetLabelSize(0.03);
  rat1->GetYaxis()->SetLabelSize(0.03);

  rat1->Draw("AP");

  rat2->SetMarkerColor(kGreen+2);
  rat2->SetMarkerStyle(20);
  rat2->SetMarkerSize(1.1);

  rat2->Draw("P");

  TLegend* leg = new TLegend(0.54,0.51,0.89,0.89);
  leg->SetHeader("variation of #chi^{2} cut:","C");
  leg->AddEntry(rat1,
      "#1 - #frac{#sigma(#chi^{2}<60)}{#sigma(#chi^{2}<80)}","P");
//       " #||{1 - #frac{#sigma(#chi^{2}<60)}{#sigma(#chi^{2}<80)}}","P");
  leg->AddEntry(rat2,
      "#1 - #frac{#sigma(#chi^{2}<100)}{#sigma(#chi^{2}<80)}","P");
//       "#||{1 - #frac{#sigma(#chi^{2}<100)}{#sigma(#chi^{2}<80)}}","P");
  leg->Draw();

  c1->Update();
  c1->Print(pdf.c_str()); // add to pdf-file

  c1->Print((pdf+"]").c_str()); // just close pdf-file
}

