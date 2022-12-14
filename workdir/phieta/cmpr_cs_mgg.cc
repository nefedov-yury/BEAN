// compare cross-sections for varying of side-band window:
//   -> cmp_cs_mgg.pdf

/*
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -> weta=2*sigma_eta
  // Sig, ErrSig - cross-sections in picobarn
  double Xisr_min=0.9;
  Xmax = 1. - Xisr_min;
  vector<double> Sig =      { // pb
        2.40626e+01, 2.82797e+01, 3.03134e+01, 2.18910e+01, 3.07289e+01,
        4.07091e+01, 2.09515e+02, 7.70231e+02, 2.26959e+03, 3.57972e+03,
        1.67886e+03, 4.39908e+02, 1.08192e+02, 6.24106e+01, 8.26878e+01,
        5.17278e+01, 2.35332e+01,
                            };
  vector<double> ErrSig =   { // pb
        5.64319e+00, 6.00827e+00, 5.53444e+00, 8.93698e+00, 6.01104e+00,
        7.15185e+00, 4.09843e+01, 8.25566e+01, 1.33736e+02, 1.69059e+02,
        1.13242e+02, 1.00463e+02, 3.80934e+01, 2.46699e+01, 2.92346e+01,
        2.98650e+01, 1.90421e+00,
                            };


  Efficiency from MCGPJ: s(after ISR)/s(ini) in [0.9,1.0]
 =========================================================

  E(GeV)    Lum.(pb-1)     Signal      Eff.(%)   Cross Section(nb)
 3.049663     14.919     20 +/-  4.7    11.39    0.0241 +/- 0.0056
 3.058707      15.06     24 +/-  5.1    11.52    0.0283 +/- 0.0060
 3.079645     17.393     30 +/-  5.5    11.64    0.0303 +/- 0.0055
 3.082510      4.769      6 +/-  2.4    11.75    0.0219 +/- 0.0089
 3.088868     15.558     28 +/-  5.5    11.98    0.0307 +/- 0.0060
 3.091774      14.91     36 +/-  6.3    12.13    0.0407 +/- 0.0072
 3.094711      2.143     28 +/-  5.5    12.75    0.2095 +/- 0.0410
 3.095444      1.816     89 +/-  9.5    13.01    0.7702 +/- 0.0826
 3.095840      2.135    305 +/- 18.0    12.87    2.2696 +/- 0.1337
 3.097227      2.069    473 +/- 22.3    13.06    3.5797 +/- 0.1691
 3.098354      2.203    233 +/- 15.7    12.88    1.6789 +/- 0.1132
 3.099056      0.756     21 +/-  4.8    12.91    0.4399 +/- 0.1005
 3.101373      1.612     11 +/-  3.9    12.90    0.1082 +/- 0.0381
 3.105594      2.106      8 +/-  3.2    12.45    0.0624 +/- 0.0247
 3.112065       1.72      8 +/-  2.8    11.50    0.0827 +/- 0.0292
 3.119892      1.264      3 +/-  1.7     9.38    0.0517 +/- 0.0299
 3.080000     126.21    169 +/- 13.7    11.64    0.0235 +/- 0.0019
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
/*
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -> weta=2.5*sigma_eta
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
/*
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -> weta=3.5*sigma(eta)
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
/*
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -> weta=4*sigma_eta
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Sig, ErrSig - cross-sections in picobarn
  double Xisr_min=0.9;
  Xmax = 1. - Xisr_min;
  vector<double> Sig =      { // pb
        2.45325e+01, 2.84122e+01, 3.17110e+01, 1.63680e+01, 3.07087e+01,
        4.19186e+01, 2.19119e+02, 7.46969e+02, 2.21324e+03, 3.58246e+03,
        1.72777e+03, 4.44681e+02, 1.26652e+02, 7.05998e+01, 8.75836e+01,
        3.83897e+01, 2.50968e+01,
                            };
  vector<double> ErrSig =   { // pb
        5.33316e+00, 5.46794e+00, 5.36014e+00, 8.66113e+00, 5.69057e+00,
        6.54658e+00, 3.99273e+01, 7.86661e+01, 1.27064e+02, 1.62294e+02,
        1.08787e+02, 9.66697e+01, 3.38493e+01, 2.44565e+01, 2.76964e+01,
        2.21643e+01, 1.81370e+00,
                            };


  Efficiency from MCGPJ: s(after ISR)/s(ini) in [0.9,1.0]
 =========================================================

  E(GeV)    Lum.(pb-1)     Signal      Eff.(%)   Cross Section(nb)
 3.049663     14.919     23 +/-  5.0    12.85    0.0245 +/- 0.0053
 3.058707      15.06     27 +/-  5.2    12.90    0.0284 +/- 0.0055
 3.079645     17.393     35 +/-  5.9    12.98    0.0317 +/- 0.0054
 3.082510      4.769      5 +/-  2.6    13.10    0.0164 +/- 0.0087
 3.088868     15.558     31 +/-  5.7    13.27    0.0307 +/- 0.0057
 3.091774      14.91     41 +/-  6.4    13.41    0.0419 +/- 0.0065
 3.094711      2.143     32 +/-  5.8    13.94    0.2191 +/- 0.0399
 3.095444      1.816     94 +/-  9.9    14.17    0.7470 +/- 0.0787
 3.095840      2.135    324 +/- 18.6    14.02    2.2132 +/- 0.1271
 3.097227      2.069    512 +/- 23.2    14.13    3.5825 +/- 0.1623
 3.098354      2.203    260 +/- 16.4    13.97    1.7278 +/- 0.1088
 3.099056      0.756     23 +/-  5.0    13.99    0.4447 +/- 0.0967
 3.101373      1.612     14 +/-  3.7    14.02    0.1267 +/- 0.0338
 3.105594      2.106     10 +/-  3.5    13.75    0.0706 +/- 0.0245
 3.112065       1.72     10 +/-  3.2    13.58    0.0876 +/- 0.0277
 3.119892      1.264      3 +/-  1.7    12.64    0.0384 +/- 0.0222
 3.080000     126.21    201 +/- 14.5    12.98    0.0251 +/- 0.0018
*/

//-------------------------------------------------------------------------
void cmpr_cs_mgg()
//-------------------------------------------------------------------------
{
  gROOT->Reset();

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLegendTextSize(0.04);

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

// standard +++
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
// weta=2*sigma_eta +
  vector<double> SigM =      { // pb
        2.40626e+01, 2.82797e+01, 3.03134e+01, 2.18910e+01, 3.07289e+01,
        4.07091e+01, 2.09515e+02, 7.70231e+02, 2.26959e+03, 3.57972e+03,
        1.67886e+03, 4.39908e+02, 1.08192e+02, 6.24106e+01, 8.26878e+01,
        5.17278e+01, 2.35332e+01,
                            };
  vector<double> ErrSigM =   { // pb
        5.64319e+00, 6.00827e+00, 5.53444e+00, 8.93698e+00, 6.01104e+00,
        7.15185e+00, 4.09843e+01, 8.25566e+01, 1.33736e+02, 1.69059e+02,
        1.13242e+02, 1.00463e+02, 3.80934e+01, 2.46699e+01, 2.92346e+01,
        2.98650e+01, 1.90421e+00,
                            };
// weta=2.5*sigma_eta
// weta=3.5*sigma_eta
// weta=4*sigma_eta
  vector<double> SigP =      { // pb
        2.45325e+01, 2.84122e+01, 3.17110e+01, 1.63680e+01, 3.07087e+01,
        4.19186e+01, 2.19119e+02, 7.46969e+02, 2.21324e+03, 3.58246e+03,
        1.72777e+03, 4.44681e+02, 1.26652e+02, 7.05998e+01, 8.75836e+01,
        3.83897e+01, 2.50968e+01,
                            };
  vector<double> ErrSigP =   { // pb
        5.33316e+00, 5.46794e+00, 5.36014e+00, 8.66113e+00, 5.69057e+00,
        6.54658e+00, 3.99273e+01, 7.86661e+01, 1.27064e+02, 1.62294e+02,
        1.08787e+02, 9.66697e+01, 3.38493e+01, 2.44565e+01, 2.76964e+01,
        2.21643e+01, 1.81370e+00,
                            };

  vector<double> RM(Nbeam), RP(Nbeam);
  string serr("  double err_mgg[]  = {");
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

  string pdf( "cmp_cs_mgg.pdf" ); // <- PDF
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

    TLegend* legD = new TLegend(0.12,0.52,0.54,0.88);
    legD->SetHeader("variation of side band for M(#gamma#gamma):","C");
    legD->AddEntry( grW1, "#sigma for normal cuts", "EP");
    legD->AddEntry( grW2, "#sigma for tight cuts", "EP");
    legD->AddEntry( grW3, "#sigma for loose cuts", "EP");
    legD->Draw();

    c1->Print(pdf.c_str()); // add to pdf-file
  } // end of if(DEBUG)

  c1->cd();
  gPad->SetGrid();
  gPad->SetLogy(false);

  TGraph* rat1 = new TGraph(Nbeam, ebeam,RM.data());
  TGraph* rat2 = new TGraph(Nbeam, ebeam,RP.data());
  string tit1 = string(";center of mass energy (GeV)")+
                string(";relative difference (%)");
  rat1->SetTitle(tit1.c_str());
  rat1->GetXaxis()->CenterTitle();
  rat1->GetYaxis()->CenterTitle();
  rat1->SetMaximum(20.);
  rat1->SetMinimum(-10.);
  rat1->SetMarkerColor(kBlue+2);
  rat1->SetMarkerStyle(20);
  rat1->SetMarkerSize(1.1);
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

//   TLegend* leg = new TLegend(0.52,0.52,0.9,0.9);
  TLegend* leg = new TLegend(0.12,0.60,0.50,0.88);
  leg->SetHeader("variation of side band for M(#gamma#gamma):","C");
  leg->AddEntry(rat1,
      "#1 - #frac{#sigma(tight)}{#sigma(norm)}","P");
//       "#||{1 - #frac{#sigma(tight)}{#sigma(norm)}}","P");
  leg->AddEntry(rat2,
      "#1 - #frac{#sigma(loose)}{#sigma(norm)}","P");
//       "#||{1 - #frac{#sigma(loose)}{#sigma(norm)}}","P");
  leg->Draw();

  c1->Update();
  c1->Print(pdf.c_str()); // add to pdf-file

  c1->Print((pdf+"]").c_str()); // just close pdf-file
}

