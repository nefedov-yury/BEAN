// compare cross-sections for varying of Mkk window:
//   -> cmp_cs_mkk.pdf

/*
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -> wphi=4*Gphi:
  // Efficiency from MCGPJ:
  double Xisr_min=0.9;
  double eff[17] = {
           0.12261,  0.12327,  0.12420,  0.12520,  0.12697,
           0.12852,  0.13391,  0.13616,  0.13490,  0.13618,
           0.13434,  0.13509,  0.13496,  0.13263,  0.12952,
           0.11712,
           0.12420 // this is 3080_new
                   };

  // Sig, ErrSig - cross-sections in picobarn
  double Xisr_min=0.9;
  Xmax = 1. - Xisr_min;
  vector<double> Sig =      { // pb
        2.57130e+01, 2.86406e+01, 3.31332e+01, 1.71250e+01, 2.89864e+01,
        4.05533e+01, 2.28037e+02, 7.69146e+02, 2.21532e+03, 3.59997e+03,
        1.67910e+03, 4.60547e+02, 1.22198e+02, 6.58921e+01, 9.17964e+01,
        4.14414e+01, 2.49178e+01,
                            };
  vector<double> ErrSig =   { // pb
        5.36154e+00, 5.61688e+00, 5.60053e+00, 9.06168e+00, 5.85614e+00,
        6.57862e+00, 4.03117e+01, 8.06097e+01, 1.28985e+02, 1.64230e+02,
        1.09036e+02, 9.60307e+01, 3.64055e+01, 2.42822e+01, 2.90286e+01,
        2.39262e+01, 1.83109e+00,
                            };


  Efficiency from MCGPJ: s(after ISR)/s(ini) in [0.9,1.0]
 =========================================================

  E(GeV)    Lum.(pb-1)     Signal      Eff.(%)   Cross Section(nb)
 3.049663     14.919     23 +/-  4.8    12.26    0.0257 +/- 0.0054
 3.058707      15.06     26 +/-  5.1    12.33    0.0286 +/- 0.0056
 3.079645     17.393     35 +/-  5.9    12.42    0.0331 +/- 0.0056
 3.082510      4.769      5 +/-  2.6    12.52    0.0171 +/- 0.0091
 3.088868     15.558     28 +/-  5.7    12.70    0.0290 +/- 0.0059
 3.091774      14.91     38 +/-  6.2    12.85    0.0406 +/- 0.0066
 3.094711      2.143     32 +/-  5.7    13.39    0.2280 +/- 0.0403
 3.095444      1.816     93 +/-  9.7    13.62    0.7691 +/- 0.0806
 3.095840      2.135    312 +/- 18.2    13.49    2.2153 +/- 0.1290
 3.097227      2.069    496 +/- 22.6    13.62    3.6000 +/- 0.1642
 3.098354      2.203    243 +/- 15.8    13.43    1.6791 +/- 0.1090
 3.099056      0.756     23 +/-  4.8    13.51    0.4605 +/- 0.0960
 3.101373      1.612     13 +/-  3.9    13.50    0.1222 +/- 0.0364
 3.105594      2.106      9 +/-  3.3    13.26    0.0659 +/- 0.0243
 3.112065       1.72     10 +/-  3.2    12.95    0.0918 +/- 0.0290
 3.119892      1.264      3 +/-  1.7    11.71    0.0414 +/- 0.0239
 3.080000     126.21    191 +/- 14.0    12.42    0.0249 +/- 0.0018
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
/*
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -> wphi=6*Gphi:
  // Efficiency from MCGPJ:
  double Xisr_min=0.9;
  double eff[17] = {
           0.12653,  0.12751,  0.12821,  0.12921,  0.13114,
           0.13267,  0.13807,  0.14071,  0.13924,  0.14020,
           0.13908,  0.13905,  0.13895,  0.13680,  0.13399,
           0.12125,
           0.12821 // this is 3080_new
                   };

  // Sig, ErrSig - cross-sections in picobarn
  double Xisr_min=0.9;
  Xmax = 1. - Xisr_min;
  vector<double> Sig =      { // pb
        2.49164e+01, 2.98181e+01, 3.20969e+01, 1.65935e+01, 3.10716e+01,
        4.03186e+01, 2.14255e+02, 7.60281e+02, 2.21506e+03, 3.67300e+03,
        1.70865e+03, 4.47431e+02, 1.00429e+02, 6.38836e+01, 9.76075e+01,
        4.00298e+01, 2.47703e+01,
                            };
  vector<double> ErrSig =   { // pb
        5.19543e+00, 5.63509e+00, 5.42536e+00, 8.78045e+00, 5.92975e+00,
        6.45614e+00, 3.97033e+01, 7.96284e+01, 1.27588e+02, 1.64279e+02,
        1.08446e+02, 9.32958e+01, 3.76437e+01, 2.35420e+01, 2.94298e+01,
        2.31112e+01, 1.79619e+00,
                            };


  Efficiency from MCGPJ: s(after ISR)/s(ini) in [0.9,1.0]
 =========================================================

  E(GeV)    Lum.(pb-1)     Signal      Eff.(%)   Cross Section(nb)
 3.049663     14.919     23 +/-  4.8    12.65    0.0249 +/- 0.0052
 3.058707      15.06     28 +/-  5.3    12.75    0.0298 +/- 0.0056
 3.079645     17.393     35 +/-  5.9    12.82    0.0321 +/- 0.0054
 3.082510      4.769      5 +/-  2.6    12.92    0.0166 +/- 0.0088
 3.088868     15.558     31 +/-  5.9    13.11    0.0311 +/- 0.0059
 3.091774      14.91     39 +/-  6.2    13.27    0.0403 +/- 0.0065
 3.094711      2.143     31 +/-  5.7    13.81    0.2143 +/- 0.0397
 3.095444      1.816     95 +/-  9.9    14.07    0.7603 +/- 0.0796
 3.095840      2.135    322 +/- 18.5    13.92    2.2151 +/- 0.1276
 3.097227      2.069    521 +/- 23.3    14.02    3.6730 +/- 0.1643
 3.098354      2.203    256 +/- 16.2    13.91    1.7086 +/- 0.1084
 3.099056      0.756     23 +/-  4.8    13.91    0.4474 +/- 0.0933
 3.101373      1.612     11 +/-  4.1    13.89    0.1004 +/- 0.0376
 3.105594      2.106      9 +/-  3.3    13.68    0.0639 +/- 0.0235
 3.112065       1.72     11 +/-  3.3    13.40    0.0976 +/- 0.0294
 3.119892      1.264      3 +/-  1.7    12.12    0.0400 +/- 0.0231
 3.080000     126.21    196 +/- 14.2    12.82    0.0248 +/- 0.0018
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
//-------------------------------------------------------------------------
void cmpr_cs_mkk()
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

// wphi=5*Gphi +++ standard
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
// wphi=4*Gphi: +
  vector<double> SigM =      { // pb
        2.57130e+01, 2.86406e+01, 3.31332e+01, 1.71250e+01, 2.89864e+01,
        4.05533e+01, 2.28037e+02, 7.69146e+02, 2.21532e+03, 3.59997e+03,
        1.67910e+03, 4.60547e+02, 1.22198e+02, 6.58921e+01, 9.17964e+01,
        4.14414e+01, 2.49178e+01,
                            };
  vector<double> ErrSigM =   { // pb
        5.36154e+00, 5.61688e+00, 5.60053e+00, 9.06168e+00, 5.85614e+00,
        6.57862e+00, 4.03117e+01, 8.06097e+01, 1.28985e+02, 1.64230e+02,
        1.09036e+02, 9.60307e+01, 3.64055e+01, 2.42822e+01, 2.90286e+01,
        2.39262e+01, 1.83109e+00,
                            };
// wphi=6*Gphi: +
  vector<double> SigP =      { // pb
        2.49164e+01, 2.98181e+01, 3.20969e+01, 1.65935e+01, 3.10716e+01,
        4.03186e+01, 2.14255e+02, 7.60281e+02, 2.21506e+03, 3.67300e+03,
        1.70865e+03, 4.47431e+02, 1.00429e+02, 6.38836e+01, 9.76075e+01,
        4.00298e+01, 2.47703e+01,
                            };
  vector<double> ErrSigP =   { // pb
        5.19543e+00, 5.63509e+00, 5.42536e+00, 8.78045e+00, 5.92975e+00,
        6.45614e+00, 3.97033e+01, 7.96284e+01, 1.27588e+02, 1.64279e+02,
        1.08446e+02, 9.32958e+01, 3.76437e+01, 2.35420e+01, 2.94298e+01,
        2.31112e+01, 1.79619e+00,
                            };

  vector<double> RM(Nbeam), RP(Nbeam);
  string serr("  double err_mkk[]  = {");
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

  string pdf( "cmp_cs_mkk.pdf" ); // <- PDF
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
    legD->SetHeader("variation of M_{KK} window:","C");
    legD->AddEntry( grW1,
        "#sigma(|M(K^{+}K^{-}) - M_{#phi}|< 5#Gamma_{#phi})", "EP");
    legD->AddEntry( grW2,
        "#sigma(|M(K^{+}K^{-}) - M_{#phi}|< 4#Gamma_{#phi})", "EP");
    legD->AddEntry( grW3,
        "#sigma(|M(K^{+}K^{-}) - M_{#phi}|< 6#Gamma_{#phi})", "EP");
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
  rat1->GetXaxis()->CenterTitle();
  rat1->GetYaxis()->CenterTitle();
  rat1->SetMaximum(15.);
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

  TLegend* leg = new TLegend(0.12,0.52,0.50,0.88);
  leg->SetHeader("variation of M_{KK} window:","C");
  leg->AddEntry(rat1,
      "#1 - #frac{#sigma(4 #times #Gamma(#phi))}"
      "{#sigma(5 #times #Gamma(#phi))}","P");
  leg->AddEntry(rat2,
      "#1 - #frac{#sigma(6 #times #Gamma(#phi))}"
      "{#sigma(5 #times #Gamma(#phi))}","P");
  leg->Draw();

  c1->Update();
  c1->Print(pdf.c_str()); // add to pdf-file

  c1->Print((pdf+"]").c_str()); // just close pdf-file
}

