// List of cuts

//-------------------------------------------------------------------------
// 1) Psi' -> J/Psi pi+ pi- 
// Side band for recoil mass of pi+ pi -
//-------------------------------------------------------------------------
  const double Mjpsi[] = {
                                3.09678, // data09
                                3.09693, // data12
                                3.09657  // MC
                         };

  double mjpsi = Mjpsi[MIDX];
  const double wjpsi = 2e-3;

  TCut c_MCtrue("dec==64");
  
  // central part
  string str_cp = "abs(Mrec-" + to_string(mjpsi) +
                      ")<" + to_string(wjpsi);
  TCut c_cp(str_cp.c_str());

  // side-band
  string str_left = "abs(Mrec-" + to_string(mjpsi-5.5*wjpsi) +
                      ")<" + to_string(wjpsi/2);
  TCut c_left(str_left.c_str());
  string str_right = "abs(Mrec-" + to_string(mjpsi+5.5*wjpsi) +
                      ")<" + to_string(wjpsi/2);
  TCut c_right(str_right.c_str());

//-------------------------------------------------------------------------
// 2) J/Psi -> phi eta 
// Side band for Meta - invariant mass of 2gamma (eta -> 2gamma)
// Window for Mphi - invariant mass of K+ K- (phi -> K+ K-)
//-------------------------------------------------------------------------
/* 
  const double mphi = 1.019461; // 1019.461 +/- 0.019   MeV
  const double meta = 0.547862; //  547.862 +/- 0.017   MeV
  const double Gphi = 4.247e-3; //   4.247 +/- 0.016 MeV
  const double seta = 0.008;

  // chi^2 cut:
  TCut c_chi2("ch2<80.");  // standard
//   TCut c_chi2("ch2<60.");  // uncertainties study
//   TCut c_chi2("ch2<100."); // uncertainties study

  // Mphi cut:
  const double wphi = 5*Gphi; // standard
//   const double wphi = 4*Gphi; // uncertainties study
//   const double wphi = 6*Gphi; // uncertainties study
  string str_cp_phi = "abs(Mkk-" + to_string(mphi) +
                      ")<" + to_string(wphi);
//   cout << " ++++ str_cp_phi= " << str_cp_phi << endl;
  TCut cp_phi(str_cp_phi.c_str());

  // side-band for Meta; central part
  const double weta = 3*seta; // standard
//   const double weta = 2*seta; // uncertainties study
//   const double weta = 4*seta; // uncertainties study
  TCutG* cg_cp = new TCutG("cp",5);
  cg_cp->SetVarX("Mgg");
  cg_cp->SetVarY("Mkk");
  cg_cp->SetPoint(0,meta-weta,mphi-wphi);
  cg_cp->SetPoint(1,meta-weta,mphi+wphi);
  cg_cp->SetPoint(2,meta+weta,mphi+wphi);
  cg_cp->SetPoint(3,meta+weta,mphi-wphi);
  cg_cp->SetPoint(4,meta-weta,mphi-wphi);
  cg_cp->SetLineColor(kRed+1);
  cg_cp->SetLineStyle(1);
  cg_cp->SetLineWidth(1);

  // side-band for Mgg only:
  TCutG* cg_etaL = new TCutG("etaL",5);
  cg_etaL->SetVarX("Mgg");
  cg_etaL->SetVarY("Mkk");
  cg_etaL->SetPoint(0,meta-seta-2*weta,mphi-wphi);
  cg_etaL->SetPoint(1,meta-seta-2*weta,mphi+wphi);
  cg_etaL->SetPoint(2,meta-seta-weta,mphi+wphi);
  cg_etaL->SetPoint(3,meta-seta-weta,mphi-wphi);
  cg_etaL->SetPoint(4,meta-seta-2*weta,mphi-wphi);
  TCutG* cg_etaR = new TCutG("etaR",5);
  cg_etaR->SetVarX("Mgg");
  cg_etaR->SetVarY("Mkk");
  cg_etaR->SetPoint(0,meta+seta+weta,mphi-wphi);
  cg_etaR->SetPoint(1,meta+seta+weta,mphi+wphi);
  cg_etaR->SetPoint(2,meta+seta+2*weta,mphi+wphi);
  cg_etaR->SetPoint(3,meta+seta+2*weta,mphi-wphi);
  cg_etaR->SetPoint(4,meta+seta+weta,mphi-wphi);
  cg_etaL->SetLineColor(kBlue+1);
  cg_etaL->SetLineStyle(1);
  cg_etaL->SetLineWidth(1);
  cg_etaR->SetLineColor(kBlue+1);
  cg_etaR->SetLineStyle(1);
  cg_etaR->SetLineWidth(1);
*/
//-------------------------------------------------------------------------
