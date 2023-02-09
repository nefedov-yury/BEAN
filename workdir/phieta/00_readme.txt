------------------------------
Scripts for phi-eta selection:
------------------------------

//--------------------------------------------------------------------
Common headers
//--------------------------------------------------------------------
masses.h   -> masses of particles (GeV) from PDG 2021

TODO: cuts.h -> List of cuts
// List of constants and cuts

Ntpls/a4c_prod2.C -> declaration of variables in prod-02
Ntpls/a4c_prod3.C -> declaration of variables in prod-03
> a4c->MakeCode("a4c_prod3.C")

//--------------------------------------------------------------------
Selection: data & MC
//--------------------------------------------------------------------
chi2_Pr.cc:
// plot chi2 of the kinematic fit
// data vs signal MC in the window for M(K+K-) and M(gg)
// -> chi2_sb_XXXXX.pdf

mass_2D.cc:
// after 4C kinematic fit and chi^2 cut:
// 1) plot M(gamma,gamma) vs M(K+K-)
//    -> Mass2D_XXXX.pdf;
//
// 2) Dalitz plots ( for >=prod-2 )
//    type = 1 M^2(K-eta) vs M^2(K+eta)
//    type = 2 M^2(K+K-) vs M^2(K+eta)
//    -> Dalitz_xxx.pdf

mass_eta.cc:
// plot M(gamma gamma) distributions
// after chi^2 cut for Mphi region
// -> mass_eta_(PR).pdf

mass_KK.cc:
// plot M(K+K-) distributions
// cuts: chi^2(4C) + Mgg
// -> mass_kk_(PR).pdf

plot_data_mc.cc
// plot Data vs MC distributions after all cuts applied
// -> data_vs_mc_XXXX.pdf

//--------------------------------------------------------------------
Background study
//--------------------------------------------------------------------
jpsi_incl.cc:
// plot eff, purity and signal/background ration as function of chi^2
// cut estimated by Monte Carlo for inclusive decays of J/Psi
// -> chi2_eff_pur.pdf
> purity for ch2<80 is 95.43%

//--------------------------------------------------------------------
Efficiency && Cross section
//--------------------------------------------------------------------
xisr_Pr.cc:
// To illustratt why X_{ISR} > 0.9 is chosen as the signal definition
// area, we plot X_{ISR} distributions (from MCGPJ) before and after
// selection
// -> Xisr_XXXX.pdf

xisr_eff.cc
// plot efficiency as a function of X_isr
// -> Xisr_eff_XXXX.pdf

eff_MC.cc:
// selection efficiency of the J/Psi -> phi eta
// as a function of true MC M(K+K-)
// fit it by a straight line
// -> eff_sig.pdf

mass_KK_fit.cc
// fit distributions M(K+K-) for data
// NOTE: see also PsipJpsiPhiEta/mass_kk_fit.cc and mkk_fitch2.cc
// -> outputs:
//    'mkk_inter/' for interference model BW with Argus
//    'mkk_noint/' for model with sum BW + Argus

cross_section.cpp:
// estimating the cross-section values from the data:
// I ) sinple subtraction of the side-band
// II) after mass_KK_fit.cc
//  *) plot efficiency of phi-eta selection on the base of MCGPJ
//     -> eff_12/18/R/all.pdf
//  *) plot/print cross-section for each energy point
//    -> cp_12/18/R/all.pdf
//    -> cs_results.h   : cpp-code with cross-section and energies
//    -> cs_results.txt : final tables

29Jan23: KK_fit, latest MC
05Feb23: just side-band subtraction, latest MC

//--------------------------------------------------------------------
Study of uncertainties
//--------------------------------------------------------------------
cmpr_cs_chi2.cc
        compare cross-sections for varying of chi^2 cut:  80 +/- 20
        -> cmp_cs_chi2.pdf
  double err_chi2[] = { 6.87, 2.57, 1.90, 18.77, 1.60,
                        1.51, 3.95, 3.07,  0.85, 1.02,
                        0.75, 1.45, 1.38,  2.09, 3.49,
                        8.14, 1.90 };
cmpr_cs_mkk.cc
        compare cross-sections for varying of Mkk window:
        -> cmp_cs_mkk.pdf
  double err_mkk[]  = { 1.83, 6.20,  1.83, 1.85, 4.76,
                        1.24, 4.33,  2.00, 0.73, 2.35,
                        2.78, 1.71, 10.41, 1.79, 8.51,
                        2.12, 0.85 };

cmpr_cs_mgg.cc
        compare cross-sections for varying of side-band window:
        -> cmp_cs_mgg.pdf
  double err_mgg[]  = { 4.71, 1.19,  6.84, 30.19, 0.97,
                        2.68, 6.45,  2.15,  1.71, 0.25,
                        2.80, 2.85, 14.44,  9.06, 8.07,
                       27.46, 5.80 };

xisr_uncer.cc - calculate the uncertainties of efficiency because of ISR
        -> xisr_un_${uncer}.pdf
A-dA
  double err_mc[] = { 0.66, 0.50, 0.67, 1.25, 1.26,
                      1.61, 0.69, 0.27, 0.19, 0.09,
                      0.20, 0.46, 1.33, 2.09, 2.03,
                      1.60, 0.67 };
A+dA
  double err_mc[] = { 0.57, 0.49, 0.71, 0.53, 1.30,
                      1.45, 0.46, 0.11, 0.08, 0.05,
                      0.09, 0.24, 0.99, 1.30, 2.02,
                      1.31, 0.71 };
phi - d phi
  double err_mc[] = { 0.82, 1.70, 2.63, 3.23, 3.51,
                      2.91, 0.54, 0.15, 0.07, 0.06,
                      0.14, 0.30, 1.25, 1.97, 2.86,
                      3.08, 2.63 };
phi + d phi
  double err_mc[] = { 0.24, 0.08, 0.31, 0.06, 0.13,
                      0.14, 0.07, 0.05, 0.02, 0.01,
                      0.01, 0.08, 0.04, 0.37, 0.30,
                      0.27, 0.31 };



-00- cmpr_cs.cc - compare cross-sections for different efficiencies
             after adjusting phase parameters in MCGPJ
             -> cmp_cs_X.pdf

//--------------------------------------------------------------------
sys_uns.py - calculate total systematic uncertainties

for cross-section:

  double err_sys[] = { 7.08, 7.23, 7.51, 7.75, 7.86,
                       7.61, 7.07, 7.05, 7.05, 7.06,
                       7.04, 7.09, 7.17, 7.31, 7.62,
                       7.70, 7.46 };

for fit of phase:

exclude all not energy dependencies
  double err_sys[] = { 1.35, 2.01, 2.85, 3.43, 3.67,
                       3.10, 1.31, 1.19, 1.20, 1.28,
                       1.17, 1.44, 1.79, 2.29, 3.12,
                       3.32, 2.72 };

energy dependencies for lum,chi2,mkk,mgg and mcgpg
  double err_sys[] = { 8.63, 7.11, 7.87, 35.76, 6.29,
                       4.54, 8.81, 4.41, 2.37, 2.87,
                       4.18, 3.90, 17.94, 9.74, 12.63,
                       28.91, 6.74 };


//--------------------------------------------------------------------
