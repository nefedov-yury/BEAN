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

chi2_Pr.cc:
// plot chi2 of the kinematic fit
// data (window for Mkk is used) for signal ans side-band events
// versus signal MC (plus another MC if any):
// -> chi2_sb_XXXXX.pdf

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
//    -> cs_results.tex or .txt : final tables in TeX format

29Jan23: KK_fit, latest MC : -> old
05Feb23: just side-band subtraction, latest MC : -> old

10Feb23: side-band subtraction, tight cut: 1.01 < Mkk < 1.03 GeV
17Mar23: mass_KK_fit (LH), more accurate parameters 

//--------------------------------------------------------------------
Study of uncertainties
//--------------------------------------------------------------------
cmpr_cs.cc
// compare cross-sections
//   -> cmpr_cs_XXXX.pdf

