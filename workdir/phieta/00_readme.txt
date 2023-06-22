------------------------------
Scripts for phi-eta selection:
------------------------------

//--------------------------------------------------------------------
Common headers
//--------------------------------------------------------------------
masses.h   -> masses of particles (GeV) from PDG 2022

cuts.h -> List of cuts
// List of constants and cuts

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

-- mass_KK.cc: --memo--
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

mcgpj_plots.cc
// plot true distributions for MCGPJ generator:
// 1) beam-spread       -> mcgpj_bs_Exxxx.pdf
// 2) cos Theta (eta)   -> mcgpj_ang_Exxxx.pdf
// 3) X_{ISR} distributions before and after selection
//    -> Xisr_XXXX.pdf

//--------------------------------------------------------------------
Background study
//--------------------------------------------------------------------
-- jpsi_incl.cc: --memo--
// plot eff, purity and signal/background ration as function of chi^2
// cut estimated by Monte Carlo for inclusive decays of J/Psi
// -> chi2_eff_pur.pdf
> purity for ch2<80 is 95.43%

//--------------------------------------------------------------------
Efficiency && Mkk fitting
//--------------------------------------------------------------------
-- xisr_eff.cc: --memo--
// plot efficiency as a function of X_isr
// NOTE: there is no cut Mkk < 1.08 in 'mcisr_xisr' histo
// -> Xisr_eff_XXXX.pdf

eff_MC.cc:
// selection efficiency of the J/Psi -> phi eta as a function
// of true M(K+K-). We fit this efficiency with constant and
// linear dependences on M(K+K-)
// -> eff_sig_XXX.pdf

mass_KK_fit.cc
// fit distributions M(K+K-) for data
// NOTE: see also PsipJpsiPhiEta/mass_kk_fit.cc and mkk_fitch2.cc
//       memo_brf_v15a: phase_angle= 0.05^+0.52_-0.63 = [-0.58,+0.58]
// -> outputs:
//    'mkk_inter/' for interference model BW with Argus
//    'mkk_noint/' for model with sum BW + Argus

get_table.py: -- memo --
# print a summary table with fitting results

//--------------------------------------------------------------------
Cross section
//--------------------------------------------------------------------
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

15Jun23: mass_KK_fit, MCGPJ_03

//--------------------------------------------------------------------
Study of uncertainties
//--------------------------------------------------------------------
sys_lumi.py
# systematic uncertainties for luminosity (pb-1)
-> sys_lumi.h

cmpr_cs.cc
// calculation tables of systematic uncertainties for
//  1) variation of the chi^2
//  2) variation of the Eta-window
//  3) variation of the mixing angle
// compare cross-sections

sys_ANGL_15Jun23.h: theta= +/- 0.58
chi^2: rat= 2.22 %, E= 3.097654
chi^2: rat= 1.75 %, E= 3.096203

Weta: rat= 3.25 %, E= 3.097654
Weta: rat= 3.17 %, E= 3.095726

sys_mcgpj.cc
// study of systematic uncertainties in efficiency due to the
// variation of parameters A and phi within their errors
// in MCGPJ generator
// DIR is the tested directory (see DIR/00_readme)
//    -> sys_mcgpj_DATE.pdf

