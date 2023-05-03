//--------------------------------------------------------------------
Scripts for Psi' -> Jpsi pi+pi-
                     ^-> phi eta     selection:
//--------------------------------------------------------------------
0) Common headers
//--------------------------------------------------------------------
masses.h   -> masses of particles (GeV) from PDG 2021

-- norm.h     -> normalization constants for inclusive MC and for
              non-resonance data. v709 ?

cuts.h     -> List of cuts for
              1) Psi' -> J/Psi pi+ pi-
              2) J/Psi -> phi eta

p10a4c.h   -> declaration of variables in prod-10 a4c ntuple
p10a4cMC.h -> declaration of variables in prod-10mc & prod-11 a4c
              ntuple

ReWeightEtaEff.h -> function to re-weight MC for eta->2gamma
                    efficiency

ReWeightTrkPid_11.h  -> function to re-weight MC for K and Pi trk+pid
                        efficiency (prod-11)
                        for other versions see archive/

//--------------------------------------------------------------------
1) Selection Psi' -> Jpsi pi+pi-
//--------------------------------------------------------------------
pipi_dataMC.cc - plot data distributions vs MC
                 for pi+ pi- in selected Psi(2S) -> Jpsi pi+pi-
               -> VarName_YEAR.pdf

MrecFitSB.cc - Calculate number of Psi(2S) -> J/Psi pi+ pi- decays
               in data using side-band method.
               MrecDraw() -> Mrec_YEAR.pdf
               DoFitSB()  -> Mrec_YEAR_fsb.pdf
// parameters for V6.6.4
if ( Npol == 3 ) {
  if ( date == 2009 ) {
    par_ini = { 1.0343, 0.0069, 0.0083, 0.0011 };
  } else if ( date == 2012 ) {
    par_ini = { 1.0478, 0.0044, 0.0076, 0.0008 };
  }
}

MrecFit.cc - FOR MC AFTER HELIX CORRECTIONS
             (archive/MrecFit.cc - version before helix corrections)
           - data fitting by correcting of MC signal and scaling of MC
             background
             -> Mrec_YEAR_fitMODEL.pdf
             -> diffYEAR_MODEL.pdf for delta(data-MC)
           - calculate number of J/Psi in Psi' -> J/Psi pi+ pi- decay
             and efficiency of selection;
           - estimate systematic associated with a fit model

//--------------------------------------------------------------------
2) Selection Jpsi -> phi eta
//           see ../phieta/ with such scripts for J/Psi scan data
//--------------------------------------------------------------------

mass_2D.cc - 1) plot M(gamma,gamma) vs M(K+K-)
                -> Mass2D_XXXX.pdf;

             2) Dalitz plots (for >=prod-9m)
                after 4C kinematic fit and chi^2 cut
                  type = 1 M^2(K-eta) vs M^2(K+eta)
                  type = 2 M^2(K+K-) vs M^2(K+eta)
                -> Dalitz_xxx.pdf

mass_eta.cc - plot M(gamma gamma) distributions
              after 4C kinematic fit and chi^2 cut
              -> mass_eta.pdf

chi2_Pr.cc - plot chi2 of the kinematic fit
             1) data vs signal MC in the window for M(K+K-) and M(gg)
                -> chi2_sb_YEAR.pdf
             2) data vs inclusive MC
                -> chi2_YEAR.pdf

kkgg_dataMC.cc - data vs MC for  Jpsi -> phi eta  selection
               -> variable_YEAR.pdf

mass_kk.cc - plot M(K+K-) for data, inclusive MC, signal MC and
             MC KKeta cuts (see cuts.h): Mrec + chi^2(4C) + Mgg
             -> mass_kk.pdf

//--------------------------------------------------------------------
// 2a) Study M(K+K-)
//--------------------------------------------------------------------
eff_mc.cc -  efficiencies for selection of the J/Psi -> phi eta
             as a function of true MC M(K+K-) and
             fit it by a straight line
             -> eff_(sig|bkg)_YYYY.pdf

mass_kk_fit.cc - unbinned LH fit of M(K+K-) distributions
                 after the cuts (see cuts.h): Mrec + chi^2(4C) + Mgg
                 *) Breit-Wigner for phi -> K+K-
                    (phi from J/Psi -> phi eta)
                 *) convolution of BW with the Gauss distribution
                 *) ARGUS function (and "reverse" Argus)
                 *) interference B-W and Argus amplitudes
                 -> mkkYY_???.pdf
                 -> mkk_cf_???.pdf (combined fit)

//--------------------------------------------------------------------
// Study systematic uncertainties:
//--------------------------------------------------------------------
mrec_4C.cc - Here we use special selection with 4C-kinematic
             constrints for data (2009,2012 and continuum) and
             inclusive MC events.
             The 'Mrec' is drawn and the contribution of
             Psi(2S) -> pi+ pi- phi eta (non J/Psi) is calculated.
             -> mrec_4Cyear.pdf

//--------------------------------------------------------------------
//      K & Pi track+PID efficiencies
//--------------------------------------------------------------------
trk_eff_sel.cc - Pictures for presentation/memo
                 Study of the track reconstruction efficiency for
                 pi and K


trk_eff_fit.cc - The reconstruction efficiency (data and MC)
                 and their ratio (data/MC) for pions and kaons.
                 Get corrections as function of P_t
                 by fitting the ratio.
                 Kolmogorov–Smirnov  && Chi2 probability tests to
                 compare the ratios for K+ and K- (pi+ and pi-)

trk_eff_wts.cc - plot weights for K+K- and pi+ pi- pairs:
                 -> wts_KK_${date}.pdf; wts_PiPi_${date}.pdf;

//--------------------------------------------------------------------
//      eta->2gamma reconstruction efficiencies
//--------------------------------------------------------------------
eta_eff_sel.cc - Pictures for presentation/memo
                 Study of the eta -> 2gamma reconstruction efficiency

eta_eff.cc -     Study of the eta -> 2gamma reconstruction efficiency
                 single photon rec.efficiency
                 -> Eff_[eta,ph]_[date]_[gamma,phi]eta.pdf

//--------------------------------------------------------------------
// Calculation of branching:
//--------------------------------------------------------------------
Lmedian.cc - calculate confidence interval

OLD:
br_aver.py - average of two periods of data taking

