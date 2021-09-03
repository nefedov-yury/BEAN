//-------------------------------------------------------------------------
Scripts for Psi' -> Jpsi pi+pi-
                     ^-> phi eta     selection:
//-------------------------------------------------------------------------
0) Cuts: cits.h
//-------------------------------------------------------------------------
 -> List of cuts for
        1) Psi' -> J/Psi pi+ pi-
        2) J/Psi -> phi eta

//-------------------------------------------------------------------------
1) Selection Psi' -> Jpsi pi+pi-
//-------------------------------------------------------------------------
num_prints.cc - print simple statistic about production:
                1) initial number of events for data and MC
                2) number of good MC events and
                   number of Psi(2S)->pi+pi-J/Psi (N64)

pipi_dataMC.cc - plot data distributions vs MC
                 for pi+ pi- in selected Psi(2S) -> Jpsi pi+pi-
               -> VarName_YEAR.pdf

MrecFitSB.cc - Calculate number of Psi(2S) -> J/Psi pi+ pi- decays
               in data using side-band method.
               MrecDraw() -> Mrec_YEAR.pdf
               DoFitSB()  -> Mrec_YEAR_fsb.pdf

MrecFit.cc - FOR MC AFTER HELIX CORRECTIONS
             (archive/MrecFit.cc - version before helix corrections)
           - fit data by de-convolution of MC inclusive signal
             with Gauss and scaling of MC inclusive background;
           - calculate number of J/Psi in Psi'->J/Psi pi+ pi- decay
           -> Mrec_YEAR_fit.pdf

//-------------------------------------------------------------------------
2) Selection Jpsi -> phi eta
//           see ../phieta/ with such scripts for J/Psi scan data
//-------------------------------------------------------------------------

mass_2D.cc - 1) plot M(gamma-gamma) vs M(K+K-)
             2) Dalitz plots ( for >=prod-9m )
                type = 1 M^2(K-eta) vs M^2(K+eta)
                type = 2 M^2(K+K-) vs M^2(K+eta)
             after 4C kinematic fit and chi^2 cut (see cuts.h)
             -> Mass2D_xxx.pdf; -> Dalitz_xxx.pdf

mass_eta.cc - plot M(gamma gamma) distributions
              after 4C kinematic fit and chi^2 cut (see cuts.h)
              -> mass_eta.pdf

mass_phi.cc - /archive/ - plot M(K+K-) distributions
              after 4C kinematic fit and chi^2 cut (see cuts.h)
              -> mass_phi.pdf

mass_kk.cc - plot M(K+K-) for data, inclusive MC, signal MC and MC KKeta
             cuts (see cuts.h): Mrec + chi^2(4C) + Mgg
             -> mass_kk.pdf

kkgg_dataMC.cc - data vs MC for  Jpsi -> phi eta  selection (see cuts.h)
               -> variable_YEAR.pdf

chi2_Pr.cc - plot chi2 of 4C kinematic fit
             1) data vs inclusive MC
                -> chi2_YEAR.pdf
             2) data vs signal MC in the window for M(K+K-) and M(gg)
                -> chi2_sb_YEAR.pdf

mrec_4C.cc - plot Mrec after 4C-kinematik fit
             we require cuts for: chi^2(4C) && |Mgg-Meta| && |Mkk-Mphi|
             Data, inclusive MC, signal MC and MC KKeta
             -> mrec_4Cyear.pdf

Mkk_tail.cc - study the tail in M(K+K-) distribution
              0) plot Mkk for data
              plot data vs MC for 1.1 < Mkk < 2.0 GeV/c^2
               or 1.08 < Mkk < 1.3 GeV/c^2 in J/Psi rest frame
              1) Mgg
              2) cos(Theta(Eta))
              3) cos(Theta(K+/-))
              4) P(K+/-)
              5) cos(K+,K-)

//-------------------------------------------------------------------------
// 2a) Study M(K+K-)
//-------------------------------------------------------------------------
eff_mc.cc -  efficiencies for selection of the J/Psi -> phi eta
             as a function of true MC M(K+K-) and
             fit it by a straight line
             -> eff_(sig|bkg)_YYYY.pdf

mass_kk_fit.cc - unbinned LH fit of M(K+K-) distributions
                 after the cuts (see cuts.h): Mrec + chi^2(4C) + Mgg
                 *) Breit-Wigner for psi -> K+K- (psi from J/Psi -> phi eta)
                 *) convolution of BW with the Gauss distribution
                 *) ARGUS function (and "reverse" Argus)
                 *) interference B-W and Argus amplitudes
                 -> mkkYY_???.pdf
                 -> mkk_cf_???.pdf (combined fit)

OLD: mkk_fitch2.cc - chi^2 fit (binning ???)

//-------------------------------------------------------------------------
// Study systematic uncertainties:
//      track efficiencies
//      photon reconstruction efficiencies
//-------------------------------------------------------------------------
trk_eff_sel.cc - Pictures for presentation/memo
        Study of the track reconstruction efficiency for pi and K:


trk_eff_fit.cc - The reconstruction efficiency (data and MC)
                 and their ratio (data/MC) for pions and kaons.
                 Get corrections as function of P_t
                 by fitting the ratio.
                 Kolmogorov–Smirnov probability test to compare
                 the ratios for K+ and K- (pi+ and pi-)

trk_eff_wts.cc -  plot weights for K+K- and pi+ pi- pairs:
                  -> wts_KK_${date}.pdf; wts_PiPi_${date}.pdf;

OLD: archive/
trk_eff.cc - Study of the track reconstruction efficiency for pi and K:
             I)  plot 1D distributions
             ->  Trkeff_[date].pdf
             II) Kolmogorov–Smirnov test to compare K+ and K- (pi+ and pi-)

trk_eff_fit.cc - Study of the track reconstruction efficiency for pi and K:
                 fit the ratio data/MC as function of cos(Theta) and P_t
                 -> trk_fit_[date].pdf
//-------------------------------------------------------------------------

eta_eff.cc - Study of the eta -> eta reconstruction efficiency &
                                 single photon rec.efficiency
           -> Eff_[eta,ph]_[date]_[gamma,phi]eta.pdf

eta_eff_pr.cc - Pictures for presentation/memo

//-------------------------------------------------------------------------
// Calculation of branching
//-------------------------------------------------------------------------
eff_mc.cc   calculate efficiencies for selection of the J/Psi -> phi eta
            > root -q -b eff_mc.cc


