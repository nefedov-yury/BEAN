// cuts.h : list of cuts
#include "masses.h"

//--------------------------------------------------------------------
// 1) Psi' -> J/Psi pi+ pi-
//--------------------------------------------------------------------
// cut on recoil mass: [3.092, 3.102]
TCut c_Mrec("abs(Mrec-3.097)<0.005");

//--------------------------------------------------------------------
// 2) J/Psi -> phi eta
//--------------------------------------------------------------------
TCut c_MCmkk("mcmkk<1.08"); // cut for MC generated events

// chi^2 cut: use helix correction ON/OFF to systematic study
// double chi2M = 40; // prod v604, chi2M = 60: very old
double chi2M = 60;
TCut c_chi2( Form("ch2<%.1f",chi2M) );

// Meta cut: central part (c_cpgg) and side-band (c_sbgg)
// Mgg is an invariant mass of 2gammas
const double seta = 0.008;
const double weta = 3*seta; // standard: 3x,
                            // uncertainties study: 2x, 4x
TCut c_cpgg( Form("abs(Mgg-%.6f)<%.6f",Meta,weta) );

// shift_eta is the starting position of the side-band
double shift_eta = 7*seta;
// double shift_eta = 6*seta; // old (prod<=10)
TCut c_sbgg( Form("abs(Mgg-%.6f)>%.6f&&abs(Mgg-%.6f)<%.6f",
         Meta,shift_eta,Meta,shift_eta+weta) );

// window for Mphi: [2*Mk, 1.08GeV], only for pictures
// Mkk is an invariant mass of K+K-
TCut c_phi( Form("Mkk>%.6f&&Mkk<1.08",2*Mk) );

//--------------------------------------------------------------------
