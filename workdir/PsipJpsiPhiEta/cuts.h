// cuts.h - list of cuts
// ATTENSION: to_string() converts only 6 digits after the point

#include "masses.h"

//--------------------------------------------------------------------
// 1) Psi' -> J/Psi pi+ pi-
//--------------------------------------------------------------------
   TCut c_MCtrue("dec==64"); // only Psi' -> J/Psi pi+ pi- decay

   // cut on recoil mass: [3.092, 3.102]
   TCut c_Mrec("abs(Mrec-3.097)<0.005");

   // cut on recoil mass the same as in BAM-0042 [3.055, 3.145]
   TCut c_MrecBAM42("abs(Mrec-3.1)<0.045");

//--------------------------------------------------------------------
// 2) J/Psi -> phi eta
// Side band for Meta - invariant mass of 2gamma (eta -> 2gamma)
// Window for Mphi - invariant mass of K+ K- (phi -> K+ K-)
//--------------------------------------------------------------------
   TCut c_MCmkk("mcmkk<1.08"); // cut for MC generated events

   // chi^2 cut: NO uncertainties study! (helix correction)
   // double chi2M = 60; // OLD standard: 60; 
   double chi2M = 40; // std: v15a
   // double chi2M = 100; // scan: 20,30,40,50,60,80,100
   string str_c_chi2 = string("ch2<") + to_string(chi2M);
   TCut c_chi2( str_c_chi2.c_str() );


   // Meta: central part
   const double seta = 0.008;
   const double weta = 3*seta; // standard: 3x,
                               // uncertainties study: 2x, 4x
   string str_c_cpgg( Form("abs(Mgg-%.6f)<%.6f",Meta,weta) );
   TCut c_cpgg(str_c_cpgg.c_str());

   // Meta: side-band
   // 'shift_eta' is the start of the side-band
   double shift_eta = 7*seta; // new for prod-11
   // double shift_eta = 6*seta; // old (prod<=10)
   string str_c_sbgg( Form("abs(Mgg-%.6f)>%.6f&&abs(Mgg-%.6f)<%.6f",
            Meta,shift_eta,Meta,shift_eta+weta) );
   TCut c_sbgg(str_c_sbgg.c_str());


   // Mphi cut: [2*Mk, 1.08GeV] -> only for pictures:
   //                              mass_eta.cc
   //                              see eff_mc.cc for efficiency
   string str_c_phi( Form("Mkk>%.6f&&Mkk<1.08",2*Mk) );
   TCut c_phi(str_c_phi.c_str());
//--------------------------------------------------------------------
