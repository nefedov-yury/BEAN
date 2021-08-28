// List of constants and cuts
//-------------------------------------------------------------------------
// 1) Psi' -> J/Psi pi+ pi-
//-------------------------------------------------------------------------
#ifdef CONSTANTS_ONLY
   const double Cto09 =  3.276; // +/- 0.021
   const double Cto12 = 10.273; // +/- 0.067
   const double Ito09 =  1.022; // +/- 0.008
   const double Ito12 =  0.853; // +/- 0.005
#else
   TCut c_MCtrue("dec==64"); // only Psi' -> J/Psi pi+ pi- decay

   // cut on recoil mass: [3.092, 3.102]
   TCut c_Mrec("abs(Mrec-3.097)<0.005");

   // cut on recoil mass the same as in BAM-0042 [3.055, 3.145]
   TCut c_MrecBAM42("abs(Mrec-3.1)<0.045");

//-------------------------------------------------------------------------
// 2) J/Psi -> phi eta
// Side band for Meta - invariant mass of 2gamma (eta -> 2gamma)
// Window for Mphi - invariant mass of K+ K- (phi -> K+ K-)
//-------------------------------------------------------------------------
   TCut c_MCmkk("mcmkk<1.08"); // cut for MC generated events

   // chi^2 cut:
   double chi2M = 80;  // standard
//    double chi2M = 40;  // uncertainties study: 60, 100
   string str_c_chi2 = string("ch2<") + to_string(chi2M);
   TCut c_chi2( str_c_chi2.c_str() );


   // Meta: central part
   const double meta = 0.547862; //  547.862 +/- 0.017 MeV
   const double seta = 0.008;
   const double weta = 3*seta; // standard
//    const double weta = 4*seta; // uncertainties study: 2x, 4x
   string str_c_cpgg = string("abs(Mgg-") + to_string(meta) + string(")<")
                     + to_string(weta);
   TCut c_cpgg(str_c_cpgg.c_str());

   // Meta: side-band
   // actually 'shft' is the gap between signal region and side-band
   double shft_eta = 6*seta - weta;
   string str_c_sbgg = string("abs(Mgg-") + to_string(meta) + string(")>")
                     + to_string(weta+shft_eta)
                     + string("&&")
                     + string("abs(Mgg-") + to_string(meta) + string(")<")
                     + to_string(2*weta+shft_eta);
   TCut c_sbgg(str_c_sbgg.c_str());


   // Mphi cut: [2*mk, 1.08GeV] -> only for pictures !
   //                              see eff_mc.cc for efficiency
   // const double mphi = 1.019461; // 1019.461 +/- 0.019 MeV
   const double mk   = 0.493677; // 493.677  +/- 0.016 MeV
   string str_c_phi = string("Mkk>") + to_string(2*mk)
                    + string("&&") + string("Mkk<1.08");
   TCut c_phi(str_c_phi.c_str());
//-------------------------------------------------------------------------
#endif
