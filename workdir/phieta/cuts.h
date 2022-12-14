// List of the cuts
// ATTENSION: to_string() converts only 6 digits after the point.
//            It's better to use the Form() function.

#include "masses.h"

//--------------------------------------------------------------------
// J/Psi -> phi eta
// Side band for Meta: invariant mass of 2gamma (eta -> 2gamma)
// Window for Mphi: invariant mass of K+ K- (phi -> K+ K-)
//--------------------------------------------------------------------
   // signal definition for MC: X_isr < 0.9
   double Xisr_min = 0.9;
   TCut c_xisr( Form("xisr>%f",Xisr_min) );
   auto f_xisr = [Xisr_min](double x)->bool{return x>Xisr_min;};

   TCut c_MCmkk("mcmkk<1.08"); // cut for MC generated events
   auto f_MCmkk = [](double mcmkk)->bool{return mcmkk<1.08;};

//--------------------------------------------------------------------
   // chi^2 cut:
   double chi2M = 80;  // std: 40 in Psi(2S)
   TCut c_chi2( Form("ch2<%f",chi2M) );
   auto f_chi2 = [chi2M](double ch2)->bool{return ch2<chi2M;};

   // Meta; central part
   const double seta = 0.008;
   const double weta = 3*seta; // standard
                               // uncertainties study: 2x, 4x
   TCut c_cpgg( Form("abs(Mgg-%.6f)<%.6f",Meta,weta) );
   auto f_cpgg = [=](double Mgg)->bool{return fabs(Mgg-Meta)<weta;};

   // Meta: side-band
   // 'shift_eta' is the start of the side-band
   double shift_eta = 6*seta; // 7 in Psi(2S)
   TCut c_sbgg( Form("abs(Mgg-%.6f)>%.6f&&abs(Mgg-%.6f)<%.6f",
            Meta,shift_eta,Meta,shift_eta+weta) );
   auto f_sbgg = [=](double Mgg) -> bool{return
      (fabs(Mgg-Meta)>shift_eta && fabs(Mgg-Meta)<shift_eta+weta);};

   // Mphi cut: [2*Mk, 1.08GeV] 
   TCut c_phi( Form("Mkk>%.6f&&Mkk<1.08",2*Mk) );
   auto f_phi = [=](double Mkk)->bool{return(Mkk>2*Mk && Mkk<1.08);};
//--------------------------------------------------------------------
