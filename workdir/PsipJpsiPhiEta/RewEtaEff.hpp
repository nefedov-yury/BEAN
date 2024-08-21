// RewEtaEff.hpp
// RewEtaEff() returns a weight to correct the reconstruction
// efficiency eta->2gamma for an MC event

//--------------------------------------------------------------------
double RewEtaEff( int DataPeriod ) {
//--------------------------------------------------------------------
   // v709n4eff: J/Psi -> gamma eta (+) J/Psi -> phi eta
   double W = 0.994; // 0.994 +/- 0.004 (error scaled 1.7 times)
   return W;
}
