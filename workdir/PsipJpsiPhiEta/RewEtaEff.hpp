// RewEtaEff.hpp - function with correction for efficiency of
// reconstruction Eta for MC events

//--------------------------------------------------------------------
double RewEtaEff( int DataPeriod ) {
//--------------------------------------------------------------------
   // v709, DelClonedTrk, helix corrections 
   // J/Psi -> gamma eta (+) J/Psi -> phi eta
   double W = 0.984; //0.9837 +/- 0.0018
   return W;
}
