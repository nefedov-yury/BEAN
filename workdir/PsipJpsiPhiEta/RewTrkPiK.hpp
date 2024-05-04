// RewTrkPiK.hpp - functions with corrections for the efficiency of
// Kaon/Pion reconstruction,including PID, for the MC events

// {{{1 RewTrkPi()
//--------------------------------------------------------------------
double RewTrkPi(int DataPeriod, double Pt, double Z) {
//--------------------------------------------------------------------
   // Corrections for the efficiency of reconstruction a pion having
   // transverse momentum Pt and sign Z.
   // The return value is the weight for the MC event.
   // v709, DelClonedTrk, helix corrections
   // Note: This function works well for MC without helix corrections

   const double Ptmin = 0.05, Ptmax = 0.4;
   Pt = max( Ptmin, Pt );
   Pt = min( Ptmax, Pt );

   // 'spline-function' rewritten from ROOT-generated function
   auto Lsp = [](double x, const double* fY,
         const double* fB, const double* fC, const double* fD) {
      const int fNp = 7;
      const double fDelta = 0.050, fXmin = 0.075, fXmax = 0.375;
      // const double fX[7] = {
         // 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375 };

      // Equidistant knots
      int klow = int( (x-fXmin)/fDelta );
      klow = max(klow,0);
      klow = min(klow,fNp-1);

      // Evaluate now
      // double dx = x - fX[klow];
      double dx = x - (fXmin + klow * fDelta);
      return fY[klow] + dx*(fB[klow] + dx*(fC[klow] + dx*fD[klow]));
   };

   double W = 1.;
   if ( DataPeriod == 2021 ) {
      if ( Z > 0 ) {
         static const double fY[7] = {
            0.981691, 0.987616, 0.985132, 0.978823, 0.981374,
            0.984894, 0.983846 };
         static const double fB[7] = {
            0.154886, 0.0457408, -0.131388, -0.0478018, 0.0971147,
            0.0236108, -0.0432394 };
         static const double fC[7] = {
            -5.55112e-16, -2.1829, -1.35968, 3.03141, -0.133076,
            -1.337, 0.05 };
         static const double fD[7] = {
            -14.5527, 5.48818, 29.2739, -21.0965, -8.02619,
            8.91336, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      } else {
         static const double fY[7] = {
            0.986859, 0.987314, 0.985509, 0.989826, 0.987409,
            0.988992, 0.981344 };
         static const double fB[7] = {
            0.0332347, -0.0391776, 0.0424775, 0.0199988, -0.00851865,
            -0.0359787, -0.211467 };
         static const double fC[7] = {
            6.93889e-17, -1.44825, 3.08135, -3.53092, 2.96057,
            -3.50977, 0.05 };
         static const double fD[7] = {
            -9.65498, 30.1973, -44.0818, 43.2766, -43.1357,
            23.3985, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      }
   } else if ( DataPeriod == 2012 ) {
      if ( Z > 0 ) {
         static const double fY[7] = {
            1.03169, 0.993442, 0.992989, 0.980393, 0.984639,
            0.982477, 0.987807 };
         static const double fB[7] = {
            -0.99227, -0.310335, -0.0883914, -0.119001, 0.0633865,
            -0.00948767, 0.164633 };
         static const double fC[7] = {
            4.44089e-15, 13.6387, -9.19985, 8.58765, -4.93989,
            3.48241, 0.05 };
         static const double fD[7] = {
            90.9248, -152.257, 118.583, -90.1836, 56.1487,
            -23.216, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      } else {
         static const double fY[7] = {
            1.03121, 0.995748, 0.994734, 0.998636, 0.985478,
            0.993828, 0.988655 };
         static const double fB[7] = {
            -0.877667, -0.372377, 0.178602, -0.16878, -0.0588284,
            0.115605, -0.212984 };
         static const double fC[7] = {
            -2.22045e-15, 10.1058, 0.9138, -7.86144, 10.0605,
            -6.57179, 0.05 };
         static const double fD[7] = {
            67.3719, -61.2799, -58.5016, 119.479, -110.882,
            43.8119, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      }
   } else if ( DataPeriod == 2009 ) {
      if ( Z > 0 ) {
         static const double fY[7] = {
            0.992996, 0.99419, 0.990618, 0.982758, 0.986737,
            0.981416, 0.995198 };
         static const double fB[7] = {
            0.0372614, -0.00290172, -0.168377, -0.00952206, -0.0263622,
            0.0344359, 0.396241 };
         static const double fC[7] = {
            -2.77556e-16, -0.803263, -2.50624, 5.68333, -6.02013,
            7.2361, 0.05 };
         static const double fD[7] = {
            -5.35509, -11.3532, 54.5971, -78.0231, 88.3749,
            -48.2407, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      } else {
         static const double fY[7] = {
            1.00415, 0.975639, 0.983931, 0.995963, 0.99217,
            0.991772, 1.00196 };
         static const double fB[7] = {
            -0.755805, -0.198795, 0.338097, 0.0658337, -0.107085,
            0.11107, 0.249998 };
         static const double fC[7] = {
            2.22045e-15, 11.1402, -0.402369, -5.04289, 1.58452,
            2.77857, 0.05 };
         static const double fD[7] = {
            74.268, -76.9504, -30.9368, 44.1827, 7.96038,
            -18.5238, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      }
   }
   return W;
}

// {{{1 RewTrk_K()
//--------------------------------------------------------------------
double RewTrk_K(int DataPeriod, double Pt, double Z) {
//--------------------------------------------------------------------
   // Corrections for the efficiency of reconstruction a kaon having
   // transverse momentum Pt and sign Z.
   // The return value is the weight for the MC event.
   // v709, DelClonedTrk, helix corrections

   const double Ptmin = 0.1, Ptmax = 1.4;
   Pt = max( Ptmin, Pt );
   Pt = min( Ptmax, Pt );

   // 'spline-function' rewritten from ROOT-generated function
   auto Lsp = [](double x, const double* fY,
         const double* fB, const double* fC, const double* fD) {
      const int fNp = 13;
      const double fDelta = 0.1, fXmin = 0.15, fXmax = 1.35;
      // const double fX[13] = {
         // 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05,
         // 1.15, 1.25, 1.35 };

      // Equidistant knots
      int klow = int( (x-fXmin)/fDelta );
      klow = max(klow,0);
      klow = min(klow,fNp-1);

      // Evaluate now
      // double dx = x - fX[klow];
      double dx = x - (fXmin + klow * fDelta);
      return fY[klow] + dx*(fB[klow] + dx*(fC[klow] + dx*fD[klow]));
   };

   double W = 1.;
   if ( DataPeriod == 2021 ) {
      if ( Z > 0 ) {
         static const double fY[13] = {
            0.953415, 0.966404, 0.963077, 0.986586, 0.983971,
            0.982765, 0.978689, 0.975506, 0.971573, 0.970367,
            0.978398, 0.978526, 0.98026 };
         static const double fB[13] = {
            0.198006, -0.00636688, 0.117306, 0.142606, -0.0609012,
            -0.0136211, -0.0430703, -0.0318821, -0.0429022, 0.0493355,
            0.0503176, -0.00584932, 0.0289374 };
         static const double fC[13] = {
            -2.77556e-16, -2.04373, 3.28046, -3.02746, 0.992394,
            -0.519593, 0.225101, -0.113219, 0.0030179, 0.919359,
            -0.909537, 0.347867, 0.1 };
         static const double fD[13] = {
            -6.81243, 17.7473, -21.0264, 13.3995, -5.03996,
            2.48231, -1.12774, 0.387457, 3.05447, -6.09632,
            4.19135, -1.15956, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      } else {
         static const double fY[13] = {
            1.00447, 0.996427, 0.967089, 0.985719, 0.989781,
            0.986789, 0.983506, 0.978769, 0.974528, 0.96308,
            0.967628, 0.976194, 0.963668 };
         static const double fB[13] = {
            0.0135562, -0.268285, -0.0617527, 0.194046, -0.0336553,
            -0.0273313, -0.0452709, -0.0321679, -0.0954092,-0.0568715,
            0.115914, -0.0133839, -0.181196 };
         static const double fC[13] = {
            1.38778e-16, -2.81841, 4.88374, -2.32575, 0.0487428,
            0.0144974, -0.193893, 0.324924, -0.957338, 1.34272,
            0.385141, -1.67812, 0.1 };
         static const double fD[13] = {
            -9.39471, 25.6738, -24.0316, 7.91499, -0.114151,
            -0.694636, 1.72939, -4.2742, 7.66684, -3.19192,
            -6.87754, 5.59374, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      }
   } else if ( DataPeriod == 2012 ) {
      if ( Z > 0 ) {
         static const double fY[13] = {
            1.0137, 0.968143, 0.967563, 0.991428, 0.996475,
            0.982962, 0.984108, 0.978996, 0.97277, 0.97812,
            0.967541, 0.974272, 0.961904 };
         static const double fB[13] = {
            -0.556072, -0.254524, 0.190091, 0.192714, -0.0935853,
            -0.072368, 0.0120715, -0.0948862, 0.0273266, -0.0407048,
            -0.0213846, 0.010797, -0.190912 };
         static const double fC[13] = {
            2.22045e-15, 3.01547, 1.43068, -1.40445, -1.45854,
            1.67071, -0.826315, -0.243261, 1.46539, -2.1457,
            2.33891, -2.01709, 0.1 };
         static const double fD[13] = {
            10.0516, -5.28265, -9.45044, -0.180282, 10.4308,
            -8.32342, 1.94351, 5.6955, -12.037, 14.9487,
            -14.52, 6.72363, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      } else {
         static const double fY[13] = {
            0.990743, 1.00858, 0.973257, 0.995647, 0.990059,
            0.993349, 0.983223, 0.985967, 0.976382, 0.971146,
            0.981296, 0.967797, 0.999896 };
         static const double fB[13] = {
            0.368389, -0.201565, -0.0867184, 0.160344, -0.0505934,
            -0.0268997, -0.0468839, -0.00703707, -0.130211, 0.0832455,
            -0.0553451, 0.0376767, 0.462645 };
         static const double fC[13] = {
            2.77556e-16, -5.69953, 6.84799, -4.37737, 2.26799,
            -2.03106, 1.83121, -1.43275, 0.201006, 1.93356,
            -3.31946, 4.24968, 0.1 };
         static const double fD[13] = {
            -18.9984, 41.8251, -37.4179, 22.1512, -14.3302,
            12.8742, -10.8799, 5.44584, 5.77518, -17.5101,
            25.2305, -14.1656, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      }
   } else if ( DataPeriod == 2009 ) {
      if ( Z > 0 ) {
         static const double fY[13] = {
            0.920735, 1.00077, 0.978848, 1.00663, 0.99453,
            1.00007, 0.980407, 0.994184, 0.987535, 0.97732,
            0.975571, 0.986549, 0.999144 };
         static const double fB[13] = {
            1.11824, 0.164452, -0.0326887, 0.142158, -0.0654663,
            -0.0771325, -0.0496844, 0.0993571, -0.133917, -0.0696089,
            0.0534208, 0.132816, 0.122511 };
         static const double fC[13] = {
            0, -9.53792, 7.56652, -5.81805, 3.74181,
            -3.85847, 4.13295, -2.64253, 0.309792, 0.33329,
            0.897007, -0.10305, 0.1 };
         static const double fD[13] = {
            -31.7931, 57.0148, -44.6152, 31.8662, -25.3343,
            26.6381, -22.5849, 9.84109, 0.0783255, 1.87906,
            -3.33352, 0.3435, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      } else {
         static const double fY[13] = {
            0.924424, 1.01926, 1.00608, 1.01769, 0.998392,
            0.989184, 0.988635, 0.981702, 0.990135, 0.977559,
            0.981579, 1.00003, 0.996868 };
         static const double fB[13] = {
            1.26186, 0.321336, -0.0974525, 0.0212766, -0.218374,
            -0.00284342, -0.0629393, 0.0301609, -0.0127175, -0.103602,
            0.17044, 0.0958909, -0.0953322 };
         static const double fC[13] = {
            2.22045e-15, -9.40519, 5.21731, -4.03002, 1.63351,
            0.521791, -1.12275, 2.05375, -2.48254, 1.57369,
            1.16674, -1.91223, 0.1 };
         static const double fD[13] = {
            -31.3506, 48.7417, -30.8244, 18.8784, -3.70574,
            -5.4818, 10.5883, -15.121, 13.5208, -1.3565,
            -10.2632, 6.3741, 1.73205 };
         W = Lsp( Pt,fY, fB,fC,fD );
      }
   }
   return W;
}

