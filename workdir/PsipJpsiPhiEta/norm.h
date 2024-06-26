// horm.h - normalization constants for inclusive MC and non-resonance
// data
// BOSS v7.0.9:
// http://docbes3.ihep.ac.cn/cgi-bin/DocDB/ShowDocument?docid=1204
// Inclusive MC:  2009 107M;  2012 341M;  2021 2300M;
// Continuum: we use only 2021 401pb-1
//-------------------------------------------------------------------
//                           Luminosity    coefficient  coefficient
//        N(ψ(3686)) δN(%)        (pb−1)   Inclusive MC  continuum
//-------------------------------------------------------------------
// 2009  107.7×10^6   0.60   161.63±0.13   1.007±0.006     0.40
// 2012  345.42×10^6  0.72   506.92±0.23   1.013±0.007     1.24
// 2021 2260.0×10^6   0.49   3208.5±1.9    0.982±0.005     7.85
//-------------------------------------------------------------------
const map<int,double> MC_Dat
   { {2009,1.007}, {2012,1.013}, {2021,0.982} };
const map<int,double> C_Dat
   { {2009,0.40},  {2012,1.24},  {2021,7.85} };
//-------------------------------------------------------------------
