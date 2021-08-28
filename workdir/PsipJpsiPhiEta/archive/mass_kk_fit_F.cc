//----------------------------------------------------------------------
vector<ValErr> CombIntErF( const ROOT::Fit::FitResult& res,
      double sl09, double sl12, unsigned int idx) {
//----------------------------------------------------------------------
// Numerical calculation of integrals and their errors
   // desired errors:
   constexpr double eps_abs = 1.e-6;
   constexpr double eps_rel = 1.e-6;
   constexpr double epsilon  = 1.e-3; // max numerical error

   // get parameters and covariation matrix of error
   int Npar = res.NPar();
   vector<double> Fp = res.Parameters();

   // parameters for 2009 and 2012 are considered separately
   int npar = Npar - 2;
   vector<double> Fpar[2];
   Fpar[0] = { Fp[0],Fp[1],Fp[4],Fp[2],Fp[6],Fp[3],sl09 }; //2009
   Fpar[1] = { Fp[0],Fp[1],Fp[5],Fp[2],Fp[7],Fp[3],sl12 }; //2012

   int covMatrStatus = res.CovMatrixStatus();
   if ( covMatrStatus != 3 ) {
      cout << " WARNING: " << __func__ << " covariance matrix"
         " status code is " << covMatrStatus << endl;
   }
   vector<double> cov_m[2];
   cov_m[0].resize(npar*npar,0);
   cov_m[1].resize(npar*npar,0);
   vector<int> I09, I12;
   I09 = {0,1,3,5,2,-1,4,-1};
   I12 = {0,1,3,5,-1,2,-1,4};
   for ( int i = 0; i < Npar; ++i ) {
      int i09 = I09[i];
      int i12 = I12[i];
      for ( int j = 0; j < Npar; ++j ) {
         int j09 = I09[j];
         if ( i09 != -1 && j09 != -1 ) {
            cov_m[0][i09*npar+j09] = res.CovMatrix(i,j);
         }
         int j12 = I12[j];
         if ( i12 != -1 && j12 != -1 ) {
            cov_m[1][i12*npar+j12] = res.CovMatrix(i,j);
         }
      }
   }

   // 1) calculate integrals of normalized component
   //    without correction of efficiency (negative Nidx)
   int Nidx = -idx;
   auto Lfc = [Nidx](const double* x, const double* p) -> double {
      return IntfrBWARGN(x[0],p,Nidx);
   };
   TF1 FC("FC", Lfc, dL, dU, npar+1); // +slope

   double Integ[2] {0.,0.}, err_Int[2] {0.,0.};
   for ( int idat = 0; idat < 2; ++idat ) {

      FC.SetParameters( Fpar[idat].data() );
      double err_num = 0;
      Integ[idat] = FC.IntegralOneDim(dL,dU,eps_rel,eps_abs, err_num);
//       cout << " err_num(" << idat << ") = " << err_num << endl;
      if ( err_num > epsilon ) {
         cout << " WARNING: " << __func__ << " numerical error"
            " of integration of FC-" << idat << " is too big "
            << err_num << endl;
      }

      // 2) calculate error of this integral
      TMatrixDSym covMatrix(npar);
      covMatrix.Use(npar,cov_m[idat].data());
//       printf("\ncovMatrix-%d : ",idat);
//       covMatrix.Print();

      TVectorD IntGrad(npar);
      double err_num2 = 0;
      for ( int i = 0; i < npar; ++i ) {
         // skip parameters with zero error
         if ( covMatrix(i,i) == 0 ) {
            continue;
         }

         auto L_dF = [FC,i](const double* x, const double* p) -> double {
            TF1 tmp(FC); // may modify 'tmp' but not FC
            return tmp.GradientPar(i,x);
         };
         TF1 dF("dF",L_dF,dL,dU,0);
         double err_num = 0;
         IntGrad[i] = dF.IntegralOneDim(dL,dU,eps_rel,eps_abs, err_num);
         err_num = covMatrix(i,i) * IntGrad[i] * err_num;
//          cout << " abs err_num(dF_" << i << ") = " << err_num << endl;
         err_num2 += SQ(err_num);
//          cout << " IntGrad-"<<idat<<"[" << i << "] = "<<IntGrad[i]<<endl;
      }

      err_Int[idat] = sqrt( covMatrix.Similarity(IntGrad) );
      err_num2 = sqrt( err_num2 / err_Int[idat] ); // abs numerical error
//       cout << " err_num(dBW-" << idat << ") = " << err_num2 << endl;
      if ( err_num2 > epsilon ) {
         cout << " WARNING: " << __func__ << " numerical error"
            " of integration of dF-" << idat << " is too big "
            << err_num2 << endl;
      }
   }

   vector<ValErr> ret { ValErr {Integ[0], err_Int[0]},
                        ValErr {Integ[1], err_Int[1]}  };
   return ret;
}

//-------------------------------------------------------------------------
void combine_Intfr_F(string fname09, string fname12, string pdf="") {
//-------------------------------------------------------------------------
   // Get un-binned data
   TH1D* hist[2];
   vector<double> mkk09 = get_mkk_hist( fname09, "mkk_09_cp", &hist[0] );
   vector<double> mkk12 = get_mkk_hist( fname12, "mkk_12_cp", &hist[1] );

   TH1D* h09 = hist[0];
   TH1D* h12 = hist[1];
   double norm09 = h09 -> Integral();
   double norm12 = h12 -> Integral();

   int n09 = mkk09.size();
   int n12 = mkk12.size();
   ROOT::Fit::DataRange dr(dL, dU);
   ROOT::Fit::UnBinData Dat(dr, n09+n12);
   for ( int i = 0; i < n09; ++i ) {
      Dat.Add(-mkk09[i]);
   }
   for ( int i = 0; i < n12; ++i ) {
      Dat.Add(mkk12[i]);
   }

   //-----------------------------------------------------------------------
   // Fit data
   // function MUST be normalized to 1 on the fit range
   const double sl09 = -1.6, sl12 = -1.9;
   auto Lsum = [sl09,sl12](const double* x,const double* p) -> double {
      double m = x[0];
      if ( m < 0 ) { // 2009
         double pp[] = { p[0],p[1],p[4],p[2],p[6],p[3],sl09 };
         return IntfrBWARGN( -m, pp, 0 );
      } else { // 2012
         double pp[] = { p[0],p[1],p[5],p[2],p[7],p[3],sl12 };
         return IntfrBWARGN( m, pp, 0 );
      }
   };
   vector<string> par_name { "M#phi", "G#phi", "A", "#vartheta",
                             "#sigma{09}", "#sigma{12}", "F09", "F12" };
//    vector<double> par_ini {Mphi,Gphi,0.,+0.7,
//                            1.5e-3,1.2e-3,0.8,1.};// neg
   vector<double> par_ini {Mphi,Gphi,0.,-0.8,
                           1.5e-3,1.2e-3,0.8,1.};// pos

   const unsigned int Npar = par_name.size(); // number of parameters
   TF1* fcom = new TF1("fcom", Lsum, dL, dU, Npar);

   ROOT::Fit::Fitter fitter;
   fitter.Config().MinimizerOptions().SetPrintLevel(3);

   ROOT::Math::WrappedTF1 WFun( *fcom );
   fitter.SetFunction( WFun, false); // false == no parameter derivatives

   fitter.Config().SetParamsSettings(Npar,par_ini.data()); // must be first
   for( unsigned int i = 0; i < Npar; ++i ) {
      fitter.Config().ParSettings(i).SetName(par_name[i]);
   }

   // fix/limit/step for parameters
//    fitter.Config().ParSettings(0).SetLimits(Mphi-0.01, Mphi+0.01);
   fitter.Config().ParSettings(0).SetValue(1.01952); // like MC-sig
   fitter.Config().ParSettings(0).Fix();                    // Mphi
   fitter.Config().ParSettings(1).SetLimits(Gphi-0.1e-3, Gphi+0.1e-3);
   fitter.Config().ParSettings(1).Fix();                    // Gphi
   fitter.Config().ParSettings(2).Fix();                    // A
   fitter.Config().ParSettings(3).SetLimits(-M_PI, M_PI);   // vartheta
   fitter.Config().ParSettings(4).SetLimits(0.2e-3, 2.e-3); // sigma09
   fitter.Config().ParSettings(5).SetLimits(0.2e-3, 2.e-3); // sigma12
   fitter.Config().ParSettings(6).SetLimits(0.01, 10.);     // F09
   fitter.Config().ParSettings(7).SetLimits(0.01, 10.);     // F12

   // == Fit
   fitter.LikelihoodFit( Dat, false ); // true == extended likelihood fit
   fitter.CalculateMinosErrors();

   ROOT::Fit::FitResult res = fitter.Result();
   res.Print(cout);
//    res.PrintCovMatrix(cout); // print error matrix and correlations

   double Lmin = res.MinFcnValue();
   ParAndErr PE(res);
   vector<double>& Fpar = PE.Fpar;

   vector<string> names { "N(KK)", "Nphi", "Nnonphi", "Nifr"  };
   vector<ValErr> Nos;
   Nos.reserve(8);
   vector<double> norm { norm09, norm12 };
   for (auto nr : norm ) {
      Nos.push_back( ValErr {nr,sqrt(nr)} );
   }
   printf("%s09 = %s %s12 = %s\n",names[0].c_str(),Nos[0].prt(".1f"),
                                  names[0].c_str(),Nos[1].prt(".1f") );
   for ( int idx = 1; idx <= 3; ++idx ) {
      vector<ValErr> Vinteg = CombIntErF( res,sl09,sl12,idx );
      for ( int j = 0; j < 2 ; ++j ) {
         const auto& Integ = Vinteg[j];
         double Num = norm[j] * Integ.val;
         double Err = fabs(Num) *
            sqrt( 1./norm[j] + SQ(Integ.err/Integ.val) );
         Nos.push_back( ValErr {Num,Err} );
      }
      printf("%s09 = %s %s12 = %s\n",
            names[idx].c_str(),Nos[2*idx].prt(".1f"),
            names[idx].c_str(),Nos[2*idx+1].prt(".1f") );
   }
   double Nphi09 = Nos[2].val, err_Nphi09 = Nos[2].err;
   double Nphi12 = Nos[3].val, err_Nphi12 = Nos[3].err;

   //-----------------------------------------------------------------------
   // "Goodness of fit" using K-S test (see goftest from ROOT-tutorial)
   // User input PDF:
   auto Lgof09 = [Fpar,sl09](double x) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[4],Fpar[2],Fpar[6],Fpar[3],sl09};
      return IntfrBWARGN( x, pp, 0 );
   };
   rmath_fun< decltype(Lgof09) > ftest09(Lgof09);
   ROOT::Math::GoFTest* goftest09 = new ROOT::Math::GoFTest(
         mkk09.size(),mkk09.data(),ftest09,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS09 = goftest09 -> KolmogorovSmirnovTest();
   cout << " pvalueKS09= " << pvalueKS09 << endl;

   auto Lgof12 = [Fpar,sl12](double x) -> double {
      double pp[] = {Fpar[0],Fpar[1],Fpar[5],Fpar[2],Fpar[7],Fpar[3],sl12};
      return IntfrBWARGN( x, pp, 0 );
   };
   rmath_fun< decltype(Lgof12) > ftest12(Lgof12);
   ROOT::Math::GoFTest* goftest12 = new ROOT::Math::GoFTest(
         mkk12.size(),mkk12.data(),ftest12,ROOT::Math::GoFTest::kPDF,dL,dU);
   double pvalueKS12 = goftest12 -> KolmogorovSmirnovTest();
   cout << " pvalueKS12= " << pvalueKS12 << endl;

   //-----------------------------------------------------------------------
   // Functions to draw
   auto Lfit09 = [norm09,Fpar,sl09](double* x,double* p) -> double {
      double m = x[0];
      double pp[] = {Fpar[0],Fpar[1],Fpar[4],Fpar[2],Fpar[6],Fpar[3],sl09};
      return bW*norm09 * IntfrBWARGN( m, pp, int(p[0]) );
   };
   TF1* ff09 = new TF1("ff09", Lfit09, dL, dU, 1);
   ff09 -> SetNpx(500);

   auto Lfit12 = [norm12,Fpar,sl12](double* x,double* p) -> double {
      double m = x[0];
      double pp[] = {Fpar[0],Fpar[1],Fpar[5],Fpar[2],Fpar[7],Fpar[3],sl12};
      return bW*norm12 * IntfrBWARGN( m, pp, int(p[0]) );
   };
   TF1* ff12 = new TF1("ff12", Lfit12, dL, dU, 1);
   ff12 -> SetNpx(500);

   //-----------------------------------------------------------------------
   // Draw results
   TCanvas* c1 = new TCanvas("c1","...",0,0,800,1000);
   c1 -> Divide(1,2);

   c1 -> cd(1);
   gPad -> SetGrid();

   SetHstFace(h09);
   h09 -> GetYaxis() -> SetTitleOffset(1.);
   h09->SetMinimum(-15.);
   h09->SetLineWidth(2);
   h09->SetLineColor(kBlack);
   h09->SetMarkerStyle(20);

   h09->Draw("EP");

   TLegend* leg = new TLegend(0.51,0.45,0.89,0.89);
   leg -> SetHeader("#bf{2009(top)   2012(bottom)}","C");
   leg -> AddEntry( h09,"Data","LEP" );

   ff09 -> SetParameter(0, 0); // Sum
   ff09 -> SetLineWidth(2);
   ff09 -> SetLineColor(kRed+1);
   ff09 -> DrawCopy("SAME");
   leg -> AddEntry( ff09 -> Clone(), "Combined fit", "L" );

   ff09 -> SetParameter(0, 1); // BW
   ff09 -> SetLineWidth(1);
   ff09 -> SetLineStyle(kDashed);
   ff09 -> SetLineColor(kGreen+2);
   ff09 -> DrawCopy("SAME");
   leg -> AddEntry( ff09 -> Clone(), "Breit-Wigner", "L");

   ff09 -> SetParameter(0, 2); // Argus
   ff09 -> SetLineColor(kBlue);
   ff09 -> DrawCopy("SAME");
   leg -> AddEntry( ff09 -> Clone(), "Argus", "L");

   ff09 -> SetParameter(0, 3); // interference
   ff09 -> SetLineColor(kMagenta+1);
   ff09 -> DrawCopy("SAME");
   leg -> AddEntry( ff09 -> Clone(), "Interference", "L");
   leg -> Draw();
   gPad -> RedrawAxis();

   c1 -> cd(2);
   gPad -> SetGrid();
   SetHstFace(h12);
   h12 -> GetYaxis() -> SetTitleOffset(1.);
   h12 -> SetMinimum(-35.);
   h12 -> SetLineWidth(2);
   h12 -> SetLineColor(kBlack);
   h12 -> SetMarkerStyle(20);

   h12 -> Draw("EP");

   ff12 -> SetParameter(0, 0); // Sum
   ff12 -> SetLineWidth(2);
   ff12 -> SetLineColor(kRed+1);
   ff12 -> DrawCopy("SAME");

   ff12 -> SetParameter(0, 1); // BW
   ff12 -> SetLineWidth(1);
   ff12 -> SetLineStyle(kDashed);
   ff12 -> SetLineColor(kGreen+2);
   ff12 -> DrawCopy("SAME");

   ff12 -> SetParameter(0, 2); // Argus
   ff12 -> SetLineColor(kBlue);
   ff12 -> DrawCopy("SAME");

   ff12 -> SetParameter(0, 3); // interference
   ff12 -> SetLineColor(kMagenta+1);
   ff12 -> DrawCopy("SAME");

   TPaveText* pt09 = new TPaveText(0.51,0.78,0.89,0.99,"NDC");
   pt09 -> SetTextAlign(12);
   pt09 -> SetTextFont(42);
   pt09 -> AddText( Form("#it{p-value(2009) = %.3f}", pvalueKS09) );
   pt09 -> AddText( Form("N_{#phi}(2009) = %.1f #pm %.1f",
                    Nphi09,err_Nphi09) );
   pt09 -> AddText( Form("#sigma(2009)= %s MeV", PE.Eform(4,".2f",1e3)) );
   pt09 -> AddText( Form("F(2009)= %s",PE.Eform(6,".2f")) );
   pt09 -> Draw();

   TPaveText* pt12 = new TPaveText(0.51,0.57,0.89,0.78,"NDC");
   pt12 -> SetTextAlign(12);
   pt12 -> SetTextFont(42);
   pt12 -> AddText( Form("#it{p-value(2012) = %.3f}", pvalueKS12) );
   pt12 -> AddText( Form("N_{#phi}(2012) = %.1f #pm %.1f",
                    Nphi12,err_Nphi12) );
   pt12 -> AddText( Form("#sigma(2012)= %s MeV", PE.Eform(5,".2f",1e3)) );
   pt12 -> AddText( Form("F(2012)= %s",PE.Eform(7,".2f")) );
   pt12 -> Draw();

   TPaveText* pt = new TPaveText(0.51,0.36,0.89,0.57,"NDC");
   pt -> SetTextAlign(12);
   pt -> SetTextFont(42);
//    pt -> AddText( Form("#it{L_{min} = %.1f}",Lmin) );
   pt -> AddText( Form("M_{#phi}= %s MeV",PE.Eform(0,".2f",1e3)) );
   pt -> AddText( Form("#Gamma_{#phi}= %s MeV",PE.Eform(1,".3f",1e3)) );
   pt -> AddText( Form("a= %s",PE.Eform(2,".2f")) );
   pt -> AddText( Form("#vartheta= %s",PE.Eform(3,".2f")) );
   pt -> Draw();

   gPad -> RedrawAxis();
   c1 -> Update();
   if ( !pdf.empty() ) {
      c1 -> Print(pdf.c_str());
   }
}

