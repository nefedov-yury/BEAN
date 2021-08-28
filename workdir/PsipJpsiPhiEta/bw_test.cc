// {{{1 Breigt Wigner for phi -> KK
//----------------------------------------------------------------------
// Breigt Wigner for phi -> KK
//----------------------------------------------------------------------

//----------------------------------------------------------------------
constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

//----------------------------------------------------------------------
double BreitWigner(double m, double mphi, double gphi) {
//----------------------------------------------------------------------
// see BAM-00117:  PhysRevD.91.112001.pdf (arXiv:1504.03194v2)
//
// we assume that mphi and gphi are the fit parameters

   static const double Mjpsi = 3.096916; // 3096.916 +/- 0.011 MeV
   static const double Meta  = 0.547862; // 547.862  +/- 0.017 MeV
   static const double Mk    = 0.493677; // 493.677  +/- 0.016 MeV

   static const double Mjpsi2 = SQ(Mjpsi);
   static const double Meta2  = SQ(Meta);
   static const double Mk2    = SQ(Mk);

   static const double R      = 3.; // GeV^{-1} for Blatt-Weisskopf ff

   static const double dL = 2*Mk; // the left cutoff = 0.987354

   if ( m < dL ) { // phase space becomes zero (see p_k)
      return 0;
   }

   double m2    = SQ(m);
   double mphi2 = SQ(mphi);
   double gphi2 = SQ(gphi);

   double gam = mphi*sqrt(mphi2+gphi2);
   double kap = 2*M_SQRT2*mphi*gphi*gam / (M_PI*sqrt(mphi2+gam));

   double p_phi0 = sqrt(SQ(Mjpsi2-mphi2-Meta2)-4*mphi2*Meta2) / (2*Mjpsi);
   double p_phi  = sqrt(SQ(Mjpsi2-m2-Meta2)-4*m2*Meta2) / (2*Mjpsi);
   double r_phi  = p_phi / p_phi0;

   double p_k0 = 0.5*sqrt(mphi2-4*Mk2);
   double p_k  = 0.5*sqrt(m2-4*Mk2);
   double r_k  = p_k / p_k0;

   // B(m)/B(m0)
   double BB_phi = sqrt( (1+SQ(R*p_phi0)) / (1+SQ(R*p_phi)) );
   double BB_k2  = (1+SQ(R*p_k0)) / (1+SQ(R*p_k));
   double BB_k = sqrt( BB_k2 );

   double fm = r_phi*r_k*BB_phi*BB_k;

   double GM = gphi*mphi*(r_k*r_k*r_k)*BB_k2; // == m*G(m)

   return (kap * SQ(fm)) / ( SQ(m2-mphi2) + SQ(GM) );
}

//-------------------------------------------------------------------------
void test_BreitWigner() {
//-------------------------------------------------------------------------

   static const double Mphi  = 1.019461; //1019.461  +/- 0.019 MeV
   static const double Gphi  = 4.247e-3; //   4.247  +/- 0.016 MeV
   static const double Mk    = 0.493677; // 493.677  +/- 0.016 MeV
   static const double dL = 2*Mk; // the left cutoff = 0.987354
   const double dU = 1.08;

   auto Lbw = [](const double* x,const double* p) -> double {
      return p[0]*BreitWigner(x[0],p[1],p[2]);
   };

   TF1* bw = new TF1("bw", Lbw, 0.98, dU, 3);
   bw->SetParNames("Norm","Mphi","Gphi");
   bw->SetParameters(1., Mphi, Gphi);
   bw->SetLineWidth(1);
   bw->SetLineColor(kBlue);
   bw->SetNpx(500);

   double norm = 1./bw->Integral(dL,dU,1e-8);
   printf("norm = %.6f\n",norm);
   bw->SetParameter(0, norm );

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad->SetGrid();
   bw->Draw();
}

// {{{1 MAIN:
//-------------------------------------------------------------------------
void bw_test() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);

   test_BreitWigner();
}
