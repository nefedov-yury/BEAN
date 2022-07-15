// {{{1

//----------------------------------------------------------------------
constexpr double SQ(double x) {
//----------------------------------------------------------------------
   return x*x;
}

//----------------------------------------------------------------------
TGraph* GrL() {
//----------------------------------------------------------------------
// #include "prod-12/Optim/Lmedian_br_scan.h"

// #include "prod-12/Optim/Lmedian_ar1_ch20.h"
// #include "prod-12/Optim/Lmedian_ar1_ch30.h"
// #include "prod-12/Optim/Lmedian_ar1_ch40.h"
// #include "prod-12/Optim/Lmedian_ar1_ch50.h"
// #include "prod-12/Optim/Lmedian_ar1_ch60.h"
// #include "prod-12/Optim/Lmedian_ar2_ch60.h"
// #include "prod-12/Optim/Lmedian_ar1_ch80.h"
// #include "prod-12/Optim/Lmedian_ar1_ch100.h"
// #include "prod-12/Optim/Lmedian_ar2_ch100.h"
// #include "prod-12/Optim/Lmedian_ar2_ch150.h"

#include "prod-12/Optim/Lmedian_br_ar1_ch40_final.h"
#include "prod-12/Optim/Lmedian_ar1_ch40_final.h"

   for( auto &b : br_phi ) {
      b *= 1e4;
   }

   // calculate rL = L / Lmax:
   vector<double> rL;
   double minL = *(min_element(Lmin.begin(),Lmin.end()));
   for( auto &l : Lmin ) {
      l -= minL;
      rL.push_back(exp(-0.5*l));
   }

   TGraph* grL = new TGraph( rL.size(),br_phi.data(),rL.data() );
   grL -> SetTitle(";Br(#phi#eta) #times 10^{4};#it{L/L_{max}}");
   grL -> GetYaxis() -> SetMaxDigits(3);
   grL -> GetYaxis() -> SetTitleOffset(1.);
   grL -> SetMarkerColor(kBlue);
   grL -> SetMarkerStyle(20);
   grL -> SetMaximum(1.1);
   grL -> SetMinimum(0.0);

   return grL;
}

//----------------------------------------------------------------------
double Median(const TGraph* gr) {
//----------------------------------------------------------------------
   int n = gr -> GetN();
   const double* x = gr -> GetX();
   const double* w = gr -> GetY();
   return TMath::Median(n, x, w);
}

// return median M and confidence interval [L,R] around M
// such that prob[LM] = prob[MR] = CL/2
//----------------------------------------------------------------------
vector<double> MConfInt(const TGraph* gr, double CL) {
//----------------------------------------------------------------------
   int n = gr -> GetN();
   const double* x = gr -> GetX();
   const double* y = gr -> GetY();
   vector<double> w(y,y+n);
   partial_sum(begin(w), end(w), begin(w));
   // normalization to one
   auto sw = w.back();
   for (auto& v : w ) {
      v /= sw;
   }

   // find median (0.5), left (0.5-cl) and right (0.5+cl)
   // ret(3+i) is value of gr->y() for ret(i)
   vector<double> p {0.5,0.5*(1-CL),0.5*(1+CL)};
   vector<double> ret(6,0);
   for( int i = 0; i < 3; i++ ) {
      auto pU = upper_bound(begin(w),end(w),p[i]);
      // linear extrapolation
      int j = distance(begin(w), pU);
      if (j > 0 && j < n ) {
         double dw = w[j] - w[j-1];
         double dx = x[j] - x[j-1];
         double dy = y[j] - y[j-1];
         ret[i]   = x[j-1] + dx*(p[i]-w[j-1])/dw;
         ret[3+i] = y[j-1] + dy*(p[i]-w[j-1])/dw;
      }
   }

   return ret;
}

// return median M and symmatric confidence interval [L,R] around M
// such that |LM| = |MR| and prob[LR] = CL
//----------------------------------------------------------------------
vector<double> MSymConfInt(const TGraph* gr, double CL) {
//----------------------------------------------------------------------
   int n = gr -> GetN();
   const double* x = gr -> GetX();
   const double* y = gr -> GetY();
   vector<double> w(y,y+n);
   partial_sum(begin(w), end(w), begin(w));
   // normalization to one
   auto sw = w.back();
   for (auto& v : w ) {
      v /= sw;
   }

   // initial estimation:
   double cl = 0.5*CL;
   double cr = 0.5*CL;
   // find median (0.5), left (0.5-cl) and right (0.5+cr)
   // ret(3+i) is value of gr->y() for ret(i)
   vector<double> ret(6,0);
   for ( int iter = 0; iter < 10; ++iter ) {
      vector<double> p {0.5,0.5-cl,0.5+cr};
      for( int i = 1 - (iter==0); i < 3; i++ ) {
         auto pU = upper_bound(begin(w),end(w),p[i]);
         // linear extrapolation
         int j = distance(begin(w), pU);
         if (j > 0 && j < n ) {
            double dw = w[j] - w[j-1];
            double dx = x[j] - x[j-1];
            double dy = y[j] - y[j-1];
            ret[i]   = x[j-1] + dx*(p[i]-w[j-1])/dw;
            ret[3+i] = y[j-1] + dy*(p[i]-w[j-1])/dw;
         }
      }
      // check that |LM| = |MR|
      double LM = ret[0]-ret[1];
      double MR = ret[2]-ret[0];
      if ( fabs(LM-MR) < 1e-3 ) {
         break;
      } else if ( LM > MR ) {
         cr *= 1 + 0.5*(LM-MR)/(LM+MR);
         cl = CL-cr;
      } else { // MR > LM
         cl *= 1 + 0.5*(MR-LM)/(LM+MR);
         cr = CL-cl;
      }
      // debug prints:
      printf("median= %.4f in [%.4f,%.4f]\n",ret[0],ret[1],ret[2]);
      printf("new cl=%.4f new cr=%.4f\n",cl,cr);
   }
   return ret;
}

// return most probable value V and confidence interval [L,R] around V
// such that prob[LV] = prob[VR] = CL/2
//----------------------------------------------------------------------
vector<double> VConfInt(const TGraph* gr, double CL) {
//----------------------------------------------------------------------
   int n = gr -> GetN();
   const double* x = gr -> GetX();
   const double* y = gr -> GetY();
   vector<double> w(y,y+n);
   partial_sum(begin(w), end(w), begin(w));
   // normalization to one
   auto sw = w.back();
   for (auto& v : w ) {
      v /= sw;
   }

   double my = 0; // max y
   double mx = 0;
   double mw = 0;
   for ( int i = 0; i < n; ++i ) {
      if (  y[i] > my ) {
         my = y[i];
         mx = x[i];
         mw = w[i];
      }
   }
   vector<double> ret {mx,0,0,my,0,0};
   printf("x-max= %.3f and y-max=%.3f,w-max=%.3f\n",mx,my,mw);

   // left (mw-cl) and right (mw+cl)
//    vector<double> p {mw-0.5*CL,mw+0.5*CL};
   vector<double> p {0.5*(1-CL),0.5*(1+CL)};
   for( int i = 0; i < 2; i++ ) {
      auto pU = upper_bound(begin(w),end(w),p[i]);
      // linear extrapolation
      int j = distance(begin(w), pU);
      if (j > 0 && j < n ) {
         double dw = w[j] - w[j-1];
         double dx = x[j] - x[j-1];
         double dy = y[j] - y[j-1];
         ret[1+i] = x[j-1] + dx*(p[i]-w[j-1])/dw;
         ret[4+i] = y[j-1] + dy*(p[i]-w[j-1])/dw;
      }
   }

   return ret;
}

// {{{1 MAIN:
//-------------------------------------------------------------------------
void Lmedian() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);

   TGraph* gr = GrL();
//    double Bmed = Median( gr );
   vector<double> Mcl = MConfInt(gr,0.69);
//    vector<double> Mcl = MSymConfInt(gr,0.683);
   printf("Br_median= %.3f and cl=[%.3f,%.3f]\n",Mcl[0],Mcl[1],Mcl[2]);

//    vector<double> Mcl = VConfInt(gr,0.69);
//    printf("Br_max= %.3f and cl=[%.3f,%.3f]\n",Mcl[0],Mcl[1],Mcl[2]);

   TLatex* tt[3];
   TLine* ll[3];
   for ( int i = 0; i < 3; i++ ) {
      auto Tcol = (i==0) ? kGreen+3 : kRed+1;
      auto& ti = tt[i];
      ti = new TLatex(Mcl[i]-0.07, 0.1, Form("%.2f",Mcl[i]));
      ti -> SetTextSize(0.03);
      ti -> SetTextAlign(12); // centered
      ti -> SetTextAngle(90);
      ti -> SetTextColor(Tcol);

      auto& li = ll[i];
      li = new TLine;
      li -> SetLineColor(Tcol);
      li -> SetLineWidth(2);
      li -> SetLineStyle((i==0) ? 7 : 9);
   }


   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1 -> cd();
   gPad -> SetGrid();

   gr -> Draw("APL");

   for ( int i = 0; i < 3; i++ ) {
      ll[i] -> DrawLine(Mcl[i],0.,Mcl[i],Mcl[3+i]);
      tt[i] -> Draw();
   }

   gPad -> RedrawAxis();
   c1 -> Update();
}
