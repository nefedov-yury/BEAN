// fit M_{K+K-} (RooFit version)

#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"

// for cuts.h
#include "TCut.h"
#include "TCutG.h"

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "Constants.h"
#include "RevArgus.h"
#include "BreitWignerPhi.h"
#include "Interference.h"

#include<iostream>
#include<string>

using namespace RooFit;
using namespace std;

//-------------------------------------------------------------------------
TH1D* plot_Mkk(string fname, string hname, int type = 0) {
//-------------------------------------------------------------------------
#include "../cuts.h"

   gStyle->SetOptStat(0);

   // name of folder with root files
   static string dir("../prod-6/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");
   TH1D* mkk = new TH1D(hname.c_str(),
                        ";M^{inv}_{K^{+}K^{-}}, GeV",
                        100, 0.98, dU  ); // dU = 1.10

   mkk->Sumw2(true);
   TCut c_here = c_Mrec+c_chi2;
   if ( type == 0 ) {   // central part
      c_here += c_cpgg;
   } else {             // side-band
      c_here += c_sbgg;
   }
   string dr = string("Mkk>>") + hname;
   a4c->Draw(dr.c_str(),c_here,"goff");

   mkk->GetXaxis()->SetTitleSize(0.04);
   mkk->GetXaxis()->SetTitleOffset(0.9);
   mkk->GetXaxis()->SetLabelFont(62);
   mkk->GetXaxis()->SetLabelSize(0.04);
   mkk->GetYaxis()->SetTitleSize(0.04);
   mkk->GetYaxis()->SetTitleOffset(0.9);
   mkk->GetYaxis()->SetLabelFont(62);
   mkk->GetYaxis()->SetLabelSize(0.04);

   mkk->SetLineWidth(2);
   mkk->SetLineColor(kBlack);
   mkk->SetMarkerStyle(20);
//    mkk->SetMarkerSize(0.5);

   return mkk;
}

//-------------------------------------------------------------------------
TH1D* Subtruct(string fname) {
//-------------------------------------------------------------------------
// Side-band subtraction
   TH1D* hst1 = plot_Mkk(fname,string("mkk_data_cp"),0);
   TH1D* hst2 = plot_Mkk(fname,string("mkk_data_sb"),1);
   TH1D* sub = (TH1D*)hst1->Clone("sub");
   sub->Add(hst1,hst2,1.,-1.);

   // errors for zero and negative bins
   for ( int ib = 1; ib <= sub->GetNbinsX(); ++ib ) {
      if ( sub->GetBinCenter(ib) < dL ) continue;
      if ( sub->GetBinContent(ib) > 1.e-3 ) {
         continue;
      }

      if( sub->GetBinError(ib) < 1e-3 ) {
         sub->SetBinError(ib,1.);
      }

      // print "not positive" bins for (hst1-hst2)
//       cout << " bin# " << ib
//            << " mkk= " << sub->GetBinCenter(ib) << " -> " << endl
//            << "\t h1: " << hst1->GetBinContent(ib)
//            << " +/- " << hst1->GetBinError(ib) << endl
//            << "\t h2: " << hst2->GetBinContent(ib)
//            << " +/- " << hst2->GetBinError(ib) << endl
//            << "\t sub:" << sub->GetBinContent(ib)
//            << " +/- " << sub->GetBinError(ib) << endl;
   }
   return sub;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void bkg_fit(string fname, string title) {
//-------------------------------------------------------------------------
   TH1D* hst = plot_Mkk(fname,string("mkk_kketa"));
   hst->SetMaximum( 2*hst->GetMaximum() );

   // import data
   RooRealVar mkk("mkk","M_{KK} (GeV)", 0.98, 1.10);
   RooDataHist hbkg("hbkg",title.c_str(),mkk,hst);

   // reverse Argus PDF for KKeta background
   RooRealVar A("A","rev. argus parameter A", 0.50, 0., 3.);
   RooRealVar B("B","rev. argus parameter B", 0., -3., 3.);
   RooConstVar L("L","the lower limit of mkk",dL);
   RooConstVar U("U","the upper limit of fit range",dU);
   RevArgus bkg("bkg","Rev_Argus",mkk,A,B,L,U);

   RooRealVar Nb("Nb","number of events in bkg.", 100., 0., 10000.);
   RooExtendPdf ebkg("ebkg","Ext.Rev.Argus",bkg,Nb);

   // Fit
//    A.setConstant(true);
   B.setConstant(true);
   RooFitResult* res = ebkg.fitTo(
                                  hbkg,
                                  Range(dL,dU),
                                  Extended(true),
                                  Save()
                                );
   assert(res);
//    res->Print();

   // Ploting
   RooPlot* xframe = mkk.frame( Title(" ") );
   hbkg.plotOn(xframe);
   ebkg.plotOn(xframe,LineColor(kBlue),LineWidth(2));

   double chi2 = xframe->chiSquare();
   ebkg.paramOn(
       xframe,
       Label(Form("#chi^{2}/ndof = %.2f",chi2)),
       ShowConstants(kTRUE),
       Layout(0.55,0.89,0.38)
              );

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();
   xframe->Draw();

   TLegend* leg = new TLegend(0.27,0.77,0.73,0.89);
   leg->AddEntry(hst,title.c_str(),"LEP");
   leg->AddEntry((TObject *)0,
         Form("#color[%i]{Reversed Argus func.}",kBlue),"");
   leg->Draw();

   c1->Update();
}

//-------------------------------------------------------------------------
void sig_fit(string fname, string title) {
//-------------------------------------------------------------------------
   TH1D* hst = plot_Mkk(fname,string("mkk_phi"));
   hst->GetYaxis()->SetMaxDigits(3);

   // import data
   RooRealVar mkk("mkk","M_{KK} (GeV)", 0.98, 1.10);
   RooDataHist hsig("hsig",title.c_str(),mkk,hst);

   // Fit signal MC by Breit-Wigner convoluted with Gauss
   RooRealVar mphi("mphi","#phi mass", Mphi, Mphi-1e-3,Mphi+1e-3);
   RooRealVar gphi("gphi","#phi width", Gphi);
//    RooRealVar gphi("gphi","#phi width", 3., 0., 10.); // this is R now
   BreitWignerPhi bwphi("bwphi","BW_phi",mkk,mphi,gphi);

   // Detector response function
   RooRealVar mg("mg","mg",0);
   RooRealVar sigma("sigma","sigma",1e-3,0.,1e-2);
   RooGaussian gauss("gauss","gauss",mkk,mg,sigma);

   // Define sampling frequency
   mkk.setBins(10000,"cache");

   // Construct convolution
   RooFFTConvPdf signal("signal","BWphi (X) gauss",mkk,bwphi,gauss);

   RooRealVar Nr("Nr","number of events in res.", 1e5,0.,1e7);
   RooExtendPdf esignal("esignal","Ext.BWphi",signal,Nr);

   // Fit
   mphi.setConstant(true);
   gphi.setConstant(true);
   RooFitResult* res = esignal.fitTo(
                                     hsig,
                                     Range(dL,1.08),
                                     Extended(true),
                                     Save()
                                   );
   assert(res);
//    res->Print();

   // Ploting
   RooPlot* xframe = mkk.frame( Title(" ") );
   hsig.plotOn(xframe);
   esignal.plotOn(xframe,LineColor(kRed),LineWidth(2));

   double chi2 = xframe->chiSquare();
   string Tchi2(Form("#chi^{2}/ndof = %.2f",chi2));
   string TNr( Form("Nr = %.0f #pm  %.0f",Nr.getValV(),Nr.getError()) );
   string Tmphi("M_{#phi} = ");
   if ( mphi.hasError() ) {
      Tmphi+=string( Form("%.3f #pm  %.3f MeV",
                          mphi.getValV()*1e3,mphi.getError()*1e3) );
   } else {
      Tmphi+=string( Form("%.3f MeV (fixed)",mphi.getValV()*1e3) );
   }
   string Tsigma( Form("#sigma = %.2f #pm  %.2f MeV",
                       sigma.getValV()*1e3,sigma.getError()*1e3) );

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);
   c1->cd();
   gPad->SetGrid();
   xframe->Draw();

   TPaveText* pt = new TPaveText(0.50,0.55,0.89,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->AddText( title.c_str() );
   TText* t1 = pt->AddText( "BW #otimes Gauss" );
   t1->SetTextColor(kRed);
   pt->AddLine(0.,0.66,1.,0.66 );
   pt->AddText( Tchi2.c_str() );
   pt->AddText( TNr.c_str() );
   pt->AddText( Tmphi.c_str() );
   pt->AddText( Tsigma.c_str() );
   pt->Draw();

   c1->Update();
}

//-------------------------------------------------------------------------
void data_fit(string fname, string title, string pdf="") {
//-------------------------------------------------------------------------
   TH1D* hst = Subtruct(fname); // subtract side-band
   hst->SetMinimum(-25.);

   // import data
   RooRealVar mkk("mkk","M_{KK} (GeV)", 0.98, 1.10);
   RooDataHist hdata("hdata",title.c_str(),mkk,hst);

   // Breit-Wigner convoluted with Gauss
   RooRealVar mphi("mphi","#phi mass", Mphi, Mphi-1e-3,Mphi+1e-3);
   RooRealVar gphi("gphi","#phi width", Gphi);
   BreitWignerPhi bwphi("bwphi","BW_phi",mkk,mphi,gphi);

   // Detector response function
   RooRealVar mg("mg","mg",0);
   RooRealVar sigma("sigma","sigma",1e-3,0.,1e-2);
   RooGaussian gauss("gauss","gauss",mkk,mg,sigma);

   // Define sampling frequency
   mkk.setBins(10000,"cache");

   // Construct convolution
   RooFFTConvPdf signal("signal","BWphi (X) gauss",mkk,bwphi,gauss);

   RooRealVar Nr("Nr","number of events in res.", 1e5,0.,1e7);
   RooExtendPdf esignal("esignal","Ext.BWphi",signal,Nr);

   // reverse Argus PDF for background
   RooRealVar A("A","rev. argus parameter A", 0.50, 0., 3.);
   RooRealVar B("B","rev. argus parameter B", 0., -3., 30.);
   RooConstVar L("L","the lower limit of mkk",dL);
   RooConstVar U("U","the upper limit of fit range",dU);
   RevArgus bkg("bkg","Rev_Argus",mkk,A,B,L,U);

   RooRealVar Nb("Nb","number of events in bkg.", 100., 0., 10000.);
   RooExtendPdf ebkg("ebkg","Ext.Rev.Argus",bkg,Nb);

   // Signal + bkg:
   RooAddPdf model("model","model",RooArgList(esignal,ebkg));
//    model.Print("t");

   // Fit
//    mphi.setConstant(true);
   gphi.setConstant(true);
//    A.setConstant(true);
//    B.setConstant(true);

//    RooFitResult* res = model.chi2FitTo( // ??? does not work ???
   RooFitResult* res = model.fitTo(
                                     hdata,
//                                      DataError(RooAbsData::SumW2),
                                     Range(dL,dU),
                                     Extended(true),
                                     Save()
                                   );
   assert(res);
//    res->Print();

   // Ploting
   RooPlot* xframe = mkk.frame( Title(" ") );
   hdata.plotOn( xframe, DataError(RooAbsData::SumW2) );
   model.plotOn(xframe,Components("bkg"),
                LineColor(kBlue),LineStyle(kDashed));
   model.plotOn(xframe,Components("signal"),
                LineColor(kGreen),LineStyle(kDashed));
   model.plotOn(xframe,LineColor(kRed),LineWidth(2));

   double chi2 = xframe->chiSquare();
   string Tchi2(Form("#chi^{2}/ndof = %.2f",chi2));
   string TNr( Form("Nr = %.0f #pm  %.0f",Nr.getValV(),Nr.getError()) );
   string TNb( Form("Nb = %.0f #pm  %.0f",Nb.getValV(),Nb.getError()) );
   string Tmphi("M(#phi) = ");
   if ( mphi.hasError() ) {
      Tmphi+=string( Form("%.3f #pm  %.3f MeV",
                          mphi.getValV()*1e3,mphi.getError()*1e3) );
   } else {
      Tmphi+=string( Form("%.3f MeV (fixed)",mphi.getValV()*1e3) );
   }
   string Tsigma( Form("#sigma = %.2f #pm  %.2f MeV",
                       sigma.getValV()*1e3,sigma.getError()*1e3) );
   string TA( Form("A = %.1f #pm  %.1f", A.getValV(),A.getError()) );
   string TB( "B = " );
   if ( B.hasError() ) {
      TB += string( Form("%.1f #pm  %.1f", B.getValV(),B.getError()) );
   } else {
      TB += string( Form("%.1f (fixed)", B.getValV()) );
   }

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad->SetGrid();
   xframe->Draw();

   TPaveText* pt = new TPaveText(0.50,0.45,0.89,0.89,"NDC");
   pt->SetTextAlign(12);
   pt->AddText( title.c_str() );
   TText* t1 = pt->AddText( "BW #otimes Gauss + Argus" );
   t1->SetTextColor(kRed);
   pt->AddLine(0.,0.775,1.,0.775 );
   pt->AddText( Tchi2.c_str() );
   pt->AddText( TNr.c_str() );
   pt->AddText( TNb.c_str() );
   pt->AddText( Tmphi.c_str() );
   pt->AddText( Tsigma.c_str() );
   pt->AddText( TA.c_str() );
   pt->AddText( TB.c_str() );
   pt->Draw();

   c1->Update();
}

//-------------------------------------------------------------------------
void data_intfr_fit(string fname, string title) {
//-------------------------------------------------------------------------
   TH1D* hst = Subtruct(fname); // subtract side-band
   hst->SetMinimum(-25.);

   // import data
   RooRealVar mkk("mkk","M_{KK} (GeV)", 0.98, 1.10);
   mkk.setBins(10000,"cache"); // Define sampling frequency
   RooDataHist hdata("hdata",title.c_str(),mkk,hst);
   hdata.printMultiline(cout,1,true);
/*
   // Detector response function
   RooRealVar mg("mg","mg",0);
   RooRealVar sigma("sigma","sigma",1e-3,0.,1e-2);
   RooGaussian gauss("gauss","gauss",mkk,mg,sigma);

   // Breit-Wigner convoluted with Gauss
   RooRealVar mphi("mphi","#phi mass", Mphi, Mphi-1e-3,Mphi+1e-3);
   RooRealVar gphi("gphi","#phi width", Gphi);
   BreitWignerPhi bwphi("bwphi","BW_phi",mkk,mphi,gphi);
   RooFFTConvPdf Res("Res","BW (X) gauss",mkk,bwphi,gauss);

   // reverse Argus PDF for background
   RooRealVar A("A","argus parameter A", 0.50, 0., 3.);
   RooRealVar B("B","argus parameter B", 0., -3., 30.);
   RooConstVar L("L","argus lower limit of mkk",dL);
   RooConstVar U("U","argus upper limit of fit range",dU);
   RevArgus argus("argus","RevArgus",mkk,A,B,L,U);
   RooFFTConvPdf Bkg("Bkg","Argus (X) gauss",mkk,argus,gauss);

   // Interference: (argus + resonance + interference)
   RooRealVar Ang("Ang","mixing angle", 0., -M_PI, M_PI);
   Interference intfr0( "intfr0","Ar + BW + Intfr",
                       mkk, mphi, gphi, A, B, L, U, Ang);
   RooFFTConvPdf Intfr("Intfr","Interference (X) gauss",mkk,intfr0,gauss);

   // Sum
   // 1) two variables responsible for number of events
//    RooRealVar Nr("Nr","number of events in res.", 3e3,0.,1e5);
   RooRealVar Fb("Fb","ratio bkg. to resonance", 1.e-1,0.,0.5);

   // 2) reparameterizing
   RooFormulaVar C1("C1","Fb-sqrt(Fb)",RooArgSet(Fb)); // bkg.
   RooFormulaVar C2("C2","1-sqrt(Fb)",RooArgSet(Fb)); // res.
   RooFormulaVar C3("C3","sqrt(Fb)",RooArgSet(Fb)); // intfr.
   RooAddPdf Model("Model","sum",
         RooArgList(Res,Bkg,Intfr),RooArgList(C1,C2,C3));

   // 3) extended pdf ????
//    RooExtendPdf Model("Model","ExtModel",model0,Nr);

   // Print GraphViz DOT file with representation of tree
   // > dot -Tpdf -o model.pdf model.dot
//    Model.graphVizTree("model.dot") ;

   // Fit
   mphi.setConstant(true);
   gphi.setConstant(true);
//    A.setConstant(true);
   B.setConstant(true);

   RooFitResult* res = Model.fitTo(
                                     hdata,
                                     Range(dL,dU),
                                     Extended(true),
                                     Save()
                                   );
   assert(res);
//    res->Print();
*/
   // Ploting
   RooPlot* xframe = mkk.frame( Title(" ") );
   hdata.plotOn( xframe, DataError(RooAbsData::SumW2) );
//    Model.plotOn(xframe,LineColor(kRed),LineWidth(2));
//    Model.plotOn(xframe,Components(Bkg),
//                 LineColor(kBlue),LineStyle(kDashed));
//    Model.plotOn(xframe,Components(Res),
//                 LineColor(kGreen),LineStyle(kDashed));

//    Model.paramOn(
//        xframe,
//        ShowConstants(kTRUE),
//        Layout(0.55,0.89,0.89)
//                 );

   TCanvas* c1 = new TCanvas("c1","...",0,0,900,900);

   c1->cd();
   gPad->SetGrid();
   xframe->Draw();

   c1->Update();
}

//----------------------------------------------------------------------
int main(int argc, char* argv[]) {
//----------------------------------------------------------------------
   vector<string> fnames = {
      "data_09psip_all.root",
      "data_12psip_all.root",
      "data_3650_all.root",
      "mcinc_09psip_all.root",
      "mcinc_12psip_all.root",
      "mcsig_kkmc_09.root",
      "mcsig_kkmc_12.root",
      "mckketa_kkmc_09.root",
      "mckketa_kkmc_12.root"
   };
   vector<string> titles = {
      "Data 2009",
      "Data 2012",
      "E=3.65 GeV",
      "MC incl. 2009",
      "MC incl. 2012",
      "MC #phi#eta 2009",
      "MC #phi#eta 2012",
      "MC KK#eta 2009",
      "MC KK#eta 2012"
   };
   // Numbers of items for date:
//    int id = 0, ii = 3, is = 5, ib = 7; // 2009
   int id = 1, ii = 4, is = 6, ib = 8; // 2012
   printf("id = %i ii = %i is = %i ib = %i\n",id,ii,is,ib);

   // root
   TApplication app("ROOT Application", &argc, argv);

//    bkg_fit(fnames.at(ib),titles.at(ib)); // fit kketa by Argus
//    sig_fit(fnames.at(is),titles.at(is)); // fit phi-eta by BW x Gauss
//    data_fit(fnames.at(id),titles.at(id));

   data_intfr_fit(fnames.at(id),titles.at(id));

   app.Run();
}

