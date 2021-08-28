// Pictures for presentation

//--------------------------------------------------------------------
// Common parameters:
struct Params {
   // names of files to use
   vector<string> fnames;

   int date;
   int slct;    // >0: g-eta; <0: phi-eta
                // 1 - photon eff; 2 - eta eff;
   int use_rew; // 0 - no weights;
                // 1 - calculate weights

   TCut Cmcsig; // for mc-signal
   TCut Cbg;    // cuts against the background
   TCut Cph;    // cuts for selection photons
   TCut Ceta;   // cuts for selection eta

   //-----------------------------------------------------------------
   Params(int dat = 2012, int slc = 2, int rew = 0) {
   //-----------------------------------------------------------------
      date = dat;
      slct = slc;
      use_rew = rew;

      // name of folder with root-files
      string dir("prod-12eff/");
      fnames = {
         "data_09psip_all.root", "mcinc_09psip_all.root",
         "data_12psip_all.root", "mcinc_12psip_all.root"
      };
      for ( auto& fn : fnames ) {
         fn = dir + fn;
      }

      // mc-signal
      if ( slct > 0 ) {    // gamma-eta
         Cmcsig += TCut("decj==22");
      } else {             // phi-eta
         Cmcsig += TCut("decj==68");
      }

      // cuts against the background
      if ( slct > 0 ) {    // gamma-eta
         Cbg += TCut("m2fr<0.002");
         Cbg += TCut("fabs(Cg0)<0.8");
         Cbg += TCut("fabs(Cg1)<0.8 && Eg1>0.2"); // ?
      } else {             // phi-eta
         Cbg += TCut("m2fr<0.002");
      }

      // cuts for selection photons
//       Cph += TCut("fabs(Cg2)<0.92");
      Cph += TCut("fabs(Cg2)<0.8||(fabs(Cg2)>0.85&&fabs(Cg2)<0.92)");
      if ( slct > 0 ) {  // gamma-eta
         Cph += TCut("Eg2>0.1&&Eg2<1.4");
      } else {             // phi-eta
         Cph += TCut("Eg2>0.05&&Eg1<1.45");
      }

      // cuts for selection eta
      Ceta += TCut("fabs(Ceta)<0.9");
//       Ceta += TCut("fabs(Cg2)<0.92");
      Ceta += TCut("fabs(Cg2)<0.8||(fabs(Cg2)>0.85&&fabs(Cg2)<0.92)");
      if ( slct > 0 ) {  // gamma-eta
         Ceta += TCut("Peta>1.3&&Peta<1.7");
      } else {             // phi-eta
         Ceta += TCut("Peta>1.15&&Peta<1.5");
      }
   }

   //-----------------------------------------------------------------
   TFile* OpenFile(int mc = 0) { // 1 for MC
   //-----------------------------------------------------------------
      // open file
      int idx = mc + ((date==2009) ? 0 : 2);
      string dfname = fnames[idx];
      TFile* froot = TFile::Open(dfname.c_str(),"READ");
      if( froot == 0 ) {
         cerr << "can not open " << dfname << endl;
         exit(0);
      }
      if ( slct > 0 ) {
         froot->cd("PsipJpsiGammaEta");
      } else {
         froot->cd("PsipJpsiPhiEta");
      }
      return froot;
   }
};
//--------------------------------------------------------------------

//--------------------------------------------------------------------
void SetHstFace(TH1* hst) {
//--------------------------------------------------------------------
   TAxis* X = hst->GetXaxis();
   if ( X ) {
      X->SetLabelFont(62);
      X->SetLabelSize(0.04);
      X->SetTitleFont(62);
      X->SetTitleSize(0.04);
   }
   TAxis* Y = hst->GetYaxis();
   if ( Y ) {
      Y->SetLabelFont(62);
      Y->SetLabelSize(0.04);
      Y->SetTitleFont(62);
      Y->SetTitleSize(0.04);
   }
   TAxis* Z = hst->GetZaxis();
   if ( Z ) {
      Z->SetLabelFont(62);
      Z->SetLabelSize(0.04);
      Z->SetTitleFont(62);
      Z->SetTitleSize(0.04);
   }
}

//--------------------------------------------------------------------
void pi0_rej(Params* p, string pdf) {
//--------------------------------------------------------------------
   TFile* froot = p->OpenFile(0); // data

   TH1D* Mg2_pi0 = (TH1D*)gDirectory->Get("Mg2_pi0");
   if ( !Mg2_pi0 ) {
      cout << __func__ << " can not find Mg2_pi0" << endl;
      exit(0);
   }
   Mg2_pi0 -> SetLineWidth(2);
   Mg2_pi0 -> SetTitle("M^{ 2}_{inv}(#gamma#gamma);(GeV/c^{2})^{2}");
   SetHstFace(Mg2_pi0);

   vector<TLine*> lines(2,nullptr);
   lines[0] = new TLine(0.013,0,0.013,Mg2_pi0 -> GetMaximum());
   lines[1] = new TLine(0.022,0,0.022,Mg2_pi0 -> GetMaximum());
   for (auto& l : lines ) {
      l->SetLineColor(kGreen+2);
      l->SetLineWidth(3);
      l->SetLineStyle(kDashed); // or kSolid
   }

   TLegend* leg = new TLegend(0.69,0.79,0.89,0.89);
//    leg -> SetHeader("data 2012", "C");
   leg -> AddEntry( Mg2_pi0, ("Data " + to_string(p->date)).c_str(), "L" );

   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,800);
   c1->cd();
   c1->Print((pdf+"[").c_str()); // open pdf-file

   Mg2_pi0->Draw();
   for (auto& l : lines ) l->Draw();
   leg->Draw();

   c1->Print(pdf.c_str()); // add to pdf-file
   c1->Print((pdf+"]").c_str()); // close pdf-file
}

//--------------------------------------------------------------------
void Mrecpipig(Params* p, string pdf) {
//--------------------------------------------------------------------
   TFile* froot = p->OpenFile(1); // MC

   TH1D* mcMrg2T = (TH1D*)gDirectory->Get("mcMrg2T");
   TH1D* mcMrg2F = (TH1D*)gDirectory->Get("mcMrg2F");
   if ( !mcMrg2T || !mcMrg2F ) {
      cout << __func__ << " can not find mcMrg2..." << endl;
      exit(0);
   }
   mcMrg2T -> SetLineWidth(2);
   mcMrg2F -> SetLineWidth(2);
   mcMrg2F -> SetLineColor(kRed+1);
   mcMrg2T -> SetTitle("M^{ 2}_{ recoil}(#pi^{+} #pi^{-} #gamma)"
                       ";(GeV/c^{2})^{2}");
   SetHstFace(mcMrg2T);

   vector<TLine*> lines(2,nullptr);
   lines[0] = new TLine(0.1,0,0.1, mcMrg2T -> GetMaximum());
   lines[1] = new TLine(0.5,0,0.5, mcMrg2T -> GetMaximum());
   for (auto& l : lines ) {
      l->SetLineColor(kGreen+2);
      l->SetLineWidth(3);
      l->SetLineStyle(kDashed); // or kSolid
   }

   TLegend* leg = new TLegend(0.69,0.69,0.89,0.89);
   leg -> SetHeader(("MC " + to_string(p->date)).c_str(), "C");
   leg -> AddEntry( mcMrg2T, "signal" , "L" );
   leg -> AddEntry( mcMrg2F, "background" , "L" );

   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,800);
   c1->cd();
   c1->Print((pdf+"[").c_str()); // open pdf-file

   mcMrg2T -> Draw();
   mcMrg2F -> Draw("SAME");
   for (auto& l : lines ) l->Draw();
   leg->Draw();

   c1->Print(pdf.c_str()); // add to pdf-file
   c1->Print((pdf+"]").c_str()); // close pdf-file
}

//--------------------------------------------------------------------
void M2mis(Params* p, string pdf) {
//--------------------------------------------------------------------
   TFile* froot = p->OpenFile(1); // MC

   TTree* eff_eta = (TTree*)gDirectory->Get("eff_eta");
   if ( !eff_eta ) {
      cout << __func__ << ": can not find eff_eta" << endl;
      exit(0);
   }

   TH1D* h1 = new TH1D("h1", "", 75,0.,0.0075);
   TH1D* h2 = new TH1D("h2", "", 75,0.,0.0075);

   TCut Cbg1; // cuts against the background without m2fr cut
   if ( p->slct > 0 ) {    // gamma-eta
      Cbg1 += TCut("fabs(Cg0)<0.8");
//       Cbg1 += TCut("fabs(Cg1)<0.8 && Eg1>0.2"); // ?
   }

   TCut cut1 = Cbg1 + p->Ceta + p->Cmcsig;
   TCut cut2 = Cbg1 + p->Ceta + !(p->Cmcsig);

   eff_eta->Draw( "m2fr>>h1",cut1,"goff");
   eff_eta->Draw( "m2fr>>h2",cut2,"goff");

   h1 -> SetLineWidth(2);
   h2 -> SetLineWidth(2);
   h2 -> SetLineColor(kRed+1);
   h1 -> SetTitle("M^{ 2}_{ missing};(GeV/c^{2})^{2}");
   SetHstFace(h1);

   vector<TLine*> lines(1,nullptr);
   lines[0] = new TLine(0.002,0,0.002, h1 -> GetMaximum());
   for (auto& l : lines ) {
      l->SetLineColor(kGreen+2);
      l->SetLineWidth(3);
      l->SetLineStyle(kDashed); // or kSolid
   }

   TLegend* leg = new TLegend(0.69,0.69,0.89,0.89);
   leg -> SetHeader(("MC " + to_string(p->date)).c_str(), "C");
   leg -> AddEntry( h1, "signal" , "L" );
   leg -> AddEntry( h2, "background" , "L" );

   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,800);
   c1->cd();
   gPad->SetLogy();
   c1->Print((pdf+"[").c_str()); // open pdf-file


   h1 -> Draw();
   h2 -> Draw("SAME");
   for (auto& l : lines ) l->Draw();
   leg->Draw();

   c1->Print(pdf.c_str()); // add to pdf-file
   c1->Print((pdf+"]").c_str()); // close pdf-file
}

//--------------------------------------------------------------------
void BgStd(Params* p, string pdf) {
//--------------------------------------------------------------------
   TFile* froot = p->OpenFile(1); // MC

   TTree* eff_eta = (TTree*)gDirectory->Get("eff_eta");
   if ( !eff_eta ) {
      cout << __func__ << ": can not find eff_eta" << endl;
      exit(0);
   }

   TH1D* he1 = new TH1D("he1", "", 68,0.,1.7);
   TH1D* he2 = new TH1D("he2", "", 68,0.,1.7);

   TCut Cbg1; // cuts against the background without Cbg1 cut
   if ( p->slct > 0 ) {    // gamma-eta
      Cbg1 += TCut("m2fr<0.002");
      Cbg1 += TCut("fabs(Cg0)<0.8");
//       Cbg1 += TCut("fabs(Cg1)<0.8 && Eg1>0.2"); // ?
   }

   TCut cut1 = Cbg1 + p->Ceta + p->Cmcsig;
   TCut cut2 = Cbg1 + p->Ceta + !(p->Cmcsig);

   eff_eta->Draw( "Eg1>>he1",cut1,"goff");
   eff_eta->Draw( "Eg1>>he2",cut2,"goff");

   he1 -> SetLineWidth(2);
   he2 -> SetLineWidth(2);
   he2 -> SetLineColor(kRed+1);
   he1 -> SetTitle("Energy of reconstructed photon;GeV");
   SetHstFace(he1);

   vector<TLine*> lines(1,nullptr);
   lines[0] = new TLine(0.2,0,0.2, he1 -> GetMaximum());
   for (auto& l : lines ) {
      l->SetLineColor(kGreen+2);
      l->SetLineWidth(3);
      l->SetLineStyle(kDashed); // or kSolid
   }

   TLegend* leg = new TLegend(0.69,0.69,0.89,0.89);
   leg -> SetHeader(("MC " + to_string(p->date)).c_str(), "C");
   leg -> AddEntry( he1, "signal" , "L" );
   leg -> AddEntry( he2, "background" , "L" );

   TCanvas* c1 = new TCanvas("c1","...",0,0,1200,800);
   c1->cd();
   c1->Print((pdf+"[").c_str()); // open pdf-file


   he1 -> Draw();
   he2 -> Draw("SAME");
   for (auto& l : lines ) l->Draw();
   leg->Draw();

   c1->Print(pdf.c_str()); // add to pdf-file
   c1->Print((pdf+"]").c_str()); // close pdf-file
}

//--------------------------------------------------------------------
void MatchPic(Params* p, int type, string pdf) {
//--------------------------------------------------------------------
   TFile* froot = p->OpenFile(1); // MC

   vector<TH1D*> mc(2,nullptr);
   if ( type == 0 ) {
      if ( p->slct > 0 ) {
         mc[0] = (TH1D*)gDirectory->Get("mcrET");
         mc[1] = (TH1D*)gDirectory->Get("mcrEF");
      } else {
         mc[0] = (TH1D*)gDirectory->Get("Emc_rET");
         mc[1] = (TH1D*)gDirectory->Get("Emc_rEF");
      }
   } else {
      if ( p->slct > 0 ) {
         mc[0] = (TH1D*)gDirectory->Get("mcdThT");
         mc[1] = (TH1D*)gDirectory->Get("mcdThF");
      } else {
         mc[0] = (TH1D*)gDirectory->Get("Emc_dThT");
         mc[1] = (TH1D*)gDirectory->Get("Emc_dThF");
      }
   }
   if ( !mc[0] || !mc[1] ) {
      cout << __func__ << " can not find histos" << endl;
      exit(0);
   }
   mc[0] -> SetLineWidth(2);
   mc[1] -> SetLineWidth(2);
   mc[1] -> SetLineColor(kRed+1);
   SetHstFace(mc[0]);
   vector<TLine*> lines;
   if ( type == 0 ) {
      mc[0] -> SetTitle("E_{#gamma}(pred)/E_{#gamma}(rec)");
      lines.push_back( new TLine(0.4,0,0.4, mc[0] -> GetMaximum()) );
      lines.push_back( new TLine(1.8,0,1.8, mc[0] -> GetMaximum()) );
   } else {
      mc[0] -> SetTitle("#delta#Theta (pred-rec); degree");
      lines.push_back( new TLine(10,0,10, mc[0] -> GetMaximum()) );
   }

   for (auto& l : lines ) {
      l->SetLineColor(kGreen+2);
      l->SetLineWidth(3);
      l->SetLineStyle(kDashed); // or kSolid
   }

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   gPad->SetLogy();
   c1->Print((pdf+"[").c_str()); // open pdf-file

   mc[0] -> Draw();
   mc[1] -> Draw("SAME");
   for (auto& l : lines ) l->Draw();

   c1->Print(pdf.c_str()); // add to pdf-file
   c1->Print((pdf+"]").c_str()); // close pdf-file
}

//--------------------------------------------------------------------
void MggFit(Params* p, int mc, string pdf) {
//--------------------------------------------------------------------
   TFile* froot = p->OpenFile(mc); // DATA or MC

   TTree* eff_eta = (TTree*)gDirectory->Get("eff_eta");
   if ( !eff_eta ) {
      cout << __func__ << ": can not find eff_eta" << endl;
      exit(0);
   }

   double meta = 0.547862;
   double weta = 3.5*0.008;
   TH1D* hf = new TH1D("hf", "", 75,meta-weta,meta+weta);
   TCut cut = p->Cbg + p->Ceta + TCut("fl>1.5");
   eff_eta->Draw( "mggf>>hf",cut,"goff");

   hf -> SetLineWidth(2);
   string title("M_{inv}(#gamma#gamma)");
   title += ( (mc == 0) ? " DATA " : " MC ");
   title += to_string(p->date);
   title += ";GeV/c^{2}";
   hf -> SetTitle(title.c_str());
   SetHstFace(hf);

   TF1* gs = (TF1*)gROOT->GetFunction("gaus");
   gs->SetLineWidth(2);
   gs->SetLineColor(kRed);
   gStyle->SetFitFormat(".3g");

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   c1->Print((pdf+"[").c_str()); // open pdf-file

   hf -> Fit(gs,"","",meta-3*0.008,meta+3*0.008);
//    hf -> Draw();

   c1->Print(pdf.c_str()); // add to pdf-file
   c1->Print((pdf+"]").c_str()); // close pdf-file
}

//--------------------------------------------------------------------
void DataMC(Params* p, int type, string pdf) {
//--------------------------------------------------------------------
   string Svar;
   vector<TH1D*> datmc(2,nullptr);

   TCut Ceta;   // cuts for selection eta without Peta & Ceta cuts
//    Ceta += TCut("fabs(Ceta)<0.9");
   Ceta += TCut("fabs(Cg2)<0.8||(fabs(Cg2)>0.85&&fabs(Cg2)<0.92)");
   TCut cut = p->Cbg + Ceta;
   if ( type >= 20 ) {
      cut += TCut("fl<1.5"); // partially
   } else if ( type >= 10 ) {
      cut += TCut("fl>1.5"); // fully
   }

   TLegend* leg = nullptr;

   TFile* froot = p->OpenFile(0); // DATA
   TTree* eff_eta = (TTree*)gDirectory->Get("eff_eta");
   if ( type%2 == 0 ) {
      Svar="Peta";
      leg = new TLegend(0.57,0.72,0.89,0.89);
//       leg -> SetHeader("P_{#eta}", "C");
      if ( p->slct > 0 ) {
         datmc[0] = new TH1D("hdat", ";P_{#eta}, GeV/c", 55,1.25,1.8);
      } else {
         datmc[0] = new TH1D("hdat", ";P_{#eta}, GeV/c", 45,1.1,1.55);
      }
   } else {
      Svar="Ceta";
      leg = new TLegend(0.32,0.20,0.68,0.40);
//       leg -> SetHeader("cos #Theta(#eta)", "C");
      if ( p->slct > 0 ) {
         datmc[0] = new TH1D("hdat", ";cos #Theta(#eta)", 40,-1.,1.);
      } else {
         datmc[0] = new TH1D("hdat", ";cos #Theta(#eta)", 30,-1.,1.);
      }
   }

   eff_eta->Draw( (Svar+">>hdat").c_str(),cut,"goff");
   cout << "data " << datmc[0]->GetEntries() << endl;

   TFile* frootMC = p->OpenFile(1); // MC
   TTree* eff_etaMC = (TTree*)gDirectory->Get("eff_eta");
   datmc[1] = (TH1D*)datmc[0]->Clone("hmc");
   eff_etaMC->Draw( (Svar+">>hmc").c_str(),cut,"goff");
   cout << "  mc " << datmc[1]->GetEntries() << endl;

   TH1D* tmp = (TH1D*)datmc[0]->Clone("tmp");
   eff_etaMC->Draw( (Svar+">>tmp").c_str(),cut+!(p->Cmcsig),"goff");
   cout << "mc-F: " << tmp->GetEntries() << endl;

   datmc[0] -> SetLineWidth(2);
   datmc[1] -> SetLineWidth(2);
   datmc[1] -> SetLineColor(kRed+1);
   double norm = datmc[0]->GetEntries()/datmc[1]->GetEntries();
   datmc[1] -> Scale(norm);
   SetHstFace(datmc[0]);
   leg -> AddEntry( datmc[0], ("Data " + to_string(p->date)).c_str(), "L" );
   leg -> AddEntry( datmc[1], ("MC " + to_string(p->date)).c_str(), "L" );

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);
   c1->cd();
   c1->Print((pdf+"[").c_str()); // open pdf-file

   datmc[0] -> Draw("E1");
   datmc[1] -> Draw("SAME HIST");
   if ( tmp ) {
      tmp -> Scale(norm);
      tmp -> SetLineWidth(2);
      tmp -> SetLineColor(kMagenta);
      tmp -> SetLineStyle(kDashed);
      leg -> AddEntry( tmp, "background", "L" );
      tmp -> Draw("SAME HIST");
   }
   leg -> Draw();
   gPad -> RedrawAxis();

   c1->Print(pdf.c_str()); // add to pdf-file
   c1->Print((pdf+"]").c_str()); // close pdf-file
}

//--------------------------------------------------------------------
void eta_eff_pr() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);
   gStyle->SetLegendFont(62);
   gStyle->SetLegendTextSize(0.04);
//    gStyle->SetStatFontSize(0.07);
//    gStyle->SetStatX(0.89);
//    gStyle->SetStatY(0.89);

//    Params* par = new Params(2012,2,0); // 2012 gamma eta
//    pi0_rej(par,"dat12_M2pi0.pdf");
//    Mrecpipig(par,"mc12_Mrecpipig.pdf");
//    M2mis(par,"mc12_M2mis.pdf");
//    BgStd(par,"mc12_Eg1.pdf");
//    MatchPic(par,0,"mc12_rE.pdf");
//    MatchPic(par,1,"mc12_dTh.pdf");
//    MggFit(par,0,"dat12_Mggf.pdf");
//    MggFit(par,1,"mc12_Mggf.pdf");
//
//    DataMC(par,0,"Peta_mcdat12.pdf");
//    DataMC(par,1,"Ceta_mcdat12.pdf");
//    DataMC(par,10,"Peta_mcdat12_f.pdf");
//    DataMC(par,11,"Ceta_mcdat12_f.pdf");
//    DataMC(par,20,"Peta_mcdat12_p.pdf");
//    DataMC(par,21,"Ceta_mcdat12_p.pdf");

   Params* par = new Params(2012,-2,0); // 2012 phi eta
   M2mis(par,"mc12phi_M2mis.pdf");
//    MatchPic(par,0,"mc12phi_rE.pdf");
//    MatchPic(par,1,"mc12phi_dTh.pdf");
//
//    DataMC(par,0,"Peta_phi_mcdat12.pdf");
//    DataMC(par,1,"Ceta_phi_mcdat12.pdf");
//    DataMC(par,10,"Peta_phi_mcdat12_f.pdf");
//    DataMC(par,11,"Ceta_phi_mcdat12_f.pdf");
//    DataMC(par,20,"Peta_phi_mcdat12_p.pdf");
//    DataMC(par,21,"Ceta_phi_mcdat12_p.pdf");
}
