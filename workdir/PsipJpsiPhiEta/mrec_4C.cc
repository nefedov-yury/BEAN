// plot Mrec after 4C-kinematik fit
// we require cuts for: chi^2(4C) && |Mgg-Meta| && |Mkk-Mphi|
// Data, inclusive MC, signal MC and MC KKeta
// -> mrec_4Cyear.pdf

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

//-------------------------------------------------------------------------
TH1D* get_Mrec(string fname, string hname, vector<double>& num) {
//-------------------------------------------------------------------------
   const double meta = 0.547862; //  547.862 +/- 0.017 MeV
   const double mphi = 1.019461; // 1019.461 +/- 0.019 MeV
   const double mk   = 0.493677; // 493.677  +/- 0.016 MeV

  // name of folder with root files
   static string dir("prod-94c/");
   fname = dir + fname;
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(0);
   }

   froot->cd("PsipJpsiPhiEta");
   TTree* a4c = (TTree*)gDirectory->Get("a4c");

   TH1D* hst = new TH1D(hname.c_str(),
                        ";M_{rec}(#pi^{+}#pi^{-}), GeV/c^{2}"
                        ";Entries / 0.001GeV/c^{2}",
                        200,3.0,3.2 // 1bin = 1MeV
//                        400,3.0,3.2 // 1bin = 0.5MeV
                       );
   hst->Sumw2(true);

   TCut c_here( Form("ch2<80"
            "&&abs(Mgg-%.6f)<0.024"
            "&&Mkk>%.6f&&Mkk<1.08", meta,2*mk));
//             "&&abs(Mgg-%.6f)<0.024"
//             "&&abs(Mkk-%.6f)<0.01", meta,mphi));
   string dr = string("Mrec>>") + hname;
   a4c->Draw(dr.c_str(),c_here,"goff");

   // number of events in central part and in side-bands:
   TCut c_Mr0("abs(Mrec-3.097)<0.005"); // [3.092, 3.102]
   TCut c_MrL("Mrec>3.00&&Mrec<3.05");  // [3.00, 3.05] left SB
   TCut c_MrR("Mrec>3.144&&Mrec<3.194");// [3.144, 3.194] right SB
   double n0 = a4c->Draw("Mrec",c_here+c_Mr0,"goff");
   double nL = a4c->Draw("Mrec",c_here+c_MrL,"goff");
   double nR = a4c->Draw("Mrec",c_here+c_MrR,"goff");
   num.push_back(n0);
   num.push_back(nL);
   num.push_back(nR);

   return hst;
}

//-------------------------------------------------------------------------
void mrec_4C() {
//-------------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetStatFont(62);
//    gStyle->SetLegendTextSize(0.05);

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
      "MC signal #phi#eta",
      "MC signal #phi#eta",
      "MC non-#phi KK#eta",
      "MC non-#phi KK#eta"
   };

//    int id = 0, ii = 3, is = 5, ib = 7; // 2009
   int id = 1, ii = 4, is = 6, ib = 8; // 2012

   string pdf = (id==0) ? string("mrec_4C09") : string("mrec_4C12");
   pdf += string(".pdf");

   TH1D* hst[20];
   vector<int> tmp { id, ii, is, 2 }; // 2 - Continuum data
   vector<double> nums;
   nums.reserve(100);
   for( auto i : tmp )  {
      string hname = string("mrec4C_") + to_string(i+1);
      hst[i] = get_Mrec(fnames.at(i),hname,nums);
   }


   for ( int it = 0; it < 1; it+=10 ) {
      hst[id+it]->SetLineWidth(2);
      hst[id+it]->SetLineColor(kBlack); // data
      hst[id+it]->SetMarkerStyle(20);
      hst[id+it]->SetMarkerSize(0.7);

      hst[is+it]->SetLineWidth(2);
      hst[is+it]->SetLineColor(kGreen+2); // MC phi eta
//
      hst[ii+it]->SetLineWidth(2);
//       hst[ii+it]->SetLineColor(kMagenta+1); // inclusive MC
      hst[ii+it]->SetLineColor(kGreen+2); // inclusive MC

      hst[2+it]->SetLineColor(kBlue+1); // Continuum data
      hst[2+it]->SetLineWidth(2);
   }


   double ninc = nums[4]+nums[5];
   double einc = sqrt(ninc);
   // normalize inclusive MC on data
   const double Ito09 =  1.022; // +/- 0.008
   const double Ito12 =  0.853; // +/- 0.005
   if ( id == 0 ) { // 2009
      hst[ii]->Scale( Ito09 );
      for ( int i = 3; i < 6; i++ ) {
         nums[i] *= Ito09;
      }
      ninc *= Ito09;
      einc *= Ito09;
   } else if ( id == 1 ) { // 2012
      hst[ii]->Scale( Ito12 );
      for ( int i = 3; i < 6; i++ ) {
         nums[i] *= Ito12;
      }
      ninc *= Ito12;
      einc *= Ito12;
   }

   double nsig = nums[7]+nums[8];
   double esig = sqrt(nsig);
   // normalize signal MC on data
   double dmax = hst[id]->GetMaximum();
   double mcs_max = hst[is]->GetMaximum();
   hst[is]->Scale( dmax/mcs_max );
   for ( int i = 6; i < 9; i++) {
         nums[i] *= dmax/mcs_max;
   }
   nsig *= dmax/mcs_max;
   esig *= dmax/mcs_max;

   // normaliza continuum data to data
   const double Cto09 =  3.276; // +/- 0.021
   const double Cto12 = 10.273; // +/- 0.067
   double n_cont= hst[2]->GetEntries(); 
   double e_cont = sqrt(n_cont);
   bool is_cont = n_cont > 0;
   if ( is_cont ) {
      if ( id == 0 ) { // 2009
         hst[2]->Scale( Cto09 );
         for ( int i = 9; i < 12; i++ ) {
            nums[i] *= Cto09;
         }
         n_cont *= Cto09;
         e_cont *= Cto09;
      } else if ( id == 1 ) { // 2012
         hst[2]->Scale( Cto12 );
         for ( int i = 9; i < 12; i++ ) {
            nums[i] *= Cto12;
         }
         n_cont *= Cto12;
         e_cont *= Cto12;
      } else {
         cout << "id= " << id << " ???\n";
         return;
      }
   }

//    cout << " DATA SB: L-R " << nums[1] << " " << nums[2] << endl;  
//    cout << " MC_I SB: L-R " << nums[4] << " " << nums[5] << endl;  
//    cout << " MC_S SB: L-R " << nums[7] << " " << nums[8] << endl;  
//    cout << " CONT SB: L-R " << nums[10] << " " << nums[11] << endl;  

   double ndat = nums[1]+nums[2];
   double edat = sqrt(ndat);

   cout << " DATA CP: " << nums[0] << endl;  
   cout << " DATA SB: " << ndat/10 << " +/- " << edat/10 << endl;  
   cout << " MC_I SB: " << ninc/10 << " +/- " << einc/10 << endl;  
   cout << " MC_S SB: " << nsig/10 << " +/- " << esig/10 << endl;  
   cout << " CONT ALL: " << n_cont/20 << " +/- " << e_cont/20 << endl;  
   cout << " CONT SB: " << nums[10] + nums[11] << endl;  


   //-----------------------------------------------------------------------

   TLine* lR = new TLine;
   lR->SetLineColor(kRed+1);
   lR->SetLineWidth(2);
   lR->SetLineStyle(7);
   TLine* lB = new TLine;
   lB->SetLineColor(kMagenta);
   lB->SetLineWidth(2);
   lB->SetLineStyle(7);

   TCanvas* c1 = new TCanvas("c1","...",0,0,800,800);

   c1->cd();
   gPad->SetGrid();

   gPad->SetLogy(true);
   double hmin = 0.5;
   hst[id]->SetMinimum(hmin);
   double hmax = hst[id]->GetMaximum();

   SetHstFace(hst[id]);
   hst[id]->GetYaxis()->SetTitleOffset(1.25);

   hst[id]->Draw("E");
   hst[ii]->Draw("HIST SAME");
//    hst[is]->Draw("HIST SAME");
   if ( is_cont ) {
      hst[2]->Draw("HIST SAME");
   }
   hst[id]->Draw("E SAME");

   lR->DrawLine(3.092,hmin,3.092,hmax);
   lR->DrawLine(3.102,hmin,3.102,hmax);
   lB->DrawLine(3.001,hmin,3.001,0.1*hmax);
   lB->DrawLine(3.050,hmin,3.050,0.1*hmax);
   lB->DrawLine(3.144,hmin,3.144,0.1*hmax);
   lB->DrawLine(3.194,hmin,3.194,0.1*hmax);

//    TLegend* leg = new TLegend(0.58,0.75,0.89,0.89);
   TLegend* leg = new TLegend(0.58,0.70,0.89,0.89);
   leg->AddEntry(hst[id],titles[id].c_str(),"PLE");
   if ( is_cont ) {
      leg->AddEntry(hst[2],"Continuum data","L");
   }
   leg->AddEntry(hst[ii],titles[ii].c_str(),"L");
//    leg->AddEntry(hst[is],titles[is].c_str(),"L");
//    leg->AddEntry(lR,"M_{rec} #in [3.092,3.102])","L");
   leg->AddEntry(lR,"signal region","L");
   leg->AddEntry(lB,"side-band","L");
   leg->Draw();

   gPad->RedrawAxis();

   c1->Update();
   c1->Print(pdf.c_str());
}

