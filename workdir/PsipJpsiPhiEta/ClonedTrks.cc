// ClonedTrks.cc - plot pictures for cloned tracks study

// {{{1 helper functions
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

// {{{1 Get histograms
//--------------------------------------------------------------------
vector<TH1D*> get_hists(string fname, vector<string> hnames) {
//--------------------------------------------------------------------
   TFile* froot = TFile::Open(fname.c_str(),"READ");
   if( froot == 0 ) {
      cerr << "can not open " << fname << endl;
      exit(EXIT_FAILURE);
   }

   froot->cd("DelClonedTrks");

   vector<TH1D*> hists;
   hists.reserve(8);

   for ( const auto& hn : hnames ) {
      auto hst = (TH1D*)gROOT->FindObject( hn.c_str() );
      if ( !hst ) {
         cerr << " can not find histo:" << hn << endl;
         exit(EXIT_FAILURE);
      }
      SetHstFace( hst );
      hst->Sumw2(true);
      hists.push_back(hst);
   }

   return hists;
}

// {{{1 plot ratio deleted tracks
//--------------------------------------------------------------------
void TrkDel(string fname, string lhead, string pdf) {
//--------------------------------------------------------------------
   vector<TH1D*> hst = get_hists(fname, {"ct_Ntrk","ct_nrm"});

   string nc = pdf.substr(4, pdf.find(".")-4);
   TCanvas* c1 = new TCanvas(nc.c_str(),nc.c_str(),0,0,800,720);
   c1->cd();
   c1->Divide(1,2,0.01,0.);

   c1->cd(1);
   gPad->SetGrid();
   gPad->SetLogy(true);

   hst[0]->SetTitle(";;Events");
   hst[0]->GetYaxis()->SetMaxDigits(3);
   hst[0]->GetYaxis()->CenterTitle(true);
   hst[0]->GetYaxis()->SetTitleSize(0.05);
   hst[0]->GetYaxis()->SetTitleOffset(0.9);
   hst[0]->SetMinimum(1.);
   hst[0]->SetLineWidth(2);
   hst[0]->Draw("HIST");

   hst[1]->SetLineColor(kRed+2);
   hst[1]->SetFillStyle(3001);
   hst[1]->SetFillColor(kRed-9);
   hst[1]->SetLineWidth(2);
   hst[1]->Draw("HIST SAME");

   c1->cd(2);
   gPad->SetLogy(true);
   gPad->SetGrid();

   TH1D* hr = (TH1D*) hst[1]->Clone( Form("rat_%s",nc.c_str()) );
   // hr->Divide(hst[1],hst[0], 100., 1.,"B");
   hr->Divide(hst[1],hst[0], 100., 1.);
   hr->SetMinimum(1e-2);
   hr->SetMaximum(20);
   // hr->SetTitle(";Negative and Positive tracks;Cloned / All (%)");
   hr->SetTitle(";Number of negative and positive tracks in event"
         ";Cloned / All (%)");
   hr->GetYaxis()->CenterTitle(true);
   hr->GetYaxis()->SetTitleSize(0.05);
   hr->GetYaxis()->SetTitleOffset(0.9);
   hr->GetXaxis()->SetNdivisions(20);
   hr->GetXaxis()->CenterTitle(true);
   hr->GetXaxis()->SetTitleSize(0.05);
   // hr->GetXaxis()->SetTitleOffset(1.);
   hr->SetLineColor(kRed+2);
   hr->SetLineWidth(2);
   hr->Draw("E1");

   // TLegend* leg = new TLegend(0.38,0.75,0.71,0.98);
   TLegend* leg = new TLegend(0.40,0.70,0.70,0.98);
   leg->SetHeader( lhead.c_str(),"C");
   leg->AddEntry(hst[0],"All tracks","L");
   leg->AddEntry(hst[1],"Cloned tracks","F");
   leg->AddEntry(hr,"Ratio Cloned/All","LE");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print( pdf.c_str() );
   }
}

// {{{1 plot angle
//--------------------------------------------------------------------
void ClnAng(string fname, string lhead, string pdf) {
//--------------------------------------------------------------------
   vector<TH1D*> hst = get_hists(fname, {"ct_ang"});

   string nc = pdf.substr(4, pdf.find(".")-4);
   TCanvas* c1 = new TCanvas(nc.c_str(),nc.c_str(),0,0,800,800);

   // box to show selection region:
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-10);
   box->SetLineWidth(1);

   c1->cd();
   gPad->SetGrid();

   hst[0]->SetTitle(";cos(#Theta);");
   hst[0]->GetYaxis()->SetMaxDigits(3);
   hst[0]->GetXaxis()->SetNdivisions(-405);
   double ymax = 1.05 * hst[0]->GetMaximum();
   hst[0]->SetMaximum(ymax);
   hst[0]->SetLineWidth(2);
   hst[0]->Draw("HIST");

   box->DrawBox(0.9998,0.,1.,ymax);

   TLegend* leg = new TLegend(0.15,0.78,0.60,0.88);
   leg->AddEntry(hst[0],lhead.c_str(),"L");
   leg->AddEntry(box,"Clone search area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 plot momentum distanse
//--------------------------------------------------------------------
void Cln_dP(string fname, string lhead, string pdf) {
//--------------------------------------------------------------------
   vector<TH1D*> hst = get_hists(fname, {"ct_dpp"});

   string nc = pdf.substr(4, pdf.find(".")-4);
   TCanvas* c1 = new TCanvas(nc.c_str(),nc.c_str(),0,0,800,800);

   // box to show selection region:
   TBox* box = new TBox;
   box->SetFillStyle(3001);
   box->SetFillColor(kRed-10);
   box->SetLineColor(kRed-10);
   box->SetLineWidth(1);

   c1->cd();
   gPad->SetGrid();

   hst[0]->SetTitle(";momentum difference, GeV/c;");
   hst[0]->GetYaxis()->SetMaxDigits(3);
   double ymax = 1.05 * hst[0]->GetMaximum();
   hst[0]->SetMaximum(ymax);
   hst[0]->SetLineWidth(2);
   hst[0]->Draw("HIST");

   box->DrawBox(-0.005,0.,0.005,ymax);

   TLegend* leg = new TLegend(0.54,0.79,0.89,0.89);
   leg->AddEntry(hst[0],lhead.c_str(),"L");
   leg->AddEntry(box,"Clone search area","F");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 plot number found in event
//--------------------------------------------------------------------
void ClnNcl(string fname, string lhead, string pdf) {
//--------------------------------------------------------------------
   vector<TH1D*> hst = get_hists(fname, {"ct_ncl"});

   string nc = pdf.substr(4, pdf.find(".")-4);
   TCanvas* c1 = new TCanvas(nc.c_str(),nc.c_str(),0,0,800,800);

   c1->cd();
   gPad->SetGrid();
   gPad->SetLogy(true);

   hst[0]->SetTitle(";number of cloned tracks in one event;");
   hst[0]->GetYaxis()->SetMaxDigits(3);
   hst[0]->SetLineWidth(2);
   hst[0]->Draw("HIST");

   TLegend* leg = new TLegend(0.54,0.79,0.89,0.89);
   leg->AddEntry(hst[0],lhead.c_str(),"L");
   leg->Draw();

   gPad->RedrawAxis();
   c1->Update();
   if ( !pdf.empty() ) {
      c1->Print(pdf.c_str());
   }
}

// {{{1 Main
//--------------------------------------------------------------------
void ClonedTrks() {
//--------------------------------------------------------------------
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetLegendFont(42);

   //========================================================
   // set the name of the folder with the root files
   const string Dir = "prod_v709n3/";
   //========================================================

   for ( int date : {2009,2012,2021} ) {
   // for ( int date : {2021} ) {
      // for ( int type : { 0, 1, 2 } ) { // data, mcinc, mcsig
      for ( int type : { 0, 1 } ) {
         string fname, lhead, pdf;
         string pdfA,pdfP,pdfN;
         switch ( type ) {
            case 0:
               fname = string(
                     Form("data_%02ipsip_all.root",date%100) );
               lhead = string( Form("Date %i",date) );
               pdf = string( Form("Cln_data_%02i.pdf",date%100) );
               pdfA = string( Form("Ang_data_%02i.pdf",date%100) );
               pdfP = string( Form("Mom_data_%02i.pdf",date%100) );
               pdfN = string( Form("Ncl_data_%02i.pdf",date%100) );
               break;
            case 1:
               fname = string(
                     Form("mcinc_%02ipsip_all.root",date%100) );
               lhead = string( Form("Inclusive MC %i",date) );
               pdf = string( Form("Cln_mcinc_%02i.pdf",date%100) );
               pdfA = string( Form("Ang_mcinc_%02i.pdf",date%100) );
               pdfP = string( Form("Mom_mcinc_%02i.pdf",date%100) );
               pdfN = string( Form("Ncl_mcinc_%02i.pdf",date%100) );
               break;
            case 2:
               fname = string(
                     Form("mcsig_kkmc_%02i.root",date%100) );
               lhead = string( Form("Signal MC %i",date) );
               pdf = string( Form("Cln_mcsig_%02i.pdf",date%100) );
               pdfA = string( Form("Ang_mcsig_%02i.pdf",date%100) );
               pdfP = string( Form("Mom_mcsig_%02i.pdf",date%100) );
               pdfN = string( Form("Ncl_mcsig_%02i.pdf",date%100) );
               break;
         }


         TrkDel(Dir+fname,lhead,pdf);
         // ClnAng(Dir+fname,lhead,pdfA);
         // Cln_dP(Dir+fname,lhead,pdfP);
         // ClnNcl(Dir+fname,lhead,pdfN);
      }
   }

}
