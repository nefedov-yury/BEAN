// calculate the uncertainties of efficiency because of ISR
//    -> xisr_un_${uncer}.pdf

//----------------------------------------------------------
void get_histo(string fname, TH1D* ret[2])
//----------------------------------------------------------
{
#include "cuts.h"

//   cout << " file: " << fname << endl;
  TFile* froot = TFile::Open(fname.c_str(),"READ");
  if( froot == 0 ) {
    cerr << "can not open " << fname << endl;
    exit(0);
  }

  froot->cd("SelectKKgg");

  TH1D* Xisr_ini = (TH1D*)gROOT->FindObject("mcisr_1x");
  Xisr_ini->Sumw2();
  if( !Xisr_ini ) {
    cerr << " can not find Xisr_ini" << endl;
    exit(0);
  }

  TTree* a4c = (TTree*)gDirectory->Get("a4c");

  // see "cuts.h"

  TCut c_tot = c_chi2 && "cp"; // use TCutG by name

  double Xisr_min = 0.9; // <----------------------------------
  // XISR cut:
  TCut c_xisr(Form("xisr>%f",Xisr_min));
  c_tot += c_xisr;

  TH1D* xf = new TH1D("xf", "after selection;1-s'/s", 100,0.,1.);
  xf->Sumw2();
  a4c->Draw("(1.-xisr)>>xf",c_tot,"goff");

  ret[0] = Xisr_ini;
  ret[1] = xf;
}

//----------------------------------------------------------
bool get_rateff(string name, string dir, TCanvas* c1, double rat[2])
//----------------------------------------------------------
{
//   cout << " name: " << name << endl;
  string fn = "ntpl_mcgpj_" + name + ".root";
  TH1D* ret[2];
  get_histo(fn, ret);
  TH1D* Xini = ret[0];
  TH1D* Xfin = ret[1];
  string tl = "ISR from MCGPJ;1 - s'/s";
  Xini->SetTitle(tl.c_str());

  string fname = dir + "/hist_phieta_" + name + "mc.root";
  TFile* fr1 = TFile::Open(fname.c_str(),"READ");
  if( fr1 == 0 ) {
    cerr << "can not open file: " << fname << endl;
    return false;
  }
  fr1->cd();
//   TH1D* x2 = (TH1D*)Xini->Clone("h_x"); // DEBUG
  TH1D* x2 = (TH1D*)gROOT->FindObject("h_x");
//   x2->Sumw2();
//   x2->Scale(0.01);
  string tit = "ratio;1 - s'/s";
  x2->SetTitle(tit.c_str());


  c1->cd(1);
  Xini->SetAxisRange(1.,2.e5,"Y");
  Xini->SetLineColor(kRed-7);
  Xini->SetFillColor(kRed-7); // for "E" options
  Xini->DrawCopy("E3");

  x2->SetLineColor(kBlue-7);
  x2->SetFillColor(kBlue-7); // for "E" options
  x2->DrawCopy("E3,SAME");

  TLatex* t1 = new TLatex(0.1,4e4, Form("0 - %s",fn.c_str()) );
  t1->SetTextColor(kRed);
  t1->Draw();
  TLatex* t2 = new TLatex(0.1,1e4, Form("1 - %s",fname.c_str()) );
  t2->SetTextColor(kBlue);
  t2->Draw();

  c1->cd(2);
  x2->Divide(x2,Xini,1.,1.);
  x2->SetAxisRange(0.,0.2,"X");
  x2->SetAxisRange(0.5,1.5,"Y");
  x2->SetLineColor(kBlue+1);
  x2->DrawCopy("E1");

  // efficiency
  int nb = Xini->GetNbinsX();
  double Nini = Xini->Integral();
//   double eff0 = Xfin->Integral()/Nini;

  double eff0 = 0, err0 = 0;
  double eff1 = 0, err1 = 0;
  for(int i = 0; i <= nb+1; i++ ) {
//     eff1 += Xfin->GetBinContent(i)*x2->GetBinContent(i);
     double a  = Xfin->GetBinContent(i);
     double b  = x2->GetBinContent(i);
     if( a == 0 || b == 0 ) continue;
     double ea = Xfin->GetBinError(i);
     double eb = x2->GetBinError(i);
     eff0 += a;
     err0 += ea*ea;
     eff1 += a*b;
     err1 += a*a*eb*eb + b*b*ea*ea;
  }
  eff0 /= Nini;
  err0  = sqrt(err0)/Nini;
  eff1 /= Nini;
  err1  = sqrt(err1)/Nini;

  double ratio = eff1/eff0;
  double err_r = ratio*sqrt( (err0*err0)/(eff0*eff0) +
                             (err1*err1)/(eff1*eff1) );

//   printf(" eff0= %.2f +/- %.2f %%", eff0*100,err0*100);
//   printf(" eff1= %.2f +/- %.2f %%\n", eff1*100,err1*100);
//   printf(" 1-eff1/eff0= %.2f +/- %.2f %%\n", (1-ratio)*100, err_r*100);

  TLatex* t3 = new TLatex(0.02,1.4,
      Form("1 - #varepsilon_{1} / #varepsilon_{0} = %.2f #pm %.2f %%",
            (1-ratio)*100, err_r*100)
                         );
  t3->SetTextColor(kMagenta+1);
  t3->Draw();

  rat[0] = ratio;
  rat[1] = err_r;

  return true;
}

//----------------------------------------------------------
 TGraphErrors* xisr_un_xx(int sw_dir, string pdf) // sw_dir = 0,1,2,3
//----------------------------------------------------------
{
  TCanvas* c1 = new TCanvas("c1","...",0,0,600,600);
  c1->Divide(1,2);
  c1->cd(1);
  gPad->SetGrid();
  gPad->SetLogy();
  c1->cd(2);
  gPad->SetGrid();

  // 3080_new MUST be last
  string names[] = {
                 "3050", "3060", "3080", "3083", "3090",
                 "3093", "3094", "3095", "3096", "3097",
                 "3098", "3099", "3102", "3106", "3112",
                 "3120",
                 "3080_new"
                   };
  double ebeam[] = {
                 3050.213, 3059.257, 3080.195, 3083.060, 3089.418,
                 3092.324, 3095.261, 3095.994, 3096.390, 3097.777,
                 3098.904, 3099.606, 3101.923, 3106.144, 3112.615,
                 3120.442,
                 3080.
                   };

  // ATTENTION: MC only
  int Nbeam = 16;
  vector<double> Ebeam(Nbeam);
  for(int i=0; i < Nbeam; ++i) {
     Ebeam[i] = ebeam[i] - 0.55;
  }

  vector<double> err_e(Nbeam,0.01);
  vector<double> ratio(Nbeam,0.0);
  vector<double> err_r(Nbeam,0.0);
  string serr("  double err_mc[] = {");

  string Mdir = "./";    // <-------------------------
//   int sw_dir = 3;

//   string Fdir[] = {"un_am", "un_ap", "un_pm", "un_pp_pi"};
  string Fdir[] = {"un_am", "un_ap", "un_pm", "un_pp"};
  string Tit[] = {"#varepsilon(A - #Delta A)",
                  "#varepsilon(A + #Delta A)",
                  "#varepsilon(#varphi - #Delta #varphi)",
                  "#varepsilon(#varphi + #Delta #varphi)" };

  string fdir = Mdir + Fdir[sw_dir];
  for(int i = 0; i < Nbeam; ++i) {

    double rat[2]={1.,0.};
    bool ok = get_rateff(names[i],fdir,c1,rat);
    ratio[i] = 100.*(1. - rat[0]); // -> %
    err_r[i] = 100.*rat[1];
    serr += string(Form(" %.2f,",fabs(ratio[i])));
    if( (i+1)%5 == 0 ) serr += string(1,'\n')+string(21,' ');
    if( !ok ) continue;

    c1->Update();
    if( pdf.size() > 0 ) c1->Print(pdf.c_str()); // add to pdf-file
  }
  // add 3080 new in serr and print out
  serr += string(Form(" %.2f };",fabs(ratio[2])));
  cout << serr << endl;

  TGraphErrors* gr = new TGraphErrors(Nbeam,
      Ebeam.data(),ratio.data(),err_e.data(),err_r.data());
  string tit = Tit[sw_dir] +
               string(";center of mass energy (MeV)") +
               string(";relative difference (%)");
  gr->SetTitle(tit.c_str());
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.8);
  gr->SetMarkerColor(kMagenta+1);
  gr->SetLineColor(kMagenta+1);
  gr->SetMinimum(-5.);
  gr->SetMaximum(+5.);
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleOffset(1.1);
  gr->GetYaxis()->SetTitleOffset(1.1);
  gr->GetXaxis()->SetLabelFont(62);
  gr->GetYaxis()->SetLabelFont(62);
  gr->GetXaxis()->SetLabelSize(0.03);
  gr->GetYaxis()->SetLabelSize(0.03);

  c1->cd();
  c1->Clear();
  gPad->SetGrid();
  gr->Draw("AP");
  c1->Update();
  if( pdf.size() > 0 ) c1->Print(pdf.c_str()); // add to pdf-file

  delete c1;

  return gr;
}

//----------------------------------------------------------
void xisr_uncer()
//----------------------------------------------------------
{
  gROOT->Reset();

  gStyle->SetOptStat(0);
//   gStyle->SetOptStat(111110);
//   gStyle->SetOptFit(1);

//   gStyle->SetStatX(0.35);
//   gStyle->SetStatY(0.85);

  int debug = 0;    // <-------------------------
  string pdf("xisr_uncer.pdf");
  string pdf_zero; // do not print intermediate results

  TCanvas* c = new TCanvas("c","...",0,0,600,600);
  c->Print((pdf+"[").c_str()); // just open pdf-file

  TGraphErrors* gr[4];
  for(int i = 0; i < 4; ++i) {
    if( debug == 1 )
      gr[i] = xisr_un_xx(i,pdf);
    else
      gr[i] = xisr_un_xx(i,pdf_zero);
  }

  c->cd();
  gPad->SetGrid();


//   gr[2]->SetTitle("");
  gr[2]->SetTitle("uncertainty in efficiency");
  gr[2]->SetMaximum(+9.);
  gr[2]->SetMarkerColor(kRed+2);
  gr[2]->SetLineColor(kRed+2);
  gr[2]->Draw("AP");

//   gr[3]->SetMarkerColor(kMagenta+2);
//   gr[3]->SetLineColor(kMagenta+2);
  gr[3]->SetMarkerColor(kRed-6);
  gr[3]->SetLineColor(kRed-6);
  gr[3]->SetMarkerStyle(21);
  gr[3]->Draw("P");

  gr[0]->SetMarkerColor(kGreen+3);
  gr[0]->SetLineColor(kGreen+3);
  gr[0]->SetMarkerStyle(22);
  gr[0]->Draw("P");

  gr[1]->SetMarkerColor(kGreen-5);
  gr[1]->SetLineColor(kGreen-5);
  gr[1]->SetMarkerStyle(23);
  gr[1]->Draw("P");

  gr[2]->Draw("P");

//   TLegend* leg = new TLegend(0.58,0.64,0.93,0.92);
  TLegend* leg = new TLegend(0.14,0.67,0.86,0.89);
  leg-> SetNColumns(2);
//   leg->SetHeader("uncertainty in efficiency due to ISR","C");
  leg->AddEntry(gr[2],"#varepsilon(#varphi - #Delta #varphi)", "PE");
  leg->AddEntry(gr[0],"#varepsilon(A - #Delta A)", "PE");
  leg->AddEntry(gr[3],"#varepsilon(#varphi + #Delta #varphi)", "PE");
  leg->AddEntry(gr[1],"#varepsilon(A + #Delta A)", "PE");
  leg->Draw();

  c->Update();
  c->Print(pdf.c_str()); // add to pdf-file
  c->Print((pdf+"]").c_str()); // just close pdf-file
}
