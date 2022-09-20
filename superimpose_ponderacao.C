#include <TFile.h>
#include <TObject.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <string>//
#include <vector>//
#include <Riostream.h>//
#include <sstream>//
#include <iostream>
#include <fstream>

#define FILES_NUMBER 4
#define CHAR_SIZE 100

void superimpose_allnpixelshistogram_transv()
{
  std::string files[]  = {"tab_chi2_7_clustering_Cs_1_300000_transv_keV", "final_transv_15_th12keV_sigmax75_range50_Hill_after", "final_transv_075_th12keV_sigmax75_range50_Hill_after", "final_transv_0254_th12keV_sigmax75_range50_Hill_after"};

  std::string legend[] = {"Tomada 2", "1,5 mm", "0,75 mm", "0,254 mm"};

  //int colors[]={1,632,416,600,920,400,616,432,800,15,215,840,7};//860 2 900
  int colors[]={1,2,3,4,5,6,7,8,9,11,41,42,49,28,32,40};

  int rebin[]={20,20,20,20,20,20,20}; // 20 allnpixels e 40 para allnpixels

  std::string hname0[]     = {"allnpixelshistogram", "singlepixelhistogram", "singlepixelhistogram", "singlepixelhistogram"};
  std::string hname1[]     = {"allnpixelshistogram", "multipixelhistogram", "multipixelhistogram", "multipixelhistogram"};
  std::string fileoutput  = "allnpixelshistogram";

  TH1F *myhist[FILES_NUMBER];
  TH1F *myhist1[FILES_NUMBER];

  TH1F *refSP[FILES_NUMBER];
  TH1F *refMP[FILES_NUMBER];

  char *histname = new char[CHAR_SIZE];
  char *legd = new char[CHAR_SIZE];

  double hs = 0.0, hd = 0.0;

  double totalSPMP = 0.0, pSP = 0.0, pMP = 0.0;

  std::string filename;
  std::string path;

  for(int m=0; m<FILES_NUMBER ; m+=1)
  {
    filename = files[m];
    path = filename+".root";

    cout << path << "\n";

    TFile *f1 = TFile::Open((path).c_str());

    sprintf(histname,"%s", (files[m]).c_str());
    sprintf(legd,"%s", (legend[m]).c_str());

    myhist[m] = (TH1F*)f1->Get( ( hname0[m] ).c_str() )->Clone();
    myhist[m]->SetName(histname);
    myhist[m]->SetTitle(legd);

    myhist1[m] = (TH1F*)f1->Get( ( hname1[m] ).c_str() )->Clone();
    myhist1[m]->SetName(histname);
    myhist1[m]->SetTitle(legd);

    refSP[m] = (TH1F*)f1->Get( "singlepixelhistogram" )->Clone();
    refSP[m]->SetName("singlepixelhistogram");
    refSP[m]->SetTitle("singlepixelhistogram");

    refMP[m] = (TH1F*)f1->Get( "multipixelhistogram" )->Clone();
    refMP[m]->SetName("multipixelhistogram");
    refMP[m]->SetTitle("multipixelhistogram");
  }

  //f1->Close();
  //delete f1;

  TFile *dfile = new TFile( ( fileoutput + ".root" ).c_str(), "RECREATE");

  TCanvas *c1 = new TCanvas("c1"," ",0,0,600,400);
  c1->GetFrame()->SetFillColor(21);
  c1->cd();
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  myhist[0]->SetLineColor(colors[0]);
  myhist[0]->SetMarkerColor(colors[0]);
  myhist[0]->GetXaxis()->SetDecimals();
  myhist[0]->GetXaxis()->SetTitleSize(0.055);
  myhist[0]->GetYaxis()->SetTitleSize(0.055);
  myhist[0]->GetXaxis()->CenterTitle();
  myhist[0]->GetYaxis()->CenterTitle();
  myhist[0]->GetXaxis()->SetTitle("Energia (keV)");
  myhist[0]->GetYaxis()->SetTitle("Frequ#hat{e}ncia");
  myhist[0]->GetYaxis()->SetTitleOffset(0.90);
  myhist[0]->GetXaxis()->SetTitleOffset(0.95);
  myhist[0]->GetXaxis()->SetLabelFont(42);
  myhist[0]->GetYaxis()->SetLabelFont(42);
  myhist[0]->GetXaxis()->SetLabelSize(0.035);
  myhist[0]->GetYaxis()->SetLabelSize(0.035);
  myhist[0]->Rebin(rebin[0]);
  hs = myhist[0]->GetBinContent(25);
  myhist[0]->Scale(1/(myhist[0]->Integral()));
  myhist[0]->GetXaxis()->SetRangeUser(0,300);
  myhist[0]->GetYaxis()->SetRangeUser(0, 0.07);
  myhist[0]->Write();
  myhist[0]->Draw();

  totalSPMP = refSP[0]->Integral() + refMP[0]->Integral();
  pSP = (refSP[0]->Integral())/totalSPMP;
  pMP = (refMP[0]->Integral())/totalSPMP;

  TLegend *leg = new TLegend(0.85,0.6,0.75,0.9); //(0.75,0.75,0.9,0.9); //(0.85,0.6,0.65,0.9); //(0.25 ou 0.75,0.4,0.1 ou 0.9,0.9);
  leg->SetTextSize(0.035);
  leg->AddEntry(myhist[0],"","lpf");

  for(int n=1; n<FILES_NUMBER ; n+=1)
  {
    myhist[n]->Scale(pSP/(myhist[n]->Integral()));
    myhist1[n]->Scale(pMP/(myhist1[n]->Integral()));
    myhist[n]->Add(myhist1[n]);
    // myhist[n]->Scale(1/(myhist[n]->Integral()));
    // myhist[n]->SetLineWidth(n/2);
    myhist[n]->SetLineColor(colors[n]);
    myhist[n]->SetMarkerColor(colors[n]);
    //myhist[n]->SetLineStyle(colors[n]);
    myhist[n]->Write();
    myhist[n]->Rebin(rebin[n]);
    hd = (myhist[n]->GetBinContent(25));//*1.25;
    myhist[n]->Scale(1/(myhist[n]->Integral()));
    // myhist[n]->Scale(hs/hd);
    myhist[n]->GetXaxis()->SetRangeUser(0,400);
    // myhist[n]->GetYaxis()->SetRangeUser(0, 0.1);
    myhist[n]->Draw("same");
    leg->AddEntry(myhist[n],"","lpf");
  }

  leg->Draw();

  c1->Print( ( fileoutput + ".pdf" ).c_str() );
  c1->Write();

  // c1->SetLogy();
  // c1->Print( ( fileoutput + "_log.pdf" ).c_str() );

  //c1->Close();

  dfile->Write();
  //dfile->Close();
}
