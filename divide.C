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

#define FILES_NUMBER 3
#define CHAR_SIZE 100

void divide()
{
  std::string files[]  = {"testando_clustering_Cs_1_300000_transv_keV_hybrid", "transv_cfilme_15_th10keV", "transv_cfilme_075_th10keV_sigmax5"};

  std::string legend[] = {"Dados", "Sim_10keV", "Sim_10keV_sigma"};

  //int colors[]={1,632,416,600,920,400,616,432,800,15,215,840,7};//860 2 900
  int colors[]={1,2,3,4,5,6,7,8,9,11,41,42,49,28,32,40};

  int rebin[]={20,20,20,20,20,20,20}; // 20 allnpixels e 40 para multipixel

  std::string hname[]     = {"multipixelhistogram", "singlepixelhistogram"};
  // std::string hname[]     = {"singlepixelhistogram", "multipixelhistogram"}; // x*log([0]*x - [1])
  std::string fileoutput  = "divide";

  TH1F *multipixelhist[FILES_NUMBER];
  TH1F *singlepixelhist[FILES_NUMBER];

  char *histname = new char[CHAR_SIZE];
  char *legd = new char[CHAR_SIZE];

  std::string filename;
  std::string path;

  for(int m=0; m<FILES_NUMBER ; m+=1)
  {
    filename = files[m];
    path = "../"+filename+".root";

    cout << path << "\n";

    TFile *f1 = TFile::Open((path).c_str());

    sprintf(histname,"%s", (files[m]).c_str());
    sprintf(legd,"%s", (legend[m]).c_str());

    multipixelhist[m] = (TH1F*)f1->Get( ( hname[0] ).c_str() )->Clone();
    multipixelhist[m]->SetName(legd);
    multipixelhist[m]->SetTitle(legd);

    singlepixelhist[m] = (TH1F*)f1->Get( ( hname[1] ).c_str() )->Clone();
    singlepixelhist[m]->SetName(legd);
    singlepixelhist[m]->SetTitle(legd);

    // TH1F* divide   = (TH1F*)f1->Get("multipixelhistogram");
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

  multipixelhist[0]->SetLineColor(colors[0]);
  multipixelhist[0]->SetMarkerColor(colors[0]);
  multipixelhist[0]->GetXaxis()->SetDecimals();
  multipixelhist[0]->GetXaxis()->SetTitleSize(0.055);
  multipixelhist[0]->GetYaxis()->SetTitleSize(0.055);
  multipixelhist[0]->GetXaxis()->CenterTitle();
  multipixelhist[0]->GetYaxis()->CenterTitle();
  multipixelhist[0]->GetXaxis()->SetTitle("keV");
  multipixelhist[0]->GetYaxis()->SetTitle("Contagens");
  multipixelhist[0]->GetYaxis()->SetTitleOffset(0.90);
  multipixelhist[0]->GetXaxis()->SetTitleOffset(0.95);
  multipixelhist[0]->GetXaxis()->SetLabelFont(42);
  multipixelhist[0]->GetYaxis()->SetLabelFont(42);
  multipixelhist[0]->GetXaxis()->SetLabelSize(0.035);
  multipixelhist[0]->GetYaxis()->SetLabelSize(0.035);

  multipixelhist[0]->Rebin(rebin[0]);
  //multipixelhist[0]->Scale(1/(multipixelhist[0]->Integral()));
  multipixelhist[0]->GetXaxis()->SetRangeUser(0,100);

  singlepixelhist[0]->Rebin(rebin[0]);
  //singlepixelhist[0]->Scale(1/(singlepixelhist[0]->Integral()));
  singlepixelhist[0]->GetXaxis()->SetRangeUser(0,100);

  multipixelhist[0]-> Divide(singlepixelhist[0]);
  multipixelhist[0]->Scale(1/(multipixelhist[0]->Integral()));
  multipixelhist[0]->Write();
  multipixelhist[0]->Draw();

  TLegend *leg = new TLegend(0.15,0.8,0.45,0.9); //(0.75,0.75,0.9,0.9); //(0.85,0.6,0.65,0.9); //(0.25 ou 0.75,0.4,0.1 ou 0.9,0.9);
  leg->SetTextSize(0.035);
  leg->AddEntry(multipixelhist[0],"","lpf");

  for(int n=1; n<FILES_NUMBER ; n+=1)
  {
    multipixelhist[n]->SetLineColor(colors[n]);
    multipixelhist[n]->SetMarkerColor(colors[n]);

    multipixelhist[n]->Rebin(rebin[n]);
    //multipixelhist[n]->Scale(1/(multipixelhist[n]->Integral()));
    multipixelhist[n]->GetXaxis()->SetRangeUser(0,100);

    singlepixelhist[n]->Rebin(rebin[0]);
    //singlepixelhist[n]->Scale(1/(singlepixelhist[n]->Integral()));
    singlepixelhist[n]->GetXaxis()->SetRangeUser(0,100);

    multipixelhist[n]-> Divide(singlepixelhist[n]);
    multipixelhist[n]->Scale(1/(multipixelhist[n]->Integral()));
    multipixelhist[n]->Write();
    multipixelhist[n]->Draw("same");

    leg->AddEntry(multipixelhist[n],"","lpf");
  }

  leg->Draw();

  c1->Print( ( fileoutput + ".pdf" ).c_str() );
  c1->Write();

  c1->SetLogy();
  c1->Print( ( fileoutput + "_log.pdf" ).c_str() );
  //c1->Close();

  dfile->Write();
  //dfile->Close();
}
