#include <TFile.h>
#include "Riostream.h"
#include "TMath.h"
#include <vector>
#include <string>

#define MAIN_FILES_NUMBER 1
#define FILES_NUMBER 4

void sum_data(){

  std::string num[] = {"1","2","3","4"};

  int file=0;

  for(int i=0; i < MAIN_FILES_NUMBER; i++)
  {

    TFile *f2                  = new TFile( "clustering_Ba133_all.root", "RECREATE" );
    TH1F *singlepixelhistogram = new TH1F("singlepixelhistogram","",10000, 0, 10000);
    TH1F *multipixelhistogram  = new TH1F("multipixelhistogram","",10000, 0, 10000);
    TH1F *allnpixelshistogram  = new TH1F("allnpixelshistogram","",10000, 0, 10000);
    TH1F *clustersizehistogram = new TH1F("clustersizehistogram","",100, 0, 100);
    TH1F *clusterperframehisto = new TH1F("clusterperframehisto","",400, 0, 400);

    while(file < FILES_NUMBER)
    {
      file +=1;
      cout << " " << i  << " " << file  << endl;

      TFile* f1 = new TFile(("clustering_Ba133_" + num[file - 1] + ".root").c_str());

      TH1F* h0   = (TH1F*)f1->Get("singlepixelhistogram");
      TH1F* h1   = (TH1F*)f1->Get("multipixelhistogram");
      TH1F* h2   = (TH1F*)f1->Get("allnpixelshistogram");
      TH1F* h3   = (TH1F*)f1->Get("clustersizehistogram");
      TH1F* h4   = (TH1F*)f1->Get("clusterperframehisto");

      singlepixelhistogram -> Add(h0,1);
      multipixelhistogram  -> Add(h1,1);
      allnpixelshistogram  -> Add(h2,1);
      clustersizehistogram -> Add(h3,1);
      clusterperframehisto -> Add(h4,1);

      f1->Close();
    }

    file = 0;

    f2 -> Write();
    f2 -> Close();

  }

}
