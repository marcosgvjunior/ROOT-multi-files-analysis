//Institute of Physics - Federal University of Rio de Janeiro
//Interdisciplinary Academic Master's Degree in Applied Physics
//Student: Marcos Vieira
//June, 2019
//Working on ROOT 5.34/36, Windows 10 x64

//Script to find several ToT peaks corresponding to different energies and sources
//Find a calibration curve and study the ralation sigma x energy

#include <TFile.h>
#include <TObject.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <string>//
#include <vector>//
#include <Riostream.h>//
#include <sstream>//
#include <iostream>
#include <fstream>
#include <algorithm>

// #include <TRandom.h>
// #include <TRandom3.h>

// Some references:
// https://root.cern.ch/root/html/tutorials/fit/multifit.C.html
// https://root.cern.ch/doc/master/FittingDemo_8C_source.html
// https://stackoverflow.com/questions/3909272/sorting-two-corresponding-arrays

//
class sort_index
{
   private:
     Float_t* mparr;
   public:
     sort_index(Float_t* parr) : mparr(parr) {}
     bool operator()(int i, int j) const { return mparr[i]<mparr[j]; }
};

// Code parameters
#define FILES_NUMBER          4

#define PEAKS_NUMBER_FILE0    5
#define PEAKS_NUMBER_FILE1    1
#define PEAKS_NUMBER_FILE2    1
#define PEAKS_NUMBER_FILE3    1

#define PEAKS_NUMBER          8

#define CHAR_SIZE             100

#define AM_MEAN_MODEL_0      13.81
#define AM_MEAN_MODEL_1      17.70
#define AM_MEAN_MODEL_2      20.70
#define AM_MEAN_MODEL_3      26.34
#define AM_MEAN_MODEL_4      59.54

#define BA_MEAN_MODEL_5       80.99

#define CO_MEAN_MODEL_6       14.41
#define CO_MEAN_MODEL_7       122.06

#define REBIN_FLAG            true
#define REBIN_VALUE           4

//Initial guess values to set the limits of a energy peak and get its maximum (bin of maximum value)
#define ANG_COEF_GUESS        10
#define LIN_COEF_GUESS        40

#define MAX_BIN_INITIAL_VALUE 3000
#define MAX_BIN_VIEW_VALUE    1000

#define ZOOM                  true

#define REPS_NUMBER          5

void fit_multi_peaks_co()
{
  // Force batch mode
  gROOT -> SetBatch( kTRUE );
  // gROOT -> ProcessLine(".x lhcbStyle.C");

  // Long64_t seed = time(NULL); //add
  // gRandom = new TRandom3(seed); // seed

  std::string rep;

  std::string repvector[]     =  {"0", "1", "2", "3", "4"};

  //input files
  std::string files[]     =  {"clustering_Am241", "clustering_Ba133", "clustering_Co57", "clustering_Co57"};
  std::string filesO[]     =  {"clustering_Am241", "clustering_Ba133", "clustering_Co57_0", "clustering_Co57_1"};

  //Histograms legends/names
  std::string legend[]    =  { "Am-241", "Ba-133", "Co-57", "Co-57" };

  // User coordinates limits for each peak
  Int_t lim[FILES_NUMBER][PEAKS_NUMBER][2]     = {{{13.7,14},{16.5,18},{20.5,22},{25,27},{50,62}},
                                                  {{62.5,95.5}},
                                                  {{14,15}},
                                                  {{100,130}}};

  //Bin number increment to left/right side of each peak
  //These values ​​are used to establish the fit ranges (given the respectives bins of maximum values)
  Int_t somalim[FILES_NUMBER][PEAKS_NUMBER][2] = {{{8,12},{8,12},{5,23},{10,21},{14,62}},
                                                  {{35,90}},
                                                  {{8,12}},
                                                  {{70,150}}};


  Int_t rebinPeaks[FILES_NUMBER][PEAKS_NUMBER] = {{REBIN_VALUE, REBIN_VALUE, REBIN_VALUE, REBIN_VALUE, REBIN_VALUE},
                                                   {8},
                                                   {REBIN_VALUE},
                                                   {40}};

  // Name of the input histogram
  std::string histogramname = "singlepixelhistogram";

  // List of colors
  Int_t colors[] = { 1, 2, 3, 6, 7, 8, 9, 11, 41, 42, 49, 28, 32, 40 }; // 4, 5,

  // Initializing a list for each histogram and its names/legends
  TH1F *mhist[FILES_NUMBER];                TH1F *mhists[PEAKS_NUMBER];

  char *histname  =  new char[CHAR_SIZE];   char *peakname  =  new char[CHAR_SIZE];
  char *legd      =  new char[CHAR_SIZE];

  // Option to rebin or not
  bool rebin      =  REBIN_FLAG;

  // Variables to store filename and its path
  std::string fileprefix  = "";   std::string filename;    std::string filenameO;   std::string path;

  // Input folder path and path/folder to save the outputs
  //**Need to create the folder before run
  std::string relativepath   = "../";
  std::string subfoldername  = "fit_peaks/";

  Double_t par[3];
  Float_t peaks[] = { AM_MEAN_MODEL_0, AM_MEAN_MODEL_1, AM_MEAN_MODEL_2, AM_MEAN_MODEL_3,
                      AM_MEAN_MODEL_4, BA_MEAN_MODEL_5, CO_MEAN_MODEL_6, CO_MEAN_MODEL_7 };

  Int_t npeaks[]  = { PEAKS_NUMBER_FILE0, PEAKS_NUMBER_FILE1, PEAKS_NUMBER_FILE2, PEAKS_NUMBER_FILE3 };

  Float_t sigma[PEAKS_NUMBER][REPS_NUMBER];      Float_t sigmaerror[PEAKS_NUMBER][REPS_NUMBER];
  Float_t mean[PEAKS_NUMBER][REPS_NUMBER];       Float_t meanerror[PEAKS_NUMBER][REPS_NUMBER];

  Float_t Mmean[PEAKS_NUMBER];  Float_t Mmeanerror[PEAKS_NUMBER];
  Float_t Msigma[PEAKS_NUMBER]; Float_t Msigmaerror[PEAKS_NUMBER];
  Float_t Sdmean[PEAKS_NUMBER];
  Float_t Sdsigma[PEAKS_NUMBER];

  Int_t m = 0, n = 0;
  Int_t partialm = 0;

  //to sort the fits
  Int_t indx[PEAKS_NUMBER];

  for( Int_t reps = 0 ; reps <= 4; reps++)
  {

    m = 0;
    n = 0;
    partialm = 0;

    rep = repvector[reps];

    ///////////////////////////////// Beginning of loop over files and peaks ////////////////////////////////////////////
    for( Int_t l = 0; l < FILES_NUMBER ; l += 1 )
    {
      filename   =  files[l] + "_" + rep;
      filenameO  =  filesO[l];
      path       =  fileprefix + filename + ".root";

      TFile *f1 =  TFile::Open( ( path ).c_str() );

      sprintf( histname, "%s", legend[l]  ); //files[l] - error for long names
      sprintf( legd    , "%s", legend[l] );

      mhist[l]  =  ( TH1F* )f1 -> Get( ( histogramname ).c_str() ) -> Clone();
      mhist[l] -> SetName( histname );
      mhist[l] -> SetTitle( legd );

      TFile *dfile  =  new TFile( ( subfoldername + rep + "_" + filenameO + ".root" ).c_str(), "RECREATE" );

      partialm = m + npeaks[l];

      indx[m] = m;

      ///////////////////////////////// Beginning of loop over peaks ////////////////////////////////////////////////////
      while( m < partialm )
      {
        // There is probably a way to put the canvas ouside the loop and optimize
        TCanvas *c1  =  new TCanvas( "c1", " ", 300, 200, 300, 200 );

        c1 -> GetFrame() -> SetFillColor( 21 );      c1 -> cd();
        c1 -> SetBorderMode( 0 );                    c1 -> SetBorderSize( 0 ); //4
        c1 -> SetFrameBorderMode( 0 );

        gStyle -> SetOptStat( "n" );        gStyle -> SetOptTitle( 0 );
        gStyle -> SetOptFit( 1111 );        gStyle -> SetStatFont( 42 );
        gStyle -> SetStatY( 0.9 );          gStyle -> SetStatX( 0.9 );
        gStyle -> SetStatH( 0.2 );          gStyle -> SetStatW( 0.18 );
        gStyle -> SetStatFontSize( 0.04 );//gStyle -> SetLabelSize( 0.04 );

        if( rebin == true ) {  mhist[l] -> Rebin( rebinPeaks[l][n] ); rebin = false; }

        mhist[l]     -> GetXaxis() -> SetRangeUser( lim[l][n][0]*ANG_COEF_GUESS - LIN_COEF_GUESS,
                                                    lim[l][n][1]*ANG_COEF_GUESS - LIN_COEF_GUESS );
        Int_t maxbin  =  mhist[l]  -> GetMaximumBin();

        if( ZOOM )
        {
          if( m < 4  ) {  mhist[l] -> GetXaxis() -> SetRangeUser( 0, 300 ); }
          if( m == 4 ) {  mhist[l] -> GetXaxis() -> SetRangeUser( 0.6*(lim[l][n][0]*ANG_COEF_GUESS - LIN_COEF_GUESS),
                                                                 1.4*(lim[l][n][1]*ANG_COEF_GUESS - LIN_COEF_GUESS)); }
          if( m > 4 ) {  mhist[l] -> GetXaxis() -> SetRangeUser( 0.6*(lim[l][n][0]*ANG_COEF_GUESS - LIN_COEF_GUESS),
                                                                 1.4*(lim[l][n][1]*ANG_COEF_GUESS - LIN_COEF_GUESS)); }
        }
        else
        {
          mhist[l]     -> GetXaxis() -> SetRangeUser( 0, MAX_BIN_VIEW_VALUE );
        }

        mhist[l] -> SetLineColor( colors[l] );             mhist[l] -> SetMarkerColor( colors[l] );
        mhist[l] -> GetXaxis() -> SetDecimals();
        mhist[l] -> GetXaxis() -> SetTitleSize( 0.055 );   mhist[l] -> GetYaxis() -> SetTitleSize( 0.055 );
        mhist[l] -> GetXaxis() -> CenterTitle();           mhist[l] -> GetYaxis() -> CenterTitle();
        mhist[l] -> GetXaxis() -> SetTitle( "ToT" );       mhist[l] -> GetYaxis() -> SetTitle( "Contagens" );
        mhist[l] -> GetYaxis() -> SetTitleOffset( 0.95 );  mhist[l] -> GetXaxis() -> SetTitleOffset( 0.95 );
        mhist[l] -> GetXaxis() -> SetLabelFont( 42 );      mhist[l] -> GetYaxis() -> SetLabelFont( 42 );
        mhist[l] -> GetXaxis() -> SetLabelSize( 0.048 );   mhist[l] -> GetYaxis() -> SetLabelSize( 0.048 );
        //mhist[l] -> GetXaxis() -> SetLabelOffset( 1.2 ); //mhist[l] -> GetYaxis() -> SetLabelOffset( 1.2 );

        //
        TF1 *g1     =  new TF1( "g1", "gaus", rebinPeaks[l][n]*( maxbin ) - somalim[l][n][0],
                                              rebinPeaks[l][n]*( maxbin ) + somalim[l][n][1] );

        g1         -> SetLineColor( colors[m + 1] );
        g1         -> SetLineWidth( 4 ); //1
        g1         -> SetMarkerColor( colors[m + 1] );

        mhist[l]   -> Fit( g1, "QR" ); //R+

        g1         -> GetParameters( &par[0] ); //not used

        mean[m][reps]      = g1 -> GetParameter( 1 );      meanerror[m][reps]     = g1 -> GetParError( 1 );
        sigma[m][reps]     = g1 -> GetParameter( 2 );      sigmaerror[m][reps]    = g1 -> GetParError( 2 );

        mhist[l] -> SetStats( 1 );
        mhist[l] -> Write();
        mhist[l] -> Draw();

        mhists[m] = (TH1F*) mhist[l] -> Clone();

        sprintf( peakname, "%d", m );
        c1 -> Print( ( subfoldername + rep + "_" + filename + "_" + peakname + "_fit.pdf" ).c_str() );
        c1 -> Write();

        c1 -> SetLogy();
        c1 -> Print( ( subfoldername + rep + "_" + filename + "_" + peakname + "_log_fit.pdf" ).c_str() );
        c1 -> Close();

        //return histogram to wanted limits
        mhist[l]     -> GetXaxis() -> SetRangeUser( 0, MAX_BIN_VIEW_VALUE );

        m = m + 1;
        n = n + 1;

        indx[m] = m;
      }
      n     = 0;
      m     = partialm;
      rebin = REBIN_FLAG;

      f1 -> Close(); dfile -> Close();
    }
  }
  ///////////////////////////////// End of loop over files and peaks.  ///////////////////////////////////////////////
  //Histograms, means, sigmas etc stored and/or printed to files

  TFile *dfile  =  new TFile( ( subfoldername + "Calibration.root" ).c_str(), "RECREATE" );

  std::sort(indx, indx+PEAKS_NUMBER, sort_index(peaks));

  for( Int_t overpeaks = 0; overpeaks < PEAKS_NUMBER; overpeaks++ )
  {
    for( Int_t meancalc = 0; meancalc < REPS_NUMBER; meancalc++ )
    {
      Mmean[overpeaks]       += mean[overpeaks][meancalc];
      Mmeanerror[overpeaks]  += meanerror[overpeaks][meancalc];
      Msigma[overpeaks]      += sigma[overpeaks][meancalc];
      Msigmaerror[overpeaks] += sigmaerror[overpeaks][meancalc];
    }
    Mmean[overpeaks]       = Mmean[overpeaks]/5.0;
    Mmeanerror[overpeaks]  = Mmeanerror[overpeaks]/5.0;
    Msigma[overpeaks]      = Msigma[overpeaks]/5.0;
    Msigmaerror[overpeaks] = Msigmaerror[overpeaks]/5.0;

    cout << Mmean[overpeaks] << '\t' << Msigma[overpeaks] << '\n';
  }

  for( Int_t overpeaks = 0; overpeaks < PEAKS_NUMBER; overpeaks++ )
  {
    for( Int_t meancalc = 0; meancalc < REPS_NUMBER; meancalc++ )
    {
      Sdmean[overpeaks]  += ( mean[overpeaks][meancalc] - Mmean[overpeaks]  )**2 ;
      Sdsigma[overpeaks] += ( sigma[overpeaks][meancalc] - Msigma[overpeaks]  )**2 ;
    }
    Sdmean[overpeaks]  = sqrt( Sdmean[overpeaks]/4.0 ) ;
    Sdsigma[overpeaks] = sqrt( Sdsigma[overpeaks]/4.0 ) ;

    cout << Sdmean[overpeaks] << '\t' << Sdsigma[overpeaks] << '\n';
  }

  Float_t sortpeaks[PEAKS_NUMBER];

  Float_t sortmean[PEAKS_NUMBER];
  Float_t sortmeanerror[PEAKS_NUMBER];
  Float_t sortmeanSD[PEAKS_NUMBER];

  Float_t sortsigma[PEAKS_NUMBER];
  Float_t sortsigmaerror[PEAKS_NUMBER];
  Float_t sortsigmaSD[PEAKS_NUMBER];

  for (Int_t s = 0; s < PEAKS_NUMBER; s++ )
  {
    sortmean[s]        = Mmean[indx[s]];
    sortmeanerror[s]   = Mmeanerror[indx[s]];
    sortmeanSD[s]      = Sdmean[indx[s]];

    sortsigma[s]       = Msigma[indx[s]];
    sortsigmaerror[s]  = Msigmaerror[indx[s]];
    sortsigmaSD[s]     = Sdsigma[indx[s]];

    sortpeaks[s]       = peaks[indx[s]];

    cout << indx[s]         << "\n";
    cout << Mmean[indx[s]]  << "\n";
    cout << peaks[indx[s]]  << "\n";
  }

  TCanvas *c3  =  new TCanvas( "c3", " ", 0, 0, 600, 400 );
  c3 -> GetFrame() -> SetFillColor( 21 );  c3 -> cd();
  c3 -> SetBorderMode( 0 );                c3 -> SetBorderSize( 4 );
  c3 -> SetFrameBorderMode( 0 );

  gStyle -> SetOptStat( "n" );        gStyle -> SetOptTitle( 0 );
  gStyle -> SetStatY( 0.9 );          gStyle -> SetStatX( 0.5 );
  gStyle -> SetOptFit( 1 );

  // TF1 *fd0     =  new TF1( "f0", "( [0]*[3] + x - [1] + sqrt( ( -1*[1] + [0]*[3] + x )^2  + 4*[0]*( -1*[2] + [1]*[3] -1*[3]*x  ) ) )/( 2*[0] )", 10, 150 );

  // TF1 *fd1     =  new TF1( "f1", "[0]*x + [1] + ( [2]/( x - [3] ) )", 10, 150 );

  // // Using the sigma as the mean error
  // TGraphErrors *gr0  =  new TGraphErrors( PEAKS_NUMBER, sortpeaks, sortmean, 0, sortsigma ); // sortmeanerror
  // //TGraphErrors *gr0  =  new TGraphErrors( PEAKS_NUMBER, sortmean, sortpeaks, sortsigma, 0 );
  // gr0 -> Fit( "f1", "MER" );           // MER
  // gr0 -> SetTitle( "" );
  // gr0 -> SetMarkerColor( 4 );                             gr0 -> SetMarkerStyle( 20 );
  // gr0 -> GetXaxis() -> SetTitleSize( 0.055 );             gr0 -> GetYaxis() -> SetTitleSize( 0.055 );
  // gr0 -> GetYaxis() -> SetTitleOffset( 0.85 );            gr0 -> GetXaxis() -> SetTitleOffset( 0.85 ); // 90 e 95
  // gr0 -> GetXaxis() -> SetLabelFont( 42 );                gr0 -> GetYaxis() -> SetLabelFont( 42 );
  // gr0 -> GetXaxis() -> SetLabelSize( 0.035 );             gr0 -> GetYaxis() -> SetLabelSize( 0.035 );
  // gr0 -> GetYaxis() -> SetTitle("E_{ medido } ( ToT )");  gr0 -> GetXaxis() -> SetTitle("E_{ teorico } ( keV )");
  // gr0 -> GetXaxis() -> CenterTitle();                     gr0 -> GetYaxis() -> CenterTitle();
  // gr0 -> Draw( "AP" );

  // c3 -> Print( ( subfoldername + "Fit_Calibration.pdf" ).c_str() );
  // c3 -> Write();
  // c3 -> Close();

  // Using the sigma as the mean error
  TGraphErrors *gr1  =  new TGraphErrors( PEAKS_NUMBER, sortmean, sortpeaks, sortmeanSD, 0 ); // sortsigma
  gr1 -> Fit( "pol1", "FQ" );                             gr1 -> SetTitle( "" );
  gr1 -> SetMarkerColor( 4 );                             gr1 -> SetMarkerStyle( 20 );
  gr1 -> GetXaxis() -> SetTitleSize( 0.055 );             gr1 -> GetYaxis() -> SetTitleSize( 0.055 );
  gr1 -> GetYaxis() -> SetTitleOffset( 0.85 );            gr1 -> GetXaxis() -> SetTitleOffset( 0.85 ); // 90 e 95
  gr1 -> GetXaxis() -> SetLabelFont( 42 );                gr1 -> GetYaxis() -> SetLabelFont( 42 );
  gr1 -> GetXaxis() -> SetLabelSize( 0.035 );             gr1 -> GetYaxis() -> SetLabelSize( 0.035 );
  gr1 -> GetYaxis() -> SetTitle("E_{ teorico } ( keV )"); gr1 -> GetXaxis() -> SetTitle("E_{ medido } ( ToT )");
  gr1 -> GetXaxis() -> CenterTitle();                     gr1 -> GetYaxis() -> CenterTitle();
  gr1 -> Draw( "AP" );

  c3 -> Print( ( subfoldername + "Fit_Calibration.pdf" ).c_str() );
  c3 -> Write();
  c3 -> Close();

  TF1 *myfunc     =  gr1    -> GetFunction( "pol1" );
  Float_t b       =  myfunc -> GetParameter( 0 );  Float_t berror  =  myfunc -> GetParError( 0 );
  Float_t a       =  myfunc -> GetParameter( 1 );  Float_t aerror  =  myfunc -> GetParError( 1 );

  Float_t sigmakev[PEAKS_NUMBER];                  Float_t meankev[PEAKS_NUMBER];
  Float_t sigmakeverror[PEAKS_NUMBER];             Float_t meankeverror[PEAKS_NUMBER];

  //Output after calibration values
  FILE *fp  =  fopen( ( subfoldername + "Calibrated_sigma_energy.txt" ).c_str(), "w" );
  fprintf( fp, "%s %s %s %s %s %s %s %s %s %s %s \n", "Peaks", "MeanToT", "MeanToTerror", "MeanToTSD",
               "SigmaToT", "SigmaToTerror", "SigmaToTSD", "Mean", "MeanError", "Sigma", "SigmaError" );

  for( Int_t q = 0; q < PEAKS_NUMBER ; q += 1 )
  {
    sigmakev[q]        = a*sortsigma[q];
    meankev[q]         = a*sortmean[q] + b;
    sigmakeverror[q]   = a*sortsigmaSD[q]; // ( ( sortsigma[q]*aerror )**2 + ( sortsigmaerror[q]*a )**2)**0.5;
    meankeverror[q]    = a*sortmeanSD[q]; // ( ( sortmean[q]*aerror )**2  + ( sortmeanerror[q]*a )**2 +  berror**2)**0.5;
    //sigmakeverror[q]   = (2**0.5)*meankeverror[q];

    fprintf( fp, "%f %f %f %f %f %f %f %f %f %f %f \n", sortpeaks[q], sortmean[q], sortmeanerror[q], sortmeanSD[q],
       sortsigma[q], sortsigmaerror[q], sortsigmaSD[q], meankev[q], meankeverror[q], sigmakev[q], sigmakeverror[q] );
  }


  TCanvas *c3  =  new TCanvas( "c3", " ", 0, 0, 600, 400 );
  c3 -> GetFrame() -> SetFillColor( 21 );        c3 -> cd();
  c3 -> SetBorderMode( 0 );                      c3 -> SetBorderSize( 4 );
  c3 -> SetFrameBorderMode( 0 );

  gStyle -> SetOptStat( "n" );              gStyle -> SetOptTitle( 0 );
  gStyle -> SetStatY( 0.9 );                gStyle -> SetStatX( 0.5 );
  gStyle -> SetOptFit( 1 );

  TF1 *fd2     =  new TF1( "f2", "[0] -(-1*[2]*x + [1])/(TMath::Exp(-1*[4]*x))", 10, 150 ); //"log( x )*pol2" 5.83 0.06 0.025
  TF1 *fd3     =  new TF1( "f3", "[0] + ( ( [1]*(x**[2]) )/( ([3]**[2]) + (x**[2]) ) )", 10, 150 );
  TF1 *fd4     =  new TF1( "f4", "1.3 + ( 17*( x**2.7 )/( (113**2.7) + (x**2.7) ) )", 10, 150 ); // c + ( (a1)/(1 + exp(-a2*(meandata - c2 ) )) )
  TF1 *fd5     =  new TF1( "f5", "[0] + ( ( [1] )/(1 + TMath::Exp(-1*[2]*( x - [3] ) ) ) )", 10, 150 ); // c + ( (a1)/(1 + exp(-a2*(meandata - c2 ) )) )
  //5.83 - (-0.06*var + 3.96)/(exp(-0.025*var))

  TGraphErrors *gr2  =  new TGraphErrors( PEAKS_NUMBER, meankev, sigmakev, 0, sigmakeverror ); //meankeverror
  //gr2 -> Fit( "pol1", "MEFQ" );
  //gr2 -> Fit( "pol3", "FQ" ); // After change the error to be the sigma peak
  gr2 -> Fit( "f5", "ME" ); //MER
  gr2 -> SetTitle( "" );
  gr2 -> SetMarkerColor( 4 );                       gr2 -> SetMarkerStyle( 20 );
  gr2 -> GetXaxis() -> SetTitleSize( 0.055 );       gr2 -> GetYaxis() -> SetTitleSize( 0.055 );
  gr2 -> GetYaxis() -> SetTitleOffset( 0.90 );      gr2 -> GetXaxis() -> SetTitleOffset( 0.95 );
  gr2 -> GetXaxis() -> SetLabelFont( 42 );          gr2 -> GetYaxis() -> SetLabelFont( 42 );
  gr2 -> GetXaxis() -> SetLabelSize( 0.035 );       gr2 -> GetYaxis() -> SetLabelSize( 0.035 );
  gr2 -> GetXaxis() -> SetTitle( "E ( keV )" );     gr2 -> GetYaxis() -> SetTitle( "Sigma ( keV )" );
  gr2 -> GetXaxis() -> CenterTitle();               gr2 -> GetYaxis() -> CenterTitle();
  // gr2 -> GetXaxis() -> SetLimits( 0, 300);
  // gr2 -> GetYaxis() -> SetLimits( 0, 30);
  gr2 -> Draw( "AP" );
  // fd4 -> Draw("same");

  c3 -> Print( ( subfoldername + "Fit_Exsigma_kev.pdf" ).c_str() );
  c3 -> Write();
  c3 -> Close();

fclose( fp );  dfile -> Write();   dfile -> Close();

}
