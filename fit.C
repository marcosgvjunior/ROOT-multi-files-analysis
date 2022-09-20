//Institute of Physics - Federal University of Rio de Janeiro
//Interdisciplinary Academic Master's Degree in Applied Physics
//Student: Marcos Vieira
//June, 2019
//Working on ROOT 5.34/36, Windows 10 x64

//Script to fit a specific peak energy for several files with different configurations

#include <TFile.h>
#include <TObject.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <string>
#include <vector>
#include <Riostream.h>
#include <sstream>
#include <iostream>
#include <fstream>

// Some references:
// https://root.cern.ch/root/html/tutorials/fit/multifit.C.html
// https://root.cern.ch/doc/master/FittingDemo_8C_source.html

#define FILES_NUMBER  10
#define CHAR_SIZE     100
#define MEAN_MODEL    59.54
#define REBIN_VALUE   1

//Initial guess values to set the limits of a energy peak and get its maximum (bin of maximum value)
#define ANG_COEF_GUESS        1
#define LIN_COEF_GUESS        0

//This flag includes/exlcudes addional analysis plots
#define ADD_PLOTS     false

void fit()
{
  // Force batch mode
  gROOT -> SetBatch( kTRUE );

  //input files
  std::string files[] = { "am241_thl0keV_noise_200_SinglePixel",
                          "am241_thl1keV_noise_200_SinglePixel",
                          "am241_thl2keV_noise_200_SinglePixel",
                          "am241_thl3keV_noise_200_SinglePixel",
                          "am241_thl4keV_noise_200_SinglePixel",
                          "am241_thl5keV_noise_200_SinglePixel",
                          "am241_thl6keV_noise_200_SinglePixel",
                          "am241_thl7keV_noise_200_SinglePixel",
                          "am241_thl8keV_noise_200_SinglePixel",
                          "am241_thl9keV_noise_200_SinglePixel" };

  // Histrograms and legend identification
  //std::string legend[] = {"Noise: 0", "Noise: 80", "Noise: 120", "Noise: 160", "Noise: 200" };
  std::string legend[] = {"TH: 0 keV",  "TH: 1 keV",  "TH: 2 keV",  "TH: 3 keV",  "TH: 4 keV",
                          "TH: 5 keV",  "TH: 6 keV",  "TH: 7 keV",  "TH: 8 keV",  "TH: 9 keV" };

  // Name of the input histogram
  std::string histogramname = "H11";

  // Inferior and superior limits for the fit
  //Int_t limi[] = {59, 58, 57, 56, 55 };
  //Int_t limf[] = {60, 62, 63, 64, 65 };
  Int_t limi[] = {54, 54, 54, 54, 54, 54, 54, 54, 54, 54 };
  Int_t limf[] = {66, 66, 66, 66, 66, 66, 66, 66, 66, 66 };

  // Increments to the limits for visualization
  Int_t somalimi[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  Int_t somalimf[] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

  // List of colors
  Int_t colors[] = {1, 2, 3, 6, 7, 8, 9, 11, 41, 42, 49, 28, 32, 40 }; //5, 4

  // Initializing a list for each histogram and its names/legends
  TH1F *mhist[FILES_NUMBER];
  char *histname = new char[CHAR_SIZE];
  char *legd     = new char[CHAR_SIZE];

  // Variables to store filename and its path
  std::string filename;                   std::string path;

  // Input folder path and path/folder to save the outputs
  //**Need to create the folder before run
  std::string relativepath   = "";        std::string subfoldername  = "fit/";

  // Output file with fit details
  FILE *fp = fopen( ( subfoldername + "fit_" + filename + ".txt" ).c_str(), "w" );

  // Variables to store parameters values
  Double_t  par[3];
  Float_t   thl[]   = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  Float_t   noise[] = {0, 80, 120, 160, 200 };
  Float_t   sigma[FILES_NUMBER];              Float_t   sigmaerror[FILES_NUMBER];
  Float_t   mean[FILES_NUMBER];               Float_t   meanerror[FILES_NUMBER];
  Float_t   sigmacalc[FILES_NUMBER];          Float_t   sigmacalcerror[FILES_NUMBER];

  // Beginning of the loop on files
  for( Int_t m = 0; m < FILES_NUMBER ; m  += 1)
  {
    filename = files[m];
    path     = relativepath + filename + ".root";

    TFile *f1 = TFile::Open( ( path ).c_str() );

    sprintf( histname, "%s", files[m]  ); // one can change to legend[m] or other list of names
    sprintf( legd    , "%s", legend[m] );

    mhist[m] = ( TH1F* )f1 -> Get( ( histogramname ).c_str() ) -> Clone();
    mhist[m] -> SetName( histname ); //
    mhist[m] -> SetTitle( legd );

    TFile *dfile = new TFile( ( subfoldername + filename + ".root" ).c_str(),  "RECREATE" );

    // There is probably a way to put the canvas ouside the loop and optimize
    TCanvas *c1 = new TCanvas( "c1", " ", 300, 200, 300, 200 ); //600, 400 ... maybe a macro value here

    c1 -> GetFrame() -> SetFillColor( 21 );         c1 -> cd();
    c1 -> SetBorderMode( 0 );                       c1 -> SetBorderSize( 0 ); //4
    c1 -> SetFrameBorderMode( 0 );

    gStyle -> SetOptStat( "n" );                  gStyle -> SetOptTitle( 0 );
    gStyle -> SetOptFit( 1111 );                  gStyle -> SetStatFont( 42 );
    gStyle -> SetStatY( 0.9 );                    gStyle -> SetStatX( 0.9 );
    gStyle -> SetStatH( 0.2 );                    gStyle -> SetStatW( 0.17 );
    //gStyle->SetStatFontSize(0.04);            //gStyle->SetLabelSize(0.04);

    // There are different possible places to rebin and set the ranges
    mhist[m] -> Rebin( REBIN_VALUE );
    mhist[m] -> GetXaxis() -> SetRangeUser( limi[m]*ANG_COEF_GUESS - LIN_COEF_GUESS,
                                            limf[m]*ANG_COEF_GUESS - LIN_COEF_GUESS );

    Int_t maxbin = mhist[m] -> GetMaximumBin();

    //mhist[m] -> GetYaxis() -> SetRangeUser( 0, 300 ); //2*( mhist[m] -> GetMaximum() )
    //mhist[m] -> GetXaxis() -> SetRange( REBIN_VALUE*maxbin-2*somalimi[m],
                                        //REBIN_VALUE*maxbin + 3*somalimf[m] );  // 0, 1000 rangeuser?

    mhist[m] -> SetLineColor( colors[m] );            mhist[m] -> SetMarkerColor( colors[m] );
    mhist[m] -> GetXaxis() -> SetDecimals();
    mhist[m] -> GetXaxis() -> SetTitleSize( 0.055 );  mhist[m] -> GetYaxis() -> SetTitleSize( 0.055 );
    mhist[m] -> GetXaxis() -> CenterTitle();          mhist[m] -> GetYaxis() -> CenterTitle();
    mhist[m] -> GetXaxis() -> SetTitle( "keV" );      mhist[m] -> GetYaxis() -> SetTitle( "Contagens" );
    mhist[m] -> GetYaxis() -> SetTitleOffset( 0.90 ); mhist[m] -> GetXaxis() -> SetTitleOffset( 0.95 );
    mhist[m] -> GetXaxis() -> SetLabelFont( 42 );     mhist[m] -> GetYaxis() -> SetLabelFont( 42 );
    mhist[m] -> GetXaxis() -> SetLabelSize( 0.045 );  mhist[m] -> GetYaxis() -> SetLabelSize( 0.045 );

    //TF1 *g1    = new TF1( "g1", "gaus", REBIN_VALUE*maxbin-somalimi[m], REBIN_VALUE*maxbin + somalimf[m] );
    TF1 *g1    = new TF1( "g1", "gaus", limi[m], limf[m] );

    g1 -> SetLineColor( colors[m + 1] );             g1 -> SetMarkerColor( colors[m + 1] );
    g1 -> SetLineWidth( 4 ); //1

    mhist[m] -> Fit( g1, "QR" );
    mhist[m] -> GetXaxis() -> SetRangeUser( limi[m] - somalimi[m], limf[m] + somalimf[m] );

    g1 -> GetParameters( &par[0] );

    mean[m]       = g1 -> GetParameter( 1 );    meanerror[m]  = g1 -> GetParError( 1 );
    sigma[m]      = g1 -> GetParameter( 2 );    sigmaerror[m] = g1 -> GetParError( 2 );

    if( ADD_PLOTS )
    {
      sigmacalc[m]      = MEAN_MODEL*sigma[m]/mean[m];
      sigmacalcerror[m] = ( ( ( MEAN_MODEL*sigma[m]*meanerror[m]/( mean[m]**2 ) )**2 )
                            + ( MEAN_MODEL*sigmaerror[m]/mean[m] )**2 )**0.5;
    }

    if ( fp! = NULL)
    {
      for ( Int_t i = 0; i < g1 -> GetNpar(); i ++ )
      {
        Float_t value = g1   -> GetParameter( i );     Float_t error = g1   -> GetParError( i );
        Float_t chi   = g1   -> GetChisquare();        Float_t max   = g1   -> GetMaximum();
        fprintf( fp, "%f %f %f %f \n", value, error, chi, max );
      }
    }
    fprintf( fp, "\n" );

    mhist[m] -> SetStats( 1 );
    mhist[m] -> Write();
    mhist[m] -> Draw( "histsame" );

    c1 -> Print( ( subfoldername + filename + "_fit.pdf" ).c_str() );
    c1 -> Write();

    c1 -> SetLogy();
    c1 -> Print( ( subfoldername + filename + "_log_fit.pdf").c_str() );
    c1 -> Close();
  }
  fclose( fp );

  TCanvas *c2 = new TCanvas( "c2", " ", 0, 0, 600, 400 );
  c2 -> GetFrame() -> SetFillColor( 21 );           c2 -> cd();
  c2 -> SetBorderMode( 0 );                         c2 -> SetBorderSize( 4 );
  c2 -> SetFrameBorderMode( 0 );

  gStyle -> SetOptStat( "n" );                  gStyle -> SetOptTitle( 0 );
  gStyle -> SetStatY( 0.9 );                    gStyle -> SetStatX( 0.9 );  //0.5

  TGraphErrors *gr = new TGraphErrors( FILES_NUMBER, thl, sigma, 0, sigmaerror );
  //TF1 *fit1 = new TF1( "fit1", "pol1",  270,  390 );
  //gr -> Fit( fit1,  "FQR" );
  gr -> Fit( "pol1", "FQ" );
  gr -> SetTitle( "THL = 0" );                      //gr -> SetName( "Sigma x THL" );
  gr -> SetFillColor( kYellow );                    gr -> SetMarkerColor( 4 );
  gr -> SetMarkerStyle( 20 );
  gr -> GetXaxis() -> SetTitleSize( 0.055 );        gr -> GetYaxis() -> SetTitleSize( 0.055 );
  gr -> GetYaxis() -> SetTitleOffset( 0.90 );       gr -> GetXaxis() -> SetTitleOffset( 0.95 );
  gr -> GetXaxis() -> SetLabelFont( 42 );           gr -> GetYaxis() -> SetLabelFont( 42 );
  gr -> GetXaxis() -> SetLabelSize( 0.035 );        gr -> GetYaxis() -> SetLabelSize( 0.035 );
  gr -> GetYaxis() -> SetTitle( "Sigma ( keV )" );  gr -> GetXaxis() -> SetTitle( "Threshold ( keV )" );
  gr -> GetXaxis() -> CenterTitle();                gr -> GetYaxis() -> CenterTitle();
  gr -> Draw( "AP" );                             //fit1 -> Draw( "E3sameAP" );

  c2 -> Print( ( subfoldername + "sigma_thl.pdf" ).c_str() );
  c2 -> Write();
  c2 -> Close();

  if( ADD_PLOTS )
  {
    TCanvas *c3 = new TCanvas( "c3", " ", 0, 0, 600, 400 );
    c3 -> GetFrame() -> SetFillColor( 21 );           c3 -> cd();
    c3 -> SetBorderMode( 0 );                         c3 -> SetBorderSize( 4 );
    c3 -> SetFrameBorderMode( 0 );

    gStyle -> SetOptStat( "n" );                      gStyle -> SetOptTitle( 0 );
    gStyle -> SetStatY( 0.9 );                        gStyle -> SetStatX( 0.5 );

    TGraphErrors *gr1 = new TGraphErrors( FILES_NUMBER, thl, mean, 0, meanerror );
    gr1 -> Fit( "pol2", "FQ" );                       gr1 -> SetTitle( "" );
    gr1 -> SetMarkerColor( 4 );                       gr1 -> SetMarkerStyle( 20 );
    gr1 -> GetXaxis() -> SetTitleSize( 0.055 );       gr1 -> GetYaxis() -> SetTitleSize( 0.055 );
    gr1 -> GetYaxis() -> SetTitleOffset( 0.90 );      gr1 -> GetXaxis() -> SetTitleOffset( 0.95 );
    gr1 -> GetXaxis() -> SetLabelFont( 42 );          gr1 -> GetYaxis() -> SetLabelFont( 42 );
    gr1 -> GetXaxis() -> SetLabelSize( 0.035 );       gr1 -> GetYaxis() -> SetLabelSize( 0.035 );
    gr1 -> GetYaxis() -> SetTitle( "Mean ( ToT )" );  gr1 -> GetXaxis() -> SetTitle( "THL" );
    gr1 -> GetXaxis() -> CenterTitle();               gr1 -> GetYaxis() -> CenterTitle();
    gr1 -> Draw( "AP" );

    c3 -> Print( ( subfoldername + "mean_THL.pdf" ).c_str() );
    c3 -> Write();
    c3 -> Close();

    TCanvas *c4 = new TCanvas( "c4", " ", 0, 0, 600, 400 );
    c4 -> GetFrame() -> SetFillColor( 21 );           c4 -> cd();
    c4 -> SetBorderMode( 0 );                         c4 -> SetBorderSize( 4 );
    c4 -> SetFrameBorderMode( 0 );

    gStyle -> SetOptStat( "n" );                    gStyle -> SetOptTitle( 0 );
    gStyle -> SetStatY( 0.9 );                      gStyle -> SetStatX( 0.5 );

    TGraphErrors *gr2 = new TGraphErrors( FILES_NUMBER, mean, sigma, meanerror, sigmaerror );
    gr2 -> Fit( "pol1", "FQ" );                           gr2 -> SetTitle( "" );
    gr2 -> SetMarkerColor( 4 );                           gr2 -> SetMarkerStyle( 20 );
    gr2 -> GetXaxis() -> SetTitleSize( 0.055 );           gr2 -> GetYaxis() -> SetTitleSize( 0.055 );
    gr2 -> GetYaxis() -> SetTitleOffset( 0.90 );          gr2 -> GetXaxis() -> SetTitleOffset( 0.95 );
    gr2 -> GetXaxis() -> SetLabelFont( 42 );              gr2 -> GetYaxis() -> SetLabelFont( 42 );
    gr2 -> GetXaxis() -> SetLabelSize( 0.035 );           gr2 -> GetYaxis() -> SetLabelSize( 0.035 );
    gr2 -> GetXaxis() -> SetTitle( "Mean ( ToT )" );      gr2 -> GetYaxis() -> SetTitle( "Sigma ( ToT )" );
    gr2 -> GetXaxis() -> CenterTitle();                   gr2 -> GetYaxis() -> CenterTitle();
    gr2 -> Draw( "APY" );

    cout << gr2 -> GetCorrelationFactor() << "\n";
    cout << gr2 -> GetCorrelationMatrix() << "\n"; //erro
    cout << gr2 -> GetCovariance()        << "\n";
    cout << gr2 -> GetCovarianceMatrix()  << "\n"; //erro

    c4 -> Print( ( subfoldername + "mean_sigma.pdf" ).c_str() );
    c4 -> Write();
    c4 -> Close();

    TCanvas *c5 = new TCanvas( "c5", " ", 0, 0, 600, 400 );
    c5 -> GetFrame() -> SetFillColor( 21 );             c5 -> cd();
    c5 -> SetBorderMode( 0 );                           c5 -> SetBorderSize( 4 );
    c5 -> SetFrameBorderMode( 0 );

    gStyle -> SetOptStat( "n" );                    gStyle -> SetOptTitle( 0 );
    gStyle -> SetStatY( 0.9 );                      gStyle -> SetStatX( 0.9 ); // 0.65 meio

    TGraphErrors *gr3 = new TGraphErrors( FILES_NUMBER, thl, sigmacalc, 0, sigmacalcerror );
    //TF1 *fit1 = new TF1( "fit1", "pol1",  270,  390 );
    //gr -> Fit( fit1,  "FQR" );
    gr3 -> Fit( "pol1", "FQ" );
    //gr -> SetName( "Sigma x THL" );
    gr3 -> SetTitle( "" );
    gr3 -> SetFillColor( kYellow );                     gr3 -> SetMarkerColor( 4 );
    gr3 -> SetMarkerStyle( 20 );
    gr3 -> GetXaxis() -> SetTitleSize( 0.055 );         gr3 -> GetYaxis() -> SetTitleSize( 0.055 );
    gr3 -> GetYaxis() -> SetTitleOffset( 0.90 );        gr3 -> GetXaxis() -> SetTitleOffset( 0.95 );
    gr3 -> GetXaxis() -> SetLabelFont( 42 );            gr3 -> GetYaxis() -> SetLabelFont( 42 );
    gr3 -> GetXaxis() -> SetLabelSize( 0.035 );         gr3 -> GetYaxis() -> SetLabelSize( 0.035 );
    gr3 -> GetYaxis() -> SetTitle( "Sigma ( keV )" );   gr3 -> GetXaxis() -> SetTitle( "THL" );
    gr3 -> GetXaxis() -> CenterTitle();                 gr3 -> GetYaxis() -> CenterTitle();
    gr3 -> Draw( "AP" );                                fit1->  Draw( "E3sameAP" );
    c5 -> Print( ( subfoldername + "sigmacalc-keV_THL.pdf" ).c_str() );
    c5 -> Write();
    c5 -> Close();
  }

  dfile -> Write();  dfile -> Close();
}
