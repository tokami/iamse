//' @name simpop
//'
//' @title Simulate an age-based population
//'
//' @details Simulate a pop
//'
//' @param logF fishing mortality
//' @param dat List with species data
//' @param set List with MSE settings
//' @param opt If 1 the function returns the yield in the last year,
//' if 2 the function returns a list with yield, TSB, SSB, and ESB over
//' the whole simulation period.
//'
//' @export

#include <iostream>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List simpop(double FM, List dat, List set) {

  // import
  int ny = as<int>(set["refYears"]);
  int amax = as<int>(dat["amax"]);
  amax = amax + 1;
  NumericVector weight = as<NumericVector>(dat["weight"]);
  NumericVector M = as<NumericVector>(dat["M"]);
  NumericVector mat = as<NumericVector>(dat["mat"]);
  NumericVector sel = as<NumericVector>(dat["sel"]);
  double logR0 = as<double>(dat["logR0"]);
  double R0 = exp(logR0);
  double h = as<double>(dat["h"]);
  std::string SR = as<std::string>(dat["SR"]);
  NumericMatrix Nage (amax, ny+1);
  NumericMatrix Bage (amax, ny+1);
  NumericMatrix SSBage (amax, ny+1);
  NumericMatrix ESBage (amax, ny+1);
  NumericMatrix CAA (amax, ny+1);
  NumericMatrix FAA (amax, ny+1);
  NumericVector CW (ny);
  std::fill( CW.begin(), CW.end(), 0);
  NumericVector SP (ny);
  std::fill( SP.begin(), SP.end(), 0);
  NumericVector SSB (ny);
  std::fill( SSB.begin(), SSB.end(), 0);
  NumericVector SSB2 (ny);
  std::fill( SSB2.begin(), SSB2.end(), 0);
  NumericVector TSB (ny);
  std::fill( TSB.begin(), TSB.end(), 0);
  NumericVector TSB1plus (ny);
  std::fill( TSB1plus.begin(), TSB1plus.end(), 0);
  NumericVector ESB (ny);
  std::fill( ESB.begin(), ESB.end(), 0);
  NumericVector SSBPR0 (ny);
  std::fill( SSBPR0.begin(), SSBPR0.end(), 0);
  NumericVector NAtmp (amax);
  NumericVector NAtmp2 (amax);
  for(int a=0; a<amax; a++){
    NAtmp(a) = R0;
    NAtmp2(a) = R0;
    Nage(a,0) = R0;
  }
  NumericVector Z (amax);
  NumericVector NAAmid (amax);
  double rec;
  double fecun = 1.0;
  NumericVector NnatM (amax);
  NumericVector survivors (amax);
  double alpha;
  double beta;

  // errors
  double sdF = as<double>(set["sigmaF"]);
  double sdR = as<double>(set["sigmaR"]);
  double sdM = as<double>(set["sigmaM"]);
  double sdH = as<double>(set["sigmaH"]);
  double sdR0 = as<double>(set["sigmaR0"]);
  double sdMat = as<double>(set["sigmaMat"]);
  NumericVector eF = rnorm(ny, 0, sdF);
  double sdF2 = pow(sdF,2);
  for(int y=0; y<ny; y++) eF(y) = eF(y) - sdF2/2;
  NumericVector eR = rnorm(ny, 0, sdR);
  double sdR2 = pow(sdR,2);
  for(int y=0; y<ny; y++) eR(y) = eR(y) - sdR2/2;
  NumericVector eM = rnorm(ny, 0, sdM);
  double sdM2 = pow(sdM,2);
  for(int y=0; y<ny; y++) eM(y) = eM(y) - sdM2/2;
  NumericVector eH = rnorm(ny, 0, sdH);
  double sdH2 = pow(sdH,2);
  for(int y=0; y<ny; y++) eH(y) = eH(y) - sdH2/2;
  NumericVector eR0 = rnorm(ny, 0, sdR0);
  double sdR02 = pow(sdR0,2);
  for(int y=0; y<ny; y++) eR0(y) = eR0(y) - sdR02/2;
  NumericVector eMat = rnorm(ny, 0, sdMat);
  double sdMat2 = pow(sdMat,2);
  for(int y=0; y<ny; y++) eMat(y) = eMat(y) - sdMat2/2;

  double hy;
  NumericVector My(M.size());
  NumericVector maty(mat.size());
  double R0y;

  // loop through years
  for(int y=0; y<ny; y++){
    hy = h * exp(eH(y));
    My = M * exp(eM(y));
    maty = mat * exp(eMat(y));
    R0y = exp(logR0) * exp(eR0(y));
    for(int a=0; a<amax; a++){
      Nage(a,y) = NAtmp(a);
      NAAmid(a) = Nage(a,y) * exp(-My(a)/2);
      FAA(a,y) = sel(a) * FM * exp(eF(y));
      Z(a) = FAA(a,y) + My(a);
      CAA(a,y) = FAA(a,y)/Z(a) * NAAmid(a) * (1 - exp(-Z(a)));
      CW(y) += CAA(a,y) * weight(a); // weightF?
      SSB(y) += NAtmp(a) * maty(a) * weight(a);
    }
    // SSB per R0
    NnatM(0) = 1;
    for(int a=1; a<amax; a++){
      NnatM(a) = NnatM(a-1) * exp(-My(a-1));
    }
    NnatM(amax-1) = NnatM(amax-1) / (1 - exp(-My(amax-2)));
    for(int a=0; a<amax; a++) SSBPR0(y) += NnatM(a) * maty(a) * fecun;
    // recruitment
    if(SR == "bevholt"){
      alpha = SSBPR0(y) * ((1-hy)/(4*hy));
      beta = (5*hy-1)/(4*hy*R0y);
      rec = SSB(y) / (alpha + beta * SSB(y));
    }else if(SR == "ricker"){
      beta = log(5 * hy) / (0.8 * R0y);
      alpha = exp(beta * R0y)/SSBPR0(y);
      rec = alpha * SSB(y) * exp(-beta * SSB(y));
    }
    NAtmp2(0) = rec * exp(eR(y));
    // survivors
    for(int a=0; a<amax; a++){
      survivors(a) = NAtmp(a) * exp(-Z(a));
    }
    for(int a=1; a<amax; a++){
      NAtmp2(a) = survivors(a-1);
    }
    // plus group
    NAtmp2(amax-1) = NAtmp2(amax-1) + survivors(amax-1);
    // next years numbers
    NAtmp = NAtmp2;
  }

  // return
  for(int y=0; y<ny; y++){
    maty = mat * exp(eMat(y));
    for(int a=0; a<amax; a++){
      Bage(a,y) = Nage(a,y) * weight(a);
      SSBage(a,y) = Nage(a,y) * weight(a) * maty(a);
      ESBage(a,y) = Nage(a,y) * weight(a) * sel(a);
      TSB(y) += Bage(a,y);
      if(a > 0) TSB1plus(y) += Bage(a,y);
      SSB2(y) += SSBage(a,y);
      ESB(y) += ESBage(a,y);
    }
  }

  for(int y=0; y<(ny-1); y++){
    SP(y) = TSB(y+1) - TSB(y) + CW(y);
  }

  List out;
  out["CW"] = CW;
  out["TSB"] = TSB;
  out["SP"] = SP;
  out["TSB1plus"] = TSB1plus;
  out["ESB"] = ESB;
  out["SSB"] = SSB2;
  return out;
}
