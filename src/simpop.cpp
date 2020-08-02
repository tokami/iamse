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
List simpop(double logFM, List dat, List set, int out) {

  double FM = exp(logFM);

  // Import
  int ny = as<int>(set["refYears"]);
  int ns = as<int>(dat["nseasons"]);
  int amax = as<int>(dat["amax"]);
  amax = amax + 1;
  double R0 = as<double>(dat["R0"]);
  double h = as<double>(dat["h"]);
  double bp = as<double>(dat["bp"]);
  double recBeta = as<double>(dat["recBeta"]);
  double recGamma = as<double>(dat["recGamma"]);
  std::string SR = as<std::string>(dat["SR"]);
  double pzbm = as<double>(dat["pzbm"]);
  double maxF = as<double>(set["maxF"]);
  NumericMatrix weights = as<NumericMatrix>(dat["weights"]);
  NumericVector weight = as<NumericVector>(dat["weight"]);
  NumericMatrix weightFs = as<NumericMatrix>(dat["weightFs"]);
  NumericMatrix Ms = as<NumericMatrix>(dat["Ms"]);
  NumericMatrix mats = as<NumericMatrix>(dat["mats"]);
  NumericVector mat = as<NumericVector>(dat["mat"]);
  NumericMatrix sels = as<NumericMatrix>(dat["sels"]);
  NumericVector sel = as<NumericVector>(dat["sel"]);
  NumericVector M = as<NumericVector>(dat["M"]);
  NumericVector initN = as<NumericVector>(dat["initN"]);
  Function getFM("getFM");

  // Containers
  NumericVector Bage (amax);
  NumericVector SSBage (amax);
  NumericVector ESBage (amax);
  NumericVector CAA (amax);
  NumericVector NAA (amax);
  NumericMatrix FAA (amax, ns);
  NumericMatrix MAA (amax, ns);
  NumericMatrix ZAA (amax, ns);
  NumericVector CW (ny);
  NumericVector SP (ny);
  NumericVector SSB (ny);
  NumericVector SSB2 (ny);
  NumericVector TSB (ny);
  NumericVector TSB1plus (ny);
  NumericVector ESB (ny);
  NumericVector SSBPR0 (ny);
  NumericVector Ntemp (amax);
  NumericMatrix maty(amax, ns);
  NumericVector NnatM (amax);
  NumericVector Myear(amax);
  NumericVector matyear(amax);
  NumericVector Mcumsum(amax);
  double rec = 0.0;
  double fecun = 1.0;
  double alpha = 0.0;
  double beta = 0.0;
  double hy = 0.0;
  double R0y = 0.0;
  double Btmp = 0.0;
  double Ctmp = 0.0;
  double Ftmp = 0.0;

  // Initialise
  std::fill( CW.begin(), CW.end(), 0);
  std::fill( SP.begin(), SP.end(), 0);
  std::fill( SSB.begin(), SSB.end(), 0);
  std::fill( SSB2.begin(), SSB2.end(), 0);
  std::fill( TSB.begin(), TSB.end(), 0);
  std::fill( TSB1plus.begin(), TSB1plus.end(), 0);
  std::fill( ESB.begin(), ESB.end(), 0);
  std::fill( SSBPR0.begin(), SSBPR0.end(), 0);
  for(int a=0; a<amax; a++) ZAA(a,0) = M(a) + FM * sel(a);
  NAA(0) = exp(initN(0)) * R0;
  for(int a=1; a<amax; a++) NAA(a) = NAA(a-1) * exp(-ZAA(a-1,0)) * exp(initN(a));

  // errors
  double sdF = as<double>(set["sigmaF"]);
  double sdM = as<double>(set["sigmaM"]);
  double sdH = as<double>(set["sigmaH"]);
  double sdR0 = as<double>(set["sigmaR0"]);
  double sdMat = as<double>(set["sigmaMat"]);
  NumericVector eF = rnorm(ny, 0, sdF);
  double sdF2 = pow(sdF,2);
  for(int y=0; y<ny; y++) eF(y) = exp(eF(y) - sdF2/2);
  NumericVector eM = rnorm(ny, 0, sdM);
  double sdM2 = pow(sdM,2);
  for(int y=0; y<ny; y++) eM(y) = exp(eM(y) - sdM2/2);
  NumericVector eH = rnorm(ny, 0, sdH);
  double sdH2 = pow(sdH,2);
  for(int y=0; y<ny; y++) eH(y) = exp(eH(y) - sdH2/2);
  NumericVector eR0 = rnorm(ny, 0, sdR0);
  double sdR02 = pow(sdR0,2);
  for(int y=0; y<ny; y++) eR0(y) = exp(eR0(y) - sdR02/2);
  NumericVector eMat = rnorm(ny, 0, sdMat);
  double sdMat2 = pow(sdMat,2);
  for(int y=0; y<ny; y++) eMat(y) = exp(eMat(y) - sdMat2/2);

  // recruitment devs
  double sdR = as<double>(set["sigmaR"]);
  double rhoR = as<double>(set["rhoR"]);
  double sdR2 = pow(sdR,2);
  double rhoR2 = pow(rhoR,2);
  NumericVector rnum = rnorm(ny, 0, sdR) - sdR2/2;
  NumericVector eR(rnum.size());
  eR(0) = rnum(0);
  for(int i=1; i<rnum.size(); i++) eR(i) = rhoR * eR(i-1) + sqrt(1-rhoR2) * rnum(i);
  eR = exp(eR);
  double eRmean;
  for(int i=0; i<eR.size(); i++) eRmean += eR(i);
  eRmean = eRmean / eR.size();
  eR = eR / eRmean;

  // years
  for(int y=0; y<ny; y++){
    // Adding noise
    hy = h * eH(y);
    MAA = Ms * eM(y);
    maty = mats * eMat(y);
    R0y = R0 * eR0(y);
    Myear = M * eM(y);
    matyear = mat * eMat(y);

    // SSB per R0
    NnatM(0) = R0y;
    Mcumsum(0) = Myear(0);
    for(int a=1; a<amax; a++) Mcumsum(a) = Myear(a) + Mcumsum(a-1);
    for(int a=1; a<(amax-1); a++){
      NnatM(a) = R0y * exp(-Mcumsum(a-1));
    }
    NnatM(amax-1) = R0y * exp(-Mcumsum(amax-2)) / (1 - exp(-Mcumsum(amax-1)));
    for(int a=0; a<amax; a++) SSBPR0(y) += NnatM(a) * matyear(a) * weight(a) * fecun;
    // SSB
    for(int a=0; a<amax; a++){
      FAA(a,0) = sels(a,0) * FM * eF(y) / ns;
      ZAA(a,0) = FAA(a,0) + MAA(a,0);
      SSB(y) += NAA(a) * maty(a,0) * weights(a,0) * exp(-pzbm * ZAA(a,0));
    }
    // Recruitment
    if(SR == "bevholt"){
      rec = 4 * hy * R0y * SSB(y) / (SSBPR0(y) * (1-hy) + SSB(y) * (5*hy-1));
    }else if(SR == "ricker"){
      beta = log(5 * hy) / (0.8 * R0y);
      alpha = exp(beta * R0y)/SSBPR0(y);
      rec = alpha * SSB(y) * exp(-beta * SSB(y));
    }else if(SR == "average"){
      rec = R0y;
    }else if(SR == "hockey-stick"){
      if(SSB(y) > bp){
        rec = R0y;
      }else{
        rec = SSB(y) * R0y/bp;
      }
    }else if(SR == "bent-hyperbola"){
      rec = recBeta * (SSB(y) + sqrt(pow(bp,2) + pow(recGamma,2)/4) -
                         sqrt(pow(SSB(y)-bp,2) + pow(recGamma,2)/4));
    }
    NAA(0) = rec * eR(y);
    //    std::cout << "NAA(0)" << NAA(0) << std::endl;

    // Seasons
    for(int s=0; s<ns; s++){
      Btmp = 0.0;
      Ctmp = 0.0;

      // ages
      for(int a=0; a<amax; a++){
        // vulnerable midyear biomass
        Btmp += NAA(a) * weights(a,s) * sels(a,s) * exp(-MAA(a,s)/2);
        // catch
        FAA(a,s) = sels(a,s) * FM * eF(y) / ns;
        ZAA(a,s) = FAA(a,s) + MAA(a,s);
        CAA(a) = FAA(a,s)/ZAA(a,s) * NAA(a) * (1 - exp(-ZAA(a,s)));
        Ctmp += CAA(a) * weightFs(a,s);
      }

      // can't catch more than what's there
      // if(Ctmp > 0.99 * Btmp){
      //   for(int a=0; a<amax; a++){
      //     Ftmp = as<double>(getFM(0.75 * Btmp, NAA, MAA(_,s), weightFs(_,s), sels(_,s)));
      //     if(Ftmp > maxF/ns){
      //       Ftmp = maxF / ns;
      //     }
      //     FAA(a,s) = sels(a,s) * Ftmp;
      //     ZAA(a,s) = FAA(a,s) + MAA(a,s);
      //     CAA(a) = FAA(a,s)/ZAA(a,s) * NAA(a) * (1 - exp(-ZAA(a,s)));
      //     CW(y) += CAA(a) * weightFs(a,s);
      //   }
      // }else{
      //   CW(y) += Ctmp;
      // }
      CW(y) += Ctmp;

      for(int a=0; a<amax; a++){
        Ntemp(a) = NAA(a) * exp(-ZAA(a,s));
        NAA(a) = Ntemp(a);
      }
      if(s == (ns-1)){
        for(int a=0; a<amax; a++){
          // biomasses
          Bage(a) = NAA(a) * weights(a,s);
          SSBage(a) = Bage(a) * maty(a,s) * exp(-pzbm * ZAA(a,s));
          ESBage(a) = Bage(a) * sels(a,s);
          TSB(y) += Bage(a);
          if(a > 0) TSB1plus(y) += Bage(a);
          SSB2(y) += SSBage(a);
          ESB(y) += ESBage(a);
        }
      }
      if(s == (ns-1)){
        NAA(amax-1) = Ntemp(amax-1) + Ntemp(amax-2);
        for(int a=1; a<(amax-1); a++){
          NAA(a) = Ntemp(a-1);
        }
      }
    }
  }

  for(int y=1; y<ny; y++){
    SP(y) = TSB(y) - TSB(y-1) + CW(y);  // SP(y) = TSB(y+1) - TSB(y) + CW(y);  // if TSB in first season
  }

  List res;
  if(out == 0){
    res["CW"] = CW;
    res["TSB"] = TSB;
    res["SP"] = SP;
    res["TSB1plus"] = TSB1plus;
    res["ESB"] = ESB;
    res["SSB"] = SSB2;
  }else if(out == 1){
    res = SP(ny-2);
  }

  return res;
}
