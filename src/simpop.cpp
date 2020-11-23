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
  std::string refmethod = as<std::string>(set["refMethod"]);
  int nyrefmsy = as<int>(set["refYearsMSY"]);
  NumericMatrix weights = as<NumericMatrix>(dat["weights"]);
  NumericVector weight = as<NumericVector>(dat["weight"]);
  NumericMatrix weightFs = as<NumericMatrix>(dat["weightFs"]);
  NumericVector Ms = as<NumericVector>(dat["Ms"]);
  NumericMatrix Msels = as<NumericMatrix>(dat["Msels"]);
  NumericMatrix mats = as<NumericMatrix>(dat["mats"]);
  NumericVector mat = as<NumericVector>(dat["mat"]);
  NumericMatrix sels = as<NumericMatrix>(dat["sels"]);
  NumericVector sel = as<NumericVector>(dat["sel"]);
  NumericVector initN = as<NumericVector>(dat["initN"]);
  int sptype = as<int>(set["spType"]);
  int recTiming = as<int>(set["recruitmentTiming"]);

  // Containers
  NumericVector Bage (amax);
  NumericVector SSBage (amax);
  NumericVector ESBage (amax);
  NumericVector CAA (amax);
  NumericVector NAA (amax);
  NumericVector Mtot (amax);
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
  NumericMatrix NnatM (amax, ns);
  NumericVector SPR (ns);
  NumericVector Myear(amax);
  NumericVector matyear(amax);
  NumericVector Mcumsum(amax);
  double rec = 0.0;
  double fecun = 1.0;
  double alpha = 0.0;
  double beta = 0.0;
  double hy = 0.0;
  double R0y = 0.0;
  double Ctmp = 0.0;

  // Initialise
  std::fill( Mtot.begin(), Mtot.end(), 0);
  for(int a=0; a<amax; a++) for(int s=0; s<ns; s++) Mtot(a) = Mtot(a) + Ms(0) * Msels(a,s);
  std::fill( SPR.begin(), SPR.end(), 0);
  std::fill( CW.begin(), CW.end(), 0);
  std::fill( SP.begin(), SP.end(), 0);
  std::fill( SSB.begin(), SSB.end(), 0);
  std::fill( SSB2.begin(), SSB2.end(), 0);
  std::fill( TSB.begin(), TSB.end(), 0);
  std::fill( TSB1plus.begin(), TSB1plus.end(), 0);
  std::fill( ESB.begin(), ESB.end(), 0);
  std::fill( SSBPR0.begin(), SSBPR0.end(), 0);
  for(int a=0; a<amax; a++) ZAA(a,0) = Mtot(a) + FM * sel(a);  // no noise on M
  NAA(0) = exp(initN(0)) * R0;
  for(int a=1; a<amax; a++) NAA(a) = NAA(a-1) * exp(-ZAA(a-1,0)) * exp(initN(a));

  // errors
  //  NumericVector eF = as<NumericVector>(set["eF"]);
  NumericVector eR = as<NumericVector>(set["eR"]);
  NumericVector eM = as<NumericVector>(set["eM"]);
  NumericVector eH = as<NumericVector>(set["eH"]);
  NumericVector eR0 = as<NumericVector>(set["eR0"]);
  NumericVector eMat = as<NumericVector>(set["eMat"]);
  //  NumericVector eImp = as<NumericVector>(set["eImp"]);

  // years
  for(int y=0; y<ny; y++){
    // Adding noise
    hy = h * eH(y);
    MAA = Ms(y) * Msels * eM(y);
    std::fill( Mtot.begin(), Mtot.end(), 0);
    for(int a=0; a<amax; a++) for(int s=0; s<ns; s++) Mtot(a) = Mtot(a) + Ms(y) * Msels(a,s);
    maty = mats * eMat(y);
    R0y = R0 * eR0(y);
    matyear = mat * eMat(y);

    // SPR for SR
    // ----------------------------------
    // set SPR to zero each year
    std::fill( SPR.begin(), SPR.end(), 0);
    // cumulative M
    Mcumsum(0) = Mtot(0);
    for(int a=1; a<amax; a++) Mcumsum(a) = Mtot(a) * eM(y) + Mcumsum(a-1);
    // pop
    NnatM(0,0) = R0y;
    for(int a=1; a<(amax-1); a++){
      NnatM(a,0) = R0y * exp(-Mcumsum(a-1));
    }
    NnatM(amax-1,0) = R0y * exp(-Mcumsum(amax-2)) / (1 - exp(-Mcumsum(amax-1)));
    // if multiple seasons
    if(ns > 1){
      for(int s=1; s<ns; s++){
        for(int a=0; a<amax; a++){
          NnatM(a,s) = NnatM(a,s-1) * exp(-MAA(a,s));
        }
      }
    }
    // SPR
    for(int s=0; s<ns; s++){
      for(int a=0; a<amax; a++){
        SPR(s) += NnatM(a,s) * maty(a,s) * weights(a,s) * fecun;
      }
    }
    SSBPR0(y) = SPR(recTiming);

    // SSB
    for(int a=0; a<amax; a++){
      FAA(a,0) = sels(a,0) * FM / ns;  // Casper 13/08: constant F for ref estimation, no noise on F // * eF(y)
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

    // seasons
    for(int s=0; s<ns; s++){
      Ctmp = 0.0;

      // ages
      for(int a=0; a<amax; a++){
        FAA(a,s) = sels(a,s) * FM / ns; // * eF(y)
        ZAA(a,s) = FAA(a,s) + MAA(a,s);
        CAA(a) = FAA(a,s)/ZAA(a,s) * NAA(a) * (1 - exp(-ZAA(a,s)));
        Ctmp += CAA(a) * weightFs(a,s);
      }
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

  if(sptype == 0){
    for(int y=1; y<ny; y++){
      SP(y) = TSB(y) - TSB(y-1) + CW(y);  // SP(y) = TSB(y+1) - TSB(y) + CW(y);  // if TSB in first season
    }
  }else if(sptype == 1){
    for(int y=1; y<ny; y++){
      SP(y) = ESB(y) - ESB(y-1) + CW(y);
    }
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
    // median SP over last 50 years (or mean)
    NumericVector sp2 = tail(SP, nyrefmsy);
    double tmp;
    if(refmethod == "mean"){
      tmp = mean(sp2);
    }else if(refmethod == "median"){
      tmp = median(sp2);
    }
    res = tmp;
    // res = SP(ny-2);
  }

  return res;
}
