
#include <iostream>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector initdist(NumericVector MAA, NumericVector FAA, double R0, NumericVector spawning, IntegerVector inds){

  int asmax = MAA.size();
  int ns = spawning.size();
  NumericMatrix NAA(asmax, ns);
  NumericVector NAAS(asmax);
  NumericVector ZAA(asmax);
  int ind1 = 0;
  int ind2 = asmax - ns + 1;
  // std::cout << asmax << "ns: " << ns << std::endl;
  NumericVector pg(ns);

  for(int a=0; a<asmax; a++){
    ZAA(a) = MAA(a) + FAA(a);
  }
  // Initial distribution
  NAA(0,_) = R0 * spawning;
  for(int a=1; a<asmax; a++){
    NAA(a,_) = NAA(a-1,_) * exp(-ZAA(a-1));
  }
  for(int s=0; s<ns; s++){
    for(int i=0; i<inds.size(); i++){
      ind1 = inds(i) + s;
//      std::cout << ind1 << std::endl;
      NAAS(ind1) += NAA(ind1,s);
    }
  }
  NAAS(0) = 0;
  // plusgroup
  if(ns == 1){
    NAAS(asmax-1) = NAAS[asmax-1] / (1 - exp(-ZAA[asmax-1]));
  }else{
    pg = NAAS[Rcpp::Range(ind2, asmax-1)] / (1 - exp(-sum(ZAA[Rcpp::Range(ind2, asmax-1)])));
    NAAS(asmax-1) = sum(pg) - sum(NAAS[Rcpp::Range(ind2, asmax-1)]);
  }

  return NAAS;

}


//' @name simpop
//'
//' @title Simulate an age-based population
//'
//' @details Simulate a pop
//'
//' @param logF fishing mortality
//' @param dat List with species data
//' @param set List with MSE settings
//' @param tvy year index for all time-variant (tv) processes (so far: Msel, Ms, sel)
//' @param opt If 1 the function returns the yield in the last year,
//' if 2 the function returns a list with yield, TSB, SSB, and ESB over
//' the whole simulation period.
//'
//' @export


// [[Rcpp::export]]
List simpop(double logFM, List dat, List set, int out) {

  double FM = exp(logFM);

  // Import
  int ny = as<int>(set["refYears"]);
  int ns = as<int>(dat["ns"]);
  int amax = as<int>(dat["amax"]);
  amax = amax + 1;
  int asmax = amax * ns;
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
  NumericMatrix weight = as<NumericMatrix>(dat["weight"]);
  NumericMatrix weightF = as<NumericMatrix>(dat["weightF"]);
  NumericVector M = as<NumericVector>(dat["M"]);
  List MselList = as<List>(dat["Msel"]);
  List selList = as<List>(dat["sel"]);
  IntegerVector s1vec = as<IntegerVector>(dat["s1vec"]);
  NumericMatrix mat = as<NumericMatrix>(dat["mat"]);
  NumericVector initN = as<NumericVector>(dat["initN"]);
  int sptype = as<int>(set["spType"]);
  NumericVector spawning = as<NumericVector>(dat["spawning"]);
  // temporary entries
  IntegerVector inds = as<IntegerVector>(dat["inds"]) - 1;
  int tvm = as<int>(set["tvm"]);
  int tvmsel = as<int>(set["tvmsel"]);
  int tvsel = as<int>(set["tvsel"]);
  // still to figure out
  IntegerVector as2a = as<IntegerVector>(dat["as2a"]) - 1;
  IntegerVector as2s = as<IntegerVector>(dat["as2s"]) - 1;

  // Containers
  NumericVector Bage (asmax);
  NumericVector SSBage (asmax);
  NumericVector ESBage (asmax);
  NumericVector CAA (asmax);
  NumericVector NAAStmp (asmax);
  NumericVector FAA (asmax);
  NumericVector MAA0 (asmax);
  NumericVector MAA (asmax);
  NumericVector ZAA (asmax);
  NumericVector CW (ny);
  NumericVector SP (ny);
  NumericVector SSB2 (ny);
  NumericVector TSB (ny);
  NumericVector ESB (ny);
  NumericVector maty(asmax);
  NumericVector sely(asmax);
  NumericVector weighty(asmax);
  NumericVector weightFy(asmax);
  double rec = 0.0;
  double fecun = 1.0;
  double beta = 0.0;
  double hy = 0.0;
  double R0y = 0.0;
  double Ctmp = 0.0;
  double fs = FM/ns;
  double SSB = 0.0;
  double SPR = 0.0;

  // errors
  //  NumericVector eF = as<NumericVector>(set["eF"]);
  NumericVector eR = as<NumericVector>(set["eR"]);
  NumericVector eM = as<NumericVector>(set["eM"]);
  NumericVector eH = as<NumericVector>(set["eH"]);
  NumericVector eR0 = as<NumericVector>(set["eR0"]);
  NumericVector eMat = as<NumericVector>(set["eMat"]);
  NumericVector eSel = as<NumericVector>(set["eSel"]);
  NumericVector eW = as<NumericVector>(set["eW"]);
  //  NumericVector eImp = as<NumericVector>(set["eImp"]);


  // Initialise
  std::fill( CW.begin(), CW.end(), 0);
  std::fill( SP.begin(), SP.end(), 0);
  std::fill( SSB2.begin(), SSB2.end(), 0);
  std::fill( TSB.begin(), TSB.end(), 0);
  std::fill( ESB.begin(), ESB.end(), 0);
  NumericMatrix Msel = as<NumericMatrix>(MselList[tvmsel-1]);
  NumericMatrix sel = as<NumericMatrix>(selList[tvsel-1]);

  for(int a=0; a<asmax; a++){
    MAA0(a) = M(as2s(a) + tvm-1) * Msel(as2a(a),as2s(a));
    //     for(int a=0; a<asmax; a++) MAA(a,_) = Mtmp(a,_) * M[Rcpp::Range(s1vec(y), s1vec(y)+ns-1)];
    FAA(a) = fs * sel(as2a(a),as2s(a));
  }

  NumericVector NAAS = initdist(MAA0 * eM(1), FAA, R0 * eR0(1), spawning, inds);

  // Years
  for(int y=0; y<ny; y++){
    // Adding noise
    hy = h * eH(y);
    R0y = R0 * eR0(y);
    for(int a=0; a<asmax; a++){
      maty(a) = mat(as2a(a),as2s(a)) * eMat(y);
      sely(a) = sel(as2a(a),as2s(a)) * eSel(y);
      weighty(a) = weight(as2a(a),as2s(a)) * eW(y);
      weightFy(a) = weightF(as2a(a),as2s(a)) * eW(y);
    }
    MAA = MAA0 * eM(y);
    ZAA = MAA + FAA;


    // Seasons
    for(int s=0; s<ns; s++){

      // Recruitment
      SSB = 0.0;
      for(int a=0; a<asmax; a++){
        SSB += NAAS(a) * maty(a) * weighty(a) * exp(-pzbm * ZAA(a));
      }
      // SR
      if(SR == "bevholt"){
        NAAStmp = initdist(MAA, FAA, 1, spawning, inds);
        SPR = 0.0;
        for(int a=0; a<asmax; a++){
          SPR += NAAStmp(a) * maty(a) * weighty(a) * fecun; // SPR(s) += NnatM(a,s) * maty(a,s) * weighty(a,s) * fecun;
        }
        rec = 4 * hy * R0y * SSB / (SPR * R0y * (1-hy) + SSB * (5*hy-1));
      }else if(SR == "ricker"){
        rec = bp * SSB * exp(-beta * SSB);
      }else if(SR == "average"){
        rec = R0y;
      }else if(SR == "hockey-stick"){
        if(SSB > bp){
          rec = R0y;
        }else{
          rec = SSB * R0y/bp;
        }
      }else if(SR == "bent-hyperbola"){
        rec = recBeta * (SSB + sqrt(pow(bp,2) + pow(recGamma,2)/4) -
                         sqrt(pow(SSB-bp,2) + pow(recGamma,2)/4));
      }

      NAAS(0) = rec * spawning(s) * eR(y);

      // Catches
      CAA = FAA/ZAA * NAAS * (1 - exp(-ZAA));
      Ctmp = 0.0;
      for(int a=0; a<asmax; a++){
        Ctmp += CAA(a) * weightFy(a);
      }
      CW(y) += Ctmp;

      // Expenential decay
      NAAS = NAAS * exp(-ZAA);
      if(s == (ns-1)){
        for(int a=0; a<asmax; a++){
          // biomasses
          Bage(a) = NAAS(a) * weighty(a);
          SSBage(a) = Bage(a) * maty(a) * exp(-pzbm * ZAA(a));
          ESBage(a) = Bage(a) * sely(a);
          TSB(y) += Bage(a);
          SSB2(y) += SSBage(a);
          ESB(y) += ESBage(a);
        }
      }
      NAAS(asmax-1) = NAAS(asmax-1) + NAAS(asmax-2);
      for(int a=(asmax-1); a>0; a--){
        NAAS(a) = NAAS(a-1);
      }
      NAAS(0) = 0;
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
