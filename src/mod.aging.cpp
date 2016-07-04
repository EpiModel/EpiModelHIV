#include <Rcpp.h>
using namespace Rcpp;

//' @title Aging Module
//'
//' @description Module for aging over time for nodes in the population.
//'
//' @param dat Master data list object of class \code{dat} containing networks,
//'        individual-level attributes, and summary statistics.
//' @param at Current time step.
//'
//' @return
//' This function returns \code{dat} after updating the nodal attribute
//' \code{age} and \code{sqrt.age}. The \code{sqrt.age} vertex attribute is also
//' updated on the three networks.
//'
//' @keywords module msm
//' @export
// [[Rcpp::export]]
List aging_msm(List dat, int at) {

  List attr = dat["attr"];
  List param = dat["param"];

  DoubleVector age = attr["age"];
  double tUnit = param["time.unit"];

  DoubleVector nage = age + (tUnit/365);

  attr["age"] = nage;
  attr["sqrt.age"] = sqrt(nage);

  dat["attr"] = attr;

  return dat;
}


//' @title Aging Module
//'
//' @description This module ages all nodes in the population by one time
//'              unit at each time step.
//'
//' @param dat Master data list object of class \code{dat} containing networks,
//'        individual-level attributes, and summary statistics.
//' @param at Current time step.
//'
//' @keywords module het
//' @export
// [[Rcpp::export]]
List aging_het(List dat, int at) {

  List attr = dat["attr"];
  List param = dat["param"];

  DoubleVector age = attr["age"];
  double tUnit = param["time.unit"];

  DoubleVector nage = age + (tUnit/365);

  attr["age"] = nage;
  dat["attr"] = attr;

  return dat;
}



