#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector vector_rcpp(int x) {

  NumericVector output(x);
  for(int i = 1; i<= x; i++){
    output[i-1] = i*i;
  }
  return output;

}


// [[Rcpp::export]]
List list_rcpp(int x) {

  List output(x);
  for(int i = 1; i<= x; i++){
    output[i-1] = i*i;
  }
  return output;

}
