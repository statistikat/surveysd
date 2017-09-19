#include <Rcpp.h>
using namespace Rcpp;


// function for rolling mean
// [[Rcpp::export]]
NumericVector rollMeanC(NumericVector x, int k, char type) {

  int n = x.size();
  int n_k = n-k+1;

  /*
  // only used für leading or trailing missings (not yet implemented)
  int i_start,i_end;
  if(type=='c'){
    i_start = ceil((double) k/2)-1;
    i_end = n - i_start + 1;
  }else if(type=='r'){
    i_start = 0;
    i_end = n-k+1;
  }else if(type=='l'){
    i_start = k-1;
    i_end = n;
  }

  return(i_start);
  */
  NumericVector x_out(n_k);


  for(int i=0;i<n_k;i++){
    x_out[i] = mean(x[seq(i,i+k-1)]);
  }

  return x_out;
}


