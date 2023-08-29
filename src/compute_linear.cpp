#include <Rcpp.h>
using namespace Rcpp;

//' @rdname computeFrac
//' @export
// [[Rcpp::export]]
NumericVector computeLinear(double curValue,
                            double target,
                            const NumericVector& x,
                            const NumericVector& w,
                            double boundLinear = 10) {
  NumericVector f(x.size());
  if(x.size()!=w.size()){
    stop("x and w of different length!");
  }
  if(x.size()==1){
    f[0] = w[0]/curValue*target;
    return f;
  }else{
    double h = 0.0;
    double j = 0.0;
    double N = 0.0;

    for(int i = 0; i < x.size(); i++){
      h += w[i]*x[i];
      j += w[i]*x[i]*x[i];
      N += w[i];
    }

    double b = (target-N*j/h)/(h-N*j/h);
    double a = (N-b*N)/h;


    for(int i = 0; i < x.size(); i++)
      f[i] = a*x[i] + b;

    //apply bounds
    if(boundLinear > 0){
      for(int i = 0; i < x.size(); i++){
        if (f[i] < 1.0/boundLinear)
          f[i] = 1.0/boundLinear;
        if (f[i] > boundLinear)
          f[i] = boundLinear;
      }
    }

    return f;
  }
}


//' @rdname computeFrac
//' @export
// [[Rcpp::export]]
NumericVector computeLinearG1_old(double curValue,
                              double target,
                              const NumericVector& x,
                              const NumericVector& w,
                              double boundLinear = 10) {
  NumericVector f(x.size());
  f=computeLinear(curValue,target,x,w,boundLinear);
  
  for(int i = 0; i < x.size(); i++){
    if (f[i]*w[i] < 1.0){
      f[i]=1/w[i];
    }
  }

  return f;
}

//' @rdname computeFrac
//' @export
// [[Rcpp::export]]
NumericVector computeLinearG1(double curValue,
                              double target,
                              const NumericVector& x,
                              const NumericVector& w,
                              double boundLinear = 10) {
  NumericVector f(x.size());
  f=computeLinear(curValue,target,x,w,boundLinear);
  
  // check if any updated weight will be smaller than 1
  LogicalVector is_smaller = f*w-1.0 < -1e-10;
  LogicalVector is_larger = f*w-1.0 > 1e-10;
  bool cond1 = is_true(any(is_smaller));
  bool adjust_larger = true;
  
  // calculate share of weights and x*w which will be added by truncating at 1
  // and redistribute difference to rest
  // int steps = 0; // number of iteration for update
  while(cond1 == true){
    
    NumericVector f_smaller = f[is_smaller==true];
    NumericVector x_smaller = x[is_smaller==true];
    NumericVector w_smaller = w[is_smaller==true];

    // all elements of x and w where f*x ist greater 1
    NumericVector x_larger = x[is_larger == true];
    
    // weighted sum of x that is added if f[i] are truncated by 1
    double target_new = sum((1.0-f_smaller*w_smaller)*x_smaller);

    // sum of w that is added if f[i] are truncated by 1
    double sum_w_smaller = sum(1.0-f_smaller*w_smaller);
    NumericVector w_new = Rcpp::rep(sum_w_smaller/x_larger.size(),x_larger.size()); 
    double curValue_new = sum(w_new*x_larger);
    NumericVector f_new = computeLinear(curValue_new,target_new,x_larger,w_new,0.0);
 
    //update vector
    // check if any weight will be negativ
    // this can happen due to instability
    NumericVector f_larger_update(f_new.size());
    int j = 0;
     for(int i=0;i<f.size();i++){
      if(is_smaller[i] == true){
        // set f such that w*f = 1
        f[i] = 1.0/w[i];
      }
      // if f*w > 1 update as well
      if(is_larger[i]==true){
        // set f to incorporate new_target and w_rest
        f_larger_update[j] = f[i] - f_new[j]*w_new[j]/w[i];
        if(f_larger_update[j]<0){
          adjust_larger = false;
        }
        j = j + 1;
      }
    }
     
 
     // if resulting weights will be negaitve break routine
     // and return current result
     if(adjust_larger == false){
       break;
     }else{
       f[is_larger == true] = f_larger_update;
     }

    // check if all f*w are >= 1
    is_smaller = f*w-1.0 < -1e-10;
    is_larger = f*w-1.0 > 1e-10;
    // update condition
    cond1 = is_true(any(is_smaller));
    // steps = steps + 1;
    // std::cout << "step: "<< steps<<"\n";
  }

  return f;
}
