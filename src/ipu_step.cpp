#include <Rcpp.h>
using namespace Rcpp;

//' Perform one step of iterative proportional updating
//'
//' C++ routines to invoke a single iteration of the Iterative proportional updating (IPU) scheme. Targets and classes
//' are assumed to be one dimensional in the `ipf_step` functions. `combine_factors` aggregates several vectors of
//' type factor into a single one to allow multidimensional ipu-steps. See examples.
//'
//' `ipf_step` returns the adjusted weights. `ipf_step_ref` does the same, but updates `w` by reference rather than
//' returning. `ipf_step_f` returns a multiplicator: adjusted weights divided by unadjusted weights. `combine_factors` is
//' designed to make `ipf_step` work with contingency tables produced by [xtabs].
//'
//' @md
//' @name ipf_step
//' @param w a numeric vector of weights. All entries should be positive.
//' @param classes a factor variable. Must have the same length as `w`.
//' @param targets key figure to target with the ipu scheme. A numeric verctor of the same length as `levels(classes)`.
//' This can also be a `table` produced by `xtabs`. See examples.
//' @examples
//'
//' ############# one-dimensional ipu ##############
//'
//' ## create random data
//' nobs <- 10
//' classLabels <- letters[1:3]
//' dat = data.frame(
//'   weight = exp(rnorm(nobs)),
//'   household = factor(sample(classLabels, nobs, replace = TRUE))
//' )
//' dat
//'
//' ## create targets (same lenght as classLabels!)
//' targets <- 3:5
//'
//' ## calculate weights
//' new_weight <- ipf_step(dat$weight, dat$household, targets)
//' cbind(dat, new_weight)
//'
//' ## check solution
//' xtabs(new_weight ~ dat$household)
//'
//' ## calculate weights "by reference"
//' ipf_step_ref(dat$weight, dat$household, targets)
//' dat
//'
//' ############# multidimensional ipu ##############
//'
//' ## load data
//' factors <- c("time", "sex", "smoker", "day")
//' tips <- data.frame(sex=c("Female","Male","Male"), day=c("Sun","Mon","Tue"),
//' time=c("Dinner","Lunch","Lunch"), smoker=c("No","Yes","No"))
//' tips <- tips[factors]
//'
//' ## combine factors
//' con <- xtabs(~., tips)
//' cf <- combine_factors(tips, con)
//' cbind(tips, cf)[sample(nrow(tips), 10, replace = TRUE),]
//'
//' ## adjust weights
//' weight <- rnorm(nrow(tips)) + 5
//' adjusted_weight <- ipf_step(weight, cf, con)
//'
//' ## check outputs
//' con2 <- xtabs(adjusted_weight ~ ., data = tips)
//' sum((con - con2)^2)
//'
//' @rdname ipf_step
//' @export
// [[Rcpp::export]]
void ipf_step_ref(NumericVector w, IntegerVector classes, NumericVector targets) {
  CharacterVector levels(classes.attr("levels"));
  int nclasses = levels.size();
  if(targets.length() != nclasses)
    stop("number of levels does not match the length of targets");
  NumericVector targets2(targets.length());
  for(int i = 0; i < w.size(); i++){
    int cl = classes[i];
    if (cl >= 0)
      targets2[cl - 1] += w[i];
  }
  for(int i = 0; i < w.size(); i++){
    int cl = classes[i];
    if (cl >= 0) {
      w[i] *= targets[cl - 1];
      w[i] /= targets2[cl - 1];
    }
  }
}

//' @rdname ipf_step
//' @export
// [[Rcpp::export]]
NumericVector ipf_step(NumericVector w, IntegerVector classes, NumericVector targets){
  NumericVector w_copy( Rcpp::clone( w ) );
  ipf_step_ref(w_copy, classes, targets);
  return w_copy;
}

//' @rdname ipf_step
//' @export
// [[Rcpp::export]]
NumericVector ipf_step_f(NumericVector w, IntegerVector classes, NumericVector targets){
  CharacterVector levels(classes.attr("levels"));
  int nclasses = levels.size();
  if(targets.length() != nclasses)
    stop("number of levels does not match the length of targets");
  NumericVector targets2(targets.length());
  for(int i = 0; i < w.size(); i++){
    int cl = classes[i];
    if (cl >= 0)
      targets2[cl - 1] += w[i];
  }
  NumericVector adjustments(classes.size());
  for(int i = 0; i < classes.size(); i++){
    int cl = classes[i];
    if(cl >= 0)
      adjustments[i] = targets[cl - 1]/targets2[cl - 1];
    else
      adjustments[i] = 1;
  }
  return adjustments;
}
