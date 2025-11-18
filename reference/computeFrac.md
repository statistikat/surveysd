# Numerical weighting functions

Customize weight-updating within factor levels in case of numerical
calibration. The functions described here serve as inputs for
[ipf](https://statistikat.github.io/surveysd/reference/ipf.md).

## Usage

``` r
computeLinear(curValue, target, x, w, boundLinear = 10)

computeLinearG1_old(curValue, target, x, w, boundLinear = 10)

computeLinearG1(curValue, target, x, w, boundLinear = 10)

computeFrac(curValue, target, x, w)
```

## Arguments

- curValue:

  Current summed up value. Same as `sum(x*w)`

- target:

  Target value. An element of `conP` in
  [ipf](https://statistikat.github.io/surveysd/reference/ipf.md)

- x:

  Vector of numeric values to be calibrated against

- w:

  Vector of weights

- boundLinear:

  The output `f` will satisfy `1/boundLinear <= f <= boundLinear`. See
  `bound` in
  [ipf](https://statistikat.github.io/surveysd/reference/ipf.md)

## Value

A weight multiplier `f`

## Details

`computeFrac` provides the "standard" IPU updating scheme given as

\$\$f = target/curValue\$\$

which means that each weight inside the level will be multtiplied by the
same factor when doing the actual update step (`w := f*w`).
`computeLinear` on the other hand calculates `f` as

f_(i) = a · x_(i) + b

where `a` and `b` are chosen, so f satisfies the following two
equations.

∑ f_(i) w_(i) x_(i) = target

∑ f_(i) w_(i) = ∑ w_(i)

`computeLinearG1` calculates `f` in the same way as `computeLinear`, but
if `f_i*w_i<1` `f_i` will be set to `1/w_i`.
