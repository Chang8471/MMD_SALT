# Version 0.0.6

## Major changes

* add functions for evaluating normalization 

# Version 0.0.5

## Major changes

* add an option to `parEst_itr` function, to either generate Z labels from a cutoff, or resample from Bernoulli(pi). Such that previous code will still run and get same result with updated package.

# Version 0.0.4

## Major changes

* change simulation setup: lower background component, center alpha_i

* change in parameter estimation: alpha_i constrain sum to 0

# Version 0.0.3

## Major changes

* add option to default calling to return z score other than p.value


## Major fix

* fixed `call_Z_logNormNegCtrl` function not using the given Ymtx


# Version 0.0.2

## Major changes

* changing monotonocity function, making low end having the same posterior likelihood as background mean


## Major fix

* fixed `call_Z_normNegCtrl` function not using the given Ymtx


# Version 0.0.1

## function added

* existing functions added (i.e. iteratively fitting model, fixing monotonically, simulate data) from 042322 file

