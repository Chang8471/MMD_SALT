
# Version 0.0.14

## Minor changes

* change `calc_S_LR` to return both the calculated test statistics and the estimated $p_1$ under the alternative

# Version 0.0.13

## Major fix

* change `calc_count_cutoff` to expand uniroot upper end if no solution

# Version 0.0.12

## Major changes

* new functions related proposed new likelihood ratio test, that takes into consideration of sample-specific sensitivities

# Version 0.0.11

## Major changes

* new function, instead of using `forceMono`, now calculate the desired count to reach a prespecified postL cutoff with `calc_count_cutoff`

* new function, `est_sensitivity` to calculate sample specific sensitivity

# Version 0.0.10

## Major changes

* change `sim_DP_sampling` function to sample according to order statistics of alpha_i

# Version 0.0.9

## Major changes

* new function `sim_DP_sampling` to simulate DP by splitting 100 samples into 2 groups

# Version 0.0.8

## Major changes

* change `sim_DP` function so that the null genes no longer have fixed prior for Z~Bernolli(p_0)

* change `sim_DP` function so that the alpha_i shift is done be moving one group down and the other up

# Version 0.0.7

## Major changes

* add function to simulate data with differential presence and unbalanced sample effect.

* add function to simulate data according to estimated parameters from real data fitting.

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

