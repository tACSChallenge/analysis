Bayesian analyses
written by Anne-Fiona Griesfeller:

The registered analysis plan includes Bayes factors for every
hypothesis, but does not specify how they are to be computed.
What was done in the scripts is therefore a proposal:

H4 to H6 are computed inside the main pipeline. Those Bayes factors compare
the registered model against the same model without the effect of interest,
using the BIC approximation (Wagenmakers, 2007).

H1 and H2 are separated because they need R and the BayesFactor package.
They rest on an estimator that is not specified anywhere in the protocol:
the observed regression coefficient minus the mean of that participant's own permutation distribution, tested against zero.

The subtraction is necessary because the coefficient is a magnitude and
cannot be negative. Subtracting each participant's own permuted mean re-centers them on their own chance
level, so that a one-sample test against zero is meaningful. This is a proposed method by AG. 
 
Note that the two are not on the same prior. The BIC approximation implies a
unit information prior; the scripts here use the JZS prior with the
BayesFactor default Cauchy scale of sqrt(2)/2.
 

How to run:

  1. Run the H1 to H3 analysis so that H1H2H3_pooled_metrics.mat exists.
  2. Run tACS_exportBayesData.m, which writes bayes_H1H2_delta.csv.
  3. Run tACS_bayes_H1H2.R on that file.


References: 

Morey, R. D., & Rouder, J. N. (2024). BayesFactor: Computation of Bayes factors for common designs. R package version 0.9.12-4.7.
Wagenmakers, E.-J. (2007). A practical solution to the pervasive problems of p values. Psychonomic Bulletin & Review, 14(5), 779–804. https://doi.org/10.3758/BF03194105