# Analytic GRIM-U

{analyticgrimu} is a forensic meta-science R package. It allows you to recalculate *p*-values and U values for Mann-Whitney U-tests and compare them against the reported values, taking into account a) the rounding of the inputs and b) various analytic choices into account.

The GRIM-U method was first proposed by Heathers and Grimes ([Medical Evidence Project, Report 1, 2025](https://medicalevidenceproject.org/grim-u-observation-establish-impossible-p-values-ranked-tests/)). Grimes' original implementation of GRIM-U is available [here](https://github.com/drg85/GRIMU).



## Installation

Currently, recalc is only available on GitHub. To install it, run the following code in your R Console:

``` r
if (!require(remotes)) install.packages("remotes")
remotes::install_github("ianhussey/analtyicgrimu")
```



## License

(c) Ian Husey (2025)

MIT license
