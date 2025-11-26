The link to Shiny App -
https://puneet-talwar.shinyapps.io/PosthocPowerShiny/

**Linear Multiple Regression: Power & Sensitivity Analysis**

-This app replicates and extends G*Power for post-hoc power and sensitivity analysis in linear multiple regression.

**Inputs**

-Sample size (N): Number of participants (must be > k + 1).

-Total predictors (k): Total predictors in the model.

-Tested predictors (u): Predictors of interest being tested (k for the full model)).

-Significance level (alpha): Type I error probability.

-Confidence level: Confidence interval for effect sizes.

-Analysis type: Choose between Target Power Analysis or Observed R² Analysis.

-Run Analysis: Press the button to compute results.

**Outputs**

-Numerical Results: Displays effect sizes (f, f², R², r, d) with confidence intervals.

-Power Curve: Shows power vs effect size (f²).

-Nomogram: Shows power curves across different sample sizes.

**How to Use**

-For sensitivity analysis: Enter N, k, u, alpha, and target power.

-For post-hoc analysis: Check 'Use observed R² or f², enter R² or f².

-Click 'Run Analysis'.

-Interpret results and plots for your study context.

**Notes**

-Ensure N > k + 1.

-Power curve assumes 1 tested predictor (u = 1).

-Confidence intervals are based on non-central F distributions.

-Use this app for replication of G*Power analyses in R, with added confidence interval reporting.-
