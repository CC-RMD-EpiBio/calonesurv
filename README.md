# Prevalence-reweighted Kaplan Meier curves

This repository is an implmentation of a test statistic found by computing the area between cohort-weighted Kaplan-Meier curves. 
The asymptotic distribution of the test statistic is used to compute p-values. The manuscript describing this statistic
is currently under review at RSOS. Here is an arXiv version: [arXiv1701.02424](https://arxiv.org/abs/1701.02424)

Here is an example to get you started with the tool

```r
library(devtools)
install_github("joshchang/calonesurv")
require(survival)
require(calonesurv)

surv_data = with(subset(survival::lung,ph.ecog %in% 0:2), 
                 data.frame(population = sex, 
                            censor = as.numeric(status==1), 
                            time = time, cohort = ph.ecog ))

out = Theta_hat(surv_data)
print(confint(out))
print(pvalue.Theta_hat(out))

plot(out,main="",ylab = "")

```
