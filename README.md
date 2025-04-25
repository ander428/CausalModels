
# CausalModels

<!-- badges: start -->
<!-- badges: end -->

The goal of CausalModels is to provide a survey of fundamental causal
inference models in one single location. While there are many packages
for these types of models, CausalModels brings them all to one place
with a simple user experience. The package uses a format that is
familiar to users using well known statistical and machine learning
models. While causal inference models require careful consideration of
variables to correctly infer a causal effect, CausalModels uses simple
code while requiring the user to make these considerations. This enables
efficient and thoughtful research using causal inference. As of May 30,
2022, the package has been [published on
CRAN](https://cran.r-project.org/package=CausalModels). \## Change Log

### Version 0.2.1 - 2025/04/27

#### Fixed

- generic treatment names bug (#7)

### Version 0.2.0 - 2022/10/27

#### Added

- one parameter g-estimation

#### Fixed

- Typo in ipweighting documentation

## Installation

You can install the development version of CausalModels from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ander428/CausalModels")
```

Since the package has been published on CRAN, the production version can
be installed with:

``` r
install.packages("CausalModels")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(CausalModels)
library(causaldata)

data(nhefs)

nhefs.nmv <- nhefs[which(!is.na(nhefs$wt82)),]
nhefs.nmv$qsmk <- as.factor(nhefs.nmv$qsmk)

confounders <- c("sex", "race", "age", "education", "smokeintensity",
                 "smokeyrs", "exercise", "active", "wt71")

# initialize package
?init_params
init_params(wt82_71, qsmk,
            covariates = confounders,
            data = nhefs.nmv, simple = F)
#> Successfully initialized!
#> 
#> Summary:
#> 
#> Outcome - wt82_71 
#> Treatment - qsmk 
#> Covariates - [ sex, race, age, education, smokeintensity, smokeyrs, exercise, active, wt71 ] 
#> 
#> Size - 1566 x 67 
#> 
#> Default formula for outcome models: 
#> wt82_71 ~ qsmk + sex + race + education + exercise + active + age + (qsmk * age) + I(age * age) + smokeintensity + (qsmk * smokeintensity) + I(smokeintensity * smokeintensity) + smokeyrs + (qsmk * smokeyrs) + I(smokeyrs * smokeyrs) + wt71 + (qsmk * wt71) + I(wt71 * wt71) 
#> 
#> Default formula for propensity models: 
#> qsmk ~ sex + race + education + exercise + active + age + I(age * age) + smokeintensity + I(smokeintensity * smokeintensity) + smokeyrs + I(smokeyrs * smokeyrs) + wt71 + I(wt71 * wt71) 

# mode the causal effect of qsmk on wt82_71
model <- standardization(nhefs.nmv)
print(model)
#> 
#> Call:  glm(formula = wt82_71 ~ qsmk + sex + race + education + exercise + 
#>     active + age + (qsmk * age) + I(age * age) + smokeintensity + 
#>     (qsmk * smokeintensity) + I(smokeintensity * smokeintensity) + 
#>     smokeyrs + (qsmk * smokeyrs) + I(smokeyrs * smokeyrs) + wt71 + 
#>     (qsmk * wt71) + I(wt71 * wt71), family = family, data = combined_data)
#> 
#> Coefficients:
#>                        (Intercept)                               qsmk1  
#>                         -0.9699812                           0.5509460  
#>                               sex1                               race1  
#>                         -1.4371844                           0.5868376  
#>                         education2                          education3  
#>                          0.8174769                           0.5824119  
#>                         education4                          education5  
#>                          1.5240890                          -0.1792422  
#>                          exercise1                           exercise2  
#>                          0.3063727                           0.3550789  
#>                            active1                             active2  
#>                         -0.9460683                          -0.2707615  
#>                                age                        I(age * age)  
#>                          0.3495673                          -0.0060652  
#>                     smokeintensity  I(smokeintensity * smokeintensity)  
#>                          0.0482197                          -0.0009597  
#>                           smokeyrs              I(smokeyrs * smokeyrs)  
#>                          0.1418662                          -0.0018076  
#>                               wt71                      I(wt71 * wt71)  
#>                          0.0393011                          -0.0009787  
#>                          qsmk1:age                qsmk1:smokeintensity  
#>                          0.0123138                           0.0448028  
#>                     qsmk1:smokeyrs                          qsmk1:wt71  
#>                         -0.0235529                           0.0291350  
#> 
#> Degrees of Freedom: 1565 Total (i.e. Null);  1542 Residual
#>   (3132 observations deleted due to missingness)
#> Null Deviance:       97180 
#> Residual Deviance: 82690     AIC: 10710
#> 
#> Average treatment effect of qsmk:
#> 3.4927 
summary(model)
#> 
#> Call:
#> glm(formula = wt82_71 ~ qsmk + sex + race + education + exercise + 
#>     active + age + (qsmk * age) + I(age * age) + smokeintensity + 
#>     (qsmk * smokeintensity) + I(smokeintensity * smokeintensity) + 
#>     smokeyrs + (qsmk * smokeyrs) + I(smokeyrs * smokeyrs) + wt71 + 
#>     (qsmk * wt71) + I(wt71 * wt71), family = family, data = combined_data)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -41.913   -4.168   -0.314    3.869   44.573  
#> 
#> Coefficients:
#>                                      Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)                        -0.9699812  4.3673208  -0.222 0.824266    
#> qsmk1                               0.5509460  2.8229123   0.195 0.845286    
#> sex1                               -1.4371844  0.4693195  -3.062 0.002235 ** 
#> race1                               0.5868376  0.5828368   1.007 0.314158    
#> education2                          0.8174769  0.6085125   1.343 0.179339    
#> education3                          0.5824119  0.5575569   1.045 0.296382    
#> education4                          1.5240890  0.8351981   1.825 0.068221 .  
#> education5                         -0.1792422  0.7462118  -0.240 0.810205    
#> exercise1                           0.3063727  0.5360193   0.572 0.567697    
#> exercise2                           0.3550789  0.5592886   0.635 0.525603    
#> active1                            -0.9460683  0.4105673  -2.304 0.021338 *  
#> active2                            -0.2707615  0.6851128  -0.395 0.692745    
#> age                                 0.3495673  0.1648157   2.121 0.034084 *  
#> I(age * age)                       -0.0060652  0.0017347  -3.496 0.000485 ***
#> smokeintensity                      0.0482197  0.0518339   0.930 0.352375    
#> I(smokeintensity * smokeintensity) -0.0009597  0.0009409  -1.020 0.307878    
#> smokeyrs                            0.1418662  0.0943836   1.503 0.133023    
#> I(smokeyrs * smokeyrs)             -0.0018076  0.0015458  -1.169 0.242437    
#> wt71                                0.0393011  0.0836422   0.470 0.638514    
#> I(wt71 * wt71)                     -0.0009787  0.0005255  -1.862 0.062751 .  
#> qsmk1:age                           0.0123138  0.0670159   0.184 0.854238    
#> qsmk1:smokeintensity                0.0448028  0.0360169   1.244 0.213712    
#> qsmk1:smokeyrs                     -0.0235529  0.0654333  -0.360 0.718931    
#> qsmk1:wt71                          0.0291350  0.0276439   1.054 0.292075    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 53.62833)
#> 
#>     Null deviance: 97176  on 1565  degrees of freedom
#> Residual deviance: 82695  on 1542  degrees of freedom
#>   (3132 observations deleted due to missingness)
#> AIC: 10706
#> 
#> Number of Fisher Scoring iterations: 2
#> 
#> Average treatment effect of qsmk:
#>                            Estimate
#> Observed effect            2.638300
#> Counterfactual (treated)   5.243703
#> Counterfactual (untreated) 1.751003
#> Risk difference            3.492700
#> Risk ratio                 2.994686
#> 
```
