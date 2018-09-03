# Overview


## Package information

- Version: 0.1.0
- Maintainer : Youjin Lee (<ylee160@jhu.edu>)
- Imports : stats

## Installation

You can download the package by:

```
install.packages("logisticRR")

# or you can directly download the development version from author's Github 
install.packages("devtools")
library(devtools)
install_github("youjin1207/logisticRR")
```


## Usage

[Here](https://github.com/youjin1207/logisticRR/blob/master/vignettes/logisticRR.Rmd) is a R vignettes for guidance. Or you can access to vignettes via:



```
install_github("youjin1207/logisticRR", build_vignettes = TRUE)
vignette("logisticRR", package = "logisticRR")
```

## Example

```
library(logisticRR)
```

### generate hypothetical data

```
n <- 500
set.seed(1234)
X <- rbinom(n, 1, 0.3)
W <- rbinom(n, 1, 0.3); W[sample(1:n, n/3)] = 2
Z <- rep(0, n)
Z[sample(1:n, n/2)] <- "female"; Z <- ifelse(Z == 0, "male", Z)
dummyZ <- ifelse(Z == "female", 1, 0)
Y <- rbinom(n, 1, plogis(X - W + 2*dummyZ))
dat <- as.data.frame(cbind(Y, X, W, Z))
dat$X <- as.numeric(dat$X); dat$X <- ifelse(dat$X == 2, 1, 0)
dat$Y <- as.numeric(dat$Y); dat$Y <- ifelse(dat$Y == 2, 1, 0)
dat$W <- as.factor(dat$W)
dat$Z <- as.factor(dat$Z)
```

```
simresult <- logisticRR(Y ~ X + W + Z, data = dat, boot = TRUE, n.boot = 200)
var(simresult$boot.rr)
simresult$delta.var

simresult$RR
```
