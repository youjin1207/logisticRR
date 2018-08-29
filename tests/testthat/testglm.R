library(logisticRR)

test_that("Test whether logisticRR captures appropriate baseline covariates",{

  data("airquality")
  ozonedat <- na.omit(airquality)
  # define binary ozone level
  ozonedat$ozone1 <- ifelse(ozonedat$Ozone < quantile(ozonedat$Ozone, prob = 0.1), 1, 0)
  ozonedat$Temp2 <- ozonedat$Temp / 10
  ozonedat$Temp2.minus <- -ozonedat$Temp2

  ozone.fit.plus <- logisticRR(ozone1 ~ Temp2 + Solar.R + Wind, data = ozonedat, basecov = min(ozonedat$Temp2),
                           fixcov = data.frame(Solar.R = mean(ozonedat$Solar.R), Wind = mean(ozonedat$Wind)), boot = FALSE)

  ozone.fit.minus <- logisticRR(ozone1 ~ Temp2.minus + Solar.R + Wind, data = ozonedat, basecov = - min(ozonedat$Temp2) - 1,
                           fixcov = data.frame(Solar.R = mean(ozonedat$Solar.R), Wind = mean(ozonedat$Wind)), boot = FALSE)

  expect_equal(1/ozone.fit.plus$RR, ozone.fit.minus$RR)

})


test_that("Test whether nominalRR captures appropriate baseline covariates",{

  data("airquality")
  ozonedat <- na.omit(airquality)
  # define binary ozone level
  ozonedat$ozone1 <- ifelse(ozonedat$Ozone < quantile(ozonedat$Ozone, prob = 0.1), 1, 0)
  ozonedat$Temp.factor <- ifelse(ozonedat$Temp <= quantile(ozonedat$Temp, prob = 0.25), "low",
                                 ifelse(ozonedat$Temp > quantile(ozonedat$Temp, prob = 0.8), "high", "medium"))
  ozonedat$Temp.factor <- as.factor(ozonedat$Temp.factor)

  ozone.fit.factor <- nominalRR(ozone1 ~ Temp.factor + Solar.R + Wind, data = ozonedat,
                                basecov = "low", comparecov = "medium",
                                fixcov = data.frame(Solar.R = mean(ozonedat$Solar.R), Wind = mean(ozonedat$Wind)), boot = FALSE)

  ozone.fit.factor2 <- nominalRR(ozone1 ~ Temp.factor + Solar.R + Wind, data = ozonedat,
                                 basecov = "medium", comparecov = "low",
                                 fixcov = data.frame(Solar.R = mean(ozonedat$Solar.R), Wind = mean(ozonedat$Wind)), boot = FALSE)


  expect_equal(ozone.fit.factor$RR, 1/ozone.fit.factor2$RR)


  ozone.fit.nobase <- nominalRR(ozone1 ~ Temp.factor + Solar.R + Wind, data = ozonedat,
                                 fixcov = data.frame(Solar.R = mean(ozonedat$Solar.R), Wind = mean(ozonedat$Wind)), boot = FALSE)

  expect_equal(ozone.fit.nobase$RR, 1)

})



test_that("Test whether logisticRR and nominalRR have proper adjusted covariates as defaults",{


  data("airquality")
  ozonedat <- na.omit(airquality)
  # define binary ozone level
  ozonedat$ozone1 <- ifelse(ozonedat$Ozone < quantile(ozonedat$Ozone, prob = 0.1), 1, 0)
  ozonedat$Temp2 <- ozonedat$Temp / 10
  ozone.fit.mean <- logisticRR(ozone1 ~ Temp2 + Solar.R + Wind, data = ozonedat, basecov = min(ozonedat$Temp2),
                               fixcov = data.frame(Solar.R = mean(ozonedat$Solar.R), Wind = mean(ozonedat$Wind)), boot = FALSE)

  ozone.fit.median <- logisticRR(ozone1 ~ Temp2 + Solar.R + Wind, data = ozonedat,
                                fixcov = data.frame(Solar.R = median(ozonedat$Solar.R), Wind = median(ozonedat$Wind)), boot = FALSE)

  expect_equal(ozone.fit.mean$fix.cov$Solar.R, mean(ozonedat$Solar.R))
  expect_equal(ozone.fit.median$fix.cov$Solar.R, median(ozonedat$Solar.R))

  expect_equal(ozone.fit.mean$fix.cov$Wind, mean(ozonedat$Wind.R))
  expect_equal(ozone.fit.median$fix.cov$Wind, median(ozonedat$Wind.R))

})



test_that("Test whether nominalRR have proper adjusted covariates as defaults",{


  data("airquality")
  ozonedat <- na.omit(airquality)
  # define binary ozone level
  ozonedat$ozone1 <- ifelse(ozonedat$Ozone < quantile(ozonedat$Ozone, prob = 0.1), 1, 0)
  ozonedat$Temp.factor <- ifelse(ozonedat$Temp <= quantile(ozonedat$Temp, prob = 0.25), "low",
                                 ifelse(ozonedat$Temp > quantile(ozonedat$Temp, prob = 0.8), "high", "medium"))
  ozonedat$Temp.factor <- as.factor(ozonedat$Temp.factor)

  ozone.fit.factor.null <- nominalRR(ozone1 ~ Temp.factor + Solar.R + Wind, data = ozonedat,
                                basecov = "low", comparecov = "medium", boot = FALSE)

  expect_equal(ozone.fit.factor.null$fix.cov$Solar.R, min(ozonedat$Solar.R))
  expect_equal(ozone.fit.factor.null$fix.cov$Wind, min(ozonedat$Wind))

  ozone.fit.factor.mean <- nominalRR(ozone1 ~ Temp.factor + Solar.R + Wind, data = ozonedat,
                                     basecov = "low", comparecov = "medium",
                                     fixcov = data.frame(Solar.R = mean(ozonedat$Solar.R), Wind = mean(ozonedat$Wind)),
                                     boot = FALSE)

  expect_equal(ozone.fit.factor.mean$fix.cov$Solar.R, mean(ozonedat$Solar.R))
  expect_equal(ozone.fit.factor.mean$fix.cov$Wind, mean(ozonedat$Wind))

})


test_that("Test whether estimated variance using Delta method and sampling variance using bootstrap reasonably match each other",{

  data("airquality")
  ozonedat <- na.omit(airquality)
  # define binary ozone level
  ozonedat$ozone1 <- ifelse(ozonedat$Ozone < quantile(ozonedat$Ozone, prob = 0.3), 1, 0)
  ozonedat$Temp2 <- ozonedat$Temp / 10
  ozone.fit20 <- logisticRR(ozone1 ~ Temp2 + Solar.R + Wind, data = ozonedat, basecov = min(ozonedat$Temp2),
                          fixcov = data.frame(Solar.R = mean(ozonedat$Solar.R), Wind = mean(ozonedat$Wind)),
                          boot = TRUE, n.boot = 20)
  ozone.fit200 <- logisticRR(ozone1 ~ Temp2 + Solar.R + Wind, data = ozonedat, basecov = min(ozonedat$Temp2),
                             fixcov = data.frame(Solar.R = mean(ozonedat$Solar.R), Wind = mean(ozonedat$Wind)),
                             boot = TRUE, n.boot = 200)

  ozone.fit20$delta.var
  var(ozone.fit20$boot.rr)
  var(ozone.fit200$boot.rr)


  expect_equal(ozone.fit20$delta.var, ozone.fit200$delta.var)
  expect_true(abs(ozone.fit20$delta.var -  var(ozone.fit20$boot.rr)) >= abs(ozone.fit200$delta.var -  var(ozone.fit200$boot.rr)))

})
