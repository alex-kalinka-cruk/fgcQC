context("test-gicc.R")

# Input dataset used in Ahrens (1976).
niris <- iris[,c("Species",colnames(iris[1:ncol(iris)-1]))]
gicc <- GICC(niris, return.cov = TRUE)

# Function to randomize subject IDs in multivariate dataset, calculate GICC.1,
#  and compare against true value.
rand.data.gicc <- function(data, num.reps){
  true_gicc <- GICC(data, return.cov = F)
  dord <- data[,1]
  i <- 1
  while(i <= num.reps){
    data[,1] <- sample(data[,1],nrow(data),replace = FALSE)
    # Make sure we don't have original data order.
    if(all(data[,1]==dord)) next
    rg <- GICC(data, return.cov=F, all.gicc = F)
    # Test both point estimate and 95% UB_rand < 95% LB_true.
    expect_lt(rg$GICC[[1]], true_gicc$GICC[[1]])
    expect_lt(rg$conf.int[2], true_gicc$conf.int[1])
    i <- i+1
  }
}


# Function to build simulated multi-variate dataset for input into GICC():
#  use random Gaussian samples with given ratio of within to between subj SD.
build.simul.mv_df <- function(num.subj, num.reps, num.var, w.sd, b.sd){
  # 'w.sd' - SD of Gaussian samples for within-subject reps.
  # 'b.sd' - SD of Gaussian samples for between-subject reps.
  subj <- as.character(rep(1:num.subj, each=num.reps))
  sm <- vector("list",num.subj)
  # Between subj means drawn from Gaussian (same across variables).
  sm <- lapply(sm, function(x) rnorm(1, 0, b.sd))
  # Replicates drawn from Gaussians with above Gaussian-sampled means.
  vars <- list()
  for(i in 1:num.var){
    vars[[i]] <- unlist(lapply(sm, function(x) rnorm(num.reps, x, w.sd)))
  }
  rd <- data.frame(subj=subj, vars)
  colnames(rd) <- c("subj",paste("V",1:num.var,sep=""))
  return(rd)
}


# Function to return GICC.1 95% CIs for simulated datasets.
simul.gicc <- function(num.subj, num.reps, num.var, N.sim, w.sd, b.sd){
  ret <- list()
  for(i in 1:N.sim){
    rd <- build.simul.mv_df(num.subj, num.reps, num.var, w.sd, b.sd)
    sg <- GICC(rd)
    # Point estimate, 95% LB, 95% UB.
    ret[[i]] <- c(sg$GICC[[1]], sg$conf.int)
  }
  return(ret)
}


### Tests ###
##
test_that("covariance matrices as expected", {
  tolerance <- 1e-10
  # Symmetric, square matrices.
  expect_equal(dim(gicc$cov.b), c(4,4))
  expect_equal(dim(gicc$cov.e), c(4,4))
  expect_true(isSymmetric(gicc$cov.b, tol=tolerance))
  expect_true(isSymmetric(gicc$cov.e, tol=tolerance))
})


test_that("sample numbers are correct", {
  expect_equal(gicc$num.var, 4)
  expect_equal(gicc$num.subj, 3)
  expect_equal(c(gicc$num.reps), c(setosa=50, versicolor=50, virginica=50))
  expect_true(gicc$balanced.design)
})

test_that("covariance matrices equal to Ahrens (1976) results for the iris dataset", {
  # Tolerance:
  # set to 1e-4:
  # 1. There are possibly small transcription differences: here, we work from the R data,
  #  but Ahrens would have worked directly with the Edgar (1935) dataset,
  #  or reproduced in Table 1 in Fisher (1936).
  # 2. Ahrens computed his values in 1976, unclear if this was by hand or not.
  tolerance <- 1e-4
  # Ahrens (1976) covariance matrices (diagonal + upper triangle, column-major) 
  #  for between and within subjects.
  # Reported in:
  # Ahrens, H. (1976). Multivariate variance-covariance components (MVCC) 
  #   and generalized intraclass correlation coefficient (GICC). Biom. Z. 18, 527-533.
  cov.bsubj <- c(0.626799,-0.201354,0.111193,1.64915,-0.573505,4.367395,0.712032,
                 -0.229954,1.866947,0.803261)
  cov.wsubj <- c(0.265034,0.092721,0.115374,0.167483,0.055238,0.18517,0.038367,0.032721,
                 0.042653,0.041905)
  expect_equal(gicc$cov.b[upper.tri(gicc$cov.b, diag=T)], cov.bsubj, tolerance=tolerance)
  expect_equal(gicc$cov.e[upper.tri(gicc$cov.e, diag=T)], cov.wsubj, tolerance=tolerance)
})


test_that("randomized iris subject IDs produce lower GICC than true data", {
  # Scramble species IDs in 'niris', calclulate GICC, and compare against true value.
  rand.data.gicc(niris, num.reps = 100)
})


test_that("GICC.1 is the same as in the Ahrens (1976) paper", {
  tolerance <- 1e-4
  # As reported in Ahrens (1976)
  gicc.1_ahrens1976 <- 0.917835
  expect_equal(gicc$GICC[[1]], gicc.1_ahrens1976, tolerance=tolerance)
})


test_that("perfect reproducibility", {
  # Input data: 3 subject IDs, 3 replicates, 2 variables, no error variance:
  # expect maximum reproducibility of 1.
  dpr <- data.frame(subj=as.character(rep(1:3,each=3)), 
                    v1=rep(c(0.5,1,1.5),each=3), 
                    v2=rep(c(9.1,21.3,11.2),each=3), stringsAsFactors = F)
  pgicc <- GICC(dpr)
  expect_equal(pgicc$GICC[[1]], 1)
})


test_that("zero reproducibility", {
  # Input data: 3 subject IDs, 3 replicates, 2 variables, no between-subject variance:
  # expect minimum reproducibility of 0.
  dzr <- data.frame(subj=as.character(rep(1:3,each=3)), 
                    v1=rep(c(11.19,2.89,7.654),3), 
                    v2=rep(c(0.96,0.25,119.2),3), stringsAsFactors = F)
  zgicc <- GICC(dzr)
  expect_equal(zgicc$GICC[[1]], 0)
})


test_that("covariances have the expected sign", {
  tolerance <- 1e-10
  # Input data: 3 subject IDs, 3 replicates, 2 variables.
  # expect positive covariance.
  dc <- data.frame(subj=as.character(rep(1:3,each=3)), 
                   v1=rep(c(1.1,7.9,21.8),each=3), 
                   v2=rep(c(0.12,4.13,5.19),each=3), stringsAsFactors = F)
  cg <- GICC(dc)
  expect_true(cg$cov.b[1,2] > 0)
  # expect negative covariance.
  dc$v2 <- rev(dc$v2)
  cg <- GICC(dc)
  expect_true(cg$cov.b[1,2] < 0)
  # expect zero covariance.
  dc$v2 <- c(1,1,-9/21.8)/sqrt(1+1+(-9/21.8)^2)
  cg <- GICC(dc)
  expect_equal(cg$cov.b[1,2], 0, tolerance=tolerance)
})


test_that("95% CI estimate contains simulated GICC.1 mean for ~95% of simulations", {
  # Simulate large-ish number of subjects/replicates 1e3 times.
  sg <- simul.gicc(30,10,3,1e3,0.95,1)
  # Mean GICC in simulations.
  msg <- mean(unlist(lapply(sg,function(x) x[1])))
  # Percent CIs bounding the simulated mean.
  pb <- round(100*sum(unlist(lapply(sg,function(x) msg>=x[2] & msg<=x[3])))/1e3)
  # <3% from 95% being contained in the estimated CIs?
  expect_gt(pb, 92)
})

