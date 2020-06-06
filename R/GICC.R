# Function to calculate Ahrens (1976) GICC.

# Confidence Interval for ICC values (using Fisher's approximation for distribution of rho).
.icc_ci <- function(rho, k, n, alpha){
  num <- 2*((1-rho)^2)*(1+(n-1)*rho)^2
  den <- k*n*(n-1)
  ci.width <- 2*qnorm(alpha/2,lower.tail = F)*sqrt(num/den)
  ci <- c(max(rho - ci.width/2,0), min(rho + ci.width/2,1))
  return(ci)
}


#' GICC
#'
#' Calculates Ahrens (1976) multivariate generalized intraclass correlation coefficient (GICC).
#'
#' @param data A matrix or data frame of multivariate data in which the first column contains between-subject identifiers and the remaining columns are the variables measured on these subjects.
#' @param return.cov A logical indicating whether to return the covariance matrix or not. Defaults to `TRUE`.
#' @param all.gicc A logical indicating whether to return all of Ahrens GICC measures. Defaults to `FALSE`.
#' @param conf.level A numeric value between 0 and 1 indicating the desired confidence level for the GICC confidence interval.
#'
#' @return A list with the following components:
#' `GICC` - A list with the various GICC measures from Ahrens original paper. Use `GICC.1` for the one advocated by Ahrens.
#' `num.var` - An integer giving the number of variables.
#' `num.samps` - An integer giving the number of samples.
#' `num.subj` - An integer giving the number of subjects.
#' `num.reps` - An integer giving the number of replicates.
#' `conf.level` - The confidence level for the CI.
#' `balanced.design` - A logical indicating if the design is balanced or not.
#' @md
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom matrixcalc direct.sum matrix.trace
#' @references Ahrens, H. (1976). Multivariate variance-covariance components (MVCC) and generalized intraclass correlation coefficient (GICC). Biom. Z. 18, 527-533.
#' @export
GICC <- function(data, return.cov = TRUE, all.gicc = FALSE, conf.level = 0.95){
  # data - columns are measured variables, rows are replicates.
  # data - first column contains between-subject identifiers.
  #
  # Y is the data: an N*p matrix (p - variables, N - total samples) - see 'iris' dataset.
  # The rest are 'relationship' matrices.
  # A is the 'direct sum' of the between-subject 1-matrices -
  #  if there are n within-subject measures, it's an n*n block matrix composed of 1/n and 0 entries only.
  # Tt is an N*N matrix of 1/N entries.
  # I is the N*N Identity Matrix.
  # H is the Mulivariate between-subject SS.
  # G is the Multivariate error SS.
  #
  ## Check input data.
  if(!inherits(data[,1],"factor") && !inherits(data[,1],"character")){
    stop("GICC: first column of input data must contain subject IDs,
         columns 2 onwards must contain numerical data")
  }
  # Is the data balanced by subject?
  bs.count <- table(data[,1])

  # Does at least one of the subjects have a minimum of 2 replicates?
  if(max(bs.count) < 2){
    warning(paste("at least one subject must have >= 2 replicates:",data[,1]))
    return(list(GICC = list(GICC.1 = NA)))
  }

  rep.by.subj <- length(unique(bs.count))
  # Data must be ordered by subject IDs.
  data <- data[order(data[,1]),]
  tryCatch(
    {
    num.subj <- length(bs.count)
    Y <- as.matrix(data[,2:ncol(data)])
    },
    error = function(e) stop(paste("GICC: error handling data:",e))
  )
  if(!is.numeric(Y)) stop("GICC: data input error: columns 2 onwards should contain numerical data")
  N <- nrow(Y)
  if(N==0) stop("GICC: no data found in input")
  num.var <- ncol(Y)
  ## Construct relationship matrices.
  ml <- list()
  for(i in 1:num.subj){
    ml[[i]] <- matrix(1/bs.count[i], bs.count[i], bs.count[i])
  }
  A <- Reduce(matrixcalc::direct.sum, ml)
  Tt <- matrix(1/N, N, N)
  I <- diag(N)
  kappa <- (1/(num.subj-1))*(N-(1/N)*sum(bs.count^2))
  H <- (1/(num.subj-1))*(t(Y) %*% (A - Tt) %*% Y)
  Cov.e <- (1/(N-num.subj))*(t(Y) %*% (I - A) %*% Y)
  Cov.a <- (1/kappa)*(H-Cov.e)
  # If any(var_w > var_b), unbiased estimate leads to negative variances:
  # poses a problem for GICC.1, as it is defined with absolute values.
  # In these cases, we do not subtract Cov.e.
  if(any(diag(H)-diag(Cov.e)<0)){
    Cov.a <- (1/kappa)*H
  }
  ## Use rho.1 from Ahrens (1976) for GICC.
  sa <- sum(abs(Cov.a))
  se <- sum(abs(Cov.e))
  GICC <- list(GICC.1 = sa/(sa+se))
  if(all.gicc){
    ## Calculate remaining GICC variants defined by Ahrens (1976).
    sa <- sum(Cov.a)
    se <- sum(Cov.e)
    GICC$GICC.2 <- sa/(sa+se)
    R <- (H - Cov.e) %*% solve(H + (kappa-1)*Cov.e)
    GICC$GICC.3 <- max(eigen(R)$val)
    GICC$GICC.4 <- matrixcalc::matrix.trace(R)
    GICC$GICC.5 <- sum(R)
  }
  if(return.cov){
    ret <- list(var.b = diag(Cov.a), var.e = diag(Cov.e),
                cov.b = Cov.a, cov.e = Cov.e,
                GICC = GICC)
  }else{
    ret <- list(var.b = diag(Cov.a), var.e = diag(Cov.e),
                GICC = GICC)
  }
  ret$num.var <- num.var
  ret$num.samps <- N
  ret$num.subj <- num.subj
  ret$num.reps <- bs.count
  ret$conf.level <- conf.level
  ret$balanced.design <- rep.by.subj == 1
  # For CI, use smallest rep size in case of unbalanced data.
  if(!ret$balanced.design){
    num.reps <- min(bs.count)
  }else{
    num.reps <- bs.count[1]
  }
  ret$conf.int <- .icc_ci(ret$GICC[[1]], num.subj, num.reps, 1-conf.level)
  class(ret) <- "GICC"
  return(ret)
}
