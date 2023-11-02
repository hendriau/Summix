
#' summix
#'
#' @description
#' Estimating mixture proportions of reference groups from large (N SNPs>10,000) genetic AF data.
#'
#' @param data A dataframe of the observed and reference allele frequencies for N genetic variants. See data formatting document at \href{https://github.com/hendriau/Summix}{https://github.com/hendriau/Summix} for more information.
#' @param reference A character vector of the column names for the reference groups.
#' @param observed A character value that is the column name for the observed group.
#' @param pi.start Length K numeric vector of the starting guess for the reference group proportions. If not specified, this defaults to 1/K where K is the number of reference groups.
#' @param goodness.of.fit Default value is TRUE. If set as FALSE, the user will override the default goodness of fit measure and return the raw objective loss from slsqp.
#' @param override_removeSmallRef Default value is FALSE. If set as TRUE, the user will override the automatic removal of reference groups with <1% global proportions - this is not recommended.
#'
#' @return A data frame with the following columns:
#' @return goodness.of.fit: scaled objective loss from slsqp() reflecting the fit of the reference data. Values between 0.5-1.5 are considered moderate fit and should be used with caution. Values greater than 1.5 indicate poor fit, and users should not perform further analyses using Summix.
#' @return iterations: number of iterations for SLSQP algorithm
#' @return time: time in seconds of SLSQP algorithm
#' @return filtered: number of genetic variants not used in the reference group mixture proportion estimation due to missing values.
#' @return K columns of mixture proportions of reference groups input into the function
#'
#'
#' @author Adelle Price, \email{adelle.price@cuanschutz.edu}
#' @author Hayley Wolff, \email{hayley.wolff@cuanschutz.edu}
#' @author Audrey Hendricks, \email{audrey.hendricks@cuanschutz.edu}
#' @references https://github.com/hendriau/Summix2
#' @keywords genetics, mixture distribution, admixture, population stratification
#'
#' @seealso \url{https://github.com/hendriau/Summix2} for further documentation. \code{\link[nloptr]{slsqp}} function in the nloptr package for further details on Sequential Quadratic Programming \url{https://www.rdocumentation.org/packages/nloptr/versions/1.2.2.2/topics/slsqp}
#'
#' @examples
#' # load the data
#' data("ancestryData")
#'
#' # Estimate 5 reference ancestry proportion values for the gnomAD African/African American group
#' # using a starting guess of .2 for each ancestry proportion.
#' summix(data = ancestryData,
#'     reference=c("reference_AF_afr",
#'         "reference_AF_eas",
#'         "reference_AF_eur",
#'         "reference_AF_iam",
#'         "reference_AF_sas"),
#'     observed="gnomad_AF_afr",
#'     pi.start = c(.2, .2, .2, .2, .2),
#'     goodness.of.fit=TRUE)
#'
#' @export

summix <- function(data, reference, observed, pi.start = NA, goodness.of.fit = TRUE, override_removeSmallRef = FALSE) {
  start_time = Sys.time()
  # get the reference group proportions using the summix_calc helper function with initial pi.start specified
  if(length(pi.start)!=0) {
    
    if(!override_removeSmallRef) {
      # first run summix and remove any reference groups <1% if 'override_removalSmallRef' = FALSE
      globTest <- invisible(summix_calc(data = data, reference = reference,
                                        observed = observed, pi.start = pi.start))
      if(sum(globTest[1,5:ncol(globTest)] < 0.01) > 0) {
        toremove <- which(globTest[1,5:ncol(globTest)] < 0.01)
        reference <- reference[-toremove]
      }
      pi.start = NA
      # get the scaled reference group proportions (after removing initial reference groups with props <1%) using the summix_calc helper function
      #with initial pi.start specified
      sum_res <- summix_calc(data = data, reference = reference,
                             observed = observed)
    }
    
    
    # get the reference group proportions for all reference groups, if 'override_removeSmallRef' = TRUE, using the summix_calc helper function
    #with pi.start specified
    sum_res <- summix_calc(data = data, reference = reference,
                           observed = observed, pi.start = pi.start)
    
    
    
    
    # get the reference group proportions using the summix_calc helper function either without initial pi.start specified
  } else {
    if(!override_removeSmallRef) {
      # first run summix and remove any reference groups <1% if 'override_removalSmallRef' = FALSE
      globTest <- invisible(summix_calc(data = data, reference = reference, observed = observed))
      if(sum(globTest[1,5:ncol(globTest)] < 0.01) > 0) {
        toremove <- which(globTest[1,5:ncol(globTest)] < 0.01)
        reference <- reference[-toremove]
      }
      # get the scaled reference group proportions (after removing initial reference groups with props <1%) using the summix_calc helper function
      #without initial pi.start specified
      sum_res <- summix_calc(data = data, reference = reference, observed = observed)
      
    }
    
    
    # get the reference group proportions for all reference groups, if 'override_removeSmallRef' = TRUE, using the summix_calc helper function
    #without pi.start specified
    sum_res <- summix_calc(data = data, reference = reference, observed = observed)
  }
  
  
  
  # get the scaled objective and replace the objective argument (first column of 'summix' function output) with the updated objective value
  # will not run if overridden by user, but is the default
  if(goodness.of.fit) {
    new_obj <- calc_scaledObj(data = data, observed = observed, reference = reference, pi.start= pi.start)
    sum_res[1] <- new_obj
  }
  end_time = Sys.time()
  ttime = end_time - start_time
  sum_res[3] <- ttime
  
  # output caution/warning messages based on the objective and determined thresholds of moderate/poor fit
  if(goodness.of.fit) {
    if(sum_res[1] >= 0.5 & sum_res[1] < 1.5) {
      print(paste0("CAUTION: Objective/1000SNPs = ", round(sum_res[1], 4),
                   " which is within the moderate fit range"))
    } else if (sum_res[1] >= 1.5) {
      print(paste0("WARNING: Objective/1000SNPs = ", round(sum_res[1], 4),
                   " which is above the poor fit threshold"))
    }
  } else {
    if(sum_res[1]/nrow(data)/1000 >= 0.5 & sum_res[1]/nrow(data)/1000 < 1.5) {
      print(paste0("CAUTION: Objective/1000SNPs = ", round(sum_res[1]/nrow(data)/1000, 4),
                   " which is within the moderate fit range"))
    } else if (sum_res[1]/nrow(data)/1000 >= 1.5) {
      print(paste0("WARNING: Objective/1000SNPs = ", round(sum_res[1]/nrow(data)/1000, 4),
                   " which is above the poor fit threshold"))
    }
  }
  
  return(sum_res)
}













#' summix_calc
#'
#' @description
#'Helper function for estimating mixture proportions of reference groups from large (N SNPs>10,000) genetic AF data, using slsqp to solve for least square difference
#'
#' @param data A dataframe of the observed and reference allele frequencies for N genetic variants. See data formatting document at \href{https://github.com/hendriau/Summix}{https://github.com/hendriau/Summix} for more information.
#' @param reference A character vector of the column names for the reference groups.
#' @param observed A character value that is the column name for the observed group.
#' @param pi.start Length K numeric vector of the starting guess for the reference group proportions. If not specified, this defaults to 1/K where K is the number of reference groups.
#'
#' @return data frame with the following columns
#' @return objective: least square value at solution
#' @return iterations: number of iterations for SLSQP algorithm
#' @return time: time in seconds of SLSQP algorithm
#' @return filtered: number of SNPs not used in estimation due to missing values
#' @return K columns of mixture proportions of reference groups input into the function
#' @importFrom methods "is"
#' @export

summix_calc = function(data, reference, observed, pi.start=NA){
  if(!is(object = data, class2 = "data.frame")){
    stop("ERROR: data must be a data.frame as described in the vignette")
  }
  if(typeof(observed)!="character"){
    stop("ERROR: 'observed' must be a character string for the column name of the observed group in data")
  }
  if(!(observed %in% names(data))){
    stop("ERROR: 'observed' must be the column name of the observed group in data")
  }
  if(typeof(reference)!="character"){
    stop("ERROR: 'reference' must be a vector of column names from data to be used in the reference")
  }
  if(all(reference %in% names(data))==FALSE){
    stop("ERROR: 'reference' must be a vector of column names from data to be used in the reference")
  }
  data <- tibble::as_tibble(data)
  #Filter NA allele frequencies out of the observed column
  filteredNA <- length(which(is.na(data[[observed]]==TRUE)))
  # The slsqp algorithm only uses the observed allele frequency vector and the reference panel
  observed.b  <- as.data.frame(data[which(is.na(data[[observed]])==FALSE), observed, exact =TRUE])
  refmatrix <- as.data.frame(data[which(is.na(data[[observed]])==FALSE), reference, exact = TRUE])
  
  
  
  ##If pi.start is specified by user, check it meets criteria and print respective errors
  if(!is.na(pi.start)[1]){
    #########################
    # Check if pi.start is numeric
    if (is.numeric(pi.start)==FALSE){
      stop("ERROR: Please make sure pi.start is a numeric vector")
    }
    # Check if length(pi.start)==length(reference)
    if (length(pi.start) != length(reference)){
      stop("ERROR: Please make sure pi.start is the same length as reference")
    }
    # Check that pi.start is positive
    if (all(pi.start>0) == FALSE){
      stop("ERROR: Please make sure pi.start is a positive numeric vector")
    }
    # Check if sum(pi.start) = 1
    if (sum(pi.start)!=1){
      stop("ERROR: Please make sure pi.start sums to one")
    }
    #########################
    ###############################
    # Set the starting guess to pi.start if criteria is met
    ###############################
    starting = pi.start
    
    
    
    #if pi.start is not specified, starting defaults to 1/K where K is the number of reference groups.
  } else{
    starting = rep( 1/ncol(refmatrix), ncol(refmatrix) )
  }
  
  
  # Here we are defining the objective function. This function is evaluated at a
  # k-dimensional set x. Each of our K reference allele frequencies are multiplied by
  # our current best guess for the reference group proportion.
  # We then subtract the allele frequency values from the observed population.
  # And finally this sum is squared to achieve a least squares form.
  #########################
  fn.refmix = function(x){
    expected=x%*%t(as.matrix(refmatrix))
    minfunc = sum((expected - observed.b)**2)
    return(minfunc)
  }
  ########################
  # Here we are defining the gradient of the objective function.
  gr.refmix <- function(x){
    gradfunc = x%*%t(as.matrix(refmatrix)) - observed.b
    gradvec <- apply(2*refmatrix*t(gradfunc), 2, sum)
    return(gradvec)
  }
  
  # H equality
  # This function returns the equality constraints for the nloptr slsqp algorithm
  # We sum up the K current proportion estimate values and subtract 1. If the estimated proportion values sum to 1 than this value should equal zero.
  heq.refmix = function(x){
    equality = sum(x)
    return(equality - 1)
  }
  
  # H inequality
  # This function returns a vector of size K
  hin.refmix <- function(x){
    h = numeric(ncol(refmatrix))
    h=x
    return(h)
  }
  
  # We use the start_time function base function
  #to record the run time for our convex optimization algorithm.
  # The output for the nloptr slsqp function is stored in the variable S.
  # This function requires 5 inputs:
  #   1. fn.refmix is the objective function.
  #   2. gr.refmix is the gradient of the objective function.
  #   3. hin.refmix defining the inequality constraints
  #   4. heq.refmix defining the equality constraints
  start_time = Sys.time()
  S = suppressMessages( nloptr::slsqp(starting,
                                      fn = fn.refmix,
                                      gr = gr.refmix,
                                      hin = hin.refmix,
                                      heq = heq.refmix)
  )
  end_time = Sys.time()
  ttime = end_time - start_time
  
  # par is the K optimal reference group propotion estimates
  # value is the minimization function evaluated at par
  # iter is the number of iterations that the algorithm took to reach the optimal solution of par
  # finally ttime is the run time for the algorithm
  
  d <- data.frame(matrix(ncol = length(reference)+4, nrow = 1))
  colnames(d) <- c("goodness.of.fit", "iterations", "time",
                   "filtered", colnames(refmatrix))
  
  d[1] <- S$value
  d[2] <- S$iter
  d[3] <- ttime
  d[4] <- filteredNA
  
  #round reference group proportion estimates to 6 decimal places
  d[5:(length(reference)+4)] <- round(S$par[1:length(reference)], 6)
  
  #return data frame with objective, iterations, time, filtered, and K reference group proportion estimate columns
  return(d)
}






#' calc_scaledObj
#'
#' @description
#' Helper function to calculate new scaled loss function using weighted AF bin objectives
#'
#' @param data a dataframe of the observed and reference allele frequencies for N genetic variants. See data formatting document at \href{https://github.com/hendriau/Summix}{https://github.com/hendriau/Summix} for more information. Uses the same input data as summix.
#' @param reference a character vector of the column names for the reference groups.
#' @param observed a string that is the column name for the observed group.
#' @param pi.start Length K numeric vector of the starting guess for the reference group proportions. If not specified, this defaults to 1/K where K is the number of reference groups.
#' @return numeric value that is the scaled objective per 1000 SNPs
#' @importFrom magrittr "%>%"
#' @export

calc_scaledObj <- function(data, reference, observed, pi.start) {
  start_time <- Sys.time()
  # Will perform summix calculation on SNPs only within these bins (only need one side of the distribution since it is symmetrical)
  bins <- c("0-0.1", "0.1-0.3", "0.3-0.5")
  # these are the weights that the bins' objective will be multiplied by
  multiplier <- c(5, 1.5, 1)
  
  data$obs <- data[[observed]]
  
  # create bins and run summix on each bin
  data <- data %>% dplyr::rowwise() %>%
    dplyr::mutate(bin = dplyr::case_when(obs <= 0.1 ~ 1,
                                         obs >0.1 & obs <= 0.3 ~ 2,
                                         obs >0.3 & obs <= 0.5 ~ 3))
  for(b in 1:3) {
    # subset data to only SNPs in given bin
    subdata <- data[which(data$bin == b),]
    
    # ensure there are SNPs in the given bin before running summix
    if(nrow(subdata) > 0) {
      res <- summix_calc(subdata, reference = reference, observed = observed, pi.start)
      
      # isolate just the objective and scale per 1000 SNPs
      res$obj_adj <-  res$goodness.of.fit/(nrow(subdata)/1000)
      res$nSNPs <- nrow(subdata)
      res$bin <- bins[b]
      if(b == 1) {
        sum_res <- res
      } else {
        sum_res <- rbind(sum_res, res)
      }
      
      # if  there are no SNPs in the bin just run summix calculation on whole dataset,
      #and replace objective/1000SNPs with NA - necessary step for properly calculating the weighted average
    } else {
      res <- summix_calc(data, reference = reference,
                         observed = observed, pi.start)
      res$obj_adj <- NA
      res$nSNPs <- 0
      res$bin <- bins[b]
      if(b == 1) {
        sum_res <- res
      } else {
        sum_res <- rbind(sum_res, res)
      }
    }
  }
  
  # calculate the weighted average for non NA bins
  objective_scaled <- 0
  nonNA <- 0
  for(b in 1:3) {
    if(!is.na(sum_res[b, ]$obj_adj)) {
      #count number of non NA bins
      nonNA <- nonNA + 1
      objective_scaled <- objective_scaled +
        (sum_res[b, ]$obj_adj*multiplier[b])
    }
  }
  #divide weighted sum (objective_scaled) by number of non NA bins to get weighted average
  objective_scaled <- objective_scaled/nonNA
  end_time <- Sys.time()
  
  #return scaled objective value for the summix reference group proportion estimates
  return(objective_scaled)
}






