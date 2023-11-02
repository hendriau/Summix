#' adjAF
#'
#' @description
#' Adjusts allele frequencies for heterogeneous populations in genetic data given proportion of reference groups
#'
#' @param data dataframe of unadjusted allele frequency for observed group, K reference group allele frequencies for N SNPs
#' @param reference character vector of the column names for K reference groups.
#' @param observed character value for the column name of observed data group
#' @param pi.target numeric vector of the mixture proportions for K reference groups in the target individual or group.
#' @param pi.observed  numeric vector of the mixture proportions for K reference groups in the observed group.
#' @param adj_method user choice of method for the allele frequency adjustment
#' @param N_reference numeric vector of the sample sizes for each of the K reference groups.
#' @param N_observed numeric value of the sample size of the observed group.
#' @param filter sets adjusted allele frequencies equal to 1 if > 1, to 0 if > -.005 and < 0, and removes adjusted allele frequencies < -.005.
#'
#'
#' @return pi: table of input reference groups, pi.observed, and pi.target
#' @return observed.data: name of the data column for the observed group from which adjusted allele frequency is estimated
#' @return Nsnps: number of SNPs for which adjusted AF is estimated
#' @return adjusted.AF: data frame of original data with an appended column of adjusted allele frequencies
#' @return effective.sample.size: The sample size of individuals effectively represented by the adjusted allele frequencies
#' @importFrom stats "na.omit"
#'
#'
#' @author Adelle Price, \email{adelle.price@cuanschutz.edu}
#' @author Hayley Wolff, \email{hayley.wolff@cuanschutz.edu}
#' @author Audrey Hendricks, \email{audrey.hendricks@cuanschutz.edu}
#' @references https://github.com/hendriau/Summix2
#' @keywords genetics, mixture distribution, admixture, population stratification
#'
#' @seealso \url{https://github.com/hendriau/Summix2} for further documentation.
#'
#'
#' @examples
#' data(ancestryData)
#' adjusted_data<-adjAF(data   = ancestryData,
#'     reference  = c("reference_AF_afr", "reference_AF_eur"),
#'     observed    = "gnomad_AF_afr",
#'     pi.target   = c(1, 0),
#'     pi.observed = c(.85, .15),
#'     adj_method = 'average',
#'     N_reference = c(704,741),
#'     N_observed = 20744,
#'     filter = TRUE)
#' adjusted_data$adjusted.AF[1:5,]
#'
#'
#'
#' @export

adjAF <- function(data,
                  reference,
                  observed,
                  pi.target,
                  pi.observed,
                  adj_method = 'average',
                  N_reference = NULL,
                  N_observed = NULL,
                  filter = TRUE) {
  
  #check for correct format of user inputs
  if(length(reference) != length(pi.target)){
    stop("ERROR: Please make sure that you are entering k groups in reference.")
  }
  if(length(N_reference) != length(reference)){
    stop("ERROR: Please make sure that the lengths of N_reference and reference match.")
  }
  if(length(pi.target) != length(pi.observed)){
    stop("ERROR: Please make sure that the lengths of pi.target and pi.observed match.")
  }
  if(!is(object = data, class2 = "data.frame")){
    stop("ERROR: data must be a data.frame as described in the vignette")
  }
  if(typeof(observed)!="character"){
    stop("ERROR: please enter the column name of the observed group")
  }
  if( !( observed %in% names(data) ) ){
    stop("ERROR: please make sure that the observed group is included in the column names of the data.")
  }
  if( !all(reference %in% names(data) ) ){
    stop("ERROR: Please make sure that all groups in reference are also column names in data")
  }
  if(typeof(reference)!="character"){
    stop("ERROR: please enter the column names for the reference groups")
  }
  if( !(all(pi.observed<=1) & all(pi.target<=1)) | !(all(pi.observed>=0) & all(pi.target>=0)) ){
    stop("ERROR: pi.observed and pi.target must contain ratios between 0 and 1.")
  }
  
  
  
  
  #Calculate average fold change of pi.observed to pi.target for across all reference groups
  ave_fold_change = mean(abs(((pi.target+.001)-(pi.observed+.001))/(pi.observed+.001)))
  print(paste0('Average fold change between observed and target group proportions is: ', round(ave_fold_change, 2)))
  #if average fold change is greater than 3 and less than or equal to 5, print warning
  if (ave_fold_change > 3 & ave_fold_change <= 5){
    print("Warning: The average fold change between the observed and target group proportions is greater than 3.")
  }
  #if average fold change is greater than 5, print warning. AF adjustment with average fold change of pi.observed to pi.target greater than 5 is not recommended
  if (ave_fold_change > 5){
    warning("The average fold change between the observed and target group proportions is greater than 5. It is not recommended to perform this AF adjustment.")
  }
  
  fold_change = (floor(sum(pmin(pi.observed, pi.target)*N_observed + N_reference*min(1,pi.target*length(pi.target))))*pi.target - (N_observed*pi.observed+N_reference))/(N_observed*pi.observed+N_reference)
  fold_change_refs = c()
  
  large_fold_change=FALSE
  #Iteratively check every pi.observed to pi.target fold change. If a given fold change is larger than 10, print warning showing which reference groups have too large of fold changes
  for (i in 1:length(N_reference)){
    if (fold_change[i] >= 10){
      fold_change_refs = append(fold_change_refs, reference[i])
      large_fold_change=TRUE
    }
  }
  if (large_fold_change==TRUE){
    warning(paste0('Based on the observed group and reference group sample sizes, the observed to target group proportion increase for the ', list(fold_change_refs),' group is too large. ', 'This group adjustment is not recommended.'))
  }
  
  
  #check the effective sample size of the adjusted AF, if less than 50% of observed sample size, print warning
  eff_samp_size = floor(sum(pmin(pi.observed, pi.target)*N_observed + N_reference*min(1,(pi.target*length(pi.target)))))
  cat('\n')
  if (eff_samp_size < .5*N_observed){
    warning('Your effective sample size of the adjusted alelle frequency is less than 50% of the observed sample size.')
  }
  
  
  
  allResults <- data.frame(matrix(nrow = nrow(data), ncol = length(reference)))
  #If user has selected 'average' method for AF adjustment (this is default)
  if (adj_method == "average"){
    for(i in 1:length(reference)) {
      #Iteratively leave each reference group out and calculate adjusted AF
      temp <- pi.target[i]
      s.targ <- pi.target[-i]
      s.targ <- c(s.targ, temp)
      
      temp <- pi.observed[i]
      s.obs <- pi.observed[-i]
      s.obs <- c(s.obs, temp)
      
      res <- invisible(adjAF_calc(data = data,
                                  reference = reference[-i],
                                  observed = observed,
                                  pi.target = s.targ,
                                  pi.observed = s.obs))
      
      allResults[,i] <- res$adjusted.AF$adjustedAF
    }
    
    #Average all adjusted AFs that leave one reference group out at a time to get final adjusted AF
    allResults$mean <- apply(allResults, 1, mean)
    
    #if user has set filter == TRUE (this is default), remove negative values less than -.005 from final adjusted AFs, set adjusted AFs less than 0 and greater than -.005 to 0, and set adjusted AFs greater than 1 to 1
    if (filter == TRUE){
      before_filt_negative = length(allResults$mean)
      before_filt_less = length(allResults$mean[allResults$mean < 0 & allResults$mean > -0.005])
      before_filt_greater = length(allResults$mean[allResults$mean > 1])
      allResults$mean[allResults$mean > 1] <- 1
      allResults$mean[allResults$mean < 0 & allResults$mean > -0.005] <- 0
      allResults$mean[allResults$mean <= -0.005] <- NA
      
      #Let user know how many adjusted AFs were set to 0, 1, or removed from the final adjusted AF dataframe
      cat('\n')
      print(paste0("Note: In this AF adjustment, ", before_filt_less, " SNPs (with adjusted AF > -.005 & < 0) were rounded to 0. ", before_filt_greater, " SNPs (with adjusted AF > 1) were rounded to 1, and ", length(allResults$mean[allResults$mean <= -0.005]), " SNPs (with adjusted AF <= -.005) were removed from the final results."))
      cat('\n')
    }
    
    
    res$adjusted.AF$adjustedAF <- allResults$mean
    #remove any NA's in dataframe (there will only be NA's if user has set filter ==TRUE and the adjusted AFs less than -.005 were set as NA in the previous code chunk)
    res$adjusted.AF <- na.omit(res$adjusted.AF)
    res$Nsnps <- length(res$adjusted.AF$adjustedAF)
    res$effective.sample.size <- eff_samp_size
    
  }
  
  
  #If user has selected 'effective' method for AF adjustment
  if (adj_method == "effective"){
    newResults <- data.frame(matrix(nrow = nrow(data), ncol = length(reference)))
    N_effective= vector()
    
    effective_new_prop <- vector()
    for(i in 1:length(reference)) {
      #Iteratively leave each reference group out and calculate adjusted AF
      temp <- pi.target[i]
      s.targ <- pi.target[-i]
      s.targ <- c(s.targ, temp)
      
      temp <- pi.observed[i]
      s.obs <- pi.observed[-i]
      s.obs <- c(s.obs, temp)
      
      res <- invisible(adjAF_calc(data = data,
                                  reference = reference[-i],
                                  observed = observed,
                                  pi.target = s.targ,
                                  pi.observed = s.obs))
      
      allResults[,i] <- res$adjusted.AF$adjustedAF
      #calculate N_effective for this iteration
      N_effective[i] <- calc_effective_N(N_reference[-i],N_observed,s.targ,s.obs)
      
    }
    
    #Create final adjusted AF as weighted sum
    #Weight all adjusted AFs that leave one reference group out at a time by N_effective proportion of that reference group
    for(m in 1:length(reference)){
      newResults[,m] <- allResults[,m]* (N_effective[m]/sum(N_effective))
    }
    #sum all N_effective-weighted adjusted AFs to get final adjusted AF
    newResults$effective_summ<- apply(newResults, 1, sum)
    
    #if user has set filter == TRUE (this is default), remove negative values less than -.005 from final adjusted AFs, set adjusted AFs less than 0 and greater than -.005 to 0, and set adjusted AFs greater than 1 to 1
    if (filter == TRUE){
      before_filt_negative = length(newResults$effective_summ)
      before_filt_less = length(newResults$effective_summ[newResults$effective_summ < 0 & newResults$effective_summ > -0.005])
      before_filt_greater = length(newResults$effective_summ[newResults$effective_summ > 1])
      newResults$effective_summ[newResults$effective_summ > 1] <- 1
      newResults$effective_summ[newResults$effective_summ < 0 & newResults$effective_summ > -0.005] <- 0
      newResults$effective_summ[newResults$effective_summ <= -0.005] <- NA
      
      #Let user know how many adjusted AFs were set to 0, 1, or removed from the final adjusted AF dataframe
      cat('\n')
      print(paste0("Note: In this AF adjustment, ", before_filt_less, " SNPs (with adjusted AF > -.005 & < 0) were rounded to 0. ", before_filt_greater, " SNPs (with adjusted AF > 1) were rounded to 1, and ", length(newResults$effective_summ[newResults$effective_summ <= -0.005]), " SNPs (with adjusted AF <= -.005) were removed from the final results."))
      cat('\n')
    }
    
    res$adjusted.AF$adjustedAF <- newResults$effective_summ
    #remove any NA's in dataframe (there will only be NA's if user has set filter ==TRUE and the adjusted AFs less than -.005 were set as NA in the previous code chunk)
    res$adjusted.AF <- na.omit(res$adjusted.AF)
    res$Nsnps <- length(res$adjusted.AF$adjustedAF)
    res$effective.sample.size <- eff_samp_size
  }
  
  print(noquote('$pi'))
  print(cbind.data.frame('ref.group' = reference, 'pi.observed' = pi.observed, 'pi.target' = pi.target))
  cat('\n')
  print(noquote('$observed.data'))
  print(paste0("observed AF data to update: ", "'", observed, "'"))
  cat('\n')
  print(noquote('$Nsnps'))
  print(length(res$adjusted.AF$adjustedAF))
  cat('\n')
  cat('\n')
  print(noquote('$effective.sample.size'))
  print(res$effective.sample.size)
  cat('\n')
  cat('\n')
  print("use $adjusted.AF$adjustedAF to see adjusted AF data")
  cat('\n')
  cat('\n')
  print("Note: The accuracy of the AF adjustment is likely lower for rare variants (< .5%).")
  
  
  return(res)
  
}






#' adjAF_calc
#'
#' @description
#' Helper function for calculating allele frequencies for heterogeneous populations in genetic data given proportion of reference groups
#'
#' @param data dataframe of unadjusted allele frequency for observed group, K-1 reference group allele frequencies for N SNPs
#' @param reference character vector of the column names for K-1 reference groups. The name of the last reference group is not included as that group is not used to estimate the adjusted allele frequencies.
#' @param observed character value for the column name of observed data group
#' @param pi.target numeric vector of the mixture proportions for K reference groups in the target sample or subject. The order must match the order of the reference columns with the last entry matching the missing reference group.
#' @param pi.observed  numeric vector of the mixture proportions for K reference groups for the observed group. The order must match the order of the reference columns with the last entry matching the missing reference group.
#'
#' @return pi: table of input reference groups, pi.observed, and pi.target
#' @return observed.data: name of the data column for the observed group from which adjusted allele frequency is estimated
#' @return Nsnps: number of SNPs for which adjusted AF is estimated
#' @return adjusted.AF: data frame of original data with an appended column of adjusted allele frequencies
#' @importFrom tidyselect "all_of"
#' @importFrom dplyr "select"
#'
#' @export

#NOTE: this is the summix v1 adjAF code, now used as a helper function
adjAF_calc <- function(
    data         ,
    reference    ,
    observed     ,
    pi.target    ,
    pi.observed  ){
  
  k = length(reference) + 1
  
  #Normalize pi. We need the pi_reference and pi.target to sum to 1
  pi.target <- pi.target/sum(pi.target)
  pi.observed <- pi.observed/sum(pi.observed)
  
  #sum the reference group allele frequencies multiplied by pi.target and name it 'hatted'
  hatted <- vector(mode = "double", length = dim(data)[1])
  #sum up the reference group allele frequencies multiplied by pi_reference and name it 'starred'
  starred <- vector(mode = "double", length = dim(data)[1])
  
  #sum the K-1 reference groups multiplied by pi.observed.
  #sum the k-1 reference groups multiplied by pi.target.
  for (j in 1:(k - 1)) {
    hatted <- hatted + (pi.observed[j] * data %>% select(reference[j]))
    starred <- starred + (pi.target[j] * data %>% select(reference[j]))
  }
  
  
  #check pi.observed for the kth reference group is greater than 0 before dividing pi.target by pi.observed
  if (pi.observed[k]>0){
    data$adjustedAF <- ((pi.target[k]/pi.observed[k]) * (data %>% select(all_of(observed)) - hatted) + starred)
    #if pi.observed is not greater than 0, set pi.target/pi.observed equal to 0
  }else{
    data$adjustedAF <- ((data %>% select(all_of(observed)) - hatted) + starred)
  }
  
  data.out <- data
  observed.data=paste("observed data to update AF: '", observed, "'", sep="")
  pi_table<-data.frame(ref.group=c(reference, "NONE"), pi.observed=pi.observed, pi.target=pi.target)
  Nsnps=nrow(data)
  tmp.out<-list("pi"=pi_table, "observed.data"=observed.data, "Nsnps"=Nsnps, "adjusted.AF"=data.out, "effective.sample.size"="0")
  
  
  return(tmp.out)
  
}






#' calc_effective_N
#'
#' @description
#' Helper function to calculate effective sample size for the group that is left out when estimating the adjusted allele frequencies in each adjAF function iteration.
#
#' @param N_reference numeric vector of the sample sizes of each K reference groups.
#' @param N_observed numeric value of the sample size of the observed group.
#' @param pi.target numeric vector of the mixture proportions for K reference groups in the target sample or subject. The order must match the order of the reference columns with the last entry matching the missing reference group.
#' @param pi.observed  numeric vector of the mixture proportions for K reference groups for the observed group. The order must match the order of the reference columns with the last entry matching the missing reference group.
#'
#' @return N_effective: effective sample size for the group that is left out when estimating the adjusted allele frequencies in each adjAF function iteration.
#'
#' @export

calc_effective_N <- function (N_reference, N_observed, pi.target, pi.observed){
  #one reference group is left out for each calculation of N_effective, must ensure k is equal to the total number of reference groups in the AF adjustment
  k = length(N_reference) + 1
  #ensure pi.target and pi.observed are normalized
  pi.target <- pi.target/sum(pi.target)
  pi.observed <- pi.observed/sum(pi.observed)
  
  #initiate and calculate the second half of the AF adjustment equation across all SNPs
  hatted <- vector(mode = "double", length = 1)
  starred <- vector(mode = "double", length = 1)
  for (j in 1:(k - 1)) {
    hatted <- hatted + (pi.observed[j] * N_reference[j])
    starred <- starred + (pi.target[j] * N_reference[j])
    
  }
  
  #calculate first half of AF adjustment equation across all SNPs
  #first check if pi.observed[k] is not equal to 0 (if pi.observed[k] is equal to 0, (pi.target[k])/(pi.observed[k]) will be undefined)
  if (pi.observed[k] == 0){
    N_effective = starred
  } else{
    N_effective <- (((pi.target[k])/(pi.observed[k])) * (N_observed - hatted) + starred)
  }
  
  #return effective sample size of this iteration of the AF adjustment leaving one reference group out
  return(N_effective)
}

