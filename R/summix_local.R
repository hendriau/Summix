#' summix_local
#'
#' @description
#' Estimates local ancestry mixture proportions in genetic summary data; Also performs a selection scan (optional) that identifies regions of selection along the given chromosome.
#'
#' @param data a dataframe of the observed and reference allele frequencies for N genetic variants of a single chromosome. Must have a column for the positions with column name "POS"
#' @param reference a character vector of the column names for the reference ancestries.
#' @param observed a character value that is the column name for the observed group.
#' @param goodness.of.fit an option to override the default scaled objective to return the raw loss from slsqp
#' @param type can be c("variants", "bp"), default "variants" to determine how the user wants to define window size (either in # variants or size in bp)
#' @param algorithm can be c("fastcatch", "windows") to use either the fastcatchup algorithm or the basic sliding windows. Usually, fastcatch is recommended for dynamic window sizes but is computationally slower
#' @param minVariants must be specified if algorithm == "fastcatch" and type == "variants" to define the smallest window size
#' @param maxVariants must be specified for either algorithm and type == "variants" to define the largest window size
#' @param minWindowSize must be specified if algorithm == "fastcatch" and type == "bp" to define the smallest window size in bp
#' @param maxWindowSize must be specified for either algorithm and if type == "bp" to define the largest window size in bp
#' @param maxStepSize default =1000bp to define the maximum gap in bp between two consecutive SNPs within a window
#' @param windowOverlap must be specified if algorithm == "windows", default = 200 to define the overlap between sliding windows
#' @param diffThreshold must be specified if algorithm == "fastcatch" default = 0.02 to define the % difference threshold to mark the end of a block
#' @param NSimRef a vector the same length as reference and in the same order with the number of individuals in each reference group - used in re-simulation for standard error estimation
#' @param override_fit default value is FALSE. An option for the user to override the autostop if the global objective is greater than 1.5 (poor fit)
#' @param override_removeSmallAnc An option for the user to override the automatic removal of reference ancestries with <2% global proportions - not recommended; default value is FALSE.
#' @param selection_scan An option to calculate local ancestry test statistic for the blocks; default value is FALSE.
#' @param position_col default value is 'POS'. Column of input data frame that contains the SNP position variable
#' @param data a data frame of the observed group and reference group allele frequencies for N genetic variants on a single chromosome. Must contain a column specifying the genetic variant positions.
#' @param reference a character vector of the column names for K reference groups.
#' @param observed a character value that is the column name for the observed group.
#' @param position_col a character value that is the column name for the genetic variants positions. Default is "POS".
#' @param maxStepSize a numeric value that defines the maximum gap in base pairs between two consecutive genetic variants within a given window. Default is 1000.
#' @param algorithm user choice of algorithm to define local ancestry blocks; options "fastcatch" and "windows" are available. "windows" uses a fixed window in a sliding windows algorithm. "fastcatch" allows dynamic window sizes. The "fastcatch" algorithm is recommended- though it is computationally slower. Default is "fastcatch".
#' @param type user choice of how to define window size; options "variants" and "bp" are available where "variants" defines window size as the number of variants in a given window and "bp" defines window size as the number of base pairs in a given window. Default is "variants".
#' @param override_fit default is FALSE. If set as TRUE, the user will override the auto-stop of summix_local() that occurs if the global goodness of fit value is greater than 1.5 (indicating a poor fit of the reference data to the observed data).
#' @param override_removeSmallAnc default is FALSE. If set as TRUE, the user will override the automatic removal of reference ancestries with <2% global proportions – this is not recommended.
#' @param selection_scan user option to perform a selection scan on the given chromosome. Default is FALSE. If set as TRUE, a test statistic will be calculated for each local ancestry block. Note: the user can expect extended computation time if this option is set as TRUE.
#' @param windowOverlap Used if algorithm = "windows". A numeric value that defines the number of variants or the number of base pairs that overlap between the given sliding windows. Default is 200.
#' @param diffThreshold Used if algorithm = "fastcatch". A numeric value that defines the percent difference threshold to mark the end of a local ancestry block. Default is 0.02.
#' @param maxVariants Used if type = "variants". A numeric value that specifies the maximum number of genetic variants allowed to define a given window.
#' @param maxWindowSize Used if type = "bp". A numeric value that defines the maximum allowed window size by the number of base pairs in a given window.
#' @param minVariants Used if algorithm = "fastcatch" and type = "variants". A numeric value that specifies the minimum number of genetic variants allowed to define a given window.
#' @param minWindowSize Used if algorithm = "fastcatch" and type = "bp". A numeric value that specifies the minimum number of base pairs allowed to define a given window.
#' @param  NSimRef Used if f selection_scan = TRUE. A numeric vector of the sample sizes for each of the K reference groups that is in the same order as the reference parameter. This is used in a simulation framework that calculates within local ancestry block standard error.
#'
#' @return data frame with a row for each local ancestry block and the following columns:
#' @return goodness.of.fit: scaled objective reflecting the fit of the reference data. Values between 0.5-1.5 are considered moderate fit and should be used with caution. Values greater than 1.5 indicate poor fit, and users should not perform further analyses using summix
#' @return iterations: number of iterations for SLSQP algorithm
#' @return time: time in seconds of SLSQP algorithm
#' @return filtered: number of SNPs not used in estimation due to missing values
#' @return K columns of mixture proportions of reference ancestry groups input into the function
#' @return nSNPs: number of SNPs in the given local ancestry block
#'
#' @author Hayley Wolff (Stoneman), \email{hayley.wolff@cuanschutz.edu}
#' @author Audrey Hendricks, \email{audrey.hendricks@cuanschutz.edu}
#'
#' @references https://github.com/hendriau/Summix2
#'
#' @keywords genetics, mixture distribution, admixture, population stratification, local ancestry
#'
#' @importFrom stats "sd" "pt"
#' @importFrom dplyr "arrange"
#'
#' @seealso \url{https://github.com/hendriau/Summix2} for further documentation.
#'
#'
#' @examples
#' data(ancestryData)
#' results <- summix_local(data = ancestryData,
#'                         reference = c("reference_AF_afr",
#'                                       "reference_AF_eas",
#'                                       "reference_AF_eur",
#'                                       "reference_AF_iam",
#'                                       "reference_AF_sas"),
#'                         NSimRef = c(704,787,741,47,545),
#'                         observed="gnomad_AF_afr",
#'                         goodness.of.fit = TRUE,
#'                         type = "variants",
#'                         algorithm = "fastcatch",
#'                         minVariants = 150,
#'                         maxVariants = 250,
#'                         maxStepSize = 1000,
#'                         diffThreshold = .02,
#'                         override_fit = FALSE,
#'                         override_removeSmallAnc = TRUE,
#'                         selection_scan = FALSE,
#'                         position_col = "POS")
#'print(results$results)
#'
#' @export

summix_local <- function(data, reference, observed, goodness.of.fit = TRUE,
                         type = "variants", algorithm = "fastcatch",
                         minVariants = 0,maxVariants = 0, maxWindowSize = 0,
                         minWindowSize = 0, windowOverlap = 200,
                         maxStepSize=1000, diffThreshold = .02, NSimRef=NULL, override_fit = FALSE,
                         override_removeSmallAnc = FALSE, selection_scan = FALSE,
                         position_col = "POS") {
  # valid input checking
  if(!(tolower(type) == "variants" | tolower(type) == "bp")) {
    stop(paste0("ERROR: type needs to be one of c('variants', 'bp'), user input: ", type))
  }
  if(!(algorithm == "fastcatch" | algorithm == "windows")) {
    stop(paste0("ERROR: algorithm needs to be one of c('fastcatch', 'windows'), user input: ",
                algorithm))
  }
  if(!is(object = data, class2 = "data.frame")){
    stop("ERROR: data must be a data.frame as described in the vignette")
  }
  if(typeof(observed)!="character"){
    stop("ERROR: 'observed' must be a character string for the column name of the observed ancestry in data")
  }
  if(!(observed %in% names(data))){
    stop("ERROR: 'observed' must be the column name of the observed ancestry in data")
  }
  if(typeof(reference)!="character"){
    stop("ERROR: 'reference' must be a vector of column names from data to be used in the reference")
  }
  if(all(reference %in% names(data))==FALSE){
    stop("ERROR: 'reference' must be a vector of column names from data to be used in the reference")
  }
  if(type == "variants" & algorithm == "fastcatch" & (maxVariants < minVariants)) {
    stop("ERROR: maxVariants must be larger than minVariants")
  }
  if(type == "bp" & algorithm == "fastcatch" & (maxWindowSize < minWindowSize)) {
    stop("ERROR: maxWindowSize must be larger than minWindowSize")
  }
  if(diffThreshold < 0 | diffThreshold > 1) {
    stop("ERROR: diffThreshold must be between 0 and 1")
  }
  if((selection_scan) & (length(NSimRef) != length(reference))) {
    stop("ERROR: NSimRef and Reference vectors must be the same length")
  }
  # first run summix global to determine if poor fit and remove any ancestries <2%
  globTest <- summix_quiet(data = data,
                           reference = reference,
                           observed = observed)
  if(globTest$goodness.of.fit > 1.5 & override_fit == F) {
    return(paste0("Reference not sufficiently well matched to observed sample: objective > 1.5 (",
                  globTest$goodness.of.fit, ") [to override stop use override = TRUE]"))
  }
  
  # remove any references with prop < 2%
  if(!override_removeSmallAnc) {
    if(sum(globTest[1,5:ncol(globTest)] < 0.02) > 0) {
      toremove <- which(globTest[1,5:ncol(globTest)] < 0.02)
      reference <- reference[-toremove]
    }
  }
  
  # initialize place to store results
  results <- data.frame()
  
  #declare vars to be used in function
  POS <- NULL
  Start_Pos <- NULL
  End_Pos <- NULL
  
  
  # set POS column in dataframe if not already there
  if(position_col != "POS") {
    data$POS <- data[,position_col]
  }
  
  POS <- data$POS
  
  # first arrange all data by position
  chrData <- data %>% arrange(POS)
  
  # local variable to store start time
  start_time <- Sys.time()
  windowSize <- 0
  
  if(algorithm == "fastcatch") {
    # algorithm for if using min/max variant window
    if (tolower(type == "variants")) {
      # absolute minimum number of variants in a window 150
      if (minVariants < 150) {
        minVariants <- 150
      }
      increment <- 1
      if(maxVariants > 10000) {increment <- floor(maxVariants/2000)}
      startPoint <- 1
      endPoint <- variantGetNext(chrData$POS, startPoint, minVariants)
      lastProportions <- data.frame(AFR = double(),
                                    EAS = double(),
                                    EUR = double(),
                                    IAM = double(),
                                    SAS = double())
      currentProportions <- data.frame(AFR = double(),
                                       EAS = double(),
                                       EUR = double(),
                                       IAM = double(),
                                       SAS = double())
      
      #record start time
      start_time <- Sys.time()
      iteration = 0
      
      # begin fast/catchup algorithm
      while(endPoint <= nrow(chrData)) {
        iteration = iteration + 1
        
        subData <- chrData[startPoint:endPoint, ]
        lastProportions <- summix_quiet(subData,
                                        reference = reference,
                                        observed = observed,
                                        goodness.of.fit = goodness.of.fit)
        
        #increment window to next point
        
        #first check if end point is at last position
        #if it is, save this block and end
        if(endPoint == nrow(chrData) | startPoint == nrow(chrData)) {
          
          results <- saveBlock(chrData, startPoint,
                               endPoint, lastProportions,
                               results)
          break
        }
        
        #check if incrementing to next step will violate maximum variants or maximum step size between variants
        if(((endPoint + 1) - startPoint - 1) > maxVariants |
           (chrData$POS[endPoint + 1] - chrData$POS[endPoint]) > maxStepSize)
        {
          #if new step violates set limits then save block and set new start and end points
          results <- saveBlock(chrData, startPoint,
                               endPoint, lastProportions,
                               results)
          startPoint <- endPoint
          endPoint <- variantGetNext(chrData$POS, startPoint, minVariants)
        }
        else {
          endPoint <- endPoint + increment
          subData <- chrData[startPoint:endPoint, ]
          currentProportions <- summix_quiet(subData,
                                             reference = reference,
                                             observed = observed,
                                             goodness.of.fit = goodness.of.fit)
          
          #determine if new ancestry proportions are different than last
          areDiff <- testDiff(lastProportions, currentProportions,
                              threshold = diffThreshold)
          if(areDiff) { # if different save new block with point prior to current end point
            
            results <- saveBlock(chrData, startPoint,
                                 endPoint - 1, lastProportions,
                                 results)
            
            startPoint <- endPoint
            endPoint <- variantGetNext(chrData$POS, startPoint, minVariants)
          }
        }
      }
    }
    
    # else algorithm for if using min/max size of window
    else if (tolower(type == "bp")) {
      startPoint <- 1
      
      endPoint <- sizeGetNext(chrData$POS, startPoint, minWindowSize)
      if (endPoint < 150) {
        endPoint <- 150
      }
      lastProportions <- data.frame(AFR = double(),
                                    EAS = double(),
                                    EUR = double(),
                                    IAM = double(),
                                    SAS = double())
      currentProportions <- data.frame(AFR = double(),
                                       EAS = double(),
                                       EUR = double(),
                                       IAM = double(),
                                       SAS = double())
      #Record start time
      start_time <- Sys.time()
      
      iteration = 0
      
      #begin fast/catchup algorithm using size instead of variants
      while(endPoint <= nrow(chrData)) {
        iteration = iteration + 1
        
        subData <- chrData[startPoint:endPoint, ]
        lastProportions <- summix_quiet(subData,
                                        reference = reference,
                                        observed = observed,
                                        goodness.of.fit = goodness.of.fit)
        
        #increment window to next point
        
        #first check if end point is at last position
        #if it is, save this block and end
        if(endPoint == nrow(chrData) | startPoint == nrow(chrData)) {
          
          results <- saveBlock(chrData, startPoint,
                               endPoint, lastProportions, results)
          break
        }
        
        #check if incrementing to next step will violate maximum variants or maximum step size between variants
        else if((chrData$POS[endPoint + 1] - chrData$POS[startPoint]) > maxWindowSize |
                (chrData$POS[endPoint + 1] - chrData$POS[endPoint]) > maxStepSize) {
          
          #if new step violates set limits then save block and set new start and end points
          results <- saveBlock(chrData, startPoint,
                               endPoint, lastProportions, results)
          startPoint <- endPoint
          endPoint <- sizeGetNext(chrData$POS, startPoint, minWindowSize)
          
        }
        
        # check if step size is too small (minimum windowsize is also 150 variants)
        else if(endPoint - startPoint < 149) {
          endPoint <- startPoint + 149
        }
        
        else {
          
          endPoint <- endPoint + 1
          subData <- chrData[startPoint:endPoint, ]
          currentProportions <- summix(subData,
                                       reference = reference,
                                       observed = observed,
                                       goodness.of.fit = goodness.of.fit)
          
          #determine if new ancestry proportions are different than last
          areDiff <- testDiff(lastProportions, currentProportions,
                              threshold = diffThreshold)
          if(areDiff) { # if different save new block with point prior to current end point
            
            results <- saveBlock(chrData, startPoint,
                                 endPoint - 1, lastProportions,
                                 results)
            startPoint <- endPoint
            endPoint <- sizeGetNext(chrData$POS, startPoint, minWindowSize)
          }
        }
      }
    }
  } else if (algorithm == "windows") {
    # initialize variables needed throughout algorithm
    startPoint <- 1
    endPoint <- 1
    start_time <- Sys.time()
    windowSize = max(maxWindowSize, maxVariants)
    
    if(type == "variants") {
      endPoint <- startPoint + windowSize - 1
    } else {
      endPoint <- getNextEndPoint(chrData, startPoint, windowSize)
    }
    
    proportions <- data.frame()
    iteration <- 0
    # begin sliding window algorithm
    while(endPoint <= nrow(chrData)) {
      iteration <- iteration + 1
      subData <- chrData[startPoint:endPoint, ]
      
      #use summix to estimate proportions
      proportions <- summix_quiet(subData, reference = reference, observed = observed,
                                  goodness.of.fit = goodness.of.fit)
      
      #save this block to results
      results <- saveBlock(chrData, startPoint, endPoint,
                           proportions, results)
      
      #increment window forward
      # first check if end point is already at end of chromosome; if yes end
      if (endPoint == nrow(chrData)) {
        break
      } else {
        if (type == "variants") {
          startPoint <- endPoint - (windowOverlap)
          endPoint <- startPoint + windowSize - 1
        } else {
          startPoint <- getNextStartPoint(chrData, startPoint, endPoint, windowOverlap)
          endPoint <- getNextEndPoint(chrData, startPoint, windowSize)
        }
      }
    }
  }
  
  print("Done getting LA proportions")
  
  if(selection_scan) {
    # calculate t and p-values if running selection scan
    sd <- (apply(results[,reference], MARGIN = 2, FUN = sd))
    # get internal simulation SE
    
    print("Running internal simulations for SE")
    
    se.2 <- doInternalSimulation(windows = results %>%
                                   select(Start_Pos, End_Pos),
                                 data = data,
                                 reference = reference,
                                 observed = observed,
                                 nRefs = NSimRef)
    
    # calculate test statistic by weighted average of SD and simulated sE
    statVals <- data.frame(matrix(ncol = length(reference),
                                  nrow = nrow(results)))
    colnames(statVals) <- paste0("t.", reference, ".avg")
    
    for(i in 1:nrow(results)) {
      for(j in 1:length(reference)) {
        statVals[i,j] <- (results[i, reference[j]] -
                            mean(results[-i,reference[j]]))/
          ((se.2$simSE[i,j] + sd[j])/2)
      }
    }
    results <- cbind(results, statVals)
    
    pvals <- data.frame(matrix(ncol = length(reference),
                               nrow = nrow(results)))
    colnames(pvals) <- paste0("p.", reference)
    for(i in 1:length(reference)){
      pvals[ , i] <- 2*pt(statVals[ , i], df = results$nSNPs, lower.tail = F)
    }
    results <- cbind(results, pvals)
    
    end_time <- Sys.time()
    print(difftime(end_time, start_time, units = "auto"))
    print(paste0("Discovered ", nrow(results), " LA blocks"))
    
    toReturn <- list(results = results,
                     sd = sd,
                     se_sim = se.2)
    return(toReturn)
  }
  return(list(results = results))
}




#' sizeGetNext
#'
#' @description
#' Helper function to get starting end point that is a minimum distance (in bases) from start point; uses indices NOT position numbers
#' @param positions list of positions of variants
#' @param start index of the current start position
#' @param minSize integer defining the minimum size in bp of the window
#' @return the new end point index
#'
#' @export

sizeGetNext <- function (positions, start, minSize) {
  size <- 0
  index <- 1
  
  for (i in 1:(length(positions) - start)) {
    
    size <- positions[start + i] - positions[start]
    
    if (size >= minSize) {
      index <- i
      return(start + index)
    } else {
      #continue
    }
  }
  return(length(positions))
}



#' variantGetNext
#'
#' @description
#' Helper function to get starting end point that is a minimum distance (in variants) from start point; uses indices NOT position numbers
#' @param positions list of positions of variants
#' @param start index of the current start position
#' @param minVariants integer defining the minimum size in number of variants of the window
#' @return the new end point index
#'
#' @export

variantGetNext <- function(positions, start, minVariants) {
  if((start + minVariants - 1) < length(positions)) {
    return(start + minVariants - 1)
  } else {
    return(length(positions))
  }
}



#' saveBlock
#'
#' @description
#' Helper function to save one block to results
#' @param data the input dataframe subsetting to just the chromosome
#' @param start index of start of block
#' @param end index of the end of block
#' @param props ancestry proportions for the block returned from summix
#' @param results current results dataframe
#'
#' @export

saveBlock <- function(data, start, end, props, results) {
  
  currResults <- c(Start_Pos = data$POS[start],
                   End_Pos = data$POS[end], props,
                   nSNPs = end-start)
  
  
  newResults <- rbind(results, currResults)
  return(newResults)
}



#' testDiff
#'
#' @description
#' Helper function to determine whether ancestry  has changed for fast/catchup window algorithm
#' @param last ancestry proportions of block returned from summix
#' @param current ancestry proportions of block returned from summix
#' @param threshold if applicable the threshold for determing change point
#' @return true if passes threshold, false if not
#'
#' @export

testDiff <- function(last, current, threshold = .01) {
  above = rep(FALSE, (length(last) - 4))
  
  for (i in 1:(length(above))) {
    
    above[i] = ifelse((abs(current[i + 4] - last [i + 4])) > threshold, TRUE, FALSE)
  }
  if (TRUE %in% above) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}



#' getNextEndPoint
#'
#' @description
#' Helper function: algorithm to get next end point in basic window algorithm; will find first point that is at least window size away from start
#' @param data the input dataframe subset to the chromosome
#' @param start index of the current start point
#' @param windowSize the window size (in bp or variants)
#' @return index of end point of window
#'
#' @export

getNextEndPoint <- function(data, start, windowSize) {
  if (start == nrow(data)) {
    return(start)
  }
  for (i in (start + 1):nrow(data)) {
    if ((data$POS[i] - data$POS[start]) >= windowSize) {
      return(i)
    }
  }
  return(nrow(data))
}


#' getNextStartPoint
#'
#' @description
#' Helper function: algorithm to get next start point; will pick the point that provides approx. the specified amount of overlap, but not more; if there are only two variants in the previous block, will jump new start point to the previous end point
#' @param data the input dataframe subset to the chromosome
#' @param start the current index of start point
#' @param end the current index of end point
#' @param overlap the desired amount of window overlap (in bp or variants)
#' @return returns index of new start point
#'
#' @export

getNextStartPoint <- function(data, start, end, overlap) {
  #if only 2 variants in window next start point needs to be current end point
  if(end - start == 1) {
    return(end)
  }
  for (i in end - 1:1) {
    if((data$POS[end] - data$POS[i]) <= overlap) {
      return(i)
    }
  }
  return(end - 1)
}

#' doInternalSimulation
#'
#' @description
#' Helper function to get the within block se using re-simulation
#' @param windows is a dataframe with the Start_Pos and End_Pos
#' @param data is the original chromosome data
#' @param reference is a list with the names of the columns with references
#' @param observed a character value that is the column name for the observed group
#' @param nRefs is a vector the same lengths as reference with the number of individuals in each reference population
#'
#' @importFrom stats "sd" "rmultinom"
#'
#' @export

doInternalSimulation <- function(windows, data, reference, observed, nRefs) {
  
  allSE <- matrix(nrow = nrow(windows), ncol = length(reference))
  colnames(allSE) <- reference
  simRefs <- paste0("S_", reference)
  allProps <- vector(mode = "list", length = nrow(windows))
  
  # for each window in the results
  for(i in 1:nrow(windows)) {
    props <- data.frame(matrix(nrow = 1000, ncol = length(reference)))
    for(boot in 1:1000) {
      # first subset data to just SNPs in window
      subdata <- data[(data$POS>=windows[i,]$Start_Pos & data$POS<=windows[i,]$End_Pos),]
      
      # simulate new reference
      for(r in 1:length(reference)) {
        newVals <- rep(0, nrow(subdata))
        vals <- subdata[[reference[r]]]
        for(s in 1:nrow(subdata)) {
          sim <- rmultinom(1, nRefs[r], c(vals[s]**2,
                                          2*vals[s]*(1-vals[s]),
                                          (1-vals[s])**2))
          newVals[s] <- ((sim[1]*2) + (sim)[2])/(2*nRefs[r])
        }
        
        subdata <- cbind.data.frame(subdata, newVals)
        names(subdata)[ncol(subdata)] <- paste0("S_", reference[r])
      }
      
      # re-estimate using summix
      sum_res <- invisible(summix_quiet(data = subdata, reference = simRefs, observed = observed))
      props[boot,] <-  sum_res[,simRefs]
    }
    allProps[[i]] <- props
    allSE[i, ] <- apply(props, 2, function(x) sd(x))
    
  }
  toReturn <- list(allProps = allProps,
                   simSE = allSE)
  return(toReturn)
}


# summix function that does not have any print outputs for local ancestry use
summix_quiet <- function(data, reference, observed, pi.start = NA, goodness.of.fit = TRUE) {
  start_time = Sys.time()
  # get the ancestry proportions using the summix_calc helper function either with or without pi.start specified
  if(!is.na(pi.start)) {
    sum_res <- summix_calc(data = data, reference = reference,
                           observed = observed, pi.start = pi.start)
  } else {
    sum_res <- summix_calc(data = data, reference = reference, observed = observed)
  }
  
  # get the scaled objective and replace the objective argument (first column) with the updated objective value
  # will not run if overridden by user, but is the default
  if(goodness.of.fit) {
    new_obj <- calc_scaledObj(data = data, observed = observed, reference = reference, pi.start= pi.start)
    sum_res[1] <- new_obj
  }
  end_time = Sys.time()
  ttime = end_time - start_time
  sum_res[3] <- ttime
  
  return(sum_res)
}
