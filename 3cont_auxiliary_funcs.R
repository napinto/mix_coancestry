# Generate allele profiles by random sampling based on cumulative allele frequencies
# alleles: vector of allele identifiers
# freqs_cum: vector of cumulative allele frequencies summing to 1
# nsims: number of profiles to simulate
generate_profiles <- function(alleles,freqs_cum,nsims){
  rand_numbers <- matrix(runif(nsims * 6),ncol = 6,byrow = T)
  
  indices <- findInterval(rand_numbers, freqs_cum,left.open = T,rightmost.closed = T) +1
  profiles <- matrix(alleles[indices], nrow = nrow(rand_numbers), ncol = ncol(rand_numbers))
  return(profiles)
}

# Count how many distinct alleles are present in each simulated profile
count_number_unique_alleles <- function(profiles){
  lvl <- sort(unique(c(profiles))) # list of alleles
  if (dim(profiles)[1] == 1){
    n_distinct_alleles <- length(c(unique(profiles)))
  } else {
    unique_counts <- vapply(lvl, function(x) rowSums(profiles == x, na.rm = TRUE), numeric(nrow(profiles))) # a matrix of size nsims by n_alleles counting how many times a given allele appears in a given row
    unique_counts_logical <- unique_counts > 0 # converts the counts into TRUE/FALSE
    n_distinct_alleles <- rowSums(unique_counts_logical)
  }
  return(n_distinct_alleles)
}

# Calculates likelihood ratio components for one genetic marker based on given allele profiles and frequencies, applying a theta correction.
# theta: population substructure correction parameter.
# profiles: matrix of allele profiles.
# freqs: allele frequencies.
# Returns calculated likelihood ratios or related statistic.
get_result_one_marker <- function(theta, profiles, freqs) {
  if(nrow(profiles)==1) {
    profiles <- rbind(profiles,profiles)
    remove_last_flag <- TRUE
  } else {
    remove_last_flag <- FALSE
  }
  
  K_theta_2 = (1+5*theta) * (1 + 6*theta)
  
  # Auxiliary function used in likelihood ratio calculations depending on parameters 
  F_theta <- function(p,c,theta){
    return((1-theta)*p + c*theta)
  }
  
  # func_1_1, func_1_2, ..., func_2_6(mixture, Ktheta2, freqs)
  # 
  # Different functions calculating likelihood ratio numerators and denominators based on the number of contributors and allele configurations in the mixture.
  # mixture: genotype mixture vector for contributors.
  # Ktheta2: precomputed constant related to theta-correction.
  # freqs: allele frequencies.
  # Each function handles specific cases in allele mixture composition for likelihood computations.
  func_1_1 <- function(mixture, K_theta_2, freqs){
    # all alleles are the same
    
    p1 = freqs[mixture[1]]
    Numerator = K_theta_2
    Denominator = F_theta(p1,6,theta) 
    LR = Numerator / Denominator
    return(LR)
  }
  
  func_1_2 <- function(mixture, K_theta_2, freqs){
    # allele_1 = allele2
    # then we may have two cases:
    # second contributor homozygous
    # second contributor heterozygous
    
    a1 = mixture[1]
    a2 = setdiff(mixture[c(3,4,5,6)],a1)[1] # the setdiff may result in more than one value (second and/or third contributor homozygous), so the [1] ensures just one allele and not a vector of more than one equal alleles
    p1 = freqs[a1]
    p2 = freqs[a2]
    # print(mixture)
    # print(a1)
    # print(freqs)
    # print(p1)
    # print(p2)
    numerator = K_theta_2 * (F_theta(p1, 2, theta)*((F_theta(p1, 3, theta)*2*F_theta(c(2,3) %*% c(p1, p2), 11, theta))+(4*F_theta(p2, 1, theta)*F_theta(p2, 2, theta)))+(F_theta(p2, 1, theta)*F_theta(p2, 2, theta)*F_theta(p2, 3, theta)))
    denominator = F_theta(p1, 2, theta)*(F_theta(p1, 3, theta)*(F_theta(p1, 4, theta)*F_theta(p1, 5, theta)*F_theta(c(6,15) %*% c(p1, p2), 51, theta)+F_theta(p2, 1, theta)*F_theta(p2, 2, theta)*5*F_theta(c(4,3) %*% c(p1, p2), 25, theta))+ 6*F_theta(p2, 1, theta)*F_theta(p2, 2, theta)*F_theta(p2, 3, theta)*F_theta(p2, 4, theta))
    LR = numerator / denominator
    return(LR)
  }
  
  func_1_3 <- function(mixture, K_theta_2, freqs){
    a1 = mixture[1]
    b = unique(setdiff(mixture[c (3, 4, 5, 6)], a1)) 
    #setdiff gives a vector with all the values different then a1 allele
    # unique(setdiff) returns the vector 'setdiff' but with duplicate elements removed.
    a2= b[1] # [1] allows selection of the first value of the vector 'unique(setdiff)'
    a3= b[2] # [2] allows selection of the second value of the vector 'unique(setdiff)'
    p1 = freqs[a1]
    p2 = freqs[a2]
    p3 = freqs[a3]
    numerator = K_theta_2 * (6*F_theta(p1, 2, theta)*F_theta(c(1,1,1) %*% c(p1, p2, p3), 5, theta)+F_theta(p2, 1, theta)*F_theta(c(2,3) %*% c(p2, p3), 7, theta)+2*F_theta(p3, 1, theta)*F_theta(p3, 2, theta))
    denominator = 15*F_theta(p1, 2, theta)*(F_theta(p1, 3, theta)*(F_theta(p1, 4, theta)*F_theta(c(1,2,2) %*% c(p1, p2, p3), 9, theta)+F_theta(p2, 1, theta)*F_theta(c(2,3) %*% c(p2, p3), 7, theta)+2*F_theta(p3, 1, theta)*F_theta(p3, 2, theta))+F_theta(p2, 1, theta)*(F_theta(p2, 2, theta)*F_theta(c(1,2) %*% c(p2, p3), 5, theta)+2*F_theta(p3, 1, theta)*F_theta(p3, 2, theta))+F_theta(p3, 1, theta)*F_theta(p3, 2, theta)*F_theta(p3, 3, theta))
    LR = numerator / denominator
    return(LR)
  }
  
  func_1_4 <- function(mixture, K_theta_2, freqs){
    a1 = mixture[1]
    b = unique(setdiff(mixture[c (3, 4, 5, 6)], a1))
    a2= b[1]
    a3= b[2]
    a4= b[3]
    p1 = freqs[a1]
    p2 = freqs[a2]
    p3 = freqs[a3]
    p4 = freqs[a4]
    numerator = K_theta_2 * (2+F_theta(c(1,1,1) %*% c(p2, p3, p4), 3, theta))
    denominator = 5*F_theta(p1, 1, theta)*F_theta(p1, 2, theta)*(F_theta(p1, 3, theta)*F_theta(c(2,3,3,3) %*% c(p1, p2, p3, p4), 17, theta)+F_theta(p2, 1, theta)*F_theta(c(2,3,3) %*% c(p2, p3, p4), 10, theta)+ F_theta(p3, 1, theta)*F_theta(c(2,3) %*% c(p3, p4), 7, theta)+ 2*F_theta(p4, 1, theta)*F_theta(p4, 2, theta))
    LR = numerator / denominator
    return(LR)
  }
  
  func_1_5 <- function(mixture, K_theta_2, freqs){
    p1 = freqs[mixture[1]]
    p2 = freqs[mixture[3]]
    p3 = freqs[mixture[4]]
    p4 = freqs[mixture[5]]
    p5 = freqs[mixture[6]]
    numerator = K_theta_2
    denominator = 30*F_theta(p1, 1, theta)*F_theta(p1, 2, theta)*F_theta(c(1,1,1,1,1) %*% c(p1, p2, p3, p4, p5), 7, theta)
    LR = numerator / denominator
    return(LR)
  }
  
  func_2_2 <- function(mixture, K_theta_2, freqs){
    p1 = freqs[mixture[1]]
    p2 = freqs[mixture[2]]
    numerator = K_theta_2 * (F_theta(p1, 1, theta)*(F_theta(p2, 1, theta)*(F_theta(p1, 2, theta)*2*F_theta(c(2,3) %*% c(p1, p2), 12, theta)+4*F_theta(p2, 2, theta)*F_theta(p2, 3, theta))+F_theta(p1, 2, theta)*F_theta(p1, 3, theta)*F_theta(p1, 4, theta))+F_theta(p2, 1, theta)*F_theta(p2, 2, theta)*F_theta(p2, 3, theta)*F_theta(p2, 4, theta))
    denominator = F_theta(p1, 1, theta)*F_theta(p2, 1, theta)*(F_theta(p1, 2, theta)*F_theta(p1, 3, theta)*(5*F_theta(p2, 2, theta)*F_theta(c(3,4) %*% c(p1, p2), 24, theta)+6*F_theta(p1, 4, theta)*F_theta(p1, 5, theta))+F_theta(p2, 2, theta)*F_theta(p2, 3, theta)*F_theta(p2, 4, theta)*3*F_theta(c(5,2) %*% c(p1, p2), 20, theta))
    LR = numerator / denominator
    return(LR)
  }
  
  func_2_3 <- function(mixture, K_theta_2, freqs){
    a1 = mixture[1]
    a2 = mixture[2]
    a3 = unique(setdiff(mixture[c (3, 4, 5, 6)], c(a1,a2)))
    p1 = freqs[a1]
    p2 = freqs[a2]
    p3 = freqs[a3]
    numerator = K_theta_2 *(F_theta(p1, 1, theta)*(2*F_theta(p1, 2, theta)*F_theta(c(2,6,3) %*% c(p1, p2, p3), 15, theta)+12*F_theta(p2, 1, theta)*F_theta(c(1,1) %*% c(p2, p3), 3, theta)+4*F_theta(p3, 1, theta)*F_theta(p3, 2, theta))+F_theta(p2, 1, theta)*(2*F_theta(p2, 2, theta)*F_theta(c(2,3) %*% c(p2, p3), 9, theta)+4*F_theta(p3, 1, theta)*F_theta(p3, 2, theta))+F_theta(p3, 1, theta)*F_theta(p3, 2, theta)*F_theta(p3, 3, theta))
    denominator = 30*F_theta(p1, 1, theta)*F_theta(p2, 1, theta)*(F_theta(p1, 2, theta)*(F_theta(p1, 3, theta)*F_theta(c(1,2,2) %*% c(p1, p2, p3), 10, theta)+F_theta(p2, 2, theta)*F_theta(c(2,3) %*% c(p2, p3), 9, theta)+2*F_theta(p3, 1, theta)*F_theta(p3, 2, theta))+F_theta(p2, 2, theta)*(F_theta(p2, 3, theta)*F_theta(c(1,2) %*% c(p2, p3), 6, theta)+2*F_theta(p3, 1, theta)*F_theta(p3, 2, theta))+F_theta(p3, 1, theta)*F_theta(p3, 2, theta)*F_theta(p3, 3, theta))
    LR = numerator / denominator
    return(LR)
  }
  
  func_2_4 <- function(mixture, K_theta_2, freqs){
    a1 = mixture[1]
    a2 = mixture[2]
    b = unique(setdiff(mixture[c (3, 4, 5, 6)], c(a1,a2)))
    a3 = b[1]
    a4 = b[2]
    p1 = freqs[a1]
    p2 = freqs[a2]
    p3 = freqs[a3]
    p4 = freqs[a4]
    numerator = K_theta_2 * (6*F_theta(p1, 1, theta)*F_theta(c(1,2,1,1) %*% c(p1, p2, p3, p4), 6, theta)+3*F_theta(p2, 1, theta)*F_theta(c(2,1,1) %*% c(p2, p3, p4), 6, theta)+F_theta(p3, 1, theta)*F_theta(c(2,3) %*% c(p3, p4), 7, theta)+2*F_theta(p4, 1, theta)*F_theta(p4, 2, theta))
    denominator = 30*F_theta(p1, 1, theta)*F_theta(p2, 1, theta)*(F_theta(p1, 2, theta)*F_theta(c(2,3,3,3) %*% c(p1, p2, p3, p4), 18, theta)+F_theta(p2, 2, theta)*F_theta(c(2,3,3) %*% c(p2, p3, p4), 12, theta)+F_theta(p3, 1, theta)*F_theta(c(2,3) %*% c(p3, p4), 7, theta)+2*F_theta(p4, 1, theta)*F_theta(p4, 2, theta))
    LR = numerator / denominator
    return(LR)
  }
  
  func_2_5 <- function(mixture, K_theta_2, freqs){
    a1 = mixture[1]
    a2 = mixture[2]
    b = unique(setdiff(mixture[c (3, 4, 5, 6)], c(a1,a2)))
    a3 = b[1]
    a4 = b[2]
    a5 = b[3]
    p1 = freqs[a1]
    p2 = freqs[a2]
    p3 = freqs[a3]
    p4 = freqs[a4]
    p5 = freqs[a5]
    numerator = K_theta_2 * F_theta(c(2,2,1,1,1) %*% c(p1, p2, p3, p4, p5), 7, theta)
    denominator = 30*F_theta(p1, 1, theta)*F_theta(p2, 1, theta)*F_theta(c(1,1,1,1,1) %*% c(p1, p2, p3, p4, p5), 7, theta)
    LR = numerator / denominator
    return(LR)
  }
  
  func_2_6 <- function(mixture, K_theta_2, freqs){
    p1 = freqs[mixture[1]]
    p2 = freqs[mixture[2]]
    numerator = K_theta_2
    denominator = 30*F_theta(p1,1,theta)*F_theta(p2,1,theta)
    LR = numerator/denominator
    return(LR)
  }
  
  # Dispatcher function that selects one of the func_1_1 to func_2_6 functions 
  # depending on the number of unique alleles (n_1) and number of contributors 
  # (n_mixture).
  # n1: number of unique alleles.
  # nmixture: number of contributors.
  # mixture, Ktheta2, freqs: passed to specific likelihood functions.
  # Returns likelihood ratio results for the given profile.
  func_split <- function(n_1,n_mixture,mixture, K_theta_2, freqs){
    stopifnot(n_1 %in% c(1,2))
    stopifnot(n_mixture %in% c(1,2,3,4,5,6))
    stopifnot((n_mixture - n_1) <= 4)
    stopifnot((n_mixture - n_1) >= 0)
    if(n_1 == 1){
      if (n_mixture == 1){
        func = func_1_1(mixture, K_theta_2, freqs)
      } else if (n_mixture == 2){
        func = func_1_2(mixture, K_theta_2, freqs)
      } else if (n_mixture == 3){
        func = func_1_3(mixture, K_theta_2, freqs)
      } else if (n_mixture == 4){
        func = func_1_4(mixture, K_theta_2, freqs)
      } else if (n_mixture == 5){
        func = func_1_5(mixture, K_theta_2, freqs)
      }
    } else if (n_1 == 2) {
      if (n_mixture == 2){
        func = func_2_2(mixture, K_theta_2, freqs)
      } else if (n_mixture == 3) {
        func = func_2_3(mixture, K_theta_2, freqs)
      } else if (n_mixture == 4) {
        func = func_2_4(mixture, K_theta_2, freqs)
      } else if (n_mixture == 5) {
        func = func_2_5(mixture, K_theta_2, freqs)
      } else if (n_mixture == 6) {
        func = func_2_6(mixture, K_theta_2, freqs)
      }
    }
  }
  
  
  n_1 <- count_number_unique_alleles(profiles = profiles[,1:2])
  n_mixture <- count_number_unique_alleles(profiles = profiles)
  
  
  #case_type = tibble(n_1, n_mixture) 
  index <- 1:nrow(profiles)
  out <- purrr::map_dbl(.x = index,
                        .f = ~{
                          n1 <- n_1[.x]
                          nMix = n_mixture[.x]
                          mixture <- profiles[.x,]
                          return(func_split(n1,nMix,mixture, K_theta_2, freqs))
                        }
  )
  if (remove_last_flag) out <- out[1]
  return (out)
}

# Generates a numerical seed from a character string by converting ASCII values.
# charstring: input string.
# Returns a numeric seed usable for random number generation.
char_to_seed <- function(char_string) {
  
  # convert each letter/number into a different ascii value
  ascii_vals <- as.integer(charToRaw(char_string))
  
  # Concatenate the ASCII values to form the seed
  seed_char <- paste(ascii_vals, collapse = "")
  
  # convert it back to a number
  seed_double <- as.double(seed_char)
  
  # depending on the size of the word, the seed_double can be huge
  # Huge numbers will not be good as seed (we can't print them), so the modulo 
  # against .Machine$integer.max (the maximum integer R can handle) will be considered
  seed <- seed_double %% .Machine$integer.max 
  
  return(seed)
}

# Calculates likelihood ratios per marker for simulated allele profiles.
# marker: target genetic marker.
# nsims: number of simulations.
# FREQ_ALELICAS: allele frequency data.
# seed: random seed.
# cuts: vector for cutting LR values into bins.
# Returns a data frame or tibble with LR values, allele counts, and summary groupings.
LR_per_marker <- function(marker,nsims = 1e6,FREQ_ALELICAS,seed = NULL,cuts = c(10^-4,10^-3,10^-2,10^-1,10^0,10^1,10^2,10^3,10^4)){
  
  freqs <- FREQ_ALELICAS[[marker]] 
  alleles <- as.character(FREQ_ALELICAS[["Allele"]])
  alleles <- alleles[!is.na(freqs)] 
  freqs <- freqs[!is.na(freqs)] 
  
  freqs <- freqs/sum(freqs)
  names(freqs) = alleles
  freqs_cum <- cumsum(freqs)
  
  if(is.null(seed)){
    set.seed(char_to_seed(marker))
  } else {
    set.seed(seed)
  }
  
  profiles <- generate_profiles(alleles,freqs_cum,nsims)
  theta = 0
  output_0 <- get_result_one_marker(theta,profiles,freqs)
  theta = 0.01
  output_001 <- get_result_one_marker(theta,profiles,freqs)
  
  output <-bind_cols(
    as_tibble(profiles) %>% setNames(c("A1","A2","A3","A4","A5", "A6")),
    tibble(
      n_1 = count_number_unique_alleles(profiles = profiles[,1:2]),
      n_mixture = count_number_unique_alleles(profiles = profiles[,1:6]),
      "LR[0]" = output_0,
      "LR[0.01]" = output_001,
      "LR[0]/LR[0.01]" = output_0 / output_001
    ) %>% 
      mutate(group = cut(`LR[0]/LR[0.01]`, cuts))
    
  )
  return(output)
} 

# Processes a list of data frames with LR results to extract summary statistics and counts by group.
# list_df: list of data frames.
# Returns a summary data frame including min, max, mean, quartiles, and group counts.
extract_summary_and_counts <- function(list_df) {
  # Apply function to each data.frame in the list
  result <- list_df %>%
    purrr::imap_dfr(function(df, marker) {
      # Calculate summary statistics
      summary_stats <- summary(df$`LR[0]/LR[0.01]`)
      # Calculate group counts
      group_counts <- table(df$group)
      group_percent <-round((group_counts %>% prop.table())*100,2)
      # Combine results into a tibble
      tibble(
        marker = marker,
        min = summary_stats["Min."],
        Q1 = summary_stats["1st Qu."],
        median = summary_stats["Median"],
        Q3 = summary_stats["3rd Qu."],
        max = summary_stats["Max."],
        mean = summary_stats["Mean"],
        !!!setNames(as.list(group_counts), paste0("count_", names(group_counts))),
        !!!setNames(as.list(group_percent), paste0("percentage_", names(group_counts)))
      )
    })
  return(result)
}

# Combines specified columns across a list of data frames by multiplication.
multiply_columns <- function(list_df, colname) {
  # Get the column from each data.frame and multiply them
  result <- Reduce(`*`, lapply(list_df, function(df) df[[colname]]))
  return(result)
}

# Extracts a specified column from each data frame in a list.
concatenate_columns <- function(list_df) {
  # Apply function to each data.frame in the list
  result <- purrr::imap_dfc(list_df, function(df, marker) {
    # Select the desired columns and rename them
    df %>%
      select(A1:A6) %>%
      rename_with(~ paste0(marker, "_", .))
  })
  return(result)
}


