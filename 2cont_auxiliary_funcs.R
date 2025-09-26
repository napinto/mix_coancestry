# Function to generate allele profiles by sampling alleles according to cumulative frequencies
# alleles: vector of allele identifiers
# freqs_cum: vector of cumulative allele frequencies (sums to 1)
# nsims: number of simulated profiles to generate
generate_profiles <- function(alleles,freqs_cum,nsims){
  # Generate uniform random values for nsims simulations with 4 allele positions
  rand_numbers <- matrix(runif(nsims * 4),ncol = 4,byrow = T)
  
  # Find interval indices where the random numbers fall within cumulative frequency ranges
  indices <- findInterval(rand_numbers, freqs_cum,left.open = T,rightmost.closed = T) +1
  # Map indices back to alleles to generate simulated profiles
  profiles <- matrix(alleles[indices], nrow = nrow(rand_numbers), ncol = ncol(rand_numbers))
  return(profiles)
}

# Function to count the number of unique alleles present in each profile (row) of the profiles matrix
count_number_unique_alleles <- function(profiles){
  # Get full list of alleles in the profiles sorted
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

# Main function to calculate likelihood ratio results for a single marker
# using allele profiles and allele frequencies with theta correction
get_result_one_marker <- function(theta, profiles, freqs) {
  if(nrow(profiles)==1) {
    # Duplicate single profile to have at least two rows (required for calculation)
    profiles <- rbind(profiles,profiles)
    remove_last_flag <- TRUE
  } else {
    remove_last_flag <- FALSE
  }
  
  # Compute theta correction constant K_theta based on population substructure parameter theta
  K_theta = (1+3*theta) * (1 + 4*theta)
  
  # Auxiliary function for theta-adjusted allele frequency combination
  F_theta <- function(p,c,theta){
    return((1-theta)*p + c*theta)
  }
  
  # Below are different likelihood ratio formulas (functions)
  # corresponding to different allele mixture scenarios.
  # The mixture vector contains allele indices for the contributors.
  func_1_1 <- function(mixture, K_theta, freqs){
    # All alleles are the same: homozygous single contributor
    p1 = freqs[mixture[1]]
    Numerator = K_theta
    Denominator = F_theta(p1,4,theta) * F_theta(p1,5,theta) 
    LR = Numerator / Denominator
    return(LR)
  }

  func_1_2 <- function(mixture, K_theta, freqs){
    # allele_1 = allele2
    # then we may have two cases:
    # second contributor homozygous
    # second contributor heterozygous
    a1 = mixture[1]
    a2 = setdiff(mixture[c(3,4)],a1)[1] # Extract unique allele from second contributor
    p1 = freqs[a1]
    p2 = freqs[a2]
    
    numerator = K_theta * F_theta(c(2,1) %*% c(p1, p2), 5, theta)
    denominator = 2 * F_theta(p1, 2, theta) * (
      F_theta(p1,3,theta) * F_theta( c(2,3) %*% c(p1,p2), 11, theta) +
        2 * F_theta(p2, 1, theta) *F_theta( p2, 2 , theta)
    )
    LR = numerator / denominator
    return(LR)
  }
  
  func_1_3 <- function(mixture, K_theta, freqs){
    # mixtures = [i,i,j,k]
    p1 = freqs[mixture[1]] 
    p2 = freqs[mixture[3]]
    p3 = freqs[mixture[4]]
    numerator = K_theta
    denominator = 6*F_theta(p1,2,theta)*F_theta(c(1,1,1)%*%c(p1,p2,p3), 5, theta)
    LR = numerator / denominator
    return(LR)
  }
  func_2_2 <- function(mixture, K_theta, freqs){
    # mixtures = [i, j, i, i]
    # mixtures = [i, j, j, j]
    # mixtures = [i, j, i, j]
    # mixtures = [i, j, j, i]
    p1 = freqs[mixture[1]]
    p2 = freqs[mixture[2]]
    numerator = K_theta*(F_theta(p1,1,theta)*F_theta(c(1,2)%*%c(p1,p2),4,theta)+F_theta(p2,1,theta)*F_theta(p2,2,theta))
    denominator = 2*(F_theta(p1,1,theta)*F_theta(p2,1,theta))*(F_theta(p1,2,theta)*F_theta(c(2,3)%*%c(p1,p2), 12, theta)+2*F_theta(p2,2,theta)*F_theta(p2,3,theta))
    LR = numerator / denominator
    return(LR)
  }
  func_2_3 <- function(mixture, K_theta, freqs){
    # mixtures = [i, j, k, k]
    # mixtures = [i, j, i, k]
    # mixtures = [i, j, k, i]
    # mixtures = [i, j, j, k]
    # mixtures = [i, j, k, j]
    p1 = freqs[mixture[1]]
    p2 = freqs[mixture[2]]
    a3 = setdiff(mixture[c(3,4)],mixture[c(1,2)])[1]
    p3 = freqs[a3]
    numerator = K_theta*F_theta(c(2,2,1)%*%c(p1,p2,p3),5,theta)
    denominator = 12*(F_theta(p1,1,theta)*F_theta(p2,1,theta))*F_theta(c(1,1,1)%*%c(p1,p2,p3),5,theta)
    LR = numerator/denominator
    return(LR)
  }
  func_2_4 <- function(mixture, K_theta, freqs){
    p1 = freqs[mixture[1]]
    p2 = freqs[mixture[2]]
    numerator = K_theta
    denominator = 12*F_theta(p1,1,theta)*F_theta(p2,1,theta)
    LR = numerator/denominator
    return(LR) 
  }
  
  
  
  # Function to select the correct LR function depending on number of alleles and contributors
  func_split <- function(n_1,n_mixture,mixture, K_theta, freqs){
    stopifnot(n_1 %in% c(1,2))
    stopifnot(n_mixture %in% c(1,2,3,4))
    stopifnot((n_mixture - n_1) <= 2)
    stopifnot((n_mixture - n_1) >= 0)
    if(n_1 == 1){
      if (n_mixture == 1){
        func = func_1_1(mixture, K_theta, freqs)
      } else if (n_mixture == 2){
        func = func_1_2(mixture, K_theta, freqs)
      } else if (n_mixture == 3) {
        func = func_1_3(mixture, K_theta, freqs)
      } 
    } else if (n_1 == 2) {
      if (n_mixture == 2){
        func = func_2_2(mixture, K_theta, freqs)
      } else if (n_mixture == 3) {
        func = func_2_3(mixture, K_theta, freqs)
      } else if (n_mixture == 4) {
        func = func_2_4(mixture, K_theta, freqs)
      }
    }
  }
  
  # Calculate the number of unique alleles in first two and all four positions
  n_1 <- count_number_unique_alleles(profiles = profiles[,1:2])
  n_mixture <- count_number_unique_alleles(profiles = profiles)
  
  
  # Use purrr to apply func_split to each profile row
  index <- 1:nrow(profiles)
  out <- purrr::map_dbl(.x = index,
                        .f = ~{
                          n1 <- n_1[.x]
                          nMix = n_mixture[.x]
                          mixture <- profiles[.x,]
                          return(func_split(n1,nMix,mixture, K_theta, freqs))
                         }
  )
  # For single input profile, return a single value; otherwise return vector
  if (remove_last_flag) out <- out[1]
  return (out)
}


# Converts a character string into a numeric seed for random number generation
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

# Main wrapper to compute likelihood ratios per marker over simulated profiles
LR_per_marker <- function(marker,nsims = 1e6,FREQ_ALELICAS,seed = NULL,cuts = c(10^-4,10^-3,10^-2,10^-1,10^0,10^1,10^2,10^3,10^4)){
  
  freqs <- FREQ_ALELICAS[[marker]] 
  alleles <- as.character(FREQ_ALELICAS[["Allele"]])
  # Remove alleles with NA frequencies and normalize to sum to 1
  alleles <- alleles[!is.na(freqs)] 
  freqs <- freqs[!is.na(freqs)] 
  
  freqs <- freqs/sum(freqs)
  names(freqs) = alleles
  
  # Get cumulative allele frequencies for sampling
  freqs_cum <- cumsum(freqs)
  
  # Set seed for reproducibility, based on marker name or explicit seed
  if(is.null(seed)){
    set.seed(char_to_seed(marker))
  } else {
    set.seed(seed)
  }
  
  # Generate simulated allele profiles
  profiles <- generate_profiles(alleles,freqs_cum,nsims)
  theta = 0
  output_0 <- get_result_one_marker(theta,profiles,freqs)
  theta = 0.01
  output_001 <- get_result_one_marker(theta,profiles,freqs)
  
  # Combine profiles and LR results into a tibble with allele counts and groups
  output <- bind_cols(
    as_tibble(profiles) %>% setNames(c("A1","A2","A3","A4")),
    tibble(
      n_1 = count_number_unique_alleles(profiles = profiles[,1:2]),
      n_mixture = count_number_unique_alleles(profiles = profiles[,1:4]),
      "LR[0]" = output_0,
      "LR[0.01]" = output_001,
      "LR[0]/LR[0.01]" = output_0 / output_001
    ) %>% 
      mutate(group = cut(`LR[0]/LR[0.01]`, cuts))
    
  )
  return(output)
} 

# Summarizes and counts occurrence groups from list of data.frames with LR results
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

# Helper function to multiply specified columns across data.frames in a list
multiply_columns <- function(list_df, colname) {
  # Get the column from each data.frame and multiply them
  result <- Reduce(`*`, lapply(list_df, function(df) df[[colname]]))
  return(result)
}

# Helper function to concatenate specified columns from list of data.frames with renamed columns
concatenate_columns <- function(list_df) {
  # Apply function to each data.frame in the list
  result <- purrr::imap_dfc(list_df, function(df, marker) {
    # Select the desired columns and rename them
    df %>%
      select(A1:A4) %>%
      rename_with(~ paste0(marker, "_", .))
  })
  return(result)
}

