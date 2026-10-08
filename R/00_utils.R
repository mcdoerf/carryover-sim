########################Michael Doerfler###########################
########################2026-10-07#################################

#hist_mat(L): returns all 2^{L} capture histories as a 0/1 matrix in reverse lexicogrpahic order.
#L: history length.  K+1 for cohort 1, K for cohort 2.

hist_mat <- function(L){
  as.matrix(expand.grid(rep(list(c(1, 0)), L)))[, L:1]
}


#hist_names(L, observable = TRUE).  The same histories as strings, e.g. "101".  Drops the all-zero
#history when observable = TRUE.

hist_names <- function(L, observable = TRUE){
  mat <-   as.matrix(expand.grid(rep(list(c(1, 0)), L)))[, L:1, drop = FALSE]
  
  result <- apply(mat, 1, paste, collapse = "")
  
  if(observable == TRUE){
    result <- result[-length(result)]
  }
  
  return(result)
}

#theta_index(j, f): returns position of theta_{j|f} in a cohort's parameter vector.

theta_index <- function(j, f = integer(0)) {
  # f: 0/1 vector of length j - 1, or a string like "01"
  if (is.character(f)) f <- as.integer(strsplit(f, "")[[1]])
  stopifnot(length(f) == j - 1, all(f %in% c(0, 1)))
  fbin <- if (j == 1) 0 else sum(f * 2^((j - 2):0))
  2^j - 1 - fbin
}


#chain_pos(K): returns positions of the chain parameters, i.e. 2^{j-1}, j=1,...,K
chain_pos <- function(K){
  2^(1:K) - 1
}


#cell_probs(theta, K, psi = NULL): returns capture-history probabilities ph, incdluing the all 0
#history.  If psi is supplied, histories have length K+1 with the anchor as the last stream.
#Arguments: theta, numeric vector of length 2^{K}-1.  psi: anchor rate or NULL.
#Returns: named numeric vector summing to 1.

cell_probs <- function(theta, K, psi = NULL){
  
  stopifnot(all(0 <= theta & theta <= 1))
  stopifnot(length(theta) == 2^{K} - 1)
  
  L <- K

  if (!is.null(psi)) {
    L <- L + 1
    
    theta <- c(theta, rep(psi, 2^{K}))
  }
  
  h <- hist_mat(L)
  
  ph = rep(NA, 2^{L})

  names(ph) <- hist_names(L, observable = FALSE)

  for(i in (1:(2^{L}))){
    hi <- h[i,]
    
    ph[i] <- 1

    for(j in (1:L)){
      
    f <- if (j == 1) integer(0) else hi[1:(j-1)]
    
    index <- theta_index(j, f)
    
    thetaj <-theta[index]
    
    ph[i] <- ph[i] * (thetaj)^{hi[j]} * (1 - thetaj)^{1 - hi[j]}
    }
  }
  
  return(ph)
}



#first_capture_counts(data).
#purpose: the sufficient statistics for the cian: bj (cohort 1, first identified by stream j)
#b_{K+1}: identified only by the anchor, aj (cohort 2).
#Arguments: a cr_data object.
#Returns: list(b = <length K + 1>, a = <length K>)

first_capture_counts <- function(data){
  K <- data$K
  
  first_capture_names1 <- rep(NA, K+1)
  
  first_capture_names2 <- rep(NA, K)
  
  first_capture_names1[1] <- "1"
  
  first_capture_names2[1] <- "1"

  for (j in 2:(K)){
  name <- paste(c(rep(0, j-1), 1), collapse = "")
  
  first_capture_names1[j] <- name
  
  first_capture_names2[j] <- name
  }
  
  first_capture_names1[K+1] <- paste(c(rep(0, K), 1), collapse = "")
  
  b <- sapply(first_capture_names1, function(p){
    matches <- startsWith(names(data$counts[[1]]), p)
    
    sum(unlist(data$counts[[1]])[matches])
  })
  
  a <- sapply(first_capture_names2, function(p){
    matches <- startsWith(names(data$counts[[2]]), p)
    
    sum(unlist(data$counts[[2]])[matches])
  })
  
 return(list(b = b, a = a))
  
}



