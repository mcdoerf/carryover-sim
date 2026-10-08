##############Michael Doerfler###############
#############2026-10-08######################
#This R script is for testing the functions in the 01_data.R file.

library(here)

for (f in sort(list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}

#Testing cr_data()-----------------------------------

#creating counts1 and counts2.  Porposely putting them out of order to see if cr_data() 
#succesfully reorders them.

counts1 <- c("111" = 12, "110" = 48, "101" = 20, "100" = 80,   # cohort 1: 2^(K+1) - 1 = 7 cells
             "011" = 15, "010" = 60, "001" = 25)

counts1_permuted <- counts1[sample(length(counts1))]

counts2 <-  c("11"  = 30, "10"  = 90, "01"  = 45) 

counts2_permuted <- counts2[sample(length(counts2))]

cr_data(counts1, counts2, psi = 0.2, N_tot1 = 5000, N_tot2 = 8000, n_dbar = 470)

cr_data(counts1_permuted, counts2_permuted, psi = 0.2, N_tot1 = 5000, N_tot2 = 8000, n_dbar = 470)

#Also, check to see if cr_data catches wrong input.

#Negative counts
counts1_bad <- c("111" = -5, "110" = 48, "101" = 20, "100" = 80,   # cohort 1: 2^(K+1) - 1 = 7 cells
             "011" = 15, "010" = 60, "001" = 25)

cr_data(counts1_bad, counts2, psi = 0.2, N_tot1 = 5000, N_tot2 = 8000, n_dbar = 470)

#Psi not between 0 and 1

cr_data(counts1, counts2, psi = 1.2, N_tot1 = 5000, N_tot2 = 8000, n_dbar = 470)

#Cell counts exceed N_tot1

cr_data(counts1, counts2, psi = 0.2, N_tot1 = 50, N_tot2 = 8000, n_dbar = 470)


# Testing print.cr_data ----------------------------------------------------------

#Creating toy cr_data.


toy <- structure(
  list(
    K   = 2,
    psi = 0.1,
    counts = list(
      c("111" = 12, "110" = 48, "101" = 20, "100" = 80,   # cohort 1: 2^(K+1) - 1 = 7 cells
        "011" = 15, "010" = 60, "001" = 25),
      c("11"  = 30, "10"  = 90, "01"  = 45)                # cohort 2: 2^K - 1 = 3 cells
    ),
    N_tot  = c(5000, 8000),
    n_dbar = 470,
    stream_names = c("pathology", "hospital")
  ),
  class = "cr_data"
)


print.cr_data(toy)
