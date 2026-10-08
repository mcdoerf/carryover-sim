######################Michael Doerfler#########################
#####################2026-10-08################################
#This script is for testing the 00_utils.R script in the R folder.

library(here)

source("./R/00_utils.R")

#First we test hist_mat.  This arranges capture histories in reverse lexicographic order.
#So, if it is working correctly, with L = 3 it should return,
#111 110 101 100 011 010 001 000

hist_mat(3)

hist_mat(2) #Should return 11 10 01 00

hist_mat(1) #Should return 1 0 

#Next we test hist_names.  Should work the same as hist_mat, except it is returning strings.

hist_names(3, observable = FALSE) #Should return 111 110 101 100 011 010 001 000

hist_names(3, observable = TRUE) #Should return 111 110 101 100 011 010 001


#Next we test theta_index(j,f).  We will do it for j = 1,2,3.  If it is working correctly, then
#theta_index(2, 1) should return 2, theta_index(2,0) should return 3, theta_index(3, 11) should return 4,
#theta_index(3, 10) should return 5, theta_index(3,01) should return 6, and theta_index(3,00) should return 7.

theta_index(2, 1); theta_index(2, "1")

theta_index(2,0); theta_index(2,"0")

theta_index(3, c(1,1)); theta_index(3, "11")

theta_index(3, c(1,0)); theta_index(3, "10")

theta_index(3, c(0,1)); theta_index(3, "01")

theta_index(3, c(0,0)); theta_index(3, "00")

#Next we test chain_pos(K).

chain_pos(3) #Should return 1, 3, 7

chain_pos(4) #Should return 1, 3, 7, 15

#Next we test cell_probs.

cell_probs(c(0.4, 0.6, 0.3), 2, psi = NULL) #I computed by hand that this should be
# 0.24, 0.16, 0.18, 0.42

cell_probs(c(0.4, 0.6, 0.3), 2, psi = 0.2) #I computed by hand that this should be
#0.048, 0.192, 0.032, 0.128, 0.036, 0.144, 0.084, 0.336

#Make sure the sum to 1

sum(cell_probs(c(0.4, 0.6, 0.3), 2, psi = NULL))

sum(cell_probs(c(0.4, 0.6, 0.3), 2, psi = 0.2))

#Make sure the function catches errors, likea mismathch between the legnth of theta and K.

cell_probs(c(0.4, 0.6, 0.3), 4, psi = 0.2) #Should return an error

cell_probs(c(1.2, 0.3, 0.2), psi = 0.2) #Should return an error


##Next I test first_capture_counts. Here is a toy cr_data object.

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

#If first_capture_counts is working correctly, then it should return a list with two vectors:
# b = (160, 75, 25) and a = (120, 45).

first_capture_counts(toy)
