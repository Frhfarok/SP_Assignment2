#  Statistical Programming Group Project 2
#  Group 18 Members and Contributions:

#  Farah Nur Jannah Binti Farok (S2891743): Sections 1 and 2
#  Marta Gorecka (S1866561): 
#  Farell Zevic (S2810226): 
############################################################

# 1. Household distribution

set.seed(0)
n <- 1000 # population test size
hmax <- 5 # max household size
h <- rep(1:(n %/% hmax + 1), sample(1:hmax, n %/% hmax + 1, replace=TRUE)) [1:n]
    # n %/% hmax is the min num of households needed if all had hmax people
    # +1 to have extra households in case division is not exact
    # replace=TRUE allows repeated household sizes
    # [1:n] trims the vector to exactly n people

# 2.  Network of regular contacts

# function to generate non-household contacts
get.net <- function(beta, h, nc = 15) {
  n <- length(beta) # total number of people
  contacts <- vector("list", n) # create empty list for each person
  
  for (i in 1:n) {
    non_household <- which(h != h[i]) # identify people NOT in the same household as person i
    contacts[[i]] <- sample(non_household, nc) # sample nc contacts from non_household population
  }
  
  # make the link bidirectional
  # ensure if person i lists person j as a contact, then person j also lists person i as a contact
  for (i in 1:n) {
    for (j in contacts[[i]] ) {
      # check if the reverse link (j -> i) already exists
      if (! (i %in% contacts[[j]]) ) { 
        contacts[[j]] <- c(contacts[[j]], i) # add i to j's contact list if bidirectional link does not exist
      }}}
  return(contacts) }

 
# 3. SEIR model implementation
# 4. Visualization of SEIR population states
# 5. Scenario-based SEIR model comparison
