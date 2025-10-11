#  Statistical Programming Group Project 2
#  Group 18 Members and Contributions:

#  Farah Nur Jannah Binti Farok (S2891743): Sections 1 and 2
#  Marta Gorecka (S1866561): Sections 2 & 3
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
beta <- runif(n, 0, 1) #n vector of βi
get.net <- function(beta, h, nc = 15) {
  # beta = 
  #
  #
  n <- length(beta) # total number of people
  bm <- mean(beta) # mean of vector beta
  contacts <- vector("list", n) # create empty list for each person
  
  for (i in 1:n) {
    non_household <- which(h != h[i]) # identify people NOT in the same household as person i
    
    probij <- nc * beta[i] * beta[non_household] / (bm^2 * (n-1)) # probability links
    isampled <- runif(length(non_household)) < probij # randomly sample which non-householders to pick
    sampled <- non_household[isampled] # sampled non-householders
    contacts[i] <- sampled
    #contacts[[i]] <- sample(non_household, nc, probij) # sample contacts from non_household population with probability probij
    
    }
  
  # make the link bidirectional
  # ensure if person i lists person j as a contact, then person j also lists person i as a contact
  for (i in 1:n) {
    for (j in contacts[[i]] ) {
      # check if the reverse link (j -> i) already exists
      if (! (i %in% contacts[[j]]) ) { 
        contacts[[j]] <- c(contacts[[j]], i) # add i to j's contact list if bidirectional link does not exist
      }}}
  return(contacts) 
  }

 
# 3. SEIR model implementation
# n = number of people
# βi = sociability parameter for ith person
# h = n vector of indicating which household the person belongs to
# beta = n vector of βi; nt = number of days to simulate 
# pinf = proportion of the initial population to *randomly* start in th eI state

## NOT FINISHED, To Be Completed
alink = get.net(beta, h, nc = 15) # non-household contacts
nseir <- function(beta, h, alink, alpha=c(.1,.01,.01), delta=.2, gamma=.4, nc=15, nt = 100, pinf = .005) {
  
  ni = round(pinf*length(beta)) # number of infected people
  iinf <- sample(n,ni) # randomly select the infectious people
  x[iinf] <- 2 # set those infected people to stage "2" (E->I)
  
  S <- E <- I <- R <- rep(0,nt) ## set up storage for pop in each state
  S[1] <- n-ni; I[1]<-ni
  
  # from LECTURE NOTES for inspiration
  for (i in 2:nt) { ## loop over days
    u <- runif(n) ## uniform random deviates
    x[x==2 & u < delta] <- 3 ## I -> R with prob delta
    x[x==1 & u < gamma] <- 2 ## E -> I with prob gamma
    x[x==0 & u < beta*I[i-1]] <- 1 ## S -> E with prob beta*I[i-1]
    S[i] <- sum(x==0); E[i] <- sum(x==1)
    I[i] <- sum(x==2); R[i] <- sum(x==3)
  }
  return(list(S=S,E=E,I=I,R=R,t=beta))
}

# 4. Visualization of SEIR population states
# 5. Scenario-based SEIR model comparison

