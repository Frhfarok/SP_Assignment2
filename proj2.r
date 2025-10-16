#  Statistical Programming Group Project 2
#  Group 18 Members and Contributions:
#  Farah Nur Jannah Binti Farok (S2891743): Sections 1 and 2
#  Marta Gorecka (S1866561): Sections 2 & 3
#  Farell Zevic (S2810226): Sections 4 and 5

# ----------------------------------------------------------------

#  Code to simulate an SEIR epidemic model on a population network.
#  The population is divided into households, and each person also has a set of non-household contacts. 
#  Transmission can occur through household, close contact, or random mixing.
#  The code implements the following components:
#  1. Generate household assignments for each person.
#  2. Generate a network of non-household contacts based on sociability.
#  3. Simulate SEIR dynamics over time.
#  4. Visualize proportions of Susceptible, Exposed, Infectious, and Recovered individuals.
#  5. Compare four scenarios combining heterogeneous/constant sociability and random/full mixing.


# ------------------------ 1. Household distribution ---------------------------
#  Each person in the population belongs to a household.
#  Household sizes are uniformly distributed between 1 and hmax.
#  The vector h indicates which household each person belongs to; people with the same household number(ID) live in the same household.

set.seed(0)
n <- 1000 # population test size
hmax <- 5 # max household size
h <- rep(1:(n %/% hmax + 1), sample(1:hmax, n %/% hmax + 1, replace=TRUE))[1:n] # generate household IDs with sizes uniformly distributed between 1 and hmax


# ------------------------- 2.  Network of regular contacts --------------------
# function to generate non-household contacts
beta <- runif(n, 0, 1) #n vector of βi
get.net <- function(beta, h, nc = 15) {
  # beta = sociability vector all individuals
  # h = household membership vector
  # nc = average number of non-household contacts per person
  n <- length(beta) # total number of people
  bm <- mean(beta) # mean of vector beta
  contacts <- vector("list", n) # create empty list for each person
  
  for (i in seq_len(n)) {
    non_household <- which(h != h[i]) # identify people NOT in the same household as person i
    probij <- nc * beta[i] * beta[non_household] / (bm^2 * (n-1)) # probability links
    isampled <- runif(length(non_household)) < probij # randomly sample which non-householders to pick
    sampled <- non_household[isampled] # sampled non-householders
    contacts[[i]] <- sampled
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

# --------------------------- 3. SEIR model implementation ---------------------
# n = number of people
# βi = sociability parameter for ith person
# h = n vector of indicating which household the person belongs to
# beta = n vector of βi; 
# nt = number of days to simulate 
# pinf = proportion of the initial population to *randomly* start in the I state'

alink <- get.net(beta, h, nc = 15) # non-household contacts
nseir <- function(beta, h, alink, alpha=c(.1,.01,.01), delta=.2, gamma=.4, nc=15, nt = 100, pinf = .005) {

# calculate combined infection rate using alpha
# Find the total pool of people that could be encountered daily.
# Assign probabilities depending on which group they are from.
# Get the total probability:
# Probability of NOT getting infected: Pnotinf = (1-alpha_h)*(1-alpha_c)*(1-alpha_r)
# Probability of getting infected: Pinf = 1 - Pnotinf = 1 - (1-alpha_h)*(1-alpha_c)*(1-alpha_r)
# Update SEIR states
  
  bm <- mean(beta)
  n <- length(beta)
  x <- rep(0,n) # vector x to store the infection state
  ni = round(pinf*n) # number of infected people
  iinf <- sample(1:n,ni) # randomly select the infectious people
  x[iinf] <- 2 # set those infected people to stage "2" (E->I)
  S <- E <- I <- R <- numeric(nt) #set up storage for pop in each state
  S[1] <- n-ni # initial susceptible people = all - infected
  I[1]<-ni # initial infected people
  
  for (t in 2:nt){ # infecting process to start the day after the initial infections happened
    for (i in iinf){
      household <- which(h == h[i])
      if (length(alink[[i]]) == 0 ) { # pretty sure NA shouldn't appear at all, but I couldn't figure out how to fix it
        next
      }
      else{
        ipluscontacts <- c(alink[[i]],i)
        randmix <- setdiff(1:n, c(household, ipluscontacts))
        randmixsampled <- sample(randmix, nc) 
        alldaily <- c(household, alink[[i]], randmixsampled) # all people person i might encounter
        probh <- 0; probc <-0; probr <-0 # reset probabilities after every run
        
        for (j in alldaily){
          if (x[j] == 0){
            if (is.na(h[j])){
              next
            }
            else if (h[i] == h[j]){
              probh <- alpha[1] # household
            }
            else if (j %in% alink[[i]]){
              probc <- alpha[2] # contacts 
            }
            probr <- alpha[3] * nc * beta[i] * beta[j] / (bm^2 * (n-1)) # random mixing
            Pinf = 1 - (1-probh)*(1-probc)*(1-probr) # calculate total probability 
            if (runif(1) < Pinf) {
              x[j] <- 1
              }
          }
        }
      }
    }
    u <- runif(n)
    x[x == 2 & u < delta] <- 3   # I -> R
    x[x == 1 & u < gamma] <- 2   # E -> I
    iinf <- which(x == 2) # update infected people
    
    # store number of people in each state
    S[t] <- sum(x == 0) 
    E[t] <- sum(x == 1)
    I[t] <- sum(x == 2)
    R[t] <- sum(x == 3)
    
  }
  return(list(S=S,E=E,I=I,R=R,t=1:nt))
}

# ------------------ 4. Visualization of SEIR population states ----------------
plot.nseir <- function(result, main_title) {
  n <- result$S[1] + result$E[1] + result$I[1] + result$R[1]
  
  plot(result$t, result$S / n, type = "l", col = "blue",
       xlab = "Day", ylab = "Proportion of Population",
       ylim = c(0, 1), main = main_title)
  lines(result$t, result$E / n, col = "orange")
  lines(result$t, result$I / n, col = "red")
  lines(result$t, result$R / n, col = "green")
  
  legend("right", legend = c("S", "E", "I", "R"),
         col = c("blue", "orange", "red", "green"), lty = 1)
}

# ------------------ 5. Scenario-based SEIR model comparison -------------------
set.seed(50)

n <- 1000
nt <- 100
beta_vector <- runif(n, 0, 1)
h_households <- rep(1:(n %/% 5 + 1), sample(1:5, n %/% 5 + 1, replace=TRUE)) [1:n]
beta_bar_vector <- rep(mean(beta_vector), n)

par(mfrow = c(2, 2))

# Scenario 1: full model with variable beta
alink_full <- get.net(beta = beta_vector, h = h_households)
result_full <- nseir(beta = beta_vector, h = h_households, alink = alink_full)
plot.nseir(result_full, "Scenario 1: Full Model (Variable Beta)")

# Scenario 2: random mixing with variable beta
alpha_random_mixing <- c(0, 0, 0.4)
result_random <- nseir(beta = beta_vector, h = h_households, alink = alink_full,
                       alpha = alpha_random_mixing)
plot.nseir(result_random, "Scenario 2: Random Mixing (Variable Beta)")

# Scenario 3: full fodel with constant beta
alink_homog <- get.net(beta = beta_bar_vector, h = h_households)
result_homog <- nseir(beta = beta_bar_vector, h = h_households, alink = alink_homog)
plot.nseir(result_homog, "Scenario 3: Full Model (Constant Beta)")

# Scenario 4: random mixing with constant beta
result_rand_homog <- nseir(beta = beta_bar_vector, h = h_households, alink = alink_homog,
                           alpha = alpha_random_mixing)
plot.nseir(result_rand_homog, "Scenario 4: Random Mixing (Constant Beta)")
par(mfrow = c(1, 1))

