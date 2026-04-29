#install.packages("deSolve") # R package to solve differential equations
#install.packages("cowplot")
#install.packages("tidyverse")

library(deSolve) # Load the deSolve package for solving differential equations
library(tidyverse)# Load the tidyverse package for data manipulation and visualization
library(cowplot) # Load the cowplot package for enhanced plot layout and customization

# Parameter List 

  Lv <- 10 # Life expectancy of mosquitoes (in days)
  Lh <- 50*365 # Life expectancy of humans (in days)
  IPh <-7 # Infectious period in humans (in days)
  IPv <-6 # Infectious period in mosquitoes (in days)
  EIP <-8.4  # Extrinsic incubation period in adult mosquitoes (in days)
  muv <- 1/Lv # Per capita mortality rate of mosquito population (1/Lv)
  muh <- 1/Lh # Per capita mortality rate of the human population (1/Lh)
  alphav <- muv # Per capita birth rate of the mosquito population. For now, we will assume that it is the same as the mortality rate.
  alphah <- muh # Per capita birth rate of the human population.  For now, we will assume that it is the same as the mortality rate
  gamma <- 1/IPh # Recovery rate in humans (1/IPh)
  delta <- 1/EIP# Extrinsic incubation rate (1/EIP)
  Nh <- 100000 # Number of humans. For this exercise, we suggest 100,000 humans. You can change this if you want according to the city you chose to model.
  m <- 2 # Density of female mosquitoes per human
  Nv <- m*Nh  # Number of mosquitoes (m * Nh)
  R0 <- 3 # Basic Reproduction Number
  ph <- 0.7 # Probability of transmission from an infectious mosquito to a susceptible human after a bite.
  pv <- 0.7 # Probability of transmission from an infectious human to a susceptible mosquito after a bite.
  b <- sqrt((R0 * muv*(muv+delta) * (muh+gamma)) /
              (m * ph * pv * delta)) # Biting rate
  betah <- ph*b # Coefficient of transmission from an infectious mosquito to a susceptible human after a bite (pH*B)
  betav <- pv*b # Coefficient of transmission from an infectious human to a susceptible mosquito after a bite (pv*b)
  TIME <- 1 # Number of years to be simulated. For this exercise, we will start with the first year of the epidemic.

  #Model equation  
  # Humans
  dSh   <-  alphah * Nh - betah * (Iv/Nh) * Sh - muh * Sh
  dIh   <- betah* (Iv/Nh) * Sh - (gamma + muh) * Ih
  dRh   <-  gamma * Ih  - muh * Rh
  
  # Mosquitoes
  dSv  <-  alphav * Nv - betav * (Ih/Nh) * Sv - muv * Sv
  dEv  <-  betav * (Ih/Nh) * Sv - (delata + muv)* Ev
  dIv  <-  delta * Ev - muv * Iv
  
  # Simple Deterministic Model (FUN)
  zika_model <- function(time, state_variable, parameters) {
    
    with(as.list(c(state_variable, parameters)), # local environment to evaluate derivatives
         {
           # Humans
           dSh   <-  alphah * Nh - betah * (Iv/Nh) * Sh - muh * Sh
           dIh   <- betah* (Iv/Nh) * Sh - (gamma + muh) * Ih
           dRh   <-  gamma * Ih  - muh * Rh
           
           # Mosquitoes
           dSv  <-  alphav * Nv - betav * (Ih/Nh) * Sv - muv * Sv
           dEv  <-  betav * (Ih/Nh) * Sv - (delta + muv)* Ev
           dIv  <-  delta * Ev - muv * Iv
           
           list(c(dSh, dIh, dRh, dSv, dEv, dIv))
         }
    )
  }
  
# Sequence of times (times)
  time <- seq(1, 365 * TIME , by = 1)

# Parameters (parms)
  parameters <- c(
    muv      = muv,
    muh      = muh,
    alphav   = alphav,
    alphah   = alphah,
    gamma    = gamma,
    delta    = delta,
    betav    = betav,
    betah    = betah,
    Nh       = Nh,
    Nv       = Nv
  )  
  
# Initial conditions of the system (y)
  
  start <- c(Sh = Nh-1 ,        
             Ih = 1 ,        # one infected human 
             Rh = 0 ,        
             Sv = Nv,        
             Ev = 0 ,        
             Iv = 0 )        
  
  # Solve the equations
out <- ode(y = start , 
             times = time ,   
             fun = zika_model ,   
             parms = parameters  
  ) %>%
    as.data.frame() # Convert to data frame

# Convert the times from days to years and weeks, respectively
out$years <- out$time / 365
out %>% View()

out$weeks <- out$time / 7
out %>% View()

# Visualise the dynamics
# Option 1: Separate plots on a grid
p1e <- ggplot(data = out, aes(y = Sh, x = weeks)) +
  geom_line(color = "royalblue", linewidth = 1) +
  labs(
    title = "Susceptible Human Population",
    x = "Weeks",
    y = "Number"
  ) +
  theme_bw()  # graph of susceptible human population

p2e <- ggplot(data = out, aes(y = Ih, x = weeks)) +
  geom_line(color = "firebrick", linewidth = 1) +
  labs(
    title = "Infectious Human Population",
    x = "Weeks",
    y = "Number"
  ) +
  theme_bw()  # graph of infectious human population

p3e <- ggplot(data = out, aes(y = Rh, x = weeks)) +
  geom_line(color = "olivedrab", linewidth = 1) +
  labs(
    title = "Recovered Human Population",
    x = "Weeks",
    y = "Number"
  ) +
  theme_bw()  # graph of recovered human population

plot_grid(p1e, p2e, p3e, ncol = 2)  # comparison graph of the susceptible, infectious human population, and

# Option 2: Show all compartments on one ggplot2 layer
# Reshape to long format
out_long <- out %>%
  pivot_longer(cols = -c(time, years, weeks),
               names_to = "state",
               values_to = "number"
  )

# Select human compartments
out_long_human <- out_long %>%
  filter(state %in% c("Sh", "Ih", "Rh"))

# Plot
sir_1_year_plot <- ggplot(data = out_long_human) +
  geom_line(aes(x = weeks, y = number, color = state), linewidth = 1) +
  labs(title = "Dynamics of the Human Population (1 year outbreak)") +
  theme_bw()

# Change the TIME parameter to 100 years and execute the model again
TIME2 <- 100  # Number of years to be simulated
# Sequence of times (times)
time2 <- seq(1, 365 * TIME2, by = 1)

out2 <- ode(
  y = start,
  times = time2,
  func = zika_model,
  parms = parameters
) %>%
  as.data.frame()  # Convert to data frame
# Convert the times from days to years and weeks, respectively
out2$weeks <- out2$time / 7
out2 %>% View()

out2$years <- out2$time / 365
#out2 %>% View()
# Examine the behavior of the model for 100 years
# Option 1. Plot as a grid
p1h <- ggplot(data = out2, aes(y = (Sh + Ih + Rh), x = years)) +
  geom_line(color = "grey68", linewidth = 1) +
  labs(
    title = "Total Human Population",
    x = "Years",
    y = "Number"
  ) +
  theme_bw()

p2h <- ggplot(data = out2, aes(y = Sh, x = years)) +
  geom_line(color = "royalblue") +
  labs(
    title = "Susceptible Human Population",
    x = "Years",
    y = "Number"
  ) +
  theme_bw()

p3h <- ggplot(data = out2, aes(y = Ih, x = years)) +
  geom_line(color = "firebrick", linewidth = 1) +
  labs(
    title = "Infectious Human Population",
    x = "Years",
    y = "Number"
  ) +
  theme_bw()
#several outbreaks within hundred years can be seen with the outbreak reducing

p4h <- ggplot(data = out2, aes(y = Rh, x = years)) +
  geom_line(color = "olivedrab", linewidth = 1) +
  labs(
    title = "Recovered Human Population",
    x = "Years",
    y = "Number"
  ) +
  theme_bw()

plot_grid(p1h, p2h, p3h, p4h, ncol = 2)

# Option 2. Plot all on one layer
# Reshape the data
out2_long <- out2 %>%
  pivot_longer(
    cols = -c(time, weeks, years),
    names_to = "state",
    values_to = "number"
  )

# Filter human compartments
out2_long_human <- out2_long %>%
  filter(state %in% c("Sh", "Ih", "Rh"))

# Make the plot
sir_100_years_plot <- ggplot(data = out2_long_human) +
  geom_line(aes(x = years, y = number, color = state), linewidth = 1) +
  labs(title = "Dynamics of the Human Population (100 years outbreak)") +
  theme_bw()
# Proportions
p1 <- ggplot(data = out2, aes(y = Sh / (Sh + Ih + Rh), x = years)) +
  geom_line(color = "royalblue", linewidth = 1) +
  labs(
    title = "Susceptible Human Population",
    x = "Years",
    y = "proportion"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1))

p2 <- ggplot(data = out2, aes(y = Ih / (Sh + Ih + Rh), x = years)) +
  geom_line(color = "firebrick", linewidth = 1) +
  labs(
    title = "Infectious Human Population",
    x = "Years",
    y = "proportion"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1))

p3 <- ggplot(data = out2, aes(y = Rh / (Sh + Ih + Rh), x = years)) +
  geom_line(color = "olivedrab", linewidth = 1) +
  labs(
    title = "Recovered Human Population",
    x = "Years",
    y = "proportion"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1))


  
