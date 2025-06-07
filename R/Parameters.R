#Case 1 ==> R0 = 9.897
# Transmission rate of 0.30 (0.30 chance of a susceptible being infected after contact with an infectious individual)
# Average infectious period around 2 weeks (1/0.07)
# Average duration of immunity around 1 month (1/0.033)
params<- list("beta"= 0.30,
              "gamma"= 0.07,
              "rho"= 0.033,
              "epsilon"= 0.01,
              "n_patches"= 9,
              "max_time"= 365*8,
              "dt"= 1,
              "amplitude"= 0.02)

#Case 2 ==> R0 = 9.897
# Transmission rate of 0.30 (0.30 chance of a susceptible being infected after contact with an infectious individual)
# Average infectious period around 2 weeks (1/0.07)
# Average duration of immunity around 2 weeks (1/0.07)
params<- list("beta"= 0.30,
              "gamma"= 0.07,
              "rho"= 0.07,
              "epsilon"= 0.01,
              "n_patches"= 9,
              "max_time"= 365*8,
              "dt"= 1,
              "amplitude"= 0.02)

#Case 3 ==> R0 = 20.99474
# Transmission rate of 0.30 (0.30 chance of a susceptible being infected after contact with an infectious individual)
# Average infectious period around 1 month (1/0.033)
# Average duration of immunity around 2 weeks (1/0.07)
params<- list("beta"= 0.30,
              "gamma"= 0.033,
              "rho"= 0.07,
              "epsilon"= 0.01,
              "n_patches"= 9,
              "max_time"= 365*8,
              "dt"= 1,
              "amplitude"= 0.02)

