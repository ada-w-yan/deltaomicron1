# the TEIV model with stages

# initial conditions
initial(T1) <- T_0
initial(L[1:n_L]) <- L_0[i]
initial(I[1:n_I]) <- I_0[i]
initial(V) <- V_0

# equations
deriv(T1) <- - max(0, beta * T1 * V)
deriv(L[1]) <- max(0, beta * T1 * V) - max(0, k1 * n_L * L[i])
deriv(L[2:n_L]) <- max(0, k1 * n_L * L[i-1]) - max(0, k1 * n_L * L[i])
deriv(I[1]) <- max(0, k1 * n_L * L[n_L]) - max(0, delta * n_I * I[i])
deriv(I[2:n_I]) <- max(0, delta * n_I * I[i-1]) - max(0, delta * n_I * I[i])
deriv(V) <- max(0, p_inf * sum(I)) - max(0, c_inf * V) - max(0, beta * T1 * V)

dim(L) <- c(n_L)
dim(I) <- c(n_I)
dim(L_0) <- c(n_L)
dim(I_0) <- c(n_I)

# parameter values
beta <- user()
k1 <- user()
delta <- user()
p_inf <- user()
c_inf <- user()
T_0 <- user()
V_0 <- user()
L_0[] <- user()
I_0[] <- user()
n_L <- user()
n_I <- user()
