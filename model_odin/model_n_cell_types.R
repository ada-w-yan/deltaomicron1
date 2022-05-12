# the TEIV model with cell types

# initial conditions
initial(T1[1:n_C]) <- T_0[i]
initial(L[1:n_C]) <- L_0[i]
initial(I[1:n_C]) <- I_0[i]
initial(V) <- V_0

# equations
deriv(T1[1:n_C]) <- - max(0, beta[i] * T1[i] * V)
deriv(L[1:n_C]) <- max(0, beta[i] * T1[i] * V) - max(0, k1 * L[i])
deriv(I[1:n_C]) <- max(0, k1 * L[i]) - max(0, delta * I[i])
deriv(V) <- max(0, sum(p_I_vec)) - max(0, c_inf * V) - max(0, sum(beta_inf_T_vec) * V)

dim(T1) <- n_C
dim(L) <- c(n_C)
dim(I) <- c(n_C)
dim(L_0) <- c(n_C)
dim(I_0) <- c(n_C)
dim(beta) <- n_C
dim(T_0) <- n_C
dim(p_I_vec) <- n_C
dim(beta_inf_T_vec) <- n_C

# parameter values
beta[] <- user()
k1 <- user()
delta <- user()
p_inf <- user()
c_inf <- user()
T_0[] <- user()
V_0 <- user()
L_0[] <- user()
I_0[] <- user()
p_I_vec[] <- p_inf * I[i]
beta_inf_T_vec[] <- beta[i] * T1[i]
n_C <- user()
