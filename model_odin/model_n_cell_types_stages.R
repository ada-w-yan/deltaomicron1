# the TEIV model with cell types

# initial conditions
initial(T1[1:n_C]) <- T_0[i]
initial(L[1:n_C, 1:n_L]) <- L_0[i,j]
initial(I[1:n_C, 1:n_I]) <- I_0[i,j]
initial(V) <- V_0

# equations
deriv(T1[1:n_C]) <- - max(0, beta[i] * T1[i] * V)
deriv(L[1:n_C,1]) <- max(0, beta[i] * T1[i] * V) - max(0, k1 * n_L * L[i,j])
deriv(L[1:n_C,2:n_L]) <- max(0, k1 * n_L * L[i,(j-1)]) - max(0, k1 * n_L * L[i,j])
deriv(I[1:n_C,1]) <- max(0, k1 * n_L * L[i,n_L]) - max(0, delta * n_I * I[i,j])
deriv(I[1:n_C,2:n_I]) <- max(0, delta * n_I * I[i,(j-1)]) - max(0, delta * n_I * I[i,j])
deriv(V) <- max(0, sum(p_I_vec)) - max(0, c_inf * V) - max(0, sum(beta_inf_T_vec) * V)

dim(T1) <- n_C
dim(L) <- c(n_C,n_L)
dim(I) <- c(n_C,n_I)
dim(L_0) <- c(n_C,n_L)
dim(I_0) <- c(n_C,n_I)
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
L_0[,] <- user()
I_0[,] <- user()

p_I_vec[] <- p_inf * sum(I[i,])
beta_inf_T_vec[] <- beta[i] * T1[i]
n_C <- user()
n_L <- user()
n_I <- user()
