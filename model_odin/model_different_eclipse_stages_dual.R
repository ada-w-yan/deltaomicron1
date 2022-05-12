# the TEIV model with cell types

# initial conditions
initial(T1[1:n_C]) <- T_0[i]
# L[1] is tmprss2+ cells infected by the endosomal pathway
# L[2] is tmprss2+ cells infected by the tmprss2 pathway
# L[3] is tmprss2- cells infected by the endosomal pathway
initial(L[1:n_C_plus_one, 1:n_L]) <- L_0[i,j]
initial(I[1:n_C, 1:n_I]) <- I_0[i,j]
initial(V) <- V_0
initial(V_tot) <- V_tot0

# equations
deriv(T1[1]) <- - max(0, (beta_endosomal + beta_tmprss2) * T1[i] * V)
deriv(T1[2]) <- - max(0, beta_endosomal * T1[i] * V)
deriv(L[1,1]) <- max(0, beta_endosomal * T1[1] * V) - max(0, k1_endosomal * n_L * L[i,j])
deriv(L[2,1]) <- max(0, beta_tmprss2 * T1[1] * V) - max(0, k1_tmprss2 * n_L * L[i,j])
deriv(L[3,1]) <- max(0, beta_endosomal * T1[2] * V) - max(0, k1_endosomal * n_L * L[i,j])
deriv(L[1,2:n_L]) <- max(0, k1_endosomal * n_L * L[i,(j-1)]) - max(0, k1_endosomal * n_L * L[i,j])
deriv(L[2,2:n_L]) <- max(0, k1_tmprss2 * n_L * L[i,(j-1)]) - max(0, k1_tmprss2 * n_L * L[i,j])
deriv(L[3,2:n_L]) <- max(0, k1_endosomal * n_L * L[i,(j-1)]) - max(0, k1_endosomal * n_L * L[i,j])
deriv(I[1,1]) <- max(0, k1_endosomal * n_L * L[1,n_L] + k1_tmprss2 * n_L * L[2,n_L]) - max(0, delta * n_I * I[i,j])
deriv(I[2,1]) <- max(0, k1_endosomal * n_L * L[3,n_L]) - max(0, delta * n_I * I[i,j])
deriv(I[1:n_C,2:n_I]) <- max(0, delta * n_I * I[i,(j-1)]) - max(0, delta * n_I * I[i,j])
deriv(V) <- max(0, sum(p_I_vec)) - max(0, c_inf * V) - max(0,  ((beta_endosomal + beta_tmprss2) * T1[1] + beta_endosomal * T1[2]) * V)
deriv(V_tot) <- max(0, sum(p_I_vec) * p_tot / p_inf) - max(0, c_tot * V_tot) - max(0,  ((beta_endosomal + beta_tmprss2) * T1[1] + beta_endosomal * T1[2]) * V)

dim(T1) <- n_C
dim(L) <- c(n_C_plus_one,n_L)
dim(I) <- c(n_C,n_I)
dim(L_0) <- c(n_C_plus_one,n_L)
dim(I_0) <- c(n_C,n_I)
dim(T_0) <- n_C
dim(p_I_vec) <- n_C

# parameter values
beta_endosomal <- user()
beta_tmprss2 <- user()
k1_endosomal <- user()
k1_tmprss2 <- user()
delta <- user()
p_inf <- user()
p_tot <- user()
c_inf <- user()
c_tot <- user()
T_0[] <- user()
V_0 <- user()
V_tot0 <- user()
L_0[,] <- user()
I_0[,] <- user()

p_I_vec[] <- p_inf * sum(I[i,])
n_C <- user()
n_L <- user()
n_I <- user()
n_C_plus_one <- n_C + 1
