# the TEIV model with cell types

# initial conditions
initial(T1) <- T_0
initial(L) <- L_0
initial(I) <- I_0
initial(V) <- V_0

# equations
deriv(T1) <- - max(0, beta * T1 * V)
deriv(L) <- max(0, beta * T1 * V) - max(0, k1 * L)
deriv(I) <- max(0, k1 * L) - max(0, delta * I)
deriv(V) <- max(0, p_inf * I) - max(0, c_inf * V) - max(0, beta * T1 * V)

# parameter values
beta <- user()
k1 <- user()
delta <- user()
p_inf <- user()
c_inf <- user()
T_0 <- user()
V_0 <- user()
L_0 <- user()
I_0 <- user()
