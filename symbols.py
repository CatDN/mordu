import sympy as sp

# constants:
# universal gas constant
R = 8.31446261815324

# symbols:
# sympy symbols used throuout all files
# they are explicitly here so its easier to keep all files consistent

#symbols for pressure, temperature, density, specific volume, compressibility factor
P, T, rho, v, Z = sp.symbols("P T rho v Z", positive=True)

#symbols for cubic equations of state
a_c, a, b = sp.symbols("a_c a b", positive = True)
alpha_T = sp.symbols("alpha_T")

#symbols for mixture composition, z: overall composition, x: liquid phase, y: vapour phase
z1, z2, x1, x2, y1, y2 = sp.symbols("z1 z2 x1 x2 y1 y2", positive = True)

# n: total number of moles, n1: moles of component 1, n2: moles of component 2
n, n1, n2 = sp.symbols("n n1 n2", positive = True)

# volume, mass, molar mass
V, m, M = sp.symbols("V m M", positive=True)

# molar density, molar volume
rho_m, v_m  = sp.symbols("rho_m v_m", positive=True)

# critical temperature, critical pressure, critical density
T_c, P_c, rho_c = sp.symbols("T_c P_c rho_c", positive=True)

# reduced variables
rho_r, T_r = sp.symbols("rho_r T_r", positive=True)
# acentric factor
omega = sp.symbols("omega")

# non dimensional density, temperature, and inverse temperature
delta, tau, T_inv = sp.symbols("delta tau T_inv", positive=True)