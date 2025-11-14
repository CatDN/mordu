# 07/11/2025

# create presets for all the pure fluid EOS relevant to the paper

from purefluid import NH3, H2

# classes
from eos import EOS
from alpha_r_cubic import AlphaRCubic
from alpha_r_helmholtz import AlphaRHelmholtz
from alpha_r_saft import AlphaRSAFT

# objects
from purefluid import NH3, H2
from cp0 import NH3_cp0_NIST, H2_cp0_NIST 
from presets_alpha_r import alpha_r_0290, alpha_r_0298, alpha_r_0300, alpha_r_0313

# algorithms and other functions and variables
from symbols import *
from other_functions import multi_root

# ============================================================ cubic EOS presets
# dictionaries for cubic EOS parameters for ammonia
vdW = {
    "alpha_r_expr":  sp.log((1/rho/(1/rho-b)))-a/(R*T*1/rho),
    "a_c_expr": 27/64*R**2*T_c**2/P_c,
    "alpha_T_expr": sp.sympify(1),
    "b_expr": 1/8*R*T_c/P_c
}

PR = {
    "alpha_r_expr": sp.log((1/rho/(1/rho-b))) + 1/(R*T)*a/(4*b)*2**0.5*sp.log(((1/rho-b*(2**0.5-1))/(1/rho+b*(2**0.5+1)))),
    "a_c_expr": 0.4572*R**2*T_c**2/P_c,
    "alpha_T_expr": (1+(0.37464+1.54226*omega - 0.2699*omega**2)*(1-(T/T_c)**0.5))**2,
    "b_expr": 0.0778*R*T_c/P_c
}

RK = {
    "alpha_r_expr": sp.log((1/rho/(1/rho-b))) + 1/(R*T)*a/b*sp.log((1/rho/(1/rho+b))),
    "a_c_expr": 0.4275*R**2*T_c**2/P_c,
    "alpha_T_expr": 1/T**0.5,
    "b_expr":  0.0867*R*T_c/P_c
}

MSRK = {
    "alpha_r_expr": sp.log(((1/rho/(1/rho-b)))) + 1/(R*T)*a/b*sp.log((1/rho/(1/rho+b))),
    "a_c_expr": 0.42748*R**2*T_c**2/P_c,
    "alpha_T_expr": (1+ (0.48503 + 1.5571*omega - 0.15613*omega**2)*(1-(T/T_c)**0.5))**2,
    "b_expr": 0.08664*R*T_c/P_c
}

# ============================================================ ideal gas law



# ============================================================ Hydrogen
H2_ideal = EOS("ideal",
                H2,
                H2_cp0_NIST,
                AlphaRHelmholtz,
                alpha_r_expr=sp.simplify(0))

# ============================================================ cubic EOS presets

H2_vdW = EOS("vdW",
              H2,
              H2_cp0_NIST,
              AlphaRCubic,
              **vdW)

H2_PR = EOS("PR",
             H2,
             H2_cp0_NIST,
             AlphaRCubic,
             **PR)

H2_RK = EOS("RK",
             H2,
             H2_cp0_NIST,
             AlphaRCubic,
             **RK)

H2_MSRK = EOS("MSRK",
               H2,
               H2_cp0_NIST,
               AlphaRCubic,
               **MSRK)

# ============================================================ Helmholtz EOS presets
H2_0313 = EOS("0313",
              H2,
              H2_cp0_NIST,
              AlphaRHelmholtz,
              alpha_r_expr = alpha_r_0313())


# ============================================================ Ammonia
NH3_ideal = EOS("ideal",
                NH3,
                NH3_cp0_NIST,
                AlphaRHelmholtz,
                alpha_r_expr=sp.simplify(0))

# ============================================================ cubic EOS presets
NH3_vdW = EOS("vdW",
              NH3,
              NH3_cp0_NIST,
              AlphaRCubic,
              **vdW)

NH3_PR = EOS("PR",
             NH3,
             NH3_cp0_NIST,
             AlphaRCubic,
             **PR)

NH3_RK = EOS("RK",
             NH3,
             NH3_cp0_NIST,
             AlphaRCubic,
             **RK)

NH3_MSRK = EOS("MSRK",
               NH3,
               NH3_cp0_NIST,
               AlphaRCubic,
               **MSRK)

# ============================================================ Helmholtz EOS presets
NH3_0290 = EOS("0290",
                NH3,
                NH3_cp0_NIST,
                AlphaRHelmholtz,
                alpha_r_expr=alpha_r_0290())

NH3_0298 = EOS("0298",
                NH3,
                NH3_cp0_NIST,
                AlphaRHelmholtz,
                alpha_r_expr=alpha_r_0298())

NH3_0300 = EOS("0300",
                NH3,
                NH3_cp0_NIST,
                AlphaRHelmholtz,
                alpha_r_expr=alpha_r_0300())

# ============================================================ SAFT EOS presets
# parameters
SAFT_0323 = {
    "epsilon": 124.3255,      #epsilon/k, in (K)
    "sigma": 2.2334,          #sigma, in (A) angstrom
    "m": 2.7078,              #segment length for ammonia, dimensionless
    "epsilon_AB": 1115.64,    #epsilon^AB/k, in (K)
    "k_AB": 0.4595,           #k^AB
    "M": 2,                    #number of association sites per molecule
    "a": [
        [0.9105631445, -0.3084016918, -0.0906148351],
        [0.6361281449, 0.1860531159, 0.4527842806],
        [2.6861347891, -2.5030047259, 0.5962700728],
        [-26.547362491, 21.419793629, -1.7241829131],
        [97.759208784, -65.255885330, -4.1302112531],
        [-159.59154087, 83.318680481, 13.776631870],
        [91.29774084, -33.746922930, -8.6728470368]
    ],
    "b": [
        [0.7240946941, -0.5755498075, 0.0976883116],
        [2.2382791861, 0.6995095521, -0.2557574982],
        [-4.0025849485, 3.8925673390, -9.1558561530],
        [-21.003576815, -17.215471648, 20.642075974],
        [26.855641363, 192.67226447, -38.804430052],
        [206.55133841, -161.82646165, 93.626774077],
        [-355.60235612, -165.20769346, -29.666905585]
    ],
    "association_scheme": "2B",
}

SAFT_0324 = {
    "epsilon": 204.63,      #epsilon/k, in (K)
    "sigma": 3.2386,          #sigma, in (A) angstrom
    "m": 1.1157,              #segment length for ammonia, dimensionless
    "epsilon_AB": 646.38,    #epsilon^AB/k, in (K)
    "k_AB": 0.00597,           #k^AB
    "M": 4,                   #number of association sites per molecule
    "association_scheme": "4B",
    "x_p": 1.3976,
    "mu": 1.469,
    "a": [
        [0.791982807, -0.623115538, -0.067775558],
        [1.071486513, 0.485734369, 0.028374114],
        [0.914746607, 1.124852696, 0.096122805],
        [-7.810606510, -2.094850156, 0.068150274],
        [25.78559770, 9.450498226, 0.059801866],
        [-56.98228765, -17.10272618, 0.286609791],
        [41.93089410, 7.776102807, -0.747016979]
    ],
    "b": [
        [0.791982807, -0.623115538, -0.067775558],
        [2.142973025, 0.971468739, 0.056748227],
        [2.744239820, 3.374558089, 0.288368414],
        [-31.24242604, -8.379400622, 0.272601097],
        [128.9279885, 47.25249113, 0.299009329],
        [-341.8937259, -102.6163571, 1.719658743],
        [ 293.5162587, 54.43271965, -5.229118852]
    ]
}

SAFT_0328 ={
    "epsilon": 283.18,      #epsilon/k, in (K)
    "sigma": 3.3476,          #sigma, in (A) angstrom
    "m": 1.503,              #segment length for ammonia, dimensionless
    "epsilon_AB": 893.1,    #epsilon^AB/k, in (K)
    "k_AB": 3.27e-2,           #k^AB
    "M": 3,                   #number of association sites per molecule
    "association_scheme": "3B",
}

# build specific alpha_disp term
tau = 0.74048
u = SAFT_0328["epsilon"]*(1+1/T)                                 #in K
d = SAFT_0328["sigma"]*(1-0.12*sp.exp(-3*SAFT_0328["epsilon"]/T))                #in angstrom

eta = pi/6*rho*SAFT_0328["m"]*d**3                        #dimensionaless


#Z_seg requires the universal constants D_ij available from [0336]
D_ij = [
    [-8.8043, 4.1646270, -48.203555, 140.43620, -195.23339, 113.51500, 0, 0, 0],
    [2.9396, -6.0865383, 40.137956, -76.230797, -133.70055, 860.25349, -1535.3224, 1221.4261, -409.10539], 
    [-2.8225, 4.7600148, 11.257177, -66.382743, 69.248785, 0, 0, 0, 0],
    [0.3400, -3.1875014, 12.231796, -12.110681, 0, 0, 0, 0, 0]
]


alpha_0_disp = sum([sum([D_ij[i][j]*(u/T)**(i+1)*(eta/tau)**(j+1) for j in range(0,9)]) for i in range(0,4)])

alpha_disp = SAFT_0328["m"]*alpha_0_disp

#fix the units
alpha_disp = alpha_disp.subs([(rho, rho*N_av*1e-30)])

SAFT_0330 ={
    "epsilon": 75.092,      #epsilon/k, in (K)
    "sigma": 2.2677,          #sigma, in (A) angstrom
    "m": 2.5485,              #segment length for ammonia, dimensionless
    "epsilon_AB": 1041.5,    #epsilon^AB/k, in (K)
    "k_AB": 0.37213,           #k^AB
    "M": 4,                   #number of association sites per molecule
    "association_scheme": "4C",
    "a": [
        [0.9105631445, -0.3084016918, -0.0906148351],
        [0.6361281449, 0.1860531159, 0.4527842806],
        [2.6861347891, -2.5030047259, 0.5962700728],
        [-26.547362491, 21.419793629, -1.7241829131],
        [97.759208784, -65.255885330, -4.1302112531],
        [-159.59154087, 83.318680481, 13.776631870],
        [91.29774084, -33.746922930, -8.6728470368]
    ],
    "b": [
        [0.7240946941, -0.5755498075, 0.0976883116],
        [2.2382791861, 0.6995095521, -0.2557574982],
        [-4.0025849485, 3.8925673390, -9.1558561530],
        [-21.003576815, -17.215471648, 20.642075974],
        [26.855641363, 192.67226447, -38.804430052],
        [206.55133841, -161.82646165, 93.626774077],
        [-355.60235612, -165.20769346, -29.666905585]
    ],
}


# SAFT
NH3_0323 = EOS("0323", NH3, NH3_cp0_NIST, AlphaRSAFT, **SAFT_0323)
NH3_0324 = EOS("0324", NH3, NH3_cp0_NIST, AlphaRSAFT, **SAFT_0324)
NH3_0328 = EOS("0328", NH3, NH3_cp0_NIST, AlphaRSAFT, **SAFT_0328,
               alpha_disp=alpha_disp)
NH3_0330 = EOS("0330", NH3, NH3_cp0_NIST, AlphaRSAFT, **SAFT_0330)