#04/03/2025
#05/03/2025
#31/07/2025
#mixture rules
import sympy as sp

def one_fluid_theory(mix, cubic_alpha_r_list, delta=0, k_ij: list =[[0,0],[0,0]]):
    z = [mix.z1, mix.z2]
    n = 2
    #a_mix
    a_ij = [[(1-k_ij[i][j])*(cubic_alpha_r_list[i].a*cubic_alpha_r_list[j].a)**0.5 for i in range(n)] for j in range(n)]

    a_mix = sum(sum(z[i]*z[j]*a_ij[i][j] for i in range(n)) for j in range(n))

    b_ij = [[(1-delta)*0.5*(cubic_alpha_r_list[i].b + cubic_alpha_r_list[j].b) for i in range(n)] for j in range(n)]

    b_mix = sum(sum(z[i]*z[j]*b_ij[i][j] for i in range(n)) for j in range(n))

    return a_mix, b_mix



