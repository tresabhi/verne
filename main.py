import pint

ur = pint.UnitRegistry()

d_tau = 2**-8

R_u = 8.3145 * ur.J / (ur.mol * ur.K)

M_C_2_H_4 = 28.05 * ur.g / ur.mol
M_O_2 = 32.00 * ur.g / ur.mol

o_f = (3 * M_O_2) / (1 * M_C_2_H_4)

M_bar = (1 * M_C_2_H_4 + 3 * M_O_2) / (1 + 3)
R_g = R_u / M_bar

# combustion
q_0 = h_f / (gamma * R_g * T_hat * (1 + o_f))
T_c_cv = T_tm + gamma * (gamma - 1) * q_0
P_c_final = P_tm  # initial approximation
P_c_cv / P_c_final = T_c_cv / T_tm
P = rho * T

# blowdown
d_rho_c = -rho_th * v_th * A_th * d_tau
d_P_c = (
    -gamma
    * (T_c * rho_th * v_th * A_th + beta * (rho_th * v_th) ** 0.8 * (T_c - T_wall))
    * d_tau
)

# refill
d_P_c = gamma * (T_tm * rho_in * v_in - T_c * rho_th * v_th * A_th) * d_tau
rho_in * v_in = (
    (
        (2 / (gamma - 1)) ** (1 / 2)
        * (P_tm / T_tm ** (1 / 2))
        * (P_c / P_tm) ** (1 / gamma)
        * (1 - (P_c / P_tm) ** ((gamma - 1) / gamma)) ** (1 / 2)
    )
    if P_c / P_tm > (2 / gamma + 1) ** (gamma / (gamma - 1))
    else (P_tm / T_tm ** (1 / 2))
    * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
)

# after convergence
w = P_c_final / (T_tm * tau_cy)
