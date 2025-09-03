p = {f, h, m, t, e, l, s, d}
α, T, β, γ, λ, θ, φ, δ, μ, η, k, r₀ = constants

P_b = 1 - ∏_{i∈p} (1 - p_i)
P_e = min(1, P_b * (1 + α * T^β))
p_eff = p_i * (1 + γ * Σ_{j≠i,j∈p} p_j)

H = 1 / (1 + exp(-k * (r - r₀)))
T_r = ∏_{s=1}^{n} (1 - q_s)
T_e = T_r * H

E = 1 - exp(-λ * T * m)
T_eff = T * exp(-δ * D)

P_y = 1 - (1 - P_e) * (1 - E) * T_e
P_30 = 1 - ∏_{y=1}^{30} (1 - P_y)

P_err = θ * (1 - exp(-φ * τ))
P_a = f + h * H + t * (1 - T_r)
P_d = s * T_eff * P_err
P_t = 1 - (1 - P_a) * (1 - P_d)

F = 1 / (1 + η * N_n)
T_reg = 1 - ∏_{r=1}^{R} (1 - T_r)
P_c = 1 - ∏_{c=1}^{C} (1 - P_c_c)
P_g = min(1, P_t * (1 - F) * T_reg * exp(γ * T_eff) * (1 - P_c))

p_n = p_i * (1 - λ_l * N_i)
q_n = q_s * exp(-μ * u)
H_n = 1 / (1 + exp(-k * (r + Δr - r₀)))
T_rn = ∏_{s=1}^{n} (1 - q_n)
T_en = T_rn * H_n
E_n = 1 - exp(-λ * (T + ΔT) * m)

P_yn = 1 - (1 - P_e) * (1 - E_n) * T_en
P_an = f + h * H_n + t * (1 - T_rn)
P_dn = s * T_eff * P_err
P_cb = 1 - (1 - P_an) * (1 - P_dn)
T_n = T + ΔT_e - ΔT_d
P_cn = 1 - ∏_{y=1}^{N} (1 - P_cb(y))

S_l = 1 / (1 + exp(-k₁ * (d - d₀)))
S_f = exp(-λ_r * t) * exp(-α_r * r)
S_i = 1 / (1 + exp(-k₂ * (I - I₀)))
S_fd = min(1, A/G)
S_m = log(1 + M) / log(1 + M_max)
S_s = exp(-β * D)
S_mv = 1 + γ * R
S = S_l * S_f * S_i * S_fd * S_m * S_s * S_mv
S_adj = S * (1 - P_cn)

T_rad = -ln(S_th) / λ_r
T_inf = I₀ / r_I
T_fd = G / A_g
T_sc = D / r_D
T_rec = max(T_rad, T_inf, T_fd, T_sc)

S_n = max(0, min(1, S * (1 + μ * (d_n - d))))
P_fb = 1 / (1 + η * N_n)
T_er = 1 - ∏_{r=1}^{R} (1 - T_r)
P_gf = min(1, P_t * (1 - P_fb) * T_er * exp(γ * T_eff))

i_u = p_i * (1 - λ_l * N_i)
q_u = q_s * exp(-μ * u_l)
T_ru = ∏_{s=1}^{n} (1 - q_u)
T_eu = T_ru * H_n
E_u = 1 - exp(-λ * (T_n + ΔT) * m)
P_f = min(1, 1 - (1 - P_e) * (1 - E_u) * T_eu)

S_au = max(0, S_n * (1 - P_f))
Ω = max(0, min(1, S_au * exp(-γ * d_f) * exp(-δ * r_f)))
