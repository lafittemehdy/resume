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

X = {X₀, X₁, X₂, ...}
Σ = {Σ₀, Σ₁, Σ₂, Σ₃, Σ₄}
κ_T, κ_W, ζ, Y_v, ρ, σ, ε, ν, ξ, τ_m, k_m, k_t, S_thresh = constants
α₁, α₂, α₃, β₁, β₂, β₃, γ₁, γ₂, δ₁, δ₂, ε₁, ε₂ = constants

Φ(y) = max(0, min(1, Φ(y-1) - (α₁ * T_L(y-1) / max(Λ_m(y-1), ε)) + (α₂ * Λ_m(y-1) * (1-Φ(y-1)))))
Φ(0) = 1

Ψ(y) = max(0, Ψ(y-1) - β₁ * T_L(y-1) * (1 + Δ(y-1)))
Ψ(0) = Ψ_max

R(y) = max(0, min(R_max, R(y-1) - δ₁*Λ_s(y-1)*(1+Δ(y-1)) + δ₂*(R_max-R(y-1))/R_max - N(y-1)*ε₁))
R(0) = R_max

Θ(y) = max(0, Θ(y-1) * (1 - Δ(y-1) * β₂) * (R(y)/R_max))
Θ(0) = 1

κ_T_eff(y) = κ_T * Λ_s(y-1) * (Ψ(y-1)/Ψ_max) * Φ(y-1) * Θ(y-1)
Λ_s(y) = max(0, Λ_s(y-1) * (1 + κ_T_eff(y)))
Λ_s(0) = 1

B(y) = max(0, B(y-1) + ε₂ * (|X(y-1)| - |X(y-2)|) * κ_T_eff(y))
B(0) = 0

T_L(y) = max(0, T_L(y-1) * (1 + κ_T_eff(y)) * (1 + B(y-1)))
T_L(0) = 1

L_mem(y) = L_mem(y-1) * exp(-1/τ_m) + N(y-1)
L_mem(0) = 0
ζ_eff(y) = ζ * L_mem(y)
Λ_m(y) = max(0, Λ_m(y-1) * (1 + κ_W * (1-Δ(y-1)) * Θ(y-1)) + ζ_eff(y) * N(y-1) * Φ(y-1))
Λ_m(0) = 1

R_d(y) = (Λ_s(y) * T_L(y)) / (1 + Λ_m(y))
P_X_new(y) = min(1, 1 - exp(-ν * R_d(y)))
X(y) = X(y-1) ∪ {X_new if rand() < P_X_new(y)}

Δ(y) = max(0, min(1, Δ(y-1) + β₂ * |X(y)| * (T_L(y)/max(Λ_m(y), ε)) - β₃ * Λ_m(y) * Θ(y)))
Δ(0) = Δ₀

M_y = T_L(y) / max(Λ_m(y), ε)
Ω_i_base = Ω(i) for i ∈ X(y)
Ω_i_eff(y) = min(1, Ω_i_base * M_y * (1 + Δ(y)))
P₀(y) = 1 - ∏_{i∈X(y)} (1 - Ω_i_eff(y))

N(y) = P₀(y) * (1 - 1/(1 + exp(-5 * (Λ_m(y) - T_L(y)))))

R_s(y) = max(ε, R_s(y-1) * (1 - ρ * T_L(y) * (1+Δ(y)) * (1+B(y))))
R_s(0) = R_s_max
P₁(y) = (1 - P₀(y)) * (1 - exp(-σ / R_s(y)))

S_c(y) = (Λ_m(y) / max(T_L(y), ε)) * Φ(y)
P₂(y) = (1 - P₀(y)) * (1 - P₁(y)) * (1 / (1 + exp(-k_m * (S_c(y) - S_thresh))))

F_h = |X(y)| * T_L(y) * (1+Δ(y))
F_w = Λ_m(y) * Φ(y) * Θ(y)
P₃(y) = (1 - P₀(y)) * (1 - P₁(y)) * (1 - P₂(y)) * (1 / (1 + exp(-k_t * (F_w - F_h))))

P₄(y) = max(0, 1 - P₀(y) - P₁(y) - P₂(y) - P₃(y))

Σ_c(y) = Σ_c(y-1)
r = rand()
if r < P₀(y): 
    Σ_c(y) = Σ₀
else if r < P₀(y) + P₁(y):
    Σ_c(y) = Σ₁
    T_L(y) = γ₁ * T_L(y)
    Λ_s(y) = γ₁ * Λ_s(y)
    Λ_m(y) = γ₂ * Λ_m(y)
else if r < P₀(y) + P₁(y) + P₂(y): 
    Σ_c(y) = Σ₂
else if r < P₀(y) + P₁(y) + P₂(y) + P₃(y): 
    Σ_c(y) = Σ₃
else: 
    Σ_c(y) = Σ₄

Π(S) = Count(Σ_c(Y_v) == S) / N_sim for S ∈ Σ