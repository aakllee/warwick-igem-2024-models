using DifferentialEquations

# Define differential equations
function stress!(du, u, p, t)
    H, L, M, s = u
    b_s, n_H1, n_H2, k_L, k_I, I_L, I_M, α, k_maxuptake, L_out, k_M, β, δ, n_M, γ, k_stress, k_recovery = p

    # Functions
    χ(s) = 1 / (1 + s / b_s)                          # Scaling function for response of dH to stress
    h(M,L) = δ*M*L                                    # ?????????
    MM(H) = (k_maxuptake * H * L_out) / (k_M + L_out) # Michaelis-Menten kinetics
    Hill(L,n,K_d) = L^n / (K_d + L^n)                 # Hill equation; L=[ligand], n=Hill coefficient, K_d=dissossiation constant
    HillR(L,n,K_d) = K_d^n / (K_d^n + L^n)

    # Differential equations
    dH = χ(s) * (HillR(L,n_H1,k_L) + Hill(I_L,n_H2,k_I)) - α*H
    dL = MM(H) - β*L - h(M,L)
    dM = 5* χ(s) * Hill(I_M,n_M,k_I) - γ*M - h(M,L)
    ds = k_stress*L + dM + dH - k_recovery*s

    # Return
    du .= [dH, dL, dM, ds]
end


# Set initial conditions
u0 = [
    1,            # H(t)        [lutH]
    0,            # L(t)        [Ln] inside cell
    0.5,            # M(t)        [unbound truncated lanM] in cytosol
    0.5             # s(t)        stress of cell
]

# Set time domain
tspan = (0.0, 500.0) # t

# Set parameter values
parameters = [
    1.0,            # b_s         base level of stress of the cell
    1.2,            # n_H1        Hill exponent for lanthanides' activation of lutH
    2.0,            # n_H2        Hill exponent for IPTG activation of lutH
    3.0,            # k_L         dissociation constant for lanthanides
    49.6,            # k_I         dissociation constant for IPTG
    100,            # I_L           [IPTG] lutH
    10000,            # I_M           IPTG LanM
    0.056,            # α           degradation rate of lutH
    0.3125,            # k_maxuptake maximum uptake of lanthanides
    1000.0,            # L_out       [Lanthanides outside the cell]
    1.0,            # k_m         Michaelis—Menton model constant
    0.1,            # β           degradation rate of lanthanides
    0.1,            # δ           formation rate of lanM—La complexes
    3.0,            # n_M         Hill exponent for IPTG activation of truncated lanM
    0.01,            # γ           degradation rate of truncated lanM
    0.1,            # k_stress    proportionality constant between [lanthanides] and stress
    0.1             # k_recovery  rate at which the cell recovers from stress
]



using Plots

# Solve problem
sol = solve(
    ODEProblem(stress!, u0, tspan, parameters),
    Tsit5())

# Plot and save to plot.png
plot(sol, label=["H" "L" "M" "s"], xlabel="time", dpi = 300)
savefig("plot.png")
