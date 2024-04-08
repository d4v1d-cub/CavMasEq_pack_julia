include("./ApproxMasEq.jl")
using .ApproxMasEq
using OrdinaryDiffEq, DiffEqCallbacks

println("Packages loaded")

# n = parse(Int64, ARGS[1])
# alpha = parse(Float64, ARGS[2])
n = 1000
alpha = 2.8
K = 3
c = K * alpha
p0 = 0.5
seed = 1

eta = 0.8
alg_str = "FMS"
rf = rate_FMS_KSAT
rargs = [eta]

t0 = 0.0
tlim = 1


saved_eners_CME = SavedValues(Float64, Float64)
cb_ener_CME = SavingCallback(save_ener_CME, saved_eners_CME)
cbs_save_CME = CallbackSet(cb_ener_CME)

tspan = [t0, tlim]

println("Running integration")
answ = CME_KSAT(rf, rargs, build_args_rate_FMS, N=n, K=K, alpha=alpha, 
                seed_g=seed, seed_l=seed, tspan=[t0, tlim], cbs_save=cbs_save_CME)

saved_eners_CDA = SavedValues(Float64, Float64)
cb_ener_CDA = SavingCallback(save_ener_CDA, saved_eners_CDA)
cbs_save_CDA = CallbackSet(cb_ener_CDA)

answ = CDA_KSAT(rf, rargs, build_args_rate_FMS, N=n, K=K, alpha=alpha, 
                seed_g=seed, seed_l=seed, tspan=[t0, tlim], cbs_save=cbs_save_CDA)


fileener = "CME_KSAT_" * alg_str * "_K_" * string(K) * "_N_" * string(n) * "_alpha_" * 
           string(alpha) * "_p0_" * string(p0) * "_eta_" * string(eta) * "_t0_" * string(t0) * 
           "_tmax_" * string(tlim) * ".txt"

print_ener(saved_eners, fileener)

println("Integration finished")