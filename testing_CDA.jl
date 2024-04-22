include("./ApproxMasEq.jl")
using .ApproxMasEq
using OrdinaryDiffEq, DiffEqCallbacks

println("Packages loaded")

# using Random
# include("./general.jl")
# include("./graphs.jl")
# include("./rates.jl")
# include("./KSAT.jl")
# include("./CME.jl")

# n = parse(Int64, ARGS[1])
# alpha = parse(Float64, ARGS[2])
N = 1000
alpha = 2.5
K = 3
c = K * alpha
p0 = 0.5
seed = 1

eta = 1.0
alg_str = "FMS"
rf = rate_FMS_KSAT
rargs = [eta]

t0 = 0.0
tlim = 10

saved_eners_CDA = SavedValues(Float64, Float64)
cb_ener_CDA = SavingCallback(save_ener_CDA, saved_eners_CDA)
cbs_save_CDA = CallbackSet(cb_ener_CDA)

tspan = [t0, tlim]

println("Running integration")
answ = CDA_KSAT(rf, rargs, build_args_rate_FMS, N=N, K=K, alpha=alpha, 
                seed_g=seed, seed_l=seed, tspan=[t0, tlim], cbs_save=cbs_save_CDA)


fileener = "CDA_KSAT_" * alg_str * "_K_" * string(K) * "_N_" * string(n) * "_alpha_" * 
           string(alpha) * "_p0_" * string(p0) * "_eta_" * string(eta) * "_t0_" * string(t0) * 
           "_tmax_" * string(tlim) * ".txt"

print_ener(saved_eners_CDA, fileener)

println("Integration finished")