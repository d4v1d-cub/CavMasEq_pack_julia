include("./CME.jl")
using .CME
using OrdinaryDiffEq, DiffEqCallbacks

println("Packages loaded")

# n = parse(Int64, ARGS[1])
# alpha = parse(Float64, ARGS[2])
n = 100
alpha = 2.3
K = 3
c = K * alpha
p0 = 0.5
seed = 1

g1 = build_ER_HGraph(n, c, K, 1)
all_l = gen_links(g1, 1)

eta = 1.0
alg_str = "FMS"
rf = rate_FMS_KSAT
rargs = [eta]

t0 = 0.0
tlim = 0.2

tspan = [t0, tlim]

eth = 1e-6

println("Running integration")

saved_eners = SavedValues(Float64, Float64)
cb_ener = SavingCallback(save_ener, saved_eners)
cbs_save = CallbackSet(cb_ener)

answ = CME_KSAT(rf, rargs, build_args_rate_FMS, N=n, K=K, alpha=alpha, 
                        seed_g=seed, seed_l=seed, tspan=[t0, tlim], cbs_save=cbs_save)

saved_eners

fileener = "CME_KSAT_" * alg_str * "_K_" * string(K) * "_N_" * string(n) * "_alpha_" * 
           string(alpha) * "_p0_" * string(p0) * "_eta_" * string(eta) * "_t0_" * string(t0) * 
           "_tmax_" * string(tlim) * "_dts_" * string(dt_sample) * "_eth_" * string(eth) * ".txt"

print_ener(e_vals, fileener)

println("Integration finished")

CallbackSet(CallbackSet(),cb_ener)