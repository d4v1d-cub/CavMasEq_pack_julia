include("./ApproxMasEq.jl")
using .ApproxMasEq
using OrdinaryDiffEq, DiffEqCallbacks, Sundials


N = parse(Int64, ARGS[1])
alpha = parse(Float64, ARGS[2])
K = parse(Int64, ARGS[3])
p0 = parse(Float64, ARGS[4])
seed = parse(Int64, ARGS[5])

eta = parse(Float64, ARGS[6])
alg_str = "FMS"
rargs = [eta]

t0 = 0.0
tlim = parse(Float64, ARGS[7])
abstol = parse(Float64, ARGS[8])
reltol = parse(Float64, ARGS[9])
custom = parse(Bool, ARGS[10])

if custom
    method=AM2_CDA

    fileener = "../Test/CDA_KSAT_custom_" * alg_str * "_K_" * string(K) * "_N_" * string(N) * "_alpha_" * 
               string(alpha) * "_p0_" * string(p0) * "_eta_" * string(eta) * "_t0_" * string(t0) * 
               "_tmax_" * string(tlim) * "_seed_" * string(seed) * ".txt"

    t_list, e_list = CDA_KSAT(rate_FMS_KSAT, rargs, build_args_rate_FMS, N=N, K=K, alpha=alpha, 
                         seed_g=seed, seed_l=seed, tspan=[t0, tlim], reltol=reltol, method=method, 
                         custom=true, fileener=fileener)
                         
else
    saved_eners_CDA = SavedValues(Float64, Float64)
    cb_ener_CDA = SavingCallback(save_ener_CDA, saved_eners_CDA)
    cbs_save_CDA = CallbackSet(cb_ener_CDA)
    
    
    method=CVODE_BDF(linear_solver = :GMRES)
    
    answ = CDA_KSAT(rate_FMS_KSAT, rargs, build_args_rate_FMS, N=N, K=K, alpha=alpha, 
                    seed_g=seed, seed_l=seed, tspan=[t0, tlim], cbs_save=cbs_save_CDA, dt_s=2.5, 
                    abstol=abstol, reltol=reltol, method=method)
    
    
    fileener = "../Test/CDA_KSAT_" * alg_str * "_K_" * string(K) * "_N_" * string(N) * "_alpha_" * 
               string(alpha) * "_p0_" * string(p0) * "_eta_" * string(eta) * "_t0_" * string(t0) * 
               "_tmax_" * string(tlim) * "_seed_" * string(seed) * ".txt"
    
    print_ener(saved_eners_CDA, fileener)
end

