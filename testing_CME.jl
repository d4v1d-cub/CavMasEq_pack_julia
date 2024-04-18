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

N = parse(Int64, ARGS[1])
alpha = parse(Float64, ARGS[2])
# N = 1000
# alpha = 2.5
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


saved_eners_CME = SavedValues(Float64, Float64)
cb_ener_CME = SavingCallback(save_ener_CME, saved_eners_CME)
cbs_save_CME = CallbackSet(cb_ener_CME)

tspan = [t0, tlim]

println("Running integration")
answ = CME_KSAT(rf, rargs, build_args_rate_FMS, N=N, K=K, alpha=alpha, 
                seed_g=seed, seed_l=seed, tspan=[t0, tlim], cbs_save=cbs_save_CME)


fileener = "CME_KSAT_" * alg_str * "_K_" * string(K) * "_N_" * string(n) * "_alpha_" * 
           string(alpha) * "_p0_" * string(p0) * "_eta_" * string(eta) * "_t0_" * string(t0) * 
           "_tmax_" * string(tlim) * ".txt"

print_ener(saved_eners_CME, fileener)

println("Integration finished")

Random.shuffle([1, 2, 3])

# efinal = 1e-4 * N
# graph = build_ER_HGraph(N, c, K, 1)
# links = gen_links(graph, 2)

# node = 10
# all_index = map(x -> get(x, node, "Error"), graph.nodes_in[graph.var_2_he[node]])

# for i in eachindex(all_index)
#     print(links[graph.var_2_he[node][i], all_index[i]])
# end

# all_lp, all_lm = all_lpm(graph, links)
# ch_u, ch_u_cond = unsat_ch(graph, links)
# p_cav = init_p_cav(graph, p0)
# probi = fill(p0, graph.N)
# len_cav = length(p_cav)

# params = graph, all_lp, all_lm, links, ch_u, ch_u_cond, rf, rargs, build_args_rate_FMS, len_cav, 
#          efinal

# u = vcat(reshape(p_cav, len_cav), probi)
# du = zeros(Float64, size(u))

# fder_KSAT_CME(du, u, params, 0.0)

# u .+= du * 0.01

# fder_KSAT_CME(du, u, params, 0.0)