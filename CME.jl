using Random
using OrdinaryDiffEq, ProgressLogging, DiffEqCallbacks

println("Packages loaded")

include("./general.jl")
include("./graphs.jl")
include("./rates.jl")
include("./KSAT.jl")
include("./export_results.jl")

# This function initializes the cavity conditional probabilities with a fully 
# independent initial distribution given by the vector 'p0'.
# It produces an array p_cav[he, i, val_cond, chain_others]
function init_p_cav(graph::HGraph, p0::Vector{Float64})
    ch_exc = graph.chains_he ÷ 2
    p_cav = zeros(Float64, (graph.M, graph.K, 2, ch_exc))
    for he in 1:graph.M
        for i in 1:graph.K
            for s in 1:2
                for ch in 0:ch_exc - 1
                    bits = digits(ch, base=2, pad=graph.K-1)
                    p_cav[he, i, s, ch + 1] = prod(bits + (1 .- 2 * bits) .* p0[graph.nodes_except[he, i, :]])
                end
            end
        end
    end
    return p_cav
end


# The user can just pass a single float 'p0' and the vector of initial conditions 
# is assumed to be homogeneous
function init_p_cav(graph::HGraph, p0::Float64)
    p0 = fill(p0, graph.N)
    init_p_cav(graph, p0)
end


# This function computes the derivative of the probabilities and feeds the julia's ODE integrator
function fder_KSAT(du::Vector{Float64}, u::Vector{Float64}, p, t::Float64)
    graph, all_lp, all_lm, links, ch_u, ch_u_cond, rfunc, rarg_cst, rarg_build, len_cav, efinal = p
    # These are the parameters of the integration
    p_cav = reshape(u[1:len_cav], (graph.M, graph.K, 2, graph.chains_he ÷ 2))
    probi = u[len_cav + 1:len_cav + graph.N]
    # The probabilities are reshaped in their original forms

    pu = comp_pu_KSAT(p_cav, graph, ch_u_cond)

    st = State(p_cav, probi, pu)  # A struct of type state is created just to pass it as an argument
                                  # for the builder of the rates' function arguments
                                  # This way, the user can choose what information to use inside the 
                                  # rate

    rates_arg = rarg_build(graph, st, ch_u, rarg_cst...)

    d_pc, d_pi = all_ders_CME_KSAT(p_cav, probi, pu, graph, all_lp, all_lm, rfunc, rates_arg, links, 
                                   ch_u_cond)

    du .= vcat(reshape(d_pc, len_cav), d_pi)
end


function save_ener(u, t, integrator)
    graph, all_lp, all_lm, links, ch_u, ch_u_cond, rfunc, rarg_cst, rarg_build,
               len_cav, efinal = integrator.p
    p_cav = reshape(u[1:len_cav], (graph.M, graph.K, 2, graph.chains_he ÷ 2))
    probi = u[len_cav + 1:len_cav + graph.N]
    pu = comp_pu_KSAT(p_cav, graph, ch_u_cond)
    e = ener(graph, probi, pu, ch_u)
    println(t, "\t", e, "\tdone")
    return e
end

function stopcond(u, t, integrator)
    graph, all_lp, all_lm, links, ch_u, ch_u_cond, rfunc, rarg_cst, rarg_build, 
              len_cav, efinal = integrator.p
    p_cav = reshape(u[1:len_cav], (graph.M, graph.K, 2, graph.chains_he ÷ 2))
    probi = u[len_cav + 1:len_cav + graph.N]
    pu = comp_pu_KSAT(p_cav, graph, ch_u_cond)
    return ener(graph, probi, pu, ch_u) - efinal
end

# This function integrates the CME's equations for a specific algorithm (given by ratefunc)
# and some boolean formula (given by graph and links)
function CME_KSAT(graph::HGraph, links::Matrix{Int8}, p0::Float64, ratefunc::Function, 
                  rargs_cst, rarg_build::Function, method, tspan::Vector{Float64}, 
                  efinal::Float64)
    all_lp, all_lm = all_lpm(graph, links)
    ch_u, ch_u_cond = unsat_ch(graph, links)
    p_cav = init_p_cav(graph, p0)
    probi = fill(p0, graph.N)
    len_cav = length(p_cav)
    params = graph, all_lp, all_lm, links, ch_u, ch_u_cond, ratefunc, rargs_cst, rarg_build, len_cav, 
             efinal
    u0 = vcat(reshape(p_cav, len_cav), probi)
    prob = ODEProblem(fder_KSAT, u0, tspan, params)

    saved_eners = SavedValues(Float64, Float64)
    cb_ener = SavingCallback(save_ener, saved_eners)
    affect!(integrator) = terminate!(integrator)
    cb_stop = ContinuousCallback(stopcond, affect!)

    cbs = CallbackSet(cb_ener, cb_stop)

    sol = solve(prob, method(), progress=true, callback=cbs)
    return sol, saved_eners
end


n = parse(Int64, ARGS[1])
alpha = parse(Int64, ARGS[2])
K = 3
c = K * alpha
p0 = 0.5

g1 = build_ER_HGraph(n, c, K, 1)
all_l = gen_links(g1, 1)

eta = 1.0
alg_str = "FMS"
rf = rate_FMS_KSAT
rargs = [eta]

t0 = 0.0
tlim = 1.0
dt_sample = 0.1

tspan = [t0, tlim]
method = VCAB4
t_save = collect(t0:dt_sample:tlim)

eth = 1e-6
efinal = eth * g1.N

println("Running integration")

answ, e_vals = CME_KSAT(g1, all_l, p0, rf, rargs, build_args_rate_FMS, method, tspan, efinal)

fileener = "CME_KSAT_" * alg_str * "_K_" * string(K) * "_N_" * string(n) * "_alpha_" * 
           string(alpha) * "_p0_" * string(p0) * "_eta_" * string(eta) * "_t0_" * string(t0) * 
           "_tmax_" * string(tlim) * "_dts_" * string(dt_sample) * "_eth_" * string(eth) * ".txt"

print_ener(e_vals, fileener)

println("Integration finished")