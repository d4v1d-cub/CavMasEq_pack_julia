module CME

using Random
using OrdinaryDiffEq, DiffEqCallbacks

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
function CME_KSAT_base(ratefunc::Function, rargs_cst, rarg_build::Function, 
                       graph::HGraph, links::Matrix{Int8}, tspan::Vector{Float64}, p0::Float64, 
                       method, eth::Float64)
    efinal = eth * graph.N
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


# Decorator for the function integrates the CME's equations for a specific algorithm (given by ratefunc)
# and some boolean formula (given by graph and links)
function CME_KSAT(ratefunc::Function, rargs_cst, rarg_build::Function;
                  graph::HGraph=build_empty_graph(), 
                  N::Int64=0, K::Int64=0, alpha::Float64=0.0, seed_g::Int64=rand(1:typemax(Int64)),
                  links::Matrix{Int8}=Matrix{Int8}(undef, 0, 0), seed_l::Int64=rand(1:typemax(Int64)), 
                  tspan::Vector{Float64}=[0.0, 1.0], 
                  p0::Float64=0.5, method=VCAB4, eth::Float64=1e-6)
    if N > 0
        if graph.N == 0
            c = K * alpha
            graph = build_ER_HGraph(N, c, K, seed_g)
        end

        if length(links) == 0
            links = gen_links(graph, seed_l)
        end

        return CME_KSAT_base(ratefunc, rargs_cst, rarg_build, graph, links, tspan, p0,
        method, eth)
    else
        throw("In CME_KSAT function: The user should provide either a graph::HGraph or valid values for N, K and alpha")
    end  
end


export CME_KSAT, HGraph, build_ER_HGraph, build_RR_HGraph, gen_links, print_ener, 
       ener, rate_FMS_KSAT, build_args_rate_FMS, build_empty_graph

end

