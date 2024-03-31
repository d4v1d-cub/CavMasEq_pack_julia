using Random
using OrdinaryDiffEq, ProgressLogging

include("./general.jl")
include("./graphs.jl")
include("./rates.jl")
include("./KSAT.jl")

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
    graph, all_lp, all_lm, links, ch_u, ch_u_cond, rfunc, rarg_cst, rarg_build, len_cav = p
    # These are the parameters of the integration
    p_cav = reshape(u[1:len_cav], (graph.M, graph.K, 2, graph.nchains_he ÷ 2))
    probi .= u[len_cav + 1:len_cav + 1 + graph.N]
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


# This function integrates the CME's equations for a specific algorithm (given by ratefunc)
# and some boolean formula (given by graph and links)
function CME_KSAT(graph::HGraph, links::Matrix{Int8}, p0::Float64, ratefunc::Function, 
                  rargs_cst, rarg_build::Function, method::Function, tspan::Vector{Float64}, 
                  t_save::Vector{Float64})
    all_lp, all_lm = all_lpm(graph, links)
    ch_u, ch_u_cond = unsat_ch(graph, links)
    p_cav = init_p_cav(graph, p0)
    probi = fill(p0, g1.N)
    len_cav = length(p_cav)
    params = graph, all_lp, all_lm, links, ch_u, ch_u_cond, ratefunc, rargs_cst, rarg_build, len_cav
    u0 = vcat(reshape(p_cav, len_cav), probi)
    prob = ODEProblem(fder_KSAT, u0, tspan, params)
    sol = solve(prob, method(), progress=true, saveat=t_save)
    return sol
end


n = 100
c = 3
K = 3
p0 = 0.5

g1 = build_ER_HGraph(n, c, K, 1)
all_l = gen_links(g1, 1)
all_lp, all_lm = all_lpm(g1, all_l)

ch_u, ch_u_cond = unsat_ch(g1, all_l)
p_cav = init_p_cav(g1, p0)
pu = comp_pu_KSAT(p_cav, g1, ch_u_cond)

p_i = fill(p0, g1.N)

rf = rate_FMS_KSAT
rargs = [1.0, 1.0, g1.K]
d_pc, d_pi = all_ders_CME_KSAT(p_cav, p_i, pu, g1, all_lp, all_lm, rf, 
                               rargs, all_l, ch_u_cond)
