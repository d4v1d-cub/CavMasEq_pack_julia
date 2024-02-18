using Random
include("graphs.jl")
include("tools.jl")

# It generates the binary links of K-SAT's boolean formula.
function gen_links(graph::HGraph, idum::Int64=rand(1:typemax(Int64)))
    Random.seed!(idum)
    links = zeros(Int8, (graph.M, graph.K))
    for he in 1:graph.M
        for i in 1:graph.K
            links[he, i] = rand([0, 1])
        end
    end
    return links
end


# This function generates two lists "lp" and "lm" for a given clause and some variable
# inside the clause
# lp contains all the clauses "a" that are positively linked with the var "i" (l_i^{a} = 1)
# lm[var] contains the negatively linked clauses (l_i^{a} = -1)
function get_lpm(node::Int64, he_index::Int64, var_2_he_loc::Vector{Int64}, 
                 nodes_in::Vector{Dict{Int64, Int64}}, links::Matrix{Int8})
    other_he = var_2_he_loc[1:end .!= he_index]
    lp = Vector{Int64}()
    lm = Vector{Int64}()
    for he in other_he
        place_in = nodes_in[he][node]
        l = links[he, place_in]
        if l == 0
            push!(lp, he)
        elseif l == 1
            push!(lm, he)
        else
            println("Error: links[" * string(he) * ", " * string(place_in) * "] is not 1 or -1")
            println("links[" * string(he) * ", " * string(place_in) * "] = " * string(l))
        end
    end
    return lp, lm
end



# This function builds the list with all lp and lm (see function get_lpm)
function all_lpm(graph::HGraph, links::Matrix{Int8})
    all_lp = [[Vector{Int64}() for _ in 1:graph.K] for _ in 1:graph.M]
    all_lm = [[Vector{Int64}() for _ in 1:graph.K] for _ in 1:graph.M]
    for he in 1:graph.M
        for i in 1:graph.K
            node = graph.he_2_var[he, i]
            lp, lm = get_lpm(node, he, graph.var_2_he[node], graph.nodes_in, links)
            append!(all_lp[he][i], lp)
            append!(all_lm[he][i], lm)
        end
    end
    return all_lp, all_lm
end


# It will be useful to have pre-computed lists of the combinations that
# unsatisfy each clause (ch_u).
# ch_u_cond[he, v, ch] is an array that contains, for each hyperedge 'he' and each
# variable 'v' in that hyperedge, the chains (coded as integers) of the rest 
# of the variables in the hyperedge that unsatisfy their links.
function unsat_ch(graph::HGraph, links::Matrix{Int8})
    ch_u = zeros(Int64, graph.M)
    ch_u_cond = zeros(Int64, (graph.M, graph.K))
    for he in 1:graph.M
        ch_u[he] = digits2int(map(x -> x ⊻ 1, links[he, :]))
        for i in 1:graph.K
            ch_u_cond[he, i] = digits2int(map(x -> x ⊻ 1, links[he, 1:end .!= i]))
        end
    end
    return ch_u_cond
end


# This function initializes the cavity conditional probabilities with a fully 
# independent initial distribution given by the vector 'p0'
# This cavity probabilities are conditioned to having an already unsatisfied link
# because those are the only probabilities that are actually needed in KSAT.
function init_p_cav_KSAT(graph::HGraph, p0::Vector{Float64})
    ch_exc = graph.chains_he ÷ 2
    p_cav = zeros(Float64, (graph.M, graph.K, ch_exc))
    for he in 1:graph.M
        for i in 1:graph.K
            for ch in 0:ch_exc - 1
                bits = digits(ch, base=2, pad=graph.K-1)
                p_cav[he, i, ch + 1] = prod(bits + (1 .- 2 * bits) .* p0[graph.nodes_except[he, i, :]])
            end
        end
    end
    return p_cav
end


# The user can just pass a single float 'p0' and the vector of initial conditions 
# is assumed to be homogeneous
function init_p_cav_KSAT(graph::HGraph, p0::Float64)
    p0 = fill(p0, graph.N)
    init_p_cav_KSAT(graph, p0)
end


# It computes the conditional cavity probabilities of having unsatisfied hyperedges.
# ch_u_cond[he, v, ch] is an array that contains, for each hyperedge 'he' and each
# variable 'v' in that hyperedge, the chains (coded as integers) of the rest 
# of the variables in the hyperedge that unsatisfy their links.
function comp_pu_KSAT(p_cav::Array{Float64, 3}, graph::HGraph, ch_u_cond::Matrix{Int64})
    pu = zeros(Float64, (graph.M, graph.K))
    for he in 1:graph.M
        for i in 1:graph.K
            pu[he, i] = p_cav[he, i, ch_u_cond[he, i] + 1]
        end
    end
    return pu
end


# This function constructs the list of cavity conditional probabilities 'pu' corresponding
# to 'node' and the clauses 'he_list'
function get_pu_lpm(node::Int64, he_list::Vector{Int64}, nodes_in::Vector{Dict{Int64, Int64}}, 
                    pu::Matrix{Float64})
    places_in = map(x -> get(x, node, "Error: Node not found in he"), nodes_in[he_list])
    iter = zip(he_list, places_in)
    return [pu[i, j] for (i, j) in iter]
end


# This function takes a variable node 'node' and a clause 'he' containing that node and
# constructs two lists with all the conditional probabilities 'pu' corresponding to the other
# clauses that contain 'node'.
# The list pu_lp contains all the probabilities that correspond to positive links
# The list pu_lm contains all the probabilities that correspond to negative links
function construct_pu_neighs(node::Int64, he::Int64, pu::Matrix{Float64}, 
                             nodes_in::Vector{Dict{Int64, Int64}}, 
                             all_lp::Vector{Vector{Vector{Int64}}}, 
                             all_lm::Vector{Vector{Vector{Int64}}})
    he_list = all_lp[he][nodes_in[he][node]]
    pu_lp = get_pu_lpm(node, he_list, nodes_in, pu)
    he_list = all_lm[he][nodes_in[he][node]]
    pu_lm = get_pu_lpm(node, he_list, nodes_in, pu)
    return pu_lp, pu_lm
end


# n = 100
# c = 3
# K = 3
# g1 = build_ER_HGraph(n, c, K)
# all_l = gen_links(g1)
# all_lp, all_lm = all_lpm(g1, all_l)

# ch_u_cond = unsat_ch(g1, all_l)
# p_cav = init_p_cav_KSAT(g1, 0.5)
# pu = comp_pu_KSAT(p_cav, g1, ch_u_cond)

# construct_pu_neighs(g1.he_2_var[1, 1], 1, pu, g1.nodes_in, all_lp, all_lm)


# This function recursively computes a vector fE[k], with k=1,..., c+1
# fE[k] is the sum of 'c' binary variable constrained to sum exactly k - 1 and weighted
# with a factorized distribution pu[i], with i = 1,..., c
function recursive_marginal(pu::Vector{Float64}, c::Int64, k::Int64, fE::Vector{Float64})
    if k <= c
        fEnew = zeros(Float64, k + 1)
        fEnew[1] = (1 - pu[k]) * fE[1]
        for i in 1:k - 1
            fEnew[i + 1] = (1 - pu[k]) * fE[i + 1] + pu[k] * fE[i]
        end
        fEnew[k + 1] = pu[k] * fE[k]
        return recursive_marginal(pu, c, k + 1, fEnew)
    else
        return fE
    end
end


# This function computes the sum-product in the CME for K-SAT.
# The rates depend on the energy for the current value of the variable (Ep)
# and on the energy after flipping the variable (Em)
function sum_product_KSAT(pu_lp::Vector{Float64}, pu_lm::Vector{Float64}, ratefunc::Function, 
                          rate_args::Vector{Any}, Eap::Int8, Eam::Int8)
    cp = length(pu_lp)
    cm = length(pu_lm)
    fEp = recursive_marginal(pu_lp, cp, 1, [1.0])
    fEm = recursive_marginal(pu_lm, cm, 1, [1.0])
    s = 0.0
    for Ep in 0:cp
        for Em in 0:cm
            s += ratefunc(Ep + Eap, Em + Eam, rate_args...) * fEp[Ep] * fEm[Em]
        end
    end
    return s
end


# This function computes the sum-product in the CME for K-SAT.
# The rates depend only on the energy for the current value of the variable (Ep)
function sum_product_KSAT(pu_lp::Vector{Float64}, ratefunc::Function, 
                          rate_args::Vector{Any}, Eap::Int8)
    cp = length(pu_lp)
    fEp = recursive_marginal(pu_lp, cp, 1, [1.0])
    s = 0.0
    for Ep in 0:cp
        s += ratefunc(Ep + Eap, rate_args...) * fEp[Ep]
    end
    return s
end