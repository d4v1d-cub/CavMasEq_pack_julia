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
    pu_lpm = Vector{Float64}()
    for (i, j) in iter
        push!(pu_lpm, pu[i, j])
    end
    return pu_lpm
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
# The rates depend on the energy for the current value of the variable (Eu)
# and on the energy after flipping the variable (Es)
# pu_lu are the probabilities for unsatisfied links
# pu_ls are the probabilities for satisfied links
# Ea is the state of the origin clause
# Ea_flip is the state when the spin flips
function sum_product_KSAT(pu_lu::Vector{Float64}, pu_ls::Vector{Float64}, ratefunc::Function, 
                          rate_args, Ea::Int, Ea_flip::Int)
    cu = length(pu_lu)
    cs = length(pu_ls)
    fE_u = recursive_marginal(pu_lu, cu, 1, [1.0])  # This corresponds to links
    fE_s = recursive_marginal(pu_ls, cs, 1, [1.0])
    cumul = 0.0
    for Eu in 0:cu
        for Es in 0:cs
            cumul += ratefunc(Eu + Ea, Es + Ea_flip, rate_args...) * fE_u[Eu + 1] * fE_s[Es + 1]
        end
    end
    return cumul
end


# This function computes the sum-product in the CME for K-SAT.
# The rates depend only on the energy for the current value of the variable (Eu)
function sum_product_KSAT(pu_lu::Vector{Float64}, ratefunc::Function, 
                          rate_args, Ea::Int)
    cu = length(pu_lu)
    fE_u = recursive_marginal(pu_lu, cu, 1, [1.0])
    cumul = 0.0
    for Eu in 0:cu
        cumul += ratefunc(Eu + Ea, rate_args...) * fE_u[Eu]
    end
    return cumul
end


# This function gives all sums needed in the computation of the derivative of the pcav
# related to he and node.
# It returns an array all_sums[val, Ea, Ea_flip]
# When val=1 (si = 1) the array of probabilities pu_lu is taken from the negative links
# because those are the unsatisfied links. At the same time, pu_ls is taken from the positive
# links, which will become unsatisfied when si flips
# When val=2 the roles are inverted.
# Ea is the value of the clause 'a' inside pcav
# Ea_flip is the value when si flips
function compute_all_sums(pu::Matrix{Float64},
                    node::Int64, he::Int64, graph::HGraph, 
                    all_lp::Vector{Vector{Vector{Int64}}}, 
                    all_lm::Vector{Vector{Vector{Int64}}}, ratefunc::Function, 
                    rate_args)
    all_sums = zeros(Float64, (2, 2, 2))
    for val in 1:2
        if val == 1
            pu_lu, pu_ls = construct_pu_neighs(node, he, pu, graph.nodes_in, all_lm, all_lp)
        else
            pu_lu, pu_ls = construct_pu_neighs(node, he, pu, graph.nodes_in, all_lp, all_lm)
        end
        all_sums[val, 1, 1] = sum_product_KSAT(pu_lu, pu_ls, ratefunc, rate_args, 0, 0)
        all_sums[val, 2, 1] = sum_product_KSAT(pu_lu, pu_ls, ratefunc, rate_args, 1, 0)
        all_sums[val, 1, 2] = sum_product_KSAT(pu_lu, pu_ls, ratefunc, rate_args, 0, 1)
        # The triad (val, 1, 1) it is not possible.
    end    
    
    return all_sums
end


# This function computes one term in the derivative of pcav. It corresponds to
# the flipping of the variable indexed as 'index_in' among the variables in the argument
# of p_cav. The array 'all_sums' is computed using the function 'compute_all_sums'
# ch_exc is an integer that codes the combination of the variables in the argument of the
# probability 'p_cav'. Remember that all p_cav are conditioned to one already unsatisfied link.
function der_pcav_node_KSAT(p_cav::Vector{Float64}, all_sums::Array{Float64, 3}, 
                            ch_exc::Int64, ch_exc_unsat::Int64, 
                            index_in::Int64)
    ch_exc_flip = (ch_exc ⊻ (2 ^ (index_in - 1)))  # The ⊻ (xor) operation flips the variable
    val = ((ch_exc >> (index_in - 1)) & 1)         # Takes the value of the variable
    Ea = (ch_exc == ch_exc_unsat)                  # As the variable in the conditional is set 
    Ea_flip = (ch_exc_flip == ch_exc_unsat)        # to always unsatisfy its link, the clause
    # will be unsatisfied only if the rest of the variables are in their unsatisfied configuration
    
    return -all_sums[val + 1, Ea + 1, Ea_flip + 1] * p_cav[ch_exc + 1] + 
            all_sums[2 - val, Ea_flip + 1, Ea + 1] * p_cav[ch_exc_flip + 1]
end


# This function computes the derivatives of all the p_cav that correspond to a clause
function der_pcav_KSAT(p_cav::Matrix{Float64}, pu::Matrix{Float64}, he::Int64, graph::HGraph, 
    all_lp::Vector{Vector{Vector{Int64}}}, 
    all_lm::Vector{Vector{Vector{Int64}}}, ratefunc::Function, 
    rate_args, links::Matrix{Int8})
    
    nch = graph.chains_he ÷ 2
    ders = zeros(Float64, (graph.K, nch))
    for i in 1:graph.K
        lvec_rev = 1 .- links[he, 1:end .!= i]  # gives the configuration that unsatisfies every link
        ch_exc_unsat = digits2int(lvec_rev)    # converting it to an integer
        for j in 1:graph.K-1
            node = graph.nodes_except[he,i,j]
            all_sums = compute_all_sums(pu, node, he, graph, all_lp, all_lm, ratefunc, rate_args)
            for ch_exc in 0:nch-1
                ders[i, ch_exc + 1] += der_pcav_node_KSAT(p_cav[i, :], all_sums, ch_exc, 
                                                         ch_exc_unsat, j)   
            end
        end
    end
    return ders
end

n = 100
c = 3
K = 3
g1 = build_ER_HGraph(n, c, K)
all_l = gen_links(g1)
all_lp, all_lm = all_lpm(g1, all_l)

ch_u_cond = unsat_ch(g1, all_l)
p_cav = init_p_cav_KSAT(g1, 0.5)
pu = comp_pu_KSAT(p_cav, g1, ch_u_cond)
g1.K

include("rates.jl")
der_pcav_KSAT(p_cav[1, :, :], pu, 1, g1, all_lp, all_lm, rate_FMS_KSAT, [1.0, 1.0, g1.K], all_l)

