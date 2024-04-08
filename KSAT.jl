# This script contains all the functions needed to compute the derivatives 
# in the integration of the CME on the K-SAT model

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


# This function generates two lists "lp" and "lm" for a given clause "a" and some variable "i"
# inside the clause
# lp[i] contains all the clauses "b != a" that are positively linked with the var "i" (l_i^{b} = 1)
# lm[i] contains the negatively linked clauses (l_i^{b} = -1)
function get_lpm(node::Int64, he_source::Int64, var_2_he_loc::Vector{Int64}, 
                 nodes_in::Vector{Dict{Int64, Int64}}, links::Matrix{Int8})
    other_he = filter(x -> x != he_source, var_2_he_loc)
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
# ch_u_cond[he, v] is an array that contains, for each hyperedge 'he' and each
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
    return ch_u, ch_u_cond
end



# It computes the conditional cavity probabilities of having partially unsatisfied hyperedges.
# ch_u_cond[he, v] is an array that contains, for each hyperedge 'he' and each
# variable 'v' in that hyperedge, the chains (coded as integers) of the rest 
# of the variables in the hyperedge that unsatisfy their links.
# pu[he, i, s] contains, for each hyperedge 'he' and node inside that hyperedge 'i',
# the cavity conditional probability of having the rest of the nodes (he \ i) in the hyperedge in their
# unsat configuration, given that 'i' is satisfying (s=1) or unsatisfying (s=2) its link
function comp_pu_KSAT(p_cav::Array{Float64, 4}, graph::HGraph, ch_u_cond::Matrix{Int64})
    pu = zeros(Float64, (graph.M, graph.K, 2))
    for he in 1:graph.M
        for i in 1:graph.K
            for s in 1:2
                pu[he, i, s] = p_cav[he, i, s, ch_u_cond[he, i] + 1]
            end
        end
    end
    return pu
end


# This function constructs the list of cavity conditional probabilities 'pu' corresponding
# to 'node' and the clauses 'he_list'. 
# The conditional in 'pu' is set from outside this function. In this way, 'pu_lpm' is already
# conditioned to a satisfied or an unsatisfied link
function get_pu_lpm(node::Int64, he_list::Vector{Int64}, nodes_in::Vector{Dict{Int64, Int64}}, 
                    pu::Matrix{Float64})
    places_in = map(x -> get(x, node, "Error: Node not found in he"), nodes_in[he_list])
    # these are the places where we can find 'node' in the hyperedges 'he_list'
    iter = zip(he_list, places_in)
    pu_lpm = Vector{Float64}(undef, length(he_list))
    counter = 1
    for (he, j) in iter
        pu_lpm[counter] = pu[he, j]
        counter += 1
    end
    return pu_lpm
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
    fE_u_given_u = recursive_marginal(pu_lu, cu, 1, [1.0])  # The currently unsat have an 
    # unsat condition (s = 2) and come from the clauses that are unsatisfied by the current value
    # of the variable
    fE_u_given_s = recursive_marginal(pu_ls, cs, 1, [1.0])  # The clauses that will be unsat
    # after flipping are among the rest of the clauses. Their probabilities are conditioned on a 
    # satisfied link
    cumul = 0.0
    for Eu in 0:cu
        for Es in 0:cs
            cumul += ratefunc(Eu + Ea, Es + Ea_flip, rate_args...) * 
                     fE_u_given_u[Eu + 1] * fE_u_given_s[Es + 1]
        end
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
function compute_all_sums(pu::Array{Float64, 3},
                    node::Int64, he::Int64, graph::HGraph, 
                    all_lp::Vector{Vector{Vector{Int64}}}, 
                    all_lm::Vector{Vector{Vector{Int64}}}, ratefunc::Function, 
                    rate_args)
    all_sums = zeros(Float64, (2, 2, 2))
    for val in 1:2
        if val == 1
            # When the variable is 1 (val = 1), the unsatisfied links are the negative ones
            he_list = all_lm[he][graph.nodes_in[he][node]]
            pu_lu = get_pu_lpm(node, he_list, graph.nodes_in, pu[:, :, 2])
            he_list = all_lp[he][graph.nodes_in[he][node]]
            pu_ls = get_pu_lpm(node, he_list, graph.nodes_in, pu[:, :, 1])
        else
            # When the variable is -1 (val = 2), the unsatisfied links are the positive ones
            he_list = all_lp[he][graph.nodes_in[he][node]]
            pu_lu = get_pu_lpm(node, he_list, graph.nodes_in, pu[:, :, 2])
            he_list = all_lm[he][graph.nodes_in[he][node]]
            pu_ls = get_pu_lpm(node, he_list, graph.nodes_in, pu[:, :, 1])
        end
        all_sums[val, 1, 1] = sum_product_KSAT(pu_lu, pu_ls, ratefunc, rate_args, 0, 0)
        all_sums[val, 2, 1] = sum_product_KSAT(pu_lu, pu_ls, ratefunc, rate_args, 1, 0)
        all_sums[val, 1, 2] = sum_product_KSAT(pu_lu, pu_ls, ratefunc, rate_args, 0, 1)
        # The triad (val, 2, 2) it is not possible.
    end    
    
    return all_sums
end