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


# This function computes one term in the derivative of pcav. It corresponds to
# the flipping of the variable indexed as 'index_in' among the variables in the argument
# of p_cav. The array 'all_sums' is computed using the function 'compute_all_sums'
# ch_exc is an integer that codes the combination of the variables in the argument of the
# probability 'p_cav'. In this case all p_cav are conditioned to one already unsatisfied link.
function der_pcav_contr_node_unsat(p_cav::Vector{Float64}, all_sums::Array{Float64, 3}, 
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



# This function computes one term in the derivative of pcav. It corresponds to
# the flipping of the variable indexed as 'index_in' among the variables in the argument
# of p_cav. The array 'all_sums' is computed using the function 'compute_all_sums'
# ch_exc is an integer that codes the combination of the variables in the argument of the
# probability 'p_cav'. In this case all p_cav are conditioned to one satisfied link.
function der_pcav_contr_node_sat(p_cav::Vector{Float64}, all_sums::Array{Float64, 3}, 
            ch_exc::Int64, index_in::Int64)
ch_exc_flip = (ch_exc ⊻ (2 ^ (index_in - 1)))  # The ⊻ (xor) operation flips the variable
val = ((ch_exc >> (index_in - 1)) & 1)         # Takes the value of the variable

return -all_sums[val + 1, 1, 1] * p_cav[ch_exc + 1] + 
    all_sums[2 - val, 1, 1] * p_cav[ch_exc_flip + 1]
end



# This function computes the part of the derivatives of the p_cav that correspond to
# the flipping of variable 'node' inside a clause 'he' conditioned on another node that is
# defined in the function all_ders_node_KSAT (see below)
# The result is cumulated in the matrix 'ders'
function der_pcav_KSAT(p_cav::Matrix{Float64}, pu::Array{Float64, 3}, he::Int64, 
    node::Int64, place_node::Int64, ch_exc_unsat::Int64, graph::HGraph, 
    all_lp::Vector{Vector{Vector{Int64}}}, 
    all_lm::Vector{Vector{Vector{Int64}}}, 
    ratefunc::Function, rate_args, ders::Matrix{Float64}, nch::Int64)
    
    all_sums = compute_all_sums(pu, node, he, graph, all_lp, all_lm, ratefunc, rate_args)
    for ch_exc in 0:nch-1
        ders[1, ch_exc + 1] += der_pcav_contr_node_sat(p_cav[1, :], all_sums, ch_exc, place_node)
        ders[2, ch_exc + 1] += der_pcav_contr_node_unsat(p_cav[2, :], all_sums, ch_exc, ch_exc_unsat, 
                                                         place_node)   
    end
    return all_sums, ders
    # It returns the sums corresponding to flipping node with fixed values of 'he' and the 
    # computed derivatives
end


# Given the sum-products corresponding to the flipping variable 'i', summed over all the clauses
# containing 'i' different from a given 'he', this function computes the derivative of the local
# probability 'probi'. It receives the cavity conditional probability of having 'he' unsatisfied
# (pu_he) and the value of the link between 'i' and 'he' (li_he)
function der_pi_KSAT(probi::Float64, pu_he::Vector{Float64}, all_sums::Array{Float64, 3}, 
                     li_he::Int8)
    cumul = 0.0
    for u_neigh in 0:1
        pr = 1 - u_neigh - (1 - 2 * u_neigh) * pu_he[li_he + 1]
        pr_flip = 1 - u_neigh - (1 - 2 * u_neigh) * pu_he[2 - li_he]
        Ea = u_neigh * (li_he == 1)             # 'u_neigh' controls the state of the rest of the variables in 'he'
        # If u_neigh = 1, the rest of the variables unsatisfy their links and u=0 otherwise 
        # However, the clause 'he' is unsat only is 'i' also unsatisfies its link
        # As 'probi' is defined for si = 0, this happens only if 'li_he=1'
        Ea_flip = u_neigh * (li_he == 0)       # On the other hand, if li_he=0 the clause will be unsatisfied
        # after flipping 'si'
        cumul += -all_sums[1, Ea + 1, Ea_flip + 1] * pr * probi + 
        all_sums[2, Ea_flip + 1, Ea + 1] * pr_flip * (1 - probi)
    end
    
    return cumul
end


# This function computes the all the derivatives related to a node 'node'
function all_ders_node_KSAT(p_cav::Array{Float64, 4}, probi::Float64, 
    pu::Array{Float64, 3}, node::Int64, graph::HGraph, all_lp::Vector{Vector{Vector{Int64}}}, 
    all_lm::Vector{Vector{Vector{Int64}}}, ratefunc::Function, rate_args, 
    links::Matrix{Int8}, ders_pcav::Array{Float64, 4}, nch::Int64, ch_u_cond::Matrix{Int64})
    
    all_sums = zeros(Float64, (2, 2, 2))

    for he in graph.var_2_he[node]           # Computing all the derivatives of cavity probabilities 
        place_in = graph.nodes_in[he][node]  # where 'node' is flipped
        for j in 1:graph.K - 1
            node_neigh = graph.nodes_except[he, place_in, j]
            place_neigh = graph.nodes_in[he][node_neigh]   
            ch_exc_unsat = ch_u_cond[he, place_neigh]  # gives the configuration that unsatisfies every link
                                                       # converting it to an integer
            place_in_exc = graph.place_there[he, place_neigh][node] # Gets the index of 'node' in the 
                                                       # array  graph.nodes_except[he, place_neigh, :]
            all_sums, ders_pcav[he, place_neigh, :, :] = 
                der_pcav_KSAT(p_cav[he, place_neigh, :, :], pu, he, node, place_in_exc, ch_exc_unsat, 
                graph, all_lp, all_lm, ratefunc, rate_args, ders_pcav[he, place_neigh, :, :], nch)   
        end
    end

    if length(graph.var_2_he[node]) > 0
        he = graph.var_2_he[node][end]
        place_in = graph.nodes_in[he][node]
        return der_pi_KSAT(probi, pu[he, place_in, :], all_sums, links[he, place_in])
    else
        return 0
    end
end


# This function computes all the derivatives for the CME in the K-SAT
function all_ders_CME_KSAT(p_cav::Array{Float64, 4}, probi::Vector{Float64}, pu::Array{Float64, 3}, 
    graph::HGraph, all_lp::Vector{Vector{Vector{Int64}}}, all_lm::Vector{Vector{Vector{Int64}}}, 
    ratefunc::Function, rate_args, links::Matrix{Int8}, ch_u_cond::Matrix{Int64})
    
    d_pcav = zeros(Float64, size(p_cav))
    d_node = zeros(Float64, size(probi))
    nch = graph.chains_he ÷ 2
    for node in 1:graph.N
        d_node[node] = all_ders_node_KSAT(p_cav, probi[node], pu, node, graph, all_lp, 
                                          all_lm, ratefunc, rate_args, links, d_pcav, nch, ch_u_cond)
    end
    return d_pcav, d_node
end