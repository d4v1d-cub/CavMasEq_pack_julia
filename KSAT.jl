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
function comp_pu_KSAT(p_cav::Array{Float64, 3}, graph::HGraph, ch_u_cond::Array{Int64, 3})
    pu = zeros(Float64, (graph.M, graph.K))
    for he in 1:graph.M
        for i in 1:graph.K
            pu[he, i] = p_cav[he, i, ch_u_cond[he, i]]
        end
    end
    return pu
end

# This function takes a variable node 'node' and a clause 'he' containing that node and
# constructs a list with all the conditional probabilities 'pu' corresponding to the other
# clauses that contain 'node'. The clause 'he' is stored in var_2_he_loc[he_index]
function construct_pu_neighs(node::Int64, he_index::Int64, var_2_he_loc::Vector{Float64}, 
                             pu::Matrix{Float64}, nodes_in::Vector{Dict{Int64, Inf64}})
    other_he = var_2_he_loc[1:end .!= he_index]
    places_in = map(x -> get(x, node, "Node not found in he"), g.nodes_in[other_he])
    indexes = [i, j in other_he, places_in]
    pu_neighs = pu[other_he, ]
end

g = build_ER_HGraph(100, 3, 3)
other_he = g.var_2_he[2][1:end .!= 1]
g.nodes_in[other_he]
places_in = map(x -> get(x, 2, "Node not found in he"), g.nodes_in[other_he])
indexes = [[i, j] for (i, j) in (other_he, places_in)]

function sum_product_KSAT(pu_lp::Vector{Int64}, pu_lm::Vector{Int64})
   println("a") 
end