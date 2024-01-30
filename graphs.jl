using Random

include("tools.jl")

mutable struct HGraph
# This structure holds the hypergraph information
    N::Int64                                   # number of variable nodes
    M::Int64                                   # number of hyperedges
    K::Int64                                   # number of nodes in each hyperedge
    var_2_he::Array{Array{Int64, 1}, 1}        # list of hyperedges for each variable
    he_2_var::Matrix{Int64}                    # list of variables per hyperedge
    degrees::Vector{Int64}                     # list of variable nodes' degrees
    nchains::Vector{Int64}                     # list equal to 2.^{degrees}
    nodes_in::Array{Dict{Int64, Int64}, 1}               # dictionary with key=node_index and val=place_in_he   
end

# Builds a hypergraph from the lists of hyperedges per variable (var_2_he),
# variables per hyperedge (he_2_var) and the list of variable nodes' degrees
function build_HGraph(var_2_he::Array{Array{Int64, 1}, 1} , he_2_var::Matrix{Int64}, 
                degrees::Vector{Int64})
    N = length(var_2_he)
    M = size(he_2_var,1)
    K = size(he_2_var,2)
    nchains = 2 .^ degrees
    nodes_in = Array{Dict{Int64, Int64}, 1}()
    for he in 1:M
        nin_he = Dict{Int64, Int64}()
        for i in 1:K
            nin_he[he_2_var[he, i]] = i
        end
        push!(nodes_in, nin_he)
    end

    return HGraph(N, M, K, var_2_he, he_2_var, degrees, nchains, nodes_in)
end


# This function creates a Random Regular Hypergaph with node connectivity 'c' and factor 
# node connectivity 'K'. The random seed is controlled by the user
function RRHyperGraph(N::Int64, c::Int64, K::Int64, idum::Int64=1)
    Random.seed!(idum)
    M = N * c / K 
    if isinteger(M)
        M = Int64(M)
        he_2_var = zeros(Int64, (M, K))
        var_2_he = zeros(Int64, (N, c))
        degrees = zeros(Int64, N)
        copynodes = repeat(1:N, c)  # We create an auxiliary array with 'c' copies of each node
        for he in 1:M
            for i in 1:K
                place = rand(1:length(copynodes))
                he_2_var[he, i] = copynodes[place]
                var_2_he[copynodes[place], degrees[copynodes[place]] + 1] = he
                degrees[copynodes[place]] += 1
                deleteat!(copynodes, place)
                # By randomly selecting the nodes in each hyperedge from the auxiliary array
                # we make sure that each node is in exactly 'c' hyperedges 
            end
        end
        return var_2_he, he_2_var, degrees        
    else
        println("The number of factor nodes 'N * c / K' must be an integer")
        return nothing
    end
end


# Builds a random regular hypergraph with parameters N, c, K and random seed 'idum'
function build_RR_HGraph(N::Int64, c::Int64, K::Int64, idum::Int64=1)
    var_2_he, he_2_var, degrees = RRHyperGraph(N, c, K, idum)
    var_2_he = slicematrix(var_2_he)             # The matrix is casted into an array of arrays
    return build_HGraph(var_2_he, he_2_var, degrees)
end


# This function creates an Erdos-Renyi Hypergaph with mean node connectivity 'c' and factor 
# node connectivity 'K'. The random seed is controlled by the user
function ERHyperGraph(N::Int64, c::Int64, K::Int64, idum::Int64=1)
    Random.seed!(idum)
    M = N * c / K 
    if isinteger(M)
        M = Int64(M)
        he_2_var = zeros(Int64, (M, K))
        var_2_he = fill(Array{Int64, 1}(), N)
        degrees = zeros(Int64, N)
        for he in 1:M
            for i in 1:K
                node = rand(1:N)
                he_2_var[he, i] = node
                push!(var_2_he[node], he)
                degrees[node] += 1
                # Each hyperedge is formed by 'K' variables chosen uniformly at random
                # This leads to a Poisson distribution of conectivities (ER Graph) 
            end
        end
        return var_2_he, he_2_var, degrees  
    else
        println("The number of factor nodes 'N * c / K' must be an integer")
        return nothing
    end      
end


# Builds a random regular hypergraph with parameters N, c, K and random seed 'idum'
function build_ER_HGraph(N::Int64, c::Int64, K::Int64, idum::Int64=1)
    var_2_he, he_2_var, degrees = ERHyperGraph(N, c, K, idum)
    return build_HGraph(var_2_he, he_2_var, degrees)
end