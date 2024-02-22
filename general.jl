# Transforms a Matrix into an Array of Arrays
function slicematrix(A::AbstractMatrix)
    return [c[:] for c in eachcol(A)]
end


# Transforms a vector of bits into an integer
function digits2int(digits::Vector{Int64})
    sum(digits .* 2 .^(0:length(digits)-1))
end


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