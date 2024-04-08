# This script contains some general functions that will be used in several other modules

# Transforms a Matrix into an Array of Arrays
function slicematrix(A::AbstractMatrix)
    return [c[:] for c in eachcol(A)]
end


# Transforms a vector of bits into an integer
function digits2int(digits::Vector{Int64})
    sum(digits .* 2 .^(0:length(digits)-1))
end



# This function recursively computes a vector fE[k], with k=1,..., c+1
# fE[k] is the sum of 'c' binary variable constrained to sum exactly k - 1 and weighted
# with a factorized distribution pu[i], with i = 1,..., c
function recursive_marginal(pu::Vector{Float64}, c::Int64, k::Int64, fE::Vector{Float64})
    if k <= c
        fEnew = Vector{Float64}(undef, k + 1)
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


# Inserts a bit ('val' - 1) in the position 'place' of the number 'ch_exc'
function convert_chain(ch_exc::Int64, val::Int64, place::Int64)
    return (val - 1) * 2 ^ (place - 1) + ch_exc % (2 ^ (place - 1)) + 
           2 ^ place * (ch_exc ÷ (2  ^ (place - 1)))
end


# It normalizes a vector with the norm 1
function normalize(vec::Vector{Float64})
    n = sum(vec)
    vec ./= n
    return vec
end