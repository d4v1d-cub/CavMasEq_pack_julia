# Transforms a Matrix into an Array of Arrays
function slicematrix(A::AbstractMatrix)
    return [c[:] for c in eachcol(A)]
end


# Transforms a vector of bits into an integer
function digits2int(digits::Vector{Int64})
    sum(digits .* 2 .^(0:length(digits)-1))
end
