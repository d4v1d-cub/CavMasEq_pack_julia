function slicematrix(A::AbstractMatrix)
    return [c[:] for c in eachcol(A)]
end