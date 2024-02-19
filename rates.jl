include("graphs.jl")

# This is the rate of FMS algorithm for KSAT
function rate_FMS_KSAT(Ep::Int64, Em::Int64, T::Float64, avE::Float64, K::Number)
    dE = Em - Ep
    if dE > 0
        return Ep / K / avE * exp(-dE / T)
    else
        return Ep / K / avE
    end
end