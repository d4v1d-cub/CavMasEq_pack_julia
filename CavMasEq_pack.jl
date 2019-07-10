using MyModule

include("testing_julia.jl")

Main.MyModule.test(2.0)

using Main.MyModule

test(2.09)
