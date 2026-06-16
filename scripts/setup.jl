using Pkg

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))


Pkg.activate(PROJECT_ROOT)



ENV["GKSwstype"] = "100"

include(joinpath(PROJECT_ROOT, "src", "CML_Model.jl"))



using .CML_Model

using Plots
using Random
using Base.Threads
using Optim
using Distributions



const OUTPUT_DIR = joinpath(PROJECT_ROOT, "output")
mkpath(OUTPUT_DIR)
