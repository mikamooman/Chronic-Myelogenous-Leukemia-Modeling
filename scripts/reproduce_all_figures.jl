project_root = normpath(joinpath(@__DIR__, ".."))

function run_figure_script(filename)
    script_path = joinpath(@__DIR__, filename)

    println()
    println("========================================")
    println("Running: $filename")
    println("========================================")

    run(`$(Base.julia_cmd()) --project=$project_root $script_path`)
end

run_figure_script("reproduce_survival_figure.jl")
run_figure_script("reproduce_chimerism_figure.jl")

println()
println("========================================")
println("All figures generated successfully.")
println("Outputs are located in:")
println(joinpath(project_root, "output"))
println("========================================")
