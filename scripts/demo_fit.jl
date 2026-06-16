include(joinpath(@__DIR__, "setup.jl"))

lower = [
    0.01, 0.01, 0.55, 0.5, 0.05,
    0.5, 0.01, 0.2, 0.5
]

upper = [
    0.1, 0.05, 0.95, 4.0, 0.5,
    4.0, 0.06, 0.25, 1.0
]

initial = (lower + upper) / 2.0

algorithm = SAMIN(
    nt = 5,
    ns = 5,
    rt = 0.95,
    neps = 10,
    f_tol = 0.5,
    x_tol = 1e-3,
    coverage_ok = true,
    verbosity = 0
)

options = Optim.Options(
    iterations = 25,
    f_calls_limit = 500,
    f_tol = 0.5,
    x_tol = 1e-3,
    store_trace = true,
    extended_trace = true,
    show_trace = true
)

const fit_seed = 1234


function objective(x)
    rng = MersenneTwister(fit_seed)
    return loss(
        x,
        15, #scale the base system size [1200, 12000] ([stem , terminal]) by 15x
        250, #number of simulations to average over
        1,  #number of niches, leave at 1
        0.025; #standard deviation for the variability added to k
        rng = rng
    )
end

result = optimize(
    objective,
    lower,
    upper,
    initial,
    algorithm,
    options
)

learned_parameters = Optim.minimizer(result)


println("Learned parameters:")
println(learned_parameters)
println("Learned loss: ", Optim.minimum(result))

figure_rng = MersenneTwister(5678)

simulation_times = [30, 60, 90, 120, 150, 180, 210, 240]
plot_times = [0; simulation_times]

println("Generating survival curves from learned parameters...")
flush(stdout)

survival_simulation = cml_survival_curve(
    simulation_times,
    [1200, 12000],
    250,                 # figure simulations
    1,
    learned_parameters,
    15,
    0.025;
    rng = figure_rng
)

simulated_survival = (
    low  = [1.0; survival_simulation[1][1]],
    mid  = [1.0; survival_simulation[2][1]],
    high = [1.0; survival_simulation[3][1]]
)

observed_survival = (
    low = [
        1.0, 1.0, 0.979031, 0.956522, 0.93453,
        0.913043, 0.892051, 0.871542, 0.851503
    ],

    mid = [
        1.0, 1.0, 0.827263, 0.572362, 0.396003,
        0.273984, 0.189563, 0.131154, 0.0907419
    ],

    high = [
        1.0, 0.900899, 0.574252, 0.36604, 0.233322,
        0.148724, 0.0947999, 0.0604274, 0.0385177
    ]
)

plot_survival_curves(
    plot_times,
    observed_survival,
    simulated_survival;
    savepath = joinpath(
        OUTPUT_DIR,
        "demo_fit_survival_curves.pdf"
    )
)

println("Demo survival plot saved.")

chimerism_rng = MersenneTwister(6789)
plot_times = [10*i for i in 1:27]
chimerism_result = prog_nonprog(
    [1, 270],
    learned_parameters,
    1,
    250,
    15,
    0.0,
    0.025;
    rng = chimerism_rng
)


low = chimerism_result[2][1]
low_prog = mean(chimerism([low[1], low[3]]))
low_nonprog = mean(chimerism([low[2], low[4]]))

mid = chimerism_result[2][2]
mid_prog = mean(chimerism([mid[1], mid[3]]))
mid_nonprog = mean(chimerism([mid[2], mid[4]]))

high = chimerism_result[2][3]
high_prog = mean(chimerism([high[1], high[3]]))
high_nonprog = mean(chimerism([high[2], high[4]]))


prog = (
    low = low_prog,
    mid = mid_prog,
    high = high_prog
)

nonprog = (
    low = low_nonprog,
    mid = mid_nonprog,
    high = high_nonprog
)



plot_chimerism_curves(
    plot_times,
    prog,
    nonprog;
    savepath = joinpath(OUTPUT_DIR, "demo_fit_chimerism_curves.pdf")
    )


    open(joinpath(OUTPUT_DIR, "demo_fit_parameters.txt"), "w") do io
        println(io, "Demo-fit parameters")
        println(io, "Minimum loss: ", Optim.minimum(result))
        println(io)

        for (i, value) in enumerate(learned_parameters)
            println(io, "parameter_$i = $value")
        end
    end
