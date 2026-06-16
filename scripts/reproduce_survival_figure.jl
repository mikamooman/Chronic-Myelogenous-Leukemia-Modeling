include(joinpath(@__DIR__, "setup.jl"))

rng = MersenneTwister(1234)

parameters = [
    0.024405393131866588,
    0.01156009633136531,
    0.7196730199636183,
    3.6675406155013843,
    0.24538071790689933,
    0.8990541058262109,
    0.018481240273735107,
    0.23071406853662568,
    0.8511500334069608
]

simulation_times = [30, 60, 90, 120, 150, 180, 210, 240]
plot_times = [0; simulation_times]


simulation = cml_survival_curve(
    simulation_times,
    [1200, 12000],
    1500,
    1,
    parameters,
    15,
    0.025;
    rng = rng
)



simulated = (
    low = [1.0; simulation[1][1]],
    mid = [1.0; simulation[2][1]],
    high = [1.0; simulation[3][1]]
)

observed = (
    low = [1.0, 1.0, 0.979031, 0.956522, 0.93453,
           0.913043, 0.892051, 0.871542, 0.851503],

    mid = [1.0, 1.0, 0.827263, 0.572362, 0.396003,
           0.273984, 0.189563, 0.131154, 0.0907419],

    high = [1.0, 0.900899, 0.574252, 0.36604, 0.233322,
            0.148724, 0.0947999, 0.0604274, 0.0385177]
)




plot_survival_curves(
    plot_times,
    observed,
    simulated;
    savepath = joinpath(OUTPUT_DIR, "cml_free_survival_curves.pdf")
)
