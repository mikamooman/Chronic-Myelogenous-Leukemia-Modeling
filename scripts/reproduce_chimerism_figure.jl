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

simulation_times = [1, 360]
plot_times = [10*i for i in 1:36]

simulation = prog_nonprog(
    simulation_times,
    parameters,
    1,      #number of niches leave at 1
    500,    #number of trials
    15,     #scale system size
    0.0,    #degree of asymmetric division events
    0.025;  #standard deviation of the variability added to k
    rng = rng
)


low = simulation[2][1]
low_prog = mean(chimerism([low[1], low[3]]))
low_nonprog = mean(chimerism([low[2], low[4]]))

mid = simulation[2][2]
mid_prog = mean(chimerism([mid[1], mid[3]]))
mid_nonprog = mean(chimerism([mid[2], mid[4]]))

high = simulation[2][3]
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
    savepath = joinpath(OUTPUT_DIR, "chimerism_curves.pdf")
    )
