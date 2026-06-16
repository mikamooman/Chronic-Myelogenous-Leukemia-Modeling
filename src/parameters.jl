Base.@kwdef struct ModelParams
    # actual CML chimerism threshold
    threshold::Float64

    # Population-dependent leukemic stem death
    α::Float64
    α2::Float64
    h::Float64
    k::Float64

    # Feedback parameters
    pbar::Float64            # p̄ and p̄ᴸ
    γ::Float64
    γ_ratio::Float64
    ϕL::Float64
    n::Float64

    # Initial chimerism distribution modes
    low_mode::Float64
    mid_mode::Float64
    high_mode::Float64

    # Fixed parameters
    d::Float64 = 0.1
    β::Float64 = 0.03
    δ::Float64 = 1.0
    ϕ::Float64 = 0.0
end

gammaL(θ::ModelParams) = θ.γ_ratio * θ.γ

function theta_from_x(x)
    k_val = (
        (2*x[3] - 1 - x[2]) /
        (1 + x[7] / (2*x[3] - 1 - x[2] - x[7]))
    )

    α2_val = x[7] / (2*x[3] - 1 - x[2] - x[7])

    γ_val = (2*x[3] - 1) / 12000

    return ModelParams(
        threshold = x[8],


        α = x[5],
        α2 = α2_val,
        h = x[6],
        k = k_val,


        pbar = x[3],
        γ = γ_val,
        γ_ratio = x[9],
        ϕL = x[1],
        n = x[4],


        low_mode = 0.029,
        mid_mode = 0.074,
        high_mode = 0.078,

    )
end

function add_variability(θ::ModelParams, σ::Float64; rng=Random.default_rng())
    var = [1 + σ * randn(rng) for _ in 1:10]

    return ModelParams(
        threshold = θ.threshold,

        α = θ.α,
        α2 = θ.α2,
        h = θ.h,
        k = θ.k*var[1],

        pbar = θ.pbar,
        γ = θ.γ,
        γ_ratio = θ.γ_ratio,
        ϕL = θ.ϕL,
        n = θ.n,

        low_mode = θ.low_mode,
        mid_mode = θ.mid_mode,
        high_mode = θ.high_mode,
    )
end
