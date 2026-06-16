#returns a vector containing the number of simulations which have developed cml at each time step
function simulate_survival_counts(
    θ::ModelParams,
    size::Vector{Int},
    times::Vector{Int},
    tr::Int64,
    niche_num::Int64,
    sc::Int64,
    σ::Float64,
    lower::Float64,
    upper::Float64,
    mode::Float64;
    rng=Random.default_rng()
)
    results1 = zeros(Int, length(times))
    ex1 = zeros(Int, length(times))

    for j in 1:tr
        input = gen_input(
            rand(rng, TriangularDist(lower, upper, mode)),
            size
        )

        temp_res = zeros(Int, length(times))
        temp_ex = zeros(Int, length(times))

        run_length = length(times) + 1

        for _ in 1:niche_num
            mut_cells, ext_cells = variability_cml_free(
                times,
                input,
                θ,
                run_length,
                sc,
                σ;
                rng = rng
            )

            temp_res .= clamp.(temp_res .+ mut_cells, 0, 1)
            temp_ex .= clamp.(temp_ex .+ ext_cells, 0, 1)

            run_length = check_esc(temp_res)

            if run_length == 1
                break
            end
        end

        results1 .+= temp_res
        ex1 .+= temp_ex
    end
    return results1, ex1
end

#function for returning cml free survival curves for each initial condition group
function cml_survival_curve(times, size, trials, niche, x, sc, σ; rng=Random.default_rng())
    θ = theta_from_x(x)

    ii = zeros(length(times))
    ii2 = zeros(length(times))
    iii = zeros(length(times))
    iii2 = zeros(length(times))
    iiii = zeros(length(times))
    iiii2 = zeros(length(times))

    for j in 1:trials
        temp_res = zeros(Int, length(times))
        temp_ex = zeros(Int, length(times))
        run_length = (length(times)+1)
        input = gen_input(rand(rng, TriangularDist(0.01,  0.03, θ.low_mode)), size*sc)
        #var = [1 + σ * randn() for i in 1:14]

        for jj in 1:niche
            mut_cells, ext = variability_cml_free(times, input, θ, run_length, sc, σ; rng = rng)
            temp_res .= clamp.(temp_res .+ mut_cells, 0, 1)
            temp_ex .= clamp.(temp_ex .+ ext, 0, 1)
        end
        ii .+= temp_res
        ii2 .+= temp_ex
    end
    for j in 1:trials
        temp_res = zeros(Int, length(times))
        temp_ex = zeros(Int, length(times))
        run_length = (length(times)+1)
        input = gen_input(rand(rng, TriangularDist(0.03,  0.075, θ.mid_mode)), size*sc)
        #var = [1 + σ * randn() for i in 1:14]

        for jj in 1:niche
            mut_cells, ext = variability_cml_free(times, input, θ, run_length, sc, σ; rng = rng)
            temp_res .= clamp.(temp_res .+ mut_cells, 0, 1)
            temp_ex .= clamp.(temp_ex .+ ext, 0, 1)
        end
        iii .+= temp_res
        iii2 .+= temp_ex
    end
    for j in 1:trials
        temp_res = zeros(Int, length(times))
        temp_ex = zeros(Int, length(times))
        run_length = (length(times)+1)
        input = gen_input(rand(rng, TriangularDist(0.075,  0.1, θ.high_mode)), size*sc)
        #var = [1 + σ * randn() for i in 1:14]

        for jj in 1:niche
            mut_cells, ext = variability_cml_free(times, input, θ, run_length, sc, σ; rng = rng)
            temp_res .= clamp.(temp_res .+ mut_cells, 0, 1)
            temp_ex .= clamp.(temp_ex .+ ext, 0, 1)
        end
        iiii .+= temp_res
        iiii2 .+= temp_ex
    end
    return [[(ii*(-1) .+ trials)/trials, ii2/trials],
            [(iii*(-1) .+ trials)/trials, iii2/trials],
            [(iiii*(-1) .+ trials)/trials, iiii2/trials]]
end
