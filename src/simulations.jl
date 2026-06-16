function sim_cycles_cml_free(
    times::Vector{Int},
    initial_condition::Vector{Int},
    θ::ModelParams,
    run_time::Int64,
    sc::Int64;
    rng=Random.default_rng()
)

    c0, c1, m0, m1 = initial_condition

    results = Int64[]
    results2 = Int64[]

    cycle_num = times[run_time - 1]

    @inbounds for i in 1:cycle_num
        c0, c1, m0, m1 = cell_cycle(
            c0, c1, m0, m1,
            θ,
            sc,
            0.0;
            rng = rng
        )

        if i in times
            push!(results, m0)
            push!(results2, m0)
        end

        if (m0 / (c0 + m0)) > θ.threshold
            if length(results) < length(times)
                for _ in 1:(length(times) - length(results))
                    push!(results, 10000000)
                    push!(results2, m0)
                end
            end
            break

        elseif m0 < 1
            if length(results) < length(times)
                for _ in 1:(length(times) - length(results))
                    push!(results, 0)
                    push!(results2, 0)
                end
            end
            break
        end
    end

    if length(results) < length(times)
        for _ in 1:(length(times) - length(results))
            push!(results, 0)
            push!(results2, 0)
        end
    end

    cml_event = Int.(results .>= 10000000)
    nonextinct = Int.(results2 .> 0)

    return cml_event, nonextinct
end

function variability_cml_free(
    times::Vector{Int},
    input::Vector{Int},
    θ_base::ModelParams,
    run_length::Int64,
    sc::Int64,
    σ::Float64;
    rng=Random.default_rng()
)

    θ_mouse = add_variability(θ_base, σ; rng=rng)

    return sim_cycles_cml_free(
        times,
        input,
        θ_mouse,
        run_length,
        sc;
        rng = rng
    )
end


function sim_cycles_cell_numbers(
    times::Vector{Int},
    cycles::Int64,
    initial_condition::Vector{Int},
    θ::ModelParams,
    sc::Int64,
    asym::Float64,
    stop::Float64;
    rng=Random.default_rng()
)

    c0, c1, m0, m1 = Int.(initial_condition * sc)

    results_m0 = Int64[]
    results_c0 = Int64[]
    results_m1 = Int64[]
    results_c1 = Int64[]

    @inbounds for i in 1:cycles
        if (m0 / (c0 + m0)) > stop
            if length(results_m0) < cycles ÷ 10
                for _ in 1:((cycles ÷ 10) - length(results_m0))
                    push!(results_m0, 100000)
                    push!(results_c0, c0)
                    push!(results_m1, m1)
                    push!(results_c1, c1)
                end
            end
            break

        elseif m0 < 1
            if length(results_m0) < cycles ÷ 10
                for _ in 1:((cycles ÷ 10) - length(results_m0))
                    push!(results_m0, m0)
                    push!(results_c0, c0)
                    push!(results_m1, m1)
                    push!(results_c1, c1)
                end
            end
            break
        end

        c0, c1, m0, m1 = cell_cycle(
            c0, c1, m0, m1,
            θ,
            sc,
            asym;
            rng = rng
        )

        if i % 10 == 0
            push!(results_m0, m0)
            push!(results_c0, c0)
            push!(results_m1, m1)
            push!(results_c1, c1)
        end
    end

    return results_m0, results_c0, results_m1, results_c1
end


function variability_cell_numbers(
    input::Vector{Int},
    times::Vector{Int},
    θ_base::ModelParams,
    niche::Int,
    sc::Int,
    asym::Float64,
    σ::Float64,
    stop::Float64;
    rng=Random.default_rng()
)

    all_m0 = []
    all_c0 = []
    all_m1 = []
    all_c1 = []

    for _ in 1:niche
        θ_mouse = add_variability(θ_base, σ; rng=rng)

        m0, c0, m1, c1 = sim_cycles_cell_numbers(
            times,
            last(times),
            input,
            θ_mouse,
            sc,
            asym,
            stop;
            rng = rng
        )

        push!(all_m0, m0)
        push!(all_c0, c0)
        push!(all_m1, m1)
        push!(all_c1, c1)
    end

    return all_m0, all_c0, all_m1, all_c1
end
