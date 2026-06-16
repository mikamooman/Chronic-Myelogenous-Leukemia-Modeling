function prog_nonprog(times, x, niche, trials, sc, asy, σ; rng=Random.default_rng())
    θ = theta_from_x(x)
    size = [1200, 12000]
    mutprog = []
    mutnon = []
    norprog = []
    nornon = []

    mutprogt = []
    mutnont = []
    norprogt = []
    nornont = []



    for i in 1:trials
        input = gen_input(rand(rng, TriangularDist(0.01,  0.03, θ.low_mode)), size)
        mut, nor, mut2, nor2 = variability_cell_numbers(input, times, θ, niche, sc, asy, σ, 0.99; rng = rng)
        tot = mut .+ nor
        l = length(tot)
        if maximum( mut[l] ./ tot[l] ) > θ.threshold
            push!(mutprog, sum(mut))
            push!(norprog, sum(nor))
            push!(mutprogt, sum(mut2))
            push!(norprogt, sum(nor2))
        else
            push!(mutnon, sum(mut))
            push!(nornon, sum(nor))
            push!(mutnont, sum(mut2))
            push!(nornont, sum(nor2))
        end
    end
    mutprog2 = []
    mutnon2 = []
    norprog2 = []
    nornon2 = []

    mutprog2t = []
    mutnon2t = []
    norprog2t = []
    nornon2t = []

    for i in 1:trials
        input = gen_input(rand(rng, TriangularDist(0.03,  0.075, θ.mid_mode)), size)
        mut, nor, mut2, nor2 = variability_cell_numbers(input, times, θ, niche, sc, asy, σ, 0.99; rng = rng)
        #maximum(last.(mut))
        tot = mut .+ nor
        l = length(tot)
        if maximum( mut[l] ./ tot[l] ) > θ.threshold
            push!(mutprog2, sum(mut))
            push!(norprog2, sum(nor))
            push!(mutprog2t, sum(mut2))
            push!(norprog2t, sum(nor2))
        else
            push!(mutnon2, sum(mut))
            push!(nornon2, sum(nor))
            push!(mutnon2t, sum(mut2))
            push!(nornon2t, sum(nor2))
        end
    end

    mutprog3 = []
    mutnon3 = []
    norprog3 = []
    nornon3 = []

    mutprog3t = []
    mutnon3t = []
    norprog3t = []
    nornon3t = []

    for i in 1:trials
        input = gen_input(rand(rng, TriangularDist(0.075,  0.10, θ.high_mode)), size)
        mut, nor, mut2, nor2 = variability_cell_numbers(input, times, θ, niche, sc, asy, σ, 0.99; rng = rng)
        tot = mut .+ nor
        l = length(tot)
        if maximum( mut[l] ./ tot[l] ) > θ.threshold
            push!(mutprog3, sum(mut))
            push!(norprog3, sum(nor))
            push!(mutprog3t, sum(mut2))
            push!(norprog3t, sum(nor2))
        else
            push!(mutnon3, sum(mut))
            push!(nornon3, sum(nor))
            push!(mutnon3t, sum(mut2))
            push!(nornon3t, sum(nor2))
        end
    end
    return [[[mutprog, mutnon, norprog, nornon],[mutprog2, mutnon2, norprog2, nornon2],[mutprog3, mutnon3, norprog3, nornon3]], [[mutprogt, mutnont, norprogt, nornont],[mutprog2t, mutnon2t, norprog2t, nornon2t],[mutprog3t, mutnon3t, norprog3t, nornon3t]]]
end
