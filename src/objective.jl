function loss(x, sc, trials, nn, σ; rng = rng)
    θ = theta_from_x(x)

    stem = 1200
    terminal = 12000
    # r = number that have developed CML, e = number that have not gone extinct (mutant stem cell number = 0)
    r1, e1 = simulate_survival_counts(
        θ, [stem * sc, terminal * sc],[30, 60, 90, 120, 150, 180, 210, 240],
        trials, nn, sc, σ, 0.01, 0.03, θ.low_mode; rng = rng
    )
    r2, e2 = simulate_survival_counts(
        θ, [stem * sc, terminal * sc],[30, 60, 90, 120, 150, 180, 210, 240],
        trials, nn, sc, σ, 0.03, 0.075, θ.mid_mode; rng = rng
    )
    r3, e3 = simulate_survival_counts(
        θ, [stem * sc, terminal * sc],[30, 60, 90, 120, 150, 180, 210, 240],
        trials, nn, sc, σ, 0.075, 0.10, θ.high_mode; rng = rng
    )
    pushfirst!(r1,0)
    pushfirst!(r2,0)
    pushfirst!(r3,0)
    e1 .-= trials
    e2 .-= trials
    e3 .-= trials
    e1 .*= -1
    e2 .*= -1
    e3 .*= -1
    pushfirst!(e1,0)
    pushfirst!(e2,0)
    pushfirst!(e3,0)
    bll1 = 0
    bll2 = 0
    bll3 = 0
    # data 1 are the number of cml free mice in the smoothed observed data for the lower initial chimerism mice
    data1 = GROUP_SIZES.low*OBSERVED_SURVIVAL.low
    for i in 2:length(r1)
        rawrisk = data1[(i-1)]
        rawevents = data1[(i-1)] - data1[i]


        risk = trials - r1[(i-1)] - e1[(i-1)]
        events = r1[i] - r1[(i-1)]
        events2 = e1[i] - e1[(i-1)]

        ##smooth it incase risk or events is 0
        prob = (events + 0.0000005)/(risk + 0.000001)

        #prob of extinction
        prob2 = (events2 + 0.0000005)/(risk + 0.000001)

        bll1 -= (rawevents*log(clamp(prob,0.00000001, 1.0)) + (rawrisk-rawevents)*log(clamp(1-prob, 0.0000001, 1.0)))
    end
    # data 2 are the number of cml free mice in the smoothed observed data for the middle initial chimerism mice
    data2 = GROUP_SIZES.mid*OBSERVED_SURVIVAL.mid
    for i in 2:length(r2)
        rawrisk = data2[(i-1)]
        rawevents = data2[(i-1)] - data2[i]


        risk = trials - r2[(i-1)] - e2[(i-1)]
        events = r2[i] - r2[(i-1)]
        events2 = e2[i] - e2[(i-1)]

        ##smooth it incase risk or events is 0
        prob = (events + 0.0000005)/(risk + 0.000001)

        #prob2 of extinction
        prob2 = (events2 + 0.0000005)/(risk + 0.000001)


        bll2 -= (rawevents*log(clamp(prob,0.00000001, 1.0)) + (rawrisk-rawevents)*log(clamp(1-prob-prob2, 0.0000001, 1.0)))
    end
    # data 3 are the number of cml free mice in the smoothed observed data for the lower initial chimerism mice
    data3 = GROUP_SIZES.high*OBSERVED_SURVIVAL.high
    for i in 2:length(r3)
        rawrisk = data3[(i-1)]
        rawevents = data3[(i-1)] - data3[i]


        risk = trials - r3[(i-1)] - e3[(i-1)]
        events = r3[i] - r3[(i-1)]
        events2 = e3[i] - e3[(i-1)]

        ##smooth it incase risk of events is 0
        prob = (events + 0.0000005)/(risk + 0.000001)

        #prob2 of extinction
        prob2 = (events2 + 0.0000005)/(risk + 0.000001)


        bll3 -= (rawevents*log(clamp(prob,0.00000001, 1.0)) + (rawrisk-rawevents)*log(clamp(1-prob-prob2, 0.0000001, 1.0)))
    end

    # now compute the loss for lower steady state.
    # we want there to be a nonzero steady state that has a chimerism below 1/300
    mut1 = []
    nor1 = []
    lk = ReentrantLock()

    @threads for i in 1:100
        input = gen_input(rand(TriangularDist(0.01,  0.03, θ.low_mode)),  [stem, terminal])
        mut, nor = variability_cell_numbers(input, [1, 400],  θ, nn, sc, 0.0, σ, θ.threshold; rng = rng)
        lock(lk) do
            if 0 < last(sum(mut)) < 100000
                push!(mut1, sum(mut))
                push!(nor1, sum(nor))
            end
        end
    end
    # we want more than half of the simulations to be at the lower steady state
    if length(mut1) < 50
        chim1 = 1
    else
        chim1 = last(mean(chimerism([mut1, nor1])))
    end
    check_lower_steady_state = 150*(((chim1*300)^5)/(1 + (chim1*300)^5))


    # we now bias the the number of cml free simulation in the high chimerism group to be below 5% after 400 cell cycles
    mut2 = []
    nor2 = []
    @threads for i in 1:200
        input = gen_input(rand(TriangularDist(0.075,  0.10, θ.high_mode)),  [stem, terminal])
        mut, nor = variability_cell_numbers(input, [1, 400], θ, nn, sc, 0.0, σ, θ.threshold; rng = rng)
        lock(lk) do
            if last(sum(mut)) < 100000
                push!(mut2, sum(mut))
                push!(nor2, sum(nor))
            end
        end
    end
    if  length(mut2) > 20
        ll3 = 2*length(mut2)
    elseif 20 > length(mut2) > 10
        ll3 = length(mut2)
    else
        ll3 = 0
    end
    #print out the loss to see how the fit is evolving
    #println([bll1, bll2, bll3, check_lower_steady_state, ll3])


    return bll1 + bll2 + bll3 + check_lower_steady_state + ll3

end
