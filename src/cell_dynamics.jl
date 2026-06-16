#given a renewal probability and the porportion of asymmetric renewal events when p = 0.5
#returns the probability of each division event outcome
function division_probs(p::Float64, asym::Float64)
    if p > 0.5
        a = (asym*2*p + 2*p - 2*asym)/(2(asym*2*p + 1 - 2*asym))
        asy = (2 - 2*p)*asym
        syex = (1 - asy)*(1-a)
        syre = (1 - asy)*(a)
    elseif p < 0.5
        a = (asym*2*p - 2*p)/(2(asym*2*p -1))
        asy = (2*p)*asym
        syex = (1 - asy)*(1-a)
        syre = (1 - asy)*(a)
    else
        a = 0.5
        asy = asym
        syex = (1 - asy)*(1-a)
        syre = (1 - asy)*(a)
    end

    return syex, asy, syre
end

#a single cell cycle
function cell_cycle(
    c0::Int64,
    c1::Int64,
    m0::Int64,
    m1::Int64,
    θ::ModelParams,
    sc::Int64,
    asym::Float64;
    rng=Random.default_rng()
)

    # Scale parameters for this simulation scale
    gamma_scaled = θ.γ / sc
    gammaL_scaled = gammaL(θ) / sc
    alpha_scaled = θ.α / sc
    phiL_scaled = θ.ϕL / sc

    phi_phiL = θ.ϕ * phiL_scaled

    cell_cycle_steps = 50
    division_chunks = allocate_counts_multinomial(c0, m0, cell_cycle_steps; rng=rng)

    @inbounds for i in 1:cell_cycle_steps
        if m0 < 0
            break
        end

        c0_div, m0_div = division_chunks[i]

        x = c1 + m1
        s = θ.β * c0 + m0

        #normal stem/progenitor cell division outcomes
        denom_inner = 1.0 / (1.0 + (phi_phiL * s)^θ.n)

        p = clamp(
            θ.pbar / (1.0 + gamma_scaled * x * denom_inner),
            0.0,
            1.0
        )

        syex, asy, syre = division_probs(p, asym)

        t_c0, t_c1 = sum(
            rand(rng, Multinomial(c0_div, [syex, asy, syre])) .*
            [[-1, 2], [0, 1], [1, 0]]
        )

        c0 += t_c0
        c1 += t_c1


        #mutant stem/progenitor cell division outcomes
        denom_inner = 1.0 / (1.0 + (phiL_scaled * s)^θ.n)
        pL = clamp(
            θ.pbar / (1.0 + gammaL_scaled * x * denom_inner),
            0.0,
            1.0
        )

        #compute the probability of each division outcome
        syex, asy, syre = division_probs(pL, asym)

        t_m0, t_m1 = sum(
            rand(rng, Multinomial(m0_div, [syex, asy, syre])) .*
            [[-1, 2], [0, 1], [1, 0]]
        )

        m0 += t_m0
        m1 += t_m1

        if m0 <= 0
            break
        end

        live = c0 + m0

        if live == 0 || live > 500000
            break
        end

        inv = (c0_div + m0_div) / live

        # Population-dependent mutant stem/progenitor death
        pos = (
            θ.α2 +
            ((alpha_scaled * m0)^θ.h) / (1 + (alpha_scaled * m0)^θ.h)
        )

        pc1 = clamp(-expm1(-θ.d * inv), 0.0, 1.0)

        pm0 = clamp(
            -expm1(-(θ.k * pos) * inv),
            1e-12,
            1.0 - 1e-12
        )

        pm1 = clamp(-expm1(-θ.d * θ.δ * inv), 0.0, 1.0)

        c1 -= rand(rng, Binomial(max(c1, 0), pc1))
        m0 -= rand(rng, Binomial(max(m0, 0), pm0))
        m1 -= rand(rng, Binomial(max(m1, 0), pm1))
    end

    return c0, c1, m0, m1
end
