#generate an input vector with a desired chimersim
function gen_input(chimerism::Float64, a::Vector)
    b = [1 - chimerism, chimerism]

    return copy(vec(vec(round.(Int, (a .* b')')')))
end


#allocate cells into groups at the beggining of each cell_cycle
function allocate_counts_multinomial(R::Int, B::Int, N::Int; rng=Random.default_rng())
    p = fill(1.0 / N, N)
    r = rand(rng, Multinomial(R, p))
    b = rand(rng, Multinomial(B, p))
    return collect(zip(r, b))
end

#returns chimerism
function chimerism(a)
    mutant, normal = a
    tot = mutant .+ normal
    chimerism = [(mutant[i] ./ tot[i]) .+ 10^(-4) for i in 1:length(tot)]
    return chimerism
end

@inline function check_esc(esc_list::Vector{Int})
    x = findfirst(>(0), esc_list)
    if x == nothing
        return (length(esc_list)+1)
    else
        return x
    end
end
