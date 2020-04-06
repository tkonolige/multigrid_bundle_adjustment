# Convert a prolongation matrix P to a tentative operator S.
# P is piecewise constant over some number of aggregates.
# S has a n x n for each aggregate where n is the size of the aggregate.

using PETScBinaryIO
using SparseArrays
using LinearAlgebra

function p_to_s(P)
    rows = rowvals(P)
    vals = nonzeros(P)
    m, n = size(P)
    is = Vector{Int}()
    js = Vector{Int}()
    ks = Vector{Float64}()
    j = 0
    bs = 16
    for i = 1:bs:n
        # get row indices of the aggregate
        inds = sort(unique(vcat(map(x -> rows[nzrange(P, i+x)], 0:(bs-1))...)))
        agg_size = length(inds)
        if agg_size > 10000
            println("Aggregate is too large ($agg_size)")
            exit(-1)
        end
        if agg_size == 1
            @warn "Ignoring aggregate of size 1"
            continue
        end
        B = Matrix(P[inds, i:(i+bs-1)])

        # I - B(BᵀB)⁻¹Bᵀ
        block = I - B * inv(B' * B) * B'
        # block = I - fill(1/agg_size, (agg_size, agg_size)) # B(BᵀB)⁻¹Bᵀ for piecewise constant

        U, S, _ = svd(block)
        l = findlast(x -> x > 1e-14, S)
        for x = 1:agg_size
            for y = 1:l
                push!(is, inds[x])
                push!(js, j + y)
                push!(ks, U[x,y])
            end
        end
        j += l

        # verify that S projects out the coarse space
        @assert norm(U[:, 1:l]' * B, Inf) < 1e14
    end
    sparse(is,js,ks,m,j)
end


if joinpath(pwd(), PROGRAM_FILE) == @__FILE__
    P = readpetsc(ARGS[1])[1]
    writepetsc(ARGS[2], p_to_s(P))
end
