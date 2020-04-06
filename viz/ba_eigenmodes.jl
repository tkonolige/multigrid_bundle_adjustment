using Makie
using HDF5
using PETScBinaryIO
using BALUtils
using ArgParse
using LinearAlgebra
using SparseArrays
using ColorSchemes
using Colors
using LightGraphs
using NetworkLayout
using SimpleWeightedGraphs
using RCall
using DataFrames
using Statistics

function poses2centers(poses)
    cameras = mapslices(x -> BALUtils.Camera(x), poses, dims=1)
    hcat(Vector.(BALUtils.center.(cameras))...)
end

function draw_diff!(scene :: Scene, initial :: Array, change :: Array, scale :: Float64; kwargs...)
    moved = initial .+ change * scale
    paired = reshape(vcat(initial, moved), size(initial,1), :)
    # linesegments!(scene, paired[1,:], paired[2,:], paired[3,:], transparency=true, show_axis=false; kwargs...)
    c = change .* scale
    arrows!(scene, initial[1,:], initial[2,:], initial[3,:], c[1,:], c[2,:], c[3,:], transparency=true, show_axis=false; kwargs...)
end

function setup_scene!(s :: Scene, points)
    plot_center = mapslices(median, points, dims=2)[:]
    old_center = mapslices(mean, points, dims=2)[:]
    plot_center = old_center
    function setcam()
        cam3d_cad!(s)
        update_cam!(s)
        off = s.data_limits[].widths[2]
        cam = cameracontrols(s)
        cam.near[] = 0.0001
        update_cam!(s, plot_center .+ [0., 170, 0.], plot_center)
    end
    lift(a -> setcam(), s.events.window_area)
    lift(a -> setcam(), s.events.window_open)
end

function gplot(A :: SparseMatrixCSC; vertexcolors=nothing, vertexcolors2=nothing)
    net = NetworkLayout.Circular.layout(A)
    s = Scene()
    is, js, ks = findnz(A)
    ls = map((i,j) -> net[i] => net[j], is, js)
    (min, max) = extrema(ks)
    linecolors = RGBA.(0.0, 0.0, 0.0, (ks .- min) ./ (max - min) .* 0.9 .+ 0.1)
    linesegments!(s, ls, color=linecolors, transparecy=true, linewidth=3)
    scatter!(s, net, color=vertexcolors, transparency=true, colormap=:pu_or, markersize=0.2)
    #scatter!(s, net, transparency=true, colormap=:blues, markersize=0.1)
    s
end

function gplot3d(A :: SparseMatrixCSC, poses :: Matrix; vertexcolors=nothing, vertexcolors2=nothing)
    centers = mapslices(x -> BALUtils.center(BALUtils.Camera(x)), poses, dims=1)
    s = Scene()
    is, js, ks = findnz(A)
    ls = mapreduce((i,j) -> hcat(centers[:,i], centers[:,j]), hcat, is, js)
    (min, max) = extrema(ks)
    linecolors = RGBA.(0.0, 0.0, 0.0, (ks .- min) ./ (max - min) .* 0.9 .+ 0.1)
    #linesegments!(s, ls[1,:], ls[2,:], ls[3,:], transparency=true, color=linecolors, linewidth=3)
    scatter!(s, centers[1,:], centers[2,:], centers[3,:], color=vertexcolors, transparency=true, colormap=:pu_or, markersize=0.02)
    #scatter!(s, net, transparency=true, colormap=:blues, markersize=0.1)
    s
end

function plot_eigen_distribution(d, index)
    bal_file = joinpath(d, "problem_noised.bal")
    bbal_file = joinpath(d, "problem_noised.bbal")
    ba = if isfile(bal_file)
        readbal(bal_file)
    else
        readbal(bbal_file)
    end
    dscale = map(x -> x[1:num_cameras(ba)*9], readpetsc(joinpath(d, "dscale.petsc")))

    r = 0:300

    x = map(r) do i
        eig = h5read(joinpath(d, "eig_$(index)_smallest.h5"), "/eigenvector$i") .* reshape(dscale[index], 9, :)
        sort(mapslices(norm, eig[1:3,:], dims=1)[:]),
        sort(mapslices(norm, eig[4:6,:], dims=1)[:]),
        sort(mapslices(norm, eig[7:9,:], dims=1)[:])
    end
    rots = getindex.(x, 1)
    trans = getindex.(x, 2)
    intrin = getindex.(x, 3)

   #  df2 = vcat(
   #             DataFrame(norm=map(norm, rots), count=map(x -> count(x .> 1e-4)/length(x), rots),
   #                       eig=1:length(rots), group="rotation"),
   #             DataFrame(norm=map(norm, trans), count=map(x -> count(x .> 1e-2)/length(x), trans),
   #                       eig=1:length(trans), group="translation"),
   #             DataFrame(norm=map(norm, intrin), count=map(x -> count(x .> 1e-3)/length(x), intrin),
   #                       eig=1:length(intrin), group="intrinsic")
   # )
   #  ggplot(df2, aes(x=:eig, y=:count, color=:group)) + geom_point()
    df = vcat(
         mapreduce(x -> DataFrame(magnitude=x[2], x=1:length(x[2]), group="rotation", eig=x[1]), vcat, enumerate(rots)),
         mapreduce(x -> DataFrame(magnitude=x[2], x=1:length(x[2]), group="translation", eig=x[1]), vcat, enumerate(trans)),
         mapreduce(x -> DataFrame(magnitude=x[2], x=1:length(x[2]), group="intrinsic", eig=x[1]), vcat, enumerate(intrin)))
    ggplot(df, aes(x=:x, y=:magnitude, color=:eig)) + geom_point(size=0.1) + scale_y_log10()
end

function visualize(d, opts, index, eigindex, scale)
    solution = readpetsc(joinpath(d, "solution.petsc"))
    bal_file = joinpath(d, "problem_noised.bal")
    bbal_file = joinpath(d, "problem_noised.bbal")
    ba = if isfile(bal_file)
        readbal(bal_file)
    else
        readbal(bbal_file)
    end
    dscale = map(x -> x[1:num_cameras(ba)*9], readpetsc(joinpath(d, "dscale.petsc")))
    poses = map(x -> reshape(x[1:num_cameras(ba)*9], 9, :), solution)[index]
    centers = mapslices(x -> BALUtils.center(BALUtils.Camera(x)), poses, dims=1)

    S = readpetsc(joinpath(d, "bamg_$(index)_$(opts)_S.petsc"))[1]
    # network = NetworkLayout.SFDP.layout(S, Point3f0)
    # centers[1:3,:] = hcat(map(x -> Vector(x), network)...)

    eig_ = reshape(h5read(joinpath(d, "cr_eigs_$(index)_$(opts).h5"), "/eigenvector$(eigindex-1)"), 9, :)
    eig = eig_ .* reshape(dscale[index][1:num_cameras(ba)*9], 9, :)
    # To visualize modes of main matrix
    # eig_ = h5read(joinpath(d, "eig_$(index)_smallest.h5"), "/eigenvector$(eigindex-1)")
    # eig = eig_ .* reshape(dscale[index], 9, :)

    # eig = reshape(h5read(joinpath(d, "cr_eigs_$(index)_$(opts).h5"), "/eigenvector$eigindex") .* dscale[index], 9, :)
    displaced = poses .+ eig
    displaced_centers = mapslices(x -> BALUtils.center(BALUtils.Camera(x)), displaced, dims=1)
    eignorms = mapslices(x -> norm(x, Inf), eig, dims=1)[:]

    # get aggregate ids for each vertex
    P = readpetsc(joinpath(d, "bamg_$(index)_$(opts)_P.petsc"))[1]
    aggs = fill(-1, num_cameras(ba))
    is, js, _ = findnz(P)
    for (i, j) in zip(is, js)
        if (i-1) % 9 == 0 && (j - 1) % 9 == 0
            aggs[(i-1)÷9+1] = (j-1)÷9+1
        end
    end
    aggs_ = copy(aggs)
    aggs_[findall(aggs_ .== -1)] = maximum(aggs_)+1:maximum(aggs_)+count(aggs_ .== -1)
    colors = distinguishable_colors(maximum(aggs_))[aggs_]

    # verify that the eigenvector is orthogonal to prolongation on the aggregate
    aggids = unique(aggs)
    for i in aggids
        if i == -1
            continue
        end
        inds = findall(aggs .== i)
        binds = mapreduce(x -> (x-1)*9+1:x*9, vcat, inds)
        n = norm(P[binds, :]' * eig_[:, inds][:], Inf)
        if length(inds) > 1 && n > 1e-10
            @warn "Aggregate is in coarse space (inner product: $n)"
        end
    end

    aggids = unique(aggs)
    diams = map(aggids) do agg
        if agg == -1
            0
        else
            inds = aggs .== agg
            diameter(SimpleWeightedGraph(1 ./ S[inds, inds]))
        end
    end
    daggs = Dict(zip(aggids, diams))

    println("Largest parts of eigenvector")
    r = sortperm(eignorms, rev=true)[1:10]
    display(eig[:, r])
    println("")
    # println("A eigenvector")
    # display(A_eig[:, r])
    # println("")
    println("Associated poses")
    display(poses[:, r])

    # i = eignorms .> 0.1
    # @show S[i, i]
    # @show findall(i)
    # @show aggs[i]
    # j = aggs .== aggs[i][1]
    # @show findall(j)
    # @show eignorms[j]
    # @show S[j, j]
    # @show sum(S[j,:], dims=2)
    # @show sum(S[j,j], dims=2)
    # @show sum(S[j,j], dims=2) ./ sum(S[j,:], dims=2)
    # @show sum(S[j,j], dims=2) ./ sum(S[j,j])
    # @show eignorms[j]
    # @show diameter(SimpleWeightedGraph(1 ./ S[j, j]))
    # net = NetworkLayout.Spring.layout(S[j,j], Point2f0, iterations=100)
    # s2 = Scene()
    # is, js, ks = findnz(S[j,j])
    # inds = vcat(is, js)[:]
    # (min, max) = extrema(ks)
    # linecolors = RGBA.(0.0, 0.0, 0.0, (ks .- min) ./ (max - min))
    # linesegments!(s2, net[inds], color = linecolors)
    # scatter!(s2, net, color=eignorms[j],transparency=true,colormap=:blues)
    # # display(s2)
    # is, js, ks = findnz(S)
    # S_ = sparse(is, js, ones(size(ks)), size(S)...)
    # k = (S_ * S_* j) .!= 0.0
    # # display(gplot(S[k,k], vertexcolors=colors[k], vertexcolors2=eignorms[k]))

    # colors = eignorms

    s = Scene()
    scatter!(s, centers[1, :], centers[2, :], centers[3, :], color=colors
            , transparency=true, markersize=1.1, colormap=:blues)
    draw_diff!(s, centers, displaced_centers - centers, scale, transparency=true, linewidth=2.0)

    # draw strength of connection
    is, js, ks = findnz(S)
    (min, max) = extrema(ks)
    linecolors = RGBA.(0.0, 0.0, 0.0, (ks .- min) ./ (max - min))
    inds = vcat(is, js)[:]
    # linesegments!(s, centers[1, inds], centers[2, inds], centers[3, inds], color=linecolors, transparency=true)

    s
end

if PROGRAM_FILE == @__FILE__
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table s begin
        "input_dir"
        help = "Bundle adjustment problem input directory"
        required = true
        "eigs"
        required = true
        "index"
        required = true
        arg_type = Int
    end
    args = parse_args(ARGS, s)

    visualize(args["input_dir"], args["eigs"], args["index"])
end
