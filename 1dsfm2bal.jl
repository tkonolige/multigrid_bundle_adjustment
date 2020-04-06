using LightGraphs
using LinearAlgebra
using BALUtils
using Rotations
using PyCall
using StaticArrays

# get camera locations
pairwise_poses = Dict()
es = []
for line in readlines(ARGS[1])
    s = split(line, " ")
    i, j = parse.(Int, s[1:2])
    push!(es, LightGraphs.SimpleEdge(i, j))
    # rotation is row major
    rot = transpose(reshape(parse.(Float64, s[3:11]), 3, 3))
    trans = parse.(Float64, s[12:14])
    pairwise_poses[(i,j)] = (rot, trans)
    pairwise_poses[(j,i)] = (transpose(rot), -(rot' * trans))
end
graph = SimpleGraphFromIterator(es)

# compute global rotations
mst = kruskal_mst(graph) |> SimpleDiGraphFromIterator
_, root = findmax(degree(mst))

poses = Dict()
poses[root] = (Matrix(I, 3, 3), zeros(3))
to_visit = [root]
while length(to_visit) > 0
    v = pop!(to_visit)
    rot_base, trans_base = poses[v]
    for n in neighbors(mst, v)

        rot, trans = pairwise_poses[v, n]
        poses[n] = (rot * rot_base, rot * trans_base + trans)

        push!(to_visit, n)
    end
end

# get point locations and image features
coords = map(x -> Vector(), 1:nv(graph))
camera_intrins = Dict()
open(ARGS[2]) do f
    while !eof(f)
        header = readline(f)
        index = parse(Int, match(r"index = ([0-9.]+)", header).captures[1])
        focal = parse(Float64, match(r"focal = ([0-9.]+)", header).captures[1])
        if haskey(poses, index)
            camera_intrins[index] = (focal, 0, 0)
        end

        px = parse(Float64, match(r"px = ([0-9.]+)", header).captures[1])
        py = parse(Float64, match(r"py = ([0-9.]+)", header).captures[1])

        keys = parse(Int, match(r"keys = ([0-9.]+)", header).captures[1])
        for i in 1:keys
            l = split(readline(f), " ")
            if haskey(poses, index)
                push!(coords[index], (parse(Int, l[1])
                                     , parse(Float64, l[2]) - px/2
                                     , parse(Float64, l[3]) - py/2
                                     ))
            end
        end
    end
end

# mapping from features to camera
features_map = Dict()
for (i, obs) in enumerate(coords)
    for o in obs
        x = o[1]
        y = (i, o[2], o[3])
        if haskey(features_map, x)
            push!(features_map[x], y)
        else
            features_map[x] = [y]
        end
    end
end
features_map = filter(x -> length(x[2]) >= 3, features_map)

renumbering = Dict(zip(1:length(poses), keys(poses)))
renumbering_inv = Dict(zip(keys(poses), 1:length(poses)))
cameras = map(1:length(poses)) do i
    rot, trans = poses[renumbering[i]]
    rv = RodriguesVec(rot)
    intrin = camera_intrins[renumbering[i]]
    Camera(SVector{3,Float64}(rv.sx, rv.sy, rv.sz), SVector{3,Float64}(trans), SVector{3,Float64}(intrin))
end

function intrinsic_matrix(intrin)
    f = intrin[1]
    [f 0 0
     0 f 0
     0 0 1]
end

function rot_matrix(x)
    RotMatrix(RodriguesVec(x[1], x[2], x[3]))
end

# triangulate points
projs = map(c -> intrinsic_matrix(c.intrinsics) * hcat(rot_matrix(c.rotation), c.pose), cameras)
opencv = pyimport("cv2")
points = map(collect(keys(features_map))) do feat
    obs = features_map[feat]
    points = map(x -> reshape([x[2], x[3]], 2, 1), obs)
    proj_mats = map(x -> projs[renumbering_inv[x[1]]], obs)
    points4d = opencv.triangulatePoints(proj_mats[1], proj_mats[2], points[1], points[2])
    points4d[1:3]/points4d[4]
end

points_map = Dict(zip(keys(features_map), 1:length(features_map)))
observations = map(1:length(cameras)) do i
    lobs = []
    for obs in coords[renumbering[i]]
        if haskey(points_map, obs[1])
            push!(lobs, (points_map[obs[1]], obs[2],  obs[3]))
        end
    end
    lobs
end

ba = BA(cameras, observations, points)
writebal(ARGS[3], ba)
