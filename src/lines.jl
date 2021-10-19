"""
Return lines and arrow heads
"""
function graphline(g::AbstractGraph{T}, locs_x, locs_y, nodemapping, nodesize::Union{R, Vector{R}}, arrowlength, angleoffset) where {T<:Integer, R<:Real}
    NE = isnothing(nodemapping) ? ne(g) : reduce(+, map(e -> ((src(e) in nodemapping) && (dst(e) in nodemapping)) ? 1 : 0, edges(g)), init=0)
    iter = isnothing(nodemapping) ? edges(g) : [e for e in edges(g) if ((src(e) in nodemapping) && (dst(e) in nodemapping))]
    lines = Array{Vector{Tuple{Float64,Float64}}}(undef, NE)
    arrows = Array{Vector{Tuple{Float64,Float64}}}(undef, NE)
    for (e_idx, e) in enumerate(iter)
        i = src(e)
        j = dst(e)
        if !isnothing(nodemapping)
            i = findfirst(nodemapping .== i)
            j = findfirst(nodemapping .== j)
        end
        Δx = locs_x[j] - locs_x[i]
        Δy = locs_y[j] - locs_y[i]
        θ = atan(Δy,Δx)
        adj = θ > 0 ? (θ > π/2 ? abs(θ - π) : θ) : (θ < -π/2 ? θ + π : abs(θ))
        adj = (1 + (0.8*adj)/π)
        startx = locs_x[i] + adj*getcycle(nodesize, i)*cos(θ)
        starty = locs_y[i] + adj*getcycle(nodesize, i)*sin(θ)
        endx = locs_x[j] + adj*getcycle(nodesize, j)*cos(θ+π)
        endy = locs_y[j] + adj*getcycle(nodesize, j)*sin(θ+π)
        lines[e_idx] = [(startx, starty), (endx, endy)]
        arr1, arr2 = arrowcoords(θ, endx, endy, arrowlength, angleoffset)
        arrows[e_idx] = [arr1, (endx, endy), arr2]
    end
    lines, arrows
end

function graphline(g::AbstractGraph{T}, locs_x, locs_y, nodemapping, nodesize::Union{R, Vector{R}}) where {T<:Integer, R<:Real}
    NE = isnothing(nodemapping) ? ne(g) : reduce(+, map(e -> ((src(e) in nodemapping) && (dst(e) in nodemapping)) ? 1 : 0, edges(g)), init=0)
    iter = isnothing(nodemapping) ? edges(g) : [e for e in edges(g) if ((src(e) in nodemapping) && (dst(e) in nodemapping))]
    lines = Array{Vector{Tuple{Float64,Float64}}}(undef, NE)
    for (e_idx, e) in enumerate(iter)
        i = src(e)
        j = dst(e)
        if !isnothing(nodemapping)
            i = findfirst(nodemapping .== i)
            j = findfirst(nodemapping .== j)
        end
        Δx = locs_x[j] - locs_x[i]
        Δy = locs_y[j] - locs_y[i]
        θ = atan(Δy,Δx)
        adj = θ > 0 ? (θ > π/2 ? abs(θ - π) : θ) : (θ < -π/2 ? θ + π : abs(θ))
        adj = (1 + (0.8*adj)/π)
        startx = locs_x[i] + adj*getcycle(nodesize, i)*cos(θ)
        starty = locs_y[i] + adj*getcycle(nodesize, i)*sin(θ)
        endx = locs_x[j] + adj*getcycle(nodesize, j)*cos(θ+π)
        endy = locs_y[j] + adj*getcycle(nodesize, j)*sin(θ+π)
        lines[e_idx] = [(startx, starty), (endx, endy)]
    end
    lines
end

function stubline(g::AbstractGraph{T}, locs_x, locs_y, nodemapping, nodesize::Union{R, Vector{R}}) where {T<:Integer, R<:Real}
    lines = Vector{Tuple{Float64, Float64}}[]
    for (i, node) in enumerate(nodemapping)
        visible = [i for (i, n) in enumerate(nodemapping) if n in neighbors(g, node)]
        Δx = [locs_x[j] - locs_x[i] for j in visible]
        Δy = [locs_y[j] - locs_y[i] for j in visible]
        visible_θ = [atan(Δy[i], Δx[i]) for i in 1:length(Δx)]
        stubs = length(neighbors(g, node)) - length(visible)
        θ = 2*π / stubs
        for k in 0:(stubs-1)
            angle = θ*k
            if angle > π
                angle = angle - 2π
            end
            if angle in visible_θ
                angle += π / (2*stubs) 
            end
            adj = angle > 0 ? (angle > π/2 ? abs(angle - π) : angle) : (angle < -π/2 ? angle + π : abs(angle))
            adj = (0.8*adj)/π
            startx = locs_x[i] + (1+adj)*getcycle(nodesize, i)*cos(angle)
            starty = locs_y[i] + (1+adj)*getcycle(nodesize, i)*sin(angle)
            endx = locs_x[i] + (2.5+adj)*getcycle(nodesize, i)*cos(angle)
            endy = locs_y[i] + (2.5+adj)*getcycle(nodesize, i)*sin(angle)
            push!(lines, [(startx, starty), (endx, endy)])
        end
    end
    lines
end

function graphcurve(g::AbstractGraph{T}, locs_x, locs_y, nodemapping, nodesize::Union{R, Vector{R}}, arrowlength, angleoffset, outangle=pi/5) where {T<:Integer, R<:Real}
    NE = isnothing(nodemapping) ? ne(g) : reduce(+, map(e -> ((src(e) in nodemapping) && (dst(e) in nodemapping)) ? 1 : 0, edges(g)), init=0)
    iter = isnothing(nodemapping) ? edges(g) : [e for e in edges(g) if ((src(e) in nodemapping) && (dst(e) in nodemapping))]
    curves = Matrix{Tuple{Float64,Float64}}(undef, ne(g), 4)
    arrows = Array{Vector{Tuple{Float64,Float64}}}(undef, NE)
    for (e_idx, e) in enumerate(iter)
        i = src(e)
        j = dst(e)
        if !isnothing(nodemapping)
            i = findfirst(nodemapping .== i)
            j = findfirst(nodemapping .== j)
        end
        Δx = locs_x[j] - locs_x[i]
        Δy = locs_y[j] - locs_y[i]
        θ = atan(Δy,Δx)
        adj = θ > 0 ? (θ > π/2 ? abs(θ - π) : θ) : (θ < -π/2 ? θ + π : abs(θ))
        adj = (1 + (0.8*adj)/π)
        startx = locs_x[i] + adj*getcycle(nodesize, i)*cos(θ+outangle)
        starty = locs_y[i] + adj*getcycle(nodesize, i)*sin(θ+outangle)
        endx = locs_x[j] + adj*getcycle(nodesize, j)*cos(θ+π-outangle)
        endy = locs_y[j] + adj*getcycle(nodesize, j)*sin(θ+π-outangle)

        d = hypot(endx-startx, endy-starty)

        if i == j
            d = 2 * π * getcycle(nodesize, i)
        end

        curves[e_idx, :] = curveedge(startx, starty, endx, endy, θ, outangle, d)

        arr1, arr2 = arrowcoords(θ-outangle, endx, endy, arrowlength, angleoffset)
        arrows[e_idx] = [arr1, (endx, endy), arr2]
    end
    return curves, arrows
end

function graphcurve(g::AbstractGraph{T}, locs_x, locs_y, nodemapping, nodesize::Union{R, Vector{R}}, outangle) where {T<:Integer, R<:Real}
    NE = isnothing(nodemapping) ? ne(g) : reduce(+, map(e -> ((src(e) in nodemapping) && (dst(e) in nodemapping)) ? 1 : 0, edges(g)), init=0)
    iter = isnothing(nodemapping) ? edges(g) : [e for e in edges(g) if ((src(e) in nodemapping) && (dst(e) in nodemapping))]
    curves = Matrix{Tuple{Float64,Float64}}(undef, NE, 4)
    for (e_idx, e) in enumerate(iter)
        i = src(e)
        j = dst(e)
        if !isnothing(nodemapping)
            i = findfirst(nodemapping .== i)
            j = findfirst(nodemapping .== j)
        end
        Δx = locs_x[j] - locs_x[i]
        Δy = locs_y[j] - locs_y[i]
        θ = atan(Δy,Δx)
        adj = θ > 0 ? (θ > π/2 ? abs(θ - π) : θ) : (θ < -π/2 ? θ + π : abs(θ))
        adj = (1 + (0.8*adj)/π)
        startx = locs_x[i] + adj*getcycle(nodesize, i)*cos(θ+outangle)
        starty = locs_y[i] + adj*getcycle(nodesize, i)*sin(θ+outangle)
        endx = locs_x[j] + adj*getcycle(nodesize, j)*cos(θ+π-outangle)
        endy = locs_y[j] + adj*getcycle(nodesize, j)*sin(θ+π-outangle)

        d = hypot(endx-startx, endy-starty)

        if i == j
            d = 2 * π * nodesize[i]
        end

        curves[e_idx, :] = curveedge(startx, starty, endx, endy, θ, outangle, d)
    end
    return curves
end

# this function is copy from [IainNZ](https://github.com/IainNZ)'s [GraphLayout.jl](https://github.com/IainNZ/GraphLayout.jl)
function arrowcoords(θ, endx, endy, arrowlength, angleoffset=20.0/180.0*π)
    arr1x = endx - arrowlength*cos(θ+angleoffset)
    arr1y = endy - arrowlength*sin(θ+angleoffset)
    arr2x = endx - arrowlength*cos(θ-angleoffset)
    arr2y = endy - arrowlength*sin(θ-angleoffset)
    return (arr1x, arr1y), (arr2x, arr2y)
end

function curveedge(x1, y1, x2, y2, θ, outangle, d; k=0.5)

    r = d * k

    # Control points for left bending curve.
    xc1 = x1 + r * cos(θ + outangle)
    yc1 = y1 + r * sin(θ + outangle)

    xc2 = x2 + r * cos(θ + π - outangle)
    yc2 = y2 + r * sin(θ + π - outangle)

    return [(x1,y1) (xc1, yc1) (xc2, yc2) (x2, y2)]
end
