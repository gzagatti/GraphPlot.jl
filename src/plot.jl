using Colors

function getcycle(x, i)
    if typeof(x) <: AbstractArray
        return x[((i-1) % length(x))+1]
    else
        return x
    end
end

# this function is copy from [GraphLayout.jl](https://github.com/IainNZ/GraphLayout.jl) and make some modifications.
"""
Given a graph and two vectors of X and Y coordinates, returns
a Compose tree of the graph layout

**Arguments**

`G`
Graph to draw

`layout`
Optional. Layout algorithm. Currently can be one of [`random_layout`,
`circular_layout`, `spring_layout`, `shell_layout`, `stressmajorize_layout`,
`spectral_layout`].
Default: `spring_layout`

`locs_x, locs_y`
Locations of the nodes. Can be any units you want,
but will be normalized and centered anyway

`nodemapping`
Optional. A vector representing the mapped graph vertices. Default: `nothing`

`NODESIZE`
Optional. Max size for the nodes. Default: `3.0/sqrt(N)`

`nodesize`
Optional. Relative size for the nodes, can be a Vector. Default: `1.0`

`nodelabel`
Optional. Labels for the vertices, a Vector or nothing. Default: `nothing`

`nodelabelc`
Optional. Color for the node labels, can be a Vector. Default: `colorant"black"`

`nodelabeldist`
Optional. Distances for the node labels from center of nodes. Default: `0.0`

`nodelabelangleoffset`
Optional. Angle offset for the node labels. Default: `π/4.0`

`NODELABELSIZE`
Optional. Largest fontsize for the vertice labels. Default: `4.0`

`nodelabelsize`
Optional. Relative fontsize for the vertice labels, can be a Vector. Default: `1.0`

`nodefillc`
Optional. Color to fill the nodes with, can be a Vector. Default: `colorant"turquoise"`

`nodestrokec`
Optional. Color for the nodes stroke, can be a Vector. Default: `nothing`

`nodestrokelw`
Optional. Line width for the nodes stroke, can be a Vector. Default: `0.0`

`edgelabel`
Optional. Labels for the edges, a Vector or nothing. Default: `[]`

`edgelabelc`
Optional. Color for the edge labels, can be a Vector. Default: `colorant"black"`

`edgelabeldistx, edgelabeldisty`
Optional. Distance for the edge label from center of edge. Default: `0.0`

`EDGELABELSIZE`
Optional. Largest fontsize for the edge labels. Default: `4.0`

`edgelabelsize`
Optional. Relative fontsize for the edge labels, can be a Vector. Default: `1.0`

`EDGELINEWIDTH`
Optional. Max line width for the edges. Default: `0.25/sqrt(N)`

`edgelinewidth`
Optional. Relative line width for the edges, can be a Vector. Default: `1.0`

`edgestrokec`
Optional. Color for the edge strokes, can be a Vector. Default: `colorant"lightgray"`

`arrowlengthfrac`
Optional. Fraction of line length to use for arrows.
Equal to 0 for undirected graphs. Default: `0.1` for the directed graphs

`arrowangleoffset`
Optional. Angular width in radians for the arrows. Default: `π/9 (20 degrees)`

"""
function gplot(g::AbstractGraph{T},
    locs_x_in::Vector{R}, locs_y_in::Vector{R};
    nodemapping = nothing,
    nodelabel = nothing,
    nodelabelc = colorant"black",
    nodelabelsize = 1.0,
    NODELABELSIZE = 4.0,
    nodelabeldist = 0.0,
    nodelabelangleoffset = π / 4.0,
    edgelabel = [],
    edgelabelc = colorant"black",
    edgelabelsize = 1.0,
    EDGELABELSIZE = 4.0,
    edgestrokec = colorant"lightgray",
    edgelinewidth = 1.0,
    EDGELINEWIDTH = 3.0 / sqrt(nv(g)),
    edgelabeldistx = 0.0,
    edgelabeldisty = 0.0,
    nodesize = 1.0,
    NODESIZE = 0.1 / sqrt(nv(g)),
    nodefillc = colorant"turquoise",
    nodestrokec = nothing,
    nodestrokelw = 0.0,
    arrowlengthfrac = is_directed(g) ? 0.1 : 0.0,
    arrowangleoffset = π / 9.0,
    linetype = "straight",
    outangle = pi/5) where {T <:Integer, R <: Real}

    length(locs_x_in) != length(locs_y_in) && error("Vectors must be same length")
    N = isnothing(nodemapping) ? nv(g) : length(nodemapping)
    NE = isnothing(nodemapping) ? ne(g) : reduce(+, map(e -> ((src(e) in nodemapping) && (dst(e) in nodemapping)) ? 1 : 0, edges(g)), init=0)
    if nodelabel != nothing && length(nodelabel) != N
        error("Must have one label per node (or none)")
    end
    if !isempty(edgelabel) && length(edgelabel) != NE
        error("Must have one label per edge (or none)")
    end

    locs_x = Float64.(locs_x_in)
    locs_y = Float64.(locs_y_in)

    # Determine sizes
    max_nodesize = NODESIZE / maximum(nodesize)
    nodesize *= max_nodesize
    max_edgelinewidth = EDGELINEWIDTH / maximum(edgelinewidth)
    edgelinewidth *= max_edgelinewidth
    max_edgelabelsize = EDGELABELSIZE / maximum(edgelabelsize)
    edgelabelsize *= max_edgelabelsize
    max_nodelabelsize = NODELABELSIZE / maximum(nodelabelsize)
    nodelabelsize *= max_nodelabelsize
    max_nodestrokelw = maximum(nodestrokelw)
    if max_nodestrokelw > 0.0
        max_nodestrokelw = EDGELINEWIDTH / max_nodestrokelw
        nodestrokelw *= max_nodestrokelw
    end

    # Scale to unit square
    min_x, max_x = extrema(locs_x)
    min_y, max_y = extrema(locs_y)
    function scaler(z, a, b)
        if (a - b) == 0.0
            return 0.5
        else
            return ((z - a) / (b - a))
        end
    end
    map!(z -> scaler(z, min_x, max_x), locs_x, locs_x)
    map!(z -> scaler(z, min_y, max_y), locs_y, locs_y)

    min_x = 0.1 -minimum([x - getcycle(nodesize, i) for (i, x) in enumerate(locs_x)])
    max_x = maximum([x + getcycle(nodesize, i) for (i, x) in enumerate(locs_x)]) - 0.9
    min_y = 0.1 -minimum([y - getcycle(nodesize, i) for (i, y) in enumerate(locs_y)])
    max_y = maximum([y + getcycle(nodesize, i) for (i, y) in enumerate(locs_y)]) - 0.9

    function rescaler(z, a, b)
        if (a - b) == 0.0
            return 0.5
        else
            r = maximum([a, b])
            return ((z + r*(1+2r)) / (1+2r))
        end
    end
    map!(z -> rescaler(z, min_x, max_x), locs_x, locs_x)
    map!(z -> rescaler(z, min_y, max_y), locs_y, locs_y)

    # Create nodes
    rs = [getcycle(nodesize, i) for i in 1:length(N)]
    nodes = circle(locs_x, locs_y, rs)

    # Create node labels if provided
    texts = []
    if !isnothing(nodelabel)
        text_locs_x = deepcopy(locs_x)
        text_locs_y = deepcopy(locs_y)
        texts = [
            (
                context(),
                text(
                    text_locs_x[i] + getcycle(nodesize, i) * (nodelabeldist * cos(nodelabelangleoffset)),
                    text_locs_y[i] - getcycle(nodesize, i) * (nodelabeldist * sin(nodelabelangleoffset)),
                    string(getcycle(nodelabel,i)), 
                    hcenter, 
                    vcenter
                ),
                fill(getcycle(nodelabelc,i)),
                stroke(getcycle(nodelabelc, i)),
                fontsize(getcycle(nodelabelsize, i)),
            )
            for i in 1:length(nodelabel)
        ]
    end

    # Create edge labels if provided
    edgetexts = []
    if !isempty(edgelabel)
        edge_locs_x = zeros(R, NE)
        edge_locs_y = zeros(R, NE)
        iter = isnothing(nodemapping) ? edges(g) : [e for e in edges(g) if ((src(e) in nodemapping) && (dst(e) in nodemapping))]
        for (e_idx, e) in enumerate(iter)
            i = src(e)
            j = dst(e)
            if !isnothing(nodemapping)
                i = findfirst(nodemapping .== i)
                j = findfirst(nodemapping .== j)
            end
            mid_x = (locs_x[i]+locs_x[j]) / 2.0
            mid_y = (locs_y[i]+locs_y[j]) / 2.0
            edge_locs_x[e_idx] = (is_directed(g) ? (mid_x+locs_x[j]) / 2.0 : mid_x) + edgelabeldistx * NODESIZE
            edge_locs_y[e_idx] = (is_directed(g) ? (mid_y+locs_y[j]) / 2.0 : mid_y) + edgelabeldisty * NODESIZE

        end
        edgetexts = [
            (
                context(),
                text(
                    edge_locs_x[i], 
                    edge_locs_y[i], 
                    string(getcycle(edgelabel, i)), 
                    hcenter, 
                    vcenter
                ),
                fill(getcycle(edgelabelc, i)),
                stroke(nothing),
                fontsize(getcycle(edgelabelsize, i)),
            )
            for i in 1:length(edgelabel)
        ]
    end

    # Create lines and arrow heads
    lines, arrows, stubs = nothing, nothing, nothing
    if linetype == "curve"
        if arrowlengthfrac > 0.0
            curves_cord, arrows_cord = graphcurve(g, locs_x, locs_y, nodemapping, nodesize, arrowlengthfrac, arrowangleoffset, outangle)
            lines = curve(curves_cord[:,1], curves_cord[:,2], curves_cord[:,3], curves_cord[:,4])
            arrows = line(arrows_cord)
        else
            curves_cord = graphcurve(g, locs_x, locs_y, nodemapping, nodesize, outangle)
            lines = curve(curves_cord[:,1], curves_cord[:,2], curves_cord[:,3], curves_cord[:,4])
        end
    else
        if arrowlengthfrac > 0.0
            lines_cord, arrows_cord = graphline(g, locs_x, locs_y, nodemapping, nodesize, arrowlengthfrac, arrowangleoffset)
            lines = line(lines_cord)
            arrows = line(arrows_cord)
        else
            lines_cord = graphline(g, locs_x, locs_y, nodemapping, nodesize)
            if !isnothing(nodemapping)
                stub_cord = stubline(g, locs_x, locs_y, nodemapping, nodesize)
                stubs = line(stub_cord)
            end
            lines = line(lines_cord)
        end
    end

    compose(
            context(),
            texts...,
            (context(), nodes, fill(nodefillc), stroke(nodestrokec), linewidth(nodestrokelw)),
            edgetexts...,
            (context(), arrows, stroke(edgestrokec), linewidth(edgelinewidth)),
            (context(), lines, stroke(edgestrokec), fill(nothing), linewidth(edgelinewidth)),
            (context(), stubs, stroke(edgestrokec), fill(nothing), linewidth(edgelinewidth))
    )
end

function gplot(g; layout::Function=spring_layout, keyargs...)
    locs_x, locs_y = layout(g)
    if haskey(keyargs, :nodemapping)
        locs_x = [locs_x[i] for i in keyargs[:nodemapping]]
        locs_y = [locs_y[i] for i in keyargs[:nodemapping]]
    end
    gplot(g, locs_x, locs_y; keyargs...)
end

# take from [Gadfly.jl](https://github.com/dcjones/Gadfly.jl)
function open_file(filename)
    if Sys.KERNEL == :Darwin
        run(`open $(filename)`)
    elseif Sys.KERNEL == :Linux || Sys.KERNEL == :FreeBSD
        run(`xdg-open $(filename)`)
    elseif Sys.KERNEL == :Windows
        run(`$(ENV["COMSPEC"]) /c start $(filename)`)
    else
        @warn("Showing plots is not supported on OS $(string(Sys.KERNEL))")
    end
end

# taken from [Gadfly.jl](https://github.com/dcjones/Gadfly.jl)
function gplothtml(g; layout::Function=spring_layout, keyargs...)
    filename = string(tempname(), ".html")
    output = open(filename, "w")

    plot_output = IOBuffer()
    draw(SVGJS(plot_output, Compose.default_graphic_width,
               Compose.default_graphic_width, false), gplot(g, layout(g)...; keyargs...))
    plotsvg = String(take!(plot_output))

    write(output,
        """
        <!DOCTYPE html>
        <html>
          <head>
            <title>GraphPlot Plot</title>
            <meta charset="utf-8">
          </head>
            <body>
            <script charset="utf-8">
                $(read(Compose.snapsvgjs, String))
            </script>
            <script charset="utf-8">
                $(read(gadflyjs, String))
            </script>
            $(plotsvg)
          </body>
        </html>
        """)
    close(output)
    open_file(filename)
end
