module plotting

using Plots

export plot_averaging_result, set_interactive_mode

# Plot backend switcher
function set_interactive_mode(interactive::Bool)
    if interactive
        try
            plotlyjs()
        catch
            println("PlotlyJS not available, using default.")
        end
    else
        gr()
    end
end

function plot_averaging_result(
    f, a, b,
    M_points::Vector{Tuple{Float64, Float64}},
    labels::Vector{String},
    plotdomain::Tuple{Float64, Float64}
)
    xs = range(plotdomain[1], plotdomain[2]; length=500)
    ys = f.(xs)

    plt = plot(xs, ys, label="f(x)", linewidth=2, legend=:bottomright)

    # Plot A and B
    scatter!([a], [f(a)], color=:red, markersize=8, label="")
    scatter!([b], [f(b)], color=:blue, markersize=8, label="")

    # Annotate A and B
    annotate!(a, f(a), text("A", :top, 15))
    annotate!(b, f(b), text("B", :bottom, 15))

    # Plot M points
    for (i, M) in enumerate(M_points)
        scatter!([M[1]], [M[2]], label=labels[i], marker=:x, markersize=10)
    end

    xlabel!("x")
    ylabel!("f(x)")
    title!("Weighted Average Points on f(x)")
    return plt
end

end
