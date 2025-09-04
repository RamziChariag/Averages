include("functions.jl")
include("methods.jl")
include("plotting.jl")

# Parameters
a, b = 0.3, 2.0
w = 0.6 # Weight on A
r = 2.0  # Too small? It will be auto-corrected
# Choose from: SqrtFunction, LogFunction, QuarterCircleFunction
function_choice = function_generator.QuarterCircleFunction

# Choose function type here
f_struct, df, plotdomain = function_generator.get_function(function_choice, a, b; r=r)
f = x -> f_struct(x)

# Compute average points
M1 = averaging_methods.average_on_tangent(f, df, a, b, w)
M2x = averaging_methods.average_by_geodesic(f, a, b, w)
M2 = (M2x, f(M2x))
M3x = averaging_methods.average_projected(f, a, b, w)
M3 = (M3x, f(M3x))

# Plotting the results
plotting.set_interactive_mode(true)   # For interactive PlotlyJS plots
plt = plotting.plot_averaging_result(f, a, b, [M1, M2, M3], ["Method 1", "Method 2", "Method 3"], plotdomain)
display(plt)

# Print results
println("Method 1: ($(round(M1[1], digits=4)), $(round(M1[2], digits=4)))")
println("Method 2: ($(round(M2[1], digits=4)), $(round(M2[2], digits=4)))")
println("Method 3: ($(round(M3[1], digits=4)), $(round(M3[2], digits=4)))")
