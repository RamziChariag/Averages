module averaging_methods

using ..function_generator
using Optim, Roots, ForwardDiff
using QuadGK

export average_on_tangent, average_by_geodesic, average_projected

"""
Method 1: Weighted average via projection onto tangent line
    average_on_tangent(f, df, a, b, w) -> (x*, y*)

Find a weighted "average" point between (a, f(a)) and (b, f(b))
using tangent-line projection.

# Arguments
- `f::Function`: Function defining the curve y = f(x).
- `df::Function`: Derivative of `f`.
- `a::Real`, `b::Real`: Endpoints on the x-axis.
- `w::Real`: Weight on point `a`. The weight on `b` is `1-w`.

# Method
1. For each candidate point M = (xM, f(xM)) on the curve:
   - Construct the tangent line at M.
   - Orthogonally project A = (a, f(a)) and B = (b, f(b)) onto that line.
   - Form the weighted average G = w*A_proj + (1-w)*B_proj on the line.
   - Compute squared distance between G and M.
2. Choose xM in [a,b] that minimizes this deviation.

# Returns
- `(x*, f(x*))`: The point on the curve that best aligns with the tangent-based weighted average.

# Notes
- This is not the same as geodesic averaging (by arc length).
- Instead, it uses a local linearization: the tangent at the solution point
is the line along which the weighted average of A and B projects back closest to M.
"""
function average_on_tangent(f, df, a, b, w)
    function deviation(xM)
        yM = f(xM)
        slope = df(xM)
        intercept = yM - slope * xM

        # Project arbitrary point (x,y) onto tangent line
        function project(x, y)
            denom = 1 + slope^2
            x_proj = (x + slope * (y - intercept)) / denom
            y_proj = (slope * x + slope^2 * y + intercept) / denom
            return (x_proj, y_proj)
        end

        A_proj = project(a, f(a))
        B_proj = project(b, f(b))
        Gx = w * A_proj[1] + (1 - w) * B_proj[1]
        Gy = w * A_proj[2] + (1 - w) * B_proj[2]

        return (Gx - xM)^2 + (Gy - yM)^2
    end

    result = optimize(deviation, a, b)
    x_opt = Optim.minimizer(result)
    return (x_opt, f(x_opt))
end


"""
Method 2: Geodesic approximation via arc-length minimization

Geodesic "average" on y=f(x) between x=a and x=b by arc-length.
Returns x* such that ∫_a^{x*} √(1 + f'(t)^2) dt = (1 - w) * ∫_a^b √(1 + f'(t)^2) dt.

Arguments:
- f :: Function      — scalar->scalar function
- a, b :: Real       — endpoints (can be in any order)
- w :: Real          — weight on point a (0≤w≤1). Weight on b is (1-w).
- df :: Function     — optional derivative; if not given, uses ForwardDiff

Notes:
- If you swap a and b, swap w ↦ 1-w for the same solution.
- Assumes the curve does not fold back in x (so arc length is monotone in x).
"""
function average_by_geodesic(f, a::Real, b::Real, w::Real; df=nothing)
    # Handle trivial / ordering
    a == b && return a
    if b < a
        # keep the objective the same by swapping weights
        x = average_by_geodesic(f, b, a, 1 - w; df=df)
        return x
    end

    # Derivative
    dfdx = df === nothing ? (x -> ForwardDiff.derivative(f, x)) : df
    integrand = x -> sqrt(1 + dfdx(x)^2)

    # Total arc length from a to b
    L, _ = quadgk(integrand, a, b)

    # Target cumulative length from a
    target = (1 - w) * L

    # Monotone cumulative arc-length from a to x
    S = x -> quadgk(integrand, a, x)[1]

    # Invert S(x) = target on [a,b] via bisection (robust)
    xstar = find_zero(x -> S(x) - target, (a, b), Bisection(); atol=1e-10, rtol=1e-10)
    return xstar
end

"""
Method 3: Projection of euclidean average onto the curve
"""
function average_projected(f, a, b, w)
    x_avg = w * a + (1 - w) * b
    y_avg = w * f(a) + (1 - w) * f(b)
    obj(x) = (x - x_avg)^2 + (f(x) - y_avg)^2
    result = optimize(obj, min(a, b), max(a, b))
    return result.minimizer
end

end
