module function_generator

export FunctionWithDerivative, SqrtFunction, LogFunction, QuarterCircleFunction, get_function

abstract type FunctionWithDerivative end

struct SqrtFunction <: FunctionWithDerivative end
struct LogFunction <: FunctionWithDerivative end
struct QuarterCircleFunction <: FunctionWithDerivative
    r::Float64
end

# Evaluation
(f::SqrtFunction)(x) = sqrt(x)
(f::LogFunction)(x) = log(x)
(f::QuarterCircleFunction)(x) = sqrt(f.r^2 - (f.r - x)^2)

# Derivatives
df(::SqrtFunction) = x -> 1 / (2 * sqrt(x))
df(::LogFunction) = x -> 1 / x
df(f::QuarterCircleFunction) = x -> (f.r - x) / sqrt(f.r^2 - (f.r - x)^2)

# Domains
domain(::SqrtFunction) = (0.0, Inf)
domain(::LogFunction) = (0.0, Inf)
domain(f::QuarterCircleFunction) = (0.0, f.r)

# Function getter
function get_function(::Type{QuarterCircleFunction}, a, b; r=nothing)
    default_r = max(a - 1, b + 1)

    if r === nothing
        r = default_r
    elseif max(a, b) > r
        println("⚠️  Provided radius r = $r is too small. Automatically adjusted to r = $default_r.")
        r = default_r
    end

    f = QuarterCircleFunction(r)
    d = domain(f)
    restricted = (max(d[1], a - 1), min(d[2], b + 1))
    return f, df(f), restricted
end

function get_function(::Type{T}, a, b; kwargs...) where T<:FunctionWithDerivative
    f = T()
    d = domain(f)
    restricted = (max(d[1], a - 1), min(d[2], b + 1))
    return f, df(f), restricted
end



end

