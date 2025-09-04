# Averages
Geometric Averages on Curves — Julia implementation of unorthodox averaging methods (tangent projection, geodesic arc-length, and projection of Euclidean averages) for points on smooth curves. Includes utilities for function generation, optimization, and interactive plotting.

# Geometric Averages on Curves

This Julia project explores different ways of defining and computing **weighted averages** of two points lying on a curve \( y = f(x) \).  
Instead of averaging their coordinates directly, the methods incorporate the curve's geometry (via tangent lines, arc length, or projections).

## 🚀 Features

- **Function generator**:
  - Provides test functions with derivatives and domains:
    - `SqrtFunction`
    - `LogFunction`
    - `QuarterCircleFunction`
  - Ensures consistency between functions and their derivatives.
  
- **Averaging methods** (`averaging_methods.jl`):
  - **Tangent-based averaging**  
    Finds the curve point whose tangent line best aligns with the weighted average of projections of the endpoints.
  - **Geodesic averaging (arc length)**  
    Chooses the point that splits the arc length between `a` and `b` in proportion to the weights.
  - **Projection of Euclidean average**  
    Projects the straight-line weighted average of `(a, f(a))` and `(b, f(b))` back onto the curve.
  
- **Plotting utilities** (`plotting.jl`):
  - Interactive (`PlotlyJS`) or static (`GR`) visualizations.
  - Highlights the original points and the computed averages.

- **Master script** (`master.jl`):
  - Example run combining the components:
    - Selects a function (e.g. quarter circle).
    - Chooses two points `a` and `b` and a weight `w`.
    - Computes the three averaging methods.
    - Displays results numerically and graphically.

## 📖 License
This project is licensed under the [Unlicense](https://unlicense.org/),  
which dedicates the work to the public domain.  
You are free to copy, modify, distribute, and use it for any purpose without restriction.
