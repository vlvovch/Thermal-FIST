# Numerical integration

`OrthogonalPolynomials.h` and `GaussianQuadrature.h` are vendored from [GaussQuad](https://github.com/vlvovch/gaussquad-cpp) v1.1.0 (commit a1dc387).

Local changes:
- namespace `gaussquad` → `thermalfist::NumericalIntegration`, include paths and guards
- `InitialGuessLaguerre`/`InitialGuessHermite` use Numerical Recipes guesses with bracketing (the upstream seeds can give duplicate roots)

To update, copy the upstream headers and reapply these changes.
