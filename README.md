# BladeCC

C port of [CCBlade.jl](https://github.com/byuflowlab/CCBlade.jl) (Andrew
Ning's blade element momentum solver for propellers/fans and wind
turbines).

## Credits and license

This project is a derivative work of two Julia packages by **Andrew Ning**
(BYU FLOW Lab):

- [CCBlade.jl](https://github.com/byuflowlab/CCBlade.jl) — the BEM solver
  ported here in full (`src/ccblade.c`, `src/airfoils.c`)
- [FLOWMath.jl](https://github.com/byuflowlab/FLOWMath.jl) — the numerical
  routines CCBlade depends on, ported as needed (`src/interp.c` and the
  Brent solver in `src/ccblade.c`): Akima splines, recursive 2D/3D
  interpolation, trapezoidal integration, and Brent root finding

Both upstream packages are MIT licensed. This port retains the same MIT
license (see [LICENSE.md](LICENSE.md)) with the original copyright notice.

## Build & test

```sh
make          # builds libbladecc.a
make test     # builds and runs the test suite (from tests/)
make check    # static analysis with cppcheck
make analyze  # static analysis with scan-build (clang)
make format   # clang-format all sources (Google style)
make format-check  # verify formatting without modifying files
```

Dependencies: a C99 compiler and libm. No external libraries.
`make check` needs cppcheck; `make analyze` needs clang's scan-build
`make format`/`format-check` need clang-format

## Layout

- `include/ccblade.h` — public API
- `src/ccblade.c` — BEM core: residual, quadrant-aware bracketing, Brent
  root solve, `simple_op`/`windturbine_op`, thrust/torque integration (incl.
  azimuth-averaged), `nondim`
- `src/airfoils.c` — airfoil data: file parsing/writing, `AlphaAF`/
  `AlphaReAF`/`AlphaMachAF`/`AlphaReMachAF`/`SimpleAF`/function callbacks,
  Prandtl–Glauert / skin-friction / DuSelig–Eggers corrections, Prandtl
  tip/hub loss, Viterna extrapolation, linear lift-curve extraction
- `src/interp.c` — Akima spline, recursive 2D/3D interpolation, trapezoid
  rule
- `tests/runtests.c` — port of the Julia test suite
- `tests/gen_reference.jl`, `tests/reference_data.h` — Julia-generated
  reference values for functions the Julia suite does not cover (viterna,
  Re/Mach airfoil evaluators, Akima extrapolation)
- `tests/extract_arrays.py`, `tests/test_arrays.h` — test data transcribed
  programmatically from `runtests.jl`
- `tests/compare_julia.jl`, `tests/compare_c.c` — end-to-end C-vs-Julia
  sweep comparison
- `tests/airfoils/` — polar data files (copied from the Julia repo)

## Port fidelity

- The Akima spline is a **direct port of FLOWMath's implementation**
  (endpoint slope estimation, clamped-index extrapolation, n=2 behavior)
- All quadrant/branch logic in `solve` (hover, reversed flow, windmill,
  empirical high-induction region) is ported verbatim.
- Verification: the three portable Julia testsets pass (normal operation
  157, qualitative check 338, correction methods — expanded per-element),
  plus ~530 checks against Julia-generated reference data, plus an
  independent 42-point thrust/torque sweep (hover to 55 m/s, 0–12000 RPM)
  agreeing with Julia to a max relative error of ~2e-11. Clean under
  ASan/UBSan/LeakSanitizer.

## Not ported (Julia-specific)

- Automatic differentiation (ImplicitAD / ForwardDiff glue) — the solver
  always evaluates numerically; call `cc_solve` for values only.
- StructArrays broadcasting, OffsetArrays/FillArrays interop — in C you
  loop over `CCSection`/`CCOperatingPoint` arrays yourself.
- The Julia `Rotor` corrections are struct fields instead of dispatched
  types; see `CCRotor` in `ccblade.h`.

## Notes

- `cc_solve` prints a warning and returns zero loads when no bracket is
  found, mirroring Julia's `@warn` behavior; use `cc_solve_ex(...,
  &success)` to detect this programmatically (relevant for batch LUT
  generation).
