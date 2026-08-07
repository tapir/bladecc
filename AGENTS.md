# AGENTS.md

Guidance for AI agents (and humans) working on BladeCC.

## What this is

C99 port of Andrew Ning's **CCBlade.jl** (blade element momentum solver for
propellers/fans/turbines) plus the **FLOWMath.jl** routines it depends on.
Upstream: https://github.com/byuflowlab/CCBlade.jl and
https://github.com/byuflowlab/FLOWMath.jl (both MIT, see LICENSE.md — retain
the attribution).

## Commands

```sh
make              # build libbladecc.a
make test         # build + run test suite (1468 checks, must all pass)
make check        # cppcheck static analysis (must be clean)
make analyze      # scan-build static analysis (must be clean)
make format       # clang-format -i (Google style, .clang-format)
make format-check # formatting gate for CI
```

No external dependencies beyond libc/libm. **Do not add library
dependencies** without confirming first.

## Prime directive: bit-fidelity to Julia

This port is verified against the Julia implementation to ~1e-11 relative
error. When changing numerical code:

- Port formulas **exactly** as written in the Julia source, including
  operation order, endpoint handling, and edge-case branches. Do not
  "improve" or mathematically rearrange expressions — small reorderings
  change floating-point results and break parity.
- The Akima spline is FLOWMath's variant (custom endpoint slopes,
  clamped-index extrapolation, n=2 zero-initialized extension slopes), not
  textbook Akima.
- The Brent solver mirrors scipy's brentq (atol=2e-12, rtol=4*eps) with its
  exact sequential-assignment swap semantics.
- `isapprox(x, y, atol=e)` in Julia with an explicit atol is a pure
  absolute-tolerance check (rtol defaults to 0) — port as `fabs(x-y) <= e`.

## Regenerating test data

- `tests/test_arrays.h` — expected values extracted from the Julia suite:
  `python3 tests/extract_arrays.py`
- `tests/reference_data.h` — Julia-generated values for viterna, Akima, and
  Re/Mach airfoil evaluators:
  `julia --project=../CCBlade.jl tests/gen_reference.jl` (from `tests/`)
- End-to-end C-vs-Julia sweep: run `tests/compare_julia.jl`, compile/run
  `tests/compare_c.c`, diff the two CSVs (both scripts print instructions;
  remove the CSVs afterwards — they are gitignored).

## Code conventions

- Style: clang-format Google (enforced by `make format-check`).
- All public API is declared in `include/ccblade.h`, prefixed `cc_`.
- Memory: every malloc must be NULL-checked and every allocation freed on
  all paths (both analyzers and LeakSanitizer enforce this). `CCAirfoil`
  owns its arrays; `CCSection.af` is a borrowed pointer.
- Julia's silent `@warn` + zero-loads failure mode is preserved, but
  `cc_solve_ex` exposes a `success` out-param — batch callers (e.g. LUT
  generation) must check it.

## Before committing

1. `make format && make test && make check && make analyze` — all clean.
2. If numerical behavior changed intentionally, regenerate reference data
   and re-run the Julia comparison sweep.
3. Keep ASan/UBSan/LeakSanitizer clean:
   `gcc -fsanitize=address,undefined,leak -Iinclude tests/runtests.c src/*.c -lm`
