# cavity-simple

A finite-volume solver for the steady incompressible Navier–Stokes equations on a
collocated Cartesian grid, applied to the lid-driven cavity benchmark.

Written in C++ as part of a self-directed study programme following Ferziger,
Perić & Street, *Computational Methods for Fluid Dynamics*. The goal of this
repository is not performance but a solver whose every discretisation choice can
be justified from first principles.

## Problem

Square cavity of side `L`. The lid translates at constant velocity `U`; no-slip
is imposed on all four walls. Steady laminar flow at `Re = ρUL/μ`.

## Numerical method

| Component | Choice |
|---|---|
| Discretisation | Finite volume, cell-centred |
| Grid | Uniform Cartesian, **collocated** |
| Convection | Central differences (2nd order) |
| Diffusion | Central differences (2nd order) |
| Pressure–velocity coupling | SIMPLE |
| Checkerboard suppression | Rhie–Chow momentum interpolation |
| Linear solver | Stone's Strongly Implicit Procedure (SIP), α = 0.92 |
| Wall treatment | Diffusive coefficient `2μ` (wall lies at `Δx/2` from the cell centre); zero convective flux |

Three points are worth stating explicitly, since they are where a collocated
SIMPLE implementation is usually either right or wrong.

**Pressure–velocity decoupling.** On a collocated grid, linear interpolation of
the cell velocities to the faces produces a mass flux that does not see the
pressure jump across the face it crosses: the two central gradients cancel and
the odd–even modes become invisible to the continuity equation. The Rhie–Chow
term subtracts the interpolated gradient and reinstates the two-point one, so the
face flux responds to `p_E − p_P` with coefficient `(V/A_P)_e`. This is
deliberately inconsistent at `O(Δx²)`, and it is the same coefficient that
appears in the pressure-correction equation — the two are consistent by
construction, which is what makes the correction drive the mass imbalance to
zero.

**Singular pressure system.** Pressure is defined up to a constant, so the
pressure-correction matrix is singular. Cell 0 is fixed as reference with
`p' = 0`. The compatibility condition is satisfied because all boundary faces
are walls, so the global mass imbalance telescopes to zero.

**Under-relaxation** is applied in implicit form, `A_P ← A_P/α_u` and
`Q ← Q + (1−α_u) A_P φ_old`, not by blending successive solutions.

Face mass fluxes are rebuilt from the Rhie–Chow expression at every outer
iteration and then corrected with `m' = −ρ(V/A_P)_f (p'_E − p'_P)`, so that
discrete continuity is satisfied at convergence.

## Build and run

```bash
g++ -O2 -o cavity cavity_simple.cpp
./cavity
```

No external dependencies. Parameters (`N`, `Re`, `U`, `ρ`, `α_u`, `α_p`, `tol`)
are compile-time constants at the top of the source file and in `main`.

Output:

- `field.csv` — one row per cell: `x, y, u, v, p` at cell centres
- `params.csv` — run parameters, iteration count, final residual

## Convergence criterion

The three residuals (`u`, `v`, mass) are normalised by `ρU²L` and `ρUL`
respectively, and the largest is compared against `tol`. The momentum residuals
are those returned by SIP before the sweep, i.e. how far the current field is
from satisfying the discrete equation.

## Status

Implemented and running:

- [x] SIP linear solver
- [x] Momentum equations with CDS convection and diffusion
- [x] Rhie–Chow interpolation
- [x] SIMPLE outer loop with pressure correction
- [x] CSV output for post-processing

Not yet done:

- [ ] Formal verification and validation (see below)
- [ ] Non-uniform grid
- [ ] Upwind and deferred-correction blending (UDS/CDS)
- [ ] Grid-convergence study with Richardson extrapolation and GCI

## Verification and validation

**Not yet performed.** The solver converges and produces a qualitatively correct
field, but no quantitative comparison is committed to this repository yet, and no
claim of validation should be read into the code as it stands.

Planned, in order:

1. Manufactured solution for the linear solver alone, to confirm that SIP
   converges to the exact solution of the discrete system.
2. Pure diffusion (`Re → 0`) against the analytical solution.
3. Velocity profiles on the centrelines against Ghia, Ghia & Shin (1982) at
   Re = 100, 400, 1000.
4. Streamfunction extrema against Ferziger & Perić §8.4.1.
5. Observed order of accuracy from three systematically refined grids, with GCI.

## Known limitations

These are current properties of the implementation, not open questions.

1. **Rhie–Chow uses the under-relaxed `A_P`.** The damping coefficient therefore
   scales with `α_u`, and the converged solution retains a weak dependence on the
   relaxation factor (Majumdar, 1988). Correcting this requires keeping an
   un-relaxed copy of `A_P` for use in the interpolation.
2. **Pressure at boundary faces uses the cell value as ghost value.** The
   resulting one-sided gradient is effectively taken over `2Δx` rather than `Δx`,
   which is first-order accurate at the walls. Linear extrapolation would restore
   second order. This is expected to degrade the observed order of accuracy in a
   refinement study.
3. **The residual is an unnormalised sum over cells.** The convergence tolerance
   is therefore grid-dependent: the same `tol` imposes a different condition at
   N = 64 and N = 256. This must be fixed before any grid-refinement study, since
   incomplete iterative convergence would contaminate the Richardson
   extrapolation.
4. Central differencing is unbounded; the solver is not expected to remain stable
   at high Reynolds number on coarse grids.

## References

- Ferziger, J. H., Perić, M. & Street, R. L., *Computational Methods for Fluid
  Dynamics*, Springer.
- Ghia, U., Ghia, K. N. & Shin, C. T. (1982), "High-Re solutions for
  incompressible flow using the Navier-Stokes equations and a multigrid method",
  *Journal of Computational Physics*, 48(3), 387–411.
- Rhie, C. M. & Chow, W. L. (1983), "Numerical study of the turbulent flow past
  an airfoil with trailing edge separation", *AIAA Journal*, 21(11), 1525–1532.
- Stone, H. L. (1968), "Iterative solution of implicit approximations of
  multidimensional partial differential equations", *SIAM Journal on Numerical
  Analysis*, 5(3), 530–558.
- Majumdar, S. (1988), "Role of underrelaxation in momentum interpolation for
  calculation of flow with nonstaggered grids", *Numerical Heat Transfer*, 13(1),
  125–132.

## Note on tooling

The solver was written by hand from the equations. An AI assistant was used for
code review and refactoring after the numerical method was working. All
discretisation choices, debugging decisions and convergence settings are the
author's own and are documented in the comments.
